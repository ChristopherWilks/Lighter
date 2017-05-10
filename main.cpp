#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sys/stat.h>
#include <map>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "utils.h"
#include "Store.h"
#include "StoreCQF.h"
#include "ErrorCorrection.h"
#include "Reads.h"
#include "KmerCode.h"
#include "GetKmers.h"
#include "pthread.h"


char LIGHTER_VERSION[] = "Torch v0.1" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

int MAX_CORRECTION ;
bool ALLOW_TRIMMING ;
int SET_NEW_QUAL ;
bool zlibVersionChecked = false ; 

struct _summary
{
	uint64_t corrCnt ;
	uint64_t trimReadsCnt ;
	uint64_t trimBaseCnt ;
	uint64_t discardReadsCnt ;
	uint64_t totalReads ;
	uint64_t errorFreeReadsCnt ;
} ;

struct _OutputThreadArg
{
	struct _summary *summary ;
	struct _Read *readBatch ;
	int batchSize ;
	bool paraDiscard ;

	Reads *reads ;
	int fileInd ;
} ;

void PrintHelp()
{
	printf( "Usage: ./lighter [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
	        "\t             The file can be fasta and fastq, and can be gzip\'ed with extension *.gz.\n"      
		"\t             When the input file is *.gz, the corresponding output file will also be gzip\'ed.\n"
		"\t-k kmer_length genome_size alpha: (see README for information on setting alpha)\n"
		"\t\t\t\t\tor\n"
		"\t-K kmer_length genom_size: in this case, the genome size should be relative accurate.\n"
		"Other parameters:\n"
		"\t-od output_file_directory: (default: ./)\n"
		"\t-t num_of_threads: number of threads to use (default: 1)\n"
		"\t-trim: allow trimming (default: false)\n"
		"\t-noQual: ignore the quality socre (default: false)\n"
		"\t-newQual ascii_quality_score: set the quality for the bases corrected to the specified score (default: not used)\n"
		"\t-h: print the help message and quit\n"
		"\t-v: print the version information and quit\n") ;
}

uint64_t StringToUint64( char *s )   
{
	int i ;
	uint64_t ret = 0 ;
	for ( i = 0 ; s[i] ; ++i )
	{
		ret = ret * 10 + s[i] - '0' ;
	}
	return ret ;
}

inline void ExtractKmer( char *s, int offset, int kmerLength, char *buffer )
{
	int i ;
	for ( i = offset ; i < offset+ kmerLength ; ++i )
		buffer[i - offset] = s[i] ;
	buffer[i - offset] = '\0' ;
}

void GetCumulativeBinomialDistribution( double F[], int l, double p )
{
	// p is the probability of getting 1.
	int i ;
	double coef = 1 ;
	double exp = pow( 1 - p, l ) ;
	F[0] = pow( 1 - p, l ) ;

	for ( i = 1 ; i <= l ; ++i )
	{
		coef = coef / i * ( l - i + 1 ) ;
		exp = exp / ( 1 - p ) * p ;
		F[i] = F[i - 1] + coef * exp ;		
	}
}

char GetGoodQuality( Reads &reads )
{
	int i ;
	int qualHisto[300] ;
	int totalCnt, cnt ;

	//Reads reads( readFile ) ;
	if ( !reads.HasQuality() )
		return 127 ;

	memset( qualHisto, 0, sizeof( qualHisto ) ) ;
	for ( i = 0 ; i < 1000000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		++qualHisto[ (int)reads.qual[0] ] ;
	}

	totalCnt = i ;
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += qualHisto[i] ;
		if ( cnt > totalCnt * 0.25 )
			break ;
	}
	return (char)i ;
}

char GetBadQuality( Reads &reads )
{
	int i ;
	int qualHisto[300], firstQualHisto[300] ;
	int totalCnt, cnt ;
	int t1, t2 ;
	//Reads reads( readFile ) ;
	if ( !reads.HasQuality() )
		return 0 ;

	memset( qualHisto, 0, sizeof( qualHisto ) ) ;
	memset( firstQualHisto, 0, sizeof( firstQualHisto )) ;
	for ( i = 0 ; i < 1000000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		++qualHisto[ (int)reads.qual[ strlen( reads.seq ) - 1 ] ] ;
		++firstQualHisto[ (int)reads.qual[0] ] ;
	}

	totalCnt = i ;
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += firstQualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t1 = i - 1 ;
	
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += qualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t2 = i ;
	return (char)( t2 < t1 ? t2 : t1 ) ;
}

double InferAlpha( Reads &reads, uint64_t genomeSize ) 
{
	uint64_t totalLen = 0 ;

	while ( reads.Next() )
		totalLen += strlen( reads.seq ) ;		
	
	return 7.0 / ( (double)totalLen / genomeSize ) ;
}

void PrintLog( const char *log ) 
{
	time_t rawtime ;
	struct tm *timeInfo ;
	char buffer[128] ;
	//FILE *fp = fopen( "lighter.log", "a" ) ;
	
	time( &rawtime ) ;
	timeInfo = localtime( &rawtime ) ;
	strftime( buffer, sizeof( buffer ), "%F %H:%M:%S", timeInfo ) ;

	fprintf( stderr, "[%s] %s\n", buffer, log ) ;

	//fclose( fp ) ;
}

void UpdateSummary( char *seq, int correction, int badSuffix, bool paraDiscard, struct _summary &summary ) 
{
	if ( correction == 0 )
		++summary.errorFreeReadsCnt ;			
	else if ( correction > 0 )
		summary.corrCnt += correction ;	
	else if ( paraDiscard ) // tmp < 0
		++summary.discardReadsCnt ;

	if ( ALLOW_TRIMMING && badSuffix > 0 )
	{
		++summary.trimReadsCnt ;
		summary.trimBaseCnt += badSuffix ;
	}
	++summary.totalReads ;	
}

void PrintSummary( const struct _summary &summary )
{
	fprintf( stderr, "Processed %" PRIu64 " reads:\n"
		"\t%" PRIu64 " are error-free\n"
		"\tCorrected %" PRIu64 " bases(%lf corrections for reads with errors)\n"
		"\tTrimmed %" PRIu64 " reads with average trimmed bases %lf\n"
		"\tDiscard %" PRIu64 " reads\n",
		summary.totalReads, summary.errorFreeReadsCnt,
		summary.corrCnt, 
		summary.totalReads == summary.errorFreeReadsCnt ? 0.0 : 
					(double)summary.corrCnt / ( summary.totalReads - summary.errorFreeReadsCnt ),
		summary.trimReadsCnt, 
		summary.trimReadsCnt == 0 ? 0.0 : (double)summary.trimBaseCnt / summary.trimReadsCnt, 
		summary.discardReadsCnt ) ;	
}


void *Output_Thread( void *arg )
{
	int i ;
	struct _OutputThreadArg *myArg = ( struct _OutputThreadArg *)arg ;
	struct _Read *readBatch = myArg->readBatch ;
	int batchSize = myArg->batchSize ;
	
	//fprintf( stderr, "hi\n"  ) ;
	for ( i = 0 ; i < batchSize ; ++i )
		UpdateSummary( readBatch[i].seq, readBatch[i].correction, readBatch[i].badSuffix, myArg->paraDiscard, *( myArg->summary ) ) ;			
	myArg->reads->OutputBatch( readBatch, batchSize, ALLOW_TRIMMING, myArg->fileInd ) ;
	pthread_exit( NULL ) ;
	return NULL ;
}

int test_kmer_counts(char* reads_file, char*** tkmers, int** counts)
{
	//if reads_file == 1M
	if(strcmp(reads_file,"../tests/lambda_reads_1.fq") == 0)
	{
		/*char* tkmers_[2]={"AACCACCAGGCCATATCTGCC","ATGGAATTAAGTCGCACACCC"};
		*tkmers = tkmers_;
		int counts_[2]={1,18};
		*counts = counts_;*/
		(*tkmers)[0]="AACCACCAGGCCATATCTGCC";
		(*tkmers)[1]="ATGGAATTAAGTCGCACACCC";
		(*counts)[0]=1;
		(*counts)[1]=18;
		return 2;
	}
	else if(strcmp(reads_file,"../tests/SRR197986_1.fastq.1m") == 0)
	{
		/*char* tkmers_[5]={"AACAGTGGCCCTTAATCAAAG","ACTGCAGGCAACAAACACAAA","ATGGGGGATTCGCGAAGAGAA","AATGGGGGATTTGCAAAGAGA","AAATCACGCGTTTTCTCTTCG"};
		*tkmers = tkmers_;
		int counts_[5]={1,11,245,2183,20391};
		*counts = counts_;*/
		(*tkmers)[0]="AACAGTGGCCCTTAATCAAAG";
		(*tkmers)[1]="ACTGCAGGCAACAAACACAAA";
		(*tkmers)[2]="ATGGGGGATTCGCGAAGAGAA";
		(*tkmers)[3]="AATGGGGGATTTGCAAAGAGA";
		(*tkmers)[4]="AAATCACGCGTTTTCTCTTCG";
		(*counts)[0]=1;
		(*counts)[1]=11;
		(*counts)[2]=245;
		(*counts)[3]=2183;
		(*counts)[4]=20391;
		return 5;
	}
	else if(strcmp(reads_file,"../tests/SRR197986_1.fastq.2m") == 0)
	{
		(*tkmers)[0]="AACAGTGGCCCTTAATCAAAG";
		(*tkmers)[1]="ACTGCAGGCAACAAACACAAA";
		(*tkmers)[2]="CCCCCATTTGACCCGAAAATC";
		(*tkmers)[3]="AATGGGGGATTTGCAAAGAGA";
		(*tkmers)[4]="CGAGATCGGTCTCGGCATTCC";
		(*tkmers)[5]="GCGGTTCAGCAGGAATGCCGA";
		(*counts)[0]=1;
		(*counts)[1]=18;
		(*counts)[2]=143;
		(*counts)[3]=4510;
		(*counts)[4]=17248;
		(*counts)[5]=119528;
		return 6;
	}
	else if(strcmp(reads_file,"../tests/SRR197986_1.fastq.4m") == 0)
	{
		(*tkmers)[0]="CGGGCCGTTGCACGCAGGTCC";
		(*tkmers)[1]="ACTGCAGGCAACAAACACAAA";
		(*tkmers)[2]="CATCCAGGGATGGTGACTCAA";
		(*tkmers)[3]="AACATCAGGCATTTTCTCTTA";
		(*tkmers)[4]="AAACGCGTGATTTTCACTTAA";
		(*tkmers)[5]="GATCGGAAGAGCGGTTCAGCA";
		(*counts)[0]=1;
		(*counts)[1]=33;
		(*counts)[2]=242;
		(*counts)[3]=1445;
		(*counts)[4]=10164;
		(*counts)[5]=307022;
		return 6;
	}
	return 0;
}



void check_known_kmers_counts(Store* kmers, StoreCQF* kmerCounters, char* reads_file, int kmerLength)
{
	//char* tkmers[2]={"AACCACCAGGCCATATCTGCC","ATGGAATTAAGTCGCACACCC"};
	//int counts[2]={1,18};
	char** tkmers = new char*[6];
	int* counts = new int[6];
	KmerCode kmerCode( kmerLength );
	int num_kmers = test_kmer_counts(reads_file, &tkmers, &counts);
	printf("%s %d %s %d\n",tkmers[0],counts[0],tkmers[1],counts[1]);

	int i;
	for(i=0; i < num_kmers; i++)
	{
		kmerCode.Restart();
		int j;
		//printf("%s %d %s %d\n",tkmers[0],counts[0],tkmers[1],counts[1]);
		//printf("%s %d\n",tkmers[i],counts[i]);
		for(j=0; j<kmerLength; j++)
		{
			kmerCode.Append( tkmers[i][j] );
		}
		uint64_t count = kmerCounters->IsIn(kmerCode);
		double diff = (double)count/counts[i];
		if(count == 0)
			diff = (double) (counts[i] - count);
		int diff_log = log10(diff);
		fprintf(stderr,"%d\t%d\t%.3f\t%d\t%s\n",counts[i],count,diff,diff_log,tkmers[i]);
		if(abs(diff_log) >= 1 && abs(count - counts[i]) > 1)
			fprintf(stderr,"ORM\t%d\t%d\t%.3f\t%d\t%s\n",counts[i],count,diff,diff_log,tkmers[i]);
	}
}

int main( int argc, char *argv[] )
{
	int kmerLength ;
	double alpha = -1 ;
	//char *readId/**, read, *qual*/ ;
	char buffer[1023] ;
	double untrustF[100][100] ;
	//double trustF[100][100] ;
	int threshold[MAX_KMER_LENGTH+1] ;
	char goodQuality = '\0', badQuality = '\0' ;
	int badPrefix, badSuffix ;
	bool paraDiscard ;
	bool ignoreQuality, inferAlpha ; //stable ;
	int zlibLevel ;
	//double bloomFilterFP = 0.0005 ;
	int i, j ;
	//uint64_t kmerCode ;
	//uint64_t mask ;
	uint64_t genomeSize = 0;
	struct _summary summary ;

	struct _SamplePattern *samplePatterns = NULL ;

	// variables for threads
	int numOfThreads ;
	pthread_attr_t pthreadAttr ;
	pthread_t *threads = NULL;
	pthread_mutex_t mutexSampleKmers, mutexStoreKmers ;

	if ( argc == 1 )
	{
		PrintHelp() ;
		exit( EXIT_FAILURE ) ;
	}

	Reads reads ;
	/*reads.AddReadFile( argv[1] ) ;
	kmerLength = atoi( argv[2] ) ;
	genomeSize = StringToUint64( argv[3] ) ;
	alpha = (double)atof( argv[4] ) ;*/

	ALLOW_TRIMMING = false ;
	SET_NEW_QUAL = -1 ;
	kmerLength = -1 ;
	numOfThreads = 1 ;
	ignoreQuality = true ; //false
	//stable = false ;
	inferAlpha = false ;
	zlibLevel = 1 ;
	paraDiscard = false;
	memset( &summary, 0, sizeof( summary ) ) ;
	
	// Parse the arguments
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			//reads.AddReadFile( argv[i + 1 ] ) ;
			++i;
			continue ; // wait to be processed after next round
		}
		else if ( !strcmp( "-k", argv[i] ) )
		{
			if(i + 1 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			kmerLength = atoi( argv[i + 1] ) ;
			if(i + 2 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			genomeSize = StringToUint64( argv[i + 2] ) ;
			if(i + 3 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			alpha = (double)atof( argv[i + 3] ) ;
			i += 3 ;
		}
		else if ( !strcmp( "-K", argv[i] ) ) 
		{
			if(i + 1 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size after -K\n");
				exit( EXIT_FAILURE );
			}
			kmerLength = atoi( argv[i + 1] ) ;
			if(i + 2 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size after -K\n");
				exit( EXIT_FAILURE );
			}
			genomeSize = StringToUint64( argv[i + 2] ) ;

			inferAlpha = true ;
			i += 2 ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			mkdir( argv[i + 1], 0700 ) ;
			reads.SetOutputDirectory( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-trim", argv[i] ) )
		{
			ALLOW_TRIMMING = true ;
		}
		else if ( !strcmp( "-t", argv[i] ) )
		{
			numOfThreads = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-noQual", argv[i] ) )
		{
			ignoreQuality = true ;
		}
		else if ( !strcmp( "-newQual", argv[i] ) )
		{
			SET_NEW_QUAL = (int)argv[i + 1][0] ;
			++i ;
		}
		else if ( !strcmp( "-h", argv[i] ) )
		{
			PrintHelp() ;
			exit( 0 ) ;
		}
		else if ( !strcmp( "-v", argv[i] ) )
		{
			printf( "%s\n", LIGHTER_VERSION ) ;
			exit( 0 ) ;
		}
		else
		{
			fprintf( stderr, "Unknown argument %s\n", argv[i] ) ;
			exit( EXIT_FAILURE ) ;
		}
	}

	// Go the second round to get the reads files
	int last_reads_file_idx = -1;
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1 ] ) ;
			++i;
			last_reads_file_idx = i;
		}
	}
	if ( kmerLength == -1 )
	{
		fprintf( stderr, "Require -k or -K parameter!\n" ) ;
		exit( EXIT_FAILURE ) ;
	}
	if ( kmerLength > MAX_KMER_LENGTH )
	{
		fprintf( stderr, "K-mer length must be no larger than %d. You can adjust the MAX_KMER_LENGTH constraints in utils.h.\n", MAX_KMER_LENGTH ) ;
		exit( EXIT_FAILURE ) ;
	}
	
	if ( alpha != -1 && inferAlpha == true )
	{
		fprintf( stderr, "Can not use both -k and -K.\n" ) ;
		exit( EXIT_FAILURE ) ;
	}

	PrintLog( "=============Start====================" ) ;
	KmerCode kmerCode( kmerLength ) ;
	reads.SetDiscard( paraDiscard ) ;	

	if ( inferAlpha )
	{
		PrintLog( "Scanning the input files to infer alpha(sampling rate)" ) ;
		alpha = InferAlpha( reads, genomeSize ) ;
		
		sprintf( buffer, "Average coverage is %.3lf and alpha is %.3lf", 7.0 / alpha, alpha ) ;
		PrintLog( buffer ) ;
		
		reads.Rewind() ;
	}

	// Prepare data structures and other data.
	//Store kmers(1000000000ull) ;
	//Store trustedKmers(1000000000ull) ;
	
	Store kmers((uint64_t)( genomeSize * 1.5 ), 0.01 ) ;
	StoreCQF kmerCounters;


	if ( ignoreQuality == false )
		badQuality = GetBadQuality( reads ) ;
	if ( badQuality != '\0' )
	{
		sprintf( buffer, "Bad quality threshold is \"%c\"", badQuality ) ;
		PrintLog( buffer ) ;
	}
	else
	{
		PrintLog( "No quality score used." ) ;
	}
	reads.Rewind() ;

	// Step 1: Sample the kmers 
	//printf( "Begin step1. \n" ) ; fflush( stdout ) ;
	srand( 17 ) ;
	size_t nuniq_kmer_count_seen = 0;
	size_t nuniq_kmer_count_added = 0;
	// It seems serialization is faster than parallel. NOT true now!
	//std::map<uint64_t, int> hash ;
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() != 0 )
		{
			SampleKmersInRead( reads.seq, reads.qual, kmerLength, alpha, kmerCode, &kmers, &kmerCounters, &nuniq_kmer_count_added, &nuniq_kmer_count_seen) ;
		}
	}

	// Update the bloom filter's false positive rate.
	// Compute the distribution of the # of sampled kmers from untrusted and trusted position
	double tableAFP = kmers.GetFP() ;
	//how many have > 1 occurrence
	uint64_t total_count = 0;
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		/*while ( reads.Next() != 0 )
		{
			total_count += CountKmers( reads.seq, reads.qual, kmerLength, kmerCode, &kmerCounters, 0) ;
		}*/
		total_count = CountKmers( reads.seq, reads.qual, kmerLength, kmerCode, &kmerCounters, 0) ;
	}
	PrintLog( "Finish counting kmers" ) ;

	//now try querying the static set of kmers per the input file
	check_known_kmers_counts(&kmers, &kmerCounters, argv[last_reads_file_idx], kmerLength);
		
	sprintf( buffer, "Bloom filter A's false positive rate: %lf\nNon-unique K-mers seen %lu\nNon-unique K-mers added %lu\ntotal # of kmers with > 1 occurrences: %lu\n", tableAFP, nuniq_kmer_count_seen, nuniq_kmer_count_added, total_count ) ;
	PrintLog( buffer ) ;

	PrintLog( "Finish counting kmers" ) ;
	//CW: start the 2nd step of using another DS here (for counting) 

	return 0 ;
}
	
	
