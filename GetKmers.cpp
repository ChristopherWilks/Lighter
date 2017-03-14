#include "GetKmers.h"

struct _SampleKmersPutThreadArg
{
	Store *kmers ;
	KmerCode *kmerCodes ;
	int kmerCodesCnt ;
} ;

//CW note: actually stores the kmers (via kmerCode) into the bloom filter
void *SampleKmers_PutThread( void *arg )
{
	struct _SampleKmersPutThreadArg *myArg = ( struct _SampleKmersPutThreadArg * )arg ;
	Store *kmers = myArg->kmers ;
	int i ;
	for ( i = 0 ; i < myArg->kmerCodesCnt ; ++i )
		kmers->Put( myArg->kmerCodes[i], true ) ;

	pthread_exit( NULL ) ;
	return NULL ;
}

void *SampleKmers_Thread( void *arg )
{
	struct _SampleKmersThreadArg *myArg = ( struct _SampleKmersThreadArg *)arg ; 	
	int i, tmp ;
	int fileInd ;
	struct _Read *readBatch = ( struct _Read *)malloc( sizeof( *readBatch ) * 128 ) ;
	struct _SamplePattern *samplePatterns = myArg->samplePatterns ;
	int kmerLength = myArg->kmerLength ;
	KmerCode kmerCode( kmerLength ) ;

	Store *kmers = myArg->kmers ;
	const int bufferSizeFactor = 9 ;
	KmerCode *kmerCodeBuffer[2] ;
	int bufferTag ;
	int kmerCodeBufferUsed = 0 ;
	
	for ( bufferTag = 0 ; bufferTag < 2 ; ++bufferTag )
	{
		kmerCodeBuffer[ bufferTag ] = ( KmerCode * )malloc( sizeof( KmerCode ) * ( bufferSizeFactor + 1 ) * MAX_READ_LENGTH ) ;
		for ( i = 0 ; i < ( bufferSizeFactor + 1 ) * MAX_READ_LENGTH ; ++i )
			kmerCodeBuffer[bufferTag][i] = kmerCode ;
	}

	void *pthreadStatus ;
	pthread_t putThread ;
	pthread_attr_t pthreadAttr ;
	struct _SampleKmersPutThreadArg putArg ;
	bool threadInit = false ;
	pthread_attr_init( &pthreadAttr ) ;
	pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;
	bufferTag = 0 ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		tmp = myArg->reads->GetBatch( readBatch, 128, fileInd, false, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;
		int k ;
		if ( tmp != 0 )
		{
			for ( k = 0 ; k < tmp ; ++k )
			{
				int len ; 
				char *id = readBatch[k].id ;
				char *read = readBatch[k].seq ;
				char *qual = readBatch[k].qual ;
				int tag = fileInd ;

				len = (int)strlen( read ) ;
				if ( read[len - 1] == '\n' )
					read[len - 1] = '\0' ;

				if ( qual[0] != '\0' )
				{
					if ( qual[len - 1] == '\n' )
						qual[len - 1] = '\0' ;
				}

				for ( i = 0 ; id[i] ; ++i )	
					tag = tag * 17  + ( id[i] - 'A' ) ;

				for ( i = 0 ; read[i] ; ++i )
					tag = tag * 7 + ( read[i] - 'A' )  ;
				
				for ( i = 0 ; qual[i] ; ++i )
					tag = tag * 17+ ( qual[i] - 'A' )  ;
					
				tag %= SAMPLE_PATTERN_COUNT ;
				if ( tag < 0 )
					tag += SAMPLE_PATTERN_COUNT ;

				if ( len - 1 < kmerLength )
					continue ;
				//CW note: reads the first kmer, then the rest, stores them as kmerCodes
				for ( i = 0 ; i < kmerLength ; ++i )
				{
					kmerCode.Append( read[i] ) ;
				}
				//printf( "%d %d %d\n", tag, (int)samplePatterns[tag].tag[2], (int)samplePatterns[tag].tag[3] ) ;	
				if ( samplePatterns[tag].tag[ ( i - 1 ) / 8 ] & ( 1 << ( ( i - 1 ) % 8 ) ) )
				{
					//kmers->Put( kmerCode ) ;
					kmerCodeBuffer[ bufferTag ][ kmerCodeBufferUsed ] = kmerCode ;
					++kmerCodeBufferUsed ;
				}

				for ( ; read[i] ; ++i )
				{
					kmerCode.Append( read[i] ) ;

					if ( samplePatterns[tag].tag[ i / 8 ] & ( 1 << ( i % 8 ) ) )
					{
						//kmers->Put( kmerCode ) ;
						//CW note: copy the current kmerCode into the current buffer position
						kmerCodeBuffer[bufferTag][ kmerCodeBufferUsed ] = kmerCode ;
						++kmerCodeBufferUsed ;
					}
				}
				
				if ( kmerCodeBufferUsed >= bufferSizeFactor * MAX_READ_LENGTH )
				{
					if ( threadInit )
						pthread_join( putThread, &pthreadStatus ) ;
					
					putArg.kmers = kmers ;
					putArg.kmerCodes = ( KmerCode *)kmerCodeBuffer[ bufferTag ] ;
					putArg.kmerCodesCnt = kmerCodeBufferUsed ;
					pthread_create( &putThread, &pthreadAttr, SampleKmers_PutThread, ( void *)&putArg ) ;

					kmerCodeBufferUsed = 0 ;
					bufferTag = 1 - bufferTag ;

					threadInit = true ;
				}
			}
		}
		else
			break ;
	}
	
	//put the remaining
	if ( threadInit )
		pthread_join( putThread, &pthreadStatus ) ;

	putArg.kmers = kmers ;
	putArg.kmerCodes = ( KmerCode *)kmerCodeBuffer[ bufferTag ] ;
	putArg.kmerCodesCnt = kmerCodeBufferUsed ;

	pthread_create( &putThread, &pthreadAttr, SampleKmers_PutThread, ( void *)&putArg ) ;
	kmerCodeBufferUsed = 0 ;
	bufferTag = 1 - bufferTag ;
	threadInit = true ;

	pthread_join( putThread, &pthreadStatus ) ;

	free( readBatch ) ;
	for ( i = 0 ; i < 2 ; ++i )
		free( kmerCodeBuffer[i] ) ;
	pthread_attr_destroy( &pthreadAttr ) ;
	pthread_exit( NULL ) ;
	return NULL ;
}

void SampleKmersInRead( char *read, char *qual, int kmerLength, double alpha, KmerCode &kmerCode, Store *kmers, size_t* kcount_added, size_t* kcount_seen )
{
	int i ;
	double p ;
	double factor = 1 ;
	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength && read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}

	if ( i < kmerLength )
		return ;

	*kcount_seen+=1;
	p = rand() / (double)RAND_MAX ;
	if ( p < alpha * factor )
	{
		*kcount_added+=1;
		kmers->Put( kmerCode ) ;
	}

	for ( ; read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;
		*kcount_seen+=1;

		p = rand() / (double)RAND_MAX ;
		if ( p < alpha * factor )
		{
			*kcount_added+=1;
			kmers->Put( kmerCode ) ;
		}
	}
}
void *StoreKmers_Thread( void *arg )
{
	struct _StoreKmersThreadArg *myArg = ( struct _StoreKmersThreadArg *)arg ; 	
	int i ;
	int tmp, fileInd ;
	int maxBatchSize = READ_BUFFER_PER_THREAD ; 
	struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;

	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		tmp = myArg->reads->GetBatch( readBatch, maxBatchSize, fileInd, false, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;

		if ( tmp == 0 )
			break ;
		
		for ( i = 0 ; i < tmp ; ++i )
		{
			char *read = readBatch[i].seq ;
			char *qual = readBatch[i].qual ;

			int len = (int)strlen( read ) ;
			if ( read[len - 1] == '\n' )
				read[len - 1] = '\0' ;

			if ( qual[0] != '\0' )
			{
				if ( qual[len - 1] == '\n' )
					qual[len - 1] = '\0' ;
			}

			StoreTrustedKmers( read, qual, myArg->kmerLength, myArg->badQuality, myArg->threshold, 
				kmerCode, myArg->kmers, myArg->trustedKmers ) ;
		}
	}

	free( readBatch ) ;

	pthread_exit( NULL ) ;
	return NULL ;
}

uint64_t CountKmers( char *read, char *qual, int kmerLength,
	KmerCode &kmerCode, Store *kmers, int cutoff )
{
	bool occur[MAX_READ_LENGTH] ;
	bool trustedPosition[MAX_READ_LENGTH] ;
	int i ;
	uint64_t total_count = 0;

	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	//CW note: do first kmer
	int count = kmers->IsIn( kmerCode);
	if(count > cutoff) {
		total_count+=1;
		kmers->decrease(kmerCode, count); 
		//trustedKmers->Put( kmerCode, true ) ;
	}

	//CW note: then the rest
	for ( ; read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;
		count = kmers->IsIn( kmerCode);
		if(count > cutoff) {
			total_count+=1;
			kmers->decrease(kmerCode, count); 
			//trustedKmers->Put( kmerCode, true ) ;
		}
	}
	return total_count;
}


void StoreTrustedKmers( char *read, char *qual, int kmerLength, char badQuality, int *threshold,
	KmerCode &kmerCode, Store *kmers, Store *trustedKmers )
{
	bool occur[MAX_READ_LENGTH] ;
	bool trustedPosition[MAX_READ_LENGTH] ;
	int i ;

	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	//CW note: do first kmer
	if ( kmers->IsIn( kmerCode) )
		occur[i - kmerLength] = true ;
	else
		occur[i - kmerLength] = false ;

	//CW note: then the rest
	for ( ; read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;

		if ( kmers->IsIn( kmerCode ) )
			occur[i - kmerLength + 1] = true ;
		else
		{
			occur[i - kmerLength + 1] = false ;
		}
	}

	int readLength = i ;
	int occurCnt = readLength - kmerLength + 1 ;
	int zeroCnt = 0, oneCnt = 0 ;


	// Set the trusted positions
	for ( i = 0 ; i < readLength ; ++i )
	{
		//CW note: only do one kmer at a time so reduce count
		//once we've moved beyond the kmer length horizon
		if ( i >= kmerLength )
		{
			if ( occur[i - kmerLength] )
				--oneCnt ;
			else
				--zeroCnt ;
		}

		if ( i < occurCnt )
		{
			if ( occur[i] )
				++oneCnt ;
			else
				++zeroCnt ;
		}
		int sum = oneCnt + zeroCnt ;	
		int adjust = 0 ;
	

		if ( qual[0] != '\0' && qual[i] <= badQuality )
		{
			trustedPosition[i] = false ;
			continue ;
		}

		if ( oneCnt > threshold[sum] + adjust )
		{
			trustedPosition[i] = true ;
		}
		else
		{
			trustedPosition[i] = false ;
		}
	}

	oneCnt = 0 ;
	kmerCode.Restart() ;


	for ( i = 0 ; i < readLength ; ++i )
	{
		if ( trustedPosition[i] )
			++oneCnt ;
		if ( i >= kmerLength )
		{
			if ( trustedPosition[i - kmerLength] ) 
				--oneCnt ;
		}

		kmerCode.Append( read[i] ) ;
		if ( oneCnt == kmerLength )
		   {
			   trustedKmers->Put( kmerCode, true ) ;
		   }
	}
}
