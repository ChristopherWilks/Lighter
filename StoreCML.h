#ifndef _MOURISL_CLASS_STORECML
#define _MOURISL_CLASS_STORECML
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "KmerCode.h"
#include "cmlsketch.h"
#include "bench_common.h"
  		

const static int SF_LENGTH=8;


class StoreCML
{
private:
	uint64_t size ;
	uint64_t uniq ;
	uint64_t uniq2 ;
	//void *cml;
	CMLSketch* cml;
	//bloom_filter bf ;
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map

	//following functions from:
	//http://stackoverflow.com/questions/9695720/how-do-i-convert-a-64bit-integer-to-a-char-array-and-back
	//needed it to convert from Lighter's uint64 to character array for SF sketch
	void int64ToChar(unsigned char a[], int64_t n) {
		  memcpy(a, &n, 8);
	}

	int64_t charTo64bitNum(unsigned char a[]) {
		  int64_t n = 0;
		  memcpy(&n, a, 8);
		  return n;
	}
	
	void increase(uint64_t val, int kmerLength, size_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
          	cml->insert(valc, SF_LENGTH, count);
	}

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
          	cml->insert(valc, SF_LENGTH, 1);
		return 1;
	}
	
	size_t IsIn( uint64_t val, int kmerLength ) 
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
       		return cml->queryPoint(valc, SF_LENGTH);
	}
	
	size_t IsInExact( uint64_t val, int kmerLength ) 
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
		if(hash.find( val ) != hash.end())
		{
			//remove so we don't count twice
			return 0;
		}
		hash[ val ] = 1;
       		return cml->queryPoint(valc, SF_LENGTH);
	}

	int numOfThreads ;
public:
	StoreCML( CMLSketch* cml_ )
	{
		numOfThreads = 1 ;
		cml = cml_;
		method = 0 ;
		uniq = 0;
		uniq2 = 0;
	}

	~StoreCML() 
	{
	}
	
	double Occupancy()
	{
		return 0;
	}
	double GetFP()
	{
		return 0;
	}
	
	int increase(KmerCode &code, size_t count) 
	{
		if ( !code.IsValid() )
			return 0 ;
		increase(code.GetCode(), code.GetKmerLength(), count) ;
	}

	int Put( KmerCode &code, bool testFirst = false ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		Put( code.GetCode(), code.GetKmerLength(), testFirst ) ;
		return 0 ;
	}
	
	int IsIn( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		return IsIn( code.GetCode(), code.GetKmerLength() ) ;
	}
	
	int IsInExact( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		return IsInExact( code.GetCode(), code.GetKmerLength() ) ;
	}

	uint64_t uniq_kmers() { return uniq; }
	uint64_t uniq_kmers_gt_2() { return uniq2; }


	void SetNumOfThreads( int in ) 
	{ 
		numOfThreads = in ;
		//bf.SetNumOfThreads( in ) ;
	}

#if MAX_KMER_LENGTH <= 32 
	uint64_t GetCanonicalKmerCode( uint64_t code, int k ) 
	{
		int i ;
		uint64_t crCode = 0ull ; // complementary code
		for ( i = 0 ; i < k ; ++i )
		{
			//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
			//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;

			uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
			crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
		}
		return crCode < code ? crCode : code ;
	}
#else
	void GetCanonicalKmerCode( uint64_t code[], int k ) 
	{
		int i, j ;
		uint64_t crCode[ K_BLOCK_NUM ] ; // complementary code
		const int largestBlock = ( k - 1 ) / ( sizeof( uint64_t ) * 4 ) ;
		for ( i = 0 ; i < K_BLOCK_NUM ; ++i )
			crCode[i] = 0 ;

		for ( i = 0, j = k - 1 ; i < k ; ++i, --j )
		{
			//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
			//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;
			
			/*int tagI = i >> 5 ;
			int tagJ = j >> 5 ;
			int residualI = i & 31 ;
			int residualJ = j & 31 ;*/
			uint64_t tmp = ( code[ i >> 5 ] >> ( 2ull * ( i & 31 ) ) ) & 3ull ;
			crCode[j >> 5] = crCode[ j >> 5 ] | ( ( 3ull - tmp ) << ( 2ull * ( j & 31 ) ) ) ;
		}

		bool crFlag = false ;
		for ( i = largestBlock ; i >= 0 ; --i )
		{
			if ( crCode[i] < code[i] )
			{
				crFlag = true ;
				break ;
			}
			else if ( crCode[i] > code[i] )
			{
				crFlag = false ;
				break ;
			}
		}
		//printf( "%llu,%llu %llu,%llu: %d\n", code[1], code[0], crCode[1], crCode[0], crFlag ) ;
		if ( crFlag )
		{
			for ( i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = crCode[i] ;
		}
	}
#endif

	void TemporaryOutput( char *file)
	{
		FILE *fp = fopen( file, "w" ) ;
		//bf.Output( fp ) ;
		fclose( fp ) ;
	}

	void TemporaryInput( char *file ) 
	{
		FILE *fp = fopen( file, "r" ) ;
		//bf.Input( fp ) ;
		fclose( fp ) ;
	}

	int Clear() 
	{
		return 0 ;
	}
} ;
#endif
