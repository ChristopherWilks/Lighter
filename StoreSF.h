#ifndef _MOURISL_CLASS_STORESF
#define _MOURISL_CLASS_STORESF
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "bloom_filter.hpp"
extern "C" {
#include "sketch.h"
}
#include "KmerCode.h"

/*extern "C" {
bool init(size_t D, size_t WL, size_t Z, size_t bits_c);
void inc(const unsigned char * str, size_t len, size_t delta);     //increase by one
size_t query(const unsigned char * str, size_t len);
}*/

const static int SF_LENGTH=8;


class StoreSF
{
private:
	uint64_t size ;
	bloom_parameters bfpara ;
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

#if MAX_KMER_LENGTH <= 32 
	void increase(uint64_t val, int kmerLength, size_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
		inc(valc, SF_LENGTH, count);
	}

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		//return 0 ;
		//printf( "%d\n", method ) ;
		//printf( "%llu\n", val ) ;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		//printf( "%llu\n", val ) ;

		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);

		//if ( numOfThreads > 1 && testFirst && bf.contains( val ) )
		if ( testFirst && query( valc, SF_LENGTH ) > 0 )
			return 0 ;
		//bf.insert( val ) ;
		inc(valc, SF_LENGTH, 1);
		//int c = query(valc, SF_LENGTH);
		//printf("count of %d\n",c);
		return 1 ;
	}
	
	int IsIn( uint64_t val, int kmerLength ) 
	{
		//printf( "1. %llu\n", val ) ;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		//printf( "2. %llu\n", val ) ;
		
		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);

		//return bf.contains( val ) ;
		return query( valc, SF_LENGTH ) ;
	}

	void decrease(uint64_t val, int kmerLength, uint64_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;

		unsigned char valc[SF_LENGTH];	
		int64ToChar(valc, val);
		
		dec(valc, SF_LENGTH, count);
	}
#else
	int Put( uint64_t code[], int kmerLength, bool testFirst )
	{
		//return 0 ;
		//printf( "%d\n", method ) ;
		//printf( "%llu\n", val ) ;
		GetCanonicalKmerCode( code, kmerLength ) ;
		//printf( "%llu\n", val ) ;
		//exit(1) ;
		if ( method == 1 )
		{
			//hash[ val ] = 1 ;
			return 1 ;
		}
		
		//printf( "1: %d\n", bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )  ;
		if ( numOfThreads > 1 && testFirst && bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )
			return 0 ;
		//printf( "2: %lld %d\n", code[0], sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) )  ;
		bf.insert( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) ) ;
		//printf( "3: %d\n", bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )  ;
		return 0 ;
	}

	int IsIn( uint64_t code[], int kmerLength ) 
	{
		//printf( "1. %llu\n", val ) ;
		GetCanonicalKmerCode( code, kmerLength ) ;
		//printf( "2. %llu\n", val ) ;
		if ( method == 1 )
		{
			//return ( hash.find( val ) != hash.end() ) ;
		}
		return bf.contains( ( char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) ) ;
	}
#endif 
	int numOfThreads ;
public:
	//uses paper defaults
	StoreSF( size_t D=5, size_t W=40000, size_t Z=3, size_t Bits_c=8*sizeof(size_t) )
	{
		//D=# of arrays in both S and F
		//W=# of buckets in S
		//W'=Z x W where W' is # of buckets in F
		//Bits_c = # of bits per counter in F
		numOfThreads = 1 ;
		method = 0 ;
		//size_t D=5,W=40000,Z=3,Bits_c=8*sizeof(size_t);
		init(D, W, Z, Bits_c);
	}

	~StoreSF() 
	{
	}
	
	double Occupancy()
	{
		return 0;
		//return bf.occupancy() ;
	}
	double GetFP()
	{
		return 0;
		
		/*if ( method == 1 )
			return 0 ;
		return bf.GetActualFP() ;*/
	}

	void decrease(KmerCode code, uint64_t count)
	{
		if ( !code.IsValid() )
			return;
#if MAX_KMER_LENGTH <= 32 
		decrease(code.GetCode(), code.GetKmerLength(), count);
#else
		uint64_t c[K_BLOCK_NUM] ;
		code.GetCode( c ) ;
		decrease(code.GetCode(), code.GetKmerLength(), count);
#endif
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
#if MAX_KMER_LENGTH <= 32 
		Put( code.GetCode(), code.GetKmerLength(), testFirst ) ;
#else
		uint64_t c[K_BLOCK_NUM] ;
		code.GetCode( c ) ;
		Put( c, code.GetKmerLength(), testFirst ) ;
#endif
		return 0 ;
	}
	
	int IsIn( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
#if MAX_KMER_LENGTH <= 32 
		return IsIn( code.GetCode(), code.GetKmerLength() ) ;
#else
		uint64_t c[K_BLOCK_NUM] ;
		code.GetCode( c ) ;
		return IsIn( c, code.GetKmerLength() ) ;
#endif
	}


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
