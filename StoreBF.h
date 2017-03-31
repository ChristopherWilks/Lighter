#ifndef _MOURISL_CLASS_STORESF
#define _MOURISL_CLASS_STORESF
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "bloom_filter.hpp"
#include "bf.h"
#include "KmerCode.h"


const static int SF_LENGTH=8;


class StoreBF
{
private:
	uint64_t size;
	static const uint64_t width=4;
	static const uint64_t num_hashes=3;
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map
	bf::spectral_mi_bloom_filter cbf;


	void increase(uint64_t val, int kmerLength, size_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		cbf.add( val );
	}

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		//printf( "%llu\n", val ) ;

		//if ( numOfThreads > 1 && testFirst && bf.contains( val ) )
		if ( testFirst && cbf.lookup( val ) > 0 )
			return 0 ;
		//bf.insert( val ) ;
		cbf.add( val );
		//int c = query(valc, SF_LENGTH);
		//printf("count of %d\n",c);
		return 1 ;
	}
	
	size_t IsIn( uint64_t val, int kmerLength ) 
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		
		return cbf.lookup( val ) ;
	}

	void decrease(uint64_t val, int kmerLength, uint64_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;

	        cbf.remove( val );	
	}

	int numOfThreads ;
public:
	StoreBF( double fprate = 0.01 ): size(10000003), cbf(bf::make_hasher(num_hashes), size, width)
	{
		size = 10000003;
		numOfThreads = 1 ;
		method = 0 ;
		//cbf = bf::spectral_mi_bloom_filter(bf::make_hasher(num_hashes), size, width);
	}

	StoreBF( uint64_t s, double fprate = 0.01 ): size(s), cbf(bf::make_hasher(num_hashes), size, width)
	{
		numOfThreads = 1 ;
		method = 0 ;
		//cbf = bf::spectral_mi_bloom_filter(bf::make_hasher(num_hashes), s, width);
	}

	~StoreBF() 
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
		decrease(code.GetCode(), code.GetKmerLength(), count);
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


	void SetNumOfThreads( int in ) 
	{ 
		numOfThreads = in ;
		//bf.SetNumOfThreads( in ) ;
	}

	uint64_t GetCanonicalKmerCode( uint64_t code, int k ) 
	{
		int i ;
		uint64_t crCode = 0ull ; // complementary code
		for ( i = 0 ; i < k ; ++i )
		{
			uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
			crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
		}
		return crCode < code ? crCode : code ;
	}

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
