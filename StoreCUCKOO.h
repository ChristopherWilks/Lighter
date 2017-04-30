#ifndef _MOURISL_CLASS_STORE_CUCKOO
#define _MOURISL_CLASS_STORE_CUCKOO
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "bloom_filter.hpp"
#include "KmerCode.h"
#include "cuckoofilter.h"

using cuckoofilter::CuckooFilter;

class StoreCUCKOO
{
private:
	uint64_t num_unique_kmers;
	static const int num_bits_per_key=12;
	CuckooFilter<uint64_t, num_bits_per_key> cf;
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		/*if ( method == 1 )
		{
			hash[ val ] = 1 ;
			return 1 ;
		}*/
		
		if ( numOfThreads > 1 && testFirst && cf.Contain(val) == cuckoofilter::Ok )
			return 0 ;
		if(cf.Add(val) != cuckoofilter::Ok)
		{
			fprintf(stderr,"couldnt insert %u into cuckoo filter\n",val);
			exit(-1);
		}
		//fprintf(stderr,"inserted %u into cuckoo filter\n",val);
		return 1;
	}

	int IsIn( uint64_t val, int kmerLength ) 
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		/*(if ( method == 1 )
		{
			return ( hash.find( val ) != hash.end() ) ;
		}*/
		return cf.Contain(val) == cuckoofilter::Ok;
	}
	int numOfThreads ;
public:
	//Store( uint64_t num_unique_kmers_, int num_bits_per_key_ ): num_unique_kmers(num_unique_kmers_), num_bits_per_key(num_bits_per_key_)
	StoreCUCKOO( uint64_t num_unique_kmers_ ): num_unique_kmers(num_unique_kmers_), cf(num_unique_kmers)
	{
		numOfThreads = 1 ;
		method = 0 ;
		//cf(num_unique_kmers);
	}

	~StoreCUCKOO() 
	{
	}
	
	double Occupancy()
	{
		return 0.0;
	}
	double GetFP()
	{
		//return bf.GetActualFP() ;
		return 0;
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
	}

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

	void TemporaryOutput( char *file)
	{
		FILE *fp = fopen( file, "w" ) ;
		fclose( fp ) ;
	}

	void TemporaryInput( char *file ) 
	{
		FILE *fp = fopen( file, "r" ) ;
		fclose( fp ) ;
	}

	int Clear() 
	{
		return 0 ;
	}
} ;
#endif
