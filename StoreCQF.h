#ifndef _MOURISL_CLASS_STORESF
#define _MOURISL_CLASS_STORESF
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "gqf.h"
#include "hashutil.h"

#include "KmerCode.h"

class StoreCQF
{
private:
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map
	QF cqf;

	void increase(uint64_t val, int kmerLength, uint64_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		val = kmercounting::HashUtil::MurmurHash64A(((void*)&val), sizeof(val), cqf.metadata->seed);
		qf_insert(&cqf, val % cqf.metadata->range, 0, count, false, false);
	}

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		// hash the kmer using murmurhash/xxHash before adding to the list
		val = kmercounting::HashUtil::MurmurHash64A(((void*)&val), sizeof(val), cqf.metadata->seed);
		//if ( numOfThreads > 1 && testFirst && bf.contains( val ) )
		if ( numOfThreads > 1 && testFirst && qf_count_key_value(&cqf, val % cqf.metadata->range, 0))
			return 0 ;
		//bf.insert( val ) ;
		//insert/update the val as a cqf key with a count of 1, don't lock, don't spin
		qf_insert(&cqf, val % cqf.metadata->range, 0, 1, false, false);
		//int c = query(valc, SF_LENGTH);
		//printf("count of %d\n",c);
		return 1 ;
	}
	
	uint64_t IsIn( uint64_t val, int kmerLength ) 
	{
		//printf( "1. %llu\n", val ) ;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		val = kmercounting::HashUtil::MurmurHash64A(((void*)&val), sizeof(val), cqf.metadata->seed);
		//printf( "2. %llu\n", val ) ;
		
		//return bf.contains( val ) ;
		//return qf_count_key_value(cqf, val, 0);
		return qf_count_key_value(&cqf, val % cqf.metadata->range, 0);
	}


	void decrease(uint64_t val, int kmerLength, uint64_t count)
	{
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		val = kmercounting::HashUtil::MurmurHash64A(((void*)&val), sizeof(val), cqf.metadata->seed);
		//qf_remove(&cqf, val, 0, count);	
	}
	
	int numOfThreads ;
public:
	StoreCQF( )
	{
		numOfThreads = 1 ;
		method = 0 ;
		int qbits = 25; //# of distinct kmers dependent //20 is example default
		int r_bits = 17; //default from squeakr //8 is example default
		int num_hash_bits = qbits + r_bits;	
		uint32_t seed = time(NULL);
		qf_init(&cqf, (1ULL<<qbits), num_hash_bits, 0, true, "", seed);
	}

	~StoreCQF() 
	{
		qf_destroy(&cqf, true);
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

	int increase(KmerCode &code, uint64_t count) 
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
	
	uint64_t count_unique_keys_gt_cutoff(int cutoff)
	{
		QFi cfi;
		uint64_t total_count = 0;
		if(qf_iterator(&cqf, &cfi, 0))
		{
			do {
				uint64_t key = 0, value = 0, count = 0;
				qfi_get(&cfi, &key, &value, &count);
				if(count > cutoff)
					total_count += 1;
			} while (!qfi_next(&cfi));
		}
		return total_count;
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
