#ifndef _MOURISL_CLASS_STORESF
#define _MOURISL_CLASS_STORESF
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>
#include <algorithm>

#include "bf.h"
#include "KmerCode.h"


class StoreBF
{
private:
	uint64_t size;
	static const uint64_t width=4;
	static const uint64_t num_hashes=3;
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map
	bf::spectral_mi_bloom_filter cbf;
	bf::spectral_mi_bloom_filter** filters;

	int* sizes;
	int* heights;
	
	int numOfThreads;
	int num_filters;

	//all Morris approx. counting related code is adapted from the SF paper/codebase for CML
	double morris_base;
	uint64_t* limits;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;


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

	bool morris_choice(int cur_count)
	{
		double r = distribution(generator);
		double lim = pow(morris_base, -cur_count);
		return r < lim;
	}

	double morris_point_query(int cur_count)
	{
		return cur_count == 0 ? 0 : pow(morris_base, cur_count - 1);
	}
	
	uint64_t tomb_query(uint64_t val, int kmerLength, int* filter_layer, bool unique=true)
	{
		if(kmerLength > 0)
			val = GetCanonicalKmerCode( val, kmerLength ) ;
		/*track for counting
		if(hash.find( val ) != hash.end())
			return 0;
		hash[ val ] = 1;*/
		int i;
		uint64_t running_count = 0;
		for(i=0;i<num_filters;i++)
		{
			//TODO: make sure this returns 0 for non-existing values
			size_t c = filters[i]->lookup( val );
			running_count += c;
			if(c < limits[i])
			{
				*filter_layer = i;
				break;
			}
		}
		if(unique && hash.find( val ) != hash.end())
		{
			hash.erase(val);
			return running_count;
		}
		else if(!unique)
			return running_count;
		return 0;
	}

	bool tomb_insert(uint64_t val, int kmerLength)
	{
		int filter_layer = -1;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		uint64_t count = tomb_query(val, 0, &filter_layer, false);
		//all filled or there's a weird error
		if(filter_layer == -1)
			return false;
		if(morris_choice(count))
			filters[filter_layer]->add( val );
		if(hash.find( val ) != hash.end())
		{
			hash[ val ] += 1;
		}
		else
		{
			hash[ val ] = 1;
		}
		return true;
	}



public:
	//StoreBF( uint64_t s, double fprate = 0.01, int num_filters, int* sizes, int* widths ): size(s), cbf(bf::make_hasher(num_hashes), size, width)
	StoreBF( uint64_t s, int num_filters_ = 4, double base = 1.08 ): size(s), num_filters(num_filters_), morris_base(base), cbf(bf::make_hasher(num_hashes), size, width)
	{
		int sizes_[4] = {1000,500,100,100};
		int heigths_[4] = {1,4,4,4};

		sizes = &sizes_[0];	
		heights = &heigths_[0];

		numOfThreads = 1 ;
		method = 0 ;
		filters = new bf::spectral_mi_bloom_filter*[num_filters];
		limits = new uint64_t[num_filters];
		int i;
		for(i=0;i<num_filters;i++)
		{
			filters[i] = new bf::spectral_mi_bloom_filter(bf::make_hasher(num_hashes), sizes[i], heights[i]);
			limits[i] = pow(morris_base, heights[i]);
		}
		//cbf = bf::spectral_mi_bloom_filter(bf::make_hasher(num_hashes), s, width);
	}

	~StoreBF() 
	{
		delete filters;
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
	
	bool TOMB_Put( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		return tomb_insert( code.GetCode(), code.GetKmerLength() ) ;
	}
	
	
	int IsIn( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		return IsIn( code.GetCode(), code.GetKmerLength() ) ;
	}

	uint64_t TOMB_query( KmerCode &code )
	{
		if ( !code.IsValid() )
			return 0 ;
		int filter_layer = -1;
		uint64_t running_count = tomb_query(code.GetCode(), code.GetKmerLength(), &filter_layer);
		if(filter_layer == -1)
			return -1;
	//for querying with Morris, use SF:CML code:
	//smallest_counter <= 1 ? morris_point_query(smallest_counter) : (int)(round((1 - morris_point_query(smallest_counter + 1)) / (1 - morris_base)));
		//TODO: need to double check this update which takes the filter_layer(s) into consideration to get accurate count since each level is an order of magnitude
		//uint64_t cur_layer_count = count <= 1 ? (uint64_t)morris_point_query(count) : (uint64_t)(round((1 - morris_point_query(count + 1)) / (1 - morris_base)));
		//if(filter_layer > 0)
		//	cur_layer_count += pow(morris_base, filter_layer);
		//THIS seems more in line with the Van Durme paper
		//the counts across all layers' slots are simply added together *then* they're used as the morris exponent
		uint64_t cur_count = running_count <= 1 ? (uint64_t)morris_point_query(running_count) : (uint64_t)(round((1 - morris_point_query(running_count + 1)) / (1 - morris_base)));
		return cur_count;
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
