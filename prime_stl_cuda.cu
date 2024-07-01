
#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>
//#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <iostream>

// $>  nvcc -std=c++11 -o prime_stl_cuda prime_stl_cuda.cu  -gencode arch=compute_61,code=sm_61


using namespace std;


struct prime_check
{

__device__ int operator()( int x )
			{ 
				for( int i=2; i < x; i++ ) 
					if( !(x%i) ) 
						return -1; 
				return x; 
			} ;
};

int main( int argc, char **argv)
{

	int nums = argc > 1 ? atoi(argv[1]) : 1000;
	thrust::device_vector<int> vint(nums-1);
	thrust::device_vector<int> pint(nums-1);

	thrust::sequence( vint.begin(), vint.end(), 2 );
	thrust::transform( vint.begin(), vint.end(), pint.begin(), prime_check() );
	auto new_end = thrust::remove( pint.begin(), pint.end(), -1 );
	thrust::copy( pint.begin(), new_end, std::ostream_iterator<int>(std::cout, "\n") );

//	for( auto i : pint )
//		if( i >= 0 )
//			cout << i << endl;	

}
