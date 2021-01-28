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



using namespace std;


int main( int argc, char **argv)
{

	int nums = argc > 1 ? atoi(argv[1]) : 1000;
	thrust::device_vector<int> vint(nums-1);
	thrust::device_vector<int> pint(nums-1);

	thrust::sequence( vint.begin(), vint.end(), 2 );

	thrust::transform( vint.begin(), vint.end(), pint.begin(), 
		[](int x) -> int 
			{ 
				for( int i=2; i < x; i++ ) 
					if( !(x%i) ) 
						return -1; 
				return x; 
			} 
		);


	for( auto i : pint )
		if( i >= 0 )
			cout << i << endl;	

}
