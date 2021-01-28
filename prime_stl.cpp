#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>


using namespace std;


int main( int argc, char **argv)
{

	int nums = argc > 1 ? atoi(argv[1]) : 1000;
	vector<int> vint(nums-1);
	vector<int> pint(nums-1);

	iota( vint.begin(), vint.end(), 2 );

	transform( vint.begin(), vint.end(), pint.begin(), 
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
