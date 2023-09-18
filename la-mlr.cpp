#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <cmath>
#include <float.h>
#include <initializer_list>
#include "la2.h"


using namespace std;

template <typename T> int VectorTest()
{

	const double PI = 3.141592653589793238463;
	const double PI2 = PI * 2.0;
	
	{
		
		cout << endl << "Linear Algebra Regression tests:" << endl;
		
		Matrix<T> mY( {
				{  9.0 },
				{ 10.0 },
				{ 13.0 },
				{ 14.0 },
				{ 16.0 }
			} );
		

		Matrix<T> mX (  {	
				{ 1.0, 1.0, 10.0 },
				{ 1.0, 3.0, 14.0 },
				{ 1.0, 4.0, 15.0 },
				{ 1.0, 6.0, 18.0 },
				{ 1.0, 7.0, 20.0 }	
			} );
		

		cout << mX << endl;

		auto mXt = mX;
		
		cout << mXt.transpose() << endl;
		


		try{
			cout << (mXt*mX) << endl << endl;
		} catch ( linear_algebra_error e )
		{
			cerr << endl << e.what() << endl;
		}
		
	}
	
	return 0;
	
	
	

}




int main()
{

	VectorTest<double>();
	
	return 0;
}


