
#include <initializer_list>
#include "la.h"

int main()
{


	Matrix<double> B( 3, 3, { 
				2, 1, 4,
				4, 4, -1,
				3, 4, 5 } );
				
        B.print() << std::endl << B.determinant() << std::endl << B << std::endl;
		
	Matrix<double> A( 4, 4, { 
				2, 1, 4, 5,
				2, 3, 4, 5,
				2, 4, 4, -1,
				5, 3, 4, 5 } );

        A.print() << std::endl << A.determinant() << std::endl << A << std::endl;
				
				
	Matrix<double> C( 5, 5, { 
				2, 1, 4, 5, 3,
				2, 3, 4, 5, -1,
				2, 4, 4, -1, 5,
				5, 3, 4, 5, 9,
				5, 8, 4, 5, 7
				 } );

        C.print() << std::endl << C.determinant() << std::endl << C << std::endl;
				
				





}


