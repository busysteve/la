
#include <initializer_list>
#include "la.h"



void print_wxMaxima_solver_query( const Matrix<double> &B1, const Matrix<double> &B2 )
{

	char alpha[] = "abcdefghijklmnopqrstuvwxyz";
	int alen = strlen( alpha );
	std::cout << "float(linsolve([";
	for( int r=0; r < B1.nr(); r++ )
	{
		if( r > 0 )
			std::cout << ",";
		for( int c=0; c < B1.nc(); c++ )
		{

			if( B1(r,c) >= 0.0 && c > 0 )
				std::cout << "+";
			std::cout << B1(r,c) << "*" << alpha[alen-(B1.nc()-c)];
		}
		
		std::cout << "=";
		std::cout << B2(r,0);
	}
	std::cout << "],[";
	
	for( int c=0; c < B1.nc(); c++ )
	{

		if( c > 0 )
			std::cout << ",";
		std::cout << alpha[alen-(B1.nc()-c)];
	}

	
	std::cout << "]));" << std::endl;

}


void print_WolframAlpha_solver_query( const Matrix<double> &B1, const Matrix<double> &B2 )
{

	char alpha[] = "abcdefghijklmnopqrstuvwxyz";
	int alen = strlen( alpha );
	std::cout << "{";
	for( int r=0; r < B1.nr(); r++ )
	{
		if( r > 0 )
			std::cout << ",";
		for( int c=0; c < B1.nc(); c++ )
		{

			if( B1(r,c) >= 0.0 && c > 0 )
				std::cout << "+";
			std::cout << B1(r,c);
			std::cout << alpha[alen-(B1.nc()-c)];
		}
		
		std::cout << "=";
		std::cout << B2(r,0);
	}
	std::cout << "}" << std::endl;

}


int main()
{

	{
		Matrix<double> B1( 2, 2, {  
					4, -1,
					4, 5 } );
					
		Matrix<double> B2( 2, 1, { 
					  8,
					  3 } );
		                          
					
		B1.print() << std::endl;
		B2.print() << std::endl;
		B1.GJ_elimination(B2).print() << std::endl << std::endl;
			
		std::cout << B1 << std::endl;
		std::cout << B2 << std::endl;

		print_WolframAlpha_solver_query( B1, B2 );
		print_wxMaxima_solver_query( B1, B2 );
	}


	{
		Matrix<double> B1( 3, 3, { 
					2, 1, 4,
					4, 4, -1,
					3, 4, 5 } );
					
		Matrix<double> B2( 3, 1, { 
					  8,
					  5,
					  3 } );
		                          
					
		B1.print() << std::endl;
		B2.print() << std::endl;
		B1.GJ_elimination(B2).print() << std::endl << std::endl;
			
		std::cout << B1 << std::endl;
		std::cout << B2 << std::endl;

		print_WolframAlpha_solver_query( B1, B2 );
		print_wxMaxima_solver_query( B1, B2 );
	}


	{
		Matrix<double> B1( 4, 4, { 
					2, 1, 4, 6,
					4, 4, -1, 2,
					7, -2, 4, 3,
					3, 4, 5, 7 } );
					
		Matrix<double> B2( 4, 1, { 
					  8,
					  5,
					  7,
					  3 } );
		                          
					
		B1.print() << std::endl;
		B2.print() << std::endl;
		B1.GJ_elimination(B2).print() << std::endl << std::endl;
			
		std::cout << B1 << std::endl;
		std::cout << B2 << std::endl;

		print_WolframAlpha_solver_query( B1, B2 );
		print_wxMaxima_solver_query( B1, B2 );
	}


	{
		Matrix<double> B1( 5, 5, { 
					7, 2, 1, 4, 6,
					7, 4, 4, -1, 2,
					-5, 3, 8, 2, 1,
					9, 7, -2, 4, 3,
					2, 3, 4, 5, 7 } );
					
		Matrix<double> B2( 5, 1, { 
					  3,
					  5,
					  9,
					  7,
					  3 } );
		                          
					
		B1.print() << std::endl;
		B2.print() << std::endl;
		B1.GJ_elimination(B2).print() << std::endl << std::endl;
			
		std::cout << B1 << std::endl;
		std::cout << B2 << std::endl;

		print_WolframAlpha_solver_query( B1, B2 );
		print_wxMaxima_solver_query( B1, B2 );
	}

}


