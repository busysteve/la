
#include "nn.h"


using namespace std;



int main()
{

	vector< Vector<double> > ins;
	vector< Vector<double> > outs;

	

	ins.push_back( { 0, 0 } );
	outs.push_back( { 0 } );

	ins.push_back( { 0, 1 } );
	outs.push_back( { 1 } );

	ins.push_back( { 1, 0 } );
	outs.push_back( { 1 } );



	ins.push_back( { 1, 1 } );
	outs.push_back( { 0 } );

	Matrix<double> w1 ( 3, 2, { .8, .2, .4, .9, .3, .5 } );
	//w1.rand( 0, 1 );

	Matrix<double> w2 ( 1, 3, { .3, .5, .9 } );
	//w2.rand( 0, 1 );

	//Vector<double> hl(3);
	//Vector<double> ol(1);
	
	for( int e=0; e < 300; e++ )
	{
		for( int i=0; i < ins.size(); i++ )
		{
			Vector<double> h1 = w1 * ins[i];
			ins[i].print();
			outs[i].print();


			h1.func( actSigmoid<double> );
			//r1.print() << endl;
			
			Vector<double> ro = w2 * h1;
			Vector<double> ao = ro;
			
			ao.func( actSigmoid<double> );
			ao.print() << "Activated outputs" << endl << endl << endl;
			
			// ===================
			
			Vector<double> dso = ao;
			dso.func( derivSigmoid<double> );
			
			Vector<double> eo = ao;
			
			double lr = .2;
			
			for( int o=0; o < eo.size(); o++ )
				eo[o] = ( outs[i][o] - ao[o] ) * lr;
				
			
			cout << eo << " - Error " << endl;
			cout << ro << " - Out Sum " << endl;
			// ===================

			dso.print() << "Deriv Outputs " << endl;
			
			
			
			dso = dso * eo;
			dso.print() << "Delta Outputs " << endl;
			
			

			w2.print() << endl;
			cout << h1 << " - Hidden Layer" << endl;
			Matrix<double> hd( dso[0] / h1 );
			cout << dso[0] << " / " << h1 << " = " << hd << endl;
			cout << hd << " - Hidden Delta Weights" << endl;
			
			auto w2n = w2 + hd;
			
			w2n.print() << endl;




			cout << endl << endl << endl;




			Matrix<double> ha( h1 );
			ha.func( derivSigmoid<double> );
			Matrix<double> id( ( dso[0] / w2 ).SchurProd( ha ) );
			cout << dso[0] << " / " << w2 << " * " << ha << " = " << id << endl;
			
			
			
			Matrix<double> w1n( id[0].outerDiv( ins[i] ) );
			cout << id << " / " << ins[i] << " = " << w1n << endl;
			
			
/*
			id.transpose();
			Matrix<double> is(ins[i]);
			Matrix<double> w1n( id * is );
			cout << id << " * " << is << " = " << w1n << endl;
*/
			//w1 = ( h1 / ins[i] ) ;
			cout << w1 << " -> " << w1n << endl;
			
			w2 = w2n;
			w1 += w1n;
			
			cout << w1 << endl;
		}
		
	}


}




