
// g++ -g -o ann ANN.cpp XMLTag/xmltag.cpp
// ./ann -w test.weights.xml -r 0.00002 -m 0.0002 -t train.txt -e 10 -i input.txt -l S2 S3 S2 S1
// or
// ./ann -w test.weights.xml -i input.txt

/* train.txt
0 0 0
0 1 1
1 0 1
1 1 0
*/

/* input.txt
0 0
0 1
1 0
1 1
*/

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
//#include "XMLTag/xmltag.h"
#include "la.h"


//#define log_verbose
//#define log_verbose   printf
#define log_verbose  if( g_verbose > 0 ) if( g_counter%g_verbose == 0 ) printf

#define log_output  if( g_output > 0 ) if( g_counter%g_output == 0 ) printf

#define MAX_NN_NAME 30

int g_output = 0;
int g_verbose = 0;
int g_counter = 0;

//const void* nullptr = NULL;

enum ActType{ linear = 0, sigmoid, tangenth, none, bias };

template<typename T>
T actNone( T n )
{
	return n;
}


template<typename T>
T actBias( T n )
{
	return (T)1.0;
}


template<typename T>
T actLinear( T n )
{
	return n;
}


template<typename T>
T actSigmoid( T n )
{
	return 1.0 / ( 1.0 + exp(-n) );
}


template<typename T>
T actTanh( T n )
{
	return tanh( n );
}


template<typename T>
T derivLinear( T n )
{
	return 1.0;
}


template<typename T>
T derivSigmoid( T n )
{
	return n * ( 1.0 - n );
}


template<typename T>
T derivTanh( T n )
{
	return 1.0 - n * n;
}


template<typename T>
T safeguard( T n )
{
	return n;
	//return n != 0.0 ? n : 0.0000000000001;
}


