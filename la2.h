#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <float.h>
#include <exception>
#include <functional>
#include <algorithm>
#include <future>
#include <thread>
#include <iterator>
#include <numeric>
#include <initializer_list>
#include <iterator>
#include <tuple>
#include <sstream>
#include <iterator>
#include <execution>

#ifndef __BS_Matrix_H
#define __BS_Matrix_H

class linear_algebra_error : public std::exception {

protected:
    std::string msg;

public:
    const char* what() const throw()
    {
        return msg.c_str();
    } //virtual const char * what () = 0;
};

class vector_algebra_error : public linear_algebra_error {
public:
    vector_algebra_error(const char* em)
    {
        msg = em;
    }
};

class matrix_algebra_error : public linear_algebra_error {
public:
    matrix_algebra_error(const char* em)
    {
        msg = em;
    }
};

const double PI = 3.141592653589793238463;
const double PI2 = PI * 2.0;

template <typename T>
class Vector;
template <typename T>
class Matrix;

template <typename From, typename To>
void convert_vector_copy(Vector<From>& f, Vector<To>& t)
{
    t.clear();
    t.resize(f.size());

    for (int i = 0; i < f.size(); i++)
        t[i] = f[i];
}

template <typename From, typename To>
void convert_matrix_copy(Matrix<From>& f, Matrix<To>& t)
{
    t.clear();
    t.resize(f.nr(), f.nc());

#pragma omp parallel for
    for (int i = 0; i < f.nr(); i++)
        for (int j = 0; j < f.nc(); j++)
            t[i][j] = f[i][j];
}




template <typename T>
class Matrix {

	


protected:

    T *_m = nullptr;

    size_t _rs;
    size_t _cs;

public:
    const double PI = 3.141592653589793238463;
    const double PI2 = PI * 2.0;



    Matrix<T>()
        : _rs(0)
        , _cs(0)
    {
    }

    Matrix<T>& resize(unsigned int rows, unsigned int cols, const T& x = (T)0.0)
    {

	delete _m;
	
	_m = nullptr;
	
	_m = new T[rows*cols];

	std::fill( _m, _m+(cols*rows), x );

        _rs = rows;
        _cs = cols;

        return (*this);
    }



    Matrix(size_t rows, size_t cols, const T& x = (T)0.0)
    {
        resize(rows, cols, x);
    }

    Matrix(size_t rows, size_t cols,  std::initializer_list<T> il)
    {
        set_rows_and_cols(rows, cols, il);
    }

    Matrix(size_t rows, size_t cols,  std::initializer_list< std::initializer_list<T> > il)
    {
        set_rows_and_cols(rows, cols, il);
    }

    Matrix( std::initializer_list< std::initializer_list<T> > il)
    {
    
    	int r=0;
    	int c=0;
    	int C=0;
    	
    	for( auto i : il )
    	{
    		if( c!=0 && C == 0 )
    			C = c;
    			
    		else if( C != c )
    			throw  matrix_algebra_error("Invalid Initializer Matrix");
    			
    		c = 0;
    			
    		r++;
    		for( auto j : i )
    			c++;
    	}
    
        resize(r, c, (T)0.0);
        
        Matrix<T>& M = (*this);

	int x=0;
	int y=0;
    	for( auto i : il )
    	{	
    		y=0;
    		for( auto j : i )
    		{
    			*( _m+(c*x+y) ) = j;
    			y++;
    		}
    		x++;
    	}

    }

    void set_cols_and_rows(size_t c, size_t r,  std::initializer_list<T> l)
    {
        resize(r, c, (T)0.0);
        
        Matrix<T>& M = (*this);

        auto k = l.begin();
        for (int i = 0; i < M.nc(); i++)
            for (int j = 0; j < M.nr(); j++, k++)
                M._m+(r*j+i) = *(k);
    }

    void set_rows_and_cols(size_t r, size_t c,  std::initializer_list<T> l)
    {
        resize(r, c, (T)0.0);
        
        Matrix<T>& M = (*this);

        auto k = l.begin();
        for (int i = 0; i < M.nc(); i++)
            for (int j = 0; j < M.nr(); j++, k++)
                *(M._m+(c*i+j)) = *(k);
    }


    void set_rows_and_cols(size_t r, size_t c,  std::initializer_list< std::initializer_list<T> > l)
    {
        resize(r, c, (T)0.0);
        
        Matrix<T>& M = (*this);

        auto k = l.begin();
        for (int i = 0; i < M.nc(); i++, k++ )
        {
            auto n = (*k).begin();
            for (int j = 0; j < M.nr(); j++, n++)
                *(M._m+(c*i+j)) = *(n);
        }
    }


    /*
    Matrix( initializer_list< initializer_list<T> > il )
    {
        for (auto l : il)
            _m.push_back(l);

        _rs = _m.size();
        _cs = 0;

        if (_rs > 0)
            _cs = _m[0].size();
    }
*/

    Matrix(const Matrix<T>& r)
    {
        (*this) = r;
    }

    virtual ~Matrix() { delete _m; }

    // recursively find determinant
    T determinant() const
    {
        const Matrix<T>& A(*this);
        T r = 0;

        if (nr() > 2 || nc() > 2) {

            auto lamSubDeterMatrix = [](const Matrix<T>& M, int o) -> Matrix<T> {
                Matrix<T> R(M.nr() - 1, M.nc() - 1);
                for (int i = 0, I = 0; i < M.nr(); i++)
                    if (i != o) {

                        for (int j = 1, J = 0; j < M.nc(); j++)
                            R(I, J++) = M(i, j);
                        I++;
                    }
                return R;
            };

            for (int i = 0; i < A.nr(); i++) {
                T p = ((i & 1) ? -1 : 1);
                Matrix<T> S(lamSubDeterMatrix(A, i));
                p *= A(i, 0) * S.determinant();
                r += p;
            }

            return r;
        }
        else if (nr() == 2 || nc() == 2) {
            return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
        }
        else if (nr() == 1 || nc() == 1) {
            return ( A(0, 0) != 0.0 ) ? 1 : 0;
        }


            throw matrix_algebra_error(
                "The matrix is not square for calculaing the determinant");
       
    }


    Matrix<T> adjugate() const
    {
        const Matrix<T>& A(*this);
        Matrix<T> R( A );
        T r = 0;

        if ( nr() == nc() && ( nr() > 2 || nc() > 2 ) ) {

            auto lamCofactorMatrix = [](const Matrix<T>& M, int o, int p ) -> Matrix<T> {
                Matrix<T> R(M.nr() - 1, M.nc() - 1);
                for (int i = 0, I = 0; i < M.nr(); i++) {
                    if (i != o ) {
                        for (int j = 0, J = 0; j < M.nc(); j++) {
                            if ( j != p) 
                            	R(I, J++) = M(i, j);
                        }
                        I++;
                    }
                }
                return R;
            };
            
            for (int i = 0; i < R.nr(); i++) {
                T p = ((i & 1) ? 1 : -1);
                for (int j = 0; j < R.nc(); j++) {
                    p *= -1;
                    Matrix<T> S(lamCofactorMatrix(A, i, j));
                    R(i, j) = p * S.determinant();
                }
            }

            R.transpose();

            return R;
        }

        throw matrix_algebra_error(
                "The matrix is not square for calculaing the adjugate");

    }

    Matrix<T> inverse() const
    {
        const Matrix<T>& A(*this);

        if ( nr() == nc() ) {
            Matrix<T> R( A.adjugate() );
            R *= ( 1.0 / A.determinant() );
            
            return R;
        }

        throw matrix_algebra_error(
                "The matrix is not square for calculaing the adjugate");
    }

    // Assignment Operator

    Matrix<T>& operator=(const Matrix<T>& rhs)
    {
        if (&rhs == this)
            return *this;

        unsigned int new_rows = rhs.nr();
        unsigned int new_cols = rhs.nc();

	delete _m;
	
	_m = nullptr;
	
        _m = new T[new_rows*new_cols];
        
        for (unsigned int i = 0; i < new_rows; i++) {
		#pragma omp parallel for
            for (unsigned int j = 0; j < new_cols; j++) {
                std::copy( rhs._m, rhs._m+(new_rows*new_cols), _m );
            }
        }
        _rs = new_rows;
        _cs = new_cols;

        return *this;
    }

    // Equality comparison
    bool operator==(const Matrix<T>& rhs)
    {
        if (&rhs == this)
            return true;

        bool ret = true;

        unsigned int new_rows = rhs.nr();
        unsigned int new_cols = rhs.nc();

	#pragma omp parallel for
        for (unsigned int i = 0; i < new_rows; i++) {
            for (unsigned int j = 0; j < new_cols; j++) {
                if ( *(_m+(new_cols*i+j)) != *(rhs._m+(new_cols*i+j)) ) {
                    ret = false;
                    break;
                }
            }
        }

        return ret;
    }

    // Addition

    Matrix<T> operator+(const Matrix<T>& rhs) const 
    {

        if (_rs != rhs._rs || _cs != rhs._cs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for addition");

        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_cs*i+j)) = *(_m+(_cs*i+j) + rhs._m+(_cs*i+j));
            }
        }

        return result;
    }


    Matrix<T> cross_add(const Matrix<T>& rhs) const 
    {

        if (_rs != rhs._rs || _cs != rhs._cs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for addition");

        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_cs*i+j)) = *(_m+(_cs*i+j) + rhs._m+(_cs*i+j));
            }
        }

        return result;
    }

    // Cumulative addition of this Matrix and another

    Matrix<T>& operator+=(const Matrix<T>& rhs)
    {
        unsigned int _rs = rhs.nr();
        unsigned int _cs = rhs.nc();

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *(_m+(_cs*i+j)) += *(rhs._m+(_cs*i+j));
            }
        }

        return *this;
    }

    // subtraction of a value from each element ofthis Matrix
    Matrix<T> operator-(T s)
    {

        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {

            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_cs*i+j)) = *(_m+(_cs*i+j)) - s;
            }
        }

        return result;
    }

    Matrix<T> operator-(const Matrix<T>& rhs) const
    {

        if (_rs != rhs._rs || _cs != rhs._cs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for addition");

        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {

            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_cs*i+j)) = *(_m+(_cs*i+j)) - *(rhs._m+(_cs*i+j));
            }
        }

        return result;
    }


    Matrix<T> operator-() const
    {

        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {

            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_cs*i+j)) = -*(_m+(_cs*i+j));
            }
        }

        return result;
    }

    // Cumulative subtraction of this Matrix and another

    Matrix<T>& operator-=(const Matrix<T>& rhs)
    {
        unsigned int _rs = rhs.nr();
        unsigned int _cs = rhs.nc();

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *(_m+(_cs*i+j)) -= *(rhs._m+(_cs*i+j));
            }
        }

        return *this;
    }

    // Left multiplication of this Matrix and scalar

    Matrix<T> mul( T s ) const
    {
        unsigned int rs = nr();
        unsigned int cs = nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                *( result._m+(_cs*i+j) ) = s * *(_m+(_cs*i+j));
            }
        }

        return result;
    }

    Matrix<T> div( T s ) const
    {
        unsigned int rs = nr();
        unsigned int cs = nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                *(result._m+(_cs*i+j)) = *(_m+(_cs*i+j)) / s;
            }
        }

        return result;
    }

    Matrix<T> divR( T s ) const
    {
        unsigned int rs = nr();
        unsigned int cs = nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                result(i, j) = s / (*this)(i, j);
            }
        }

        return result;
    }

    class mxIterator 
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        mxIterator(pointer ptr) : m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }
        mxIterator& operator++() { m_ptr++; return *this; }  
        mxIterator operator++(int) { mxIterator tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const mxIterator& a, const mxIterator& b) { return a.m_ptr == b.m_ptr; };
        friend bool operator!= (const mxIterator& a, const mxIterator& b) { return a.m_ptr != b.m_ptr; };  

    protected:
        pointer m_ptr;
    };

    // Define an iterator for the container
    class rowIterator : public mxIterator {
    public:
        rowIterator(T* ptr) : mxIterator(ptr) {  }

        T& operator*() {
            return *mxIterator::m_ptr;
        }

        rowIterator& operator++() {
            mxIterator::m_ptr++;
            return *this;
        }

        bool operator!=(const rowIterator& other) {
            return mxIterator::m_ptr != other.m_ptr;
        }

    };

    rowIterator row_begin(size_t r) const {
        return rowIterator(_m+(_cs*r));
    }

    rowIterator row_end(size_t r) const {
        return rowIterator(_m+(_cs*r+_cs));
    }



    // Define an iterator for the container
    class colIterator : public mxIterator {
    public:
        colIterator(size_t c, T* ptr) : mxIterator(ptr), step_( c ) {}

        T& operator*() {
            return *mxIterator::m_ptr;
        }

        colIterator& operator++() {
            mxIterator::m_ptr += step_;
            return *this;
        }

        bool operator!=(const colIterator& other) {
            return mxIterator::m_ptr != other.m_ptr;
        }

    private:
        size_t step_;
    };

    colIterator col_begin(size_t c) const {
        return colIterator(_cs, _m+(c));
    }

    colIterator col_end(size_t c) const {
        return colIterator(_cs, _m+(_cs*_rs+c));
    }



    Matrix<T> operator*(const Matrix<T>& rhs) const
    {
    
        if (_rs != rhs._cs || _cs != rhs._rs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for multiplication");
    
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        Matrix result(rs, cs, 0.0);

	size_t ra = nr();
	size_t ca = nc();
	size_t rb = rhs.nr();
	size_t cb = rhs.nc();
	

	#pragma omp parallel for
	for ( size_t i = 0; i < ra; i++ )
	{
		#pragma omp parallel for
		for ( size_t j = 0; j < cb; j++ )
		{
			//#pragma omp parallel for
			for ( size_t k = 0; k < ca; k++ )
			{
				auto x = *(_m+(_cs*i+k));
				auto y = *(rhs._m+(rhs._cs*k+j));
				
				*(result._m+(result._cs*i+j)) += x * y;
			}
		}
	}
        return result;
    }


    Matrix<T> mmul(const Matrix<T>& rhs) const
    {
    
        if (_rs != rhs._cs || _cs != rhs._rs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for multiplication");
    
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        Matrix result(rs, cs, 0.0);

	size_t ra = nr();
	size_t ca = nc();
	size_t rb = rhs.nr();
	size_t cb = rhs.nc();
	

	#pragma omp parallel for
	for ( size_t i = 0; i < ra; i++ )
	{
		#pragma omp parallel for
		for ( size_t j = 0; j < cb; j++ )
		{
			//#pragma omp parallel for
			for ( size_t k = 0; k < ca; k++ )
			{
				auto x = *(_m+(_cs*i+k));
				auto y = *(rhs._m+(rhs._cs*k+j));
				
				*(result._m+(result._cs*i+j)) += x * y;
			}
		}
	}
        return result;
    }


    Matrix<T> mul(const Matrix<T>& rhs) const
    {
    
        if (_rs != rhs._rs || _cs != rhs._cs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for cross multiplication");
    
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        Matrix result(rs, cs, 0.0);

	size_t ra = nr();
	size_t ca = nc();
	size_t rb = rhs.nr();
	size_t cb = rhs.nc();
	

	#pragma omp parallel for
	for ( size_t i = 0; i < ra; i++ )
	{
		#pragma omp parallel for
		for ( size_t j = 0; j < ca; j++ )
		{
			auto x = *(_m+(_cs*i+j));
			auto y = *(rhs._m+(rhs._cs*i+j));
			
			*(result._m+(result._cs*i+j)) = x * y;
		}
	}
        return result;
    }


    Matrix<T> it_mul( const Matrix<T>& rhs) // const
    {
        if (_rs != rhs._cs || _cs != rhs._rs)
            //return Matrix<T>();
            throw matrix_algebra_error("Matrix dimensions do not match as needed for multiplication");
    
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        Matrix result(rs, cs, 0.0);

	size_t ra = nr();
	size_t ca = nc();
	size_t rb = rhs.nr();
	size_t cb = rhs.nc();
	

	#pragma omp parallel for
	for ( size_t i = 0; i < ra; i++ )
	{
		#pragma omp parallel for
		for ( size_t j = 0; j < cb; j++ )
		{
		    	*(result._m+(result._cs*i+j)) = 
		    	std::transform_reduce( std::execution::par,
			row_begin(j), row_end(j), rhs.col_begin(i), 0,
				[i,j](const T& a, const T& b ) {
			    		return a + b;
				},
				[i,j](const T& a, const T& b ) {
			    		return a * b;
				}
			);
		}
	}
        return result;
    }


    Matrix<T>& rand( T start=0, T end=0, int seed=0 )
    {
	T diff = end - start;

	if( seed == 0 )
		std::srand( std::time(nullptr) );
	else
		std::srand(seed);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                _m+(_cs*i+j) = start + ((T)std::rand()) / ( ( (T)RAND_MAX  ) / diff );
            }
        }

        return (*this);
    }

    // Matrix/scalar addition

    Matrix<T>& operator+(const T& rhs)
    {
        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *( result._m+(_cs*i+j) ) = *(_m+(_cs*i+j) ) + rhs;
            }
        }

        return result;
    }


    Matrix<T> operator*(const T& rhs) const 
    {
        Matrix result(_rs, _cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *( result._m+(_cs*i+j) ) = *(_m+(_cs*i+j) ) + rhs;
            }
        }

        return result;
    }


    // Matrix transpose
    Matrix<T>& transpose()
    {
        auto cs = _rs;
        auto rs = _cs;
        Matrix result(rs, cs, 0.0);

	#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                *(result._m+(_rs*j+i) ) = *(_m+(_cs*i+j) );
            }
        }

        //delete _m;
        
        //_m = nullptr;
        
        //_m = result._m;
        (*this) = result;
        //_rs = rs;
        //_cs = cs;

        return *this;
    }


    // Access the individual elements
    T& operator()(const unsigned int& row, const unsigned int& col)
    {
        return *(_m+(_cs*row+col) );
    }

    // Access the individual elements (const)
    const T& operator()(const unsigned int& row, const unsigned int& col) const
    {
        return *(_m+(_cs*row+col) );
    }

    // Get the number of rows of the Matrix
    unsigned int nr() const
    {
        return this->_rs;
    }

    // Get the number of columns of the Matrix
    unsigned int nc() const
    {
        return this->_cs;
    }



    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
    {

        //os << "{";

        if (m.nr() > 0) {
            for (int i = 0; i < m.nr(); i++)
            {
                os << std::endl << "{";
		    for (int j = 0; j < m.nc(); j++)
		        os << "," << m(i,j);
                os << "}";
            }
        }
        else
            os << "{}";

        //os << "}" << std::endl;
        os << std::endl;

        return os;
    }


};

#endif
