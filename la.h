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

template <class T>
class Vector : public std::vector<T> {

public:
    //static const double PI = 3.141592653589793238463;
    //static const double PI2 = PI * 2.0;

    Vector<T>() {}

    Vector<T>(int s, T x = (T)(0.0))
    {
        std::vector<T>::resize(s, x);
    }

    Vector<T>(const Vector<T>& r)
    {
        for (auto k : r)
            std::vector<T>::push_back(k);
    }

    Vector<T>(std::initializer_list<T> l)
    {
        for (auto k : l)
            std::vector<T>::push_back(k);
    }

    Vector<T>(std::initializer_list<std::initializer_list<T> >& ll)
    {
        for (auto l : ll)
            std::vector<T>::push_back(T(l));
    }

    Vector<T>(Vector<Vector<T> > ll)
    {
        for (auto l : ll)
            std::vector<T>::push_back(T(l));
    }

    virtual ~Vector<T>() {}

    // Vector Outer Products
    Matrix<T> outer(const Vector<T>& v) const
    {
        Matrix<T> m((*this).size(), v.size(), 0.0);
        const Vector<T>& t(*this);

#pragma omp parallel for
        for (unsigned int i = 0; i < m.nr(); i++) {
            //#pragma omp parallel for
            for (unsigned int j = 0; j < m.nc(); j++)
                m(i, j) = t[i] * v[j];
        }

        return m;
    }

    Matrix<T> outerDiv(const Vector<T>& v) const
    {
        Matrix<T> m((*this).size(), v.size(), 0.0);
        const Vector<T>& t(*this);

#pragma omp parallel for
        for (unsigned int i = 0; i < m.nr(); i++) {
            //#pragma omp parallel for
            for (unsigned int j = 0; j < m.nc(); j++)
                m(i, j) = (v[j] != 0) ? ( t[i] / v[j] ) : ( (T)0.000000000000000001 );
        }

        return m;
    }

    // Vector | Vector Arithmetic section
    Vector<T>& func(std::function<T(T)> f)
    {

// #pragma omp parallel for
        for ( auto &i : (*this) ) {
                i = f(i);
        }

        return (*this);
    }

    Vector<T>& rand( T start=0, T end=0, int seed=0 )
    {
	T diff = end - start;

	if( seed == 0 )
		std::srand( std::time(nullptr) );
	else
		std::srand(seed);

//#pragma omp parallel for
	for ( auto &i : (*this) ) {
		i = start + ((T)std::rand()) / ( ( (T)RAND_MAX  ) / diff );
	}

        return (*this);
    }


    T dot(const Vector<T>& v) const
    {
        if (v.size() == 0)
            throw vector_algebra_error("Passed Vector is zero length for dot product");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for dot product");

        if (this->size() != v.size() || v.size() == 0)
            throw vector_algebra_error("Vector sizes do not match as needed for dot product");

        T ret = 0;
        for (int i = 0; i < this->size(); i++)
            ret += (*this)[i] * v[i];

        return ret;
    }

    Vector<T> mul(const Vector<T>& v) const
    {
        if (v.size() == 0)
            throw vector_algebra_error("Passed Vector is zero length for product");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for product");

        if (this->size() != v.size() || v.size() == 0)
            throw vector_algebra_error("Vector sizes do not match as needed for product");

        Vector<T> vret(v.size(), (T)0.0);
        for (int i = 0; i < v.size(); i++)
            vret[i] = (*this)[i] * v[i];

        return vret;
    }

    Vector<T> div(const Vector<T>& v) const
    {
        if (v.size() == 0)
            throw vector_algebra_error("Passed Vector is zero length for divsion");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for divsion");

        if (this->size() != v.size() || v.size() == 0)
            throw vector_algebra_error("Vector sizes do not match as needed for divsion");

        Vector<T> vret(v.size(), (T)0.0);
        for (int i = 0; i < v.size(); i++)
            vret[i] = (*this)[i] / v[i];

        return vret;
    }

    Vector<T> add(const Vector<T>& v) const
    {
        if (v.size() == 0)
            throw vector_algebra_error("Passed Vector is zero length for addition");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for addition");

        if (this->size() != v.size() || v.size() == 0)
            throw vector_algebra_error("Vector sizes do not match as needed for addition");

        Vector<T> vret(v.size(), (T)0.0);
        for (int i = 0; i < v.size(); i++)
            vret[i] = (*this)[i] + v[i];

        return vret;
    }

    Vector<T> sub(const Vector<T>& v) const
    {
        if (v.size() == 0)
            throw vector_algebra_error("Passed Vector is zero length for subtraction");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for subtraction");

        if (this->size() != v.size() || v.size() == 0)
            throw vector_algebra_error("Vector sizes do not match as needed for subtraction");

        Vector<T> vret(v.size(), (T)0.0);
        for (int i = 0; i < v.size(); i++)
            vret[i] = (*this)[i] - v[i];

        return vret;
    }

    // Vector | Scalar Arithmetic section

    Vector<T> mul(T s) const
    {

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for product");

        Vector<T> vret(std::vector<T>::size(), (T)0.0);
        for (int i = 0; i < std::vector<T>::size(); i++)
            vret[i] = (*this)[i] * s;

        return vret;
    }

    Vector<T> div(T s) const
    {
        if (s == (T)0.0)
            throw vector_algebra_error("Passed in scalar is zero for divsion");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for divsion");

        Vector<T> vret(std::vector<T>::size(), (T)0.0);
        for (int i = 0; i < std::vector<T>::size(); i++)
            vret[i] = (*this)[i] / s;

        return vret;
    }

    Vector<T> divR(T s) const
    {
        if (s == (T)0.0)
            throw vector_algebra_error("Passed in scalar is zero for divsion");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for divsion");

        Vector<T> vret(std::vector<T>::size(), (T)0.0);
        for (int i = 0; i < std::vector<T>::size(); i++)
            vret[i] = s / (*this)[i];

        return vret;
    }

/*
    Vector<T>& div(T s)
    {
        if (s == (T)0.0)
            throw vector_algebra_error("Passed in scalar is zero for divsion");

        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for divsion");

        //Vector<T> vret( std::vector<T>::size(), (T)0.0 );
        for (int i = 0; i < std::vector<T>::size(); i++)
            (*this)[i] /= s;

        return (*this);
    }
*/
    Vector<T> add(T s) const
    {
        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for addition");

        Vector<T> vret(std::vector<T>::size(), (T)0.0);
        for (int i = 0; i < std::vector<T>::size(); i++)
            vret[i] = (*this)[i] + s;

        return vret;
    }

    Vector<T> sub(T s) const
    {
        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for subtraction");

        Vector<T> vret(std::vector<T>::size(), (T)0.0);
        for (int i = 0; i < std::vector<T>::size(); i++)
            vret[i] = (*this)[i] - s;

        return vret;
    }

    // Vector Operator Overloading
    Vector<T> operator-() const
    {
        if (this->size() == 0)
            throw vector_algebra_error("This Vector is zero length for negation");

        Vector<T> ret((*this).mul((T)-1.0));

        return ret;
    }

    Vector<T> operator-(const Vector<T>& r) const
    {
        return sub(r);
    }

    Vector<T> operator-(const T& r) const
    {
        return sub(r);
    }

    Vector<T> operator+(const Vector<T>& r) const
    {
        return add(r);
    }

    Vector<T> operator+(const T& r) const
    {
        return add(r);
    }

    Vector<T> operator*(const Vector<T>& r) const
    {
        return mul(r);
    }

    Vector<T> operator*(const T& r) const
    {
        return mul(r);
    }

    Vector<T> operator/(const Vector<T>& r) const
    {
        return div(r);
    }
/*
    Vector<T> operator/(const T& r) const
    {
        return div(r);
    }
*/
    Vector<T>& operator/=(const T& r)
    {
        return div(r);
    }

    Vector<T>& operator=(const T& r)
    {
        std::vector<T>::clear();

        std::vector<T>::resize(r.size());

        for (int a = 0; a < r.size(); a++)
            r[a] = (*this)[a];

        return (*this);
    }

    Vector<T>& range(double begin, double end, double step = 1.0)
    {
        this->clear();

        for (double t = begin; t < end; t += step)
            this->push_back(t);

        return (*this);
    }

    Vector<T>& range(double end)
    {
        return range(0.0, end);
    }

    /*
	Vector<T> range( double begin, double end, double step = 0.0 ) const
	{
		Vector<T> v;
		
		v.range( begin, end, step );
			 
		return v;
	}
	*/

    friend Vector<T> operator-(const T& l, const Vector<T>& v)
    {
        return (-v).sub(l);
    }

    friend Vector<T> operator/(const T l, const Vector<T>& v)
    {
    	//std::cout << __func__ << std::endl;
        return ((v).divR((T)l));
    }

    friend Vector<T> operator/( const Vector<T>& v, const T l )
    {
    	//std::cout << __func__ << std::endl;
        return ((v).div((T)l));
    }

    friend Vector<T> operator*(const T& l, const Vector<T>& v)
    {
        return (v).mul(l);
    }

    friend Vector<T> operator+(const T& l, const Vector<T>& v)
    {
        return (v).add(l);
    }

    //T& operator[](unsigned int i)  { return (*this)[i]; }

    friend std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
    {

        os << "{";

        if (v.size() > 0) {
            os << v[0];
            for (int i = 1; i < v.size(); i++)
                os << "," << v[i];
        }

        os << "}";

        return os;
    }

    std::ostream& print(std::ostream& os = std::cout)
    {
        os << "|";

        if ((*this).size() > 0) {
            os << (*this)[0];
            for (int i = 1; i < (*this).size(); i++)
                os << "\t" << (*this)[i];
        }

        os << "|" << std::endl;

        return os;
    }
};

template <typename T>
class Matrix {

    Vector<Vector<T> > _m;

    unsigned int _rs;
    unsigned int _cs;

public:
    const double PI = 3.141592653589793238463;
    const double PI2 = PI * 2.0;

    Matrix<T>()
        : _rs(0)
        , _cs(0)
    {
    }

    Matrix(unsigned int rows, unsigned int cols, const T& x = (T)0.0)
    {
        resize(rows, cols, x);
    }

    Matrix(unsigned int rows, unsigned int cols, const std::initializer_list<T>& il)
    {
        set_rows_and_cols(rows, cols, il);
    }

    Matrix<T>& resize(unsigned int rows, unsigned int cols, const T& x = (T)0.0)
    {
        _m.resize(rows);

        for (unsigned int i = 0; i < _m.size(); i++)
            _m[i].resize(cols, x);

        _rs = rows;
        _cs = cols;

        return (*this);
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

    Matrix(Vector<Vector<T> > vl)
    {
        for (auto l : vl)
            _m.push_back(l);

        _rs = _m.size();
        _cs = 0;

        if (_rs > 0)
            _cs = _m[0].size();
    }

    Matrix(Vector<T> v)
    {
    	_m.push_back(v);

        _rs = _m.size();
        _cs = 0;

        if (_rs > 0)
            _cs = _m[0].size();
    }

    Matrix(const Matrix<T>& r)
    {
        (*this) = r;
    }

    virtual ~Matrix() {}

    void clear()
    {
        _m.clear();
    }

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

        _m.resize(new_rows);
        for (unsigned int i = 0; i < _m.size(); i++) {
            _m[i].resize(new_cols);
        }

        for (unsigned int i = 0; i < new_rows; i++) {
#pragma omp parallel for
            for (unsigned int j = 0; j < new_cols; j++) {
                _m[i][j] = rhs(i, j);
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
                if (_m[i][j] != rhs(i, j)) {
                    ret = false;
                    i = j = new_cols;
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
                result(i, j) = this->_m[i][j] + rhs(i, j);
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
                this->_m[i][j] += rhs(i, j);
            }
        }

        return *this;
    }

    // subtraction of this Matrix and another
    Matrix<T> operator-(T s)
    {

        Matrix result(_rs, _cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {

            for (unsigned int j = 0; j < _cs; j++) {
                result(i, j) = this->_m[i][j] - s;
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
                result(i, j) = this->_m[i][j] - rhs(i, j);
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
                result(i, j) = -(this->_m[i][j]);
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
                this->_m[i][j] -= rhs(i, j);
            }
        }

        return *this;
    }

    // Left multiplication of this Matrix and another

    Matrix<T> mul(const Matrix<T>& rhs)
    {
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                for (unsigned int k = 0; k < rhs.nr(); k++) {
                    result(i, j) += this->_m[i][k] * rhs(k, j);
                }
            }
        }

        return result;
    }

    Matrix<T> mul( T s ) const
    {
        unsigned int rs = nr();
        unsigned int cs = nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                result(i, j) = s * (*this)(i, j);
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
                result(i, j) = (*this)(i, j) / s;
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

    Matrix<T> operator*(const Matrix<T>& rhs) const
    {
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                for (unsigned int k = 0; k < rhs.nr(); k++) {
                    result(i, j) += this->_m[i][k] * rhs(k, j);
                }
            }
        }

        return result;
    }

    Matrix<T> SchurProd(const Matrix<T>& rhs) const
    {
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                result(i, j) = this->_m[i][j] * rhs(i, j);
            }
        }

        return result;
    }

    Matrix<T> SchurDiv(const Matrix<T>& rhs) const
    {
        unsigned int rs = nr();
        unsigned int cs = rhs.nc();

        //cout << rs << "x" << cs << endl;

        Matrix result(rs, cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < rs; i++) {
            for (unsigned int j = 0; j < cs; j++) {
                result(i, j) = this->_m[i][j] / rhs(i, j);
            }
        }

        return result;
    }

    // Cumulative left multiplication of this Matrix and another

    Matrix<T>& operator*=(const Matrix& rhs)
    {
        Matrix result = (*this) * rhs;
        (*this) = result;
        return *this;
    }

    Matrix<T>& func(std::function<T(T)> f)
    {
//Matrix result(_rs, _cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                (*this)(i, j) = f(this->_m[i][j]);
            }
        }

        return (*this);
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
                (*this)(i, j) = start + ((T)std::rand()) / ( ( (T)RAND_MAX  ) / diff );
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
                result(i, j) = this->_m[i][j] + rhs;
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
                result(i, j) = this->_m[i][j] * rhs;
            }
        }

        return result;
    }

    Matrix<T>& operator*=(const T& rhs)
    {
        //Matrix result(_rs, _cs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                this->_m[i][j] *= rhs;
            }
        }

        return *this;
    }

    Matrix<T>& operator=(const T& r)
    {
        _m.clear();

        _m.resize(r.nr());

        for (int a = 0; a < r.nr(); a++)
            r[a] = (*this)[a];

        return (*this);
    }

    // Multiply a Matrix with a Vector

    Vector<T> operator*(const Vector<T>& r) const
    {

        if (r.size() != nc())
        {
            std::stringstream ss;
            ss << "The vector of " << r.size() << " dimension does not match this matrix column dimension " << nc();
            throw matrix_algebra_error( ss.str().c_str() );
		}
		
        Vector<T> result(_rs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                result[i] += this->_m[i][j] * r[j];
            }
        }

        return result;
    }

    Vector<T> pmul(const Vector<T>& r) const
    {

        if (r.size() != nc())
        {
            std::stringstream ss;
            ss << "The vector of " << r.size() << " dimension does not match this matrix column dimension " << nc();
            throw matrix_algebra_error( ss.str().c_str() );
		}
		
        Vector<T> result(_rs, 0.0);
		
		/*
#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                result[i] += this->_m[i][j] * r[j];
            }
        }
		*/
		
		std::transform( std::execution::par, _m.begin(), _m.end(), result.begin(), 
				[&](auto c) -> T { return std::transform_reduce( std::execution::seq, c.begin(), c.end(), r.begin(), 0.0 ); }
			);

        return result;
    }

    Vector<T> operator+(const Vector<T>& r) const
    {
        Vector<T> result(_rs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            for (unsigned int j = 0; j < _cs; j++) {
                result[i] += this->_m[i][j] + r[j];
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
                result(j, i) = this->_m[i][j];
            }
        }

        _m.clear();
        //_m = result._m;
        (*this) = result;
        _rs = rs;
        _cs = cs;

        return *this;
    }

    // Matrix transpose
    Matrix<T> GJ_elimination( const Matrix<T> &B ) const 
    {
        const Matrix<T> &A( *this );
        if( !A.square() )	
        {
             throw matrix_algebra_error(
                 "The source matrix must be square to perform elimination");
        }
        else if( A.nr() != B.nr() )	
        {
             throw matrix_algebra_error(
                 "The augment matrix must have the same number of rows as the source matrix to perform elimination");
        }

        Matrix<T> G( A );
        Matrix<T> J( B );


	for( int c = 0; c < G.nc(); c++ )
	{			
		//T f = 1.0 / G(c,c);
		for( int r=c; r < G.nr(); r++ )
		for( int s=r; s < G.nr(); s++ )
		{
			if( r != s && G(s,c) != 0.0 )
			{
				T f = G(s,c) / G(r,c);
				//G[r] = G[r] - ( G[c]*f );
				G[s] = G[s] - ( G[r] * f );
				J[s] = J[s] - ( J[r] * f );
			}
			
		}
	}
        
	for( int c = G.nc()-1; c > 0; c-- )
	{			
		//T f = 1.0 / G(c,c);
		for( int r=c; r >= 0; r-- )
		for( int s=r; s >= 0; s-- )
		{
			if( r != s && G(s,c) != 0.0 )
			{
				T f = G(s,c) / G(r,c);
				//G[r] = G[r] - ( G[c]*f );
				G[s] = G[s] - ( G[r] * f );
				J[s] = J[s] - ( J[r] * f );
			}
			
		}
	}
        
	for( int c = 0; c < G.nc(); c++ )
	{			
		T f = 1.0 / G(c,c);
		if( f != 1.0  )
		{
			G[c] = ( G[c] * f );
			J[c] = ( J[c] * f );
		}
		
	}
        
        return J;
    }


    // Matrix transpose
    inline bool square() const 
    {
        return ( _rs == _cs );        
    }

    // Obtain a Vector of the diagonal elements

    Vector<T> diag_vec()
    {
        Vector<T> result(_rs, 0.0);

#pragma omp parallel for
        for (unsigned int i = 0; i < _rs; i++) {
            result[i] = this->_m[i][i];
        }

        return result;
    }

    void set_cols_and_rows(int c, int r, std::initializer_list<T>& l)
    {
        resize(r, c, (T)0.0);
        Matrix<T>& M = (*this);

        auto k = l.begin();
        for (int i = 0; i < M.nc(); i++)
            for (int j = 0; j < M.nr(); j++, k++)
                M[i][j] = *(k);
    }

    void set_rows_and_cols(int r, int c, std::initializer_list<T> l)
    {
        resize(r, c, (T)0.0);
        Matrix<T>& M = (*this);

        auto k = l.begin();
        for (int i = 0; i < M.nr(); i++)
            for (int j = 0; j < M.nc(); j++, k++)
                M[i][j] = *(k);
    }

    // Access the individual elements
    T& operator()(const unsigned int& row, const unsigned int& col)
    {
        return this->_m[row][col];
    }

    // Access the individual elements (const)
    const T& operator()(const unsigned int& row, const unsigned int& col) const
    {
        return this->_m[row][col];
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


    Vector<T> operator[](unsigned int i) const
    {
        return _m[i];
    }

    Vector<T>& operator[](unsigned int i)
    {
        return _m[i];
    }

    friend Matrix<T> operator/(const T l, const Matrix<T>& v)
    {
    	//std::cout << __func__ << std::endl;
        return ((v).divR((T)l));
    }

    friend Matrix<T> operator/( const Matrix<T>& v, const T l )
    {
    	//std::cout << __func__ << std::endl;
        return ((v).div((T)l));
    }



    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
    {

        os << "{";

        if (m.nr() > 0) {
            os << m[0];
            for (int i = 1; i < m.nr(); i++)
                os << "," << m[i];
        }
        else
            os << "{}";

        os << "}";

        return os;
    }

    std::ostream& print(std::ostream& os = std::cout)
    {
        if (nr() > 0) {
            //os << "|";
            //os << _m[0];
            for (int i = 0; i < nr(); i++)
                _m[i].print(os);
            //os << "|";
        }

        return os;
    }


};

#endif

