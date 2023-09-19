# la.h - a C++ template class set for supporting the Vector and Matrix operations of Linear Algebra

I never really "learned" linear algebra.  So I have decided to learn it by writing a C++ template class that supports float, double, and std::complex\<T\> numeric types.  As I learn more the class functionality and features will grow and improve.  I am open to comments and critiques as they are great tools for learning.  And a great way to get a fixed or improved version of the class.


# Files

The only file needed is **la.h**.  The other files are for test cases and usage examples.

Don't forget to set **export OMP_NUM_THREADS=8** (or some other cpu core count number)

## Current features and capabilities

### Vector arithmetic and operations
The **Vector\<T\>** class supports the following:
* adding
* subtracting
* multiplying
* dividing
* outer product
* dot product
* as well as operations with scalars values


### Matrix arithmetic and operations
The **Martix\<T\>** class supports the following:
* adding
* subtracting
* multiplying
* dividing
* determinant
* adjugate
* transpose
* inverse
* Gauss-Jordan elimination
* as well as operations with Vector\<T\> types and scalars values



