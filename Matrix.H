//
//  Interface to 2-D Matrix Class
//

#ifndef _Matrix_H
#define _Matrix_H

#include <iostream.h>

#define TINY 	1.0e-20		// some small number not equal to zero

class Matrix 
{
public:
    unsigned short myRows, myColumns;
    float **m;

// Constructors and Destructors
    Matrix( unsigned short = 1, unsigned short = 1);
    Matrix( Matrix & );
    ~Matrix( void );

// Type Conversion Operator

// Matrix Information methods
    inline unsigned short rows( void ){ return myRows;    }
    inline unsigned short cols( void ){ return myColumns; }
    inline float* operator [] ( unsigned short i ){ return m[ i ]; }

// Matrix - Matrix Operations
    friend Matrix& operator + (Matrix &, Matrix &);
    friend Matrix& operator - (Matrix &, Matrix &);
    friend Matrix& operator * (Matrix &, Matrix &);
    friend Matrix& operator / (Matrix &, Matrix &);

// Scalar - Matrix Operations 
    friend Matrix& operator * (float, Matrix &);
    friend inline Matrix& operator * (Matrix &matA, float scalar){
        return (scalar * matA);
    }

    friend Matrix& operator + ( Matrix &, float );
    friend inline Matrix& operator + ( float scalar, Matrix &matA){
        return ( matA + scalar );
    }
    friend inline Matrix& operator - ( Matrix &matA, float scalar){
        return ( matA + ( - scalar ) );
    }
    friend inline Matrix& operator - ( float scalar, Matrix &matA){
        return -( matA + scalar );
    }
    friend inline Matrix& operator - (Matrix &matA ){
        return (-1.0 * matA);
    }
    friend inline Matrix& operator / (Matrix &matA, float scalar){
        return ( ( 1.0/scalar ) * matA );
    }
    friend inline Matrix& operator / (float scalar, Matrix &matA){
        return ( scalar * inverse( matA ) );
    }

// Assignment Operations
    Matrix&  operator = (Matrix &);

// Transformations
    friend Matrix& transpose( Matrix &); 
    friend Matrix& inverse( Matrix &);

// Matrix elementwise functional mappings
    friend Matrix& operator ^ ( Matrix &, Matrix &);  // elementwise multiply
    friend Matrix& operator % ( Matrix &, Matrix &);  // elementwise divide
#ifndef GCC
    friend Matrix& map ( float (*)( float, float ), Matrix &, Matrix &);
    friend Matrix& map ( float (*)( float ), Matrix &);
#else
    friend Matrix& map ( float (*)(), Matrix &, Matrix &);
    friend Matrix& map ( float (*)(), Matrix &);
#endif

// Input Output
    friend ostream& operator << (ostream &, Matrix &);
    friend istream& operator >> (istream &, Matrix &);
};

class LU_Decomp
{
protected:
    unsigned short *row_permute;
    unsigned short size;
    float permute_parity;  // even # of row permutations = +1.0, odd = -1.0
    float **lumat;

public:

    LU_Decomp( Matrix & );
    ~LU_Decomp();

    Matrix& L( void );
    Matrix& U( void );
    Matrix& solve_for( Matrix & );

};

#endif // _Matrix_H

