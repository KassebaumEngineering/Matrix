//
//  Interface to 2D Matrix Class
//
//  Author: John Kassebaum
//  
// 
//  This file and it's contents are my unique creation, and you
//  can't have it.  I retain all copy rights.  I may from time to time 
//  grant license to others to use it, but I retain ownership of
//  MY software.  I hereby give credit my classes in Linear Algenbra and 
//  and materials which I have read (especially Numerical Recipes in C) for 
//  concepts and direction. 
//
//  Revision:
/*  $Id: Matrix.H,v 1.3 1993/11/15 20:29:41 jak Exp $
 */
//  History:
/*  $Log: Matrix.H,v $
/*  Revision 1.3  1993/11/15 20:29:41  jak
/*  Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
/**/
// =====================================

#ifndef _Matrix_H
#define _Matrix_H
static char rcsid_H[] =  "$Id: Matrix.H,v 1.3 1993/11/15 20:29:41 jak Exp $";

#include <iostream.h>
#include <sys/types.h>
#include <sys/stat.h>

const float TINY = 1.0e-20; 	// some small number not equal to zero

class Matrix 
{
public:
    unsigned short myRows, myColumns;
    float **m;

// Constructors and Destructors
    Matrix( unsigned short = 1, unsigned short = 1);
    Matrix( const Matrix& );
    ~Matrix( void );

// Type Conversion Operator

// Matrix Information methods
    inline unsigned short rows( void ) const { return myRows;    }
    inline unsigned short cols( void ) const { return myColumns; }
    inline float* operator [] ( unsigned short i ) const { return m[ i ]; }

// Matrix - Matrix Operations
    friend Matrix& operator + ( const Matrix&, const Matrix& );
    friend Matrix& operator - ( const Matrix&, const Matrix& );
    friend Matrix& operator * ( const Matrix&, const Matrix& );
    friend Matrix& operator / ( const Matrix&, const Matrix& );

// Scalar - Matrix Operations 
    friend Matrix& operator * (float, const Matrix&);
    friend inline Matrix& operator * (const Matrix &matA, float scalar){
        return (scalar * matA);
    }

    friend Matrix& operator + ( const Matrix &, float );
    friend inline Matrix& operator + (float scalar, const Matrix &matA){
        return ( matA + scalar );
    }
    friend inline Matrix& operator - (const Matrix &matA, float scalar){
        return ( matA + ( - scalar ) );
    }
    friend inline Matrix& operator - (float scalar, const Matrix &matA){
        return -( matA + scalar );
    }
    friend inline Matrix& operator - (const Matrix &matA ){
        return (-1.0 * matA);
    }
    friend inline Matrix& operator / (const Matrix &matA, float scalar){
        return ( ( 1.0/scalar ) * matA );
    }
    friend inline Matrix& operator / (float scalar, const Matrix &matA){
        return ( scalar * inverse( matA ) );
    }

// Assignment Operations
    Matrix&  operator = (const Matrix &);

// Transformations
    friend Matrix& transpose( const Matrix &); 
    friend Matrix& inverse( const Matrix &);

// Matrix elementwise functional mappings
    friend Matrix& operator ^ ( const Matrix &, const Matrix &);  // elementwise multiply
    friend Matrix& operator % ( const Matrix &, const Matrix &);  // elementwise divide

    friend Matrix& map ( float(*)( float ), const Matrix& );
    friend Matrix& map ( float(*)( float, float ), const Matrix&, const Matrix& );

// Input Output
    friend ostream & operator << (ostream &, const Matrix&);
    friend istream & operator >> (istream &, const Matrix&);
};

class LU_Decomp
{
protected:
    unsigned short *row_permute;
    unsigned short size;
    float permute_parity;  // even # of row permutations = +1.0, odd = -1.0
    float **lumat;

public:

    LU_Decomp( const Matrix & );
    ~LU_Decomp();

    Matrix& L( void ) const;
    Matrix& U( void ) const;
    Matrix& solve_for( const Matrix & );

};

#endif // _Matrix_H

