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
/*  $Id: Matrix.H,v 1.4 1993/11/18 07:29:25 jak Exp $
 */
//  History:
/*  $Log: Matrix.H,v $
/*  Revision 1.4  1993/11/18 07:29:25  jak
/*  Added alot of increased functionality, including support for
/*  non-zero aligned matrices.  This supports dealing with
/*  arbitrary matrix partitions.  Also, LU decompositions are
/*  stored with the matrices the derived from, and are recovered
/*  rather than re-computed if a matrix is re-used.   -jak
/*
 * Revision 1.3  1993/11/15  20:29:41  jak
 * Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
 **/
// =====================================

#ifndef _Matrix_H
#define _Matrix_H
static char rcsid_Matrix_H[] =  "$Id: Matrix.H,v 1.4 1993/11/18 07:29:25 jak Exp $";

#include <iostream.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
#include "Linked_List_Template.H"

const float TINY = 1.0e-20; 	// some small number not equal to zero

class Matrix;
class Composition;
class LU_Decomposition;
class IntArray;

//
// Decomposition Components (keys for linked list)
//
#define LU_DECOMP  0
#define P_OF_LU    1
#define L_OF_LU    2
#define U_OF_LU    3


// template <class NUMTYPE>
// Range Template
// ===
// User must define:
//    <NUMTYPE>::operator += ()
//    <NUMTYPE>::operator > ()
// note: these should be already defined for 
// builtin numeric types (e.g. int, char, float)
template <class NUMTYPE>
class Range {
    NUMTYPE start;
    NUMTYPE end;
    NUMTYPE step;
    NUMTYPE current;

public:
    Range(NUMTYPE s, NUMTYPE e, NUMTYPE st): 
         start(s), end(e), step(st), current(s) {};
    ~Range(){};

    inline NUMTYPE start( void ) const { return start; }
    inline void setStart( NUMTYPE s ) const { start = s; }

    inline NUMTYPE end( void ) const { return end; }
    inline void setEnd( NUMTYPE e ) const { end = e; }

    inline NUMTYPE step( void ) const{ return step; }
    inline void setStep( NUMTYPE st ) const{ step = st; }

    inline NUMTYPE value() const{ return current; }
    inline NUMTYPE next(){ 
        register NUMTYPE rtn;
        rtn = current;
        current += step;
        if (current > end) current = start;
        return rtn; 
    };
    inline NUMTYPE reset(){ 
        current = start; 
        return current; 
    };
};

#define MAXNAMELEN 31

class Matrix 
{
protected:
    char name[ MAXNAMELEN ];
    int first_row, first_col;
    unsigned int number_of_rows, number_of_cols;
    float **m;
    Linked_List< Composition > decomp_list;
    void allocate(void);
    void deallocate(void);

public:

// Constructors and Destructors
    Matrix( unsigned int = 0, unsigned int = 0, const char* = "");
    Matrix( unsigned int, unsigned int, int, int, const char* = "");
    Matrix( const Matrix& );
    ~Matrix( void );

// Assignment Operations
    Matrix&  operator = (const Matrix &);

// Type Conversion Operator

//  Submatrix Access Methods
    Matrix operator() ( int, int, int, int ) const;
    Matrix&  shift_to ( int, int );

// Matrix Information methods
    inline char* getName( void ) const { return name; };
    inline void setName( const char* aname ) {
        strncpy(name, aname, MAXNAMELEN-1); 
    };
    inline unsigned int rows( void ) const { return number_of_rows;    };
    inline int firstrow( void ) const { return first_row;    };
    inline unsigned int cols( void ) const { return number_of_cols; };
    inline int firstcol( void ) const { return first_col; };
    inline float* operator [] ( int i ) const { return m[ i ]; };
    inline int has_element( int i, int j) const {
        if (     (i <  first_row) 
              || (i >= first_row + number_of_rows ) 
              || (j <  first_col) 
              || (j >= first_col + number_of_cols ))
            return 0;
        else
            return 1;
    };
    inline int is_empty() const {
        return m?0:1;
    };

// Matrix Operators
    Matrix & operator += ( const Matrix& );
    Matrix & operator += ( float  );
    Matrix & operator -= ( const Matrix& );
    Matrix & operator -= ( float );
    Matrix & operator *= ( const Matrix& );
    Matrix & operator *= ( float );
    Matrix & operator /= ( const Matrix& );
    inline Matrix & operator /= ( float scalar ){
        return (this->operator *= (1.0/scalar));
    };
    Matrix & operator &= ( const Matrix & );  // elementwise multiply
    Matrix & operator %= ( const Matrix & );  // elementwise divide

    Matrix & apply ( float(*)( float ) );
    Matrix & apply ( float(*)( float, float ), const Matrix& );

// Matrix Norms
// default Euclidean norm p = 2
// p = 1 -> 1-norm, p=0 -> infinity norm
    friend float norm ( const Matrix& , float = 2.0 ); 

// Comparison Operations
    friend int operator == ( const Matrix&, const Matrix& );

// Matrix - Matrix Operations
    friend Matrix operator + ( const Matrix&, const Matrix& );
    friend Matrix operator - ( const Matrix&, const Matrix& );
    friend Matrix operator * ( const Matrix&, const Matrix& );
    friend Matrix operator / ( const Matrix&, const Matrix& );

// Scalar - Matrix Operations 
    friend Matrix operator * (float, const Matrix&);
    friend inline Matrix operator * (const Matrix &matA, float scalar){
        return (scalar * matA);
    }

    friend Matrix operator + ( const Matrix &, float );
    friend inline Matrix operator + (float scalar, const Matrix &matA){
        return ( matA + scalar );
    }
    friend inline Matrix operator - (const Matrix &matA, float scalar){
        return ( matA + ( - scalar ) );
    }
    friend inline Matrix operator - (float scalar, const Matrix &matA){
        return -( matA + scalar );
    }
    friend inline Matrix operator - (const Matrix &matA ){
        return (-1.0 * matA);
    }
    friend inline Matrix operator / (const Matrix &matA, float scalar){
        return ( ( 1.0/scalar ) * matA );
    }
    friend inline Matrix operator / (float scalar, const Matrix &matA){
        return ( scalar * inverse( matA ) );
    }

// Transformations
    friend Matrix transpose( const Matrix &); 
    friend Matrix inverse( const Matrix &);

// Matrix elementwise functions (for signal processing)
    friend Matrix operator & ( const Matrix &, const Matrix &); // elementwise multiply
    friend Matrix operator % ( const Matrix &, const Matrix &); // elementwise divide

    friend Matrix map ( float(*)( float ), const Matrix& );
    friend Matrix map ( float(*)( float, float ), const Matrix&, const Matrix& );

// Input Output
    friend ostream & operator << (ostream &, const Matrix&);
    friend istream & operator >> (istream &, const Matrix&);

// Decompositions
    void         setDecompose ( Composition *, int   );
    Composition *getDecompose ( int   ) const ;
    void         delDecompose ( int   );

// Special Matrices
    friend Matrix Ones( int, int );
    friend Matrix Identity( int );
    friend Matrix UpperTriangle( int );
    friend Matrix LowerTriangle( int );
};

class Composition {
public:
    Composition( void ){};
    virtual ~Composition( void ){};

    virtual Matrix operator[] ( int  key ) const  = 0 ;
    virtual Matrix solve_for ( const Matrix & )   = 0 ;
};

class IntArray {
    int first_index;
    unsigned int length;
    int *array;

    void allocate(void);
    void deallocate(void);

public:

// Constructors and Destructors
    IntArray( unsigned int = 0 );
    IntArray( int, unsigned int);
    IntArray( const IntArray & );
    ~IntArray( void );

// Assignment Operations
    IntArray&  operator = (const IntArray &);

//  Sub array Access Methods
    IntArray operator() ( int, unsigned int ) const;
    IntArray&  shift_to ( int );

    inline unsigned int getLength( void ) const { return length;    };
    inline int first( void ) const { return first_index;    };

    inline int& operator [] ( int i ) const { return array[ i ]; };
    inline int has_element( int i ) const {
        if ( (i <  first_index) || (i >= first_index + length ))
            return 0;
        else
            return 1;
    };
    inline int is_empty() const {
        return array?0:1;
    };

// Input Output
    friend ostream & operator << (ostream &, const IntArray&);
    friend istream & operator >> (istream &, const IntArray&);
};

class LU_Decomposition : public Composition
{
public:
    IntArray row_permute;
    Matrix   lumat;

    LU_Decomposition( Matrix& = Matrix() );
    LU_Decomposition( const LU_Decomposition& );
    virtual ~LU_Decomposition( void );

// Assignment Operations
    LU_Decomposition&  operator = (const LU_Decomposition &);

    inline int is_empty() const {
        return lumat.is_empty();
    };
    inline LU_Decomposition&  shift_to (int fr, int fc){
        row_permute.shift_to( fr );
        lumat.shift_to( fr, fc );
        return *this;
    };

    virtual Matrix operator[] ( int  ) const;
    virtual Matrix solve_for ( const Matrix & );

};


#endif // _Matrix_H

