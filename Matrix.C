// 
// Implementation of 2D Matrix Class
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
/*  $Id: Matrix.C,v 1.3 1993/11/15 20:29:39 jak Exp $
 */
//  History:
/*  $Log: Matrix.C,v $
/*  Revision 1.3  1993/11/15 20:29:39  jak
/*  Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
/**/
// =====================================

static char rcsid_C[] =  "$Id: Matrix.C,v 1.3 1993/11/15 20:29:39 jak Exp $";


#ifdef LIBGpp
#include <new.h>
#endif

#include <stdlib.h>
#include <math.h>
#include "Matrix.H"

void Abort( char *s ){
    cerr << "Matrix::" << s;
    exit( -1 );
};

Matrix::Matrix(unsigned short r, unsigned short c)
{
    int i;
    float *newmat;

    myRows = r;
    myColumns = c;

    newmat = (float *) new float[ r * c ];
    m = (float **) new float*[ sizeof(float*) * r ];

    if (newmat == (float  *)0 ) cerr << "newmat == 0 !\n";	
    if (m      == (float **)0 ) cerr << "m == 0 !\n";

    for (i = 0; i < (int)r; i++){
        m[i] = &(newmat[i*c]);
    }

    for( i=0; i < (int)r*(int)c; i++){
        *newmat++ = 0.0;
    }
};

Matrix::Matrix(const Matrix &matA)    // Copy Constructor
{
    int i;
    float *newmat, *ptr;

    myRows = matA.rows();
    myColumns = matA.cols();

    newmat = (float *) new float[ myRows * myColumns ];
    m = (float **) new float*[ sizeof(float*) * myRows];

    if (newmat == (float  *)0 ) cerr << "newmat == 0 !\n";	
    if (m      == (float **)0 ) cerr << "m == 0 !\n";
	
    for (i = 0; i< (int)myRows; i++){
        m[i] = &(newmat[ i * myColumns ]);
    }

    ptr = matA[0];
    for( i=0; i < (int)myRows*(int)myColumns; i++){
        *newmat++ = *ptr++;
    }
};

Matrix::~Matrix( void )
{
    delete  *m;
    delete  m;
};

Matrix& Matrix::operator=(const Matrix &matB)   // performs a copy
{
    register short int r,c;
    register float *newmat, *ptr;

    if ((matB.rows() != myRows) || (matB.cols() != myColumns)){
        delete *m;
        delete m;
        myRows    = matB.rows();
        myColumns = matB.cols();
		newmat = (float *) new float[ myRows * myColumns ];
		m = (float **) new float*[ sizeof(float*) * myRows];
		for (r = 0; r< (int)myRows; r++){
			m[r] = &(newmat[ r * myColumns ]);
		}
    }

    ptr = matB[0];
    for( r=0; r < (int)myRows*(int)myColumns; r++){
        *newmat++ = *ptr++;
    }

    return (*this);
};

Matrix& transpose( const Matrix &matA )
{
    register float *ptr_r;
    register unsigned short r,c;
    unsigned short Rows, Cols;
    Matrix  *result;

    result = new Matrix( Rows = matA.cols(), Cols = matA.rows() );

    for( r=0; r< Rows; r++ ){
        ptr_r = (*result)[r];
        for( c=0; c<Cols; c++ ){
            *ptr_r++ = matA[c][r];
        }
    }

    return (*result);
};

Matrix& operator+( const Matrix &matA, const Matrix &matB )
{
    Matrix *result;
    register short int r,c, Rows, Cols;
    register float *ptr_r, *ptr_a, *ptr_b;

    if ((Rows = matA.rows()) != matB.rows()) 
        Abort("Incompatible Row Dimensions");
    if ((Cols = matA.cols()) != matB.cols()) 
        Abort("Incompatible Col Dimensions");

    result = new Matrix( Rows, Cols );
    
    for( r=0; r<Rows; r++ ){  
        ptr_r = (*result)[r];
        ptr_a = matA[r];
        ptr_b = matB[r];
        for( c=0; c<Cols ; c++)
            *ptr_r++ = *ptr_a++ + *ptr_b++;
    }

    return (*result);
};

Matrix& operator-( const Matrix &matA, const Matrix &matB )
{
    Matrix *result;
    register short int r,c, Rows, Cols;
    register float *ptr_r, *ptr_a, *ptr_b;

    if ((Rows = matA.rows()) != matB.rows()) 
        Abort("Incompatible Row Dimensions");
    if ((Cols = matA.cols()) != matB.cols()) 
        Abort("Incompatible Col Dimensions");

    result = new Matrix( Rows, Cols );
    
    for( r=0; r<Rows; r++ ){  
        ptr_r = (*result)[r];
        ptr_a = matA[r];
        ptr_b = matB[r];
        for( c=0; c<Cols ; c++)
            *ptr_r++ = *ptr_a++ - *ptr_b++;
    }

    return (*result);
};

Matrix& operator*( const Matrix &matA, const Matrix &matB )
{    
    Matrix *result;
    register unsigned short r,c,q;
    unsigned short Rows_A, Cols_A, Rows_B, Cols_B;
    register float *ptr_r, *ptr_a, *ptr_b;

    Cols_A = matA.cols();  Cols_B = matB.cols();
    Rows_A = matA.rows();  Rows_B = matB.rows();

    result = new Matrix( Rows_A, Cols_B );

// MxQ times QxN
    if ( Cols_A == Rows_B ){
        Matrix B_bycols( transpose( matB ) );

        for (r=0; r< Rows_A; r++){
            ptr_r = (*result)[r];
            for (c=0; c<Cols_B; c++, ptr_r++){
                ptr_a = matA[r];
                ptr_b = (B_bycols)[c];
                *ptr_r = 0.0;
                for (q=0; q<Cols_A; q++)
                    *ptr_r += *ptr_a++ * *ptr_b++; 
            }
        }
        return (*result);
     }

// 1x1 times NxM
    else if ((Cols_A == 1) && (Rows_A == 1)){
        return matA[0][0] * matB;
    } else if ((Rows_B == 1) && (Cols_B == 1)){
        return matA * matB[0][0];
    }

// Error 
    else {
        Abort("Incompatible Matrix Dimensions");
    }

    return (*result * 0.0);
};

Matrix& operator*( float scalar, const Matrix &matA )
{
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] * scalar;

    return (*result);
};

Matrix& operator+( const Matrix &matA, float scalar){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] + scalar;

    return (*result);
};

Matrix& map( float (*f)(float), const Matrix &matA){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = (*f)( matA[r][c] );

    return (*result);
};

Matrix& map( float (*f)(float,float), const Matrix &matA, const Matrix &matB){

    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = (*f)( matA[r][c] , matB[r][c]);

    return (*result);
};

Matrix& operator^( const Matrix &matA, const Matrix &matB){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] * matB[r][c];

    return (*result);
};

Matrix& operator%( const Matrix &matA, const Matrix &matB){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] / matB[r][c];

    return (*result);
};

ostream & operator << (ostream &cbuf, const Matrix &matA)
{
    register float *ptr;
    register unsigned short r,c;

    for (r=0; r<matA.rows(); r++){
        ptr = matA[r];
        for (c=0; c< matA.cols(); c++){
            cbuf << *ptr++ << "\t";
        }
        cbuf << "\n";
    }
    cbuf << "\n";
    return cbuf;
};

istream & operator >> (istream &cbuf, const Matrix &matA)
{
    register float *ptr;
    register unsigned short r,c;

    for (r=0; r<matA.rows(); r++){
        ptr = matA[r];
        for (c=0; c<matA.cols(); c++){
            cbuf >> *ptr++;
        }
    }
    return cbuf;
};

Matrix& operator / (const Matrix &matA, const Matrix &matB)
{
  // solve X = matA / matB
  // means :
  //  matB * X = matA
  //         X = inverse( matB ) * matA
  //  NOT
  //  X * matB = matA
  //         X = matA * inverse( matB );
  // if the second case occurs, then the alg will try to
  // recognize it, and transpose everything to fix it.
  // this is not guaranteed to succeed - but I think it will
  //
    if ( matA.rows() < matA.cols() ){
        LU_Decomp LU( transpose(matB) );
        return  transpose( LU.solve_for( transpose( matA ) ) );
    } else {
        LU_Decomp LU( matB );
        return  LU.solve_for( matA );
    }
};

Matrix& inverse( const Matrix &matA )
{
    register unsigned short r;
    unsigned short size;

    if ((size = matA.rows()) != matA.cols())
        Abort("Singular - Inverse requires square matrix!");
    // else
    {    
        LU_Decomp LU( matA );
        Matrix identity( size, size );
        for( r=0; r< size; r++){
            identity[r][r] = 1.0;
        }
        return LU.solve_for( identity );
    }    
};

LU_Decomp::LU_Decomp( const Matrix &matA )
{
// This algorithm is a distinct implementation of an algorithm found
// in Numerical Recipes in C.  The ideas expressed as code here, were 
// derived from that work, which this programmer highly recommends.
// This routine is based on Crout's method with implicit partial pivoting 
// as described by Numerical Recipes in C for the ludcmp() routine.

    register unsigned short row, col, k;
    float *row_scaling;
    float *m,*ptr;

// LU routine
//
    if ( matA.rows() != matA.cols() ) 
        Abort("LU_Decomp::LU_Decomp - Square Matrix Required!");

// Constructor
//
    size = matA.rows();
    m = (float *) new (float[ size * size ]);
    lumat = (float **) new float*[ sizeof(float*) * size];
    
    for ( row=0; row<size; row++){
        lumat[row] = &( m[ row*size ] );
    }

    row_permute = new unsigned short[ size ];
    if (row_permute  == (unsigned short *)0 ) cerr << "row_permute == 0 !\n";
	
    for ( row=0; row< size; row++){
        row_permute[ row ] = row;
    }

    permute_parity = 1.0;

    ptr = matA[0];
    for( k=0; (int)k<(int)size*(int)size; k++){
        *m++ = *ptr++;
    }
//
// End Construction

   row_scaling = new float[ size ];

// Find the Maximum element value in each row and record a scaling factor
//
    for( row=0; row<size; row++){
        float max, tempval;

        max = 0.0;
        for( col=0; col<size; col++ ){
            tempval = fabs( lumat[row][col] );
            if(tempval > max) max = tempval;
        }
        if (max = 0.0) cerr << "Matrix::LU_Decomp - Singular Matrix\n";
        row_scaling[ row ] = 1.0 / max;
    }

// Crout's Method of Gaussian Elimination
//
    for( col=0; col<size; col++){
        float sum, max, tempval;
        unsigned int max_index;

        for( row=0; row<col; row++) {
            sum = lumat[row][col];
            for( k=0; k<row; k++)
                sum -= lumat[row][k] * lumat[k][col];
            lumat[row][col] = sum;
        }
        max = 0.0;
        max_index = col;
        for( row=col; row<size; row++){
            sum = lumat[row][col];
            for( k=0; k<col; k++)
                sum -= lumat[row][k] * lumat[k][col];
            lumat[row][col] = sum;
            if((tempval = row_scaling[row] * fabs( sum )) >= max ){
                max_index = row;
                max = tempval;
            }
        }
        if ( col != max_index ){
            for( k=0; k<size; k++ ){
                tempval = lumat[max_index][k];
                lumat[max_index][k] = lumat[col][k];
                lumat[col][k] = tempval;
            }
            permute_parity = -permute_parity;
            row_scaling[ max_index ] = row_scaling[ col ];
        }
        row_permute[ col ] = max_index ;
        if ( lumat[col][col] == 0.0 ) lumat[col][col] = TINY;
        if ( col != size-1 ){
            tempval = 1.0 / (lumat[col][col]);
            for( row=col+1; row<size; row++) lumat[row][col] *= tempval;
        }
    }

// End LU routine
// 

    delete row_scaling;

};

LU_Decomp::~LU_Decomp()
{
    delete  *lumat;
    delete  lumat;
    delete  row_permute;
};

Matrix& LU_Decomp::L() const
{    
    Matrix *result;
    register short int row, col;

    result = new Matrix( size, size );

    for( row=0; (int)row<(int)size; row++ )
        for( col=0; (int)col<(int)size; col++ )
            if( row > col )
                (*result)[row][col] = lumat[row][col];
            else if ( row == col )
                (*result)[row][col] = 1.0;
            else // ( row < col )
                (*result)[row][col] = 0.0;
//
// Return the L matrix in its permuted row order.
//
    for( row=0; (int)row<(int)size; row++){
        float tempval;
        if (row != row_permute[ row ])
            for( col=0; (int)col<(int)size; col++){
                tempval = (*result)[row][col];
                (*result)[row][col] = (*result)[row_permute[ row ]][col];
                (*result)[row_permute[ row ]][col] = tempval;
            }
    }

    return (*result);
};

Matrix& LU_Decomp::U() const
{
    Matrix *result;
    register short int row, col;

    result = new Matrix( size, size );

    for( row=0; (int)row<(int)size; row++ )
        for( col=0; (int)col<(int)size; col++ )
            if( row <= col )
                (*result)[row][col] = lumat[row][col];
            else  // ( row > col )
                (*result)[row][col] = 0.0;
//   
// Only the L matrix is returned in its permuted row order.
//
    return (*result);
};

Matrix& LU_Decomp::solve_for( const Matrix &matA )
{
    Matrix *result;
    register short row, col;
    unsigned short Rows, Cols;

    Rows = matA.rows();
    Cols = matA.cols();
    result = new Matrix( matA );

    for( col=0; (int)col<(int)Cols; col++ ){
        int index;
        register unsigned short k;
        unsigned short row_p;
        float sum;
//
// forward substitution
//
        index = -1;
        for( row=0; (int)row<(int)Rows; row++ ){
            row_p = row_permute[row];
// --
            sum = (*result)[ row_p ][ col ];
            (*result)[ row_p ][ col ] = (*result)[ row ][ col ];
            if (index >= 0)
                for( k=index; (int)k<(int)row;  k++) 
                    sum -= lumat[row][k] * (*result)[k][col];
            else if (sum != 0.0)
                index = row;
            (*result)[row][col] = sum;
        }
//
// backsubstitution
//
        for( row=Rows-1; (int)row>=0; row--){
            sum = (*result)[row][col];
            for( k=row+1; k<Rows; k++){
                sum -= lumat[row][k] * (*result)[k][col];
            }
            (*result)[row][col] = sum / lumat[row][row];
        }
    }

    return (*result);
};
