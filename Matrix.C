// 
// Implementation of 2-D MAtrix Class
//

#include <stdlib.h>
#include <math.h>
#include "Matrix.H"

void Abort( char *s ){
    cerr << "Matrix::" << s;
    exit( -1 );
};

Matrix::Matrix(unsigned short r, unsigned short c)
{
    unsigned short i;
    float *newmat;

    myRows = r;
    myColumns = c;

    newmat = (float *) new float[ r ][ c ];
    m = new float (*)[ r ];

    if (newmat == (float  *)0 ) cerr << "newmat == 0 !\n";	
    if (m      == (float **)0 ) cerr << "m == 0 !\n";

    for (i = 0; i< r; i++){
        m[i] = &(newmat[i*c]);
    }

    for( i=0; i<myRows*myColumns; i++){
        *newmat++ = 0.0;
    }
};

Matrix::Matrix(Matrix &matA)    // Copy Constructor
{
    unsigned short i, r, c;
    float *newmat, *ptr;

    myRows = matA.rows();
    myColumns = matA.cols();

    newmat = (float *) new float[ myRows ][ myColumns ];
    m = new float (*)[ myRows ];

    if (newmat == (float  *)0 ) cerr << "newmat == 0 !\n";	
    if (m      == (float **)0 ) cerr << "m == 0 !\n";
	
    for (i = 0; i< myRows; i++){
        m[i] = &(newmat[ i * myColumns ]);
    }

    ptr = matA[0];
    for( i=0; i<myRows*myColumns; i++){
        *newmat++ = *ptr++;
    }
};

Matrix::~Matrix( void )
{
    unsigned short i;

    delete  *m;
    delete  m;
};

Matrix operator=(Matrix &matA, Matrix &matB)   // performs a copy
{
    Matrix *result;
    register short int r,c, Rows, Cols;
    float *ptr_A, *ptr_B;

    result = &matA;
    if ((matB.rows() != matA.rows()) || (matB.cols() != matA.cols())){
        delete result;
        result = new Matrix( matB.rows(), matB.cols() );
    }
    for( r=0; r<matA.rows(); r++){
        ptr_A = (*result)[r];
        ptr_B = matB[r];
        for( c=0; c<matB.cols(); c++)
            *ptr_A++ = *ptr_B++;
    }

    return (*result);
};

Matrix *operator=(Matrix *matA, Matrix &matB)  // returns a reference
{
    return (matA = &matB);
};

Matrix transpose( Matrix &matA )
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

Matrix operator+( Matrix &matA, Matrix &matB )
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

Matrix operator-( Matrix &matA, Matrix &matB )
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

Matrix operator*( Matrix &matA, Matrix &matB )
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

Matrix operator*( float scalar, Matrix &matA )
{
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] * scalar;

    return (*result);
};

Matrix operator+( Matrix &matA, float scalar){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] + scalar;

    return (*result);
};

Matrix map( float (*f)( ), Matrix &matA){
    Matrix *result;
    register short int r,c, Rows, Cols;
    float (*func)( float ) = (float (*)( float ))f;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = func( matA[r][c] );

    return (*result);
};

Matrix map( float (*f)( ), Matrix &matA, Matrix &matB){
    Matrix *result;
    register short int r,c, Rows, Cols;
    float (*func)( float, float ) = (float (*)( float, float ))f;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = func( matA[r][c] , matB[r][c]);

    return (*result);
};

Matrix operator^( Matrix &matA, Matrix &matB){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] * matB[r][c];

    return (*result);
};

Matrix operator%( Matrix &matA, Matrix &matB){
    Matrix *result;
    register short int r,c, Rows, Cols;

    result = new Matrix( Rows = matA.rows(), Cols = matA.cols() );

    for( r=0; r<Rows; r++)
        for( c=0; c<Cols; c++)
            (*result)[r][c] = matA[r][c] / matB[r][c];

    return (*result);
};

ostream& operator << (ostream &cbuf, Matrix &matA)
{
    register float *ptr;
    register unsigned short r,c;

    for (r=0; r<matA.rows(); r++){
        ptr = matA[r];
        for (c=0; c<matA.cols(); c++){
            cbuf << *ptr++;   cbuf << "\t";
        }
        cbuf << "\n";
    }
    cbuf << "\n";
    return cbuf;
};

istream& operator >> (istream &cbuf, Matrix &matA)
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

Matrix operator / (Matrix &matA, Matrix &matB)
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

Matrix inverse( Matrix &matA)
{
    register unsigned short r,c;
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

LU_Decomp::LU_Decomp( Matrix &matA )
{
// This algorithm is a distinct implementation of an algorithm found
// in Numerical Recipes in C.  The ideas expressed as code here, were 
// derived from that work, which this programmer highly recommends.
// This routine is based on Crout's method with implicit partial pivoting 
// as described by Numerical Recipes in C for the ludcmp() routine.

    register unsigned short row, col, k;
    unsigned short myRows, myColumns;
    float *row_scaling;
    float *m,*ptr;

// LU routine
//
    if ( matA.rows() != matA.cols() ) 
        Abort("LU_Decomp::LU_Decomp - Square Matrix Required!");

// Constructor
//
    size = matA.rows();
    m = (float *) new float[ size ][ size ];
    lumat = new float (*)[ size ];
    
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
    for( k=0; k<size*size; k++){
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
            tempval = lumat[row][col];
            tempval = (-tempval > tempval)? -tempval : tempval;
            if(tempval > max) max = tempval;
        }
        if (max = 0.0) cerr << "Matrix::LU_Decomp - Singular Matrix\n";
        row_scaling[ row ] = 1.0 / max;
    }

// Crout's Method
//
    for( col=0; col<size; col++){
        float sum, max, tempval;
        unsigned short max_index;

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
            if((tempval = row_scaling[row] * (-sum > sum)? -sum:sum) >= max ){
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

Matrix LU_Decomp::L()
{    
    Matrix *result;
    register short int row, col, Rows, Cols;

    result = new Matrix( size, size );

    for( row=0; row<size; row++ )
        for( col=0; col<size; col++ )
            if( row > col )
                (*result)[row][col] = lumat[row][col];
            else if ( row == col )
                (*result)[row][col] = 1.0;
            else // ( row < col )
                (*result)[row][col] = 0.0;
//
// Return the L matrix in its permuted row order.
//
    for( row=0; row<size; row++){
        float tempval;
        if (row != row_permute[ row ])
            for( col=0; col<size; col++){
                tempval = (*result)[row][col];
                (*result)[row][col] = (*result)[row_permute[ row ]][col];
                (*result)[row_permute[ row ]][col] = tempval;
            }
    }

    return (*result);
};

Matrix LU_Decomp::U()
{
    Matrix *result;
    register short int row, col, Rows, Cols;

    result = new Matrix( size, size );

    for( row=0; row<size; row++ )
        for( col=0; col<size; col++ )
            if( row <= col )
                (*result)[row][col] = lumat[row][col];
            else  // ( row > col )
                (*result)[row][col] = 0.0;
//   
// Only the L matrix is returned in its permuted row order.
//
    return (*result);
};

Matrix LU_Decomp::solve_for( Matrix &matA )
{
    Matrix *result;
    register short row, col;
    unsigned short Rows, Cols;

    Rows = matA.rows();
    Cols = matA.cols();
    result = new Matrix( matA );

    for( col=0; col<Cols; col++ ){
        int index;
        register unsigned short k;
        unsigned short row_p;
        float sum;
//
// forward substitution
//
        index = -1;
        for( row=0; row<Rows; row++ ){
            row_p = row_permute[row];
// --
            sum = (*result)[ row_p ][ col ];
            (*result)[ row_p ][ col ] = (*result)[ row ][ col ];
            if (index >= 0)
                for( k=index; k<row;  k++) 
                    sum -= lumat[row][k] * (*result)[k][col];
            else if (sum != 0.0)
                index = row;
            (*result)[row][col] = sum;
        }
//
// backsubstitution
//
        for( row=Rows-1; row>=0; row--){
            sum = (*result)[row][col];
            for( k=row+1; k<Rows; k++){
                sum -= lumat[row][k] * (*result)[k][col];
            }
            (*result)[row][col] = sum / lumat[row][row];
        }
    }

    return (*result);
};