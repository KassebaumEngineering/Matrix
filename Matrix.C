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
/*  $Id: Matrix.C,v 1.11 1993/11/27 00:20:22 jak Exp $
 */
//  History:
/*  $Log: Matrix.C,v $
/*  Revision 1.11  1993/11/27 00:20:22  jak
/*  Matrix Class has been ported for use with the AT&T cfront compiler version 3
/*  (with templates).   -jak
/*
 * Revision 1.10  1993/11/24  00:13:56  jak
 * Major Bug Fixed.  A new'ed pointer with offset was being incorrectly
 * deleted in LU_Decmposition. -jak
 *
 * Revision 1.9  1993/11/23  21:06:53  jak
 * Bug Fixes especially for rare cases of NUll Matrices.  -jak
 *
 * Revision 1.8  1993/11/21  08:45:45  jak
 * Changes to inverse and divide to handle the case of a single element
 * matrix as a scalar.  Also fixed the Matrix operator [] to return the
 * reference to the edge vector pointers so they can be used as lvalues
 * in an assignment.  -jak
 *
 * Revision 1.7  1993/11/20  21:53:14  jak
 * Fixed a bug in the Linked_List_Template to allow it to be correctly
 * included and used in a library situation.  -jak
 *
 * Revision 1.6  1993/11/20  03:18:43  jak
 * Added the matrix determinant function.  -jak
 *
 * Revision 1.5  1993/11/20  02:19:39  jak
 * Added Time and resource usage programs.  Also, the class is now
 * built into a library (libMatrix.a).  The Linked_List now has
 * reference counts and is correctly copied and deleted by the new
 * inc and dec ,methods for the reference count.  -jak
 *
 * Revision 1.4  1993/11/18  07:29:23  jak
 * Added alot of increased functionality, including support for
 * non-zero aligned matrices.  This supports dealing with
 * arbitrary matrix partitions.  Also, LU decompositions are
 * stored with the matrices the derived from, and are recovered
 * rather than re-computed if a matrix is re-used.   -jak
 *
 * Revision 1.3  1993/11/15  20:29:39  jak
 * Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
 **/
// =====================================

static char rcsid_MATRIX_C[] =  "$Id: Matrix.C,v 1.11 1993/11/27 00:20:22 jak Exp $";

#include <new.h>
#include <stdlib.h>
#include <math.h>

#pragma implementation
#include "Matrix.H"
#include "Linked_List_Template.H"

#define ABS(x) ((x)>=0?(x):-(x))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void Abort( const char *s ){
    cerr <<  s << "\n";
    abort();
};

Matrix::Matrix(unsigned int r, unsigned int c, const char*  aname) : 
first_row(0), number_of_rows(r), first_col(0), number_of_cols(c), m(0), decomp_list( new Linked_List< Composition > )
{
    strncpy(name, aname, MAXNAMELEN-1); 
    allocate();
};

Matrix::Matrix(unsigned int fr, unsigned int fc, int r, int c, const char* aname): 
first_row(fr), number_of_rows(r), first_col(fc), number_of_cols(c), m(0), decomp_list( new Linked_List< Composition > )
{
    strncpy(name, aname, MAXNAMELEN-1);
    allocate();
};

// Copy Constructor
Matrix::Matrix(const Matrix &matA):  first_row(matA.first_row), number_of_rows(matA.number_of_rows), first_col(matA.first_col), number_of_cols(matA.number_of_cols), m(0)
{
    register int r,c;
    register float *ptrA, *ptrB;

    decomp_list = matA.decomp_list;
	decomp_list->inc_refcount();
	
    strncpy(name, matA.name, MAXNAMELEN-1);
 
    allocate();

    if( &(m[first_row]) != (float **)0){
		for( r=first_row; r < first_row+number_of_rows; r++){
			ptrA = &(m[r][first_col]);
			ptrB = &(matA[r][first_col]);
			for( c=0; c< number_of_cols; c++){
				*ptrA++ = *ptrB++;
			}
		}
    }
};

Matrix::~Matrix( void )
{
    decomp_list->dec_refcount();
    deallocate();
};

void Matrix::deallocate()
{
	if ( &(m[first_row]) != (float **)0 ){
		delete [] &(m[first_row][first_col]);
		delete [] &(m[first_row]);
	}
};

void Matrix::allocate()
{
    register int r, size;
    register float *newmat, *ptr;

    if( (size = number_of_rows * number_of_cols) > 0 ){

        newmat = (float *) new float [ size ];
        m = (float **) new float* [ number_of_rows ];
		if ((newmat == (float  *)0 )||(m == (float **)0 )) {
			cerr << "Matrix::allocate(): couln't allocate space for a ";
			cerr << "Matrix( " << first_row << ","  ;
			cerr << number_of_rows << "," << first_col << "," ;
			cerr << number_of_cols << " )\n" ;
			abort();
		}

		for( ptr = newmat, r=0; r < size; r++){
			*ptr++ = 0.0;
		}
	
		for (ptr = newmat, r = 0; r < number_of_rows; r++){
			m[r] = ptr - first_col;
			ptr += number_of_cols;
		}
		m = m - first_row;

    } else {

        newmat = (float *)0;
        m = (float **)0;

    }
};


Matrix& Matrix::operator=(const Matrix &matB)   // performs a copy
{
    register int r,c;
    register float *ptrB, *ptrA;

	decomp_list->dec_refcount();
    decomp_list = matB.decomp_list;
	decomp_list->inc_refcount();

    if ((matB.number_of_rows != number_of_rows) || (matB.number_of_cols != number_of_cols)){
 		deallocate();
        number_of_rows = matB.number_of_rows;
        number_of_cols = matB.number_of_cols;
        first_row = matB.first_row;
        first_col = matB.first_col;
        allocate();
    } else 
        this->shift_to(matB.first_row, matB.first_col);


    if( &(m[first_row]) != (float **)0){
		for( r=first_row; r < first_row+number_of_rows; r++){
			ptrA = &(m[r][first_col]);
			ptrB = &(matB[r][first_col]);
			for( c=first_col; c< first_col+number_of_cols; c++){
				*ptrA++ = *ptrB++;
			}
		}
    }
    return (*this);
};

//  Submatrix Access Methods
Matrix Matrix:: operator() ( int fr, int fc, int rows, int cols ) const
{
    Matrix result;

    result = Ones(rows,cols);
    result.shift_to(fr,fc);

    return (result & *this);
};

Matrix& Matrix:: shift_to ( int new_fr, int new_fc )
{
    register int r;

// adjust collumns
    if( &(m[first_row]) != (float **)0 ){
		for(r = first_row; r< first_row + number_of_rows; r++)
				m[r] = m[r] + first_col - new_fc;
    }
    first_col = new_fc;

// adjust rows
    m = m + first_row - new_fr;
    first_row = new_fr;

    return *this;
};

// Matrix Operators
Matrix & Matrix::operator += ( const Matrix& matA )
{
    int fr,fc,rows,cols;
    register int r,c;
    register float *ptr_r, *ptr_a;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;
	
    fr = MIN( first_row, matA.first_row );
    fc = MIN( first_col, matA.first_col );
    rows = MAX( first_row + number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( first_col + number_of_cols, matA.first_col + matA.number_of_cols );
    rows -= fr;
    cols -= fc;

    if((rows != number_of_rows)||(cols != number_of_cols)){
        Matrix Temp( *this );
 		deallocate();
        number_of_rows = rows;
        number_of_cols = cols;
        first_row = fr;
        first_col = fc;
        allocate();
		for( r=Temp.first_row; r<Temp.first_row+Temp.number_of_rows; r++ ){  
			ptr_r = &(m[r][Temp.first_col]);
			ptr_a = &(Temp[r][Temp.first_col]);
			for( c=0; c < Temp.number_of_cols; c++)
				*ptr_r++ = *ptr_a++;
		}
    }

    for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++ ){  
        ptr_r = &(m[r][matA.first_col]);
        ptr_a = &(matA[r][matA.first_col]);
        for( c=0; c < matA.number_of_cols; c++)
            *ptr_r++ += *ptr_a++;
    }

    return (*this);
};

Matrix & Matrix::operator += ( float scalar )
{
    register int r,c;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    for( r=first_row; r<number_of_rows; r++)
        for( c=first_col; c<number_of_cols; c++)
            m[r][c] += scalar;

    return *this;
};

Matrix & Matrix::operator -= ( const Matrix& matA)
{
    int fr,fc,rows,cols;
    register int r,c;
    register float *ptr_r, *ptr_a;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    fr = MIN( first_row, matA.first_row );
    fc = MIN( first_col, matA.first_col );
    rows = MAX( first_row + number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( first_col + number_of_cols, matA.first_col + matA.number_of_cols );
    rows -= fr;
    cols -= fc;

    if((rows != number_of_rows)||(cols != number_of_cols)){
        Matrix Temp( *this );
 		deallocate();
        number_of_rows = rows;
        number_of_cols = cols;
        first_row = fr;
        first_col = fc;
        allocate();
		for( r=Temp.first_row; r<Temp.first_row+Temp.number_of_rows; r++ ){  
			ptr_r = &(m[r][Temp.first_col]);
			ptr_a = &(Temp[r][Temp.first_col]);
			for( c=0; c < Temp.number_of_cols; c++)
				*ptr_r++ = *ptr_a++;
		}
    }

    for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++ ){  
        ptr_r = &(m[r][matA.first_col]);
        ptr_a = &(matA[r][matA.first_col]);
        for( c=0; c < matA.number_of_cols; c++)
            *ptr_r++ -= *ptr_a++;
    }

    return (*this);
};

Matrix & Matrix::operator -= ( float scalar)
{
    register int r,c;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    for( r=first_row; r<number_of_rows; r++)
        for( c=first_col; c<number_of_cols; c++)
            m[r][c] -= scalar;

    return *this;
};

Matrix & Matrix::operator *= ( const Matrix& matA)
{
    register int r,c,q;
    register float *ptr_r, *ptr_a, *ptr_b;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

// MxQ times QxN
    if ( number_of_cols == matA.number_of_rows ){
        Matrix Temp( *this );
        Matrix A_bycols;
        A_bycols = transpose( matA );
 		deallocate();
        number_of_cols = matA.number_of_cols;
        first_col = matA.first_col;
        allocate();
		for( r=Temp.first_row; r<Temp.first_row+Temp.number_of_rows; r++ ){  
//			ptr_r = &(m[r][Temp.first_col]);
			ptr_r = &(m[r][first_col]);
            for (c=first_col; c<first_col+ number_of_cols; c++, ptr_r++){
                ptr_a = &(Temp[r][Temp.first_col]);
                ptr_b = &(A_bycols[c][matA.first_row]);
                *ptr_r = 0.0;
                for (q=0; q<Temp.number_of_cols; q++)
                    *ptr_r += *ptr_a++ * *ptr_b++; 
            }
		}
     }

// 1x1 times NxM
    else if ((number_of_rows == 1) && (number_of_cols == 1)){
        *this = m[first_row][first_col] * matA;

    } else if ((matA.number_of_rows == 1) && (matA.number_of_cols == 1)){
         *this  *= matA[matA.first_row][matA.first_col];

    }

// Error 
    else {
        Abort("operator*=(const Matrix& ):Incompatible Matrix Dimensions");
    }

    return *this;
};

Matrix & Matrix::operator *= ( float scalar )
{
    register int r,c;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    for( r=first_row; r<number_of_rows; r++)
        for( c=first_col; c<number_of_cols; c++)
            m[r][c] *= scalar;

    return *this;
};

Matrix & Matrix::operator /= ( Matrix& matA)
{
	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    *this = *this / matA ;
	
	return *this;
};

Matrix & Matrix::operator &= ( const Matrix& matA)
{
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    fr = MIN( first_row, matA.first_row );
    fc = MIN( first_col, matA.first_col );
    rows = MAX( first_row + number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( first_col + number_of_cols, matA.first_col + matA.number_of_cols );
    rows = rows - fr;
    cols = cols - fc;

    if((rows > 0) && (cols > 0)){
        Matrix Temp( *this );
 		deallocate();
        first_row = fr;
        first_col = fc;
        number_of_rows = rows;
        number_of_cols = cols;
        allocate();
		for( r=fr; r < fr + rows; r++){
			ptr_r = &(m[r][fc]);
			ptr_a = &(Temp[r][fc]);
			ptr_b = &(matA[r][fc]);
			for( c=0; c < cols; c++)
				*ptr_r++ = *ptr_a++ * *ptr_b++ ;
		}
    } else {
		ptr_r = &(m[first_row][first_col]);
        for( r=0; r<number_of_rows; r++)
            for( c=0; c<number_of_cols; c++)
                *ptr_r++ = 0.0;
    }
    return *this;
};

Matrix & Matrix::operator %= ( const Matrix & matA)
{
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    fr = MIN( first_row, matA.first_row );
    fc = MIN( first_col, matA.first_col );
    rows = MAX( first_row + number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( first_col + number_of_cols, matA.first_col + matA.number_of_cols );
    rows = rows - fr;
    cols = cols - fc;

    if((rows > 0) && (cols > 0)){
        Matrix Temp( *this );
 		deallocate();
        first_row = fr;
        first_col = fc;
        number_of_rows = rows;
        number_of_cols = cols;
        allocate();
		for( r=fr; r < fr + rows; r++){
			ptr_r = &(m[r][fc]);
			ptr_a = &(Temp[r][fc]);
			ptr_b = &(matA[r][fc]);
			for( c=0; c < cols; c++)
				*ptr_r++ = *ptr_a++ / *ptr_b++ ;
		}
    } else {
		ptr_r = &(m[first_row][first_col]);
        for( r=0; r<number_of_rows; r++)
            for( c=0; c<number_of_cols; c++)
                *ptr_r++ = 0.0;
    }
    return *this;
};

Matrix & Matrix::apply ( float(*f)( float ) )
{
    *this = map( f, *this );
	 
	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    return *this;
};

Matrix & Matrix::apply ( float(*f)( float, float ), const Matrix& matA)
{
    *this = map( f, *this , matA);

	decomp_list->dec_refcount();
    decomp_list = new Linked_List< Composition >;

    return *this;
};

// Matrix Norms
// default Euclidean norm p = 2
// p = 1 -> 1-norm, p=0 -> infinity norm
float norm ( const Matrix& matA, float p )
{
    register int r,c;
    register double temp, max_val;

    if ( p == 0.0 ){        // infinity norm (max row-sum norm)
        for( max_val=0.0, r=matA.first_row;
            r< matA.first_row+matA.number_of_rows; r++){
			for(temp=0.0,c=matA.first_col;
                c<matA.first_col+matA.number_of_cols; c++){
                temp += fabs( matA[r][c] );
            }
            max_val = MAX( max_val, temp );
        }
    } else if ( p == 1.0 ){ // 1 - norm      (max col-sum norm)
		for(max_val=0.0,c=matA.first_col;
			c<matA.first_col+matA.number_of_cols; c++){
			for(temp=0.0, r=matA.first_row;
				r< matA.first_row+matA.number_of_rows; r++){
                temp += fabs( matA[r][c] );
            }
            max_val = MAX( max_val, temp );
        }
    } else {              // p - norm => p=2 euclidean norm
		for(max_val=0.0,c=matA.first_col;
			c<matA.first_col+matA.number_of_cols; c++){
			for(temp=0.0, r=matA.first_row;
				r< matA.first_row+matA.number_of_rows; r++){
                temp += fabs( pow((double)matA[r][c],(double)p ));
            }
            temp = pow( temp , (1.0 / p) );
            max_val = MAX( max_val, temp );
        }
    }

    return (float)max_val;
};

float determinant ( Matrix& matA )
{
    LU_Decomposition lu;
    unsigned int size;

    if ((size = matA.rows()) != matA.cols())
        Abort("determinant( const Matrix& ):Singular! - Matrix Argument is not square!");

    lu = LU_Decomposition( matA );
    
    return lu.determinant();
}; 

// Comparison Operations
int operator == ( const Matrix& matA, const Matrix&matB )
{
    register int r,c;
    register float result;

    if ( matA.number_of_rows != matB.number_of_rows ) return 0;
    if ( matA.first_row != matB.first_row )  return 0;
    if ( matA.number_of_cols != matB.number_of_cols )  return 0;
    if ( matA.first_col != matB.first_col )  return 0;

    result = 0.0;
    for( r=matA.first_row; r<matA.first_row + matA.number_of_rows; r++)
        for( c=matA.first_col; c<matA.first_col + matA.number_of_cols; c++)
            result += fabs( matA[r][c] - matB[r][c] );

    return (result == 0.0);
};

Matrix transpose( const Matrix &matA )
{
    register float *ptr_dest;
    register int r,c;
    Matrix   result;

    result = Matrix( matA.first_col, matA.first_row, matA.number_of_cols, matA.number_of_rows );

    if ( !result.is_empty() ){
		for( r=result.first_row; r < result.first_row + result.number_of_rows; r++ ){
			ptr_dest = &(result[r][result.first_col]);
			for( c=result.first_col; c < result.first_col + result.number_of_cols; c++ ){
				*ptr_dest++ = matA[c][r];
			}
		}
    }
    return result;
};

Matrix operator+( const Matrix &matA, const Matrix &matB )
{
    Matrix   result;
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

    fr = MIN( matB.first_row, matA.first_row );
    fc = MIN( matB.first_col, matA.first_col );
    rows = MAX( matB.first_row + matB.number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( matB.first_col + matB.number_of_cols, matA.first_col + matA.number_of_cols );
    rows -= fr;
    cols -= fc;

    result = Matrix( fr, fc, rows,  cols );
    
    if ( !result.is_empty() ){
		for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++ ){  
			ptr_r = &(result[r][matA.first_col]);
			ptr_a = &(matA[r][matA.first_col]);
			for( c=0; c < matA.number_of_cols; c++)
				*ptr_r++ = *ptr_a++;
		}
		for( r=matB.first_row; r<matB.first_row+matB.number_of_rows; r++ ){  
			ptr_r = &(result[r][matB.first_col]);
			ptr_b = &(matB[r][matB.first_col]);
			for( c=0; c < matB.number_of_cols; c++)
				*ptr_r++ += *ptr_b++;
		}
    }

    return result;
};

Matrix operator-( const Matrix &matA, const Matrix &matB )
{
    Matrix   result;
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

    fr = MIN( matB.first_row, matA.first_row );
    fc = MIN( matB.first_col, matA.first_col );
    rows = MAX( matB.first_row + matB.number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MAX( matB.first_col + matB.number_of_cols, matA.first_col + matA.number_of_cols );
    rows -= fr;
    cols -= fc;

    result = Matrix( fr, fc, rows, cols );
    
    if ( !result.is_empty() ){
		for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++ ){  
			ptr_r = &(result[r][matA.first_col]);
			ptr_a = &(matA[r][matA.first_col]);
			for( c=0; c < matA.number_of_cols; c++)
				*ptr_r++ = *ptr_a++;
		}
		for( r=matB.first_row; r<matB.first_row+matB.number_of_rows; r++ ){  
			ptr_r = &(result[r][matB.first_col]);
			ptr_b = &(matB[r][matB.first_col]);
			for( c=0; c < matB.number_of_cols; c++)
				*ptr_r++ -= *ptr_b++;
		}
    }

    return result;
};

// row_to_col = first_row - first_col;
// col_to_row = first_col - first_row;

Matrix operator*( const Matrix &matA, const Matrix &matB )
{    
    Matrix   result;
    int fr, fc, rows, cols;
    register int r,c,q;
    register float *ptr_r, *ptr_a, *ptr_b;

// NUll Matrix in multiply
    if( matA.is_empty() || matB.is_empty()) {
        result = Matrix(0, 0, 0, 0);

// MxQ times QxN
    } else if ( matA.number_of_cols == matB.number_of_rows ){
        Matrix B_bycols;
        B_bycols = transpose( matB );
        fr = matA.first_row;
        rows = matA.number_of_rows;
        fc = matB.first_col;
        cols = matB.number_of_cols;
        result = Matrix(fr, fc, rows,  cols );

		for (r=fr; r< fr + rows; r++){
			ptr_r = &(result[r][fc]);
			for (c=fc; c<fc + cols; c++, ptr_r++){
				ptr_a = &(matA[r][matA.first_col]);
				ptr_b = &(B_bycols[c][matB.first_row]);
				*ptr_r = 0.0;
				for (q=0; q<matA.number_of_cols; q++,ptr_a++,ptr_b++)
					*ptr_r += (*ptr_a * *ptr_b); 
			}
		}

// 1x1 times NxM
    } else if ((matA.number_of_rows == 1) && (matA.number_of_cols == 1)){
        result = matA[matA.first_row][matA.first_col] * matB;

    } else if ((matB.number_of_rows == 1) && (matB.number_of_cols == 1)){
        result = matA * matB[matB.first_row][matB.first_col];

// Error 
    } else {
        Abort("operator*(const Matrix&, const Matrix& ):Incompatible Matrix Dimensions");
    }

    return result;
};

Matrix operator*( float scalar, const Matrix &matA )
{
    Matrix result( matA );
    register int r,c;

	for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++)
		for( c=matA.first_col; c<matA.first_col+matA.number_of_cols; c++)
			result[r][c] = matA[r][c] * scalar;

    return result;
};

Matrix operator+( const Matrix &matA, float scalar){
    Matrix result( matA );
    register int r,c;

	for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++)
		for( c=matA.first_col; c<matA.first_col+matA.number_of_cols; c++)
			result[r][c] = matA[r][c] + scalar;

    return result;
};

Matrix map( float (*f)(float), const Matrix &matA){
    Matrix result( matA );
    register int r,c;

	for( r=matA.first_row; r<matA.first_row+matA.number_of_rows; r++)
		for( c=matA.first_col; c<matA.first_col+matA.number_of_cols; c++)
			result[r][c] = (*f)(matA[r][c]);
    return result;
};

Matrix map( float (*f)(float,float), const Matrix &matA, const Matrix &matB){
    Matrix result;
    register int r,c;

    if ( matA.number_of_rows != matB.number_of_rows ) 
        Abort("map(float(*)(float,float),const Matrix&,const Matrix&):Incompatible Row Dimensions. \n");
    if ( matA.first_row != matB.first_row ) 
        Abort("map(float(*)(float,float),const Matrix&,const Matrix&):Incompatible Row Indices. \n");
    if ( matA.number_of_cols != matB.number_of_cols ) 
        Abort("map(float(*)(float,float),const Matrix&,const Matrix&):Incompatible Column Dimensions.\n");
    if ( matA.first_col != matB.first_col ) 
        Abort("map(float(*)(float,float),const Matrix&,const Matrix&):Incompatible Column Indices. \n");

    result = Matrix( matA );
	for( r=matA.first_row; r<matA.first_row + matA.number_of_rows; r++)
		for( c=matA.first_col; c<matA.first_col + matA.number_of_cols; c++)
			result[r][c] = (*f)( matA[r][c], matB[r][c] );
    return result;
};

Matrix operator & ( const Matrix &matA, const Matrix &matB){
    Matrix result;
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

    fr = MAX( matB.first_row, matA.first_row );
    fc = MAX( matB.first_col, matA.first_col );
    rows = MIN( matB.first_row + matB.number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MIN( matB.first_col + matB.number_of_cols, matA.first_col + matA.number_of_cols );
    rows = rows - fr;
    cols = cols - fc;

    if((rows > 0) && (cols > 0)){
		result = Matrix( fr, fc, rows,  cols );
		for( r=fr; r < fr + rows; r++){
			ptr_r = &(result[r][fc]);
			ptr_a = &(matA[r][fc]);
			ptr_b = &(matB[r][fc]);
			for( c=0; c < cols; c++)
				*ptr_r++ = *ptr_a++ * *ptr_b++ ;
		}
    } else {
		result = Matrix(matA.rows(),matA.cols(),"Empty");
    }
    return result;
};

Matrix operator % ( const Matrix &matA, const Matrix &matB){
    Matrix result;
    int fr, fc, rows, cols;
    register int r,c;
    register float *ptr_r, *ptr_a, *ptr_b;

    fr = MAX( matB.first_row, matA.first_row );
    fc = MAX( matB.first_col, matA.first_col );
    rows = MIN( matB.first_row + matB.number_of_rows, matA.first_row + matA.number_of_rows );
    cols = MIN( matB.first_col + matB.number_of_cols, matA.first_col + matA.number_of_cols );
    rows = rows - fr;
    cols = cols - fc;

    if((rows > 0) && (cols > 0)){
		result = Matrix( fr, fc, rows, cols );
		for( r=fr; r < fr + rows; r++){
			ptr_r = &(result[r][fc]);
			ptr_a = &(matA[r][fc]);
			ptr_b = &(matB[r][fc]);
			for( c=0; c < cols; c++)
				*ptr_r++ = *ptr_a++ / *ptr_b++ ;
		}
    } else {
		result = Matrix(matA.rows(),matA.cols(),"Empty");
    }
    return result;
};

ostream & operator << (ostream &cbuf, const Matrix &matA)
{
    register float *ptr;
    register unsigned r,c;

    if ( strlen( matA.name ) > 0)
        cbuf << matA.name << " ";
    cbuf << "at ( " << matA.first_row ;
    cbuf << " , " << matA.first_col << " )\n";
    if( matA.is_empty() )
        cbuf << "Null Matrix\n";
    else {
		for (r=matA.first_row; r < matA.first_row + matA.number_of_rows; r++){
			ptr = &(matA[r][matA.first_col]);
			for ( c=0; c < matA.number_of_cols; c++){
				cbuf << *ptr++ << "\t";
			}
			cbuf << "\n";
		}
		cbuf << "\n";
    }
    return cbuf;
};

istream & operator >> (istream &cbuf, const Matrix &matA)
{
    register float *ptr;
    register unsigned short r,c;

    for (r=matA.first_row; r < matA.first_row + matA.number_of_rows; r++){
        ptr = &(matA[r][matA.first_col]);
        for (c=0; c < matA.number_of_cols; c++){
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
  // if the second case occurs, then the alg cannot
  // recognize it. But you can transpose the input to fix it.
  //
    if ((matB.rows() == 1)&&(matB.cols() == 1)){
        return (matA * (1.0 / matB[matB.firstrow()][matB.firstcol()] ));
    } else {
        LU_Decomposition LU( matB );
        return  LU.solve_for( matA );
    }
};

Matrix inverse( Matrix &matA )
{
    Matrix b, soln;
    LU_Decomposition lu;
    unsigned int size;

    if ((size = matA.rows()) != matA.cols())
        Abort("inverse( const Matrix& ):Singular! - Matrix Argument is not square!");
    // else
    if ((matA.rows() == 1)&&(matA.cols() == 1)){
        soln = matA;
        soln[matA.firstcol()][matA.firstrow()] = 1.0 / matA[matA.firstcol()][matA.firstrow()];
    } else {
		b = Identity(size);
		b.shift_to( matA.firstcol(), matA.firstrow() );
		lu = LU_Decomposition( matA );
		soln = ( lu.solve_for( b ) );
    }
	return soln;
};

// Decompositions
void    Matrix:: setDecompose ( Composition *comA, int  key )
{
    decomp_list->add( comA, key );
};

Composition* Matrix:: getDecompose ( int key ) const
{
    return decomp_list->find( key );
};

void    Matrix:: delDecompose ( int  key )
{
    decomp_list->del( key );
};

// ==============================================
// Special Matrices
//

Matrix Ones( int rows, int cols ){
	register int i,j;
	Matrix id(0,0,rows,cols,"Ones");
	for(i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			id[i][j] = 1.0;
	return id;
};
Matrix Identity( int size ){
	register int i;
	Matrix id(0,0,size,size,"Identity");
	for(i=0;i<size;i++)
		id[i][i] = 1.0;
	return id;
};
Matrix UpperTriangle( int  size ){
	register int i,j;
	Matrix id(0,0,size,size,"UpperTriangular");
	for(i=0;i<size;i++)
		for(j=i;j<size;j++)
			id[i][j] = 1.0;
	return id;
};
Matrix LowerTriangle( int  size ){
	register int i,j;
	Matrix id(0,0,size,size,"LowerTriangular");
	for(i=0;i<size;i++)
		for(j=0;j<=i;j++)
			id[i][j] = 1.0;
	return id;
};

// ==============================================
// Decomposition Implementations
//

LU_Decomposition::LU_Decomposition( void )
{
};

// copy constructor
LU_Decomposition::LU_Decomposition( const LU_Decomposition &lu )
{
    if ( ! lu.is_empty() ){
		lumat = lu.lumat;
		row_permute = lu.row_permute;
    }
};

// Assignment Operations
LU_Decomposition&  LU_Decomposition::operator = (const LU_Decomposition &lu)
{
    if ( ! lu.is_empty() ){
		lumat = lu.lumat;
		row_permute = lu.row_permute;
    }
	return *this;
};

//
// LU Decomposition of a Matrix
// 
LU_Decomposition::LU_Decomposition( Matrix &matA )
{
// This algorithm is a distinct implementation of an algorithm found
// in Numerical Recipes in C.  The ideas expressed as code here, were 
// derived from that work, which this programmer highly recommends.
// This routine is based on Crout's method with implicit partial pivoting 
// as described by Numerical Recipes in C.

    int first_row, first_col, col_to_row, size;
    register int row, col, k;
    float *row_scaling;
    float permute_parity;  // even # of row permutations = +1.0, odd = -1.0
    LU_Decomposition *temp;

    permute_parity = 1.0;
	
    if ( matA.is_empty() ){
        return;
    }

// Get MatA's stored LU Decomposition if available
    temp = (LU_Decomposition *)( matA.getDecompose( LU_DECOMP ) );
    if ( temp != (LU_Decomposition *)0 ) {
        *this = *temp;
        return;
    }

//
// LU routine initialization
//
    if ( matA.rows() != matA.cols() ) 
        Abort("LU_Decompose( const Matrix& ): Square Matrix Required!");
    size =  matA.rows(); // == matA.cols()
    first_row = matA.firstrow();
    first_col = matA.firstcol();

    lumat = matA;  
    row_permute = IntArray( first_row,  matA.rows());

    row_scaling = new float[ size ];
    if (row_scaling  == (float *)0 ){
        Abort("LU_Decompose( const Matrix& ): Could not allocate memory for row_scaling[]!\n");
    }
    row_scaling = row_scaling - first_row;

//
// Find the Maximum element value in each row and record a scaling factor
//
    for( row=first_row; row < first_row + size; row++){
        float max, tempval;

        max = TINY;
        for( col=first_col; col < first_col + size; col++ ){
            tempval = fabs( lumat[row][col] );
            if(tempval > max) max = tempval;
        }
        if (max = 0.0) cerr << "Matrix::LU_Decomp - Singular Matrix\n";
        row_scaling[ row ] = 1.0 / max;
    }

//
// Crout's Method of Gaussian Elimination
//
    col_to_row = first_row - first_col;
    for( col=first_col; col < first_col + size; col++){
        float sum, max, tempval;
        unsigned int max_index;

        for( row=first_row; row < col+col_to_row; row++) {
            sum = lumat[row][col];
            for( k=0; k<row-first_row ; k++)
                sum -= lumat[row][k+first_col] * lumat[k+first_row][col];
            lumat[row][col] = sum;
        }
        max = TINY;
        max_index = col;
        for( row = col+col_to_row; row < first_row + size; row++){
            sum = lumat[row][col];
            for( k=0; k<col-first_col; k++)
                sum -= lumat[row][k+first_col] * lumat[k+first_row][col];
            lumat[row][col] = sum;
            if((tempval = row_scaling[row] * fabs( sum )) >= max ){
                max_index = row;
                max = tempval;
            }
        }
        if ( col != max_index ){
            for( k=0; k<size; k++ ){
                tempval = lumat[max_index][k+first_col];
                lumat[max_index][k+first_col] 
                    = lumat[ col+col_to_row ][k+first_col];
                lumat[ col+col_to_row ][k+first_col] = tempval;
            }
            permute_parity = -permute_parity;
            row_scaling[ max_index ] = row_scaling[ col+col_to_row ];
        }
        row_permute[ col+col_to_row ] = max_index ;
        if ( lumat[ col+col_to_row ][col] == 0.0 )
            lumat[ col+col_to_row ][col] = TINY;
        if ( col != size-1 ){
            tempval = 1.0 / (lumat[ col+col_to_row ][col]);
            for( row=col+col_to_row+1; row < first_row + size; row++) 
                lumat[row][col] *= tempval;
        }
    }

//
// End LU routine
// 

	matA.setDecompose( new LU_Decomposition( *this ), LU_DECOMP );

    delete &(row_scaling[first_row]); 
};

LU_Decomposition::~LU_Decomposition()
{
};


Matrix LU_Decomposition::solve_for( const Matrix &matA )
{
    Matrix result;
    register int row, col;
    int first_row, first_col;
    int Rows, Cols;

    if( matA.rows() != lumat.rows() ){
        Abort("LU_Decomposition::solve_for(const Matrix&): Incompatible Row Dimensions!\n");
    }
    if((matA.firstrow() != lumat.firstrow())
           ||(matA.firstcol() != lumat.firstcol())){
        this->shift_to( matA.firstrow(), matA.firstcol());
    }

    first_row = matA.firstrow();
    Rows = matA.rows();
    first_col = matA.firstcol();
    Cols = matA.cols();

    result = matA;
    for( col=first_col; col < first_col + Cols; col++ ){
        int flag, index;
        register unsigned short k;
        unsigned short row_p;
        float sum;
//
// forward substitution
//
        flag = 0; index = -100000;
        for( row=first_row; row < first_row + Rows; row++ ){
            row_p = row_permute[row];
// --
            sum = result[ row_p ][ col ];
            result[ row_p ][ col ] = result[ row ][ col ];
            if (flag)
                for( k=index; k < row;  k++) 
                    sum -= lumat[row][k+first_col-first_row] * result[k][col];
            else if (sum != 0.0){
                flag = 1; index = row;
            }
            result[row][col] = sum;
        }
//
// backsubstitution
//
        for( row=first_row + Rows - 1; row >= first_row; row--){
            sum = result[row][col];
            for( k=row+1; k < first_row + Rows; k++){
                sum -= lumat[row][k+first_col-first_row] * result[k][col];
            }
            result[row][col] = sum / lumat[row][row+first_col-first_row];
        }
    }

    result.setName("Solution");
    return result;
};

float  LU_Decomposition:: determinant() 
{
    register int i, col_to_row;
    float result;

    col_to_row = lumat.firstrow() - lumat.firstcol();
    result = 1.0;
    for(i=lumat.firstrow(); i< lumat.firstrow()+lumat.rows();i++)
        result *= lumat[i+col_to_row][i];

    return result;
};

Matrix LU_Decomposition:: operator[] ( int  key ) const
{
    register int i,j;
    Matrix rtn( lumat );
    Matrix temp;
    float value;

    if ( key == L_OF_LU ){
        temp = LowerTriangle( rtn.rows() );
        temp.shift_to(rtn.firstrow(), rtn.firstcol());
        rtn = rtn & temp;
        for(i=rtn.firstrow(); i<rtn.firstrow()+rtn.rows(); i++)
            rtn[i][i] = 1.0;
    } else if (key == U_OF_LU) {
        temp = UpperTriangle( rtn.rows() );
        temp.shift_to(rtn.firstrow(), rtn.firstcol());
        rtn = rtn & temp;
    } else if (key == P_OF_LU) {
        temp = Identity( rtn.rows() );
        temp.shift_to(rtn.firstrow(), rtn.firstcol());
        for(i=rtn.firstrow()+rtn.rows()-1; i>=rtn.firstrow(); i--){
            for(j=rtn.firstcol(); j< rtn.firstcol()+rtn.cols(); j++){
                value = temp[i][j]; 
                temp[i][j] = temp[ row_permute[i] ][j]; 
                temp[ row_permute[i] ][j] = value;
            }
        }
        rtn = temp;
    } else {
        rtn = Matrix();
    }

    return rtn;
};

//
//  IntArray definitiions
//

void IntArray:: allocate()
{
    register int r;
    register int *ptr;

    if( length > 0 ){
        array = (int *) new int [ length ];
		if (array == (int *)0 ) {
			cerr << "IntArray::allocate(): couln't allocate space for an ";
			cerr << "IntArray( " << length  << " )\n" ;
			abort();
		}
		for( ptr = array, r=0; r < length; r++){
			*ptr++ = 0;
		}	
		array = array - first_index;
    } else {
        array = (int *)0;
    }
};
void IntArray:: deallocate()
{
	if ( array != (int *)0 ){
		delete [] &(array[first_index]);
	}
};


// Constructors and Destructors
IntArray:: IntArray( unsigned int len ): 
length(len), first_index(0)
{
    allocate();
};
IntArray:: IntArray( int fi, unsigned int len): 
length(len), first_index(fi)
{
    allocate();
};
IntArray:: IntArray( const IntArray &ar ): 
length(ar.length), first_index(ar.first_index)
{
    register int r;
    register int *ptrA, *ptrB;

    allocate();

    ptrA = &(array[first_index]);
    ptrB = &(ar[first_index]);
    for( r=first_index; r < first_index+length; r++){
         *ptrA++ = *ptrB++;
    }
};
IntArray:: ~IntArray( void )
{
    deallocate();
};

// Assignment Operations
IntArray& IntArray:: operator = (const IntArray &ar)
{
    register int r;
    register int *ptrA, *ptrB;

    if (length != ar.length){
 		deallocate();
        length = ar.length;
        first_index = ar.first_index;
        allocate();
    } else if(ar.first_index != first_index) {
        this->shift_to(ar.first_index);
    }

    ptrA = &(array[first_index]);
    ptrB = &(ar[first_index]);
    for( r=first_index; r < first_index+length; r++){
         *ptrA++ = *ptrB++;
    }
   
    return *this;
};

//  Sub array Access Methods
IntArray  IntArray:: operator() ( int fi, unsigned int len ) const 
{
    IntArray result( fi, len );
    register int r, start, stop;
    register int *ptrA, *ptrB;

    start = MAX( fi, first_index );
    stop  = MIN( fi+len, first_index+length );

    ptrA = &(result[first_index]);
    ptrB = &(array[first_index]);
    for( r=start; r < stop; r++){
         *ptrA++ = *ptrB++;
    }

    return result;
};

IntArray& IntArray:: shift_to ( int fi )
{
    register int i;

    for(i=first_index; i<first_index+length; i++)
        array[i] += ( fi - first_index );

	array -= ( fi - first_index );
    first_index = fi;

    return *this;
};

ostream & operator << (ostream &cbuf, const IntArray &Array)
{
    register int r;

    cbuf << "at ( " << Array.first_index << " )\n";
    if( Array.is_empty() )
        cbuf << "Null Array\n";
    else {
		for (r=Array.first_index; r < Array.first_index + Array.length; r++){
			cbuf << Array[r] << "  ";
		}
		cbuf << "\n";
    }
    return cbuf;
};

istream & operator >> (istream &cbuf, const IntArray &Array)
{
    register int r;

    for (r=Array.first_index; r < Array.first_index + Array.length; r++){
        cbuf >> Array[r];
    }
    return cbuf;
};
