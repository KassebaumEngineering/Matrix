//
// Test program for the Matrix Class
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
/*  $Id: matrixtest.C,v 1.6 1993/11/20 06:10:05 jak Exp $
 */
//  History:
/*  $Log: matrixtest.C,v $
/*  Revision 1.6  1993/11/20 06:10:05  jak
/*  Bug fixes and optimization turned on.   -jak
/*
 * Revision 1.5  1993/11/20  03:18:46  jak
 * Added the matrix determinant function.  -jak
 *
 * Revision 1.4  1993/11/18  07:29:26  jak
 * Added alot of increased functionality, including support for
 * non-zero aligned matrices.  This supports dealing with
 * arbitrary matrix partitions.  Also, LU decompositions are
 * stored with the matrices the derived from, and are recovered
 * rather than re-computed if a matrix is re-used.   -jak
 *
 * Revision 1.3  1993/11/15  20:29:42  jak
 * Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
 **/
// =====================================

static char rcsid_MAIN_C[] =  "$Id: matrixtest.C,v 1.6 1993/11/20 06:10:05 jak Exp $";

#ifdef LIBGpp
#include <std.h>
#endif

#include <iostream.h>
#include "Matrix.H"
#include "TimeUse.H"

float square( float );

float divide( float , float );

main()
{
    Matrix A(3,3), B(2,3), C(2,2), D(2,1), E(3,1), F(2,2);
    Matrix *Big;
    LU_Decomposition *anLU;
    TimeUse stopwatch;

    A[0][0] = 1; A[0][1] = 2; A[0][2] = -3;
    A[1][0] = 2; A[1][1] = -1; A[1][2] = 4;
    A[2][0] = 4; A[2][1] = 3; A[2][2] = 2;

    B[0][0] = 1; B[0][1] = 2; B[0][2] = 3;
    B[1][0] = 4; B[1][1] = 5; B[1][2] = 6;

    C[0][0] = 2; C[0][1] = 3;
    C[1][0] = 5; C[1][1] = 7;

    D[0][0] = 1; // D[0][1] = 1;
    D[1][0] = 3; // D[1][1] = 3;

    E[0][0] = 6; //E[0][1] = 6;  E[0][2] = 6;
    E[1][0] = 2; //E[1][1] = 2;  E[1][2] = 2;
    E[2][0] = 14;//E[2][1] = 14; E[2][2] = 14;

    F[0][0] = 2; F[0][1] = 5;
    F[1][0] = 1; F[1][1] = 3;

    cout << "A = \n";
    cout << A;

    cout << "B = \n";
    cout << B;

    cout << "transpose( A ) = \n";
    cout << (transpose( A ));

    cout << "B * A = \n";
    cout << (B * A);

    cout << "A + A = \n";
    cout << (A + A);

    cout << "A - A = \n";
    cout << (A - A);

    cout << "5.0 * A \n";
    cout << (5.0 * A);

    cout << "A * 5.0\n";
    cout << (A * 5.0);

    cout << "map( square, A ) = \n";
    cout << ( map( square, A ) );

    cout << "map( divide, A, A) = \n";
    cout << ( map( divide, A, A ) );

    cout << "A & A = \n";
    cout << (A & A);

    cout << "A % A = \n";
    cout << (A % A);

    cout << "Identity(3) = \n";
    cout <<  Identity(3);

    cout << "UpperTriangle(3) = \n";
    cout <<  UpperTriangle(3);

    cout << "LowerTriangle(3) = \n";
    cout <<  LowerTriangle(3);

    cout << "LU_Decomposition( A ) = \n";
    stopwatch.start( 1 );
    anLU =  new LU_Decomposition( A ) ;

    cout << "LU decomposition of A took " << stopwatch.stop( 1 ) << "\n";
    cout << (*anLU)[ P_OF_LU ]  << "*\n\n"<< ((*anLU)[ L_OF_LU ]) << "*\n\n" << ((*anLU)[ U_OF_LU ]);
    cout << "P * L * U = \n";
    cout << ( (*anLU)[ P_OF_LU ] * (*anLU)[ L_OF_LU ] * (*anLU)[ U_OF_LU ] );

    cout << "anLU->solve_for( E ) = \n";
    cout << (anLU->solve_for( E ));

    cout <<  "E =" << "\n";
    cout <<  E  <<"\n";

    cout << "A * (anLU->solve_for( E )) =" <<"\n";
    cout << (A * (anLU->solve_for( E ))) <<"\n";

    cout << "D / C  = \n";
    cout << ( D/C );

    cout << "inverse(C) * D = \n";
    cout << inverse(C) * D;

    cout << "1.0 / F = \n";
    cout << (1.0 / F);

    cout << "inverse(F) = \n";
    cout << inverse(F);

    cout << "inverse(A) = \n";
    cout << inverse(A);

    cout << "1.0 / A = \n";
    cout << (1.0 / A);

    cout << "inverse(A) * A = \n";
    cout << inverse(A) * A;

    cout << "A.shift_to(1,2) = \n";
    A.shift_to(1,2);
    cout << A;

    cout << "inverse(A.shift_to(1,2)) = \n";
    cout << inverse(A);

    cout << "inverse(A.shift_to(1,2)) * A.shift_to(1,2) = \n";
    cout << inverse(A) * A;
    {
        Matrix Ainv( inverse(A) );
        cout << "Ainv *= A =\n";
        cout << (Ainv *= A);
    }

    cout << "norm( A , 0.0 ) = " << norm( A, 0.0 ) << "\n" ;
    cout << "norm( A , 1.0 ) = " << norm( A, 1.0 ) << "\n" ;
    cout << "norm( A , 2.0 ) = " << norm( A, 2.0 ) << "\n" ;
    cout << "norm( A , 3.0 ) = " << norm( A, 3.0 ) << "\n" ;
    cout << "\n";
    cout << "norm( inverse(A), 0.0 ) = " << norm(inverse(A), 0.0 ) << "\n" ;
    cout << "norm( inverse(A), 1.0 ) = " << norm(inverse(A), 1.0 ) << "\n" ;
    cout << "norm( inverse(A), 2.0 ) = " << norm(inverse(A), 2.0 ) << "\n" ;
    cout << "norm( inverse(A), 3.0 ) = " << norm(inverse(A), 3.0 ) << "\n" ;

    cout << "determinant( A ) =  ";
    cout << determinant( A ) << "\n\n" ;

    cout << "Ashift_to(1,2) + B = \n";
    cout <<  A + B ;

    cout << "Ashift_to(1,2) - B = \n";
    cout <<  A - B ;

    cout << "Ashift_to(1,2) & B = \n";
    cout <<  (A & B);

    cout << "Ashift_to(1,2) % B = \n";
    cout <<  (A % B);

    cout << "A(1,2,2,2) = \n";
    cout <<  A(1,2,2,2);


    Big = new Matrix( UpperTriangle( 100 ) );
    cout << "Big = \n";

//    cout << "Big = \n";
//    cout << *Big;

//    cout << "transpose( *Big ) = \n";
//    cout << transpose( *Big );

//    cout << "LU_Decomp( *Big ) = \n";
//    LU = new LU_Decomp( *Big );
//    cout << (LU->L()) << "*\n\n" << (LU->U());
    stopwatch.start( 2 );
//    cout << "1.0 / *Big = \n";
//    cout << ( 1.0 / (*Big) );
    inverse (*Big);
    cout << "Big Inversion took " << stopwatch.stop( 2 ) << "\n";
    delete Big;

    cout << "Finished after " << stopwatch.stop() << "\n";
};

float square( float a ){
    return (a*a);
};

float divide( float a, float b){
    return a/b;
};
