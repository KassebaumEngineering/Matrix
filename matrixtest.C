
#ifdef LIBGpp
#include <std.h>
#endif

#include <iostream.h>
#include "Matrix.H"

float square( float );
float divide( float , float );

main()
{
    Matrix A(3,3), B(2,3), C(2,2), D(2,1), E(3,1), F(2,2);
    LU_Decomp *LU;

    A[0][0] = 1; A[0][1] = 2; A[0][2] = -3;
    A[1][0] = 2; A[1][1] = -1; A[1][2] = 4;
    A[2][0] = 4; A[2][1] = 3; A[2][2] = -2;

    B[0][0] = 1; B[0][1] = 2; B[0][2] = 3;
    B[1][0] = 4; B[1][1] = 5; B[1][2] = 6;

    C[0][0] = 2; C[0][1] = 3;
    C[1][0] = 5; C[1][1] = 7;

    D[0][0] = 1;
    D[1][0] = 3;

    E[0][0] = 6;
    E[1][0] = 2;
    E[2][0] = 14;

    F[0][0] = 2; F[0][1] = 5;
    F[1][0] = 1; F[1][1] = 3;

    cout << "A = \n";
    cout << A;

    cout << "B = \n";
    cout << B;

    cout << "transpose( A ) = \n";
    cout << transpose( A );

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

    cout << "map( &square, A ) = \n";
    cout << ( map( &square, A ) );

    cout << "map( &divide, A ) = \n";
    cout << ( map( &divide, A, A ) );

    cout << "A ^ A = \n";
    cout << (A ^ A);

    cout << "A % A = \n";
    cout << (A % A);

    cout << "LU_Decomp( A ) = \n";
    LU = new LU_Decomp( A );
    cout << LU->L() << "*\n\n" << LU->U();

    cout << "L * U = \n";
    cout << ( LU->L() * LU->U() );

    cout << "LU->solve_for( E ) = \n";
    cout << LU->solve_for( E );

    cout << "D / C  = \n";
    cout << ( D/C );

    cout << "1.0 / F = \n";
    cout << 1.0 / F;
};

float square( float a ){
    return (a*a);
};

float divide( float a, float b){
    return a/b;
};
