
#include <std.h>
#include <iostream.h>
#include <String.h>
#include <Uniform.h>
#include <ACG.h>
#include <math.h>
#include "Matrix.H"

ACG gen(100, 20);
Uniform rnd(-10.0,10.0,&gen);

Matrix X1(1,100), X2(1,100), X3(1,100), Y(1,100);
Matrix Z11(1,100), Z22(1,100), Z33(1,100), Z12(1,100), Z13(1,100), Z23(1,100);

Matrix W12(1,6), W13(1,6), W23(1,6);
   
Matrix X12( 6, 100 ), X13( 6, 100), X23( 6, 100 );
Matrix XY12( 1, 6 ), XY13( 1, 6 ), XY23( 1, 6 );
Matrix XX12( 6 , 6 ), XX13( 6 , 6 ), XX23( 6 , 6 );

Matrix R12(1, 100), R13(1, 100), R23(1, 100);
Matrix Y12(1, 100), Y13(1, 100), Y23(1, 100);

float sum_of_squares;
float _sumsq( float a )
{
    float val;

    val = a * a ;
    sum_of_squares += val;

    return val;
}

float sumsq( Matrix &matA )
{
   sum_of_squares = 0.0;
   map( &_sumsq, matA );
   return sum_of_squares;
}

main()
{
    Matrix W1(1,6), W2(1,6), W3(1,6);
    float E, Elast;
    int i,j;

    for (i=0; i<100; i++) {
        X1[0][i] = rnd();
        X2[0][i] = rnd();
        X3[0][i] = rnd();
    }
   
    Y = X1*3.0 + X2*2.0 + (X2^X3)*4.0 + X3*6.0 + 5.0;

// Form the 2nd degree sequences
    Z11 =  X1 ^ X1 ;
    Z12 =  X1 ^ X2 ;
    Z13 =  X1 ^ X3 ;
    Z22 =  X2 ^ X2 ;
    Z23 =  X2 ^ X3 ;
    Z33 =  X3 ^ X3 ;

// Form the parts of the equations
//
// node 1-2:
    for (i=0; i<100; i++) X12[0][i] = 1.0;
    for (i=0; i<100; i++) X12[1][i] = X1[0][i];
    for (i=0; i<100; i++) X12[2][i] = X2[0][i];
    for (i=0; i<100; i++) X12[3][i] = Z12[0][i];
    for (i=0; i<100; i++) X12[4][i] = Z11[0][i];
    for (i=0; i<100; i++) X12[5][i] = Z22[0][i];

    XX12 = X12 * transpose( X12 );

//
// node 1-3:
    for (i=0; i<100; i++) X13[0][i] = 1.0;
    for (i=0; i<100; i++) X13[1][i] = X1[0][i];
    for (i=0; i<100; i++) X13[2][i] = X3[0][i];
    for (i=0; i<100; i++) X13[3][i] = Z13[0][i];
    for (i=0; i<100; i++) X13[4][i] = Z11[0][i];
    for (i=0; i<100; i++) X13[5][i] = Z23[0][i];

    XX13 = X13 * transpose( X13 );

//
// node 2-3:
    for (i=0; i<100; i++) X23[0][i] = 1.0;
    for (i=0; i<100; i++) X23[1][i] = X2[0][i];
    for (i=0; i<100; i++) X23[2][i] = X3[0][i];
    for (i=0; i<100; i++) X23[3][i] = Z23[0][i];
    for (i=0; i<100; i++) X23[4][i] = Z22[0][i];
    for (i=0; i<100; i++) X23[5][i] = Z33[0][i];

    XX23 = X23 * transpose( X23 );

//
// Training
//

// First step:
    XY12 = Y * transpose( X12 );
    W12 = XY12 / XX12 ;
    Y12 = W12 * X12;
    R12 = Y - Y12;

    XY13 = R12 * transpose( X13 );
    W13 = XY13 / XX13 ;
    Y13 = W13 * X13;
    R13 = R12 - Y13;
 
    XY23 = R13 * transpose( X23 );
    W23 = XY23 / XX23 ;
    Y23 = W23 * X23;
    R23 = R13 - Y23;

    cout << "Sum of Square Errors initial step:\n";

//    cout << "R12 = " << sumsq( R12 ) << "\n";
//    cout << "R13 = " << sumsq( R13 ) << "\n";
//    cout << "R23 = " << sumsq( R23 ) << "\n";
    cout << "E = " << (E = sumsq( Y - (W12*X12) - (W13*X13) - (W23*X23) )) << "\n";

// Iteration step:
//
    for( Elast=10e20, i=0; i<100, fabs((double)(Elast-E))>10e-6 ; i++){
        Matrix T12(1,100), T13(1,100), T23(1,100);

        W1 = W12;
        W2 = W13;
        W3 = W23;

        T12 = Y - Y13 - Y23;
        XY12 = T12 * transpose( X12 );
        W12  = XY12 / XX12 ;
        Y12 = W12 * X12;
        R12 = T12 - Y12;

        T13 = Y - Y12 - Y23;
        XY13 = T13  * transpose( X13 );
        W13 = XY13 / XX13 ;
        Y13 = W13 * X13;
        R13 = T13 - Y13;

        T23 = Y - Y12 - Y13;
        XY23 = T23 * transpose( X23 );
        W23 = XY23 / XX23 ;
        Y23 = W23 * X23;
        R23 = T23 - Y23;
       
        cout << "\nSum of Square Errors step " << i << "\n";
//        cout << "R12 = " << sumsq( R12 ) << "\n";
//        cout << "R13 = " << sumsq( R13 ) << "\n";
//        cout << "R23 = " << sumsq( R23 ) << "\n";
        Elast = E;
        cout << "E = " << (E = sumsq( Y - (W12*X12) - (W13*X13) - (W23*X23))) << "\n";
    }
// 
//   Weights:
//
    cout << "W12 = \n";
    cout <<  W1 ;

    cout << "W13 = \n";
    cout <<  W2 ;

    cout << "W23 = \n";
    cout <<  W3 ;

    cout << "E = " << Elast << "\n";
cout << "Error vector = \n";
    cout << transpose(Y - ( W1 * X12 ) - ( W2 * X13 ) - ( W3 * X23 ) );
};
