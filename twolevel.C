
#include <std.h>
#include <iostream.h>
#include <String.h>
#include <Uniform.h>
#include <ACG.h>
#include <math.h>
#include "Matrix.H"

ACG gen(100, 20);
Uniform rnd(-1.0,1.0,&gen);

Matrix X1(1,100), X2(1,100), X3(1,100), X4(1,100), Y(1,100);
Matrix Z11(1,100), Z22(1,100), Z33(1,100), Z44(1,100);
Matrix Z12(1,100), Z13(1,100), Z14(1,100);
Matrix Z23(1,100), Z24(1,100);
Matrix Z34(1,100);

Matrix X5(1,100), X6(1,100);
Matrix Z55(1,100), Z66(1,100);
Matrix Z56(1,100);

Matrix W12(1,6), W34(1,6), W56(1,6);
   
Matrix X12( 6, 100 ), X34( 6, 100), X56( 6, 100 );
Matrix XY12( 1, 6 ), XY34( 1, 6 ), XY56( 1, 6 );
Matrix XX12( 6 , 6 ), XX34( 6 , 6 ), XX56( 6 , 6 );

Matrix R12(1, 100), R34(1, 100), R56(1, 100);
Matrix Y12(1, 100), Y34(1, 100), Y56(1, 100);

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

float sgn( float a )
{
    if (a > 0.0) {
        return 1.0;
    } else if (a == 0.0 ){
        return 0.0;
    } else {
        return -1.0;
    }
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
        X4[0][i] = rnd();
    }
   
    Y = X1*3.0 + X2*2.0 + X1^X1*1.0
        + X3*6.0 + X4*4.0 + X3^X4*2.0;

// Form the 2nd degree sequences
    Z11 =  X1 ^ X1 ;
    Z22 =  X2 ^ X2 ;
    Z33 =  X3 ^ X3 ;
    Z12 =  X1 ^ X2 ;
    Z13 =  X1 ^ X3 ;
    Z14 =  X1 ^ X4 ;
    Z23 =  X2 ^ X3 ;
    Z24 =  X2 ^ X4 ;

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
// node 3-4:
    for (i=0; i<100; i++) X34[0][i] = 1.0;
    for (i=0; i<100; i++) X34[1][i] = X3[0][i];
    for (i=0; i<100; i++) X34[2][i] = X4[0][i];
    for (i=0; i<100; i++) X34[3][i] = Z34[0][i];
    for (i=0; i<100; i++) X34[4][i] = Z33[0][i];
    for (i=0; i<100; i++) X34[5][i] = Z44[0][i];

    XX34 = X34 * transpose( X34 );

//
// Training
//

// First step:
    XY12 = Y * transpose( X12 );
    W12 = XY12 / XX12 ;
    Y12 = W12 * X12;
    R12 = Y - Y12;

    XY34 = R12 * transpose( X34 );
    W34 = XY34 / XX34 ;
    Y34 = W34 * X34;
    R34 = R12 - Y34;
 
//
// node 5-6:
    X5 = Y12;
    X6 = Y34;
    Z56 = X5 ^ X6;
    Z55 = X5 ^ X5;
    Z66 = X6 ^ X6;
    for (i=0; i<100; i++) X56[0][i] = 1.0;
    for (i=0; i<100; i++) X56[1][i] = X5[0][i];
    for (i=0; i<100; i++) X56[2][i] = X6[0][i];
    for (i=0; i<100; i++) X56[3][i] = Z56[0][i];
    for (i=0; i<100; i++) X56[4][i] = Z55[0][i];
    for (i=0; i<100; i++) X56[5][i] = Z66[0][i];

    XX56 = X56 * transpose( X56 );

//

    XY56 = Y * transpose( X56 );
    W56 = XY56 / XX56 ;
    Y56 = W56 * X56;
    R56 = Y - Y56;

    cout << "Sum of Square Errors initial step:\n";
    cout << "E = " << (E = sumsq( R56 )) << "\n";

// Iteration step:
//
    for( Elast=10e20, i=0; i<100, fabs((double)(Elast-E))>10e-6 ; i++){
        Matrix T12(1,100), T34(1,100), T56(1,100);
        int c;

        W1 = W12;
        W2 = W34;
        W3 = W56;

    // Calculate new desired T12;
    //    Y56 = W0 + W1X5 + W2X6 + W3X5X6 + W4X5X5 + W5X6X6
        {    
             Matrix WA(1,100), WB(1,100), WC(1,100), q(1,100);
             int i;

             for (i=0; i<100; i++){
                 WC[0][i] = W56[0][0] + W56[0][2]*X6[0][i] + W56[0][5]*X6[0][i]*X6[0][i] - R56[0][i];
                 WB[0][i] = W56[0][1] + W56[0][3]*X6[0][i];
                 WA[0][i] = W56[0][4];
             }

             for( i=0; i< 100; i++){
                 float det;
                 det = WB[0][i]*WB[0][i] - 4.0*(WA[0][i]*WC[0][i]);
                 if (det < 0)
                     q[0][i] = -0.5 * WB[0][i];
                 else 
                     q[0][i] = -0.5*( WB[0][i] + sgn( WB[0][i] )*sqrt(det) );
             }
             
             T12 = WC % q + Y12;
        }

        XY12 = T12 * transpose( X12 );
        W12  = XY12 / XX12 ;
        Y12 = W12 * X12;
        R12 = T12 - Y12;

    //
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 ^ X6;
        Z55 = X5 ^ X5;
        Z66 = X6 ^ X6;
        for (c=0; c<100; c++) X56[0][c] = 1.0;
        for (c=0; c<100; c++) X56[1][c] = X5[0][c];
        for (c=0; c<100; c++) X56[2][c] = X6[0][c];
        for (c=0; c<100; c++) X56[3][c] = Z56[0][c];
        for (c=0; c<100; c++) X56[4][c] = Z55[0][c];
        for (c=0; c<100; c++) X56[5][c] = Z66[0][c];

        XX56 = X56 * transpose( X56 );
    
    //
    
        XY56 = Y * transpose( X56 );
        W56 = XY56 / XX56 ;
        Y56 = W56 * X56;
        R56 = Y - Y56;

    // Calculate new desired T34;
    //    Y56 = W0 + W1X5 + W2X6 + W3X5X6 + W4X5X5 + W5X6X6
        {    
             Matrix WA(1,100), WB(1,100), WC(1,100), q(1,100);
             int i;

             for (i=0; i<100; i++){
                 WC[0][i] = W56[0][0] + W56[0][1]*X5[0][i] + W56[0][4]*X5[0][i]*X5[0][i] - R56[0][i];
                 WB[0][i] = W56[0][2] + W56[0][3]*X5[0][i];
                 WA[0][i] = W56[0][5];
             }
             
             //q = -0.5 * (WB + map( &sgn, WB ) ^ map( &sqrt, (WB^WB - 4.0*(WA^WC))));
             for( i=0; i< 100; i++){
                 float det;
                 det = WB[0][i]*WB[0][i] - 4.0*(WA[0][i]*WC[0][i]);
                 if (det < 0)
                     q[0][i] = -0.5 * WB[0][i];
                 else 
                     q[0][i] = -0.5*( WB[0][i] + sgn( WB[0][i] )*sqrt(det) );
             }
             
             T34 = WC % q + Y34;
        }

        XY34 = T34  * transpose( X34 );
        W34 = XY34 / XX34 ;
        Y34 = W34 * X34;
        R34 = T34 - Y34;

    //
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 ^ X6;
        Z55 = X5 ^ X5;
        Z66 = X6 ^ X6;
        for (c=0; c<100; c++) X56[0][c] = 1.0;
        for (c=0; c<100; c++) X56[1][c] = X5[0][c];
        for (c=0; c<100; c++) X56[2][c] = X6[0][c];
        for (c=0; c<100; c++) X56[3][c] = Z56[0][c];
        for (c=0; c<100; c++) X56[4][c] = Z55[0][c];
        for (c=0; c<100; c++) X56[5][c] = Z66[0][c];

        XX56 = X56 * transpose( X56 );
    
    //
    
        XY56 = Y * transpose( X56 );
        W56 = XY56 / XX56 ;
        Y56 = W56 * X56;
        R56 = Y - Y56;
       
        cout << "\nSum of Square Errors step " << i << "\n";
        Elast = E;
        cout << "E = " << (E = sumsq( Y - Y56 )) << "\n";

    }
// 
//   Weights:
//
    cout << "W12 = \n";
    cout <<  W1 ;

    cout << "W34 = \n";
    cout <<  W2 ;

    cout << "W56 = \n";
    cout <<  W3 ;

    cout << "E = " << Elast << "\n";
    cout << "Error vector = \n";
        Y12 = W1 * X12;
        Y34 = W2 * X34;
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 ^ X6;
        Z55 = X5 ^ X5;
        Z66 = X6 ^ X6;
        for (i=0; i<100; i++) X56[0][i] = 1.0;
        for (i=0; i<100; i++) X56[1][i] = X5[0][i];
        for (i=0; i<100; i++) X56[2][i] = X6[0][i];
        for (i=0; i<100; i++) X56[3][i] = Z56[0][i];
        for (i=0; i<100; i++) X56[4][i] = Z55[0][i];
        for (i=0; i<100; i++) X56[5][i] = Z66[0][i];
        Y56 = W3 * X56;
    cout << transpose(Y - Y56);
};
