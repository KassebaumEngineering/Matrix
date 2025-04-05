
#include <iostream> // Use modern C++ header
#include <cmath>    // Use modern C++ header
#include <random>   // For modern C++ random numbers
#include "Matrix.H"

using namespace std; // For cout

// Setup for modern C++ random number generation
std::mt19937 gen(10134); // Mersenne Twister engine seeded with 10134 (like original srand48)
std::uniform_real_distribution<float> rnd(-10.0f, 10.0f);

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

float sumsq( const Matrix &matA ) // Accept const reference
{
   sum_of_squares = 0.0;
   map( &_sumsq, matA );
   return sum_of_squares;
}

int main() // Standard C++ main signature
{
    Matrix W1(1,6), W2(1,6), W3(1,6);
    float E, Elast;
    int i;

    // Seeding is done at initialization of 'gen' above

    for (i=0; i<100; i++) {
        X1[0][i] = rnd(gen); // Generate random number using engine and distribution
        X2[0][i] = rnd(gen);
        X3[0][i] = rnd(gen);
    }
   
    // NOTE: The original code uses operator^ for element-wise multiplication.
    // This is non-standard. Assuming Matrix class overloads it this way.
    // If compilation fails here, Matrix.C/Matrix.H might need adjustment
    // or this line needs to use the element-wise multiply operator (e.g., operator&).
    // Let's assume operator^ is overloaded correctly for now.
    Y = X1*3.0f + X2*2.0f + (X1&X2)*4.0f + X3*6.0f + 5.0f; // Using operator& for element-wise mult, and float literals

// Form the 2nd degree sequences (using element-wise multiply operator &)
    Z11 =  X1 & X1 ;
    Z12 =  X1 & X2 ;
    Z13 =  X1 & X3 ;
    Z22 =  X2 & X2 ;
    Z23 =  X2 & X3 ;
    Z33 =  X3 & X3 ;

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
    Matrix XY12_T = transpose(XY12); // Transpose XY12 (1x6) to 6x1
    Matrix W12_T = XY12_T / XX12;    // Solve XX12 * W12_T = XY12_T. Result W12_T is 6x1
    W12 = transpose(W12_T);          // Transpose result back to W12 (1x6)
    Y12 = W12 * X12;
    R12 = Y - Y12;

    XY13 = R12 * transpose( X13 );
    Matrix XY13_T = transpose(XY13); // Transpose XY13 (1x6) to 6x1
    Matrix W13_T = XY13_T / XX13;    // Solve XX13 * W13_T = XY13_T. Result W13_T is 6x1
    W13 = transpose(W13_T);          // Transpose result back to W13 (1x6)
    Y13 = W13 * X13;
    R13 = R12 - Y13;
 
    XY23 = R13 * transpose( X23 );
    Matrix XY23_T = transpose(XY23); // Transpose XY23 (1x6) to 6x1
    Matrix W23_T = XY23_T / XX23;    // Solve XX23 * W23_T = XY23_T. Result W23_T is 6x1
    W23 = transpose(W23_T);          // Transpose result back to W23 (1x6)
    Y23 = W23 * X23;
    R23 = R13 - Y23;

    cout << "Sum of Square Errors initial step:\n";

//    cout << "R12 = " << sumsq( R12 ) << "\n";
//    cout << "R13 = " << sumsq( R13 ) << "\n";
//    cout << "R23 = " << sumsq( R23 ) << "\n";
    cout << "E = " << (E = sumsq( Y - (W12*X12) - (W13*X13) - (W23*X23) )) << "\n";

// Iteration step:
//
    // Using && in loop condition instead of comma operator
    for( Elast=1.0e21f, i=0; i<100 && std::fabs(Elast-E) > 1.0e-6f ; i++){
        Matrix T12(1,100), T13(1,100), T23(1,100);

        W1 = W12;
        W2 = W13;
        W3 = W23;

        T12 = Y - Y13 - Y23;
        XY12 = T12 * transpose( X12 );
        Matrix XY12_T_iter = transpose(XY12); // Transpose XY12 (1x6) to 6x1
        Matrix W12_T_iter = XY12_T_iter / XX12;    // Solve XX12 * W12_T = XY12_T. Result W12_T is 6x1
        W12 = transpose(W12_T_iter);          // Transpose result back to W12 (1x6)
        Y12 = W12 * X12;
        R12 = T12 - Y12;

        T13 = Y - Y12 - Y23;
        XY13 = T13  * transpose( X13 );
        Matrix XY13_T_iter = transpose(XY13); // Transpose XY13 (1x6) to 6x1
        Matrix W13_T_iter = XY13_T_iter / XX13;    // Solve XX13 * W13_T = XY13_T. Result W13_T is 6x1
        W13 = transpose(W13_T_iter);          // Transpose result back to W13 (1x6)
        Y13 = W13 * X13;
        R13 = T13 - Y13;

        T23 = Y - Y12 - Y13;
        XY23 = T23 * transpose( X23 );
        Matrix XY23_T_iter = transpose(XY23); // Transpose XY23 (1x6) to 6x1
        Matrix W23_T_iter = XY23_T_iter / XX23;    // Solve XX23 * W23_T = XY23_T. Result W23_T is 6x1
        W23 = transpose(W23_T_iter);          // Transpose result back to W23 (1x6)
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
    cout <<  W12 ;

    cout << "W13 = \n";
    cout <<  W13 ;

    cout << "W23 = \n";
    cout <<  W23 ;

    cout << "E = " << E << "\n";
    cout << "Error vector = \n";
    cout << transpose(Y - ( W12 * X12 ) - ( W13 * X13 ) - ( W23 * X23 ) );

    return 0; // Return 0 from main
};

