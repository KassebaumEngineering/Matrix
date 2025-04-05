
#include <iostream> // Use modern C++ header
#include <cmath>    // Use modern C++ header
#include <random>   // For modern C++ random numbers
#include "Matrix.H"

using namespace std; // For cout

// Setup for modern C++ random number generation
std::mt19937 gen(10348); // Mersenne Twister engine seeded with 10348 (like original srand48)
std::uniform_real_distribution<float> rnd(-1.0f, 1.0f); // Reduced range back to -1.0 to 1.0


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

// Removed the functor, sumsq will use a loop

float sumsq( const Matrix &matA ) // Accept const reference
{
   float current_sum_of_squares = 0.0f;
   // Assuming Matrix provides access to dimensions and elements
   // Adjust accessors (e.g., .rows(), .cols(), .get(r, c)) if different
   for (size_t r = 0; r < matA.rows(); ++r) { // Use size_t for indices
       for (size_t c = 0; c < matA.cols(); ++c) { // Use size_t for indices
           float val = matA[r][c]; // Or matA.get(r, c)
           current_sum_of_squares += val * val;
       }
   }
   return current_sum_of_squares;
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
int main() // Standard C++ main signature
{
    Matrix W1(1,6), W2(1,6), W3(1,6);
    float E, Elast;
    const float lambda = 1e-6f; // Regularization parameter (Reverted from 1.0f)
    int i;

    // Seeding is done at initialization of 'gen' above

    for (i=0; i<100; i++) {
        X1[0][i] = rnd(gen); // Generate random number using engine and distribution
        X2[0][i] = rnd(gen);
        X3[0][i] = rnd(gen);
        X4[0][i] = rnd(gen);
    }
   
    // Using operator& for element-wise multiplication and float literals
    Y = X1*3.0f + X2*2.0f + (X1&X1)*1.0f
        + X3*6.0f + X4*4.0f + (X3&X4)*2.0f;

// Form the 2nd degree sequences (using element-wise multiply operator &)
    Z11 =  X1 & X1 ;
    Z22 =  X2 & X2 ;
    Z33 =  X3 & X3 ;
    Z44 =  X4 & X4 ; // Added missing Z44 calculation
    Z12 =  X1 & X2 ;
    Z13 =  X1 & X3 ;
    Z14 =  X1 & X4 ;
    Z23 =  X2 & X3 ;
    Z24 =  X2 & X4 ;
    Z34 =  X3 & X4 ; // Added missing Z34 calculation

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
    Matrix XY12_T = transpose(XY12); // Transpose XY12 (1x6) to 6x1
    // Add regularization
    Matrix I6_1 = Identity(6); // Use free function Identity()
    Matrix XX12_reg = XX12 + I6_1 * lambda;
    Matrix W12_T = XY12_T / XX12_reg; // Solve (XX12 + lambda*I) * W12_T = XY12_T
    W12 = transpose(W12_T);          // Transpose result back to W12 (1x6)
    Y12 = W12 * X12;
    R12 = Y - Y12;

    XY34 = R12 * transpose( X34 );
    Matrix XY34_T = transpose(XY34); // Transpose XY34 (1x6) to 6x1
    // Add regularization
    Matrix I6_2 = Identity(6); // Use free function Identity()
    Matrix XX34_reg = XX34 + I6_2 * lambda;
    Matrix W34_T = XY34_T / XX34_reg; // Solve (XX34 + lambda*I) * W34_T = XY34_T
    W34 = transpose(W34_T);          // Transpose result back to W34 (1x6)
    Y34 = W34 * X34;
    R34 = R12 - Y34;
 
//
// node 5-6:
    X5 = Y12;
    X6 = Y34;
    Z56 = X5 & X6; // Use element-wise multiply
    Z55 = X5 & X5; // Use element-wise multiply
    Z66 = X6 & X6; // Use element-wise multiply
    for (i=0; i<100; i++) X56[0][i] = 1.0;
    for (i=0; i<100; i++) X56[1][i] = X5[0][i];
    for (i=0; i<100; i++) X56[2][i] = X6[0][i];
    for (i=0; i<100; i++) X56[3][i] = Z56[0][i];
    for (i=0; i<100; i++) X56[4][i] = Z55[0][i];
    for (i=0; i<100; i++) X56[5][i] = Z66[0][i];

    XX56 = X56 * transpose( X56 );

//

    XY56 = Y * transpose( X56 );
    Matrix XY56_T = transpose(XY56); // Transpose XY56 (1x6) to 6x1
    // Add regularization
    Matrix I6_3 = Identity(6); // Use free function Identity()
    Matrix XX56_reg = XX56 + I6_3 * lambda;
    Matrix W56_T = XY56_T / XX56_reg; // Solve (XX56 + lambda*I) * W56_T = XY56_T
    W56 = transpose(W56_T);          // Transpose result back to W56 (1x6)
    Y56 = W56 * X56;
    R56 = Y - Y56;

    cout << "Sum of Square Errors initial step:\n";
    cout << "E = " << (E = sumsq( R56 )) << "\n";

// Iteration step:
//
    // Using && in loop condition instead of comma operator, and float literals
    for( Elast=1.0e21f, i=0; i<100 && std::fabs(Elast-E) > 1.0e-6f ; i++){
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

             Matrix T12_temp(1, 100); // Temporary matrix for division result
             for( i=0; i< 100; i++){
                 float det = WB[0][i]*WB[0][i] - 4.0f*(WA[0][i]*WC[0][i]);

                 // Check det
                 if (std::isnan(det) || std::isinf(det)) {
                     cerr << "Warning: det is nan/inf at i=" << i << " in T12 calc. Setting q=0.\n";
                     q[0][i] = 0.0f; // Assign safe value
                 } else {
                     // Clamp det to avoid sqrt(negative) -> NaN
                     float safe_det = std::max(0.0f, det);
                     // Use stable quadratic formula part
                     q[0][i] = -0.5f * (WB[0][i] + sgn(WB[0][i]) * std::sqrt(safe_det));
                 }

                 // Check q
                 if (std::isnan(q[0][i]) || std::isinf(q[0][i])) {
                     cerr << "Warning: q is nan/inf at i=" << i << " in T12 calc. Setting T12_temp=0.\n";
                     T12_temp[0][i] = 0.0f; // Assign safe value
                 } else if (std::fabs(q[0][i]) < 1e-20f) { // Check for division by zero/small number
                     T12_temp[0][i] = 0.0f; // Avoid division by zero, assign a default (e.g., 0)
                 } else {
                     T12_temp[0][i] = WC[0][i] / q[0][i]; // Standard float division
                     // Check result of division
                     if (std::isnan(T12_temp[0][i]) || std::isinf(T12_temp[0][i])) {
                         cerr << "Warning: T12_temp is nan/inf after division at i=" << i << ". Setting T12_temp=0.\n";
                         T12_temp[0][i] = 0.0f; // Assign safe value
                     }
                 }
             }
             T12 = T12_temp + Y12; // Add Y12 after handling division
        }

        XY12 = T12 * transpose( X12 );
        Matrix XY12_T_iter = transpose(XY12); // Transpose XY12 (1x6) to 6x1
        // Add regularization (re-use XX12_reg if XX12 hasn't changed, otherwise recalculate)
        // Assuming XX12 is constant within the loop, we can reuse XX12_reg from line ~125
        // If XX12 *can* change, recalculate:
        // Matrix I6_iter1 = Identity(6); // Use free function Identity()
        // Matrix XX12_reg_iter = XX12 + I6_iter1 * lambda;
        // Matrix W12_T_iter = XY12_T_iter / XX12_reg_iter;
        Matrix W12_T_iter = XY12_T_iter / XX12_reg; // Using XX12_reg from initial step
        W12 = transpose(W12_T_iter);          // Transpose result back to W12 (1x6)
        Y12 = W12 * X12;
        // Clamp Y12 values to prevent explosion
        for(int k=0; k<100; ++k) {
            Y12[0][k] = std::max(-1.0e6f, std::min(1.0e6f, Y12[0][k]));
            if (std::isnan(Y12[0][k])) Y12[0][k] = 0.0f; // Also handle potential NaNs
        }
        R12 = T12 - Y12;

    //
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 & X6; // Use element-wise multiply
        Z55 = X5 & X5; // Use element-wise multiply
        Z66 = X6 & X6; // Use element-wise multiply
        for (c=0; c<100; c++) X56[0][c] = 1.0;
        for (c=0; c<100; c++) X56[1][c] = X5[0][c];
        for (c=0; c<100; c++) X56[2][c] = X6[0][c];
        for (c=0; c<100; c++) X56[3][c] = Z56[0][c];
        for (c=0; c<100; c++) X56[4][c] = Z55[0][c];
        for (c=0; c<100; c++) X56[5][c] = Z66[0][c];

        XX56 = X56 * transpose( X56 );
    
    //
    
        XY56 = Y * transpose( X56 );
        Matrix XY56_T_iter1 = transpose(XY56); // Transpose XY56 (1x6) to 6x1
        // Add regularization (XX56 is recalculated in the loop)
        Matrix I6_iter2 = Identity(6); // Use free function Identity()
        Matrix XX56_reg_iter1 = XX56 + I6_iter2 * lambda;
        Matrix W56_T_iter1 = XY56_T_iter1 / XX56_reg_iter1; // Solve (XX56 + lambda*I) * W56_T = XY56_T
        W56 = transpose(W56_T_iter1);          // Transpose result back to W56 (1x6)
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
             
             Matrix T34_temp(1, 100); // Temporary matrix for division result
             for( i=0; i< 100; i++){
                 float det = WB[0][i]*WB[0][i] - 4.0f*(WA[0][i]*WC[0][i]);

                 // Check det
                      // --- DEBUG PRINTS START ---
                      // Use ::i for outer loop iteration count (already named 'i' here)
                      // Removed impossible condition (i >= 25 && i < 5)
                      // --- DEBUG PRINTS END ---
                 if (std::isnan(det) || std::isinf(det)) {
                     cerr << "Warning: det is nan/inf at i=" << i << " in T34 calc. Setting q=0.\n";
                     q[0][i] = 0.0f; // Assign safe value
                 } else {
                     // Clamp det to avoid sqrt(negative) -> NaN
                     float safe_det = std::max(0.0f, det);
                     // Use stable quadratic formula part
                     q[0][i] = -0.5f * (WB[0][i] + sgn(WB[0][i]) * std::sqrt(safe_det));
                 }

                 // Check q
                      // --- DEBUG PRINTS START ---
                      // Removed impossible condition (i >= 25 && i < 5)
                      // --- DEBUG PRINTS END ---
                 if (std::isnan(q[0][i]) || std::isinf(q[0][i])) {
                     cerr << "Warning: q is nan/inf at i=" << i << " in T34 calc. Setting T34_temp=0.\n";
                     T34_temp[0][i] = 0.0f; // Assign safe value
                 } else if (std::fabs(q[0][i]) < 1e-20f) { // Check for division by zero/small number
                     T34_temp[0][i] = 0.0f; // Avoid division by zero, assign a default (e.g., 0)
                 } else {
                     T34_temp[0][i] = WC[0][i] / q[0][i]; // Standard float division
                     // Check result of division
                     if (std::isnan(T34_temp[0][i]) || std::isinf(T34_temp[0][i])) {
                         cerr << "Warning: T34_temp is nan/inf after division at i=" << i << ". Setting T34_temp=0.\n";
                         T34_temp[0][i] = 0.0f; // Assign safe value
                     }
                 }
             }
             T34 = T34_temp + Y34; // Add Y34 after handling division
        }

        XY34 = T34  * transpose( X34 );
        Matrix XY34_T_iter = transpose(XY34); // Transpose XY34 (1x6) to 6x1
        // Add regularization (re-use XX34_reg if XX34 hasn't changed, otherwise recalculate)
        // Assuming XX34 is constant within the loop, we can reuse XX34_reg from line ~135
        // If XX34 *can* change, recalculate:
        // Matrix I6_iter3 = Identity(6); // Use free function Identity()
        // Matrix XX34_reg_iter = XX34 + I6_iter3 * lambda;
        // Matrix W34_T_iter = XY34_T_iter / XX34_reg_iter;
        Matrix W34_T_iter = XY34_T_iter / XX34_reg; // Using XX34_reg from initial step
        W34 = transpose(W34_T_iter);          // Transpose result back to W34 (1x6)
        Y34 = W34 * X34;
        // Clamp Y34 values to prevent explosion
        for(int k=0; k<100; ++k) {
            Y34[0][k] = std::max(-1.0e6f, std::min(1.0e6f, Y34[0][k]));
             if (std::isnan(Y34[0][k])) Y34[0][k] = 0.0f; // Also handle potential NaNs
       }
        R34 = T34 - Y34;

    //
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 & X6; // Use element-wise multiply
        Z55 = X5 & X5; // Use element-wise multiply
        Z66 = X6 & X6; // Use element-wise multiply
        for (c=0; c<100; c++) X56[0][c] = 1.0;
        for (c=0; c<100; c++) X56[1][c] = X5[0][c];
        for (c=0; c<100; c++) X56[2][c] = X6[0][c];
        for (c=0; c<100; c++) X56[3][c] = Z56[0][c];
        for (c=0; c<100; c++) X56[4][c] = Z55[0][c];
        for (c=0; c<100; c++) X56[5][c] = Z66[0][c];

        XX56 = X56 * transpose( X56 );
    
    //
    
        XY56 = Y * transpose( X56 );
        Matrix XY56_T_iter2 = transpose(XY56); // Transpose XY56 (1x6) to 6x1
        // Add regularization (XX56 is recalculated in the loop)
        Matrix I6_iter4 = Identity(6); // Use free function Identity()
        Matrix XX56_reg_iter2 = XX56 + I6_iter4 * lambda;
        Matrix W56_T_iter2 = XY56_T_iter2 / XX56_reg_iter2; // Solve (XX56 + lambda*I) * W56_T = XY56_T
        W56 = transpose(W56_T_iter2);          // Transpose result back to W56 (1x6)
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
    cout <<  W12 ;

    cout << "W34 = \n";
    cout <<  W34 ;

    cout << "W56 = \n";
    cout <<  W56 ;

    cout << "E = " << E << "\n";
    cout << "Error vector = \n";
        Y12 = W12 * X12;
        Y34 = W34 * X34;
    // node 5-6:   SECOND LAYER NODE
        X5 = Y12;
        X6 = Y34;
        Z56 = X5 & X6; // Use element-wise multiply
        Z55 = X5 & X5; // Use element-wise multiply
        Z66 = X6 & X6; // Use element-wise multiply
        for (i=0; i<100; i++) X56[0][i] = 1.0;
        for (i=0; i<100; i++) X56[1][i] = X5[0][i];
        for (i=0; i<100; i++) X56[2][i] = X6[0][i];
        for (i=0; i<100; i++) X56[3][i] = Z56[0][i];
        for (i=0; i<100; i++) X56[4][i] = Z55[0][i];
        for (i=0; i<100; i++) X56[5][i] = Z66[0][i];
        Y56 = W56 * X56;
    cout << transpose(Y - Y56);

    return 0; // Return 0 from main
};
