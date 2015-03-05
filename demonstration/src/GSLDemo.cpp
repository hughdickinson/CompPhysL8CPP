// GSLDemo.cpp - Demonstrates the use of the GNU Scientific Library

// The STL <iostream> header provides terminal I/O functionality
#include <iostream>
// The STL <iomanip> header provides advanced formatting of textual I/O
#include <iomanip>
/* The STL <algorithm> header provides the std::copy alogrithm that is
 * used to initialize matrix values.
 */
#include <algorithm>

// Include header files provided by the GSL

// The "gsl_vector.h" header provides a gsl_vector type
#include "gsl/gsl_vector.h"
// The "gsl_matrix.h" header provides a gsl_matrix type
#include "gsl/gsl_matrix.h"
/* The "gsl_blas.h" header provides an interface to the
 * Basic Linear Algebra SubPrograms (BLAS) framework.
 */
 #include "gsl/gsl_blas.h"

/* A helper function to output a summary of a GSL matrix
 * to the terminal.
 *
 * The GSL also provides output functions that can interface
 * with the C I/O API that is provided by the <cstdio> header.
 * However, that output is not as pretty!
 */
void printMatrix(gsl_matrix * mat){
    std::cout << "[ ";
    for(int i = 0; i < mat->size1; ++i){
        if(i > 0){
            std::cout << "  ";
        }
        for(int j = 0; j < mat->size2; ++j){
            std::cout << "(" << i << ", " << j <<  ") : " << std::setw(4) << gsl_matrix_get(mat, i, j);
            if(j < mat->size2 - 1){
                std::cout << ",\t";
            }
            else if(i*j < (mat->size1-1)*(mat->size2-1)){
                std::cout << ",";
            }
        }
        if (i == mat->size1 -1 ){
            std::cout << " ]";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

/* A helper function to output a summary of a GSL vector
 * to the terminal.
 *
 * The GSL also provides output functions that can interface
 * with the C I/O API that is provided by the <cstdio> header.
 */
void printVector(gsl_vector * vec){
    std::cout << "[ ";
    for(int i = 0; i < vec->size; ++i){
        if(i > 0){
            std::cout << "  ";
        }
        std::cout << i << " :\t" << gsl_vector_get(vec, i);
        if(i < vec->size - 1){
            std::cout << ",";
        }
        else{
            std::cout << " ]";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

/* The program entry point. 
 *
 * This program demonstrates the use of routines provided by
 * the GNU Scientific Library. It demonstrates how to allocate
 * GSL-provided matrix and vector "objects" as well as how to 
 * invoke Linear Algebra functions that are provided by the GSL
 * in a separate Basic Linear Algebra SubPrograms (BLAS) library.
 *
 * 1) Allocates GSL vector and matrix "objects"
 * 2) Invokes the GSL BLAS functions to compute:
 *      a) The scalar product of two vectors.
 *      b) The product of a vector and a matrix.
 *      c) The product of two matrices.
 */
int main(int argc, const char * argv[]) {

    // BASIC INITIALIZATION
    
    /* Declare and initialize basic arrays to contain
     * local copies of the test vector and matrix elements.
     */
    double vectorComponents1[3] = {1.0, 2.0, 3.0};
    double vectorComponents2[3] = {4.0, 5.0, 6.0};
    
    unsigned int matrixDimensionality1[2] = {2, 3}; // 2 columns, 3 rows
    double matrixElements1[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    
    unsigned int matrixDimensionality2[2] = {3, 2}; // 3 columns, 2 rows
    double matrixElements2[6] = {7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

    // START OF LINEAR ALGEBRA COMPUTATIONS
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL vector with the identifier "testVector1" with three components.
     */
    gsl_vector * testVector1 = gsl_vector_calloc(3);
   
    /* Loop over all the components of testVector1 and use the
     * gsl_vector_set(...) function to assign the values from
     * vectorComponents1 to the corresponding components of
     * testVector1
     */
    for(int i = 0; i < testVector1->size; ++i){
        gsl_vector_set(testVector1, i, vectorComponents1[i]);
    }
    
    /* Output a summary of the contents of testVector1 to the
     * terminal.
     */
    std::cout << "testVector1:\n";
    printVector(testVector1);
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL vector with the identifier "testVector2" with three components.
     */
    gsl_vector * testVector2 = gsl_vector_calloc(3);
    
    /* Loop over all the components of testVector2 and use the
     * gsl_vector_set(...) function to assign the values from
     * vectorComponents2 to the corresponding components of
     * testVector2
     */
    for(int i = 0; i < testVector2->size; ++i){
        gsl_vector_set(testVector2, i, vectorComponents2[i]);
    }
    
    /* Output a summary of the contents of testVector2 to the
     * terminal.
     */
    std::cout << "testVector2:\n";
    printVector(testVector2);
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL matrix with the identifier "testMatrix1" with two dimensions.
     * The number of components for each dimension are set using the elements
     * of the matrixDimensionality1 array.
     */
    gsl_matrix * testMatrix1 = gsl_matrix_calloc(matrixDimensionality1[1], matrixDimensionality1[0]);

    /* Initialize the elements of testMatrix1.
     * 1) Obtain a pointer to the first element of the internal data 
     * storage array of "testMatrix1"
     * 2) Use the std::copy algorithm to initialize the values of the
     * internal data storage array of "testMatrix1" to be equal to the
     * elements of "matrixElements1"
     */
    double * matrixElementPtr1 = gsl_matrix_ptr(testMatrix1, 0, 0);
    std::copy(matrixElements1, matrixElements1 + (testMatrix1->size1*testMatrix1->size2), matrixElementPtr1);
    
    /* Output a summary of the contents of testMatrix1 to the
     * terminal.
     */
    std::cout << "testMatrix1:\n";
    printMatrix(testMatrix1);
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL matrix with the identifier "testMatrix2" with two dimensions.
     * The number of components for each dimension are set using the elements
     * of the matrixDimensionality2 array.
     */
    gsl_matrix * testMatrix2 = gsl_matrix_calloc(matrixDimensionality2[1], matrixDimensionality2[0]);
    
    /* Initialize the elements of testMatrix1.
     * 1) Obtain a pointer to the first element of the internal data
     * storage array of "testMatrix1"
     * 2) Use the std::copy algorithm to initialize the values of the
     * internal data storage array of "testMatrix1" to be equal to the
     * elements of "matrixElements1"
     */
    double * matrixElementPtr2 = gsl_matrix_ptr(testMatrix2, 0, 0);
    std::copy(matrixElements2, matrixElements2 + 6, matrixElementPtr2);
    
    /* Output a summary of the contents of testMatrix1 to the
     * terminal.
     */
    std::cout << "testMatrix2:\n";
    printMatrix(testMatrix2);
    
    // SCALAR PRODUCT:
    
    /* Declare and zero-initialize a double precision variable to
     * receive the computed scalar product.
     */
    double scalarProduct(0.0);
    /* Declare and zero-initialize an integer variable that will
     * will be used to capture the exectution status of each GSL 
     * BLAS function that is called.
     */
    int state(0);
    
    /* Call the appropriate GSL BLAS function to perform a Scalar
     * product. Inspecting the gsl/gsl_blas.h header file
     * (or online documentation) reveals that the required function
     * signature is:
     *
     * int gsl_blas_ddot(const gsl_vector * X,
     *                   const gsl_vector * Y,
     *                   double * result
     *                   );
     *
     * This function performs the following operation:
     *
     * result = X . Y
     *
     * The function returns an integer that is defined to be zero on
     * success.
     */
    state = gsl_blas_ddot(testVector1, testVector2, &scalarProduct);
    
    /* Output a summary of the scalar product result to the
     * terminal.
     */
    std::cout << "Scalar product:\n";
    std::cout << "State = " << state << "\n";
    std::cout << scalarProduct << "\n" << std::endl;
    
    // MATRIX-VECTOR PRODUCT:
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL vector to contain the result of a Matrix-Vector multiplication.
     */
    gsl_vector * resultVector = gsl_vector_calloc(testMatrix2->size1);
    
    /* Call the appropriate GSL BLAS function to perform a Matrix
     * multiplication. Inspecting the gsl/gsl_blas.h header file 
     * (or online documentation) reveals that the required function
     * signature is:
     *
     * int gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA,
     *                    double alpha,
     *                    const gsl_matrix * A,
     *                    const gsl_vector * X,
     *                    double beta,
     *                    gsl_vector * Y);
     *
     * In fact, this function is more general than we need.
     * It performs the following operation:
     *
     * Y = alpha*A*X + beta*Y
     *
     * The mysterious CBLAS_TRANSPOSE_t type is actually an
     * ENUMERATION that can be one of three named values:
     *
     * enum CBLAS_TRANSPOSE {
     *   CblasNoTrans=111,
     *   CblasTrans=112,
     *   CblasConjTrans=113
     * };
     *
     * and specifies whether the input matrix should be transposed
     * prior to performing the multiplication.
     *
     * The alpha parameters can be used to apply a scale factors
     * to the result of the multiplication and the beta parameter
     * can be used to add a scaled copy of the original contents of
     * the Y (result) vector.
     *
     * The function returns an integer that is defined to be zero on
     * success.
     */
    state = gsl_blas_dgemv(CblasNoTrans, // do not transpose matrix
                           1.0, // do not scale multiplication result
                           testMatrix2, // input matrix
                           testVector1, // input vector
                           0.0, // do not add a copy of the "previous result"
                           resultVector // result of multiplication
                           );
    
    /* Output a summary of the Matrix-vector multiplication result 
     * to the terminal.
     */
    std::cout << "Matrix-vector multiplication:\n";
    std::cout << "State = " << state << "\n";
    printVector(resultVector);
    
    // MATRIX-MATRIX MULTIPLICATION
    
    /* Call the appropriate GSL function to allocate and zero-initialize
     * a GSL vector to contain the result of a Matrix-Vector multiplication.
     * Use the appropriate size* member datum of the input GSL matrices to
     * specify the required matrix dimension sizes.
     */
    gsl_matrix * resultMatrix = gsl_matrix_calloc(testMatrix1->size1, testMatrix2->size2);
    
    /* Call the appropriate GSL BLAS function to perform a Matrix
     * multiplication. Inspecting the gsl/gsl_blas.h header file
     * (or online documentation) reveals that the required function
     * signature is:
     *
     * int gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA,
     *                    CBLAS_TRANSPOSE_t TransB,
     *                    double alpha,
     *                    const gsl_matrix * A,
     *                    const gsl_matrix * B,
     *                    double beta,
     *                    gsl_matrix * C);
     *
     * In fact, this function is more general than we need.
     * It performs the following operation:
     *
     * Y = alpha*A*B + beta*C
     *
     * The mysterious CBLAS_TRANSPOSE_t type is actually an
     * ENUMERATION that can be one of three named values:
     *
     * enum CBLAS_TRANSPOSE {
     *   CblasNoTrans=111,
     *   CblasTrans=112,
     *   CblasConjTrans=113
     * };
     *
     * and specifies whether the input matrix should be transposed
     * prior to performing the multiplication.
     *
     * The alpha parameters can be used to apply a scale factors
     * to the result of the multiplication and the beta parameter
     * can be used to add a scaled copy of the original contents of
     * the C (result) matrix.
     *
     * The function returns an integer that is defined to be zero on
     * success.
     */
    state = gsl_blas_dgemm(CblasNoTrans,
                           CblasNoTrans,
                           1.0,
                           testMatrix1,
                           testMatrix2,
                           0.0,
                           resultMatrix);

    /* Output a summary of the Matrix-Matrix multiplication result
     * to the terminal.
     */
    std::cout << "Matrix-Matrix multiplication:\n";
    std::cout << "State = " << state << "\n";
    printMatrix(resultMatrix);
    
    /* As usual, it is important to deallocate any dynamically allocated 
     * memory. Since the GSL is a C library, it is NOT possible to use
     * the "delete" or "delete[]" operators to do so.
     *
     * Instead, it is neccessary to use the library-provided deallocation
     * functions. These are detailed in the online documentation.
     */
    gsl_vector_free(testVector1);
    gsl_vector_free(testVector2);
    gsl_matrix_free(testMatrix1);
    gsl_matrix_free(testMatrix2);
    
    /* Since "state" will be zero if all the Linar Algebra operations
     * succeeded, it makes sense to return it as the overall program
     * result.
     */
    return state;
}
