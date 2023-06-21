#include <stdlib.h>
#include <time.h>
#include "la.h"

int main() {
    srand(0);

    size_t a_rows = 5;
    size_t a_cols = 5;

    MatrixD a = newMatrixD(a_rows, a_cols);
    for (size_t i = 0; i < a_rows; ++i) {
        for (size_t j = 0; j < a_cols; ++j) {
            a.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    size_t b_rows = 5;
    size_t b_cols = 5;

    MatrixD b = newMatrixD(b_rows, b_cols);
    for (size_t i = 0; i < b_rows; ++i) {
        for (size_t j = 0; j < b_cols; ++j) {
            b.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    printf("a =\n");
    displayMatrixD(a);

    printf("b =\n");
    displayMatrixD(b);

    MatrixD m = newMatrixD(a.rows, a.cols);
    addMatrixD(&m, &a, &b);
    printf("a + b =\n");
    displayMatrixD(m);

    multiplyMatrixD(&m, &a, &b);
    printf("a * b =\n");
    displayMatrixD(m);

    double da2 = fastDeterminantMatrixD(a);
    printf("|a| = %f\n\n", da2);

    double db2 = fastDeterminantMatrixD(b);
    printf("|b| = %f\n\n", db2);

    MatrixD a_echelon = copyMatrixD(a);
    echelonMatrixD(&a_echelon);
    printf("Echelon form of a:\n");
    displayMatrixD(a_echelon);

    MatrixD ainv = newMatrixD(a.rows, a.cols);
    inverseMatrixD(&ainv, &a);
    printf("Inverse of a:\n");
    displayMatrixD(ainv);

    MatrixD at = transposeMatrixD(a);
    printf("Transpose of a:\n");
    displayMatrixD(at);

    MatrixD *q = malloc(sizeof(MatrixD));
    MatrixD *r = malloc(sizeof(MatrixD));
    qrDecompositionMatrixD(a, q, r);
    printf("QR Decomposition of a:\n");
    printf("Q = \n");
    displayMatrixD(*q);
    printf("R = \n");
    displayMatrixD(*r);


    MatrixD *p = malloc(sizeof(MatrixD));
    MatrixD *d = malloc(sizeof(MatrixD));
    qrAlgorithmMatrixD(a, 100, d, p);
    printf("QR Algorithm:\n");
    printf("A is similar to:\n");
    displayMatrixD(*d);
    printf("by orthogonal matrix:\n");
    displayMatrixD(*p);


    MatrixD lambda = newMatrixD(a.rows, 2);
    MatrixD v = newMatrixD(a.rows, a.cols);
    eigensMatrixD(&a, &lambda, &v);
    printf("Eigenvalues of A:\n");
    displayMatrixD(lambda);
    printf("Corresponding Eigenvectors:\n");
    displayMatrixD(v);

    MatrixD vinv = newMatrixD(v.rows, v.cols);
    inverseMatrixD(&vinv, &v);
    MatrixD D = newMatrixD(a.rows, a.cols);
    multiplyMatrixD(&D, &a, &v);
    multiplyMatrixD(&D, &vinv, &D);
    printf("Computed Jordan Form of A:\n");
    displayMatrixD(D);

    freeMatrixD(&vinv);
    freeMatrixD(&D);
    
    freeMatrixD(&lambda);
    freeMatrixD(&v);


    printf("----------------------------------------\n\n");

    MatrixD *x = malloc(sizeof(MatrixD));
    MatrixD *N = malloc(sizeof(MatrixD));
    MatrixD A = appendMatrixD(a, b);
    printf("A = \n");
    displayMatrixD(A);

    MatrixD e1 = newMatrixD(A.rows, 1);
    e1.matrix[0][0] = 1;

    printf("Solving Ax = e1:\n");

    int dims = solveLinearMatrixD(A, e1, x, N);

    printf("Solution space has dimension %d\n\n", dims);
    if (dims >= 0) {
        printf("Particular Solution:\n");
        displayMatrixD(*x);

        freeMatrixD(x);
        free(x);
    }
    if (dims > 0) {
        printf("Homogeneous Solution Space:\n");
        displayMatrixD(*N);

        freeMatrixD(N);
        free(N);
    }
    if (dims == -1) {
        printf("No solution\n");
    }


    printf("----------------------------------------\n\n");

    printf("Regular Multiplication vs. Strassen Multiplication\n");
    printf("--------------------\n");

    freeMatrixD(&A);
    MatrixD B, C1, C2, Cdiff;

    size_t maxDim = 128; // Set to 4096 takes ~10 minutes
    for (size_t dim = 64; dim <= maxDim; dim *= 2) {
        A = newMatrixD(dim, dim);
        for (size_t i = 0; i < A.rows; ++i) {
            for (size_t j = 0; j < A.cols; ++j) {
                A.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
            }
        }
        
        B = newMatrixD(dim, dim);
        for (size_t i = 0; i < B.rows; ++i) {
            for (size_t j = 0; j < B.cols; ++j) {
                B.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
            }
        }

        clock_t t1, t2, t3;
        t1 = clock();

        C1 = newMatrixD(dim, dim);
        multiplyMatrixD(&C1, &A, &B);
        
        t2 = clock();
        
        C2 = newMatrixD(dim, dim);
        strassenMultiplyMatrixD(&C2, &A, &B);
        
        t3 = clock();

        Cdiff = newMatrixD(dim, dim);
        subtractMatrixD(&Cdiff, &C1, &C2);
        double max = 0;
        for (size_t i = 0; i < Cdiff.rows; ++i) {
            for (size_t j = 0; j < Cdiff.cols; ++j) {
                max = (fabs(Cdiff.matrix[i][j]) > max) ? fabs(Cdiff.matrix[i][j]) : max;
            }
        }
        printf("%lu x %lu\n", dim, dim);
        printf("Difference: %e\n", max);
        printf("%f | %f\n", (double) (t2 - t1) / CLOCKS_PER_SEC, (double) (t3 - t2) / CLOCKS_PER_SEC);
        printf("--------------------\n");

        freeMatrixD(&A);
        freeMatrixD(&B);
        freeMatrixD(&C1);
        freeMatrixD(&C2);
        freeMatrixD(&Cdiff);
    }

    printf("----------------------------------------\n\n");

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&m);
    freeMatrixD(&a_echelon);
    freeMatrixD(&ainv);
    freeMatrixD(&at);

    freeMatrixD(q);
    freeMatrixD(r);
    freeMatrixD(p);
    freeMatrixD(d);
    free(q);
    free(r);
    free(p);
    free(d);

    freeMatrixD(&e1);

    return 0;
}