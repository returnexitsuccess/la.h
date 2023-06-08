#include <stdlib.h>
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

    MatrixD m = addMatrixD(a, b);
    printf("a + b =\n");
    displayMatrixD(m);

    m = multiplyMatrixD(a, b);
    printf("a * b =\n");
    displayMatrixD(m);

    double da2 = fastDeterminantMatrixD(a);
    printf("|a| = %f\n\n", da2);

    double db2 = fastDeterminantMatrixD(b);
    printf("|b| = %f\n\n", db2);

    MatrixD a_echelon = echelonMatrixD(a);
    printf("Echelon form of a:\n");
    displayMatrixD(a_echelon);

    MatrixD ainv = inverseMatrixD(a);
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
    qrAlgorithmMatrixD(a, 100, d, p); // not producing upper triangular d
    printf("QR Algorithm:\n");
    printf("A is similar to:\n");
    displayMatrixD(*d);
    printf("by orthogonal matrix:\n");
    displayMatrixD(*p);


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
    }
    if (dims > 0) {
        printf("Homogeneous Solution Space:\n");
        displayMatrixD(*N);
    }
    if (dims == -1) {
        printf("No solution\n");
    }

    return 0;
}