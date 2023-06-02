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
    qrAlgorithmMatrixD(a, 10, d, p); // not producing upper triangular d
    displayMatrixD(*d);
    displayMatrixD(*p);

    return 0;
}