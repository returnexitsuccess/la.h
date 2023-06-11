#include <assert.h>
#include "la.h"

int Test_newMatrixD();
int Test_identityMatrixD();
int Test_copyMatrixD();

int main() {
    int status = 0;

    status |= Test_newMatrixD();
    status |= Test_identityMatrixD();
    status |= Test_copyMatrixD();

    return status;
}

int Test_newMatrixD() {
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    if (a.rows != dim || a.cols != dim) {
        printf("FAILED Test_newMatrixD (row and column size)\n");
        return 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (a.matrix[i][j] != 0) {
                printf("FAILED Test_newMatrixD (not all values zero)\n");
                return 1;
            }
        }
    }

    freeMatrixD(&a);

    printf("PASSED Test_newMatrixD\n");
    return 0;
}

int Test_identityMatrixD() {
    size_t dim = 10;

    MatrixD a = identityMatrixD(dim);
    if (a.rows != dim || a.cols != dim) {
        printf("FAILED Test_identityMatrixD (row and column size)\n");
        return 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (i != j && a.matrix[i][j] != 0) {
                printf("FAILED Test_identityMatrixD (non-zero off-diagonal value)\n");
                return 1;
            } else if (i == j && a.matrix[i][j] != 1) {
                printf("FAILED Test_identityMatrixD (diagonal value not one)\n");
                return 1;
            }
        }
    }

    freeMatrixD(&a);

    printf("PASSED Test_identityMatrixD\n");
    return 0;
}

int Test_copyMatrixD() {
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
        }
    }

    MatrixD b = copyMatrixD(a);
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("FAILED Test_copyMatrixD (different size)\n");
        return 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (a.matrix[i][j] != b.matrix[i][j]) {
                printf("FAILED Test_copyMatrixD (copy not equal)\n");
                return 1;
            }
            a.matrix[i][j] += 1;
            if (a.matrix[i][j] == b.matrix[i][j]) {
                printf("FAILED Test_copyMatrixD (copy not independent)\n");
                return 1;
            }
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&b);

    printf("PASSED Test_copyMatrixD\n");
    return 0;
}