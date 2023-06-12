#include <assert.h>
#include "la.h"

static int Test_newMatrixD();
static int Test_identityMatrixD();
static int Test_copyMatrixD();
static int Test_equalsMatrixD();
static int Test_addMatrixD();

int main() {
    int status = 0;

    status |= Test_newMatrixD();
    status |= Test_identityMatrixD();
    status |= Test_copyMatrixD();
    status |= Test_equalsMatrixD();
    status |= Test_addMatrixD();

    return status;
}

static int Test_newMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    if (a.rows != dim || a.cols != dim) {
        printf("FAILED Test_newMatrixD (row and column size)\n");
        status = 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (a.matrix[i][j] != 0) {
                printf("FAILED Test_newMatrixD (not all values zero)\n");
                status = 1;
            }
        }
    }

    freeMatrixD(&a);

    if (status == 0) printf("PASSED Test_newMatrixD\n");

    return status;
}

static int Test_identityMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = identityMatrixD(dim);
    if (a.rows != dim || a.cols != dim) {
        printf("FAILED Test_identityMatrixD (row and column size)\n");
        status = 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (i != j && a.matrix[i][j] != 0) {
                printf("FAILED Test_identityMatrixD (non-zero off-diagonal value)\n");
                status = 1;
            } else if (i == j && a.matrix[i][j] != 1) {
                printf("FAILED Test_identityMatrixD (diagonal value not one)\n");
                status = 1;
            }
        }
    }

    freeMatrixD(&a);

    if (status == 0) printf("PASSED Test_identityMatrixD\n");

    return status;
}

static int Test_copyMatrixD() {
    int status = 0;
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
        status = 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (a.matrix[i][j] != b.matrix[i][j]) {
                printf("FAILED Test_copyMatrixD (copy not equal)\n");
                status = 1;
            }
            a.matrix[i][j] += 1;
            if (a.matrix[i][j] == b.matrix[i][j]) {
                printf("FAILED Test_copyMatrixD (copy not independent)\n");
                status = 1;
            }
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&b);

    if (status == 0) printf("PASSED Test_copyMatrixD\n");

    return status;
}

static int Test_equalsMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    MatrixD c = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
            b.matrix[i][j] = a.matrix[i][j] + ((double) ((i + j) % 2) - 0.5);
            c.matrix[i][j] = a.matrix[i][j] + ((double) ((i + j) % 2) - 0.5) * 2;
        }
    }

    if (!equalsMatrixD(a, a, 0)) {
        printf("FAILED Test_equalsMatrixD (a, a, 0)\n");
        status = 1;
    }

    if (!equalsMatrixD(a, a, 1)) {
        printf("FAILED Test_equalsMatrixD (a, a, 1)\n");
        status = 1;
    }

    if (equalsMatrixD(a, b, 0)) {
        printf("FAILED Test_equalsMatrixD (a, b, 0)\n");
        status = 1;
    }

    if (equalsMatrixD(a, b, 0.1)) {
        printf("FAILED Test_equalsMatrixD (a, b, 0.1)\n");
        status = 1;
    }

    if (!equalsMatrixD(a, b, 0.5)) {
        printf("FAILED Test_equalsMatrixD (a, b, 0.5)\n");
        status = 1;
    }

    if (!equalsMatrixD(a, b, 1)) {
        printf("FAILED Test_equalsMatrixD (a, b, 1)\n");
        status = 1;
    }

    if (equalsMatrixD(a, c, 0)) {
        printf("FAILED Test_equalsMatrixD (a, c, 0)\n");
        status = 1;
    }

    if (equalsMatrixD(a, c, 0.5)) {
        printf("FAILED Test_equalsMatrixD (a, c, 0.5)\n");
        status = 1;
    }

    if (!equalsMatrixD(a, c, 1)) {
        printf("FAILED Test_equalsMatrixD (a, c, 1)\n");
        status = 1;
    }

    if (!equalsMatrixD(a, c, 2)) {
        printf("FAILED Test_equalsMatrixD (a, c, 2)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_equalsMatrixD\n");

    return status;
}

static int Test_addMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * j);
            b.matrix[i][j] = (double) (i + j);
        }
    }

    MatrixD c = addMatrixD(a, b);
    if (c.rows != dim || c.cols != dim) {
        printf("FAILED Test_addMatrixD (size not equal)\n");
        status = 1;
    }
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (c.matrix[i][j] != (double) (i * j + i + j)) {
                printf("FAILED Test_addMatrixD (values not equal)\n");
                status = 1;
            }
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_addMatrixD\n");

    return status;
}