#include <assert.h>
#include "la.h"

static int Test_newMatrixD();
static int Test_identityMatrixD();
static int Test_copyMatrixD();

static int Test_equalsMatrixD();

static int Test_addMatrixD();
static int Test_subtractMatrixD();
static int Test_multiplyMatrixD();
static int Test_addMultiplyMatrixD();
static int Test_scaleMatrixD();
static int Test_innerProductMatrixD();
static int Test_outerProductMatrixD();

static int Test_slowDeterminantMatrixD();
static int Test_fastDeterminantMatrixD();

static int Test_appendMatrixD();

static int Test_inverseMatrixD();
static int Test_transposeMatrixD();
static int Test_normMatrixD();
static int Test_solveLinearMatrixD();

static int Test_strassenMultiplyMatrixD();

int main() {
    int status = 0;

    status |= Test_newMatrixD();
    status |= Test_identityMatrixD();
    status |= Test_copyMatrixD();

    status |= Test_equalsMatrixD();

    status |= Test_addMatrixD();
    status |= Test_subtractMatrixD();
    status |= Test_multiplyMatrixD();
    status |= Test_addMultiplyMatrixD();
    status |= Test_scaleMatrixD();
    status |= Test_innerProductMatrixD();
    status |= Test_outerProductMatrixD();

    status |= Test_slowDeterminantMatrixD();
    status |= Test_fastDeterminantMatrixD();

    status |= Test_appendMatrixD();

    status |= Test_inverseMatrixD();
    status |= Test_transposeMatrixD();
    status |= Test_normMatrixD();
    status |= Test_solveLinearMatrixD();

    status |= Test_strassenMultiplyMatrixD();

    printf("----------------------------------------\n");

    if (status == 0) {
        printf("ALL TESTS PASSED\n");
    } else if (status == 1) {
        printf("SOME TESTS FAILED\n");
    }

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
                break;
            }
        }
        if (status == 1) break;
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
                break;
            } else if (i == j && a.matrix[i][j] != 1) {
                printf("FAILED Test_identityMatrixD (diagonal value not one)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
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
                break;
            }
            a.matrix[i][j] += 1;
            if (a.matrix[i][j] == b.matrix[i][j]) {
                printf("FAILED Test_copyMatrixD (copy not independent)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
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

    MatrixD c = newMatrixD(dim, dim);
    addMatrixD(&c, &a, &b);

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (c.matrix[i][j] != (double) (i * j + i + j)) {
                printf("FAILED Test_addMatrixD (values not equal)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_addMatrixD\n");

    return status;
}

static int Test_subtractMatrixD() {
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

    MatrixD c = newMatrixD(dim, dim);
    subtractMatrixD(&c, &a, &b);

    for (int i = 0; i < (int) dim; ++i) {
        for (int j = 0; j < (int) dim; ++j) {
            if (c.matrix[i][j] != (double) (i * j - i - j)) {
                printf("FAILED Test_subtractMatrixD (values not equal)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_subtractMatrixD\n");

    return status;
}

static int Test_multiplyMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = identityMatrixD(dim);
    MatrixD c = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
        }
    }
    
    multiplyMatrixD(&c, &a, &b);
    if (!equalsMatrixD(a, c, 0)) {
        printf("FAILED Test_multiplyMatrixD (a * I == a)\n");
        status = 1;
    }

    multiplyMatrixD(&c, &b, &a);
    if (!equalsMatrixD(a, c, 0)) {
        printf("FAILED Test_multiplyMatrixD (I * a == a)\n");
        status = 1;
    }

    multiplyMatrixD(&c, &b, &b);
    if (!equalsMatrixD(b, c, 0)) {
        printf("FAILED Test_multiplyMatrixD (I * I == I)\n");
        status = 1;
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 1;
        }
    }

    multiplyMatrixD(&a, &a, &a);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (a.matrix[i][j] != dim) {
                printf("FAILED Test_multiplyMatrixD (multiply in-place)\n");
                status = 1;
            }
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_multiplyMatrixD\n");

    return status;
}

static int Test_addMultiplyMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = identityMatrixD(dim);
    MatrixD c;
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
        }
    }
    
    c = newMatrixD(dim, dim);
    addMultiplyMatrixD(&c, &a, &b);
    if (!equalsMatrixD(a, c, 0)) {
        printf("FAILED Test_addMultiplyMatrixD (a * I == a)\n");
        status = 1;
    }
    freeMatrixD(&c);

    c = newMatrixD(dim, dim);
    addMultiplyMatrixD(&c, &b, &a);
    if (!equalsMatrixD(a, c, 0)) {
        printf("FAILED Test_addMultiplyMatrixD (I * a == a)\n");
        status = 1;
    }
    freeMatrixD(&c);

    c = newMatrixD(dim, dim);
    addMultiplyMatrixD(&c, &b, &b);
    if (!equalsMatrixD(b, c, 0)) {
        printf("FAILED Test_addMultiplyMatrixD (I * I == I)\n");
        status = 1;
    }
    freeMatrixD(&c);

    c = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 1;
        }
    }

    addMultiplyMatrixD(&c, &a, &a);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (c.matrix[i][j] != dim) {
                printf("FAILED Test_addMultiplyMatrixD (a * a)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
    }

    addMultiplyMatrixD(&c, &a, &a);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            if (c.matrix[i][j] != 2 * dim) {
                printf("FAILED Test_addMultiplyMatrixD (a * a + a * a)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_multiplyMatrixD\n");

    return status;
}

static int Test_scaleMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
        }
    }

    MatrixD y = newMatrixD(dim, dim);
    for (int k = -5; k <= 5; ++k) {
        scaleMatrixD(&y, k, &a);
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                if (y.matrix[i][j] != (double) (k * (int) (i * i * j))) {
                    printf("FAILED Test_scaleMatrixD (incorrect value)\n");
                    status = 1;
                    break;
                }
                if (a.matrix[i][j] != (double) (i * i * j)) {
                    printf("FAILED Test_scaleMatrixD (original matrix changed)\n");
                    status = 1;
                    break;
                }
            }
            if (status == 1) break;
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&y);

    if (status == 0) printf("PASSED Test_scaleMatrixD\n");

    return status;
}

static int Test_innerProductMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    MatrixD c1 = newMatrixD(dim, dim);
    MatrixD c2 = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
            b.matrix[i][j] = (double) (2 * i + j);
        }
    }

    MatrixD at = transposeMatrixD(a);
    innerProductMatrixD(&c1, &a, &b);
    multiplyMatrixD(&c2, &at, &b);
    if (!equalsMatrixD(c1, c2, DBL_EPSILON)) {
        printf("FAILED Test_innerProductMatrixD (incorrect value)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c1);
    freeMatrixD(&c2);
    freeMatrixD(&at);

    if (status == 0) printf("PASSED Test_innerProductMatrixD\n");

    return status;
}

static int Test_outerProductMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    MatrixD c1 = newMatrixD(dim, dim);
    MatrixD c2 = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
            b.matrix[i][j] = (double) (2 * i + j);
        }
    }

    MatrixD bt = transposeMatrixD(b);
    outerProductMatrixD(&c1, &a, &b);
    multiplyMatrixD(&c2, &a, &bt);
    if (!equalsMatrixD(c1, c2, DBL_EPSILON)) {
        printf("FAILED Test_outerProductMatrixD (incorrect value)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c1);
    freeMatrixD(&c2);
    freeMatrixD(&bt);

    if (status == 0) printf("PASSED Test_outerProductMatrixD\n");

    return status;
}

static int Test_slowDeterminantMatrixD() {
    int status = 0;
    size_t dim = 5;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = identityMatrixD(dim);
    MatrixD c = identityMatrixD(dim);
    scaleMatrixD(&c, 2, &c);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 1;
        }
    }

    if (slowDeterminantMatrixD(a) != 0) {
        printf("FAILED Test_slowDeterminantMatrixD (det(a) == 0)\n");
        status = 1;
    }

    if (slowDeterminantMatrixD(b) != 1) {
        printf("FAILED Test_slowDeterminantMatrixD (det(b) == 1)\n");
        status = 1;
    }

    if (slowDeterminantMatrixD(c) != pow(2, (double) dim)) {
        printf("FAILED Test_slowDeterminantMatrixD (det(c) == 2^dim)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_slowDeterminantMatrixD\n");

    return status;
}

static int Test_fastDeterminantMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = identityMatrixD(dim);
    MatrixD c = identityMatrixD(dim);
    scaleMatrixD(&c, 2, &c);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 1;
        }
    }

    if (fastDeterminantMatrixD(a) != 0) {
        printf("FAILED Test_fastDeterminantMatrixD (det(a) == 0)\n");
        status = 1;
    }

    if (fastDeterminantMatrixD(b) != 1) {
        printf("FAILED Test_fastDeterminantMatrixD (det(b) == 1)\n");
        status = 1;
    }

    if (fastDeterminantMatrixD(c) != pow(2, (double) dim)) {
        printf("FAILED Test_fastDeterminantMatrixD (det(c) == 2^dim)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_fastDeterminantMatrixD\n");

    return status;
}

static int Test_appendMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
            b.matrix[i][j] = (double) (2 * i + j);
        }
    }

    MatrixD c = appendMatrixD(a, b);

    if (c.rows != dim || c.cols != 2 * dim) {
        printf("FAILED Test_appendMatrixD (wrong size)\n");
        status = 1;
    }
    for (size_t i = 0; i < c.rows; ++i) {
        for (size_t j = 0; j < c.cols; ++j) {
            if ((j < dim && c.matrix[i][j] != (double) (i * i * j)) || (j >= dim && c.matrix[i][j] != (double) (2 * i + j - dim))) {
                printf("FAILED Test_appendMatrixD (incorrect value)\n");
                status = 1;
                break;
            }
            if (j < dim && (a.matrix[i][j] != (double) (i * i * j) || b.matrix[i][j] != (double) (2 * i + j))) {
                printf("FAILED Test_appendMatrixD (inputs changed)\n");
                status = 1;
                break;
            }
        }
        if (status == 1) break;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    if (status == 0) printf("PASSED Test_appendMatrixD\n");

    return status;
}

static int Test_inverseMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD I = identityMatrixD(dim);
    MatrixD a = newMatrixD(dim, dim);
    for (int i = 0; i < (int) dim; ++i) {
        for (int j = 0; j < (int) dim; ++j) {
            a.matrix[i][j] = exp(-((double) (i - j) * (double) (i - j)));
        }
    }

    MatrixD b = inverseMatrixD(I);
    if (!equalsMatrixD(I, b, 0)) {
        printf("FAILED Test_inverseMatrixD (inverse of identity)\n");
        status = 1;
    }

    freeMatrixD(&b);
    b = inverseMatrixD(a);
    MatrixD c = newMatrixD(dim, dim);
    MatrixD d = newMatrixD(dim, dim);
    multiplyMatrixD(&c, &a, &b);
    multiplyMatrixD(&d, &b, &a);
    if (!equalsMatrixD(I, c, 1e-10)) {
        printf("FAILED Test_inverseMatrixD (right inverse)\n");
        status = 1;
    }
    if (!equalsMatrixD(I, d, 1e-10)) {
        printf("FAILED Test_inverseMatrixD (left inverse)\n");
        status = 1;
    }

    freeMatrixD(&I);
    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);
    freeMatrixD(&d);

    if (status == 0) printf("PASSED Test_inverseMatrixD\n");

    return status;
}

static int Test_transposeMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, 1);
    for (size_t i = 0; i < dim; ++i) {
        a.matrix[i][0] = (double) (i * i + i + 1);
    }

    MatrixD b = transposeMatrixD(a);

    if (b.rows != 1 || b.cols != dim) {
        printf("FAILED Test_transposeMatrixD (wrong size)\n");
        status = 1;
    }

    for (size_t j = 0; j < dim; ++j) {
        if (b.matrix[0][j] != (double) (j * j + j + 1)) {
            printf("FAILED Test_transposeMatrixD (wrong values)\n");
            status = 1;
            break;
        }
    }

    freeMatrixD(&a);
    freeMatrixD(&b);

    if (status == 0) printf("PASSED Test_transposeMatrixD\n");

    return status;
}

static int Test_normMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = ((i + j) % 2 == 0) ? 1 : -1;
        }
    }

    if (normMatrixD(a, ONE_NORM) != dim * dim) {
        printf("FAILED Test_normMatrixD (one_norm)\n");
        status = 1;
    }

    if (normMatrixD(a, TWO_NORM) != dim) {
        printf("FAILED Test_normMatrixD (two_norm)\n");
        status = 1;
    }

    if (normMatrixD(a, INF_NORM) != 1) {
        printf("FAILED Test_normMatrixD (inf_norm)\n");
        status = 1;
    }

    freeMatrixD(&a);

    if (status == 0) printf("PASSED Test_normMatrixD\n");

    return status;
}

static int Test_solveLinearMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 1;
        }
    }
    MatrixD e1 = newMatrixD(dim, 1);
    e1.matrix[0][0] = 1;

    MatrixD *x = malloc(sizeof(MatrixD));
    MatrixD *N = malloc(sizeof(MatrixD));
    int result;

    result = solveLinearMatrixD(a, e1, x, N);
    if (result != -1) {
        printf("FAILED Test_solveLinearMatrixD (no solution)\n");
        status = 1;
    }

    MatrixD I = identityMatrixD(dim);
    result = solveLinearMatrixD(I, e1, x, N);
    if (result != 0) {
        printf("FAILED Test_solveLinearMatrixD (one solution)\n");
        status = 1;
    }
    if (!equalsMatrixD(e1, *x, 0)) {
        printf("FAILED Test_solveLinearMatrixD (wrong solution)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&e1);
    freeMatrixD(&I);
    freeMatrixD(x);
    freeMatrixD(N);
    free(x);
    free(N);

    if (status == 0) printf("PASSED Test_solveLinearMatrixD\n");

    return status;
}

static int Test_strassenMultiplyMatrixD() {
    int status = 0;
    size_t dim = 10;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = (double) (i * i * j);
            b.matrix[i][j] = (double) (2 * i + j);
        }
    }

    MatrixD c1 = newMatrixD(dim, dim);
    MatrixD c2 = newMatrixD(dim, dim);
    multiplyMatrixD(&c1, &a, &b);
    strassenMultiplyMatrixD(&c2, &a, &b);

    if (!equalsMatrixD(c1, c2, DBL_EPSILON)) {
        printf("FAILED Test_strassenMultiplyMatrixD (incorrect value)\n");
        status = 1;
    }

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c1);
    freeMatrixD(&c2);

    if (status == 0) printf("PASSED Test_strassenMultiplyMatrixD\n");

    return status;
}