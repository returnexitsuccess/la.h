#include <stdlib.h>
#include <time.h>
#include "la.h"

double benchmarkMultiply(size_t dim, size_t iterations);
double benchmarkStrassenMultiply(size_t dim, size_t iterations);
double benchmarkDeterminant(size_t dim, size_t iterations);
double benchmarkInverse(size_t dim, size_t iterations);


int main() {
    srand(0);
    double t;

    printf("# Benchmarks\n\n");

    printf("| Test                                     | Time Taken |\n");
    printf("| ---------------------------------------- | ---------- |\n");

    t = benchmarkMultiply(10, 100000);
    printf("| 10x10 multiply, 100,000 iterations       | %fs  |\n", t);

    t = benchmarkMultiply(1000, 1);
    printf("| 1000x1000 multiply, 1 iteration          | %fs  |\n", t);

    t = benchmarkStrassenMultiply(1000, 1);
    printf("| 1000x1000 strassen multiply, 1 iteration | %fs  |\n", t);

    t = benchmarkDeterminant(100, 1000);
    printf("| 100x100 determinant, 1000 iterations     | %fs  |\n", t);

    t = benchmarkDeterminant(1000, 1);
    printf("| 1000x1000 determinant, 1 iteration       | %fs  |\n", t);

    t = benchmarkInverse(100, 100);
    printf("| 100x100 inverse, 100 iterations          | %fs  |\n", t);

    t = benchmarkInverse(1000, 1);
    printf("| 1000x1000 inverse, 1 iteration           | %fs  |\n", t);

    return 0;
}

double benchmarkMultiply(size_t dim, size_t iterations) {
    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    MatrixD c = newMatrixD(dim, dim);

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
            b.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    clock_t t0, t1;

    t0 = clock();

    for (size_t k = 0; k < iterations; ++k) {
        multiplyMatrixD(&c, &a, &b);
    }

    t1 = clock();

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    return (double) (t1 - t0) / CLOCKS_PER_SEC;
}

double benchmarkStrassenMultiply(size_t dim, size_t iterations) {
    MatrixD a = newMatrixD(dim, dim);
    MatrixD b = newMatrixD(dim, dim);
    MatrixD c = newMatrixD(dim, dim);

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
            b.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    clock_t t0, t1;

    t0 = clock();

    for (size_t k = 0; k < iterations; ++k) {
        strassenMultiplyMatrixD(&c, &a, &b);
    }

    t1 = clock();

    freeMatrixD(&a);
    freeMatrixD(&b);
    freeMatrixD(&c);

    return (double) (t1 - t0) / CLOCKS_PER_SEC;
}

double benchmarkDeterminant(size_t dim, size_t iterations) {
    MatrixD a = newMatrixD(dim, dim);

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    clock_t t0, t1;

    t0 = clock();

    for (size_t k = 0; k < iterations; ++k) {
        fastDeterminantMatrixD(a);
    }

    t1 = clock();

    freeMatrixD(&a);

    return (double) (t1 - t0) / CLOCKS_PER_SEC;
}

double benchmarkInverse(size_t dim, size_t iterations) {
    double sum = 0;

    MatrixD a = newMatrixD(dim, dim);
    MatrixD b;

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            a.matrix[i][j] = 2 * ((double) rand() / RAND_MAX) - 1;
        }
    }

    clock_t t0, t1;

    for (size_t k = 0; k < iterations; ++k) {
        t0 = clock();
        b = inverseMatrixD(a);
        t1 = clock();
        sum += (double) (t1 - t0);
        freeMatrixD(&b);
    }

    freeMatrixD(&a);

    return sum / CLOCKS_PER_SEC;
}