#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

typedef struct {
    size_t rows;
    size_t cols;
    double** matrix;
} MatrixD;


MatrixD newMatrixD(size_t rows, size_t cols);
MatrixD identityMatrixD(size_t rows);
MatrixD copyMatrixD(MatrixD a);
void displayMatrixD(MatrixD m);
bool equalsMatrixD(MatrixD a, MatrixD b);
MatrixD addMatrixD(MatrixD a, MatrixD b);
MatrixD multiplyMatrixD(MatrixD a, MatrixD b);
MatrixD scaleMatrixD(double k, MatrixD a);
double slowDeterminantMatrixD(MatrixD a);
MatrixD appendMatrixD(MatrixD a, MatrixD b);
MatrixD echelonMatrixD(MatrixD a);
void _swapRowMatrixD(MatrixD *m, size_t row1, size_t row2);
void _rowAdditionMatrixD(MatrixD *m, double k, size_t row1, size_t row2);
double fastDeterminantMatrixD(MatrixD m);
MatrixD reducedEchelonMatrixD(MatrixD m);
void _scaleRowMatrixD(MatrixD *m, double k, size_t row);
MatrixD inverseMatrixD(MatrixD m);
MatrixD _submatrixMatrixD(MatrixD m, size_t row, size_t col, size_t height, size_t width);
MatrixD transposeMatrixD(MatrixD m);
void qrDecompositionMatrixD(MatrixD m, MatrixD *q, MatrixD *r);
void qrAlgorithmMatrixD(MatrixD m, size_t iterations, MatrixD *d, MatrixD *p);


MatrixD newMatrixD(size_t rows, size_t cols) {
    double** mat = malloc(sizeof(double*) * rows);
    for (size_t i = 0; i < rows; ++i) {
        mat[i] = calloc(cols, sizeof(double));
    }
    MatrixD m = {rows, cols, mat};
    return m;
}

MatrixD identityMatrixD(size_t rows) {
    MatrixD m = newMatrixD(rows, rows);
    for (size_t i = 0; i < rows; ++i) {
        m.matrix[i][i] = 1;
    }
    return m;
}

MatrixD copyMatrixD(MatrixD a) {
    MatrixD b = newMatrixD(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            b.matrix[i][j] = a.matrix[i][j];
        }
    }
    return b;
}

void displayMatrixD(MatrixD m) {
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            printf("%f ", m.matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

bool equalsMatrixD(MatrixD a, MatrixD b) {
    if (a.rows != b.rows || a.cols != b.cols) return false;

    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            if (a.matrix[i][j] - b.matrix[i][j] <= DBL_EPSILON) return false;
        }
    }

    return true;
}

MatrixD addMatrixD(MatrixD a, MatrixD b) {
    // Check rows and columns match
    if (a.rows != b.rows || a.cols != b.cols) {
        fprintf(stderr, "Error: Could not add matrices of shape (%lu, %lu) and (%lu, %lu)\n", a.rows, a.cols, b.rows, b.cols);
        exit(1);
    }

    MatrixD m = newMatrixD(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            m.matrix[i][j] = a.matrix[i][j] + b.matrix[i][j];
        }
    }

    return m;
}

MatrixD multiplyMatrixD(MatrixD a, MatrixD b) {
    //Check a.cols == b.rows
    if (a.cols != b.rows) {
        fprintf(stderr, "Error: Could not multiply matrices of shape (%lu, %lu) and (%lu, %lu)\n", a.rows, a.cols, b.rows, b.cols);
        exit(1);
    }

    MatrixD m = newMatrixD(a.rows, b.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < b.cols; ++j) {
            for (size_t k = 0; k < a.cols; ++k) {
                m.matrix[i][j] += a.matrix[i][k] * b.matrix[k][j];
            }
        }
    }

    return m;
}

MatrixD scaleMatrixD(double k, MatrixD a) {
    MatrixD m = newMatrixD(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            m.matrix[i][j] = k * a.matrix[i][j];
        }
    }

    return m;
}

double slowDeterminantMatrixD(MatrixD a) {
    if (a.rows != a.cols) {
        fprintf(stderr, "Error: Cannot compute determinant of non-square matrix with shape (%lu, %lu)\n", a.rows, a.cols);
        exit(1);
    }

    if (a.rows == 1) return a.matrix[0][0];

    if (a.rows == 2) return a.matrix[0][0] * a.matrix[1][1] - a.matrix[0][1] * a.matrix[1][0];

    // Expand along the first column
    double sum = 0;
    int unit_factor = 1;
    MatrixD m = newMatrixD(a.rows - 1, a.rows - 1);
    for (size_t i = 0; i < a.rows; ++i) {
        size_t m_row = 0;
        for (size_t a_row = 0; a_row < a.rows; ++a_row) {
            if (i == a_row) continue;
            for (size_t col = 1; col < a.cols; ++col) {
                m.matrix[m_row][col - 1] = a.matrix[a_row][col];
            }
            m_row++;
        }
        double f = slowDeterminantMatrixD(m);
        sum += unit_factor * a.matrix[i][0] * f;
        unit_factor *= -1;
    }

    return sum;
}

MatrixD appendMatrixD(MatrixD a, MatrixD b) {
    if (a.rows != b.rows) {
        fprintf(stderr, "Error: Cannot append matrices of shape (%lu, %lu) and (%lu, %lu)\n", a.rows, a.cols, b.rows, b.cols);
        exit(1);
    }

    MatrixD m = newMatrixD(a.rows, a.cols + b.cols);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            m.matrix[i][j] = a.matrix[i][j];
        }
        for (size_t j = 0; j < b.cols; ++j) {
            m.matrix[i][j + a.cols] = b.matrix[i][j];
        }
    }

    return m;
}

MatrixD echelonMatrixD(MatrixD m) {
    MatrixD a = copyMatrixD(m);

    size_t pivots = 0;
    for (size_t col = 0; col < a.cols; ++col) {
        size_t row = pivots;
        while (a.matrix[row][col] == 0) {
            row++;
            if (row >= a.rows) {
                row = a.rows - 1;
                break;
            }
        }
        if (a.matrix[row][col] != 0) {
            _swapRowMatrixD(&a, pivots, row);
            for (size_t i = pivots + 1; i < a.rows; ++i) {
                _rowAdditionMatrixD(&a, -a.matrix[i][col] / a.matrix[pivots][col], pivots, i);
            }
            pivots++;
        }
        
        if (pivots == a.rows) break;
    }

    return a;
}

void _swapRowMatrixD(MatrixD *m, size_t row1, size_t row2) {
    if (row1 >= m->rows || row2 >= m->rows) {
        fprintf(stderr, "Error: Could not swap row %lu with row %lu in matrix of shape (%lu, %lu)\n", row1, row2, m->rows, m->cols);
        exit(1);
    }

    if (row1 == row2) return;

    double temp;
    for (size_t j = 0; j < m->cols; ++j) {
        temp = m->matrix[row1][j];
        m->matrix[row1][j] = m->matrix[row2][j];
        m->matrix[row2][j] = temp;
    }
}

// replace row2 with k * row1 + row2
void _rowAdditionMatrixD(MatrixD *m, double k, size_t row1, size_t row2) {
    if (row1 >= m->rows || row2 >= m->rows) {
        fprintf(stderr, "Error: Could not add row %lu to row %lu in matrix of shape (%lu, %lu)\n", row1, row2, m->rows, m->cols);
        exit(1);
    }

    for (size_t j = 0; j < m->cols; ++j) {
        m->matrix[row2][j] = m->matrix[row2][j] + k * m->matrix[row1][j];
    }
}

double fastDeterminantMatrixD(MatrixD m) {
    if (m.rows != m.cols) {
        fprintf(stderr, "Error: Cannot compute determinant of non-square matrix with shape (%lu, %lu)\n", m.rows, m.cols);
        exit(1);
    }

    if (m.rows == 1) return m.matrix[0][0];

    if (m.rows == 2) return m.matrix[0][0] * m.matrix[1][1] - m.matrix[0][1] * m.matrix[1][0];

    // compute echelon form of matrix
    int swapSign = 1;
    MatrixD a = copyMatrixD(m);

    size_t pivots = 0;
    for (size_t col = 0; col < a.cols; ++col) {
        size_t row = pivots;
        while (a.matrix[row][col] == 0) {
            row++;
            if (row >= a.rows) {
                row = a.rows - 1;
                break;
            }
        }
        if (a.matrix[row][col] != 0) {
            if (pivots != row) {
                _swapRowMatrixD(&a, pivots, row);
                swapSign *= -1;
            }
            for (size_t i = pivots + 1; i < a.rows; ++i) {
                _rowAdditionMatrixD(&a, -a.matrix[i][col] / a.matrix[pivots][col], pivots, i);
            }
            pivots++;
        }
        
        if (pivots == a.rows) break;
    }

    double det = swapSign;
    for (size_t i = 0; i < a.rows; ++i) {
        det *= a.matrix[i][i];
    }

    return det;
}

MatrixD reducedEchelonMatrixD(MatrixD m) {
    MatrixD a = copyMatrixD(m);

    size_t pivots = 0;
    for (size_t col = 0; col < a.cols; ++col) {
        size_t row = pivots;
        while (a.matrix[row][col] == 0) {
            row++;
            if (row >= a.rows) {
                row = a.rows - 1;
                break;
            }
        }
        if (a.matrix[row][col] != 0) {
            _swapRowMatrixD(&a, pivots, row);
            _scaleRowMatrixD(&a, 1 / a.matrix[pivots][col], pivots);
            for (size_t i = 0; i < a.rows; ++i) {
                if (i == pivots) continue;
                _rowAdditionMatrixD(&a, -a.matrix[i][col], pivots, i);
            }
            pivots++;
        }
        
        if (pivots == a.rows) break;
    }

    return a;
}

void _scaleRowMatrixD(MatrixD *m, double k, size_t row) {
    if (row >= m->rows) {
        fprintf(stderr, "Error: Cannot scale row %lu in matrix of shape (%lu, %lu)\n", row, m->rows, m->cols);
        exit(1);
    }

    for (size_t j = 0; j < m->cols; ++j) {
        m->matrix[row][j] *= k;
    }
}

MatrixD inverseMatrixD(MatrixD m) {
    if (m.rows != m.cols) {
        fprintf(stderr, "Error: Cannot compute inverse of non-square matrix with shape (%lu, %lu)\n", m.rows, m.cols);
        exit(1);
    }

    MatrixD a = appendMatrixD(m, identityMatrixD(m.rows));
    a = reducedEchelonMatrixD(a);
    return _submatrixMatrixD(a, 0, m.rows, m.rows, m.rows);
}

MatrixD _submatrixMatrixD(MatrixD m, size_t row, size_t col, size_t height, size_t width) {
    if (row + height - 1 >= m.rows || col + width - 1 >= m.cols) {
        fprintf(stderr, "Error: Cannot take submatrix with entry [%lu][%lu] which exceeds bounds of shape (%lu, %lu)\n", row + height - 1, col + width - 1, m.rows, m.cols);
        exit(1);
    }

    MatrixD a = newMatrixD(height, width);
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            a.matrix[i][j] = m.matrix[row + i][col + j];
        }
    }

    return a;
}

MatrixD transposeMatrixD(MatrixD m) {
    MatrixD a = newMatrixD(m.cols, m.rows);

    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            a.matrix[i][j] = m.matrix[j][i];
        }
    }

    return a;
}

double normMatrixD(MatrixD m) {
    double sum = 0;
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            sum += m.matrix[i][j] * m.matrix[i][j];
        }
    }

    return sqrt(sum);
}

void qrDecompositionMatrixD(MatrixD m, MatrixD *q, MatrixD *r) {
    size_t t = (m.rows > m.cols) ? m.cols : m.rows; // minimum

    MatrixD mprime = copyMatrixD(m);
    *q = identityMatrixD(m.rows);
    *r = m;
    for (size_t i = 0; i < t; ++i) {
        MatrixD x = _submatrixMatrixD(mprime, 0, 0, mprime.rows, 1);

        double alpha = normMatrixD(x);
        if (x.matrix[0][0] >= 0) alpha *= -1;

        MatrixD e1 = newMatrixD(mprime.rows, 1);
        e1.matrix[0][0] = 1;

        MatrixD u = addMatrixD(x, scaleMatrixD(-alpha, e1));
        MatrixD v = scaleMatrixD(1 / normMatrixD(u), u);

        MatrixD qprime = addMatrixD(identityMatrixD(mprime.rows), scaleMatrixD(-2, multiplyMatrixD(v, transposeMatrixD(v))));

        MatrixD Qi = newMatrixD(m.rows, m.rows);
        for (size_t j = 0; j < i; ++j) {
            Qi.matrix[j][j] = 1;
        }
        for (size_t j = i; j < m.rows; ++j) {
            for (size_t k = i; k < m.rows; ++k) {
                Qi.matrix[j][k] = qprime.matrix[j-i][k-i];
            }
        }

        *q = multiplyMatrixD(*q, Qi);
        *r = multiplyMatrixD(Qi, *r);

        mprime = _submatrixMatrixD(*r, i + 1, i + 1, m.rows - i - 1, m.rows - i - 1);
    }
}

// TODO: Add order reduction after convergence check
void qrAlgorithmMatrixD(MatrixD m, size_t iterations, MatrixD *d, MatrixD *p) {
    if (m.rows != m.cols) {
        fprintf(stderr, "Error: Cannot perform QR Algorithm on non-square matrix of shape (%lu, %lu)\n", m.rows, m.cols);
        exit(1);
    }

    MatrixD *q = malloc(sizeof(MatrixD));
    MatrixD *r = malloc(sizeof(MatrixD));

    *d = m;
    *p = identityMatrixD(m.rows);
    for (size_t i = 0; i < iterations; ++i) {
        double lambda = d->matrix[d->rows - 1][d->cols - 1]; // Rayleigh shift
        
        qrDecompositionMatrixD(addMatrixD(*d, scaleMatrixD(-lambda, identityMatrixD(d->rows))), q, r);
        *d = addMatrixD(multiplyMatrixD(*r, *q), scaleMatrixD(lambda, identityMatrixD(d->rows)));
        *p = multiplyMatrixD(*p, *q);
    }

    free(q);
    free(r);
}

// Solve Ax = b, the columns of N form a basis for the null space of A
// Function returns the dimension of the solution space, -1 if there is no solution
int solveLinearMatrixD(MatrixD A, MatrixD b, MatrixD *x, MatrixD *N) {
    if (A.rows != b.rows || b.cols != 1) {
        fprintf(stderr, "Error: Cannot solve Ax=b with A of shape (%lu, %lu) and b of shape (%lu, %lu)\n", A.rows, A.cols, b.rows, b.cols);
        exit(1);
    }

    MatrixD augmented = appendMatrixD(A, b);
    augmented = reducedEchelonMatrixD(augmented);

    MatrixD sol = newMatrixD(A.cols, 1);
    MatrixD null;
    size_t pivots = augmented.rows;
    bool foundPivot;
    size_t lastNullColumn;
    size_t lastPivotColumn = augmented.cols - 1;
    size_t *freeColumns = malloc(sizeof(size_t) * (A.cols));
    for (int i = augmented.rows - 1; i >= 0; --i) {
        // look for pivot
        foundPivot = false;
        for (size_t j = i; j < augmented.cols; ++j) {
            if (j == augmented.cols - 1) {
                if (augmented.matrix[i][j] != 0 && !foundPivot) {
                    // inconsistent
                    return -1;
                } else if (!foundPivot) {
                    // row of all zeros
                    pivots--;
                    break;
                } else {
                    break;
                }
            }
            if (augmented.matrix[i][j] != 0) {
                if (!foundPivot) {
                    // found pivot
                    foundPivot = true;
                    if (i == pivots - 1) {
                        // first row from the bottom with a pivot
                        // null space is A.cols - pivots
                        null = newMatrixD(A.cols, A.cols - pivots);
                        lastNullColumn = A.cols - pivots;
                    }
                    sol.matrix[j][0] = augmented.matrix[i][augmented.cols - 1];
                    for (size_t k = lastPivotColumn - 1; k > j; --k) {
                        lastNullColumn--;
                        null.matrix[k][lastNullColumn] = 1;
                        freeColumns[k] = lastNullColumn;
                    }
                    lastPivotColumn = j;
                } else {
                    null.matrix[lastPivotColumn][freeColumns[j]] = -augmented.matrix[i][j];
                }
            }
        }
    }

    *x = sol;
    *N = null;

    return null.cols;
}