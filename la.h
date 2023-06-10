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

typedef struct {
    size_t rows;
    size_t cols;
    size_t rowOffset;
    size_t colOffset;
    double** matrix;
} _BlockMatrixD;


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
int solveLinearMatrixD(MatrixD A, MatrixD b, MatrixD *x, MatrixD *N);
MatrixD strassenMultiplyMatrixD(MatrixD a, MatrixD b);
_BlockMatrixD _toBlockMatrixD(MatrixD m, size_t rows, size_t cols, size_t rowOffset, size_t colOffset);
MatrixD _fromBlockMatrixD(_BlockMatrixD mBlock);
_BlockMatrixD _addBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
_BlockMatrixD _subtractBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
void _updateAddBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b);
void _updateSubtractBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b);
_BlockMatrixD _strassenMultiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
_BlockMatrixD _multiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);


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
        MatrixD qprime;
        if (fabs(alpha) > 1e-100) {
            if (x.matrix[0][0] >= 0) alpha *= -1;

            MatrixD e1 = newMatrixD(mprime.rows, 1);
            e1.matrix[0][0] = 1;

            MatrixD u = addMatrixD(x, scaleMatrixD(-alpha, e1));
            MatrixD v = scaleMatrixD(1 / normMatrixD(u), u);

            qprime = addMatrixD(identityMatrixD(mprime.rows), scaleMatrixD(-2, multiplyMatrixD(v, transposeMatrixD(v))));
        } else {
            qprime = identityMatrixD(mprime.rows);
        }

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

void qrAlgorithmMatrixD(MatrixD m, size_t iterations, MatrixD *d, MatrixD *p) {
    if (m.rows != m.cols) {
        fprintf(stderr, "Error: Cannot perform QR Algorithm on non-square matrix of shape (%lu, %lu)\n", m.rows, m.cols);
        exit(1);
    }

    MatrixD *q = malloc(sizeof(MatrixD));
    MatrixD *r = malloc(sizeof(MatrixD));

    *d = m;
    *p = identityMatrixD(m.rows);

    size_t blocksize = 1;
    for (size_t i = 0; i < iterations; ++i) {

        bool converged = true;

        // check if complex eigenvalues in lower right 2x2 block
        double trace = d->matrix[d->rows - 2][d->cols - 2] + d->matrix[d->rows - 1][d->cols - 1];
        double det = d->matrix[d->rows - 2][d->cols - 2] * d->matrix[d->rows - 1][d->cols - 1] - d->matrix[d->rows - 2][d->cols - 1] * d->matrix[d->rows - 1][d->cols - 2];
        if (trace * trace < 4 * det) {
            // complex eigenvalues
            // double implicit shift strategy
            MatrixD M = addMatrixD(multiplyMatrixD(*d, *d), addMatrixD(scaleMatrixD(-trace, *d), scaleMatrixD(det, identityMatrixD(d->rows))));
            qrDecompositionMatrixD(M, q, r);
            *d = multiplyMatrixD(transposeMatrixD(*q), multiplyMatrixD(*d, *q));
            *p = multiplyMatrixD(*p, *q);

            blocksize = 2;

            for (size_t j = 0; j < d->cols - 2; ++j) {
                converged = converged && (fabs(d->matrix[d->rows - 2][j]) < 1e-100);
                converged = converged && (fabs(d->matrix[d->rows - 1][j]) < 1e-100);
            }

            if (converged) {
                //printf("2: Exited early at iteration %lu for stage %lu\n", i, d->rows);
                break;
            }
        } else {
            double lambda = d->matrix[d->rows - 1][d->cols - 1]; // Rayleigh shift
            
            qrDecompositionMatrixD(addMatrixD(*d, scaleMatrixD(-lambda, identityMatrixD(d->rows))), q, r);
            *d = addMatrixD(multiplyMatrixD(*r, *q), scaleMatrixD(lambda, identityMatrixD(d->rows)));
            *p = multiplyMatrixD(*p, *q);

            blocksize = 1;

            for (size_t j = 0; j < d->cols - 1; ++j) {
                converged = converged && (fabs(d->matrix[d->rows - 1][j]) < 1e-100);
            }

            if (converged) {
                //printf("1: Exited early at iteration %lu for stage %lu\n", i, d->rows);
                break;
            }
        }
    }

    if (m.rows > 1 + blocksize) {
        MatrixD mi = _submatrixMatrixD(*d, 0, 0, m.rows - blocksize, m.cols - blocksize);
        MatrixD *di = malloc(sizeof(MatrixD));
        MatrixD *pi = malloc(sizeof(MatrixD));
        qrAlgorithmMatrixD(mi, iterations, di, pi);

        MatrixD piFull = newMatrixD(m.rows, m.cols);
        for (size_t i = 0; i < m.rows - blocksize; ++i) {
            for (size_t j = 0; j < m.cols - blocksize; ++j) {
                piFull.matrix[i][j] = pi->matrix[i][j];
            }
        }
        for (size_t i = m.rows - blocksize; i < m.rows; ++i) {
            piFull.matrix[i][i] = 1;
        }
        *d = multiplyMatrixD(transposeMatrixD(piFull), multiplyMatrixD(*d, piFull));
        *p = multiplyMatrixD(*p, piFull);

        free(di);
        free(pi);
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

// Strassen Algorithm for square matrices (for now)
MatrixD strassenMultiplyMatrixD(MatrixD a, MatrixD b) {
    if (a.rows != a.cols || b.rows != b.cols) {
        fprintf(stderr, "Error: Cannot perform Strassen algorithm on non-square matrices of shape (%lu, %lu) and (%lu, %lu)\n", a.rows, a.cols, b.rows, b.cols);
        exit(1);
    }

    if (a.cols != b.rows) {
        fprintf(stderr, "Error: Cannot perform Strassen algorithm on matrices of shape (%lu, %lu) and (%lu, %lu)\n", a.rows, a.cols, b.rows, b.cols);
        exit(1);
    }

    _BlockMatrixD A11 = _toBlockMatrixD(a, (a.rows + 1) / 2, (a.cols + 1) / 2, 0, 0);
    _BlockMatrixD B11 = _toBlockMatrixD(b, (b.rows + 1) / 2, (b.cols + 1) / 2, 0, 0);

    _BlockMatrixD A12 = _toBlockMatrixD(a, (a.rows + 1) / 2, a.cols / 2, 0, (a.cols + 1) / 2);
    _BlockMatrixD B12 = _toBlockMatrixD(b, (b.rows + 1) / 2, b.cols / 2, 0, (b.cols + 1) / 2);

    _BlockMatrixD A21 = _toBlockMatrixD(a, a.rows / 2, (a.cols + 1) / 2, (a.rows + 1) / 2, 0);
    _BlockMatrixD B21 = _toBlockMatrixD(b, b.rows / 2, (b.cols + 1) / 2, (b.rows + 1) / 2, 0);

    _BlockMatrixD A22 = _toBlockMatrixD(a, a.rows / 2, a.cols / 2, (a.rows + 1) / 2, (a.cols + 1) / 2);
    _BlockMatrixD B22 = _toBlockMatrixD(b, b.rows / 2, b.cols / 2, (b.rows + 1) / 2, (b.cols + 1) / 2);

    printf("Initialized blocks\n");

    _BlockMatrixD M1 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A11, A22), _addBlockMatrixD(B11, B22));
    _BlockMatrixD M2 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A21, A22), B11);
    _BlockMatrixD M3 = _strassenMultiplyBlockMatrixD(A11, _subtractBlockMatrixD(B12, B22));
    _BlockMatrixD M4 = _strassenMultiplyBlockMatrixD(A22, _subtractBlockMatrixD(B21, B11));
    _BlockMatrixD M5 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A11, A12), B22);
    _BlockMatrixD M6 = _strassenMultiplyBlockMatrixD(_subtractBlockMatrixD(A21, A11), _addBlockMatrixD(B11, B12));
    _BlockMatrixD M7 = _strassenMultiplyBlockMatrixD(_subtractBlockMatrixD(A12, A22), _addBlockMatrixD(B21, B22));

    printf("Computed Ms\n");

    // C11 = M1 + M4 - M5 + M7
    _BlockMatrixD C11 = _addBlockMatrixD(M1, M4);
    _updateSubtractBlockMatrixD(&C11, M5);
    _updateAddBlockMatrixD(&C11, M7);

    // C12 = M3 + M5
    _BlockMatrixD C12 = _addBlockMatrixD(M3, M5);

    // C21 = M2 + M4
    _BlockMatrixD C21 = _addBlockMatrixD(M2, M4);

    // C22 = M1 - M2 + M3 + M6
    _BlockMatrixD C22 = _subtractBlockMatrixD(M1, M2);
    _updateAddBlockMatrixD(&C22, M3);
    _updateAddBlockMatrixD(&C22, M6);

    printf("Computed Cs\n");

    MatrixD C = newMatrixD(a.rows, a.cols);

    for (size_t i = 0; i < C.rows; ++i) {
        for (size_t j = 0; j < C.cols; ++j) {
            if (i < C11.rows) {
                if (j < C11.cols) {
                    C.matrix[i][j] = C11.matrix[i + C11.rowOffset][j + C11.colOffset];
                } else {
                    C.matrix[i][j] = C12.matrix[i + C12.rowOffset][j - C11.cols + C12.colOffset];
                }
            } else {
                if (j < C21.cols) {
                    C.matrix[i][j] = C21.matrix[i - C11.rows + C21.rowOffset][j + C21.colOffset];
                } else {
                    C.matrix[i][j] = C22.matrix[i - C11.rows + C22.rowOffset][j - C21.cols + C22.colOffset];
                }
            }
        }
    }

    printf("Computed C\n");

    return C;
}

_BlockMatrixD _toBlockMatrixD(MatrixD m, size_t rows, size_t cols, size_t rowOffset, size_t colOffset) {
    if (rowOffset + rows > m.rows || colOffset + cols > m.cols) {
        fprintf(stderr, "Error: Resulting _BlockMatrixD would overrun bounds of the array\n");
        exit(1);
    }
    _BlockMatrixD mBlock = { rows, cols, rowOffset, colOffset, m.matrix };
    return mBlock;
}

MatrixD _fromBlockMatrixD(_BlockMatrixD mBlock) {
    double **matrix = malloc(sizeof(double*) * mBlock.rows);
    for (size_t i = 0; i < mBlock.rows; ++i) {
        matrix[i] = mBlock.matrix[i + mBlock.rowOffset] + mBlock.colOffset;
    }

    MatrixD m = { mBlock.rows, mBlock.cols, matrix };
    return m;
}

_BlockMatrixD _addBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    printf("_addBlockMatrixD\n");
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum

    printf("a\n");
    printf("%lu, %lu, %lu, %lu\n", a.rows, a.cols, a.rowOffset, a.colOffset);
    displayMatrixD(_fromBlockMatrixD(a));
    printf("b\n");
    displayMatrixD(_fromBlockMatrixD(b));

    double **matrix = malloc(sizeof(double*) * rows);
    for (size_t i = 0; i < rows; ++i) {
        matrix[i] = malloc(sizeof(double) * cols);
        for (size_t j = 0; j < cols; ++j) {
            if (i >= a.rows || j >= a.cols) {
                matrix[i][j] = b.matrix[i + b.rowOffset][j + b.colOffset];
            } else if (i >= b.rows || j >= b.cols) {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset];
            } else {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset] + b.matrix[i + b.rowOffset][j + b.colOffset];
            }
        }
    }

    _BlockMatrixD c = { rows, cols, 0, 0, matrix };

    printf("Exit _addBlockMatrixD\n");

    return c;
}

_BlockMatrixD _subtractBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    printf("_subtractBlockMatrixD\n");
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum

    double **matrix = malloc(sizeof(double*) * rows);
    for (size_t i = 0; i < rows; ++i) {
        matrix[i] = malloc(sizeof(double) * cols);
        for (size_t j = 0; j < cols; ++j) {
            if (i >= a.rows || j >= a.cols) {
                matrix[i][j] = -b.matrix[i + b.rowOffset][j + b.colOffset];
            } else if (i >= b.rows || j >= b.cols) {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset];
            } else {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset] - b.matrix[i + b.rowOffset][j + b.colOffset];
            }
        }
    }

    _BlockMatrixD c = { rows, cols, 0, 0, matrix };
    return c;
}

// a = a + b
void _updateAddBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b) {
    if (a->rows < b.rows || a->cols < b.cols) {
        fprintf(stderr, "Error: Cannot update add matrix of shape (%lu, %lu) with (%lu, %lu)\n", a->rows, a->cols, b.rows, b.cols);
        exit(1);
    }

    for (size_t i = 0; i < b.rows; ++i) {
        for (size_t j = 0; j < b.cols; ++j) {
            a->matrix[i + a->rowOffset][j + a->colOffset] += b.matrix[i + b.rowOffset][j + b.colOffset];
        }
    }
}

// a = a - b
void _updateSubtractBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b) {
    if (a->rows < b.rows || a->cols < b.cols) {
        fprintf(stderr, "Error: Cannot update add matrix of shape (%lu, %lu) with (%lu, %lu)\n", a->rows, a->cols, b.rows, b.cols);
        exit(1);
    }

    for (size_t i = 0; i < b.rows; ++i) {
        for (size_t j = 0; j < b.cols; ++j) {
            a->matrix[i + a->rowOffset][j + a->colOffset] -= b.matrix[i + b.rowOffset][j + b.colOffset];
        }
    }
}

_BlockMatrixD _strassenMultiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    printf("_strassenMultiplyBlockMatrixD\n");
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum
    size_t dim = (rows > cols) ? rows : cols;

    if (rows < 4 || cols < 4) {
        return _multiplyBlockMatrixD(a, b);
    }

    printf("%lu, %lu, %lu, %lu\n", a.rows, a.cols, a.rowOffset, a.colOffset);

    _BlockMatrixD A11 = { (dim + 1) / 2, (dim + 1) / 2, a.rowOffset, a.colOffset, a.matrix };
    _BlockMatrixD B11 = { (dim + 1) / 2, (dim + 1) / 2, b.rowOffset, b.colOffset, b.matrix };

    _BlockMatrixD A12 = { (dim + 1) / 2, a.cols - (dim + 1) / 2, a.rowOffset, a.colOffset + (dim + 1) / 2, a.matrix };
    _BlockMatrixD B12 = { (dim + 1) / 2, b.cols - (dim + 1) / 2, b.rowOffset, b.colOffset + (dim + 1) / 2, b.matrix };

    _BlockMatrixD A21 = { a.rows - (dim + 1) / 2, (dim + 1) / 2, a.rowOffset + (dim + 1) / 2, a.colOffset, a.matrix };
    _BlockMatrixD B21 = { b.rows - (dim + 1) / 2, (dim + 1) / 2, b.rowOffset + (dim + 1) / 2, b.colOffset, b.matrix };

    _BlockMatrixD A22 = { a.rows - (dim + 1) / 2, a.cols - (dim + 1) / 2, a.rowOffset + (dim + 1) / 2, a.colOffset + (dim + 1) / 2, a.matrix };
    _BlockMatrixD B22 = { b.rows - (dim + 1) / 2, b.cols - (dim + 1) / 2, b.rowOffset + (dim + 1) / 2, b.colOffset + (dim + 1) / 2, b.matrix };

    _BlockMatrixD M1 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A11, A22), _addBlockMatrixD(B11, B22));
    _BlockMatrixD M2 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A21, A22), B11);
    _BlockMatrixD M3 = _strassenMultiplyBlockMatrixD(A11, _subtractBlockMatrixD(B12, B22));
    _BlockMatrixD M4 = _strassenMultiplyBlockMatrixD(A22, _subtractBlockMatrixD(B21, B11));
    _BlockMatrixD M5 = _strassenMultiplyBlockMatrixD(_addBlockMatrixD(A11, A12), B22);
    _BlockMatrixD M6 = _strassenMultiplyBlockMatrixD(_subtractBlockMatrixD(A21, A11), _addBlockMatrixD(B11, B12));
    _BlockMatrixD M7 = _strassenMultiplyBlockMatrixD(_subtractBlockMatrixD(A12, A22), _addBlockMatrixD(B21, B22));

    // C11 = M1 + M4 - M5 + M7
    _BlockMatrixD C11 = _addBlockMatrixD(M1, M4);
    _updateSubtractBlockMatrixD(&C11, M5);
    _updateAddBlockMatrixD(&C11, M7);

    // C12 = M3 + M5
    _BlockMatrixD C12 = _addBlockMatrixD(M3, M5);

    // C21 = M2 + M4
    _BlockMatrixD C21 = _addBlockMatrixD(M2, M4);

    // C22 = M1 - M2 + M3 + M6
    _BlockMatrixD C22 = _subtractBlockMatrixD(M1, M2);
    _updateAddBlockMatrixD(&C22, M3);
    _updateAddBlockMatrixD(&C22, M6);

    double **matrix = malloc(sizeof(double*) * dim);
    for (size_t i = 0; i < dim; ++i) {
        matrix[i] = malloc(sizeof(double) * dim);
        for (size_t j = 0; j < dim; ++j) {
            if (i < C11.rows) {
                if (j < C11.cols) {
                    matrix[i][j] = C11.matrix[i + C11.rowOffset][j + C11.colOffset];
                } else {
                    matrix[i][j] = C12.matrix[i + C12.rowOffset][j - C11.cols + C12.colOffset];
                }
            } else {
                if (j < C21.cols) {
                    matrix[i][j] = C21.matrix[i - C11.rows + C21.rowOffset][j + C21.colOffset];
                } else {
                    matrix[i][j] = C22.matrix[i - C11.rows + C22.rowOffset][j - C21.cols + C22.colOffset];
                }
            }
        }
    }

    _BlockMatrixD C = { dim, dim, 0, 0, matrix };

    _BlockMatrixD Calt = _multiplyBlockMatrixD(a, b);
    double max = 0;
    for (size_t i = 0; i < C.rows; ++i) {
        for (size_t j = 0; j < C.cols; ++j) {
            max = (fabs(C.matrix[i][j] - Calt.matrix[i][j]) > max) ? fabs(C.matrix[i][j] - Calt.matrix[i][j]) : max;
        }
    }
    printf("(%lu, %lu) difference: %e\n", C.rows, C.cols, max);

    if (max > 1) {
        displayMatrixD(_fromBlockMatrixD(a));
        displayMatrixD(_fromBlockMatrixD(A11));
        displayMatrixD(_fromBlockMatrixD(A12));
        displayMatrixD(_fromBlockMatrixD(A21));
        displayMatrixD(_fromBlockMatrixD(A22));
        displayMatrixD(_fromBlockMatrixD(C));
        displayMatrixD(_fromBlockMatrixD(Calt));
        exit(1);
    }

    return C;
}

_BlockMatrixD _multiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    printf("_multiplyBlockMatrixD\n");
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum
    size_t dim = (rows > cols) ? rows : cols;

    double **matrix = malloc(sizeof(double*) * dim);
    for (size_t i = 0; i < dim; ++i) {
        matrix[i] = malloc(sizeof(double) * dim);
        for (size_t j = 0; j < dim; ++j) {
            if (i >= a.rows) {
                matrix[i][j] = 0;
            } else if (j >= b.cols) {
                matrix[i][j] = 0;
            } else {
                double sum = 0;
                for (size_t k = 0; k < dim; ++k) {
                    if (k < a.cols && k < b.rows) {
                        sum += a.matrix[i + a.rowOffset][k + a.colOffset] * b.matrix[k + b.rowOffset][j + b.colOffset];
                    }
                }
                matrix[i][j] = sum;
            }
        }
    }

    _BlockMatrixD c = { dim, dim, 0, 0, matrix };
    return c;
}