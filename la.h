#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#define ONE_NORM 0
#define TWO_NORM 1
#define INF_NORM 2

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
bool equalsMatrixD(MatrixD a, MatrixD b, double tolerance);

void addMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
void subtractMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
void multiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
void addMultiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
void scaleMatrixD(MatrixD *y, double k, MatrixD *a);
void innerProductMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
void outerProductMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);

double slowDeterminantMatrixD(MatrixD a);

MatrixD appendMatrixD(MatrixD a, MatrixD b);

void echelonMatrixD(MatrixD *m);
void _swapRowMatrixD(MatrixD *m, size_t row1, size_t row2);
void _rowAdditionMatrixD(MatrixD *m, double k, size_t row1, size_t row2);
double fastDeterminantMatrixD(MatrixD m);
void reducedEchelonMatrixD(MatrixD *m);
void _scaleRowMatrixD(MatrixD *m, double k, size_t row);

MatrixD inverseMatrixD(MatrixD m);
MatrixD _submatrixMatrixD(MatrixD m, size_t row, size_t col, size_t height, size_t width);
MatrixD transposeMatrixD(MatrixD m);
double normMatrixD(MatrixD m, int norm_type);
void qrDecompositionMatrixD(MatrixD m, MatrixD *q, MatrixD *r);
void qrAlgorithmMatrixD(MatrixD m, size_t iterations, MatrixD *d, MatrixD *p);
int solveLinearMatrixD(MatrixD A, MatrixD b, MatrixD *x, MatrixD *N);

void strassenMultiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b);
_BlockMatrixD _toBlockMatrixD(MatrixD m, size_t rows, size_t cols, size_t rowOffset, size_t colOffset);
MatrixD _fromBlockMatrixD(_BlockMatrixD mBlock);
void _freeBlockMatrixD(_BlockMatrixD *m);
_BlockMatrixD _addBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
_BlockMatrixD _subtractBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
void _updateAddBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b);
void _updateSubtractBlockMatrixD(_BlockMatrixD *a, _BlockMatrixD b);
_BlockMatrixD _strassenMultiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);
_BlockMatrixD _multiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b);

void eigensMatrixD(MatrixD *a, MatrixD *lambda, MatrixD *v);
void _twoByTwoEigenvaluesMatrixD(MatrixD *m, MatrixD *lambda);
double _eigenvalueIterationMatrixD(MatrixD *a, MatrixD *vi, double lambda);
void _complexEigenvectorsMatrixD(MatrixD *a, MatrixD *v, double real, double imag);


MatrixD newMatrixD(size_t rows, size_t cols) {
    double** mat = malloc(sizeof(double*) * rows);
    for (size_t i = 0; i < rows; ++i) {
        mat[i] = calloc(cols, sizeof(double));
    }
    MatrixD m = {rows, cols, mat};
    return m;
}

void freeMatrixD(MatrixD *a) {
    for (size_t i = 0; i < a->rows; ++i) {
        free(a->matrix[i]);
    }
    free(a->matrix);
    a->matrix = NULL;
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

bool equalsMatrixD(MatrixD a, MatrixD b, double tolerance) {
    if (a.rows != b.rows || a.cols != b.cols) return false;

    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < a.cols; ++j) {
            if (fabs(a.matrix[i][j] - b.matrix[i][j]) > tolerance) return false;
        }
    }

    return true;
}

// c = a + b
void addMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    // Check rows and columns match
    if (a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "Error: Could not add matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->rows != c->rows || a->cols != c->cols) {
        fprintf(stderr, "Error: Could not store sum of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, a->cols, c->rows, c->cols);
        exit(1);
    }

    for (size_t i = 0; i < c->rows; ++i) {
        for (size_t j = 0; j < c->cols; ++j) {
            c->matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];
        }
    }
}

// c = a - b
void subtractMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    // Check rows and columns match
    if (a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "Error: Could not subtract matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->rows != c->rows || a->cols != c->cols) {
        fprintf(stderr, "Error: Could not store difference of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, a->cols, c->rows, c->cols);
        exit(1);
    }

    for (size_t i = 0; i < c->rows; ++i) {
        for (size_t j = 0; j < c->cols; ++j) {
            c->matrix[i][j] = a->matrix[i][j] - b->matrix[i][j];
        }
    }
}

// c = a * b
void multiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Could not multiply matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->rows != c->rows || b->cols != c->cols) {
        fprintf(stderr, "Error: Could not store product of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, b->cols, c->rows, c->cols);
        exit(1);
    }

    if (c == a || c == b) {
        MatrixD m = newMatrixD(c->rows, c->cols);
        for (size_t i = 0; i < a->rows; ++i) {
            for (size_t j = 0; j < b->cols; ++j) {
                for (size_t k = 0; k < a->cols; ++k) {
                    m.matrix[i][j] += a->matrix[i][k] * b->matrix[k][j];
                }
            }
        }
        for (size_t i = 0; i < a->rows; ++i) {
            for (size_t j = 0; j < b->cols; ++j) {
                c->matrix[i][j] = m.matrix[i][j];
            }
        }
        freeMatrixD(&m);
    } else {
        for (size_t i = 0; i < a->rows; ++i) {
            for (size_t j = 0; j < b->cols; ++j) {
                c->matrix[i][j] = 0;
                for (size_t k = 0; k < a->cols; ++k) {
                    c->matrix[i][j] += a->matrix[i][k] * b->matrix[k][j];
                }
            }
        }
    }
}

// c = c + a * b
void addMultiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Could not multiply matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->rows != c->rows || b->cols != c->cols) {
        fprintf(stderr, "Error: Could not store product of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, b->cols, c->rows, c->cols);
        exit(1);
    }

    for (size_t i = 0; i < a->rows; ++i) {
        for (size_t j = 0; j < b->cols; ++j) {
            for (size_t k = 0; k < a->cols; ++k) {
                c->matrix[i][j] += a->matrix[i][k] * b->matrix[k][j];
            }
        }
    }
}

// y = k * a
void scaleMatrixD(MatrixD *y, double k, MatrixD *a) {
    if (y->rows != a->rows || y->cols != a->cols) {
        fprintf(stderr, "Error: Cannot store scalar multiple of shape (%lu, %lu) in matrix of shape (%lu, %lu)\n", a->rows, a->cols, y->rows, y->cols);
        exit(1);
    }

    for (size_t i = 0; i < a->rows; ++i) {
        for (size_t j = 0; j < a->cols; ++j) {
            y->matrix[i][j] = k * a->matrix[i][j];
        }
    }
}

// c = a^T * b
void innerProductMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    if (a->rows != b->rows) {
        fprintf(stderr, "Error: Could not compute inner product for matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->cols != c->rows || b->cols != c->cols) {
        fprintf(stderr, "Error: Could not store product of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->cols, b->cols, c->rows, c->cols);
        exit(1);
    }

    for (size_t i = 0; i < c->rows; ++i) {
        for (size_t j = 0; j < c->cols; ++j) {
            c->matrix[i][j] = 0;
            for (size_t k = 0; k < a->rows; ++k) {
                c->matrix[i][j] += a->matrix[k][i] * b->matrix[k][j];
            }
        }
    }
}

// c = a * b^T
void outerProductMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    if (a->cols != b->cols) {
        fprintf(stderr, "Error: Could not multiply matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->rows != c->rows || b->rows != c->cols) {
        fprintf(stderr, "Error: Could not store product of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, b->rows, c->rows, c->cols);
        exit(1);
    }

    for (size_t i = 0; i < c->rows; ++i) {
        for (size_t j = 0; j < c->cols; ++j) {
            c->matrix[i][j] = 0;
            for (size_t k = 0; k < a->cols; ++k) {
                c->matrix[i][j] += a->matrix[i][k] * b->matrix[j][k];
            }
        }
    }
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
    freeMatrixD(&m);

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

// Computes the echelon matrix of m in-place
void echelonMatrixD(MatrixD *m) {
    size_t pivots = 0;
    for (size_t col = 0; col < m->cols; ++col) {
        size_t row = pivots;
        while (m->matrix[row][col] == 0) {
            row++;
            if (row >= m->rows) {
                row = m->rows - 1;
                break;
            }
        }
        if (m->matrix[row][col] != 0) {
            _swapRowMatrixD(m, pivots, row);
            for (size_t i = pivots + 1; i < m->rows; ++i) {
                _rowAdditionMatrixD(m, -m->matrix[i][col] / m->matrix[pivots][col], pivots, i);
            }
            pivots++;
        }
        
        if (pivots == m->rows) break;
    }
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

    freeMatrixD(&a);

    return det;
}

void reducedEchelonMatrixD(MatrixD *m) {
    size_t pivots = 0;
    for (size_t col = 0; col < m->cols; ++col) {
        size_t row = pivots;
        for (size_t i = pivots + 1; i < m->rows; ++i) {
            if (fabs(m->matrix[i][col]) > fabs(m->matrix[row][col])) row = i;
        }
        if (fabs(m->matrix[row][col]) > 1e-10) {
            _swapRowMatrixD(m, pivots, row);
            _scaleRowMatrixD(m, 1 / m->matrix[pivots][col], pivots);
            for (size_t i = 0; i < m->rows; ++i) {
                if (i == pivots) continue;
                _rowAdditionMatrixD(m, -m->matrix[i][col], pivots, i);
                m->matrix[i][col] = 0; // manually set to avoid error
            }
            pivots++;
        } else {
            for (size_t i = pivots; i < m->rows; ++i) {
                m->matrix[i][col] = 0; // manually set to avoid error
            }
        }
        
        if (pivots == m->rows) break;
    }
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

    MatrixD I = identityMatrixD(m.rows);
    MatrixD a = appendMatrixD(m, I);
    reducedEchelonMatrixD(&a);
    MatrixD inv = _submatrixMatrixD(a, 0, m.rows, m.rows, m.rows);

    freeMatrixD(&I);
    freeMatrixD(&a);

    return inv;
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

// Computes the vector 1-norm, 2-norm, or inf-norm of m
double normMatrixD(MatrixD m, int norm_type) {
    double sum = 0;
    if (norm_type == ONE_NORM) {
        for (size_t i = 0; i < m.rows; ++i) {
            for (size_t j = 0; j < m.cols; ++j) {
                sum += fabs(m.matrix[i][j]);
            }
        }
    } else if (norm_type == TWO_NORM) {
        for (size_t i = 0; i < m.rows; ++i) {
            for (size_t j = 0; j < m.cols; ++j) {
                sum += m.matrix[i][j] * m.matrix[i][j];
            }
        }
        sum = sqrt(sum);
    } else if (norm_type == INF_NORM) {
        for (size_t i = 0; i < m.rows; ++i) {
            for (size_t j = 0; j < m.cols; ++j) {
                if (fabs(m.matrix[i][j]) > sum) sum = fabs(m.matrix[i][j]);
            }
        }
    }

    return sum;
}

void qrDecompositionMatrixD(MatrixD m, MatrixD *q, MatrixD *r) {
    size_t t = (m.rows > m.cols) ? m.cols : m.rows; // minimum

    MatrixD mprime = copyMatrixD(m);
    *q = identityMatrixD(m.rows);
    *r = copyMatrixD(m);

    MatrixD qprime, x;
    MatrixD Qi = newMatrixD(m.rows, m.rows);

    for (size_t i = 0; i < t; ++i) {
        x = _submatrixMatrixD(mprime, 0, 0, mprime.rows, 1);

        double alpha = normMatrixD(x, TWO_NORM);
        if (fabs(alpha) > 1e-100) {
            qprime = newMatrixD(mprime.rows, mprime.rows);
            if (x.matrix[0][0] >= 0) alpha *= -1;

            MatrixD e1 = newMatrixD(mprime.rows, 1);
            e1.matrix[0][0] = 1;
            scaleMatrixD(&e1, -alpha, &e1);

            MatrixD u = newMatrixD(mprime.rows, 1);
            addMatrixD(&u, &x, &e1);
            scaleMatrixD(&u, 1 / normMatrixD(u, TWO_NORM), &u);

            outerProductMatrixD(&qprime, &u, &u);
            scaleMatrixD(&qprime, -2, &qprime);
            MatrixD I = identityMatrixD(mprime.rows);
            addMatrixD(&qprime, &I, &qprime);

            freeMatrixD(&e1);
            freeMatrixD(&u);
            freeMatrixD(&I);
        } else {
            qprime = identityMatrixD(mprime.rows);
        }

        for (size_t j = 0; j < m.rows; ++j) {
            for (size_t k = 0; k < m.rows; ++k) {
                if (j < i) {
                    if (j == k) {
                        Qi.matrix[j][k] = 1;
                    } else {
                        Qi.matrix[j][k] = 0;
                    }
                } else if (k < i) {
                    Qi.matrix[j][k] = 0;
                } else {
                    Qi.matrix[j][k] = qprime.matrix[j-i][k-i];
                }
            }
        }

        multiplyMatrixD(q, q, &Qi);
        multiplyMatrixD(r, &Qi, r);

        freeMatrixD(&mprime);
        freeMatrixD(&qprime);
        freeMatrixD(&x);

        mprime = _submatrixMatrixD(*r, i + 1, i + 1, m.rows - i - 1, m.rows - i - 1);
    }

    freeMatrixD(&mprime);
    freeMatrixD(&Qi);
}

void qrAlgorithmMatrixD(MatrixD m, size_t iterations, MatrixD *d, MatrixD *p) {
    if (m.rows != m.cols) {
        fprintf(stderr, "Error: Cannot perform QR Algorithm on non-square matrix of shape (%lu, %lu)\n", m.rows, m.cols);
        exit(1);
    }

    MatrixD *q = malloc(sizeof(MatrixD));
    MatrixD *r = malloc(sizeof(MatrixD));

    *d = copyMatrixD(m);
    *p = identityMatrixD(m.rows);

    MatrixD M, I;

    size_t blocksize = 1;
    for (size_t i = 0; i < iterations; ++i) {

        bool converged = true;

        M = newMatrixD(d->rows, d->cols);
        I = identityMatrixD(d->rows);

        // check if complex eigenvalues in lower right 2x2 block
        double trace = d->matrix[d->rows - 2][d->cols - 2] + d->matrix[d->rows - 1][d->cols - 1];
        double det = d->matrix[d->rows - 2][d->cols - 2] * d->matrix[d->rows - 1][d->cols - 1] - d->matrix[d->rows - 2][d->cols - 1] * d->matrix[d->rows - 1][d->cols - 2];
        if (trace * trace < 4 * det) {
            // complex eigenvalues
            // double implicit shift strategy

            // Compute M = d^2 - trace * d + det * I
            scaleMatrixD(&I, det, &I);
            scaleMatrixD(&M, -trace, d);
            addMultiplyMatrixD(&M, d, d);
            addMatrixD(&M, &M, &I);

            qrDecompositionMatrixD(M, q, r);
            MatrixD qt = transposeMatrixD(*q);
            multiplyMatrixD(d, d, q);
            multiplyMatrixD(d, &qt, d); // TODO: Replace with Ut A V function
            multiplyMatrixD(p, p, q);

            freeMatrixD(&I);
            freeMatrixD(&M);
            freeMatrixD(&qt);

            freeMatrixD(q);
            freeMatrixD(r);

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
            
            // Compute M = d - lambda * I
            scaleMatrixD(&I, lambda, &I);
            subtractMatrixD(&M, d, &I);

            qrDecompositionMatrixD(M, q, r);

            // d = r * q + lambda * I
            multiplyMatrixD(d, r, q);
            addMatrixD(d, d, &I);
            multiplyMatrixD(p, p, q);

            freeMatrixD(&I);
            freeMatrixD(&M);

            freeMatrixD(q);
            freeMatrixD(r);

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
        MatrixD piFullT = transposeMatrixD(piFull);
        multiplyMatrixD(d, d, &piFull);
        multiplyMatrixD(d, &piFullT, d); // TODO: Replace with Ut A V function
        multiplyMatrixD(p, p, &piFull);

        freeMatrixD(&mi);
        freeMatrixD(&piFull);
        freeMatrixD(&piFullT);
        freeMatrixD(di);
        freeMatrixD(pi);
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
    reducedEchelonMatrixD(&augmented);

    MatrixD sol = newMatrixD(A.cols, 1);
    MatrixD null;
    size_t pivots = augmented.rows;
    bool foundPivot;
    size_t lastNullColumn;
    size_t lastPivotColumn = augmented.cols - 1;
    size_t *freeColumns = malloc(sizeof(size_t) * (A.cols));
    for (int i = (int) augmented.rows - 1; i >= 0; --i) {
        // look for pivot
        foundPivot = false;
        for (size_t j = (size_t) i; j < augmented.cols; ++j) {
            if (j == augmented.cols - 1) {
                if ((fabs(augmented.matrix[i][j]) >= DBL_EPSILON) && !foundPivot) {
                    // inconsistent
                    freeMatrixD(&sol);
                    freeMatrixD(&augmented);
                    free(freeColumns);
                    return -1;
                } else if (!foundPivot) {
                    // row of all zeros
                    pivots--;
                    break;
                } else {
                    break;
                }
            }
            if (fabs(augmented.matrix[i][j]) >= DBL_EPSILON) {
                if (!foundPivot) {
                    // found pivot
                    foundPivot = true;
                    if (i == (int) pivots - 1) {
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

    freeMatrixD(&augmented);
    free(freeColumns);

    return (int) null.cols;
}

// Strassen Algorithm for square matrices (for now)
void strassenMultiplyMatrixD(MatrixD *c, MatrixD *a, MatrixD *b) {
    if (a->rows != a->cols || b->rows != b->cols) {
        fprintf(stderr, "Error: Cannot perform Strassen algorithm on non-square matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Cannot perform Strassen algorithm on matrices of shape (%lu, %lu) and (%lu, %lu)\n", a->rows, a->cols, b->rows, b->cols);
        exit(1);
    }

    if (c->rows != a->rows || c->cols != b->cols) {
        fprintf(stderr, "Error: Cannot store matrix of shape (%lu, %lu) into matrix of shape (%lu, %lu)\n", a->rows, b->cols, c->rows, c->cols);
        exit(1);
    }

    _BlockMatrixD A11 = _toBlockMatrixD(*a, (a->rows + 1) / 2, (a->cols + 1) / 2, 0, 0);
    _BlockMatrixD B11 = _toBlockMatrixD(*b, (b->rows + 1) / 2, (b->cols + 1) / 2, 0, 0);

    _BlockMatrixD A12 = _toBlockMatrixD(*a, (a->rows + 1) / 2, a->cols / 2, 0, (a->cols + 1) / 2);
    _BlockMatrixD B12 = _toBlockMatrixD(*b, (b->rows + 1) / 2, b->cols / 2, 0, (b->cols + 1) / 2);

    _BlockMatrixD A21 = _toBlockMatrixD(*a, a->rows / 2, (a->cols + 1) / 2, (a->rows + 1) / 2, 0);
    _BlockMatrixD B21 = _toBlockMatrixD(*b, b->rows / 2, (b->cols + 1) / 2, (b->rows + 1) / 2, 0);

    _BlockMatrixD A22 = _toBlockMatrixD(*a, a->rows / 2, a->cols / 2, (a->rows + 1) / 2, (a->cols + 1) / 2);
    _BlockMatrixD B22 = _toBlockMatrixD(*b, b->rows / 2, b->cols / 2, (b->rows + 1) / 2, (b->cols + 1) / 2);

    _BlockMatrixD M1a = _addBlockMatrixD(A11, A22);
    _BlockMatrixD M1b = _addBlockMatrixD(B11, B22);
    _BlockMatrixD M1 = _strassenMultiplyBlockMatrixD(M1a, M1b);
    _freeBlockMatrixD(&M1a);
    _freeBlockMatrixD(&M1b);

    _BlockMatrixD M2a = _addBlockMatrixD(A21, A22);
    _BlockMatrixD M2 = _strassenMultiplyBlockMatrixD(M2a, B11);
    _freeBlockMatrixD(&M2a);

    _BlockMatrixD M3a = _subtractBlockMatrixD(B12, B22);
    _BlockMatrixD M3 = _strassenMultiplyBlockMatrixD(A11, M3a);
    _freeBlockMatrixD(&M3a);

    _BlockMatrixD M4a = _subtractBlockMatrixD(B21, B11);
    _BlockMatrixD M4 = _strassenMultiplyBlockMatrixD(A22, M4a);
    _freeBlockMatrixD(&M4a);

    _BlockMatrixD M5a = _addBlockMatrixD(A11, A12);
    _BlockMatrixD M5 = _strassenMultiplyBlockMatrixD(M5a, B22);
    _freeBlockMatrixD(&M5a);

    _BlockMatrixD M6a = _subtractBlockMatrixD(A21, A11);
    _BlockMatrixD M6b = _addBlockMatrixD(B11, B12);
    _BlockMatrixD M6 = _strassenMultiplyBlockMatrixD(M6a, M6b);
    _freeBlockMatrixD(&M6a);
    _freeBlockMatrixD(&M6b);

    _BlockMatrixD M7a = _subtractBlockMatrixD(A12, A22);
    _BlockMatrixD M7b = _addBlockMatrixD(B21, B22);
    _BlockMatrixD M7 = _strassenMultiplyBlockMatrixD(M7a, M7b);
    _freeBlockMatrixD(&M7a);
    _freeBlockMatrixD(&M7b);

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

    for (size_t i = 0; i < c->rows; ++i) {
        for (size_t j = 0; j < c->cols; ++j) {
            if (i < C11.rows) {
                if (j < C11.cols) {
                    c->matrix[i][j] = C11.matrix[i + C11.rowOffset][j + C11.colOffset];
                } else {
                    c->matrix[i][j] = C12.matrix[i + C12.rowOffset][j - C11.cols + C12.colOffset];
                }
            } else {
                if (j < C21.cols) {
                    c->matrix[i][j] = C21.matrix[i - C11.rows + C21.rowOffset][j + C21.colOffset];
                } else {
                    c->matrix[i][j] = C22.matrix[i - C11.rows + C22.rowOffset][j - C21.cols + C22.colOffset];
                }
            }
        }
    }

    _freeBlockMatrixD(&M1);
    _freeBlockMatrixD(&M2);
    _freeBlockMatrixD(&M3);
    _freeBlockMatrixD(&M4);
    _freeBlockMatrixD(&M5);
    _freeBlockMatrixD(&M6);
    _freeBlockMatrixD(&M7);

    _freeBlockMatrixD(&C11);
    _freeBlockMatrixD(&C12);
    _freeBlockMatrixD(&C21);
    _freeBlockMatrixD(&C22);
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

void _freeBlockMatrixD(_BlockMatrixD *m) {
    if (m->rowOffset != 0 || m->colOffset != 0) {
        fprintf(stderr, "Error: Cannot free Block Matrix that is a child of another Matrix / Block Matrix\n");
        exit(1);
    }

    for (size_t i = 0; i < m->rows; ++i) {
        free(m->matrix[i]);
    }
    free(m->matrix);
    m->matrix = NULL;
}

_BlockMatrixD _addBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum
    size_t dim = (rows > cols) ? rows : cols;

    double **matrix = malloc(sizeof(double*) * dim);
    for (size_t i = 0; i < dim; ++i) {
        matrix[i] = malloc(sizeof(double) * dim);
        for (size_t j = 0; j < dim; ++j) {
            if ((i >= a.rows || j >= a.cols) && (i >= b.rows || j >= b.cols)) {
                matrix[i][j] = 0;
            } else if (i >= a.rows || j >= a.cols) {
                matrix[i][j] = b.matrix[i + b.rowOffset][j + b.colOffset];
            } else if (i >= b.rows || j >= b.cols) {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset];
            } else {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset] + b.matrix[i + b.rowOffset][j + b.colOffset];
            }
        }
    }

    _BlockMatrixD c = { dim, dim, 0, 0, matrix };
    return c;
}

_BlockMatrixD _subtractBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum
    size_t dim = (rows > cols) ? rows : cols;

    double **matrix = malloc(sizeof(double*) * dim);
    for (size_t i = 0; i < dim; ++i) {
        matrix[i] = malloc(sizeof(double) * dim);
        for (size_t j = 0; j < dim; ++j) {
            if ((i >= a.rows || j >= a.cols) && (i >= b.rows || j >= b.cols)) {
                matrix[i][j] = 0;
            } else if (i >= a.rows || j >= a.cols) {
                matrix[i][j] = -b.matrix[i + b.rowOffset][j + b.colOffset];
            } else if (i >= b.rows || j >= b.cols) {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset];
            } else {
                matrix[i][j] = a.matrix[i + a.rowOffset][j + a.colOffset] - b.matrix[i + b.rowOffset][j + b.colOffset];
            }
        }
    }

    _BlockMatrixD c = { dim, dim, 0, 0, matrix };
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
    size_t maxSize = 128;
    if (a.rows < maxSize || a.cols < maxSize || b.rows < maxSize || b.cols < maxSize) {
        return _multiplyBlockMatrixD(a, b);
    }

    size_t rows = (a.rows > b.rows) ? a.rows : b.rows; // maximum
    size_t cols = (a.cols > b.cols) ? a.cols : b.cols; // maximum
    size_t dim = (rows > cols) ? rows : cols;

    _BlockMatrixD A11 = { (dim + 1) / 2, (dim + 1) / 2, a.rowOffset, a.colOffset, a.matrix };
    _BlockMatrixD B11 = { (dim + 1) / 2, (dim + 1) / 2, b.rowOffset, b.colOffset, b.matrix };

    _BlockMatrixD A12 = { (dim + 1) / 2, a.cols - (dim + 1) / 2, a.rowOffset, a.colOffset + (dim + 1) / 2, a.matrix };
    _BlockMatrixD B12 = { (dim + 1) / 2, b.cols - (dim + 1) / 2, b.rowOffset, b.colOffset + (dim + 1) / 2, b.matrix };

    _BlockMatrixD A21 = { a.rows - (dim + 1) / 2, (dim + 1) / 2, a.rowOffset + (dim + 1) / 2, a.colOffset, a.matrix };
    _BlockMatrixD B21 = { b.rows - (dim + 1) / 2, (dim + 1) / 2, b.rowOffset + (dim + 1) / 2, b.colOffset, b.matrix };

    _BlockMatrixD A22 = { a.rows - (dim + 1) / 2, a.cols - (dim + 1) / 2, a.rowOffset + (dim + 1) / 2, a.colOffset + (dim + 1) / 2, a.matrix };
    _BlockMatrixD B22 = { b.rows - (dim + 1) / 2, b.cols - (dim + 1) / 2, b.rowOffset + (dim + 1) / 2, b.colOffset + (dim + 1) / 2, b.matrix };

    _BlockMatrixD M1a = _addBlockMatrixD(A11, A22);
    _BlockMatrixD M1b = _addBlockMatrixD(B11, B22);
    _BlockMatrixD M1 = _strassenMultiplyBlockMatrixD(M1a, M1b);
    _freeBlockMatrixD(&M1a);
    _freeBlockMatrixD(&M1b);

    _BlockMatrixD M2a = _addBlockMatrixD(A21, A22);
    _BlockMatrixD M2 = _strassenMultiplyBlockMatrixD(M2a, B11);
    _freeBlockMatrixD(&M2a);

    _BlockMatrixD M3a = _subtractBlockMatrixD(B12, B22);
    _BlockMatrixD M3 = _strassenMultiplyBlockMatrixD(A11, M3a);
    _freeBlockMatrixD(&M3a);

    _BlockMatrixD M4a = _subtractBlockMatrixD(B21, B11);
    _BlockMatrixD M4 = _strassenMultiplyBlockMatrixD(A22, M4a);
    _freeBlockMatrixD(&M4a);

    _BlockMatrixD M5a = _addBlockMatrixD(A11, A12);
    _BlockMatrixD M5 = _strassenMultiplyBlockMatrixD(M5a, B22);
    _freeBlockMatrixD(&M5a);

    _BlockMatrixD M6a = _subtractBlockMatrixD(A21, A11);
    _BlockMatrixD M6b = _addBlockMatrixD(B11, B12);
    _BlockMatrixD M6 = _strassenMultiplyBlockMatrixD(M6a, M6b);
    _freeBlockMatrixD(&M6a);
    _freeBlockMatrixD(&M6b);

    _BlockMatrixD M7a = _subtractBlockMatrixD(A12, A22);
    _BlockMatrixD M7b = _addBlockMatrixD(B21, B22);
    _BlockMatrixD M7 = _strassenMultiplyBlockMatrixD(M7a, M7b);
    _freeBlockMatrixD(&M7a);
    _freeBlockMatrixD(&M7b);

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

    _freeBlockMatrixD(&M1);
    _freeBlockMatrixD(&M2);
    _freeBlockMatrixD(&M3);
    _freeBlockMatrixD(&M4);
    _freeBlockMatrixD(&M5);
    _freeBlockMatrixD(&M6);
    _freeBlockMatrixD(&M7);

    _freeBlockMatrixD(&C11);
    _freeBlockMatrixD(&C12);
    _freeBlockMatrixD(&C21);
    _freeBlockMatrixD(&C22);

    return C;
}

_BlockMatrixD _multiplyBlockMatrixD(_BlockMatrixD a, _BlockMatrixD b) {
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

// Computes eigenvalues of nxn matrix a and stores them in nx2 matrix lambda
// Computes eigenvectors of nxn matrix a and stores them in nxn matrix v
// The first and second column of lambda are the real and imaginary parts of the eigenvalues
void eigensMatrixD(MatrixD *a, MatrixD *lambda, MatrixD *v) {
    if (a->rows != a->cols) {
        fprintf(stderr, "Error: Cannot compute eigenvalues of non-square matrix of shape (%lu, %lu)\n", a->rows, a->cols);
        exit(1);
    }

    if (a->rows != lambda->rows || lambda->cols != 2) {
        fprintf(stderr, "Error: Cannot store eigenvalue matrix of shape (%lu, 2) in matrix of shape (%lu, %lu)\n", a->rows, lambda->rows, lambda->cols);
    }

    if (v->rows != a->rows || v->cols != a->cols) {
        fprintf(stderr, "Error: Cannot store eigenvector matrix of shape (%lu, %lu) in matrix of shape (%lu, %lu)\n", a->rows, a->cols, v->rows, v->cols);
        exit(1);
    }

    MatrixD *d = malloc(sizeof(MatrixD));
    MatrixD *p = malloc(sizeof(MatrixD));
    qrAlgorithmMatrixD(*a, 100, d, p);

    MatrixD vi;
    MatrixD b;
    MatrixD lambda2 = newMatrixD(2, 2);

    size_t i = 0;
    while (i < a->rows) {
        if (i == a->rows - 1 || fabs(d->matrix[i+1][i]) < 1e-10) {
            // real eigenvalue
            vi = _submatrixMatrixD(*p, 0, i, p->rows, 1);
            
            lambda->matrix[i][0] = _eigenvalueIterationMatrixD(a, &vi, d->matrix[i][i]);
            lambda->matrix[i][1] = 0;

            for (size_t j = 0; j < vi.rows; ++j) {
                v->matrix[j][i] = vi.matrix[j][0];
            }

            freeMatrixD(&vi);
            i += 1;
        } else {
            // complex eigenvalue
            b = _submatrixMatrixD(*d, i, i, 2, 2);
            _twoByTwoEigenvaluesMatrixD(&b, &lambda2);

            for (size_t row = 0; row < 2; ++row) {
                for (size_t col = 0; col < 2; ++col) {
                    lambda->matrix[i + row][col] = lambda2.matrix[row][col];
                }
            }

            vi = newMatrixD(a->rows, 2);
            _complexEigenvectorsMatrixD(a, &vi, lambda2.matrix[0][0], lambda2.matrix[0][1]);

            for (size_t j = 0; j < vi.rows; ++j) {
                v->matrix[j][i] = vi.matrix[j][0];
                v->matrix[j][i + 1] = vi.matrix[j][1];
            }

            freeMatrixD(&b);
            freeMatrixD(&vi);
            i += 2;
        }
    }

    freeMatrixD(&lambda2);
    freeMatrixD(d);
    freeMatrixD(p);
    free(d);
    free(p);
}

void _twoByTwoEigenvaluesMatrixD(MatrixD *m, MatrixD *lambda) {
    if (m->rows != 2 || m->cols != 2) {
        fprintf(stderr, "Error: Cannot compute eigenvalues of matrix with shape (%lu, %lu) which is not 2x2\n", m->rows, m->cols);
        exit(1);
    }

    if (lambda->rows != 2 || lambda->cols != 2) {
        fprintf(stderr, "Error: Cannot store eigenvalues in matrix of shape (%lu, %lu) which is not 2x2\n", lambda->rows, lambda->cols);
        exit(1);
    }

    // x^2 - trace * x + det
    double trace = m->matrix[0][0] + m->matrix[1][1];
    double det = m->matrix[0][0] * m->matrix[1][1] - m->matrix[0][1] * m->matrix[1][0];
    double disc = trace * trace - 4 * det;
    if (disc >= 0) {
        // real eigenvalues
        lambda->matrix[0][0] = (trace + sqrt(disc)) / 2;
        lambda->matrix[0][1] = 0;
        lambda->matrix[1][0] = (trace - sqrt(disc)) / 2;
        lambda->matrix[1][1] = 0;
    } else {
        // complex eigenvalues
        lambda->matrix[0][0] = trace / 2;
        lambda->matrix[0][1] = sqrt(-disc) / 2;
        lambda->matrix[1][0] = trace / 2;
        lambda->matrix[1][1] = -sqrt(-disc) / 2;
    }
}

double _eigenvalueIterationMatrixD(MatrixD *a, MatrixD *vi, double lambda) {
    if (a->rows != a->cols) {
        fprintf(stderr, "Error: Cannot iterate eigenvalue of non-square matrix of shape (%lu, %lu)\n", a->rows, a->cols);
        exit(1);
    }

    if (vi->rows != a->rows || vi->cols != 1) {
        fprintf(stderr, "Error: Cannot iterate on eigenvector of shape (%lu, %lu)\n", vi->rows, vi->cols);
        exit(1);
    }

    size_t iterations = 100;
    MatrixD A = newMatrixD(a->rows, a->cols);
    MatrixD I = identityMatrixD(a->rows);
    MatrixD Ainv;

    scaleMatrixD(&A, -lambda, &I);
    addMatrixD(&A, a, &A);
    MatrixD *x = malloc(sizeof(MatrixD));
    MatrixD *N = malloc(sizeof(MatrixD));
    MatrixD zero = newMatrixD(a->rows, 1);
    int dims = solveLinearMatrixD(A, zero, x, N);
    if (dims >= 1) {
        for (size_t i = 0; i < vi->rows; ++i) {
            vi->matrix[i][0] = N->matrix[i][0]; // TODO: Handle eigenspace of dimension >1
        }
        scaleMatrixD(vi, 1 / normMatrixD(*vi, TWO_NORM), vi);
    } else {
        Ainv = inverseMatrixD(A);
        for (size_t k = 0; k < iterations; ++k) {
            multiplyMatrixD(vi, &Ainv, vi);
            scaleMatrixD(vi, 1 / normMatrixD(*vi, TWO_NORM), vi);
        }
        freeMatrixD(&Ainv);

        MatrixD avi = newMatrixD(a->rows, 1);
        MatrixD lambdaM = newMatrixD(1, 1);

        multiplyMatrixD(&avi, a, vi);
        innerProductMatrixD(&lambdaM, vi, &avi);
        lambda = lambdaM.matrix[0][0];

        freeMatrixD(&avi);
        freeMatrixD(&lambdaM);
    }

    freeMatrixD(&A);
    freeMatrixD(&I);
    freeMatrixD(&zero);
    freeMatrixD(x);
    freeMatrixD(N);
    free(x);
    free(N);

    return lambda;
}

void _complexEigenvectorsMatrixD(MatrixD *a, MatrixD *v, double real, double imag) {
    if (a->rows != a->cols) {
        fprintf(stderr, "Error: Cannot compute eigenvectors of non-square matrix of shape (%lu, %lu)\n", a->rows, a->cols);
        exit(1);
    }

    if (v->rows != a->rows || v->cols != 2) {
        fprintf(stderr, "Error: Cannot store eigenvector matrix of shape (%lu, 2) in matrix of shape (%lu, %lu)\n", a->rows, v->rows, v->cols);
        exit(1);
    }

    MatrixD I = identityMatrixD(a->rows);
    scaleMatrixD(&I, real * real + imag * imag, &I);
    MatrixD A2 = newMatrixD(a->rows, a->cols);
    scaleMatrixD(&A2, -2 * real, a);
    addMatrixD(&A2, &A2, &I);
    addMultiplyMatrixD(&A2, a, a);

    MatrixD *x = malloc(sizeof(MatrixD));
    MatrixD *U = malloc(sizeof(MatrixD));
    MatrixD zero = newMatrixD(a->rows, 1);

    int dims = solveLinearMatrixD(A2, zero, x, U);
    if (dims != 2) {
        fprintf(stderr, "Error: A2 has null space of dimension %d instead of dimension 2\n", dims);
        exit(1);
    }

    MatrixD ip = newMatrixD(2, 2);
    innerProductMatrixD(&ip, U, U);
    MatrixD B = inverseMatrixD(ip);

    MatrixD aU = newMatrixD(U->rows, U->cols);
    multiplyMatrixD(&aU, a, U);
    MatrixD UtaU = newMatrixD(2, 2);
    innerProductMatrixD(&UtaU, U, &aU);
    multiplyMatrixD(&B, &B, &UtaU);

    MatrixD P = newMatrixD(2, 2);
    P.matrix[0][0] = real - B.matrix[1][1];
    P.matrix[0][1] = imag;
    P.matrix[1][0] = B.matrix[1][0];
    P.matrix[1][1] = 0;

    multiplyMatrixD(v, U, &P);
    scaleMatrixD(v, 1 / normMatrixD(*v, TWO_NORM), v);

    freeMatrixD(&I);
    freeMatrixD(&A2);
    freeMatrixD(x);
    freeMatrixD(U);
    free(x);
    free(U);
    freeMatrixD(&zero);
    freeMatrixD(&ip);
    freeMatrixD(&B);
    freeMatrixD(&aU);
    freeMatrixD(&UtaU);
    freeMatrixD(&P);
}