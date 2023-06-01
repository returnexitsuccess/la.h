#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

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