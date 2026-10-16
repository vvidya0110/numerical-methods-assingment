//solving a system of linear equations using both gaussian elimination as well as LU decomposition
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 1e-12

double **alloc_matrix(int rows, int cols) {
    double **M = malloc(rows * sizeof(double *));
    if (!M) return NULL;
    for (int i = 0; i < rows; i++) {
        M[i] = malloc(cols * sizeof(double));
        if (!M[i]) {                       // free what was already allocated
            for (int k = 0; k < i; k++) free(M[k]);
            free(M);
            return NULL;
        }
    }
    return M;
}

void free_matrix(double **matrix, int rows) {
    if (!matrix) return;
    for (int i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

void swap_matrix_rows(double **A, int row1, int row2, int cols) {
    for (int j = 0; j < cols; j++) {
        double t = A[row1][j];
        A[row1][j] = A[row2][j];
        A[row2][j] = t;
    }
}

void print_solution(const double *x, int N) {
    printf("Solution:\n");
    for (int i = 0; i < N; i++)
        printf("  x[%d] = %.8f\n", i + 1, x[i]);
}

double **read_matrix_from_file(const char *filename, int N) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { printf("Error opening file: %s\n", filename); return NULL; }

    double **A = alloc_matrix(N, N + 1);
    if (!A) { fclose(fp); return NULL; }

    for (int i = 0; i < N; i++)
        for (int j = 0; j <= N; j++)
            if (fscanf(fp, "%lf", &A[i][j]) != 1) {
                printf("Error reading matrix.\n");
                fclose(fp);
                free_matrix(A, N);
                return NULL;
            }
    fclose(fp);
    return A;
}

// Gaussian elimination with partial pivoting.
// x holds b on entry and the solution on exit. Returns 0 on success, -1 if singular.
int gaussian_elimination(double **A, double *x, int N) {
    double **M = alloc_matrix(N, N + 1);
    if (!M) return -1;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) M[i][j] = A[i][j];
        M[i][N] = x[i];
    }

    for (int i = 0; i < N; i++) {                  // includes the last pivot
        int maxRow = i;
        for (int k = i + 1; k < N; k++)
            if (fabs(M[k][i]) > fabs(M[maxRow][i])) maxRow = k;
        if (fabs(M[maxRow][i]) < EPSILON) {
            free_matrix(M, N);
            return -1;                             // x left untouched
        }
        if (maxRow != i) swap_matrix_rows(M, i, maxRow, N + 1);

        for (int k = i + 1; k < N; k++) {
            double factor = M[k][i] / M[i][i];
            for (int j = i; j <= N; j++)
                M[k][j] -= factor * M[i][j];
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        double sum = M[i][N];
        for (int j = i + 1; j < N; j++) sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }

    free_matrix(M, N);
    return 0;
}

// Doolittle LU (unit diagonal L) with row pivoting.
// x holds b on entry and the solution on exit. Returns 0 on success, -1 if singular.
int lu_decomposition(double **A, double *x, int N) {
    double **LU = alloc_matrix(N, N);
    double *y = malloc(N * sizeof(double));
    int *perm = malloc(N * sizeof(int));
    if (!LU || !y || !perm) { free_matrix(LU, N); free(y); free(perm); return -1; }

    for (int i = 0; i < N; i++) {
        perm[i] = i;
        for (int j = 0; j < N; j++) LU[i][j] = A[i][j];
    }

    for (int i = 0; i < N; i++) {
        int maxRow = i;
        for (int k = i + 1; k < N; k++)
            if (fabs(LU[k][i]) > fabs(LU[maxRow][i])) maxRow = k;
        if (fabs(LU[maxRow][i]) < EPSILON) {
            free_matrix(LU, N); free(y); free(perm);
            return -1;
        }
        if (maxRow != i) {
            swap_matrix_rows(LU, i, maxRow, N);
            int t = perm[i]; perm[i] = perm[maxRow]; perm[maxRow] = t;
        }
        for (int k = i + 1; k < N; k++) {
            LU[k][i] /= LU[i][i];                  // L entry
            for (int j = i + 1; j < N; j++)
                LU[k][j] -= LU[k][i] * LU[i][j];   // U entries
        }
    }

    // forward substitution: L y = P b
    for (int i = 0; i < N; i++) {
        double sum = x[perm[i]];
        for (int j = 0; j < i; j++) sum -= LU[i][j] * y[j];
        y[i] = sum;
    }
    // backward substitution: U x = y
    for (int i = N - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < N; j++) sum -= LU[i][j] * x[j];
        x[i] = sum / LU[i][i];
    }

    free_matrix(LU, N); free(y); free(perm);
    return 0;
}

// largest |A x - b| entry, to check a solution
static double residual(double **A, const double *x, const double *b, int N) {
    double worst = 0.0;
    for (int i = 0; i < N; i++) {
        double s = -b[i];
        for (int j = 0; j < N; j++) s += A[i][j] * x[j];
        if (fabs(s) > worst) worst = fabs(s);
    }
    return worst;
}

static void solve_both(const char *title, double **A, const double *b, int N) {
    printf("=== %s (N = %d) ===\n", title, N);
    double *x = malloc(N * sizeof(double));
    if (!x) { printf("Memory allocation error.\n"); return; }

    memcpy(x, b, N * sizeof(double));
    printf("Gaussian elimination:\n");
    if (gaussian_elimination(A, x, N) == 0) {
        print_solution(x, N);
        printf("  max residual = %.2e\n", residual(A, x, b, N));
    } else printf("  Matrix is singular.\n");

    memcpy(x, b, N * sizeof(double));
    printf("LU decomposition:\n");
    if (lu_decomposition(A, x, N) == 0) {
        print_solution(x, N);
        printf("  max residual = %.2e\n", residual(A, x, b, N));
    } else printf("  Matrix is singular.\n");
    printf("\n");
    free(x);
}

int main(int argc, char **argv) {
    if (argc == 3) {
        int N = atoi(argv[2]);
        if (N < 1) { fprintf(stderr, "N must be >= 1\n"); return 1; }
        double **aug = read_matrix_from_file(argv[1], N);
        if (!aug) return 1;

        double **A = alloc_matrix(N, N);
        double *b = malloc(N * sizeof(double));
        if (!A || !b) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) A[i][j] = aug[i][j];
            b[i] = aug[i][N];
        }
        solve_both(argv[1], A, b, N);
        free_matrix(aug, N); free_matrix(A, N); free(b);
        return 0;
    }

    // demo 1: ordinary 3x3 system, solution (2, 3, -1)
    double d1[3][3] = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    double b1[3] = {8, -11, -3};
    double **A1 = alloc_matrix(3, 3);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) A1[i][j] = d1[i][j];
    solve_both("Demo 1: 3x3 system", A1, b1, 3);
    free_matrix(A1, 3);

    // demo 2: zero on the diagonal, solution (5, 3); the old LU divided by zero here
    double d2[2][2] = {{0, 1}, {1, 0}};
    double b2[2] = {3, 5};
    double **A2 = alloc_matrix(2, 2);
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) A2[i][j] = d2[i][j];
    solve_both("Demo 2: zero pivot", A2, b2, 2);
    free_matrix(A2, 2);

    // demo 3: singular matrix, should be reported rather than returning garbage
    double d3[2][2] = {{1, 2}, {2, 4}};
    double b3[2] = {3, 6};
    double **A3 = alloc_matrix(2, 2);
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) A3[i][j] = d3[i][j];
    solve_both("Demo 3: singular matrix", A3, b3, 2);
    free_matrix(A3, 2);

    return 0;
}
