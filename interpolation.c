//comparing different types of graph interpolation to find the final value
//uses that to find (x,y) at a specific location, when data is noisy 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ---------- Lagrange ----------
double lagrange_interpolate(const double *x, const double *y, int M, double xp) {
    double result = 0.0;
    for (int i = 0; i < M; i++) {
        double term = y[i];
        for (int j = 0; j < M; j++) {
            if (i != j)
                term *= (xp - x[j]) / (x[i] - x[j]);
        }
        result += term;
    }
    return result;
}

// ---------- Newton ----------
void newton_build_table(const double *x, const double *y, int M, double **table) {
    for (int i = 0; i < M; i++)
        table[i][0] = y[i];
    for (int j = 1; j < M; j++)
        for (int i = 0; i < M - j; i++)
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
}

double newton_interpolate(const double *x, double **table, int M, double xp) {
    double result = table[0][0];
    for (int i = 1; i < M; i++) {
        double term = table[0][i];
        for (int j = 0; j < i; j++)
            term *= (xp - x[j]);
        result += term;
    }
    return result;
}

// ---------- Natural cubic spline (tridiagonal solve) ----------
double cubic_spline_interpolate(const double *x, const double *y, int M, double xp) {
    if (M < 2) return 0.0;
    double *h    = malloc((M - 1) * sizeof(double));
    double *diag = malloc(M * sizeof(double));
    double *m    = malloc(M * sizeof(double));
    if (!h || !diag || !m) { free(h); free(diag); free(m); return 0.0; }

    for (int i = 0; i < M - 1; i++) {
        h[i] = x[i + 1] - x[i];
        if (h[i] <= 0.0) { free(h); free(diag); free(m); return 0.0; }
    }
    m[0] = 0.0;
    m[M - 1] = 0.0;
    for (int i = 1; i < M - 1; i++) {
        diag[i] = 2.0 * (h[i - 1] + h[i]);
        m[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }
    for (int i = 2; i < M - 1; i++) {
        double w = h[i - 1] / diag[i - 1];
        diag[i] -= w * h[i - 1];
        m[i]    -= w * m[i - 1];
    }
    if (M > 2) {
        m[M - 2] /= diag[M - 2];
        for (int i = M - 3; i >= 1; i--)
            m[i] = (m[i] - h[i] * m[i + 1]) / diag[i];
    }

    int lo = 0, hi = M - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (xp >= x[mid]) lo = mid; else hi = mid;
    }
    double hh = h[lo];
    double a = (x[lo + 1] - xp) / hh;
    double b = (xp - x[lo]) / hh;
    double S = a * y[lo] + b * y[lo + 1] +
               ((a * a * a - a) * m[lo] + (b * b * b - b) * m[lo + 1]) * (hh * hh) / 6.0;

    free(h); free(diag); free(m);
    return S;
}

double polynomial_eval(double x, const double *coeffs, int N) {
    double result = 0.0;
    double xn = 1.0;
    for (int i = 0; i <= N; i++) {
        result += coeffs[i] * xn;
        xn *= x;
    }
    return result;
}

double gaussian_noise(double stdev) {
    static int hasSpare = 0;
    static double spare;
    if (hasSpare) {
        hasSpare = 0;
        return stdev * spare;
    }
    hasSpare = 1;
    double u, v, s;
    do {
        u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1 || s == 0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return stdev * (u * s);
}

double calculate_rmse(const double *y_true, const double *y_pred, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double diff = y_true[i] - y_pred[i];
        sum += diff * diff;
    }
    return sqrt(sum / N);
}


int main(int argc, char **argv) {
    int M = (argc > 1) ? atoi(argv[1]) : 8;           
    double noise = (argc > 2) ? atof(argv[2]) : 0.1;  
    if (M < 3 || noise < 0) {
        fprintf(stderr, "Usage: %s [M >= 3] [noise_stdev >= 0]\n", argv[0]);
        return 1;
    }
    srand((unsigned)time(NULL));


    double coeffs[] = {1.0, -2.0, 0.5, 0.3};
    int deg = 3;
    double x_lo = -3.0, x_hi = 3.0;


    double *x = malloc(M * sizeof(double));
    double *y = malloc(M * sizeof(double));
    double **table = malloc(M * sizeof(double *));
    if (!x || !y || !table) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
    for (int i = 0; i < M; i++) {
        table[i] = malloc(M * sizeof(double));
        if (!table[i]) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
        x[i] = x_lo + (x_hi - x_lo) * i / (M - 1);
        y[i] = polynomial_eval(x[i], coeffs, deg) + gaussian_noise(noise);
    }
    newton_build_table(x, y, M, table);


    int T = 100;
    double *y_true = malloc(T * sizeof(double));
    double *y_lag  = malloc(T * sizeof(double));
    double *y_new  = malloc(T * sizeof(double));
    double *y_spl  = malloc(T * sizeof(double));
    if (!y_true || !y_lag || !y_new || !y_spl) { fprintf(stderr, "Memory allocation error.\n"); return 1; }

    printf("M = %d samples, noise stdev = %.3f\n\n", M, noise);
    printf("%-10s %-12s %-12s %-12s %-12s\n", "x", "true", "lagrange", "newton", "spline");
    for (int t = 0; t < T; t++) {
        double xp = x_lo + (x_hi - x_lo) * t / (T - 1);
        y_true[t] = polynomial_eval(xp, coeffs, deg);
        y_lag[t]  = lagrange_interpolate(x, y, M, xp);
        y_new[t]  = newton_interpolate(x, table, M, xp);
        y_spl[t]  = cubic_spline_interpolate(x, y, M, xp);
        if (t % 10 == 0)
            printf("%-10.4f %-12.5f %-12.5f %-12.5f %-12.5f\n",
                   xp, y_true[t], y_lag[t], y_new[t], y_spl[t]);
    }

    printf("\nRMSE vs true polynomial:\n");
    printf("  Lagrange : %.6f\n", calculate_rmse(y_true, y_lag, T));
    printf("  Newton   : %.6f\n", calculate_rmse(y_true, y_new, T));
    printf("  Spline   : %.6f\n", calculate_rmse(y_true, y_spl, T));

    for (int i = 0; i < M; i++) free(table[i]);
    free(table); free(x); free(y);
    free(y_true); free(y_lag); free(y_new); free(y_spl);
    return 0;
}
