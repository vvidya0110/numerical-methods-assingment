//fitting data using cubic spline 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Cubic spline interpolation with natural boundary conditions.
// x must be strictly increasing. Returns 0.0 on invalid input or allocation failure.
double cubic_spline_interpolate(const double *x, const double *y, int M, double xp) {
    if (M < 2) return 0.0;

    double *h    = malloc((M - 1) * sizeof(double));
    double *diag = malloc(M * sizeof(double));
    double *m    = malloc(M * sizeof(double));   // second derivatives
    if (!h || !diag || !m) { free(h); free(diag); free(m); return 0.0; }

    for (int i = 0; i < M - 1; i++) {
        h[i] = x[i + 1] - x[i];
        if (h[i] <= 0.0) { free(h); free(diag); free(m); return 0.0; }
    }

    // natural boundary conditions
    m[0] = 0.0;
    m[M - 1] = 0.0;

    // interior system: h[i-1]*m[i-1] + 2(h[i-1]+h[i])*m[i] + h[i]*m[i+1] = rhs[i]
    for (int i = 1; i < M - 1; i++) {
        diag[i] = 2.0 * (h[i - 1] + h[i]);
        m[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }

    // Thomas algorithm: forward elimination
    for (int i = 2; i < M - 1; i++) {
        double w = h[i - 1] / diag[i - 1];
        diag[i] -= w * h[i - 1];
        m[i]    -= w * m[i - 1];
    }
    // back substitution (only if there are interior points)
    if (M > 2) {
        m[M - 2] /= diag[M - 2];
        for (int i = M - 3; i >= 1; i--)
            m[i] = (m[i] - h[i] * m[i + 1]) / diag[i];
    }

    // find the interval containing xp (binary search, clamped at the ends)
    int lo = 0, hi = M - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (xp >= x[mid]) lo = mid; else hi = mid;
    }
    int idx = lo;

    double hi_ = h[idx];
    double a = (x[idx + 1] - xp) / hi_;
    double b = (xp - x[idx]) / hi_;

    double S = a * y[idx] + b * y[idx + 1] +
               ((a * a * a - a) * m[idx] + (b * b * b - b) * m[idx + 1]) * (hi_ * hi_) / 6.0;

    free(h); free(diag); free(m);
    return S;
}

int main(void) {
    // sample sin(x) at 7 points on [0, pi]
    double pi = 3.14159265358979323846;
    int M = 7;
    double x[7], y[7];
    for (int i = 0; i < M; i++) {
        x[i] = pi * i / (M - 1);
        y[i] = sin(x[i]);
    }

    // spline should match the data points exactly
    printf("At the data points (error should be ~0):\n");
    for (int i = 0; i < M; i++)
        printf("  x=%.4f  spline=%.6f  sin=%.6f\n",
               x[i], cubic_spline_interpolate(x, y, M, x[i]), y[i]);

    // and stay close to sin between the points
    printf("Between the data points:\n");
    for (double xp = 0.1; xp < pi; xp += 0.5)
        printf("  x=%.4f  spline=%.6f  sin=%.6f  err=%.2e\n",
               xp, cubic_spline_interpolate(x, y, M, xp), sin(xp),
               fabs(cubic_spline_interpolate(x, y, M, xp) - sin(xp)));
    return 0;
}
