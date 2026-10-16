//command line input four arguments: mean1, mean2, stdev and N
//Generate 2 arrays of N uniformly distributed random numbers 
//Box-Muller transformation to convert it into 2 arrays of numbers in a normal distribution
//The distribution is shifted so they have their averages at mean1 and mean2
//Student t- test and value of t is printed.




#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define PI 3.141592653589793

static double mean_of(double *a, int n) {
    double s = 0; for (int i = 0; i < n; i++) s += a[i]; return s / n;
}
static double var_of(double *a, int n, double m) {
    double s = 0; for (int i = 0; i < n; i++) s += (a[i]-m)*(a[i]-m); return s / (n-1);
}
static double unif(void) { return (rand() + 1.0) / (RAND_MAX + 2.0); } 
int main(int argc, char *argv[]) {
    if (argc != 5) { fprintf(stderr, "Usage: %s mean1 mean2 stdev N\n", argv[0]); return 1; }
    double mean1 = atof(argv[1]), mean2 = atof(argv[2]), stdev = atof(argv[3]);
    int N = atoi(argv[4]);
    if (N < 2) { fprintf(stderr, "N must be >= 2\n"); return 1; }

    srand((unsigned)time(NULL) ^ (unsigned)N * 2654435761u);

    double *u1 = malloc(N * sizeof(double)), *u2 = malloc(N * sizeof(double));
    double *x = malloc(N * sizeof(double)),  *y = malloc(N * sizeof(double));
    if (!u1 || !u2 || !x || !y) { fprintf(stderr, "malloc failed\n"); return 1; }

    for (int i = 0; i < N; i++) { u1[i] = unif(); u2[i] = unif(); }

    for (int i = 0; i < N; i++) {              
        double r = sqrt(-2.0 * log(u1[i]));
        double th = 2.0 * PI * u2[i];
        x[i] = mean1 + stdev * r * cos(th);    
        y[i] = mean2 + stdev * r * sin(th);    
    }

    double m1 = mean_of(x, N), m2 = mean_of(y, N);
    double v1 = var_of(x, N, m1), v2 = var_of(y, N, m2);

    double sp2 = ((N-1)*v1 + (N-1)*v2) / (2.0*N - 2.0);
    double t = (m1 - m2) / sqrt(sp2 * (2.0 / N));
    printf("t = %f\n", t);

    FILE *fp = fopen("t_vs_N.dat", "a");
    if (fp) { fprintf(fp, "%d %.6f\n", N, t); fclose(fp); }

    free(u1); free(u2); free(x); free(y);
    return 0;
}
