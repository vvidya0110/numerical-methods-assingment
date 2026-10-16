//for a particular temperature, the average of all the neighbouring spins is taken and based on probability (from the energy calculation) flips the atom. Done till the magnetisation remains constant for a temperature. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int **alloc_ising_grid(int N) {
    int **grid = malloc(N * sizeof(int *));
    if (!grid) return NULL;
    for (int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(int));
        if (!grid[i]) {
            for (int k = 0; k < i; k++) free(grid[k]);
            free(grid);
            return NULL;
        }
    }
    return grid;
}

void free_ising_grid(int **grid, int N) {
    if (!grid) return;
    for (int i = 0; i < N; i++) free(grid[i]);
    free(grid);
}

void initialize_random_spins(int **M, int N) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            M[i][j] = (rand() % 2) * 2 - 1;   // +1 or -1
}

// sum of the 4 neighbours, with wrap-around at the edges
int neighbor_sum(int **M, int i, int j, int N) {
    int top    = (i == 0) ? N - 1 : i - 1;
    int bottom = (i == N - 1) ? 0 : i + 1;
    int left   = (j == 0) ? N - 1 : j - 1;
    int right  = (j == N - 1) ? 0 : j + 1;
    return M[top][j] + M[bottom][j] + M[i][left] + M[i][right];
}

double calculate_magnetization(int **M, int N) {
    long sum = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            sum += M[i][j];
    return (double)sum / ((double)N * N);
}

// one sweep, in place: each spin sees its neighbours' latest values
void metropolis_sweep(int **M, int N, double J, double T) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double deltaE = 2.0 * J * M[i][j] * neighbor_sum(M, i, j, N);
            if (deltaE <= 0 || rand() / (RAND_MAX + 1.0) < exp(-deltaE / T))
                M[i][j] = -M[i][j];
        }
}

int ising_simulation(int N, double J, double Tstart, double Tstop, double dT,
                     double threshold, int patience, int max_iter) {
    int **M = alloc_ising_grid(N);
    if (!M) { fprintf(stderr, "Memory allocation error.\n"); return -1; }
    initialize_random_spins(M, N);

    FILE *fptr = fopen("data.dat", "w");
    if (!fptr) {
        fprintf(stderr, "Error opening data.dat\n");
        free_ising_grid(M, N);
        return -1;
    }

    int nT = (int)((Tstop - Tstart) / dT + 0.5);
    printf("%-12s %-12s %s\n", "Temperature", "<|m|>", "sweeps");

    for (int k = 0; k <= nT; k++) {
        double T = Tstart + k * dT;

        double prev = fabs(calculate_magnetization(M, N));
        int stable = 0, cnt = 0, iter;
        double acc = 0.0;

        for (iter = 1; iter <= max_iter; iter++) {
            metropolis_sweep(M, N, J, T);
            double m = fabs(calculate_magnetization(M, N));

            if (fabs(m - prev) < threshold) {
                stable++;                 // magnetisation unchanged this sweep
                acc += m; cnt++;
            } else {
                stable = 0;               // changed: start counting again
                acc = 0.0; cnt = 0;
            }
            prev = m;
            if (stable >= patience) break; // constant for `patience` sweeps in a row
        }

        // average over the stable window (or the last value if it never settled)
        double avg_mag = (cnt > 0) ? acc / cnt : prev;

        fprintf(fptr, "%.6f,%.6f\n", T, avg_mag);
        printf("%-12.4f %-12.4f %d%s\n", T, avg_mag, iter > max_iter ? max_iter : iter,
               iter > max_iter ? "  (hit max_iter)" : "");
    }
    fclose(fptr);
    free_ising_grid(M, N);

    FILE *plot = popen("gnuplot", "w");
    if (!plot) {
        fprintf(stderr, "gnuplot not found; data is in data.dat\n");
        return 0;
    }
    fprintf(plot, "set terminal jpeg size 1200,800\n");
    fprintf(plot, "set output 'ising_simulation.jpg'\n");
    fprintf(plot, "set datafile separator ','\n");
    fprintf(plot, "set xlabel 'Temperature'\n");
    fprintf(plot, "set ylabel 'Magnetization |m|'\n");
    fprintf(plot, "unset key\n");
    fprintf(plot, "plot 'data.dat' using 1:2 with linespoints\n");
    fprintf(plot, "set output\n");
    pclose(plot);
    return 0;
}

int main(int argc, char **argv) {
    int    N         = (argc > 1) ? atoi(argv[1]) : 32;
    double J         = (argc > 2) ? atof(argv[2]) : 1.0;
    double Tstart    = (argc > 3) ? atof(argv[3]) : 1.0;
    double Tstop     = (argc > 4) ? atof(argv[4]) : 4.0;
    double dT        = (argc > 5) ? atof(argv[5]) : 0.1;
    double threshold = (argc > 6) ? atof(argv[6]) : 0.005;
    int    patience  = (argc > 7) ? atoi(argv[7]) : 10;
    int    max_iter  = (argc > 8) ? atoi(argv[8]) : 5000;

    if (N < 2 || Tstart <= 0 || Tstop < Tstart || dT <= 0 ||
        threshold <= 0 || patience < 1 || max_iter < 1) {
        fprintf(stderr, "Usage: %s [N>=2] [J] [Tstart>0] [Tstop>=Tstart] [dT>0] "
                        "[threshold>0] [patience>=1] [max_iter>=1]\n", argv[0]);
        return 1;
    }

    srand((unsigned)time(NULL));
    return ising_simulation(N, J, Tstart, Tstop, dT, threshold, patience, max_iter) == 0 ? 0 : 1;
}
