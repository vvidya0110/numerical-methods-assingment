/*
 * Balls hopping randomly between n places in a line.
 * Start: places 1..n/2 hold p balls each, the rest hold 0.
 * Every step, each ball moves left or right with probability 0.5
 * (at the two ends, a ball moves inward with probability 0.5, else stays).
 * Stops when the "difference" measure drops below w.
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double rnd(void) { return rand() / (RAND_MAX + 1.0); }

static double avg_change(const int *a, const int *b, int n) {
    long long sum = 0;
    for (int m = 0; m < n; m++)
        sum += abs(a[m] - b[m]);
    return (double)sum / n;
}

static double pair_difference(const int *a, int n) {
    long long sum = 0;
    for (int s = 0; s < n - 1; s++)
        for (int i = s + 1; i < n; i++)
            sum += abs(a[s] - a[i]);
    return 2.0 * sum / ((double)(n - 1) * (n - 2));
}

int main(int argc, char **argv) {
    int    n         = argc > 1 ? atoi(argv[1]) : 20;
    int    p         = argc > 2 ? atoi(argv[2]) : 1000;
    double w         = argc > 3 ? atof(argv[3]) : 20.0;
    long   max_steps = argc > 4 ? atol(argv[4]) : 1000000;

    if (n < 3) { fprintf(stderr, "n must be at least 3\n"); return 1; }

    srand((unsigned)time(NULL));

    int *cur  = calloc(n, sizeof(int));   
    int *next = calloc(n, sizeof(int));  
    int *prev = calloc(n, sizeof(int));   
    if (!cur || !next || !prev) { fprintf(stderr, "out of memory\n"); return 1; }

    for (int i = 0; i < n / 2; i++) cur[i] = p;   
    double difference = 100.0;
    long   step = 0;

    while (difference >= w && step < max_steps) {
        for (int i = 0; i < n; i++) prev[i] = cur[i];   
        for (int i = 0; i < n; i++) next[i] = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < cur[i]; j++) {
                double y = rnd();
                if (i == 0) {                    
                    if (y > 0.5) next[i + 1]++; else next[i]++;
                } else if (i == n - 1) {         
                    if (y < 0.5) next[i - 1]++; else next[i]++;
                } else {                         
                    if (y > 0.5) next[i + 1]++; else next[i - 1]++;
                }
            }
        }
        for (int i = 0; i < n; i++) cur[i] = next[i];


        double change = avg_change(cur, prev, n);   
        (void)change;
        difference = pair_difference(cur, n);       
        step++;
    }

    printf("Stopped after %ld steps, difference = %.4f\n", step, difference);
    for (int i = 0; i < n; i++) printf("x%d = %d\n", i + 1, cur[i]);

    free(cur); free(next); free(prev);
    return 0;
}
