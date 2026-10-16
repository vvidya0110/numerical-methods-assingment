/*
Droplet detection: peak search, FWHM width, spacing from previous droplet
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_LEN     8192
#define THRESH_FRAC  0.5   
#define VALLEY_FRAC  0.5   

static double lerp(double x1, double y1, double x2, double y2, double yt) {
    if (y2 == y1) return x1;
    return x1 + (yt - y1) * (x2 - x1) / (y2 - y1);
}

static int parse_two(const char *s, double *a, double *b) {
    char *e1, *e2;
    *a = strtod(s, &e1);
    if (e1 == s) return 0;
    while (*e1 == ',' || *e1 == ' ' || *e1 == '\t' || *e1 == ';') e1++;
    *b = strtod(e1, &e2);
    return e2 != e1;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <datafile>\n", argv[0]);
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) {
        fprintf(stderr, "Failed to open %s\n", argv[1]);
        return 1;
    }

    long cap = 100000, n = 0;
    double *x = malloc(cap * sizeof(double));
    double *y = malloc(cap * sizeof(double));
    if (!x || !y) { fprintf(stderr, "Memory allocation error.\n"); return 1; }

    char line[LINE_LEN];
    while (fgets(line, sizeof line, f)) {
        if (!strchr(line, '\n') && !feof(f)) {          
            int c;
            while ((c = fgetc(f)) != '\n' && c != EOF) ;
        }
        double a, b;
        if (!parse_two(line, &a, &b)) continue;          
        if (n == cap) {
            cap *= 2;
            double *nx = realloc(x, cap * sizeof(double));
            double *ny = realloc(y, cap * sizeof(double));
            if (!nx || !ny) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
            x = nx; y = ny;
        }
        x[n] = a; y[n] = b; n++;
    }
    fclose(f);
    if (n < 5) {
        fprintf(stderr, "Insufficient data points.\n");
        return 1;
    }

    double ymin = y[0], ymax = y[0];
    for (long i = 0; i < n; i++) {
        if (y[i] < ymin) ymin = y[i];
        if (y[i] > ymax) ymax = y[i];
    }
    double base = ymin;
    double thr  = base + THRESH_FRAC * (ymax - base);

    long *cand = malloc((n / 2 + 1) * sizeof(long));
    if (!cand) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
    long m = 0;
    for (long i = 1; i < n - 1; i++)
        if (y[i] > y[i - 1] && y[i] >= y[i + 1] && y[i] > thr)
            cand[m++] = i;
    if (m == 0) {
        printf("Number of droplets: 0\n");
        return 0;
    }

    long k = 1;
    for (long c = 1; c < m; c++) {
        long a = cand[k - 1], b = cand[c];
        double valley = y[a];
        for (long i = a; i <= b; i++)
            if (y[i] < valley) valley = y[i];
        double lowH = (y[a] < y[b] ? y[a] : y[b]) - base;
        if (valley - base > VALLEY_FRAC * lowH) {       
            if (y[b] > y[a]) cand[k - 1] = b;            
        } else {
            cand[k++] = b;                              
        }
    }

    double *pos = malloc(k * sizeof(double));
    double *wid = malloc(k * sizeof(double));
    if (!pos || !wid) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
    long cnt = 0;

    for (long j = 0; j < k; j++) {
        long p = cand[j];
        double half = base + 0.5 * (y[p] - base);

        long l = p;
        while (l > 0 && y[l] >= half) l--;
        int hasL = (y[l] < half);
        double xl = hasL ? lerp(x[l], y[l], x[l + 1], y[l + 1], half) : 0.0;

        long r = p;
        while (r < n - 1 && y[r] >= half) r++;
        int hasR = (y[r] < half);
        double xr = hasR ? lerp(x[r - 1], y[r - 1], x[r], y[r], half) : 0.0;

        double w;
        if (hasL && hasR) w = xr - xl;
        else if (hasL)    w = 2.0 * (x[p] - xl);         /* droplet cut off at right edge */
        else if (hasR)    w = 2.0 * (xr - x[p]);         /* droplet cut off at left edge */
        else continue;

        pos[cnt] = x[p];
        wid[cnt] = w;
        cnt++;
    }

    printf("Number of droplets: %ld\n", cnt);
    printf("%-14s %-14s %-14s\n", "Location", "Width", "Spacing");
    for (long i = 0; i < cnt; i++) {
        printf("%-14.4f %-14.4f ", pos[i], wid[i]);
        if (i == 0) printf("-\n");
        else        printf("%.4f\n", pos[i] - pos[i - 1]);
    }

    free(x); free(y); free(cand); free(pos); free(wid);
    return 0;
}
