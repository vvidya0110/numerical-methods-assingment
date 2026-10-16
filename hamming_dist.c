//haming distance
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// packs a string of '0'/'1' into bytes, MSB first; returns 0 if a bad character is found
static int pack(const char *s, unsigned char *out) {
    for (int i = 0; s[i]; i++) {
        if (s[i] != '0' && s[i] != '1') return 0;
        if (s[i] == '1') out[i / 8] |= 1 << (7 - i % 8);
    }
    return 1;
}

static int get_bit(const unsigned char *a, int i) {
    return (a[i / 8] >> (7 - i % 8)) & 1;
}

int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Usage: %s <stream_bits> <codeword_bits> [max_distance]\n", argv[0]);
        return 1;
    }

    int n = strlen(argv[1]), m = strlen(argv[2]);
    int maxd = (argc == 4) ? atoi(argv[3]) : 0;
    if (m == 0 || m > n) {
        fprintf(stderr, "Codeword must be non-empty and not longer than the stream.\n");
        return 1;
    }

    unsigned char *stream = calloc((n + 7) / 8, 1);
    unsigned char *code   = calloc((m + 7) / 8, 1);
    if (!stream || !code) { fprintf(stderr, "Memory allocation error.\n"); return 1; }
    if (!pack(argv[1], stream) || !pack(argv[2], code)) {
        fprintf(stderr, "Inputs must contain only 0 and 1.\n");
        return 1;
    }

    printf("%-10s %s\n", "Location", "Hamming distance");
    int best = m + 1, best_loc = -1, shown = 0;
    for (int s = 0; s + m <= n; s++) {
        int d = 0;
        for (int k = 0; k < m; k++)
            if (get_bit(stream, s + k) != get_bit(code, k)) d++;
        if (d < best) { best = d; best_loc = s; }
        if (d <= maxd) { printf("%-10d %d\n", s, d); shown++; }
    }

    if (shown == 0) printf("(no windows within distance %d)\n", maxd);
    printf("\nBest match: location %d, distance %d\n", best_loc, best);

    free(stream); free(code);
    return 0;
}
