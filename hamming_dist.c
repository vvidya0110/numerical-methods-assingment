//finding the hamming distance in a array of 0's and 1's 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

// -------------------------------
// GET BIT
// -------------------------------
int get_bit(const unsigned char *arr, int bit_index) {
    int byte = bit_index / 8;
    int bit  = bit_index % 8;
    return (arr[byte] >> (7 - bit)) & 1;
}

// -------------------------------
// SET BIT
// -------------------------------
void set_bit(unsigned char *arr, int bit_index, int value) {
    int byte = bit_index / 8;
    int bit  = bit_index % 8;
    if (value)
        arr[byte] |= (1 << (7 - bit));
    else
        arr[byte] &= ~(1 << (7 - bit));
}

// -------------------------------
// FLIP BIT
// -------------------------------
void flip_bit(unsigned char *arr, int bit_index) {
    int byte = bit_index / 8;
    int bit  = bit_index % 8;
    arr[byte] ^= (1 << (7 - bit));
}

// -------------------------------
// GENERATE RANDOM BITSTREAM
// -------------------------------
void generate_random_bits(unsigned char *arr, int num_bits) {
    int num_bytes = (num_bits + 7) / 8;
    for (int i = 0; i < num_bytes; i++)
        arr[i] = rand() & 0xFF;

    // Clear extra bits at end
    int extra = num_bytes * 8 - num_bits;
    if (extra > 0) {
        arr[num_bytes - 1] &= (0xFF << extra);
    }
}

// -------------------------------
// FLIP WITH PROBABILITY p
// -------------------------------
void flip_bits_with_probability(unsigned char *arr, int num_bits, double probability) {
    for (int i = 0; i < num_bits; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < probability)
            flip_bit(arr, i);
    }
}

// -------------------------------
// STRING → BIT ARRAY
// -------------------------------
void string_to_bits(const char *str, unsigned char *arr, int length) {
    for (int i = 0; i < length; i++) {
        set_bit(arr, i, str[i] == '1');
    }
}

// -------------------------------
// WRITE BITS TO FILE
// -------------------------------
void write_bits_to_file(const char *filename,
                        const unsigned char *arr,
                        int num_bits) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("File write error");
        return;
    }
    int num_bytes = (num_bits + 7) / 8;
    fwrite(arr, 1, num_bytes, fp);
    fclose(fp);
}

// -------------------------------
// READ BITS FROM FILE
// -------------------------------
void read_bits_from_file(const char *filename,
                         unsigned char *arr,
                         int num_bits) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("File read error");
        return;
    }
    int num_bytes = (num_bits + 7) / 8;
    fread(arr, 1, num_bytes, fp);
    fclose(fp);
}

// -------------------------------
// HAMMING DISTANCE
// -------------------------------
int hamming_distance(const unsigned char *bits1,
                     const unsigned char *bits2,
                     int num_bits)
{
    int dist = 0;
    for (int i = 0; i < num_bits; i++) {
        if (get_bit(bits1, i) != get_bit(bits2, i))
            dist++;
    }
    return dist;
}

// -------------------------------
// SLIDING WINDOW CODEWORD SEARCH
// -------------------------------
int find_codeword(const unsigned char *bitstream,
                  int stream_length,
                  const unsigned char *codeword,
                  int codeword_length,
                  HammingResult *results,
                  int max_results)
{
    int num_windows = stream_length - codeword_length + 1;
    if (num_windows <= 0) return 0;

    int count = 0;

    for (int start = 0; start < num_windows; start++) {

        int hd = 0;
        for (int k = 0; k < codeword_length; k++) {
            if (get_bit(bitstream, start + k) != get_bit(codeword, k))
                hd++;
        }

        if (count < max_results) {
            results[count].location = start;
            results[count].distance = hd;
        }
        count++;
    }

    return (count < max_results ? count : max_results);
}

