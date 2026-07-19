#ifndef GRANDPRE_H
#define GRANDPRE_H
#include "colossus.h"

// =====================================================================
//  Grandpre cipher primitive (grandpre.c)
// =====================================================================
//
// An ACA "Grandpre" (150-200 plaintext letters). An N x N square (N in 6..10, 8x8
// standard and preferred) is filled with N-letter words written horizontally so that
// EVERY ROW is a word and the FIRST COLUMN spells a further (N+1'th) word; ALL 26
// letters appear at least once. Each plaintext letter is enciphered by a 2-DIGIT
// number -- the (row, column) coordinates of ANY cell holding that letter. Because a
// letter usually occupies several cells, it has several homophones; the encoder picks
// one freely (frequency flattening). Rows/columns are numbered 1..N (0..9 for N=10),
// so the ciphertext is a stream of 2-digit codes.
//
//   Square (ACA 8x8 worked example; first column spells LACQUERS):
//        1 2 3 4 5 6 7 8
//     1  L A D Y B U G S
//     2  A Z I M U T H S
//     3  C A L F S K I N
//     4  Q U A C K I S H
//     5  U N J O V I A L
//     6  E V U L S I O N
//     7  R O W D Y I S M
//     8  S E X T U P L Y
//
//   pt "thefirstcolumn" -> 84 27 82 34 56 71 77 26 44 54 64 63 78
//
// DECODING IS UNIQUE (code -> the single letter in that cell); ENCODING is polyphonic
// (a letter -> any of its cells). So cryptanalytically Grandpre is a HOMOPHONIC
// substitution over <= N^2 numeric codes -> 26 letters (see grandpre_solver.c). Row/col
// labels are fixed (the grid indices), so -- unlike the Checkerboard / Syllabary -- the
// actual square is directly recoverable, not merely up to a label permutation. The
// word structure of the rows/first-column is a constraint the ACA solver exploits but
// is NOT required for a blind n-gram attack, so this primitive does not enforce it.
//
// The square is built once by grandpre_build() (or grandpre_build_from_words());
// encrypt / the generator / the tests share it, so the convention can never drift.

#define GRANDPRE_MIN_SIDE 6
#define GRANDPRE_MAX_SIDE 10
#define GRANDPRE_MAX_GRID 100        // GRANDPRE_MAX_SIDE^2

// A built square: the letter grid plus the reverse index (letter -> its homophone codes).
typedef struct {
    int N;                                   // side (6..10)
    int cell[GRANDPRE_MAX_GRID];             // row-major: cell[r*N + c] = letter 0..25
    int ncodes[ALPHABET_SIZE];               // homophone count per letter (>= 1 for all 26)
    int codes[ALPHABET_SIZE][GRANDPRE_MAX_GRID]; // letter -> its 2-digit codes (row-major order)
} GrandpreSquare;

// The digit label of a 0-based row/column index, and its inverse. Labels are 1..N for a
// 6..9 square and 0..9 for the 10x10 (ACA convention). grandpre_index returns -1 for a
// label outside the square.
static inline int grandpre_label(int idx, int N) { return (N >= 10) ? idx : idx + 1; }
static inline int grandpre_index(int label, int N) {
    int idx = (N >= 10) ? label : label - 1;
    return (idx >= 0 && idx < N) ? idx : -1;
}
// The 2-digit ciphertext code of cell (r, c): the two labels concatenated as base-10 digits.
static inline int grandpre_code(int r, int c, int N) {
    return grandpre_label(r, N) * 10 + grandpre_label(c, N);
}

// Build the square from the row-major letter arrangement `letters[0..N*N-1]` (each 0..25).
// Requires N in [6,10] and ALL 26 letters present. Fills the reverse index. Returns 0 on
// success, -1 on a bad side, an out-of-range letter, or a missing letter.
int grandpre_build(GrandpreSquare *g, int N, const int letters[]);

// Convenience builder from N words of length N (letters only, case-insensitive), one per
// row. Returns 0 on success, -1 if a word is not exactly N letters or a letter is missing.
int grandpre_build_from_words(GrandpreSquare *g, int N, const char *const words[]);

// Encipher plaintext `plain[0..n-1]` into the 2-digit code stream out_codes[]. Each letter
// (index 0..25) becomes the code of one of its cells, chosen UNIFORMLY AT RANDOM among the
// homophones (rand_int -- seed for reproducibility); any non-letter symbol is skipped.
// Returns the number of codes written (== the number of letters in plain[]).
int grandpre_encrypt(const int plain[], int n, const GrandpreSquare *g, int out_codes[]);

// Decipher the 2-digit code stream `codes[0..ncodes-1]` into out_pt[] (unique: code -> the
// letter in that cell). An out-of-range code (a label outside the square) is skipped.
// Returns the number of plaintext letters written.
int grandpre_decrypt(const int codes[], int ncodes, const GrandpreSquare *g, int out_pt[]);

#endif
