#ifndef TRIDIGITAL_H
#define TRIDIGITAL_H
#include "colossus.h"

// =====================================================================
//  Tridigital cipher primitive (tridigital.c)
// =====================================================================
//
// An ACA "Tridigital" (75-100 plaintext letters). A 10-letter keyword produces a
// numerical key placed above a block 10 columns wide; a 26-letter keyed alphabet is
// written ROW-MAJOR into the block leaving the LAST column blank. Each plaintext letter
// is enciphered by the DIGIT heading its column; the digit above the (blank) last column
// is the WORD SEPARATOR. So the ciphertext is a stream of digits 0-9.
//
//   Keyword A (10 distinct letters) -> the 10 column digit-labels by alphabetical rank
//     (rank 1..10, with 10 written as 0):  NOVELCRAFT -> 6 7 0 3 5 2 8 1 4 9
//   Keyword B -> the 26-letter keyed alphabet (keyword deduped, then the unused letters
//     A..Z in order; NO J->I merge):       DRAGONFLY -> D R A G O N F L Y B C E H I ...
//
//   Block (ACA worked example):
//        6 7 0 3 5 2 8 1 4 9      column digit-labels (col 9 = separator = 9)
//        D R A G O N F L Y -      row 0  (cols 0..8)
//        B C E H I J K M P -      row 1  (cols 0..8)
//        Q S T U V W X Z - -      row 2  (cols 0..7)
//
//   pt "the ides of march" -> 0 3 0 9 5 6 0 7 9 5 8 9 1 0 7 7 3   (9 = word separator)
//
// The keyed alphabet fills the 9 letter-columns (0..8) row-major: row0 cols0..8, row1
// cols0..8, row2 cols0..7 (26 letters; row2 col8 is blank). So columns 0..7 hold 3
// letters and column 8 holds 2. The 10th physical column (index 9) holds NO letter -- its
// digit is the separator. Because a letter's ROW does not affect its digit, the cipher is
// deliberately AMBIGUOUS: a single digit stands for any of the (up to 3) letters in its
// column -- decryption is a language-model problem, not a unique inverse (the solver's job).
//
// The grid is built once by tridigital_build_grid() (or tridigital_build_from_keywords());
// encrypt / the generator / the tests share it, so the convention can never drift.

#define TRI_NALPHA       26    // full A..Z (no J->I merge)
#define TRI_COLS         10    // 10 physical columns; column 9 is the blank separator column
#define TRI_LETTER_COLS  9     // columns 0..8 bear letters (row-major fill 9 + 9 + 8 = 26)
#define TRI_MAXCOLLET    3     // max letters in a column (cols 0..7 hold 3, col 8 holds 2)

// The space / word-break sentinel in an index stream: index_to_char(TRI_SPACE) == ' '
// (see colossus.h). A negative value, so ngram_score() treats it as transparent.
#define TRI_SPACE  (-((int) ' ') - 1)

// A built grid: everything encrypt / the solver need, derived once from (keyed_alpha, col_digit).
typedef struct {
    int digit_of_col[TRI_COLS];        // column 0..9 -> its digit label 0..9 (a permutation)
    int col_of_digit[10];              // inverse: digit 0..9 -> column 0..9
    int sep_digit;                     // = digit_of_col[9] (the word-break digit)
    int col_of_letter[TRI_NALPHA];     // letter 0..25 -> its letter-column 0..8
    int digit_of_letter[TRI_NALPHA];   // letter 0..25 -> its ciphertext digit
    int letters_of_col[TRI_LETTER_COLS][TRI_MAXCOLLET]; // column -> its letters (in fill order)
    int ncol[TRI_LETTER_COLS];         // letters in each column (3 for cols 0..7, 2 for col 8)
} TridigitalGrid;

// Build the grid from the 26-letter arrangement `keyed_alpha` (a permutation of 0..25 in
// row-major reading order: row0 cols0..8, row1 cols0..8, row2 cols0..7) and the 10 column
// digit-labels `col_digit` (a permutation of 0..9; col_digit[9] is the separator digit).
// Returns 0 on success, -1 if either input is not a valid permutation.
int tridigital_build_grid(TridigitalGrid *g, const int keyed_alpha[TRI_NALPHA],
                          const int col_digit[TRI_COLS]);

// Convenience builder from the two ACA keywords (letters only, case-insensitive):
//   kw_cols  : 10 DISTINCT letters -> the column digit-labels (alphabetical rank, 10->0).
//   kw_alpha : any keyword -> the 26-letter keyed alphabet (deduped keyword + A..Z tail).
// Returns 0 on success, -1 if kw_cols does not yield 10 distinct letters.
int tridigital_build_from_keywords(TridigitalGrid *g, const char *kw_cols, const char *kw_alpha);

// Encipher plaintext `plain[0..n-1]` into the digit stream out[] (values 0..9). A letter
// (index 0..25) becomes the digit heading its column; ANY other symbol (a space / word-break
// sentinel, punctuation) becomes the separator digit. 1:1 -- returns n (out[] holds n digits).
int tridigital_encrypt(const int plain[], int n, const TridigitalGrid *g, int out[]);

// Candidate plaintext letters for a ciphertext digit: writes the (up to 3) letters sharing
// that digit's column into out_letters[] and returns the count. Returns 0 if `digit` is the
// separator (or out of range) -- i.e. it stands for a word break, not a letter.
int tridigital_candidates(int digit, const TridigitalGrid *g, int out_letters[TRI_MAXCOLLET]);

#endif
