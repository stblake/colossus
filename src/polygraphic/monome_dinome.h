#ifndef MONOME_DINOME_H
#define MONOME_DINOME_H
#include "colossus.h"

// =====================================================================
//  Monome-Dinome cipher primitive (monome_dinome.c)
// =====================================================================
//
// A variable-length fractionation cipher (ACA "Monome-Dinome", 60-120 plaintext
// letters). A 3-row by 8-column keyed box enciphers plaintext into a DIGIT stream:
//
//   - The box holds a 24-LETTER keyed alphabet (keyword deduped, then the unused
//     letters in order). Two letters SHARE cells so 26 -> 24: J is written as I and
//     Z as Y (the ACA convention). So the plaintext alphabet is A..Z minus {J,Z}.
//   - The 10 digits 0..9 split into TWO ROW INDICATORS (heading the 2nd and 3rd rows)
//     and EIGHT COLUMN LABELS (over the 8 columns). The top row has NO indicator.
//   - A top-row letter enciphers as ONE digit  (its column label)          -- a "monome".
//     A 2nd/3rd-row letter enciphers as TWO digits (row indicator + column label)
//                                                                          -- a "dinome".
//
//   Board (ACA worked example, keyword RMASTERTON -> 6318927054: rows 6,3; cols 1 8 9 2 7 0 5 4):
//        1 8 9 2 7 0 5 4        column labels
//        N O T A R I E S        top row    (1-digit / monome codes)
//     6  B C D F G H K L        row indicator 6 (codes 6x)
//     3  M P Q U V W X Y        row indicator 3 (codes 3x)
//   HIGHFREQUENCY... -> 60 0 67 60 62 7 5 39 32 5 1 68 34 65 ...
//
// Unlike the Straddling Checkerboard, EVERY cell holds a letter (24 cells, 24 letters):
// there is no NULL cell and no figure-shift -- the plaintext is letters only. The two
// row-indicator digits never appear as a column label, so a dinome's SECOND digit is
// always a column-label digit; a decode whose second digit is an indicator is invalid
// (a strong signal a wrong indicator pair was assumed -- the solver's validity reward).
//
// The box is built by monome_dinome_build_board() from the 24-letter arrangement
// `letters` (reading order: top row cols 0..7, indicator-r0 row cols 0..7, indicator-r1
// row cols 0..7), the 8 column-label digits, and the ordered indicator pair (r0, r1);
// the solver / generator / tests share it, so the convention can never drift.

#define MD_NALPHA  24    // letters in the keyed box (A..Z with J->I, Z->Y)
#define MD_COLS    8     // columns
#define MD_CELLS   24    // total box cells (8 top + 8 + 8), one letter each

// A built box: everything encrypt/decrypt need, derived once from (letters, col_label, r0, r1).
typedef struct {
    int letter_at_cell[MD_CELLS];   // cell (reading order) -> letter 0..23
    int cell_of_letter[MD_NALPHA];  // inverse: letter 0..23 -> cell

    int ind[2];                     // row indicators: ind[0] heads row 1, ind[1] heads row 2
    int col_label[MD_COLS];         // column c (0..7) -> its digit label
    int col_of_label[10];           // digit -> column (0..7), or -1 if the digit is an indicator
    int is_ind[10];                 // is_ind[d] == 1 iff digit d is a row indicator

    int clen[MD_CELLS];             // code length of each cell (1 = top, 2 = lower)
    int cd[MD_CELLS][2];            // the code digits of each cell
} MonomeDinomeBoard;

// Build the box from the 24-letter arrangement `letters` (perm of 0..23, reading order:
// top row cols 0..7, then the r0 row, then the r1 row), the 8 column-label digits
// `col_label` (8 distinct digits in 0..9, none equal to r0/r1), and the two indicator
// digits r0, r1 (distinct, not among the column labels). Returns 0 on success, -1 if the
// inputs are inconsistent (bad permutation, r0 == r1, overlap of labels and indicators).
int monome_dinome_build_board(MonomeDinomeBoard *b, const int letters[MD_NALPHA],
                              const int col_label[MD_COLS], int r0, int r1);

// Encipher plaintext letters `plain[0..n-1]` (each 0..23) into the digit stream out[]
// (values 0..9). Returns the digit-stream length. out[] must hold up to 2*n digits.
// A symbol outside 0..23 is skipped defensively.
int monome_dinome_encrypt(const int plain[], int n, const MonomeDinomeBoard *b, int out[]);

// Decipher the digit stream `digits[0..clen-1]` under box `b`. Writes recovered plaintext
// letters into out[] (0..23) and returns the count. *n_tokens / *n_valid (either may be
// NULL) receive the total / legal token counts (a truncated final dinome or a dinome whose
// second digit is a row indicator is invalid; a `filler` letter is emitted for it). On the
// true box this returns exactly the original plaintext with every token valid.
int monome_dinome_decrypt(const int digits[], int clen, const MonomeDinomeBoard *b,
                          int out[], int filler, int *n_tokens, int *n_valid);

#endif
