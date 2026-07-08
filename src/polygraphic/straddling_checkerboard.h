#ifndef STRADDLING_CHECKERBOARD_H
#define STRADDLING_CHECKERBOARD_H
#include "colossus.h"

// =====================================================================
//  Straddling Checkerboard cipher primitive (straddling_checkerboard.c)
// =====================================================================
//
// A variable-length fractionation cipher (the classic VIC-cipher building block). A
// 3-row board over 10 columns enciphers plaintext into a DIGIT stream:
//
//   - 10 physical columns carry KEYED digit labels `labels[c]` (a permutation of 0..9,
//     0..9-in-order in the basic cipher, scrambled by a key in the VIC form).
//   - TWO physical columns are BLANK in the top row; their labels are the two ROW
//     INDICATORS r0 < r1. The 8 non-blank top cells encode as ONE digit (`labels[c]`);
//     each of the 20 lower-row cells (rows headed by r0 and r1) encodes as TWO digits
//     (`10*indicator + labels[c]`).
//
//   Board (Wikipedia example, labels 0..9, blanks at columns 2 & 6):
//        0 1 2 3 4 5 6 7 8 9
//          E T   A O N   R I S      top row  (1-digit codes)
//     2  B C D F G H J K L M        row indicator 2 (codes 20..29)
//     6  P Q / U V W X Y Z .        row indicator 6 (codes 60..69)
//   `ATTACK AT DAWN` -> 3113212731223655.
//
// The 28 code-bearing cells hold 27 board SYMBOLS: the 26 letters (0..25) plus a
// FIGURE-SHIFT marker (SC_FIG). One cell is a NULL (unused; a token landing there is
// invalid). Board symbols live in [0, SC_FIG]; the '/' cell of the ACA/Wikipedia board
// is the FIG marker, the '.' cell is the NULL.
//
// FIGURE SHIFT (numeric plaintext). PLAINTEXT symbols are 0..25 (letters) and 26..35
// (the digits 0..9), so digit-bearing text can be enciphered. FIG toggles letter/figure
// mode: in figure mode each token decodes to the DIGIT equal to its cell's column label
// (`labels[c]`), and a second FIG toggles back. Self-inverting; figure-mode boundaries
// need only FIG's cell (part of the board key). Encryption emits a FIG code at each
// letter<->digit transition. The scheme is entirely inside this file.
//
// The board is built by straddling_build_board() from the 27-symbol arrangement `seq`
// (a keyed-alphabet permutation of 0..SC_FIG), the label permutation, and the indicator
// pair; the solver / generator / tests share it, so the convention can never drift.

#define SC_FIG    26     // board-symbol index of the figure-shift marker (letters are 0..25)
#define SC_NSYM   27     // board symbols placed in cells: 26 letters + FIG
#define SC_CELLS  28     // total board cells (8 top + 10 + 10 lower); one is the NULL cell
#define SC_NULL   (-1)   // sym_at_cell value for the single unused (NULL) cell

// A built board: everything encrypt/decrypt need, derived once from (seq, labels, r0, r1).
typedef struct {
    int labels[10];              // column c -> digit label (permutation of 0..9)
    int col_of_label[10];        // inverse: digit -> physical column
    int ind[2];                  // the two row indicators, ind[0] < ind[1]
    int is_ind[10];              // is_ind[d] == 1 iff digit d is a row indicator

    int sym_at_cell[SC_CELLS];   // board symbol at each cell (0..25, SC_FIG, or SC_NULL)
    int cell_of_sym[SC_NSYM];    // inverse for the 27 real symbols (letters + FIG)

    int top_cell_of_col[10];     // cell index of a column's top cell, or -1 if blank
    int lower_cell[2][10];       // cell index of (row ind[0]/ind[1], column c)
    int fig_cell;                // cell holding the FIG marker
    int null_cell;               // the single NULL cell (always the last reading-order cell)
    int digit_cell[10];          // a cell used to emit digit d in figure mode (not FIG/NULL),
                                 //   or -1 if this board cannot encipher digit d

    int clen[SC_CELLS];          // code length of each cell (1 = top, 2 = lower)
    int cd[SC_CELLS][2];         // the code digits of each cell
} StraddlingBoard;

// Build the board from the 27-symbol keyed arrangement `seq` (seq[k] = board symbol at
// reading-order cell k, a permutation of 0..SC_FIG), the label permutation `labels`
// (perm of 0..9), and the two indicator digits r0, r1 (order-independent). Reading order:
// the 8 non-blank top cells left-to-right, then the ind[0] row (cols 0..9), then the
// ind[1] row (cols 0..9); the 28th cell is the NULL. Returns 0 on success, -1 if the
// inputs are not valid permutations / r0 == r1.
int straddling_build_board(StraddlingBoard *b, const int seq[SC_NSYM],
                           const int labels[10], int r0, int r1);

// Encipher plaintext `plain[0..n-1]` (symbols 0..25 letters, 26..35 the digits 0..9)
// into the digit stream out[] (values 0..9). Returns the digit-stream length. out[] must
// hold up to ~2*n + a few figure-shift codes.
int straddling_encrypt(const int plain[], int n, const StraddlingBoard *b, int out[]);

// Decipher the digit stream `digits[0..clen-1]` under board `b`. Writes recovered
// plaintext symbols into out[] (0..25 letters, 26..35 digits) and returns the count.
// *n_tokens / *n_valid (either may be NULL) receive the total / legal token counts (a
// NULL-cell hit or a truncated final token is invalid; a `filler` letter is emitted for
// it). On the true board this returns exactly the original plaintext (all tokens valid).
int straddling_decrypt(const int digits[], int clen, const StraddlingBoard *b,
                       int out[], int filler, int *n_tokens, int *n_valid);

#endif
