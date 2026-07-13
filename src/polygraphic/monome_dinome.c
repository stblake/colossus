//
//  Monome-Dinome cipher primitives (see monome_dinome.h for the full convention). A keyed
//  3x8 box over a 24-letter alphabet (J->I, Z->Y) enciphers plaintext into a variable-length
//  DIGIT stream (0..9). Top-row letters take one digit (the column label); 2nd/3rd-row
//  letters take two (a row-indicator digit + a column-label digit). Everything is derived
//  once by monome_dinome_build_board(); encrypt/decrypt share it with the solver and the
//  generator, so the convention can never drift.
//

#include "monome_dinome.h"

// Validate that p[0..n-1] is a permutation of 0..n-1.
static int md_is_perm(const int p[], int n) {
    int seen[MD_NALPHA];
    if (n > MD_NALPHA) return 0;
    for (int i = 0; i < n; i++) seen[i] = 0;
    for (int i = 0; i < n; i++) {
        if (p[i] < 0 || p[i] >= n || seen[p[i]]) return 0;
        seen[p[i]] = 1;
    }
    return 1;
}

int monome_dinome_build_board(MonomeDinomeBoard *b, const int letters[MD_NALPHA],
                              const int col_label[MD_COLS], int r0, int r1) {
    if (r0 < 0 || r0 > 9 || r1 < 0 || r1 > 9 || r0 == r1) return -1;
    if (!md_is_perm(letters, MD_NALPHA)) return -1;

    // The 8 column labels must be distinct digits, none equal to a row indicator.
    int seen[10]; for (int d = 0; d < 10; d++) seen[d] = 0;
    seen[r0] = seen[r1] = 1;
    for (int c = 0; c < MD_COLS; c++) {
        int d = col_label[c];
        if (d < 0 || d > 9 || seen[d]) return -1;   // out of range, an indicator, or repeated
        seen[d] = 1;
    }

    b->ind[0] = r0; b->ind[1] = r1;
    for (int d = 0; d < 10; d++) { b->is_ind[d] = (d == r0 || d == r1); b->col_of_label[d] = -1; }
    for (int c = 0; c < MD_COLS; c++) { b->col_label[c] = col_label[c]; b->col_of_label[col_label[c]] = c; }

    // Fill cells in reading order: top row (cols 0..7, 1-digit), then the r0 row and the r1
    // row (cols 0..7, 2-digit each). Cell k holds letter letters[k].
    int k = 0;
    for (int c = 0; c < MD_COLS; c++) {                 // top row: monome
        b->letter_at_cell[k] = letters[k];
        b->clen[k] = 1; b->cd[k][0] = col_label[c];
        k++;
    }
    for (int row = 0; row < 2; row++) {                 // the two indicator rows: dinome
        for (int c = 0; c < MD_COLS; c++) {
            b->letter_at_cell[k] = letters[k];
            b->clen[k] = 2; b->cd[k][0] = b->ind[row]; b->cd[k][1] = col_label[c];
            k++;
        }
    }
    for (int s = 0; s < MD_NALPHA; s++) b->cell_of_letter[s] = -1;
    for (int cell = 0; cell < MD_CELLS; cell++) b->cell_of_letter[b->letter_at_cell[cell]] = cell;
    return 0;
}

int monome_dinome_encrypt(const int plain[], int n, const MonomeDinomeBoard *b, int out[]) {
    int o = 0;
    for (int i = 0; i < n; i++) {
        int x = plain[i];
        if (x < 0 || x >= MD_NALPHA) continue;          // skip anything outside the 24-letter box
        int cell = b->cell_of_letter[x];
        for (int j = 0; j < b->clen[cell]; j++) out[o++] = b->cd[cell][j];
    }
    return o;
}

int monome_dinome_decrypt(const int digits[], int clen, const MonomeDinomeBoard *b,
                          int out[], int filler, int *n_tokens, int *n_valid) {
    int o = 0, nt = 0, nv = 0, i = 0;
    while (i < clen) {
        int g = digits[i], cell;
        if (g < 0 || g > 9) { i++; continue; }          // stray non-digit (shouldn't occur)
        if (b->is_ind[g]) {                             // a 2-digit token (dinome)
            if (i + 1 >= clen) { out[o++] = filler; nt++; i++; continue; }   // truncated -> invalid
            int g2 = digits[i + 1];
            i += 2;
            if (g2 < 0 || g2 > 9 || b->col_of_label[g2] < 0) {              // 2nd digit not a label
                out[o++] = filler; nt++; continue;                         //   (indicator) -> invalid
            }
            int row = (g == b->ind[1]) ? 1 : 0;
            cell = MD_COLS + row * MD_COLS + b->col_of_label[g2];           // rows 1,2 start at 8,16
        } else {                                        // a 1-digit token (monome)
            i += 1;
            cell = b->col_of_label[g];                                     // top row, always a label
        }
        nt++; nv++;
        out[o++] = b->letter_at_cell[cell];
    }
    if (n_tokens) *n_tokens = nt;
    if (n_valid)  *n_valid  = nv;
    return o;
}
