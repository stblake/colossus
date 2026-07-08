//
//  Straddling Checkerboard cipher primitives (see straddling_checkerboard.h for the full
//  convention). A keyed 3-row board over 10 columns enciphers plaintext (letters 0..25 and
//  the digits 26..35) into a variable-length DIGIT stream (0..9). High-frequency letters on
//  the top row take one digit; the rest take two (an indicator digit + a column-label
//  digit). A FIGURE-SHIFT board cell toggles letter/figure mode so numeric plaintext can be
//  sent. Everything is derived once by straddling_build_board(); encrypt/decrypt share it
//  with the solver and the generator, so the convention can never drift.
//

#include "straddling_checkerboard.h"

// Validate that p[0..n-1] is a permutation of 0..n-1.
static int sc_is_perm(const int p[], int n) {
    int seen[SC_NSYM];
    if (n > SC_NSYM) return 0;
    for (int i = 0; i < n; i++) seen[i] = 0;
    for (int i = 0; i < n; i++) {
        if (p[i] < 0 || p[i] >= n || seen[p[i]]) return 0;
        seen[p[i]] = 1;
    }
    return 1;
}

int straddling_build_board(StraddlingBoard *b, const int seq[SC_NSYM],
                           const int labels[10], int r0, int r1) {
    if (r0 < 0 || r0 > 9 || r1 < 0 || r1 > 9 || r0 == r1) return -1;
    if (!sc_is_perm(labels, 10) || !sc_is_perm(seq, SC_NSYM)) return -1;

    int lo = (r0 < r1) ? r0 : r1;
    int hi = (r0 < r1) ? r1 : r0;
    b->ind[0] = lo; b->ind[1] = hi;
    for (int d = 0; d < 10; d++) b->is_ind[d] = (d == lo || d == hi);
    for (int c = 0; c < 10; c++) { b->labels[c] = labels[c]; b->col_of_label[labels[c]] = c; }

    int blank0 = b->col_of_label[lo], blank1 = b->col_of_label[hi];

    // Fill cells in reading order: 8 non-blank top cells, then the ind[0] row (cols 0..9),
    // then the ind[1] row (cols 0..9). Cells 0..SC_FIG take seq[]; the last cell is the NULL.
    for (int c = 0; c < 10; c++) b->top_cell_of_col[c] = -1;
    int k = 0;
    for (int c = 0; c < 10; c++) {
        if (c == blank0 || c == blank1) continue;          // blank top column
        b->top_cell_of_col[c] = k;
        b->sym_at_cell[k] = seq[k];
        b->clen[k] = 1; b->cd[k][0] = labels[c];
        k++;
    }
    for (int row = 0; row < 2; row++) {
        for (int c = 0; c < 10; c++) {
            b->lower_cell[row][c] = k;
            b->sym_at_cell[k] = (k < SC_NSYM) ? seq[k] : SC_NULL;   // the last cell is NULL
            b->clen[k] = 2; b->cd[k][0] = b->ind[row]; b->cd[k][1] = labels[c];
            k++;
        }
    }
    b->null_cell = SC_CELLS - 1;                            // cell 27 (row ind[1], col 9)
    b->sym_at_cell[b->null_cell] = SC_NULL;

    for (int s = 0; s < SC_NSYM; s++) b->cell_of_sym[s] = -1;
    for (int cell = 0; cell < SC_CELLS; cell++)
        if (b->sym_at_cell[cell] >= 0) b->cell_of_sym[b->sym_at_cell[cell]] = cell;
    b->fig_cell = b->cell_of_sym[SC_FIG];

    // digit_cell[d]: a cell in column col_of_label[d] (so it decodes to digit d in figure
    // mode) that is neither the FIG nor the NULL cell. Prefer the top cell, then the lower
    // rows. -1 marks the pathological board that cannot encipher digit d.
    for (int d = 0; d < 10; d++) {
        int col = b->col_of_label[d], chosen = -1;
        int cand[3], nc = 0;
        if (b->top_cell_of_col[col] >= 0) cand[nc++] = b->top_cell_of_col[col];
        cand[nc++] = b->lower_cell[0][col];
        cand[nc++] = b->lower_cell[1][col];
        for (int j = 0; j < nc; j++)
            if (cand[j] != b->fig_cell && cand[j] != b->null_cell) { chosen = cand[j]; break; }
        b->digit_cell[d] = chosen;
    }
    return 0;
}

// Append cell `cell`'s code digits to out[o..], returning the new length.
static int sc_emit(const StraddlingBoard *b, int out[], int o, int cell) {
    for (int j = 0; j < b->clen[cell]; j++) out[o++] = b->cd[cell][j];
    return o;
}

int straddling_encrypt(const int plain[], int n, const StraddlingBoard *b, int out[]) {
    int o = 0, mode = 0;                                   // mode 0 = letters, 1 = figures
    for (int i = 0; i < n; i++) {
        int x = plain[i], cell;
        if (x >= 0 && x < 26) {                            // a letter
            if (mode == 1) { o = sc_emit(b, out, o, b->fig_cell); mode = 0; }
            cell = b->cell_of_sym[x];
        } else if (x >= 26 && x < 36) {                    // a digit 0..9
            int d = x - 26;
            if (mode == 0) { o = sc_emit(b, out, o, b->fig_cell); mode = 1; }
            cell = b->digit_cell[d];
        } else {
            continue;                                      // skip anything else defensively
        }
        if (cell < 0) return -1;                           // board can't encipher this symbol
        o = sc_emit(b, out, o, cell);
    }
    return o;
}

int straddling_decrypt(const int digits[], int clen, const StraddlingBoard *b,
                       int out[], int filler, int *n_tokens, int *n_valid) {
    int o = 0, nt = 0, nv = 0, mode = 0, i = 0;
    while (i < clen) {
        int g = digits[i], cell, coldigit;
        if (g < 0 || g > 9) { i++; continue; }             // stray non-digit (shouldn't occur)
        if (b->is_ind[g]) {                                // a 2-digit token
            if (i + 1 >= clen) { out[o++] = filler; nt++; i++; continue; }   // truncated -> invalid
            int g2 = digits[i + 1];
            i += 2;
            if (g2 < 0 || g2 > 9) { out[o++] = filler; nt++; continue; }
            int row = (g == b->ind[1]) ? 1 : 0;
            cell = b->lower_cell[row][b->col_of_label[g2]];
            coldigit = g2;
        } else {                                           // a 1-digit token
            i += 1;
            cell = b->top_cell_of_col[b->col_of_label[g]];
            coldigit = g;
        }
        int s = (cell >= 0) ? b->sym_at_cell[cell] : SC_NULL;
        if (s == SC_NULL) { out[o++] = filler; nt++; continue; }             // NULL cell -> invalid
        if (s == SC_FIG)  { mode ^= 1; nt++; nv++; continue; }               // figure-shift toggle
        nt++; nv++;
        out[o++] = (mode == 0) ? s : (26 + coldigit);      // letter, or figure-mode digit
    }
    if (n_tokens) *n_tokens = nt;
    if (n_valid)  *n_valid  = nv;
    return o;
}
