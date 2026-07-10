//
//  Digrafid cipher primitives (digraphic fractionation over two keyed 27-symbol alphabets).
//
//  The tableau is TWO mixed alphabets over the 27 symbols A..Z + '#':
//    - the HORIZONTAL grid H, 3 rows x 9 cols (cell = hr*9 + hc, hr 0..2, hc 0..8), and
//    - the VERTICAL   grid V, 9 rows x 3 cols (cell = vr*3 + vc, vr 0..8, vc 0..2),
//  each a permutation of 0..26. gridH[p] / gridV[p] is the symbol at cell p.
//
//  A plaintext DIGRAPH (a, b) maps to a 3-digit number (top, mid, bot), each digit 0..8
//  (the ACA tableau labels them 1..9):
//    top = a's column in H (hc),  bot = b's row in V (vr),  mid = (a's row in H)*3 + (b's
//    col in V) = hr*3 + vc.  This is a bijection between the 729 digraphs and 729 triples.
//  The inverse (triple -> digraph): hr = mid/3, vc = mid%3, a = gridH[hr*9 + top],
//  b = gridV[bot*3 + vc].
//
//  Fractionation: the plaintext is taken in digraphs grouped into blocks of `period`
//  digraphs.  Within a block of g digraphs the g triples are stacked as 3 rows
//  (tops / mids / bots), read row-major into a 3g-digit stream, re-split into g consecutive
//  triples, and each new triple is mapped back through the tableau to one ciphertext
//  digraph -- exactly the Trifid reshape, but over digraphs.  Decryption reverses the
//  reshape.  A short final block (g < period) is handled in place; a lone trailing letter
//  (odd length) passes through unchanged.  Period 1 is the identity (a single-digraph
//  block's reshape is a no-op, and digraph -> triple -> same digraph).
//
//  The solver needs only digrafid_decrypt(); encrypt / grid_from_keyword serve the
//  test-data generator and the unit tests.
//

#include "colossus.h"

// Single-threaded scratch for the per-block digit stream (3 digits per digraph, at most
// 3*period <= 3*(len/2) entries). Kept off the stack like bifid's stream scratch.
static _Thread_local int g_digrafid_stream[3 * MAX_CIPHER_LENGTH];

// Plaintext/ciphertext digraph (a, b) -> 3-digit number (top, mid, bot), each 0..8, using
// the precomputed inverse tables posH (symbol -> H cell) and posV (symbol -> V cell).
static inline void digrafid_digraph_to_triple(int a, int b, const int posH[], const int posV[],
                                              int *top, int *mid, int *bot) {
    int cellH = posH[a];
    int hr = cellH / DIGRAFID_HCOLS, hc = cellH % DIGRAFID_HCOLS;   // 0..2, 0..8
    int cellV = posV[b];
    int vr = cellV / DIGRAFID_VCOLS, vc = cellV % DIGRAFID_VCOLS;   // 0..8, 0..2
    *top = hc;
    *bot = vr;
    *mid = hr * DIGRAFID_VCOLS + vc;                                // hr*3 + vc, 0..8
}

// 3-digit number (top, mid, bot) -> digraph (a, b) via the grids.
static inline void digrafid_triple_to_digraph(int top, int mid, int bot,
                                              const int gridH[], const int gridV[],
                                              int *a, int *b) {
    int hr = mid / DIGRAFID_VCOLS, vc = mid % DIGRAFID_VCOLS;       // 0..2, 0..2
    *a = gridH[hr * DIGRAFID_HCOLS + top];
    *b = gridV[bot * DIGRAFID_VCOLS + vc];
}

// Encipher plaintext (any length) under the two keyed grids, in blocks of `period`
// digraphs. A lone trailing letter (odd len) passes through unchanged.
void digrafid_encrypt(const int plain[], int len, const int gridH[], const int gridV[],
                      int period, int out[]) {
    int posH[DIGRAFID_GRID], posV[DIGRAFID_GRID];
    bifid_build_inverse(gridH, posH, DIGRAFID_GRID);
    bifid_build_inverse(gridV, posV, DIGRAFID_GRID);
    int *s = g_digrafid_stream;
    int ndig = len / 2;                                            // complete digraphs
    if (period < 1) period = 1;

    for (int doff = 0; doff < ndig; doff += period) {
        int g = (doff + period <= ndig) ? period : (ndig - doff);
        // digraph -> triple, written row-major: tops | mids | bots.
        for (int j = 0; j < g; j++) {
            int a = plain[2 * (doff + j)], b = plain[2 * (doff + j) + 1];
            int top, mid, bot;
            digrafid_digraph_to_triple(a, b, posH, posV, &top, &mid, &bot);
            s[0 * g + j] = top;
            s[1 * g + j] = mid;
            s[2 * g + j] = bot;
        }
        // re-split into consecutive triples -> output digraphs.
        for (int j = 0; j < g; j++) {
            int a, b;
            digrafid_triple_to_digraph(s[3 * j], s[3 * j + 1], s[3 * j + 2], gridH, gridV, &a, &b);
            out[2 * (doff + j)]     = a;
            out[2 * (doff + j) + 1] = b;
        }
    }
    if (len & 1) out[len - 1] = plain[len - 1];                    // lone trailing letter
}

// Decipher ciphertext (any length) into out[]. Inverse of digrafid_encrypt.
void digrafid_decrypt(const int cipher[], int len, const int gridH[], const int gridV[],
                      int period, int out[]) {
    int posH[DIGRAFID_GRID], posV[DIGRAFID_GRID];
    bifid_build_inverse(gridH, posH, DIGRAFID_GRID);
    bifid_build_inverse(gridV, posV, DIGRAFID_GRID);
    int *s = g_digrafid_stream;
    int ndig = len / 2;
    if (period < 1) period = 1;

    for (int doff = 0; doff < ndig; doff += period) {
        int g = (doff + period <= ndig) ? period : (ndig - doff);
        // ciphertext digraph -> its (new) triple, written as consecutive triples.
        for (int j = 0; j < g; j++) {
            int a = cipher[2 * (doff + j)], b = cipher[2 * (doff + j) + 1];
            int top, mid, bot;
            digrafid_digraph_to_triple(a, b, posH, posV, &top, &mid, &bot);
            s[3 * j]     = top;
            s[3 * j + 1] = mid;
            s[3 * j + 2] = bot;
        }
        // read the stream row-major to recover the original per-digraph triples.
        for (int c = 0; c < g; c++) {
            int a, b;
            digrafid_triple_to_digraph(s[0 * g + c], s[1 * g + c], s[2 * g + c], gridH, gridV, &a, &b);
            out[2 * (doff + c)]     = a;
            out[2 * (doff + c) + 1] = b;
        }
    }
    if (len & 1) out[len - 1] = cipher[len - 1];                   // lone trailing letter
}

// Build a grid from a keyword (alphabet indices): the keyword letters in order with
// duplicates removed, then the remaining alphabet letters ascending (the keyed sequence),
// placed into a rows x cols grid (n = rows*cols cells, stored row-major). column_major == 0
// fills row by row (the horizontal alphabet); != 0 fills column by column (the vertical
// alphabet, e.g. keyword VERTICAL read down the first column of the 9x3 grid).
void digrafid_grid_from_keyword(const int keyword[], int kwlen, int grid[],
                                int rows, int cols, int column_major) {
    int n = rows * cols;
    int seq[DIGRAFID_GRID];
    char used[DIGRAFID_GRID];
    for (int l = 0; l < n; l++) used[l] = 0;
    int m = 0;
    for (int i = 0; i < kwlen && m < n; i++) {
        int l = keyword[i];
        if (l < 0 || l >= n || used[l]) continue;
        used[l] = 1;
        seq[m++] = l;
    }
    for (int l = 0; l < n && m < n; l++)
        if (!used[l]) { used[l] = 1; seq[m++] = l; }

    if (!column_major) {
        for (int k = 0; k < n; k++) grid[k] = seq[k];              // row-major
    } else {
        for (int k = 0; k < n; k++) grid[(k % rows) * cols + (k / rows)] = seq[k];
    }
}
