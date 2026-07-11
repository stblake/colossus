
// Transpositions

#include "colossus.h"

void transperoffset(int plaintext[], int len, int d, int n) {

    if (d == 1 && n == 0) return; // Identity transformation.

    int indx, temp[MAX_CIPHER_LENGTH];
    
    // Periodic decimation.
    for (int i = 0; i < len; i++) {
        temp[i] = plaintext[(d * i) % len];
    }

    // Rotation (Offset.)
    for (int i = 0; i < len; i++) {
    	indx = (i + n) % len;
    	if (indx < 0) {
    		indx += len;
    	}
        plaintext[i] = temp[indx];
    }
    return ;
}



void matrix_rotate(int text[], int len, int width, int clockwise) {
    if (width <= 1 || width >= len) return; // Identity or 1D matrix

    int R = (len + width - 1) / width; // Ceiling division for rows
    int W = width;
    int temp[MAX_CIPHER_LENGTH];
    int idx = 0;

    if (clockwise) {
        // Read columns left-to-right, but from bottom row to top row
        for (int c = 0; c < W; c++) {
            for (int r = R - 1; r >= 0; r--) {
                int old_idx = r * W + c;
                // Only read if the cell is valid (handles incomplete final rows)
                if (old_idx < len) {
                    temp[idx++] = text[old_idx];
                }
            }
        }
    } else {
        // Anti-clockwise: Read columns right-to-left, top to bottom
        for (int c = W - 1; c >= 0; c--) {
            for (int r = 0; r < R; r++) {
                int old_idx = r * W + c;
                if (old_idx < len) {
                    temp[idx++] = text[old_idx];
                }
            }
        }
    }

    // Copy back to original array
    for (int i = 0; i < len; i++) {
        text[i] = temp[i];
    }
}

void transmatrix(int text[], int len, int w1, int w2, int clockwise) {
    // Perform a K3-like double rotation.
    matrix_rotate(text, len, w1, clockwise);
    matrix_rotate(text, len, w2, clockwise);
}


// Invert one columnar transposition stage.
//
// Encryption writes the plaintext into a grid of K columns, row by row
// left-to-right, then reads the columns off in order `order[0..K-1]` (each column
// top-to-bottom for COL_READ_TB, bottom-to-top for COL_READ_BT). The grid's last
// row is short when len % K != 0: the leftmost (len % K) columns are one cell
// taller than the rest, so each column's height is known up front.
//
// To decrypt we slice the ciphertext back into those columns in read order,
// refill the grid, and read it row-major. out[] must not alias cipher[].
void decrypt_columnar(int cipher[], int len, int K, int order[], int dir, int out[]) {

    if (K <= 1 || K > len) {            // degenerate: a single column is the identity
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int grid[MAX_CIPHER_LENGTH];
    int R = (len + K - 1) / K;          // number of rows (ceiling)
    int rem = len % K;                  // tall columns are 0..rem-1 (all K if rem==0)

    // Refill the grid one column at a time, consuming the ciphertext in read order.
    int pos = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        int h = (rem == 0 || c < rem) ? R : R - 1;   // height of grid column c
        if (dir == COL_READ_BT) {
            for (int r = h - 1; r >= 0; r--) grid[r * K + c] = cipher[pos++];
        } else { // COL_READ_TB
            for (int r = 0; r < h; r++) grid[r * K + c] = cipher[pos++];
        }
    }

    // Read the grid row-major to recover the plaintext (skip missing short-row cells).
    int o = 0;
    for (int r = 0; r < R; r++) {
        for (int c = 0; c < K; c++) {
            int h = (rem == 0 || c < rem) ? R : R - 1;
            if (r < h) out[o++] = grid[r * K + c];
        }
    }
}


// Columnar transposition with a uniform within-column row permutation L (the
// "track order"). The grid is refilled from the K columns exactly as
// decrypt_columnar (column read per `dir`), but the recovered plaintext reads the
// R grid rows in the order L[0..R-1] -- output row i is grid row L[i] -- each row
// left-to-right (row_rev == 0) or right-to-left (row_rev == 1). Because a single
// permutation is applied inside every column, it is equivalent to permuting whole
// grid rows, so the tracked read needs a COMPLETE grid (len % K == 0, R = len/K);
// a ragged grid falls back to plain decrypt_columnar (L/row_rev ignored). With
// L == identity and row_rev == 0 this is exactly decrypt_columnar on a complete
// grid. out[] must not alias cipher[].
void decrypt_columnar_tracked(int cipher[], int len, int K, int order[], int dir,
                              int L[], int row_rev, int out[]) {
    if (K <= 1 || K > len || len % K != 0) {        // tracked read needs a complete grid
        decrypt_columnar(cipher, len, K, order, dir, out);
        return;
    }
    int R = len / K;
    int grid[MAX_CIPHER_LENGTH];
    int pos = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        if (dir == COL_READ_BT) for (int r = R - 1; r >= 0; r--) grid[r * K + c] = cipher[pos++];
        else                    for (int r = 0; r < R; r++)     grid[r * K + c] = cipher[pos++];
    }
    int o = 0;
    for (int i = 0; i < R; i++) {
        int r = L[i];
        if (row_rev) for (int c = K - 1; c >= 0; c--) out[o++] = grid[r * K + c];
        else         for (int c = 0; c < K; c++)     out[o++] = grid[r * K + c];
    }
}


// Sub-grid / tile transposition. The ciphertext is W columns each R tall, laid out
// contiguously (block c = cipher[c*R .. c*R+R-1]); the global key `order` permutes
// the columns into an R x W grid (G[r][c] = block order[c], row r). Within that grid
// every complete h x w tile (tiled from the origin) had the SAME cell permutation
// `perm` (length h*w, cell s at tile offset (s/w, s%w)) applied; leftover rows
// (R % h) and columns (W % w) pass through untouched. Decryption undoes the tile
// permutation (grid cell at tile offset s takes the value at offset perm[s]) and
// reads the grid row-major. Complete grid only (len % W == 0, R = len/W); a ragged
// or degenerate W is the identity. out[] must not alias cipher[].
void decrypt_tile(int cipher[], int len, int W, int order[], int h, int w,
                  int perm[], int out[]) {
    if (W <= 1 || W > len || len % W != 0 || h < 1 || w < 1) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }
    int R = len / W;
    int G[MAX_CIPHER_LENGTH], P[MAX_CIPHER_LENGTH];
    for (int c = 0; c < W; c++)
        for (int r = 0; r < R; r++) G[r * W + c] = cipher[order[c] * R + r];
    for (int i = 0; i < len; i++) P[i] = G[i];           // leftover cells pass through
    int m = h * w;
    for (int ti = 0; ti < R / h; ti++)
        for (int tj = 0; tj < W / w; tj++) {
            int r0 = ti * h, c0 = tj * w;
            for (int s = 0; s < m; s++) {
                int dr = s / w, dc = s % w;
                int sr = perm[s] / w, sc = perm[s] % w;
                P[(r0 + dr) * W + (c0 + dc)] = G[(r0 + sr) * W + (c0 + sc)];
            }
        }
    for (int i = 0; i < len; i++) out[i] = P[i];
}


// Apply a ciphertext->plaintext position map. pt_of_ct[k] is the plaintext
// position that received ciphertext position k, so standard decryption scatters
// the ciphertext (out[pt_of_ct[k]] = cipher[k]) while -variant (encryption swapped
// for decryption) gathers it (out[k] = cipher[pt_of_ct[k]]). Shared by the
// permutation-style transposition primitives below. out[] must not alias cipher[].
static void apply_perm(const int cipher[], const int pt_of_ct[], int len, int variant, int out[]) {
    if (variant) for (int k = 0; k < len; k++) out[k] = cipher[pt_of_ct[k]];
    else         for (int k = 0; k < len; k++) out[pt_of_ct[k]] = cipher[k];
}


// =====================================================================
//  Rail fence (and Variant Rail fence)
// =====================================================================
//
// Encryption writes the plaintext along a zigzag that runs down then up across
// `rails` rows -- position i lands on rail r(i) where the rail index sweeps
// 0,1,...,rails-1,rails-2,...,1,0,1,... with period P = 2*(rails-1). The
// ciphertext is then the rows read off top to bottom, each left to right.
//
// `offset` shifts the starting phase of the zigzag: a standard rail fence has
// offset == 0 (begins on the top rail going down), while the ACA "Variant Rail
// fence" begins partway through the sweep. Sweeping offset over [0, P-1] lets the
// one primitive cover both, so a single solver attacks standard and variant rails.
//
// To decrypt we recompute r(i) for every position, then walk the ciphertext rail
// by rail (the order it was written out) handing each cipher letter back to its
// original position. out[] must not alias cipher[].
//
// `variant` swaps encryption for decryption (the global -variant convention): the
// cipher was produced by reading the plaintext along the rails and laying it back
// out along the zigzag -- the inverse layout -- so recovery applies the forward map
// (out[pos] = cipher[i]) instead of the inverse (out[i] = cipher[pos]).
void decrypt_railfence(int cipher[], int len, int rails, int offset, int variant, int out[]) {

    if (rails <= 1 || rails >= len) {   // degenerate: a single rail is the identity
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int P = 2 * (rails - 1);            // length of one down-and-up sweep
    int pos = 0;                        // cursor into the rail-by-rail read order

    // The ciphertext is rail 0's letters, then rail 1's, ... Walk the rails in
    // that order and, within each, the positions left to right, to invert it.
    for (int r = 0; r < rails; r++) {
        for (int i = 0; i < len; i++) {
            int ph = (i + offset) % P;
            int rail = (ph < rails) ? ph : P - ph;
            if (rail != r) continue;
            if (variant) out[pos++] = cipher[i];   // apply the forward permutation
            else         out[i] = cipher[pos++];   // apply the inverse (standard)
        }
    }
}


// =====================================================================
//  Route transposition
// =====================================================================
//
// Encryption writes the plaintext row-major into an R x C grid, then reads the
// cells out along a geometric "route". The grid may be ragged: with a short final
// row ((R-1)*C < len < R*C) the cells whose row-major index r*C + c is >= len do
// not exist, so `route_cells` skips them and walks the route over the surviving
// cells only. `route_cells` returns that route as the reading order: cells[k] is
// the row-major index (r*C + c) of the k-th cell emitted, and the return value is
// the number of cells emitted (== len for a valid grid). N_ROUTES routes are defined:
//
//   0  rows, boustrophedon (snake) -- row 0 left->right, row 1 right->left, ...
//   1  columns, boustrophedon      -- col 0 top->bottom, col 1 bottom->top, ...
//   2  spiral inward, clockwise from the top-left corner (start heading right)
//   3  spiral inward, counter-clockwise from the top-left corner (start heading down)
//   4  anti-diagonals (r+c constant), boustrophedon (alternate scan direction)
//   5  anti-diagonals (r+c constant), constant direction (always r ascending)
//
// Plain row-major and plain column reads are already reachable through the
// columnar solver, so the route set focuses on the snake, spiral and diagonal
// patterns that the dedicated columnar primitive cannot express.
int route_cells(int R, int C, int len, int route_id, int cells[]) {
    int k = 0;
    // EMIT records cell (idx) only if it exists in the (possibly ragged) grid, i.e.
    // its row-major index is within the len letters actually written. For a complete
    // grid (R*C == len) the guard is always true and every cell is emitted.
    #define EMIT(idx) do { int _i = (idx); if (_i >= 0 && _i < len) cells[k++] = _i; } while (0)
    switch (route_id) {
        case 0: // rows, snake
            for (int r = 0; r < R; r++) {
                if ((r & 1) == 0) for (int c = 0;     c < C; c++) EMIT(r * C + c);
                else              for (int c = C - 1; c >= 0; c--) EMIT(r * C + c);
            }
            break;
        case 1: // columns, snake
            for (int c = 0; c < C; c++) {
                if ((c & 1) == 0) for (int r = 0;     r < R; r++) EMIT(r * C + c);
                else              for (int r = R - 1; r >= 0; r--) EMIT(r * C + c);
            }
            break;
        case 2: { // clockwise spiral inward from top-left
            int top = 0, bottom = R - 1, left = 0, right = C - 1;
            while (top <= bottom && left <= right) {
                for (int c = left; c <= right; c++)  EMIT(top * C + c);
                top++;
                for (int r = top; r <= bottom; r++)  EMIT(r * C + right);
                right--;
                if (top <= bottom)
                    for (int c = right; c >= left; c--) EMIT(bottom * C + c);
                bottom--;
                if (left <= right)
                    for (int r = bottom; r >= top; r--) EMIT(r * C + left);
                left++;
            }
            break;
        }
        case 3: { // counter-clockwise spiral inward from top-left
            int top = 0, bottom = R - 1, left = 0, right = C - 1;
            while (top <= bottom && left <= right) {
                for (int r = top; r <= bottom; r++)  EMIT(r * C + left);
                left++;
                for (int c = left; c <= right; c++)  EMIT(bottom * C + c);
                bottom--;
                if (left <= right)
                    for (int r = bottom; r >= top; r--) EMIT(r * C + right);
                right--;
                if (top <= bottom)
                    for (int c = right; c >= left; c--) EMIT(top * C + c);
                top++;
            }
            break;
        }
        case 4:   // anti-diagonals, boustrophedon
        case 5: { // anti-diagonals, constant direction
            for (int d = 0; d <= R + C - 2; d++) {
                int r_lo = (d - (C - 1) > 0) ? d - (C - 1) : 0;  // c = d - r <= C-1
                int r_hi = (d < R - 1) ? d : R - 1;              // r <= R-1
                if (route_id == 4 && (d & 1)) {                  // reverse on odd diagonals
                    for (int r = r_hi; r >= r_lo; r--) EMIT(r * C + (d - r));
                } else {
                    for (int r = r_lo; r <= r_hi; r++) EMIT(r * C + (d - r));
                }
            }
            break;
        }
        default: // unknown route -> identity (row-major)
            for (int i = 0; i < len; i++) cells[k++] = i;
            break;
    }
    #undef EMIT
    return k;
}

// Invert one route transposition. The plaintext was written row-major, so the
// recovered plaintext is just the grid read row-major; cell cells[k] received the
// k-th ciphertext letter, hence out[cells[k]] = cipher[k]. out[] must not alias
// cipher[]. The grid may be ragged (a short final row): a valid grid holds exactly
// len cells, i.e. (R-1)*C < len <= R*C; anything else is rejected as the identity so
// the caller can skip that shape.
//
// `variant` swaps encryption for decryption (the global -variant convention): the
// cipher was produced by writing along the route and reading row-major, so recovery
// applies the forward map (out[k] = cipher[cells[k]]) instead of the inverse.
void decrypt_route(int cipher[], int len, int R, int C, int route_id, int variant, int out[]) {
    if (R < 2 || C < 2 || (long)(R - 1) * C >= len || (long)R * C < len) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }
    int cells[MAX_CIPHER_LENGTH];
    int n = route_cells(R, C, len, route_id, cells);
    if (n != len) {   // route could not cover the grid exactly: leave unchanged
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }
    if (variant) for (int k = 0; k < len; k++) out[k] = cipher[cells[k]];
    else         for (int k = 0; k < len; k++) out[cells[k]] = cipher[k];
}


// =====================================================================
//  Amsco
// =====================================================================
//
// Encryption fills a grid of K columns row by row, each cell taking alternately
// one then two letters (`start` = the size of the very first cell, 1 or 2), with
// the 1/2 alternation continuing without break from one row into the next. The
// columns are then read off in key order `order[0..K-1]`, each top to bottom,
// concatenating each cell's one-or-two letters.
//
// Because the alternation runs continuously by cell index, the size of cell m
// (row m/K, column m%K, filled row-major) is fixed up front, so we can rebuild the
// exact per-cell letter layout, slice the ciphertext back into columns in key
// order, and scatter each cell's letters to their original plaintext offsets.
// `variant` swaps read/write. out[] must not alias cipher[].
void decrypt_amsco(int cipher[], int len, int K, int order[], int start, int variant, int out[]) {

    if (K <= 1 || K > len) {            // degenerate: a single column is the identity
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    // 1. Walk the cells in row-major fill order, recording each cell's letter
    //    count and the plaintext offset where its letters begin, until the text
    //    is exhausted (the final cell may be clipped to the letters that remain).
    int cell_size[MAX_CIPHER_LENGTH];  // 1 or 2 (last cell possibly fewer)
    int cell_off[MAX_CIPHER_LENGTH];   // plaintext offset of the cell's first letter
    int sz_even = start;               // size of even-indexed cells
    int sz_odd  = (start == 1) ? 2 : 1;
    int n_cells = 0, placed = 0;
    while (placed < len) {
        int nominal = ((n_cells & 1) == 0) ? sz_even : sz_odd;
        int sz = (len - placed < nominal) ? (len - placed) : nominal;
        cell_off[n_cells]  = placed;
        cell_size[n_cells] = sz;
        placed += sz;
        n_cells++;
    }

    // 2. Read the columns in key order, each top-to-bottom (cells c, c+K, c+2K,
    //    ...), emitting each cell's letters to build the ciphertext->plaintext map.
    int pt_of_ct[MAX_CIPHER_LENGTH];
    int k = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        for (int m = c; m < n_cells; m += K)
            for (int i = 0; i < cell_size[m]; i++)
                pt_of_ct[k++] = cell_off[m] + i;
    }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Myszkowski
// =====================================================================
//
// A columnar transposition whose keyword may contain repeated letters. `rank[c]`
// is the numeric rank of column c (equal ranks for equal keyword letters). The
// grid is filled row-major into K columns (incomplete final row when len % K != 0).
// Columns are taken off by ascending rank; a rank shared by several columns is
// read **row by row across those columns together** (left to right) rather than
// column by column -- the feature that distinguishes Myszkowski from plain columnar
// (to which it reduces when all ranks are distinct).
//
// `variant` swaps read/write. out[] must not alias cipher[]. Grid cell (r,c) holds
// plaintext position r*K + c and exists iff that position is < len.
void decrypt_myszkowski(int cipher[], int len, int K, int rank[], int variant, int out[]) {

    if (K <= 1 || K > len) {            // degenerate: a single column is the identity
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int R = (len + K - 1) / K;          // number of rows (ceiling)
    int pt_of_ct[MAX_CIPHER_LENGTH];
    int k = 0;

    // Process ranks in ascending order. To avoid sorting, sweep candidate rank
    // values and, for each, gather the columns carrying it (left to right).
    bool done[MAX_COLS] = { false };
    int processed = 0;
    while (processed < K) {
        // Find the smallest not-yet-processed rank value.
        int v = 0; bool have_v = false;
        for (int c = 0; c < K; c++) {
            if (done[c]) continue;
            if (!have_v || rank[c] < v) { v = rank[c]; have_v = true; }
        }

        // Columns sharing rank v, in left-to-right order.
        int group[MAX_COLS], g = 0;
        for (int c = 0; c < K; c++) if (!done[c] && rank[c] == v) { group[g++] = c; done[c] = true; }
        processed += g;

        if (g == 1) {
            // Single column: read it top-to-bottom.
            int c = group[0];
            for (int r = 0; r < R; r++) {
                int pos = r * K + c;
                if (pos < len) pt_of_ct[k++] = pos;
            }
        } else {
            // Tied columns: read row by row across the group, left to right.
            for (int r = 0; r < R; r++)
                for (int gi = 0; gi < g; gi++) {
                    int pos = r * K + group[gi];
                    if (pos < len) pt_of_ct[k++] = pos;
                }
        }
    }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Redefence (keyed rail fence)
// =====================================================================
//
// A rail fence (zigzag over `rails` rows, starting phase `offset`) whose rails are
// read off not top-to-bottom but in the keyed order `order[0..rails-1]`. With
// order == identity this is exactly decrypt_railfence; the solver climbs `order`.
// `variant` swaps read/write. out[] must not alias cipher[].
void decrypt_redefence(int cipher[], int len, int rails, int offset, int order[], int variant, int out[]) {

    if (rails <= 1 || rails >= len) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int P = 2 * (rails - 1);
    int pt_of_ct[MAX_CIPHER_LENGTH];
    int k = 0;

    // Emit the rails in keyed order; within each rail, positions left to right.
    for (int j = 0; j < rails; j++) {
        int target = order[j];
        for (int i = 0; i < len; i++) {
            int ph = (i + offset) % P;
            int rail = (ph < rails) ? ph : P - ph;
            if (rail == target) pt_of_ct[k++] = i;
        }
    }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Cadenus
// =====================================================================
//
// The plaintext is written row-major into a grid of `rows` rows (25 in the classic
// cipher) and K = len/rows columns. Each column is cyclically rotated upward, then
// the columns are reordered, and the grid is read row-major. The keyword supplies
// both the rotation amount and the column order; the solver decouples them and
// searches `order` (the read-order permutation, length K) and `rot` (the per-column
// upward rotation, 0..rows-1) freely, which subsumes any keyword/alphabet
// convention. `variant` swaps read/write. out[] must not alias cipher[].
//
// Read column p (0..K-1) comes from original column c = order[p]; output row r of
// that column comes from original grid row (r + rot[c]) mod rows.
void decrypt_cadenus(int cipher[], int len, int K, int order[], int rot[], int variant, int out[]) {

    if (K <= 1 || K > len || len % K != 0) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int rows = len / K;
    int pt_of_ct[MAX_CIPHER_LENGTH];

    for (int r = 0; r < rows; r++)
        for (int p = 0; p < K; p++) {
            int c = order[p];
            int src_row = (r + rot[c]) % rows;
            pt_of_ct[r * K + p] = src_row * K + c;
        }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Nihilist transposition
// =====================================================================
//
// The plaintext fills an N x N grid row-major; `rowperm`/`colperm` permute the rows
// and columns (the classic Nihilist transposition ties them to one keyword, but we
// keep them independent so the solver covers that and the two-key variant), and the
// grid is read off either row-major (readmode 0) or column-major (readmode 1).
// N = sqrt(len). Grid cell (r,c) holds plaintext rowperm[r]*N + colperm[c].
// `variant` swaps read/write. out[] must not alias cipher[].
void decrypt_nihilist(int cipher[], int len, int N, int rowperm[], int colperm[],
                      int readmode, int variant, int out[]) {

    if (N < 2 || N * N != len) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int pt_of_ct[MAX_CIPHER_LENGTH];
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int ct_pos = (readmode == 1) ? (c * N + r) : (r * N + c);
            pt_of_ct[ct_pos] = rowperm[r] * N + colperm[c];
        }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Swagman
// =====================================================================
//
// The plaintext fills an N-row x W-column grid row-major (W = len/N). An N x N key
// square (each column a permutation of 0..N-1) is applied column-wise and repeats
// every N columns: plaintext column jj is sorted vertically by the digits of square
// column (jj mod N) -- output row i takes the plaintext row whose square digit is i.
// `square[r*N + j]` is the digit at row r, column j of the key square. `readmode`
// selects how the rearranged grid is read off: 0 = row-major, 1 = column-major
// (sources differ, so the solver sweeps it). `variant` swaps read/write. out[] must
// not alias cipher[].
void decrypt_swagman(int cipher[], int len, int N, int square[], int readmode, int variant, int out[]) {

    if (N < 2 || N > 7 || len % N != 0) {   // Swagman squares are 3x3..7x7
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    int W = len / N;
    int pt_of_ct[MAX_CIPHER_LENGTH];

    // pos_of_digit for each of the N square columns: pod[j][d] = row r with
    // square[r][j] == d.
    int pod[7][7];                      // N <= 7 for Swagman
    for (int j = 0; j < N; j++)
        for (int r = 0; r < N; r++)
            pod[j][square[r * N + j]] = r;

    for (int i = 0; i < N; i++) {        // output row within a column
        for (int jj = 0; jj < W; jj++) { // plaintext/cipher column
            int src_row = pod[jj % N][i];
            int pt_pos  = src_row * W + jj;
            int ct_pos  = (readmode == 1) ? (jj * N + i) : (i * W + jj);
            pt_of_ct[ct_pos] = pt_pos;
        }
    }

    apply_perm(cipher, pt_of_ct, len, variant, out);
}


// =====================================================================
//  Turning grille
// =====================================================================
//
// An N x N turning grille has, for each rotation orbit of cells (the 4 cells a cell
// visits under repeated 90 degrees turns), exactly one open hole. The plaintext is
// written through the holes, the grille is turned 90 degrees, and so on for four
// turns, filling every cell once; the ciphertext is the grid read row-major.
//
// `key[orbit]` in {0,1,2,3} chooses, for each orbit, which turn exposes that
// orbit's hole. Orbits are numbered by their representative cell (the smallest
// row-major index among the four), in increasing representative order, so the key
// indexing is stable across calls for a given N. The number of orbits is written to
// *n_orbits when non-NULL (the solver needs it as the key length). `variant` swaps
// read/write. out[] must not alias cipher[].
//
// A cell at rotation index `ridx` within its orbit is exposed at turn t iff
// ridx == (key[orbit] + t) mod 4, so across the four turns every orbit cell is
// written exactly once. Odd N leaves a single centre cell (an orbit of size one),
// handled by the same rule. The grille's grid read row-major IS the ciphertext, so
// pt_of_ct[cell] = the plaintext index written into that cell.
void decrypt_grille(int cipher[], int len, int N, int key[], int variant, int out[], int *n_orbits) {

    int side = N;
    if (side < 2 || side * side != len) {
        if (n_orbits) *n_orbits = 0;
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        return;
    }

    // For each cell compute its orbit representative (min row-major index over the
    // four 90-degree CW rotations (r,c)->(c, N-1-r)) and its rotation index.
    int rep[MAX_CIPHER_LENGTH], ridx[MAX_CIPHER_LENGTH], orbit[MAX_CIPHER_LENGTH];
    for (int idx = 0; idx < len; idx++) {
        int r = idx / side, c = idx % side;
        int best = idx, bestk = 0;
        for (int t = 1; t < 4; t++) {
            int nr = c, nc = side - 1 - r;     // rotate 90 CW
            r = nr; c = nc;
            int j = r * side + c;
            if (j < best) { best = j; bestk = t; }
        }
        rep[idx] = best;
        // rotation index of idx relative to its representative: turns from rep to idx
        ridx[idx] = (4 - bestk) % 4;
    }

    // Assign stable orbit ids in representative-index order.
    int norb = 0;
    for (int idx = 0; idx < len; idx++) {
        if (rep[idx] == idx) { orbit[idx] = norb++; }
    }
    for (int idx = 0; idx < len; idx++) orbit[idx] = orbit[rep[idx]];
    if (n_orbits) *n_orbits = norb;

    // Walk the four turns; at each, scan cells row-major and assign the next
    // plaintext index to every exposed cell.
    int pt_of_ct[MAX_CIPHER_LENGTH];
    int m = 0;
    for (int t = 0; t < 4; t++)
        for (int idx = 0; idx < len; idx++)
            if (ridx[idx] == ((key[orbit[idx]] + t) & 3))
                pt_of_ct[idx] = m++;

    apply_perm(cipher, pt_of_ct, len, variant, out);
}

