//
//  Unit + stress tests for the Checkerboard primitives (keyed 5x5 square; plaintext letter ->
//  (row label, column label) digraph; simple = 1 label/axis, complex = 2 labels/axis homophonic).
//
//  Framework-free: build with `make test`, which links this against checkerboard.c + bifid.c +
//  utils.c. Exits non-zero if any check fails.
//
//  Strategy:
//   - The two ACA worked examples (p.39) are pinned as known-answer vectors, verified against the
//     printed square: the SIMPLE encrypt reproduces ct_a EXACTLY (it is RNG-free); the COMPLEX
//     case is pinned decrypt-only (encrypt randomizes the 2x2 choice). detect() is pinned on both
//     ACA ciphertexts, including the spec's own complex example that shows only 7 (not 10) column
//     labels -- so detection must threshold at > side, never == 2*side.
//   - The spiral route is pinned: square_from_keyword(KNIGHT, spiral) reproduces the printed grid.
//   - Stress round-trips (>=5000 each) over random squares x labels x lengths assert
//     decrypt(encrypt(pt)) == pt and n_valid == len, for simple / complex / both mixed / 6x6.
//   - Choice randomization: all 4 codes of a letter decrypt to it; a 4-bin chi-squared over 200k
//     single-letter encrypts catches a biased pick; a no-canonical-collapse guard catches a
//     deterministic pick (which would also collapse the observed label set).
//   - Label-order folding: permuting the label->line assignment and correspondingly permuting the
//     square's rows/columns leaves the decrypt unchanged (the "order is not identifiable" fact the
//     solver relies on) -- for both simple and complex.
//

#include "colossus.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

static int arrays_equal(const int a[], const int b[], int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return 0;
    return 1;
}

// A..Z string -> alphabet indices, merging J into I (the 25-letter convention). Returns count.
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c == 'J') c = 'I';
        if (c < 'A' || c > 'Z') continue;
        out[n++] = g_char_to_idx[c];
    }
    return n;
}

// Build the flat encrypt-label array (lbl[line*n_lbl + k]) from label index arrays.
static void enc_labels(int lbl[], int n_lbl, const int l1[], const int l2[], int side) {
    for (int r = 0; r < side; r++) {
        lbl[r * n_lbl + 0] = l1[r];
        if (n_lbl > 1) lbl[r * n_lbl + 1] = l2[r];
    }
}

// Build the decrypt label->line map (map[label letter] = line, -1 elsewhere) from label indices.
static void label_map(int map[], const int l1[], const int l2[], int n_lbl, int side) {
    for (int i = 0; i < MAX_ALPHABET_SIZE; i++) map[i] = -1;
    for (int r = 0; r < side; r++) {
        map[l1[r]] = r;
        if (n_lbl > 1) map[l2[r]] = r;
    }
}

// --- Known-answer vectors: the two ACA worked examples (p.39) ------------------
//
// Printed square (a), row-major: KNIGH / PQRST / OYZUA / MXWVB / LFEDC.
// Simple:  row labels BLACK, col labels WHITE.       pt NUMBERSCANALSOBEU.
// Complex: row labels HORSE(1)+BLACK(2), col GRAYS(1)+WHITE(2); same square, same pt.
static const char *KAT_SQUARE = "KNIGHPQRSTOYZUAMXWVBLFEDC";
static const char *KAT_PT     = "NUMBERSCANALSOBEU";
static const char *KAT_CT_A   = "BHATCWCEKILILTKEAEBHAEKWLTAWCEKIAT";
static const char *KAT_CT_B   = "HRRYCGSSEALAOTKSRSBRASEGLYAGCSEIAT";  // one valid ACA choice

static void test_known_answer(void) {
    int side = CHECKERBOARD_SIDE;
    int grid[CHECKERBOARD_GRID];
    int gn = str_to_idx(KAT_SQUARE, grid);
    CHECK(gn == CHECKERBOARD_GRID, "KAT square is not 25 letters (%d)", gn);

    int pt[64];  int ptn = str_to_idx(KAT_PT, pt);
    int cta[128]; int ctan = str_to_idx(KAT_CT_A, cta);
    int ctb[128]; int ctbn = str_to_idx(KAT_CT_B, ctb);
    CHECK(ctan == 2 * ptn && ctbn == 2 * ptn, "KAT ct lengths wrong (%d %d vs 2*%d)", ctan, ctbn, ptn);

    int blk[8], whi[8], hor[8], gry[8];
    str_to_idx("BLACK", blk); str_to_idx("WHITE", whi);
    str_to_idx("HORSE", hor); str_to_idx("GRAYS", gry);

    // KAT-1: simple encrypt is RNG-free -> reproduces ct_a EXACTLY.
    int rl[8], cl[8];
    enc_labels(rl, 1, blk, NULL, side);
    enc_labels(cl, 1, whi, NULL, side);
    int enc[128];
    checkerboard_encrypt(pt, ptn, grid, side, rl, 1, cl, 1, enc);
    CHECK(arrays_equal(enc, cta, ctan), "KAT-1: simple encrypt != ACA ct_a");

    // KAT-2: simple decrypt of ct_a recovers pt, all positions legal.
    int smap_r[MAX_ALPHABET_SIZE], smap_c[MAX_ALPHABET_SIZE];
    label_map(smap_r, blk, NULL, 1, side);
    label_map(smap_c, whi, NULL, 1, side);
    int dec[64];
    int nv = checkerboard_decrypt(cta, ctan, grid, side, smap_r, smap_c, dec);
    CHECK(arrays_equal(dec, pt, ptn) && nv == ptn, "KAT-2: simple decrypt != pt (nv=%d)", nv);

    // KAT-3: complex decrypt of the ACA ct_b recovers pt (encrypt is not pinned -- it randomizes).
    int cmap_r[MAX_ALPHABET_SIZE], cmap_c[MAX_ALPHABET_SIZE];
    label_map(cmap_r, hor, blk, 2, side);
    label_map(cmap_c, gry, whi, 2, side);
    nv = checkerboard_decrypt(ctb, ctbn, grid, side, cmap_r, cmap_c, dec);
    CHECK(arrays_equal(dec, pt, ptn) && nv == ptn, "KAT-3: complex decrypt != pt (nv=%d)", nv);

    // KAT-4: ct_a is the "always-2nd-label" complex encoding (BLACK/WHITE are the 2nd keywords),
    // so it decrypts under the COMPLEX label maps too.
    nv = checkerboard_decrypt(cta, ctan, grid, side, cmap_r, cmap_c, dec);
    CHECK(arrays_equal(dec, pt, ptn) && nv == ptn, "KAT-4: ct_a under complex maps != pt (nv=%d)", nv);

    // KAT-5: cross -- ct_b under the SIMPLE labels does NOT recover pt (its H/O/R/G/Y letters are
    // not simple labels), so some positions are illegal.
    nv = checkerboard_decrypt(ctb, ctbn, grid, side, smap_r, smap_c, dec);
    CHECK(!(arrays_equal(dec, pt, ptn)) && nv < ptn, "KAT-5: ct_b wrongly legal under simple labels (nv=%d)", nv);

    // KAT-6/7: detection. The simple example shows 5/5 labels; the ACA complex example shows
    // 10 row labels but only 7 column labels -- pins the ">side" rule, not "==2*side".
    int nR, nC;
    checkerboard_detect(cta, ctan, &nR, &nC);
    CHECK(nR == 5 && nC == 5, "KAT-6: detect(ct_a) = (%d,%d), expected (5,5)", nR, nC);
    checkerboard_detect(ctb, ctbn, &nR, &nC);
    CHECK(nR == 10 && nC == 7, "KAT-7: detect(ct_b) = (%d,%d), expected (10,7)", nR, nC);

    // KAT-8: all four homophone codes of a single letter decrypt to that letter. 'N' is at cell
    // (0,1): row labels {H,B}, col labels {R,H} -> codes HR,HH,BR,BH.
    int letter_N = g_char_to_idx['N'];
    int cells_r[2] = { hor[0], blk[0] };   // H, B (row 0)
    int cells_c[2] = { gry[1], whi[1] };   // R, H (col 1)
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++) {
            int digraph[2] = { cells_r[a], cells_c[b] };
            int out1[1];
            int v = checkerboard_decrypt(digraph, 2, grid, side, cmap_r, cmap_c, out1);
            CHECK(v == 1 && out1[0] == letter_N, "KAT-8: code (%d,%d) did not decrypt to N", a, b);
        }
}

// --- spiral route pins the ACA square ------------------------------------------
static void test_spiral_route(void) {
    int want[CHECKERBOARD_GRID];
    str_to_idx(KAT_SQUARE, want);
    int kw[8]; int kwn = str_to_idx("KNIGHT", kw);

    int grid[CHECKERBOARD_GRID];
    checkerboard_square_from_keyword(kw, kwn, CB_ROUTE_SPIRAL_CW, grid, CHECKERBOARD_GRID);
    CHECK(arrays_equal(grid, want, CHECKERBOARD_GRID), "spiral(KNIGHT) != ACA square");

    // Row-major fill is a different square (the sequence laid straight into rows).
    int rm[CHECKERBOARD_GRID];
    checkerboard_square_from_keyword(kw, kwn, CB_ROUTE_ROWMAJOR, rm, CHECKERBOARD_GRID);
    CHECK(!arrays_equal(rm, want, CHECKERBOARD_GRID), "row-major fill coincides with spiral");
}

// --- random square (a permutation of 0..n-1) -----------------------------------
static void random_square(int grid[], int n) {
    for (int i = 0; i < n; i++) grid[i] = i;
    shuffle(grid, n);
}

// --- random distinct labels: two disjoint sets of `side` letters drawn from 0..alpha-1 ---------
static void random_labels(int l1[], int l2[], int n_lbl, int side, int alpha) {
    int pool[CHECKERBOARD_MAX_GRID];
    for (int i = 0; i < alpha; i++) pool[i] = i;
    shuffle(pool, alpha);
    for (int r = 0; r < side; r++) l1[r] = pool[r];
    if (n_lbl > 1) for (int r = 0; r < side; r++) l2[r] = pool[side + r];
}

// --- stress round-trip for a given (n_row_lbl, n_col_lbl) combo -----------------
static void test_roundtrip(int n_row_lbl, int n_col_lbl, int side, int trials, const char *tag) {
    int alpha = side * side, ncell = alpha;
    for (int t = 0; t < trials; t++) {
        int grid[CHECKERBOARD_MAX_GRID];
        random_square(grid, ncell);
        int rl1[CHECKERBOARD_MAX_SIDE], rl2[CHECKERBOARD_MAX_SIDE];
        int cl1[CHECKERBOARD_MAX_SIDE], cl2[CHECKERBOARD_MAX_SIDE];
        random_labels(rl1, rl2, n_row_lbl, side, alpha);
        random_labels(cl1, cl2, n_col_lbl, side, alpha);
        int rl[CHECKERBOARD_MAX_LABELS], cl[CHECKERBOARD_MAX_LABELS];
        enc_labels(rl, n_row_lbl, rl1, rl2, side);
        enc_labels(cl, n_col_lbl, cl1, cl2, side);
        int mr[MAX_ALPHABET_SIZE], mc[MAX_ALPHABET_SIZE];
        label_map(mr, rl1, rl2, n_row_lbl, side);
        label_map(mc, cl1, cl2, n_col_lbl, side);

        int len = 1 + rand_int(0, 500);
        int pt[512], ct[1024], back[512];
        for (int i = 0; i < len; i++) pt[i] = rand_int(0, ncell);
        checkerboard_encrypt(pt, len, grid, side, rl, n_row_lbl, cl, n_col_lbl, ct);
        int nv = checkerboard_decrypt(ct, 2 * len, grid, side, mr, mc, back);
        CHECK(arrays_equal(back, pt, len) && nv == len,
            "%s round-trip mismatch (len=%d nv=%d)", tag, len, nv);
    }
}

// --- detection over the four combos x several lengths --------------------------
static void test_detect(void) {
    int side = CHECKERBOARD_SIDE, alpha = side * side;
    int combos[4][2] = { {1,1}, {2,2}, {2,1}, {1,2} };
    int lens[3] = { 90, 150, 300 };
    for (int ci = 0; ci < 4; ci++) {
        int nrl = combos[ci][0], ncl = combos[ci][1];
        for (int li = 0; li < 3; li++) {
            int len = lens[li], misses = 0;
            for (int t = 0; t < 400; t++) {
                int grid[CHECKERBOARD_GRID]; random_square(grid, alpha);
                int rl1[8], rl2[8], cl1[8], cl2[8];
                random_labels(rl1, rl2, nrl, side, alpha);
                random_labels(cl1, cl2, ncl, side, alpha);
                int rl[CHECKERBOARD_MAX_LABELS], cl[CHECKERBOARD_MAX_LABELS];
                enc_labels(rl, nrl, rl1, rl2, side);
                enc_labels(cl, ncl, cl1, cl2, side);
                int pt[512], ct[1024];
                for (int i = 0; i < len; i++) pt[i] = rand_int(0, alpha);
                checkerboard_encrypt(pt, len, grid, side, rl, nrl, cl, ncl, ct);
                int nR, nC;
                checkerboard_detect(ct, 2 * len, &nR, &nC);
                // A simple axis can NEVER exceed `side`; a complex axis is > side once enough
                // labels appear. Assert the classification (>side <=> complex) at len >= 150.
                int row_complex = (nR > side), col_complex = (nC > side);
                if (len >= 150 && (row_complex != (nrl > 1) || col_complex != (ncl > 1))) misses++;
            }
            if (len >= 150)
                CHECK(misses == 0, "detect misclassified %d/400 at combo(%d,%d) len=%d",
                    misses, nrl, ncl, len);
        }
    }
}

// --- choice randomization: uniform 2x2 pick (chi-squared) + all-4-valid ---------
static void test_choice(void) {
    int side = CHECKERBOARD_SIDE, alpha = side * side;
    int grid[CHECKERBOARD_GRID]; random_square(grid, alpha);
    int rl1[8], rl2[8], cl1[8], cl2[8];
    random_labels(rl1, rl2, 2, side, alpha);
    random_labels(cl1, cl2, 2, side, alpha);
    int rl[CHECKERBOARD_MAX_LABELS], cl[CHECKERBOARD_MAX_LABELS];
    enc_labels(rl, 2, rl1, rl2, side);
    enc_labels(cl, 2, cl1, cl2, side);

    // Pick one plaintext letter; find its cell (r,c) so we know its 4 legal codes.
    int letter = rand_int(0, alpha), cell = 0;
    for (int p = 0; p < alpha; p++) if (grid[p] == letter) cell = p;
    int r = cell / side, c = cell % side;

    const int M = 200000;
    int count[4] = {0,0,0,0};
    int pt1[1] = { letter }, out[2];
    for (int m = 0; m < M; m++) {
        checkerboard_encrypt(pt1, 1, grid, side, rl, 2, cl, 2, out);
        int kr = (out[0] == rl1[r]) ? 0 : 1;
        int kc = (out[1] == cl1[c]) ? 0 : 1;
        count[kr * 2 + kc]++;
    }
    double exp = M / 4.0, chi2 = 0.0;
    for (int b = 0; b < 4; b++) chi2 += (count[b] - exp) * (count[b] - exp) / exp;
    CHECK(chi2 < 21.0, "choice not uniform: chi2=%.2f counts={%d %d %d %d}",
        chi2, count[0], count[1], count[2], count[3]);
}

// --- no canonical collapse: a randomized complex ct differs from all 4 canonical picks and uses
//     all 10 row/col labels ------------------------------------------------------
static void test_no_canonical_collapse(void) {
    int side = CHECKERBOARD_SIDE, alpha = side * side;
    int grid[CHECKERBOARD_GRID]; random_square(grid, alpha);
    int rl1[8], rl2[8], cl1[8], cl2[8];
    random_labels(rl1, rl2, 2, side, alpha);
    random_labels(cl1, cl2, 2, side, alpha);
    int rl[CHECKERBOARD_MAX_LABELS], cl[CHECKERBOARD_MAX_LABELS];
    enc_labels(rl, 2, rl1, rl2, side);
    enc_labels(cl, 2, cl1, cl2, side);
    int pos[CHECKERBOARD_GRID];
    for (int p = 0; p < alpha; p++) pos[grid[p]] = p;

    int len = 400;
    int pt[400], ct[800];
    for (int i = 0; i < len; i++) pt[i] = rand_int(0, alpha);
    checkerboard_encrypt(pt, len, grid, side, rl, 2, cl, 2, ct);

    // Build the 4 canonical (fixed choice i,j) encodings and confirm the random ct matches none.
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            int canon[800];
            for (int m = 0; m < len; m++) {
                int cell = pos[pt[m]], rr = cell / side, cc = cell % side;
                canon[2*m]   = (i == 0) ? rl1[rr] : rl2[rr];
                canon[2*m+1] = (j == 0) ? cl1[cc] : cl2[cc];
            }
            CHECK(!arrays_equal(ct, canon, 2 * len),
                "randomized ct coincides with canonical pick (%d,%d)", i, j);
        }
    int nR, nC;
    checkerboard_detect(ct, 2 * len, &nR, &nC);
    CHECK(nR == 10 && nC == 10, "randomized 400-letter ct used only (%d,%d) labels", nR, nC);
}

// --- label-order folding: relabel the line assignment + permute the square's rows/cols, the
//     decrypt is unchanged (order is not identifiable) ---------------------------
static void test_order_folds(int n_row_lbl, int n_col_lbl, const char *tag) {
    int side = CHECKERBOARD_SIDE, alpha = side * side;
    for (int t = 0; t < 2000; t++) {
        int grid[CHECKERBOARD_GRID]; random_square(grid, alpha);
        int rl1[8], rl2[8], cl1[8], cl2[8];
        random_labels(rl1, rl2, n_row_lbl, side, alpha);
        random_labels(cl1, cl2, n_col_lbl, side, alpha);
        int rl[CHECKERBOARD_MAX_LABELS], cl[CHECKERBOARD_MAX_LABELS];
        enc_labels(rl, n_row_lbl, rl1, rl2, side);
        enc_labels(cl, n_col_lbl, cl1, cl2, side);
        int mr[MAX_ALPHABET_SIZE], mc[MAX_ALPHABET_SIZE];
        label_map(mr, rl1, rl2, n_row_lbl, side);
        label_map(mc, cl1, cl2, n_col_lbl, side);

        int len = 1 + rand_int(0, 200);
        int pt[256], ct[512], dec1[256], dec2[256];
        for (int i = 0; i < len; i++) pt[i] = rand_int(0, alpha);
        checkerboard_encrypt(pt, len, grid, side, rl, n_row_lbl, cl, n_col_lbl, ct);
        checkerboard_decrypt(ct, 2 * len, grid, side, mr, mc, dec1);

        // Random row perm pi and col perm sigma. mr'(letter) = pi(mr(letter)); grid'(r,c) =
        // grid(pi^{-1}(r), sigma^{-1}(c)). Decrypt must be identical.
        int pi[8], si[8], pinv[8], sinv[8];
        for (int i = 0; i < side; i++) { pi[i] = i; si[i] = i; }
        shuffle(pi, side); shuffle(si, side);
        for (int i = 0; i < side; i++) { pinv[pi[i]] = i; sinv[si[i]] = i; }

        int mr2[MAX_ALPHABET_SIZE], mc2[MAX_ALPHABET_SIZE];
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++) { mr2[i] = -1; mc2[i] = -1; }
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++) {
            if (mr[i] >= 0) mr2[i] = pi[mr[i]];
            if (mc[i] >= 0) mc2[i] = si[mc[i]];
        }
        int grid2[CHECKERBOARD_GRID];
        for (int r = 0; r < side; r++)
            for (int c = 0; c < side; c++)
                grid2[r * side + c] = grid[pinv[r] * side + sinv[c]];

        checkerboard_decrypt(ct, 2 * len, grid2, side, mr2, mc2, dec2);
        CHECK(arrays_equal(dec1, dec2, len), "%s order-fold changed the decrypt (len=%d)", tag, len);
    }
}

int main(void) {
    seed_rand(20240716u);
    init_alphabet("J");
    CHECK(g_alpha == CHECKERBOARD_GRID, "alphabet size %d, expected %d", g_alpha, CHECKERBOARD_GRID);

    test_known_answer();
    test_spiral_route();

    test_roundtrip(1, 1, CHECKERBOARD_SIDE, 5000, "simple");
    test_roundtrip(2, 2, CHECKERBOARD_SIDE, 5000, "complex");
    test_roundtrip(2, 1, CHECKERBOARD_SIDE, 2000, "mixed-rc");
    test_roundtrip(1, 2, CHECKERBOARD_SIDE, 2000, "mixed-cr");
    test_roundtrip(1, 1, 6, 2000, "6x6-simple");
    test_roundtrip(2, 2, 6, 2000, "6x6-complex");

    test_detect();
    test_choice();
    test_no_canonical_collapse();
    test_order_folds(1, 1, "simple");
    test_order_folds(2, 2, "complex");

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
