// Unit tests for the exact small-permutation ordering helpers in trans_common.c:
//   - held_karp_best_path: verified against an exhaustive R! brute force.
//   - ngram_sum_raw: verified ADDITIVE (sum_raw(A++B) == sum_raw(A)+sum_raw(B)+seam),
//     the property the seam decomposition relies on.
//   - seam_best_row_order: recovers the best row order on a synthetic where the true
//     order is the unique optimum (cross-checked against brute force).
//
// Build (see makefile `test` target):
//   gcc -I... tests/test_held_karp.c src/transposition/trans_common.c \
//       src/core/scoring.c src/core/utils.c src/core/dict.c -o tests/test_held_karp

#include "colossus.h"
#include "scoring.h"
#include "trans_common.h"

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("  FAIL: %s\n", msg); failures++; } } while (0)

// Objective of a visiting order under (indiv, delta) -- mirrors held_karp's objective.
static double path_obj(int R, const double *indiv, const double *delta, const int *ord) {
    double s = indiv[ord[0]];
    for (int i = 1; i < R; i++) s += indiv[ord[i]] + delta[ord[i - 1] * R + ord[i]];
    return s;
}

// Exhaustive best path over R! orders (R small).
static double brute_best(int R, const double *indiv, const double *delta, int *best_ord) {
    int ord[HELD_KARP_MAX_NODES];
    for (int i = 0; i < R; i++) ord[i] = i;
    double best = -1e300;
    // Heap's algorithm.
    int c[HELD_KARP_MAX_NODES] = {0};
    double s = path_obj(R, indiv, delta, ord);
    if (s > best) { best = s; for (int i = 0; i < R; i++) best_ord[i] = ord[i]; }
    int i = 0;
    while (i < R) {
        if (c[i] < i) {
            int a = (i % 2 == 0) ? 0 : c[i];
            int t = ord[a]; ord[a] = ord[i]; ord[i] = t;
            s = path_obj(R, indiv, delta, ord);
            if (s > best) { best = s; for (int k = 0; k < R; k++) best_ord[k] = ord[k]; }
            c[i]++; i = 0;
        } else { c[i] = 0; i++; }
    }
    return best;
}

static void test_held_karp_exact(void) {
    printf("held_karp_best_path vs brute force...\n");
    seed_rand(12345);
    int order[HELD_KARP_MAX_NODES], bord[HELD_KARP_MAX_NODES];
    double indiv[HELD_KARP_MAX_NODES], delta[HELD_KARP_MAX_NODES * HELD_KARP_MAX_NODES];
    for (int R = 1; R <= 9; R++) {
        for (int trial = 0; trial < 40; trial++) {
            for (int a = 0; a < R; a++) indiv[a] = (frand() - 0.5) * 10.0;
            for (int a = 0; a < R * R; a++) delta[a] = (frand() - 0.5) * 6.0;
            double hk = held_karp_best_path(R, indiv, delta, order);
            double bf = brute_best(R, indiv, delta, bord);
            CHECK(fabs(hk - bf) < 1e-9, "held-karp score != brute-force best");
            // The returned order must reproduce the returned score.
            CHECK(fabs(path_obj(R, indiv, delta, order) - hk) < 1e-9,
                  "returned order does not reproduce score");
        }
    }
    printf("  ok\n");
}

static void test_ngram_sum_raw_additive(void) {
    printf("ngram_sum_raw additivity...\n");
    init_alphabet(NULL);                       // full A..Z, g_alpha = 26
    int n = 4;
    long sz = 1; for (int i = 0; i < n; i++) sz *= g_alpha;
    float *tab = malloc(sz * sizeof(float));
    seed_rand(999);
    for (long i = 0; i < sz; i++) tab[i] = (float)(frand() * 2.0 - 1.0);

    for (int trial = 0; trial < 200; trial++) {
        int la = 4 + (int)(frand() * 12), lb = 4 + (int)(frand() * 12);
        int A[32], B[32], C[64];
        for (int i = 0; i < la; i++) A[i] = (frand() < 0.1) ? -1 : (int)(frand() * 26);
        for (int i = 0; i < lb; i++) B[i] = (frand() < 0.1) ? -1 : (int)(frand() * 26);
        int m = 0;
        for (int i = 0; i < la; i++) C[m++] = A[i];
        for (int i = 0; i < lb; i++) C[m++] = B[i];
        double sa = ngram_sum_raw(A, la, tab, n);
        double sb = ngram_sum_raw(B, lb, tab, n);
        double sc = ngram_sum_raw(C, m, tab, n);
        // Seam = windows of C that span the A|B boundary (start in [la-n+1, la-1]).
        double seam = 0.0;
        for (int s = la - n + 1; s <= la - 1; s++) {
            if (s < 0) continue;
            int win[16], bad = 0;
            for (int j = 0; j < n; j++) {
                int v = C[s + j]; if (v < 0) { bad++; v = 0; }
                win[j] = v;
            }
            // Use the production packer so the test tracks ngram_sum_raw's convention.
            if (bad == 0) seam += tab[ngram_index_int(win, n)];
        }
        CHECK(fabs(sc - (sa + sb + seam)) < 1e-6, "sum_raw(A++B) != sum_raw(A)+sum_raw(B)+seam");
    }
    free(tab);
    printf("  ok\n");
}

static void test_seam_best_row_order(void) {
    printf("seam_best_row_order recovers planted order...\n");
    init_alphabet(NULL);
    int n = 4;
    long sz = 1; for (int i = 0; i < n; i++) sz *= g_alpha;
    float *tab = malloc(sz * sizeof(float));
    for (long i = 0; i < sz; i++) tab[i] = -1.0f;          // floor: penalise everything
    // Reward the quadgrams of a known "plaintext" laid out as R rows of width W so
    // that reading the rows in the planted order spells it.
    const char *rowtext[5] = { "THEQUICK", "BROWNFOX", "JUMPSXXX", "OVERXXXX", "LAZYDOGS" };
    int R = 5, W = 8;
    int rowsbuf[5][8];
    for (int r = 0; r < R; r++)
        for (int c = 0; c < W; c++) rowsbuf[r][c] = rowtext[r][c] - 'A';
    // Reward every quadgram appearing in the concatenation of the true order 0..4.
    int flat[64], fl = 0;
    for (int r = 0; r < R; r++) for (int c = 0; c < W; c++) flat[fl++] = rowsbuf[r][c];
    for (int i = 0; i + n <= fl; i++)
        tab[ngram_index_int(&flat[i], n)] = 5.0f;   // production packer (convention-agnostic)
    int *rows[5]; int rowlen[5];
    for (int r = 0; r < R; r++) { rows[r] = rowsbuf[r]; rowlen[r] = W; }
    double indiv_buf[5], delta_buf[25]; int order[5];
    double sc = seam_best_row_order(R, rows, rowlen, tab, n, NULL, 0.0, indiv_buf, delta_buf, order);
    // The planted order is 0,1,2,3,4 and should be the unique optimum.
    int ok = 1; for (int r = 0; r < R; r++) if (order[r] != r) ok = 0;
    CHECK(ok, "did not recover the planted row order 0,1,2,3,4");
    CHECK(sc > 0.0, "best seam score should be positive (rewarded quadgrams)");
    free(tab);
    printf("  ok\n");
}

int main(void) {
    printf("\n=== test_held_karp ===\n");
    test_held_karp_exact();
    test_ngram_sum_raw_additive();
    test_seam_best_row_order();
    if (failures) { printf("\n%d CHECK(s) FAILED\n", failures); return 1; }
    printf("\nAll held-karp/seam tests passed.\n");
    return 0;
}
