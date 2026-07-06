// End-to-end regression for the double columnar transposition divide-and-conquer
// solver (dct_solve_core). Plants double-transposition ciphers of several key-length
// pairs and lengths at fixed seeds, runs the full screen -> refine -> finish pipeline
// (with trimmed budgets so the test stays a few seconds), and asserts the recovered
// plaintext matches the planted plaintext. Also prints a length/key-length capability
// characterization -- the regime where the IDP divide-and-conquer succeeds, well past
// the ~key-length-12-15 ceiling of parallel hill climbing (TRANSCOL2).
//
// Built with -DCOLOSSUS_NO_MAIN and linked against the full solver, like the other
// *_solver regressions.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "colossus.h"
#include "double_transposition_solver.h"
#include "scoring.h"

static int failures = 0;

// A long English plaintext source (trimmed per test).
static const char *SRC =
    "THEDOUBLETRANSPOSITIONCIPHERWASCONSIDEREDTOBEONEOFTHEMOSTSECURETYPESOFMAN"
    "UALCIPHERSITWASEXTENSIVELYUSEDINBOTHWORLDWARSANDDURINGTHECOLDWARINNINETEE"
    "NNINETYNINEOTTOLEIBERICHTHEFORMERHEADOFTHEGERMANFEDERALOFFICEFORINFORMATI"
    "ONSECURITYSUGGESTEDTHATADOUBLETRANSPOSITIONCHALLENGEBEPUBLISHEDWITHSPECIF"
    "ICPARAMETERSDESIGNEDTOENSUREITSSECURITYSUCHAWASPUBLISHEDBYKLAUSSCHMEHINTW"
    "OTHOUSANDSEVENHECHOSEANENGLISHPLAINTEXTANDENCRYPTEDITUSINGTWOTRANSPOSITIO"
    "NKEYSBOTHDERIVEDFROMENGLISHKEYPHRASESLONGERTHANTWENTYTHELENGTHOFTHEPLAINT"
    "EXTWASFIVEHUNDREDNINETYNINETHECHALLENGEWASRANKEDNUMBERFIVEINALISTOFUNSOLV"
    "EDCIPHERSANDWASCONSIDEREDTOBEUNBREAKABLEBYMANYEXPERTSINTHEFIELDOFCLASSICA"
    "LCRYPTANALYSISFORMORETHANADECADEUNTILITWASFINALLYSOLVED";

// An English plaintext WITH spaces and punctuation, for the space-significant case
// (the transposition grid carries the spaces/punctuation as real cells).
static const char *SRC_SPACES =
    "THE DOUBLE TRANSPOSITION CIPHER, LONG THOUGHT SECURE, MOVES EVERY CELL OF THE "
    "MESSAGE -- LETTERS AND SPACES ALIKE -- THROUGH TWO KEYED COLUMNAR STAGES. WHEN "
    "THE SPACES ARE PART OF THE TEXT, THEY RIDE THROUGH THE GRID JUST LIKE LETTERS.";

static int keyword_order(const char *kw, int order[]) {
    int lets[MAX_COLS], K = 0;
    for (const char *k = kw; *k && K < MAX_COLS; k++) {
        int c = toupper((unsigned char) *k);
        if (c >= 'A' && c <= 'Z') lets[K++] = c - 'A';
    }
    int used[MAX_COLS];
    for (int i = 0; i < K; i++) used[i] = 0;
    for (int j = 0; j < K; j++) {
        int best = -1;
        for (int c = 0; c < K; c++)
            if (!used[c] && (best < 0 || lets[c] < lets[best])) best = c;
        used[best] = 1; order[j] = best;
    }
    return K;
}

static void columnar_encrypt(const int in[], int len, int K, const int order[], int out[]) {
    static int grid[MAX_CIPHER_LENGTH];
    int R = (len + K - 1) / K, rem = len % K, pos = 0;
    for (int r = 0; r < R; r++)
        for (int c = 0; c < K; c++) {
            int h = (rem == 0 || c < rem) ? R : R - 1;
            if (r < h) grid[r * K + c] = in[pos++];
        }
    int o = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j], h = (rem == 0 || c < rem) ? R : R - 1;
        for (int r = 0; r < h; r++) out[o++] = grid[r * K + c];
    }
}

// Plant a double-transposition cipher (first `len` letters of SRC), solve it with the
// given budget, and return the fraction of plaintext positions recovered. If crib_cover
// > 0, supply the first crib_cover plaintext letters as a crib (over plaintext positions).
static double run_case(const char *kw1, const char *kw2, int len, unsigned long long seed,
                       int lo, int hi, int restarts, int iters, int crib_cover,
                       int *out_L1, int *out_L2) {
    static int P[MAX_CIPHER_LENGTH], mid[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH];
    int srclen = (int) strlen(SRC);
    if (len > srclen) len = srclen;
    for (int i = 0; i < len; i++) P[i] = SRC[i] - 'A';

    int o1[MAX_COLS], o2[MAX_COLS];
    int L1 = keyword_order(kw1, o1), L2 = keyword_order(kw2, o2);
    columnar_encrypt(P,   len, L1, o1, mid);
    columnar_encrypt(mid, len, L2, o2, C);

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TRANSCOL2_DC;
    cfg.ngram_size = 4;
    strcpy(cfg.ngram_file, "english_quadgrams.txt");
    cfg.min_cols = lo; cfg.max_cols = hi;
    cfg.read_direction = COL_READ_TB;
    cfg.n_restarts = restarts; cfg.n_hill_climbs = iters;

    SharedData shared;
    memset(&shared, 0, sizeof shared);
    g_ngram_logprob = true;                     // -logprob fitness (letters preserved anyway)
    shared.ngram_data = load_ngrams(cfg.ngram_file, cfg.ngram_size, false);
    double *bg = dct_load_bigrams(cfg.ngram_file, cfg.ngram_size);

    // Optional crib: the first crib_cover plaintext letters at their plaintext positions.
    int crib_idx[MAX_CIPHER_LENGTH], crib_pos[MAX_CIPHER_LENGTH], n_cribs = 0;
    for (int i = 0; i < crib_cover && i < len; i++) { crib_pos[n_cribs] = i; crib_idx[n_cribs] = P[i]; n_cribs++; }

    rng_state = seed;
    DCTResult res;
    dct_solve_core(C, len, &cfg, shared.ngram_data, bg,
                   n_cribs ? crib_idx : NULL, n_cribs ? crib_pos : NULL, n_cribs, &res);

    int match = 0;
    for (int i = 0; i < len; i++) if (res.plaintext[i] == P[i]) match++;
    if (out_L1) *out_L1 = res.L1;
    if (out_L2) *out_L2 = res.L2;
    free(shared.ngram_data); free(bg);
    return (double) match / len;
}

// Space-significant variant: plant a double columnar over the first `len` chars of
// SRC_SPACES with the NON-LETTERS carried as reversible sentinels (letters -> 0..25,
// other bytes -> -(byte)-1, exactly like ord()), so the spaces/punctuation ride the grid
// as real cells. Solve and return the fraction of ALL positions (letters and non-letters)
// recovered -- the solver must reproduce the punctuation in place.
static double run_case_spaces(const char *kw1, const char *kw2, int len,
                              unsigned long long seed, int lo, int hi,
                              int restarts, int iters, int *out_L1, int *out_L2) {
    static int P[MAX_CIPHER_LENGTH], mid[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH];
    int srclen = (int) strlen(SRC_SPACES);
    if (len > srclen) len = srclen;
    for (int i = 0; i < len; i++) {
        unsigned char ch = (unsigned char) SRC_SPACES[i];
        int c = toupper(ch);
        P[i] = (c >= 'A' && c <= 'Z') ? (c - 'A') : (-(int) ch - 1);
    }

    int o1[MAX_COLS], o2[MAX_COLS];
    int L1 = keyword_order(kw1, o1), L2 = keyword_order(kw2, o2);
    columnar_encrypt(P,   len, L1, o1, mid);
    columnar_encrypt(mid, len, L2, o2, C);

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TRANSCOL2_DC;
    cfg.ngram_size = 4;
    strcpy(cfg.ngram_file, "english_quadgrams.txt");
    cfg.min_cols = lo; cfg.max_cols = hi;
    cfg.read_direction = COL_READ_TB;
    cfg.n_restarts = restarts; cfg.n_hill_climbs = iters;

    SharedData shared;
    memset(&shared, 0, sizeof shared);
    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(cfg.ngram_file, cfg.ngram_size, false);
    double *bg = dct_load_bigrams(cfg.ngram_file, cfg.ngram_size);

    rng_state = seed;
    DCTResult res;
    dct_solve_core(C, len, &cfg, shared.ngram_data, bg, NULL, NULL, 0, &res);

    int match = 0;
    for (int i = 0; i < len; i++) if (res.plaintext[i] == P[i]) match++;
    if (out_L1) *out_L1 = res.L1;
    if (out_L2) *out_L2 = res.L2;
    free(shared.ngram_data); free(bg);
    return (double) match / len;
}

int main(void) {
    printf("test_double_transposition_solver:\n");
    init_alphabet(NULL);

    // --- capability: recover several key-length pairs (tight column ranges,
    // lean budgets so the test stays a few seconds). The divide-and-conquer solver
    // needs K2 recovered EXACTLY for the rotation-resolved finish to lock K1, so
    // these assert full recovery -- key lengths (6-10) well past TRANSCOL2's ceiling.
    // Lengths are chosen NOT to be multiples of a key length (incomplete rectangles,
    // the realistic / harder case the IDP start-windows handle).
    struct { const char *k1, *k2; int len, lo, hi; } cases[] = {
        { "PLANET",     "OCEANS",     320, 5,  7 },   //  6 x 6
        { "SUBMARINE",  "WATERFALL",  350, 8, 10 },   //  9 x 9
        { "HELICOPTER", "BLACKSMITH", 410, 9, 11 },   // 10 x 10 (distinct positions)
    };
    int nc = (int) (sizeof cases / sizeof cases[0]);
    for (int i = 0; i < nc; i++) {
        int L1, L2;
        double frac = run_case(cases[i].k1, cases[i].k2, cases[i].len,
                               12345ULL + 1000ULL * i, cases[i].lo, cases[i].hi,
                               14, 4500, 0, &L1, &L2);
        printf("  case %-11s x %-11s len=%d -> %.1f%% (recovered K1=%d K2=%d)\n",
            cases[i].k1, cases[i].k2, cases[i].len, 100.0 * frac, L1, L2);
        if (frac < 0.99) { printf("    FAIL: expected full recovery\n"); failures++; }
    }

    // --- length cliff characterization (K1=SUBMARINE 9 x K2=WATERFALL 9) -------
    printf("  length cliff (K1=SUBMARINE 9 x K2=WATERFALL 9):\n");
    int lens[] = { 250, 350, 440 };
    for (int i = 0; i < (int)(sizeof lens/sizeof lens[0]); i++) {
        int L1, L2;
        double frac = run_case("SUBMARINE", "WATERFALL", lens[i], 4242ULL,
                               8, 10, 14, 4500, 0, &L1, &L2);
        printf("    len=%d -> %.1f%%\n", lens[i], 100.0 * frac);
        if (lens[i] == 440 && frac < 0.99) {
            printf("    FAIL: expected full recovery at len=440\n"); failures++;
        }
    }

    // --- crib support: a crib (known leading plaintext) GUIDES the K2 climb (via the
    // rotation-resolved crib probe), so it rescues a solve at a budget too small for
    // pure IDP. Contrast no-crib vs a 40-letter crib at a tiny 4x2000 budget on the same
    // cipher/seed; only the crib run is asserted to recover (no-crib is printed for
    // contrast and may occasionally succeed on luck).
    {
        int L1, L2;
        double no_crib   = run_case("SUBMARINE", "WATERFALL", 350, 20260706ULL,
                                    8, 10, 4, 2000, 0,  &L1, &L2);
        double with_crib = run_case("SUBMARINE", "WATERFALL", 350, 20260706ULL,
                                    8, 10, 4, 2000, 40, &L1, &L2);
        printf("  crib-guided K2 (tiny 4x2000 budget): no-crib=%.1f%%  with-crib=%.1f%% "
               "(K1=%d K2=%d)\n", 100.0 * no_crib, 100.0 * with_crib, L1, L2);
        if (with_crib < 0.99) { printf("    FAIL: a 40-letter crib should rescue the solve\n"); failures++; }
    }

    // --- space-significant: the transposition grid carries spaces/punctuation as real
    // cells, and the solver must recover them in place (letters-only IDP + n-gram scoring,
    // sentinels reversed on output). Full recovery is asserted at a moderate budget.
    {
        int L1, L2;
        double frac = run_case_spaces("CIPHERTEXT", "TRANSPOSE",
                                      (int) strlen(SRC_SPACES), 20260706ULL,
                                      9, 10, 30, 9000, &L1, &L2);
        printf("  spaces/punctuation carried through grid: %.1f%% (K1=%d K2=%d, len=%d w/ spaces)\n",
               100.0 * frac, L1, L2, (int) strlen(SRC_SPACES));
        if (frac < 0.99) { printf("    FAIL: expected full recovery incl. spaces/punctuation\n"); failures++; }
    }

    if (failures == 0) printf("  all double-transposition solver tests passed.\n");
    else               printf("  %d CHECK(S) FAILED.\n", failures);
    return failures ? 1 : 0;
}
