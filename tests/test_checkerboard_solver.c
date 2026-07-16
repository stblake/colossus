//
//  In-process stress / limits tests for the Checkerboard solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  The Checkerboard maps each plaintext letter to a (row label, column label) digraph over a keyed
//  5x5 square. Once the per-axis label grouping is fixed the search is a free 25-code -> 25-letter
//  bijection (an Aristocrat over 25 symbols) on the homophonic incremental fast path. The SIMPLE
//  case (one label/axis) recovers like a 25-letter Aristocrat; the COMPLEX case (two labels/axis,
//  homophonic) additionally needs the chi-squared pre-pass to rank the true label PAIRING first,
//  which takes far more text than the ACA 60-90 range. This suite (also the basis for tuning the
//  SearchDefaults schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. blind case detection (checkerboard_detect) across the four simple/complex combos;
//    3. the pre-pass RANK curve -- rank of the true pairing vs length (the complex viability test);
//    4. the SIMPLE length cliff (assert at the reliable length, print the ACA-range margins);
//    5. a SIMPLE multi-keyword sweep (mean/worst);
//    6. the COMPLEX capability at a long length (bounded budget; breakage catch);
//    7. per-scheme calibration (anneal / shotgun / pso) on the simple case.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"                // apply_cipher_defaults
#include "scoring.h"               // load_ngrams
#include "checkerboard_solver.h"   // checkerboard_pairing_rank
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// ~1000 letters of natural English (Dickens / Austen), letters + spaces (spaces are stripped).
static const char *TEXT =
    "it was the best of times it was the worst of times it was the age of wisdom it was the age of "
    "foolishness it was the epoch of belief it was the epoch of incredulity it was the season of light "
    "it was the season of darkness it was the spring of hope it was the winter of despair we had "
    "everything before us we had nothing before us we were all going direct to heaven we were all going "
    "direct the other way it is a truth universally acknowledged that a single man in possession of a "
    "good fortune must be in want of a wife however little known the feelings or views of such a man may "
    "be on his first entering a neighbourhood this truth is so well fixed in the minds of the surrounding "
    "families that he is considered the rightful property of some one or other of their daughters";

// A..Z string -> alphabet indices, merging J into I. Returns count.
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

// Plant a Checkerboard cipher: encipher pt_len letters of TEXT under the keyed square + labels.
// rk2/ck2 NULL => that axis is simple. Fills prepared[] (expected plaintext, 0..24) and cipher_str.
// If ct_out non-NULL, also writes the 2*n label-letter indices. If true_row/true_col non-NULL, writes
// label letter -> true line. choice_seed makes the complex 2x2 choice reproducible. Returns n.
static int plant(const char *sq_kw, const char *rk1, const char *ck1, const char *rk2, const char *ck2,
                 int pt_len, uint32_t choice_seed, int prepared[], char cipher_str[],
                 int ct_out[], int true_row[], int true_col[]) {
    int side = CHECKERBOARD_SIDE;
    static int raw[MAX_CIPHER_LENGTH];
    int n = 0;
    for (int i = 0; TEXT[i] && n < pt_len; i++) {
        int c = toupper((unsigned char) TEXT[i]);
        if (c == 'J') c = 'I';
        if (c < 'A' || c > 'Z') continue;
        raw[n] = g_char_to_idx[c]; prepared[n] = raw[n]; n++;
    }

    int kw[64]; int kwn = str_to_idx(sq_kw, kw);
    int grid[CHECKERBOARD_GRID];
    checkerboard_square_from_keyword(kw, kwn, CB_ROUTE_SPIRAL_CW, grid, CHECKERBOARD_GRID);

    int nrl = rk2 ? 2 : 1, ncl = ck2 ? 2 : 1;
    int r1[8], r2[8], c1[8], c2[8];
    str_to_idx(rk1, r1); str_to_idx(ck1, c1);
    if (rk2) str_to_idx(rk2, r2);
    if (ck2) str_to_idx(ck2, c2);
    int rowlbl[CHECKERBOARD_MAX_LABELS], collbl[CHECKERBOARD_MAX_LABELS];
    for (int r = 0; r < side; r++) { rowlbl[r * nrl] = r1[r]; if (nrl > 1) rowlbl[r * nrl + 1] = r2[r]; }
    for (int c = 0; c < side; c++) { collbl[c * ncl] = c1[c]; if (ncl > 1) collbl[c * ncl + 1] = c2[c]; }

    if (true_row) {
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++) true_row[i] = -1;
        for (int r = 0; r < side; r++) { true_row[r1[r]] = r; if (rk2) true_row[r2[r]] = r; }
    }
    if (true_col) {
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++) true_col[i] = -1;
        for (int c = 0; c < side; c++) { true_col[c1[c]] = c; if (ck2) true_col[c2[c]] = c; }
    }

    static int ct[2 * MAX_CIPHER_LENGTH];
    seed_rand(choice_seed);
    checkerboard_encrypt(raw, n, grid, side, rowlbl, nrl, collbl, ncl, ct);
    for (int i = 0; i < 2 * n; i++) {
        cipher_str[i] = index_to_char(ct[i]);
        if (ct_out) ct_out[i] = ct[i];
    }
    cipher_str[2 * n] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the budget; keep > 0 bounds the complex
// config count via -nperiods. `method` overrides the scheme.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int keep, int method, uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = CHECKERBOARD;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) cfg.n_restarts = r;
    if (h > 0) cfg.n_hill_climbs = h;
    if (keep > 0) cfg.n_periods = keep;   // bounds the complex config count (top-K pairings per axis)
    if (method == METHOD_PSO) { cfg.n_particles = 12; cfg.refine_steps = 4; }

    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* still proceed */ }
    seed_rand(seed);
    solve_cipher((char *) cipher_str, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);
    if (secs_out) *secs_out = ((double) clock() - t0) / CLOCKS_PER_SEC;

    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = CHECKERBOARD; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "checkerboard registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 200000,
        "checkerboard anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.1499 && cfg.init_temp < 0.1501,
        "checkerboard anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = CHECKERBOARD; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "checkerboard registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 60 && cfg.n_hill_climbs == 200000,
        "checkerboard shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. blind case detection --------------------------------------------------

static void test_detect(void) {
    struct { const char *rk2, *ck2; int wr, wc; const char *tag; } combos[] = {
        { NULL,   NULL,   0, 0, "simple/simple" },
        { "BLACK","WHITE",1, 1, "complex/complex" },
        { "BLACK",NULL,   1, 0, "complex/simple" },
        { NULL,   "WHITE",0, 1, "simple/complex" },
    };
    printf("\n[blind case detection @ ~250 letters]\n");
    int side = CHECKERBOARD_SIDE;
    for (int ci = 0; ci < 4; ci++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        static int ct[2 * MAX_CIPHER_LENGTH];
        int n = plant("KNIGHTSTEMPLAR", "HORSE", "GRAYS", combos[ci].rk2, combos[ci].ck2,
                      250, 0x1000u + ci, prepared, cs, ct, NULL, NULL);
        int nR, nC;
        checkerboard_detect(ct, 2 * n, &nR, &nC);
        int rc = (nR > side), cc = (nC > side);
        printf("  %-16s : rows=%d cols=%d -> %s/%s\n", combos[ci].tag, nR, nC,
            rc ? "complex" : "simple", cc ? "complex" : "simple");
        CHECK(rc == combos[ci].wr && cc == combos[ci].wc,
            "detect %s wrong: (%d,%d)", combos[ci].tag, nR, nC);
    }
}

// --- 3. pre-pass rank curve (the complex-viability calibration) ----------------

static void test_rank_curve(void) {
    struct { const char *sq, *rk1, *rk2, *ck1, *ck2; } sets[] = {
        { "KNIGHTSTEMPLAR", "HORSE", "BLACK", "GRAYS", "WHITE" },
        { "KRYPTOS",        "MOUTH", "BLAZE", "DUNCE", "WIGHT" },
        { "CIPHERMACHINE",  "PLUMB", "WRECK", "FIGHT", "DOZEN" },
    };
    int lens[] = { 90, 150, 250, 400, 600 };
    printf("\n[pre-pass rank of the TRUE pairing among 945 (0 = ranked first); complex case]\n");
    printf("  set  N=90         150        250        400        600\n");
    int rank600_worst = 0;
    for (int s = 0; s < 3; s++) {
        printf("  %-3d ", s);
        for (int li = 0; li < 5; li++) {
            static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
            static int ct[2 * MAX_CIPHER_LENGTH]; static int trow[MAX_ALPHABET_SIZE], tcol[MAX_ALPHABET_SIZE];
            int n = plant(sets[s].sq, sets[s].rk1, sets[s].ck1, sets[s].rk2, sets[s].ck2,
                          lens[li], 0x2000u + 31 * s + li, prepared, cs, ct, trow, tcol);
            int rr = checkerboard_pairing_rank(ct, 2 * n, 0, trow);
            int rc = checkerboard_pairing_rank(ct, 2 * n, 1, tcol);
            printf(" r%-3d/c%-3d ", rr, rc);
            if (lens[li] == 600 && (rr > rank600_worst)) rank600_worst = rr;
            if (lens[li] == 600 && (rc > rank600_worst)) rank600_worst = rc;
        }
        printf("\n");
    }
    // At 600 letters the pre-pass must rank the true pairing FIRST on both axes for every set.
    CHECK(rank600_worst == 0, "pre-pass did not rank the true pairing first at N=600 (worst rank %d)",
        rank600_worst);
}

// --- 4. simple length cliff ---------------------------------------------------

static void test_simple_cliff(void) {
    int lens[] = { 60, 90, 130, 200, 300 };
    printf("\n[simple length cliff: square KNIGHTSTEMPLAR, BLACK/WHITE, registry anneal]\n");
    double frac200 = 0.0, frac300 = 0.0;
    for (int li = 0; li < 5; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("KNIGHTSTEMPLAR", "BLACK", "WHITE", NULL, NULL,
                      lens[li], 0, prepared, cs, NULL, NULL, NULL);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 0, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (lens[li] == 200) frac200 = frac;
        if (lens[li] == 300) frac300 = frac;
    }
    CHECK(frac300 >= 0.98, "simple 300-letter recovered only %.3f", frac300);
    CHECK(frac200 >= 0.90, "simple 200-letter recovered only %.3f", frac200);
}

// --- 5. simple multi-keyword sweep --------------------------------------------

static void test_simple_multi(void) {
    struct { const char *sq, *rk, *ck; } sets[] = {
        { "KRYPTOS",       "BLACK", "WHITE" },
        { "ZEBRAS",        "GRAYS", "HORSE" },
        { "CIPHERMACHINE", "MOUTH", "RIGHT" },
        { "PALMERSTON",    "FIGHT", "BLOCK" },
    };
    int plen = 200;
    double sum = 0, worst = 1.0;
    printf("\n[simple multi-keyword @ %d letters, registry schedule]\n", plen);
    for (int k = 0; k < 4; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(sets[k].sq, sets[k].rk, sets[k].ck, NULL, NULL, plen, 0, prepared, cs, NULL, NULL, NULL);
        double frac = solve_and_frac(cs, prepared, n, 0, 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL);
        printf("  %-14s : %.1f%%\n", sets[k].sq, 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / 4, 100.0 * worst);
    CHECK(sum / 4 > 0.95, "simple multi-keyword mean too low: %.1f%%", 100.0 * sum / 4);
    CHECK(worst > 0.90, "simple multi-keyword worst too low: %.1f%%", 100.0 * worst);
}

// --- 6. complex capability (bounded budget; below the ACA range, a breakage catch) -------------

static void test_complex_capability(void) {
    int lens[] = { 300, 600 };
    printf("\n[complex length: set A (KNIGHTSTEMPLAR, HORSE/BLACK, GRAYS/WHITE), keep=2, 8x40000]\n");
    double frac600 = 0.0;
    for (int li = 0; li < 2; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("KNIGHTSTEMPLAR", "HORSE", "GRAYS", "BLACK", "WHITE",
                      lens[li], 0xB0A2Du + li, prepared, cs, NULL, NULL, NULL);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 8, 40000, 2, METHOD_DEFAULT, 0xC0FFEEu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (lens[li] == 600) frac600 = frac;
    }
    // Below the reliable length this is high-variance; assert a breakage catch, not accuracy.
    CHECK(frac600 > 0.40, "complex 600-letter recovered only %.1f%% (breakage?)", 100.0 * frac600);
}

// --- 7. per-scheme calibration (simple case) ----------------------------------

static void test_per_scheme(void) {
    int plen = 200;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHERMACHINE", "BLACK", "WHITE", NULL, NULL, plen, 0, prepared, cs, NULL, NULL, NULL);
    struct { int method; const char *name; int r, h; } M[] = {
        { METHOD_DEFAULT, "anneal ", 0, 0 },
        { METHOD_SHOTGUN, "shotgun", 60, 200000 },
        { METHOD_PSO,     "pso    ",  3,    400 },
    };
    printf("\n[per-scheme @ %d letters, keyword CIPHERMACHINE (pso bounded, not asserted)]\n", plen);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, M[m].r, M[m].h, 0, M[m].method, 0x5C8E0u + m, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet("J");                         // 25-letter (J->I) alphabet, as -type checkerboard forces
    CHECK(g_alpha == CHECKERBOARD_GRID, "alphabet size %d, expected %d", g_alpha, CHECKERBOARD_GRID);

    // A free 25-code bijection rides the log-probability quadgram fitness (-logprob).
    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_detect();
    test_rank_curve();
    test_simple_cliff();
    test_simple_multi();
    test_complex_capability();
    test_per_scheme();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
