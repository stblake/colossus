//
//  In-process stress / limits tests for the Aristocrat / Patristocrat solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  The Aristocrat is a simple monoalphabetic substitution (26-letter bijection pt = key[ct]) with
//  the word divisions preserved; the Patristocrat is the identical cipher with the spaces removed.
//  Both run ONE solver: a free 26-permutation climbed by n-gram score with the homophonic-style
//  incremental fast path (each swap scored as a delta). A free 26-sub is weaker signal than
//  Ragbaby's known-shift alphabet, so it rides the log-probability quadgram fitness (-logprob) and
//  recovers reliably across the ACA ~80-150-letter range. This suite (also the basis for tuning the
//  SearchDefaults 12x200000 / inittemp-0.15 schedule) checks:
//    1. registry validation (apply_cipher_defaults) for BOTH types + a non-registry type untouched;
//    2. the substitution move INVARIANT (aristocrat_move keeps a 26-permutation);
//    3. a capability floor (recovery %) across keywords in the ACA range;
//    4. a length cliff (recovery vs length);
//    5. a multi-keyword sweep (mean/worst recovery);
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso;
//    7. a Patristocrat smoke test (same core, no word divisions).
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"             // apply_cipher_defaults
#include "scoring.h"            // load_ngrams
#include "aristocrat_solver.h"  // aristocrat_move
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

// Word-divided natural English (letters + single spaces).
static const char *PLAINTEXT =
    "it is a truth universally acknowledged that a single man in possession of a good fortune must "
    "be in want of a wife however little known the feelings or views of such a man may be on his "
    "first entering a neighbourhood this truth is so well fixed in the minds of the surrounding "
    "families that he is considered the rightful property of some one or other of their daughters "
    "my dear mister bennet said his lady to him one day have you heard that netherfield park is let "
    "at last mister bennet replied that he had not but it is returned she for missus long has just "
    "been here and she told me all about it mister bennet made no answer do you not want to know "
    "who has taken it cried his wife impatiently you want to tell me and i have no objection to hear";

// Plant an Aristocrat cipher: take whole words of PLAINTEXT up to ~pt_len letters, encipher under
// the keyword's K2 keyed alphabet, and build the SPACED cipher string the solver decodes. Fills
// prepared[] (the expected letters-only plaintext) and cipher_str; returns the letter count n.
static int plant(const char *kw, int pt_len, int prepared[], char cipher_str[]) {
    char pt[MAX_CIPHER_LENGTH]; int m = 0, letters_so_far = 0;
    for (int i = 0; PLAINTEXT[i] && letters_so_far < pt_len; i++) {
        char c = PLAINTEXT[i];
        if (c == ' ') { if (m > 0 && pt[m - 1] != ' ') pt[m++] = ' '; continue; }
        pt[m++] = c;
        letters_so_far++;
    }
    while (m > 0 && pt[m - 1] == ' ') m--;
    pt[m] = '\0';

    int cmap[26];
    aristocrat_build_map(ARIST_K2, kw, kw, 5, cmap);

    int n = 0, cn = 0; bool prev_space = false, started = false;
    for (int i = 0; pt[i]; i++) {
        if (pt[i] == ' ') { prev_space = true; continue; }
        int idx = g_char_to_idx[(int) toupper((unsigned char) pt[i])];
        if (idx < 0 || idx >= 26) continue;
        prepared[n] = idx;
        if (prev_space && started) cipher_str[cn++] = ' ';
        cipher_str[cn++] = index_to_char(cmap[idx]);
        prev_space = false; started = true;
        n++;
    }
    cipher_str[cn] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the budget (else the registry schedule
// is used). `method` overrides the search scheme; `type` selects Aristocrat vs Patristocrat.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen, int type,
        int r, int h, int method, uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) { cfg.n_restarts = r; }
    if (h > 0) { cfg.n_hill_climbs = h; }
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
    int types[2] = { ARISTOCRAT, PATRISTOCRAT };
    for (int ti = 0; ti < 2; ti++) {
        init_config(&cfg); cfg.cipher_type = types[ti]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "type %d registry: no entry applied", types[ti]);
        CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 200000,
            "type %d anneal defaults wrong: %dx%d", types[ti], cfg.n_restarts, cfg.n_hill_climbs);
        CHECK(cfg.init_temp > 0.1499 && cfg.init_temp < 0.1501,
            "type %d anneal inittemp wrong: %.4f", types[ti], cfg.init_temp);

        init_config(&cfg); cfg.cipher_type = types[ti]; cfg.method = METHOD_SHOTGUN;
        CHECK(apply_cipher_defaults(&cfg, false), "type %d registry (shotgun): no entry", types[ti]);
        CHECK(cfg.n_restarts == 60 && cfg.n_hill_climbs == 200000,
            "type %d shotgun defaults wrong: %dx%d", types[ti], cfg.n_restarts, cfg.n_hill_climbs);
    }

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. substitution move invariant -------------------------------------------

static void test_move_invariant(void) {
    int bad_perm = 0;
    seed_rand(0xA11CEu);
    for (int trial = 0; trial < 400; trial++) {
        int key[26];
        for (int i = 0; i < 26; i++) key[i] = i;
        for (int step = 0; step < 500; step++) {
            aristocrat_move(key, 26);
            int seen[26] = {0};
            for (int i = 0; i < 26; i++) {
                if (key[i] < 0 || key[i] >= 26 || seen[key[i]]) { bad_perm++; break; }
                seen[key[i]] = 1;
            }
        }
    }
    printf("[move invariant] 400 chains x 500 moves: perm-violations=%d\n", bad_perm);
    CHECK(bad_perm == 0, "aristocrat_move broke the permutation %d times", bad_perm);
}

// --- 3. capability floor (registry schedule) ----------------------------------

static void test_capability_floor(void) {
    const char *kws[] = { "KRYPTOS", "CRYPTOGRAM", "ZEBRA" };
    int plen = 120;
    printf("\n[capability floor @ ~%d chars (ACA range), registry schedule]\n", plen);
    double best = 0.0;
    for (int k = 0; k < 3; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, ARISTOCRAT, 0, 0, METHOD_DEFAULT, 0xC0FFEEu + 17 * k, &secs);
        printf("  %-10s : %.1f%%  [%.1fs]\n", kws[k], 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.98, "aristocrat capability floor: best only %.1f%% at ~%d chars", 100.0 * best, plen);
}

// --- 4. length cliff (registry schedule) --------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 60, 100, 200 };
    printf("\n[length cliff: keyword ZEBRA, registry anneal]\n");
    double frac100 = 0.0, frac200 = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("ZEBRA", lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, ARISTOCRAT, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 1) frac100 = frac;
        if (li == 2) frac200 = frac;
    }
    CHECK(frac200 >= 0.98, "aristocrat 200ch recovered only %.3f", frac200);
    CHECK(frac100 >= 0.90, "aristocrat 100ch recovered only %.3f", frac100);
}

// --- 5. multi-keyword sweep ---------------------------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "MONARCHY", "GADGETRY", "CIPHERTEXT" };
    int plen = 130, nk = 4;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars (ACA range), registry schedule]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, ARISTOCRAT, 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL);
        printf("  %-10s : %.1f%%\n", kws[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int plen = 130;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", plen, prepared, cs);
    struct { int method; const char *name; int r, h; } M[] = {
        { METHOD_DEFAULT, "anneal ", 0, 0 },
        { METHOD_SHOTGUN, "shotgun", 60, 200000 },
        { METHOD_PSO,     "pso    ",  3,    400 },
    };
    printf("\n[per-scheme @ ~%d chars, keyword CIPHER (pso bounded, not asserted)]\n", plen);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, ARISTOCRAT, M[m].r, M[m].h, M[m].method, 0x5C8E0u + m, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

// --- 7. Patristocrat smoke (same core, no word divisions) ----------------------

static void test_patristocrat(void) {
    int plen = 130;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("KRYPTOS", plen, prepared, cs);   // spaced planting; the solver strips spaces anyway
    double secs;
    double frac = solve_and_frac(cs, prepared, n, PATRISTOCRAT, 0, 0, METHOD_DEFAULT, 0xBEEFu, &secs);
    printf("\n[patristocrat @ ~%d chars, keyword KRYPTOS]\n  %.1f%%  [%.1fs]\n", plen, 100.0 * frac, secs);
    CHECK(frac > 0.95, "patristocrat recovery only %.1f%%", 100.0 * frac);
}

int main(void) {
    init_alphabet(NULL);                        // full 26-letter alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    // A free 26-letter substitution rides the log-probability quadgram fitness (-logprob).
    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_move_invariant();
    test_capability_floor();
    test_length_cliff();
    test_multi_keyword();
    test_per_scheme();
    test_patristocrat();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
