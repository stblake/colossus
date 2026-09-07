//
//  In-process stress / limits tests for the Seriated Playfair solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Seriated Playfair (ACA) is plain Playfair over a single 5x5 keyed square, with the digraphs
//  the vertical pairs of a two-row seriated layout of period P. There is NO per-column
//  decoupling -- one square enciphers every pair -- so the attack is Playfair's single-grid
//  anneal, with the seriation period SWEPT (one config per P; the n-gram score picks it). It
//  effectively needs the log-probability fitness, so this suite enables g_ngram_logprob and
//  loads quadgrams in that mode. The suite (also the basis for tuning the SearchDefaults
//  schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. a capability floor (recovery %) across several keywords, period pinned;
//    3. a length cliff (recovery vs length), period pinned;
//    4. a multi-keyword sweep (mean/worst recovery), period pinned;
//    5. a BLIND period solve (P swept over a bounded range) -- the reported P must match;
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

static const char *PLAINTEXT =
    "WHENINTHECOURSEOFHUMANEVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICAL"
    "BANDSWHICHHAVECONNECTEDTHEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATE"
    "ANDEQUALSTATIONTOWHICHTHELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHE"
    "OPINIONSOFMANKINDREQUIRESTHATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATION"
    "WEHOLDTHESETRUTHSTOBESELFEVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIR"
    "CREATORWITHCERTAINUNALIENABLERIGHTSTHATAMONGTHESEARELIFELIBERTYANDTHEPURSUITOFHAPPINESS";

// A..Z char -> 0..24 alphabet index, merging J into I (the 25-letter convention).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Seriated Playfair cipher: take the first pt_len letters of PLAINTEXT, build the
// grid from `keyword`, prepare into the seriated period-`period` layout (nulls + padding) and
// encrypt. Fills prepared[] (the expected solution) and cipher_str. Returns the prepared len.
static int plant(const char *keyword, int period, int pt_len, int prepared[], char cipher_str[]) {
    int raw[MAX_CIPHER_LENGTH], n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) raw[n++] = idx;
    }
    int kw[64], kwn = 0;
    for (int i = 0; keyword[i] && kwn < 64; i++) {
        int idx = letter_to_index((unsigned char) keyword[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int grid[PLAYFAIR_GRID];
    playfair_grid_from_keyword(kw, kwn, grid);

    int filler = g_char_to_idx['X'], alt = g_char_to_idx['Q'];
    int plen = seriated_playfair_prepare(raw, n, period, filler, alt, prepared, MAX_CIPHER_LENGTH);

    int cipher[MAX_CIPHER_LENGTH];
    seriated_playfair_encrypt(prepared, plen, grid, period, cipher);
    for (int i = 0; i < plen; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[plen] = '\0';
    return plen;
}

// Solve and return the recovered fraction (vs the planted prepared plaintext). period>0 pins
// P; minc/maxc>0 bound the blind sweep; method overrides the search scheme. Always applies the
// tuned per-type registry schedule.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int period, int minc, int maxc, int method, uint32_t seed,
        int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = SERIATED_PLAYFAIR;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (minc > 0) cfg.min_cols = minc;
    if (maxc > 0) cfg.max_cols = maxc;

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

    if (period_out) *period_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = SERIATED_PLAYFAIR; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "seriated-playfair registry: no entry applied");
    CHECK(cfg.n_restarts == 6 && cfg.n_hill_climbs == 400000,
        "seriated-playfair anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "seriated-playfair anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = SERIATED_PLAYFAIR; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "seriated-playfair registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 20 && cfg.n_hill_climbs == 300000,
        "seriated-playfair shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor (period pinned) --------------------------------------

static void test_capability_floor(void) {
    // Keywords/periods chosen so the quadgram global optimum IS the true grid. (Some grids
    // hit the inherent rare-letter ambiguity of a square attack: the filler X appears only as
    // nulls, so it is weakly constrained, and for an unlucky grid geometry the quadgram optimum
    // swaps X with a neighbour -- a deterministic ~92% ceiling no budget escapes, as in Playfair.
    // These four land 100%; the cliff/sweep below cover the same plaintext at other lengths.)
    const char *kws[] = { "CIPHER", "KRYPTOS", "MONARCHYBDF", "CIPHERKEYWORD" };
    int periods[] =     { 6,        7,         9,             8 };
    int plen = 500;
    printf("\n[capability floor @ ~%d chars, period pinned]\n", plen);
    for (int k = 0; k < (int) (sizeof kws / sizeof kws[0]); k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], periods[k], plen, prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, periods[k], 0, 0, METHOD_DEFAULT,
            0xC0FFEEu + 17 * k, NULL, &secs);
        printf("  %-13s P=%d : %.1f%%  [%.1fs]\n", kws[k], periods[k], 100.0 * frac, secs);
        CHECK(frac > 0.95, "seriated(%s) floor: only %.1f%% at %d chars", kws[k], 100.0 * frac, n);
    }
}

// --- 3. length cliff (period pinned) ------------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 150, 250, 350, 500 };
    const char *kw = "KRYPTOS"; int P = 7;
    printf("\n[length cliff: keyword=%s, P=%d, registry anneal]\n", kw, P);
    double frac_longest = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kw, P, lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == (int) (sizeof lens / sizeof lens[0]) - 1) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.95, "seriated 500ch (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 4. multi-keyword sweep (period pinned) -----------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "MONARCHY", "CIPHER", "KRYPTOS", "PALMERSTON" };
    int periods[] =     { 5,       8,          6,        7,         10 };
    int plen = 500, nk = (int) (sizeof kws / sizeof kws[0]);
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars, period pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], periods[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, periods[k], 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  %-11s P=%2d : %.1f%%\n", kws[k], periods[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 5. blind period solve (P swept over a bounded range) ---------------------

static void test_blind_period(void) {
    const char *kw = "KRYPTOS"; int P = 7, plen = 500;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant(kw, P, plen, prepared, cs);
    double secs; int pout;
    // Bound the sweep (5..9) so the multi-config full-grid anneal stays tractable in CI.
    double frac = solve_and_frac(cs, prepared, n, 0, 5, 9, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind P in 5..9, true keyword=%s P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        kw, P, pout, 100.0 * frac, secs);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
    CHECK(frac > 0.90, "blind-P recovery only %.1f%%", 100.0 * frac);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    const char *kw = "MONARCHYBDF"; int P = 9, plen = 500;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant(kw, P, plen, prepared, cs);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ ~%d chars, keyword=%s, P=%d pinned]\n", plen, kw, P);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, 0, 0, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet("J");                 // 25-letter Seriated Playfair alphabet (J merged into I)
    CHECK(g_alpha == PLAYFAIR_GRID, "alphabet size %d, expected %d", g_alpha, PLAYFAIR_GRID);

    g_ngram_logprob = true;             // Seriated Playfair needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability_floor();
    test_length_cliff();
    test_multi_keyword();
    test_blind_period();
    test_per_scheme();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
