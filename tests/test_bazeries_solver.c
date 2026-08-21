//
//  In-process stress / limits tests for the Bazeries solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Bazeries is keyed by ONE number N < 10^6, climbed as its decimal digits (one engine config
//  per digit count D), with a square-quality monogram reward decoupling the square from the
//  digit-grouped transposition. This suite, also the basis for tuning the SearchDefaults
//  schedule, checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. a capability floor over the ACA 150-250-letter band, several numbers, D pinned;
//    3. a length cliff (recovery vs length, D pinned);
//    4. a multi-number sweep (mean/worst recovery, D pinned);
//    5. a BLIND digit-count solve (D swept) -- the reported digit count must match the true one;
//    6. per-scheme calibration: the same cipher solved under -method anneal / shotgun / pso,
//       reporting recovery + time for each (the data the schedule is tuned against).
//
//  Bazeries runs on the 25-letter J->I alphabet and needs the log-probability fitness. Run from
//  the source directory (loads ngram_data/english/english_quadgrams.txt).
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "bazeries.h"             // primitives, for planting ciphers
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
    "WEHOLDTHESETRUTHSTOBESELFEVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIR";

// Count the decimal digits of N (== the digit-count config the solver reports as cycleword_len).
static int ndigits_of(long N) { int d = 0; if (N == 0) return 1; while (N > 0) { d++; N /= 10; } return d; }

// Plant a Bazeries cipher: first pt_len plaintext letters (J->I folded) under key number N.
// Fills prepared[] (expected solution, J->I indices) and cipher_str (bare letters).
static int plant(long N, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c == 'J') c = 'I';
        if (c >= 'A' && c <= 'Z') {
            int idx = g_char_to_idx[c];
            if (idx >= 0) prepared[n++] = idx;
        }
    }
    static int cipher[MAX_CIPHER_LENGTH];
    bazeries_encrypt(prepared, n, N, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. digits>0 pins the digit count D; method overrides
// the search scheme (METHOD_DEFAULT keeps the registry/anneal default).
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int digits, int method, uint32_t seed, int *digits_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = BAZERIES;
    cfg.ngram_size = NGRAM_SIZE;        // log-prob table selected by g_ngram_logprob (set in main)
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (digits > 0) { cfg.period_present = true; cfg.period = digits; }

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

    if (digits_out) *digits_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = BAZERIES; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "bazeries registry: no entry applied");
    CHECK(cfg.n_restarts == 40 && cfg.n_hill_climbs == 20000,
        "bazeries anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. capability floor (D pinned) -------------------------------------------

static void test_capability_floor(void) {
    long nums[] = { 3752, 81257, 246813, 50407 };   // 4..6 digits, incl. internal zeros
    int plen = 200;
    printf("\n[capability floor @ %d chars, D pinned]\n", plen);
    for (int k = 0; k < (int) (sizeof nums / sizeof nums[0]); k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(nums[k], plen, prepared, cs);
        double secs; int dout;
        double frac = solve_and_frac(cs, prepared, n, ndigits_of(nums[k]),
            METHOD_DEFAULT, 0xC0FFEEu + k, &dout, &secs);
        printf("  N=%-7ld D=%d : %.1f%%  [%.1fs]\n", nums[k], ndigits_of(nums[k]), 100.0 * frac, secs);
        CHECK(frac > 0.95, "N=%ld capability floor: only %.1f%% at %d chars", nums[k], 100.0 * frac, n);
    }
}

// --- 3. length cliff (D pinned) -----------------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 100, 130, 160, 200, 250 };
    long N = 81257;
    printf("\n[length cliff: N=%ld D=%d, D pinned]\n", N, ndigits_of(N));
    double best = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(N, lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, ndigits_of(N), METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 4. multi-number sweep (D pinned) -----------------------------------------

static void test_multi_number(void) {
    long nums[] = { 837, 3752, 81257, 50407, 246813, 909090 };
    int plen = 200, nk = (int) (sizeof nums / sizeof nums[0]);
    double sum = 0, worst = 1.0;
    printf("\n[multi-number sweep @ %d chars, D pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(nums[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, ndigits_of(nums[k]), METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  N=%-7ld D=%d : %.1f%%\n", nums[k], ndigits_of(nums[k]), 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-number mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 5. blind digit-count solve (D swept) -------------------------------------

static void test_blind_digits(void) {
    long N = 81257; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
    int n = plant(N, plen, prepared, cs);
    double secs; int dout;
    double frac = solve_and_frac(cs, prepared, n, 0, METHOD_DEFAULT, 0xB11Du, &dout, &secs);
    printf("\n[blind D, true N=%ld D=%d]: reported D=%d, %.1f%%  [%.1fs]\n",
        N, ndigits_of(N), dout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-D recovery only %.1f%%", 100.0 * frac);
    CHECK(dout == ndigits_of(N), "blind-D reported D=%d (true %d)", dout, ndigits_of(N));
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    long N = 3752; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
    int n = plant(N, plen, prepared, cs);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d chars, N=%ld, D pinned]\n", plen, N);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, ndigits_of(N), M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    g_ngram_logprob = true;                       // AZDecrypt-style penalising n-gram fitness
    init_alphabet("J");                           // 25-letter J->I alphabet (as the binary forces)
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("cannot load %s\n", NGRAM_FILE); return 1; }

    test_registry();
    test_capability_floor();
    test_length_cliff();
    test_multi_number();
    test_blind_digits();
    test_per_scheme();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
