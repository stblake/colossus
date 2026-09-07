//
//  In-process stress / limits tests for the Portax solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Portax (ACA "periodic digraphic Porta") enciphers vertical pairs over a Porta slide; the key
//  is P per-column Porta shifts (0..12), one engine config per swept period P, with a per-column
//  monogram-fit warm start decoupling each column's shift. This suite, also the basis for tuning
//  the SearchDefaults schedule, checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. a capability floor (recovery %, several keywords, period pinned);
//    3. a length cliff (recovery vs length, period pinned);
//    4. a multi-keyword sweep (mean/worst recovery, period pinned);
//    5. a BLIND period solve (P swept) -- the reported period must match the true one;
//    6. per-scheme calibration: the same cipher solved under -method anneal / shotgun / pso,
//       reporting recovery + time for each (the data the schedule is tuned against).
//
//  Portax runs on the full 26-letter alphabet and -- like the other Porta-family ciphers -- rides
//  the reward-only quadgram table (no -logprob). Run from the source directory (loads the table).
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "portax.h"               // primitives, for planting ciphers
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

// Keyword string -> per-column key letters (indices 0..25); returns the period P.
static int kw_to_key(const char *kw, int key[]) {
    int P = 0;
    for (int i = 0; kw[i]; i++) key[P++] = toupper((unsigned char) kw[i]) - 'A';
    return P;
}

// Plant a Portax cipher: the first pt_len plaintext letters under `kw`, padded with X to a full
// number of row-pairs (a multiple of 2*P). Fills prepared[] (expected solution) and cipher_str.
static int plant(const char *kw, int pt_len, int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN];
    int P = kw_to_key(kw, key);
    if (P_out) *P_out = P;
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int block = 2 * P;
    while (n % block != 0) prepared[n++] = 'X' - 'A';       // fill an even number of rows
    static int cipher[MAX_CIPHER_LENGTH];
    portax_encrypt(cipher, prepared, n, key, P);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. period>0 pins P; method overrides the search scheme.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int period, int method, uint32_t seed, int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = PORTAX;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0) { cfg.period_present = true; cfg.period = period; }

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
    init_config(&cfg); cfg.cipher_type = PORTAX; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "portax registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 20000,
        "portax anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. capability floor (period pinned) --------------------------------------

static void test_capability_floor(void) {
    const char *kws[] = { "PORTAX", "KRYPTOS", "CIPHER", "MONARCHY" };
    int plen = 180;
    printf("\n[capability floor @ %d chars, period pinned]\n", plen);
    for (int k = 0; k < (int) (sizeof kws / sizeof kws[0]); k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kws[k], plen, prepared, cs, &P);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, METHOD_DEFAULT, 0xC0FFEEu + k, NULL, &secs);
        printf("  %-9s P=%d : %.1f%%  [%.1fs]\n", kws[k], P, 100.0 * frac, secs);
        CHECK(frac > 0.95, "keyword %s capability floor: only %.1f%% at %d chars",
            kws[k], 100.0 * frac, n);
    }
}

// --- 3. length cliff (period pinned) ------------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 60, 80, 100, 140, 200 };
    const char *kw = "KRYPTOS";
    printf("\n[length cliff: keyword=%s, period pinned]\n", kw);
    double best = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kw, lens[li], prepared, cs, &P);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 4. multi-keyword sweep (period pinned) -----------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "PORTAX", "CIPHER", "KRYPTOS", "MONARCHY", "PALMERSTON" };
    int plen = 200, nk = (int) (sizeof kws / sizeof kws[0]);
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ %d chars, period pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kws[k], plen, prepared, cs, &P);
        double frac = solve_and_frac(cs, prepared, n, P, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  %-11s P=%2d : %.1f%%\n", kws[k], P, 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 5. blind period solve (P swept) ------------------------------------------

static void test_blind_period(void) {
    const char *kw = "KRYPTOS"; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, plen, prepared, cs, &P);
    double secs; int pout;
    double frac = solve_and_frac(cs, prepared, n, 0, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind P, true keyword=%s P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        kw, P, pout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-P recovery only %.1f%%", 100.0 * frac);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    const char *kw = "PORTAX"; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, plen, prepared, cs, &P);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d chars, keyword=%s, period pinned]\n", plen, kw);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet(NULL);                          // full 26-letter alphabet
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("cannot load %s\n", NGRAM_FILE); return 1; }

    test_registry();
    test_capability_floor();
    test_length_cliff();
    test_multi_keyword();
    test_blind_period();
    test_per_scheme();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
