//
//  In-process stress / limits tests for the Slidefair solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Slidefair (ACA "periodic digraphic Vigenere/Variant/Beaufort") enciphers consecutive digraphs
//  over a two-row slide; the key is P per-column key letters (0..25), one engine config per swept
//  period P, with a per-column monogram-fit warm start decoupling each column's key. This suite,
//  also the basis for tuning the SearchDefaults schedule, checks:
//    1. registry validation (apply_cipher_defaults) for all three codes + a non-registry type;
//    2. a capability floor (recovery %) across all three variants, period pinned;
//    3. a length cliff (recovery vs length), period pinned;
//    4. a multi-keyword sweep (mean/worst recovery), period pinned;
//    5. a BLIND period solve (P swept) -- the reported period must match the true one;
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Slidefair runs on the full 26-letter alphabet and -- like the rest of the Vigenere family --
//  rides the reward-only quadgram table (no -logprob). Run from the source directory.
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "slidefair.h"            // primitives, for planting ciphers
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

// Plant a Slidefair cipher: the first pt_len plaintext letters under `kw`/`type`, padded with X to
// an even length (a whole number of digraphs). Fills prepared[] (expected solution) and cipher_str.
static int plant(const char *kw, int type, int pt_len, int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN];
    int P = kw_to_key(kw, key);
    if (P_out) *P_out = P;
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    if (n % 2 != 0) prepared[n++] = 'X' - 'A';              // even number of letters (digraphs)
    static int cipher[MAX_CIPHER_LENGTH];
    slidefair_encrypt(cipher, prepared, n, key, P, type);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. period>0 pins P; method overrides the search scheme.
static double solve_and_frac(const char *cipher_str, int type, const int prepared[], int plen,
        int period, int method, uint32_t seed, int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
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

static const char *vname(int type) {
    return type == SLIDEFAIR_VAR ? "var" : type == SLIDEFAIR_BEAU ? "beau" : "vig ";
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    int types[3] = { SLIDEFAIR, SLIDEFAIR_VAR, SLIDEFAIR_BEAU };
    for (int t = 0; t < 3; t++) {
        ColossusConfig cfg;
        init_config(&cfg); cfg.cipher_type = types[t]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "slidefair(%s) registry: no entry applied", vname(types[t]));
        CHECK(cfg.n_restarts == 8 && cfg.n_hill_climbs == 10000,
            "slidefair(%s) anneal defaults wrong: %dx%d", vname(types[t]), cfg.n_restarts, cfg.n_hill_climbs);
    }
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. capability floor (all three variants, period pinned) ------------------

static void test_capability_floor(void) {
    int types[3] = { SLIDEFAIR, SLIDEFAIR_VAR, SLIDEFAIR_BEAU };
    const char *kws[] = { "SLIDE", "KRYPTOS", "CIPHER", "MONARCHY" };
    int plen = 180;
    printf("\n[capability floor @ %d chars, period pinned]\n", plen);
    for (int t = 0; t < 3; t++) {
        for (int k = 0; k < (int) (sizeof kws / sizeof kws[0]); k++) {
            int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
            int n = plant(kws[k], types[t], plen, prepared, cs, &P);
            double secs;
            double frac = solve_and_frac(cs, types[t], prepared, n, P, METHOD_DEFAULT,
                0xC0FFEEu + 100*t + k, NULL, &secs);
            printf("  %s %-9s P=%d : %.1f%%  [%.1fs]\n", vname(types[t]), kws[k], P, 100.0 * frac, secs);
            CHECK(frac > 0.95, "slidefair(%s) keyword %s floor: only %.1f%% at %d chars",
                vname(types[t]), kws[k], 100.0 * frac, n);
        }
    }
}

// --- 3. length cliff (period pinned) ------------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 50, 70, 100, 140, 200 };
    const char *kw = "KRYPTOS";
    printf("\n[length cliff: keyword=%s, Vigenere, period pinned]\n", kw);
    double best = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kw, SLIDEFAIR, lens[li], prepared, cs, &P);
        double secs;
        double frac = solve_and_frac(cs, SLIDEFAIR, prepared, n, P, METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 4. multi-keyword sweep (period pinned) -----------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "SLIDE", "CIPHER", "KRYPTOS", "MONARCHY", "PALMERSTON" };
    int plen = 200, nk = (int) (sizeof kws / sizeof kws[0]);
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ %d chars, Vigenere, period pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kws[k], SLIDEFAIR, plen, prepared, cs, &P);
        double frac = solve_and_frac(cs, SLIDEFAIR, prepared, n, P, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
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
    int n = plant(kw, SLIDEFAIR, plen, prepared, cs, &P);
    double secs; int pout;
    double frac = solve_and_frac(cs, SLIDEFAIR, prepared, n, 0, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind P, true keyword=%s P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        kw, P, pout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-P recovery only %.1f%%", 100.0 * frac);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    const char *kw = "SLIDE"; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, SLIDEFAIR, plen, prepared, cs, &P);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d chars, keyword=%s, period pinned]\n", plen, kw);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, SLIDEFAIR, prepared, n, P, M[m].method, 0x5C8E0u + m, NULL, &secs);
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
