//
//  In-process stress / limits tests for the Progressive Key solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  The Progressive Key cipher is a periodic base cipher (Vigenere / Variant / Beaufort) under a
//  letter keyword, composed with a per-group constant key drift (the progression index). The
//  whole key is P per-column base shifts PLUS the progression 0..25; IoC period estimation fails
//  through the drift, so the solver brute-forces the period and enumerates the progression (one
//  engine config per (P, prog) pair), with a per-column monogram-fit warm start decoupling each
//  column's shift on the de-progressed ciphertext. This suite, also the basis for tuning the
//  SearchDefaults schedule, checks:
//    1. registry validation (apply_cipher_defaults) for all three codes + a non-registry type;
//    2. a capability floor over the ACA ~150-letter band, per base (P/prog pinned);
//    3. a length cliff (recovery vs length, P/prog pinned);
//    4. a multi-keyword sweep (mean/worst recovery, P/prog pinned);
//    5. a BLIND period solve (P swept, prog pinned) -- the reported period must match the true one;
//    6. a BLIND progression solve (prog swept, P pinned) -- the reported prog must match the true one;
//    7. per-scheme calibration: the same cipher under -method anneal / shotgun / pso (recovery + time).
//
//  Progressive Key runs on the full 26-letter alphabet and -- like the rest of the Vigenere family
//  -- rides the reward-only quadgram table (no -logprob). Run from the source directory.
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
    "WEHOLDTHESETRUTHSTOBESELFEVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIR";

static const char *BASE_NAME[3] = { "vig", "var", "beau" };
static int base_type(int base) {
    return base == PROGKEY_BASE_VAR ? PROGKEY_VAR
         : base == PROGKEY_BASE_BEAU ? PROGKEY_BEAU : PROGKEY;
}

// Keyword string -> per-column base shifts (indices 0..25); returns the period P.
static int kw_to_key(const char *kw, int key[]) {
    int P = 0;
    for (int i = 0; kw[i]; i++) key[P++] = toupper((unsigned char) kw[i]) - 'A';
    return P;
}

// Plant a Progressive Key cipher: the first pt_len plaintext letters under (kw, prog, base).
// The cipher is positional, so no padding is needed. Fills prepared[] + cipher_str.
static int plant(const char *kw, int pt_len, int prog, int base,
                 int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN];
    int P = kw_to_key(kw, key);
    if (P_out) *P_out = P;
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    static int cipher[MAX_CIPHER_LENGTH];
    progkey_encrypt(cipher, prepared, n, key, P, prog, base);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. period>0 pins P; progression>=0 pins prog; method
// overrides the search scheme. Reports the recovered period / progression when requested.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int cipher_type, int period, int progression, int method, uint32_t seed,
        int *period_out, int *prog_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = cipher_type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (progression >= 0) { cfg.progression_present = true; cfg.progression = progression; }

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
    if (prog_out)   *prog_out   = res.solved ? res.progression : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    int codes[3] = { PROGKEY, PROGKEY_VAR, PROGKEY_BEAU };
    for (int i = 0; i < 3; i++) {
        ColossusConfig cfg;
        init_config(&cfg); cfg.cipher_type = codes[i]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "progkey %s registry: no entry applied", BASE_NAME[i]);
        CHECK(cfg.n_restarts == 3 && cfg.n_hill_climbs == 2500,
            "progkey %s anneal defaults wrong: %dx%d", BASE_NAME[i], cfg.n_restarts, cfg.n_hill_climbs);
    }
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. capability floor per base (P/prog pinned) -----------------------------

static void test_capability_floor(void) {
    int bases[3] = { PROGKEY_BASE_VIG, PROGKEY_BASE_VAR, PROGKEY_BASE_BEAU };
    const char *kws[] = { "GRAPEFRUIT", "KRYPTOS", "PALMERSTON" };
    int progs[]       = { 1,            2,         3 };
    int plen = 150;
    printf("\n[capability floor @ %d chars, P/prog pinned, per base]\n", plen);
    for (int b = 0; b < 3; b++) {
        for (int k = 0; k < 3; k++) {
            int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
            int n = plant(kws[k], plen, progs[k], bases[b], prepared, cs, &P);
            double secs;
            double frac = solve_and_frac(cs, prepared, n, base_type(bases[b]),
                P, progs[k], METHOD_DEFAULT, 0xC0FFEEu + 7u * b + k, NULL, NULL, &secs);
            printf("  %-4s %-11s P=%2d prog=%d : %.1f%%  [%.1fs]\n",
                BASE_NAME[b], kws[k], P, progs[k], 100.0 * frac, secs);
            CHECK(frac > 0.95, "progkey %s/%s capability floor: only %.1f%% at %d chars",
                BASE_NAME[b], kws[k], 100.0 * frac, n);
        }
    }
}

// --- 3. length cliff (P/prog pinned) ------------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 24, 32, 48, 72, 120 };
    const char *kw = "GRAPEFRUIT"; int prog = 1;
    printf("\n[length cliff: keyword=%s prog=%d (Vigenere), P/prog pinned]\n", kw, prog);
    double best = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kw, lens[li], prog, PROGKEY_BASE_VIG, prepared, cs, &P);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, PROGKEY, P, prog,
            METHOD_DEFAULT, 0x5EEDu + li, NULL, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 4. multi-keyword sweep (P/prog pinned) -----------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "CIPHER", "KRYPTOS", "MONARCHY", "GRAPEFRUIT", "PALMERSTON" };
    int plen = 150, nk = (int) (sizeof kws / sizeof kws[0]);
    int prog = 4;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ %d chars, prog=%d (Vigenere), P/prog pinned]\n", plen, prog);
    for (int k = 0; k < nk; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant(kws[k], plen, prog, PROGKEY_BASE_VIG, prepared, cs, &P);
        double frac = solve_and_frac(cs, prepared, n, PROGKEY, P, prog,
            METHOD_DEFAULT, 0xABCDu + k, NULL, NULL, NULL);
        printf("  %-11s P=%2d : %.1f%%\n", kws[k], P, 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 5. blind period solve (P swept, prog pinned) -----------------------------

static void test_blind_period(void) {
    const char *kw = "KRYPTOS"; int plen = 180, prog = 3;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, plen, prog, PROGKEY_BASE_VIG, prepared, cs, &P);
    double secs; int pout;
    double frac = solve_and_frac(cs, prepared, n, PROGKEY, 0, prog,
        METHOD_DEFAULT, 0xB11Du, &pout, NULL, &secs);
    printf("\n[blind P, true keyword=%s P=%d prog=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        kw, P, prog, pout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-P recovery only %.1f%%", 100.0 * frac);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
}

// --- 6. blind progression solve (prog swept, P pinned) ------------------------

static void test_blind_progression(void) {
    const char *kw = "GRAPEFRUIT"; int plen = 180, prog = 5;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, plen, prog, PROGKEY_BASE_VIG, prepared, cs, &P);
    double secs; int progout;
    double frac = solve_and_frac(cs, prepared, n, PROGKEY, P, -1,
        METHOD_DEFAULT, 0xB12Eu, NULL, &progout, &secs);
    printf("\n[blind prog, true keyword=%s P=%d prog=%d]: reported prog=%d, %.1f%%  [%.1fs]\n",
        kw, P, prog, progout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-prog recovery only %.1f%%", 100.0 * frac);
    CHECK(progout == prog, "blind-prog reported prog=%d (true %d)", progout, prog);
}

// --- 7. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    const char *kw = "KRYPTOS"; int plen = 150, prog = 2;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant(kw, plen, prog, PROGKEY_BASE_VIG, prepared, cs, &P);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d chars, keyword=%s prog=%d, P/prog pinned]\n", plen, kw, prog);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, PROGKEY, P, prog,
            M[m].method, 0x5C8E0u + m, NULL, NULL, &secs);
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
    test_blind_progression();
    test_per_scheme();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
