//
//  In-process stress / limits tests for the Nicodemus solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Nicodemus is a COUPLED substitution+transposition cipher solved by annealing the per-block
//  COLUMN ORDER alone and DERIVING the per-column shifts by monogram fit (the -optimalcycle-style
//  decoupling). This suite, also the basis for tuning the SearchDefaults schedule, checks:
//    1. registry validation (apply_cipher_defaults) for all three codes + a non-registry type;
//    2. a capability floor (recovery fraction at ~250+ letters, P/H pinned) per variant;
//    3. a length cliff (recovery vs length, P/H pinned);
//    4. a multi-keyword sweep (mean/worst recovery, P/H pinned);
//    5. blind solves -- P swept (H pinned) and H swept (P pinned) -- validating each sweep axis.
//
//  Nicodemus runs on the full 26-letter alphabet and needs the log-probability fitness. Run from
//  the source directory (loads ngram_data/english/english_quadgrams.txt).
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "nicodemus.h"            // primitives, for planting ciphers
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
    "THEUNANIMOUSDECLARATIONOFTHETHIRTEENUNITEDSTATESOFAMERICAWHENINTHECOURSEOFHUMAN"
    "EVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICALBANDSWHICHHAVECONNECTED"
    "THEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATEANDEQUALSTATIONTOWHICH"
    "THELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHEOPINIONSOFMANKINDREQUIRES"
    "THATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATIONWEHOLDTHESETRUTHSTOBESELF"
    "EVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIRCREATORWITHCERTAINUNALIENABLE";

// Map a cipher-type code to the primitive's NICO_* substitution variant (mirrors the
// solver's nico_variant_of). Used so the planted cipher matches the type being solved.
static int nico_var(int type) {
    if (type == NICODEMUS_VARIANT)  return NICO_VARIANT;
    if (type == NICODEMUS_BEAUFORT) return NICO_BEAU;
    return NICO_VIG;
}

// Plant a Nicodemus cipher of cipher-type `type`: first pt_len plaintext letters under
// `keyword`, block height H. Fills prepared[] (expected solution) and cipher_str (bare letters).
static int plant(int type, const char *keyword, int H, int pt_len,
                 int prepared[], char cipher_str[]) {
    int variant = nico_var(type);
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int kw[MAX_COLS], P = 0;
    for (int i = 0; keyword[i] && P < MAX_COLS; i++)
        if (keyword[i] >= 'A' && keyword[i] <= 'Z') kw[P++] = keyword[i] - 'A';
    int order[MAX_COLS], shifts[MAX_COLS];
    nicodemus_key_from_keyword(kw, P, order, shifts);
    static int cipher[MAX_CIPHER_LENGTH];
    nicodemus_encrypt(prepared, n, P, H, order, shifts, variant, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. period>0 pins P; block_h>0 pins H.
static double solve_and_frac(int type, const char *cipher_str, const int prepared[], int plen,
        int period, int block_h, uint32_t seed, int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
    cfg.ngram_size = NGRAM_SIZE;        // log-prob table is selected by g_ngram_logprob (set in main)
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0)  { cfg.period_present = true; cfg.period = period; }
    if (block_h > 0) { cfg.block_height = block_h; }

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
    int codes[3] = { NICODEMUS, NICODEMUS_VARIANT, NICODEMUS_BEAUFORT };
    const char *names[3] = { "nicodemus", "nicodemus-variant", "nicodemus-beaufort" };
    for (int i = 0; i < 3; i++) {
        init_config(&cfg); cfg.cipher_type = codes[i]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "%s registry: no entry applied", names[i]);
        CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 20000,
            "%s anneal defaults wrong: %dx%d", names[i], cfg.n_restarts, cfg.n_hill_climbs);
    }
    // A non-registry type must be left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. capability floor (P/H pinned) -----------------------------------------

static void test_capability_floor(void) {
    struct { int variant; const char *name; const char *kw; } V[] = {
        { NICODEMUS,          "vigenere", "SPHINX"  },
        { NICODEMUS_VARIANT,  "variant",  "PLANET"  },
        { NICODEMUS_BEAUFORT, "beaufort", "FALCON"  },
    };
    int H = 5, plen = 300;
    printf("\n[capability floor @ %d chars, H=%d, P/H pinned]\n", plen, H);
    for (int v = 0; v < 3; v++) {
        int P = (int) strlen(V[v].kw);
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(V[v].variant, V[v].kw, H, plen, prepared, cs);
        double secs; int pout;
        double frac = solve_and_frac(V[v].variant, cs, prepared, n, P, H, 0xC0FFEEu + v, &pout, &secs);
        printf("  %-9s kw=%-7s P=%d : %.1f%%  [%.1fs]\n", V[v].name, V[v].kw, P, 100.0 * frac, secs);
        CHECK(frac > 0.95, "%s capability floor: only %.1f%% at %d chars", V[v].name, 100.0 * frac, n);
    }
}

// --- 3. length cliff (Vigenere, P/H pinned) -----------------------------------

static void test_length_cliff(void) {
    int lens[] = { 120, 160, 200, 260, 340, 460 };
    int H = 5; const char *kw = "SPHINX"; int P = 6;
    printf("\n[length cliff: Vigenere kw=%s P=%d H=%d, P/H pinned]\n", kw, P, H);
    double best = 0.0;
    for (int li = 0; li < 6; li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(NICODEMUS, kw, H, lens[li], prepared, cs);
        double secs; double frac = solve_and_frac(NICODEMUS, cs, prepared, n, P, H, 0x5EEDu + li, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 4. multi-keyword sweep (P/H pinned) --------------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "SPHINX", "JUPITER", "MONARCH", "DIAGRAM", "VICTORY", "TANGENT" };
    int H = 5, plen = 320;
    double sum = 0, worst = 1.0; int nk = 6;
    printf("\n[multi-keyword sweep @ %d chars, Vigenere, H=%d, P/H pinned]\n", plen, H);
    for (int k = 0; k < nk; k++) {
        int P = (int) strlen(kws[k]);
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
        int n = plant(NICODEMUS, kws[k], H, plen, prepared, cs);
        double frac = solve_and_frac(NICODEMUS, cs, prepared, n, P, H, 0xABCDu + k, NULL, NULL);
        printf("  kw=%-8s P=%d : %.1f%%\n", kws[k], P, 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 5. blind solves (one axis swept at a time) -------------------------------

static void test_blind_period(void) {
    // H pinned (the ACA-standard 5), P swept: the solver must report the true P and recover.
    const char *kw = "JUPITER"; int P = 7, H = 5, plen = 320;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
    int n = plant(NICODEMUS, kw, H, plen, prepared, cs);
    double secs; int pout;
    double frac = solve_and_frac(NICODEMUS, cs, prepared, n, 0, H, 0xB11Du, &pout, &secs);
    printf("\n[blind P (H=%d pinned), true P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        H, P, pout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-P recovery only %.1f%%", 100.0 * frac);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
}

static void test_blind_blockheight(void) {
    // P pinned, H swept: a wrong block height de-transposes to garbage, so a high recovery
    // fraction confirms the true H was found.
    const char *kw = "MONARCH"; int P = 7, H = 4, plen = 320;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH];
    int n = plant(NICODEMUS, kw, H, plen, prepared, cs);
    double secs;
    double frac = solve_and_frac(NICODEMUS, cs, prepared, n, P, 0, 0xB10Cu, NULL, &secs);
    printf("[blind H (P=%d pinned), true H=%d]: %.1f%%  [%.1fs]\n", P, H, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-H recovery only %.1f%%", 100.0 * frac);
}

int main(void) {
    g_ngram_logprob = true;                       // AZDecrypt-style penalising n-gram fitness
    init_alphabet(NULL);                          // full 26-letter alphabet
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("cannot load %s\n", NGRAM_FILE); return 1; }

    test_registry();
    test_capability_floor();
    test_length_cliff();
    test_multi_keyword();
    test_blind_period();
    test_blind_blockheight();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
