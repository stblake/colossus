//
//  In-process stress / limits tests for the ADFGX / ADFGVX solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each stochastic solve
//  deterministic.
//
//  ADFGVX is a COUPLED keyed-square + keyed-columnar cipher, the hardest of the
//  polygraphic family. This suite is deliberately thorough:
//    1. registry validation (apply_cipher_defaults) for ADFGX and ADFGVX, plus a
//       non-registry type left untouched (regression safety);
//    2. ADFGX capability FLOOR with the column count K pinned (isolates the joint
//       square+order recovery), recovered ~100% and the reported K correct;
//    3. ADFGX BLIND K-selection -- K swept over a range, the solver must still recover
//       and report the true K (this exercises the structural IoC reward that decouples
//       the column-order search from the square);
//    4. ADFGX length cliff -- recovery vs ciphertext length printed (K pinned);
//    5. ADFGX multi-keyword sweep -- mean / worst recovery over several square keywords;
//    6. ADFGVX (6x6, 36-symbol) capability floor -- the side-generic / digit-alphabet
//       path recovered ~100% end to end.
//
//  ADFGX (25 letters) and ADFGVX (36 symbols) need different alphabets AND different
//  n-gram-table packing bases, so the two are run as separate phases, each re-initing the
//  alphabet and reloading the quadgram table. Both need the log-probability fitness, so
//  g_ngram_logprob is enabled (in this mode the n-gram scale is 1 regardless of g_alpha,
//  so switching alphabets mid-process is safe). Run from the source directory.
//

#include "colossus.h"
#include "engine.h"        // apply_cipher_defaults
#include "scoring.h"       // load_ngrams

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A long chunk of natural English (US Declaration of Independence, opening), letters only.
static const char *PLAINTEXT =
    "THEUNANIMOUSDECLARATIONOFTHETHIRTEENUNITEDSTATESOFAMERICAWHENINTHECOURSEOFHUMAN"
    "EVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICALBANDSWHICHHAVECONNECTED"
    "THEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATEANDEQUALSTATIONTOWHICH"
    "THELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHEOPINIONSOFMANKINDREQUIRES"
    "THATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATIONWEHOLDTHESETRUTHSTOBESELF"
    "EVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIRCREATORWITHCERTAINUNALIENABLE"
    "RIGHTSTHATAMONGTHESEARELIFELIBERTYANDTHEPURSUITOFHAPPINESS";

// A..Z (and digits for ADFGVX) char -> alphabet index. ADFGX (side 5) merges J into I.
static int sym_to_index(int c, int side) {
    c = toupper(c);
    if (side == ADFGX_SIDE && c == 'J') c = 'I';
    if (c >= 'A' && c <= 'Z') return g_char_to_idx[c];
    if (side == ADFGVX_SIDE && c >= '0' && c <= '9') return g_char_to_idx[c];
    return -1;
}

// Columnar read order from a transposition keyword (columns in keyword-alphabetical
// order, ties left to right). Returns K.
static int build_order(const char *kw, int side, int order[]) {
    int idx[MAX_COLS], K = 0;
    for (int i = 0; kw[i] && K < MAX_COLS; i++) {
        int v = sym_to_index((unsigned char) kw[i], side);
        if (v >= 0) idx[K++] = v;
    }
    char used[MAX_COLS];
    for (int c = 0; c < K; c++) used[c] = 0;
    for (int j = 0; j < K; j++) {
        int best = -1;
        for (int c = 0; c < K; c++)
            if (!used[c] && (best < 0 || idx[c] < idx[best])) best = c;
        used[best] = 1; order[j] = best;
    }
    return K;
}

// Plant an ADFGX/ADFGVX cipher: first pt_len plaintext symbols under the keyed square
// (from sqkw) and the columnar order (from trkw). Fills prepared[] (the expected
// solution, length n) and cipher_str (the label ciphertext, length 2n). Returns n.
static int plant(const char *sqkw, const char *trkw, int side, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = sym_to_index((unsigned char) PLAINTEXT[i], side);
        if (idx >= 0) prepared[n++] = idx;
    }
    int kw[64], kwn = 0;
    for (int i = 0; sqkw[i] && kwn < 64; i++) {
        int idx = sym_to_index((unsigned char) sqkw[i], side);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int square[SQUARE_MAX_GRID];
    bifid_grid_from_keyword(kw, kwn, square, g_alpha);

    int order[MAX_COLS];
    int K = build_order(trkw, side, order);

    static int cipher[2 * MAX_CIPHER_LENGTH];
    adfgvx_encrypt(prepared, n, square, side, K, order, COL_READ_TB, cipher);
    const char *labels = adfgvx_labels(side);
    for (int i = 0; i < 2 * n; i++) cipher_str[i] = labels[cipher[i]];
    cipher_str[2 * n] = '\0';
    return n;
}

// Run one solve (K swept over [klo,khi]) and return the recovered fraction; K_out gets
// the reported column count.
static double solve_and_frac(int type, const char *cipher_str, const int prepared[], int plen,
        int klo, int khi, int nr, int nh, double it, double bt, uint32_t seed,
        int *K_out, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    cfg.min_cols = klo; cfg.max_cols = khi;
    cfg.n_restarts = nr; cfg.n_hill_climbs = nh;
    cfg.init_temp = it; cfg.backtracking_probability = bt;

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

    if (K_out) *K_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    init_config(&cfg);
    cfg.cipher_type = ADFGX; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "adfgx registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 600000,
        "adfgx anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "adfgx anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg);
    cfg.cipher_type = ADFGVX; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "adfgvx registry: no entry applied");
    CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 800000,
        "adfgvx anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. ADFGX capability floor (K pinned) -------------------------------------

static void test_capability_adfgx(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant("KRYPTOSABCWORLD", "BERLIN", ADFGX_SIDE, 480, prepared, cipher_str);  // K=6
    int K = -1; double secs = 0.0;
    double frac = solve_and_frac(ADFGX, cipher_str, prepared, plen, 6, 6,
        8, 200000, 0.08, 0.30, 1u, &K, &secs);
    printf("[ADFGX capability  %dch pt, K=6 pinned, 8x200k]  frac=%.3f K=%d  %.1fs\n",
        plen, frac, K, secs);
    CHECK(frac >= 0.97, "ADFGX %dch (K pinned) recovered only %.3f", plen, frac);
    CHECK(K == 6, "ADFGX reported K=%d, planted 6", K);
}

// --- 3. ADFGX blind K-selection (K swept) -------------------------------------

static void test_kselect_adfgx(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant("KRYPTOSABCWORLD", "BERLIN", ADFGX_SIDE, 360, prepared, cipher_str);  // K=6
    int K = -1; double secs = 0.0;
    double frac = solve_and_frac(ADFGX, cipher_str, prepared, plen, 2, 8,
        6, 120000, 0.08, 0.30, 7u, &K, &secs);
    printf("[ADFGX K-select    %dch pt, K swept 2..8, 6x120k] frac=%.3f K=%d  %.1fs\n",
        plen, frac, K, secs);
    CHECK(frac >= 0.95, "ADFGX blind K-sweep recovered only %.3f", frac);
    CHECK(K == 6, "ADFGX blind sweep picked K=%d, planted 6", K);
}

// --- 4. ADFGX length cliff (K pinned) -----------------------------------------

static void test_length_cliff_adfgx(void) {
    int lengths[] = {160, 300, 480};
    printf("ADFGX recovery vs plaintext length (keyword KRYPTOSABCWORLD, K=6 pinned, 8x180k, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("KRYPTOSABCWORLD", "BERLIN", ADFGX_SIDE, lengths[li], prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(ADFGX, cipher_str, prepared, plen, 6, 6,
            8, 180000, 0.08, 0.30, 1u, NULL, &secs);
        printf("    pt~%-4d  frac=%.3f  %.1fs\n", plen, frac, secs);
        if (li == 2) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.97, "ADFGX longest (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 5. ADFGX multi-keyword sweep ---------------------------------------------

static void test_multikeyword_adfgx(void) {
    const char *kws[] = {"KRYPTOSABCWORLD", "PLAYFAIREXAMPLE", "ZEBRASCRAMBLED"};
    double sum = 0.0, worst = 1.0;
    printf("ADFGX multi-keyword (360ch pt, K=6 pinned, 10x200k, seed 2):\n");
    for (int i = 0; i < 3; i++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(kws[i], "BERLIN", ADFGX_SIDE, 360, prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(ADFGX, cipher_str, prepared, plen, 6, 6,
            10, 200000, 0.08, 0.30, 2u, NULL, &secs);
        printf("    keyword %-16s frac=%.3f  %.1fs\n", kws[i], frac, secs);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("    mean=%.3f worst=%.3f\n", sum / 3.0, worst);
    CHECK(worst >= 0.95, "ADFGX multi-keyword worst recovery only %.3f", worst);
}

// --- 6. ADFGVX (6x6, 36-symbol) capability floor ------------------------------

static void test_capability_adfgvx(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant("KRYPTOSABCWORLD", "BERLIN", ADFGVX_SIDE, 360, prepared, cipher_str);  // K=6
    int K = -1; double secs = 0.0;
    double frac = solve_and_frac(ADFGVX, cipher_str, prepared, plen, 6, 6,
        12, 300000, 0.08, 0.30, 3u, &K, &secs);
    printf("[ADFGVX capability %dch pt, K=6 pinned, 12x300k] frac=%.3f K=%d  %.1fs\n",
        plen, frac, K, secs);
    CHECK(frac >= 0.97, "ADFGVX %dch (K pinned) recovered only %.3f", plen, frac);
    CHECK(K == 6, "ADFGVX reported K=%d, planted 6", K);
}

int main(void) {
    g_ngram_logprob = true;             // ADFGVX needs the log-probability fitness

    // ---- ADFGX phase (25-letter alphabet, base-25 quadgrams) ----
    init_alphabet("J");
    CHECK(g_alpha == 25, "ADFGX alphabet size %d, expected 25", g_alpha);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }
    test_registry();
    test_capability_adfgx();
    test_kselect_adfgx();
    test_length_cliff_adfgx();
    test_multikeyword_adfgx();
    free(shared.ngram_data);

    // ---- ADFGVX phase (36-symbol alphabet, base-36 quadgrams) ----
    init_alphabet_adfgvx();
    CHECK(g_alpha == 36, "ADFGVX alphabet size %d, expected 36", g_alpha);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s for ADFGVX phase\n", NGRAM_FILE);
        return 1;
    }
    test_capability_adfgvx();
    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
