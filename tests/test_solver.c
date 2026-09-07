//
//  In-process regression tests for the polyalphabetic optimizer (solve_cipher).
//
//  Framework-free: build with `make test`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be
//  driven directly and its SolveResult inspected -- no stdout scraping. A fixed
//  -seed makes each stochastic solve deterministic.
//
//  Strategy (planted-cipher recovery): encrypt a known English plaintext under a
//  known key, run solve_cipher with a fixed seed and a bounded (restarts x
//  hill-climbs) budget, and assert it recovers the plaintext.
//
//  The budget per case is the effectiveness contract: it is set comfortably
//  above what is needed today, but an optimizer change that makes the climber
//  need *more* iterations to reach the same instance will fail recovery at this
//  fixed budget -- catching effectiveness regressions, not just correctness.
//
//  Run from the source directory (as `make test` does) so the n-gram table is
//  found in the cwd.
//

#include "colossus.h"
#include "scoring.h"          // load_ngrams, crib_score (moved out of colossus.h)
#include "polyalpha_solver.h" // constrain_cycleword

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A chunk of natural English (Declaration of Independence), letters only.
static const char *PLAINTEXT =
    "WHENINTHECOURSEOFHUMANEVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVE"
    "THEPOLITICALBANDSWHICHHAVECONNECTEDTHEMWITHANOTHERANDTOASSUMEAMONGTHE"
    "POWERSOFTHEEARTHTHESEPARATEANDEQUALSTATIONTOWHICHTHELAWSOFNATUREENTITLE";

static double recovered_fraction(SolveResult *res, int P[], int len) {
    int ok = 0;
    for (int i = 0; i < len; i++) if (res->decrypted[i] == P[i]) ok++;
    return (double) ok / (double) len;
}

// One planted-cipher recovery case. pt_kw/ct_kw are keyword strings for the
// keyed alphabets (NULL -> straight). pt_len/ct_len pin the keyword-length
// search (0 -> leave at the cipher-type default). crib_in may be NULL.
static void run_case(const char *name, int cipher_type, int variant,
    const char *pt_kw, const char *ct_kw, int *planted_cw, int cwlen,
    int pt_len, int ct_len, const char *crib_in,
    int n_restarts, int n_hillclimbs, uint32_t seed, double min_fraction) {

    int len = (int) strlen(PLAINTEXT);
    int P[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH];
    int ptkw[ALPHABET_SIZE], ctkw[ALPHABET_SIZE];
    char cipher_str[MAX_CIPHER_LENGTH];

    ord((char *) PLAINTEXT, P);

    if (pt_kw) make_keyed_alphabet((char *) pt_kw, ptkw); else straight_alphabet(ptkw, ALPHABET_SIZE);
    if (ct_kw) make_keyed_alphabet((char *) ct_kw, ctkw); else straight_alphabet(ctkw, ALPHABET_SIZE);

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = cipher_type;
    cfg.variant = variant;

    bool is_autokey = (cipher_type >= AUTOKEY_0 && cipher_type <= AUTOKEY_PORTA);

    if (cipher_type == VIGENERE)      vigenere_encrypt(C, P, len, planted_cw, cwlen, variant);
    else if (cipher_type == BEAUFORT) beaufort_encrypt(C, P, len, planted_cw, cwlen);
    else if (cipher_type == PORTA)    porta_encrypt(C, P, len, planted_cw, cwlen);
    else if (is_autokey)              autokey_encrypt(&cfg, C, P, len, ptkw, ctkw, planted_cw, cwlen);
    else                              quagmire_encrypt(C, P, len, ptkw, ctkw, planted_cw, cwlen, variant);

    for (int i = 0; i < len; i++) cipher_str[i] = C[i] + 'A';
    cipher_str[len] = '\0';

    cfg.ngram_size = NGRAM_SIZE;
    cfg.n_restarts = n_restarts;
    cfg.n_hill_climbs = n_hillclimbs;
    cfg.cycleword_len_present = true;     // pin the period (bypass IoC estimation)
    cfg.cycleword_len = cwlen;
    cfg.max_cycleword_len = cwlen + 1;
    cfg.dictionary_present = false;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (cipher_type == BEAUFORT) cfg.beaufort = true;
    if (pt_len > 0) { cfg.plaintext_keyword_len_present = true; cfg.plaintext_keyword_len = pt_len; }
    if (ct_len > 0) { cfg.ciphertext_keyword_len_present = true; cfg.ciphertext_keyword_len = ct_len; }

    const char *crib = crib_in ? crib_in : "";

    // Silence solve_cipher's own (verbose) reporting during the run; we assert
    // on the returned SolveResult, not stdout. Restore the real stdout after.
    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved_stdout = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* fall through, just noisy */ }
    seed_rand(seed);
    solve_cipher(cipher_str, (char *) crib, &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved_stdout, fileno(stdout));
    close(saved_stdout);
    clearerr(stdout);
    double secs = ((double) clock() - t0) / CLOCKS_PER_SEC;

    double frac = res.solved ? recovered_fraction(&res, P, len) : 0.0;
    printf("[%-20s] solved=%d frac=%.3f score=%.2f %.1fs (budget %dx%d)\n",
        name, res.solved, frac, res.solved ? res.score : 0.0, secs, n_restarts, n_hillclimbs);

    CHECK(res.solved, "%s: solver reported no solution", name);
    CHECK(frac >= min_fraction, "%s: recovered %.3f of plaintext (< %.3f)", name, frac, min_fraction);
}

// Build a crib string: the plaintext at every `stride`-th position, '_' elsewhere.
static void make_crib(char *out, int stride) {
    int len = (int) strlen(PLAINTEXT);
    for (int i = 0; i < len; i++) out[i] = (i % stride == 0) ? PLAINTEXT[i] : '_';
    out[len] = '\0';
}

int main(void) {
    // Build the runtime alphabet (g_monograms, g_char_to_idx, ...) exactly as the
    // real binary's main does. derive_optimal_cycleword scores columns against
    // g_monograms, so without this the optimal-cycleword ciphers (porta, ...) fail.
    init_alphabet(NULL);

    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    int cw7[]  = {10, 17, 24, 15, 19, 14, 18};       // KRYPTOS
    int cw6[]  = {18, 4, 2, 17, 4, 19};              // SECRET
    int cw5[]  = {2, 14, 5, 5, 4};                   // COFFE

    // --- straight-keyword ciphers (cycleword derived deterministically) -------
    run_case("vigenere",          VIGENERE, 0, NULL, NULL, cw7, 7, 0, 0, NULL,   60, 1500, 12345u, 1.0);
    run_case("vigenere-variant",  VIGENERE, 1, NULL, NULL, cw5, 5, 0, 0, NULL,   60, 1500, 12345u, 1.0);
    run_case("beaufort",          BEAUFORT, 0, NULL, NULL, cw7, 7, 0, 0, NULL,   60, 1500, 12345u, 1.0);
    run_case("porta",             PORTA,    0, NULL, NULL, cw6, 6, 0, 0, NULL,   60, 1500, 12345u, 1.0);

    // --- keyword-search ciphers (the stochastic hill climber) -----------------
    // Budgets are set ~2x above the reliable floor measured across seeds
    // (Q1/Q2/Q3 solve at 120x1500, Q4 at 500x2500); a climber regression that
    // needs more iterations to reach the same instance fails recovery here.
    run_case("quagmire1",         QUAGMIRE_1, 0, "CODE", NULL,   cw6, 6, 4, 0, NULL,  250, 1500, 7u, 1.0);
    run_case("quagmire2",         QUAGMIRE_2, 0, NULL,   "CODE", cw6, 6, 0, 4, NULL,  250, 1500, 7u, 1.0);
    run_case("quagmire3",         QUAGMIRE_3, 0, "CODE", "CODE", cw6, 6, 4, 4, NULL,  250, 1500, 7u, 1.0);
    run_case("quagmire4",         QUAGMIRE_4, 0, "CAT",  "DOG",  cw5, 5, 3, 3, NULL,  700, 2500, 7u, 1.0);

    // --- crib-assisted (exercises crib_score + constrain_cycleword) ----------
    static char crib[MAX_CIPHER_LENGTH];
    make_crib(crib, 5);   // ~20% of positions revealed
    run_case("quagmire4-crib",    QUAGMIRE_4, 0, "CAT", "DOG", cw5, 5, 3, 3, crib, 300, 2500, 7u, 1.0);

    // --- autokey (straight-alphabet tableaus) --------------------------------
    // Here the "cycleword" is the autokey primer; the keystream is self-extending,
    // so optimal-cycleword derivation does not apply and the primer is searched
    // stochastically. cw5 (length 5) is the planted primer. Budgets are ~2x the
    // reliable floor (all three solve 8/8 across seeds at 60x1000).
    //
    // The KEYED-alphabet autokeys (auto1-4) are deliberately NOT in this suite:
    // their joint keyword+primer search does not recover reliably at any sane
    // budget (auto1 0/8 at 400x2000 and only 3/8 at 1200x3000; auto3 6/8 at
    // 600x2500), which would make for flaky regression tests. Their encrypt/
    // decrypt correctness is pinned by the KAT + round-trip checks in
    // test_ciphers.c instead.
    run_case("autokey0",          AUTOKEY_0,     0, NULL, NULL, cw5, 5, 0, 0, NULL, 120, 1500, 7u, 1.0);
    run_case("autokey-beaufort",  AUTOKEY_BEAU,  0, NULL, NULL, cw5, 5, 0, 0, NULL, 120, 1500, 7u, 1.0);
    run_case("autokey-porta",     AUTOKEY_PORTA, 0, NULL, NULL, cw5, 5, 0, 0, NULL, 120, 1500, 7u, 1.0);

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
