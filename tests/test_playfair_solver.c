//
//  In-process stress / limits tests for the Playfair solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be
//  driven directly and its SolveResult inspected. A fixed -seed makes each
//  stochastic solve deterministic.
//
//  Strategy (planted-cipher recovery): encipher a known English plaintext under a
//  known Playfair keyword, run solve_cipher with a fixed seed and a bounded budget,
//  and measure the recovered fraction. The suite does three things:
//    1. validates the per-type schedule registry (apply_cipher_defaults), including
//       that a non-registry type is left untouched (regression safety);
//    2. asserts the capability FLOOR -- a long (800-char) cipher is recovered to
//       ~100% at the tuned schedule, under both the registry default and an explicit
//       budget, and that annealing beats shotgun;
//    3. characterizes the LIMIT -- recovery vs ciphertext length is printed across a
//       sweep so the short-text cliff is visible (Playfair is genuinely near the edge
//       of a quadgram attack below a few hundred characters). The marginal lengths are
//       reported, not asserted, so the suite stays deterministic rather than flaky.
//
//  Run from the source directory (as `make testopt` does) so the n-gram table is
//  found in the cwd. Playfair needs the discriminating log-probability scoring, so
//  this enables g_ngram_logprob and loads the quadgram table in that mode.
//

#include "colossus.h"
#include "engine.h"   // apply_cipher_defaults (moved out of colossus.h)
#include "scoring.h"  // load_ngrams

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A long chunk of natural English (Pride and Prejudice, opening), letters only.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNE"
    "MUSTBEINWANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHIS"
    "FIRSTENTERINGANEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDING"
    "FAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGHTERS"
    "MYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLETATLAST"
    "MRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASJUSTBEENHEREANDSHE"
    "TOLDMEALLABOUTITMRBENNETMADENOANSWERDOYOUNOTWANTTOKNOWWHOHASTAKENITCRIEDHISWIFE"
    "IMPATIENTLYYOUWANTTOTELLMEANDIHAVENOOBJECTIONTOHEARINGITTHISWASINVITATIONENOUGH"
    "WHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATNETHERFIELDISTAKENBYAYOUNGMANOFLARGEFORTUNE"
    "FROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEETHEPLACE";

// A..Z char -> 0..24 alphabet index, merging J into I (the 25-letter convention).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Playfair cipher: take the first pt_len letters of PLAINTEXT, build the grid
// from `keyword`, prepare + encrypt. Fills prepared[] (the expected solution) and
// cipher_str (the A..Z ciphertext to hand to solve_cipher). Returns the prepared len.
static int plant(const char *keyword, int pt_len, int prepared[], char cipher_str[]) {
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
    int plen = playfair_prepare(raw, n, filler, alt, prepared, MAX_CIPHER_LENGTH);

    int cipher[MAX_CIPHER_LENGTH];
    playfair_encrypt(prepared, plen, grid, cipher);
    for (int i = 0; i < plen; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[plen] = '\0';
    return plen;
}

// Run one solve and return the recovered fraction (vs the planted prepared plaintext).
// `use_registry` applies the tuned per-type schedule instead of the passed budget.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int method, bool use_registry, int n_restarts, int n_hillclimbs,
        double init_temp, double backtrack, uint32_t seed, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = PLAYFAIR;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");

    if (use_registry) {
        apply_cipher_defaults(&cfg, false);
    } else {
        cfg.n_restarts = n_restarts;
        cfg.n_hill_climbs = n_hillclimbs;
        cfg.init_temp = init_temp;
        cfg.backtracking_probability = backtrack;
    }

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

    init_config(&cfg);
    cfg.cipher_type = PLAYFAIR; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "playfair registry: no entry applied");
    CHECK(cfg.n_restarts == 6 && cfg.n_hill_climbs == 400000,
        "playfair anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "playfair anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg);
    cfg.cipher_type = PLAYFAIR; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "playfair registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 30 && cfg.n_hill_climbs == 300000,
        "playfair shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor ------------------------------------------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant("CIPHERKEYWORD", 800, prepared, cipher_str);

    // The tuned anneal schedule (via the per-type registry -- this also exercises
    // apply_cipher_defaults through a real solve) recovers an 800-char cipher to
    // ~100% (3/3 seeds in calibration).
    double secs = 0.0;
    double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, true,
        0, 0, 0.0, 0.0, 1u, &secs);
    printf("[capability registry-anneal 800ch] frac=%.3f %.1fs\n", frac, secs);
    CHECK(frac >= 0.99, "playfair 800ch registry-default recovered only %.3f", frac);

    // Annealing beats shotgun on Playfair. Shotgun is run at a comparable budget and
    // is expected to trail; we assert annealing wins, and only report shotgun's score.
    double sg = solve_and_frac(cipher_str, prepared, plen, METHOD_SHOTGUN, false,
        10, 200000, 0.08, 0.20, 1u, &secs);
    printf("[shotgun 800ch 10x200k] frac=%.3f %.1fs\n", sg, secs);
    CHECK(frac >= sg, "annealing (%.3f) did not beat shotgun (%.3f)", frac, sg);
}

// --- 3. length cliff ----------------------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {150, 300, 500, 800};
    printf("recovery vs ciphertext length (keyword PLAYFAIRGRID, 4x250k anneal, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("PLAYFAIRGRID", lengths[li], prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
            4, 250000, 0.08, 0.30, 1u, &secs);
        printf("    len~%-4d (%d digraphs)  frac=%.3f  %.1fs\n",
            lengths[li], plen / 2, frac, secs);
        if (li == 3) frac_longest = frac;
    }
    // Only the longest length is a hard assertion (the short ones document the cliff).
    CHECK(frac_longest >= 0.99, "playfair 800ch (cliff sweep) recovered only %.3f", frac_longest);
}

int main(void) {
    init_alphabet("J");                 // 25-letter Playfair alphabet (J merged into I)
    CHECK(g_alpha == PLAYFAIR_GRID, "alphabet size %d, expected %d", g_alpha, PLAYFAIR_GRID);

    g_ngram_logprob = true;             // Playfair needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability();
    test_length_cliff();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
