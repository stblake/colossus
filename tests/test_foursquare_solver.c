//
//  In-process stress / limits tests for the Four-Square solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each solve deterministic.
//
//  Four-Square is the hardest of the square family here: its state is TWO INDEPENDENT keyed
//  squares (the upper-right and lower-left, 50 cells), with the plaintext squares fixed. The
//  suite (planted-cipher recovery) does:
//    1. validates the per-type schedule registry (and that a non-registry type is untouched);
//    2. asserts the capability FLOOR -- a long cipher is recovered to ~100% at a tuned budget;
//    3. characterizes the LIMIT -- recovery vs ciphertext length is printed across a sweep
//       (the short-text cliff is visible; the marginal lengths are reported, not asserted);
//    4. extra-thorough: a PER-SQUARE recovery breakdown. In a Four-Square decrypt the first
//       letter of every digraph (even positions) is read through the upper-right square and
//       the second (odd positions) through the lower-left, so even- vs odd-position recovery
//       is a direct proxy for how well each keyed square was recovered; plus a MULTI-KEYWORD
//       sweep over several random key pairs reporting mean / worst recovery.
//
//  Run from the source directory so the n-gram table is found. Four-Square needs the
//  discriminating log-probability fitness, so this enables g_ngram_logprob.
//

#include "colossus.h"
#include "engine.h"
#include "scoring.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A long chunk of natural English (Pride and Prejudice, opening), letters only (~900).
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

static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Build the two keyed squares (upper-right, lower-left) from two keywords (J->I).
static void build_squares(const char *urk, const char *llk, int ur[], int ll[]) {
    int kw[64], kwn;
    kwn = 0;
    for (int i = 0; urk[i] && kwn < 64; i++) { int x = letter_to_index((unsigned char) urk[i]); if (x >= 0) kw[kwn++] = x; }
    playfair_grid_from_keyword(kw, kwn, ur);
    kwn = 0;
    for (int i = 0; llk[i] && kwn < 64; i++) { int x = letter_to_index((unsigned char) llk[i]); if (x >= 0) kw[kwn++] = x; }
    playfair_grid_from_keyword(kw, kwn, ll);
}

// Plant a Four-Square cipher: first pt_len letters of PLAINTEXT (padded to even with X),
// enciphered under the two keyword squares. Fills prepared[] (the expected solution) and
// cipher_str. Returns the prepared (even) length.
static int plant(const char *urk, const char *llk, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    if (n % 2 != 0) prepared[n++] = letter_to_index('X');
    int ur[SQUARE_GRID], ll[SQUARE_GRID];
    build_squares(urk, llk, ur, ll);
    int cipher[MAX_CIPHER_LENGTH];
    foursquare_encrypt(prepared, n, ur, ll, SQUARE_SIDE, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Run one solve; return the recovered fraction. If recovered != NULL, copy the decrypted
// plaintext into it (for the per-square even/odd breakdown).
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int method, bool use_registry, int n_restarts, int n_hillclimbs,
        double init_temp, double backtrack, uint32_t seed, double *secs_out, int *recovered) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = FOUR_SQUARE;
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
    for (int i = 0; i < plen; i++) {
        if (recovered) recovered[i] = res.decrypted[i];
        if (res.decrypted[i] == prepared[i]) ok++;
    }
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    init_config(&cfg);
    cfg.cipher_type = FOUR_SQUARE; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "foursquare registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 700000,
        "foursquare anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "foursquare anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg);
    cfg.cipher_type = FOUR_SQUARE; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "foursquare registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 40 && cfg.n_hill_climbs == 600000,
        "foursquare shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor + per-square breakdown -------------------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    static int recovered[MAX_CIPHER_LENGTH];
    double secs;

    int plen = plant("CIPHERKEYWORD", "PALIMPSEST", 900, prepared, cipher_str);
    double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
        8, 350000, 0.08, 0.30, 1u, &secs, recovered);
    printf("[capability anneal ~900ch 8x350k] frac=%.3f %.1fs\n", frac, secs);
    CHECK(frac >= 0.99, "foursquare ~900ch recovered only %.3f", frac);

    // Per-square breakdown: even positions decrypt through the upper-right square, odd through
    // the lower-left. Equal high recovery on both confirms BOTH independent squares were found.
    int ev_ok = 0, ev = 0, od_ok = 0, od = 0;
    for (int i = 0; i < plen; i++) {
        if (i % 2 == 0) { ev++; if (recovered[i] == prepared[i]) ev_ok++; }
        else            { od++; if (recovered[i] == prepared[i]) od_ok++; }
    }
    double fe = ev ? (double) ev_ok / ev : 0.0, fo = od ? (double) od_ok / od : 0.0;
    printf("[per-square] upper-right (even pos)=%.3f  lower-left (odd pos)=%.3f\n", fe, fo);
    CHECK(fe >= 0.99 && fo >= 0.99, "foursquare per-square recovery UR=%.3f LL=%.3f", fe, fo);
}

// --- 3. length cliff ----------------------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {300, 500, 700, 900};
    printf("recovery vs ciphertext length (foursquare, keys FOURSQUARE/CRYPTO, 8x350k, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("FOURSQUARE", "CRYPTO", lengths[li], prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
            8, 350000, 0.08, 0.30, 1u, &secs, NULL);
        printf("    len~%-4d (%d digraphs)  frac=%.3f  %.1fs\n", lengths[li], plen / 2, frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.99, "foursquare 900ch (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 4. multi-keyword sweep ---------------------------------------------------

static void test_multi_keyword(void) {
    const char *keys[][2] = {
        {"KRYPTOSABCDEF", "PALIMPSEST"},
        {"NORTHERLY",     "ABANDONED"},
        {"SHADOWFORCE",   "IODINEGAS"},
    };
    int nk = 3;
    double sum = 0.0, worst = 1.0;
    printf("multi-keyword sweep (foursquare, ~900ch, 8x350k, seed 1):\n");
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(keys[k][0], keys[k][1], 900, prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
            8, 350000, 0.08, 0.30, 1u, &secs, NULL);
        printf("    keys=%-13s/%-11s frac=%.3f %.1fs\n", keys[k][0], keys[k][1], frac, secs);
        sum += frac;
        if (frac < worst) worst = frac;
    }
    printf("    mean=%.3f worst=%.3f\n", sum / nk, worst);
    CHECK(sum / nk >= 0.95, "foursquare multi-keyword mean recovery only %.3f", sum / nk);
}

int main(void) {
    init_alphabet("J");                 // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == SQUARE_GRID, "alphabet size %d, expected %d", g_alpha, SQUARE_GRID);

    g_ngram_logprob = true;             // Four-Square needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability();
    test_length_cliff();
    test_multi_keyword();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
