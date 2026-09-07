//
//  In-process stress / limits tests for the Two-Square solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each stochastic solve
//  deterministic.
//
//  Two-Square carries a PAIR of 5x5 squares (50 cells, double Playfair's), so it needs
//  more text and a bigger budget than Playfair. The suite (planted-cipher recovery) does:
//    1. validates the per-type schedule registry for BOTH arrangements (horizontal +
//       vertical), incl. that a non-registry type is left untouched;
//    2. asserts the capability FLOOR -- a long cipher is recovered to ~100% for each
//       arrangement at the tuned schedule;
//    3. characterizes the LIMIT -- recovery vs ciphertext length is printed across a sweep
//       so the short-text cliff is visible (the marginal lengths are reported, not asserted);
//    4. extra-thorough: measures the digraph TRANSPARENCY rate (~20% of digraphs leak as
//       same-row reversals / same-column fixed pairs -- the cipher's documented weakness)
//       and confirms recovery survives it, and runs a MULTI-KEYWORD sweep (several random
//       key pairs at the capability length) reporting mean / worst recovery.
//
//  Run from the source directory so the n-gram table is found. Two-Square needs the
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

// Build two keyed squares from two keywords (J->I), reusing the Playfair keyword build.
static void build_squares(const char *kw1, const char *kw2, int sq1[], int sq2[]) {
    int kw[64], kwn;
    kwn = 0;
    for (int i = 0; kw1[i] && kwn < 64; i++) { int x = letter_to_index((unsigned char) kw1[i]); if (x >= 0) kw[kwn++] = x; }
    playfair_grid_from_keyword(kw, kwn, sq1);
    kwn = 0;
    for (int i = 0; kw2[i] && kwn < 64; i++) { int x = letter_to_index((unsigned char) kw2[i]); if (x >= 0) kw[kwn++] = x; }
    playfair_grid_from_keyword(kw, kwn, sq2);
}

// Plant a Two-Square cipher: take the first pt_len letters of PLAINTEXT (padded to an even
// length with X), encipher under the two keywords + arrangement. Fills prepared[] (the
// expected solution) and cipher_str. Returns the prepared (even) length.
static int plant(const char *kw1, const char *kw2, int variant, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    if (n % 2 != 0) prepared[n++] = letter_to_index('X');
    int sq1[SQUARE_GRID], sq2[SQUARE_GRID];
    build_squares(kw1, kw2, sq1, sq2);
    int cipher[MAX_CIPHER_LENGTH];
    twosquare_encrypt(prepared, n, sq1, sq2, SQUARE_SIDE, variant, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Run one solve and return the recovered fraction (vs the planted prepared plaintext).
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int cipher_type, int method, bool use_registry, int n_restarts, int n_hillclimbs,
        double init_temp, double backtrack, uint32_t seed, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = cipher_type;
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
    int types[2] = { TWO_SQUARE, TWO_SQUARE_V };
    for (int ti = 0; ti < 2; ti++) {
        init_config(&cfg);
        cfg.cipher_type = types[ti]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "twosquare registry: no entry for type %d", types[ti]);
        CHECK(cfg.n_restarts == 8 && cfg.n_hill_climbs == 600000,
            "twosquare anneal defaults wrong (type %d): %dx%d", types[ti], cfg.n_restarts, cfg.n_hill_climbs);
        CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
            "twosquare anneal inittemp wrong: %.4f", cfg.init_temp);

        init_config(&cfg);
        cfg.cipher_type = types[ti]; cfg.method = METHOD_SHOTGUN;
        CHECK(apply_cipher_defaults(&cfg, false), "twosquare registry (shotgun): no entry");
        CHECK(cfg.n_restarts == 30 && cfg.n_hill_climbs == 500000,
            "twosquare shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    }

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor (both arrangements) ----------------------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    double secs;

    // Horizontal (ACA): exercise apply_cipher_defaults end-to-end via the registry default.
    int plen = plant("CIPHERKEYWORD", "PALIMPSEST", TWO_SQ_HORIZONTAL, 900, prepared, cipher_str);
    double fh = solve_and_frac(cipher_str, prepared, plen, TWO_SQUARE, METHOD_DEFAULT, true,
        0, 0, 0.0, 0.0, 1u, &secs);
    printf("[capability registry-anneal horizontal ~900ch] frac=%.3f %.1fs\n", fh, secs);
    CHECK(fh >= 0.99, "twosquare horizontal ~900ch registry-default recovered only %.3f", fh);

    // Vertical (self-inverse): a trimmed explicit budget (keeps the suite quick).
    plen = plant("CIPHERKEYWORD", "PALIMPSEST", TWO_SQ_VERTICAL, 900, prepared, cipher_str);
    double fv = solve_and_frac(cipher_str, prepared, plen, TWO_SQUARE_V, METHOD_DEFAULT, false,
        6, 300000, 0.08, 0.30, 1u, &secs);
    printf("[capability anneal vertical ~900ch 6x300k] frac=%.3f %.1fs\n", fv, secs);
    CHECK(fv >= 0.99, "twosquare vertical ~900ch recovered only %.3f", fv);
}

// --- 3. length cliff (horizontal) ---------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {300, 500, 700, 900};
    printf("recovery vs ciphertext length (twosquare horizontal, keys TWOSQUARE/CRYPTO, 6x300k, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("TWOSQUARE", "CRYPTO", TWO_SQ_HORIZONTAL, lengths[li], prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, TWO_SQUARE, METHOD_DEFAULT, false,
            6, 300000, 0.08, 0.30, 1u, &secs);
        printf("    len~%-4d (%d digraphs)  frac=%.3f  %.1fs\n", lengths[li], plen / 2, frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.99, "twosquare 900ch (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 4a. transparency rate ----------------------------------------------------
//
// Measure the fraction of plaintext digraphs that are transparencies (horizontal: the two
// letters share a row across the two squares -> the cipher reverses them). ~20% is the
// documented weakness; we report it and confirm the capability solve survived it.

static void test_transparency_rate(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant("TWOSQUARE", "CRYPTO", TWO_SQ_HORIZONTAL, 900, prepared, cipher_str);
    int sq1[SQUARE_GRID], sq2[SQUARE_GRID];
    build_squares("TWOSQUARE", "CRYPTO", sq1, sq2);
    int pos1[SQUARE_GRID], pos2[SQUARE_GRID];
    for (int p = 0; p < SQUARE_GRID; p++) { pos1[sq1[p]] = p; pos2[sq2[p]] = p; }
    int s = SQUARE_SIDE, transparencies = 0, pairs = 0;
    for (int i = 0; i + 1 < plen; i += 2, pairs++)
        if (pos1[prepared[i]] / s == pos2[prepared[i + 1]] / s) transparencies++;
    double rate = pairs ? (double) transparencies / pairs : 0.0;
    printf("[transparency] horizontal same-row digraphs: %d/%d (%.1f%%)\n",
        transparencies, pairs, 100.0 * rate);
    // Loose sanity band (English over these squares lands near the ~1/5 prior).
    CHECK(rate > 0.05 && rate < 0.40, "twosquare transparency rate %.3f outside expected band", rate);
}

// --- 4b. multi-keyword sweep --------------------------------------------------

static void test_multi_keyword(void) {
    const char *keys[][2] = {
        {"KRYPTOSABCDEF", "PALIMPSEST"},
        {"NORTHERLY",     "ABANDONED"},
        {"SHADOWFORCE",   "IODINEGAS"},
    };
    int nk = 3;
    double sum = 0.0, worst = 1.0;
    printf("multi-keyword sweep (twosquare horizontal, ~900ch, 6x300k, seed 1):\n");
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(keys[k][0], keys[k][1], TWO_SQ_HORIZONTAL, 900, prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, TWO_SQUARE, METHOD_DEFAULT, false,
            6, 300000, 0.08, 0.30, 1u, &secs);
        printf("    keys=%-13s/%-11s frac=%.3f %.1fs\n", keys[k][0], keys[k][1], frac, secs);
        sum += frac;
        if (frac < worst) worst = frac;
    }
    printf("    mean=%.3f worst=%.3f\n", sum / nk, worst);
    CHECK(sum / nk >= 0.95, "twosquare multi-keyword mean recovery only %.3f", sum / nk);
}

int main(void) {
    init_alphabet("J");                 // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == SQUARE_GRID, "alphabet size %d, expected %d", g_alpha, SQUARE_GRID);

    g_ngram_logprob = true;             // Two-Square needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability();
    test_length_cliff();
    test_transparency_rate();
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
