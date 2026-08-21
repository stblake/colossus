//
//  In-process stress / limits tests for the Phillips solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each stochastic solve
//  deterministic.
//
//  Strategy (planted-cipher recovery): encipher a known English plaintext under a known
//  base square (built from a keyword) for each variant, run solve_cipher with a fixed seed
//  and a bounded budget, and measure the recovered fraction. The suite does three things:
//    1. validates the per-type schedule registry (apply_cipher_defaults) for all three
//       Phillips types, including that a non-registry type is left untouched;
//    2. asserts the capability FLOOR -- a long (~760-char) cipher is recovered to ~100%
//       for EACH variant (Row via the registry default, which also exercises
//       apply_cipher_defaults through a real solve; Column and Row-Column via an explicit
//       budget);
//    3. characterizes the LIMIT -- recovery vs ciphertext length is printed across a sweep
//       so the short-text cliff is visible (Phillips, being monographic, recovers from
//       shorter text than digraphic Playfair -- reliably from ~200 characters). Only the
//       longest length is asserted, so the suite stays deterministic rather than flaky.
//
//  Run from the source directory (as `make testopt` does) so the n-gram table is found in
//  the cwd. Like Playfair/Bifid, Phillips needs the discriminating log-probability scoring,
//  so this enables g_ngram_logprob and loads the quadgram table in that mode.
//

#include "colossus.h"
#include "engine.h"   // apply_cipher_defaults
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

// Plant a Phillips cipher of the given variant: take the first pt_len letters of PLAINTEXT,
// build the base square from `keyword`, encipher. Fills plain[] (the expected solution) and
// cipher_str (the A..Z ciphertext to hand to solve_cipher). Returns the length.
static int plant(const char *keyword, int variant, int pt_len, int plain[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) plain[n++] = idx;
    }
    int kw[64], kwn = 0;
    for (int i = 0; keyword[i] && kwn < 64; i++) {
        int idx = letter_to_index((unsigned char) keyword[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int grid[PHILLIPS_GRID];
    phillips_grid_from_keyword(kw, kwn, grid, PHILLIPS_GRID);

    int cipher[MAX_CIPHER_LENGTH];
    phillips_encrypt(plain, n, grid, PHILLIPS_SIDE, variant, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Map a Phillips square-generation variant to its solver cipher type.
static int type_for_variant(int variant) {
    if (variant == PHILLIPS_COL) return PHILLIPS_C;
    if (variant == PHILLIPS_ROWCOL) return PHILLIPS_RC;
    return PHILLIPS;
}

// Run one solve and return the recovered fraction (vs the planted plaintext). `use_registry`
// applies the tuned per-type schedule instead of the passed budget.
static double solve_and_frac(int cipher_type, const char *cipher_str, const int plain[], int plen,
        int method, bool use_registry, int n_restarts, int n_hillclimbs,
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
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == plain[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    const int types[3] = { PHILLIPS, PHILLIPS_C, PHILLIPS_RC };
    const char *names[3] = { "phillips", "phillips-c", "phillips-rc" };
    for (int i = 0; i < 3; i++) {
        ColossusConfig cfg;
        init_config(&cfg);
        cfg.cipher_type = types[i]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "%s registry: no entry applied", names[i]);
        CHECK(cfg.n_restarts == 4 && cfg.n_hill_climbs == 250000,
            "%s anneal defaults wrong: %dx%d", names[i], cfg.n_restarts, cfg.n_hill_climbs);
        CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
            "%s anneal inittemp wrong: %.4f", names[i], cfg.init_temp);

        init_config(&cfg);
        cfg.cipher_type = types[i]; cfg.method = METHOD_SHOTGUN;
        CHECK(apply_cipher_defaults(&cfg, false), "%s registry (shotgun): no entry", names[i]);
        CHECK(cfg.n_restarts == 20 && cfg.n_hill_climbs == 250000,
            "%s shotgun defaults wrong: %dx%d", names[i], cfg.n_restarts, cfg.n_hill_climbs);
    }

    // Regression safety: a type with no registry entry is left untouched.
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor (each variant) ---------------------------------------

static void test_capability(void) {
    static int plain[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];

    // Row via the registry default -- this also exercises apply_cipher_defaults through a
    // real solve. ~760 chars recovers to ~100%.
    int plen = plant("KRYPTOSPHILLIPS", PHILLIPS_ROW, 760, plain, cipher_str);
    double secs = 0.0;
    double frac = solve_and_frac(PHILLIPS, cipher_str, plain, plen, METHOD_DEFAULT, true,
        0, 0, 0.0, 0.0, 1u, &secs);
    printf("[capability registry-anneal row 760ch] frac=%.3f %.1fs\n", frac, secs);
    CHECK(frac >= 0.99, "phillips(row) 760ch registry-default recovered only %.3f", frac);

    // Column and Row-Column at an explicit budget (4x200000 lands ~100% in ~15s).
    const int variants[2] = { PHILLIPS_COL, PHILLIPS_ROWCOL };
    const char *names[2] = { "phillips-c", "phillips-rc" };
    for (int v = 0; v < 2; v++) {
        plen = plant("KRYPTOSPHILLIPS", variants[v], 760, plain, cipher_str);
        frac = solve_and_frac(type_for_variant(variants[v]), cipher_str, plain, plen,
            METHOD_DEFAULT, false, 4, 200000, 0.08, 0.30, 1u, &secs);
        printf("[capability %s 760ch 4x200k] frac=%.3f %.1fs\n", names[v], frac, secs);
        CHECK(frac >= 0.99, "%s 760ch recovered only %.3f", names[v], frac);
    }
}

// --- 3. length cliff (Row) ----------------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {150, 200, 400, 760};
    printf("recovery vs ciphertext length (Row, keyword PHILLIPSGRID, 4x200k anneal, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int plain[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("PHILLIPSGRID", PHILLIPS_ROW, lengths[li], plain, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(PHILLIPS, cipher_str, plain, plen, METHOD_DEFAULT, false,
            4, 200000, 0.08, 0.30, 1u, &secs);
        printf("    len~%-4d  frac=%.3f  %.1fs\n", lengths[li], frac, secs);
        if (li == 3) frac_longest = frac;
    }
    // Only the longest length is a hard assertion (the short ones document the cliff).
    CHECK(frac_longest >= 0.99, "phillips 760ch (cliff sweep) recovered only %.3f", frac_longest);
}

int main(void) {
    init_alphabet("J");                 // 25-letter Phillips alphabet (J merged into I)
    CHECK(g_alpha == PHILLIPS_GRID, "alphabet size %d, expected %d", g_alpha, PHILLIPS_GRID);

    g_ngram_logprob = true;             // Phillips needs the log-probability fitness
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
