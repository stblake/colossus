//
//  In-process stress / limits tests for the Hill solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each stochastic solve
//  deterministic.
//
//  Strategy (planted-cipher recovery): encipher a known English plaintext under a known
//  invertible key matrix and block size k, then attack it ciphertext-only. The suite does
//  four things:
//    1. validates the per-type schedule registry (apply_cipher_defaults), including that a
//       non-registry type is left untouched (regression safety);
//    2. tests BLOCK-SIZE SELECTION -- with k swept (not pinned) the solver must pick the
//       planted k (a wrong k decrypts to gibberish and loses on n-gram score) and recover
//       the plaintext;
//    3. asserts the capability FLOOR -- a long k=2 and a long k=3 cipher are recovered to
//       ~100%;
//    4. characterizes the LIMIT -- recovery vs ciphertext length is printed (k=2 pinned)
//       so the short-text cliff is visible.
//
//  k=2 (26^4 keys) and k=3 are reliably breakable ciphertext-only; k=4/5 are exercised
//  only by the primitive round-trip/inverse tests (test_hill.c), not asserted here.
//
//  Run from the source directory so the n-gram table is found in the cwd. Hill, like
//  Playfair/Bifid/Trifid, needs the discriminating log-probability scoring, so this
//  enables g_ngram_logprob and loads the quadgram table in that mode.
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

// A long chunk of natural English (Pride and Prejudice, opening), letters only. Long
// enough that even k=3 (a 9-entry matrix) has the material it needs.
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
    "FROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEETHEPLACEANDWAS"
    "SOMUCHDELIGHTEDWITHITTHATHEAGREEDWITHMRMORRISIMMEDIATELYTHATHEISTOTAKEPOSSESSION"
    "BEFOREMICHAELMASANDSOMEOFHISSERVANTSARETOBEINTHEHOUSEBYTHEENDOFNEXTWEEKWHATISHIS"
    "NAMEBINGLEYISHEMARRIEDORSINGLEOHSINGLEMYDEARTOBESUREASINGLEMANOFLARGEFORTUNEFOUR"
    "ORFIVETHOUSANDAYEARWHATAFINETHINGFOROURGIRLSHOWSOHOWCANITAFFECTTHEMMYDEARMRBENNET"
    "REPLIEDHISWIFEHOWCANYOUBESOTIRESOMEYOUMUSTKNOWTHATIAMTHINKINGOFHISMARRYINGONEOFTHEM";

// A..Z char -> 0..25 alphabet index (full 26-letter alphabet, no J merge).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Hill cipher: take the first pt_len letters of PLAINTEXT (truncated to a whole
// number of k-blocks), build an invertible key matrix from `keyword` (the same
// deterministic invertibility retry the generator uses), and encipher. Fills prepared[]
// (the expected solution) and cipher_str (the A..Z ciphertext). Returns the length.
static int plant(const char *keyword, int pt_len, int k, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    n -= n % k;                                  // whole number of k-blocks

    int kw[64], kwn = 0;
    for (int i = 0; keyword[i] && kwn < 64; i++) {
        int idx = letter_to_index((unsigned char) keyword[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }

    // Invertible-mod-26 key from the keyword: offset entries by a base-26 odometer indexed
    // by `attempt` until invertible (same scheme as tools/hill_gen.c, so the planted ciphers
    // match generator output). attempt 0 is the untweaked keyword matrix.
    int km = k * k, mat[HILL_MAX_KEY], inv[HILL_MAX_KEY], found = 0;
    for (int attempt = 0; attempt < 1000000 && !found; attempt++) {
        hill_matrix_from_keyword(kw, kwn, mat, k);
        for (int a = attempt, i = 0; a > 0 && i < km; a /= ALPHABET_SIZE, i++)
            mat[i] = (mat[i] + a % ALPHABET_SIZE) % ALPHABET_SIZE;
        found = hill_mat_inverse(mat, k, inv);
    }

    int cipher[MAX_CIPHER_LENGTH];
    hill_encrypt(prepared, n, mat, k, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Run one solve and return the recovered fraction. `period` == 0 sweeps block sizes (up
// to max_period); >0 pins k. `use_registry` applies the tuned per-type schedule. If
// k_out is non-NULL it receives the block size the solver reported.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        bool use_registry, int n_restarts, int n_hillclimbs, double init_temp,
        double backtrack, int period, int max_period, uint32_t seed,
        int *k_out, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = HILL;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (max_period > 0) cfg.max_period = max_period;

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

    if (k_out) *k_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    init_config(&cfg);
    cfg.cipher_type = HILL; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "hill registry: no entry applied");
    CHECK(cfg.n_restarts == 250 && cfg.n_hill_climbs == 8000,
        "hill anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0999 && cfg.init_temp < 0.1001,
        "hill anneal inittemp wrong: %.4f", cfg.init_temp);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. block-size selection (k swept, not pinned) ----------------------------

static void test_block_size_selection(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];

    // Plant k=2; sweep k in [2..3]; the solver must pick k=2 and recover the plaintext.
    int plen = plant("KRYPTOSHILLCIPHER", 600, 2, prepared, cipher_str);
    int got_k = -1; double secs = 0.0;
    double frac = solve_and_frac(cipher_str, prepared, plen, false,
        20, 40000, 0.10, 0.25, 0, 3, 1u, &got_k, &secs);
    printf("[selection k=2 swept 2..3, 600ch] frac=%.3f k=%d %.1fs\n", frac, got_k, secs);
    CHECK(got_k == 2, "hill selection chose k=%d, planted 2", got_k);
    CHECK(frac >= 0.95, "hill k=2 (swept) recovered only %.3f", frac);
}

// --- 3. capability floor ------------------------------------------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];

    // k=2, long cipher, registry anneal (k pinned so only the k=2 config runs -- this
    // exercises the registry schedule end-to-end without the full 2..5 sweep cost).
    int plen2 = plant("KRYPTOSHILLCIPHER", 800, 2, prepared, cipher_str);
    int got_k = -1; double secs = 0.0;
    double frac2 = solve_and_frac(cipher_str, prepared, plen2, true,
        0, 0, 0.0, 0.0, 2, 0, 1u, &got_k, &secs);
    printf("[capability registry-anneal k=2 ~800ch] frac=%.3f k=%d %.1fs\n",
        frac2, got_k, secs);
    CHECK(frac2 >= 0.95, "hill k=2 (~800ch, registry) recovered only %.3f", frac2);
    CHECK(got_k == 2, "hill recovered k=%d, planted 2", got_k);

    // k=3, long cipher, k pinned. k=3 (26^9 keys) is the hard ciphertext-only case; many
    // short restarts (the registry philosophy) crack it from ~1100+ characters.
    int plen3 = plant("KRYPTOSHILLCIPHER", 1200, 3, prepared, cipher_str);
    secs = 0.0;
    double frac3 = solve_and_frac(cipher_str, prepared, plen3, false,
        500, 8000, 0.10, 0.25, 3, 0, 1u, NULL, &secs);
    printf("[capability k=3 pinned ~1200ch] frac=%.3f %.1fs\n", frac3, secs);
    CHECK(frac3 >= 0.95, "hill k=3 (~1200ch, pinned) recovered only %.3f", frac3);
}

// --- 4. length cliff (k=2 pinned) ---------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {100, 200, 300, 500};
    printf("recovery vs ciphertext length (keyword KRYPTOSHILLCIPHER, k=2 pinned, 10x80k anneal, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("KRYPTOSHILLCIPHER", lengths[li], 2, prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(cipher_str, prepared, plen, false,
            10, 80000, 0.10, 0.25, 2, 0, 1u, NULL, &secs);
        printf("    len~%-4d  frac=%.3f  %.1fs\n", lengths[li], frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.95, "hill k=2 500ch (cliff sweep) recovered only %.3f", frac_longest);
}

int main(void) {
    init_alphabet(NULL);                 // full 26-letter Hill alphabet
    CHECK(g_alpha == ALPHABET_SIZE, "alphabet size %d, expected 26", g_alpha);

    g_ngram_logprob = true;              // Hill needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_block_size_selection();
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
