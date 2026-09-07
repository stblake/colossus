//
//  In-process stress / limits tests for the Trifid solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be
//  driven directly and its SolveResult inspected. A fixed -seed makes each stochastic
//  solve deterministic.
//
//  Strategy (planted-cipher recovery): encipher a known English plaintext under a known
//  keyed cube and period, then attack it. The suite does four things:
//    1. validates the per-type schedule registry (apply_cipher_defaults), including that
//       a non-registry type is left untouched (regression safety);
//    2. tests the PERIOD ESTIMATOR in isolation (trifid_estimate_periods) -- the true
//       period must land in the top-K columnar-IoC candidates across periods and lengths
//       (the IoC also peaks at multiples of the true period, so rank 1 is not required,
//       only top-K membership, which is all the solver needs);
//    3. asserts the capability FLOOR end-to-end -- a long cipher is recovered to ~100%
//       with the period ESTIMATED (not pinned), and the reported period is correct;
//    4. characterizes the LIMIT -- recovery vs ciphertext length is printed (period
//       pinned to isolate cube recovery) so the short-text cliff is visible.
//
//  Run from the source directory so the n-gram table is found in the cwd. Trifid, like
//  Bifid and Playfair, needs the discriminating log-probability scoring, so this enables
//  g_ngram_logprob and loads the quadgram table in that mode.
//

#include "colossus.h"
#include "engine.h"        // apply_cipher_defaults
#include "scoring.h"       // load_ngrams
#include "trifid_solver.h" // trifid_estimate_periods

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

// A..Z char -> 0..25 alphabet index over the 27-symbol Trifid alphabet.
static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Trifid cipher: take the first pt_len letters of PLAINTEXT, build the cube
// from `keyword`, encipher with `period`. Fills prepared[] (the expected solution, bare
// A..Z) and cipher_str (the ciphertext over the 27-symbol alphabet). Returns the length.
static int plant(const char *keyword, int pt_len, int period, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    int kw[64], kwn = 0;
    for (int i = 0; keyword[i] && kwn < 64; i++) {
        int idx = letter_to_index((unsigned char) keyword[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int cube[TRIFID_CELLS];
    trifid_cube_from_keyword(kw, kwn, cube, g_alpha);

    int cipher[MAX_CIPHER_LENGTH];
    trifid_encrypt(prepared, n, cube, TRIFID_SIDE, period, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Run one solve and return the recovered fraction. `period` == 0 estimates the period;
// >0 pins it. `use_registry` applies the tuned per-type schedule. If period_out is
// non-NULL it receives the period the solver reported.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        bool use_registry, int n_restarts, int n_hillclimbs, double init_temp,
        double backtrack, int period, int n_periods, uint32_t seed,
        int *period_out, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TRIFID;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (n_periods > 0) cfg.n_periods = n_periods;

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

    if (period_out) *period_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    init_config(&cfg);
    cfg.cipher_type = TRIFID; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "trifid registry: no entry applied");
    CHECK(cfg.n_restarts == 6 && cfg.n_hill_climbs == 300000,
        "trifid anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "trifid anneal inittemp wrong: %.4f", cfg.init_temp);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. period estimator ------------------------------------------------------

static void test_period_estimator(void) {
    int periods[] = {5, 7, 9, 11, 13};
    int lengths[] = {760, 450};
    int trials = 0, hits = 0, rank1 = 0;
    for (int pi = 0; pi < 5; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cipher_str[MAX_CIPHER_LENGTH];
            int p = periods[pi];
            int plen = plant("KRYPTOSABCDEF", lengths[li], p, prepared, cipher_str);

            int cidx[MAX_CIPHER_LENGTH];
            for (int i = 0; i < plen; i++) {
                int c = toupper((unsigned char) cipher_str[i]);
                cidx[i] = (c < 128) ? g_char_to_idx[c] : -1;
            }

            int out[8];
            int n = trifid_estimate_periods(cidx, plen, 2, 20, 5, out, false);
            trials++;
            int found = 0;
            for (int i = 0; i < n; i++) if (out[i] == p) { found = 1; if (i == 0) rank1++; }
            if (found) hits++;
            else printf("  [estimator MISS] period %d len %d -> {%d %d %d %d %d}\n",
                p, lengths[li], out[0], out[1], out[2], out[3], out[4]);
        }
    }
    printf("[period estimator] true period in top-5 for %d/%d cases (rank-1 for %d)\n",
        hits, trials, rank1);
    // The true period must be a top-5 candidate every time -- that is all the solver
    // needs (it anneals all five and the n-gram score discards multiples-of-period).
    CHECK(hits == trials, "period estimator missed the true period in %d/%d cases",
        trials - hits, trials);
}

// --- 3. capability floor (period estimated, not pinned) -----------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int period = 7;
    int plen = plant("KRYPTOSABCDEF", 760, period, prepared, cipher_str);

    // Estimate the period (period == 0) and break the cube at the registry default.
    // n_periods is trimmed to 3 (the true period is rank-1 ~70% of the time and always
    // top-3) so the end-to-end estimated-period solve stays under ~1 minute.
    int got_period = -1;
    double secs = 0.0;
    double frac = solve_and_frac(cipher_str, prepared, plen, true,
        0, 0, 0.0, 0.0, 0, 3, 1u, &got_period, &secs);
    printf("[capability registry-anneal 760ch, period estimated, top-3] frac=%.3f period=%d %.1fs\n",
        frac, got_period, secs);
    CHECK(frac >= 0.95, "trifid 760ch (period estimated) recovered only %.3f", frac);
    CHECK(got_period == period, "trifid recovered period %d, planted %d", got_period, period);
}

// --- 4. length cliff (period pinned to isolate cube recovery) -----------------

static void test_length_cliff(void) {
    int lengths[] = {300, 500, 760};
    int nlen = 3;
    printf("recovery vs ciphertext length (keyword KRYPTOSABCDEF, period 7 pinned, 6x300k anneal, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < nlen; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("KRYPTOSABCDEF", lengths[li], 7, prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(cipher_str, prepared, plen, false,
            6, 300000, 0.08, 0.30, 7, 1, 1u, NULL, &secs);
        printf("    len~%-4d  frac=%.3f  %.1fs\n", lengths[li], frac, secs);
        if (li == nlen - 1) frac_longest = frac;
    }
    // Only the longest length is a hard assertion (the short ones document the cliff).
    CHECK(frac_longest >= 0.95, "trifid 760ch (cliff sweep) recovered only %.3f", frac_longest);
}

int main(void) {
    init_alphabet_trifid();             // 27-symbol Trifid alphabet (A..Z + '+')
    CHECK(g_alpha == TRIFID_CELLS, "alphabet size %d, expected %d", g_alpha, TRIFID_CELLS);

    g_ngram_logprob = true;             // Trifid needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_period_estimator();
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
