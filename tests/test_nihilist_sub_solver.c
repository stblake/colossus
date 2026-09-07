//
//  In-process stress / limits tests for the Nihilist Substitution solver (solve_cipher),
//  run for EACH addition sub-type (carry / no-carry / mod-100).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Nihilist Substitution is a COUPLED keyed-square + periodic-additive cipher (the twin of
//  ADFGVX); the additive key is recovered by a square-independent validity reward, the square
//  by n-grams. This suite is deliberately thorough and per sub-type:
//    1. registry validation (apply_cipher_defaults) for all three codes + a non-registry type;
//    2. period-estimator top-K hit rate, measured separately per convention (the columnar-IoC
//       peak is the additive's, present under every add rule);
//    3. capability floor + length cliff PER CONVENTION (period pinned): recovery printed vs
//       length, the longest asserted >= 0.95;
//    4. BLIND period (carry): period estimated end-to-end, recovery high and the reported
//       period a multiple of the true one (IoC peaks at multiples too);
//    5. multi-keyword sweep (carry): mean / worst recovery over several square+additive keys;
//    6. keyed-label end-to-end (carry): a label-keyed cipher recovered as the relabelled square.
//
//  All three conventions share the 25-letter (J->I) alphabet and base-25 quadgrams, and need
//  the log-probability fitness (g_ngram_logprob). Run from the source directory.
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "nihilist_sub_solver.h"  // nihilist_sub_estimate_periods

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

static int sym_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    return (c >= 'A' && c <= 'Z') ? g_char_to_idx[c] : -1;
}

static const char *conv_name(int type) {
    return type == NIHILIST_SUB_NC ? "no-carry" : type == NIHILIST_SUB_M100 ? "mod-100" : "carry";
}

// Plant a Nihilist Substitution cipher: first pt_len plaintext letters under the keyed square
// (from sqkw) and the periodic additive key (from addkw, its letters' cells). conv selects the
// addition convention; rowlbl/collbl (NULL => fixed 1..5) select the labels. Fills prepared[]
// (expected solution, length n), cipher_nums[] (the numbers), and cipher_str (space-separated
// numbers). Returns n.
static int plant(const char *sqkw, const char *addkw, int conv,
                 const int *rowlbl_in, const int *collbl_in,
                 int pt_len, int prepared[], int cipher_nums[], char cipher_str[]) {
    int side = NIHILIST_SUB_SIDE, gs = NIHILIST_SUB_GRID;
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = sym_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    int kw[64], kwn = 0;
    for (int i = 0; sqkw[i] && kwn < 64; i++) {
        int idx = sym_to_index((unsigned char) sqkw[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int grid[NIHILIST_SUB_GRID];
    bifid_grid_from_keyword(kw, kwn, grid, gs);
    int pos[NIHILIST_SUB_GRID];
    bifid_build_inverse(grid, pos, gs);

    int key_cells[64], period = 0;
    for (int i = 0; addkw[i] && period < 64; i++) {
        int idx = sym_to_index((unsigned char) addkw[i]);
        if (idx >= 0) key_cells[period++] = pos[idx];
    }

    int rowlbl[NIHILIST_SUB_SIDE], collbl[NIHILIST_SUB_SIDE];
    nihilist_sub_fixed_labels(rowlbl, collbl, side);
    if (rowlbl_in) for (int i = 0; i < side; i++) rowlbl[i] = rowlbl_in[i];
    if (collbl_in) for (int i = 0; i < side; i++) collbl[i] = collbl_in[i];

    nihilist_sub_encrypt(prepared, n, grid, rowlbl, collbl, side, key_cells, period, conv, cipher_nums);

    int w = 0;
    for (int i = 0; i < n; i++) w += sprintf(cipher_str + w, "%s%d", i ? " " : "", cipher_nums[i]);
    return n;
}

// Run one solve (period pinned if period>0, else estimated with n_periods candidates); returns
// the recovered fraction. period_out gets the reported period.
static double solve_and_frac(int type, const char *cipher_str, const int prepared[], int plen,
        int period, int n_periods, int nr, int nh, double it, double bt, uint32_t seed,
        int *period_out, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (n_periods > 0) cfg.n_periods = n_periods;
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

    if (period_out) *period_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    int codes[3] = { NIHILIST_SUB, NIHILIST_SUB_NC, NIHILIST_SUB_M100 };
    for (int i = 0; i < 3; i++) {
        init_config(&cfg);
        cfg.cipher_type = codes[i]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "%s registry: no entry applied", conv_name(codes[i]));
        CHECK(cfg.n_restarts == 8 && cfg.n_hill_climbs == 300000,
            "%s anneal defaults wrong: %dx%d", conv_name(codes[i]), cfg.n_restarts, cfg.n_hill_climbs);
        CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
            "%s anneal inittemp wrong: %.4f", conv_name(codes[i]), cfg.init_temp);
    }
    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. period estimator top-K hit rate (per convention) ----------------------

static void test_period_estimator(int type) {
    int periods[] = {5, 6, 7, 9, 11};
    int lengths[] = {600, 350};
    int trials = 0, hits = 0, rank1 = 0;
    const char *addkw[] = {"ABCDE", "BERLIN", "MERCURY", "KEYWORDXY", "ABCDEFGHIKL"};
    for (int pi = 0; pi < 5; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH], cnums[MAX_CIPHER_LENGTH];
            static char cipher_str[8 * MAX_CIPHER_LENGTH];
            int p = periods[pi];
            int plen = plant("KRYPTOSABCD", addkw[pi], type, NULL, NULL, lengths[li],
                prepared, cnums, cipher_str);
            int out[8];
            int nn = nihilist_sub_estimate_periods(cnums, plen, 1, 15, 5, out, false);
            trials++;
            for (int i = 0; i < nn; i++) if (out[i] == p) { hits++; if (i == 0) rank1++; break; }
        }
    }
    printf("[%s period estimator] true period in top-5 for %d/%d cases (rank-1 for %d)\n",
        conv_name(type), hits, trials, rank1);
    CHECK(hits == trials, "%s estimator missed the true period in %d/%d cases",
        conv_name(type), trials - hits, trials);
}

// --- 3. capability floor + length cliff (per convention, period pinned) -------

static void test_capability_cliff(int type) {
    int lengths[] = {250, 400};
    printf("%s recovery vs length (sq=KRYPTOSABCD, key=BERLIN, period 6 pinned, 5x150k, seed 1):\n",
        conv_name(type));
    double frac_longest = 0.0;
    for (int li = 0; li < 2; li++) {
        static int prepared[MAX_CIPHER_LENGTH], cnums[MAX_CIPHER_LENGTH];
        static char cipher_str[8 * MAX_CIPHER_LENGTH];
        int plen = plant("KRYPTOSABCD", "BERLIN", type, NULL, NULL, lengths[li],
            prepared, cnums, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(type, cipher_str, prepared, plen, 6, 0,
            5, 150000, 0.08, 0.30, 1u, NULL, &secs);
        printf("    pt~%-4d  frac=%.3f  %.1fs\n", plen, frac, secs);
        if (li == 1) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.95, "%s 400ch (period pinned) recovered only %.3f",
        conv_name(type), frac_longest);
}

// --- 4. blind period (carry): estimated end to end ----------------------------

static void test_blind_period(void) {
    static int prepared[MAX_CIPHER_LENGTH], cnums[MAX_CIPHER_LENGTH];
    static char cipher_str[8 * MAX_CIPHER_LENGTH];
    int truep = 6;
    int plen = plant("KRYPTOSABCD", "BERLIN", NIHILIST_SUB, NULL, NULL, 480,
        prepared, cnums, cipher_str);
    int period = -1; double secs = 0.0;
    double frac = solve_and_frac(NIHILIST_SUB, cipher_str, prepared, plen, 0, 5,
        6, 150000, 0.08, 0.30, 1u, &period, &secs);
    printf("[carry blind period  %dch, estimated]  frac=%.3f period=%d (true %d)  %.1fs\n",
        plen, frac, period, truep, secs);
    CHECK(frac >= 0.95, "carry blind-period recovered only %.3f", frac);
    // IoC peaks at multiples of the true period too, so a multiple (e.g. 12, key BERLINBERLIN)
    // is an equally correct solve -- accept any multiple of the true period.
    CHECK(period > 0 && period % truep == 0, "carry blind period %d not a multiple of %d", period, truep);
}

// --- 5. multi-keyword sweep (carry) -------------------------------------------

static void test_multikeyword(void) {
    const char *sqkw[] = {"KRYPTOSABCD", "PLAYFAIREXAM", "ZEBRASCMBLD"};
    const char *adkw[] = {"BERLIN", "MERCURY", "PARISXY"};
    double sum = 0.0, worst = 1.0;
    printf("carry multi-keyword (400ch, period varies, pinned, 6x180k, seed 2):\n");
    for (int i = 0; i < 3; i++) {
        static int prepared[MAX_CIPHER_LENGTH], cnums[MAX_CIPHER_LENGTH];
        static char cipher_str[8 * MAX_CIPHER_LENGTH];
        int plen = plant(sqkw[i], adkw[i], NIHILIST_SUB, NULL, NULL, 400, prepared, cnums, cipher_str);
        int per = (int) strlen(adkw[i]);
        double secs = 0.0;
        double frac = solve_and_frac(NIHILIST_SUB, cipher_str, prepared, plen, per, 0,
            6, 180000, 0.08, 0.30, 2u, NULL, &secs);
        printf("    sq=%-13s key=%-8s frac=%.3f  %.1fs\n", sqkw[i], adkw[i], frac, secs);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("    mean=%.3f worst=%.3f\n", sum / 3.0, worst);
    CHECK(worst >= 0.95, "carry multi-keyword worst recovery only %.3f", worst);
}

// --- 6. keyed-label end to end (carry) ----------------------------------------

static void test_keyed_labels(void) {
    int rowlbl[5] = {3, 1, 4, 5, 2};
    int collbl[5] = {2, 5, 1, 3, 4};
    static int prepared[MAX_CIPHER_LENGTH], cnums[MAX_CIPHER_LENGTH];
    static char cipher_str[8 * MAX_CIPHER_LENGTH];
    int plen = plant("KRYPTOSABCD", "BERLIN", NIHILIST_SUB, rowlbl, collbl, 400,
        prepared, cnums, cipher_str);
    double secs = 0.0;
    double frac = solve_and_frac(NIHILIST_SUB, cipher_str, prepared, plen, 6, 0,
        5, 150000, 0.08, 0.30, 1u, NULL, &secs);
    printf("[carry keyed-label   %dch, period 6 pinned] frac=%.3f  %.1fs (solver assumes fixed labels)\n",
        plen, frac, secs);
    CHECK(frac >= 0.95, "carry keyed-label recovered only %.3f (as relabelled square)", frac);
}

int main(void) {
    g_ngram_logprob = true;             // Nihilist Substitution needs the log-probability fitness
    init_alphabet("J");
    CHECK(g_alpha == NIHILIST_SUB_GRID, "alphabet size %d, expected %d", g_alpha, NIHILIST_SUB_GRID);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    int codes[3] = { NIHILIST_SUB, NIHILIST_SUB_NC, NIHILIST_SUB_M100 };

    test_registry();
    for (int i = 0; i < 3; i++) test_period_estimator(codes[i]);
    for (int i = 0; i < 3; i++) test_capability_cliff(codes[i]);
    test_blind_period();
    test_multikeyword();
    test_keyed_labels();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
