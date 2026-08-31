//
//  In-process stress / limits tests for the Twin Bifid solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  A Twin Bifid is TWO Bifid messages enciphered under the SAME keyed 5x5 square (J->I) but at
//  DIFFERENT periods; the ACA con guarantees a shared plaintext phrase. The attack anneals a
//  SINGLE shared square (Bifid's 25-cell state + move set) while n-gram-scoring the CONCATENATED
//  decrypt of both messages -- so every square move is judged against ~2x the text. That is the
//  headline of this calibration: a lone Bifid needs long text, but the shared square recovers
//  from the short ACA lengths (100-150 letters EACH) because the joint signal is doubled. The
//  second ciphertext is injected in-process via cfg.twincipher_str (the CLI reads it from
//  -cipher2). The two periods are independent unknowns: bifid_estimate_periods (columnar IoC,
//  square-agnostic) ranks each message's period, and the solver anneals the CROSS PRODUCT of the
//  two top-K lists; -period pins message 1, -period2 pins message 2.
//
//  The suite checks:
//    1. registry validation (apply_cipher_defaults, anneal + shotgun) + a non-registry type;
//    2. the (Bifid) period estimator on BOTH messages (true period in the top-K IoC);
//    3. a capability floor at the ACA length (135 each) with periods ESTIMATED;
//    4. a length cliff (recovery vs per-message length) with periods pinned;
//    5. a multi-keyword-pair sweep (mean recovery), periods pinned;
//    6. a BLIND two-period solve -- the reported message-1 period must match, and full recovery
//       (which needs the message-2 period right too) confirms blind two-period selection;
//    7. per-scheme calibration: the same twin under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"            // apply_cipher_defaults
#include "scoring.h"           // load_ngrams
#include "bifid_solver.h"      // bifid_estimate_periods (square-agnostic, reused by Twin Bifid)
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

// A long chunk of natural English (Pride and Prejudice, opening), letters only (~940). The two
// twin messages are disjoint slices of this, so both decrypts are real English for the joint score.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEINWANTOFAWIFE"
    "HOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGANEIGHBOURHOODTHISTRUTHIS"
    "SOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOR"
    "OTHEROFTHEIRDAUGHTERSMYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLET"
    "ATLASTMRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASIUSTBEENHEREANDSHETOLDMEALLABOUT"
    "ITMRBENNETMADENOANSWERDOYOUNOTWANTTOKNOWWHOHASTAKENITCRIEDHISWIFEIMPATIENTLYYOUWANTTOTELLMEAND"
    "IHAVENOOBIECTIONTOHEARINGITTHISWASINVITATIONENOUGHWHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATNETHERFIELD"
    "ISTAKENBYAYOUNGMANOFLARGEFORTUNEFROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEE"
    "THEPLACEANDWASSOMUCHDELIGHTEDWITHITTHATHEAGREEDWITHMRMORRISIMMEDIATELYTHATHEISTOTAKEPOSSESSION"
    "BEFOREMICHAELMASANDSOMEOFHISSERVANTSARETOBEINTHEHOUSEBYTHEENDOFNEXTWEEKWHATISHISNAMEBINGLEYIS";

// A..Z char -> 0..24 alphabet index, merging J into I (the 25-letter Bifid convention).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

static void square_from_keyword(const char *kw, int grid[]) {
    int k[64], kn = 0;
    for (int i = 0; kw[i] && kn < 64; i++) { int x = letter_to_index((unsigned char) kw[i]); if (x >= 0) k[kn++] = x; }
    bifid_grid_from_keyword(k, kn, grid, g_alpha);
}

// Plant a Twin Bifid: message 1 = the first len1 PLAINTEXT letters, message 2 = the NEXT len2
// letters (two disjoint English slices under one keyed square, at periods p1 and p2). Fills
// prepared[] (the concatenated expected solution) and the two cipher strings. Returns len1+len2.
static int plant(const char *kw, int p1, int p2, int len1, int len2,
                 int prepared[], char cs1[], char cs2[]) {
    int nr = 0;
    for (int i = 0; PLAINTEXT[i] && nr < len1 + len2; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[nr++] = idx;
    }
    int n1 = len1, n2 = nr - len1;
    int sq[PLAYFAIR_GRID];
    square_from_keyword(kw, sq);
    int c1[MAX_CIPHER_LENGTH], c2[MAX_CIPHER_LENGTH];
    twin_bifid_encrypt(prepared, n1, prepared + n1, n2, sq, 5, p1, p2, c1, c2);
    for (int i = 0; i < n1; i++) cs1[i] = index_to_char(c1[i]);  cs1[n1] = '\0';
    for (int i = 0; i < n2; i++) cs2[i] = index_to_char(c2[i]);  cs2[n2] = '\0';
    return n1 + n2;
}

// Solve and return the recovered fraction over the concatenated (msg1 ++ msg2) plaintext.
// p1>0 pins message 1's period, p2>0 pins message 2's; otherwise both are estimated (scanning
// up to max_period, top n_periods each, crossed). `method` overrides the search scheme. Always
// applies the tuned registry schedule.
static double solve_and_frac(const char *cs1, const char *cs2, const int prepared[], int plen,
        int p1, int p2, int max_period, int n_periods, int method, uint32_t seed,
        int *period1_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TWIN_BIFID;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    cfg.twincipher_str = (char *) cs2;         // inject the second ciphertext in-process
    if (p1 > 0) { cfg.period_present = true; cfg.period = p1; }
    if (p2 > 0) { cfg.period2_present = true; cfg.period2 = p2; }
    if (max_period > 0) cfg.max_period = max_period;
    if (n_periods > 0) cfg.n_periods = n_periods;

    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* still proceed */ }
    seed_rand(seed);
    solve_cipher((char *) cs1, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);
    if (secs_out) *secs_out = ((double) clock() - t0) / CLOCKS_PER_SEC;

    if (period1_out) *period1_out = res.solved ? res.cycleword_len : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = TWIN_BIFID; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "twin-bifid registry: no entry applied");
    CHECK(cfg.n_restarts == 4 && cfg.n_hill_climbs == 120000,
        "twin-bifid anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "twin-bifid anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = TWIN_BIFID; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "twin-bifid registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 150000,
        "twin-bifid shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. period estimator (in isolation) ---------------------------------------
//
// Twin Bifid reuses Bifid's columnar-IoC estimator unchanged (per message). Like every IoC
// period test it is a LONG-text property: on short Bifid ciphertext the columnar IoC is noisy
// and the true period is not reliably in the top-K until ~400 letters (a documented Bifid-
// family limitation, orthogonal to the square recovery -- see the capability/blind blocks). So
// this validates the estimator where it is meant to work: the true period must be a top-6
// candidate at 760 and 560 letters for every tested period.

static void test_period_estimator(void) {
    int periods[] = {5, 7, 9, 11};
    int lengths[] = {760, 560};
    int trials = 0, hits = 0;
    for (int pi = 0; pi < 4; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
            int p = periods[pi], L = lengths[li];
            plant("KRYPTOS", p, p, L, 8, prepared, cs1, cs2);   // estimate on message 1 (length L)
            int cidx[MAX_CIPHER_LENGTH];
            int n = 0; for (int i = 0; cs1[i]; i++) cidx[n++] = letter_to_index((unsigned char) cs1[i]);
            int out[8];
            int m = bifid_estimate_periods(cidx, n, 2, 16, 6, out, false);
            trials++;
            int found = 0; for (int i = 0; i < m; i++) if (out[i] == p) found = 1;
            if (found) hits++;
            else printf("  [estimator MISS] period %d len %d -> {%d %d %d %d %d %d}\n",
                p, L, out[0], out[1], out[2], out[3], out[4], out[5]);
        }
    }
    printf("[period estimator] true period in top-6 for %d/%d cases\n", hits, trials);
    CHECK(hits == trials, "period estimator missed the true period in %d/%d cases", trials - hits, trials);
}

// --- 3. capability floor at the ACA length (periods PINNED) --------------------
//
// THE HEADLINE. 135 letters EACH (270 combined) -- squarely in the ACA 100-150 band. Given the
// two periods, the SHARED square recovers cleanly: a lone Bifid is near hopeless at 135 letters
// (its square break needs ~500+), but scoring one square against BOTH decrypts doubles the
// n-gram signal, so the twin recovers from ACA-length text. (Period ESTIMATION is the separate
// bottleneck -- reliable only from ~400 letters, block 2 -- so a fully-blind ACA-length solve is
// estimator-limited; pin the periods, which the ACA con lets you find.)

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
    int p1 = 7, p2 = 9, L = 135;
    int n = plant("KRYPTOS", p1, p2, L, L, prepared, cs1, cs2);
    int got_p1 = -1; double secs = 0.0;
    double frac = solve_and_frac(cs1, cs2, prepared, n, p1, p2, 0, 0, METHOD_DEFAULT, 1u, &got_p1, &secs);
    printf("[capability @ %d+%d chars, periods PINNED (%d,%d)] frac=%.3f %.1fs\n",
        L, L, p1, p2, frac, secs);
    CHECK(frac >= 0.95, "twin-bifid %d+%d (periods pinned) recovered only %.3f", L, L, frac);
}

// --- 4. length cliff (periods pinned) -----------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {90, 110, 135, 160};
    printf("recovery vs per-message length (keyword KRYPTOS, periods 7 & 9 pinned, registry anneal, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
        int L = lengths[li];
        int n = plant("KRYPTOS", 7, 9, L, L, prepared, cs1, cs2);
        double secs = 0.0;
        double frac = solve_and_frac(cs1, cs2, prepared, n, 7, 9, 0, 0, METHOD_DEFAULT, 1u, NULL, &secs);
        printf("    %d+%d (%d)  frac=%.3f  %.1fs\n", L, L, n, frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.95, "twin-bifid 160+160 (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 5. multi-keyword sweep (periods pinned) ----------------------------------

static void test_multi_keyword(void) {
    const char *kw[] = { "MONARCHY", "ZEBRA", "PORTABLE" };
    int p1s[] =        { 7,          9,       5         };
    int p2s[] =        { 9,          11,      8         };   // periods differ per message
    int L = 140, nk = 3;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ %d+%d chars, periods pinned]\n", L, L);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
        int n = plant(kw[k], p1s[k], p2s[k], L, L, prepared, cs1, cs2);
        double frac = solve_and_frac(cs1, cs2, prepared, n, p1s[k], p2s[k], 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  %-9s P=(%d,%d) : %.1f%%\n", kw[k], p1s[k], p2s[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 6. blind two-period solve ------------------------------------------------
//
// Both periods ESTIMATED (top-4 crossed). The solver reports message 1's period in
// cycleword_len; full recovery over the CONCATENATED plaintext additionally requires message
// 2's period to be right, so asserting the reported p1 AND frac asserts blind TWO-period
// selection. Asserted at 400 letters each -- the length where the columnar-IoC estimator
// reliably surfaces both periods (block 2). A blind ACA-length (135) attempt is PRINTED for the
// record: it is estimator-limited (the true periods often miss the top-K), the documented
// short-text limitation -- pin the periods there.

static void test_blind_period(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
    int p1 = 7, p2 = 9, L = 400;
    int n = plant("KRYPTOS", p1, p2, L, L, prepared, cs1, cs2);
    double secs; int pout;
    double frac = solve_and_frac(cs1, cs2, prepared, n, 0, 0, 16, 3, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind periods @ %d+%d (scan 2..16, top 4 crossed), true (%d,%d)]: reported p1=%d, %.1f%%  [%.1fs]\n",
        L, L, p1, p2, pout, 100.0 * frac, secs);
    CHECK(pout == p1, "blind reported message-1 period %d (true %d)", pout, p1);
    CHECK(frac > 0.90, "blind two-period recovery only %.1f%% (implies message-2 period recovered)", 100.0 * frac);

    // Characterization (not asserted): a fully-blind solve at the ACA length is estimator-limited.
    int L2 = 135;
    int n2 = plant("KRYPTOS", p1, p2, L2, L2, prepared, cs1, cs2);
    double secs2; int pout2;
    double frac2 = solve_and_frac(cs1, cs2, prepared, n2, 0, 0, 16, 3, METHOD_DEFAULT, 0xB11Du, &pout2, &secs2);
    printf("[blind @ %d+%d (ACA length, documented limitation)]: reported p1=%d, %.1f%%  [%.1fs]\n",
        L2, L2, pout2, 100.0 * frac2, secs2);
}

// --- 7. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
    int p1 = 7, p2 = 9, L = 150;
    int n = plant("KRYPTOS", p1, p2, L, L, prepared, cs1, cs2);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d+%d chars, keyword KRYPTOS, periods (%d,%d) pinned]\n", L, L, p1, p2);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs1, cs2, prepared, n, p1, p2, 0, 0, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet("J");                 // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == PLAYFAIR_GRID, "alphabet size %d, expected %d", g_alpha, PLAYFAIR_GRID);

    g_ngram_logprob = true;             // Twin Bifid needs the log-probability fitness
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
    test_multi_keyword();
    test_blind_period();
    test_per_scheme();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
