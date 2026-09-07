//
//  In-process stress / limits tests for the CM Bifid solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  CM Bifid (Conjugated Matrix Bifid) is Bifid fractionation over TWO independently keyed 5x5
//  squares: square 1 fractionates the plaintext into coords and square 2 recombines the
//  re-paired coordinate stream into ciphertext letters. The attack is a JOINT two-square
//  anneal (the proven Two/Four-Square break) -- the state is the pair of squares packed back-
//  to-back, each move perturbing ONE square. There is NO square-independent decoupling reward
//  (both squares are entangled in the n-gram fitness), so it is a genuine joint search and,
//  like every square type, effectively needs the log-probability fitness -- this suite enables
//  g_ngram_logprob and loads quadgrams in that mode. The period is recovered by the SAME
//  columnar-IoC estimator as Bifid (bifid_estimate_periods is square-agnostic, reused), swept
//  on top (one engine config per candidate period; the n-gram score picks it).
//
//  *** ODD-PERIOD PROPERTY (the headline of the calibration). *** With an EVEN period P the
//  rows-then-cols re-paired stream splits cleanly -- the first P/2 output pairs are pure-ROW
//  coordinates and the last P/2 pure-COLUMN -- so a fractionation row and column never share an
//  output pair, leaving a transpose-like symmetry that makes the squares ambiguous CIPHERTEXT-
//  ONLY (a different square pair decrypts to equally-English text). Recovery of the PLANTED key
//  then fails even with a huge budget / long text -- not a search shortfall but a genuine
//  degeneracy (measured: even P recovers ~5%, the noise floor). With an ODD period the boundary
//  pair MIXES a row and a column, breaking the symmetry, and recovery is clean (~100%). So the
//  suite asserts on ODD periods and documents the even-period degeneracy explicitly. Empirically
//  the odd-period cliff sits near ~400 letters and recovery is reliable from ~480+; the longer
//  per-restart climb is the critical lever (more restarts with shorter climbs does WORSE near
//  the cliff), which is why the registry schedule is 8x400000. The suite checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the (Bifid) PERIOD ESTIMATOR in isolation (true period in the top-K columnar-IoC);
//    3. a capability floor at ~520 chars (odd P) + the ODD-vs-EVEN contrast on one cipher;
//    4. a length cliff (recovery vs length) showing the ~400 cliff / ~480 floor (odd P);
//    5. a multi-keyword-pair sweep (mean/worst recovery), odd periods pinned;
//    6. a BLIND period solve (P estimated over a bounded range) -- the reported P must match;
//    7. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"            // apply_cipher_defaults
#include "scoring.h"           // load_ngrams
#include "bifid_solver.h"      // bifid_estimate_periods (square-agnostic, reused by CM Bifid)
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

// A long chunk of natural English (Pride and Prejudice, opening), letters only (~940).
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

// A..Z char -> 0..24 alphabet index, merging J into I (the 25-letter CM Bifid convention).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Build a keyed square (alphabet indices) from a keyword string, dropping non-letters.
static void square_from_keyword(const char *kw, int grid[]) {
    int k[64], kn = 0;
    for (int i = 0; kw[i] && kn < 64; i++) { int x = letter_to_index((unsigned char) kw[i]); if (x >= 0) k[kn++] = x; }
    bifid_grid_from_keyword(k, kn, grid, g_alpha);
}

// Plant a CM Bifid cipher: take the first pt_len letters of PLAINTEXT (J->I), build square 1
// from kw1 and square 2 from kw2, and encipher with `period`. Fills prepared[] (the expected
// solution) and cipher_str. Returns the prepared length.
static int plant(const char *kw1, const char *kw2, int period, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    int sq1[PLAYFAIR_GRID], sq2[PLAYFAIR_GRID];
    square_from_keyword(kw1, sq1);
    square_from_keyword(kw2, sq2);

    int cipher[MAX_CIPHER_LENGTH];
    cm_bifid_encrypt(prepared, n, sq1, sq2, 5, period, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Solve and return the recovered fraction. period>0 pins P; otherwise the period is
// estimated, scanning 2..max_period and annealing the top n_periods candidates. `method`
// overrides the search scheme. Always applies the tuned per-type registry schedule.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int period, int max_period, int n_periods, int method, uint32_t seed,
        int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = CM_BIFID;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0) { cfg.period_present = true; cfg.period = period; }
    if (max_period > 0) cfg.max_period = max_period;
    if (n_periods > 0) cfg.n_periods = n_periods;

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
    init_config(&cfg); cfg.cipher_type = CM_BIFID; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "cm-bifid registry: no entry applied");
    CHECK(cfg.n_restarts == 8 && cfg.n_hill_climbs == 400000,
        "cm-bifid anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "cm-bifid anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = CM_BIFID; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "cm-bifid registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 20 && cfg.n_hill_climbs == 300000,
        "cm-bifid shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. period estimator (in isolation) ---------------------------------------
//
// CM Bifid reuses Bifid's columnar-IoC estimator unchanged -- square 2 only relabels the
// coordinate pairs and the IoC is relabel-invariant, so the peak still lands on the true
// period (and its multiples, which the n-gram score discards). The estimator works for BOTH
// parities (it ranks the period; the EVEN-period weakness is in key RECOVERY, not estimation),
// so this checks odd and even alike. The true period must be a top-5 candidate every time.

static void test_period_estimator(void) {
    int periods[] = {5, 6, 7, 9};
    int lengths[] = {760, 560};
    int trials = 0, hits = 0, rank1 = 0;
    for (int pi = 0; pi < 4; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cs[MAX_CIPHER_LENGTH];
            int p = periods[pi];
            int plen = plant("KRYPTOS", "MACHINE", p, lengths[li], prepared, cs);
            int cidx[MAX_CIPHER_LENGTH];
            for (int i = 0; i < plen; i++) cidx[i] = letter_to_index((unsigned char) cs[i]);
            int out[8];
            int n = bifid_estimate_periods(cidx, plen, 2, 14, 5, out, false);
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
    CHECK(hits == trials, "period estimator missed the true period in %d/%d cases",
        trials - hits, trials);
}

// --- 3. capability floor (odd P) + the ODD-vs-EVEN contrast --------------------
//
// Run at ~520 chars -- comfortably above the ~400 odd-period cliff. Then re-encipher the SAME
// plaintext/keys at an EVEN period and show recovery collapses to the noise floor: the
// documented even-period square ambiguity (NOT asserted, since an occasional even-period key
// still resolves -- it is fragile, not impossible -- but printed so the property is on record).

static void test_capability_floor(void) {
    int plen = 520;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    double secs;

    printf("\n[capability floor @ ~%d chars, KRYPTOS/MACHINE, ODD period 7]\n", plen);
    int n = plant("KRYPTOS", "MACHINE", 7, plen, prepared, cs);
    double frac = solve_and_frac(cs, prepared, n, 7, 0, 0, METHOD_DEFAULT, 0xC0FFEEu, NULL, &secs);
    printf("  P=7 (odd)  : %.1f%%  [%.1fs]\n", 100.0 * frac, secs);
    CHECK(frac > 0.95, "cm-bifid floor (odd P): only %.1f%% at %d chars", 100.0 * frac, n);

    int n6 = plant("KRYPTOS", "MACHINE", 6, plen, prepared, cs);
    double frac6 = solve_and_frac(cs, prepared, n6, 6, 0, 0, METHOD_DEFAULT, 0xC0FFEEu, NULL, &secs);
    printf("  P=6 (even) : %.1f%%  [%.1fs]  <- even-period square ambiguity (documented, not asserted)\n",
        100.0 * frac6, secs);
}

// --- 4. length cliff (odd P) --------------------------------------------------
//
// Show recovery climbing across 320 -> 480 -> 640 (odd period 7): ~320 is below the cliff
// (documented), the floor is reliable from ~480.

static void test_length_cliff(void) {
    int lens[] = { 320, 480, 640 };
    printf("\n[length cliff: keyword KRYPTOS/MACHINE, P=7 (odd), registry anneal -- cliff ~400, floor ~480]\n");
    double frac480 = 0.0, frac640 = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("KRYPTOS", "MACHINE", 7, lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 7, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 1) frac480 = frac;
        if (li == 2) frac640 = frac;
    }
    CHECK(frac640 >= 0.95, "cm-bifid 640ch recovered only %.3f", frac640);
    CHECK(frac480 >= 0.95, "cm-bifid 480ch (above the cliff) recovered only %.3f", frac480);
}

// --- 5. multi-keyword sweep (odd periods pinned) ------------------------------

static void test_multi_keyword(void) {
    const char *kw1[] = { "MONARCHY", "ZEBRA",      "PORTABLE" };
    const char *kw2[] = { "SHADOW",   "PALMERSTON", "GADGETRY" };
    int periods[] =     { 7,          5,            9          };   // all ODD
    int plen = 600, nk = 3;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars, odd periods pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kw1[k], kw2[k], periods[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, periods[k], 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  %-9s/%-11s P=%d : %.1f%%\n", kw1[k], kw2[k], periods[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 6. blind period solve (period estimated over a bounded range) ------------

static void test_blind_period(void) {
    int P = 7, plen = 600;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("KRYPTOS", "MACHINE", P, plen, prepared, cs);
    double secs; int pout;
    // Bound the scan to 2..9 (so 2*P=14 is excluded) and anneal the top 3 candidates so the
    // multi-config two-square anneal stays tractable in CI; the true period reliably tops the
    // columnar IoC there and the n-gram score discards the degenerate even-period peaks.
    double frac = solve_and_frac(cs, prepared, n, 0, 9, 3, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind P (scan 2..9, top 3), true P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        P, pout, 100.0 * frac, secs);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
    CHECK(frac > 0.90, "blind-P recovery only %.1f%%", 100.0 * frac);
}

// --- 7. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int P = 7, plen = 520;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("KRYPTOS", "MACHINE", P, plen, prepared, cs);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ ~%d chars, keyword KRYPTOS/MACHINE, P=%d (odd) pinned]\n", plen, P);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, 0, 0, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet("J");                 // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == PLAYFAIR_GRID, "alphabet size %d, expected %d", g_alpha, PLAYFAIR_GRID);

    g_ngram_logprob = true;             // CM Bifid needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_period_estimator();
    test_capability_floor();
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
