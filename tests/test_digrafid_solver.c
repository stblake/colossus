//
//  In-process stress / limits tests for the Digrafid solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Digrafid is a digraphic fractionation cipher over TWO independently keyed 27-symbol
//  alphabets (A..Z + '#'): a horizontal 3x9 grid and a vertical 9x3 grid. ACA Digrafid grids
//  are KEYED ALPHABETS (a keyword, duplicates dropped, then the rest of the alphabet
//  ascending), so the solver searches that structured space -- the state is two keyed-alphabet
//  SEQUENCES, each "keyword prefix of length kw + sorted tail", perturbed by structure-
//  preserving keyword moves (digrafid_move_seq) rather than the old free 54-cell square break.
//  The keyed-alphabet prior collapses the keyspace by orders of magnitude, so recovery no
//  longer needs ~700-800 letters (the free-permutation cliff) -- it tracks the keyword and
//  recovers reliably from ~300 letters, the way an ACA solver does. A fractionation period is
//  SWEPT on top (one engine config per candidate period; the n-gram score picks it). It still
//  effectively needs the log-probability fitness, so this suite enables g_ngram_logprob and
//  loads quadgrams in that mode. The suite (also the basis for tuning the SearchDefaults
//  48x150000 / inittemp-0.30 schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the keyed-alphabet move INVARIANT (digrafid_move_seq keeps a permutation + sorted tail);
//    3. the PERIOD ESTIMATOR in isolation (true period in the top-K per-lane-IoC candidates);
//    4. a capability floor (recovery %) across keywords at ~420 chars -- BELOW the old cliff;
//    5. a length cliff (recovery vs length) showing the ~300-char floor (vs the old ~700);
//    6. a multi-keyword sweep (mean/worst recovery), period pinned;
//    7. a BLIND period solve (P estimated over a bounded range) -- the reported P must match;
//    8. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"            // apply_cipher_defaults
#include "scoring.h"           // load_ngrams
#include "digrafid_solver.h"   // digrafid_estimate_periods
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
    "ATLASTMRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASJUSTBEENHEREANDSHETOLDMEALLABOUT"
    "ITMRBENNETMADENOANSWERDOYOUNOTWANTTOKNOWWHOHASTAKENITCRIEDHISWIFEIMPATIENTLYYOUWANTTOTELLMEAND"
    "IHAVENOOBJECTIONTOHEARINGITTHISWASINVITATIONENOUGHWHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATNETHERFIELD"
    "ISTAKENBYAYOUNGMANOFLARGEFORTUNEFROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEE"
    "THEPLACEANDWASSOMUCHDELIGHTEDWITHITTHATHEAGREEDWITHMRMORRISIMMEDIATELYTHATHEISTOTAKEPOSSESSION"
    "BEFOREMICHAELMASANDSOMEOFHISSERVANTSARETOBEINTHEHOUSEBYTHEENDOFNEXTWEEKWHATISHISNAMEBINGLEYIS";

// A..Z char -> 0..25 alphabet index (Digrafid keeps the full 26 letters; no J->I merge).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Digrafid cipher: take the first pt_len (rounded down to even) letters of PLAINTEXT,
// build the horizontal grid from kwH (row-major) and the vertical from kwV (column-major),
// and encipher with `period`. Fills prepared[] (the expected solution) and cipher_str.
// Returns the prepared length.
static int plant(const char *kwH, const char *kwV, int period, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    n &= ~1;                                  // even length (whole digraphs)

    int kw[64], kn;
    int gridH[DIGRAFID_GRID], gridV[DIGRAFID_GRID];
    kn = 0; for (int i = 0; kwH[i] && kn < 64; i++) { int x = letter_to_index((unsigned char) kwH[i]); if (x >= 0) kw[kn++] = x; }
    digrafid_grid_from_keyword(kw, kn, gridH, DIGRAFID_HROWS, DIGRAFID_HCOLS, 0);
    kn = 0; for (int i = 0; kwV[i] && kn < 64; i++) { int x = letter_to_index((unsigned char) kwV[i]); if (x >= 0) kw[kn++] = x; }
    digrafid_grid_from_keyword(kw, kn, gridV, DIGRAFID_VROWS, DIGRAFID_VCOLS, 1);

    int cipher[MAX_CIPHER_LENGTH];
    digrafid_encrypt(prepared, n, gridH, gridV, period, cipher);
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
    cfg.cipher_type = DIGRAFID;
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
    init_config(&cfg); cfg.cipher_type = DIGRAFID; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "digrafid registry: no entry applied");
    CHECK(cfg.n_restarts == 48 && cfg.n_hill_climbs == 150000,
        "digrafid anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "digrafid anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = DIGRAFID; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "digrafid registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 200 && cfg.n_hill_climbs == 150000,
        "digrafid shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. keyed-alphabet move invariant -----------------------------------------
//
// The whole keyed-alphabet search rests on digrafid_move_seq preserving the invariant
// "seq is a permutation of 0..26 AND seq[kw..26] is strictly ascending" (and kw stays in
// [DIGRAFID_KW_MIN .. DIGRAFID_KW_MAX]). Stress it: from many random keyed seeds, apply a
// long chain of moves and assert the invariant holds at every step.

static void test_move_invariant(void) {
    int bad_perm = 0, bad_tail = 0, bad_kw = 0;
    seed_rand(0xA11CEu);
    for (int trial = 0; trial < 400; trial++) {
        int seq[DIGRAFID_GRID];
        int kw = rand_int(3, 14);                       // DIGRAFID_KW_MIN..MAX
        random_keyword(seq, DIGRAFID_GRID, kw);
        for (int step = 0; step < 500; step++) {
            digrafid_move_seq(seq, DIGRAFID_GRID, &kw);
            // permutation of 0..26?
            int seen[DIGRAFID_GRID] = {0};
            for (int i = 0; i < DIGRAFID_GRID; i++) {
                if (seq[i] < 0 || seq[i] >= DIGRAFID_GRID || seen[seq[i]]) { bad_perm++; break; }
                seen[seq[i]] = 1;
            }
            // tail strictly ascending past kw?
            for (int i = kw; i + 1 < DIGRAFID_GRID; i++)
                if (seq[i] >= seq[i + 1]) { bad_tail++; break; }
            // kw in range?
            if (kw < 3 || kw > 13) bad_kw++;
        }
    }
    printf("[move invariant] 400 chains x 500 moves: perm-violations=%d tail-violations=%d kw-violations=%d\n",
        bad_perm, bad_tail, bad_kw);
    CHECK(bad_perm == 0, "digrafid_move_seq broke the permutation %d times", bad_perm);
    CHECK(bad_tail == 0, "digrafid_move_seq broke the sorted tail %d times", bad_tail);
    CHECK(bad_kw == 0, "digrafid_move_seq drove kw out of [3,13] %d times", bad_kw);
}

// --- 3. period estimator (in isolation) ---------------------------------------

static void test_period_estimator(void) {
    int periods[] = {4, 5, 6, 7};
    int lengths[] = {760, 560};
    int trials = 0, hits = 0, rank1 = 0;
    for (int pi = 0; pi < 4; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cs[MAX_CIPHER_LENGTH];
            int p = periods[pi];
            int plen = plant("CIPHER", "MACHINE", p, lengths[li], prepared, cs);
            int cidx[MAX_CIPHER_LENGTH];
            for (int i = 0; i < plen; i++) cidx[i] = letter_to_index((unsigned char) cs[i]);
            int out[8];
            int n = digrafid_estimate_periods(cidx, plen, 2, 12, 5, out, false);
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
    // The true period must be a top-5 candidate every time -- that is all the solver needs
    // (it anneals all five and the n-gram score discards the multiples-of-period peaks).
    CHECK(hits == trials, "period estimator missed the true period in %d/%d cases",
        trials - hits, trials);
}

// --- 4. capability floor (period pinned) --------------------------------------
//
// Run at ~420 chars -- comfortably BELOW the old free-permutation cliff (~700-800), where the
// keyed-alphabet search is the whole reason a solve is possible. Two keyword pairs spanning
// short (CIPHER=6) and long (BERLINCLOCK=11, near DIGRAFID_KW_MAX) keywords.

static void test_capability_floor(void) {
    const char *kwH[] = { "CIPHER",  "KRYPTOS" };
    const char *kwV[] = { "MACHINE", "BERLINCLOCK" };
    int periods[] =     { 5,         7 };
    int plen = 420;
    printf("\n[capability floor @ ~%d chars (below the old ~700 cliff), period pinned]\n", plen);
    for (int k = 0; k < 2; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kwH[k], kwV[k], periods[k], plen, prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, periods[k], 0, 0, METHOD_DEFAULT,
            0xC0FFEEu + 17 * k, NULL, &secs);
        printf("  %-8s/%-12s P=%d : %.1f%%  [%.1fs]\n", kwH[k], kwV[k], periods[k], 100.0 * frac, secs);
        CHECK(frac > 0.95, "digrafid(%s/%s) floor: only %.1f%% at %d chars",
            kwH[k], kwV[k], 100.0 * frac, n);
    }
}

// --- 5. length cliff (period pinned) ------------------------------------------
//
// The headline of the keyed-alphabet rewrite: the cliff moved from ~700 (free permutation) to
// ~300. Show recovery climbing across 240 -> 320 -> 480 and assert a clean solve at 320+.

static void test_length_cliff(void) {
    int lens[] = { 240, 320, 480 };
    printf("\n[length cliff: keyword CIPHER/MACHINE, P=5, registry anneal -- floor ~300, was ~700]\n");
    double frac320 = 0.0, frac480 = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("CIPHER", "MACHINE", 5, lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 5, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, NULL, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 1) frac320 = frac;
        if (li == 2) frac480 = frac;
    }
    CHECK(frac480 >= 0.95, "digrafid 480ch recovered only %.3f", frac480);
    CHECK(frac320 >= 0.95, "digrafid 320ch (just above the new cliff) recovered only %.3f", frac320);
}

// --- 6. multi-keyword sweep (period pinned) -----------------------------------

static void test_multi_keyword(void) {
    const char *kwH[] = { "ZEBRA",      "MONARCHY", "PORTABLE" };
    const char *kwV[] = { "PALMERSTON", "SHADOW",   "GADGETRY" };
    int periods[] =     { 5,            6,          8 };
    int plen = 480, nk = 3;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars, period pinned]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kwH[k], kwV[k], periods[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, periods[k], 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL, NULL);
        printf("  %-8s/%-11s P=%d : %.1f%%\n", kwH[k], kwV[k], periods[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 7. blind period solve (period estimated over a bounded range) ------------

static void test_blind_period(void) {
    int P = 5, plen = 480;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", "MACHINE", P, plen, prepared, cs);
    double secs; int pout;
    // Bound the scan to 2..6 and anneal the top 3 candidates so the multi-config full-grid
    // anneal stays tractable in CI (the true period reliably tops the per-lane IoC there).
    double frac = solve_and_frac(cs, prepared, n, 0, 6, 3, METHOD_DEFAULT, 0xB11Du, &pout, &secs);
    printf("\n[blind P (scan 2..6, top 3), true P=%d]: reported P=%d, %.1f%%  [%.1fs]\n",
        P, pout, 100.0 * frac, secs);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
    CHECK(frac > 0.90, "blind-P recovery only %.1f%%", 100.0 * frac);
}

// --- 8. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int P = 5, plen = 360;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", "MACHINE", P, plen, prepared, cs);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ ~%d chars, keyword CIPHER/MACHINE, P=%d pinned]\n", plen, P);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, P, 0, 0, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet_digrafid();           // 27-symbol Digrafid alphabet (A..Z + '#')
    CHECK(g_alpha == DIGRAFID_GRID, "alphabet size %d, expected %d", g_alpha, DIGRAFID_GRID);

    g_ngram_logprob = true;             // Digrafid needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_move_invariant();
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
