//
//  In-process stress / limits tests for the Twin Trifid solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  A Twin Trifid is TWO Trifid messages enciphered under the SAME keyed 3x3x3 cube (27 symbols,
//  A..Z + '+') but at DIFFERENT periods; the ACA con guarantees a shared plaintext phrase. The
//  attack anneals a SINGLE shared cube (Trifid's 27-cell state + move set) while n-gram-scoring
//  the CONCATENATED decrypt of both messages -- so every cube move is judged against ~2x the
//  text. The second ciphertext is injected in-process via cfg.twincipher_str (the CLI reads it
//  from -cipher2).
//
//  *** QUINTGRAMS + LENGTH (the headline of the calibration). *** The 27-cell cube is a much
//  rougher landscape than Twin Bifid's 25-cell square (single-message Trifid needs ~500+
//  letters), so even with the shared cube's doubled signal it EFFECTIVELY NEEDS QUINTGRAMS
//  (quadgrams do not carry enough signal at the ACA lengths -- verified below) and is reliable
//  from ~210 letters EACH. That is still well below a lone Trifid's floor -- the twin advantage
//  -- but above Twin Bifid's ~110. So this suite loads QUINTGRAMS and pins the periods; blind
//  two-period recovery is a documented limitation (the columnar-IoC estimator is unreliable on
//  short Trifid ciphertext -- same as Bifid, block 2 -- and a blind cube sweep multiplies the
//  already-heavy per-config budget over every period pair), so pin the periods, which the ACA
//  con lets you find. The suite checks:
//    1. registry validation (apply_cipher_defaults, anneal + shotgun) + a non-registry type;
//    2. the (Trifid) period estimator on a message (true period in the top-K IoC, long text);
//    3. a capability floor at 210 letters each (periods pinned);
//    4. a length cliff (recovery vs per-message length) showing the ~210 floor;
//    5. per-scheme calibration: the same twin under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram tables are found in the cwd.
//

#include "colossus.h"
#include "engine.h"            // apply_cipher_defaults
#include "scoring.h"           // load_ngrams
#include "trifid_solver.h"     // trifid_estimate_periods
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define QUINT_FILE "ngram_data/english/english_quintgrams.txt"

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

// A..Z char -> alphabet index (no J merge; the Trifid alphabet is A..Z + '+' and '+' is a cube
// symbol, not used in plaintext).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

static void cube_from_keyword(const char *kw, int cube[]) {
    int k[64], kn = 0;
    for (int i = 0; kw[i] && kn < 64; i++) { int x = letter_to_index((unsigned char) kw[i]); if (x >= 0) k[kn++] = x; }
    trifid_cube_from_keyword(k, kn, cube, g_alpha);
}

// Plant a Twin Trifid: message 1 = the first len1 PLAINTEXT letters, message 2 = the NEXT len2
// letters (two disjoint English slices under one keyed cube, at periods p1 and p2). Fills
// prepared[] (the concatenated expected solution) and the two cipher strings. Returns len1+len2.
static int plant(const char *kw, int p1, int p2, int len1, int len2,
                 int prepared[], char cs1[], char cs2[]) {
    int nr = 0;
    for (int i = 0; PLAINTEXT[i] && nr < len1 + len2; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[nr++] = idx;
    }
    int n1 = len1, n2 = nr - len1;
    int cube[TRIFID_CELLS];
    cube_from_keyword(kw, cube);
    int c1[MAX_CIPHER_LENGTH], c2[MAX_CIPHER_LENGTH];
    twin_trifid_encrypt(prepared, n1, prepared + n1, n2, cube, TRIFID_SIDE, p1, p2, c1, c2);
    for (int i = 0; i < n1; i++) cs1[i] = index_to_char(c1[i]);  cs1[n1] = '\0';
    for (int i = 0; i < n2; i++) cs2[i] = index_to_char(c2[i]);  cs2[n2] = '\0';
    return n1 + n2;
}

// Solve and return the recovered fraction over the concatenated (msg1 ++ msg2) plaintext.
// p1>0 pins message 1's period, p2>0 pins message 2's. `method` overrides the scheme. Always
// applies the tuned registry schedule. Scores with the currently-loaded n-gram table.
static double solve_and_frac(const char *cs1, const char *cs2, const int prepared[], int plen,
        int p1, int p2, int ngram_size, int method, uint32_t seed, int *period1_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TWIN_TRIFID;
    cfg.ngram_size = ngram_size;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    cfg.twincipher_str = (char *) cs2;
    if (p1 > 0) { cfg.period_present = true; cfg.period = p1; }
    if (p2 > 0) { cfg.period2_present = true; cfg.period2 = p2; }

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
    init_config(&cfg); cfg.cipher_type = TWIN_TRIFID; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "twin-trifid registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 300000,
        "twin-trifid anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "twin-trifid anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = TWIN_TRIFID; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "twin-trifid registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 20 && cfg.n_hill_climbs == 200000,
        "twin-trifid shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. period estimator (in isolation, long text) ----------------------------

static void test_period_estimator(void) {
    int periods[] = {5, 7, 8, 11};
    int lengths[] = {760, 560};
    int trials = 0, hits = 0;
    for (int pi = 0; pi < 4; pi++) {
        for (int li = 0; li < 2; li++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
            int p = periods[pi], L = lengths[li];
            plant("KRYPTOS", p, p, L, 8, prepared, cs1, cs2);
            int cidx[MAX_CIPHER_LENGTH];
            int n = 0; for (int i = 0; cs1[i]; i++) cidx[n++] = letter_to_index((unsigned char) cs1[i]);
            int out[8];
            int m = trifid_estimate_periods(cidx, n, 2, 16, 6, out, false);
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

// --- 3. capability floor at 210 letters each (periods pinned) ------------------
//
// 210 letters EACH (420 combined), periods pinned. The shared-cube joint score recovers here --
// well below a lone Trifid's ~500+ floor (the twin advantage), though above Twin Bifid's ~110
// (the 27-cell cube is harder). Quintgrams.

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
    int p1 = 7, p2 = 8, L = 210;
    int n = plant("KRYPTOS", p1, p2, L, L, prepared, cs1, cs2);
    int got_p1 = -1; double secs = 0.0;
    double frac = solve_and_frac(cs1, cs2, prepared, n, p1, p2, 5, METHOD_DEFAULT, 1u, &got_p1, &secs);
    printf("[capability @ %d+%d chars, periods PINNED (%d,%d), quintgrams] frac=%.3f %.1fs\n",
        L, L, p1, p2, frac, secs);
    CHECK(frac >= 0.95, "twin-trifid %d+%d (periods pinned) recovered only %.3f", L, L, frac);
}

// --- 4. length cliff (periods pinned) -----------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {150, 180, 210};
    printf("recovery vs per-message length (keyword KRYPTOS, periods 7 & 8 pinned, registry anneal, quints, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
        int L = lengths[li];
        int n = plant("KRYPTOS", 7, 8, L, L, prepared, cs1, cs2);
        double secs = 0.0;
        double frac = solve_and_frac(cs1, cs2, prepared, n, 7, 8, 5, METHOD_DEFAULT, 1u, NULL, &secs);
        printf("    %d+%d (%d)  frac=%.3f  %.1fs\n", L, L, n, frac, secs);
        if (li == 2) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.95, "twin-trifid 210+210 (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 5. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cs1[MAX_CIPHER_LENGTH], cs2[MAX_CIPHER_LENGTH];
    int p1 = 7, p2 = 8, L = 210;
    int n = plant("KRYPTOS", p1, p2, L, L, prepared, cs1, cs2);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d+%d chars, keyword KRYPTOS, periods (%d,%d) pinned, quints]\n", L, L, p1, p2);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs1, cs2, prepared, n, p1, p2, 5, M[m].method, 0x5C8E0u + m, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet_trifid();             // 27-symbol alphabet (A..Z + '+')
    CHECK(g_alpha == TRIFID_CELLS, "alphabet size %d, expected %d", g_alpha, TRIFID_CELLS);

    g_ngram_logprob = true;             // Twin Trifid needs the log-probability fitness
    shared.ngram_data = load_ngrams(QUINT_FILE, 5, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", QUINT_FILE);
        return 1;
    }

    test_registry();
    test_period_estimator();
    test_capability();
    test_length_cliff();
    test_per_scheme();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
