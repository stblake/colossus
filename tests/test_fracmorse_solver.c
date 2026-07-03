//
//  In-process stress / limits tests for the Fractionated Morse solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Fractionated Morse writes the plaintext in Morse with single 'x' separators, groups the
//  {DOT,DASH,X} stream into trigraphs (rank 9a+3b+c, xxx excluded), and maps each trigraph to a
//  ciphertext letter through a keyed 26-letter alphabet. The ONLY unknown is that alphabet, and
//  it is an ACA KEYED alphabet (keyword + ascending tail), so the solver searches that structured
//  space -- the state is a keyed-alphabet SEQUENCE perturbed by structure-preserving keyword moves
//  (fracmorse_move_seq), the length-26 twin of Digrafid's. There is NO period. The decode length
//  varies per key, so the decrypt hook tiles the decode to the fixed ciphertext length C and folds
//  a Morse-validity reward into score_adjust; it effectively needs the log-probability fitness, so
//  this suite enables g_ngram_logprob and loads quadgrams. The keyed-alphabet prior makes recovery
//  reliable well into the short ACA range (~110-150 letters). This suite (also the basis for tuning
//  the SearchDefaults 16x120000 / inittemp-0.30 schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the keyed-alphabet move INVARIANT (fracmorse_move_seq keeps a permutation + sorted tail);
//    3. a capability floor (recovery %) across keywords at ~150 chars (the ACA range);
//    4. a length cliff (recovery vs length) -- strong from ~120 up;
//    5. a multi-keyword sweep (mean/worst recovery);
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"             // apply_cipher_defaults
#include "scoring.h"            // load_ngrams
#include "fracmorse_solver.h"   // fracmorse_move_seq / fracmorse_canonicalize
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
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
    "ISTAKENBYAYOUNGMANOFLARGEFORTUNEFROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEE";

static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Keyed alphabet sigma (rank -> letter): keyword letters (dedup) then the rest ascending.
static void build_sigma(const char *kw, int sigma[]) {
    char used[26];
    for (int i = 0; i < 26; i++) used[i] = 0;
    int m = 0;
    for (int i = 0; kw[i]; i++) {
        int l = letter_to_index((unsigned char) kw[i]);
        if (l < 0 || used[l]) continue;
        used[l] = 1; sigma[m++] = l;
    }
    for (int l = 0; l < 26 && m < 26; l++) if (!used[l]) { used[l] = 1; sigma[m++] = l; }
}

// Plant a Fractionated Morse cipher: take the first pt_len letters of PLAINTEXT, build sigma from
// the keyword, and encipher. Fills prepared[] (the plaintext == the expected solution) and
// cipher_str; returns the plaintext length n.
static int plant(const char *kw, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    int sigma[26];
    build_sigma(kw, sigma);
    static int cipher[MAX_CIPHER_LENGTH];
    int clen = fracmorse_encrypt(prepared, n, sigma, cipher);
    for (int i = 0; i < clen; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[clen] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the budget (else the registry
// schedule is used). `method` overrides the search scheme.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int method, uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = FRAC_MORSE;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) { cfg.n_restarts = r; }
    if (h > 0) { cfg.n_hill_climbs = h; }
    // PSO cost is restarts x iterations x particles x refine; the global 30x50 particle/refine
    // defaults make a large iteration budget pathological, so bound the swarm for the per-scheme
    // characterization (PSO is not the recommended scheme -- only anneal's recovery is asserted).
    if (method == METHOD_PSO) { cfg.n_particles = 12; cfg.refine_steps = 4; }

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
    init_config(&cfg); cfg.cipher_type = FRAC_MORSE; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "fracmorse registry: no entry applied");
    CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 120000,
        "fracmorse anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "fracmorse anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = FRAC_MORSE; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "fracmorse registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 120 && cfg.n_hill_climbs == 120000,
        "fracmorse shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. keyed-alphabet move invariant -----------------------------------------

static void test_move_invariant(void) {
    int bad_perm = 0, bad_tail = 0, bad_kw = 0;
    seed_rand(0xA11CEu);
    for (int trial = 0; trial < 400; trial++) {
        int seq[26];
        int kw = rand_int(3, 14);                       // FM_KW_MIN..MAX
        random_keyword(seq, 26, kw);
        for (int step = 0; step < 500; step++) {
            fracmorse_move_seq(seq, &kw);
            int seen[26] = {0};
            for (int i = 0; i < 26; i++) {
                if (seq[i] < 0 || seq[i] >= 26 || seen[seq[i]]) { bad_perm++; break; }
                seen[seq[i]] = 1;
            }
            for (int i = kw; i + 1 < 26; i++)
                if (seq[i] >= seq[i + 1]) { bad_tail++; break; }
            if (kw < 3 || kw > 13) bad_kw++;
        }
    }
    printf("[move invariant] 400 chains x 500 moves: perm-violations=%d tail-violations=%d kw-violations=%d\n",
        bad_perm, bad_tail, bad_kw);
    CHECK(bad_perm == 0, "fracmorse_move_seq broke the permutation %d times", bad_perm);
    CHECK(bad_tail == 0, "fracmorse_move_seq broke the sorted tail %d times", bad_tail);
    CHECK(bad_kw == 0, "fracmorse_move_seq drove kw out of [3,13] %d times", bad_kw);
}

// --- 3. capability floor (registry schedule) ----------------------------------
//
// The ACA range is ~110-150 plaintext letters. Run several keywords at ~150 through the tuned
// registry schedule; the keyed-alphabet prior recovers these cleanly (a persistent ~2% shortfall
// on some short texts is a rare-letter ambiguity, so the floor asserts > 0.90).

static void test_capability_floor(void) {
    const char *kws[] = { "KRYPTOS", "PORTABLE", "MACHINE" };
    int plen = 150;
    printf("\n[capability floor @ ~%d chars (ACA range), registry schedule]\n", plen);
    for (int k = 0; k < 3; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 0, 0, METHOD_DEFAULT, 0xC0FFEEu + 17 * k, &secs);
        printf("  %-10s : %.1f%%  [%.1fs]\n", kws[k], 100.0 * frac, secs);
        CHECK(frac > 0.90, "fracmorse(%s) floor: only %.1f%% at %d chars", kws[k], 100.0 * frac, n);
    }
}

// --- 4. length cliff (registry schedule) --------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 90, 150, 300 };
    printf("\n[length cliff: keyword SHADOW, registry anneal]\n");
    double frac150 = 0.0, frac300 = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("SHADOW", lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 1) frac150 = frac;
        if (li == 2) frac300 = frac;
    }
    CHECK(frac300 >= 0.95, "fracmorse 300ch recovered only %.3f", frac300);
    CHECK(frac150 >= 0.90, "fracmorse 150ch recovered only %.3f", frac150);
}

// --- 5. multi-keyword sweep ---------------------------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "MONARCHY", "GADGETRY", "CIPHERTEXT" };
    int plen = 250, nk = 4;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars, budget 12x100000]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, 12, 100000, METHOD_DEFAULT, 0xABCDu + k, NULL);
        printf("  %-10s : %.1f%%\n", kws[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int plen = 200;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", plen, prepared, cs);
    // Per-method budgets: anneal/shotgun get a real budget; PSO (weak + expensive here) a bounded
    // characterization run (see the swarm clamp in solve_and_frac). Only anneal is asserted.
    struct { int method; const char *name; int r, h; } M[] = {
        { METHOD_DEFAULT, "anneal ", 12, 80000 },
        { METHOD_SHOTGUN, "shotgun", 12, 80000 },
        { METHOD_PSO,     "pso    ",  3,   400 },
    };
    printf("\n[per-scheme @ ~%d chars, keyword CIPHER (pso bounded, not asserted)]\n", plen);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, M[m].r, M[m].h, M[m].method, 0x5C8E0u + m, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.90, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet(NULL);                        // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    g_ngram_logprob = true;                     // Fractionated Morse needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_move_invariant();
    test_capability_floor();
    test_length_cliff();
    test_multi_keyword();
    test_per_scheme();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
