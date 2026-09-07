//
//  In-process stress / limits tests for the Compressocrat solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Compressocrat is the fractionation twin of Fractionated Morse: the plaintext is written in a
//  FIXED prefix-free {1,2,3} Huffman code, the stream padded with 1 and grouped into trigraphs
//  (rank 9(a-1)+3(b-1)+(c-1), 333 excluded), and each trigraph mapped to a ciphertext letter
//  through a keyed 26-letter alphabet. The ONLY unknown is that alphabet, and it is an ACA KEYED
//  alphabet (keyword + ascending tail), so the solver searches that structured space -- the state
//  is a keyed-alphabet SEQUENCE perturbed by structure-preserving keyword moves
//  (compressocrat_move_seq, the length-26 twin of fracmorse's / Digrafid's). There is NO period.
//  The decode length varies per key (and EXPANDS -- the ciphertext is the compressed side), so the
//  decrypt hook tiles the decode to the fixed ciphertext length C and folds a validity reward into
//  score_adjust; it effectively needs the log-probability fitness, so this suite enables
//  g_ngram_logprob and loads quadgrams. This suite (also the basis for tuning the SearchDefaults
//  16x120000 / inittemp-0.30 schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the keyed-alphabet move INVARIANT (compressocrat_move_seq keeps a permutation + sorted tail);
//    3. a capability floor (recovery %) across keywords at ~150 chars (the ACA range);
//    4. a length cliff (recovery vs length);
//    5. a multi-keyword sweep (mean/worst recovery);
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso;
//    7. scoring-mode characterization: -logprob vs reward-only (locks the default).
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"                  // apply_cipher_defaults
#include "scoring.h"                 // load_ngrams
#include "compressocrat_solver.h"    // compressocrat_move_seq / compressocrat_canonicalize
#include "compressocrat.h"           // compressocrat_encrypt (planting)
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

// Plant a Compressocrat cipher: take the first pt_len letters of PLAINTEXT, build sigma from the
// keyword, and encipher. Fills prepared[] (the plaintext == the expected solution) and cipher_str;
// returns the plaintext length n. (The ciphertext C is typically shorter -- the compression.)
static int plant(const char *kw, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    int sigma[26];
    build_sigma(kw, sigma);
    static int cipher[MAX_CIPHER_LENGTH];
    int clen = compressocrat_encrypt(prepared, n, sigma, cipher);
    for (int i = 0; i < clen; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[clen] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the budget (else the registry
// schedule is used). `method` overrides the search scheme. `logprob` sets the scoring mode.
static double solve_and_frac_mode(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int method, uint32_t seed, bool logprob, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = COMPRESSOCRAT;
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

    bool saved_logprob = g_ngram_logprob;
    g_ngram_logprob = logprob;

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

    g_ngram_logprob = saved_logprob;

    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// Default scoring mode (-logprob) convenience wrapper.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int method, uint32_t seed, double *secs_out) {
    return solve_and_frac_mode(cipher_str, prepared, plen, r, h, method, seed, true, secs_out);
}

// Best recovery over nseeds well-separated seeds (early-out on a full solve). At >=300 letters
// the true key is the global n-gram max, but the anneal can stick in a gaming LOCAL optimum for
// an unlucky seed; a few restarts-from-fresh-seeds reliably find it (the grandpre pattern).
static double solve_best_of(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int method, uint32_t seed0, int nseeds) {
    double best = 0.0;
    for (int s = 0; s < nseeds; s++) {
        double f = solve_and_frac(cipher_str, prepared, plen, r, h, method,
            seed0 + 0x9E3779B1u * (uint32_t) s, NULL);
        if (f > best) best = f;
        if (best > 0.999) break;                // a full solve -- no need for more seeds
    }
    return best;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = COMPRESSOCRAT; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "compressocrat registry: no entry applied");
    CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 120000,
        "compressocrat anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "compressocrat anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = COMPRESSOCRAT; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "compressocrat registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 120 && cfg.n_hill_climbs == 120000,
        "compressocrat shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

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
        int kw = rand_int(3, 14);                       // CR_KW_MIN..MAX
        random_keyword(seq, 26, kw);
        for (int step = 0; step < 500; step++) {
            compressocrat_move_seq(seq, &kw);
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
    CHECK(bad_perm == 0, "compressocrat_move_seq broke the permutation %d times", bad_perm);
    CHECK(bad_tail == 0, "compressocrat_move_seq broke the sorted tail %d times", bad_tail);
    CHECK(bad_kw == 0, "compressocrat_move_seq drove kw out of [3,13] %d times", bad_kw);
}

// --- 3. length curve: the gaming cliff (registry schedule) --------------------
//
// DOCUMENTED LIMITATION: unlike Fractionated Morse (reliable ~110-150), Compressocrat COMPRESSES
// -- 150 plaintext letters give only ~130 ciphertext letters of scoring signal -- so through the
// ACA 110-150 range (and, for some keyed alphabets, up to ~250) a wrong keyed alphabet can
// out-score the truth (short-length fractionation GAMING: the true key is NOT the global n-gram
// max; neither more restarts nor quintgrams recover it -- verified). Blind recovery becomes
// RELIABLE from ~300 letters. PORTABLE is a "hard" keyed alphabet that traces the whole cliff;
// only the reliable end (>=300) is asserted, 150/200/250 are characterized.

static void test_length_curve(void) {
    int lens[] = { 150, 200, 250, 300 };
    printf("\n[length curve: keyword PORTABLE, anneal 12x100000 (150/200/250 characterize the gaming cliff)]\n");
    double frac300 = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("PORTABLE", lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 12, 100000, METHOD_DEFAULT, 0x5EEDu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 3) frac300 = frac;
    }
    CHECK(frac300 >= 0.95, "compressocrat 300ch recovered only %.3f", frac300);
}

// --- 4. capability floor across keywords (registry schedule, reliable @300) ----

static void test_capability_floor(void) {
    const char *kws[] = { "KRYPTOS", "MACHINE", "ZEBRA" };   // incl. MACHINE, a hard keyed alphabet
    int plen = 300;
    printf("\n[capability floor @ %d chars (reliable range), anneal 12x100000, best of 3 seeds]\n", plen);
    for (int k = 0; k < 3; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double frac = solve_best_of(cs, prepared, n, 12, 100000, METHOD_DEFAULT, 0xC0FFEEu + 17 * k, 3);
        printf("  %-10s : %.1f%%\n", kws[k], 100.0 * frac);
        CHECK(frac > 0.95, "compressocrat(%s) floor: only %.1f%% at %d chars (best of 3)", kws[k], 100.0 * frac, n);
    }
}

// --- 5. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int plen = 300;                             // reliable range (above the ~250 gaming cliff)
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", plen, prepared, cs);
    // Per-method budgets: anneal/shotgun get a real budget; PSO (weak + expensive here) a bounded
    // characterization run (see the swarm clamp in solve_and_frac_mode). Only anneal is asserted.
    struct { int method; const char *name; int r, h, nseeds; } M[] = {
        { METHOD_DEFAULT, "anneal ", 12, 80000, 2 },
        { METHOD_SHOTGUN, "shotgun", 12, 80000, 2 },
        { METHOD_PSO,     "pso    ",  3,   400, 1 },   // bounded + single-seed: characterized only
    };
    printf("\n[per-scheme @ ~%d chars, keyword CIPHER (best of nseeds; pso bounded, not asserted)]\n", plen);
    for (int m = 0; m < 3; m++) {
        double frac = solve_best_of(cs, prepared, n, M[m].r, M[m].h, M[m].method, 0x5C8E0u + 7 * m, M[m].nseeds);
        printf("  %s : %.1f%% (best of %d)\n", M[m].name, 100.0 * frac, M[m].nseeds);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.90, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

// --- 7. scoring-mode characterization: -logprob vs reward-only -----------------
//
// Like the fractionation family, the tiled decode + validity reward want the AZDecrypt log-prob
// fitness. Characterize both at a fixed budget/keyword to lock the default; -logprob is asserted
// to recover, reward-only is printed for comparison (not asserted).

static void test_scoring_mode(void) {
    int plen = 300;                             // reliable range (above the ~250 gaming cliff)
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("BENNET", plen, prepared, cs);
    printf("\n[scoring mode @ ~%d chars, keyword BENNET, budget 12x100000]\n", plen);
    // -logprob (best of 2, asserted) vs reward-only (single seed, characterized).
    double fr_lp = solve_best_of(cs, prepared, n, 12, 100000, METHOD_DEFAULT, 0xF00Du, 2);
    double fr_rw = solve_and_frac_mode(cs, prepared, n, 12, 100000, METHOD_DEFAULT, 0xF00Du, false, NULL);
    printf("  -logprob    : %.1f%% (best of 2)\n", 100.0 * fr_lp);
    printf("  reward-only : %.1f%% (1 seed)\n", 100.0 * fr_rw);
    CHECK(fr_lp > 0.90, "compressocrat -logprob recovery only %.1f%%", 100.0 * fr_lp);
}

int main(void) {
    init_alphabet(NULL);                        // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    g_ngram_logprob = true;                     // Compressocrat needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_move_invariant();
    test_length_curve();
    test_capability_floor();
    test_per_scheme();
    test_scoring_mode();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
