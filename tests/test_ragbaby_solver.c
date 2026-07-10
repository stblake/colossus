//
//  In-process stress / limits tests for the Ragbaby solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Ragbaby shifts each plaintext LETTER forward by its word-position number (mod 24) in a keyed
//  24-letter alphabet (I/J and W/X paired). The per-letter numbering is FIXED by the (spaced)
//  ciphertext, so the ONLY unknown is the keyed alphabet, and it is an ACA KEYED alphabet (keyword
//  + ascending tail) -- the solver searches that structured space (ragbaby_move_seq, the 24-cell
//  twin of Digrafid/Fractionated-Morse). There is NO period. Because the per-letter shift is KNOWN
//  the keyed alphabet is heavily constrained, so Ragbaby rides the REWARD-ONLY quadgram table (no
//  -logprob, like the Vigenere/Porta family; the true alphabet scores highest) and recovers
//  reliably across the ACA ~100-150-letter range. At the very short (~80) end recovery is
//  seed-sensitive (a keyed-alphabet cliff, like every square/alphabet solver). This suite (also the
//  basis for tuning the SearchDefaults 16x120000 / inittemp-0.30 schedule) checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the keyed-alphabet move INVARIANT (ragbaby_move_seq keeps a permutation + sorted tail);
//    3. a capability floor (recovery %) across keywords in the ACA range;
//    4. a length cliff (recovery vs length) -- strong from ~110 up;
//    5. a multi-keyword sweep (mean/worst recovery);
//    6. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"           // apply_cipher_defaults
#include "scoring.h"          // load_ngrams
#include "ragbaby_solver.h"   // ragbaby_move_seq / ragbaby_canonicalize
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4
#define RAG_ALPHA  24

static SharedData shared;

// Word-divided natural English (letters + single spaces). Word boundaries drive the numbering, so
// this MUST keep its spaces (unlike the other solver tests' solid blocks).
static const char *PLAINTEXT =
    "it is a truth universally acknowledged that a single man in possession of a good fortune must "
    "be in want of a wife however little known the feelings or views of such a man may be on his "
    "first entering a neighbourhood this truth is so well fixed in the minds of the surrounding "
    "families that he is considered the rightful property of some one or other of their daughters "
    "my dear mister bennet said his lady to him one day have you heard that netherfield park is let "
    "at last mister bennet replied that he had not but it is returned she for missus long has just "
    "been here and she told me all about it mister bennet made no answer do you not want to know "
    "who has taken it cried his wife impatiently you want to tell me and i have no objection to hear";

// Build the keyed alphabet ka[] from a keyword (folded through the live 24-letter alphabet).
static void build_ka(const char *kw, int ka[]) {
    ragbaby_build_keyed_alphabet(ka, kw, RAG_ALPHA);
}

// Plant a Ragbaby cipher: take whole words of PLAINTEXT up to ~pt_len letters, encipher under the
// keyword's keyed alphabet, and build the SPACED cipher string the solver re-parses. Fills
// prepared[] (the folded plaintext letters == the expected solution) and cipher_str; returns the
// letter count n.
static int plant(const char *kw, int pt_len, int prepared[], char cipher_str[]) {
    // Copy whole words of PLAINTEXT until we have >= pt_len letters.
    char pt[MAX_CIPHER_LENGTH]; int m = 0, letters_so_far = 0; bool in_word = false;
    for (int i = 0; PLAINTEXT[i] && letters_so_far < pt_len; i++) {
        char c = PLAINTEXT[i];
        if (c == ' ') { in_word = false; if (m > 0) pt[m++] = ' '; continue; }
        pt[m++] = c;
        if (!in_word) in_word = true;
        letters_so_far++;
    }
    // trim a trailing space
    while (m > 0 && pt[m - 1] == ' ') m--;
    pt[m] = '\0';

    static int num[MAX_CIPHER_LENGTH];
    int n = 0;
    ragbaby_number_stream(pt, RAG_ALPHA, prepared, num, &n);   // prepared[] = folded plaintext letters

    int ka[RAG_ALPHA], ka_inv[RAG_ALPHA];
    build_ka(kw, ka);
    for (int p = 0; p < RAG_ALPHA; p++) ka_inv[ka[p]] = p;
    static int ct[MAX_CIPHER_LENGTH];
    ragbaby_encrypt(prepared, num, n, ka, ka_inv, RAG_ALPHA, ct);

    // Emit the spaced ciphertext: replace each plaintext letter with its cipher letter, keep single
    // spaces at word boundaries (exactly what tools/ragbaby_gen and the solver's parser expect).
    int idx = 0, cn = 0; bool prev_space = false, started = false;
    for (int i = 0; pt[i]; i++) {
        if (pt[i] == ' ') { prev_space = true; continue; }
        if (prev_space && started) cipher_str[cn++] = ' ';
        cipher_str[cn++] = index_to_char(ct[idx++]);
        prev_space = false; started = true;
    }
    cipher_str[cn] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the budget (else the registry schedule
// is used). `method` overrides the search scheme.
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int r, int h, int method, uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = RAGBABY;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) { cfg.n_restarts = r; }
    if (h > 0) { cfg.n_hill_climbs = h; }
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
    init_config(&cfg); cfg.cipher_type = RAGBABY; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "ragbaby registry: no entry applied");
    CHECK(cfg.n_restarts == 16 && cfg.n_hill_climbs == 120000,
        "ragbaby anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "ragbaby anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = RAGBABY; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "ragbaby registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 120 && cfg.n_hill_climbs == 120000,
        "ragbaby shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

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
        int seq[RAG_ALPHA];
        int kw = rand_int(3, 13);                        // RAG_KW_MIN..MAX
        random_keyword(seq, RAG_ALPHA, kw);
        for (int step = 0; step < 500; step++) {
            ragbaby_move_seq(seq, &kw);
            int seen[RAG_ALPHA] = {0};
            for (int i = 0; i < RAG_ALPHA; i++) {
                if (seq[i] < 0 || seq[i] >= RAG_ALPHA || seen[seq[i]]) { bad_perm++; break; }
                seen[seq[i]] = 1;
            }
            for (int i = kw; i + 1 < RAG_ALPHA; i++)
                if (seq[i] >= seq[i + 1]) { bad_tail++; break; }
            if (kw < 3 || kw > 12) bad_kw++;
        }
    }
    printf("[move invariant] 400 chains x 500 moves: perm-violations=%d tail-violations=%d kw-violations=%d\n",
        bad_perm, bad_tail, bad_kw);
    CHECK(bad_perm == 0, "ragbaby_move_seq broke the permutation %d times", bad_perm);
    CHECK(bad_tail == 0, "ragbaby_move_seq broke the sorted tail %d times", bad_tail);
    CHECK(bad_kw == 0, "ragbaby_move_seq drove kw out of [3,12] %d times", bad_kw);
}

// --- 3. capability floor (registry schedule) ----------------------------------

static void test_capability_floor(void) {
    const char *kws[] = { "GROSBEAK", "CRYPTOGRAM", "ZEBRA" };
    int plen = 150;
    printf("\n[capability floor @ ~%d chars (ACA range), registry schedule]\n", plen);
    double best = 0.0;
    for (int k = 0; k < 3; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 0, 0, METHOD_DEFAULT, 0xC0FFEEu + 17 * k, &secs);
        printf("  %-10s : %.1f%%  [%.1fs]\n", kws[k], 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    // At the ACA ~150 range at least one keyword recovers essentially fully (some keyed alphabets
    // hit short-text near-miss basins at a fixed seed -- characterized by the sweep below).
    CHECK(best > 0.98, "ragbaby capability floor: best only %.1f%% at ~%d chars", 100.0 * best, plen);
}

// --- 4. length cliff (registry schedule) --------------------------------------

static void test_length_cliff(void) {
    int lens[] = { 80, 150, 300 };
    printf("\n[length cliff: keyword ZEBRA, registry anneal]\n");
    double frac150 = 0.0, frac300 = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant("ZEBRA", lens[li], prepared, cs);
        double secs;
        double frac = solve_and_frac(cs, prepared, n, 0, 0, METHOD_DEFAULT, 0x5EEDu + li, &secs);
        printf("  len~%-4d (%d) : %.1f%%  [%.1fs]\n", lens[li], n, 100.0 * frac, secs);
        if (li == 1) frac150 = frac;
        if (li == 2) frac300 = frac;
    }
    CHECK(frac300 >= 0.98, "ragbaby 300ch recovered only %.3f", frac300);
    CHECK(frac150 >= 0.95, "ragbaby 150ch recovered only %.3f", frac150);
}

// --- 5. multi-keyword sweep ---------------------------------------------------

static void test_multi_keyword(void) {
    // Ragbaby is an 80-150-letter ACA cipher; sweep several keywords across that range.
    const char *kws[] = { "ZEBRA", "MONARCHY", "GADGETRY", "CIPHERTEXT" };
    int plen = 150, nk = 4;
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ ~%d chars (ACA range), registry schedule]\n", plen);
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
        int n = plant(kws[k], plen, prepared, cs);
        double frac = solve_and_frac(cs, prepared, n, 0, 0, METHOD_DEFAULT, 0xABCDu + k, NULL);
        printf("  %-10s : %.1f%%\n", kws[k], 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    // A keyed alphabet can hit a short-text near-miss basin at a fixed seed (a few swapped cells);
    // the mean stays high. The floor/cliff above are the strong capability guards.
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 6. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    int plen = 150;
    static int prepared[MAX_CIPHER_LENGTH]; static char cs[MAX_CIPHER_LENGTH];
    int n = plant("CIPHER", plen, prepared, cs);
    struct { int method; const char *name; int r, h; } M[] = {
        { METHOD_DEFAULT, "anneal ", 0, 0 },
        { METHOD_SHOTGUN, "shotgun", 60, 120000 },
        { METHOD_PSO,     "pso    ",  3,    400 },
    };
    printf("\n[per-scheme @ ~%d chars, keyword CIPHER (pso bounded, not asserted)]\n", plen);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cs, prepared, n, M[m].r, M[m].h, M[m].method, 0x5C8E0u + m, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet_ragbaby();                    // 24-letter Ragbaby alphabet (J->I, X->W)
    CHECK(g_alpha == RAG_ALPHA, "alphabet size %d, expected 24", g_alpha);

    // Reward-only quadgrams (NO -logprob): the known per-letter shift constrains the keyed alphabet
    // so the true alphabet scores highest under the reward-only table (like the Vigenere/Porta family).
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
