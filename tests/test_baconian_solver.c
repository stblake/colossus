//
//  In-process stress / limits tests for the Baconian solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Baconian biliteral-encodes the plaintext (five a/b symbols per letter, fixed 24-letter table)
//  and CONCEALS the a/b stream in cover text under a hidden CLASSIFIER (per-letter or per-word).
//  The solver searches that 26-letter a/b classifier: canonical rules (A-M/N-Z, vowel/consonant x
//  polarity) are single SWEEP cells, and one FREE-CLIMB config per mode anneals a general
//  labelling, guided by a biliteral-VALIDITY reward. It effectively needs -logprob, so this suite
//  enables g_ngram_logprob and loads quadgrams. It checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. EXACT recovery of a canonical (A-M/N-Z) cover text -- the sweep decodes the truth exactly;
//    3. a length curve (CHARACTERIZED at the ACA <=25-letter maximum, ASSERTED where reliable);
//    4. per-scheme calibration: the same cipher under -method anneal / shotgun / pso;
//    5. per-word grouping + auto mode-selection;
//    6. reward-only vs -logprob characterization.
//
//  Blind classifier recovery at the ACA <=25-letter maximum is a documented gaming limitation
//  (n-gram is weak there); the reliable case -- a canonical classifier caught by the sweep -- is
//  what is asserted. Synthetic cover text picks random letters within the correct group (an
//  optimistic bound vs English-crafted cover). Run from the source directory (n-gram table in cwd).
//

#include "colossus.h"
#include "engine.h"             // apply_cipher_defaults
#include "scoring.h"            // load_ngrams
#include "baconian.h"           // baconian_encrypt
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

static int letter_to_index(int c) { c = toupper(c); return (c < 'A' || c > 'Z') ? -1 : g_char_to_idx[c]; }
static int fold(int l) { return (l == 9) ? 8 : (l == 21) ? 20 : l; }   // J->I, V->U (the 24-letter decode)
static int is_vowel(int l) { return (l == 0 || l == 4 || l == 8 || l == 14 || l == 20); }
static int classify(int l, int rule_am, int pol) {
    int c = rule_am ? ((l >= 13) ? 1 : 0) : (is_vowel(l) ? 0 : 1);
    return pol ? (c ^ 1) : c;
}

// Plant a Baconian cover text: biliteral-encode the first pt_len letters of PLAINTEXT, then
// conceal the a/b bits under (rule_am, pol) in the given grouping mode. Fills expected[] (=
// fold(plaintext), what the solver recovers) and cover_str; returns the plaintext length n.
static int plant(int pt_len, int word_mode, int rule_am, int pol, uint32_t seed,
                 int expected[], char cover_str[]) {
    int pt[512], n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) { pt[n] = idx; expected[n] = fold(idx); n++; }
    }
    static int bits[8 * 512];
    int nbits = baconian_encrypt(pt, n, bits);

    int grp[2][26], ng[2] = {0, 0};
    for (int l = 0; l < 26; l++) { int g = classify(l, rule_am, pol); grp[g][ng[g]++] = l; }

    seed_rand(seed);
    int p = 0;
    for (int i = 0; i < nbits; i++) {
        int b = bits[i];
        if (word_mode) {
            int wlen = 2 + (int) rand_bounded(5);
            cover_str[p++] = index_to_char(grp[b][rand_bounded(ng[b])]);
            for (int k = 1; k < wlen; k++) cover_str[p++] = index_to_char(rand_bounded(26));
            if (i + 1 < nbits) cover_str[p++] = ' ';
        } else {
            cover_str[p++] = index_to_char(grp[b][rand_bounded(ng[b])]);
            if ((i + 1) % 5 == 0 && i + 1 < nbits) cover_str[p++] = ' ';
        }
    }
    cover_str[p] = '\0';
    return n;
}

// Solve and return the recovered fraction. r/h > 0 override the registry budget; `method` the
// scheme; `mode_cfg` the -baconmode. Compares res.decrypted to expected[] (the folded plaintext).
static double solve_and_frac(const char *cover_str, const int expected[], int plen,
        int mode_cfg, int r, int h, int method, uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = BACONIAN;
    cfg.bacon_mode = mode_cfg;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) cfg.n_restarts = r;
    if (h > 0) cfg.n_hill_climbs = h;
    if (method == METHOD_PSO) { cfg.n_particles = 12; cfg.refine_steps = 4; }

    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* still proceed */ }
    seed_rand(seed);
    solve_cipher((char *) cover_str, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);
    if (secs_out) *secs_out = ((double) clock() - t0) / CLOCKS_PER_SEC;

    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == expected[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = BACONIAN; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "baconian registry: no entry applied");
    CHECK(cfg.n_restarts == 24 && cfg.n_hill_climbs == 60000,
        "baconian anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "baconian anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = BACONIAN; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "baconian registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 120 && cfg.n_hill_climbs == 60000,
        "baconian shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs; double t0 = cfg.init_temp;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0 && cfg.init_temp == t0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. exact recovery of a canonical (A-M/N-Z) cover text ---------------------
//
// The canonical classifier is a SWEEP cell, so it decodes the truth EXACTLY -- at a comfortable
// length the true plaintext (all-valid, real English) is the global n-gram max. Per-letter mode,
// pinned, so the solve is fast; assert ~100%.

static void test_exact_canonical(void) {
    static int expected[512]; static char cover[20000];
    printf("\n[exact recovery, canonical A-M/N-Z, per-letter]\n");
    int lens[] = { 120, 250 };
    for (int li = 0; li < 2; li++) {
        int n = plant(lens[li], 0, 1, 0, 0xBAC0u + li, expected, cover);
        double secs;
        double frac = solve_and_frac(cover, expected, n, BAC_MODE_LETTER, 0, 0, METHOD_DEFAULT,
            0x5EEDu + li, &secs);
        printf("  len %3d : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        CHECK(frac > 0.99, "canonical exact recovery only %.1f%% at %d letters", 100.0 * frac, n);
    }
}

// --- 3. length curve (characterize the ACA range, assert where reliable) -------

static void test_length_curve(void) {
    int lens[] = { 15, 25, 50, 100 };
    static int expected[512]; static char cover[20000];
    printf("\n[length curve: canonical A-M/N-Z per-letter, registry anneal]\n");
    double frac100 = 0.0;
    for (int li = 0; li < 4; li++) {
        int n = plant(lens[li], 0, 1, 0, 0xC0DEu + li, expected, cover);
        double secs;
        double frac = solve_and_frac(cover, expected, n, BAC_MODE_LETTER, 0, 0, METHOD_DEFAULT,
            0xA5A5u + li, &secs);
        printf("  len %3d : %.1f%%  [%.1fs]%s\n", n, 100.0 * frac, secs,
            n <= 25 ? "   (ACA range, characterized)" : "");
        if (li == 3) frac100 = frac;
    }
    CHECK(frac100 > 0.95, "baconian 100-letter recovery only %.1f%%", 100.0 * frac100);
}

// --- 4. per-scheme calibration (anneal / shotgun / pso) -----------------------

static void test_per_scheme(void) {
    static int expected[512]; static char cover[20000];
    int n = plant(120, 0, 1, 0, 0xF00Du, expected, cover);
    struct { int method; const char *name; int r, h; } M[] = {
        { METHOD_DEFAULT, "anneal ", 0, 0 },
        { METHOD_SHOTGUN, "shotgun", 0, 0 },
        { METHOD_PSO,     "pso    ", 3, 400 },
    };
    printf("\n[per-scheme @ 120 letters, canonical per-letter (pso bounded, not asserted)]\n");
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_and_frac(cover, expected, n, BAC_MODE_LETTER, M[m].r, M[m].h,
            M[m].method, 0x5C8E0u + m, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

// --- 5. per-word grouping + auto mode-selection -------------------------------

static void test_word_mode_and_auto(void) {
    static int expected[512]; static char cover[40000];
    printf("\n[per-word grouping + auto mode-selection, canonical A-M/N-Z]\n");

    // Pinned per-word (at a reliable length -- like per-letter, ~60 letters is still marginal).
    int n = plant(100, 1, 1, 0, 0x7717u, expected, cover);
    double secs;
    double frac = solve_and_frac(cover, expected, n, BAC_MODE_WORD, 0, 0, METHOD_DEFAULT, 0x1234u, &secs);
    printf("  per-word (pinned)  len %d : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
    CHECK(frac > 0.95, "per-word canonical recovery only %.1f%%", 100.0 * frac);

    // Auto mode: the same cover text solved without pinning -- the solver must pick per-word.
    frac = solve_and_frac(cover, expected, n, BAC_MODE_AUTO, 0, 0, METHOD_DEFAULT, 0x1234u, &secs);
    printf("  per-word (auto)    len %d : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
    CHECK(frac > 0.95, "auto-mode did not recover the per-word cover text (%.1f%%)", 100.0 * frac);
}

// --- 6. reward-only vs -logprob characterization ------------------------------

static void test_logprob_vs_reward(void) {
    static int expected[512]; static char cover[20000];
    int n = plant(120, 0, 1, 0, 0x9001u, expected, cover);
    printf("\n[reward-only vs -logprob @ %d letters, canonical per-letter]\n", n);
    double secs;
    double f_lp = solve_and_frac(cover, expected, n, BAC_MODE_LETTER, 0, 0, METHOD_DEFAULT, 0x2222u, &secs);
    printf("  -logprob    : %.1f%%\n", 100.0 * f_lp);
    g_ngram_logprob = false;
    double f_rw = solve_and_frac(cover, expected, n, BAC_MODE_LETTER, 0, 0, METHOD_DEFAULT, 0x2222u, &secs);
    g_ngram_logprob = true;
    printf("  reward-only : %.1f%%   (characterization only)\n", 100.0 * f_rw);
    CHECK(f_lp > 0.95, "logprob recovery only %.1f%%", 100.0 * f_lp);
}

int main(void) {
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_exact_canonical();
    test_length_curve();
    test_per_scheme();
    test_word_mode_and_auto();
    test_logprob_vs_reward();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
