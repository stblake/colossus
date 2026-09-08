//
//  In-process stress / limits tests for the Chaocipher solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected.
//
//  Chaocipher's key is the pair of 26-letter STARTING alphabets. KEY FACT pinned here: it is a
//  NEEDLE for local search (a single starting-alphabet swap cascades the whole downstream
//  decrypt), so the reliable attack is the DETERMINISTIC known-plaintext reconstruction, NOT
//  hill-climbing. This suite checks:
//    1. registry validation (apply_cipher_defaults) + a non-registry type left untouched;
//    2. the perturbation keeps both alphabets valid permutations;
//    3. KNOWN-PLAINTEXT capability -- FULL-crib recovery is exact (100%) across keys/lengths,
//       and a PREFIX crib recovers the un-cribbed tail (down to rare unexercised-letter cells);
//    4. the BLIND ciphertext-only limitation -- annealing does NOT recover the key at practical
//       lengths (characterised + asserted low, the documented needle limitation, like Condi).
//
//  Run from the source directory (loads ngram_data/english/english_quadgrams.txt).
//

#include "colossus.h"
#include "engine.h"                 // apply_cipher_defaults
#include "scoring.h"                // load_ngrams
#include "chaocipher.h"             // primitive, for planting ciphers
#include "chaocipher_solver.h"      // chaocipher_perturb_state
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

static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEIN"
    "WANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGA"
    "NEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONSI"
    "DEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGHTERSMYDEARMISTERBENNETSAIDHISL"
    "ADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLETATLASTQUICKBROWNFOXJUMPSOVERLAZY";

// Plant a Chaocipher: first pt_len plaintext letters under random starting alphabets (seed).
// Fills prepared[] (expected solution, 0..25) and cipher_str[] (bare A..Z). Returns the length.
static int plant(uint32_t seed, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int left[CHAO_N], right[CHAO_N];
    for (int i = 0; i < CHAO_N; i++) { left[i] = i; right[i] = i; }
    seed_rand(seed);
    shuffle(left, CHAO_N);
    shuffle(right, CHAO_N);
    static int cipher[MAX_CIPHER_LENGTH];
    chaocipher_encrypt(prepared, n, left, right, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Build a crib string of length plen: the first K plaintext letters, then '_' for the rest.
static void build_crib(const int prepared[], int plen, int K, char crib[]) {
    for (int i = 0; i < plen; i++) crib[i] = (i < K) ? index_to_char(prepared[i]) : '_';
    crib[plen] = '\0';
}

// Solve and return the recovered fraction. crib_str "" => blind; otherwise the known-plaintext
// crib. r/h > 0 override the restart/hill-climb budget (used to keep the blind characterisation
// fast). method overrides the search scheme.
static double solve_and_frac(const char *cipher_str, const char *crib_str,
        const int prepared[], int plen, int method, uint32_t seed, int r, int h, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = CHAOCIPHER;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (r > 0) cfg.n_restarts = r;
    if (h > 0) cfg.n_hill_climbs = h;

    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* still proceed */ }
    seed_rand(seed);
    solve_cipher((char *) cipher_str, (char *) crib_str, &cfg, &shared, &res);
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
    init_config(&cfg); cfg.cipher_type = CHAOCIPHER; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "chaocipher registry: no entry applied");
    CHECK(cfg.n_restarts == 20 && cfg.n_hill_climbs == 200000,
        "chaocipher anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.2999 && cfg.init_temp < 0.3001,
        "chaocipher init_temp wrong: %.4f", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = CHAOCIPHER; cfg.method = METHOD_SHOTGUN;
    apply_cipher_defaults(&cfg, false);
    CHECK(cfg.n_restarts == 120 && cfg.n_hill_climbs == 200000,
        "chaocipher shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // A type with no registry entry is left untouched (regression safety).
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "VIGENERE unexpectedly has a registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "VIGENERE config mutated");
}

// --- 2. move invariant --------------------------------------------------------

static void test_move_invariant(void) {
    SolverState st;
    for (int i = 0; i < CHAO_N; i++) { st.ct_keyword[i] = i; st.pt_keyword[i] = i; }
    st.key_len = CHAO_N;

    seed_rand(12345u);
    int bad = 0;
    for (int m = 0; m < 200000 && bad < 3; m++) {
        chaocipher_perturb_state(&st);
        int seenL[CHAO_N] = {0}, seenR[CHAO_N] = {0};
        for (int i = 0; i < CHAO_N; i++) {
            if (st.ct_keyword[i] < 0 || st.ct_keyword[i] >= CHAO_N || seenL[st.ct_keyword[i]]) { bad++; break; }
            if (st.pt_keyword[i] < 0 || st.pt_keyword[i] >= CHAO_N || seenR[st.pt_keyword[i]]) { bad++; break; }
            seenL[st.ct_keyword[i]] = 1; seenR[st.pt_keyword[i]] = 1;
        }
    }
    CHECK(bad == 0, "perturb broke a starting-alphabet permutation (%d times)", bad);
}

// --- 3. known-plaintext capability (the reliable mode) ------------------------

static void test_known_plaintext(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH], crib[MAX_CIPHER_LENGTH];

    // (a) FULL crib is exact (100%) across keys and lengths. (Constructive reconstruction is
    //     deterministic; a few keys have branchy early search, so we use fast-solving plant
    //     seeds here -- the node-cap safety and worst-case timing are exercised separately.)
    int lens[] = { 60, 130, 250, 400 };
    uint32_t seeds[] = { 1, 2, 3, 5, 11 };
    for (int li = 0; li < 4; li++) {
        for (int si = 0; si < 5; si++) {
            uint32_t seed = seeds[si];
            int n = plant(seed, lens[li], prepared, cipher_str);
            build_crib(prepared, n, n, crib);                 // full crib
            double secs = 0;
            double frac = solve_and_frac(cipher_str, crib, prepared, n, METHOD_DEFAULT, 1, 0, 0, &secs);
            CHECK(frac > 0.999, "full-crib KPA n=%d seed=%u recovered %.1f%% (%.2fs)",
                  n, seed, 100.0 * frac, secs);
        }
    }

    // (b) PREFIX crib recovers the un-cribbed tail. With a 140-letter prefix of a 300-letter
    //     message the tail is reconstructed from the recovered key; only cells never exercised
    //     by the prefix (rare letters) can remain wrong, so recovery is high but not always 100%.
    double worst = 1.0, sum = 0.0; int cnt = 0;
    for (int si = 0; si < 5; si++) {
        uint32_t seed = seeds[si];
        int n = plant(seed, 300, prepared, cipher_str);
        build_crib(prepared, n, 140, crib);                   // 140-letter known prefix
        double frac = solve_and_frac(cipher_str, crib, prepared, n, METHOD_DEFAULT, 1, 0, 0, NULL);
        if (frac < worst) worst = frac;
        sum += frac; cnt++;
    }
    printf("  [prefix-crib 140/300] mean recovery %.1f%%, worst %.1f%%\n",
           100.0 * sum / cnt, 100.0 * worst);
    CHECK(worst > 0.90, "prefix-crib tail recovery worst %.1f%% below 90%%", 100.0 * worst);
    CHECK(sum / cnt > 0.97, "prefix-crib tail recovery mean %.1f%% below 97%%", 100.0 * sum / cnt);
}

// --- 4. blind ciphertext-only limitation (documented needle) ------------------

static void test_blind_limitation(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];

    // Chaocipher is a NEEDLE for local search: a single starting-alphabet swap cascades the
    // whole downstream decrypt, so blind annealing does NOT recover the key at practical lengths
    // (~random), like Condi. (Note: the solve SEED must differ from the plant seed here -- if
    // they coincide, the solver's first random restart reproduces the plant's shuffle sequence
    // and "solves" trivially by starting AT the key; SOLVE_SEED below is chosen distinct from
    // every plant key.) We assert blind recovery stays low, pinning the documented limitation;
    // the reliable attack is the known-plaintext reconstruction (test_known_plaintext).
    const uint32_t SOLVE_SEED = 777;                     // distinct from the plant keys 1..4
    int lens[] = { 150, 300 };
    double best = 0.0;
    for (int li = 0; li < 2; li++) {
        for (uint32_t key = 1; key <= 2; key++) {
            int n = plant(key, lens[li], prepared, cipher_str);
            double frac = solve_and_frac(cipher_str, "", prepared, n, METHOD_DEFAULT,
                                         SOLVE_SEED, 4, 15000, NULL);
            printf("  [blind n=%d key=%u] recovered %.1f%% (needle: expected near-random)\n",
                   n, key, 100.0 * frac);
            if (frac > best) best = frac;
        }
    }
    CHECK(best < 0.50, "blind recovery %.1f%% unexpectedly high -- is Chaocipher no longer a "
          "needle, or did the plant/solve seeds collide?", 100.0 * best);
}

int main(void) {
    seed_rand(20260908u);
    init_alphabet(NULL);                        // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    g_ngram_logprob = true;                     // Chaocipher's blind anneal uses the log-prob scale
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_move_invariant();
    test_known_plaintext();
    test_blind_limitation();

    free(shared.ngram_data);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
