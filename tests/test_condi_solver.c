//
//  In-process stress / limits tests for the Condi solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Condi is a plaintext-feedback substitution over a keyed alphabet sigma. The KEY FINDING this
//  suite pins is a STRUCTURAL limitation (analogous to the CM-Bifid even-period degeneracy): the
//  feedback makes the true sigma an ISOLATED NEEDLE -- one sigma cell swap re-derives the whole
//  downstream plaintext (an alternating-sign cascade), so there is no basin and hence no gradient
//  for ANY local search. This suite:
//    1. validates the SearchDefaults registry entry (+ a non-registry type left untouched);
//    2. ASSERTS the needle: the true sigma out-scores its best single-swap neighbour by a clear
//       margin AND its neighbours sit near the random floor (a regression guard on the finding);
//    3. characterizes (printed, NOT asserted) that a blind solve does not recover, that pinning the
//       starter + a large budget still does not recover (no budget escapes the needle), and the
//       per-scheme (anneal / shotgun / pso) behaviour.
//
//  Condi effectively needs -logprob, so ngrams are loaded with the log-probability table. Run from
//  the source directory.
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams, ngram_score
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

static const char *PLAINTEXT =
    "WHENINTHECOURSEOFHUMANEVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICAL"
    "BANDSWHICHHAVECONNECTEDTHEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATE"
    "ANDEQUALSTATIONTOWHICHTHELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHE"
    "OPINIONSOFMANKINDREQUIRESTHATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATION"
    "WEHOLDTHESETRUTHSTOBESELFEVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIR";

static int take_plaintext(int prepared[], int pt_len) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    return n;
}

static void sigma_from_keyword(const char *kw, int shift, int sigma[]) {
    int keyed[ALPHABET_SIZE];
    char buf[64]; strncpy(buf, kw, 63); buf[63] = '\0';
    make_keyed_alphabet(buf, keyed);
    int s = ((shift % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    for (int k = 0; k < ALPHABET_SIZE; k++) sigma[k] = keyed[(k + s) % ALPHABET_SIZE];
}

// Plant a Condi cipher; returns n, fills prepared[], cipher_str, and the true sigma.
static int plant_condi(const char *kw, int shift, int starter, int pt_len,
                       int prepared[], char cipher_str[], int sigma_out[]) {
    int sigma[ALPHABET_SIZE], sinv[ALPHABET_SIZE];
    sigma_from_keyword(kw, shift, sigma);
    for (int k = 0; k < ALPHABET_SIZE; k++) sinv[sigma[k]] = k;
    int n = take_plaintext(prepared, pt_len);
    static int cipher[MAX_CIPHER_LENGTH];
    condi_encrypt(prepared, n, sigma, sinv, starter, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    if (sigma_out) for (int k = 0; k < ALPHABET_SIZE; k++) sigma_out[k] = sigma[k];
    return n;
}

// Solve harness. startkey>=0 pins the starter; nrestarts/nhill>0 override the schedule.
static double solve_ex(const char *cipher_str, const int prepared[], int plen,
        int startkey, int method, int nrestarts, int nhill, uint32_t seed,
        int *start_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = CONDI;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (startkey >= 0) { cfg.startkey_present = true; cfg.startkey = startkey; }
    if (nrestarts > 0) cfg.n_restarts = nrestarts;
    if (nhill > 0)     cfg.n_hill_climbs = nhill;

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
    if (secs_out)  *secs_out = ((double) clock() - t0) / CLOCKS_PER_SEC;
    if (start_out) *start_out = res.solved ? res.condi_start : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg; init_config(&cfg);
    cfg.cipher_type = CONDI;
    bool hit = apply_cipher_defaults(&cfg, false);
    CHECK(hit, "apply_cipher_defaults should have a CONDI entry");
    CHECK(cfg.n_restarts == 6 && cfg.n_hill_climbs == 60000,
          "CONDI schedule expected 6x60000, got %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // A non-registry type is left untouched (VIGENERE rides the global defaults).
    ColossusConfig a, b; init_config(&a); init_config(&b);
    a.cipher_type = VIGENERE;
    bool hit2 = apply_cipher_defaults(&a, false);
    CHECK(!hit2, "VIGENERE should have no registry entry");
    CHECK(a.n_restarts == b.n_restarts && a.n_hill_climbs == b.n_hill_climbs,
          "non-registry type must be unmodified");
    printf("[registry] CONDI -> %dx%d (anneal)\n", cfg.n_restarts, cfg.n_hill_climbs);
}

// --- 2. the needle: assert the structural no-gradient finding ------------------

static double score_sigma(const int cipher[], int n, const int sigma[], int starter) {
    int sinv[ALPHABET_SIZE], dec[MAX_CIPHER_LENGTH];
    for (int k = 0; k < ALPHABET_SIZE; k++) sinv[sigma[k]] = k;
    condi_decrypt(cipher, n, sigma, sinv, starter, dec);
    return ngram_score(dec, n, shared.ngram_data, NGRAM_SIZE);
}

static void test_needle(void) {
    int prepared[MAX_CIPHER_LENGTH], sigma[ALPHABET_SIZE];
    char cstr[MAX_CIPHER_LENGTH];
    int starter = 7;
    int n = plant_condi("CIPHER", 0, starter, 360, prepared, cstr, sigma);
    int cipher[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) cipher[i] = cstr[i] - 'A';

    double s_true = score_sigma(cipher, n, sigma, starter);

    // All single-swap neighbours of the true sigma.
    double best = -1e9, sum = 0; int cnt = 0;
    for (int a = 0; a < ALPHABET_SIZE; a++)
        for (int b = a + 1; b < ALPHABET_SIZE; b++) {
            int s[ALPHABET_SIZE];
            for (int k = 0; k < ALPHABET_SIZE; k++) s[k] = sigma[k];
            int t = s[a]; s[a] = s[b]; s[b] = t;
            double v = score_sigma(cipher, n, s, starter);
            if (v > best) best = v;
            sum += v; cnt++;
        }
    double mean_nb = sum / cnt;

    // Random-sigma baseline.
    double rsum = 0; int rn = 60;
    for (int r = 0; r < rn; r++) {
        int s[ALPHABET_SIZE]; for (int k = 0; k < ALPHABET_SIZE; k++) s[k] = k;
        seed_rand(r + 1); shuffle(s, ALPHABET_SIZE);
        rsum += score_sigma(cipher, n, s, starter);
    }
    double rand_floor = rsum / rn;

    printf("[needle] n=%d true=%.4f best-1swap=%.4f mean-1swap=%.4f random=%.4f\n",
           n, s_true, best, mean_nb, rand_floor);

    // The cliff: the true key beats EVERY single-swap neighbour by a clear margin (no smooth ascent).
    CHECK(s_true - best > 0.3,
          "expected a needle cliff (true - best_neighbour > 0.3), got %.4f", s_true - best);
    // Neighbours sit far below the true key, near the random floor (no basin).
    CHECK((s_true - mean_nb) > 0.6 * (s_true - rand_floor),
          "expected neighbours near the random floor (no basin)");
    CHECK(mean_nb < s_true - 1.0, "expected mean neighbour >1.0 below the true key");
}

// --- 3. characterization (printed, NOT asserted) -------------------------------

static void test_blind_characterization(void) {
    int prepared[MAX_CIPHER_LENGTH]; char cstr[MAX_CIPHER_LENGTH];
    int n = plant_condi("CIPHER", 0, 7, 180, prepared, cstr, NULL);
    double secs; int sout;
    double frac = solve_ex(cstr, prepared, n, -1, METHOD_DEFAULT, 0, 0, 1, &sout, &secs);
    printf("[blind] n=%d recovery=%.1f%% (starter %s) [%.1fs] -- expected LOW (needle)\n",
           n, 100.0 * frac, sout == 7 ? "hit" : "miss", secs);
}

static void test_pinned_big_budget(void) {
    int prepared[MAX_CIPHER_LENGTH]; char cstr[MAX_CIPHER_LENGTH];
    int n = plant_condi("CIPHER", 0, 7, 360, prepared, cstr, NULL);
    double secs;
    // Pin the true starter AND grant a large budget: demonstrates no budget escapes the needle.
    double frac = solve_ex(cstr, prepared, n, 7, METHOD_DEFAULT, 12, 300000, 1, NULL, &secs);
    printf("[pinned+big] n=%d starter=7 12x300000 recovery=%.1f%% [%.1fs] -- expected LOW\n",
           n, 100.0 * frac, secs);
}

static void test_per_scheme(void) {
    int prepared[MAX_CIPHER_LENGTH]; char cstr[MAX_CIPHER_LENGTH];
    int n = plant_condi("CIPHER", 0, 7, 180, prepared, cstr, NULL);
    struct { int m; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal" }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso" } };
    for (int i = 0; i < 3; i++) {
        double secs;
        double frac = solve_ex(cstr, prepared, n, 7, M[i].m, 0, 0, 1, NULL, &secs);
        printf("[per-scheme] %-7s starter=7 recovery=%.1f%% [%.1fs]\n",
               M[i].name, 100.0 * frac, secs);
    }
}

int main(void) {
    init_alphabet(NULL);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, true);   // -logprob table
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    test_registry();
    test_needle();
    test_blind_characterization();
    test_pinned_big_budget();
    test_per_scheme();

    printf("%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
