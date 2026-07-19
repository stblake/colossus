//
//  Solver calibration / capability tests for the Grandpre cipher (type 84). Built with
//  -DCOLOSSUS_NO_MAIN so solve_cipher() links directly; fixed seeds -> deterministic. Covers
//  registry validation and a letter-capability PROBE across lengths.
//
//  Grandpre decoding is code -> letter (unique), so cryptanalytically it is a HOMOPHONIC map
//  over <= N^2 numeric codes -> 26 letters. An ~64-code (8x8) map over the ACA 150-200 range is
//  undersampled (~2-3 samples/code), so recovery is PARTIAL and rare-letter homophones (Q/X/Z)
//  are the last to resolve -- a documented characterization, not a 99% bar. This test asserts a
//  lenient best-of-seeds capability floor (breakage catch) and PRINTS the accuracy-vs-length
//  curve. Effectively needs -logprob with QUINTGRAMS.
//
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "grandpre.h"
#include "grandpre_solver.h"
#include <fcntl.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE (char *) "english_quintgrams.txt"
#define NGRAM_SIZE 5

static SharedData shared;

// ~450 letters of English (A Tale of Two Cities opening), letters only.
static const char *PLAINTEXT =
    "ITWASTHEBESTOFTIMESITWASTHEWORSTOFTIMESITWASTHEAGEOFWISDOMITWASTHEAGEOFFOOLISHNESS"
    "ITWASTHEEPOCHOFBELIEFITWASTHEEPOCHOFINCREDULITYITWASTHESEASONOFLIGHTITWASTHESEASON"
    "OFDARKNESSITWASTHESPRINGOFHOPEITWASTHEWINTEROFDESPAIRWEHADEVERYTHINGBEFOREUSWEHADNO"
    "THINGBEFOREUSWEWEREALLGOINGDIRECTTOHEAVENWEWEREALLGOINGDIRECTTHEOTHERWAY";

// First L letters into out[] (0..25); returns the count.
static int text_of(const char *s, int L, int out[]) {
    int n = 0;
    for (int i = 0; s[i] && n < L; i++)
        if (s[i] >= 'A' && s[i] <= 'Z') out[n++] = s[i] - 'A';
    return n;
}

// A random valid N x N square (all 26 letters present), reproducible from `sqseed`.
static void make_square(GrandpreSquare *g, int N, uint32_t sqseed) {
    seed_rand(sqseed);
    int M = N * N, letters[GRANDPRE_MAX_GRID];
    for (int i = 0; i < 26; i++) letters[i] = i;
    for (int i = 26; i < M; i++) letters[i] = rand_int(0, 26);
    for (int i = M - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = letters[i]; letters[i] = letters[j]; letters[j] = t; }
    grandpre_build(g, N, letters);
}

// Plant a Grandpre cipher (random homophones under `plantseed`) and solve it in-process; return
// the letter match fraction.
static double plant_and_solve(const int *pt, int n, int N, uint32_t sqseed,
                              uint32_t plantseed, uint32_t solveseed) {
    GrandpreSquare g; make_square(&g, N, sqseed);
    static int codes[MAX_CIPHER_LENGTH];
    seed_rand(plantseed);
    int clen = grandpre_encrypt(pt, n, &g, codes);
    static char cs[2 * MAX_CIPHER_LENGTH + 1];
    int p = 0;
    for (int i = 0; i < clen; i++) { cs[p++] = '0' + codes[i] / 10; cs[p++] = '0' + codes[i] % 10; }
    cs[p] = '\0';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = GRANDPRE;
    cfg.ngram_size = NGRAM_SIZE;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    cfg.n_restarts = 8; cfg.n_hill_climbs = 90000;   // reduced budget: keeps the test fast (~30s)

    SolveResult res;
    res.solved = false; res.decrypted_len = 0;
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* proceed anyway */ }
    seed_rand(solveseed);
    solve_cipher(cs, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);

    int lim = (res.decrypted_len < n) ? res.decrypted_len : n;
    int ok = 0;
    for (int i = 0; i < lim; i++) if (res.decrypted[i] == pt[i]) ok++;
    return (n > 0) ? (double) ok / n : 0.0;
}

int main(void) {
    init_alphabet(NULL);                        // full 26-letter A..Z, as -type grandpre expects
    CHECK(g_alpha == DEFAULT_ALPHABET_SIZE, "alphabet %d, expected %d", g_alpha, DEFAULT_ALPHABET_SIZE);
    g_ngram_logprob = true;                     // Grandpre effectively needs -logprob (quintgrams)
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    static int pt[MAX_CIPHER_LENGTH];

    // 1) Registry validation: the tuned annealed-square schedule is applied; a non-registry
    //    type (VIGENERE) is left untouched.
    {
        ColossusConfig cfg; init_config(&cfg);
        cfg.cipher_type = GRANDPRE;
        bool applied = apply_cipher_defaults(&cfg, false);
        CHECK(applied, "no SearchDefaults entry for GRANDPRE");
        CHECK(cfg.init_temp > 0.1499 && cfg.init_temp < 0.1501, "unexpected inittemp %.4f", cfg.init_temp);
        CHECK(cfg.n_restarts >= 12, "unexpected a_n_restarts %d", cfg.n_restarts);
        ColossusConfig v; init_config(&v); v.cipher_type = VIGENERE;
        CHECK(!apply_cipher_defaults(&v, false), "VIGENERE unexpectedly has a registry entry");
    }

    // 2) LETTER capability PROBE across lengths (best of a few seeds -- recovery is high-variance
    //    at the short end). Asserts only a lenient breakage floor; prints the curve. Rare-letter
    //    homophones are the last to resolve, so exact recovery is capped below 99% at ACA lengths.
    printf("\n  Grandpre letter recovery vs length (random 8x8 square, best of 2 seeds, quintgrams):\n");
    {
        int lens[] = {90, 150, 200, 285};
        double best_at_long = 0.0;
        for (int li = 0; li < 4; li++) {
            int n = text_of(PLAINTEXT, lens[li], pt);
            double best = 0.0;
            for (uint32_t s = 1; s <= 2; s++) {
                double f = plant_and_solve(pt, n, 8, 11u, 7u, s);
                if (f > best) best = f;
            }
            printf("    L=%3d letters: best=%.3f\n", n, best);
            if (lens[li] >= 200 && best > best_at_long) best_at_long = best;
        }
        CHECK(best_at_long > 0.60, "capability floor: best recovery %.3f at L>=200 (want > 0.60)", best_at_long);
        printf("    (partial + rare-letter-capped by design: homophonic over ~64 codes, below the 99%% floor at ACA lengths)\n");
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(shared.ngram_data);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
