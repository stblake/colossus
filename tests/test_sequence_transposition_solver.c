//
//  In-process capability + calibration tests for the Sequence Transposition solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN and
//  this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Sequence Transposition is a COLUMNAR transposition over 10 digit columns whose fill is driven
//  by a Gromark chain-addition primer. Two regimes:
//    * PRIMER GIVEN (-primer, the ACA convention transmits the primer): the only unknown is the
//      10-bucket read-order permutation -- a small, reliable transposition climb. This is the
//      calibrated path: all three -method schedules recover it. Needs -logprob (the interleaving
//      makes reward-only quadgrams gameable at practical budgets, like the square/Morse types).
//    * BLIND (no primer): a pre-pass ranks the 10^5 primer space. This is a DOCUMENTED LIMITATION
//      -- every wrong primer gets to pick its own best-scoring (gamed) read order, so with 10^5
//      alternatives the true primer does not reliably achieve the global n-gram maximum even at
//      300 letters. Characterized here, not asserted.
//
//  Checks:
//    1. no per-type SearchDefaults entry (uses the global transposition anneal schedule);
//    2. primer-given exact recovery across keywords at a solid length (asserted 100%);
//    3. -method calibration (default/shotgun/anneal/pso) on the primer-given path (asserted);
//    4. primer-given recovery-vs-length cliff (printed; asserted floor);
//    5. blind primer recovery characterization (printed; not asserted -- the gaming limitation);
//    6. determinism.
//
//  Runs with -logprob on the quadgram table; run from the source directory.
//

#include "colossus.h"
#include "engine.h"                      // apply_cipher_defaults
#include "scoring.h"                     // load_ngrams
#include "sequence_transposition_solver.h"

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
    "WHENINTHECOURSEOFHUMANEVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICALBANDS"
    "WHICHHAVECONNECTEDTHEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATEANDEQUAL"
    "STATIONTOWHICHTHELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHEOPINIONSOF"
    "MANKINDREQUIRESTHATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATIONWEHOLDTHESE";

// 10-letter keywords (any 10 letters work; ranks of 10 items always form a permutation of 1..10)
// paired with 5-digit primers.
static const char *KEYWORDS[] = { "BACKGROUND", "CRYPTOGRAM", "GUMMYBEARS", "NIGHTMARES" };
static const char *PRIMERS[]  = { "31415",      "92653",      "69315",      "58979" };

// Encipher the first pt_len letters of PLAINTEXT under keyword+primer into prepared[] (expected
// solution) and cipher_str (bare-letter ciphertext). Returns the length.
static int plant(const char *keyword, const char *primerstr, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int pi[SEQ_TRANS_BUCKETS];
    sequence_transposition_pi_from_keyword(keyword, pi);
    int primer[SEQ_TRANS_MAX_PRIMER], P = 0;
    for (int i = 0; primerstr[i]; i++)
        if (primerstr[i] >= '0' && primerstr[i] <= '9') primer[P++] = primerstr[i] - '0';
    static int cipher[MAX_CIPHER_LENGTH];
    sequence_transposition_encrypt(prepared, n, primer, P, pi, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Drive solve_cipher and return the recovered fraction. primer_str != NULL pins the primer
// (-primer); NULL runs the blind pre-pass.
static double solve_frac(const char *cipher_str, const int prepared[], int plen,
        const char *primer_str, int method, int nrestarts, int nhill, int nparticles,
        uint32_t seed, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = SEQUENCE_TRANSPOSITION;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    cfg.n_restarts = nrestarts;
    cfg.n_hill_climbs = nhill;
    if (nparticles > 0) cfg.n_particles = nparticles;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (primer_str) {
        int P = 0;
        for (int i = 0; primer_str[i] && P < SEQ_TRANS_MAX_PRIMER; i++)
            if (primer_str[i] >= '0' && primer_str[i] <= '9') cfg.seq_primer[P++] = primer_str[i] - '0';
        cfg.seq_primer_len = P;
    }

    SolveResult res;
    res.solved = false;
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

// --- 1. registry: no per-type entry (uses the global transposition anneal schedule) ----------

static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = SEQUENCE_TRANSPOSITION; cfg.method = METHOD_DEFAULT;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false),
          "sequence-transposition should have no registry entry (uses globals)");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "globals were modified");
}

// --- 2. primer-given exact recovery across keywords at L=140 ----------------------------------

static int nkw(void) { return (int) (sizeof KEYWORDS / sizeof KEYWORDS[0]); }

static void test_primer_given_recovery(void) {
    int ok = 0, total = 0;
    for (int ki = 0; ki < nkw(); ki++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(KEYWORDS[ki], PRIMERS[ki], 140, prepared, cipher_str);
        double f = solve_frac(cipher_str, prepared, plen, PRIMERS[ki],
                              METHOD_DEFAULT, 40, 20000, 0, 1u, NULL);
        total++; ok += (f > 0.999);
    }
    CHECK(ok == total, "primer-given exact recovery at L=140: %d/%d", ok, total);
}

// --- 3. -method calibration on the primer-given path (L=140) ----------------------------------

static void test_method_calibration(void) {
    struct { const char *name; int method, nr, nh, np; } M[] = {
        { "default", METHOD_DEFAULT, 40, 20000, 0  },
        { "shotgun", METHOD_SHOTGUN, 40, 20000, 0  },
        { "anneal",  METHOD_ANNEAL,  40, 20000, 0  },
        { "pso",     METHOD_PSO,      3,  1500, 16 },
    };
    int nm = (int) (sizeof M / sizeof M[0]);
    printf("\n  -method calibration on the primer-given path (mean exact over %d keywords, L=140):\n",
           nkw());
    double def_acc = -1.0;
    for (int mi = 0; mi < nm; mi++) {
        double sum = 0.0; int cnt = 0;
        for (int ki = 0; ki < nkw(); ki++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cipher_str[MAX_CIPHER_LENGTH];
            int plen = plant(KEYWORDS[ki], PRIMERS[ki], 140, prepared, cipher_str);
            double f = solve_frac(cipher_str, prepared, plen, PRIMERS[ki],
                                  M[mi].method, M[mi].nr, M[mi].nh, M[mi].np, 1u, NULL);
            sum += (f > 0.999) ? 1.0 : 0.0; cnt++;
        }
        double mean = sum / cnt;
        printf("    %-8s exact %.2f\n", M[mi].name, mean);
        if (M[mi].method == METHOD_DEFAULT) def_acc = mean;
        // Every schedule is expected to solve this small permutation; assert each recovers all.
        CHECK(mean > 0.999, "-method %s recovered only %.2f at L=140", M[mi].name, mean);
    }
    CHECK(def_acc > 0.999, "default method recovered only %.2f", def_acc);
}

// --- 4. primer-given recovery-vs-length cliff -------------------------------------------------

static void test_length_cliff(void) {
    int lens[] = {60, 80, 100, 120, 140};
    int nlens = (int) (sizeof lens / sizeof lens[0]);
    printf("\n  primer-given recovery vs length (mean exact over %d keywords, -method default):\n",
           nkw());
    double acc120 = 0.0, acc140 = 0.0;
    for (int li = 0; li < nlens; li++) {
        int L = lens[li];
        double sum = 0.0; int cnt = 0;
        for (int ki = 0; ki < nkw(); ki++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cipher_str[MAX_CIPHER_LENGTH];
            int plen = plant(KEYWORDS[ki], PRIMERS[ki], L, prepared, cipher_str);
            double f = solve_frac(cipher_str, prepared, plen, PRIMERS[ki],
                                  METHOD_DEFAULT, 60, 20000, 0, 1u, NULL);
            sum += (f > 0.999) ? 1.0 : 0.0; cnt++;
        }
        double mean = sum / cnt;
        printf("    len %4d: exact %.2f\n", L, mean);
        if (L == 120) acc120 = mean;
        if (L == 140) acc140 = mean;
    }
    CHECK(acc120 > 0.999, "primer-given recovery at L=120: %.2f (want all)", acc120);
    CHECK(acc140 > 0.999, "primer-given recovery at L=140: %.2f (want all)", acc140);
}

// --- 5. blind primer recovery characterization (NOT asserted -- gaming limitation) ------------

static void test_blind_characterization(void) {
    printf("\n  blind (no -primer) recovery -- DOCUMENTED LIMITATION (10^5 primers each game the\n"
           "  n-gram, so the true primer rarely wins the global max; recovery is noisy and\n"
           "  non-monotonic in length/budget). Characterized on 2 instances, NOT asserted:\n");
    int lens[] = {140, 290};                             // two lengths, one keyword each (blind ~15s)
    int nlens = (int) (sizeof lens / sizeof lens[0]);
    for (int li = 0; li < nlens; li++) {
        int L = lens[li];
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(KEYWORDS[li], PRIMERS[li], L, prepared, cipher_str);
        double secs = 0.0;
        double f = solve_frac(cipher_str, prepared, plen, NULL, METHOD_DEFAULT, 6, 12000, 0, 1u, &secs);
        printf("    %-10s len ~%-4d exact %s  (%.1fs)\n",
               KEYWORDS[li], L, f > 0.999 ? "yes" : "no", secs);
    }
}

// --- 6. determinism ---------------------------------------------------------------------------

static void test_determinism(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant(KEYWORDS[0], PRIMERS[0], 140, prepared, cipher_str);
    double a = solve_frac(cipher_str, prepared, plen, PRIMERS[0], METHOD_DEFAULT, 40, 20000, 0, 7u, NULL);
    double b = solve_frac(cipher_str, prepared, plen, PRIMERS[0], METHOD_DEFAULT, 40, 20000, 0, 7u, NULL);
    CHECK(a == b, "not deterministic at fixed seed: %.4f vs %.4f", a, b);
}

int main(void) {
    g_ngram_logprob = true;                      // sequence transposition needs -logprob
    init_alphabet(NULL);                         // full 26-letter alphabet
    CHECK(g_alpha == ALPHABET_SIZE, "alphabet size %d, expected 26", g_alpha);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_primer_given_recovery();
    test_method_calibration();
    test_length_cliff();
    test_blind_characterization();
    test_determinism();

    free(shared.ngram_data);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
