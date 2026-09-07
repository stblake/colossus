//
//  In-process stress / limits tests for the Gromark / Periodic Gromark solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  Gromark is a COUPLED keyed-alphabet + chain-running-key cipher with NO clean decoupling
//  reward, so it is attacked by a PRIMER PRE-PASS (an assignment-based frequency attack ranks
//  the finite primer space; the top-K primers are annealed). This suite checks:
//    1. registry validation (apply_cipher_defaults) for both codes + a non-registry type;
//    2. pre-pass hit rate: gromark_rank_primers lands the true primer in the top-K, vs length;
//    3. capability floor + length cliff (basic, blind end to end);
//    4. a Periodic Gromark blind solve (period swept).
//
//  Gromark runs on the full 26-letter alphabet and needs the log-probability fitness. Run from
//  the source directory (loads ngram_data/english/english_quadgrams.txt).
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include "gromark_solver.h"       // gromark_rank_primers

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
    "THEUNANIMOUSDECLARATIONOFTHETHIRTEENUNITEDSTATESOFAMERICAWHENINTHECOURSEOFHUMAN"
    "EVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICALBANDSWHICHHAVECONNECTED"
    "THEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATEANDEQUALSTATIONTOWHICH"
    "THELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHEOPINIONSOFMANKINDREQUIRES"
    "THATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATIONWEHOLDTHESETRUTHSTOBESELF"
    "EVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIRCREATORWITHCERTAINUNALIENABLE";

// Plant a basic Gromark cipher: first pt_len plaintext letters under the K2M alphabet of
// `keyword` and the digit `primer`. Fills prepared[] (expected solution) and cipher_str
// (the bare-letter ciphertext). Returns the length.
static int plant(const char *keyword, const char *primerstr, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int sigma[26];
    gromark_mixed_alphabet(keyword, sigma);
    int primer[26], P = 0;
    for (int i = 0; primerstr[i]; i++)
        if (primerstr[i] >= '0' && primerstr[i] <= '9') primer[P++] = primerstr[i] - '0';
    static int cipher[MAX_CIPHER_LENGTH];
    gromark_encrypt(prepared, n, sigma, primer, P, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

static double solve_and_frac(int type, const char *cipher_str, const int prepared[], int plen,
        int nprimers, int period, uint32_t seed, int *period_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (nprimers > 0) cfg.n_primers = nprimers;
    if (period > 0) { cfg.period_present = true; cfg.period = period; }

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
    init_config(&cfg); cfg.cipher_type = GROMARK; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "gromark registry: no entry applied");
    CHECK(cfg.n_restarts == 3 && cfg.n_hill_climbs == 120000,
        "gromark anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = GROMARK_PERIODIC; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "gromark-periodic registry: no entry applied");
    CHECK(cfg.n_restarts == 4 && cfg.n_hill_climbs == 160000,
        "gromark-periodic anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. pre-pass hit rate (basic): true primer in the top-K, vs length --------

static void test_prepass_hitrate(void) {
    const char *keywords[] = {"KRYPTOSABCDEF", "MERCURYXZ", "PALIMPSEST", "BERLINCLOCK"};
    const char *primers[]  = {"31415", "92653", "58979", "32384"};
    int lengths[] = {200, 150, 120};
    printf("Gromark pre-pass: true primer's presence in top-K (basic, K=64) vs length:\n");
    for (int li = 0; li < 3; li++) {
        int hits = 0, trials = 0;
        for (int ki = 0; ki < 4; ki++) {
            static int prepared[MAX_CIPHER_LENGTH];
            static char cipher_str[MAX_CIPHER_LENGTH];
            int plen = plant(keywords[ki], primers[ki], lengths[li], prepared, cipher_str);
            static int cidx[MAX_CIPHER_LENGTH];
            for (int i = 0; i < plen; i++) cidx[i] = cipher_str[i] - 'A';
            static int outp[64 * GROMARK_MAX_PRIMER], outper[64], outw[64 * (26 + GROMARK_MAX_PRIMER)];
            int K = gromark_rank_primers(cidx, plen, GROMARK, 5, 5, 5,
                shared.ngram_data, NGRAM_SIZE, 64, outper, outp, outw, false);
            int truep[5]; for (int i = 0; i < 5; i++) truep[i] = primers[ki][i] - '0';
            int found = 0;
            for (int k = 0; k < K && !found; k++) {
                int same = 1;
                for (int i = 0; i < 5; i++) if (outp[k * GROMARK_MAX_PRIMER + i] != truep[i]) same = 0;
                if (same) found = 1;
            }
            hits += found; trials++;
        }
        printf("    len~%-4d  %d/%d primers in top-64\n", lengths[li], hits, trials);
        if (lengths[li] >= 150) CHECK(hits == trials, "pre-pass missed %d/%d at len %d",
            trials - hits, trials, lengths[li]);
    }
}

// --- 3. capability floor + length cliff (basic, blind) ------------------------

static void test_capability_cliff(void) {
    int lengths[] = {120, 150, 200};
    printf("basic Gromark blind recovery vs length (sq=KRYPTOSABCDEF, primer=31415, seed 1):\n");
    double frac_long = 0.0;
    for (int li = 0; li < 3; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("KRYPTOSABCDEF", "31415", lengths[li], prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(GROMARK, cipher_str, prepared, plen, 0, 0, 1u, NULL, &secs);
        printf("    pt~%-4d  frac=%.3f  %.1fs\n", plen, frac, secs);
        if (li == 2) frac_long = frac;
    }
    CHECK(frac_long >= 0.95, "basic Gromark 200ch blind recovered only %.3f", frac_long);
}

// --- 4. periodic blind solve --------------------------------------------------

static int plant_periodic(const char *keyword, int pt_len, int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    int sigma[26];
    gromark_mixed_alphabet(keyword, sigma);
    // primer = alphabetical ranks of the distinct keyword letters; offsets = their positions.
    int keyed[26], seen[26] = {0}, P = 0;
    make_keyed_alphabet((char *) keyword, keyed);
    for (int i = 0; keyword[i]; i++) { int x = keyword[i] - 'A'; if (x >= 0 && x < 26 && !seen[x]) { seen[x] = 1; P++; } }
    int sinv[26]; for (int i = 0; i < 26; i++) sinv[sigma[i]] = i;
    int primer[26], offsets[26];
    for (int g = 0; g < P; g++) {
        int rank = 1; for (int h = 0; h < P; h++) if (keyed[h] < keyed[g]) rank++;
        primer[g] = rank; offsets[g] = sinv[keyed[g]];
    }
    static int cipher[MAX_CIPHER_LENGTH];
    gromark_periodic_encrypt(prepared, n, sigma, primer, P, offsets, cipher);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

static void test_periodic_blind(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int plen = plant_periodic("MERCURY", 200, prepared, cipher_str);   // 6 distinct -> period 6
    double secs = 0.0; int period = -1;
    double frac = solve_and_frac(GROMARK_PERIODIC, cipher_str, prepared, plen, 0, 0, 1u, &period, &secs);
    printf("[periodic blind %dch, period swept]  frac=%.3f period=%d (true 6)  %.1fs\n",
        plen, frac, period, secs);
    CHECK(frac >= 0.90, "periodic Gromark blind recovered only %.3f", frac);
}

int main(void) {
    g_ngram_logprob = true;
    init_alphabet(NULL);                // full 26-letter alphabet
    CHECK(g_alpha == ALPHABET_SIZE, "alphabet size %d, expected 26", g_alpha);
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_prepass_hitrate();
    test_capability_cliff();
    test_periodic_blind();

    free(shared.ngram_data);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
