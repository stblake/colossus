//
//  In-process stress / limits tests for the Gronsfeld solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each stochastic solve
//  deterministic.
//
//  Strategy (planted-cipher recovery): encipher a known English plaintext under a known
//  numeric Gronsfeld key, then attack it ciphertext-only. The suite does three things:
//    1. confirms Gronsfeld inherits the polyalphabetic search defaults -- it has NO
//       per-type registry entry, so apply_cipher_defaults must leave the config untouched
//       (regression safety, exactly like Vigenere);
//    2. asserts the capability FLOOR -- a Gronsfeld with its period ESTIMATED end-to-end
//       (the realistic ciphertext-only attack) is recovered to ~100%, and again with the
//       period pinned;
//    3. characterizes the LIMIT -- recovery vs ciphertext length is printed (period pinned)
//       so the short-text behaviour is visible.
//
//  Gronsfeld is solved by the same optimal-cycleword frequency attack as Vigenere, with the
//  per-column shift search bounded to the digit domain 0..9; like Vigenere it does NOT need
//  the log-probability fitness, so this runs on the default reward-only quadgram table.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "engine.h"        // apply_cipher_defaults
#include "scoring.h"       // load_ngrams

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A long chunk of natural English (Pride and Prejudice, opening), letters only.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNE"
    "MUSTBEINWANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHIS"
    "FIRSTENTERINGANEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDING"
    "FAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGHTERS"
    "MYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLETATLAST"
    "MRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASJUSTBEENHEREANDSHE"
    "TOLDMEALLABOUTITMRBENNETMADENOANSWERDOYOUNOTWANTTOKNOWWHOHASTAKENITCRIEDHISWIFE"
    "IMPATIENTLYYOUWANTTOTELLMEANDIHAVENOOBJECTIONTOHEARINGITTHISWASINVITATIONENOUGH"
    "WHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATNETHERFIELDISTAKENBYAYOUNGMANOFLARGEFORTUNE"
    "FROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAISEANDFOURTOSEETHEPLACEANDWAS";

// A..Z char -> 0..25 alphabet index (full 26-letter alphabet).
static int letter_to_index(int c) {
    c = toupper(c);
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Plant a Gronsfeld cipher: take the first pt_len letters of PLAINTEXT, encipher with the
// numeric key `digitkey` (a string of digits 0..9). Fills prepared[] (the expected
// solution) and cipher_str (the A..Z ciphertext). Returns the length; *period_out gets the
// true key length.
static int plant(const char *digitkey, int pt_len, int *period_out,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }

    int key[64], keylen = 0;
    for (int i = 0; digitkey[i] && keylen < 64; i++)
        if (digitkey[i] >= '0' && digitkey[i] <= '9') key[keylen++] = digitkey[i] - '0';

    int cipher[MAX_CIPHER_LENGTH];
    gronsfeld_encrypt(cipher, prepared, n, key, keylen);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    if (period_out) *period_out = keylen;
    return n;
}

// Run one solve and return the recovered fraction. `period` == 0 ESTIMATES the cycleword
// length (IoC); >0 pins it. If period_out is non-NULL it receives the period the solver
// reported (the recovered cycleword length).
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int n_restarts, int n_hillclimbs, int period, uint32_t seed,
        int *period_out, double *secs_out) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = GRONSFELD;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    strcpy(cfg.ciphertext_file, "in-process-test");
    if (period > 0) { cfg.cycleword_len_present = true; cfg.cycleword_len = period; }
    cfg.n_restarts = n_restarts;
    cfg.n_hill_climbs = n_hillclimbs;

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

// --- 1. registry / defaults ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    // Gronsfeld rides the polyalphabetic pipeline and has no tuned schedule, so (like
    // Vigenere) apply_cipher_defaults must report "no entry" and leave the config alone.
    init_config(&cfg);
    cfg.cipher_type = GRONSFELD;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "gronsfeld should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
    CHECK(cfg.optimal_cycleword, "gronsfeld should default to the optimal-cycleword attack");
}

// --- 2. capability floor ------------------------------------------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    int true_period = 0;

    // Period ESTIMATED end-to-end (the realistic attack): plant a period-8 key over ~400
    // characters; the solver must estimate the period (or a multiple) and recover the text.
    int plen = plant("31415926", 400, &true_period, prepared, cipher_str);
    int got_p = -1; double secs = 0.0;
    double frac = solve_and_frac(cipher_str, prepared, plen, 200, 500, 0, 1u, &got_p, &secs);
    printf("[capability period-estimated ~400ch, true period %d] frac=%.3f period=%d %.1fs\n",
        true_period, frac, got_p, secs);
    CHECK(frac >= 0.99, "gronsfeld (~400ch, period estimated) recovered only %.3f", frac);
    CHECK(got_p > 0 && got_p % true_period == 0,
        "gronsfeld recovered period %d is not a multiple of the true %d", got_p, true_period);

    // Same cipher with the period PINNED: must also recover ~fully.
    double frac_pin = solve_and_frac(cipher_str, prepared, plen, 200, 500, true_period, 1u, NULL, &secs);
    printf("[capability period-pinned ~400ch] frac=%.3f %.1fs\n", frac_pin, secs);
    CHECK(frac_pin >= 0.99, "gronsfeld (~400ch, period pinned) recovered only %.3f", frac_pin);
}

// --- 3. length cliff (period pinned) ------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {80, 120, 200, 300};
    printf("recovery vs ciphertext length (key 31415926, period 8 pinned, 200x500, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int true_period = 0;
        int plen = plant("31415926", lengths[li], &true_period, prepared, cipher_str);
        double secs = 0.0;
        double frac = solve_and_frac(cipher_str, prepared, plen, 200, 500, true_period, 1u, NULL, &secs);
        printf("    len~%-4d  frac=%.3f  %.1fs\n", lengths[li], frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.99, "gronsfeld 300ch (cliff sweep) recovered only %.3f", frac_longest);
}

int main(void) {
    init_alphabet(NULL);                 // full 26-letter alphabet
    CHECK(g_alpha == ALPHABET_SIZE, "alphabet size %d, expected 26", g_alpha);

    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability();
    test_length_cliff();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
