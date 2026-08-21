//
//  In-process stress / limits tests for the Tri-Square solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with
//  -DCOLOSSUS_NO_MAIN and this file supplies its own main, so solve_cipher can be driven
//  directly and its SolveResult inspected. A fixed -seed makes each solve deterministic.
//
//  Tri-Square is the largest square state of the family: THREE independent keyed 5x5 squares
//  (75 cells), jointly annealed with no decoupling reward. But its polyphonic cipher letters
//  (c0/c2, chosen at random from a column/row) spread the full alphabet over every position,
//  so the n-gram gradient is sharp and it recovers MORE easily than Four-Square -- reliable
//  from ~500 plaintext letters (a 300-400 cliff). The suite (planted-cipher recovery) does:
//    1. validates the per-type schedule registry (and that a non-registry type is untouched);
//    2. asserts the capability FLOOR -- a long cipher is recovered to ~100% at the registry budget;
//    3. characterizes the LIMIT -- recovery vs plaintext length is printed across a sweep
//       (the short-text cliff is visible; the marginal lengths are reported, not asserted);
//    4. extra-thorough: a first-of-pair / second-of-pair recovery breakdown (even positions are
//       the first plaintext letter of each digraph, decrypted through sq1; odd positions the
//       second, through sq2), plus a MULTI-KEYWORD sweep reporting mean / worst recovery.
//
//  Run from the source directory so the n-gram table is found. Tri-Square needs the
//  discriminating log-probability fitness, so this enables g_ngram_logprob. Because the encode
//  is polyphonic (RNG-driven), plant() seeds the RNG to a fixed value so each cipher is
//  deterministic regardless of the preceding solve's RNG consumption.
//

#include "colossus.h"
#include "engine.h"
#include "scoring.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// A long chunk of natural English (Pride and Prejudice, opening), letters only (~1369) -- long
// enough for the length sweep to reach well past the ~500-char capability floor.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEM"
    "USTBEINWANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHIS"
    "FIRSTENTERINGANEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDIN"
    "GFAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGH"
    "TERSMYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISL"
    "ETATLASTMRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASJUSTBEEN"
    "HEREANDSHETOLDMEALLABOUTITMRBENNETMADENOANSWERDOYOUNOTWANTTOKNOWWHOHASTAKE"
    "NITCRIEDHISWIFEIMPATIENTLYYOUWANTTOTELLMEANDIHAVENOOBJECTIONTOHEARINGITTHI"
    "SWASINVITATIONENOUGHWHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATNETHERFIELDISTAKENB"
    "YAYOUNGMANOFLARGEFORTUNEFROMTHENORTHOFENGLANDTHATHECAMEDOWNONMONDAYINACHAI"
    "SEANDFOURTOSEETHEPLACEANDWASSOMUCHDELIGHTEDWITHITTHATHEAGREEDWITHMRMORRISI"
    "MMEDIATELYTHATHEISTOTAKEPOSSESSIONBEFOREMICHAELMASANDSOMEOFHISSERVANTSARET"
    "OBEINTHEHOUSEBYTHEENDOFNEXTWEEKWHATISHISNAMEBINGLEYISHEMARRIEDORSINGLEOHSI"
    "NGLEMYDEARTOBESUREASINGLEMANOFLARGEFORTUNEFOURORFIVETHOUSANDAYEARWHATAFINE"
    "THINGFOROURGIRLSHOWSOHOWCANITAFFECTTHEMMRBENNETHOWCANYOUBESOTIRESOMEYOUMUS"
    "TKNOWTHATIAMTHINKINGOFHISMARRYINGONEOFTHEMISTHATHISDESIGNINSETTLINGHEREDES"
    "IGNNONSENSEHOWCANYOUTALKSOBUTITISVERYLIKELYTHATHEMAYFALLINLOVEWITHONEOFTHE"
    "MANDTHEREFOREYOUMUSTVISITHIMASSOONASHECOMESISEENOOCCASIONFORTHATYOUANDTHEG"
    "IRLSMAYGOORYOUMAYSENDTHEMBYTHEMSELVES";

static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Build one keyed square from a keyword (J->I) via the shared Playfair grid builder.
static void build_square(const char *kwstr, int sq[]) {
    int kw[64], kwn = 0;
    for (int i = 0; kwstr[i] && kwn < 64; i++) { int x = letter_to_index((unsigned char) kwstr[i]); if (x >= 0) kw[kwn++] = x; }
    playfair_grid_from_keyword(kw, kwn, sq);
}

// Plant a Tri-Square cipher: first pt_len letters of PLAINTEXT (padded to even with X),
// enciphered under three keyword squares (polyphonic, so seed the RNG for reproducibility).
// Fills prepared[] (the expected solution, length 2M) and cipher_str (length 3M). Returns the
// prepared (plaintext / scoring) length.
static int plant(const char *k1, const char *k2, const char *k3, int pt_len,
                 int prepared[], char cipher_str[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int idx = letter_to_index((unsigned char) PLAINTEXT[i]);
        if (idx >= 0) prepared[n++] = idx;
    }
    if (n % 2 != 0) prepared[n++] = letter_to_index('X');
    int sq1[SQUARE_GRID], sq2[SQUARE_GRID], sq3[SQUARE_GRID];
    build_square(k1, sq1); build_square(k2, sq2); build_square(k3, sq3);
    int cipher[MAX_CIPHER_LENGTH];
    seed_rand(0x77213u);                     // deterministic polyphonic choices
    int clen = trisquare_encrypt(prepared, n, sq1, sq2, sq3, SQUARE_SIDE, cipher);
    for (int i = 0; i < clen; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[clen] = '\0';
    return n;
}

// Run one solve; return the recovered fraction. If recovered != NULL, copy the decrypted
// plaintext into it (for the first/second-of-pair breakdown).
static double solve_and_frac(const char *cipher_str, const int prepared[], int plen,
        int method, bool use_registry, int n_restarts, int n_hillclimbs,
        double init_temp, double backtrack, uint32_t seed, double *secs_out, int *recovered) {

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TRI_SQUARE;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");

    if (use_registry) {
        apply_cipher_defaults(&cfg, false);
    } else {
        cfg.n_restarts = n_restarts;
        cfg.n_hill_climbs = n_hillclimbs;
        cfg.init_temp = init_temp;
        cfg.backtracking_probability = backtrack;
    }

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
    for (int i = 0; i < plen; i++) {
        if (recovered) recovered[i] = res.decrypted[i];
        if (res.decrypted[i] == prepared[i]) ok++;
    }
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    ColossusConfig cfg;

    init_config(&cfg);
    cfg.cipher_type = TRI_SQUARE; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "trisquare registry: no entry applied");
    CHECK(cfg.n_restarts == 12 && cfg.n_hill_climbs == 500000,
        "trisquare anneal defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0799 && cfg.init_temp < 0.0801,
        "trisquare anneal inittemp wrong: %.4f", cfg.init_temp);

    init_config(&cfg);
    cfg.cipher_type = TRI_SQUARE; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "trisquare registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 40 && cfg.n_hill_climbs == 400000,
        "trisquare shotgun defaults wrong: %dx%d", cfg.n_restarts, cfg.n_hill_climbs);

    // Regression safety: a type with no registry entry is left untouched.
    init_config(&cfg);
    cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0,
        "non-registry type was modified by apply_cipher_defaults");
}

// --- 2. capability floor + first/second-of-pair breakdown ---------------------

static void test_capability(void) {
    static int prepared[MAX_CIPHER_LENGTH];
    static char cipher_str[MAX_CIPHER_LENGTH];
    static int recovered[MAX_CIPHER_LENGTH];
    double secs;

    int plen = plant("CIPHERKEYWORD", "PALIMPSEST", "NORTHERLY", 900, prepared, cipher_str);
    double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, true,
        0, 0, 0.0, 0.0, 1u, &secs, recovered);
    printf("[capability anneal ~900pt (registry 12x500k)] frac=%.3f %.1fs\n", frac, secs);
    CHECK(frac >= 0.99, "trisquare ~900pt recovered only %.3f", frac);

    // Breakdown: even positions are the first plaintext letter of each digraph (decrypted
    // through sq1), odd positions the second (through sq2). High recovery on both confirms
    // both the sq1-side and sq2-side of the joint state were found.
    int ev_ok = 0, ev = 0, od_ok = 0, od = 0;
    for (int i = 0; i < plen; i++) {
        if (i % 2 == 0) { ev++; if (recovered[i] == prepared[i]) ev_ok++; }
        else            { od++; if (recovered[i] == prepared[i]) od_ok++; }
    }
    double fe = ev ? (double) ev_ok / ev : 0.0, fo = od ? (double) od_ok / od : 0.0;
    printf("[per-side] first-of-pair (sq1, even pos)=%.3f  second-of-pair (sq2, odd pos)=%.3f\n", fe, fo);
    CHECK(fe >= 0.99 && fo >= 0.99, "trisquare per-side recovery sq1=%.3f sq2=%.3f", fe, fo);
}

// --- 3. length cliff ----------------------------------------------------------

static void test_length_cliff(void) {
    int lengths[] = {400, 600, 900, 1300};
    printf("recovery vs plaintext length (trisquare, keys TRISQUARE/CRYPTO/BALTIMORE, 8x400k, seed 1):\n");
    double frac_longest = 0.0;
    for (int li = 0; li < 4; li++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant("TRISQUARE", "CRYPTO", "BALTIMORE", lengths[li], prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
            8, 400000, 0.08, 0.30, 1u, &secs, NULL);
        printf("    pt~%-4d (%d trigraphs)  frac=%.3f  %.1fs\n", lengths[li], plen / 2, frac, secs);
        if (li == 3) frac_longest = frac;
    }
    CHECK(frac_longest >= 0.99, "trisquare 1300pt (cliff sweep) recovered only %.3f", frac_longest);
}

// --- 4. multi-keyword sweep ---------------------------------------------------

static void test_multi_keyword(void) {
    const char *keys[][3] = {
        {"KRYPTOSABC", "PALIMPSEST", "NORTHERLY"},
        {"SHADOWFORCE", "IODINEGAS", "ABANDONED"},
        {"NORTHERLY",   "BERLINCLOCK", "EASTWARD"},
    };
    int nk = 3;
    double sum = 0.0, worst = 1.0;
    printf("multi-keyword sweep (trisquare, ~900pt, 8x400k, seed 1):\n");
    for (int k = 0; k < nk; k++) {
        static int prepared[MAX_CIPHER_LENGTH];
        static char cipher_str[MAX_CIPHER_LENGTH];
        int plen = plant(keys[k][0], keys[k][1], keys[k][2], 900, prepared, cipher_str);
        double secs;
        double frac = solve_and_frac(cipher_str, prepared, plen, METHOD_DEFAULT, false,
            8, 400000, 0.08, 0.30, 1u, &secs, NULL);
        printf("    keys=%-11s/%-11s/%-11s frac=%.3f %.1fs\n",
            keys[k][0], keys[k][1], keys[k][2], frac, secs);
        sum += frac;
        if (frac < worst) worst = frac;
    }
    printf("    mean=%.3f worst=%.3f\n", sum / nk, worst);
    CHECK(sum / nk >= 0.95, "trisquare multi-keyword mean recovery only %.3f", sum / nk);
}

int main(void) {
    init_alphabet("J");                 // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == SQUARE_GRID, "alphabet size %d, expected %d", g_alpha, SQUARE_GRID);

    g_ngram_logprob = true;             // Tri-Square needs the log-probability fitness
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) {
        printf("FAIL: could not load %s (run from the source directory)\n", NGRAM_FILE);
        return 1;
    }

    test_registry();
    test_capability();
    test_length_cliff();
    test_multi_keyword();

    free(shared.ngram_data);

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}
