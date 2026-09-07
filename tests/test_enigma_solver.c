// In-process solver regression for the Enigma type. Built with -DCOLOSSUS_NO_MAIN and linked
// against the whole solver so it calls solve_cipher() directly. Validates the SearchDefaults
// registry entry, characterises the ciphertext-only capability (recovery vs length/plugs with
// the wheel order pinned -- the realistic fast path), and checks the Bombe recovers a planted
// crib exactly. Run from the source directory so the n-gram table is found in the cwd.

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "enigma_solver.h"
#include "enigma.h"

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

static SharedData shared;

// A long English passage (Pride & Prejudice), letters only, used as planted plaintext. Must
// be at least as long as the largest planted `len` below (currently 450) plus margin -- plant()
// reads PLAINTEXT[0..len-1], so a short string would run past the terminator (a NUL enciphered
// as an out-of-range letter). Kept comfortably over 600 (the plant() pt[]/ct[] buffer size).
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEIN"
    "WANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGA"
    "NEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONS"
    "IDEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGHTERSMYDEARMRBENNETSAIDHISLAD"
    "YTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLETATLASTMRBENNETREPLIEDTHATHEHADNOT"
    "BUTITISRETURNEDSHEFORMRSLONGHASJUSTBEENHEREANDSHETOLDMEALLABOUTITMRBENNETMADENOAN"
    "SWERDOYOUNOTWANTTOKNOWWHOHASTAKENITCRIEDHISWIFEIMPATIENTLYYOUWANTTOTELLMEANDIHAVE"
    "NOOBJECTIONTOHEARINGITTHISWASINVITATIONENOUGHWHYMYDEARYOUMUSTKNOWMRSLONGSAYSTHATN"
    "ETHERFIELDISTAKENBYAYOUNGMANOFLARGEFORTUNEFROMTHENORTHOFENGLANDTHATHECAMEDOWNONMO";

// Build a random 3-rotor key with a pinned wheel order and `nplugs` disjoint steckers.
static void make_key(EnigmaKey *k, const int rotors[3], int nplugs, unsigned seed) {
    seed_rand(seed);
    k->reflector = ENIGMA_UKW_B;
    k->n_wheels = 3;
    for (int i = 0; i < 3; i++) { k->rotor[i] = rotors[i]; k->ring[i] = rand_int(0, 26); k->pos[i] = rand_int(0, 26); }
    enigma_plug_identity(k->plug);
    int avail[26]; for (int i = 0; i < 26; i++) avail[i] = i; int na = 26;
    for (int p = 0; p < nplugs && na >= 2; p++) {
        int i = rand_int(0, na); int a = avail[i]; avail[i] = avail[--na];
        int j = rand_int(0, na); int b = avail[j]; avail[j] = avail[--na];
        enigma_plug_set_pair(k->plug, a, b);
    }
}

// Plant: encrypt the first `len` letters of PLAINTEXT with key k.
static void plant(const EnigmaKey *k, int len, int *prepared, char *cipher_str) {
    int pt[600], ct[600];
    for (int i = 0; i < len; i++) pt[i] = PLAINTEXT[i] - 'A';
    for (int i = 0; i < len; i++) prepared[i] = pt[i];
    enigma_encrypt(pt, len, k, ct);
    for (int i = 0; i < len; i++) cipher_str[i] = 'A' + ct[i];
    cipher_str[len] = '\0';
}

// One in-process solve with the wheel order pinned; returns the plaintext recovery fraction.
// bombe => run the crib attack with a crib_len-letter crib at the message start.
static double solve_once(const int rotors[3], int len, int nplugs,
                         unsigned plant_seed, unsigned solve_seed, int bombe, int crib_len) {
    EnigmaKey k;
    make_key(&k, rotors, nplugs, plant_seed);
    int prepared[600];
    char cipher_str[601];
    plant(&k, len, prepared, cipher_str);

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = ENIGMA;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    cfg.enigma_rotors_present = true;
    for (int i = 0; i < 3; i++) cfg.enigma_rotors[i] = rotors[i];
    cfg.enigma_bombe = bombe;
    cfg.n_threads = 4;               // the phase-1 / Bombe searches parallelise over threads
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);

    char crib[601]; crib[0] = '\0';
    if (crib_len > 0) {
        for (int i = 0; i < len; i++) crib[i] = (i < crib_len) ? (char)('A' + prepared[i]) : '_';
        crib[len] = '\0';
    }

    SolveResult res;
    res.solved = false;
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* proceed anyway */ }
    seed_rand(solve_seed);
    solve_cipher(cipher_str, (crib_len > 0) ? crib : (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout)); close(saved); clearerr(stdout);

    if (!res.solved || res.decrypted_len != len) return 0.0;
    int ok = 0;
    for (int i = 0; i < len; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) len;
}

// Best recovery over `ntry` different planted KEYS (phase-1 IoC is deterministic, so varying
// the solve seed alone cannot rescue a phase-1 miss). This characterises "the attack recovers
// at least one of N random keys at this length/plug count" -- honest for a probabilistic,
// length/plug-limited ciphertext-only method (Gillogly).
static double best_over_keys(const int rotors[3], int len, int nplugs, unsigned base, int ntry) {
    double best = 0.0;
    for (int s = 0; s < ntry; s++) {
        double f = solve_once(rotors, len, nplugs, base + 101u * s, 1u, 0, 0);
        if (f > best) best = f;
        if (best > 0.999) break;
    }
    return best;
}

// --- SearchDefaults registry entry --------------------------------------------------
static void test_registry(void) {
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = ENIGMA; cfg.method = METHOD_DEFAULT;
    CHECK(apply_cipher_defaults(&cfg, false), "enigma registry: no entry applied");
    CHECK(cfg.n_restarts == 2 && cfg.n_hill_climbs == 3000,
          "enigma anneal budget %dx%d, expected 2x3000", cfg.n_restarts, cfg.n_hill_climbs);
    CHECK(cfg.init_temp > 0.0599 && cfg.init_temp < 0.0601, "enigma init_temp %.4f, expected 0.06", cfg.init_temp);

    init_config(&cfg); cfg.cipher_type = ENIGMA; cfg.method = METHOD_SHOTGUN;
    CHECK(apply_cipher_defaults(&cfg, false), "enigma registry (shotgun): no entry");
    CHECK(cfg.n_restarts == 10 && cfg.n_hill_climbs == 3000,
          "enigma shotgun budget %dx%d, expected 10x3000", cfg.n_restarts, cfg.n_hill_climbs);

    init_config(&cfg); cfg.cipher_type = VIGENERE;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
}

// --- Ciphertext-only capability: recovery vs length x plugs (wheel order pinned) ----
// The IoC attack is length-limited (Gillogly): reliable from a few hundred letters, marginal
// below. Recovery is characterised (printed); floors are asserted only where reliable.
static void test_capability(void) {
    const int rotors[3] = { ENIGMA_II, ENIGMA_I, ENIGMA_III };
    struct { int len, plugs; double floor; } cases[] = {
        { 250, 3, 0.85 },   // short: recovers (best-of-keys ~100%)
        { 450, 3, 0.85 },   // longer, few plugs (~90%: usually all but one stecker)
        { 450, 6, 0.80 },   // longer, more plugs (~90%; a near-solution, refinable)
    };
    printf("\nEnigma ciphertext-only recovery (rotors pinned II I III, best of 4 keys):\n");
    for (int c = 0; c < 3; c++) {
        double f = best_over_keys(rotors, cases[c].len, cases[c].plugs, 424242u + 7000u * c, 4);
        printf("  len %-4d plugs %d : %.1f%%  (floor %.0f%%)\n",
               cases[c].len, cases[c].plugs, 100.0 * f, 100.0 * cases[c].floor);
        CHECK(f > cases[c].floor, "len%d %d-plug recovery %.2f, expected > %.2f",
              cases[c].len, cases[c].plugs, f, cases[c].floor);
    }
}

// --- Bombe (crib) recovery: exact rotor/pos + full plugboard from a short crib -------
static void test_bombe(void) {
    const int rotors[3] = { ENIGMA_III, ENIGMA_I, ENIGMA_II };
    printf("\nEnigma Bombe (crib) recovery (rotors pinned, 30-letter crib):\n");
    double f = solve_once(rotors, 250, 6, 20260907u, 0, /*bombe*/ 1, /*crib*/ 30);
    printf("  len 250, 6 plugs, crib 30 : %.1f%%\n", 100.0 * f);
    CHECK(f > 0.90, "Bombe recovery %.2f, expected > 0.90 from a crib", f);
}

int main(void) {
    init_alphabet(NULL);
    enigma_init();
    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("FAIL: could not load %s (run from the source dir)\n", NGRAM_FILE); return 1; }

    test_registry();
    test_capability();
    test_bombe();

    free(shared.ngram_data);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
