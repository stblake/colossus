// In-process solver regression for the Enigma type. Built with -DCOLOSSUS_NO_MAIN and linked
// against the whole solver so it calls solve_cipher() directly. Validates the SearchDefaults
// registry entry, characterises the ciphertext-only capability (recovery vs length/plugs with
// the wheel order pinned -- the realistic fast path), checks the Bombe recovers a planted
// crib exactly, and checks the position-free crib DRAG (-bombe -cribdrag) lands a stop on a
// true crib, none on a false one, and honours -cribdragmaxoffset. Run from the source
// directory so the n-gram table is found in the cwd.

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
// adaptive => -enigmaadaptive (the Ostwald-Weierud completed-n-gram config ranking).
static double solve_once(const int rotors[3], int len, int nplugs,
                         unsigned plant_seed, unsigned solve_seed, int bombe, int crib_len,
                         int adaptive) {
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
    cfg.enigma_adaptive = adaptive;
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
        double f = solve_once(rotors, len, nplugs, base + 101u * s, 1u, 0, 0, /*adaptive*/ 0);
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
    double f = solve_once(rotors, 250, 6, 20260907u, 0, /*bombe*/ 1, /*crib*/ 30, /*adaptive*/ 0);
    printf("  len 250, 6 plugs, crib 30 : %.1f%%\n", 100.0 * f);
    CHECK(f > 0.90, "Bombe recovery %.2f, expected > 0.90 from a crib", f);
}

// --- Bombe CRIB DRAG: a position-free crib recovers the config + plugboard -----------
// One in-process drag: plant with the wheel order AND ring pinned (so the stop decrypts
// exactly and the run stays quick -- the fast-ring sweep is already covered by test_bombe),
// then call solve_enigma_bombe_drag() with a single dragged `crib_word`. Returns the
// plaintext recovery fraction and sets *found to the drag's "a stop was completed" flag.
static double drag_run(const int rotors[3], int len, int nplugs, unsigned seed,
                       const char *crib_word, int maxoffset, int *found) {
    EnigmaKey k;
    make_key(&k, rotors, nplugs, seed);
    for (int i = 0; i < 3; i++) k.ring[i] = 0;        // ring AAA, pinned below -> exact decrypt
    int prepared[600];
    char cipher_str[601];
    plant(&k, len, prepared, cipher_str);
    int cipher_idx[600];
    for (int i = 0; i < len; i++) cipher_idx[i] = cipher_str[i] - 'A';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = ENIGMA;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = METHOD_DEFAULT;
    cfg.enigma_rotors_present = true;
    for (int i = 0; i < 3; i++) cfg.enigma_rotors[i] = rotors[i];
    cfg.enigma_ring_present = true;
    for (int i = 0; i < 3; i++) cfg.enigma_ring[i] = 0;
    cfg.enigma_bombe = true;
    cfg.cribdrag_present = true;
    cfg.cribdrag.nwords = 1;
    int L = (int) strlen(crib_word);
    cfg.cribdrag.wordlen[0] = L;
    for (int i = 0; i < L; i++) cfg.cribdrag.words[0][i] = crib_word[i] - 'A';
    cfg.cribdrag_max_offset = maxoffset;
    cfg.n_threads = 4;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);

    SolveResult res;
    res.solved = false;
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* proceed anyway */ }
    seed_rand(1u);
    bool f = solve_enigma_bombe_drag(&cfg, &shared, cipher_idx, len, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout)); close(saved); clearerr(stdout);

    if (found) *found = f ? 1 : 0;
    if (!res.solved || res.decrypted_len != len) return 0.0;
    int ok = 0;
    for (int i = 0; i < len; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) len;
}

static void test_bombe_drag(void) {
    const int rotors[3] = { ENIGMA_IV, ENIGMA_II, ENIGMA_V };
    printf("\nEnigma Bombe CRIB DRAG (position-free crib, rotors+ring pinned):\n");

    // 1) A long TRUE crib (the first 30 plaintext letters, so its true offset is 0) lands a
    //    stop that completes to the whole key and decrypts the message.
    char ctrue[64]; memcpy(ctrue, PLAINTEXT, 30); ctrue[30] = '\0';
    int f1 = 0;
    double r1 = drag_run(rotors, 150, 3, 20260922u, ctrue, -1, &f1);
    printf("  true  crib (30, off 0,  cap -1): found=%d  %.1f%%\n", f1, 100.0 * r1);
    CHECK(f1 && r1 > 0.90, "drag true crib: expected a stop recovering > 0.90 (found=%d, %.2f)", f1, r1);

    // 2) A long FALSE crib (not in the plaintext) lands NO stop anywhere.
    int f2 = 0;
    double r2 = drag_run(rotors, 150, 3, 20260922u, "THEQUICKBROWNFOXJUMPSOVERLAZYD", -1, &f2);
    printf("  false crib (30,          cap -1): found=%d  %.1f%%\n", f2, 100.0 * r2);
    CHECK(!f2 && r2 < 0.5, "drag false crib: expected no stop (found=%d, %.2f)", f2, r2);

    // 3/4) -cribdragmaxoffset gates a crib whose only true offset is 40: cap 40 keeps it,
    //      cap 39 excludes it (so no stop).
    char cmid[64]; memcpy(cmid, PLAINTEXT + 40, 26); cmid[26] = '\0';
    int f3 = 0, f4 = 0;
    double r3 = drag_run(rotors, 150, 3, 20260922u, cmid, 40, &f3);
    double r4 = drag_run(rotors, 150, 3, 20260922u, cmid, 39, &f4);
    printf("  mid   crib (26, off 40, cap 40): found=%d  %.1f%%\n", f3, 100.0 * r3);
    printf("  mid   crib (26, off 40, cap 39): found=%d  %.1f%%\n", f4, 100.0 * r4);
    CHECK(f3 && r3 > 0.90, "drag mid crib (cap includes true offset): expected > 0.90 (found=%d, %.2f)", f3, r3);
    CHECK(!f4, "drag mid crib (cap 39 excludes true offset 40): expected no stop (found=%d)", f4);
}

// --- -enigmaadaptive: short-message config ranking (Ostwald-Weierud) -----------------
// Colossus's default pipeline ranks rotor configs by EMPTY-plugboard IoC, which drops the true
// config on short/many-plug messages before any plugboard climb (the real short-message floor --
// selection, not plugboard-climb quality). -enigmaadaptive reranks the top configs by a plugboard-
// COMPLETED n-gram (Tier 1a), plus a forced-E-Stecker partial exhaustion below ~300 letters (Tier
// 1b), which rescues the true config. Recovery is probabilistic per key (Gillogly), so the gain is
// a MEAN over a fixed key set, not a per-key guarantee. Two regimes:
//   130 letters / 3 plugs -- the tractable case: adaptive reliably lifts mean recovery well clear.
//   150 letters / 6 plugs -- HARD: adaptive lifts mean recovery several-fold (the E-Stecker win)
//                            but full solves stay rare; the >=6-plug short floor is near the
//                            fundamental limit (documented, matches the literature).
static void test_adaptive(void) {
    const int rotors[3] = { ENIGMA_II, ENIGMA_I, ENIGMA_III };
    const int nk = 10;
    double moff3 = 0, mon3 = 0, moff6 = 0, mon6 = 0;
    for (int s = 0; s < nk; s++) {
        unsigned ks = 500000u + 1009u * s;
        moff3 += solve_once(rotors, 130, 3, ks, 1u, 0, 0, /*adaptive*/ 0);
        mon3  += solve_once(rotors, 130, 3, ks, 1u, 0, 0, /*adaptive*/ 1);
        moff6 += solve_once(rotors, 150, 6, ks, 1u, 0, 0, /*adaptive*/ 0);
        mon6  += solve_once(rotors, 150, 6, ks, 1u, 0, 0, /*adaptive*/ 1);
    }
    moff3 /= nk; mon3 /= nk; moff6 /= nk; mon6 /= nk;
    printf("\nEnigma -enigmaadaptive short-message gain (mean over %d keys):\n", nk);
    printf("  len 130 / 3 plugs : OFF %.1f%%   ON %.1f%%\n", 100.0 * moff3, 100.0 * mon3);
    printf("  len 150 / 6 plugs : OFF %.1f%%   ON %.1f%%  (hard regime)\n", 100.0 * moff6, 100.0 * mon6);
    // Tier 1a: the tractable few-plug case improves substantially and clears a real floor.
    CHECK(mon3 > moff3, "adaptive mean %.2f should beat default %.2f at 130/3", mon3, moff3);
    CHECK(mon3 > 0.45, "adaptive mean %.2f should clear the short-message floor at 130/3", mon3);
    // Tier 1b: the E-Stecker exhaustion lifts the hard many-plug case (mean, not necessarily solves).
    CHECK(mon6 >= moff6, "adaptive mean %.2f should not underperform default %.2f at 150/6", mon6, moff6);
    CHECK(mon6 > moff6, "adaptive (E-Stecker) mean %.2f should beat default %.2f at 150/6", mon6, moff6);
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
    test_bombe_drag();
    test_adaptive();

    free(shared.ngram_data);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
