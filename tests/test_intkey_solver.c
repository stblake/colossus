//
//  In-process stress / limits tests for the Interrupted Key solver (solve_cipher).
//
//  Framework-free: build with `make testopt`. colossus.c is compiled with -DCOLOSSUS_NO_MAIN
//  and this file supplies its own main, so solve_cipher is driven directly and its SolveResult
//  inspected. A fixed -seed makes each stochastic solve deterministic.
//
//  The Interrupted Key cipher is a periodic base cipher (Vig/Var/Beau) under a keyword whose key
//  index RESETS at break points. The solver sweeps the period and enumerates the interruption
//  STRATEGY: ciphertext-interruptor (ct, breaks known from the cipher -> decouples like Vigenere),
//  plaintext-interruptor (pt, breaks causal -> EM warm start), supplied breaks (random / word-
//  division), and a blind joint keyword+break-mask anneal. This suite -- also the basis for tuning
//  the SearchDefaults schedule -- checks:
//    1. registry validation (apply_cipher_defaults) for all three codes + a non-registry type;
//    2. a CT capability floor (blind interruptor, P pinned) + a length cliff;
//    3. a PT capability floor (blind interruptor, P pinned);
//    4. a BLIND period solve (P swept) -- the reported period must match the true one;
//    5. a BLIND interruptor+scheme solve -- the reported letter AND scheme must match;
//    6. a multi-keyword sweep (mean/worst, CT);
//    7. a SUPPLIED-BREAKS (random / word-division) solve via -breaks (asserted);
//    8. a JOINT blind-random characterization (printed, not asserted -- convergence is text-length
//       dependent, like the CM-Bifid even-period finding), with/without a crib;
//    9. per-scheme calibration: the same cipher under -method anneal / shotgun / pso.
//
//  Interrupted Key runs on the full 26-letter alphabet and -- like the rest of the Vigenere family
//  -- rides the reward-only quadgram table (no -logprob). Run from the source directory.
//

#include "colossus.h"
#include "engine.h"               // apply_cipher_defaults
#include "scoring.h"              // load_ngrams
#include <unistd.h>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4
#define BREAKS_TMP "test_intkey_breaks.tmp"

static SharedData shared;

static const char *PLAINTEXT =
    "WHENINTHECOURSEOFHUMANEVENTSITBECOMESNECESSARYFORONEPEOPLETODISSOLVETHEPOLITICAL"
    "BANDSWHICHHAVECONNECTEDTHEMWITHANOTHERANDTOASSUMEAMONGTHEPOWERSOFTHEEARTHTHESEPARATE"
    "ANDEQUALSTATIONTOWHICHTHELAWSOFNATUREANDOFNATURESGODENTITLETHEMADECENTRESPECTTOTHE"
    "OPINIONSOFMANKINDREQUIRESTHATTHEYSHOULDDECLARETHECAUSESWHICHIMPELTHEMTOTHESEPARATION"
    "WEHOLDTHESETRUTHSTOBESELFEVIDENTTHATALLMENARECREATEDEQUALTHATTHEYAREENDOWEDBYTHEIR";

static const char *BASE_NAME[3] = { "vig", "var", "beau" };
static int base_type(int base) {
    return base == IK_BASE_VAR ? INTERRUPTED_KEY_VAR
         : base == IK_BASE_BEAU ? INTERRUPTED_KEY_BEAU : INTERRUPTED_KEY;
}

static int kw_to_key(const char *kw, int key[]) {
    int P = 0;
    for (int i = 0; kw[i]; i++) key[P++] = toupper((unsigned char) kw[i]) - 'A';
    return P;
}

static int take_plaintext(int prepared[], int pt_len) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < pt_len; i++) {
        int c = PLAINTEXT[i];
        if (c >= 'A' && c <= 'Z') prepared[n++] = c - 'A';
    }
    return n;
}

// Plant a ciphertext-interruptor cipher; returns n, fills prepared[] + cipher_str, *P_out.
static int plant_ct(const char *kw, int pt_len, int base, int interruptor,
                    int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN]; int P = kw_to_key(kw, key); if (P_out) *P_out = P;
    int n = take_plaintext(prepared, pt_len);
    static int cipher[MAX_CIPHER_LENGTH];
    intkey_encrypt_ctint(cipher, prepared, n, key, P, base, interruptor);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

static int plant_pt(const char *kw, int pt_len, int base, int interruptor,
                    int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN]; int P = kw_to_key(kw, key); if (P_out) *P_out = P;
    int n = take_plaintext(prepared, pt_len);
    static int cipher[MAX_CIPHER_LENGTH];
    intkey_encrypt_ptint(cipher, prepared, n, key, P, base, interruptor);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    return n;
}

// Plant a random / word-division cipher and write its break positions to BREAKS_TMP.
static int plant_random(const char *kw, int pt_len, int base, uint32_t rseed,
                        int prepared[], char cipher_str[], int *P_out) {
    int key[MAX_KEYWORD_LEN]; int P = kw_to_key(kw, key); if (P_out) *P_out = P;
    int n = take_plaintext(prepared, pt_len);
    static int is_break[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) is_break[i] = 0;
    seed_rand(rseed);
    int pos = 0, first = 1;
    while (pos < n) { is_break[pos] = 1; int L = first ? P : rand_int(1, P + 1); first = 0; pos += L; }
    static int cipher[MAX_CIPHER_LENGTH];
    intkey_encrypt_mask(cipher, prepared, n, key, P, base, is_break);
    for (int i = 0; i < n; i++) cipher_str[i] = index_to_char(cipher[i]);
    cipher_str[n] = '\0';
    FILE *fp = fopen(BREAKS_TMP, "w");
    if (fp) { for (int i = 0; i < n; i++) if (is_break[i]) fprintf(fp, "%d ", i); fclose(fp); }
    return n;
}

// Solve harness. period>0 pins P; intscheme>=0 pins the strategy; interruptor>=0 pins the letter;
// breaks_file != NULL supplies -breaks; crib != NULL is a cipher-length crib string ('_' unknown).
static double solve_ex(const char *cipher_str, const int prepared[], int plen,
        int cipher_type, int period, int intscheme, int interruptor,
        const char *breaks_file, const char *crib, int method, uint32_t seed,
        int *period_out, int *intr_out, int *scheme_out, double *secs_out) {
    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = cipher_type;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.method = method;
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    if (period > 0)     { cfg.period_present = true; cfg.period = period; }
    if (intscheme >= 0) { cfg.intscheme_present = true; cfg.intscheme = intscheme; }
    if (interruptor >= 0) { cfg.interruptor_present = true; cfg.interruptor = interruptor; }
    if (breaks_file)    { cfg.breaks_present = true; strcpy(cfg.breaks_file, breaks_file); }

    SolveResult res;
    clock_t t0 = clock();
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* still proceed */ }
    seed_rand(seed);
    solve_cipher((char *) cipher_str, (char *) (crib ? crib : ""), &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);
    if (secs_out)   *secs_out = ((double) clock() - t0) / CLOCKS_PER_SEC;
    if (period_out) *period_out = res.solved ? res.cycleword_len : -1;
    if (intr_out)   *intr_out   = res.solved ? res.interruptor : -1;
    if (scheme_out) *scheme_out = res.solved ? res.intscheme : -1;
    if (!res.solved || res.decrypted_len != plen) return 0.0;
    int ok = 0;
    for (int i = 0; i < plen; i++) if (res.decrypted[i] == prepared[i]) ok++;
    return (double) ok / (double) plen;
}

// --- 1. registry validation ---------------------------------------------------

static void test_registry(void) {
    int codes[3] = { INTERRUPTED_KEY, INTERRUPTED_KEY_VAR, INTERRUPTED_KEY_BEAU };
    for (int i = 0; i < 3; i++) {
        ColossusConfig cfg;
        init_config(&cfg); cfg.cipher_type = codes[i]; cfg.method = METHOD_DEFAULT;
        CHECK(apply_cipher_defaults(&cfg, false), "intkey %s registry: no entry applied", BASE_NAME[i]);
        CHECK(cfg.n_restarts == 6 && cfg.n_hill_climbs == 8000,
            "intkey %s anneal defaults wrong: %dx%d", BASE_NAME[i], cfg.n_restarts, cfg.n_hill_climbs);
    }
    ColossusConfig cfg;
    init_config(&cfg); cfg.cipher_type = VIGENERE;
    int r0 = cfg.n_restarts, h0 = cfg.n_hill_climbs;
    CHECK(!apply_cipher_defaults(&cfg, false), "vigenere should have no registry entry");
    CHECK(cfg.n_restarts == r0 && cfg.n_hill_climbs == h0, "non-registry type was modified");
}

// --- 2. CT capability floor (blind interruptor, P pinned) + length cliff --------

static void test_ct_capability_floor(void) {
    int bases[3] = { IK_BASE_VIG, IK_BASE_VAR, IK_BASE_BEAU };
    const char *kws[] = { "GRAPEFRUIT", "KRYPTOS", "PALMERSTON" };
    int intr[]        = { 4,            17,        8 };          // interruptor letters E / R / I
    int plen = 160;
    printf("\n[CT capability floor @ %d chars, P pinned, blind interruptor, per base]\n", plen);
    for (int b = 0; b < 3; b++) {
        for (int k = 0; k < 3; k++) {
            int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
            int n = plant_ct(kws[k], plen, bases[b], intr[k], prepared, cs, &P);
            double secs;
            double frac = solve_ex(cs, prepared, n, base_type(bases[b]),
                P, IK_STRAT_CT, -1, NULL, NULL, METHOD_DEFAULT, 0xC0FFEEu + 7u * b + k,
                NULL, NULL, NULL, &secs);
            printf("  %-4s %-11s P=%2d intr=%c : %.1f%%  [%.1fs]\n",
                BASE_NAME[b], kws[k], P, 'A' + intr[k], 100.0 * frac, secs);
            CHECK(frac > 0.95, "intkey CT %s/%s capability floor: only %.1f%% at %d chars",
                BASE_NAME[b], kws[k], 100.0 * frac, n);
        }
    }
}

static void test_ct_length_cliff(void) {
    int lens[] = { 40, 60, 90, 130, 180 };
    const char *kw = "GRAPEFRUIT"; int intr = 4;
    printf("\n[CT length cliff: keyword=%s intr=E (Vigenere), P pinned, blind interruptor]\n", kw);
    double best = 0.0;
    for (int li = 0; li < (int) (sizeof lens / sizeof lens[0]); li++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant_ct(kw, lens[li], IK_BASE_VIG, intr, prepared, cs, &P);
        double secs;
        double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_CT, -1, NULL, NULL,
            METHOD_DEFAULT, 0x5EEDu + li, NULL, NULL, NULL, &secs);
        printf("  %3d chars : %.1f%%  [%.1fs]\n", n, 100.0 * frac, secs);
        if (frac > best) best = frac;
    }
    CHECK(best > 0.95, "CT length cliff: never recovered (best %.1f%%)", 100.0 * best);
}

// --- 3. PT capability floor + fragility characterization ------------------------
//
// The plaintext-interruptor is the weaker (blind) mode: because the reset trigger is a plaintext
// letter, a wrong (keyword, interruptor) pair can still de-cohere into fairly English text, and a
// very common interruptor relative to P starves the high keyword columns. So we ASSERT the PT
// machinery on reliable pinned-interruptor cases (proving the EM warm-start + causal decrypt +
// anneal recover the keyword when the interruptor is known), and CHARACTERIZE the blind fragility
// (printed): blind PT is reliable on longer text but can mispick the interruptor on short text,
// and some interruptor letters are degenerate even when pinned.

static void test_pt_capability_floor(void) {
    // Reliable pinned-interruptor cases (interruptor known): assert the PT machinery recovers --
    // the EM warm-start + causal decrypt + anneal do work when the interruptor is not guessed.
    struct { const char *kw; int intr; int plen; uint32_t seed; } ok[] = {
        { "KRYPTOS",    19, 220, 0x9A1Du },   // T
        { "GRAPEFRUIT", 17, 220, 0x9A1Eu },   // R
    };
    printf("\n[PT capability floor, interruptor PINNED (known), P pinned]\n");
    for (int k = 0; k < 2; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int n = plant_pt(ok[k].kw, ok[k].plen, IK_BASE_VIG, ok[k].intr, prepared, cs, &P);
        double secs;
        double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_PT, ok[k].intr, NULL, NULL,
            METHOD_DEFAULT, ok[k].seed, NULL, NULL, NULL, &secs);
        printf("  %-11s P=%2d intr=%c : %.1f%%  [%.1fs]\n", ok[k].kw, P, 'A' + ok[k].intr, 100.0 * frac, secs);
        CHECK(frac > 0.90, "intkey PT %s/intr=%c pinned: only %.1f%% at %d chars",
            ok[k].kw, 'A' + ok[k].intr, 100.0 * frac, n);
    }

    // Characterization (not asserted): PT is the FRAGILE mode. Because the reset trigger is a
    // plaintext letter, a wrong (keyword, interruptor) pair can still de-cohere into fairly English
    // text, so blind PT can mispick the interruptor and even a pinned interruptor can hit a rugged
    // multi-modal basin (seed/content sensitive). CT is the reliable blind workhorse; PT is a
    // best-effort second lens. The sweep below documents the variance.
    printf("  [PT fragility, GRAPEFRUIT @ 220, interruptor swept -- characterization, not asserted]\n");
    int sweep[] = { 4, 13, 17, 18, 7 };   // E N R S H
    for (int i = 0; i < 5; i++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P; int iout;
        int n = plant_pt("GRAPEFRUIT", 220, IK_BASE_VIG, sweep[i], prepared, cs, &P);
        double blind = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_PT, -1, NULL, NULL,
            METHOD_DEFAULT, 0x9A1Du, NULL, &iout, NULL, NULL);
        double pinned = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_PT, sweep[i], NULL, NULL,
            METHOD_DEFAULT, 0x9A1Du, NULL, NULL, NULL, NULL);
        printf("    intr=%c : blind %.1f%% (reported %c)  pinned %.1f%%\n", 'A' + sweep[i],
            100.0 * blind, iout >= 0 ? 'A' + iout : '?', 100.0 * pinned);
    }
}

// --- 4. blind period solve (P swept, CT) ----------------------------------------

static void test_blind_period(void) {
    const char *kw = "KRYPTOS"; int plen = 200, intr = 4;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant_ct(kw, plen, IK_BASE_VIG, intr, prepared, cs, &P);
    double secs; int pout;
    double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, 0, IK_STRAT_CT, -1, NULL, NULL,
        METHOD_DEFAULT, 0xB11Du, &pout, NULL, NULL, &secs);
    printf("\n[blind P (CT), true keyword=%s P=%d intr=E]: reported P=%d, %.1f%%  [%.1fs]\n",
        kw, P, pout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-P recovery only %.1f%%", 100.0 * frac);
    CHECK(pout == P, "blind-P reported P=%d (true %d)", pout, P);
}

// --- 5. blind interruptor + scheme solve (fully blind) --------------------------

static void test_blind_interruptor(void) {
    const char *kw = "GRAPEFRUIT"; int plen = 200, intr = 8;      // interruptor 'I'
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant_ct(kw, plen, IK_BASE_VIG, intr, prepared, cs, &P);
    double secs; int iout, sout;
    // Fully blind: P pinned (period estimation is out of scope), interruptor + scheme enumerated.
    double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, -1, -1, NULL, NULL,
        METHOD_DEFAULT, 0xB1FFu, NULL, &iout, &sout, &secs);
    printf("\n[blind interruptor+scheme, true intr=%c scheme=ct]: reported intr=%c scheme=%d, %.1f%%  [%.1fs]\n",
        'A' + intr, iout >= 0 ? 'A' + iout : '?', sout, 100.0 * frac, secs);
    CHECK(frac > 0.95, "blind-interruptor recovery only %.1f%%", 100.0 * frac);
    CHECK(iout == intr, "blind-interruptor reported %d (true %d)", iout, intr);
    CHECK(sout == IK_STRAT_CT, "blind-scheme reported %d (true CT=%d)", sout, IK_STRAT_CT);
}

// --- 6. multi-keyword sweep (CT) ------------------------------------------------

static void test_multi_keyword(void) {
    const char *kws[] = { "ZEBRA", "CIPHER", "KRYPTOS", "MONARCHY", "GRAPEFRUIT", "PALMERSTON" };
    int plen = 160, nk = (int) (sizeof kws / sizeof kws[0]);
    double sum = 0, worst = 1.0;
    printf("\n[multi-keyword sweep @ %d chars (CT, Vigenere), P pinned, blind interruptor]\n", plen);
    for (int k = 0; k < nk; k++) {
        int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
        int intr = 4 + (k % 5);                                  // vary the interruptor letter
        int n = plant_ct(kws[k], plen, IK_BASE_VIG, intr, prepared, cs, &P);
        double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_CT, -1, NULL, NULL,
            METHOD_DEFAULT, 0xABCDu + k, NULL, NULL, NULL, NULL);
        printf("  %-11s P=%2d intr=%c : %.1f%%\n", kws[k], P, 'A' + intr, 100.0 * frac);
        sum += frac; if (frac < worst) worst = frac;
    }
    printf("  mean=%.1f%%  worst=%.1f%%\n", 100.0 * sum / nk, 100.0 * worst);
    CHECK(sum / nk > 0.90, "multi-keyword mean too low: %.1f%%", 100.0 * sum / nk);
}

// --- 7. supplied-breaks (random / word-division) solve --------------------------

static void test_supplied_breaks(void) {
    const char *kw = "KRYPTOS"; int plen = 200;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant_random(kw, plen, IK_BASE_VIG, 0x1234u, prepared, cs, &P);
    double secs;
    double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_BREAKS, -1, BREAKS_TMP, NULL,
        METHOD_DEFAULT, 0x0B00u, NULL, NULL, NULL, &secs);
    printf("\n[supplied-breaks (random), keyword=%s P=%d]: %.1f%%  [%.1fs]\n", kw, P, 100.0 * frac, secs);
    CHECK(frac > 0.95, "supplied-breaks recovery only %.1f%%", 100.0 * frac);
}

// --- 8. JOINT blind-random characterization (printed, not asserted) -------------

static void test_joint_characterize(void) {
    const char *kw = "KRYPTOS"; int plen = 300;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant_random(kw, plen, IK_BASE_VIG, 0x1234u, prepared, cs, &P);

    // Build a crib revealing the first 20 plaintext letters (positional: decrypted[i]==pt[i]).
    static char crib[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < n; i++) crib[i] = (i < 20) ? index_to_char(prepared[i]) : '_';
    crib[n] = '\0';

    double s1, s2;
    double blind = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_JOINT, -1, NULL, NULL,
        METHOD_DEFAULT, 0x30117u, NULL, NULL, NULL, &s1);
    double cribbed = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_JOINT, -1, NULL, crib,
        METHOD_DEFAULT, 0x30117u, NULL, NULL, NULL, &s2);
    printf("\n[JOINT blind-random characterization, keyword=%s P=%d @ %d chars] (not asserted)\n"
           "  blind        : %.1f%%  [%.1fs]\n"
           "  +20-char crib : %.1f%%  [%.1fs]\n",
        kw, P, n, 100.0 * blind, s1, 100.0 * cribbed, s2);
}

// --- 9. per-scheme calibration (anneal / shotgun / pso) -------------------------

static void test_per_scheme(void) {
    const char *kw = "KRYPTOS"; int plen = 160, intr = 4;
    int prepared[MAX_CIPHER_LENGTH]; char cs[MAX_CIPHER_LENGTH]; int P;
    int n = plant_ct(kw, plen, IK_BASE_VIG, intr, prepared, cs, &P);
    struct { int method; const char *name; } M[] = {
        { METHOD_DEFAULT, "anneal " }, { METHOD_SHOTGUN, "shotgun" }, { METHOD_PSO, "pso    " },
    };
    printf("\n[per-scheme @ %d chars (CT), keyword=%s intr=E, P pinned]\n", plen, kw);
    for (int m = 0; m < 3; m++) {
        double secs;
        double frac = solve_ex(cs, prepared, n, INTERRUPTED_KEY, P, IK_STRAT_CT, -1, NULL, NULL,
            M[m].method, 0x5C8E0u + m, NULL, NULL, NULL, &secs);
        printf("  %s : %.1f%%  [%.1fs]\n", M[m].name, 100.0 * frac, secs);
        if (M[m].method == METHOD_DEFAULT)
            CHECK(frac > 0.95, "default (anneal) scheme recovery only %.1f%%", 100.0 * frac);
    }
}

int main(void) {
    init_alphabet(NULL);                          // full 26-letter alphabet
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("cannot load %s\n", NGRAM_FILE); return 1; }

    test_registry();
    test_ct_capability_floor();
    test_ct_length_cliff();
    test_pt_capability_floor();
    test_blind_period();
    test_blind_interruptor();
    test_multi_keyword();
    test_supplied_breaks();
    test_joint_characterize();
    test_per_scheme();

    remove(BREAKS_TMP);
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
