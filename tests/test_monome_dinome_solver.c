//
//  Solver calibration / capability tests for the Monome-Dinome cipher (type 78). Built with
//  -DCOLOSSUS_NO_MAIN so solve_cipher() links directly; fixed seeds -> deterministic. Covers:
//  registry validation, the STRUCTURAL config-discrimination invariant (the true indicator
//  pair is always exactly 100% valid, so validity retains it), a letter capability floor
//  (fixed AND keyed boards), and a length cliff. Monome-Dinome effectively NEEDS -logprob with
//  QUINTGRAMS and a DICTIONARY (the cross-config selector is dictionary coverage -- see
//  monome_dinome_solver.h), so this test loads both, mirroring a real -type md run.
//
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "monome_dinome.h"
#include "monome_dinome_solver.h"
#include <fcntl.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE (char *) "english_quintgrams.txt"
#define NGRAM_SIZE 5
#define DICT_FILE  (char *) "OxfordEnglishWords.txt"

static SharedData shared;

// ~260 letters of English (Pride & Prejudice opening); J/Z fold to I/Y in the 24-letter set.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEIN"
    "WANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGA"
    "NEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONS"
    "IDEREDTHERIGHTFULPROPERTYOFSOMEONEOFTHEIRDAUGHTERS";

// Build a Monome-Dinome board from a keyword (24-letter keyed alphabet), the 8 column-label
// digits, and the two indicators. Uses g_char_to_idx (J->I, Z->Y already folded).
static void build_board_kw(MonomeDinomeBoard *b, const char *kw, const int col_label[8], int r0, int r1) {
    int letters[MD_NALPHA], used[MD_NALPHA] = {0}, m = 0;
    for (const char *k = kw; *k; k++) {
        int c = g_char_to_idx[(int) *k];
        if (c >= 0 && c < MD_NALPHA && !used[c]) { used[c] = 1; letters[m++] = c; }
    }
    for (int c = 0; c < MD_NALPHA; c++) if (!used[c]) letters[m++] = c;
    monome_dinome_build_board(b, letters, col_label, r0, r1);
}

// Structural token validity of config {a,b} (a copy of the solver's, for the invariant test).
static double config_validity(const int *dg, int n, int a, int b) {
    int i = 0, nt = 0, nv = 0;
    while (i < n) {
        int g = dg[i];
        if (g == a || g == b) { if (i + 1 >= n) { nt++; i++; continue; } int g2 = dg[i + 1]; i += 2;
            if (g2 == a || g2 == b) { nt++; continue; } nt++; nv++; }
        else { i++; nt++; nv++; }
    }
    return (nt > 0) ? (double) nv / nt : 0.0;
}

// Plant a cipher from plaintext letters and solve it in-process; return the letter match count.
static int plant_and_solve(const int *pt, int n, const char *kw, const int col_label[8],
                           int r0, int r1, uint32_t seed, int *matched) {
    MonomeDinomeBoard board;
    build_board_kw(&board, kw, col_label, r0, r1);
    static int cipher[2 * MAX_CIPHER_LENGTH];
    int clen = monome_dinome_encrypt(pt, n, &board, cipher);
    static char cs[2 * MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < clen; i++) cs[i] = '0' + cipher[i];
    cs[clen] = '\0';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = MONOME_DINOME;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.dictionary_present = (shared.dict != NULL);
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);

    SolveResult res;
    res.solved = false; res.decrypted_len = 0;
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* proceed anyway */ }
    seed_rand(seed);
    solve_cipher(cs, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);

    int ok = 0, lim = (res.decrypted_len < n) ? res.decrypted_len : n;
    for (int i = 0; i < lim; i++) if (res.decrypted[i] == pt[i]) ok++;
    if (matched) *matched = ok;
    return n;
}

// Fill pt[] with the first L letters of PLAINTEXT mapped into the 24-letter alphabet.
static int letters_pt(int L, int pt[]) {
    for (int i = 0; i < L; i++) pt[i] = g_char_to_idx[(int) PLAINTEXT[i]];
    return L;
}

int main(void) {
    init_alphabet_monome_dinome();              // 24-letter A..Z with J->I, Z->Y, as -type md forces
    CHECK(g_alpha == MD_NALPHA, "alphabet %d, expected %d", g_alpha, MD_NALPHA);
    g_ngram_logprob = true;                     // Monome-Dinome effectively needs -logprob (quintgrams)
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }
    load_dictionary(DICT_FILE, &shared.dict, &shared.n_dict_words, &shared.max_dict_word_len, false);
    CHECK(shared.dict != NULL, "could not load %s (the coverage selector needs it)", DICT_FILE);

    static int pt[MAX_CIPHER_LENGTH];
    int col_fixed[8] = {1, 8, 9, 2, 7, 0, 5, 4};   // the ACA example column labels (indicators 6,3)
    int col_keyed[8] = {8, 1, 4, 0, 2, 7, 9, 5};   // a scrambled label order

    // 1) Registry validation: the tuned schedule is applied.
    {
        ColossusConfig cfg; init_config(&cfg);
        cfg.cipher_type = MONOME_DINOME;
        bool applied = apply_cipher_defaults(&cfg, false);
        CHECK(applied, "no SearchDefaults entry for MONOME_DINOME");
        CHECK(cfg.init_temp > 0.2 && cfg.init_temp < 0.4, "unexpected inittemp %.2f", cfg.init_temp);
        CHECK(cfg.n_restarts >= 12, "unexpected a_n_restarts %d (want the restart-heavy profile)", cfg.n_restarts);
        ColossusConfig v; init_config(&v); v.cipher_type = VIGENERE;
        CHECK(!apply_cipher_defaults(&v, false), "VIGENERE unexpectedly has a registry entry");
    }

    // 2) Config-discrimination INVARIANT: the true indicator pair is always exactly 100% valid
    //    (a structural guarantee), so the validity pre-filter can never drop it. Verified over
    //    several lengths; also reports how many pairs tie (the residual n-gram/coverage job).
    {
        int allok = 1;
        for (int L = 80; L <= 200; L += 40) {
            int n = letters_pt(L, pt);
            MonomeDinomeBoard board; build_board_kw(&board, "SECRETKEYWORD", col_fixed, 6, 3);
            static int cipher[2 * MAX_CIPHER_LENGTH];
            int clen = monome_dinome_encrypt(pt, n, &board, cipher);
            double vtrue = config_validity(cipher, clen, 3, 6);
            double maxv = 0; int nties = 0;
            for (int a = 0; a < 10; a++) for (int b = a + 1; b < 10; b++) {
                double v = config_validity(cipher, clen, a, b);
                if (v > maxv) maxv = v;
            }
            for (int a = 0; a < 10; a++) for (int b = a + 1; b < 10; b++)
                if (config_validity(cipher, clen, a, b) >= maxv - 1e-9) nties++;
            printf("  L=%3d: true-config validity=%.4f, max=%.4f, fully-valid pairs=%d\n",
                   L, vtrue, maxv, nties);
            if (vtrue < 0.99999 || vtrue < maxv - 1e-9) allok = 0;
        }
        CHECK(allok, "true indicator pair is not always the (tied) validity maximum");
    }

    // 3) LETTER capability floor at 140 chars -- fixed labels (2 keywords) and keyed labels.
    {
        int n = letters_pt(140, pt);
        struct { const char *kw; const int *lab; int r0, r1; } cases[] = {
            {"SECRET", col_fixed, 6, 3}, {"CIPHERKEY", col_fixed, 6, 3}, {"WARTHOG", col_keyed, 4, 1},
        };
        int allok = 1;
        for (int c = 0; c < 3; c++) {
            int ok; plant_and_solve(pt, n, cases[c].kw, cases[c].lab, cases[c].r0, cases[c].r1, 1u + c, &ok);
            double f = (double) ok / n;
            printf("  floor L=140 kw=%s: %.3f (%d/%d)\n", cases[c].kw, f, ok, n);
            if (f < 0.95) allok = 0;
        }
        CHECK(allok, "letter capability floor below 0.95 at L=140");
    }

    // 4) Length cliff (keyword SECRETKEYWORD, keyed labels, indicators 6,3).
    printf("\n  Monome-Dinome letter recovery vs length (keyword SECRETKEYWORD, keyed labels):\n");
    {
        int lens[] = {100, 120, 140, 160};
        double f120 = 0, f140 = 0, f160 = 0;
        for (int li = 0; li < 4; li++) {
            int n = letters_pt(lens[li], pt);
            int ok; plant_and_solve(pt, n, "SECRETKEYWORD", col_keyed, 6, 3, 7u, &ok);
            double f = (double) ok / n;
            printf("    L=%4d: %.3f (%d/%d)\n", lens[li], f, ok, n);
            if (lens[li] == 120) f120 = f;
            if (lens[li] == 140) f140 = f;
            if (lens[li] == 160) f160 = f;
        }
        CHECK(f140 > 0.95, "length cliff: L=140 recovery %.3f (want > 0.95)", f140);
        CHECK(f160 > 0.95, "length cliff: L=160 recovery %.3f (want > 0.95)", f160);
        CHECK(f120 > 0.90, "length cliff: L=120 recovery %.3f (want > 0.90)", f120);
        printf("    (the ACA 60-100 low end is below the blind floor -- a documented limitation)\n");
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(shared.ngram_data);
    if (shared.dict) free_dictionary(shared.dict, shared.n_dict_words);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
