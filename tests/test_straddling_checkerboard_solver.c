//
//  Solver calibration / capability tests for the Straddling Checkerboard (type 76).
//  Built with -DCOLOSSUS_NO_MAIN so solve_cipher() is linked directly. Fixed seeds ->
//  deterministic. Covers: registry validation, a LETTER capability floor (fixed AND keyed
//  labels), a length cliff, and a FIGURE-SHIFT (numeric) characterization -- letters recover,
//  the digit VALUES are a documented limitation (weak FIG-placement gradient), so digit
//  recovery is printed, not asserted.
//
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "straddling_checkerboard.h"
#include "straddling_checkerboard_solver.h"
#include <fcntl.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE (char *) "english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

// ~260 letters of English (Pride & Prejudice opening), the plaintext source.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEIN"
    "WANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGA"
    "NEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONS"
    "IDEREDTHERIGHTFULPROPERTYOFSOMEONEOFTHEIRDAUGHTERS";

// Build a board from an arrangement keyword (keyed alphabet + FIG), optional keyed labels
// (NULL => identity 0..9), and the two blank columns.
static void build_board_kw(StraddlingBoard *b, const char *kw, const int *labels_in, int b0, int b1) {
    int seq[SC_NSYM], used[26] = {0}, m = 0;
    for (const char *k = kw; *k; k++) {
        int c = *k - 'A';
        if (c >= 0 && c < 26 && !used[c]) { used[c] = 1; seq[m++] = c; }
    }
    for (int c = 0; c < 26; c++) if (!used[c]) seq[m++] = c;
    seq[26] = SC_FIG;
    int labels[10];
    if (labels_in) for (int c = 0; c < 10; c++) labels[c] = labels_in[c];
    else           for (int c = 0; c < 10; c++) labels[c] = c;
    straddling_build_board(b, seq, labels, labels[b0], labels[b1]);
}

// Plant a cipher from plaintext syms and solve it in-process; return the letter- and
// digit-position match counts (via the out-params). Stdout is muted during the solve.
static void plant_and_solve(const int *pt, int n, const char *kw, const int *labels,
                            int b0, int b1, uint32_t seed,
                            int *letters_ok, int *letters_tot, int *digits_ok, int *digits_tot) {
    StraddlingBoard board;
    build_board_kw(&board, kw, labels, b0, b1);
    static int cipher[3 * MAX_CIPHER_LENGTH];
    int clen = straddling_encrypt(pt, n, &board, cipher);
    static char cs[3 * MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < clen; i++) cs[i] = '0' + cipher[i];
    cs[clen] = '\0';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = STRADDLING_CHECKERBOARD;
    cfg.ngram_size = NGRAM_SIZE;
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

    int lo = 0, lt = 0, dok = 0, dt = 0;
    int lim = (res.decrypted_len < n) ? res.decrypted_len : n;
    for (int i = 0; i < n; i++) {
        if (pt[i] < 26) { lt++; if (i < lim && res.decrypted[i] == pt[i]) lo++; }
        else            { dt++; if (i < lim && res.decrypted[i] == pt[i]) dok++; }
    }
    if (letters_ok) *letters_ok = lo;
    if (letters_tot) *letters_tot = lt;
    if (digits_ok) *digits_ok = dok;
    if (digits_tot) *digits_tot = dt;
}

// Fill pt[] with the first L letters of PLAINTEXT (all letters, 0..25). Returns L.
static int letters_pt(int L, int pt[]) {
    for (int i = 0; i < L; i++) pt[i] = PLAINTEXT[i] - 'A';
    return L;
}

int main(void) {
    init_alphabet_adfgvx();                     // 36-symbol A..Z + 0..9, as -type sc forces
    CHECK(g_alpha == 36, "alphabet %d, expected 36", g_alpha);
    g_ngram_logprob = true;                     // straddling effectively needs -logprob
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    shared.dict = NULL; shared.n_dict_words = 0; shared.max_dict_word_len = 0;
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    static int pt[MAX_CIPHER_LENGTH];

    // 1) Registry validation: the tuned schedule is applied.
    {
        ColossusConfig cfg; init_config(&cfg);
        cfg.cipher_type = STRADDLING_CHECKERBOARD;
        bool applied = apply_cipher_defaults(&cfg, false);
        CHECK(applied, "no SearchDefaults entry for STRADDLING_CHECKERBOARD");
        CHECK(cfg.init_temp > 0.2 && cfg.init_temp < 0.4, "unexpected inittemp %.2f", cfg.init_temp);
        // A non-registered type must be left untouched.
        ColossusConfig v; init_config(&v); v.cipher_type = VIGENERE;
        CHECK(!apply_cipher_defaults(&v, false), "VIGENERE unexpectedly has a registry entry");
    }

    // 2) LETTER capability floor at 220 chars -- fixed 0..9 labels, two keywords.
    {
        int n = letters_pt(220, pt);
        const char *kws[] = {"SECRET", "CIPHERKEY"};
        int allok = 1;
        for (int k = 0; k < 2; k++) {
            int lo, lt;
            plant_and_solve(pt, n, kws[k], NULL, 2, 6, 1u + k, &lo, &lt, NULL, NULL);
            double f = (double) lo / lt;
            printf("  letter floor L=220 kw=%s: %.3f (%d/%d)\n", kws[k], f, lo, lt);
            if (f < 0.99) allok = 0;
        }
        CHECK(allok, "letter capability floor (fixed labels) below 0.99 at L=220");
    }

    // 3) KEYED-LABEL capability floor at 220 chars (a scrambled column-heading permutation).
    {
        int n = letters_pt(220, pt);
        int labels[10] = {7, 2, 9, 0, 5, 1, 8, 3, 6, 4};
        int lo, lt;
        plant_and_solve(pt, n, "WARTHOG", labels, 0, 6, 3u, &lo, &lt, NULL, NULL);
        double f = (double) lo / lt;
        printf("  keyed-label floor L=220: %.3f (%d/%d)\n", f, lo, lt);
        CHECK(f > 0.99, "keyed-label capability floor below 0.99 at L=220 (%.3f)", f);
    }

    // 4) Length cliff (fixed labels, keyword SECRET, blanks 2&6).
    printf("\n  Straddling letter recovery vs length (keyword SECRET, fixed labels):\n");
    {
        int lens[] = {100, 140, 180, 220};
        double f180 = 0, f220 = 0;
        for (int li = 0; li < 4; li++) {
            int n = letters_pt(lens[li], pt);
            int lo, lt;
            plant_and_solve(pt, n, "SECRET", NULL, 2, 6, 7u, &lo, &lt, NULL, NULL);
            double f = (double) lo / lt;
            printf("    L=%4d: %.3f (%d/%d)\n", lens[li], f, lo, lt);
            if (lens[li] == 180) f180 = f;
            if (lens[li] == 220) f220 = f;
        }
        CHECK(f220 > 0.99, "length cliff: L=220 recovery %.3f (want > 0.99)", f220);
        CHECK(f180 > 0.90, "length cliff: L=180 recovery %.3f (want > 0.90)", f180);
    }

    // 5) FIGURE-SHIFT (numeric) characterization -- a DOCUMENTED LIMITATION, printed, NOT
    //    asserted. Figure-mode digits carry ~0 n-gram weight so FIG placement has too weak a
    //    gradient to recover blind; worse, an unrecognized FIG code decodes as an extra letter,
    //    SHIFTING alignment, so mid-stream digit runs corrupt the following letters. The
    //    primitive + solver fully support figure-shift (a known/recovered board decodes the
    //    digits), but blind numeric-bearing solves are the limitation, characterized here.
    printf("\n  Figure-shift (numeric) characterization -- DOCUMENTED LIMITATION (not asserted):\n");
    {
        // (a) trailing digit run: letters ahead of the run still recover well.
        int n = letters_pt(200, pt);
        int run[] = {1, 8, 1, 3}; for (int i = 0; i < 4; i++) pt[n++] = 26 + run[i];
        int lo, lt, dok, dt;
        plant_and_solve(pt, n, "SECRET", NULL, 2, 6, 1u, &lo, &lt, &dok, &dt);
        printf("    trailing digits  (L=%d): letters %.3f (%d/%d), digit values %d/%d\n",
               n, (double) lo / lt, lo, lt, dok, dt);
        // (b) mid-stream digit run: the FIG-insertion alignment shift corrupts trailing letters.
        n = 0;
        for (int i = 0; i < 120; i++) pt[n++] = PLAINTEXT[i] - 'A';
        for (int i = 0; i < 4; i++) pt[n++] = 26 + run[i];
        for (int i = 120; i < 200; i++) pt[n++] = PLAINTEXT[i] - 'A';
        plant_and_solve(pt, n, "SECRET", NULL, 2, 6, 1u, &lo, &lt, &dok, &dt);
        printf("    mid-stream digits(L=%d): letters %.3f (%d/%d), digit values %d/%d\n",
               n, (double) lo / lt, lo, lt, dok, dt);
        printf("    (letter recovery of NON-numeric ciphers is the asserted capability above)\n");
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(shared.ngram_data);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
