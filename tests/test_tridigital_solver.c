//
//  Solver calibration / capability tests for the Tridigital cipher (type 81). Built with
//  -DCOLOSSUS_NO_MAIN so solve_cipher() links directly; fixed seeds -> deterministic. Covers:
//  registry validation, the map-independent SEPARATOR-DISCRIMINATION property (the true
//  separator digit's word-length histogram fits English, so the pre-filter keeps it), and a
//  letter-capability PROBE across lengths.
//
//  Tridigital is a DENSE polyphonic cipher (9 symbols carry 26 letters) whose decode is a free
//  max-over-disambiguations, so n-gram gaming is severe and the max-likelihood decode does not
//  equal the truth (the n-gram often prefers the wrong one of a column's 3 letters). Recovery is
//  therefore PARTIAL and high-variance, and the ACA 75-100 range sits well below the blind floor
//  -- a documented limitation, like Tri-Square / CM-Bifid. This test asserts the reliable
//  structural properties and a lenient best-of-seeds capability floor (to catch total breakage),
//  and PRINTS the accuracy-vs-length curve for characterization. Effectively needs -logprob with
//  QUINTGRAMS and a DICTIONARY (the cross-config selector is whole-word coverage).
//
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "tridigital.h"
#include "tridigital_solver.h"
#include <fcntl.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE (char *) "english_quintgrams.txt"
#define NGRAM_SIZE 5
#define DICT_FILE  (char *) "OxfordEnglishWords.txt"

static SharedData shared;

// ~450 letters of spaced English (A Tale of Two Cities opening) -- word breaks drive the separator.
static const char *PLAINTEXT =
    "IT WAS THE BEST OF TIMES IT WAS THE WORST OF TIMES IT WAS THE AGE OF WISDOM IT WAS "
    "THE AGE OF FOOLISHNESS IT WAS THE EPOCH OF BELIEF IT WAS THE EPOCH OF INCREDULITY IT "
    "WAS THE SEASON OF LIGHT IT WAS THE SEASON OF DARKNESS IT WAS THE SPRING OF HOPE IT WAS "
    "THE WINTER OF DESPAIR WE HAD EVERYTHING BEFORE US WE HAD NOTHING BEFORE US WE WERE ALL "
    "GOING DIRECT TO HEAVEN WE WERE ALL GOING DIRECT THE OTHER WAY";

// Map a spaced string into an index stream (letters 0..25, word breaks -> TRI_SPACE), up to a
// target letter count L; stops at the next word boundary so no word is cut. Returns the stream
// length (letters + separators); *nletters receives the letter count.
static int text_of(const char *s, int L, int out[], int *nletters) {
    int n = 0, lc = 0;
    for (int i = 0; s[i]; i++) {
        char c = s[i];
        if (c >= 'A' && c <= 'Z') {
            if (lc >= L) break;                 // reached the target at a word boundary
            out[n++] = c - 'A'; lc++;
        } else {
            if (lc >= L) break;
            if (n > 0 && out[n - 1] != TRI_SPACE) out[n++] = TRI_SPACE;  // one break, no leading/doubled
        }
    }
    while (n > 0 && out[n - 1] == TRI_SPACE) n--;    // no trailing separator
    *nletters = lc;
    return n;
}

// Map-independent word-length fit of a candidate separator digit (a copy of the solver's, for
// the discrimination test): split on the digit, histogram run lengths, score fit to English.
static double sep_fit(const int *digits, int clen, int sep) {
    double hist[25]; for (int b = 0; b < 25; b++) hist[b] = 0.0;
    int nruns = 0, run = 0;
    for (int i = 0; i < clen; i++) {
        if (digits[i] == sep) { if (run > 0) { int b = (run > 25) ? 25 : run; hist[b - 1] += 1.0; nruns++; run = 0; } }
        else run++;
    }
    if (run > 0) { int b = (run > 25) ? 25 : run; hist[b - 1] += 1.0; nruns++; }
    if (nruns == 0) return -1e9;
    double l1 = 0.0;
    for (int b = 0; b < 25; b++) { hist[b] /= nruns; l1 += fabs(hist[b] - english_word_length_frequencies[b]); }
    return -l1;
}

// Rank of the true separator digit among all 10 by word-length fit (0 = best).
static int sep_rank(const int *digits, int clen, int true_sep) {
    double tf = sep_fit(digits, clen, true_sep);
    int rank = 0;
    for (int d = 0; d < 10; d++) if (d != true_sep && sep_fit(digits, clen, d) > tf) rank++;
    return rank;
}

// Plant a cipher from a spaced plaintext index stream and solve it in-process; return the letter
// match count (letters only; separator positions are deterministic and excluded).
static double plant_and_solve(const int *pt, int n, const char *kw_cols, const char *kw_alpha,
                              uint32_t seed) {
    TridigitalGrid g;
    tridigital_build_from_keywords(&g, kw_cols, kw_alpha);
    static int cipher[MAX_CIPHER_LENGTH];
    int clen = tridigital_encrypt(pt, n, &g, cipher);
    static char cs[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < clen; i++) cs[i] = '0' + cipher[i];
    cs[clen] = '\0';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = TRIDIGITAL;
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

    int lim = (res.decrypted_len < n) ? res.decrypted_len : n;
    int letters = 0, ok = 0;
    for (int i = 0; i < lim; i++) {
        if (pt[i] == TRI_SPACE) continue;
        letters++;
        if (res.decrypted[i] == pt[i]) ok++;
    }
    return (letters > 0) ? (double) ok / letters : 0.0;
}

int main(void) {
    init_alphabet(NULL);                        // full 26-letter A..Z, as -type td expects
    CHECK(g_alpha == DEFAULT_ALPHABET_SIZE, "alphabet %d, expected %d", g_alpha, DEFAULT_ALPHABET_SIZE);
    g_ngram_logprob = true;                     // Tridigital effectively needs -logprob (quintgrams)
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }
    load_dictionary(DICT_FILE, &shared.dict, &shared.n_dict_words, &shared.max_dict_word_len, false);
    CHECK(shared.dict != NULL, "could not load %s (the whole-word selector needs it)", DICT_FILE);

    static int pt[MAX_CIPHER_LENGTH];

    // 1) Registry validation: the tuned schedule is applied.
    {
        ColossusConfig cfg; init_config(&cfg);
        cfg.cipher_type = TRIDIGITAL;
        bool applied = apply_cipher_defaults(&cfg, false);
        CHECK(applied, "no SearchDefaults entry for TRIDIGITAL");
        CHECK(cfg.init_temp > 0.2 && cfg.init_temp < 0.4, "unexpected inittemp %.2f", cfg.init_temp);
        CHECK(cfg.n_restarts >= 8, "unexpected a_n_restarts %d", cfg.n_restarts);
        ColossusConfig v; init_config(&v); v.cipher_type = VIGENERE;
        CHECK(!apply_cipher_defaults(&v, false), "VIGENERE unexpectedly has a registry entry");
    }

    // 2) SEPARATOR-DISCRIMINATION property: the map-independent word-length fit keeps the true
    //    separator well within the kept top-K (TRI_KEEP). Unlike Monome-Dinome's exact validity,
    //    this is statistical, so we assert a lenient bound and print the rank.
    {
        printf("\n  separator word-length rank (true sep must be kept in top %d):\n", TRI_KEEP);
        int allok = 1;
        struct { const char *cols, *alpha; } keys[] = {
            {"NOVELCRAFT", "DRAGONFLY"}, {"BLACKHORSE", "CIPHERKEY"}, {"DUMPTRUCKS", "SECRETWORD"},
        };
        for (int c = 0; c < 3; c++) {
            int nl, n = text_of(PLAINTEXT, 300, pt, &nl);
            TridigitalGrid g; tridigital_build_from_keywords(&g, keys[c].cols, keys[c].alpha);
            static int cipher[MAX_CIPHER_LENGTH];
            int clen = tridigital_encrypt(pt, n, &g, cipher);
            int rank = sep_rank(cipher, clen, g.sep_digit);
            printf("    cols=%-11s sep=%d  word-length rank=%d\n", keys[c].cols, g.sep_digit, rank);
            if (rank >= TRI_KEEP) allok = 0;    // must survive the pre-filter (kept = top TRI_KEEP)
        }
        CHECK(allok, "true separator ranked outside the kept top-%d by word-length fit", TRI_KEEP);
    }

    // 3) LETTER capability PROBE across lengths (best of a few seeds -- the recovery is
    //    high-variance). Asserts only a lenient breakage floor; prints the curve. Above the floor
    //    the separator is recovered and much of the text is readable, but exact recovery is capped
    //    by the polyphonic ambiguity (below the 99% bar at every practical length).
    printf("\n  Tridigital letter recovery vs length (keywords NOVELCRAFT / DRAGONFLY, best of 3 seeds):\n");
    {
        int lens[] = {250, 400};
        double best_overall = 0.0;
        for (int li = 0; li < 2; li++) {
            int nl, n = text_of(PLAINTEXT, lens[li], pt, &nl);
            double best = 0.0;
            for (uint32_t s = 1; s <= 3; s++) {
                double f = plant_and_solve(pt, n, "NOVELCRAFT", "DRAGONFLY", s);
                if (f > best) best = f;
            }
            printf("    L=%3d letters: best=%.3f\n", nl, best);
            if (best > best_overall) best_overall = best;
        }
        CHECK(best_overall > 0.40, "capability floor: best recovery %.3f (want > 0.40, breakage catch)", best_overall);
        printf("    (partial + high-variance by design: dense polyphonic decode, below the 99%% floor)\n");
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(shared.ngram_data);
    if (shared.dict) free_dictionary(shared.dict, shared.n_dict_words);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
