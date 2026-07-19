//
//  Solver calibration / capability tests for the Syllabary cipher (type 85). Built with
//  -DCOLOSSUS_NO_MAIN so solve_cipher() links directly; fixed seeds -> deterministic.
//
//  Syllabary is a SUBSTITUTION over 100 numeric codes -> 100 KNOWN 1-3 letter tokens (a
//  length-changing decode). Blind recovery is a DOCUMENTED LIMITATION well below the ACA range:
//  the 100-token bijection is severely undersampled (~1-1.5 samples/code) and, critically, the
//  token set carries common syllables (THE, AND, ING, RE, ...) that let the search assemble
//  high-scoring English REGARDLESS of the true map -- so the n-gram global max GAMES the true
//  plaintext (like the Checkerboard complex case / jarl5). This test therefore validates the
//  registry and the solve MACHINERY (the tiled length-fair scoring, the composite-map search,
//  and the report all function -- a broken pipeline yields gibberish with no words), and PRINTS
//  the true-letter recovery to characterize the gaming. The DECODE correctness itself is pinned
//  by test_syllabary.c (all four ACA isolog vectors). A small search budget keeps this fast.
//
#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "syllabary.h"
#include "syllabary_solver.h"
#include <fcntl.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE (char *) "english_quintgrams.txt"
#define NGRAM_SIZE 5
#define DICT_FILE  (char *) "OxfordEnglishWords.txt"

static SharedData shared;

static const char *PLAINTEXT =
    "ITWASTHEBESTOFTIMESITWASTHEWORSTOFTIMESITWASTHEAGEOFWISDOMITWASTHEAGEOFFOOLISHNESS"
    "ITWASTHEEPOCHOFBELIEFITWASTHEEPOCHOFINCREDULITYITWASTHESEASONOFLIGHTITWASTHESEASON"
    "OFDARKNESSITWASTHESPRINGOFHOPEITWASTHEWINTEROFDESPAIRWEHADEVERYTHINGBEFOREUSWEHADNO"
    "THINGBEFOREUSWEWEREALLGOINGDIRECTTOHEAVENWEWEREALLGOINGDIRECTTHEOTHERWAY";

static int text_of(const char *s, int L, int out[]) {
    int n = 0;
    for (int i = 0; s[i] && n < L; i++)
        if (s[i] >= 'A' && s[i] <= 'Z') out[n++] = s[i] - 'A';
    return n;
}

// Plant a Syllabary cipher and solve it in-process with a SMALL budget; return the true-letter
// recovery fraction and (via *nwords) the dictionary-word count of the recovered decode.
static double plant_and_solve(const int *pt, int n, uint32_t sqseed, int mode,
                              uint32_t solveseed, int *nwords, int *declen) {
    SyllabarySquare sq;
    syllabary_build_random(&sq, sqseed);
    static int codes[MAX_CIPHER_LENGTH];
    int clen = syllabary_encrypt(pt, n, &sq, mode, codes);
    static char cs[2 * MAX_CIPHER_LENGTH + 1];
    int p = 0;
    for (int i = 0; i < clen; i++) { cs[p++] = '0' + codes[i] / 10; cs[p++] = '0' + codes[i] % 10; }
    cs[p] = '\0';

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = SYLLABARY;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.dictionary_present = (shared.dict != NULL);
    strcpy(cfg.ciphertext_file, "in-process-test");
    apply_cipher_defaults(&cfg, false);
    cfg.n_restarts = 4; cfg.n_hill_climbs = 50000;   // small budget: this test checks health, not capability

    SolveResult res;
    res.solved = false; res.decrypted_len = 0; res.n_words = 0;
    fflush(stdout);
    int saved = dup(fileno(stdout));
    if (freopen("/dev/null", "w", stdout) == NULL) { /* proceed anyway */ }
    seed_rand(solveseed);
    solve_cipher(cs, (char *) "", &cfg, &shared, &res);
    fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    clearerr(stdout);

    *nwords = res.n_words;
    *declen = res.decrypted_len;
    int lim = (res.decrypted_len < n) ? res.decrypted_len : n;
    int ok = 0;
    for (int i = 0; i < lim; i++) if (res.decrypted[i] == pt[i]) ok++;
    return (n > 0) ? (double) ok / n : 0.0;
}

int main(void) {
    init_alphabet(NULL);
    CHECK(g_alpha == DEFAULT_ALPHABET_SIZE, "alphabet %d, expected %d", g_alpha, DEFAULT_ALPHABET_SIZE);
    g_ngram_logprob = true;
    shared.ngram_data = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }
    load_dictionary(DICT_FILE, &shared.dict, &shared.n_dict_words, &shared.max_dict_word_len, false);

    static int pt[MAX_CIPHER_LENGTH];

    // 1) Registry validation: the tuned schedule is applied; a non-registry type is untouched.
    {
        ColossusConfig cfg; init_config(&cfg);
        cfg.cipher_type = SYLLABARY;
        bool applied = apply_cipher_defaults(&cfg, false);
        CHECK(applied, "no SearchDefaults entry for SYLLABARY");
        CHECK(cfg.init_temp > 0.1499 && cfg.init_temp < 0.1501, "unexpected inittemp %.4f", cfg.init_temp);
        CHECK(cfg.n_restarts >= 16, "unexpected a_n_restarts %d", cfg.n_restarts);
        ColossusConfig v; init_config(&v); v.cipher_type = VIGENERE;
        CHECK(!apply_cipher_defaults(&v, false), "VIGENERE unexpectedly has a registry entry");
    }

    // 2) MACHINERY HEALTH + gaming characterization. The composite-map search + tiled length-fair
    //    scoring + report must all function: the decode is a plausible length and (because the
    //    n-gram optimizer assembles fluent English) yields real dictionary words. True-letter
    //    recovery stays near chance -- the documented gaming -- and is printed, not asserted.
    printf("\n  Syllabary machinery + gaming characterization (random square, greedy isolog, quintgrams):\n");
    {
        int lens[] = {120, 200};
        int health_ok = 1;
        for (int li = 0; li < 2; li++) {
            int n = text_of(PLAINTEXT, lens[li], pt);
            int nwords = 0, declen = 0;
            double rec = plant_and_solve(pt, n, 3u, SYL_TOK_GREEDY, 1u, &nwords, &declen);
            printf("    L=%3d letters: decode=%3d letters, words=%3d, true-letter recovery=%.3f (gamed)\n",
                   n, declen, nwords, rec);
            if (declen < n / 2 || declen > 3 * n) health_ok = 0;   // tiling/decode produced a sane length
            if (shared.dict && nwords < 3) health_ok = 0;          // engine reached English-scoring output
        }
        CHECK(health_ok, "solve machinery health check failed (decode length or word yield implausible)");
        printf("    (blind recovery is a documented limitation: 100-token set with syllable magnets games the n-gram)\n");
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(shared.ngram_data);
    if (shared.dict) free_dictionary(shared.dict, shared.n_dict_words);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}
