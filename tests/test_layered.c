//
//  In-process unit tests for the layered (Paradigm) solvers' shared engine.
//
//  Framework-free: build with `make test`. The two new types -- QUAG_TRANS (quagtrans)
//  and HILL_QUAG (hillquag) -- both rest on quag_strip_and_refine(): strip a plain
//  period-P Quagmire III of English by the monogram statistic, then polish the cycleword
//  by n-gram coordinate ascent over its COMPONENTS. This file plants known Quagmire
//  ciphers (KRYPTOS keyed alphabet) and asserts that engine recovers them, which is the
//  core the layered attacks depend on. The FULL layered capability (outer Quagmire/Hill
//  over an inner transposition/Quagmire) is asserted end-to-end by the fixed-seed
//  ciphers/tests/run_tests.sh cases paradigm_pk{3,4,6,7} (all real puzzles, ~100%).
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "scoring.h"
#include "layered_solver.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

#define NGRAM_FILE "ngram_data/english/english_quadgrams.txt"
#define NGRAM_SIZE 4

static SharedData shared;

static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNE"
    "MUSTBEINWANTOFAWIFEHOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHIS"
    "FIRSTENTERINGANEIGHBOURHOODTHISTRUTHISSOWELLFIXEDINTHEMINDSOFTHESURROUNDING"
    "FAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOROTHEROFTHEIRDAUGHTERS"
    "MYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLETATLAST";

// Plant a plain Quagmire III (KRYPTOS keyed alphabet) of the first `plen` plaintext
// letters under the cycleword `cw`. Fills prepared[] (expected plaintext, 0..25) and
// cipher[] (0..25). Returns the length used.
static int plant_quag(const char *cw, int plen, int pt_kw[], int ct_kw[],
                      int prepared[], int cipher[]) {
    int n = 0;
    for (int i = 0; PLAINTEXT[i] && n < plen; i++) {
        int idx = g_char_to_idx[toupper((unsigned char) PLAINTEXT[i]) & 127];
        if (idx >= 0 && idx < 26) prepared[n++] = idx;
    }
    int P = 0, cwi[MAX_CYCLEWORD_LEN];
    for (int i = 0; cw[i]; i++) cwi[P++] = g_char_to_idx[toupper((unsigned char) cw[i]) & 127];
    // Quagmire III encrypt: pt_kw == ct_kw; C = ct_kw[(posn_pt(P) + posn_ct(K)) mod 26].
    int pt_inv[ALPHABET_SIZE], ct_inv[ALPHABET_SIZE];
    for (int i = 0; i < 26; i++) { pt_inv[pt_kw[i]] = i; ct_inv[ct_kw[i]] = i; }
    for (int i = 0; i < n; i++) {
        int s = ct_inv[cwi[i % P]];
        cipher[i] = ct_kw[(pt_inv[prepared[i]] + s) % 26];
    }
    return n;
}

static double strip_frac(const char *cw, int plen) {
    int pt_kw[ALPHABET_SIZE], ct_kw[ALPHABET_SIZE];
    make_keyed_alphabet("KRYPTOS", pt_kw);
    for (int i = 0; i < 26; i++) ct_kw[i] = pt_kw[i];

    int prepared[MAX_CIPHER_LENGTH], cipher[MAX_CIPHER_LENGTH];
    int n = plant_quag(cw, plen, pt_kw, ct_kw, prepared, cipher);
    int P = (int) strlen(cw);

    ColossusConfig cfg;
    init_config(&cfg);
    cfg.cipher_type = QUAG_TRANS;
    cfg.ngram_size = NGRAM_SIZE;
    cfg.weight_ngram = 1.0f;

    int cw_out[MAX_CYCLEWORD_LEN], pt_out[MAX_CIPHER_LENGTH];
    quag_strip_and_refine(&cfg, cipher, n, pt_kw, ct_kw, P, shared.ngram_data,
                          NULL, NULL, 0, cw_out, pt_out);
    int ok = 0;
    for (int i = 0; i < n; i++) if (pt_out[i] == prepared[i]) ok++;
    return (double) ok / (double) n;
}

int main(void) {
    ColossusConfig cfg; init_config(&cfg);          // defaults g_alpha=26, fills g_char_to_idx
    g_ngram_logprob = true;                          // strip+refine uses the discriminating fitness
    shared.ngram_data = load_ngrams((char *) NGRAM_FILE, NGRAM_SIZE, false);
    if (!shared.ngram_data) { printf("FAIL: could not load %s (run from the src dir)\n", NGRAM_FILE); return 1; }

    // The monogram strip + n-gram component refine recovers a plain Quagmire III at the
    // periods the layered attacks strip (PK6/PK7 use 6; PK3/PK4 collapse to 40/45).
    double f6  = strip_frac("PORTAL", 300);
    CHECK(f6  > 0.99, "quag_strip_and_refine period 6 recovered only %.1f%%",  100 * f6);
    double f12 = strip_frac("ABSCISSA", 360);        // period 8
    CHECK(f12 > 0.99, "quag_strip_and_refine period 8 recovered only %.1f%%", 100 * f12);

    printf("test_layered: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
