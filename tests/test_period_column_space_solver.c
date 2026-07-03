//
//  In-process capability tests for the space-robust Period column order solver
//  (period_column_space_search).
//
//  Framework-free: build with `make testopt`. The solver keeps the observed
//  ciphertext (spaces included -- they are structurally significant) and searches
//  FROM SCRATCH over small EDITS to it, modelling a cipher in which a few
//  characters were accidentally dropped OR spuriously added:
//    * INSERT a blank "gap" cell (repairs a dropped char -- observed too short);
//    * DELETE an observed cell (repairs an added char -- observed too long).
//  Either edit re-aligns the columns and lets the edited length reach a
//  complete-grid factorisation; the edit positions are annealed to maximise the
//  exhaustive period-column n-gram fitness.
//
//  Checks:
//    1. (0,0) superset: on a CLEAN cipher the solver recovers the plaintext
//       exactly, i.e. it is a strict superset of the deterministic solver.
//    2. dropped-character repair: delete one cipher cell, and the solver must
//       INSERT a gap that restores the length and recover all but the lost letter.
//    3. two dropped characters (printed + a loose floor).
//    4. added-character repair: insert a spurious letter, and the solver must
//       DELETE it, restoring the exact cipher and recovering the plaintext.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "scoring.h"                    // load_ngrams
#include "period_column_space_solver.h" // period_column_space_search
#include "period_column_solver.h"       // period_column_transform, PCStage

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4

static float *ngrams = NULL;

// Encode a spaced A..Z string the way decode_cipher/ord does: letters -> 0..25,
// anything else -> a reversible negative sentinel.
static int encode_spaced(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++) {
        unsigned char c = (unsigned char) s[i];
        int v = g_char_to_idx[toupper(c)];
        out[n++] = (v >= 0) ? v : (-(int) c - 1);
    }
    return n;
}

static int count_matches(const int *a, const int *b, int len) {
    int m = 0;
    for (int i = 0; i < len; i++) if (a[i] == b[i]) m++;
    return m;
}

int main(void) {
    init_alphabet(NULL);                    // full 26-letter alphabet (no J->I)
    seed_rand(0x9E3779B9u);                 // deterministic RNG for reproducibility
    ngrams = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!ngrams) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    // A natural spaced plaintext, truncated to L = 168 so the grid factorises richly
    // (168 = 4x42 = 6x28 = 8x21 = 24x7 = 56x3 ...). Spaces are real grid cells.
    const char *SPACED =
        "IT IS A TRUTH UNIVERSALLY ACKNOWLEDGED THAT A SINGLE MAN IN POSSESSION OF A "
        "GOOD FORTUNE MUST BE IN WANT OF A WIFE HOWEVER LITTLE KNOWN THE FEELINGS OR "
        "VIEWS OF SUCH A MAN MAYBE";
    static int pt[MAX_CIPHER_LENGTH];
    int full = encode_spaced(SPACED, pt, MAX_CIPHER_LENGTH);
    const int L = 168;
    if (full < L) { printf("FAIL: SPACED too short (%d < %d)\n", full, L); return 1; }

    // Encrypt with a single forward stage.
    const int DX = 8, P = 3, UTP = 0;
    static int cipher[MAX_CIPHER_LENGTH];
    period_column_transform(pt, cipher, L, DX, P, UTP);

    // Output buffers (need L + max_ins room).
    static int out_edited[MAX_CIPHER_LENGTH], out_dec[MAX_CIPHER_LENGTH];
    static int gap_pos[MAX_CIPHER_LENGTH], del_pos[MAX_CIPHER_LENGTH];
    PCStage stages[2];
    int nstg = 0, ngaps = 0, ndels = 0, out_len = 0;

    // --- 1. (0,0) superset: clean cipher recovers exactly ----------------------------
    {
        double s = period_column_space_search(cipher, L,
            /*max_ins=*/0, /*max_dels=*/0, /*inner_depth=*/1, /*final_depth=*/1,
            /*n_restarts=*/1, /*n_hillclimbs=*/1, 0.10, 0.001,
            ngrams, NGRAM_SIZE, 1.0f, 0.0f, NULL, false,
            out_edited, &out_len, out_dec, stages, &nstg, gap_pos, &ngaps, del_pos, &ndels, false);
        (void) s;
        CHECK(out_len == L && ngaps == 0 && ndels == 0,
              "(0,0) keeps length/no edits (len=%d gaps=%d drops=%d)", out_len, ngaps, ndels);
        CHECK(count_matches(out_dec, pt, L) == L, "(0,0) clean cipher recovers exactly (%d/%d)",
              count_matches(out_dec, pt, L), L);
    }

    // --- 2. dropped-character repair (INSERT a gap) ----------------------------------
    // Delete one LETTER cell from the ciphertext; the solver must insert one gap to
    // restore length 168 and recover all but that letter.
    {
        int drop = 80;
        while (drop < L && cipher[drop] < 0) drop++;   // ensure we drop a letter, not a space
        static int cdrop[MAX_CIPHER_LENGTH];
        int m = 0;
        for (int i = 0; i < L; i++) if (i != drop) cdrop[m++] = cipher[i];   // length L-1

        double s = period_column_space_search(cdrop, L - 1,
            /*max_ins=*/2, /*max_dels=*/0, /*inner_depth=*/1, /*final_depth=*/1,
            /*n_restarts=*/3, /*n_hillclimbs=*/1000, 0.10, 0.001,
            ngrams, NGRAM_SIZE, 1.0f, 0.0f, NULL, false,
            out_edited, &out_len, out_dec, stages, &nstg, gap_pos, &ngaps, del_pos, &ndels, false);

        int matched = count_matches(out_dec, pt, L);
        printf("\n  dropped-1 (insert): len %d->%d, gaps=%d drops=%d, match %d/%d, score %.2f\n",
               L - 1, out_len, ngaps, ndels, matched, L, s);
        CHECK(out_len == L, "dropped-1 restores length to %d (got %d)", L, out_len);
        CHECK(ngaps >= 1, "dropped-1 inserts at least one gap (got %d)", ngaps);
        CHECK(matched >= L - 2, "dropped-1 recovers all but the lost letter (%d/%d)", matched, L);
    }

    // --- 3. two dropped characters (printed + a loose floor) -------------------------
    {
        int d1 = 55, d2 = 120;
        while (d1 < L && cipher[d1] < 0) d1++;
        while (d2 < L && cipher[d2] < 0) d2++;
        static int cdrop2[MAX_CIPHER_LENGTH];
        int m = 0;
        for (int i = 0; i < L; i++) if (i != d1 && i != d2) cdrop2[m++] = cipher[i];  // length L-2

        double s = period_column_space_search(cdrop2, L - 2,
            /*max_ins=*/2, /*max_dels=*/0, /*inner_depth=*/1, /*final_depth=*/1,
            /*n_restarts=*/4, /*n_hillclimbs=*/1500, 0.10, 0.001,
            ngrams, NGRAM_SIZE, 1.0f, 0.0f, NULL, false,
            out_edited, &out_len, out_dec, stages, &nstg, gap_pos, &ngaps, del_pos, &ndels, false);

        int matched = count_matches(out_dec, pt, L);
        printf("\n  dropped-2 (insert): len %d->%d, gaps=%d drops=%d, match %d/%d, score %.2f\n",
               L - 2, out_len, ngaps, ndels, matched, L, s);
        CHECK(matched >= (int) (0.85 * L), "dropped-2 recovers >= 85%% (%d/%d)", matched, L);
    }

    // --- 4. added-character repair (DELETE a cell) -----------------------------------
    // Splice a spurious letter (Q) into the ciphertext; deleting exactly that cell
    // restores the original cipher, so the plaintext recovers cleanly.
    {
        int add = 90;
        static int cadd[MAX_CIPHER_LENGTH];
        int m = 0;
        for (int i = 0; i < L; i++) {
            if (i == add) cadd[m++] = g_char_to_idx['Q'];   // spurious inserted letter
            cadd[m++] = cipher[i];
        }
        // cadd length is L + 1.

        double s = period_column_space_search(cadd, L + 1,
            /*max_ins=*/0, /*max_dels=*/2, /*inner_depth=*/1, /*final_depth=*/1,
            /*n_restarts=*/3, /*n_hillclimbs=*/1000, 0.10, 0.001,
            ngrams, NGRAM_SIZE, 1.0f, 0.0f, NULL, false,
            out_edited, &out_len, out_dec, stages, &nstg, gap_pos, &ngaps, del_pos, &ndels, false);

        int matched = count_matches(out_dec, pt, L);
        printf("\n  added-1 (delete): len %d->%d, gaps=%d drops=%d, match %d/%d, score %.2f\n",
               L + 1, out_len, ngaps, ndels, matched, L, s);
        CHECK(out_len == L, "added-1 restores length to %d (got %d)", L, out_len);
        CHECK(ndels >= 1, "added-1 deletes at least one cell (got %d)", ndels);
        CHECK(matched >= L - 1, "added-1 recovers the plaintext (%d/%d)", matched, L);
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(ngrams);
    return failures ? 1 : 0;
}
