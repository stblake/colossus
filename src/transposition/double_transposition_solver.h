#ifndef DOUBLE_TRANSPOSITION_SOLVER_H
#define DOUBLE_TRANSPOSITION_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Double columnar transposition -- DIVIDE AND CONQUER (TRANSCOL2_DC)
// =====================================================================
//
//  A dedicated solver for the DOUBLE COLUMNAR TRANSPOSITION cipher that,
//  unlike TRANSCOL2 (which hill-climbs K1 and K2 *in parallel* over the
//  combined K1!*K2! keyspace), attacks the two keys SEPARATELY, following
//  Lasry, Kopal & Wacker, "Solving the Double Transposition Challenge with
//  a Divide and Conquer Approach" (Cryptologia, 2014).
//
//  The pivot is the Index of Digraphic Potential (IDP): a fitness for a
//  putative SECOND key K2 that can be computed WITHOUT knowing the first
//  key K1. Once the second transposition is undone with the correct K2 the
//  residue is a *single* columnar transposition of the plaintext under K1,
//  which is easy to break. So the O(K1!*K2!) joint search collapses to an
//  O(K2!) search for K2 (scored by the IDP) followed by an O(K1!) single-
//  columnar finish -- the reason this cracks key lengths (20-25) at which
//  parallel hill-climbing (TRANSCOL2 / the paper's "Step 1") fails.
//
//  Convention: the cipher is C = colenc(colenc(P, K1), K2); each stage is
//  a standard columnar transposition in colossus's convention (see
//  decrypt_columnar: plaintext written row-major, leftmost `len % K`
//  columns one cell taller, columns read off in key order, top-to-bottom
//  for COL_READ_TB). Both stages use cfg->read_direction.

// Build a 26x26 English bigram log10-probability table by marginalizing an
// n-gram COUNT file (e.g. english_quadgrams.txt): every internal bigram of each
// listed n-gram contributes its count, the result is normalized to log10
// probabilities and unseen bigrams get a rare-but-not-impossible floor. Returns a
// malloc'd double[26*26] (index a*26 + b for the ordered pair (a,b)), or NULL on
// failure. ngram_size must be >= 2. (The IDP needs digraph log-frequencies, which
// the loaded quadgram/quintgram fitness table does not expose directly.)
double *dct_load_bigrams(const char *ngram_file, int ngram_size);

// Index of Digraphic Potential of a putative second key K2 (order2[0..L2-1]),
// scored WITHOUT knowing K1. Undoes the second columnar transposition (dir2 =
// COL_READ_TB/BT) to obtain the intermediate text, hypothesizes a first key of
// length L1, and returns the total weight of a greedy digraph-adjacency
// reconstruction of the intermediate's L1 columns -- higher when K2 is closer to
// correct. `work` is caller scratch of length >= len (receives the intermediate).
// If `chain_out` is non-NULL it receives the reconstructed left-to-right column
// read-order (a candidate inverse-K1, used to warm-start the single-columnar
// finish). bg is the table from dct_load_bigrams.
double dct_idp(const int *cipher, int len, int L1, int L2, int dir2,
               const int *order2, const double *bg, int *work, int *chain_out);

// Recovered solution of the divide-and-conquer search (filled by dct_solve_core).
typedef struct {
    double score;                       // n-gram score of the recovered plaintext
    int L1, L2;                         // recovered key lengths (columns)
    int dir1, dir2;                     // recovered read directions (COL_READ_TB/BT)
    int order1[MAX_COLS], order2[MAX_COLS];   // recovered column read orders
    int plaintext[MAX_CIPHER_LENGTH];   // recovered plaintext (0..25 indices)
    int len;
} DCTResult;

// Divide-and-conquer solve core (no I/O): screen every (L1, L2, dir) length hypothesis
// by a short IDP climb, refine the top few, break the residual single columnar, and
// return the best plaintext in *res. cfg supplies the column range (min_cols/max_cols),
// read direction, n-gram size, the refinement budget (n_restarts/n_hill_climbs), and the
// crib weights (weight_ngram/weight_crib). ngram_data is the loaded fitness table; bg is
// from dct_load_bigrams. Cribs (over PLAINTEXT positions) are blended into the finish and
// rotation scoring via state_score -- they guide the K1 single-columnar finish and the
// selection across configs (n_cribs==0 => pure n-gram, unchanged). Shared by the CLI
// entry (which then reports) and the regression test.
void dct_solve_core(const int *cipher, int len, ColossusConfig *cfg,
                    float *ngram_data, const double *bg,
                    int *crib_indices, int *crib_positions, int n_cribs,
                    DCTResult *res);

// Full solver entry, dispatched from solve_cipher().
void solve_double_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs);
#endif
