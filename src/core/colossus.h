#ifndef COLOSSUS_H
#define COLOSSUS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h> 
#include <string.h>
#include <math.h>
#include <time.h>
#include <strings.h>
#include <unistd.h>
#include <stdint.h>

#define KRYPTOS 0
#define CRIB_CHECK 0
#define PARTIAL_CRIB_MATCH 1

#define VIGENERE 0
#define QUAGMIRE_1 1
#define QUAGMIRE_2 2
#define QUAGMIRE_3 3
#define QUAGMIRE_4 4
#define BEAUFORT   5
#define PORTA      6
#define AUTOKEY_0  7 
#define AUTOKEY_1  8 
#define AUTOKEY_2  9 
#define AUTOKEY_3  10
#define AUTOKEY_4  11
#define AUTOKEY_BEAU 12
#define AUTOKEY_PORTA 13
#define TRANSMATRIX    14
#define TRANSPEROFFSET 15
#define TRANSPOSITION  16
#define TRANSCOL       17   // single columnar transposition (dedicated solver)
#define TRANSCOL2      18   // double (nested) columnar transposition
#define RAILFENCE      19   // rail fence + variant rail fence (phase-offset sweep)
#define ROUTE          20   // route transposition (snake / spiral over an R x C grid)
#define AMSCO          21   // Amsco (alternating 1/2-letter columnar)
#define MYSZKOWSKI     22   // Myszkowski (columnar with tied keyword ranks)
#define REDEFENCE      23   // Redefence (keyed rail fence)
#define CADENUS        24   // Cadenus (rotated-column transposition, 25 rows)
#define NIHILIST       25   // Nihilist transposition (single perm on rows + columns)
#define SWAGMAN        26   // Swagman (N x N Latin-square column transposition)
#define GRILLE         27   // Turning grille
#define INDEP_PERIODIC 28   // period-P substitution with P INDEPENDENT mixed alphabets
#define HOMOPHONIC     29   // homophonic substitution (ciphertext alphabet > plaintext)
#define PLAYFAIR       30   // Playfair (digraphic substitution over a 5x5 keyed grid)
#define BIFID          31   // Bifid (Delastelle): fractionation over a keyed Polybius square
#define TRIFID         32   // Trifid (Delastelle): fractionation over a 3x3x3 keyed cube
#define HILL           33   // Hill (polygraphic substitution by a k x k matrix mod 26)
#define GRONSFELD      34   // Gronsfeld (Vigenere with a numeric key: per-column shifts 0..9)
#define PHILLIPS       35   // Phillips ("Row" type): 8 keyed-square substitution, period 40
#define PHILLIPS_C     36   // Phillips-C: the column-shift dual of the Row square generation
#define PHILLIPS_RC    37   // Phillips-RC: rows for squares 2-5, columns for squares 6-8
#define TWO_SQUARE     38   // Two-Square horizontal (ACA): two keyed 5x5 squares side by side
#define TWO_SQUARE_V   39   // Two-Square vertical (Wikipedia): two keyed squares stacked (self-inverse)
#define FOUR_SQUARE    40   // Four-Square: two keyed ciphertext squares + two standard plaintext squares
#define TRANSCOL_L     41   // columnar with a within-column row permutation L (exact seam best-L)
#define TRANSROUTECOL  42   // two-stage chain: fixed read-route global + searched column key (seam best-L)
#define TRANSTILE      43   // sub-grid / tile (h x w) transposition: uniform tile cell permutation
#define ADFGX          44   // ADFGX: 5x5 keyed-square fractionation + keyed columnar transposition
#define ADFGVX         45   // ADFGVX: 6x6 (36-symbol) keyed-square fractionation + keyed columnar
#define NIHILIST_SUB      46   // Nihilist Substitution (carry): periodic additive over a keyed Polybius square
#define NIHILIST_SUB_NC   47   // Nihilist Substitution, no-carry (per-digit add mod 10)
#define NIHILIST_SUB_M100 48   // Nihilist Substitution, add mod 100
#define GROMARK           49   // Gromark: keyed-alphabet substitution + chain-addition running key
#define GROMARK_PERIODIC  50   // Periodic Gromark: + per-group offset (period = keyword length)
#define NICODEMUS          51  // Nicodemus: periodic Vigenere + per-block columnar (one keyword)
#define NICODEMUS_VARIANT  52  // Nicodemus, Variant substitution (C = P - k)
#define NICODEMUS_BEAUFORT 53  // Nicodemus, Beaufort substitution (C = k - P, reciprocal)
#define BAZERIES           54  // Bazeries: keyed-square substitution + digit-grouped reversal (one number key)
#define PORTAX             55  // Portax: periodic digraphic Porta (vertical pairs over a Porta slide)
#define PROGKEY            56  // Progressive Key (Vigenere base): periodic key + per-group constant drift
#define PROGKEY_VAR        57  // Progressive Key, Variant base (C = P - k per pass)
#define PROGKEY_BEAU       58  // Progressive Key, Beaufort base (C = k - P per pass, reciprocal)
#define SLIDEFAIR          59  // Slidefair: periodic digraphic Vigenere (rectangle over a shift slide)
#define SLIDEFAIR_VAR      60  // Slidefair, Variant slide (bottom = col - k)
#define SLIDEFAIR_BEAU     61  // Slidefair, Beaufort slide (bottom = k - col)
#define SERIATED_PLAYFAIR  62  // Seriated Playfair: digraphic Playfair over vertical pairs of a two-row seriated layout (period swept)
#define DIGRAFID           63  // Digrafid: digraphic fractionation over two keyed 27-symbol alphabets (period swept)
#define CM_BIFID           64  // CM Bifid (Conjugated Matrix Bifid): Bifid fractionation over two keyed Polybius squares (period swept)
#define TRI_SQUARE         65  // Tri-Square: digraphic substitution over three keyed 5x5 squares (digraph -> trigraph)
#define INTERRUPTED_KEY      66  // Interrupted Key (Vigenere base): periodic keyword that resets to key letter 1 at break points
#define INTERRUPTED_KEY_VAR  67  // Interrupted Key, Variant base (C = P - k)
#define INTERRUPTED_KEY_BEAU 68  // Interrupted Key, Beaufort base (C = k - P, reciprocal)
#define CONDI              69  // Condi: plaintext-feedback substitution over a keyed alphabet (shift = position of the previous plaintext letter)
#define FRAC_MORSE         70  // Fractionated Morse: Morse fractionation over a keyed 26-letter alphabet (trigraph substitution)
#define PERIOD_COLUMN      71  // Period column order: periodic-column-permutation transposition (AZdecrypt), composable
#define PERIOD_COLUMN_SPACE 72 // Period column order, space-robust: inserts searched blank/gap cells (dropped-char repair)

#define GRONSFELD_DIGITS 10     // Gronsfeld key digits are 0..9 (the shift domain, vs 26)

// Progressive Key base type (which periodic cipher both the keyword and the
// progressive-drift passes use) and its shift math, shared by the primitive
// (progkey.c) and the solver (progkey_solver.c). The progression index ranges
// 0..PROGKEY_MAX_PROG (0 = no drift -> a plain periodic Vigenere/Variant cipher;
// for the Beaufort base the group-0 pass is still a reflection, not identity).
#define PROGKEY_MAX_PROG 25
enum { PROGKEY_BASE_VIG = 0, PROGKEY_BASE_VAR = 1, PROGKEY_BASE_BEAU = 2 };

static inline int progkey_base(int cipher_type) {
    if (cipher_type == PROGKEY_VAR)  return PROGKEY_BASE_VAR;
    if (cipher_type == PROGKEY_BEAU) return PROGKEY_BASE_BEAU;
    return PROGKEY_BASE_VIG;
}
// progkey_base_encrypt / progkey_base_decrypt are defined just below the
// ALPHABET_SIZE definition (they use it as the mod base).

// Interrupted Key (intkey.c / intkey_solver.c). A periodic Vig/Var/Beau keyword whose key index
// RESETS to the first key letter at break points (then increments mod P again). The base shift
// math is identical to Progressive Key's, so IK_BASE_* reuse the PROGKEY_BASE_* values and the
// primitives call progkey_base_encrypt/decrypt. `intkey_base` maps the three type codes to a base;
// the interruption STRATEGY (how break points arise) is enumerated by the solver (IK_STRAT_*).
enum { IK_BASE_VIG = PROGKEY_BASE_VIG, IK_BASE_VAR = PROGKEY_BASE_VAR, IK_BASE_BEAU = PROGKEY_BASE_BEAU };
enum {
    IK_STRAT_CT     = 0,   // reset AFTER a chosen ciphertext letter (breaks known from cipher alone)
    IK_STRAT_PT     = 1,   // reset AFTER a chosen plaintext letter (breaks causal, keyword-dependent)
    IK_STRAT_BREAKS = 2,   // reset at user-supplied break positions (random / word-division, known)
    IK_STRAT_JOINT  = 3    // blind: keyword + break-mask searched jointly (random, characterized)
};

static inline int intkey_base(int cipher_type) {
    if (cipher_type == INTERRUPTED_KEY_VAR)  return IK_BASE_VAR;
    if (cipher_type == INTERRUPTED_KEY_BEAU) return IK_BASE_BEAU;
    return IK_BASE_VIG;
}

#define GROMARK_PRIMER_LEN 5    // basic Gromark standard primer length (ACA convention)
#define GROMARK_MAX_PRIMER 26   // max primer length == max period (periodic: keyword len <= 26)

#define ADFGX_SIDE  5           // ADFGX Polybius square side (5x5, 25 letters, J->I)
#define ADFGVX_SIDE 6           // ADFGVX Polybius square side (6x6, 36 symbols A..Z + 0..9)

#define HILL_MAX_K   5          // largest Hill block size (matrix dimension) supported
#define HILL_MAX_KEY (HILL_MAX_K * HILL_MAX_K)  // largest k*k matrix (=25), fits the key lane

#define PLAYFAIR_SIDE 5         // Playfair grid side (the classic 5x5)
#define PLAYFAIR_GRID 25        // Playfair grid size (PLAYFAIR_SIDE * PLAYFAIR_SIDE)

#define TWO_SQ_HORIZONTAL 0     // Two-Square arrangement: squares side by side (ACA)
#define TWO_SQ_VERTICAL   1     // Two-Square arrangement: squares stacked (Wikipedia, self-inverse)
#define SQUARE_SIDE 5           // Two/Four-Square square side (the classic 5x5)
#define SQUARE_GRID 25          // one square's cell count (SQUARE_SIDE * SQUARE_SIDE)
#define SQUARE_MAX_GRID 36      // largest single square (6x6) the side-generic primitives handle

#define BIFID_MAX_SIDE 6        // largest Bifid square side supported (6x6, 36 cells)
#define BIFID_MAX_GRID 36       // BIFID_MAX_SIDE * BIFID_MAX_SIDE

#define NIHILIST_SUB_SIDE     5   // Nihilist Substitution Polybius square side (5x5, 25 letters, J->I)
#define NIHILIST_SUB_GRID     25  // NIHILIST_SUB_SIDE^2
#define NIHILIST_SUB_MAX_SIDE 6   // largest side the side-generic primitive/tests support
#define NIHILIST_SUB_MAX_GRID 36  // NIHILIST_SUB_MAX_SIDE^2
// Nihilist Substitution addition conventions (one per sub-type code; chosen by cipher_type).
enum { NIH_ADD_CARRY = 0,    // integer add with carry (ACA standard): cipher 22..110
       NIH_ADD_NOCARRY = 1,  // per-digit add mod 10, no carry: cipher 00..99 (2-digit)
       NIH_ADD_MOD100 = 2 }; // 2-digit add mod 100: cipher 00..99

#define TRIFID_SIDE 3           // Trifid cube side (the classic 3x3x3)
#define TRIFID_CELLS 27         // Trifid cube size (TRIFID_SIDE^3) == the 27-symbol alphabet
#define TRIFID_MAX_SIDE 4       // largest cube side the side-generic primitive supports
#define TRIFID_MAX_CELLS 64     // TRIFID_MAX_SIDE^3 (headroom for primitive stress tests)
#define TRIFID_EXTRA_CHAR '+'   // the 27th symbol completing A..Z into a 27-cell cube

// Digrafid: a digraphic fractionation cipher over TWO independently keyed 27-symbol
// alphabets (A..Z + '#'). The horizontal alphabet is a 3x9 grid, the vertical a 9x3 grid;
// each grid is a permutation of 0..26. The two grids are packed back-to-back in the key
// lane (H = key[0..26], V = key[27..53]), so the search state is DIGRAFID_STATE cells.
#define DIGRAFID_GRID       27   // cells per grid == the 27-symbol alphabet
#define DIGRAFID_STATE      54   // both grids packed back-to-back (2 * DIGRAFID_GRID)
#define DIGRAFID_HROWS      3    // horizontal alphabet: 3 rows x 9 cols (cell = hr*9 + hc)
#define DIGRAFID_HCOLS      9
#define DIGRAFID_VROWS      9    // vertical alphabet: 9 rows x 3 cols (cell = vr*3 + vc)
#define DIGRAFID_VCOLS      3
#define DIGRAFID_EXTRA_CHAR '#'  // the 27th symbol completing A..Z into a 27-cell alphabet

#define PHILLIPS_SIDE 5         // Phillips grid side (the classic 5x5)
#define PHILLIPS_GRID 25        // Phillips grid size (PHILLIPS_SIDE * PHILLIPS_SIDE)
#define PHILLIPS_MAX_SIDE 6     // largest side the side-generic primitive/tests support
#define PHILLIPS_MAX_GRID 36    // PHILLIPS_MAX_SIDE * PHILLIPS_MAX_SIDE
#define PHILLIPS_MAX_SQUARES 10 // 2*PHILLIPS_MAX_SIDE - 2 (squares derived per base)
// Phillips square-generation variants (which axis the derived squares permute).
enum { PHILLIPS_ROW, PHILLIPS_COL, PHILLIPS_ROWCOL };

#define ALPHABET_SIZE 26        // the classical 26-letter alphabet. This is BOTH the size
                                // of the polyalphabetic keyword lanes AND the hardcoded
                                // mod base of the Vigenere/Beaufort/Porta/Quagmire/Autokey
                                // primitives ((x + ALPHABET_SIZE) % ALPHABET_SIZE), so it
                                // must stay 26 -- do not repurpose it as a max-alphabet.
#define DEFAULT_ALPHABET_SIZE 26 // the standard full A..Z alphabet size (init_alphabet(NULL))
#define MAX_ALPHABET_SIZE 36    // largest RUNTIME alphabet g_alpha can take: the ADFGVX
                                // 6x6 square (A..Z + 0..9, 36 symbols); the Trifid 27-symbol
                                // cube (A..Z + '+') also fits. Sizes only the runtime
                                // alphabet maps (g_idx_to_char_arr, g_monograms) and the
                                // ciphertext-facing IoC scratch -- everything that may be
                                // indexed by a live symbol id in [0, g_alpha). All other
                                // (polyalphabetic) arrays stay ALPHABET_SIZE.
#define MAX_CIPHER_LENGTH 10000

// Progressive Key per-letter base shift math (mod ALPHABET_SIZE). Shared by the
// primitive (progkey.c) and the solver's per-column monogram warm start.
static inline int progkey_base_encrypt(int p, int k, int base) {
    if (base == PROGKEY_BASE_VAR)  return (p - k + ALPHABET_SIZE) % ALPHABET_SIZE;
    if (base == PROGKEY_BASE_BEAU) return (k - p + ALPHABET_SIZE) % ALPHABET_SIZE;
    return (p + k) % ALPHABET_SIZE;                                   // VIGENERE
}
static inline int progkey_base_decrypt(int c, int k, int base) {
    if (base == PROGKEY_BASE_VAR)  return (c + k) % ALPHABET_SIZE;
    if (base == PROGKEY_BASE_BEAU) return (k - c + ALPHABET_SIZE) % ALPHABET_SIZE;
    return (c - k + ALPHABET_SIZE) % ALPHABET_SIZE;                   // VIGENERE
}

#define MAX_SYMBOLS    512      // distinct ciphertext symbols a homophonic cipher may use
#define MAX_TOKEN_LEN  16       // longest surface form of one ciphertext symbol token
#define MAX_FILENAME_LEN 4096   // must hold absolute paths; main() strcpy's CLI args in unbounded
#define MAX_KEYWORD_LEN 26
#define MAX_CYCLEWORD_LEN 300
#define MAX_NGRAM_SIZE 8
#define MAX_DICT_WORD_LEN 30
#define MAX_COLS 100            // upper bound on columnar column count (sizes order arrays)
#define MAX_TRANS_KEY 1024      // upper bound on a climbed transposition key length
                                // (grille orbit map, packed Cadenus order+rotation, etc.)

// Columnar read direction (cfg->read_direction).
#define COL_READ_TB   0        // read each column top-to-bottom (canonical)
#define COL_READ_BT   1        // read each column bottom-to-top
#define COL_READ_BOTH 2        // search both directions

// Columnar ROW read direction (cfg->read_row_direction), used by the track-perm
// columnar / chain solvers (transcol-L, transroutecol). Default LR is the
// canonical row read, so every other solver is unaffected.
#define ROW_READ_LR   0        // read each grid row left-to-right (canonical)
#define ROW_READ_RL   1        // read each grid row right-to-left
#define ROW_READ_BOTH 2        // search both row directions

// Optimization method (cfg->method): which acceptance strategy the engine uses.
// METHOD_DEFAULT keeps each cipher model's built-in SearchShape; the others force
// one strategy on EVERY cipher type (cipher-agnostic), set via -method.
#define METHOD_DEFAULT 0       // per-model default (shotgun, or anneal for transposition)
#define METHOD_SHOTGUN 1       // accept-worse with slip_probability
#define METHOD_ANNEAL  2       // geometric-cooling Metropolis (simulated annealing)
#define METHOD_PSO     3       // memetic particle swarm (discrete, cipher-agnostic)

#define FREQUENCY_WEIGHTED_SELECTION 1
#define DICTIONARY 1
#define INACTIVE -9999

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

// --- Tokenized ciphertext I/O ---
//
// A SymbolTable interns the distinct surface tokens of a ciphertext so that ciphers
// whose ciphertext alphabet is larger than (or disjoint from) A..Z -- homophonic
// substitution, in particular -- can be entered and displayed consistently. The
// ciphertext is decoded into a sequence of integer symbol ids (0..n-1) indexing this
// table. The single-character / 0..25 letter path (every other cipher type) does NOT
// build a table: decode_cipher reproduces ord() byte-for-byte when tab is unused.
typedef struct {
    int  n;                                    // distinct symbols interned so far
    char tokens[MAX_SYMBOLS][MAX_TOKEN_LEN];   // surface form of each symbol id
    int  freq[MAX_SYMBOLS];                    // occurrence count of each symbol id
    char delimiter;                            // 0 => per-character; else the field separator
} SymbolTable;

// --- Configuration Structs ---

typedef struct {
    int cipher_type;
    int ngram_size;
    int n_hill_climbs;
    int n_restarts;
    
    // Keyword/Cycleword constraints
    int plaintext_keyword_len;
    int ciphertext_keyword_len;
    int min_keyword_len;
    int plaintext_max_keyword_len;
    int ciphertext_max_keyword_len;
    int max_cycleword_len;
    int cycleword_len;

    // Bifid period (fractionation block size). period_present pins a single period;
    // otherwise the solver estimates and anneals the top n_periods candidates in
    // [2 .. max_period] (see bifid_solver.c).
    int period;
    bool period_present;
    int max_period;
    int n_periods;
    int n_primers;      // Gromark: top-K primers the pre-pass keeps to anneal (0 => auto by length)

    // Progressive Key: pin the per-group progression index (else the solver sweeps 0..25).
    bool progression_present;
    int  progression;

    // Interrupted Key: pin the interruptor letter / strategy (else the solver enumerates them);
    // -breaks supplies known group-start positions (random / word-division scheme, decoded in
    // solve_intkey from a whitespace-separated file of 0-based positions).
    bool interruptor_present;
    int  interruptor;             // 0..25
    bool intscheme_present;
    int  intscheme;               // IK_STRAT_*
    bool breaks_present;
    char breaks_file[MAX_FILENAME_LEN];

    // Condi: pin the starter offset (else the solver enumerates the 26 starter values 0..25).
    bool startkey_present;
    int  startkey;                // 0..25

    // Input Flags for lengths
    bool plaintext_keyword_len_present;
    bool ciphertext_keyword_len_present;
    bool cycleword_len_present;

    // Explicit User Keywords (Strings)
    char user_plaintext_keyword[ALPHABET_SIZE + 1];
    char user_ciphertext_keyword[ALPHABET_SIZE + 1];
    bool user_plaintext_keyword_present;
    bool user_ciphertext_keyword_present;

    // Probabilities & Thresholds
    double n_sigma_threshold;
    double ioc_threshold;
    double backtracking_probability;
    double keyword_permutation_probability;
    double slip_probability;

    // Weights
    float weight_ngram;
    float weight_crib;
    float weight_ioc;
    float weight_entropy;
    float weight_structure; // general transposition: reward regular (columnar) key steps
    float weight_monogram;  // homophonic: penalise decrypted letter-frequency deviation
                            // from English (chi-squared) -- stops the map collapsing
                            // many symbols onto a few common letters

    // Files
    char ciphertext_file[MAX_FILENAME_LEN];
    char batch_file[MAX_FILENAME_LEN];
    char crib_file[MAX_FILENAME_LEN];
    char dictionary_file[MAX_FILENAME_LEN];
    char ngram_file[MAX_FILENAME_LEN];

    // Flags
    bool verbose;
    bool skip_spaces;   // strip spaces/punctuation from the ciphertext entirely
    bool multiline;     // read the entire cipher file, not just its first line
    bool cipher_present;
    bool batch_present;
    bool crib_present;
    bool dictionary_present;
    bool variant;
    bool beaufort;

    bool optimal_cycleword;
    bool same_key_cycle;

    int method;          // METHOD_DEFAULT / METHOD_SHOTGUN / METHOD_ANNEAL (-method)

    // Simulated-annealing schedule (SHAPE_ANNEAL only). Backtracking is shared with
    // the shotgun path via backtracking_probability above. cooling_rate <= 0 means
    // "derive a geometric init_temp -> min_temp schedule over n_hill_climbs steps".
    double init_temp;    // starting temperature (-inittemp), default 0.10
    double min_temp;     // floor temperature for the derived schedule (-mintemp), default 0.001
    double cooling_rate; // per-iteration multiplier (-coolingrate); <= 0 => derive, default 0

    // Particle swarm schedule (SHAPE_PSO only, via -method pso). The swarm runs
    // n_particles particles for n_hill_climbs iterations, relaunched n_restarts
    // times; each iteration pulls a particle toward its personal best (cognitive)
    // and the global best (social) with random momentum (inertia), then refines it
    // with refine_steps local hill-climb steps (memetic).
    int    n_particles;  // swarm size (-nparticles), default 30
    double inertia;      // momentum: random-perturb scale (-inertia), default 0.7
    double cognitive;    // pull-toward-pbest scale (-cognitive), default 1.5
    double social;       // pull-toward-gbest scale (-social), default 1.5
    int    refine_steps; // per-particle local hill-climb steps (-refine), default 50

    // Transpositions.
    bool transperoffset_present;
    int trans_offset;
    int trans_period;

    bool transmatrix_present;
    int trans_w1;
    int trans_w2;
    int trans_clockwise; // 1 for clockwise, 0 for anti-clockwise

    // Columnar transposition (TRANSCOL / TRANSCOL2).
    int min_cols;        // smallest column count to search
    int max_cols;        // largest column count to search
    int read_direction;  // COL_READ_TB / COL_READ_BT / COL_READ_BOTH
    int read_row_direction; // ROW_READ_LR / ROW_READ_RL / ROW_READ_BOTH (transcol-L / chain)
    bool crib_anchored;  // -cribanchored: pin crib-implied positions in the transposition search
    float weight_word;   // optional dictionary word-fraction reward (default 0 => unused)
    int tile_h, tile_w;  // sub-grid tile shape for TRANSTILE (default 2x2)
    int trans_depth;     // PERIOD_COLUMN: max composed stages to search (default 2)
    int max_gaps;        // PERIOD_COLUMN_SPACE: max blank/gap cells to INSERT (dropped-char repair; default 4)
    int max_dels;        // PERIOD_COLUMN_SPACE: max observed cells to DELETE (added-char repair; default 2)

    // Nicodemus block height (rows per block). block_height > 0 pins a single height;
    // 0 => sweep [2 .. max_block_height]. max_block_height 0 => default top of the sweep.
    int block_height;
    int max_block_height;

    // Tokenized I/O. delimiter == 0 keeps the historical per-character / 0..25 letter
    // decode (bit-identical to ord()); a non-zero delimiter (default ',' for
    // HOMOPHONIC) splits the ciphertext into multi-character symbol tokens.
    char delimiter;
    bool delimiter_present;

} ColossusConfig;

typedef struct {
    float *ngram_data;
    char **dict;
    int n_dict_words;
    int max_dict_word_len;
} SharedData;

// Outcome of a polyalphabetic solve. solve_cipher fills this (when a non-NULL
// pointer is passed) so callers -- in particular the regression tests -- can
// inspect the recovered solution instead of scraping stdout. solved is false if
// no valid configuration was found, or for the transposition cipher types
// (which are dispatched to their own solvers and report separately).
typedef struct {
    bool solved;
    int cipher_type;
    double score;
    int n_words;
    int cycleword_len;
    int plaintext_keyword[ALPHABET_SIZE];
    int ciphertext_keyword[ALPHABET_SIZE];
    int cycleword[MAX_CYCLEWORD_LEN];
    int progression;            // Progressive Key: recovered progression index (else 0)
    int interruptor;            // Interrupted Key: recovered interruptor letter 0..25 (else -1)
    int intscheme;              // Interrupted Key: recovered strategy IK_STRAT_* (else -1)
    int condi_start;            // Condi: recovered starter offset 0..25 (else 0)
    int decrypted[MAX_CIPHER_LENGTH];
    int decrypted_len;
} SolveResult;

// =====================================================================
//  Cipher-type-agnostic search engine (run_solver) -- see CLAUDE.md.
// =====================================================================
//
// Every cipher is solved by the same skeleton: enumerate the outer search
// configurations, then for each either evaluate a single candidate (a pure
// parameter sweep) or shotgun-restart + hill-climb a candidate state, scoring
// each decrypt with the shared n-gram (+ crib) state_score and keeping the best.
// The skeleton lives in run_solver()/run_one_config(); each cipher type supplies
// a CipherModel (a vtable of hooks) describing how to seed/perturb/decrypt/report
// its own state. The engine owns the loops, acceptance, backtracking, best
// tracking, and the single state_score call site; the model owns the cipher math.

#define MAX_SOLVER_CONFIGS 65536   // cap on enumerated outer configs per solve

// Invariant problem instance + engine-owned scratch, passed to every hook.
typedef struct {
    ColossusConfig *cfg;
    SharedData     *shared;
    int            *cipher;          // cipher_indices[0..cipher_len-1]
    int             cipher_len;
    int            *crib_indices;
    int            *crib_positions;
    int             n_cribs;
    char           *cribtext;        // raw crib string ('_' = no crib) for the report hooks
    float          *ngram_data;      // == shared->ngram_data, hoisted for the hot path
    int            *hist_by_col;     // engine scratch: optimal-cycleword per-column histogram
                                     // (hist_by_col[col*ALPHABET_SIZE + c]); built when model->needs_hist
    void           *model_scratch;   // model-private per-config cache (e.g. indep_periodic seed)
    SolveResult    *result;          // polyalpha: report hook fills it (may be NULL)
} SolverCtx;

// One outer enumeration point. Field meanings are per-model: `period` is the
// cycleword length / column count K / rail count / permutation length; `j`,`k`
// are the polyalphabetic pt/ct keyword lengths; `aux` holds any remaining fixed
// parameters (amsco start, columnar read direction, route id, rail offset, ...).
typedef struct {
    int period;
    int j, k;
    int aux[2];
} SolverConfig;

// Superset candidate state; a model uses only the lane(s) it needs and its
// copy_state hook copies only those, so the per-iteration copy stays cheap.
typedef struct {
    int pt_keyword[ALPHABET_SIZE];   // polyalpha lane
    int ct_keyword[ALPHABET_SIZE];
    int cycleword[MAX_CYCLEWORD_LEN];
    int key[MAX_CIPHER_LENGTH];       // transposition lane (perm / order / rank / rot / square / turns)
    int key_len;
    int aux[8];                       // small per-state scratch (stage K's, dirs, carried period, ...)
} SolverState;

typedef enum {
    SHAPE_SHOTGUN,        // accept-worse with slip_probability (polyalpha, transmatrix)
    SHAPE_ANNEAL,         // geometric-cooling Metropolis acceptance (columnar, permutation)
    SHAPE_DETERMINISTIC,  // shotgun + first-improvement break (Vig/Beau/Porta under -optimalcycle)
    SHAPE_PSO             // memetic particle swarm; engine-only, forced via -method pso
} SearchShape;

// Live engine counters, filled by the engine for the verbose-report hook.
typedef struct {
    int    n_iterations;
    int    n_restarts;
    int    n_backtracks;
    int    n_slips;
    int    n_contradictions;
    clock_t start_time;
} EngineStats;

// The per-cipher-type model. Optional hooks may be NULL.
typedef struct CipherModel {
    const char *name;
    SearchShape shape;
    bool        needs_hist;          // engine builds ctx->hist_by_col before each config

    // Fill out[] with up to cap outer configs; return the count (<= cap).
    int  (*enumerate_configs)(const SolverCtx *ctx, SolverConfig *out, int cap);

    // Return the climbed key length for this config; 0 => this config is a single
    // SWEEP candidate (seed+decrypt+score once, no restarts). NULL => always climb.
    int  (*key_len)(const SolverCtx *ctx, const SolverConfig *cfg_c);

    // Produce a fresh restart state (random, or deterministic for a SWEEP cell).
    void (*seed)(const SolverCtx *ctx, const SolverConfig *cfg_c, SolverState *st);

    // One neighbour move on `st` (already a copy of current). *force_primary is the
    // cross-iteration "must perturb the primary lane" flag (polyalpha's
    // perturbate_keyword_p): the engine sets it true at the start of each restart;
    // the model reads it and writes back the value for the next iteration. Models
    // with no primary/secondary distinction ignore it. Any deterministic refine
    // (optimal cycleword), crib constrain, and -samekey coupling that the original
    // per-cipher climber did around the move live inside this hook.
    void (*perturb)(const SolverCtx *ctx, const SolverConfig *cfg_c, SolverState *st,
                    bool *force_primary);

    // Copy only the live lane(s) of `src` into `dst`.
    void (*copy_state)(const SolverConfig *cfg_c, const SolverState *src, SolverState *dst);

    // Decrypt `st` into out[0..cipher_len-1]. *score_adjust (init 0 by the engine)
    // receives any additive score term the model wants folded into state_score
    // (e.g. the general-transposition structure score); st may be mutated to carry
    // state forward (e.g. the detected period). Most models leave *score_adjust 0.
    void (*decrypt)(const SolverCtx *ctx, const SolverConfig *cfg_c, SolverState *st,
                    int *out, double *score_adjust);

    // Final report (human block + ">>>" CSV) for the global-best state.
    void (*report)(const SolverCtx *ctx, const SolverConfig *cfg_c, const SolverState *st,
                   double score, int *decrypted);
    // Optional: live display on each best-improvement when -verbose.
    void (*report_verbose)(const SolverCtx *ctx, const SolverConfig *cfg_c, const SolverState *st,
                           double score, int *decrypted, const EngineStats *stats);

    // --- Optional incremental fast path (all three NULL => generic path) ---
    // When a model supplies these three hooks the engine drives an incremental
    // inner loop instead of the generic copy/perturb/decrypt/full-rescore one: it
    // keeps the CURRENT (accepted) state's decryption live and scores each
    // neighbour as a cheap delta over only the cipher positions the move touched,
    // avoiding the per-iteration O(cipher_len) re-decrypt + n-gram rescan. The
    // returned score must equal what decrypt()+state_score()+score_adjust would
    // produce for the neighbour (verified by the regression suite). Used by the
    // homophonic model; every other type leaves these NULL and is byte-for-byte
    // unaffected.
    //
    // sync_caches: rebuild the model's incremental caches from `dec` (the current
    //   state's full decryption). Called once per restart/backtrack after a full
    //   decrypt resets the current state.
    void (*sync_caches)(const SolverCtx *ctx, const SolverConfig *cfg_c, const int *dec);
    // score_neighbor: score `loc` (= `cur` after one perturb) incrementally from
    //   the caches (synced to `cur`, whose decryption is `cur_dec`, score
    //   `cur_score`). Stashes the pending delta for a possible commit; must not
    //   mutate the caches or cur_dec.
    double (*score_neighbor)(const SolverCtx *ctx, const SolverConfig *cfg_c,
                             const SolverState *cur, const SolverState *loc,
                             const int *cur_dec, double cur_score);
    // commit_neighbor: apply the delta stashed by the most recent score_neighbor
    //   to the caches and to cur_dec, advancing the current state's decryption to
    //   the just-accepted neighbour. Called only when the engine accepts the move.
    void (*commit_neighbor)(const SolverCtx *ctx, const SolverConfig *cfg_c, int *cur_dec);
} CipherModel;


// --- Statistics Data ---

static int n_english_word_length_frequency_letters = 25;
static double english_word_length_frequencies[] = {
	0.0316, 0.16975, 0.21192, 0.15678, 0.10852, 0.08524, 0.07724, 
	0.05623, 0.04032, 0.02766, 0.01582, 0.00917, 0.00483, 0.00262, 
	0.00099, 0.0005, 0.00027, 0.00022, 0.00011, 0.00006, 0.00005, 
	0.00002, 0.00001, 0.00001, 0.00001};

static double english_monograms[] = {
	0.085517, 0.016048, 0.031644, 0.038712, 0.120965, 0.021815, 
	0.020863, 0.049557, 0.073251, 0.002198, 0.008087, 0.042065, 
	0.025263, 0.071722, 0.074673, 0.020662, 0.001040, 0.063327, 
	0.067282, 0.089381, 0.026816, 0.010593, 0.018254, 0.001914, 
	0.017214, 0.001138
};

// Per-cipher-type tuned search defaults. The global init_config() values suit the
// polyalphabetic / transposition score scale; substitution types with a very
// different scale (Playfair's mean log-probability, in particular) need their own
// schedule. Each entry carries a profile for BOTH search shapes, so -method
// shotgun|anneal each get sane defaults on that type. Precedence is:
//   init_config globals  <  this registry  <  explicit CLI flags
// (the registry only fills fields the user did not set). Types with no entry keep
// the global defaults exactly, so existing solves stay bit-identical.
typedef struct {
    int         cipher_type;     // which type this profile is for
    SearchShape default_shape;   // the model's own shape (which profile applies when
                                 // -method is not given)
    // Simulated-annealing profile (used when the effective shape is SHAPE_ANNEAL).
    int    a_n_restarts;
    int    a_n_hill_climbs;
    double a_init_temp;
    double a_min_temp;
    double a_cooling_rate;       // <= 0 => derive the geometric schedule
    double a_backtracking_probability;
    // Shotgun profile (used when the effective shape is SHAPE_SHOTGUN).
    int    s_n_restarts;
    int    s_n_hill_climbs;
    double s_slip_probability;
    double s_backtracking_probability;
    // Particle-swarm profile (used when the effective shape is SHAPE_PSO, i.e.
    // -method pso). p_n_particles == 0 means "no PSO profile for this type" -- the
    // global init_config defaults are kept, exactly like a type with no entry.
    int    p_n_restarts;
    int    p_n_hill_climbs;
    int    p_n_particles;
    double p_inertia;
    double p_cognitive;
    double p_social;
    int    p_refine_steps;
} SearchDefaults;


// Core Logic
void init_config(ColossusConfig *cfg);
// result may be NULL (CLI use); when supplied it receives the recovered solution.
void solve_cipher(char *ciphertext_str, char *cribtext_str, ColossusConfig *cfg,
    SharedData *shared, SolveResult *result);


// Porta cipher
void porta_decrypt(int output[], int input[], int len, int cycleword_indices[], int cycleword_len);
void porta_encrypt(int output[], int input[], int len, int cycleword_indices[], int cycleword_len);

// Portax cipher (periodic digraphic Porta). cycleword_indices are key letters 0..25.
void portax_encrypt(int output[], int input[], int len, int cycleword_indices[], int cycleword_len);
void portax_decrypt(int output[], int input[], int len, int cycleword_indices[], int cycleword_len);

// Vigenere cipher
void vigenere_decrypt(int decrypted[], int cipher_indices[], int cipher_len,
    int cycleword_indices[], int cycleword_len, bool variant);
void vigenere_encrypt(int encrypted[], int plaintext_indices[], int cipher_len,
    int cycleword_indices[], int cycleword_len, bool variant);

// Gronsfeld cipher (Vigenere with a numeric key; key_digits[] are per-column shifts 0..9)
void gronsfeld_decrypt(int decrypted[], int cipher_indices[], int cipher_len,
    int key_digits[], int key_len);
void gronsfeld_encrypt(int encrypted[], int plaintext_indices[], int plaintext_len,
    int key_digits[], int key_len);

// Progressive Key cipher (progkey.c). A periodic base cipher (Vig/Var/Beau, `base`) under a
// P-letter keyword (per-column shifts 0..25, in `keyword[]`), composed with a per-GROUP
// constant drift: group g = i/P adds (g*prog) mod 26 via a second base pass. So the per-
// position shift is keyword[i%P] (+/-) (g*prog). `progkey_deprogress` undoes only the drift
// pass, leaving a pure periodic base cipher (C1) the solver attacks per column.
void progkey_encrypt(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int prog, int base);
void progkey_decrypt(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int prog, int base);
void progkey_deprogress(int out[], int cipher_indices[], int cipher_len,
    int P, int prog, int base);

// Interrupted Key cipher (intkey.c). A periodic base cipher (Vig/Var/Beau, `base`) under a
// P-letter keyword (per-column shifts 0..25, in `keyword[]`) whose key index k starts at 0,
// increments (k+1)%P each position, and RESETS to 0 at each break point. The *-mask forms take
// an is_break[] array (is_break[i]==1 => position i starts a new group; is_break[0] implicit);
// the *-ptint forms reset causally after the plaintext letter equals `interruptor`.
void intkey_encrypt_mask(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, const int is_break[]);
void intkey_decrypt_mask(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, const int is_break[]);
void intkey_encrypt_ptint(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, int interruptor);
void intkey_decrypt_ptint(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, int interruptor);
// Ciphertext-interruptor: reset causally AFTER the produced/consumed ciphertext letter equals
// `interruptor`. Encryption must be causal (the reset depends on the ciphertext being produced);
// decryption is equivalent to intkey_decrypt_mask over the mask built by intkey_build_mask_ct,
// which is what the solver uses (build the mask once, then anneal the keyword).
void intkey_encrypt_ctint(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, int interruptor);
void intkey_decrypt_ctint(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, int interruptor);
// Build the ciphertext-interruptor break mask: is_break[i]==1 iff cipher[i-1]==interruptor
// (reset AFTER the interruptor CT letter). Keyword-independent, so precomputable by the solver.
void intkey_build_mask_ct(int is_break[], const int cipher_indices[], int cipher_len, int interruptor);

// Condi cipher (condi.c). A plaintext-feedback substitution over a keyed alphabet sigma (a
// 26-permutation; sigma[k] = letter at position k, sigma_inv[x] = position of letter x, supplied by
// the caller). The shift for each letter is the 1-indexed position of the PREVIOUS plaintext letter
// in sigma; the first letter uses `starter` (0..25): ct_i = sigma[(idx(pt_i) + off_i) mod 26] with
// off_1 = starter, off_{i+1} = (idx(pt_i) + 1) mod 26. Decryption is causal.
void condi_encrypt(const int plaintext_indices[], int plaintext_len,
    const int sigma[], const int sigma_inv[], int starter, int encrypted[]);
void condi_decrypt(const int cipher_indices[], int cipher_len,
    const int sigma[], const int sigma_inv[], int starter, int decrypted[]);

// Gromark / Periodic Gromark cipher (gromark.c). A keyed 26-letter substitution sigma (a
// permutation of A..Z) composed with a chain-addition running key from a P-digit primer
// (d[i] = primer[i] for i<P, else (d[i-P] + d[i-P+1]) mod 10; one digit per letter):
//   BASIC     C[i] = sigma[(p[i] + d[i]) mod 26]
//   PERIODIC  C[i] = sigma[(p[i] + d[i] + offset[(i/P) mod P]) mod 26], P = period, one
//             offset per group of P letters (offset[g] = position of keyword[g] in sigma).
// gromark_chain_key generates the running key; gromark_mixed_alphabet builds sigma from a
// keyword via the K2M transposition block (generator + unit tests only -- the solver
// hill-climbs sigma as a free permutation). gromark_decrypt_core is the solver hot path
// (precomputed running key d[] + inverse alphabet; offsets == NULL => basic).
void gromark_chain_key(const int primer[], int primer_len, int n, int out_digits[]);
void gromark_mixed_alphabet(const char *keyword, int sigma[]);
void gromark_build_from_keyword_idx(const int kw[], int P, int sigma[],
                                    int primer[], int offsets[]);
void gromark_encrypt(const int plain[], int len, const int sigma[],
                     const int primer[], int primer_len, int out[]);
void gromark_decrypt(const int cipher[], int len, const int sigma[],
                     const int primer[], int primer_len, int out[]);
void gromark_periodic_encrypt(const int plain[], int len, const int sigma[],
                     const int primer[], int period, const int offsets[], int out[]);
void gromark_periodic_decrypt(const int cipher[], int len, const int sigma[],
                     const int primer[], int period, const int offsets[], int out[]);
void gromark_decrypt_core(const int cipher[], int len, const int sigma_inv[],
                     const int d[], const int offsets[], int period, int out[]);

// Beaufort cipher
void beaufort_decrypt(int decrypted[], int cipher_indices[], int cipher_len, 
    int cycleword_indices[], int cycleword_len);
void beaufort_encrypt(int encrypted[], int plaintext_indices[], int cipher_len, 
    int cycleword_indices[], int cycleword_len);

// Quagmire I - IV ciphers
void quagmire_decrypt(int decrypted[], int cipher_indices[], int cipher_len, 
    int plaintext_keyword_indices[], int ciphertext_keyword_indices[], 
    int cycleword_indices[], int cycleword_len, bool variant);
void quagmire_encrypt(int encrypted[], int plaintext_indices[], int cipher_len, 
    int plaintext_keyword_indices[], int ciphertext_keyword_indices[], 
    int cycleword_indices[], int cycleword_len, bool variant);

// Autokey
void autokey_encrypt(ColossusConfig *cfg, int ciphertext[], int plaintext_indices[],
    int plaintext_len, int plaintext_keyword[], int ciphertext_keyword[],
    int primer_indices[], int primer_len);
void autokey_decrypt(ColossusConfig *cfg, int decrypted[], int cipher_indices[],
    int cipher_len, int plaintext_keyword[], int ciphertext_keyword[],
    int key_indices[], int key_len);

// Transposition
void transperoffset(int plaintext[], int len, int d, int n);
void matrix_rotate(int text[], int len, int width, int clockwise);
void transmatrix(int text[], int len, int w1, int w2, int clockwise);

// Columnar transposition primitive (one stage). Inverts a columnar encryption:
// the ciphertext is `K` columns concatenated in read order `order[0..K-1]`, each
// column read top-to-bottom (dir == COL_READ_TB) or bottom-to-top (COL_READ_BT).
// Incomplete grids (len % K != 0) are handled via per-column heights. Writes the
// recovered row-major plaintext to out[0..len-1] (out must differ from cipher).
void decrypt_columnar(int cipher[], int len, int K, int order[], int dir, int out[]);

// Columnar with a uniform within-column row permutation L (the "track order"):
// like decrypt_columnar, but the grid's R rows are read in the order L[0..R-1]
// (output row i = grid row L[i]), each left-to-right (row_rev 0) or reversed
// (row_rev 1). Needs a complete grid (len % K == 0, R = len/K); a ragged grid
// falls back to decrypt_columnar. out[] must not alias cipher[].
void decrypt_columnar_tracked(int cipher[], int len, int K, int order[], int dir,
                              int L[], int row_rev, int out[]);

// Sub-grid / tile transposition primitive: the W contiguous R-tall columns are
// reordered by `order`, then every complete h x w tile of the resulting R x W grid
// is inverted by the cell permutation `perm` (length h*w), and the grid is read
// row-major. Complete grid only (len % W == 0). out[] must not alias cipher[].
void decrypt_tile(int cipher[], int len, int W, int order[], int h, int w,
                  int perm[], int out[]);

// Period column order transposition primitive (one stage; AZdecrypt's "Period column
// order"). The input is laid row-major into a `dx`-wide grid of `dy = ceil(len/dx)`
// rows; its columns are reordered by a periodic schedule of period `p` (visit columns
// 0, p, 2p, ..., then 1, 1+p, ..., i.e. `for i in 0..p-1: for j = i; j < dx; j += p`),
// then the grid is read back row-major, skipping the (dx*dy - len) padding cells. `utp`
// selects the transform (0) or its exact inverse (1). Values are carried opaquely
// (letters 0..g_alpha-1 AND negative space/punctuation sentinels ride through), so it
// composes with the transposition pipeline. Length-preserving; `out` must not alias
// `in`. A complete grid (len % dx == 0) makes utp==1 the bit-exact inverse of utp==0.
void period_column_transform(const int *in, int *out, int len, int dx, int p, int utp);

// Number of geometric routes recognised by route_cells()/decrypt_route().
#define N_ROUTES 6

// Rail fence (and variant rail fence) primitive: invert a zigzag over `rails`
// rows with starting phase `offset` in [0, 2*(rails-1)). `variant` swaps the read
// and write directions (the -variant convention). out[] must not alias cipher[].
void decrypt_railfence(int cipher[], int len, int rails, int offset, int variant, int out[]);

// Route transposition primitives. route_cells() fills cells[] with the row-major
// cell indices in reading order for route_id in [0, N_ROUTES), emitting only cells
// that exist (row-major index < len, so a ragged final row is handled), and returns
// the number of cells emitted. decrypt_route() inverts one route over an R x C grid
// (complete or ragged: (R-1)*C < len <= R*C); `variant` swaps read/write directions.
int route_cells(int R, int C, int len, int route_id, int cells[]);
void decrypt_route(int cipher[], int len, int R, int C, int route_id, int variant, int out[]);

// Amsco primitive: invert an alternating 1/2-letter columnar transposition over K
// columns read in `order`, with `start` (1 or 2) the size of the first cell.
// `variant` swaps read/write. out[] must not alias cipher[].
void decrypt_amsco(int cipher[], int len, int K, int order[], int start, int variant, int out[]);

// Myszkowski primitive: invert a columnar transposition over K columns whose
// keyword ranks `rank[0..K-1]` may tie; tied columns are read row-by-row together.
// `variant` swaps read/write. out[] must not alias cipher[].
void decrypt_myszkowski(int cipher[], int len, int K, int rank[], int variant, int out[]);

// Redefence primitive: rail fence over `rails` rows (phase `offset`) whose rails
// are read in keyed order `order[0..rails-1]`. `variant` swaps read/write.
void decrypt_redefence(int cipher[], int len, int rails, int offset, int order[], int variant, int out[]);

// Cadenus primitive: K = len/rows columns rotated vertically by rot[c] and reordered
// by order[p] (read col p -> original col order[p]), grid read row-major.
void decrypt_cadenus(int cipher[], int len, int K, int order[], int rot[], int variant, int out[]);

// Nihilist transposition primitive: independent rowperm/colperm on the N x N grid
// (N = sqrt(len)), read row-major (readmode 0) or column-major (readmode 1).
// `variant` swaps read/write.
void decrypt_nihilist(int cipher[], int len, int N, int rowperm[], int colperm[],
                      int readmode, int variant, int out[]);

// Swagman primitive: N x N key square (column j a permutation of 0..N-1) applied
// column-wise over an N x (len/N) grid; readmode selects row/column-major read-off.
void decrypt_swagman(int cipher[], int len, int N, int square[], int readmode, int variant, int out[]);

// Turning-grille primitive: key[orbit] in {0..3} chooses the turn that exposes each
// rotation orbit of the N x N grid (N = sqrt(len)). Writes the orbit count to
// *n_orbits (the climbed key length) when non-NULL. `variant` swaps read/write.
void decrypt_grille(int cipher[], int len, int N, int key[], int variant, int out[], int *n_orbits);


// Playfair cipher (playfair.c). Primitives operate on 0..g_alpha-1 alphabet indices,
// with the grid a permutation of the active 25-letter alphabet (g_alpha == 25):
// grid[p] is the letter at grid position p (row p/5, col p%5); pos[] is its inverse.
// The solver needs only playfair_decrypt(); encrypt/prepare/grid_from_keyword serve
// the test-data generator and the round-trip + known-answer unit tests.
void playfair_build_inverse(const int grid[], int pos[]);
void playfair_encrypt(const int plain[], int len, const int grid[], int out[]);
void playfair_decrypt(const int cipher[], int len, const int grid[], int out[]);
int  playfair_prepare(const int raw[], int len, int filler, int alt, int out[], int out_cap);
void playfair_grid_from_keyword(const int keyword[], int kwlen, int grid[]);

// Seriated Playfair cipher (seriated_playfair.c). Plain Playfair over the same 5x5 keyed
// grid, but digraphs are vertical pairs of a two-row seriated layout of period P: in each
// 2P block, pair j couples block-letter j with j+P. Reuses playfair_build_inverse and
// playfair_grid_from_keyword. The solver needs only seriated_playfair_decrypt(); encrypt
// + prepare (null insertion for doubled vertical pairs) serve the generator + unit tests.
void seriated_playfair_encrypt(const int plain[], int len, const int grid[], int period, int out[]);
void seriated_playfair_decrypt(const int cipher[], int len, const int grid[], int period, int out[]);
int  seriated_playfair_prepare(const int raw[], int len, int period, int filler, int alt, int out[], int out_cap);


// Bifid cipher (bifid.c). Side-generic over a side x side keyed Polybius square (a
// permutation of the active n = side*side letter alphabet carried in 0..n-1 indices):
// grid[p] is the letter at cell p (row p/side, col p%side); pos[] is its inverse. The
// fractionation works in blocks of `period`: encryption writes the block's row
// coordinates then its column coordinates as one 2*L stream, then re-pairs that stream
// consecutively into ciphertext cells; decryption is the inverse. An incomplete final
// block (L < period) is handled in place. The solver needs only bifid_decrypt();
// encrypt/grid_from_keyword serve the test-data generator and the unit tests.
void bifid_build_inverse(const int grid[], int pos[], int n);
void bifid_encrypt(const int plain[], int len, const int grid[], int side, int period, int out[]);
void bifid_decrypt(const int cipher[], int len, const int grid[], int side, int period, int out[]);
void bifid_grid_from_keyword(const int keyword[], int kwlen, int grid[], int n);

// CM Bifid (Conjugated Matrix Bifid) cipher (cm_bifid.c). Plain Bifid with TWO keyed
// squares: fractionate plaintext into coords via sq1, then map the re-paired coordinate
// stream to ciphertext letters via a different sq2 (decrypt mirrors: expand via sq2,
// recombine via sq1). sq1 == sq2 reduces to Bifid. Side-generic; reuses bifid_build_inverse.
void cm_bifid_encrypt(const int plain[], int len, const int sq1[], const int sq2[],
                      int side, int period, int out[]);
void cm_bifid_decrypt(const int cipher[], int len, const int sq1[], const int sq2[],
                      int side, int period, int out[]);


// Nihilist Substitution cipher (nihilist_sub.c). A periodic ADDITIVE cipher over a keyed
// side x side Polybius square (a permutation of the active n = side*side alphabet, carried
// in 0..n-1 indices like Bifid). Each plaintext letter -> its cell -> a 2-digit coordinate
// number rowlbl[row]*10 + collbl[col] (fixed labels are 1..side); a periodic key (its own
// cells -> numbers) is added per position under one of the NIH_ADD_* conventions, so the
// ciphertext is a stream of decimal NUMBERS (not letters). The primitive is side-generic and
// label-aware (rowlbl/collbl are a permutation of 1..side; nihilist_sub_fixed_labels fills the
// standard 1..side). The solver always assumes fixed labels and recovers any keyed-label
// cipher as the equivalent relabeled square. nihilist_sub_decrypt returns the count of
// positions that decrypted to a legal coordinate (the square-independent validity signal the
// solver rewards). key_cells[] are coordinate indices 0..n-1 (the additive's own cells).
void nihilist_sub_fixed_labels(int rowlbl[], int collbl[], int side);
int  nihilist_sub_num_valid(int num, int side);     // both digits of num in [1..side]?
void nihilist_sub_encrypt(const int plain[], int n, const int grid[],
        const int rowlbl[], const int collbl[], int side,
        const int key_cells[], int period, int conv, int out_nums[]);
int  nihilist_sub_decrypt(const int nums[], int n, const int grid[],
        const int rowlbl[], const int collbl[], int side,
        const int key_cells[], int period, int conv, int out_letters[]);


// ADFGVX / ADFGX cipher (adfgvx.c). Fractionation over a keyed side x side Polybius
// square (built/inverted with the shared bifid helpers) composed with a keyed columnar
// transposition. Works in COORDINATE space: a length-2N ciphertext stream of
// coordinates 0..side-1 (the solver maps the label characters A/D/F/G/(V/)X to
// coordinates up front). adfgvx_encrypt maps N plaintext symbols (0..side*side-1) to a
// 2N coordinate stream; adfgvx_decrypt inverts (undo columnar via decrypt_columnar,
// then pair consecutive coordinates and look up the square). `order[0..K-1]` is the
// columnar read order, `dir` the column read direction (COL_READ_TB / COL_READ_BT).
// adfgvx_labels(side) returns the label string ("ADFGX" for side 5, "ADFGVX" for 6).
extern const char ADFGX_LABELS[];
extern const char ADFGVX_LABELS[];
const char *adfgvx_labels(int side);
void adfgvx_encrypt(const int plain[], int n, const int square[], int side, int K,
                    const int order[], int dir, int out[]);
void adfgvx_decrypt(const int cipher[], int len2, const int square[], int side, int K,
                    const int order[], int dir, int out[]);


// Phillips cipher (phillips.c). A periodic monographic substitution over `nsq = 2*side-2`
// keyed Polybius squares derived from one base square (a permutation of the active
// n = side*side letter alphabet, carried in 0..n-1 indices: base[p] is the letter at cell
// p, row p/side col p%side). The plaintext is split into blocks of `side` letters and
// block b is enciphered with square (b mod nsq); each letter is replaced by the one
// diagonally down-right (with wrap), so the overall period is nsq*side (40 for the 5x5).
// The `variant` selects how the derived squares are built (PHILLIPS_ROW / _COL / _ROWCOL,
// see phillips.c). The solver needs only phillips_decrypt(); build_squares/encrypt/
// grid_from_keyword serve the test-data generator and the unit tests.
void phillips_build_squares(const int base[], int side, int variant, int squares[]);
void phillips_encrypt(const int plain[], int len, const int base[], int side, int variant, int out[]);
void phillips_decrypt(const int cipher[], int len, const int base[], int side, int variant, int out[]);
void phillips_grid_from_keyword(const int keyword[], int kwlen, int grid[], int n);


// Two-Square cipher (twosquare.c). A digraphic substitution over TWO keyed side x side
// Polybius squares (each a permutation of the active n = side*side letter alphabet in
// 0..n-1 indices: sq[p] is the letter at cell p, row p/side col p%side). For each
// plaintext digraph (a, b): a is located in square 1, b in square 2, and the cipher pair
// is the opposite corners of the rectangle they span. `variant` is the arrangement:
//   TWO_SQ_HORIZONTAL (ACA): squares side by side -- out = (sq2[r1][c2], sq1[r2][c1]);
//     a same-row digraph maps to the reversed pair (a "transparency").
//   TWO_SQ_VERTICAL (Wikipedia): squares stacked -- out = (sq1[r1][c2], sq2[r2][c1]),
//     which is self-inverse (decrypt == encrypt); a same-column digraph maps to itself.
// No doubled-letter handling and no padding beyond an even length (an odd trailing letter
// passes through). The solver needs only twosquare_decrypt(); encrypt + the shared
// playfair_grid_from_keyword serve the generator and the unit tests.
void twosquare_encrypt(const int plain[], int len, const int sq1[], const int sq2[],
                       int side, int variant, int out[]);
void twosquare_decrypt(const int cipher[], int len, const int sq1[], const int sq2[],
                       int side, int variant, int out[]);


// Four-Square cipher (foursquare.c). A digraphic substitution over a 2x2 layout of four
// side x side squares: the upper-left and lower-right are the FIXED standard square (cell
// p holds letter p), the upper-right (ur) and lower-left (ll) are the keyed unknowns (each
// a permutation of 0..n-1). For each plaintext digraph (p1, p2): p1 sits at (r1,c1) of the
// standard UL, p2 at (r2,c2) of the standard LR, and the cipher pair is (ur[r1][c2],
// ll[r2][c1]); decryption is the inverse (find c1 in ur, c2 in ll, read the standard
// corners). An odd trailing letter passes through. The solver needs only foursquare_decrypt();
// encrypt + standard_square serve the generator and the unit tests.
void foursquare_standard_square(int sq[], int n);
void foursquare_encrypt(const int plain[], int len, const int ur[], const int ll[],
                        int side, int out[]);
void foursquare_decrypt(const int cipher[], int len, const int ur[], const int ll[],
                        int side, int out[]);


// Tri-Square cipher (trisquare.c). A digraphic substitution over THREE independent keyed
// side x side squares (each a permutation of 0..n-1). A plaintext digraph (p1 in sq1 at
// (r1,c1), p2 in sq2 at (r2,c2)) becomes a ciphertext TRIGRAPH (a 3:2 expansion):
//   c0 = any letter in p1's COLUMN in sq1 (decodes to c1);  c1 = sq3[r1][c2] (deterministic);
//   c2 = any letter in p2's ROW in sq2 (decodes to r2). Encryption is polyphonic; decryption
// is exact for any clerk choice. trisquare_encrypt picks a canonical representative (top of
// column / left of row); a lone trailing letter passes through. The solver needs only
// trisquare_decrypt(); encrypt + the shared playfair_grid_from_keyword serve the generator
// and the unit tests. Side-generic; reuses bifid_build_inverse. Both return the length written.
int trisquare_encrypt(const int plain[], int len, const int sq1[], const int sq2[],
                      const int sq3[], int side, int out[]);
int trisquare_decrypt(const int cipher[], int len3, const int sq1[], const int sq2[],
                      const int sq3[], int side, int out[]);


// Trifid cipher (trifid.c). The 3D generalization of Bifid: side-generic over a
// side x side x side keyed cube (a permutation of the active n = side^3 letter
// alphabet carried in 0..n-1 indices). cube[p] is the letter at cell p, whose three
// coordinates are (c0, c1, c2) = (p/(side*side), (p/side)%side, p%side) -- "layer",
// "row", "column" -- each in 0..side-1; pos[] is the inverse. The fractionation works
// in blocks of `period`: encryption lays out the block's layer coordinates, then its
// row coordinates, then its column coordinates as one 3*L stream, then re-groups that
// stream into consecutive triples that index new cube cells; decryption is the inverse.
// An incomplete final block (L < period) is handled in place. The default cube is the
// 3x3x3, 27-symbol (A..Z + '+') cube. The solver needs only trifid_decrypt(); the rest
// serve the test-data generator and the unit tests.
void trifid_build_inverse(const int cube[], int pos[], int n);
void trifid_encrypt(const int plain[], int len, const int cube[], int side, int period, int out[]);
void trifid_decrypt(const int cipher[], int len, const int cube[], int side, int period, int out[]);
void trifid_cube_from_keyword(const int keyword[], int kwlen, int cube[], int n);


// Digrafid cipher (digrafid.c). A digraphic fractionation cipher over two keyed 27-symbol
// alphabets: the horizontal grid H (3x9, cell = hr*9 + hc) and the vertical grid V (9x3,
// cell = vr*3 + vc), each a permutation of 0..26. A plaintext digraph (a,b) -> a 3-digit
// number (top = a's column 0..8, bot = b's row 0..8, mid = a's row*3 + b's col, 0..8).
// The digraphs of each PERIOD-long group have their 3-digit numbers stacked as 3 rows,
// read row-major into a 3*g stream, re-split into consecutive triples, and mapped back
// through the tableau to the output digraphs (Trifid-style fractionation over digraphs).
// The solver needs only digrafid_decrypt(); encrypt/grid_from_keyword serve the test-data
// generator and the unit tests. digrafid_grid_from_keyword places the keyed sequence
// row-major (column_major == 0, for H) or column-major (!= 0, for V) into a rows x cols grid.
void digrafid_encrypt(const int plain[], int len, const int gridH[], const int gridV[], int period, int out[]);
void digrafid_decrypt(const int cipher[], int len, const int gridH[], const int gridV[], int period, int out[]);
void digrafid_grid_from_keyword(const int keyword[], int kwlen, int grid[], int rows, int cols, int column_major);


// Fractionated Morse (fracmorse.c). Plaintext -> Morse with single 'x' separators -> trigraphs
// of {DOT,DASH,X} (trigraph (a,b,c) has rank 9a+3b+c, and xxx never occurs) -> ciphertext via a
// keyed alphabet sigma (rank -> letter). fracmorse_encrypt returns the ciphertext length C (a
// length change from N plaintext letters). fracmorse_decrypt inverts sigma, expands to the fixed
// 3C-symbol Morse stream, splits on x, decodes each run: a legal codeword -> its letter (a valid
// token), else `filler` (an invalid token); it returns the letter count and reports the total /
// legal-codeword token counts via n_tokens/n_valid (the solver's Morse-validity reward). The
// solver needs only fracmorse_decrypt(); encrypt serves the generator + unit tests.
int fracmorse_encrypt(const int plain[], int n, const int sigma[], int out[]);
int fracmorse_decrypt(const int cipher[], int clen, const int sigma[], int out[],
                      int filler, int *n_tokens, int *n_valid);


// Hill cipher (hill.c). A polygraphic substitution: a block of k plaintext letters
// (a column vector p) is enciphered c = K*p mod 26 with a k x k key matrix K (row-major
// in mat[]); deciphering is p = K^-1*c mod 26, which exists iff gcd(det K, 26) == 1. The
// mod base is ALPHABET_SIZE (26), so -type hill keeps the full 26-letter alphabet. The
// one hot-path primitive is hill_mat_mul_blocks (the solver hill-climbs the DECRYPTION
// matrix and applies it straight to the ciphertext); a trailing partial block (len % k)
// is copied through unchanged. The determinant / inverse / keyword build serve the
// generator, the unit tests, and the report hook (inverting the recovered decryption
// matrix to display the encryption key).
void hill_mat_mul_blocks(const int mat[], int k, const int in[], int len, int out[]);
void hill_encrypt(const int plain[], int len, const int key[], int k, int out[]);
void hill_decrypt(const int cipher[], int len, const int key[], int k, int out[]);
int  hill_det_mod(const int mat[], int k);
int  hill_mod_inverse(int a, int m);
int  hill_mat_inverse(const int mat[], int k, int out[]);
void hill_matrix_from_keyword(const int keyword[], int kwlen, int mat[], int k);


// hist_by_col, when non-NULL, is a caller-supplied per-column ciphertext
// histogram laid out as hist_by_col[col*ALPHABET_SIZE + c] (counts of cipher
// char c in column col for this cycleword_len). It depends only on the (fixed)
// ciphertext and cycleword_len, so the caller computes it once per
// shotgun_hill_climber call instead of every derive. Pass NULL to have the
// routine build it locally (standalone/test use).
void derive_optimal_cycleword(
    ColossusConfig *cfg,
    int cipher_indices[], int cipher_len,
    int plaintext_keyword_indices[], int ciphertext_keyword_indices[],
    int cycleword_state[], int cycleword_len, int *hist_by_col);


void straight_alphabet(int keyword[], int len);
void make_keyed_alphabet(char *keyword_str, int *output_indices); // NEW
double entropy(int text[], int len);
double chi_squared(int plaintext[], int len);

// I/O & Data
void load_dictionary(char *filename, char ***dict, int *n_dict_words, int *max_dict_word_len, bool verbose);
void free_dictionary(char **dict, int n_dict_words);
int find_dictionary_words(char *plaintext, char **dict, int n_dict_words, int max_dict_word_len);

// Fast word-set membership + the length-weighted dictionary-coverage objective
// (Rec 2), used by the space-preserving transposition solvers. Built on demand from
// the loaded dictionary; default solves never touch it.
typedef struct WordSet WordSet;
WordSet *word_set_build(char **dict, int n_words);
void word_set_free(WordSet *ws);
int word_set_contains(const WordSet *ws, const char *upper_word);
double word_coverage(const int *text, int len, const WordSet *ws);


void shuffle(int *array, size_t n);

// Stats
double mean_ioc(int text[], int len, int len_cycleword, int *caesar_column);
void estimate_cycleword_lengths(int text[], int len, int max_cycleword_len, 
	double n_sigma_threshold, double ioc_threshold, 
	int *n_cycleword_lengths, int cycleword_lengths[], bool verbose);
double vec_mean(double vec[], int len);
double vec_stddev(double vec[], int len);

// Utils
static inline uint32_t fast_rand(void);
static inline void seed_fast_rand(uint32_t seed);
static inline uint32_t fast_rand_bounded(uint32_t range);
int gcd(int a, int b);
int str_eq(const char *a, const char *b);
int parse_cipher_type(const char *arg);
int unique_len(char *str);
void vec_print(int vec[], int len);
void print_text(int indices[], int len);
void ord(char *text, int indices[]);

// Decode a raw ciphertext string into integer indices, returning the index count.
// LETTER mode (cfg->cipher_type != HOMOPHONIC and no -delimiter): byte-for-byte ord()
// -- one index per character, A..Z -> 0..g_alpha-1, everything else a negative
// sentinel; tab is left untouched (may be NULL). SYMBOL mode (HOMOPHONIC, or any type
// with -delimiter): split into tokens and intern them into tab, writing symbol ids.
int decode_cipher(const char *text, const ColossusConfig *cfg, int indices[], SymbolTable *tab);

// Echo a decoded ciphertext. tab == NULL reproduces print_text() exactly (letter
// mode); otherwise the symbol tokens are printed joined by tab->delimiter.
void print_cipher(const int indices[], int len, const SymbolTable *tab);

// Non-alphabetic bytes (spaces, punctuation) in a ciphertext are carried through
// the integer-index arrays as negative values that reversibly encode the original
// byte: index = -(unsigned char)c - 1, so 'A'..'Z' stay 0..25 and everything else
// is < 0. Scoring and frequency stats skip these positions; printing restores the
// original character. This lets pure-transposition solves permute the spaces along
// with the letters and reveal word boundaries in the recovered plaintext.
// Runtime alphabet. Defaults to the full 26-letter A..Z (g_alpha == 26,
// g_idx_to_char_arr == "ABC..Z", g_char_to_idx the identity) so the historical
// behaviour is bit-identical. -excludeletter / -alphabet shrink it (e.g. the
// 25-letter A..Z-minus-P, mod 25) for ciphers built on a reduced alphabet.
// Array sizes everywhere stay ALPHABET_SIZE (26, the max); only loop bounds,
// modular arithmetic, and the n-gram packing base use g_alpha.
extern bool g_ngram_logprob;      // n-gram scoring mode (see utils.c); false = legacy
extern int  g_alpha;              // runtime alphabet size (<= ALPHABET_SIZE)
extern int  g_char_to_idx[128];   // ASCII (upper) -> alphabet index, or -1 if absent
extern char g_idx_to_char_arr[MAX_ALPHABET_SIZE + 1];  // alphabet index -> char (room for 27)
extern double g_monograms[MAX_ALPHABET_SIZE];          // English monogram freqs, reindexed to runtime alphabet
void init_alphabet(const char *excluded);          // (re)build the maps; NULL => full A..Z
// Build the 27-symbol Trifid alphabet: A..Z (0..25) plus TRIFID_EXTRA_CHAR ('+') at
// index 26, so a 3x3x3 cube has exactly 27 cells. Unlike init_alphabet this registers a
// non-letter in g_char_to_idx (so '+' decodes from the ciphertext) and gives it a
// negligible English monogram weight (the cube attack uses no monogram penalty).
void init_alphabet_trifid(void);
// Build the 36-symbol ADFGVX alphabet: A..Z (0..25) plus the digits 0..9 (26..35), so a
// 6x6 Polybius square has exactly 36 cells. Like init_alphabet_trifid this registers the
// non-letter digits in g_char_to_idx (so a digit-bearing plaintext decodes) and gives
// them a negligible English monogram weight (the square attack uses no monogram penalty,
// and the A..Z n-gram table never indexes a digit). ADFGX uses init_alphabet("J") (25).
void init_alphabet_adfgvx(void);
// Build the 27-symbol Digrafid alphabet: A..Z (0..25) plus DIGRAFID_EXTRA_CHAR ('#') at
// index 26, so each 3x9 / 9x3 grid has exactly 27 cells. Like init_alphabet_trifid this
// registers the non-letter in g_char_to_idx (so '#' decodes from the ciphertext) and gives
// it a negligible English monogram weight (the grid attack uses no monogram penalty).
void init_alphabet_digrafid(void);

static inline int index_to_char(int idx) {
    return (idx >= 0) ? (unsigned char) g_idx_to_char_arr[idx] : (-(idx + 1));
}

float index_of_coincidence(int plaintext[], int len);
void tally(int plaintext[], int len, int frequencies[], int n_frequencies);
bool file_exists(const char * filename);
void vec_copy(int src[], int dest[], int len);
int int_pow(int base, int exp);

extern uint32_t rng_state;

// Fast inline Xorshift32 generator
static inline uint32_t fast_rand(void) {
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return rng_state = x;
}

// Function to seed the PRNG
static inline void seed_rand(uint32_t seed) {
    if (seed == 0) seed = 1; // State cannot be 0
    rng_state = seed;
}

// Lemire's method to map the 32-bit random integer to a specific range [0, range)
static inline uint32_t rand_bounded(uint32_t range) {
    return (uint32_t)(((uint64_t)fast_rand() * range) >> 32);
}

static inline int rand_int(int min, int max) {
    if (min >= max) return min; 
    uint32_t range = (uint32_t)(max - min);
    return min + (int)rand_bounded(range);
}

static inline double frand(void) {
    // Divided by 2^32 to ensure the result is strictly < 1.0
    return (double)fast_rand() / 4294967296.0; 
}

#endif
