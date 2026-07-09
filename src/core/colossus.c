//
//  Colossus - a classical cipher solver
//

// A stochastic, slippery shotgun-restarted hill climber with backtracking for solving
// Vigenere, Beaufort, Porta, Quagmire I - IV, and Autokey ciphers with variants.

// One cipher-agnostic search engine -- run_solver() / run_one_config() -- drives EVERY
// cipher type. Each type supplies a CipherModel (colossus.h): a vtable of
// seed / perturb / decrypt / enumerate / report hooks plus a SearchShape (SHOTGUN slip,
// ANNEAL Metropolis, or DETERMINISTIC first-improvement). The engine owns the restart /
// hill-climb / accept / backtrack / best-tracking loop and the single state_score call;
// the model owns the cipher math. The polyalphabetic family (Vigenere/Quagmire/Beaufort/
// Porta/Autokey) is one model whose seed/perturb keep the explicit per-type
// switch(cipher_type) ladders. The transposition families are their own models:
//   - transmatrix / transperoffset : climb the transform's small parameter vector.
//   - transposition                : AZDecrypt-style full-permutation-key climb with a
//                                    periodic-redundancy structure term (key_structure_score,
//                                    weight -weightstructure) folded into the decrypt score.
//   - transcol / transcol2         : climb the per-stage column-order permutation(s).
//   - railfence / route            : exhaustive parameter sweeps (key_len 0 => no climb).
//   - amsco / myszkowski / redefence / cadenus / nihilist / swagman / grille : anneal a
//                                    short integer key, with a shared TransKeyOps seed/move.
// (indep_periodic keeps its own coordinate-ascent iterated-local-search climber, a search
// shape that does not fit the shotgun/anneal engine.)
//
// The -type values above are distinct from the -transmatrix / -transperoffset post-decrypt
// STAGE flags, which apply a fixed, user-supplied transposition after a polyalphabetic solve.

// Written by Sam Blake, started 14 July 2023.

/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
*/

// Reference for n-gram data: http://practicalcryptography.com/cryptanalysis/letter-frequencies-various-languages/english-letter-frequencies/

/* Usage
    -----
    $ ./colossus [options]

    Parameters
    ----------
    Input/Output:
        -cipher <file> : str
            Path to the ciphertext file. The cipher should be on the first line.
        -batch <file> : str
            Path to a file containing multiple ciphers (one per line) for batch processing.
        -crib <file> : str
            Path to the crib file. Use "_" for unknown characters. Must match cipher length.
        -dictionary <file> : str, optional
            Path to a dictionary file (one word per line). Defaults to "OxfordEnglishWords.txt".
        -verbose : flag
            Enable detailed output during execution.
        -seed <int> : int, optional
            Fix the PRNG seed for reproducible runs. Defaults to the current
            Unix time. Used by the regression tests to make stochastic solves
            deterministic.

    Cipher Configuration:
        -type <int> : int
            The cipher algorithm to solve:
            all, any, sweep   : run every plausible cipher type (a subprocess per type;
                                types the ciphertext plainly cannot be are skipped) and
                                print a best-score leaderboard.
            vigenere, vig, 0  : Vigenere
            quagmire1, quag1, q1, 1  : Quagmire I
            quagmire2, quag2, q2, 2  : Quagmire II
            quagmire3, quag3, q3, 3  : Quagmire III
            quagmire4, quag4, q4, 4  : Quagmire IV
            beaufort, beau, 5        : Beaufort
            porta, 6                 : Porta
            auto, autokey, 7         : Autokey (Vigenere tableau)
            auto1, autokey1, 8       : Autokey (Quagmire I tableau)
            auto2, autokey2, 9       : Autokey (Quagmire II tableau)
            auto3, autokey3, 10      : Autokey (Quagmire III tableau)
            auto4, autokey4, 11      : Autokey (Quagmire IV tableau)
            auto5, autobeau          : Autokey (Beaufort tableau)
            auto6, autoporta         : Autokey (Porta tableau)
            transmatrix, tmatrix, 14 : Transposition - double grid rotation (K3-style),
                                       solved by optimizing (w1, w2, direction).
            transperoffset, tpo, 15  : Transposition - periodic decimation + rotation,
                                       solved by optimizing (period d, offset n).
            transposition, trans, 16 : General transposition (columnar / route) - hill
                                       climbs the full permutation key (AZDecrypt-style),
                                       guarded by a periodic-structure score
                                       (-weightstructure). Restarts seeded from columnar
                                       layouts; column-swap move reorders whole columns.
            transcol, 17             : Dedicated single columnar - hill climbs only the
                                       column-order permutation (length K) via
                                       decrypt_columnar(), sweeping K over
                                       -mincols..-maxcols. Read direction via -readdir.
            transcol2, 18            : Dedicated double columnar - randomises the two
                                       column counts (K1, K2) per restart and anneals
                                       both column-order permutations together.
        -transperiodoffset <int> <int> : int, int
            Applies a periodic decimation and rotation to the decrypted text.
            The first integer specifies the offset (rotation), and the second 
            integer specifies the period (decimation step). 
            Aliases: -transperoffset, -transperoff
        -transmatrix <int> <int> <str> : int, int, str
            Applies a double matrix transposition (like Kryptos K3).
            The first integer is the initial grid width (w1), the second is 
            the subsequent grid width (w2). The string specifies the rotation 
            direction: 'cw' (clockwise) or 'ccw' (anti-clockwise).
        -variant : flag
            Enable the Quagmire variant (which swaps decryption for encryption.)
        -samekey : flag
            Forces the cycleword (indicator word) to be the same as the plaintext / ciphertext 
            keyword.

    Optimization Strategy:
        -optimalcycle : flag (default true)
            Enables hybrid deterministic solving. The cycleword is mathematically derived 
            (using Chi-squared/Dot-product) for every keyword candidate, rather than being 
            perturbed stochastically. Highly recommended for ciphers without cribs.
        -stochasticcycle : flag (default false)
            Enables stochastic solving for cycleword (indicator word.) The cycle is not 
            derived using the Chi-squared/Dot-product algorithm, rather it is perturbed 
            stochastically. 
        -nhillclimbs <int> : int
            Number of iterations per restart in the hill climber.
        -nrestarts <int> : int
            Number of times to restart the hill climber from a random state.

    Constraints (Lengths & Fixed Keywords):
        -plaintextkeyword <str> : str
            Fixes the plaintext keyword to a specific string.
        -ciphertextkeyword <str> : str
            Fixes the ciphertext keyword to a specific string.
        -cyclewordlen <int> : int
            Fixes the cycleword (period) length.
        -maxcyclewordlen <int> : int
            Maximum allowable length for the cycleword (if not fixed).
        -plaintextkeywordlen <int> : int
            Fixed length for the plaintext keyword.
        -ciphertextkeywordlen <int> : int
            Fixed length for the ciphertext keyword.
        -maxkeywordlen <int> : int
            Sets both plaintext and ciphertext max keyword lengths.
        -keywordlen <int> : int
            Sets both plaintext and ciphertext fixed lengths.

    Statistics & Resources:
        -ngramfile <file> : str
            Path to the n-gram statistics file.
        -ngramsize <int> : int
            The size of the n-grams (e.g., 3 for trigrams, 4 for quadgrams).

    Tuning Probabilities & Thresholds:
        -backtrackprob <float> : float
            Probability (0.0 - 1.0) of backtracking to the best known solution 
            instead of a random state during a restart.
        -keywordpermprob <float> : float
            Probability (0.0 - 1.0) of permuting the keyword vs the cycleword 
            (Ignored if -optimalcycle is enabled).
        -slipprob <float> : float
            Probability (0.0 - 1.0) of accepting a worse score to escape local maxima.
        -nsigmathreshold <float> : float
            Sigma threshold for cycleword length estimation via IoC.
        -iocthreshold <float> : float
            Minimum IoC required to consider a cycleword length valid.

    Scoring Weights:
        -weightngram <float> : float
            Weight of the n-gram score component.
        -weightcrib <float> : float
            Weight of the crib match component.
        -weightioc <float> : float
            Weight of the Index of Coincidence component.
        -weightentropy <float> : float
            Weight of the entropy component.
        -weightstructure <float> : float
            (General transposition only.) Weight of the columnar-structure reward
            that biases the permutation key toward regular (periodic) layouts.

    Columnar transposition (transcol / transcol2):
        -mincols <int> / -maxcols <int> : int
            Column-count search range (default 2..30, clamped to len/2). Set
            -mincols == -maxcols to target a single, known column count.
        -readdir <tb|bt|both> :
            Column read direction: top-to-bottom (default), bottom-to-top, or
            search both. Variants are only tried when explicitly requested.

    Notes
    -----
    Cipher Types:
        Ciphers are as defined by the American Cryptogram Association. 
        https://www.cryptogram.org/resource-area/cipher-types/


*/

#include "colossus.h"
#include "engine.h"
#include "scoring.h"
#include "polyalpha_solver.h"
#include "transmatrix_solver.h"
#include "permutation_solver.h"
#include "columnar_solver.h"
#include "columnar_track_solver.h"
#include "route_chain_solver.h"
#include "tile_solver.h"
#include "period_column_solver.h"
#include "period_column_space_solver.h"
#include "railfence_solver.h"
#include "route_solver.h"
#include "amsco_solver.h"
#include "myszkowski_solver.h"
#include "redefence_solver.h"
#include "cadenus_solver.h"
#include "nihilist_solver.h"
#include "swagman_solver.h"
#include "grille_solver.h"
#include "indep_solver.h"
#include "homophonic_solver.h"
#include "playfair_solver.h"
#include "bifid_solver.h"
#include "phillips_solver.h"
#include "twosquare_solver.h"
#include "foursquare_solver.h"
#include "trisquare_solver.h"
#include "trifid_solver.h"
#include "hill_solver.h"
#include "adfgvx_solver.h"
#include "nihilist_sub_solver.h"
#include "gromark_solver.h"
#include "nicodemus_solver.h"
#include "bazeries_solver.h"
#include "portax_solver.h"
#include "progkey_solver.h"
#include "slidefair_solver.h"
#include "seriated_playfair_solver.h"
#include "digrafid_solver.h"
#include "cm_bifid_solver.h"
#include "intkey_solver.h"
#include "condi_solver.h"
#include "fracmorse_solver.h"
#include "double_transposition_solver.h"
#include "pollux_solver.h"
#include "morbit_solver.h"
#include "straddling_checkerboard_solver.h"

#include <sys/wait.h>   // waitpid() for the "-type all" subprocess sweep

void init_config(ColossusConfig *cfg) {
    // Set Defaults
    cfg->cipher_type = -1;
    cfg->ngram_size = 0;
    cfg->n_hill_climbs = 1000;
    cfg->n_restarts = 1;

    cfg->ciphertext_keyword_len = 5;
    cfg->plaintext_keyword_len = 5;
    cfg->ciphertext_max_keyword_len = 12;
    cfg->min_keyword_len = 5;
    cfg->plaintext_max_keyword_len = 12;
    cfg->max_cycleword_len = 20;
    cfg->cycleword_len = 0; // 0 implies not set by user

    cfg->period = 0;            // bifid: 0 => estimate (not pinned by user)
    cfg->period_present = false;
    cfg->progression_present = false;   // progkey: 0 => sweep the progression index 0..25
    cfg->progression = 0;
    cfg->interruptor_present = false;   // intkey: else enumerate the 26 interruptor letters
    cfg->interruptor = 0;
    cfg->intscheme_present = false;     // intkey: else blind default (ct + pt)
    cfg->intscheme = IK_STRAT_CT;
    cfg->breaks_present = false;        // intkey: -breaks supplies known group-start positions
    cfg->breaks_file[0] = '\0';
    cfg->startkey_present = false;      // condi: else enumerate the 26 starter offsets 0..25
    cfg->startkey = 0;
    cfg->max_period = 0;        // 0 => derive from ciphertext length (min(20, len/2))
    cfg->n_periods = 5;         // anneal the estimator's top-K candidate periods
    cfg->n_primers = 0;         // Gromark pre-pass top-K (0 => auto by ciphertext length)

    cfg->plaintext_keyword_len_present = false;
    cfg->ciphertext_keyword_len_present = false;
    cfg->cycleword_len_present = false;
    cfg->user_plaintext_keyword_present = false;
    cfg->user_ciphertext_keyword_present = false;

    cfg->cipher_present = false;
    cfg->batch_present = false;
    cfg->crib_present = false;
    cfg->dictionary_present = false;
    cfg->verbose = false;
    cfg->skip_spaces = false;
    cfg->multiline = false;
    cfg->variant = false;
    cfg->beaufort = false;

    cfg->n_sigma_threshold = 1.0;
    cfg->ioc_threshold = 0.047;
    cfg->backtracking_probability = 0.15;
    cfg->keyword_permutation_probability = 0.95;
    cfg->slip_probability = 0.001;

    cfg->weight_ngram = 12.0;
    cfg->weight_crib = 36.0;
    cfg->weight_ioc = 0.0;
    cfg->weight_entropy = 0.0;
    cfg->weight_structure = 4.0;
    cfg->weight_monogram = 1.0;   // homophonic anti-collapse penalty (chi-squared vs English)

    cfg->optimal_cycleword = true;
    cfg->same_key_cycle = false;

    cfg->method = METHOD_DEFAULT;
    cfg->init_temp = 0.10;     // matches the previously hardcoded annealing schedule
    cfg->min_temp = 0.001;
    cfg->cooling_rate = 0.0;   // 0 => derive the geometric schedule over n_hill_climbs

    cfg->n_particles = 30;     // particle swarm (SHAPE_PSO, -method pso)
    cfg->inertia = 0.7;
    cfg->cognitive = 1.5;
    cfg->social = 1.5;
    cfg->refine_steps = 50;

    cfg->transperoffset_present = false;
    cfg->trans_offset = 0;
    cfg->trans_period = 1;

    cfg->transmatrix_present = false;
    cfg->trans_w1 = 0;
    cfg->trans_w2 = 0;
    cfg->trans_clockwise = 1; // Default to clockwise

    cfg->min_cols = 2;
    cfg->max_cols = 30;
    cfg->block_height = 0;       // Nicodemus: 0 => sweep block heights
    cfg->max_block_height = 0;   // 0 => default top of the sweep
    cfg->read_direction = COL_READ_TB; // canonical only; bottom-to-top is opt-in
    cfg->read_row_direction = ROW_READ_LR; // canonical only; right-to-left is opt-in (transcol-L/chain)
    cfg->crib_anchored = false;        // -cribanchored: structural crib constraint, off by default
    cfg->weight_word = 0.0;            // dictionary word-fraction reward off => bit-identical scoring
    cfg->tile_h = 2; cfg->tile_w = 2;  // default sub-grid tile shape for TRANSTILE
    cfg->trans_depth = 2;              // PERIOD_COLUMN: search compositions up to 2 stages
    cfg->max_gaps = 4;                 // PERIOD_COLUMN_SPACE: default blank/gap insertion budget
    cfg->max_dels = 2;                 // PERIOD_COLUMN_SPACE: default observed-cell deletion budget

    cfg->delimiter = 0;                 // 0 => per-character / 0..25 letter decode (ord())
    cfg->delimiter_present = false;
}


// Guarded so the regression tests can link the whole solver (solve_cipher and
// its dependencies live in this file) while supplying their own main:
// compile this translation unit with -DCOLOSSUS_NO_MAIN.
#ifndef COLOSSUS_NO_MAIN

// --- "-type all": sweep every plausible cipher type ------------------------------
//
// Every cipher type sets up its own alphabet size and n-gram base (25/26/27/36) and,
// for the digit ciphers, its own ciphertext parser -- all keyed on the type and fixed
// once, before load_ngrams(), in main(). Rather than tear that one-shot global setup
// apart to loop in-process, "-type all" re-runs this very binary once per candidate
// type (a subprocess inherits ALL the same flags with only the -type value replaced),
// so each solve gets its correct, untouched setup. A cheap pre-filter over the
// ciphertext skips types the input plainly cannot be (see cipher_type_plausible).

// Read the ciphertext sample used ONLY for the pre-filter, mirroring main's reader:
// the first line of -cipher, or the first non-trivial line of a -batch file. The real
// solve is done by each child on the untouched -cipher/-batch argument.
static void read_cipher_sample(int argc, char **argv, char *buf, int cap) {
    const char *cipher_file = NULL, *batch_file = NULL;
    bool multiline = false;
    buf[0] = '\0';
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-cipher") == 0 && i + 1 < argc) cipher_file = argv[i + 1];
        else if (strcmp(argv[i], "-batch") == 0 && i + 1 < argc) batch_file = argv[i + 1];
        else if (strcmp(argv[i], "-multiline") == 0) multiline = true;
    }
    const char *path = cipher_file ? cipher_file : batch_file;
    if (!path || !file_exists(path)) return;
    FILE *fp = fopen(path, "r");
    if (!fp) return;
    if (batch_file && !cipher_file) {
        char line[MAX_CIPHER_LENGTH];
        while (fgets(line, sizeof line, fp)) {
            line[strcspn(line, "\r\n")] = 0;
            if (strlen(line) >= 5) { strncpy(buf, line, cap - 1); buf[cap - 1] = '\0'; break; }
        }
    } else {
        int ci = 0, ch;
        while ((ch = fgetc(fp)) != EOF && (multiline || ch != '\n') && ci < cap - 1) {
            if (ch == '\r' || ch == '\n') continue;
            buf[ci++] = (char) ch;
        }
        buf[ci] = '\0';
    }
    fclose(fp);
}

// Fast, conservative "could this ciphertext be of this type?" test. Only returns
// false when the input is IMMEDIATELY, structurally incompatible with the type -- it
// never gambles on statistics, so a plausible type is always run.
//   - The six digit-stream ciphers (Pollux, Morbit, Straddling, Nihilist-Sub x3) need
//     an all-digit ciphertext; every letter-based type needs letters.
//   - Morbit's digits are 1..9, so a '0' rules it out.
//   - ADFGX / ADFGVX ciphertext is written only in their label letters.
static bool cipher_type_plausible(int type, const char *cipher) {
    bool has_letters = false, has_digits = false, has_zero = false;
    bool seen[26]; memset(seen, 0, sizeof seen);
    for (const char *p = cipher; *p; p++) {
        unsigned char c = (unsigned char) *p;
        if (isalpha(c)) { has_letters = true; seen[toupper(c) - 'A'] = true; }
        else if (isdigit(c)) { has_digits = true; if (c == '0') has_zero = true; }
    }
    if (!has_letters && !has_digits) return true;   // degenerate input: don't skip

    bool digit_family = (type == POLLUX || type == MORBIT ||
                         type == STRADDLING_CHECKERBOARD ||
                         type == NIHILIST_SUB || type == NIHILIST_SUB_NC ||
                         type == NIHILIST_SUB_M100);

    if (has_digits && !has_letters) {
        // Pure-digit ciphertext: only the digit-stream ciphers are possible.
        if (!digit_family) return false;
        if (type == MORBIT && has_zero) return false;   // Morbit emits digits 1..9
        return true;
    }

    // Ciphertext contains letters: the digit-stream ciphers are ruled out.
    if (digit_family) return false;

    // ADFGX / ADFGVX ciphertext uses only its fractionation-label letters.
    if (type == ADFGX || type == ADFGVX) {
        const char *labels = (type == ADFGX) ? "ADFGX" : "ADFGVX";
        for (int c = 0; c < 26; c++)
            if (seen[c] && !strchr(labels, 'A' + c)) return false;
    }
    return true;
}

// Run the "-type all" sweep: for every plausible cipher type, re-exec this binary
// with the -type value overridden to that type's integer code, tee the child's output
// through to our stdout, and record its best ">>>" score for a final leaderboard.
static int run_all_types(int argc, char **argv) {
    // Find the -type value slot so each child can override it in place.
    int type_slot = -1;
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-type") == 0 && i + 1 < argc) { type_slot = i + 1; break; }
    if (type_slot < 0) {
        printf("\n\nERROR: -type all needs a -type argument.\n\n");
        return 0;
    }

    char sample[MAX_CIPHER_LENGTH];
    read_cipher_sample(argc, argv, sample, sizeof sample);

    printf("\n=== -type all: sweeping every plausible cipher type ===\n");
    if (sample[0])
        printf("Ciphertext sample (%zu chars): %.72s%s\n",
               strlen(sample), sample, strlen(sample) > 72 ? " ..." : "");
    printf("Each type is solved in a subprocess with all your other flags unchanged.\n");

    // A mutable copy of argv, NUL-terminated for execvp.
    char **child_argv = malloc((size_t)(argc + 1) * sizeof *child_argv);
    for (int i = 0; i < argc; i++) child_argv[i] = argv[i];
    child_argv[argc] = NULL;

    typedef struct { int type; double score; bool ran; } SweepResult;
    SweepResult board[N_CIPHER_TYPES];
    int nboard = 0, n_run = 0, n_skipped = 0;

    for (int t = 0; t < N_CIPHER_TYPES; t++) {
        const char *name = cipher_type_name(t);
        if (!name) continue;   // not a real type code

        if (!cipher_type_plausible(t, sample)) {
            printf("  SKIP  type %2d  %-38s (ciphertext not of this form)\n", t, name);
            n_skipped++;
            continue;
        }

        printf("\n\n#################### -type %d  (%s) ####################\n", t, name);
        fflush(stdout);

        char numbuf[16];
        snprintf(numbuf, sizeof numbuf, "%d", t);
        child_argv[type_slot] = numbuf;

        int pipefd[2];
        if (pipe(pipefd) != 0) { perror("pipe"); continue; }
        pid_t pid = fork();
        if (pid < 0) { perror("fork"); close(pipefd[0]); close(pipefd[1]); continue; }
        if (pid == 0) {
            // Child: send stdout+stderr up the pipe, then become the concrete solve.
            dup2(pipefd[1], STDOUT_FILENO);
            dup2(pipefd[1], STDERR_FILENO);
            close(pipefd[0]);
            close(pipefd[1]);
            execvp(argv[0], child_argv);
            _exit(127);   // exec failed
        }

        // Parent: tee the child's output through, capturing its best ">>>" score.
        close(pipefd[1]);
        FILE *cout = fdopen(pipefd[0], "r");
        double best = 0.0; bool got = false;
        if (cout) {
            char line[MAX_CIPHER_LENGTH];
            while (fgets(line, sizeof line, cout)) {
                fputs(line, stdout);
                if (strncmp(line, ">>> ", 4) == 0) {
                    double sc;
                    if (sscanf(line + 4, "%lf", &sc) == 1 && (!got || sc > best)) {
                        best = sc; got = true;
                    }
                }
            }
            fclose(cout);
        }
        int status;
        waitpid(pid, &status, 0);

        board[nboard].type = t;
        board[nboard].score = got ? best : -1e30;
        board[nboard].ran = got;
        nboard++;
        n_run++;
    }
    free(child_argv);

    // Leaderboard, best score first.
    printf("\n\n=== -type all leaderboard  (%d run, %d skipped) ===\n", n_run, n_skipped);
    printf("NOTE: scores are NOT directly comparable across types with different\n");
    printf("      alphabets (25/26/27/36 symbols) or -logprob settings -- use as a\n");
    printf("      rough guide and read the decrypts of the top few candidates.\n\n");
    for (int a = 0; a < nboard; a++) {
        int bi = a;
        for (int b = a + 1; b < nboard; b++)
            if (board[b].score > board[bi].score) bi = b;
        if (bi != a) {
            SweepResult tmp = board[a]; board[a] = board[bi]; board[bi] = tmp;
        }
        if (board[a].ran)
            printf("  %8.2f   type %2d   %s\n",
                   board[a].score, board[a].type, cipher_type_name(board[a].type));
        else
            printf("     (no score)   type %2d   %s\n",
                   board[a].type, cipher_type_name(board[a].type));
    }
    printf("\n");
    return 1;
}

// --- "-help": full usage / option / cipher-type reference -------------------------
//
// Printed on -help / -h / --help (or when the binary is run with no arguments). The
// cipher-type table is driven by cipher_type_name() so it stays authoritative as
// types are added; the option list mirrors the parse loop in main() and the defaults
// in init_config().
static void print_help(const char *prog) {
    printf(
"COLOSSUS -- a polyalphabetic / polygraphic / transposition cipher solver.\n"
"A stochastic, slippery, shotgun-restarted hill climber with backtracking\n"
"(plus simulated annealing and particle-swarm search), built to crack the\n"
"Kryptos sculpture's K1-K4. Cipher conventions follow the ACA\n"
"(https://www.cryptogram.org/resource-area/cipher-types/).\n"
"\n"
"USAGE\n"
"  %s -type <type> -cipher <file> -ngramsize <n> -ngramfile <file> [options]\n"
"  %s -type <type> -batch  <file> -ngramsize <n> -ngramfile <file> [options]\n"
"  %s -help\n"
"\n"
"  Run from a directory holding the n-gram table, dictionary and ciphertext\n"
"  (paths are resolved relative to the current working directory).\n"
"\n"
"REQUIRED\n"
"  -type <type>            Cipher type: an alias or an integer code (see CIPHER TYPES).\n"
"                          'all' sweeps every plausible type in a subprocess each.\n"
"  -cipher <file>          Ciphertext file (only the FIRST line is read unless\n"
"                          -multiline). Mutually exclusive with -batch.\n"
"  -batch <file>           Solve every ciphertext line in the file in turn.\n"
"  -ngramsize <n>          N-gram order for scoring (typically 4, or 5 with -logprob).\n"
"  -ngramfile <file>       N-gram frequency table (e.g. english_quadgrams.txt).\n"
"\n"
"INPUT / OUTPUT\n"
"  -multiline              Read the whole -cipher file, concatenating lines into one\n"
"                          symbol stream (e.g. a homophonic grid).           [off]\n"
"  -delimiter <c|space>    Field separator for tokenized input; 'space'/'char'/'none'\n"
"                          selects per-character decode, else the char (e.g. ',').\n"
"                          [auto: comma for homophonic input, else per-character]\n"
"  -crib <file>            Known-plaintext crib (aligned to plaintext positions).\n"
"  -cribs <file>           Alias for -crib.\n"
"  -dictionary <file>      Word list for the readability report (-dict).\n"
"                          [auto-loads OxfordEnglishWords.txt if present]\n"
"  -excludeletter <L>      Drop a letter from the alphabet (e.g. J for 25-letter\n"
"                          square ciphers). Must precede n-gram load internally.\n"
"  -verbose                Stream a live search dialog.                        [off]\n"
"  -seed <uint>            Fix the PRNG seed for reproducible runs.   [time-based]\n"
"\n"
"SCORING / N-GRAMS\n"
"  -logprob                AZDecrypt-style log-probability n-grams with an unseen-\n"
"                          n-gram floor penalty (recommended for hard substitution\n"
"                          and square ciphers). Alias -azdecrypt.              [off]\n"
"  -reversengrams          Reversal-invariant table (each n-gram and its reverse\n"
"                          share the max weight). Alias -revngrams.            [off]\n"
"  -weightngram <f>        N-gram score weight.                              [12.0]\n"
"  -weightcrib <f>         Crib-match score weight.                          [36.0]\n"
"  -weightioc <f>          Index-of-coincidence score weight.                 [0.0]\n"
"  -weightentropy <f>      Entropy score weight.                              [0.0]\n"
"  -weightstructure <f>    Periodic-redundancy guard (general transposition). [4.0]\n"
"  -weightmono <f>         Homophonic monogram chi-squared anti-collapse.     [1.0]\n"
"  -weightword <f>         Dictionary word-coverage reward (transposition).   [0.0]\n"
"\n"
"SEARCH CONTROL\n"
"  -method <m>             shotgun | anneal (sa) | pso.        [per-type default]\n"
"  -nrestarts <n>          Random restarts (the robustness lever).             [1]\n"
"  -nhillclimbs <n>        Iterations per restart.                          [1000]\n"
"  -inittemp <f>           Annealing start temperature (-initialtemp).      [0.10]\n"
"  -mintemp <f>            Annealing floor temperature.                    [0.001]\n"
"  -coolingrate <f>        Geometric cooling rate; 0 => derive from budget.   [0.0]\n"
"  -backtrackprob <f>      Backtrack-to-best probability.                    [0.15]\n"
"  -slipprob <f>           Accept-worse (slip) probability (shotgun).       [0.001]\n"
"  -keywordpermprob <f>    Per-iteration keyword perturbation probability.  [0.95]\n"
"  -optimalcycle           Derive each column key deterministically (default).  [on]\n"
"  -stochasticcycle        Perturb the cycleword randomly instead.           [off]\n"
"\n"
"PARTICLE SWARM (-method pso)\n"
"  -nparticles <n>         Swarm size (-npart).                               [30]\n"
"  -inertia <f>            Inertia weight.                                   [0.7]\n"
"  -cognitive <f>          Cognitive (pbest) pull (-cog).                    [1.5]\n"
"  -social <f>             Social (gbest) pull (-soc).                       [1.5]\n"
"  -refine <n>             Greedy local-refinement steps per particle (-psorefine). [50]\n"
"\n"
"PERIODIC / KEYWORD\n"
"  -period <n>             Pin the fractionation/pairing period (else estimate/sweep).\n"
"  -maxperiod <n>          Upper bound for the period estimator. [min(20, len/2)]\n"
"  -nperiods <n>           Top-K candidate periods to anneal.                  [5]\n"
"  -cyclewordlen <n>       Pin the cycleword length.\n"
"  -maxcyclewordlen <n>    Upper bound for cycleword-length search.           [20]\n"
"  -keywordlen <n>         Pin the keyword length.\n"
"  -maxkeywordlen <n>      Upper bound for keyword-length search.\n"
"  -plaintextkeyword <w>   Pin the plaintext keyed alphabet.\n"
"  -ciphertextkeyword <w>  Pin the ciphertext keyed alphabet.\n"
"  -plaintextkeywordlen <n> / -ciphertextkeywordlen <n>   Pin those lengths.\n"
"  -nsigmathreshold <f>    IoC period-estimation Z-score threshold.          [1.0]\n"
"  -iocthreshold <f>       IoC acceptance threshold.                       [0.047]\n"
"  -samekey                Tie the keyword and cycleword together.           [off]\n"
"  -variant                Swap encrypt/decrypt in the Quagmire/Vigenere math.[off]\n"
"  -mincols <n>            Bottom of a swept column/period range.              [2]\n"
"  -maxcols <n>            Top of a swept column/period range.                [30]\n"
"\n"
"PER-CIPHER OPTIONS\n"
"  -nprimers <n>           Gromark primer pre-pass top-K.        [auto by length]\n"
"  -progression <n>        Progressive Key: pin the drift index (else sweep 0..25).\n"
"  -interruptor <A-Z>      Interrupted Key: pin the interruptor letter (else 26).\n"
"  -intscheme <s>          Interrupted Key strategy: ct | pt | breaks | joint.\n"
"  -breaks <file>          Interrupted Key: supplied group-start positions.\n"
"  -startkey <n>           Condi: pin the starter offset 0..25 (else enumerate).\n"
"  -blockheight <n>        Nicodemus: pin rows per block (else sweep).\n"
"  -maxblockheight <n>     Nicodemus: top of the block-height sweep.\n"
"  -skipspaces             Strip spaces before solving.                       [off]\n"
"\n"
"TRANSPOSITION\n"
"  -readdir <tb|bt|both>   Columnar read direction.                          [tb]\n"
"  -readrowdir <lr|rl|both> Within-row read direction (transcol-L/chain).     [lr]\n"
"  -cribanchored           Use the crib as a structural column-order constraint.[off]\n"
"  -tile <h> <w>           Sub-grid tile shape for transtile.                [2 2]\n"
"  -depth <1|2>            Period-column: max composed stages searched.         [2]\n"
"  -maxgaps <n>            Period-column-space: max inserted gap cells.         [4]\n"
"  -maxdels <n>            Period-column-space: max deleted observed cells.     [2]\n"
"  -transperoffset <o> <p> Post-decrypt periodic-decimation transposition stage.\n"
"  -transmatrix <w1> <w2> <cw|ccw>  Post-decrypt double-rotation transposition.\n"
"\n", prog, prog, prog);

    printf("CIPHER TYPES  (-type accepts the integer code or any listed alias)\n");
    printf("  code  name                                aliases\n");
    for (int t = 0; t <= STRADDLING_CHECKERBOARD; t++) {
        const char *name = cipher_type_name(t);
        if (name)
            printf("  %3d   %-35s %s\n", t, name, cipher_type_aliases(t));
    }
    printf("   --   %-35s %s\n", "sweep every plausible type", "all, any, sweep");
    printf(
"\n"
"EXAMPLES\n"
"  colossus -type q3 -cipher cipher.txt -ngramsize 4 -ngramfile english_quadgrams.txt\n"
"  colossus -type playfair -cipher pf.txt -ngramsize 5 -ngramfile english_quintgrams.txt -logprob\n"
"  colossus -type transcol -cipher ct.txt -ngramsize 4 -ngramfile english_quadgrams.txt -mincols 2 -maxcols 15\n"
"  colossus -type all -cipher ct.txt -ngramsize 4 -ngramfile english_quadgrams.txt\n"
"\n"
"See README.md and CLAUDE.md for the full writeup and per-type notes.\n"
"\n");
}

int main(int argc, char **argv) {
    ColossusConfig cfg;
    SharedData shared;
    int i;
    char single_ciphertext_buffer[MAX_CIPHER_LENGTH];
    char cribtext[MAX_CIPHER_LENGTH];

    printf("\n\nCOLOSSUS Cipher Solver\n\n");
    printf("Written by Sam Blake, started 14 July 2023.\n\n");

    // -help / -h / --help (or no arguments): print the full usage reference and exit.
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-h") == 0 ||
            strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "help") == 0) {
            print_help(argv[0]);
            return 0;
        }
    }
    if (argc == 1) {
        print_help(argv[0]);
        return 0;
    }

    // Seed the PRNG with the current Unix time (seconds since Epoch).
    // A -seed <uint> argument (parsed below) overrides this for reproducible runs.
    uint32_t rng_seed = (uint32_t)time(NULL);
    seed_rand(rng_seed);

    // "-type all": sweep every plausible cipher type in a subprocess each (each type
    // needs its own alphabet / n-gram setup, so a fresh process per type is the clean
    // way). Detected here so the normal single-type setup below is skipped entirely.
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-type") == 0 && i + 1 < argc &&
            parse_cipher_type(argv[i + 1]) == TYPE_ALL) {
            return run_all_types(argc, argv) ? 0 : 1;
        }
    }

    init_config(&cfg);

    // Default to the full 26-letter A..Z alphabet; -excludeletter (parsed below)
    // shrinks it (must happen before load_ngrams and any ord() call).
    init_alphabet(NULL);

    // Pre-scan -type and -method so a tuned per-type schedule (apply_cipher_defaults)
    // can be overlaid onto the init_config globals BEFORE the main parse loop runs --
    // any explicit -nrestarts/-inittemp/... below then overrides it (globals <
    // registry < CLI). Reading these two flags twice is harmless (the loop re-parses
    // and echoes them); the registry is a no-op for types without an entry.
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-type") == 0 && i + 1 < argc) {
            cfg.cipher_type = parse_cipher_type(argv[i + 1]);
        } else if (strcmp(argv[i], "-method") == 0 && i + 1 < argc) {
            char *m = argv[i + 1];
            if (strcasecmp(m, "shotgun") == 0) cfg.method = METHOD_SHOTGUN;
            else if (strcasecmp(m, "sa") == 0 || strcasecmp(m, "anneal") == 0 ||
                     strcasecmp(m, "simanneal") == 0 || strcasecmp(m, "simulatedannealing") == 0)
                cfg.method = METHOD_ANNEAL;
            else if (strcasecmp(m, "pso") == 0 || strcasecmp(m, "particle") == 0 ||
                     strcasecmp(m, "particleswarm") == 0 || strcasecmp(m, "swarm") == 0)
                cfg.method = METHOD_PSO;
        }
    }
    apply_cipher_defaults(&cfg, true);

    // Initialize shared data pointers.
    shared.ngram_data = NULL;
    shared.dict = NULL;
    shared.n_dict_words = 0;
    shared.max_dict_word_len = 0;

    // --- Argument Parsing ---
    for(i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-type") == 0) {
            char *type_arg = argv[++i];            
            cfg.cipher_type = parse_cipher_type(type_arg); 
            printf("-type %s\n", type_arg);      
        } else if (strcmp(argv[i], "-cipher") == 0) {
            cfg.cipher_present = true;
            strcpy(cfg.ciphertext_file, argv[++i]);
            printf("-cipher %s\n", cfg.ciphertext_file);
        } else if (strcmp(argv[i], "-batch") == 0) {
            cfg.batch_present = true;
            strcpy(cfg.batch_file, argv[++i]);
            printf("-batch %s\n", cfg.batch_file);
        } else if (strcmp(argv[i], "-crib") == 0 || strcmp(argv[i], "-cribs") == 0) {
            cfg.crib_present = true;
            strcpy(cfg.crib_file, argv[++i]);
            printf("-crib %s\n", cfg.crib_file);
        } else if (strcmp(argv[i], "-ngramsize") == 0) {
            cfg.ngram_size = atoi(argv[++i]);
            printf("-ngramsize %d\n", cfg.ngram_size);
        } else if (strcmp(argv[i], "-ngramfile") == 0) {
            strcpy(cfg.ngram_file, argv[++i]);
            printf("-ngramfile %s\n", cfg.ngram_file);
        } else if (strcmp(argv[i], "-excludeletter") == 0) {
            // Drop one (or more) letters from the alphabet, shrinking it to an
            // N<26 letter alphabet with mod-N arithmetic. E.g. -excludeletter P
            // gives the 25-letter A..Z-minus-P alphabet (mod 25). Must be set
            // before the ngram table is loaded and before any ciphertext is read.
            init_alphabet(argv[++i]);
            printf("-excludeletter %s  (alphabet size now %d: %s)\n",
                argv[i], g_alpha, g_idx_to_char_arr);
        } else if (strcmp(argv[i], "-maxkeywordlen") == 0) {
            cfg.plaintext_keyword_len = atoi(argv[++i]);
            cfg.ciphertext_keyword_len = cfg.plaintext_keyword_len;
            printf("-maxkeywordlen %d\n", cfg.plaintext_keyword_len);
        } else if (strcmp(argv[i], "-keywordlen") == 0) {
            cfg.plaintext_keyword_len_present = true;
            cfg.ciphertext_keyword_len_present = true;
            cfg.plaintext_keyword_len = atoi(argv[++i]);
            cfg.ciphertext_keyword_len = cfg.plaintext_keyword_len;
            cfg.plaintext_max_keyword_len = max(cfg.plaintext_max_keyword_len, 1 + cfg.plaintext_keyword_len);
            cfg.ciphertext_max_keyword_len = max(cfg.ciphertext_max_keyword_len, 1 + cfg.ciphertext_keyword_len);
            cfg.min_keyword_len = cfg.plaintext_keyword_len;
            printf("-keywordlen %d\n", cfg.plaintext_keyword_len);
        } else if (strcmp(argv[i], "-plaintextkeywordlen") == 0) {
            cfg.plaintext_keyword_len_present = true;
            cfg.plaintext_keyword_len = atoi(argv[++i]);
            cfg.plaintext_max_keyword_len = max(cfg.plaintext_max_keyword_len, 1 + cfg.plaintext_keyword_len);
            cfg.min_keyword_len = cfg.plaintext_keyword_len;
            printf("-plaintextkeywordlen %d\n", cfg.plaintext_keyword_len);
        } else if (strcmp(argv[i], "-ciphertextkeywordlen") == 0) {
            cfg.ciphertext_keyword_len_present = true;
            cfg.ciphertext_keyword_len = atoi(argv[++i]);
            cfg.ciphertext_max_keyword_len = max(cfg.ciphertext_max_keyword_len, 1 + cfg.ciphertext_keyword_len);
            cfg.min_keyword_len = cfg.ciphertext_keyword_len;
            printf("-ciphertextkeywordlen %d\n", cfg.ciphertext_keyword_len);
        } else if (strcmp(argv[i], "-plaintextkeyword") == 0) {
            // Explicit Plaintext Keyword
            cfg.user_plaintext_keyword_present = true;
            strcpy(cfg.user_plaintext_keyword, argv[++i]);
            int len = unique_len(cfg.user_plaintext_keyword);
            cfg.plaintext_keyword_len = len;
            cfg.plaintext_max_keyword_len = len + 1;
            cfg.plaintext_keyword_len_present = true;
            printf("-plaintextkeyword %s\n", cfg.user_plaintext_keyword);
        } else if (strcmp(argv[i], "-ciphertextkeyword") == 0) {
            // Explicit Ciphertext Keyword
            cfg.user_ciphertext_keyword_present = true;
            strcpy(cfg.user_ciphertext_keyword, argv[++i]);
            int len = unique_len(cfg.user_ciphertext_keyword);
            cfg.ciphertext_keyword_len = len;
            cfg.ciphertext_max_keyword_len = len + 1;
            cfg.ciphertext_keyword_len_present = true;
            printf("-ciphertextkeyword %s\n", cfg.user_ciphertext_keyword);
        } else if (strcmp(argv[i], "-maxcyclewordlen") == 0) {
            cfg.max_cycleword_len = atoi(argv[++i]);
            printf("-maxcyclewordlen %d\n", cfg.max_cycleword_len);
        } else if (strcmp(argv[i], "-cyclewordlen") == 0) {
            cfg.cycleword_len_present = true;
            cfg.cycleword_len = atoi(argv[++i]);
            cfg.max_cycleword_len = max(cfg.max_cycleword_len, 1 + cfg.cycleword_len);
            printf("-cyclewordlen %d\n", cfg.cycleword_len);
        } else if (strcmp(argv[i], "-nsigmathreshold") == 0) {
            cfg.n_sigma_threshold = atof(argv[++i]);
            printf("-nsigmathreshold %.4f\n", cfg.n_sigma_threshold);
        } else if (strcmp(argv[i], "-nhillclimbs") == 0) {
            cfg.n_hill_climbs = atoi(argv[++i]);
            printf("-nhillclimbs %d\n", cfg.n_hill_climbs);
        } else if (strcmp(argv[i], "-nrestarts") == 0) {
            cfg.n_restarts = atoi(argv[++i]);
            printf("-nrestarts %d\n", cfg.n_restarts);
        } else if (strcmp(argv[i], "-backtrackprob") == 0) {
            cfg.backtracking_probability = atof(argv[++i]);
            printf("-backtrackprob %.6f\n", cfg.backtracking_probability);
        } else if (strcmp(argv[i], "-keywordpermprob") == 0) {
            cfg.keyword_permutation_probability = atof(argv[++i]);
            printf("-keywordpermprob %.4f\n", cfg.keyword_permutation_probability);
        } else if (strcmp(argv[i], "-slipprob") == 0) {
            cfg.slip_probability = atof(argv[++i]);
            printf("-slipprob %.6f\n", cfg.slip_probability);
        } else if (strcmp(argv[i], "-method") == 0) {
            // Cipher-agnostic optimization strategy override.
            char *m = argv[++i];
            if (strcasecmp(m, "shotgun") == 0) {
                cfg.method = METHOD_SHOTGUN;
            } else if (strcasecmp(m, "sa") == 0 || strcasecmp(m, "anneal") == 0 ||
                       strcasecmp(m, "simanneal") == 0 || strcasecmp(m, "simulatedannealing") == 0) {
                cfg.method = METHOD_ANNEAL;
            } else if (strcasecmp(m, "pso") == 0 || strcasecmp(m, "particle") == 0 ||
                       strcasecmp(m, "particleswarm") == 0 || strcasecmp(m, "swarm") == 0) {
                cfg.method = METHOD_PSO;
            } else {
                printf("Unknown -method '%s' (expected shotgun | sa | anneal | simanneal | simulatedannealing | pso).\n", m);
                return 1;
            }
            printf("-method %s\n", cfg.method == METHOD_SHOTGUN ? "shotgun" :
                                   cfg.method == METHOD_PSO ? "particle-swarm" : "simulated-annealing");
        } else if (strcmp(argv[i], "-nparticles") == 0 || strcmp(argv[i], "-npart") == 0) {
            cfg.n_particles = atoi(argv[++i]);
            printf("-nparticles %d\n", cfg.n_particles);
        } else if (strcmp(argv[i], "-inertia") == 0) {
            cfg.inertia = atof(argv[++i]);
            printf("-inertia %.6f\n", cfg.inertia);
        } else if (strcmp(argv[i], "-cognitive") == 0 || strcmp(argv[i], "-cog") == 0) {
            cfg.cognitive = atof(argv[++i]);
            printf("-cognitive %.6f\n", cfg.cognitive);
        } else if (strcmp(argv[i], "-social") == 0 || strcmp(argv[i], "-soc") == 0) {
            cfg.social = atof(argv[++i]);
            printf("-social %.6f\n", cfg.social);
        } else if (strcmp(argv[i], "-refine") == 0 || strcmp(argv[i], "-psorefine") == 0) {
            cfg.refine_steps = atoi(argv[++i]);
            printf("-refine %d\n", cfg.refine_steps);
        } else if (strcmp(argv[i], "-inittemp") == 0 || strcmp(argv[i], "-initialtemp") == 0 ||
                   strcmp(argv[i], "-inittemperature") == 0 || strcmp(argv[i], "-initialtemperature") == 0) {
            cfg.init_temp = atof(argv[++i]);
            printf("-inittemp %.6f\n", cfg.init_temp);
        } else if (strcmp(argv[i], "-mintemp") == 0 || strcmp(argv[i], "-mintemperature") == 0) {
            cfg.min_temp = atof(argv[++i]);
            printf("-mintemp %.6f\n", cfg.min_temp);
        } else if (strcmp(argv[i], "-coolingrate") == 0 || strcmp(argv[i], "-cooling") == 0) {
            cfg.cooling_rate = atof(argv[++i]);
            printf("-coolingrate %.6f\n", cfg.cooling_rate);
        } else if (strcmp(argv[i], "-iocthreshold") == 0) {
            cfg.ioc_threshold = atof(argv[++i]);
            printf("-iocthreshold %.4f\n", cfg.ioc_threshold);
        } else if (strcmp(argv[i], "-dictionary") == 0 || strcmp(argv[i], "-dict") == 0) {
            cfg.dictionary_present = true;
            strcpy(cfg.dictionary_file, argv[++i]);
            printf("-dictionary %s\n", cfg.dictionary_file);
        } else if (strcmp(argv[i], "-weightngram") == 0) { 
            cfg.weight_ngram = atof(argv[++i]);
            printf("-weightngram %.4f\n", cfg.weight_ngram);
        } else if (strcmp(argv[i], "-weightcrib") == 0) { 
            cfg.weight_crib = atof(argv[++i]);
            printf("-weightcrib %.4f\n", cfg.weight_crib);
        } else if (strcmp(argv[i], "-weightioc") == 0) { 
            cfg.weight_ioc = atof(argv[++i]);
            printf("-weightioc %.4f\n", cfg.weight_ioc);
        } else if (strcmp(argv[i], "-weightentropy") == 0) {
            cfg.weight_entropy = atof(argv[++i]);
            printf("-weightentropy %.4f\n", cfg.weight_entropy);
        } else if (strcmp(argv[i], "-weightstructure") == 0) {
            cfg.weight_structure = atof(argv[++i]);
            printf("-weightstructure %.4f\n", cfg.weight_structure);
        } else if (strcmp(argv[i], "-variant") == 0) { 
            cfg.variant = true;
            printf("-variant\n");
        } else if (strcmp(argv[i], "-seed") == 0) {
            // Fix the PRNG seed for reproducible runs (regression tests, debugging).
            rng_seed = (uint32_t)strtoul(argv[++i], NULL, 10);
            seed_rand(rng_seed);
            printf("-seed %u\n", rng_seed);
        } else if (strcmp(argv[i], "-verbose") == 0) {
            cfg.verbose = true;
            printf("-verbose\n");
        } else if (strcmp(argv[i], "-skipspaces") == 0) {
            cfg.skip_spaces = true;
            printf("-skipspaces\n");
        } else if (strcmp(argv[i], "-multiline") == 0) {
            cfg.multiline = true;
            printf("-multiline\n");
        } else if (strcmp(argv[i], "-logprob") == 0 || strcmp(argv[i], "-azdecrypt") == 0) {
            g_ngram_logprob = true;
            printf("-logprob (AZDecrypt-style n-gram fitness: log-probabilities with an unseen-n-gram floor)\n");
        } else if (strcmp(argv[i], "-reversengrams") == 0 || strcmp(argv[i], "-revngrams") == 0) {
            g_ngram_reverse = true;
            printf("-reversengrams (reversal-invariant scoring: each n-gram and its reverse share the max weight; reads text with any words/segments reversed)\n");
        } else if (strcmp(argv[i], "-weightmono") == 0) {
            cfg.weight_monogram = atof(argv[++i]);
            printf("-weightmono %.3f\n", cfg.weight_monogram);
        } else if (strcmp(argv[i], "-delimiter") == 0) {
            // Field separator for tokenized input. The literal word "space" / "char"
            // (or an empty arg) selects per-character tokenization; otherwise the
            // first character of the argument is the delimiter (e.g. -delimiter ,).
            const char *d = argv[++i];
            if (str_eq(d, "space") || str_eq(d, "char") || str_eq(d, "none") || d[0] == '\0')
                cfg.delimiter = 0;
            else
                cfg.delimiter = d[0];
            cfg.delimiter_present = true;
            printf("-delimiter '%c' (code %d)\n", cfg.delimiter ? cfg.delimiter : ' ', cfg.delimiter);
        } else if (strcmp(argv[i], "-optimalcycle") == 0) {
            cfg.optimal_cycleword = true;
            printf("-optimalcycle\n");
        } else if (strcmp(argv[i], "-stochasticcycle") == 0) {
            cfg.optimal_cycleword = false;
            printf("-stochasticcycle\n");
        } else if (strcmp(argv[i], "-samekey") == 0) {
            cfg.same_key_cycle = true;
        } else if (strcmp(argv[i], "-transperiodoffset") == 0 || 
                strcmp(argv[i], "-transperoffset") == 0 || 
                strcmp(argv[i], "-transperoff") == 0) {
            cfg.trans_offset = atoi(argv[++i]);
            cfg.trans_period = atoi(argv[++i]);
            cfg.transperoffset_present = true;
            printf("-transperiodoffset %d %d\n", cfg.trans_offset, cfg.trans_period);
        } else if (strcmp(argv[i], "-transmatrix") == 0) {
            cfg.transmatrix_present = true;
            cfg.trans_w1 = atoi(argv[++i]);
            cfg.trans_w2 = atoi(argv[++i]);
            
            char *dir_arg = argv[++i];
            // Flexible parsing for clockwise vs anti-clockwise
            if (strcasecmp(dir_arg, "cw") == 0 || strcasecmp(dir_arg, "clockwise") == 0 || strcmp(dir_arg, "1") == 0) {
                cfg.trans_clockwise = 1;
            } else if (strcasecmp(dir_arg, "ccw") == 0 || strcasecmp(dir_arg, "anticlockwise") == 0 || strcmp(dir_arg, "0") == 0) {
                cfg.trans_clockwise = 0;
            } else {
                printf("\n\nERROR: Invalid direction '%s' for -transmatrix. Use 'cw' or 'ccw'.\n\n", dir_arg);
                return 0;
            }
            printf("-transmatrix %d %d %s\n", cfg.trans_w1, cfg.trans_w2, cfg.trans_clockwise ? "cw" : "ccw");
        } else if (strcmp(argv[i], "-mincols") == 0) {
            cfg.min_cols = atoi(argv[++i]);
            printf("-mincols %d\n", cfg.min_cols);
        } else if (strcmp(argv[i], "-maxcols") == 0) {
            cfg.max_cols = atoi(argv[++i]);
            printf("-maxcols %d\n", cfg.max_cols);
        } else if (strcmp(argv[i], "-blockheight") == 0) {
            cfg.block_height = atoi(argv[++i]);   // Nicodemus: pin the rows-per-block
            printf("-blockheight %d\n", cfg.block_height);
        } else if (strcmp(argv[i], "-maxblockheight") == 0) {
            cfg.max_block_height = atoi(argv[++i]); // Nicodemus: top of the block-height sweep
            printf("-maxblockheight %d\n", cfg.max_block_height);
        } else if (strcmp(argv[i], "-readdir") == 0) {
            char *dir_arg = argv[++i];
            // Flexible parsing of the columnar read direction.
            if (strcasecmp(dir_arg, "tb") == 0 || strcasecmp(dir_arg, "topbottom") == 0 || strcmp(dir_arg, "0") == 0) {
                cfg.read_direction = COL_READ_TB;
            } else if (strcasecmp(dir_arg, "bt") == 0 || strcasecmp(dir_arg, "bottomtop") == 0 || strcmp(dir_arg, "1") == 0) {
                cfg.read_direction = COL_READ_BT;
            } else if (strcasecmp(dir_arg, "both") == 0 || strcmp(dir_arg, "2") == 0) {
                cfg.read_direction = COL_READ_BOTH;
            } else {
                printf("\n\nERROR: Invalid direction '%s' for -readdir. Use 'tb', 'bt' or 'both'.\n\n", dir_arg);
                return 0;
            }
            printf("-readdir %s\n", dir_arg);
        } else if (strcmp(argv[i], "-readrowdir") == 0) {
            // Row read direction for transcol-L / transroutecol (default lr).
            char *dir_arg = argv[++i];
            if (strcasecmp(dir_arg, "lr") == 0 || strcasecmp(dir_arg, "leftright") == 0 || strcmp(dir_arg, "0") == 0) {
                cfg.read_row_direction = ROW_READ_LR;
            } else if (strcasecmp(dir_arg, "rl") == 0 || strcasecmp(dir_arg, "rightleft") == 0 || strcmp(dir_arg, "1") == 0) {
                cfg.read_row_direction = ROW_READ_RL;
            } else if (strcasecmp(dir_arg, "both") == 0 || strcmp(dir_arg, "2") == 0) {
                cfg.read_row_direction = ROW_READ_BOTH;
            } else {
                printf("\n\nERROR: Invalid direction '%s' for -readrowdir. Use 'lr', 'rl' or 'both'.\n\n", dir_arg);
                return 0;
            }
            printf("-readrowdir %s\n", dir_arg);
        } else if (strcmp(argv[i], "-cribanchored") == 0) {
            // Use the crib as a STRUCTURAL constraint on the transposition column
            // order (transcol / transcol-L), not just a score term (see Rec 3).
            cfg.crib_anchored = true;
            printf("-cribanchored\n");
        } else if (strcmp(argv[i], "-weightword") == 0) {
            // Optional dictionary word-fraction reward for the space-preserving
            // transposition solvers (default 0 => scoring is bit-identical).
            cfg.weight_word = atof(argv[++i]);
            printf("-weightword %.4f\n", cfg.weight_word);
        } else if (strcmp(argv[i], "-tile") == 0) {
            // Sub-grid tile shape h x w for TRANSTILE (e.g. -tile 2 2).
            cfg.tile_h = atoi(argv[++i]);
            cfg.tile_w = atoi(argv[++i]);
            printf("-tile %dx%d\n", cfg.tile_h, cfg.tile_w);
        } else if (strcmp(argv[i], "-depth") == 0) {
            // PERIOD_COLUMN: max number of composed period-column stages to search
            // (1 or 2; exhaustive). Default 2.
            cfg.trans_depth = atoi(argv[++i]);
            printf("-depth %d\n", cfg.trans_depth);
        } else if (strcmp(argv[i], "-maxgaps") == 0) {
            // PERIOD_COLUMN_SPACE: max number of blank/gap cells the solver may
            // insert (dropped-character repair / grid re-factorisation). Default 4.
            cfg.max_gaps = atoi(argv[++i]);
            printf("-maxgaps %d\n", cfg.max_gaps);
        } else if (strcmp(argv[i], "-maxdels") == 0) {
            // PERIOD_COLUMN_SPACE: max number of observed cells the solver may
            // delete (spuriously-added-character repair). Default 2.
            cfg.max_dels = atoi(argv[++i]);
            printf("-maxdels %d\n", cfg.max_dels);
        } else if (strcmp(argv[i], "-period") == 0) {
            // Bifid/Trifid: pin the fractionation period (block size) vs estimating it.
            cfg.period_present = true;
            cfg.period = atoi(argv[++i]);
            printf("-period %d\n", cfg.period);
        } else if (strcmp(argv[i], "-maxperiod") == 0) {
            // Bifid/Trifid: largest period the IoC estimator scans (default min(20, len/2)).
            cfg.max_period = atoi(argv[++i]);
            printf("-maxperiod %d\n", cfg.max_period);
        } else if (strcmp(argv[i], "-nperiods") == 0) {
            // Bifid/Trifid: how many top-IoC candidate periods to anneal (default 5).
            cfg.n_periods = atoi(argv[++i]);
            printf("-nperiods %d\n", cfg.n_periods);
        } else if (strcmp(argv[i], "-nprimers") == 0) {
            cfg.n_primers = atoi(argv[++i]);
            printf("-nprimers %d\n", cfg.n_primers);
        } else if (strcmp(argv[i], "-progression") == 0) {
            // Progressive Key: pin the per-group progression index (else swept 0..25).
            cfg.progression_present = true;
            cfg.progression = atoi(argv[++i]);
            printf("-progression %d\n", cfg.progression);
        } else if (strcmp(argv[i], "-interruptor") == 0) {
            // Interrupted Key: pin the interruptor letter (else the 26 letters are enumerated).
            char c = toupper((unsigned char) argv[++i][0]);
            if (c < 'A' || c > 'Z') { printf("\n\nERROR: -interruptor needs a letter A..Z\n\n"); return 0; }
            cfg.interruptor_present = true;
            cfg.interruptor = c - 'A';
            printf("-interruptor %c\n", c);
        } else if (strcmp(argv[i], "-intscheme") == 0) {
            // Interrupted Key: pin the interruption strategy (else blind ct+pt, or breaks if -breaks).
            const char *s = argv[++i];
            if (strcmp(s, "ct") == 0)          cfg.intscheme = IK_STRAT_CT;
            else if (strcmp(s, "pt") == 0)     cfg.intscheme = IK_STRAT_PT;
            else if (strcmp(s, "breaks") == 0) cfg.intscheme = IK_STRAT_BREAKS;
            else if (strcmp(s, "joint") == 0)  cfg.intscheme = IK_STRAT_JOINT;
            else { printf("\n\nERROR: -intscheme must be ct|pt|breaks|joint\n\n"); return 0; }
            cfg.intscheme_present = true;
            printf("-intscheme %s\n", s);
        } else if (strcmp(argv[i], "-breaks") == 0) {
            // Interrupted Key: file of whitespace-separated 0-based group-start positions.
            cfg.breaks_present = true;
            strncpy(cfg.breaks_file, argv[++i], MAX_FILENAME_LEN - 1);
            cfg.breaks_file[MAX_FILENAME_LEN - 1] = '\0';
            printf("-breaks %s\n", cfg.breaks_file);
        } else if (strcmp(argv[i], "-startkey") == 0) {
            // Condi: pin the starter offset 0..25 (else the 26 values are enumerated).
            cfg.startkey_present = true;
            cfg.startkey = atoi(argv[++i]);
            printf("-startkey %d\n", cfg.startkey);
        } else {
            printf("\n\nERROR: unknown command line arg: \'%s\'\n\n", argv[i]);
            return 0;
        }
    }

    printf("\n\n");

    if (cfg.cipher_type == VIGENERE) {
        printf("\nAttacking a Vigenere cipher.\n\n");
    } else if (cfg.cipher_type == GRONSFELD) {
        printf("\nAttacking a Gronsfeld cipher (Vigenere with a numeric key, shifts 0-9).\n\n");
    } else if (cfg.cipher_type == QUAGMIRE_1) {
        printf("\nAttacking a Quagmire I cipher.\n\n");
    } else if (cfg.cipher_type == QUAGMIRE_2) {
        printf("\nAttacking a Quagmire II cipher.\n\n");
    } else if (cfg.cipher_type == QUAGMIRE_3) {
        printf("\nAttacking a Quagmire III cipher.\n\n");
    } else if (cfg.cipher_type == QUAGMIRE_4) {
        printf("\nAttacking a Quagmire IV cipher.\n\n");
    } else if (cfg.cipher_type == BEAUFORT) {
        printf("\nAttacking a Beaufort cipher.\n\n");
    } else if (cfg.cipher_type == PORTA) {
        printf("\nAttacking a Porta cipher.\n\n");
    } else if (cfg.cipher_type == AUTOKEY_0) {
        printf("\nAttacking a Autokey cipher (Vigenere tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_1) {
        printf("\nAttacking a Autokey cipher (Quagmire I tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_2) {
        printf("\nAttacking a Autokey cipher (Quagmire II tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_3) {
        printf("\nAttacking a Autokey cipher (Quagmire III tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_4) {
        printf("\nAttacking a Autokey cipher (Quagmire IV tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_BEAU) {
        printf("\nAttacking a Autokey cipher (Beaufort tableau.)\n\n");
    } else if (cfg.cipher_type == AUTOKEY_PORTA) {
        printf("\nAttacking a Autokey cipher (Porta tableau.)\n\n");
    } else if (cfg.cipher_type == TRANSMATRIX) {
        printf("\nAttacking a transmatrix (double grid rotation) transposition cipher.\n\n");
    } else if (cfg.cipher_type == TRANSPEROFFSET) {
        printf("\nAttacking a transperoffset (periodic decimation + rotation) transposition cipher.\n\n");
    } else if (cfg.cipher_type == TRANSPOSITION) {
        printf("\nAttacking a general transposition cipher (permutation-key hill climber).\n\n");
    } else if (cfg.cipher_type == TRANSCOL) {
        printf("\nAttacking a columnar transposition cipher (column-order hill climber).\n\n");
    } else if (cfg.cipher_type == TRANSCOL2) {
        printf("\nAttacking a double columnar transposition cipher (column-order hill climber).\n\n");
    } else if (cfg.cipher_type == TRANSCOL_L) {
        printf("\nAttacking a columnar transposition with a within-column track permutation (exact seam best-L).\n\n");
    } else if (cfg.cipher_type == TRANSROUTECOL) {
        printf("\nAttacking a route + column-key two-stage transposition chain (seam best-L reading).\n\n");
    } else if (cfg.cipher_type == TRANSTILE) {
        printf("\nAttacking a sub-grid / tile transposition (uniform h x w tile cell permutation).\n\n");
    } else if (cfg.cipher_type == PERIOD_COLUMN) {
        printf("\nAttacking a period column order transposition (periodic column permutation, composable to 2 stages).\n\n");
    } else if (cfg.cipher_type == PERIOD_COLUMN_SPACE) {
        printf("\nAttacking a period column order transposition, space-robust (inserts searched blank/gap cells for dropped-character repair).\n\n");
    } else if (cfg.cipher_type == TRANSCOL2_DC) {
        printf("\nAttacking a double columnar transposition (divide & conquer: Index of Digraphic Potential scores K2 independently of K1).\n\n");
    } else if (cfg.cipher_type == RAILFENCE) {
        printf("\nAttacking a rail fence transposition cipher (rail-count + phase enumeration).\n\n");
    } else if (cfg.cipher_type == ROUTE) {
        printf("\nAttacking a route transposition cipher (grid + route enumeration).\n\n");
    } else if (cfg.cipher_type == AMSCO) {
        printf("\nAttacking an Amsco transposition cipher (column-order hill climber).\n\n");
    } else if (cfg.cipher_type == MYSZKOWSKI) {
        printf("\nAttacking a Myszkowski transposition cipher (rank-vector hill climber).\n\n");
    } else if (cfg.cipher_type == REDEFENCE) {
        printf("\nAttacking a redefence (keyed rail fence) cipher (rail-order hill climber).\n\n");
    } else if (cfg.cipher_type == CADENUS) {
        printf("\nAttacking a Cadenus transposition cipher (order + rotation hill climber).\n\n");
    } else if (cfg.cipher_type == NIHILIST) {
        printf("\nAttacking a Nihilist transposition cipher (single-permutation hill climber).\n\n");
    } else if (cfg.cipher_type == SWAGMAN) {
        printf("\nAttacking a Swagman transposition cipher (key-square hill climber).\n\n");
    } else if (cfg.cipher_type == GRILLE) {
        printf("\nAttacking a turning grille transposition cipher (orbit-assignment hill climber).\n\n");
    } else if (cfg.cipher_type == INDEP_PERIODIC) {
        printf("\nAttacking an independent-periodic substitution (P independent mixed alphabets, joint hill climber).\n\n");
    } else if (cfg.cipher_type == HOMOPHONIC) {
        printf("\nAttacking a homophonic substitution (ciphertext alphabet larger than the plaintext alphabet).\n\n");
    } else if (cfg.cipher_type == PLAYFAIR) {
        printf("\nAttacking a Playfair cipher (digraphic substitution over a 5x5 keyed grid).\n\n");
    } else if (cfg.cipher_type == SERIATED_PLAYFAIR) {
        printf("\nAttacking a Seriated Playfair cipher (digraphic Playfair over vertical pairs of a two-row seriated layout).\n\n");
    } else if (cfg.cipher_type == BIFID) {
        printf("\nAttacking a Bifid cipher (fractionation over a keyed Polybius square).\n\n");
    } else if (cfg.cipher_type == TRIFID) {
        printf("\nAttacking a Trifid cipher (fractionation over a keyed 3x3x3 cube).\n\n");
    } else if (cfg.cipher_type == HILL) {
        printf("\nAttacking a Hill cipher (polygraphic substitution by a k x k matrix mod 26).\n\n");
    } else if (cfg.cipher_type == PHILLIPS || cfg.cipher_type == PHILLIPS_C ||
               cfg.cipher_type == PHILLIPS_RC) {
        printf("\nAttacking a Phillips cipher (8-square keyed-Polybius substitution, period 40).\n\n");
    } else if (cfg.cipher_type == TWO_SQUARE || cfg.cipher_type == TWO_SQUARE_V) {
        printf("\nAttacking a Two-Square cipher (digraphic substitution over two keyed 5x5 squares).\n\n");
    } else if (cfg.cipher_type == FOUR_SQUARE) {
        printf("\nAttacking a Four-Square cipher (digraphic substitution over four 5x5 squares).\n\n");
    } else if (cfg.cipher_type == TRI_SQUARE) {
        printf("\nAttacking a Tri-Square cipher (digraphic substitution over three keyed 5x5 squares; digraph -> trigraph).\n\n");
    } else if (cfg.cipher_type == ADFGX) {
        printf("\nAttacking an ADFGX cipher (5x5 keyed-square fractionation + keyed columnar transposition).\n\n");
    } else if (cfg.cipher_type == ADFGVX) {
        printf("\nAttacking an ADFGVX cipher (6x6 keyed-square fractionation + keyed columnar transposition).\n\n");
    } else if (cfg.cipher_type == NIHILIST_SUB || cfg.cipher_type == NIHILIST_SUB_NC ||
               cfg.cipher_type == NIHILIST_SUB_M100) {
        printf("\nAttacking a Nihilist Substitution cipher (periodic additive over a keyed Polybius square, %s).\n\n",
            cfg.cipher_type == NIHILIST_SUB_NC ? "no-carry mod 10" :
            cfg.cipher_type == NIHILIST_SUB_M100 ? "add mod 100" : "integer add with carry");
    } else if (cfg.cipher_type == GROMARK || cfg.cipher_type == GROMARK_PERIODIC) {
        printf("\nAttacking a %s Gromark cipher (keyed-alphabet substitution + chain-addition running key%s).\n\n",
            cfg.cipher_type == GROMARK_PERIODIC ? "Periodic" : "basic",
            cfg.cipher_type == GROMARK_PERIODIC ? " + per-group offset" : "");
    } else if (cfg.cipher_type == NICODEMUS || cfg.cipher_type == NICODEMUS_VARIANT ||
               cfg.cipher_type == NICODEMUS_BEAUFORT) {
        printf("\nAttacking a Nicodemus cipher (periodic %s substitution + per-block columnar transposition).\n\n",
            cfg.cipher_type == NICODEMUS_VARIANT ? "Variant"
            : cfg.cipher_type == NICODEMUS_BEAUFORT ? "Beaufort" : "Vigenere");
    } else if (cfg.cipher_type == BAZERIES) {
        printf("\nAttacking a Bazeries cipher (keyed-square substitution + digit-grouped reversal, one number key).\n\n");
    } else if (cfg.cipher_type == PORTAX) {
        printf("\nAttacking a Portax cipher (periodic digraphic Porta: vertical pairs over a Porta slide).\n\n");
    } else if (cfg.cipher_type == PROGKEY || cfg.cipher_type == PROGKEY_VAR ||
               cfg.cipher_type == PROGKEY_BEAU) {
        printf("\nAttacking a Progressive Key cipher (periodic %s + per-group constant key drift).\n\n",
            cfg.cipher_type == PROGKEY_VAR ? "Variant"
            : cfg.cipher_type == PROGKEY_BEAU ? "Beaufort" : "Vigenere");
    } else if (cfg.cipher_type == SLIDEFAIR || cfg.cipher_type == SLIDEFAIR_VAR ||
               cfg.cipher_type == SLIDEFAIR_BEAU) {
        printf("\nAttacking a Slidefair cipher (periodic digraphic %s: rectangle over a shift slide).\n\n",
            cfg.cipher_type == SLIDEFAIR_VAR ? "Variant"
            : cfg.cipher_type == SLIDEFAIR_BEAU ? "Beaufort" : "Vigenere");
    } else if (cfg.cipher_type == DIGRAFID) {
        printf("\nAttacking a Digrafid cipher (digraphic fractionation over two keyed 27-symbol alphabets).\n\n");
    } else if (cfg.cipher_type == CM_BIFID) {
        printf("\nAttacking a CM Bifid cipher (Bifid fractionation over two keyed Polybius squares).\n\n");
    } else if (cfg.cipher_type == INTERRUPTED_KEY || cfg.cipher_type == INTERRUPTED_KEY_VAR ||
               cfg.cipher_type == INTERRUPTED_KEY_BEAU) {
        printf("\nAttacking an Interrupted Key cipher (periodic %s keyword reset at break points).\n\n",
            cfg.cipher_type == INTERRUPTED_KEY_VAR ? "Variant"
            : cfg.cipher_type == INTERRUPTED_KEY_BEAU ? "Beaufort" : "Vigenere");
    } else if (cfg.cipher_type == CONDI) {
        printf("\nAttacking a Condi cipher (plaintext-feedback substitution over a keyed alphabet).\n\n");
    } else if (cfg.cipher_type == FRAC_MORSE) {
        printf("\nAttacking a Fractionated Morse cipher (Morse fractionation over a keyed 26-letter alphabet).\n\n");
    } else if (cfg.cipher_type == POLLUX) {
        printf("\nAttacking a Pollux cipher (Morse over a digit -> dot/dash/divider map; deterministic exhaustive 3^10 search).\n\n");
    } else if (cfg.cipher_type == MORBIT) {
        printf("\nAttacking a Morbit cipher (Morse taken in pairs over a pair <-> digit map; deterministic exhaustive 9! search).\n\n");
    } else if (cfg.cipher_type == STRADDLING_CHECKERBOARD) {
        printf("\nAttacking a Straddling Checkerboard cipher (keyed-board digit fractionation; keyed labels + figure-shift).\n\n");
    } else {
        printf("\n\nERROR: Unknown cipher type %d.\n\n", cfg.cipher_type);
        return 0;
    }


    if (cfg.cipher_type == BEAUFORT) {
        cfg.beaufort = true;
    }

    // --- Validation ---
    if (cfg.cipher_type == -1) {
        printf("\n\nERROR: missing cipher type. Use -type /name or integer code/. \n\n");
        return 0;
    }

    if (!cfg.cipher_present && !cfg.batch_present) {
        printf("\n\nERROR: No cipher input specified. Use -cipher or -batch.\n\n");
        return 0;
    }

    if (cfg.ngram_size == 0) {
        printf("\n\nERROR: -ngramsize missing.\n\n");
        return 0;
    }
    if (!file_exists(cfg.ngram_file)) {
        printf("\nERROR: missing file '%s'\n", cfg.ngram_file);
        return 0;
    }

    // Default Dictionary Check
    char oxford_english_words[] = "OxfordEnglishWords.txt";
    if (!cfg.dictionary_present && file_exists(oxford_english_words)) {
        cfg.dictionary_present = true;
        strcpy(cfg.dictionary_file, oxford_english_words);
        if (cfg.verbose) printf("\nDefault dictionary = %s\n", cfg.dictionary_file);
    }

    // Playfair runs on a 25-letter grid: force a 25-letter alphabet (J merged into I
    // by ACA convention) unless the user has already shrunk it with -excludeletter.
    // Must happen before load_ngrams so the n-gram table is built over the same 25
    // letters (mod-25 packing), and before any ciphertext is decoded.
    if (cfg.cipher_type == PLAYFAIR && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type playfair: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Seriated Playfair runs on the same 5x5 (25-letter, J->I) grid as Playfair.
    if (cfg.cipher_type == SERIATED_PLAYFAIR && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type seriated-playfair: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Bifid defaults to the same 5x5 (25-letter, J->I) square as Playfair. Force it
    // here -- before load_ngrams -- unless the user already shrank the alphabet (e.g.
    // -excludeletter for a different excluded letter, or a 36-letter 6x6 alphabet).
    if (cfg.cipher_type == BIFID && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type bifid: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // CM Bifid runs on the same 5x5 (25-letter, J->I) grid as Bifid (two such squares).
    // Force it here -- before load_ngrams -- unless the user already shrank the alphabet.
    if (cfg.cipher_type == CM_BIFID && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type cm-bifid: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Phillips runs on the same 5x5 (25-letter, J->I) grid as Playfair/Bifid. Force it
    // here -- before load_ngrams -- unless the user already shrank the alphabet.
    if ((cfg.cipher_type == PHILLIPS || cfg.cipher_type == PHILLIPS_C ||
         cfg.cipher_type == PHILLIPS_RC) && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type phillips: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Two-Square, Four-Square and Tri-Square run on the same 5x5 (25-letter, J->I) squares as
    // Playfair. Force it here -- before load_ngrams -- unless the user already shrank the alphabet.
    if ((cfg.cipher_type == TWO_SQUARE || cfg.cipher_type == TWO_SQUARE_V ||
         cfg.cipher_type == FOUR_SQUARE || cfg.cipher_type == TRI_SQUARE) &&
        g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type %s: alphabet forced to %d letters (J->I): %s\n",
            (cfg.cipher_type == FOUR_SQUARE) ? "foursquare" :
            (cfg.cipher_type == TRI_SQUARE)  ? "trisquare" : "twosquare",
            g_alpha, g_idx_to_char_arr);
    }

    // Trifid runs on a 27-symbol cube (A..Z + '+'): force that alphabet here -- before
    // load_ngrams, so the n-gram table is built over the same 27 symbols (base-27
    // packing) and the ciphertext '+' decodes -- unless the user already changed it.
    if (cfg.cipher_type == TRIFID && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet_trifid();
        printf("-type trifid: alphabet forced to %d symbols (A..Z + '%c'): %s\n",
            g_alpha, TRIFID_EXTRA_CHAR, g_idx_to_char_arr);
    }

    // Digrafid runs on a 27-symbol alphabet (A..Z + '#') over two grids (3x9 and 9x3):
    // force that alphabet here -- before load_ngrams, so the n-gram table is built over the
    // same 27 symbols (base-27 packing) and the ciphertext '#' decodes -- unless the user
    // already changed it.
    if (cfg.cipher_type == DIGRAFID && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet_digrafid();
        printf("-type digrafid: alphabet forced to %d symbols (A..Z + '%c'): %s\n",
            g_alpha, DIGRAFID_EXTRA_CHAR, g_idx_to_char_arr);
    }

    // ADFGX runs on the same 5x5 (25-letter, J->I) square as Playfair/Bifid; ADFGVX on a
    // 6x6 (36-symbol, A..Z + 0..9) square. Force the alphabet here -- before load_ngrams,
    // so the n-gram table is built over the same symbols and the ciphertext labels decode
    // -- unless the user already shrank/changed it.
    if (cfg.cipher_type == ADFGX && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type adfgx: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }
    if (cfg.cipher_type == ADFGVX && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet_adfgvx();
        printf("-type adfgvx: alphabet forced to %d symbols (A..Z + 0..9): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Nihilist Substitution runs on the same 5x5 (25-letter, J->I) square as Playfair/Bifid.
    // Force it here -- before load_ngrams -- unless the user already shrank the alphabet.
    if ((cfg.cipher_type == NIHILIST_SUB || cfg.cipher_type == NIHILIST_SUB_NC ||
         cfg.cipher_type == NIHILIST_SUB_M100) && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type nihilist-sub: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Bazeries runs on the same 5x5 (25-letter, J->I) square as Playfair/Bifid. Force it
    // here -- before load_ngrams -- unless the user already shrank the alphabet.
    if (cfg.cipher_type == BAZERIES && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet("J");
        printf("-type bazeries: alphabet forced to %d letters (J->I): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // Straddling Checkerboard reuses the 36-symbol ADFGVX alphabet (A..Z + 0..9) so that
    // figure-shifted numeric plaintext is representable and scored (at negligible weight).
    // Force it here -- before load_ngrams -- unless the user already changed the alphabet.
    if (cfg.cipher_type == STRADDLING_CHECKERBOARD && g_alpha == DEFAULT_ALPHABET_SIZE) {
        init_alphabet_adfgvx();
        printf("-type sc: alphabet forced to %d symbols (A..Z + 0..9): %s\n",
            g_alpha, g_idx_to_char_arr);
    }

    // --- Resource Loading ---

    shared.ngram_data = load_ngrams(cfg.ngram_file, cfg.ngram_size, cfg.verbose);

    if (cfg.dictionary_present) {
        load_dictionary(cfg.dictionary_file, &shared.dict, &shared.n_dict_words, &shared.max_dict_word_len, cfg.verbose);
    }

    cribtext[0] = '\0';
    if (cfg.crib_present) {
        if (file_exists(cfg.crib_file)) {
            FILE *fp_crib = fopen(cfg.crib_file, "r");
            fscanf(fp_crib, "%s", cribtext);
            fclose(fp_crib);
            if (cfg.verbose) printf("cribtext = \n\'%s\'\n\n", cribtext);
        } else {
            printf("\nERROR: missing file '%s'\n", cfg.crib_file);
            return 0;
        }
    }


    // --- Execution Flow ---

    printf("\nRNG seed = %u (override with -seed)\n", rng_seed);

    if (cfg.batch_present) {
        if (!file_exists(cfg.batch_file)) {
             printf("\nERROR: missing batch file '%s'\n", cfg.batch_file);
             return 0;
        }

        FILE *fp_batch = fopen(cfg.batch_file, "r");
        char line_buffer[MAX_CIPHER_LENGTH];
        
        printf("\n--- Starting Batch Processing ---\n");

        int n_ciphers = 0;
        while (fgets(line_buffer, sizeof(line_buffer), fp_batch)) {
            // Strip newline
            line_buffer[strcspn(line_buffer, "\r\n")] = 0;
            
            // Skip empty lines
            if (strlen(line_buffer) < 5) continue; 

            n_ciphers ++;
            printf("\nProcessing %d: %s\n", n_ciphers, line_buffer);
            solve_cipher(line_buffer, cribtext, &cfg, &shared, NULL);
        }
        fclose(fp_batch);

    } else {
        // Single Cipher Mode
        if (!file_exists(cfg.ciphertext_file)) {
             printf("\nERROR: missing cipher file '%s'\n", cfg.ciphertext_file);
             return 0;
        }

        // Read the cipher file as the ciphertext, preserving any internal spaces and
        // punctuation (unlike fscanf("%s"), which stops at the first whitespace). By
        // default only the first line is read: stopping at the newline keeps the
        // historical behaviour of ignoring trailing lines (e.g. a "plaintext = ..."
        // annotation). With -multiline the whole file is read and newlines are dropped
        // (not turned into cipher positions), so a ciphertext laid out over several
        // lines -- e.g. a homophonic grid -- is concatenated into one symbol stream.
        // These non-alphabetic characters are kept as positions and carried through the
        // decryption; scoring skips them. Use -skipspaces to drop them entirely.
        FILE *fp_cipher = fopen(cfg.ciphertext_file, "r");
        int ci = 0, ch;
        while ((ch = fgetc(fp_cipher)) != EOF && (cfg.multiline || ch != '\n')
               && ci < MAX_CIPHER_LENGTH - 1) {
            if (ch == '\r' || ch == '\n') continue;
            single_ciphertext_buffer[ci++] = (char) ch;
        }
        single_ciphertext_buffer[ci] = '\0';
        // Trim trailing whitespace so a stray space at the end of the line does not
        // become an extra cipher position. EXCEPTION: the pure-transposition types
        // carry spaces/periods as genuine grid cells, so a trailing space is a real
        // position (e.g. the last column of a columnar grid may end in a space). For
        // those, trimming would shorten the grid and misalign a full-length crib, so
        // it is skipped. (All existing transposition test ciphers are trailing-space
        // free, so this is bit-identical for them.)
        int t = cfg.cipher_type;
        bool space_significant = (t >= TRANSMATRIX && t <= GRILLE) ||
                                 (t >= TRANSCOL_L && t <= TRANSTILE) ||
                                 (t == PERIOD_COLUMN) ||
                                 (t == PERIOD_COLUMN_SPACE) ||
                                 (t == TRANSCOL2_DC);
        if (!space_significant)
            while (ci > 0 && isspace((unsigned char) single_ciphertext_buffer[ci - 1]))
                single_ciphertext_buffer[--ci] = '\0';
        fclose(fp_cipher);

        if (cfg.verbose) printf("ciphertext = \n\'%s\'\n\n", single_ciphertext_buffer);

        solve_cipher(single_ciphertext_buffer, cribtext, &cfg, &shared, NULL);
    }

    // --- Cleanup ---
    free(shared.ngram_data);
    if (shared.dict != NULL) {
        free_dictionary(shared.dict, shared.n_dict_words);
    }

    return 1;
}
#endif // COLOSSUS_NO_MAIN
void solve_cipher(char *ciphertext_str, char *cribtext_str, ColossusConfig *cfg,
    SharedData *shared, SolveResult *result) {

    // Default to "not solved": every early return (transposition dispatch, no
    // periodicities found, no valid configuration) then leaves a correct result.
    if (result) result->solved = false;

    // Tokenized (symbol) input: HOMOPHONIC always, or any cipher run with -delimiter.
    // In this mode delimiters and symbol tokens are significant, so -skipspaces (which
    // would strip them) is suppressed and the ciphertext is decoded into symbol ids.
    bool symbol_mode = (cfg->cipher_type == HOMOPHONIC) || cfg->delimiter_present;

    // -skipspaces: drop spaces/punctuation from the ciphertext entirely, so they
    // are not even carried as transposition positions. Default (flag off) keeps
    // them -- ord() encodes them as negative sentinels that ride through the
    // decryption and are skipped only by scoring.
    if (cfg->skip_spaces && !symbol_mode) {
        int w = 0;
        for (int r = 0; ciphertext_str[r] != '\0'; r++) {
            unsigned char c = (unsigned char) ciphertext_str[r];
            // Keep letters and any other char registered in the active alphabet (the
            // Trifid '+'); drop spaces/punctuation. For the default A..Z alphabet this
            // is byte-for-byte the historical isalpha() filter.
            if (isalpha(c) || (c < 128 && g_char_to_idx[toupper(c)] >= 0))
                ciphertext_str[w++] = ciphertext_str[r];
        }
        ciphertext_str[w] = '\0';
    }

    int cipher_indices[MAX_CIPHER_LENGTH];
    int n_cribs = 0;
    int crib_positions[MAX_CIPHER_LENGTH];
    int crib_indices[MAX_CIPHER_LENGTH];
    SymbolTable symtab;

    // Prepare indices. Letter ciphers get the historical 0..25/sentinel encoding
    // (decode_cipher reproduces ord() byte-for-byte); HOMOPHONIC fills symtab with the
    // distinct ciphertext symbols and emits one symbol id per position.
    int cipher_len = decode_cipher(ciphertext_str, cfg, cipher_indices, &symtab);

    // Process Cribs (Local to this cipher)
    if (strlen(cribtext_str) > 0) {
        if ((int)strlen(cribtext_str) != cipher_len) {
            if (cfg->verbose) printf("Crib length mismatch (Crib: %lu, Cipher: %d). Ignoring crib.\n", strlen(cribtext_str), cipher_len);
            n_cribs = 0;
        } else {
            for (int i = 0; i < cipher_len; i++) {
                if (cribtext_str[i] != '_') {
                    crib_positions[n_cribs] = i;
                    crib_indices[n_cribs] = g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127];
                    n_cribs++;
                }
            }
        }
    }


    // --- TRANSPOSITION CIPHERS ---
    // These are pure transpositions solved by optimization over the transform
    // parameters, not via the keyword/cycleword/period machinery below.
    if (cfg->cipher_type == TRANSMATRIX || cfg->cipher_type == TRANSPEROFFSET) {
        solve_transposition(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSPOSITION) {
        solve_general_transposition(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSCOL || cfg->cipher_type == TRANSCOL2) {
        solve_columnar(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSCOL_L) {
        solve_columnar_track(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSROUTECOL) {
        solve_route_chain(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSTILE) {
        solve_tile(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == PERIOD_COLUMN) {
        solve_period_column(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == PERIOD_COLUMN_SPACE) {
        solve_period_column_space(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == TRANSCOL2_DC) {
        solve_double_transposition(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == POLLUX) {
        solve_pollux(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }
    if (cfg->cipher_type == MORBIT) {
        solve_morbit(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }
    if (cfg->cipher_type == STRADDLING_CHECKERBOARD) {
        solve_straddling_checkerboard(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }
    if (cfg->cipher_type == RAILFENCE) {
        solve_railfence(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == ROUTE) {
        solve_route(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == AMSCO) {
        solve_amsco(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == MYSZKOWSKI) {
        solve_myszkowski(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == REDEFENCE) {
        solve_redefence(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == CADENUS) {
        solve_cadenus(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == NIHILIST) {
        solve_nihilist(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == SWAGMAN) {
        solve_swagman(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == GRILLE) {
        solve_grille(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == INDEP_PERIODIC) {
        solve_indep_periodic(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
        return ;
    }
    if (cfg->cipher_type == HOMOPHONIC) {
        solve_homophonic(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, &symtab);
        return ;
    }
    if (cfg->cipher_type == PLAYFAIR) {
        solve_playfair(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }
    if (cfg->cipher_type == SERIATED_PLAYFAIR) {
        solve_seriated_playfair(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == BIFID) {
        solve_bifid(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == TRIFID) {
        solve_trifid(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == DIGRAFID) {
        solve_digrafid(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == CM_BIFID) {
        solve_cm_bifid(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == HILL) {
        solve_hill(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == PHILLIPS || cfg->cipher_type == PHILLIPS_C ||
        cfg->cipher_type == PHILLIPS_RC) {
        solve_phillips(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == TWO_SQUARE || cfg->cipher_type == TWO_SQUARE_V) {
        solve_twosquare(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == FOUR_SQUARE) {
        solve_foursquare(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == TRI_SQUARE) {
        solve_trisquare(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == ADFGX || cfg->cipher_type == ADFGVX) {
        solve_adfgvx(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == NIHILIST_SUB || cfg->cipher_type == NIHILIST_SUB_NC ||
        cfg->cipher_type == NIHILIST_SUB_M100) {
        solve_nihilist_sub(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == GROMARK || cfg->cipher_type == GROMARK_PERIODIC) {
        solve_gromark(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == NICODEMUS || cfg->cipher_type == NICODEMUS_VARIANT ||
        cfg->cipher_type == NICODEMUS_BEAUFORT) {
        solve_nicodemus(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == BAZERIES) {
        solve_bazeries(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == PORTAX) {
        solve_portax(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == PROGKEY || cfg->cipher_type == PROGKEY_VAR ||
        cfg->cipher_type == PROGKEY_BEAU) {
        solve_progkey(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == SLIDEFAIR || cfg->cipher_type == SLIDEFAIR_VAR ||
        cfg->cipher_type == SLIDEFAIR_BEAU) {
        solve_slidefair(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == INTERRUPTED_KEY || cfg->cipher_type == INTERRUPTED_KEY_VAR ||
        cfg->cipher_type == INTERRUPTED_KEY_BEAU) {
        solve_intkey(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == CONDI) {
        solve_condi(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }

    if (cfg->cipher_type == FRAC_MORSE) {
        solve_fracmorse(ciphertext_str, cribtext_str, cfg, shared,
            cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
        return ;
    }



    // --- POLYALPHABETIC CIPHERS ---
    // Vigenere / Quagmire I-IV / Beaufort / Porta / Autokey* share one model;
    // solve_polyalpha (polyalpha_solver.c) owns POLYALPHA_MODEL and the engine call.
    solve_polyalpha(ciphertext_str, cribtext_str, cfg, shared,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs, result);
}
