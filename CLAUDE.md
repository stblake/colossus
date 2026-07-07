# CLAUDE.md

Guidance for working in this repository.

## Scope

This directory is the **entire project** and the git root. It tracks
`https://github.com/stblake/colossus` (branch `main`). Everything outside this
directory is out of scope — git can't see it, and neither should you. (The parent
folder holds unrelated experiment runs, logs, and candidate dumps; ignore it.)

## What this is

Colossus is a polyalphabetic substitution cipher solver in C by Sam Blake (started 14 July 2023).
It attacks **Vigenère, Gronsfeld, Beaufort, Porta, Quagmire I–IV, Autokey, Progressive
Key, Interrupted Key, and Condi** ciphers (plus
their variants and Beaufort/Porta autokey tableaus), optionally composed with a
transposition stage. The engine is a **stochastic, slippery, shotgun-restarted hill
climber with backtracking**. Cipher conventions follow the American Cryptogram
Association (https://www.cryptogram.org/resource-area/cipher-types/). It exists to
crack the Kryptos sculpture's K1–K4. See `README.md` for the author's full writeup.

## Layout (sources under `src/<cipher-class>/`)

The C sources are grouped by cipher class under `src/`; everything else (build,
tests, tools, data, ciphers) stays at the repo root. All local `#include`s are
**flat** (`#include "foo.h"`, no directory prefix) — the makefile resolves them
with one `-I` per `src/` subdir (the `INCLUDES` variable), so a header is found
regardless of which subdirectory it lives in. Add a new `src/` subdir → add it
to `INCLUDES`.

```
src/core/        # cipher-agnostic engine + shared infrastructure
  colossus.c       # main(): arg parsing, init_config(), solve_cipher() dispatcher
  colossus.h       # shared CORE header: config/ctx/model structs, constants, cipher-type
                   #   codes, globals, inline RNG, and the cipher-PRIMITIVE prototypes
  engine.c/.h      # cipher-agnostic search engine: run_solver(), run_one_config(),
                   #   make_solver_ctx(), the SearchDefaults registry + apply_cipher_defaults()
  scoring.c/.h     # state_score / ngram_score / crib_score, load_ngrams, keyword/cycleword RNG
  parse.c          # parse_cipher_type(): string/int aliases -> cipher-type code
  perioc.c         # estimate_cycleword_lengths(): IoC period estimation (Z-score + threshold)
  optimal_cycleword.c  # derive_optimal_cycleword(): deterministic per-column frequency attack
  dict.c           # dictionary load + word-finding (scores plaintext readability)
  utils.c          # ord/print, decode_cipher/print_cipher (symbol I/O), IoC, chi-squared, etc.

src/polyalphabetic/   # Vigenère family — searched inside POLYALPHA_MODEL
  polyalpha_solver.c/.h  # POLYALPHA_MODEL (vig/quag/beau/porta/autokey) + crib/cycleword helpers
                         #   + solve_polyalpha(); solve_cipher() dispatches the polyalpha types here
  vigenere.c gronsfeld.c beaufort.c porta.c quagmire.c autokey.c   # per-cipher encrypt/decrypt primitives
  gromark.c gromark_solver.c/.h  # Gromark + Periodic Gromark: keyed-alphabet substitution + chain-
                         #   addition running key. Basic = primer pre-pass (10^5 space) then a sigma
                         #   anneal; Periodic = anneal the KEYWORD directly (it derives sigma/primer/
                         #   offsets). Has its own CipherModels, not part of POLYALPHA_MODEL.
  nicodemus.c nicodemus_solver.c/.h  # Nicodemus: periodic Vigenere/Variant/Beaufort substitution
                         #   (one shift per column) composed with a per-block columnar transposition,
                         #   both keyed by one keyword. Solver anneals the COLUMN ORDER and derives the
                         #   per-column shifts by monogram fit (decoupling). Own CipherModel; sweeps
                         #   (period P, block height H).
  progkey.c progkey_solver.c/.h  # Progressive Key (Vig/Var/Beau base): a periodic base cipher
                         #   under a keyword + a per-GROUP constant key drift (the progression
                         #   index). Own CipherModel (NOT in POLYALPHA_MODEL); period brute-forced
                         #   x progression 0..25 enumerated (IoC fails through the drift). For a
                         #   fixed prog, DE-PROGRESSING decouples the keyword -> per-column monogram
                         #   warm start, then n-gram anneal. 3 type codes share the solver.
  intkey.c intkey_solver.c/.h  # Interrupted Key (Vig/Var/Beau base): a periodic keyword whose key
                         #   index RESETS to key letter 1 at break points. Own CipherModel (NOT in
                         #   POLYALPHA_MODEL); IoC fails, so period swept + interruption STRATEGY
                         #   enumerated. ct-interruptor (reset after a ct letter -> breaks known from
                         #   cipher -> decouples like Vigenere, the reliable blind mode), pt-interruptor
                         #   (reset after a pt letter -> causal, EM warm start, FRAGILE), supplied
                         #   -breaks (random/word-division, reliable), joint keyword+break-mask anneal
                         #   (blind random, characterized). Gromark-style pre-pass ranks (period,
                         #   strategy, interruptor); top-K annealed. 3 type codes share the solver.
  condi.c condi_solver.c/.h  # Condi: PLAINTEXT-FEEDBACK substitution over a keyed alphabet sigma --
                         #   idx(ct_i)=(idx(pt_i)+idx(pt_{i-1})+1) mod 26, first letter uses a starter.
                         #   Own CipherModel (NOT in POLYALPHA_MODEL); no period. Shipped attack is a
                         #   free-permutation sigma anneal with the 26 starters enumerated as configs.
                         #   KEY FINDING: the feedback makes the true sigma an ISOLATED NEEDLE (one cell
                         #   swap re-derives the whole downstream decrypt -> no basin/gradient), so NO
                         #   local search cracks it blind at any budget (documented, asserted in the
                         #   solver test). Cribs passed through positionally (substrate for the untapped
                         #   tractable attack: crib-anchored constraint solving of pos[]).

src/transposition/    # pure-transposition solvers + shared helpers
  trans_common.c/.h    # shared transposition-solver helpers: report_transposition(),
                       #   TransKeyOps seed/move, perm_move/seed, sweep no-ops, exact_isqrt
  transpositions.c     # transperoffset() (periodic decimation), transmatrix() (K3-style double rotation)
  transmatrix_solver.c/.h permutation_solver.c/.h columnar_solver.c/.h
  railfence_solver.c/.h route_solver.c/.h amsco_solver.c/.h myszkowski_solver.c/.h
  redefence_solver.c/.h cadenus_solver.c/.h nihilist_solver.c/.h swagman_solver.c/.h grille_solver.c/.h
  columnar_track_solver.c/.h   # transcol-L: columnar + within-column row permutation L (seam best-L);
                               #   also the structural -cribanchored block<->column matcher
  route_chain_solver.c/.h      # transroutecol: fixed read-route global + searched column key (seam best-L)
  tile_solver.c/.h             # transtile: sub-grid h x w tile transposition (joint column-order + tile perm)
  period_column.c              # period_column_transform(): AZdecrypt's "Period column order" -- one
                               #   periodic column-permutation transposition stage (complete + incomplete grids)
  period_column_solver.c/.h    # period-column: DETERMINISTIC EXHAUSTIVE solver -- tries every complete-grid
                               #   stage (depth 1) and every ordered pair (depth 2), n-gram-scores each, keeps
                               #   the global best (reproduces AZdecrypt's stacked "Period column order" solutions);
                               #   honours -readdir tb|bt|both (forward / reversed cipher stream)
  period_column_space_solver.c/.h # period-column-space: space-robust variant for POORLY enciphered text --
                               #   keeps the observed cipher (spaces intact) but searches FROM SCRATCH over
                               #   INSERTING blank/gap cells (dropped-char repair) and DELETING observed cells
                               #   (added-char repair); a restart+anneal edit search whose fitness is the
                               #   exhaustive period-column search (depth-1 inner, depth-<=2 final); -readdir too
  double_transposition_solver.c/.h # transcol2-dc: DOUBLE columnar, DIVIDE-AND-CONQUER (Lasry/Kopal/Wacker
                               #   2014). Scores the SECOND key K2 INDEPENDENTLY of K1 via the Index of
                               #   Digraphic Potential (IDP: undo K2, then greedily reconstruct the residual
                               #   single columnar's column adjacency by digraph fit) -- so the O(K1!*K2!)
                               #   joint search of transcol2 collapses to an O(K2!) IDP climb + an O(K1!)
                               #   single-columnar finish (rotation-resolved from the IDP's own greedy chain).
                               #   Cracks key lengths (20-25) at which transcol2's parallel hill-climb fails.
  # trans_common.c also carries the shared exact-ordering helpers: held_karp_best_path()
  #   (max-weight Hamiltonian path) and seam_best_row_order() (exact best within-column
  #   track order L via a per-row + seam-delta decomposition), plus trans_word_set().

src/polygraphic/      # square/cube/matrix ciphers — each: primitive + a CipherModel solver
  playfair.c playfair_solver.c/.h   # Playfair: 5x5 keyed grid; grid build / prepare / encrypt / decrypt
  bifid.c    bifid_solver.c/.h      # Bifid: side-generic keyed Polybius square build / encrypt / decrypt
  trifid.c   trifid_solver.c/.h     # Trifid: side-generic keyed 3x3x3 cube build / encrypt / decrypt
  hill.c     hill_solver.c/.h       # Hill: matrix multiply / encrypt / decrypt / det+inverse mod 26 (generic k x k)
  phillips.c phillips_solver.c/.h   # Phillips: derive 8 squares from a base / encrypt / decrypt (side-generic, 3 variants)
  twosquare.c twosquare_solver.c/.h # Two-Square: two keyed squares, rectangle rule; horizontal (ACA) + vertical (self-inverse), one solver
  foursquare.c foursquare_solver.c/.h # Four-Square: 2 keyed + 2 fixed-standard squares / encrypt / decrypt (side-generic)
  adfgvx.c   adfgvx_solver.c/.h     # ADFGVX/ADFGX: keyed-square fractionation (reuses bifid square build/inverse) +
                                     #   keyed columnar (reuses decrypt_columnar); coordinate-space, side-generic 5x5/6x6
  nihilist_sub.c nihilist_sub_solver.c/.h # Nihilist Substitution: periodic ADDITIVE over a keyed Polybius square
                                     #   (numeric ciphertext); 3 add conventions (carry/no-carry/mod-100), one solver;
                                     #   square-independent validity reward decouples the additive key (à la ADFGVX)
  bazeries.c bazeries_solver.c/.h    # Bazeries: keyed-square substitution (N spelled out -> ct square,
                                     #   fixed column-major pt square) composed with a digit-grouped reversal
                                     #   transposition, both keyed by ONE number N < 10^6. Solver CLIMBS N's
                                     #   decimal digits (one config per digit count D); a transposition-independent
                                     #   monogram reward decouples the square from the digits. Reuses bifid square build.
  portax.c   portax_solver.c/.h      # Portax: PERIODIC DIGRAPHIC Porta. Plaintext written row-major at width P
                                     #   (= keyword len), rows taken in PAIRS; the vertical pair in column c is
                                     #   enciphered as a unit over a Porta slide by keyword[c] (only its shift key/2
                                     #   matters). Self-reciprocal. Solver anneals the P per-column shifts (cycleword
                                     #   lane, one config per swept period P); a per-column monogram-fit warm start
                                     #   decouples each column (every pair uses ONE column key). Cribs supported.
  slidefair.c slidefair_solver.c/.h  # Slidefair: PERIODIC DIGRAPHIC Vigenere/Variant/Beaufort. Plaintext taken in
                                     #   consecutive DIGRAPHS; digraph i is keyed by keyword[i mod P] over a two-row
                                     #   slide (top = standard alphabet, bottom = a shift alphabet). The pair forms
                                     #   diagonal corners of a rectangle; substitutes are the other two corners (top
                                     #   first), vertical pair -> the pair one column right (decrypt: left). Self-
                                     #   reciprocal but for that vertical step. 3 type codes share the solver. Solver
                                     #   anneals the P per-column key letters (0..25); a per-column monogram-fit warm
                                     #   start decouples each column (every digraph uses ONE column key). Cribs supported.
  seriated_playfair.c seriated_playfair_solver.c/.h  # Seriated Playfair: plain Playfair over a
                                     #   SINGLE 5x5 keyed square, but digraphs are the VERTICAL PAIRS of a two-row
                                     #   seriated layout of period P: in each 2P block, pair j couples block-letter j
                                     #   (top) with j+P (bottom); cipher taken off block by block (top row then bottom).
                                     #   Reuses playfair.c's grid/inverse + keyword build. NO per-column decoupling (one
                                     #   square), so the attack IS Playfair's single-grid anneal with the seriation period
                                     #   SWEPT (one config per P; n-gram score picks it). Needs -logprob; no cribs.
  digrafid.c digrafid_solver.c/.h    # Digrafid: digraphic fractionation over TWO keyed 27-symbol
                                     #   alphabets (A..Z + #) -- a horizontal 3x9 grid + a vertical 9x3 grid. A
                                     #   plaintext digraph -> a 3-digit number (top=col in H, bot=row in V, mid from
                                     #   the 3x3 intersection); each PERIOD-long group's numbers are stacked + Trifid-
                                     #   fractionated, then mapped back to ciphertext digraphs. The two grids are KEYED
                                     #   ALPHABETS (keyword + ascending tail), searched AS SUCH (anneal two short keywords,
                                     #   not free 54-cell permutations) -- drops the blind cliff from ~700 to ~300 letters;
                                     #   period SWEPT (one config per P; n-gram score picks it). Reuses bifid.c's
                                     #   build-inverse. Needs -logprob; no cribs.
  cm_bifid.c cm_bifid_solver.c/.h    # CM Bifid (Conjugated Matrix Bifid): plain Bifid, but the
                                     #   coordinate-pair -> letter recombination uses a SECOND keyed 5x5 square
                                     #   (sq1 fractionates, sq2 recombines; sq1==sq2 reduces to Bifid). JOINT two-
                                     #   square anneal (state = both squares packed back-to-back, perturb one per
                                     #   move) -- no decoupling reward, like Four-Square; period SWEPT via Bifid's
                                     #   estimator (reused, square-agnostic). ODD periods recover (~480 letters);
                                     #   EVEN periods are degenerate ciphertext-only (rows/cols never share an output
                                     #   pair -> transpose-like square ambiguity). Reuses bifid.c. Needs -logprob; no cribs.
  trisquare.c trisquare_solver.c/.h  # Tri-Square: digraphic substitution over THREE independent
                                     #   keyed 5x5 squares. A plaintext digraph (p1 in sq1, p2 in sq2) -> a ciphertext
                                     #   TRIGRAPH (3:2 expansion): c0 = any letter in p1's COLUMN of sq1, c1 =
                                     #   sq3[row(p1)][col(p2)], c2 = any letter in p2's ROW of sq2. Encryption is
                                     #   POLYPHONIC (c0/c2 random members); decryption is exact for any choice. JOINT
                                     #   three-square anneal (state = the 3 squares packed back-to-back, perturb one per
                                     #   move) -- no decoupling reward, like Four-Square/CM-Bifid. NO period. The 3:2
                                     #   length change is handled like ADFGVX (cipher in scratch, plaintext length n=2M
                                     #   passed to the engine). Reuses bifid_build_inverse. Needs -logprob; no cribs;
                                     #   recovers from ~500 plaintext letters (the polyphonic letters spread the full
                                     #   alphabet over every position, so it is EASIER than Four-Square despite 75 cells).
  fracmorse.c fracmorse_solver.c/.h  # Fractionated Morse: plaintext -> Morse with single 'x'
                                     #   letter separators -> the {DOT,DASH,X} stream grouped into TRIGRAPHS (rank
                                     #   9a+3b+c, xxx excluded) -> ciphertext via a keyed 26-letter alphabet sigma
                                     #   (rank -> letter). A length CHANGE (N plaintext -> C ciphertext letters). NO
                                     #   period. The alphabet is an ACA KEYED alphabet, searched AS SUCH (keyword prefix
                                     #   + ascending tail, fracmorse_move_seq, the length-26 twin of digrafid's) -- not a
                                     #   free 26! permutation -- so it tracks the keyword and cracks the short (~110-150-
                                     #   letter) ACA range. The decode length VARIES per key but the engine scores a fixed
                                     #   length, so the decrypt hook decodes to plaintext (a filler for each invalid Morse
                                     #   token) and CYCLICALLY TILES it to the ciphertext length C; a MORSE-VALIDITY reward
                                     #   (fraction of tokens that are legal codewords) is folded into score_adjust (à la
                                     #   ADFGVX's IoC term). Needs -logprob; no cribs.
  pollux.c pollux_solver.c/.h        # Pollux: plaintext -> Morse with single 'x' letter separators
                                     #   (xx between words) -> the {DOT,DASH,X} stream mapped to a DIGIT stream by a
                                     #   KEY assigning each digit 0..9 to one of the three symbols (polyphonic encode:
                                     #   a random assigned digit per symbol; deterministic decode). A length CHANGE (C
                                     #   digits -> N letters). NO period. The ENTIRE keyspace is 3^10 = 59049 maps, so
                                     #   -- like period_column -- the solver is DETERMINISTIC and EXHAUSTIVE (a needle:
                                     #   one digit re-assignment re-parses the whole stream, no gradient), NOT a
                                     #   CipherModel/anneal: pollux_search decodes+scores every key and keeps the global
                                     #   best. Fitness = mean n-gram (length-fair) + a Morse-VALIDITY reward (à la
                                     #   fracmorse; illegal xxx runs / >4-symbol codewords count against it). Reuses a
                                     #   small self-contained Morse table (fracmorse.c left byte-identical). Rides the
                                     #   reward-only quadgram table (NO -logprob -- the validity weight is tuned to that
                                     #   scale; -logprob degrades short-text recovery). Recovers from ~24 letters, far
                                     #   below the ACA 80-100 range. No cribs; digit input parsed from ciphertext_str.

src/substitution/     # monoalphabetic / homophonic substitution solvers
  indep_solver.c/.h homophonic_solver.c/.h   # each: a CipherModel + solve_<type>()

makefile
README.md  LICENSE
example.sh           # canonical usage example
cipher.txt  crib.txt # sample ciphertext + crib
tools/homophonic_gen.c       # standalone homophonic-cipher test-data generator (make homophonic_gen)
tools/playfair_gen.c         # standalone Playfair test-data generator (make playfair_gen)
tools/bifid_gen.c            # standalone Bifid test-data generator (make bifid_gen)
tools/trifid_gen.c           # standalone Trifid test-data generator (make trifid_gen)
tools/hill_gen.c             # standalone Hill test-data generator (make hill_gen)
tools/gronsfeld_gen.c        # standalone Gronsfeld test-data generator (make gronsfeld_gen)
tools/phillips_gen.c         # standalone Phillips test-data generator (make phillips_gen)
tools/twosquare_gen.c        # standalone Two-Square test-data generator (make twosquare_gen)
tools/foursquare_gen.c       # standalone Four-Square test-data generator (make foursquare_gen)
tools/nihilist_sub_gen.c     # standalone Nihilist Substitution generator (make nihilist_sub_gen)
tools/gromark_gen.c          # standalone Gromark / Periodic Gromark generator (make gromark_gen)
tools/nicodemus_gen.c        # standalone Nicodemus generator (make nicodemus_gen)
tools/bazeries_gen.c         # standalone Bazeries generator (make bazeries_gen)
tools/portax_gen.c           # standalone Portax generator (make portax_gen)
tools/slidefair_gen.c        # standalone Slidefair generator (make slidefair_gen)
tools/seriated_playfair_gen.c # standalone Seriated Playfair generator (make seriated_playfair_gen)
tools/digrafid_gen.c         # standalone Digrafid generator (make digrafid_gen)
tools/cm_bifid_gen.c         # standalone CM Bifid generator (make cm_bifid_gen)
tools/trisquare_gen.c        # standalone Tri-Square generator (make trisquare_gen)
tools/progkey_gen.c          # standalone Progressive Key generator (make progkey_gen)
tools/intkey_gen.c           # standalone Interrupted Key generator (make intkey_gen)
tools/condi_gen.c            # standalone Condi generator (make condi_gen)
tools/fracmorse_gen.c        # standalone Fractionated Morse generator (make fracmorse_gen)
tools/period_column_gen.c    # standalone Period column order generator (make period_column_gen; space-aware)
tools/pollux_gen.c           # standalone Pollux generator (make pollux_gen; assignment string + RNG seed)
tools/double_transposition_gen.c # standalone double columnar transposition generator (make double_transposition_gen)
english_quadgrams.txt        # n-gram table (quadgrams); english_quintgrams.txt (5-grams) optional, with -logprob
OxfordEnglishWords.txt       # default dictionary (auto-loaded if present in cwd)
ciphers/kryptos/     # K1–K4 ciphertexts + run scripts
ciphers/tests/       # per-cipher test cases (cipher + expected solution)
```

## Build

```bash
make            # gcc -Wall -O3; builds ./colossus
make clean
```

Two caveats:
- The active `CC` line does **not** include `-lm`. Links on macOS (clang folds libm
  into libc) but fails on Linux — add `-lm` there.
- `make` also runs `cp colossus ..` (and `../quagmire`), copying the binary
  *outside* this directory. That predates the isolation of this repo; the in-tree
  `./colossus` is the one that matters here.
- The translation-unit list lives in makefile variables, each path prefixed with the
  source-class dir vars (`$(CORE)`, `$(POLY)`, `$(TRANS)`, `$(GRAPH)`, `$(SUBST)`):
  `PRIMITIVES` (the cipher decrypt math + utils), `SOLVERS` (the engine/scoring/trans_common
  core + the per-cipher-type solver modules), and
  `SOLVER_SRC = $(PRIMITIVES) $(SOLVERS) $(CORE)/colossus.c` (used by both `all` and the
  `testopt` harnesses). Every compile also passes `$(INCLUDES)` (one `-I` per `src/` subdir).
  Add a new solver module to `SOLVERS` with its `src/<class>/` prefix.

`make test` builds and runs the framework-free unit tests in
`tests/test_transpositions.c` (the transposition primitives, including the columnar
`decrypt_columnar`: known-answer, round-trip across complete/incomplete grids and both
read directions, double-columnar composition), `tests/test_ciphers.c`,
`tests/test_optimal_cycleword.c`, and `tests/test_playfair.c` (the Playfair primitives:
the Wikipedia known-answer vector, the prepare/X-insertion rules, encrypt/decrypt
round-trips over random grids, and the cyclic row/column key-square equivalence), and
`tests/test_bifid.c` (the Bifid primitives: the Wikipedia known-answer vector, the keyed-
square build, encrypt/decrypt round-trips over random 5x5 and 6x6 squares × random
lengths/periods incl. incomplete blocks, and the period-1 identity), and
`tests/test_trifid.c` (the Trifid primitives: the Wikipedia known-answer vector — two
groups, AIDET→FMJFV and OILEC→OISSU, over the 27-symbol cube — the keyed-cube build,
encrypt/decrypt round-trips over random 3x3x3 cubes and side-generic 2x2x2/4x4x4 cubes ×
random lengths/periods incl. incomplete blocks, and the period-1 identity), and
`tests/test_hill.c` (the Hill primitives: the Wikipedia known-answer vector — key
GYBNQKURP, ACT→POH and CAT→FIN — the modular inverse over all residues mod 26, the
matrix inverse by M·inv==I and inv(inv)==M plus singular rejection, encrypt/decrypt
round-trips over random invertible keys for block sizes k=2..5 × random lengths incl.
non-multiples of k, and the k=1 edge cases), and
`tests/test_gronsfeld.c` (the Gronsfeld primitives: a known-answer vector pinning the
numeric-key convention — HELLOWORLD + key 12345 → IGOPTXQUPI — plus the zero-shift
identity and mod-26 wrap, encrypt/decrypt round-trips over random digit keys × random
lengths, exact agreement with `vigenere_*` fed the same digits as its cycleword, and the
keylen-1 / over-long-key edge cases), and
`tests/test_phillips.c` (the Phillips primitives: the full ACA `DIAGONALS` worked-example
known-answer vector — the 81-letter `squares one…is forty` → `KZWLY…GREYXO` — plus an
assertion pinning `phillips_build_squares`' 8-square Row table cell-for-cell against the
ACA's printed squares #1–#8, encrypt/decrypt round-trips over random base squares × random
lengths for ALL THREE variants and a side-generic 6x6, and the documented structural facts
— the column-rotation symmetry, its row-rotation dual, the #1≡#5 / #2≡#8 cyclic equivalences,
and that the 8 derived squares are distinct grids), and
`tests/test_twosquare.c` (the Two-Square primitives: the ACA worked-example known-answer
vector for the horizontal type — the two printed squares, `anothe…px` → `IRRTEHMKGIMEQGRUNMMZSV`
— and the Wikipedia worked-example for the vertical type — EXAMPLE/KEYWORD → `HEDLXW…` — plus
encrypt/decrypt round-trips over random squares × random lengths for BOTH arrangements and a
side-generic 6x6, the vertical type asserted SELF-INVERSE, and the documented transparencies
asserted directly — horizontal same-row → reversed pair, vertical same-column → unchanged), and
`tests/test_foursquare.c` (the Four-Square primitives: the Wikipedia worked-example known-answer
vector — keyed squares EXAMPLE/KEYWORD over the fixed standard squares, `HELP…` → `FYGMKY…` —
encrypt/decrypt round-trips over random keyed squares × random lengths and a side-generic 6x6,
and a degenerate identity-square check pinning the exact coordinate algebra), and
`tests/test_adfgvx.c` (the ADFGVX/ADFGX primitives: a HAND-COMPUTED known-answer vector pinning
the whole convention end to end — identity 5x5 square, `ATTACK` + keyword `KEY` → `AGADAGAFGGAX`
— the label tables, the K=1 (identity-columnar) edge case, and encrypt/decrypt round-trips over
random squares × random column counts × random lengths incl. ragged grids and both read
directions, for both the 5x5 (ADFGX) and the side-generic 6x6 (ADFGVX, 36 cells) squares), and
`tests/test_nihilist_sub.c` (the Nihilist Substitution primitives: a HAND-COMPUTED known-answer
vector pinning ALL THREE conventions at once on the same plaintext/key/square — identity 5x5,
`ZAEZ` + key `ZK` → carry `110 36 70 80` / no-carry `0 36 60 70` / mod-100 `10 36 70 80`, the
carry-triggering positions making the conventions diverge so a mix-up is caught — the validity
predicate asserted vs the legal set, and per-convention stress: encrypt/decrypt round-trips over
random squares × keys × lengths × periods (incl. p=1, p>len, incomplete), a keyed-label round-trip
that also asserts the label-keyed cipher equals the relabelled fixed-label cipher, and a
side-generic 6x6), and
`tests/test_gromark.c` (the Gromark / Periodic Gromark primitives: the TWO ACA worked-example
known-answer vectors pinned cell-for-cell — keyword `ENIGMA`, basic primer `23452`,
`thereare…` → `NFYCK…`, and periodic period-6 `wintry…` → `RHNAAX…` — the K2M mixed-alphabet builder
(`ENIGMA` → `AJRXEBKSYGFPVIDOUMHQWNCLTZ`), the chain-addition rule for P=5 and P=6, the periodic
primer/offset derivation pinned (ranks `264351`, offsets `4 21 13 9 17 0`), encrypt/decrypt round-
trips over random alphabets × primers × lengths for both variants, the identity-alphabet reduction
to a pure chain-shift, and periodic-with-zero-offsets == basic), and
`tests/test_nicodemus.c` (the Nicodemus primitives: a HAND-COMPUTED known-answer vector pinning all
three substitution conventions at once on the same plaintext/keyword — `ATTACKATDAWN`, keyword `KEY`,
H=2 → VIG `XGKKRIXAKKBL` / VARIANT `PYQQVMPSQQFP` / BEAU `LCKKFOLIKKVL` — `nicodemus_key_from_keyword`
(shifts + stable-argsort order incl. repeated keyword letters), encrypt/decrypt round-trips over
random orders/shifts × lengths × P × block heights — incl. ragged final blocks, H=1 and the
single-block degenerate — for all three variants, and agreement of the per-column substitution with
`vigenere_decrypt`/`beaufort_decrypt` fed the same shift as a length-1 cycleword), and
`tests/test_bazeries.c` (the Bazeries primitives: the ACA worked-example known-answer vector pinning
the whole convention end to end — plaintext `simplesubstitution…`, `N=3752` → ciphertext `ACYYU…GQGCI`,
the spelled-out keyed square `THREOUSANDVFIYWBCGKLMPQXZ`, the digits `3,7,5,2`, and the intermediate
reversed-groups string `missbuselp…` — the column-major pt vs row-major ct substitution convention
(identity-square check + fsub/invsub mutual inverse), the transposition asserted an involution over
random digit patterns, encrypt/decrypt round-trips over random N × lengths (incl. the 150–250 band),
and edge cases — a 0 digit, a 1-digit key, ragged final groups), and
`tests/test_portax.c` (the Portax primitives: the ACA worked-example known-answer vectors — the
mini pairs `IN→JL`/`NO→UA`/`NA→DB` (key U/V) and `TA→NM`/`BG→QH` (key E, same-column), plus the
keyword `EASY` end-to-end `THEEARLYBIRDGETSTHEWORMX → NIJAMPBGQCWKHQJEUIKYMPAT` — the pair operation
asserted an involution over all (s, a, b), the self-reciprocal `decrypt == encrypt` and the
shift/key-letter forms agreeing over random keys × lengths, per-column independence, a ragged final
block (lone top letters pass through), and edge cases — `P=1`, `P>len`, a single pair), and
`tests/test_progkey.c` (the Progressive Key primitives: the ACA worked-example known-answer vector —
Vigenère, key `GRAPEFRUIT`, P=10, prog=1, `THISCIPHERCANBEUSEDWITHANYOFTH → ZYIHGNGBMKJSORJAKZMQQMJRTFHBDC`
with the per-group drift A,B,C — encrypt/decrypt round-trips over random keyword × prog × length ×
period for all three bases (incl. ragged final group, P=1, over-long key), agreement of
`progkey_encrypt` with an INDEPENDENT two-pass reference built from `vigenere_*`/`beaufort_*` (keyword
pass then a full-length drift-key pass), the prog=0 degeneration to a plain periodic Vigenère/Variant,
and `progkey_deprogress` inverting the drift pass exactly), and
`tests/test_slidefair.c` (the Slidefair primitives: the ACA worked-example known-answer vectors — the
mini pairs for all three variants (key B: `ca→ZD/BB/BZ`, `de→EF/FC/XY`, `de` being the vertical/same-
column case), plus the full Vigenère example keyword `DIGRAPH`,
`THESLIDEFAIRCANBEUSEDWITHVIGENEREVARIANTORBEAUFORT → EWKMCRNUAFCXTJ…` — each encrypt pair asserted to
invert under decrypt over every (key, variant, p1, p2), decrypt(encrypt) == identity over random keys ×
lengths × all three variants (incl. odd length → lone final letter passthrough), per-column independence,
and edge cases — `P=1`, `P>ndigraphs`, a single digraph), and
`tests/test_seriated_playfair.c` (the Seriated Playfair primitives: the ACA worked-example known-answer
vector — square `LOGARITHMBCDEFKNPQSUVWXYZ`, period 6,
`comequicklyweneedhelpimmediatelytom` → prepared `COMEQU…XELPIM…` (one null splits the `e/e` vertical
pair) → cipher `NLBCSPCDFGXZQQCDCMGCGQTBHCFTRHFGWHGB` — the prepare/encrypt/decrypt verified cell for cell;
a targeted null-insertion check (filler + the `X→Q` alt rule); a `P=1` equivalence asserting
`seriated_decrypt == playfair_decrypt` (block 2 = consecutive pairs); and encrypt/decrypt round-trips over
random grids × plaintexts × periods incl. ragged buffers (lone top letters pass through) and `P>length`), and
`tests/test_digrafid.c` (the Digrafid primitives: the TWO ACA worked-example known-answer vectors — grids
from keywords `KEYWORD` (H, 3x9) / `VERTICAL` (V, 9x3), `THISISTHEFORESTPRI` → period 3 `HJMXWSWJADWGFCSPYI`
and period 4 `HJTKVHYUFFWDSQYPRI` (the period-4 case a ragged final group of one digraph) — the two grids
asserted cell-for-cell, the `P=1` identity, an odd-length lone-trailing-letter passthrough, and
decrypt/encrypt round-trips over random independent grids × even lengths × periods incl. ragged blocks and
`P>digraph count`), and
`tests/test_cm_bifid.c` (the CM Bifid primitives: the ACA worked-example known-answer vector — squares
`EXTRAKLMPOHWZQDGVUSIFCBYN` (pt) / `NCDRSOBFQUVAGPWEYHMXLTIKZ` (CT), `ODDPERIODSAREPOPULAR` period 7 →
`FANXZEXFENUKKRBYNKAK` — the **`sq1==sq2` ⇒ equals `bifid_encrypt`/`bifid_decrypt`** invariant over random
squares/lengths/periods (anchors the new primitive to the proven Bifid one), decrypt/encrypt round-trips over
random INDEPENDENT square pairs × lengths × periods incl. ragged final blocks and `P>len`, a period-1
monoalphabetic-map check (`x → sq2[pos1(x)]`, identity when `sq1==sq2`), and a 6x6 side-generic round-trip), and
`tests/test_trisquare.c` (the Tri-Square primitives: the ACA worked-example known-answer vector pinned as a
DECRYPT-only test — squares `NSFMUOAGPWVBHQXECIRYLDKTZ` / `READINGBCFHKLMOPQSTUVWXYZ` / `PASTINOQRMLYZUEKXWVBHGFDC`,
`RHLQXR…AAABFZ` → `THREEKEYSQUARESUSEDX` (the ACA CT's polyphonic clerk choices are not reproduced, but
decryption is deterministic) — the decrypt/encrypt round-trip over random independent square triples × lengths
(incl. odd, lone-trailing passthrough) and a side-generic 6x6, a cipher-length check (`3*(len/2)+len%2`), and the
**POLYPHONIC-INVARIANCE** structural test: for a digraph it enumerates all 25 `(c0,c2)` column/row alternatives
and asserts every one decrypts to the SAME digraph — the Tri-Square analogue of the Two-Square transparency /
Four-Square identity-algebra structural check), and
`tests/test_intkey.c` (the Interrupted Key primitives: the ACA worked-example known-answer vector via the mask
form — keyword `ORANGE`, the word-division break lengths 4,6,2,3,4,3,1,1,5,1,2,3,5, `THISCIPHER…PERIODICS` →
`HYIFQZPUKV…ICUIPY` — encrypt asserted cell-for-cell + decrypt round-trip; encrypt/decrypt round-trips for the
ct / pt / mask forms over random keyword × length × period × interruptor and all three bases (incl. ragged, P=1,
len=1); agreement with an INDEPENDENT reference (raw shift formulas + a hand-rolled reset walk, sharing no code
with intkey.c) for all three forms; the ct consistency check (`decrypt_ctint` == `decrypt_mask` over
`build_mask_ct`); and edge cases — an interruptor that never occurs (pt → plain periodic base cipher), one that
occurs at every position (→ monoalphabetic on the first key letter), and P=1), and
`tests/test_condi.c` (the Condi primitives: the ACA worked-example known-answer vector — keyword `STRANGE`,
shifted σ `VWXYZSTRANGEBCDFHIJKLMOPQU`, starter 25, `OURS… → MORC…` — encrypt asserted cell-for-cell +
decrypt round-trip; agreement with an INDEPENDENT reference following the ACA verbal description via a
linear alphabet search (no shared code) over random keyword × starter × length; encrypt/decrypt round-trips
over random σ × starter × length incl. len=1; and edge cases — the starter-only len=1, `starter=0` self-
encipherment, the offset-wrap self-encipherment (a previous letter at the last σ position → shift 0), and a
cyclically shifted keyed alphabet), and
`tests/test_fracmorse.c` (the Fractionated Morse primitives: three HAND-COMPUTED known-answer vectors under
the single-`x` convention — keyword `MORSE`, `SOS → MWLR`, `EE → B`, `ET → C` — pinning the Morse table, the
trigraph rank `9a+3b+c`, the padding and the keyed-alphabet map; agreement with an INDEPENDENT in-test
reference encrypt (a fresh implementation, own Morse table) over random plaintexts/keys; encrypt→decrypt
round-trips == identity over random permutations σ × lengths (incl. the three padding residues, single-letter,
and all-`E`/all-`T` extremes); and the invalid-token path — a run of >4 dots is not a legal codeword → a
filler, counted as not valid), and
`tests/test_period_column.c` (the Period column order primitive: a tiny hand-computed KAT pinning the
column-visit/fill/readout convention (`[0..5]` dx3 p2 → `0,2,1,3,5,4`) and its `utp` inverse; the REAL
AZdecrypt worked example as a two-stage composition KAT — the 168-char spaced ciphertext decrypts via
`[4x42,TP,P:3]` then `[56x3,UTP,P:2]` to `I LIKE KILLING PEOPLE …` and the inverse re-encrypts it, spaces
included; exact `utp==0`→`utp==1` round-trip over random COMPLETE grids × periods (incl. period 1, `dx==len`,
and negative space sentinels); and multiset-preservation over random INCOMPLETE grids — with a printed note
that per-stage `utp` inversion holds for only ~16% of incomplete cases, the reason the solver restricts to
complete grids), and
`tests/test_double_transposition.c` (the double-transposition divide-and-conquer core: the bigram table, a
double encrypt→decrypt round-trip, and the IDP's paper properties — SELECTIVITY (the true `K2` out-scores 400
random keys and peaks at the true `L1` under the per-edge normalization), MONOTONICITY (the IDP degrades as
`K2` is progressively perturbed, Fig 2), and the rotation-resolved greedy `K1` reconstruction from the true
`K2`), and
`tests/test_pollux.c` (the Pollux primitives: the ACA worked-example DECODE KAT — the published 39-digit
ciphertext decodes to `LUCKHELPS` under the ACA key — a MORSE-STREAM KAT that encodes `LUCK HELPS` with
the word space and maps the polyphonic digits back through the key to reproduce the exact 39-symbol Morse
stream (exercising the `xx` word divider), encode→decode round-trips == identity over random keys ×
plaintexts × lengths, and the edge paths — an illegal `xxx` run flagged invalid, a codeword > 4 symbols →
filler, a key missing an element rejected by encode, and single-letter recovery).
`make testopt` additionally runs the in-process solver regressions
`tests/test_solver.c` (polyalphabetic), `tests/test_playfair_solver.c` (Playfair:
validates the per-type schedule registry, asserts an 800-char capability floor, and
prints recovery vs ciphertext length to characterize the short-text cliff),
`tests/test_bifid_solver.c` (Bifid: registry validation, the period-estimator top-K hit
rate, a capability floor with the period *estimated* end-to-end, and the length cliff),
`tests/test_trifid_solver.c` (Trifid: the same four checks over the 27-symbol cube),
`tests/test_hill_solver.c` (Hill: registry validation, block-size *selection* with
k swept, a k=2 and a k=3 capability floor, and the k=2 length cliff), and
`tests/test_gronsfeld_solver.c` (Gronsfeld: confirms it has *no* registry entry and rides
the polyalpha defaults, a capability floor with the period *estimated* end-to-end and again
pinned, and the length cliff), and
`tests/test_phillips_solver.c` (Phillips: registry validation for all three variant types
plus a non-registry type left untouched, a ~760-char capability floor for EACH variant — Row
through the registry default, Column/Row-Column at an explicit budget — and a Row length
cliff showing recovery from ~200 characters),
`tests/test_twosquare_solver.c` (Two-Square: registry validation for both arrangements plus a
non-registry type left untouched, a ~900-char capability floor — horizontal through the
registry default, vertical at an explicit budget — a horizontal length cliff, and the
extra-thorough additions: a transparency-rate measurement (~20% same-row digraphs) and a
multi-keyword sweep reporting mean/worst recovery), and
`tests/test_foursquare_solver.c` (Four-Square: registry validation, a ~900-char capability
floor with a PER-SQUARE recovery breakdown — even positions decrypt through the upper-right
square, odd through the lower-left — a length cliff, and a multi-keyword sweep), and
`tests/test_adfgvx_solver.c` (ADFGX/ADFGVX: registry validation for both types plus a
non-registry type left untouched; an ADFGX capability floor with the column count K pinned;
a BLIND K-selection test — K swept, the solver must report the true K, exercising the IoC
decoupling term; an ADFGX length cliff; an ADFGX multi-keyword sweep (mean/worst); and an
ADFGVX 6x6/36-symbol capability floor over the digit-bearing alphabet), and
`tests/test_nihilist_sub_solver.c` (Nihilist Substitution, run PER CONVENTION: registry validation
for all three codes plus a non-registry type; a period-estimator top-K hit rate measured separately
per convention; a capability floor + length cliff per convention (period pinned, ~250/400 chars);
a BLIND-period carry solve — period estimated end-to-end, the reported period asserted a multiple
of the true one; a carry multi-keyword sweep (mean/worst); and a carry keyed-label end-to-end solve
recovered as the relabelled square), and
`tests/test_gromark_solver.c` (Gromark / Periodic Gromark: registry validation for both codes plus a
non-registry type; a basic-Gromark primer pre-pass hit-rate (the true primer in the top-K vs length,
over several keyword/primer pairs); a basic capability floor + length cliff (blind, the true primer
recovered end-to-end at ~120/150/200 chars); and a Periodic Gromark blind solve with the period
swept, the reported period asserted == the true one), and
`tests/test_nicodemus_solver.c` (Nicodemus: registry validation for all three codes plus a
non-registry type; a per-variant capability floor (~300 chars, P/H pinned); a length cliff
(recovery from ~120 chars); a multi-keyword sweep (mean/worst); and two blind solves — P swept
with H pinned (the reported P asserted == the true one) and H swept with P pinned — validating each
sweep axis. Also the basis the `SearchDefaults` schedule was tuned against), and
`tests/test_bazeries_solver.c` (Bazeries: registry validation plus a non-registry type; a capability
floor over the ACA 150–250-letter band across several numbers (D pinned); a length cliff; a
multi-number sweep (mean/worst); a BLIND digit-count solve — D swept, the reported digit count
asserted == the true one; and a per-scheme pass running the same cipher under `-method`
anneal / shotgun / pso, reporting recovery + time for each (the data the schedule is tuned against)), and
`tests/test_portax_solver.c` (Portax: registry validation plus a non-registry type; a capability floor
across several keywords (period pinned); a length cliff (recovers cleanly from ~70 chars — the
per-column monogram warm start makes it strong); a multi-keyword sweep (mean/worst); a BLIND period
solve — P swept, the reported period asserted == the true one; and a per-scheme pass under `-method`
anneal / shotgun / pso. Rides the reward-only quadgram table — no `-logprob` needed), and
`tests/test_progkey_solver.c` (Progressive Key: registry validation for all three codes plus a
non-registry type; a capability floor over the ACA ~150-letter band for EACH base (P/prog pinned); a
length cliff; a multi-keyword sweep (mean/worst); a BLIND period solve (P swept, prog pinned, the
reported P asserted == the true one); a BLIND progression solve (prog swept, P pinned, the reported
prog asserted == the true one); and a per-scheme pass under `-method` anneal / shotgun / pso. Rides
the reward-only quadgram table — no `-logprob` needed; recovers cleanly from very short text since the
de-progressed columns are a pure Vigenère), and
`tests/test_slidefair_solver.c` (Slidefair: registry validation for all three codes plus a non-registry
type; a capability floor across all three variants (period pinned); a length cliff (recovers ~100% from
~50 chars — the per-column monogram warm start makes it strong); a multi-keyword sweep (mean/worst); a
BLIND period solve (P swept, the reported P asserted == the true one); and a per-scheme pass under
`-method` anneal / shotgun / pso. Rides the reward-only quadgram table — no `-logprob` needed), and
`tests/test_seriated_playfair_solver.c` (Seriated Playfair: registry validation + a non-registry type;
a capability floor across several keywords (period pinned, ~500 chars, `-logprob` — keywords chosen to
avoid the rare-letter X ambiguity that pins an unlucky grid at ~92%); a length cliff (recovers from
~250 chars); a multi-keyword sweep (mean/worst); a BLIND period solve (P swept over a bounded range, the
reported P asserted == the true one); and a per-scheme pass under `-method` anneal / shotgun / pso. The
basis the `SearchDefaults` `6x400000` schedule is tuned against), and
`tests/test_digrafid_solver.c` (Digrafid: registry validation + a non-registry type; the keyed-alphabet move
INVARIANT (`digrafid_move_seq` keeps a permutation + sorted tail over 400×500 random move chains); the PERIOD
ESTIMATOR in isolation (true period in the top-K per-lane IoC across periods × lengths); a capability floor
across keywords (period pinned, **~420 chars — BELOW the old free-permutation cliff**, `-logprob`); a length
cliff showing the **new ~300-char floor** (~9% at 240, ~100% by 320, vs the old ~700); a multi-keyword sweep
(mean/worst); a BLIND period solve (period estimated over a bounded scan, the reported P asserted == the true
one); and a per-scheme pass under `-method` anneal / shotgun / pso. The basis the `SearchDefaults`
`48x150000` / `inittemp 0.30` keyed-alphabet schedule is tuned against), and
`tests/test_cm_bifid_solver.c` (CM Bifid: registry validation + a non-registry type; the PERIOD ESTIMATOR in
isolation (Bifid's columnar-IoC, reused unchanged — true period in the top-5 for BOTH parities, estimation
being parity-independent); a capability floor at ~520 chars (odd period 7, `-logprob`) shown alongside the
**ODD-vs-EVEN contrast** on the same cipher (P=7 ~100% vs P=6 ~noise-floor — the documented even-period square
ambiguity, printed not asserted since it is fragile-not-impossible); a length cliff (~400 odd-period cliff,
reliable from ~480); a multi-keyword sweep over ODD periods (mean/worst); a BLIND period solve (P estimated
over a bounded scan 2..9, the reported P asserted == the true one); and a per-scheme pass under `-method`
anneal / shotgun / pso. The basis the `SearchDefaults` `8x400000` two-square schedule is tuned against — and
the source of the ODD-PERIOD finding: with EVEN P the rows-then-cols re-pairing splits into pure-row / pure-col
output pairs, so a fractionation row and column never share a pair, leaving a transpose-like square ambiguity
no budget escapes; ODD P mixes a row/col at the boundary pair and recovers cleanly), and
`tests/test_trisquare_solver.c` (Tri-Square: registry validation + a non-registry type; a capability floor at
~900 plaintext solved through the **registry default** (`12x500000`) with a first-of-pair/second-of-pair
recovery breakdown (even positions decrypt through sq1, odd through sq2 — both must clear 0.99); a length cliff
(`8x400000`, lengths 400/600/900/1300 — the ~500-char floor visible, only the longest asserted); and a
multi-keyword sweep (mean/worst). The polyphonic c0/c2 letters spread the full alphabet over every position, so
recovery is EASIER than Four-Square despite the 75-cell state (reliable from ~500 letters, a 300-400 cliff) — the
basis the `SearchDefaults` `12x500000` three-square schedule is tuned against; because plant()'s encode is
polyphonic (RNG-driven) it seeds the RNG for a reproducible cipher), and
`tests/test_intkey_solver.c` (Interrupted Key: registry validation for all three codes + a non-registry type; a
**CT capability floor** (blind interruptor, P pinned) across bases/keywords — ~100% from ~40 letters, the reliable
workhorse — + a length cliff; a **BLIND period** solve (P swept, reported P asserted == true); a **BLIND
interruptor+scheme** solve (reported letter AND scheme asserted == true); a multi-keyword sweep (mean/worst, CT); a
**supplied-breaks** (random/word-division) solve via `-breaks` (asserted); a **PT** capability floor asserted only
with the interruptor PINNED, plus a printed **fragility characterization** (blind pt mispicks the interruptor / hits
rugged basins — seed/content sensitive, not asserted); a **JOINT** blind-random characterization (printed, a 20-char
crib lifting it ~29%→~64%); and per-scheme calibration under `-method anneal/shotgun/pso`. The basis the
`SearchDefaults` `6x8000` schedule is tuned against, and the source of the CT-reliable / PT-fragile finding),
and `tests/test_condi_solver.c` (Condi: registry validation (`6x60000`) + a non-registry type left untouched;
**ASSERTS the NEEDLE** — on a real 360-letter solve the true σ out-scores its best single-swap neighbour by a
clear margin AND its neighbours sit near the random floor (the regression guard on the finding that the
plaintext feedback leaves no gradient); and characterizes, printed-not-asserted, that a blind solve does not
recover (~6%), that pinning the starter + a 12×300000 budget still does not (~4% — no budget escapes the
needle), and the per-scheme (anneal / shotgun / pso) behaviour. The source of the Condi needle finding), and
`tests/test_fracmorse_solver.c` (Fractionated Morse: registry validation (`16x120000`) + a non-registry type
left untouched; the keyed-alphabet move INVARIANT (`fracmorse_move_seq` keeps a permutation + sorted tail over
long move chains); a capability floor across keywords at ~150 chars (the ACA range, recovers ~100%); a length
cliff (~97% at 90, ~100% from ~150); a multi-keyword sweep (mean/worst); and per-scheme calibration under
`-method anneal/shotgun/pso` (PSO bounded — its swarm × particle × refine cost makes a large iteration budget
pathological, so the test clamps particles/refine + a tiny budget; only anneal is asserted). The basis the
`SearchDefaults` `16x120000` / `inittemp 0.30` schedule is tuned against), and
`tests/test_period_column_solver.c` (Period column order: depth-1 exact recovery across widths/periods/
directions; depth-2 exact recovery (the AZdecrypt two-stage scenario); the REAL AZdecrypt worked example
asserted **end to end** — the 168-char spaced ciphertext → `I LIKE KILLING PEOPLE …` with the recovered
stages asserted `== [4x42,TP,P:3][56x3,UTP,P:2]`; and a length cliff (recovery vs length, clean from ~60).
Since the solver is deterministic-exhaustive there is no schedule to tune — the test is a pure correctness /
capability guard), and
`tests/test_period_column_space_solver.c` (Period column order, space-robust: the `(0,0)` superset check —
a CLEAN cipher recovers exactly, so the solver is a strict superset of the deterministic one; a **dropped-
character** repair — delete one cipher cell, the solver INSERTS a gap restoring the length and recovers all
but the lost letter; a two-dropped-character floor; and an **added-character** repair — splice a spurious
letter, the solver DELETES exactly that cell and recovers the plaintext 100%. All over a 168-cell spaced
plaintext at a fixed seed, so recovery is deterministic; trimmed restart/climb budgets keep it ~25s), and
`tests/test_double_transposition_solver.c` (Double columnar divide & conquer: end-to-end recovery via
`dct_solve_core` of several key-length pairs (6×6, 9×9, 10×10) at incomplete-rectangle lengths, plus a
length cliff — all asserted at 100% at fixed seeds with lean budgets, exercising the regime well past
`transcol2`'s ~12-15 key-length ceiling; ~45s), and
`tests/test_pollux_solver.c` (Pollux: exact recovery across assignments × RNG seeds at a solid length
(asserted 100%); a recovery-vs-length cliff printed down to the limit with an asserted floor (the cliff
sits ~24 letters — 0% at 12, ~31% at 16, ~97% at 20, 100% from 24 — far below the ACA 80-100 band); a
validity-weight sweep in the marginal regime showing the reward is essential there (weight 0 → 0% vs
weight 6 → ~98% at 16 letters); a reward-only-vs-`-logprob` characterization (reward-only wins;
`-logprob` degrades short-text recovery); and determinism. Since the solver is deterministic-exhaustive
there is no schedule to tune — the sweeps CALIBRATE the single knob, `POLLUX_VALID_WEIGHT`; ~5s).
`ciphers/tests/` additionally holds
end-to-end cases (ciphertext + `*_solution.txt`, plus `*_solve.sh` runners — e.g. the
`transcol_*_solve.sh` columnar recovery tests and `playfair_solve.sh`) you can run by hand.

`ciphers/tests/run_tests.sh` is the **accuracy regression suite**: a manifest of
74 end-to-end cases (Vigenère, Gronsfeld, Beaufort, Porta, Quagmire I–IV, autokey, the ACA
`q*_p1xx` puzzles, pure-transposition types, a homophonic substitution, a Playfair
cipher, a Bifid cipher, a Trifid cipher, a Hill cipher, the three Phillips
variants — Row / Column / Row-Column — a Two-Square (horizontal + vertical), a
Four-Square cipher, an ADFGX cipher (`adfgx_decl`, K pinned), the Nihilist
Substitution family (carry / no-carry / mod-100, plus a keyed-label cipher solved as a
relabelled square), a Gromark (`gromark_decl`, blind) + Periodic Gromark
(`gromark_periodic_decl`, blind, period swept), the three Nicodemus variants
(`nicodemus_decl` / `nicodemus_variant_decl` / `nicodemus_beaufort_decl`, P/H pinned), a
Bazeries cipher (`bazeries_decl`, digit count pinned), a Portax cipher (`portax_decl`, period
pinned), the three Progressive Key bases (`progkey_decl` / `progkey_var_decl` /
`progkey_beau_decl`, period + progression pinned), and the three Slidefair bases (`slidefair_decl` /
`slidefair_var_decl` / `slidefair_beau_decl`, period pinned), a Seriated Playfair cipher
(`seriated_playfair_decl`, period pinned, `-logprob`), a Digrafid cipher
(`digrafid_pride`, period pinned, `-logprob`), a CM Bifid cipher
(`cm_bifid_pride`, ODD period 7 pinned, `-logprob`), a Tri-Square cipher
(`trisquare_pride`, three keyed squares, `-logprob`), the Interrupted Key family
(`intkey_decl` / `intkey_var_decl` / `intkey_beau_decl` — ct-interruptor, blind interruptor;
`intkey_ptint_decl` — pt-interruptor, interruptor pinned; `intkey_breaks_decl` — supplied `-breaks`;
all period pinned), a Fractionated Morse cipher (`fracmorse_pride`, keyed-alphabet anneal,
`-logprob`), a Period column order cipher (`period_column_pp`, a two-stage 168-letter
composition solved deterministically, `-depth 2`), and a Pollux cipher (`pollux_pp`, a ~90-letter digit
cipher solved by the deterministic exhaustive 3^10 search, no budget args)) that each
solve to ~100% with a **fixed `-seed`** and quadgrams. It runs the solver, pulls the
recovered plaintext from the last field of the `>>>` CSV line, compares it
character-for-character to a sibling `<name>.solution` (bare A–Z plaintext), and prints
per-test accuracy + time + mean, exiting non-zero if any test drops below the threshold
(default 99%; the homophonic case lands ~99.9%). Because the seed is fixed, a
bit-identical refactor keeps every score at 100% and any behavioural regression shows up
immediately. Each test's `-nrestarts`/`-nhillclimbs` are trimmed to the smallest that
still lands on the solution at the seed, so the full run is ~2 min (was ~45 before
trimming). The manifest tags each case `fast` or `slow`:
`./run_tests.sh --fast` runs the 38-case fast tier in ~62s (use while iterating; incl. the three
~0s Progressive Key bases, the three ~0s Slidefair bases, the five ~1s Interrupted Key cases, the
~2s deterministic Period column order case, and the ~0s deterministic Pollux case),
`--slow` the 32 heavier ciphers (incl. the ~24s Playfair, ~6s Bifid, ~18s Trifid, the
three ~13s Phillips solves, the two ~10s Two-Square solves, the ~17s Four-Square, the ~20s Tri-Square, the
~10s ADFGX, the four ~6–8s Nihilist Substitution solves, the three ~1s Nicodemus
solves, the ~1s Bazeries solve, the ~3s Seriated Playfair solve, the ~26s Digrafid solve, and the ~3s
Fractionated Morse solve), no flag runs both.
Add a case by appending a
`tier|name|type|cipher|args` line and running `./run_tests.sh --generate <name>`
once the recovered text is verified correct.

## Run

Run from this directory — the binary loads its n-gram table, dictionary, and
ciphertext from the current working directory.

```bash
./example.sh
# or, minimally:
./colossus -type q3 -cipher cipher.txt -ngramsize 4 -ngramfile english_quadgrams.txt
```

Required flags: `-type`, a cipher source (`-cipher <file>` or `-batch <file>`),
`-ngramsize`, and `-ngramfile`. Everything else has defaults (see `init_config`).
`-type` accepts aliases or integer codes: `vig`/`0`, `q1`..`q4`/`1`..`4`, `beau`/`5`,
`porta`/`6`, `auto`/`7`, `auto1`..`auto4`/`8`..`11`, `autobeau`, `autoporta`,
`transmatrix`/`14`, `transperoffset`/`15`, `transposition`/`16`, `transcol`/`17`,
`transcol2`/`18`, `indep`/`28`, `homophonic`/`29`, `playfair`/`pf`/`30`, `bifid`/`bf`/`31`,
`trifid`/`tf`/`tri`/`32`, `hill`/`33`, `gronsfeld`/`gron`/`34`,
`phillips`/`phil`/`35`, `phillips-c`/`36`, `phillips-rc`/`37`,
`twosquare`/`ts`/`38`, `twosquare-v`/`tsv`/`39`, `foursquare`/`fs`/`40`,
`transcol-l`/`coltrack`/`41`, `transroutecol`/`routecol`/`42`, `transtile`/`tile`/`43`,
`adfgx`/`44`, `adfgvx`/`adfg`/`45`,
`nihilist-sub`/`nihsub`/`46`, `nihilist-sub-nc`/`47`, `nihilist-sub-m100`/`48`,
`gromark`/`gm`/`49`, `gromark-periodic`/`pgromark`/`50`,
`nicodemus`/`nico`/`51`, `nicodemus-variant`/`nicov`/`52`, `nicodemus-beaufort`/`nicob`/`53`,
`bazeries`/`baz`/`54`, `portax`/`ptx`/`55`,
`progkey`/`pk`/`56`, `progkey-var`/`pkv`/`57`, `progkey-beau`/`pkb`/`58`,
`slidefair`/`sf`/`59`, `slidefair-var`/`sfv`/`60`, `slidefair-beau`/`sfb`/`61`,
`seriated-playfair`/`serpf`/`spf`/`62`,
`digrafid`/`df`/`dgf`/`63`,
`cm-bifid`/`cmbifid`/`cmb`/`64`,
`trisquare`/`tri-square`/`3square`/`3sq`/`trisq`/`65`,
`interrupted-key`/`intkey`/`ik`/`66`, `interrupted-key-var`/`intkey-var`/`ikv`/`67`,
`interrupted-key-beau`/`intkey-beau`/`ikb`/`68`,
`condi`/`cond`/`69`,
`fractionated-morse`/`fracmorse`/`fmorse`/`fm`/`70`,
`period-column`/`periodcol`/`pcol`/`transpercol`/`71`,
`period-column-space`/`pcol-space`/`pcolspace`/`pcolsp`/`transpercolspace`/`72`,
`transcol2-dc`/`transcol2dc`/`dctrans`/`doublecol-dc`/`dcol`/`73`,
`pollux`/`pol`/`74`
(full list in `parse.c`; codes in `colossus.h`). Output is a human-readable block followed by a
`>>> ...` one-line CSV summary that batch runs grep/sort.

By default only the **first line** of the `-cipher` file is read (the rest is ignored,
e.g. a trailing `plaintext = ...` annotation). Pass **`-multiline`** to read the whole
file, dropping newlines so a ciphertext laid out over several lines (e.g. a homophonic
grid like Zodiac-408) is concatenated into one symbol stream.

**Tokenized symbol I/O.** Ciphertext is decoded by `decode_cipher()` (utils.c), not the
bare `ord()`. For every type except homophonic with no `-delimiter`, this is byte-for-byte
the historical per-character / 0..25-letter encoding (so the regression suite stays
bit-identical). For `homophonic` (or any type run with `-delimiter <char>`) it tokenizes
the input into a `SymbolTable` of distinct surface tokens and emits one **symbol id** per
position, so a ciphertext alphabet larger than A..Z -- comma-separated numbers
(`12,5,99,12`) or arbitrary ASCII symbols -- can be entered and displayed consistently
(`print_cipher()`). Default delimiter: auto (comma if the homophonic input contains one,
else per-character).

**`-logprob` (a.k.a. `-azdecrypt`).** Opt-in AZDecrypt / Practical-Cryptography n-gram
fitness: `load_ngrams` builds log10-probabilities with a floor that **penalises unseen
n-grams**, instead of the default reward-only normalized `log(1+count)` table (which
leaves unseen n-grams at 0). `ngram_score` keeps the score at scale 1 in this mode (a
mean log-probability). Default off => the table and every existing solve are unchanged.
Recommended with higher-order n-grams (e.g. `-ngramsize 5 -ngramfile english_quintgrams.txt`)
for hard substitution attacks; quintgrams take a homophonic solve from ~98% (quadgrams)
to ~100%.

**`-reversengrams` (a.k.a. `-revngrams`).** Opt-in **reversal-invariant** n-gram fitness:
after `load_ngrams` builds and normalizes the table (either scoring mode above), it
symmetrizes it so every n-gram and its **digit-reversed twin** carry the same (`max`)
weight, controlled by the global `g_ngram_reverse`. Because a reversed word's n-grams are
the reverses of the forward word's, a plaintext written with any words/segments **reversed**
then scores like clean English (verified: reversed English → the *same* score as forward,
vs the random floor without it; forward text unchanged). Intended for transposition attacks
where the pre-transposition plaintext may have reversed runs (e.g. the W168 alternate-word-
reversal hypothesis). Trade-off: the symmetric table can no longer tell a text from its
reverse, so it roughly doubles the acceptable solution set (mild extra gaming room at large
transposition keys). Default off => the table and every existing solve are **bit-identical**.

Five **pure transposition** cipher types bypass the keyword/cycleword/period machinery
and are solved by optimization instead (all isolated from the polyalphabetic pipeline by
an early branch in `solve_cipher`):
- `transmatrix`/`transperoffset` → `solve_transposition()` +
  `shotgun_transposition_climber()`: optimize the transform's own small parameter vector
  (`transmatrix` → `w1,w2,direction`; `transperoffset` → `period d, offset n`) with the
  shotgun/slip hill climber and n-gram scoring.
- `transposition` → `solve_general_transposition()` + `shotgun_permutation_climber()`: an
  AZDecrypt-style solver that hill-climbs the **full permutation key** (`decrypted[i] =
  cipher[key[i]]`). Restarts are seeded from columnar layouts; a periodic-redundancy
  structure term (`key_structure_score`, weight `-weightstructure`, default 4) guards
  against n-gram-gaming; simulated-annealing acceptance; and a period-targeted column-swap
  move reorders whole columns. Stochastic — run more restarts/iterations for hard ciphers.
- `transcol`/`transcol2` → `solve_columnar()` + `shotgun_columnar_climber()`: a
  **dedicated columnar** solver that, unlike the general one, optimizes only the small
  per-stage **column-order permutation** (length `K` = column count) via the
  `decrypt_columnar()` primitive (`transpositions.c`). Single (`transcol`) sweeps the
  column count over `-mincols..-maxcols` (default 2..30); double (`transcol2`) randomises
  `(K1,K2)` per restart and anneals both keys. Read direction is opt-in via
  `-readdir tb|bt|both` (default `tb`); incomplete final rows are handled automatically.
  No structure-score guard is needed — every candidate is a genuine columnar layout. Move
  set is column swaps (dominant) + short reverses/block-moves with the same Metropolis
  annealing as the permutation climber.

These `-type` values are distinct from the `-transmatrix`/`-transperoffset` *post-decrypt
stage* flags, which apply a fixed, user-supplied transposition after a polyalphabetic solve.

**Three additional transposition solvers** (added from the W168 toolkit review; see
`COLOSSUS_ADDITIONS.md`) extend the columnar family. All preserve spaces/periods as grid
cells and are isolated by their own early branch in `solve_cipher`:
- `transcol-l`/`coltrack`/`41` → `solve_columnar_track()` (`columnar_track_solver.c`): a
  **columnar with a within-column row permutation `L`** (the jarl / "dave transposition"
  scheme). The engine anneals the column order under the cheap **identity-L** reading (the
  within-row n-grams discriminate the order); the exact best `L` is recovered **once at
  report** via the Held-Karp **seam decomposition** (`seam_best_row_order`, `trans_common.c`)
  -- nesting best-L per eval is both slower and worse (it flattens the column-order contrast).
  Needs a complete grid (`len % K == 0`) for best-L; a ragged grid degrades to a plain
  columnar. Sweeps `K` over `-mincols..-maxcols`, optionally crossed with `-readdir` (column
  dir) and `-readrowdir lr|rl|both` (row dir, Rec 4). The reward-only quadgram table makes the
  seam ambiguous, so it effectively needs **`-logprob`** and/or the dictionary term.
  **`-cribanchored`** switches to a STRUCTURAL crib attack (a soft crib gives no gradient on a
  shallow many-column grid): the cipher's R-char blocks ARE the grid columns, a crib fixes some
  cells of each column, and a backtracking **block<->column matching** (most-constrained
  column first, n-gram + word-coverage tie-break) collapses the otherwise-intractable K-column
  search -- the only reliable attack on a shallow keyed columnar (assumes `L = identity`). This
  cracks all ten length-28 keyed columnars in `ciphers/W168/dave_tests/` from a 3-row crib
  (`ciphers/W168/dave_tests/solve_all.sh`, 10/10).
- `transroutecol`/`routecol`/`42` → `solve_route_chain()` (`route_chain_solver.c`): a two-stage
  **chain** -- a fixed geometric read-route global (`route_cells`, the 6 colossus routes)
  composed with a searched **column key**, read with the seam best-L. Sweeps the complete-grid
  rectangles and all routes; anneals the column key. Blind recovery is hard (the route adds a
  layer); best on short / favorable ciphers.
- `transtile`/`tile`/`43` → `solve_tile()` (`tile_solver.c`): a **sub-grid / tile
  transposition** -- every `h x w` tile of the grid is permuted by the same cell permutation
  (`-tile h w`, default 2x2), composed with a columnar column-order global. The engine
  **jointly** anneals the column order and the tile permutation; primitive `decrypt_tile`
  (`transpositions.c`). Complete grid only.

**Shared scoring additions (Rec 2, `dict.c`):** an optional `-weightword <f>` (default 0 =>
bit-identical) folds a length-weighted dictionary **word-coverage** reward into the
space-preserving transposition decrypt score (`word_set_build`/`word_coverage` + a fast hash
set); the same coverage breaks ties in the seam best-L. New raw additive n-gram sum
`ngram_sum_raw` (`scoring.c`) makes the seam decomposition exact. These are used only by the
new solvers; every existing solve is byte-for-byte unchanged. Note: `main()` no longer trims
trailing whitespace for the pure-transposition types (a trailing space is a real grid cell);
all existing transposition test ciphers are trailing-space free, so this is bit-identical.

The **homophonic** type (`solve_homophonic()`) is a `CipherModel` plugged into the shared
`run_solver()` engine just like every other type (`HOMOPHONIC_MODEL`, `SHAPE_ANNEAL`). Its
state is the many-to-one map `symbol_id -> plaintext letter` (carried in the `key` lane);
`decrypted[i] = key[cipher[i]]`. Unlike a 26->26 substitution (a bijection), a homophonic
map is free to fold many symbols onto E/T/A... to tile common n-grams -- a fixed point that
out-scores the true plaintext on raw n-grams. Two things prevent that collapse: (1) a
**monogram chi-squared penalty** on the decrypted letter distribution (`-weightmono`,
default 1.0), folded in via the decrypt hook's `score_adjust` and the greedy move; and (2) a
move set built around a **greedy coordinate step** (best plaintext letter for one symbol)
plus a **letter-class swap** (exchange the whole homophone classes of two letters, to cross
equal-frequency ambiguities like W<->M that single-symbol moves cannot). Seeds draw each
symbol's letter from the English monogram distribution. With quadgrams it recovers
~98%/~99%+ depending on homophone density; `-logprob` + quintgrams take it to ~100%.
Generate test ciphers with `tools/homophonic_gen.c` (`make homophonic_gen`).

The **Playfair** type (`solve_playfair()`, `PLAYFAIR_MODEL`, `SHAPE_ANNEAL`) is a
digraphic substitution over a **5x5 keyed grid** of 25 letters. The binary forces a
25-letter alphabet for `-type playfair` (J merged into I by ACA convention) via
`init_alphabet("J")` *before* `load_ngrams`, so the n-gram table is built over the same
25 letters (base-25 packing) and the grid is simply a permutation of `0..24` carried in
the `key` lane. The primitives live in `playfair.c` (`playfair_decrypt` is all the solver
needs; `encrypt`/`prepare`/`grid_from_keyword` serve the generator + unit tests). The
attack hill-climbs / anneals the grid with n-gram scoring (the classic SA Playfair
break): the move set is a single **cell swap** (dominant) plus **row/column swaps** and
**grid reflections** — the larger moves jump the local optima a cell swap can't escape;
cyclic row/column *rotations* are deliberately excluded (they re-encipher identically, so
the recovered grid is unique only up to such a rotation, but the recovered plaintext is
not). No anti-collapse penalty is needed (a grid is a bijection). Playfair is genuinely
near the limit of a quadgram attack: `-logprob` is effectively required, and recovery is
reliable from ~600+ characters and falls off a cliff below a few hundred (see
`tests/test_playfair_solver.c`, which prints the curve). Like every type it honours
`-method shotgun|anneal`; annealing is the default and far stronger here. Generate test
ciphers with `tools/playfair_gen.c` (`make playfair_gen`).

The **Bifid** type (`solve_bifid()`, `BIFID_MODEL`, `SHAPE_ANNEAL`) is Delastelle
fractionation over a **keyed Polybius square**: each letter splits into (row, col)
coordinates, the block's rows-then-cols coordinate stream is re-paired into ciphertext
letters, and the block size is the **period**. It defaults to the same 5x5, 25-letter,
J->I square as Playfair (`init_alphabet("J")` before `load_ngrams`), with the square a
permutation of `0..24` in the `key` lane; the primitives (`bifid.c`) are **side-generic**
(`side`/`n = side*side` parameters) so a 6x6 (36-cell) square works once a 36-letter
alphabet is active. The square attack is identical to Playfair's (cell-swap-dominated
anneal + row/column swaps and reflections, no anti-collapse penalty — a square is a
bijection). The extra dimension is the period: `bifid_estimate_periods()` ranks trial
periods by the **columnar Index of Coincidence** (the existing `mean_ioc()`), which peaks
at the true period because each within-block position becomes a coordinate-pure column.
The IoC also peaks at *multiples* of the true period, so the true period is not always
rank 1 — but it is reliably in the top-K, and the solver anneals the top `-nperiods` (default
5) candidates as separate engine configs and lets the n-gram score pick the winner (a wrong
period decrypts to gibberish). `-period N` pins a single period; `-maxperiod` bounds the
estimator's scan (default `min(20, len/2)`). Like Playfair it effectively needs `-logprob`
and recovers reliably from ~350+ characters (see `tests/test_bifid_solver.c`, which prints
the period-estimator hit rate and the length cliff). Generate test ciphers with
`tools/bifid_gen.c` (`make bifid_gen`).

The **Trifid** type (`solve_trifid()`, `TRIFID_MODEL`, `SHAPE_ANNEAL`) is Bifid lifted
into **three dimensions**: each letter splits into (layer, row, col) coordinates over a
**keyed 3x3x3 cube**, the block's layers-then-rows-then-cols coordinate stream is re-
grouped into consecutive *triples* that index new cube cells, and the block size is the
**period**. A 3x3x3 cube has **27 cells**, one more than the 26-letter alphabet, so Trifid
runs on a **27-symbol alphabet — A..Z plus a 27th symbol `+`** (`init_alphabet_trifid()`
before `load_ngrams`; the cube is a permutation of `0..26` in the `key` lane). This is the
one type whose runtime alphabet `g_alpha` exceeds 26: `ALPHABET_SIZE` (26) stays the
hardcoded mod base of the polyalphabetic primitives, and a separate `MAX_ALPHABET_SIZE`
(27) sizes only the runtime alphabet maps (`g_idx_to_char_arr`, `g_monograms`) and the
ciphertext-facing IoC scratch — everything indexable by a live symbol id. The `+` decodes
because `ord()`/`char_to_index()` consult `g_char_to_idx` for any registered ASCII char,
not just letters (bit-identical for A..Z input). The primitives (`trifid.c`) are
**side-generic** (`side`/`n = side^3`). The cube attack is the same anneal as Bifid/Playfair
(cell-swap-dominated + structured plane-swap/reflection moves along the cube's three axes,
no anti-collapse penalty — a cube is a bijection). The period is recovered exactly as
Bifid's (`trifid_estimate_periods()` over `mean_ioc()`, top-`-nperiods` annealed). Like
Bifid it effectively needs `-logprob` and recovers reliably from ~450+ characters (see
`tests/test_trifid_solver.c`). Generate test ciphers with `tools/trifid_gen.c`
(`make trifid_gen`).

The **Hill** type (`solve_hill()`, `HILL_MODEL`, `SHAPE_ANNEAL`) is a **polygraphic
substitution**: a block of `k` plaintext letters (a column vector) is multiplied by a
`k×k` key matrix **mod 26**, so it runs on the **full 26-letter alphabet unchanged**
(`ALPHABET_SIZE` is already the mod base — no `init_alphabet` forcing). The crucial design
choice: the state carried in the `key` lane **is the decryption matrix `D`**, applied
straight to the ciphertext (`plain = D·cipher mod 26`), so the hot path never inverts a
matrix; the true plaintext came from an invertible encryption key `K`, so the climb
converges on `D = K⁻¹`, and the matrix is inverted **only at report time**
(`hill_mat_inverse`, via a cofactor determinant + adjugate) to display the recovered
encryption key. The block size `k` has no IoC-style estimator, so the solver simply
**sweeps `k = 2..5`** (one engine config each, `-period` pins one; `cipher_len >= 2*k`
required) and the n-gram score picks the winner — a wrong `k` decrypts to gibberish. The
attack is a matrix anneal: ~85% change one element to a different random value (the
dominant fine move), ~10% randomize a whole row, ~5% add a random multiple of one row to
another mod 26 (a coarse jump). **A singular decryption matrix (det not coprime to 26) is
penalised** (`HILL_SINGULAR_PENALTY`, via the decrypt hook's `score_adjust`): unlike
Playfair/Bifid/Trifid a Hill matrix is a bijection *only* when invertible, and a singular
one folds the ciphertext onto a sub-lattice, decrypting to a low-entropy repetitive string
(the zero matrix → all `A`s) that out-scores real plaintext on n-grams and would otherwise
attract the climb into a collapse. Like the other near-the-limit types it effectively needs
`-logprob`. **The
search lever is restarts, not iterations**: the matrices are small (k=2 is only 26⁴ keys),
greedy climbs converge fast, and the landscape is rugged, so the schedule favours **many
short restarts**. k=2 and k=3 are reliably breakable ciphertext-only; k≥4 is exercised
only by the primitive round-trip/inverse tests (`tests/test_hill.c`), not asserted as a
ciphertext-only solve (`tests/test_hill_solver.c` characterizes the k=2/k=3 capability).
The primitives (`hill.c`) are **generic in `k`** (the determinant/inverse use cofactor
expansion, fine for `k <= HILL_MAX_K = 5`). Generate test ciphers with `tools/hill_gen.c`
(`make hill_gen`).

The **Gronsfeld** type (`gronsfeld`/`gron`/`34`) is a **Vigenère cipher with a numeric
key**: the per-column shift is a key digit `0..9` (`C = P + d`, `P = C − d`, mod 26), so
it is exactly Vigenère restricted to the 10 smallest shifts. Unlike Hill/Bifid/… it is
**not** a separate `CipherModel` — it is a type *inside* the shared `POLYALPHA_MODEL`
(`polyalpha_solver.c`), wired to behave like `VIGENERE` (straight pt/ct alphabets, only the
cycleword is searched) with one change: the **cycleword/shift domain is bounded to `0..9`**
(`GRONSFELD_DIGITS`) in the seed, perturb, and — crucially — the `derive_optimal_cycleword`
column search (`smax` in `optimal_cycleword.c`). So it reuses the whole polyalpha pipeline
(IoC period estimation, the deterministic optimal-cycleword frequency attack, crib
handling); the digit bound is a strong prior that makes recovery faster and reliable from
shorter text than an unconstrained Vigenère solve. The cycleword *is* the key, so the
decrypt path calls a tight direct primitive (`gronsfeld_decrypt`, `gronsfeld.c`) rather
than the Quagmire indirection, and the report prints the recovered key as digits. It rides
the polyalphabetic search defaults (no registry entry) and the default reward-only quadgram
table (no `-logprob` needed, like Vigenère). Generate test ciphers with
`tools/gronsfeld_gen.c` (`make gronsfeld_gen`; the key arg is a digit string, e.g. `31415926`).

The **Phillips** type (`solve_phillips()`, `PHILLIPS_MODEL`, `SHAPE_ANNEAL`) is a **periodic
monographic substitution over 8 keyed Polybius squares** derived from one base 5x5 square. The
plaintext is split into blocks of 5 letters; block `b` is enciphered with square `(b mod 8)`,
each letter going to the one diagonally **down-right with wrap** (`cipher = sq[(r+1)%5][(c+1)%5]`),
so the overall period is 40. The 8 squares are derived from the base by a fixed row-reinsertion
table (squares 1–5: base row 0 reinserted at row positions 0–4; squares 6–8: then row 1 reinserted
at positions 1–3) — verified cell-for-cell against the ACA worked example. It runs on the same
25-letter (J→I) grid as Playfair/Bifid (`init_alphabet("J")` before `load_ngrams`), with the base
square a permutation of `0..24` in the `key` lane. The only unknown is that base square, so the
attack is **identical to Playfair's** (cell-swap-dominated anneal + row/column swaps and grid
reflections, no anti-collapse penalty — every derived square, and hence the whole map, is a
bijection) — but *monographic* (no digraph prep) and with **no period to estimate** (block size 5
and the 8-square cycle are fixed → a single engine config). Because it carries more signal per
character than digraphic Playfair, it recovers reliably from **~200+ characters** (see
`tests/test_phillips_solver.c`) and rides a leaner registry budget than Playfair (`4×250000`).
Like the other near-the-limit square types it effectively needs `-logprob`. The cyclic-column
rotation of the base re-enciphers identically, so the recovered square is unique only up to that
(the plaintext is unique). **Three variants** select how the 8 squares are built, all sharing the
solver/primitive: `phillips` (the ACA-standard **Row** type), `phillips-c` (**Column** — the same
reinsertion applied to columns), and `phillips-rc` (**Row-Column** — rows for squares 2–5, columns
for 6–8). The primitives (`phillips.c`) are **side-generic** (`2·side−2` squares of `side` letters).
Generate test ciphers with `tools/phillips_gen.c` (`make phillips_gen`; the variant arg is
`row`/`col`/`rowcol`).

The **Two-Square** (`solve_twosquare()`, `TWOSQUARE_MODEL`) and **Four-Square**
(`solve_foursquare()`, `FOURSQUARE_MODEL`) types are digraphic substitutions over **two keyed
5x5 squares** (`SHAPE_ANNEAL`). They run on the same 25-letter (J→I) grid as Playfair
(`init_alphabet("J")` before `load_ngrams`), but the state is a **pair** of squares packed
back-to-back in the `key` lane (sq1 = `key[0..24]`, sq2 = `key[25..49]`) — **50 cells, double
Playfair's**. For each plaintext digraph the two letters span a rectangle and the cipher pair is
the rectangle's other two corners (`twosquare.c` / `foursquare.c`; verified cell-for-cell against
the ACA Two-Square and Wikipedia Four-Square worked examples). The attack is the **same SA square
break as Playfair** — cell-swap-dominated moves plus row/column swaps and reflections, no
anti-collapse penalty (every square is a bijection) — but each move perturbs **one of the two
squares** (chosen uniformly) and there is **no period to estimate** (a single engine config). The
only prep is padding an odd plaintext with one `X` (no doubled-letter handling, unlike Playfair).
The larger 50-cell state needs more text and a bigger budget; like the other square types it
effectively needs `-logprob`. **Two-Square has two arrangements** sharing the solver/primitive,
selected by `cfg->cipher_type`: `twosquare` (the ACA **horizontal**, squares side by side — a
same-ROW digraph maps to the reversed pair) and `twosquare-v` (**vertical**, squares stacked — a
same-COLUMN digraph maps to itself, and the whole cipher is **self-inverse** so decrypt ==
encrypt). About **20% of digraphs are such transparencies**, a documented weakness that leaks
plaintext and makes Two-Square recover from short text. **Four-Square** keeps the upper-left and
lower-right squares as the **fixed standard alphabet** (`foursquare_standard_square()`); only the
upper-right (UR) and lower-left (LL) squares are keyed, and the two are independent — making it
the hardest of the family (its report breaks recovery down per-square: even cipher positions
decrypt through UR, odd through LL). The primitives are **side-generic** (`side`/`n = side*side`).
Generate test ciphers with `tools/twosquare_gen.c` / `tools/foursquare_gen.c`
(`make twosquare_gen foursquare_gen`; the Two-Square variant arg is `h`/`v`).

The **ADFGVX** (`adfgvx`/`adfg`/`45`) and **ADFGX** (`adfgx`/`44`) types (`solve_adfgvx()`,
`ADFGVX_MODEL`, `SHAPE_ANNEAL`) are the classic WWI cipher: a **keyed Polybius-square
fractionation composed with a keyed columnar transposition** — the first
*fractionation-then-transposition* type, and the hardest of the polygraphic family because
the square and the column order must be recovered **jointly**. Each plaintext symbol becomes
its two cell coordinates, each a label from `{A,D,F,G,X}` (ADFGX, 5x5, 25 letters J→I) or
`{A,D,F,G,V,X}` (ADFGVX, 6x6, **36 symbols** A..Z + 0..9), so the ciphertext is **2N labels**
for an N-symbol plaintext; that 2N coordinate stream is then columnar-transposed under a
keyword of length K. ADFGX runs on the same 25-letter J→I alphabet as Playfair/Bifid;
ADFGVX forces a **36-symbol alphabet** (`init_alphabet_adfgvx()`, A..Z + digits, with the
digits at negligible monogram weight like the Trifid `+`), which is why `MAX_ALPHABET_SIZE`
is **36** (the largest runtime alphabet, sizing only the runtime alphabet maps + the IoC
scratch; the n-gram table is then base-36 — quadgrams 36⁴≈6.7M floats are fine, but **avoid
quintgrams** at base 36, ~240MB). The primitive (`adfgvx.c`) works in **coordinate space**
(the solver maps the label characters to coordinates 0..side-1 up front) and is **side-generic**,
reusing `bifid_grid_from_keyword`/`bifid_build_inverse` for the square and the exposed
`decrypt_columnar` for the transposition stage — so it adds almost no new cipher math.
The state carries the square (`key[0..n-1]`, a permutation) **and** the column order
(`key[n..n+K-1]`); one engine config is enumerated per column count K in `-mincols..-maxcols`,
the square move set is Bifid's, the column-order move set the columnar solver's, and a move
perturbs one or the other. **The key to making the coupled search tractable** is a
**structural IoC reward folded into `score_adjust`**: after undoing the *correct* columnar
the paired cell ids are a monoalphabetic image of the plaintext, so the decrypt's Index of
Coincidence is English (~0.066) — and that IoC depends ONLY on the column order, not the
square (a square just relabels cells). The reward therefore gives the column-order search a
gradient independent of the square, decoupling the two halves (the climb locks the column
order by IoC, then the n-gram score recovers the square). Like the other near-the-limit
square types it effectively needs `-logprob`; ADFGX recovers reliably from ~200+ characters,
ADFGVX (the 36-cell square) needs more text and a bigger budget. Cribs are not used (the crib
positions are over the ciphertext, which is 2x the plaintext). Generate test ciphers with
`tools/adfgvx_gen.c` (`make adfgvx_gen`; args are a square keyword, a transposition keyword,
and `adfgx`/`adfgvx`).

The **Nihilist Substitution** family (`nihilist-sub`/`46`, `nihilist-sub-nc`/`47`,
`nihilist-sub-m100`/`48`; `solve_nihilist_sub()`, `NIHILIST_SUB_MODEL`, `SHAPE_ANNEAL`) is a
**periodic ADDITIVE** cipher over a keyed 5x5 Polybius square (25 letters, J→I, same alphabet as
Bifid). Each plaintext letter → its 2-digit coordinate **number** `rowlbl[row]*10+collbl[col]`
(fixed labels are `1..5`, so the legal set is `V = {11..15,21..25,…,51..55}`); a periodic additive
key (its own coordinate numbers) is added per position, so the ciphertext is a stream of decimal
**NUMBERS** (not letters) — `solve_nihilist_sub` parses them from the raw string itself (numbers
are space/comma/any-separated; not the per-character decode). **Three addition conventions** are
distinct `-type` codes sharing one primitive (`nihilist_sub.c`) + one solver, branched on
`cfg->cipher_type`: `nihilist-sub` (integer add **with carry**, ACA standard, cipher 22–110),
`nihilist-sub-nc` (per-digit add **mod 10, no carry**, cipher 00–99) and `nihilist-sub-m100` (2-digit
add **mod 100**). It is structurally the twin of ADFGVX — a COUPLED square + periodic-key search —
and uses the **same decoupling trick**: `pt_num = cipher_num − key_num` is INDEPENDENT of the
square, so the fraction of positions decrypting to a legal coordinate (the **validity** reward,
`NIH_VALID_WEIGHT * n_valid/n` folded into `score_adjust`) gives the additive-key search a gradient
flat in the square dimension — the climb locks the additive by validity, then the n-gram score
recovers the square (the additive cells live in `key[grid_size..]`, the square in `key[0..grid_size-1]`,
exactly like ADFGVX). The period is recovered by columnar IoC over the ciphertext **numbers**
(`nihilist_sub_estimate_periods`, top-`-nperiods` annealed; IoC also peaks at multiples, so a
multiple of the true period — with the keyword repeated — is an equally correct solve). Cribs are not
used (positions are over numbers, not plaintext). Like the other square types it effectively needs
`-logprob`; it recovers reliably from ~250+ characters. The primitive is **side-generic** and
**label-aware**: the standard fixed `1..5` labels can be replaced by a **keyed-label** permutation
(`-labels` in the generator) — but since a label permutation only permutes which cell gets which
number (leaving `V` unchanged), the labels are **not separately identifiable ciphertext-only**: they
fold into the recovered square, so the (fixed-label) solver cracks a label-keyed cipher as the
equivalent **relabelled square** (proven by the `nihilist_sub_kl` regression case). Generate test
ciphers with `tools/nihilist_sub_gen.c` (`make nihilist_sub_gen`; args are a square keyword, an
additive keyword, `carry`/`nc`/`m100`, and optional `-labels <rowkey> <colkey>`).

The **Gromark** (`gromark`/`gm`/`49`) and **Periodic Gromark** (`gromark-periodic`/`pgromark`/`50`)
types (`solve_gromark()`, `gromark.c` + `gromark_solver.c`) are the ACA "GROnsfeld with Mixed
Alphabet and Running Key" cipher (DUMBO 1969) and its periodic variant (1973). Both run on the
**full 26-letter alphabet** (no J-merge). A keyed cipher alphabet **σ** (a permutation of A..Z,
built by the **K2M transposition block** of a keyword — the simple keyed alphabet laid row-major at
keyword-width and read off by columns in heading-letter order) is composed with a **chain-addition
running key**: from a P-digit primer, `d[i]=primer[i]` for `i<P` and `d[i]=(d[i-P]+d[i-P+1]) mod 10`
otherwise (one digit per letter). **Basic Gromark** (5-digit primer): `C[i]=σ[(p[i]+d[i]) mod 26]`.
**Periodic Gromark**: the plaintext is split into consecutive groups of P (= the keyword length);
group `g` (cycling `mod P`) adds a per-group **offset** = the keyword letter's position in σ, so
`C[i]=σ[(p[i]+d[i]+offset[(i/P) mod P]) mod 26]`, and the primer is the keyword letters' alphabetical
ranks. The primitives (`gromark.c`) are hand-verified cell-for-cell against the two ACA worked
examples (keyword ENIGMA). **The two are attacked by DIFFERENT models** (both `SHAPE_ANNEAL`,
effectively need `-logprob`), because there is **no σ-independent decoupling reward** (the additive
shift sits inside the permutation σ, and chain addition amplifies a wrong primer digit into a
fully-wrong key, so naive joint annealing of (σ, primer) never accepts a primer move):
- **Basic Gromark** uses a **primer PRE-PASS** (the analog of `bifid_estimate_periods`): for each of
  the 10⁵ primers the running key is known, so recovering σ collapses to a monoalphabetic-
  substitution-with-known-per-position-shift — a 26×26 **max-weight assignment** (Hungarian) that
  picks σ to maximize an English monogram fit — and the primer is ranked by the n-gram score of the
  provisional decrypt. The top-K primers (length-adaptive, `-nprimers` overrides) become engine
  configs, each annealing σ from a random restart warm-started by the provisional σ; the n-gram
  score across configs picks the winner. Recovers reliably from ~120+ letters and reports the
  recovered primer.
- **Periodic Gromark** instead anneals the **KEYWORD directly** — its *entire* key is one keyword
  of P distinct letters (~28 bits), which derives σ, the primer, AND the offsets together, so
  treating those as independent unknowns blows a tiny key into an intractable coupled space the
  primer cannot be ranked in. The state is the keyword (`key[0..P-1]`); each decrypt rebuilds
  σ/primer/offsets via `gromark_build_from_keyword_idx`; one engine config per swept period P (like
  a fractionation period), no pre-pass. This is both far faster and far more reliable than a free-σ
  search (recovers ~150+ letters in seconds, reporting the recovered keyword). Cribs are not used.
  Generate test ciphers with `tools/gromark_gen.c` (`make gromark_gen`; args are a plaintext, a
  keyword, and a `<primer-digits>` for basic or the literal `periodic` for periodic).

The **Nicodemus** family (`nicodemus`/`51`, `nicodemus-variant`/`52`, `nicodemus-beaufort`/`53`;
`solve_nicodemus()`, `NICODEMUS_MODEL`, `SHAPE_ANNEAL`) is the ACA **substitution + transposition**
composite — the first such type besides ADFGVX, and the most plausibly K4-relevant (a masking
construction). A single **keyword** of length P drives two stages over blocks of **H rows × P
columns** (row-major fill, final block ragged; H is the ACA-standard 5 but swept here): (1) each
grid **column** is enciphered by its keyword letter — three ACA substitution conventions are
distinct `-type` codes sharing one primitive (`nicodemus.c`) + solver, `nicodemus` (Vigenère,
`C=P+k`), `nicodemus-variant` (`C=P−k`) and `nicodemus-beaufort` (`C=k−P`, reciprocal); (2) the
block's columns are read off top-to-bottom in the keyword's **alphabetical rank order** (a per-block
columnar transposition, the same incomplete-grid rule as `decrypt_columnar` but applied per H·P
block — its `K>len` guard would misfire on a final block narrower than P). Full 26-letter alphabet
(no J-merge); cribs are not used (the per-block transposition scrambles plaintext positions).
**The key to making the composite tractable** is the same decoupling the default `-optimalcycle`
path uses for the polyalphabetic ciphers: the annealed state is the **COLUMN ORDER alone** (a
permutation of `0..P-1` in `key[0..P-1]`, perturbed with the columnar `perm_move`/`perm_seed`), and
the P per-column **shifts are DERIVED deterministically** for each candidate order — after
de-transposing, grid-column `g = i % P` is a Caesar sample, and the shift maximising its monogram
fit against `g_monograms` is chosen (specialising `derive_optimal_cycleword`). Because every column
is fit to English *monograms* regardless of order, the **n-gram (quadgram) score drives the order**
search — cross-column digraphs only form at the true order — so it effectively needs `-logprob`.
Solving this general (order, shifts) form also cracks the ACA cipher (true key = the special case
`order == argsort(shifts)`); the recovered plaintext is what's checked. Note Vigenère and Variant
are **not separately identifiable** (a free derived shift absorbs the sign), so either solver cracks
a shift-substitution Nicodemus; only Beaufort (a reflection) is distinct. The period P **and** block
height H are **swept** — one engine config per `(P, H)` pair (`-period`/`-blockheight` pin them,
`-mincols`/`-maxcols`/`-maxblockheight` bound the sweep) — IoC period estimation is useless through
a transposition. It recovers reliably from ~120+ letters (P/H known) — see
`tests/test_nicodemus_solver.c`, which prints the length cliff and tunes the schedule. Generate test
ciphers with `tools/nicodemus_gen.c` (`make nicodemus_gen`; args are a plaintext, a keyword, a
`vig`/`variant`/`beau` substitution, and an optional block height).

The **Bazeries** type (`bazeries`/`baz`/`54`; `solve_bazeries()`, `BAZERIES_MODEL`, `SHAPE_ANNEAL`)
is the ACA **"simple substitution plus transposition"** cipher, keyed by ONE number **N < 1,000,000**
that drives BOTH stages over the 25-letter J→I alphabet (so `init_alphabet("J")` is forced before
`load_ngrams`, same as Playfair/Bifid). (1) **Transposition:** the plaintext is split into groups
whose sizes cycle through N's **decimal digits** (e.g. `3752` → 3,7,5,2,3,7,5,2,…) and each group is
**reversed** (an involution; a 0 digit is a zero-length group, skipped). (2) **Substitution:** a fixed
monoalphabetic map between two 5×5 squares — the PLAINTEXT square is the alphabet entered **column-major**
(fixed), the CIPHERTEXT square is **N spelled out** ("three thousand seven hundred fifty two") used as
the keyword of a keyed square entered **row-major** (so the square build reuses `bifid_grid_from_keyword`
/ `bifid_build_inverse`; the substitution is `pt_square[(r,c)] → ct_square[(r,c)]`). Encryption =
transpose then substitute; decryption = inverse-substitute then un-transpose. The primitive
(`bazeries.c`) is hand-verified cell-for-cell against the ACA worked example (`N=3752`, the square
`THREOUSANDVFIYW…`, ciphertext `ACYYU…`). **Because Colossus is optimisation-only and the whole key is
a number, the solver CLIMBS N's decimal digits rather than adding an exhaustive driver:** the state is
the D digit values (`key[0..D-1]`, leading digit 1..9), with **one engine config per digit count D in
1..6** (`-period` pins a single D; the union over D covers all of 1..999999). The rugged < 10⁶ digit
landscape is made navigable by a **square-quality monogram reward folded into `score_adjust`** (the
analog of ADFGVX's structural IoC term): the inverse substitution is monoalphabetic and a transposition
leaves the monogram multiset unchanged, so the decrypt's **mean English-monogram fit depends only on the
square (which N is spelled), not on the transposition digits** — rewarding it pulls the digit climb toward
the correct square, after which the n-gram score discriminates the exact number. Like the other square
types it effectively needs `-logprob`; **RESTARTS are the robustness lever** (each reseeds a fresh random
number), and it recovers reliably blind from ~150+ letters. Because the model also implements
`seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it (calibrated in
`tests/test_bazeries_solver.c`). Cribs are not used (the transposition scrambles plaintext positions).
Generate test ciphers with `tools/bazeries_gen.c` (`make bazeries_gen`; args are a plaintext and the
number N).

The **Portax** type (`portax`/`ptx`/`55`; `solve_portax()`, `PORTAX_MODEL`, `SHAPE_ANNEAL`) is the
ACA **periodic DIGRAPHIC Porta** — the first digraphic polyalphabetic type. The plaintext is written
row-major into a block of width **P** (= keyword length) and rows are taken in **vertical PAIRS**
(rows 2g, 2g+1); the pair in column c — (top, bottom) — is enciphered **as a unit** over a Porta
**slide** keyed by `keyword[c]`, of which **only the Porta shift `key/2` (0..12) matters** (key
letters U and V are identical). The slide has a fixed upper half (A–M), a sliding upper half (N–Z),
and a two-row lower alphabet (even/odd letters); the two plaintext letters are diagonally opposite
corners of a rectangle whose **other two corners are the substitutes** (with a same-vertical-line
special case taking the other two cells of that line). The map is **self-reciprocal** (decrypt ==
encrypt). Full 26-letter alphabet (**no J→I merge**); the primitive (`portax.c`) is hand-verified
cell-for-cell against the ACA worked examples (key U/V `IN→JL`/`NO→UA`/`NA→DB`; key E `TA→NM`/`BG→QH`;
keyword `EASY`, `the early bird gets the worm` → `NIJAMPBGQCWKHQJEUIKYMPAT`). **The whole key is P
Porta shifts (a 13^P space)**, carried in the cycleword lane, with **one engine config per swept
period P** (IoC period estimation is useless through the digraphic pairing, so P is swept and the
n-gram score picks it — the rigid pairing makes a wrong P decrypt to gibberish; `-period` pins one,
`-mincols`/`-maxcols` bound the sweep). **The key to efficiency is per-column independence:** a
vertical pair is enciphered ENTIRELY by its column key, so every pair in column c decrypts from
`shift[c]` alone — the **per-column monogram-fit shift** (the analog of `derive_optimal_cycleword`)
**warm-starts the seed**, and the n-gram (quadgram) score then drives the anneal (cross-column
digraphs only form at the true shifts) and corrects any column the monogram fit mis-set. No
`score_adjust` is needed (every cycleword is a valid bijective decrypt). Unlike the other square /
fractionation types it **rides the reward-only quadgram table (no `-logprob`)**, like the rest of the
Porta family, and recovers reliably from **~70+ letters** (see `tests/test_portax_solver.c`). **Cribs
are supported** (the cipher is positional: `decrypted[i]` is plaintext position i). Because the model
implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it. Generate test ciphers
with `tools/portax_gen.c` (`make portax_gen`; args are a plaintext and a keyword).

The **Progressive Key** family (`progkey`/`pk`/`56`, `progkey-var`/`pkv`/`57`, `progkey-beau`/`pkb`/`58`;
`solve_progkey()`, `PROGKEY_MODEL`, `SHAPE_ANNEAL`) is the ACA **periodic key with a constant per-group
drift** — the first **non-stationary periodic key** type (a leading K4 hypothesis). The plaintext is
set in groups of **P** (= keyword length); a **primary** periodic encipherment of base type T under the
keyword yields `C1`, then a **second** encipherment **of the same type T** adds the **progressive key
letter** `Kp[g] = (g·prog) mod 26` per group `g = i/P` (the **progression index** `prog`: 1→A,B,C…,
2→A,C,E…). So the per-position shift is `keyword[i%P]` drifting by `(g·prog)`. **Three base types** are
distinct `-type` codes sharing one primitive (`progkey.c`) + solver, branched on `cfg->cipher_type`:
`progkey` (**Vigenère**, `C=P+k`), `progkey-var` (**Variant**, `C=P−k`) and `progkey-beau` (**Beaufort**,
`C=k−P`). Full 26-letter alphabet (**no J→I merge**); the primitive is hand-verified cell-for-cell
against the ACA worked example (Vigenère, key `GRAPEFRUIT`, P=10, prog=1, `thiscipher… → ZYIHG NGBMK
JSORJ AKZMQ QMJRT FHBDC…`). (`prog=0` is no drift — a plain periodic Vigenère/Variant; the Beaufort
group-0 pass is still a reflection, not identity.) **The whole key is P base shifts (0..25) plus the
progression**, carried in the cycleword lane. **IoC period estimation FAILS** (within a column each
group carries a different drifted shift, so columns are not monoalphabetic — like autokey), so the
**period is brute-forced** and the **progression enumerated 0..25**: one engine config per `(P, prog)`
pair (`cc->aux[0] = prog`; `-period`/`-cyclewordlen` pin P, `-progression` pins prog, `-mincols`/
`-maxcols` bound the period sweep). **The key to efficiency** is the same decoupling the polyalphabetic
`-optimalcycle` path uses: for a fixed `prog`, **DE-PROGRESSING** the ciphertext (undoing only the
drift pass via `progkey_deprogress`, leaving the primary base cipher `C1`) makes every column an
independent base sample under its own keyword shift — so the **per-column monogram-fit shift**
**warm-starts the seed**, and the n-gram (quadgram) score then drives the anneal AND, across all
`(P, prog)` configs, selects the true period and progression (a wrong P or prog leaves columns drifted
→ gibberish → low score). No `score_adjust` is needed (every cycleword is a valid bijective decrypt).
Unlike the square/fractionation types it **rides the reward-only quadgram table (no `-logprob`)**, like
the rest of the Vigenère family, and **recovers from very short text** (the de-progressed columns are a
pure Vigenère; see `tests/test_progkey_solver.c`). It is a **dedicated `CipherModel`, NOT inside
`POLYALPHA_MODEL`** (the engine's `engine_build_hist` keys per-column histograms on the raw cipher +
period, but progkey needs them on the de-progressed buffer, which depends on `prog` — so the solver
does its own per-column derivation, keeping the cipher-agnostic engine and shared `optimal_cycleword.c`
untouched). **Cribs are supported** (the cipher is positional: `decrypted[i]` is plaintext position i).
Because the model implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it.
Generate test ciphers with `tools/progkey_gen.c` (`make progkey_gen`; args are a plaintext, a keyword,
a progression index, and `vig`/`var`/`beau`).

The **Interrupted Key** family (`interrupted-key`/`intkey`/`ik`/`66`, `intkey-var`/`67`,
`intkey-beau`/`68`; `solve_intkey()`, `INTKEY_MODEL`, `SHAPE_ANNEAL`) is the ACA **periodic keyword
that is INTERRUPTED** — a leading K4 hypothesis and the second **non-stationary periodic key** type
(after Progressive Key). A base cipher (Vig/Var/Beau) enciphers under a P-letter keyword whose key
index **resets to the first key letter** at *break* points and otherwise increments `(k+1) mod P`
(the ACA: "enciphered with 1, 2, 3 or more letters of the keyword which is interrupted ... Return to
the first key letter each time ... The entire keyword must be used at least once"). The base shift
math is Progressive Key's, so the primitive (`intkey.c`) **reuses** `progkey_base_encrypt/decrypt`;
it is hand-verified cell-for-cell against the ACA worked example (Vigenère, keyword `ORANGE`, the
word-division break lengths 4,6,2,3,4,3,1,1,5,1,2,3,5, `THISCIPHER…PERIODICS` →
`HYIFQZPUKV…ICUIPY`). **The unifying model:** every scheme reduces to *assigning each position a key
index that resets to 0 at break points*; schemes differ only in how breaks arise (the interruption
**STRATEGY**, `IK_STRAT_*`). IoC period estimation fails through the interruption (columns are not
monoalphabetic, like autokey/progkey), so the **period is SWEPT** and the strategy is **enumerated**,
the n-gram score picking the winner:
- **ct-interruptor** (`ct`): reset AFTER a chosen **ciphertext** letter → the break positions are a
  function of the ciphertext ALONE (keyword-independent), so the per-position key index is KNOWN;
  grouping positions by index makes each column an independent Caesar/Beaufort sample recovered by a
  **monogram fit** (the `intkey_build_mask_ct` mask + the analog of `derive_optimal_cycleword`). This
  is the **reliable blind workhorse** — recovers like a plain Vigenère, ~100% from ~40 letters.
- **pt-interruptor** (`pt`): reset AFTER a chosen **plaintext** letter → breaks are causal (keyword-
  dependent), warm-started by an **EM** loop (fit keyword under the plain-periodic assumption, causal-
  decrypt to get the actual columns, refit). **FRAGILE** — because the reset trigger is a plaintext
  letter, a wrong (keyword, interruptor) pair can still de-cohere into fairly English text, so blind
  `pt` can mispick the interruptor on short text and even a *pinned* interruptor can hit a rugged
  multi-modal basin (seed/content sensitive). Characterized, not guaranteed (see the solver test).
- **supplied breaks** (`breaks`, `-breaks <file>`): user-supplied group-start positions (random /
  word-division). Mask known → decouples exactly like `ct`; **reliable**. `P` = max group length.
- **joint** (`joint`): blind random — the keyword AND an N-bit **break-mask** (in the `key` lane) are
  annealed together; the only genuinely blind attack on random breaks but reliable only on longer text
  / with cribs (characterized — a crib lifts it markedly, e.g. ~29% → ~64%).
**Efficiency**: a **Gromark-style pre-pass** scores every `(period, strategy, interruptor)` by the
n-gram of its monogram-derived decrypt and keeps only the **top-K** (`-nprimers`, default 12) to
anneal — `ct` typically solves on the warm seed. `score_adjust` stays 0 (every keyword is a bijective
decrypt; n-grams discriminate). Default blind enumeration is **ct + pt**; `-breaks` selects `breaks`,
`-intscheme ct|pt|breaks|joint` pins a strategy, `-interruptor <A-Z>` pins the letter, `-period`/
`-mincols`/`-maxcols` bound the period sweep. Like the rest of the Vigenère family it **rides the
reward-only quadgram table (no `-logprob`)**. **Cribs are supported** (positional: `decrypted[i]` is
plaintext position i). It is a **dedicated `CipherModel`, NOT inside `POLYALPHA_MODEL`** (the break
structure is not a stationary period the engine's histograms model). Because the model implements
`seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it. Generate test ciphers with
`tools/intkey_gen.c` (`make intkey_gen`; args are a plaintext, a keyword, `vig`/`var`/`beau`,
`ct`/`pt`/`random`, and the interruptor letter (ct/pt) or an RNG seed (random)).

The **Condi** type (`condi`/`cond`/`69`; `solve_condi()`, `CONDI_MODEL`, `SHAPE_ANNEAL`) is the ACA
**"Condi"** — a **plaintext-feedback substitution over a keyed alphabet σ** (a 26-permutation; full
26-letter alphabet, no J→I). The shift applied to each plaintext letter is the **1-indexed position of
the PREVIOUS plaintext letter in σ**; the first letter uses a **starter** offset (0..25). So
`idx(ct_i) = (idx(pt_i) + idx(pt_{i-1}) + 1) mod 26` for i≥2, and decryption is **causal** (each offset
needs the already-recovered previous plaintext letter). Only alphabetic letters participate (the ACA
keeps word divisions / punctuation, which Colossus strips, so the feedback is over consecutive letters).
The primitive (`condi.c`, hand-verified cell-for-cell against the ACA worked example — keyword `STRANGE`,
shifted σ `VWXYZSTRANGEBCDFHIJKLMOPQU`, starter 25, `OURS… → MORC…`) takes σ and its inverse (caller
precomputes the inverse; the hook caches it per decrypt). It is a **dedicated `CipherModel`, NOT inside
`POLYALPHA_MODEL`** — there is **no period** (a feedback stream cipher) and **no per-column /
monogram decoupling** (the running key is plaintext-derived, entangled with σ, unlike Gromark's
σ-independent primer). The shipped attack is a **free-permutation σ anneal** (random shuffle seed,
cell-swap moves, `score_adjust`=0 — every σ is a bijection) with the **26 starter values enumerated as
engine configs** (a wrong starter cascades into gibberish → low n-gram score; `-startkey` pins one). σ
is a free permutation (not a keyword+tail) because the ACA alphabet may be cyclically **shifted**, and a
rotation genuinely changes the plaintext (idx(ct) shifts by r, the feedback sum by 2r — no cancellation),
so the free permutation absorbs the rotation. **THE KEY PROPERTY (a documented structural limitation, the
analog of CM Bifid's even-period degeneracy):** the plaintext feedback makes the true σ an **ISOLATED
NEEDLE**. Because `idx(pt_i)` feeds the next offset, one σ cell swap re-derives the **entire downstream
plaintext** (an alternating-sign shift cascades from the first changed position), so there is **no basin
and no gradient** — measured on a real 360-letter solve the true σ scores ~**3.40** while its **best
single-swap neighbour scores only ~2.48** (a ~0.9 cliff) and the mean neighbour ~**1.37** ≈ the random
floor ~**1.01**. Hence **NO local search (anneal / shotgun / PSO) recovers Condi blind at any budget,
even with the starter pinned or a crib supplied** (empirically ~4–6%); the solver is a bounded honest
attempt, and `tests/test_condi_solver.c` **asserts the needle** (true ≫ best-neighbour, neighbours near
random) as a regression guard while characterizing the blind/pinned/per-scheme failure. It effectively
needs `-logprob`. **Cribs are passed through** (positional: `decrypted[i]` is plaintext position i) — the
substrate for the **tractable attack that is NOT yet built**: crib-anchored constraint solving (a crib
gives linear congruences `pos[ct_i] = pos[pt_i] + pos[pt_{i-1}] + 1 (mod 26)` that pin σ directly, cf.
the `transcol-l` `-cribanchored` structural attack; also the K4-relevant form, since K4 has cribs). Its
registry schedule is a **modest `6x60000`** (a large budget is pointless against the needle). There is
**no `run_tests.sh` end-to-end recovery case** (blind Condi cannot reach the ~99% threshold). Because the
model implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it. Generate test
ciphers with `tools/condi_gen.c` (`make condi_gen`; args are a plaintext, a keyword, a starter, and an
optional cyclic shift of the keyed alphabet).

The **Fractionated Morse** type (`fractionated-morse`/`fracmorse`/`fmorse`/`fm`/`70`;
`solve_fracmorse()`, `FRACMORSE_MODEL`, `SHAPE_ANNEAL`) is the ACA **"Fractionated Morse"** — the first
**Morse-fractionation** type. The plaintext is written in Morse with a **single `x` separator between
letters** (Colossus works on bare A–Z, so the ACA `xx`-between-words convention does not apply — external
worked examples that use word divisions will NOT match byte-for-byte), the `{DOT,DASH,X}` stream is padded
with trailing `x` to a multiple of 3 (the body has only single-`x` separators, so the run `xxx` can never
occur) and grouped into **trigraphs**; each trigraph `(a,b,c)` has **rank `r = 9a+3b+c`** (base 3,
`DOT<DASH<X`), and since `xxx` = rank 26 never occurs every group has `r ∈ 0..25`, mapped to a ciphertext
letter through a **keyed 26-letter alphabet σ** (rank → letter). This is a **length CHANGE** (N plaintext
→ C ciphertext letters) and there is **NO period** (the grouping is fixed at 3; full 26-letter alphabet,
no J→I). The primitive (`fracmorse.c`, hand-verified against KATs computed under the single-`x`
convention — keyword `MORSE`, `SOS → MWLR`, `EE → B`, `ET → C`) is table-free for the rank math and keeps
the Morse stream in a file-static scratch. **The key design points:** (1) the alphabet is an ACA **KEYED
alphabet**, searched **as such** — the state `st->key` **is σ** (a keyed-alphabet sequence, keyword prefix
`st->aux[0]` + ascending tail) perturbed by `fracmorse_move_seq`, the length-26 twin of `digrafid_move_seq`
(4% grow/shrink kw, 48% keyword↔tail + re-sort tail, 48% in-keyword reorder), seeded by `random_keyword`
with `kw` sampled per restart from `[FM_KW_MIN..FM_KW_MAX]=[3..13]` — **not** a free 26! permutation, which
collapses the keyspace and (like Digrafid) tracks the keyword so it cracks the **short ~110–150-letter ACA
range** reliably (the free permutation is still reachable at kw=26). (2) **Length handling:** the decoded
plaintext length varies per key but the engine scores a fixed length, so the decrypt hook decodes to
plaintext (a **filler** letter for each **invalid** Morse token — a run that is not a legal codeword) and
**cyclically TILES** it to the ciphertext length `C` so the mean n-gram score is length-fair across keys;
a **MORSE-VALIDITY reward** (`FM_VALID_WEIGHT * n_valid/n_tokens`) is folded into `score_adjust` (the
analog of ADFGVX's IoC term), giving the anneal a structural gradient toward clean-Morse keys. The
report/`>>>` CSV/`result` carry a **clean re-decode** (the real N letters), not the tiled buffer, so
`run_tests.sh`'s char-compare against `.solution` works. It effectively needs **`-logprob`**; **cribs are
NOT used** (the length change + tiling break the positional mapping). Its registry schedule is
`16x120000` at a **warm `inittemp 0.30`** (the keyed-alphabet moves are coarse, RESTARTS are the lever;
single config, no sweep), tuned against `tests/test_fracmorse_solver.c` (which recovers ~100% at the ACA
~150-char range, ~97% at 90). Because the model implements `seed`/`perturb`/`copy`, `-method
anneal|shotgun|pso` all run on it (PSO is bounded/weak here — the swarm × particle × refine cost makes a
large iteration budget pathological, so the test clamps it). Generate test ciphers with
`tools/fracmorse_gen.c` (`make fracmorse_gen`; args are a plaintext and a keyword).

The **Pollux** type (`pollux`/`pol`/`74`; `solve_pollux()` + `pollux_search()` in `pollux_solver.c`)
is the ACA **"Pollux"** — the second **Morse** type. The plaintext is written in International Morse
with a single `x` separator between letters and `xx` between words (a run of `xxx` is impossible), and
the three stream symbols `{DOT, DASH, X}` are mapped to a **DIGIT stream** by a **KEY that assigns each
of the ten digits 0..9 to one of the three symbols** (usually 4 dots / 3 dashes / 3 x, but any split
with ≥1 of each is legal). So the ciphertext is a run of digits 0-9, **one digit per Morse symbol** — a
length CHANGE (C digits → N letters), **no period**, full 26-letter alphabet. Encryption is
**POLYPHONIC** (each Morse symbol becomes a *random* digit among those the key assigns to it);
decryption is fully **deterministic** given the key. The primitive (`pollux.c`, hand-verified against
the ACA worked example — key `1→x 2→- 3→. 4→. 5→x 6→. 7→- 8→- 9→x 0→.`, `LUCK HELPS` → CT
`08639 34257 02417 68596 30414 56234 90874 5360`) reuses **a small self-contained Morse table** (the
`fracmorse.c` reverse-token trick, copied so `fracmorse.c` stays byte-identical, the seriated_playfair
precedent). **The decisive design point:** the entire keyspace is **3^10 = 59,049** digit→symbol maps,
so — like the **Period column order** solver — Pollux is a **DETERMINISTIC EXHAUSTIVE** solver, **NOT a
hill-climbing `CipherModel`** (no engine schedule / `SearchDefaults` / `-method`): `pollux_search`
decodes + scores **every** key and keeps the global best, guaranteed optimal in well under a second (~30M
ops). This is the right shape because the landscape is a **needle** (re-assigning one digit re-parses the
whole Morse stream → no local gradient), exactly where a stochastic climber flails but enumeration is
certain. Fitness is the shared **mean** n-gram score of the decoded plaintext (length-fair for free, so
variable-length decodes compare directly — **no tiling**, unlike fracmorse which tiles only because it
runs inside the fixed-length engine) plus a **Morse-VALIDITY reward** (`POLLUX_VALID_WEIGHT * n_valid/
n_tokens`; illegal `xxx` runs and >4-symbol codewords count against it), the fracmorse analogue — it
breaks near-ties and, at the very short end, is essential (measured: at 16 letters, validity weight 0 →
0% recovery, weight 3 → 31%, weight 6 → 98%). It **rides the reward-only quadgram table (NO `-logprob`)**
like period_column — and unlike the fracmorse family, `-logprob` **degrades** short-text recovery (the
validity weight is tuned to the reward-only scale). It recovers reliably **from ~24 letters**, far below
the ACA 80-100 range (see `tests/test_pollux_solver.c`, which prints the length cliff and the
validity-weight / `-logprob` calibration sweeps). **Cribs are not used** (the length change breaks the
positional map, like fracmorse). Digit input is parsed from `ciphertext_str` directly (mirroring
`solve_nihilist_sub`'s number parser, but one digit per symbol), not the A-Z decode; the report emits a
fracmorse-style `>>>` CSV with the recovered `key=` (a 10-char `.-x` string) and `valid=`. Dispatched by
its own early branch in `solve_cipher`. Generate test ciphers with `tools/pollux_gen.c`
(`make pollux_gen`; args are a plaintext, a 10-char `./-/x` assignment indexed by digit 0..9, and an
optional RNG seed). Unit tests: `tests/test_pollux.c` (primitive — the ACA decode KAT, the Morse-stream
KAT via the `xx` word divider, round-trips, and the illegal-`xxx` / long-codeword / missing-element
edges) and `tests/test_pollux_solver.c` (exact recovery, the length cliff, and the validity-weight +
`-logprob` calibration sweeps); `ciphers/tests/run_tests.sh` carries `pollux_pp` (a ~90-letter cipher,
~100%).

The **Period column order** type (`period-column`/`periodcol`/`pcol`/`transpercol`/`71`;
`solve_period_column()` + `period_column_search()` in `period_column_solver.c`) is AZdecrypt's
**"Period column order"** transposition — a **columnar transposition whose column permutation is fixed
by a PERIOD, not a keyword**. The message is laid row-major into a `dx`-wide grid of `dy = ceil(len/dx)`
rows; the columns are visited in the periodic order `0, p, 2p, … , 1, 1+p, … , 2, …` (`for i in 0..p-1:
for j = i; j < dx; j += p`) and the k-th visited column receives (`utp==0`, "TP") / supplies (`utp==1`,
"UTP") the natural k-th column; the grid is read back row-major, **skipping the padding cells**, so it is
exactly length-preserving. The primitive (`period_column.c`, a faithful port of AZdecrypt's case 13,
hand-verified cell-for-cell — the 168-char worked example `I LIKE KILLING PEOPLE …` decrypts via
`[4x42,TP,P:3]` then `[56x3,UTP,P:2]`) carries **spaces/punctuation as genuine grid cells** (they ride
through the permutation like letters), so it composes with the transposition pipeline; `PERIOD_COLUMN` is
therefore in the **space-significant** set (a trailing space is a real cell). Because a single stage's key
is tiny — a **complete-grid width `dx | len` in `[3, len]`**, a **period `p ∈ [2, dx-1]`**, and the TP/UTP
flag — the whole single-stage keyspace is enumerable and a **composition of two stages** is enumerable as
the product; so the solver is **DETERMINISTIC and EXHAUSTIVE**, NOT a hill climber: it tries every stage
(depth 1) and every ordered pair (depth 2), scores each decrypt with the shared reward-only quadgram
fitness (`state_score`, **no `-logprob`** — a transposition preserves letter frequencies), and keeps the
**global best** (`-depth` pins 1 or 2, default 2 — the depth-2 search includes all depth-1 solutions). This
is the right shape because a wrong stage yields **no n-gram gradient** (a needle), so a stochastic climber
would flail where exhaustive enumeration is guaranteed and, for a few-hundred-char cipher, sub-second
(depth 2 is `O(K²)` decrypt+score over `K` ≈ a few hundred to a couple thousand stages; `PC_MAX_STAGE_LIST`
clamps pathologically composite lengths). It is a **dedicated solver, NOT a `CipherModel`** (no engine
schedule / `SearchDefaults` entry — there is nothing to anneal), dispatched by its own early branch in
`solve_cipher`. **Widths are complete-grid divisors only** (`dx | len`): the primitive handles incomplete
grids for faithfulness, but per-stage inversion (`utp==1` == inverse of `utp==0`) is exact **only** on a
complete grid, so the solver restricts to those (the vast majority of real puzzles; the motivating cipher is
`6×28 = 168`). Cribs are supported positionally (blended via `state_score` when supplied) but not required.
This is the solver that **reproduces AZdecrypt end-to-end** on the transposed Zodiac-style plaintext that
colossus's free-permutation and keyword-columnar solvers could not (they don't model the composed
period-column structure). Generate test ciphers with `tools/period_column_gen.c` (`make period_column_gen`;
args are a plaintext then one or more `<dx> <period> <tp|utp>` stage triples; it carries spaces, so the
plaintext may contain them). Unit tests: `tests/test_period_column.c` (primitive — the tiny hand KAT, the
full AZ 168-char composition KAT with spaces, exact round-trip over random complete grids, multiset
preservation over incomplete grids) and `tests/test_period_column_solver.c` (depth-1/depth-2 exact
recovery, the AZ worked example asserted end-to-end with the recovered stages, and a length cliff);
`ciphers/tests/run_tests.sh` carries `period_column_pp` (a two-stage 168-letter case, ~100%). Both this
solver and the space-robust one below honour **`-readdir tb|bt|both`** — the period-column read applies to
the whole cipher **stream** (forwards = `tb`, reversed = `bt`), `both` runs each orientation and keeps the
higher n-gram score (the winning direction is tagged `dir=bt` in the report/CSV params). Default `tb` is one
forward search, so it is **bit-identical** to before (`period_column_pp` still 100%).

The **Period column order, space-robust** type (`period-column-space`/`pcol-space`/`pcolspace`/`pcolsp`/
`transpercolspace`/`72`; `solve_period_column_space()` + `period_column_space_search()` in
`period_column_space_solver.c`) generalises the deterministic solver above to **POORLY enciphered**
ciphertext — one where a few characters were accidentally **dropped** or **added**. The observed cipher is
kept **exactly, spaces and all** (the spaces are structurally significant — they ride the transposition as
real grid cells and a single misplaced space re-shuffles every downstream letter — so they are **never
stripped or moved**). Instead the solver applies two **symmetric** edits, both searched **from scratch**:
**INSERT** up to `-maxgaps` blank/gap cells (a synthesised space; repairs a **dropped** char — the observed
stream is too short — by restoring the grid coordinates of every following cell) and **DELETE** up to
`-maxdels` observed cells (repairs a spuriously **added** char — the stream is too long). A gap is
transparent to the letters-only n-gram fitness but its POSITION changes the grid layout and hence the
decrypt; both edit kinds also let the edited length reach a richer complete-grid factorisation. Because the
edit-position space is combinatorial, this is a **stochastic** solver (unlike the deterministic base): for
each `(n_ins, n_dels)` pair in `[0,maxgaps] × [0,maxdels]` — `(0,0)` is exactly the plain period-column
search, so it never scores below it — the edit positions are annealed by a restart+anneal hill climb whose
**fitness is the fast exhaustive period-column search** (depth-1 inner), and the single best edited cipher
gets a full **`-depth` (≤2) exhaustive finish**. Rides the reward-only quadgram table (no `-logprob`);
optional `-weightword` blends dictionary word-coverage. Honours `-nrestarts`/`-nhillclimbs`/`-inittemp`/
`-mintemp` (the stochastic knobs) and `-readdir` (as above). **Cribs are not used** (inserted/deleted cells
shift plaintext positions). It is a **dedicated solver, NOT a `CipherModel`**, with its own early branch in
`solve_cipher`. Verified on synthetically corrupted ciphers: a dropped letter recovers ~all-but-one position
by inserting one gap; a spuriously added letter recovers **100%** by deleting exactly that cell (restoring
the original cipher). Unit tests: `tests/test_period_column_space_solver.c` (the `(0,0)` clean-cipher
superset check, dropped-1/dropped-2 insertion repair, and added-1 deletion repair — all over a 168-cell
spaced plaintext at a fixed seed). No generator of its own (corrupt a `period_column_gen` cipher). No
`run_tests.sh` case (recovery of a corrupted cipher lands just under the char-for-char ~99% threshold — the
deterministic solver test is the guard).

The **Double columnar transposition, divide & conquer** type
(`transcol2-dc`/`transcol2dc`/`dctrans`/`doublecol-dc`/`dcol`/`73`; `solve_double_transposition()` +
`dct_solve_core()` in `double_transposition_solver.c`) is a dedicated solver for the **double columnar
transposition** cipher — `C = colenc(colenc(P, K1), K2)` — that attacks the two keys **SEPARATELY**,
following Lasry, Kopal & Wacker, *"Solving the Double Transposition Challenge with a Divide and Conquer
Approach"* (Cryptologia, 2014). Where `transcol2` (type 18) hill-climbs `K1` and `K2` **in parallel** over
the combined `K1!·K2!` keyspace — the paper's "Step 1", which they measured to work only for key lengths
**≤ ~15** and which fails on the famous 599-letter, key-length-21/23 challenge — this solver uses the
**Index of Digraphic Potential (IDP)** to score the SECOND key `K2` **without knowing `K1`**: undo the
second transposition with a putative `K2`, and the residue is a *single* columnar transposition of the
plaintext under `K1`, whose "reconstructability" is scored by greedily chaining the residue's `K1` columns
into a left-to-right order by **digraph fit** (a max-weight Hamiltonian path over a column-adjacency matrix
`B[i][j]`, whose entries are the best per-row bigram log-frequency of placing column `j` right of column
`i`, maximised over the incomplete-rectangle **start-position windows**). Higher IDP ⇒ `K2` closer to
correct. So the `O(K1!·K2!)` joint search collapses to an `O(K2!)` **IDP hill-climb for `K2`** followed by
an `O(K1!)` **single-columnar finish for `K1`**. The primitive is `decrypt_columnar` (`transpositions.c`),
so the solver adds no new cipher math; the IDP needs a **26×26 bigram log-frequency table**, built by
`dct_load_bigrams` by marginalising the `-ngramfile` (no separate bigram file). **Key design points:**
(1) the IDP is returned **per-edge normalized** (÷ `L1−1`) so it is comparable across `L1` hypotheses — a
raw sum trivially favours fewer columns and mis-ranks the length sweep; (2) a columnar transposition is
recoverable only **up to a cyclic rotation** (the Hamiltonian path has free endpoints and the last column's
spurious "best right neighbour" closes a near-cycle), so the finish **tries all `L1` rotations of the IDP's
own greedy chain** (and its reversal) as columnar orders — the correct rotation cuts the near-cycle at the
true first/last-column boundary and recovers `K1` almost exactly when `K2` is correct, giving the finish a
near-deterministic warm start that only needs a light n-gram polish; (3) the anneal temperature is matched
to the per-edge IDP scale (~-2.3, single-swap deltas ~0.02-0.2) — a schedule tuned to the un-normalized sum
is ~`L1`× too hot and random-walks the budget. The pipeline is **SCREEN** (a short IDP climb of every
`(L1, L2, dir)` config in `-mincols..-maxcols²`, ranked by best IDP) → **REFINE** (a deep IDP climb of the
top `keep`=8) → **FINISH** (rotation-resolved single-columnar break of each, kept by the real n-gram score,
which — not the screening IDP — makes the final call across the kept set). It is a **dedicated solver, NOT a
`CipherModel`** (the two-phase K2/K1 structure has no single annealed state), dispatched by its own early
branch in `solve_cipher`; `-nrestarts`/`-nhillclimbs` scale the refine budget (screening depth is a fraction
of it) and `-readdir tb|bt|both` sets the read direction (both stages share it). The K2 hill-climb uses the
paper's **Step-4 refinements**: segment-wise moves (segment slide / segment swap / 3-partite rotation, more
effective than single swaps on transposition keys) plus a **left-to-right best-improvement heuristic**
applied to each restart's seed and as a final polish. **Cribs ARE supported** (`-crib`, over **plaintext**
positions) and — crucially — **guide the K2 climb itself**: since undoing K2 leaves a single columnar the
IDP's greedy chain reconstructs up to a rotation, a crib-aware K2 fitness cheaply rotation-resolves K1 from
the chain and folds the best crib match into the objective (`dct_k2_fitness`, `crib_w 1.0`). Far from the
solution every K2's crib fraction sits at the floor (the IDP drives), but as K2 nears the truth the crib
fraction climbs steeply and sharpens onto the exact key — so a known-plaintext fragment lets the solver
crack a cipher at a budget too small for pure IDP (verified: the true K2 rotation-resolves to crib fraction
1.0 vs ~0.22 for random keys). The crib also blends into the **finish + rotation scoring** (`state_score`,
default `weight_ngram 12` / `weight_crib 36`) to lock K1 and select across configs. It rides the reward-only
quadgram table but **`-logprob`**
sharpens the finish. **Spaces/punctuation are carried as REAL grid cells** (like the period-column solvers,
`TRANSCOL2_DC` is in the **space-significant** set so a trailing space is a real position): when the
transposition grid includes the spaces/punctuation (they ride both columnar stages), the observed cipher
keeps them as reversible negative sentinels, `decrypt_columnar` permutes them opaquely, and they are printed
back in place by `index_to_char`. Both scorers are **letters-only**: the IDP skips any column-pair row touching
a non-letter (mean bigram log-freq over the LETTER pairs only — `nv==q` for an all-letter cipher, so it is
bit-identical to the old sum/q) and `ngram_score` already projects onto letters. So a double transposition whose
plaintext contains spaces (e.g. the W168 hypothesis) is attacked with the spaces intact (no `-skipspaces`);
recovery is exact including the punctuation. Under **`-verbose`** it streams a live dialog in the polyalphabetic solver's
`value⇥[label]` convention: per-config screening lines, refine phase markers, live best-IDP improvements
during the `K2` climb, and a metrics block + the full decrypted plaintext on each new global best. It cracks
key lengths **(20-25)** at which `transcol2`'s parallel hill-climb fails; the IDP is **selective on the real
599-letter challenge** (the true `K2` is the global IDP optimum, `L1=21` identifiable), so given enough
search budget it recovers the challenge plaintext. Generate test ciphers with `tools/double_transposition_gen.c`
(`make double_transposition_gen`; args are a plaintext and two keywords). Unit tests:
`tests/test_double_transposition.c` (the IDP's paper properties — SELECTIVITY: true `K2` > random and peaks
at the true `L1`; MONOTONICITY: IDP degrades as `K2` is perturbed, Fig 2; plus the round-trip and the
rotation-resolved greedy `K1` reconstruction) and `tests/test_double_transposition_solver.c` (end-to-end
recovery of several key-length pairs 6-10, a length cliff, a crib-guided-K2 case, and a **space-significant
case** — a double transposition over text with spaces/punctuation recovered 100% incl. the punctuation, all at
fixed seeds — the regime well past
`transcol2`'s ~12-15 ceiling). No `run_tests.sh` case (the solver is stochastic and each solve is tens of
seconds; the two unit tests are the guard).

The **Slidefair** family (`slidefair`/`sf`/`59`, `slidefair-var`/`sfv`/`60`, `slidefair-beau`/`sfb`/`61`;
`solve_slidefair()`, `SLIDEFAIR_MODEL`, `SHAPE_ANNEAL`) is the ACA **periodic DIGRAPHIC
Vigenère/Variant/Beaufort** — the digraphic sibling of Portax (rectangle-over-a-slide applied to the
three shift tableaus instead of the Porta tableau). The plaintext is taken in consecutive **digraphs**;
a **keyword** fixes the period **P**, and digraph i (cipher letters 2i, 2i+1) is keyed by `keyword[i%P]`
over a **two-row slide**: the TOP row is the standard alphabet (`top[col]=col`), the BOTTOM row a shift
alphabet — Vigenère `bottom[col]=(col+k)`, Variant `(col−k)`, Beaufort `(k−col)` (the three are distinct
`-type` codes sharing one primitive `slidefair.c` + solver, branched on `cfg->cipher_type`). A plaintext
pair (p1 in the top row at column p1, p2 in the bottom row at the column where `bottom==p2`) forms diagonal
corners of a 2-row rectangle; the **substitutes are the other two corners, the TOP one first**
(`c1=top[col2]`, `c2=bottom[col1]`); a **vertical pair** (both in one column) maps to the vertical pair
**one column to the right** (decrypt: to the **left** — the only place decrypt ≠ encrypt, the rectangle
case being self-reciprocal). Full 26-letter alphabet (**no J→I merge**); the primitive is hand-verified
cell-for-cell against the ACA worked examples (key B: `ca→ZD/BB/BZ`, `de→EF/FC/XY`; keyword `DIGRAPH`,
`the slidefair can be used... → EWKMCRNUAFCXTJ...`). **The whole key is P key letters (0..25), a 26^P
space**, carried in the cycleword lane, with **one engine config per swept period P** (IoC period
estimation is useless through the digraphic pairing, so P is swept and the n-gram score picks it — a wrong
P decrypts to gibberish; `-period` pins one, `-mincols`/`-maxcols` bound the sweep). **The key to
efficiency is per-column independence:** a digraph is enciphered ENTIRELY by its column key, so every
digraph in column c decrypts from `key[c]` alone — the **per-column monogram-fit key** (the analog of
`derive_optimal_cycleword`) **warm-starts the seed**, and the n-gram (quadgram) score then drives the
anneal (cross-column digraphs only form at the true keys) and corrects any column the monogram fit mis-set.
No `score_adjust` is needed (every cycleword is a valid bijective decrypt). Vigenère and Variant are
**not separately identifiable** (a free per-column key absorbs the sign), so either solver cracks a
shift-Slidefair; only Beaufort is distinct. Like the rest of the Vigenère/Porta family it **rides the
reward-only quadgram table (no `-logprob`)** and **recovers cleanly from very short text** (~100% from
~50 letters — the 26-value per-column monogram fit is a strong warm start; see
`tests/test_slidefair_solver.c`). **Cribs are supported** (the cipher is positional: `decrypted[i]` is
plaintext position i). Because the model implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso`
all run on it. Generate test ciphers with `tools/slidefair_gen.c` (`make slidefair_gen`; args are a
plaintext, a keyword, and `vig`/`var`/`beau`).

The **Seriated Playfair** type (`seriated-playfair`/`serpf`/`spf`/`62`; `solve_seriated_playfair()`,
`SERIATED_PLAYFAIR_MODEL`, `SHAPE_ANNEAL`) is the ACA **"Seriated Playfair"** — plain Playfair over a
**single 5x5 keyed square**, but the digraphs are the **VERTICAL PAIRS of a two-row seriated layout** of
period **P** instead of consecutive horizontal pairs. The plaintext is laid into blocks of **2P** letters
(first P = top row, next P = bottom row); within a block, vertical pair `j` couples block-letter `j` (top)
with `j+P` (bottom) for `j = 0..P-1`, each pair enciphered by the three standard Playfair rules and written
**back to the same positions**, and the cipher is the blocks serialized left-to-right (top row then bottom —
the "taken off horizontally" readout). The whole cipher is thus
`for each 2P block, for j: (out[j], out[j+P]) = playfair_pair(in[j], in[j+P], dir)`. Runs on the same
25-letter (J→I) grid as Playfair (`init_alphabet("J")` before `load_ngrams`); the primitive
(`seriated_playfair.c`, hand-verified cell-for-cell against the ACA worked example — square
`LOGARITHMBCDEFKNPQSUVWXYZ`, period 6, `comequickly… → NLBCSP…`) **reuses** `playfair.c`'s
`playfair_build_inverse`/`playfair_grid_from_keyword` (only the 3-rule pair is a small inlined copy so
`playfair.c` stays byte-identical). **The key design point is what it is NOT:** unlike Portax/Slidefair there
is **NO per-column independence** — one square enciphers every pair — so there is **no
`derive_optimal_cycleword`-style decoupling**. The attack is therefore **exactly Playfair's single-grid
anneal** (random-grid seed; cell-swap-dominated moves + row/column swaps + reflections; bijection, no
anti-collapse penalty) with **one addition: the seriation period P is SWEPT** (IoC is useless through the
pairing, so one engine config per P and the n-gram score picks the true one — only the exact P pairs
correctly, a multiple does not; `-period` pins, `-mincols`/`-maxcols` bound the sweep, default top 15). Each
config is a full Playfair-scale grid anneal, so the blind sweep **multiplies the Playfair cost**; like every
square type it effectively **needs `-logprob`** and recovers reliably from **~250+ letters** (see
`tests/test_seriated_playfair_solver.c`, which prints the length cliff). On rare grids the recovered solution
hits the inherent **rare-letter ambiguity** of a square attack (the filler `X` appears only as nulls, so it
is weakly constrained and a few X-positions may flip — a deterministic ~92% ceiling no budget escapes, as in
Playfair). The grid is unique only up to a cyclic row/column rotation (all re-decrypt identically). **Cribs
are not used** (the prepare null-insertion shifts plaintext positions, like Playfair). Because the model
implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it. The encrypt/prepare path
(null insertion on a doubled vertical pair → a filler `X`, or `Q` if the letter is `X`; final block padded)
serves only the generator + unit tests; the solver only ever decrypts. Generate test ciphers with
`tools/seriated_playfair_gen.c` (`make seriated_playfair_gen`; args are a plaintext, a keyword, and a period).

The **Digrafid** type (`digrafid`/`df`/`dgf`/`63`; `solve_digrafid()`, `DIGRAFID_MODEL`, `SHAPE_ANNEAL`)
is the ACA **"Digrafid"** — a **digraphic fractionation cipher over TWO independently keyed 27-symbol
alphabets** (A..Z + `#`, so `init_alphabet_digrafid()` is forced before `load_ngrams`, like Trifid's
A..Z+`+`; **no J→I merge**). The tableau is a **horizontal** grid `H` (3 rows × 9 cols, keyed alphabet
entered **row-major**) and a **vertical** grid `V` (9 rows × 3 cols, entered **column-major**). A plaintext
**digraph** (a, b) → a **3-digit number** (top, mid, bot), each 0..8: `top` = a's column in `H`, `bot` =
b's row in `V`, `mid` = (a's row in `H`)·3 + (b's col in `V`) (the 3×3 intersection block) — a bijection
between the 729 digraphs and 729 triples. **Fractionation** (period = digraphs per group): within each
group of `g` digraphs the `g` triples are stacked as 3 rows (tops/mids/bots), read **row-major** into a
`3g`-digit stream, re-split into `g` consecutive triples, and each new triple mapped back through the
tableau to one ciphertext digraph — exactly the Trifid reshape, but over digraphs. The primitive
(`digrafid.c`, hand-verified cell-for-cell against the ACA worked example — keywords `KEYWORD`/`VERTICAL`,
`THISISTHEFORESTPRI` → period 3 `HJMXWSWJADWGFCSPYI`, period 4 `HJTKVHYUFFWDSQYPRI`) **reuses** `bifid.c`'s
`bifid_build_inverse` and is O(len). **The key design point: the two grids are KEYED ALPHABETS, and the
solver searches them AS SUCH** — not as free 54-cell permutations. An ACA Digrafid grid is a keyword
(duplicates dropped) followed by the rest of the alphabet ascending (exactly what `digrafid_grid_from_keyword`
builds), so the state is **two keyed-alphabet SEQUENCES** (`H = key[0..26]`, `V = key[27..53]`), each kept as
**"keyword prefix of length `kw` + sorted tail"** (the keyword lengths live in `st->aux[0..1]`). The move
(`digrafid_move_seq`, the same keyed-alphabet structure the Quagmire solver uses via `random_keyword`/
`perturbate_keyword`, but re-weighted for the coarser fractionation landscape) perturbs ONE sequence: ~4% a
keyword-length change (re-canonicalising the tail with `digrafid_canonicalize`), ~48% a keyword↔tail swap
(coarse set search, re-sorts the tail), ~48% an in-keyword swap (smooth reorder); the grids are rebuilt from
the sequences at decrypt (`digrafid_grid_from_keyword` with the whole sequence = a pure sequence→grid map,
H row-major / V column-major). Every keyed alphabet is a bijection, so there is **no anti-collapse penalty**
(`score_adjust` stays 0). **Why this matters:** collapsing each grid from a free 27! permutation to a short
keyword shrinks the keyspace by orders of magnitude and makes the sorted tail a strong prior, so recovery no
longer needs ~700–800 letters (the old free-permutation square break's cliff) — it **tracks the keyword and
recovers reliably from ~300 letters** (the way an ACA solver does; on the ACA index corpus this flips the
shorter puzzles from unsolvable to solved). A keyword longer than `DIGRAFID_KW_MAX` (13), or a non-keyed/
fully-shuffled grid, is not representable — but every ACA Digrafid, and the generator, uses a keyed alphabet.
A **fractionation period is SWEPT** on top (IoC is useless through the digraphic pairing, so one engine config
per candidate period and the n-gram score picks the true one — a wrong period regroups → gibberish). The
period is recovered by `digrafid_estimate_periods` — the mean per-lane **Index of Coincidence** over the 2P
lanes (each (digraph-position-in-group, first/second-letter role)) peaks at the true period (and at its
multiples, which lose on the n-gram score), exactly like Bifid's columnar-IoC estimator; top-`-nperiods`
(default 5) annealed, `-period` pins one, `-maxperiod` bounds the scan. Like the other square types it
effectively **needs `-logprob`**; **RESTARTS are the lever** (each samples/refines a keyword length, and the
coarse keyword moves want a **warm temperature**, so the schedule is many restarts at `inittemp 0.30`).
**Cribs are not used** (the fractionation couples a crib to its whole group → weak gradient, as in
Bifid/Trifid). Because the model implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on
it. Generate test ciphers with `tools/digrafid_gen.c` (`make digrafid_gen`; args are a plaintext, a
horizontal keyword, a vertical keyword, and a period).

The **CM Bifid** type (`cm-bifid`/`cmbifid`/`cmb`/`64`; `solve_cm_bifid()`, `CM_BIFID_MODEL`, `SHAPE_ANNEAL`)
is the ACA **"Conjugated Matrix Bifid"** — plain Bifid (Delastelle fractionation, block size = period) but the
two coordinate lookups use **TWO different keyed 5x5 squares**: square 1 fractionates each plaintext letter
into its (row, col) coords, and after the standard rows-then-cols reshape and consecutive re-pairing each
coordinate **pair** is mapped to a ciphertext letter through a **second** square (decrypt mirrors: expand the
cipher letters via sq2's inverse, recombine via sq1). When `sq1 == sq2` it is **exactly Bifid** (a property the
unit tests pin). Runs on the same 25-letter (J→I) grid as Playfair/Bifid (`init_alphabet("J")` before
`load_ngrams`); the primitive (`cm_bifid.c`, hand-verified cell-for-cell against the ACA worked example —
squares `EXTRAKLMPOHWZQDGVUSIFCBYN` (pt) / `NCDRSOBFQUVAGPWEYHMXLTIKZ` (CT), `ODDPERIODSAREPOPULAR` period 7 →
`FANXZEXFENUKKRBYNKAK`) **reuses** `bifid.c`'s `bifid_build_inverse` and is side-generic (a 6x6/36-cell pair
works once a 36-letter alphabet is active). **The key design point: there is NO square-independent decoupling
reward** (both squares are entangled in the n-gram fitness — unlike ADFGVX's transposition-only IoC or
Nihilist-sub's additive-only validity), so the attack is the proven **JOINT two-square anneal** (cf.
Two/Four-Square): the state is the pair of squares packed back-to-back in the `key` lane (sq1 = `key[0..24]`,
sq2 = `key[25..49]`), each move perturbing **one** square (chosen uniformly) with the Bifid/Playfair move set
(cell-swap-dominant + row/column swaps + reflections), `score_adjust` stays 0 (every square is a bijection), and
like every square type it effectively **needs `-logprob`**. The **period is recovered exactly as Bifid's** —
`bifid_estimate_periods` is square-AGNOSTIC (square 2 only relabels the coordinate pairs and columnar IoC is
relabel-invariant), so it is **reused unchanged** and the period is SWEPT (top-`-nperiods` annealed, `-period`
pins, `-maxperiod` bounds; the n-gram score picks the winner). **The headline property is ODD-vs-EVEN period:**
with an **EVEN** period the rows-then-cols re-paired stream splits cleanly — the first P/2 output pairs are
pure-ROW coordinates and the last P/2 pure-COLUMN — so a fractionation row and column **never share an output
pair**, leaving a transpose-like square ambiguity that makes the squares **degenerate ciphertext-only** (a
different square pair decrypts to equally-English text); recovery of the planted key then collapses to the
noise floor (~5%) **no budget or text length escapes it**. With an **ODD** period the boundary pair **mixes** a
row and a column, breaking the symmetry, and recovery is clean — reliable from **~480 letters** (a ~400 cliff;
the long per-restart climb is the critical lever, so the schedule is `8x400000`, not many short restarts). This
parity effect is the analog of Playfair's rare-letter-X ceiling: an inherent property, not a solver weakness.
**Cribs are not used** (the fractionation couples a crib to its whole group → weak gradient, as in
Bifid/Trifid). Because the model implements `seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on
it. Generate test ciphers with `tools/cm_bifid_gen.c` (`make cm_bifid_gen`; args are a plaintext, a square-1
keyword, a square-2 keyword, and a period).

The **Tri-Square** type (`trisquare`/`tri-square`/`3square`/`3sq`/`trisq`/`65`; `solve_trisquare()`,
`TRISQUARE_MODEL`, `SHAPE_ANNEAL`) is the ACA **"Tri-Square"** — a digraphic substitution over **THREE
independent keyed 5x5 squares** (25 letters, J→I, same alphabet as Playfair/Two/Four-Square, so
`init_alphabet("J")` is forced before `load_ngrams`). Plaintext is taken in pairs: p1 is located in square 1
at (r1,c1), p2 in square 2 at (r2,c2), and the digraph enciphers to a ciphertext **TRIGRAPH** (a **3:2 length
expansion**, so an N-letter plaintext yields 3N/2 ciphertext letters): `c0` = **any** letter in p1's **column**
in sq1, `c1` = `sq3[r1][c2]` (the deterministic middle letter), `c2` = **any** letter in p2's **row** in sq2.
The first and third letters are **polyphonic** on encode (the ACA lets the clerk pick any of the 5 column/row
members — the primitive `trisquare.c` picks a **random** representative under the RNG, faithful to the ACA and
far better-conditioned for the solver than any fixed canonical choice, which would starve the gradient by
concentrating c0/c2 into a size-5 subset); **decryption is exact for any choice** (a square maps any column
member back to its column, any row member back to its row — `col1 = pos1[c0]%side`, `(row1,col2) = pos3[c1]`,
`row2 = pos2[c2]/side`, then `p1 = sq1[row1][col1]`, `p2 = sq2[row2][col2]`). The primitive is hand-verified
cell-for-cell against the ACA worked example (squares `NSFMUOAGPWVBHQXECIRYLDKTZ` / `READINGBCFHKLMOPQSTUVWXYZ`
/ `PASTINOQRMLYZUEKXWVBHGFDC`, `RHLQXR…AAABFZ` → `THREEKEYSQUARESUSEDX`) and **reuses** `bifid_build_inverse`;
it is side-generic (a 6x6/36-cell triple works once a 36-letter alphabet is active). **The search state is the
three squares packed back-to-back** in `st->key` (sq1 = `key[0..24]`, sq2 = `key[25..49]`, sq3 = `key[50..74]` —
75 cells, the **largest square state of the family**); there is **no square-independent decoupling reward**
(every square is a bijection, so any triple decrypts to a permuted-but-valid stream and only JOINT correctness
yields English n-grams), so the attack is the proven **joint multi-square anneal** (cf. Four-Square / CM-Bifid) —
each move perturbs **one** of the three squares (chosen uniformly) with the classic Playfair move set (cell-swap
dominant + row/column swaps + reflections, no cyclic rotation), `score_adjust` stays 0, one engine config (no
period to estimate), effectively **needs `-logprob`**. **The 3:2 length change is handled exactly like ADFGVX:**
the raw trigraph stream lives in the scratch, the solver passes the **plaintext/scoring length `n = 2M`** (not
`3M`) to `make_solver_ctx`, and the decrypt hook emits `n` plaintext symbols the engine n-gram-scores. Cribs are
not used (positions are over the trigraph stream and the encode is polyphonic). The squares are recoverable only
up to the cipher's structural symmetries (sq1's columns / sq2's rows are unconstrained by the c0/c2 redundancy),
but the recovered **plaintext is unique**. **Counter-intuitively it recovers MORE easily than Four-Square**
despite the bigger state: the polyphonic c0/c2 letters spread the full alphabet over every ciphertext position,
giving the n-gram gradient sharp signal, so it is **reliable from ~500 plaintext letters (~750 cipher)** with a
**300-400 cliff** — comfortably covering the ACA "100-125 groups" range. Because the model implements
`seed`/`perturb`/`copy`, `-method anneal|shotgun|pso` all run on it. Generate test ciphers with
`tools/trisquare_gen.c` (`make trisquare_gen`; args are a plaintext and three keywords; the generator seeds the
RNG so the polyphonic ciphertext is reproducible).

**Per-cipher-type search schedules (`SearchDefaults`, `apply_cipher_defaults`).** The
`init_config()` globals (`inittemp 0.10`, `1x1000`, ...) suit the polyalphabetic /
transposition reward-score scale; a type whose score lives on a very different scale
needs its own schedule. A small compiled-in registry (`g_search_defaults[]` in
`colossus.c`) keyed by cipher type carries a tuned profile for **all three** search shapes
(anneal `a_*` + shotgun `s_*` + particle-swarm `p_*`; a zero `p_n_particles` means "no PSO
profile, keep globals"). `main()` pre-scans `-type`/`-method` and overlays the matching profile
*before* the main arg loop, so precedence is **globals < registry < explicit CLI flags**.
Types with no entry keep the global defaults bit-for-bit (so the regression suite is
unaffected) — currently Playfair (`SHAPE_ANNEAL`, `6x400000`, `inittemp 0.08`,
`backtrack 0.30`), Bifid (`SHAPE_ANNEAL`, `4x200000` per period, `inittemp 0.08`,
`backtrack 0.30`), Trifid (`SHAPE_ANNEAL`, `6x300000` per period, `inittemp 0.08`,
`backtrack 0.30` — a larger budget for the 27-cell cube), Hill (`SHAPE_ANNEAL`,
`250x8000` per swept block size, `inittemp 0.10`, `backtrack 0.25` — many short restarts,
since the small matrix climbs converge fast), Phillips and its two variants
(`SHAPE_ANNEAL`, `4x250000`, `inittemp 0.08`, `backtrack 0.30` — a single config, leaner
than Playfair since monographic Phillips recovers from shorter text), Two-Square and its
vertical variant (`SHAPE_ANNEAL`, `8x600000`, `inittemp 0.08`, `backtrack 0.30`) and
Four-Square (`SHAPE_ANNEAL`, `12x700000` — the biggest budget, for its two independent
keyed squares), ADFGX (`SHAPE_ANNEAL`, `12x600000` per swept column count K) and ADFGVX
(`SHAPE_ANNEAL`, `16x800000` per K — more, for the 36-cell square; both anneal at
`inittemp 0.08`, `backtrack 0.30`) have tuned entries (the two/four-square budgets are
larger than Playfair's for the 50-cell two-square state; the ADFGVX budgets are the largest,
for the coupled square+columnar search), the three Nihilist Substitution conventions
(`SHAPE_ANNEAL`, `8x300000` per period, `inittemp 0.08`, `backtrack 0.30` — between Bifid and
ADFGX, for the coupled square+additive search), Gromark (`SHAPE_ANNEAL`, `3x120000` per top-K
primer config — a lean per-config 26-letter substitution anneal, since the primer pre-pass and the
provisional-σ warm start do most of the work) and Periodic Gromark (`SHAPE_ANNEAL`, `4x160000` per
swept period — the keyword anneal over a ~28-bit key; both at `inittemp 0.08`, `backtrack 0.30`),
and the three Nicodemus codes (`SHAPE_ANNEAL`, `16x20000` per swept `(P, H)` pair, `inittemp 0.08`,
`backtrack 0.30` — many short restarts, since the climbed state is just a short column-order
permutation, so restarts are the robustness lever, not climbs), and Bazeries
(`SHAPE_ANNEAL`, `40x20000` per swept digit count D, `inittemp 0.08`, `backtrack 0.30` — many
restarts, since the climbed state is a short digit string over a rugged < 10⁶ keyspace, so each
restart reseeds a fresh random number and restarts carry the robustness), and Portax
(`SHAPE_ANNEAL`, `12x20000` per swept period P, `inittemp 0.08`, `backtrack 0.30` — a lean budget,
since the monogram-fit warm start gets most of the short per-column-shift state right on seed and
the anneal/n-gram pass only corrects a few columns), and the three Progressive Key bases
(`SHAPE_ANNEAL`, `3x2500` per swept `(P, prog)` pair, `inittemp 0.08`, `backtrack 0.30` — a very
lean per-config budget, since MANY `(P, prog)` configs are enumerated and the per-column monogram
warm start already gets most columns right, so a few short restarts suffice), and the three Slidefair
bases (`SHAPE_ANNEAL`, `8x10000` per swept period P, `inittemp 0.08`, `backtrack 0.30` — a very lean
budget, since the 26-value per-column monogram warm start recovers ~100% from ~50 letters, so a few
short restarts and modest climbs suffice), and Seriated Playfair
(`SHAPE_ANNEAL`, **`6x400000` per swept seriation period P**, `inittemp 0.08`, `backtrack 0.30` — Playfair's
own proven square-anneal budget, since each period is a full single-grid anneal with no per-column
decoupling; the 400000-climb cooling schedule is what reliably reaches the optimum, and the budget is
*per P* so the blind sweep multiplies it), and Digrafid
(`SHAPE_ANNEAL`, **`48x150000` per swept period P**, **`inittemp 0.30`**, `backtrack 0.30` — the two-grid
KEYED-ALPHABET search (not the free 54-cell square break) per period; the coarse keyword moves want a warm
temperature and short ciphers want many basins tried, so it is MANY warm restarts, which recovers from ~300
letters — vs the old free-permutation ~700-800 cliff — and the budget is *per P* so the blind sweep multiplies it),
and CM Bifid
(`SHAPE_ANNEAL`, **`8x400000` per swept period P**, `inittemp 0.08`, `backtrack 0.30` — the JOINT two-square
anneal (no decoupling reward) per period; the long per-restart climb is the critical lever (more restarts with
shorter climbs does WORSE near the cliff), so it is a few LONG climbs, which recovers from ~480 letters at an ODD
period — EVEN periods are a documented ciphertext-only degeneracy no budget escapes — and the budget is *per P*
so the blind sweep multiplies it), and the three Interrupted Key bases
(`SHAPE_ANNEAL`, **`6x8000` per kept pre-pass config**, `inittemp 0.08`, `backtrack 0.30` — the climbed
state is just the P per-column key letters, and the Gromark-style pre-pass warm-starts each kept
`(period, strategy, interruptor)` config close to the answer, so a few short restarts x modest climbs
suffice; the ct strategy usually solves on the warm seed. Same lean profile as Progressive Key / Slidefair),
and Condi (`SHAPE_ANNEAL`, **`6x60000` per swept starter (26 configs)**, `inittemp 0.08`, `backtrack 0.30` —
kept deliberately MODEST because the plaintext-feedback cascade makes the true σ an isolated needle with no
basin, so no local-search budget cracks it blind; the entry is a bounded honest attempt, not a tuning
target), and Fractionated Morse (`SHAPE_ANNEAL`, **`16x120000`**, **`inittemp 0.30`**, `backtrack 0.30` — the
single KEYED-alphabet anneal (keyword + tail, `fracmorse_move_seq`, not a free permutation); like Digrafid the
coarse keyword moves want a WARM temperature and RESTARTS are the lever, but there is NO period so the budget
is a single config, not multiplied by a sweep. Tuned against `test_fracmorse_solver`).
This is the mechanism for moving the magic
per-type budgets out of the run scripts and into the binary; add tuned entries for other
types incrementally. The registry is validated end-to-end in `tests/test_playfair_solver.c`.

**Optimisation methods (`-method`, cipher-agnostic).** The engine offers three search
methods, all driven by the *same* `run_solver`/`run_one_config` skeleton over the model
hooks — none of them know anything about the cipher representation:
- **Shotgun hill-climbing** (`SHAPE_SHOTGUN`): greedy uphill + accept-worse with a flat
  `slip_probability`; escape is external (many restarts + backtracking).
- **Simulated annealing** (`SHAPE_ANNEAL`): greedy uphill + Metropolis accept-worse
  `exp((Δ)/temp)` on a geometric `inittemp → mintemp` schedule. Each model declares its
  default shape; `-method shotgun|anneal` overrides it on every type.
- **Particle swarm** (`SHAPE_PSO`, `-method pso`): a memetic, **discrete swap-sequence**
  PSO in `run_one_config_pso` (`engine.c`). It is *only* reachable via `-method pso` (never
  a model's default), and is **completely cipher-agnostic**: a particle's position *is* a
  `SolverState`; "pull toward an attractor (pbest/gbest)" applies the model's own
  `perturb()` moves and keeps only those that reduce a generic Hamming distance over the raw
  state lanes (`state_distance`) — so a permutation stays a permutation, a keyword a keyword,
  a homophone map a map, with no per-cipher code. "Inertia" is a few random `perturb()`
  moves; each particle then does a short greedy local refinement (`-refine`) before its
  decrypt+score updates pbest/gbest. It reuses the budget knobs (`-nhillclimbs` = swarm
  iterations, `-nrestarts` = swarm relaunches) plus `-nparticles`/`-inertia`/`-cognitive`/
  `-social`/`-refine` (defaults in `init_config`; per-type `p_*` registry overrides). Works
  on **every** cipher type (verified end-to-end on polyalpha, square/period, and
  transposition lanes); whether it beats annealing on a given type is a tuning/benchmark
  question, not assumed. Because everything is gated behind `-method pso`, the
  `METHOD_DEFAULT` path — and the whole regression suite — stays byte-for-byte identical.

## How the solver works (mental model)

`solve_cipher()` (in `colossus.c`) is the pipeline:

1. **Period estimation.** For periodic ciphers, `estimate_cycleword_lengths`
   (`perioc.c`) picks candidate cycleword lengths via columnar IoC Z-scores. For
   **autokey** and **transposition-composed** ciphers IoC is useless, so it
   brute-forces lengths `1..max_cycleword_len`.
2. **Shotgun loop.** Nested loops over `(cycleword_len, pt_keyword_len, ct_keyword_len)`.
   Per-cipher-type rules constrain which `(j,k)` keyword-length pairs are valid (e.g.
   Vigenère/Beaufort/Porta force straight alphabets → length 1; Q3/A3 force `j==k`).
   These are the dense `if (...) continue;` blocks.
3. **`shotgun_hill_climber()`.** Random restarts, per-iteration keyword perturbation,
   optional slip (accept-worse to escape local maxima), and backtracking to best.
   Two cycleword strategies:
   - **`-optimalcycle` (default)**: the cycleword is *not* perturbed;
     `derive_optimal_cycleword()` solves each column's key char deterministically by
     maximizing the dot product of decrypted-column letter frequencies vs English
     monograms. Preferred for crib-free attacks.
   - **`-stochasticcycle`**: the cycleword is perturbed randomly like the keyword.
4. **Scoring** (`state_score`): n-gram log-prob is the backbone; with cribs it blends
   in a partial-match `crib_score`. `weight_ioc`/`weight_entropy` default to 0. The
   n-gram table itself has two modes (`load_ngrams`): the default reward-only normalized
   `log(1+count)` (unseen n-grams contribute 0), or, under `-logprob`, AZDecrypt-style
   log10-probabilities with an unseen-n-gram floor penalty (`g_ngram_logprob`).
5. **Reporting**: re-decrypts the best state, applies any transposition, counts
   dictionary words, prints results.

Text is carried internally as **0–25 integer index arrays**, not chars (`ord()` in,
`+ 'A'` out). A "keyword" is a 26-entry keyed-alphabet permutation; a "cycleword" is
the periodic key (sequence of shifts).

## Conventions & gotchas

- **Shared core header + thin per-module headers.** `colossus.h` is the shared
  *core*: the config/ctx/model structs, constants, cipher-type codes, globals, inline
  RNG, and the cipher-*primitive* prototypes (`vigenere_decrypt`, `decrypt_columnar`,
  `playfair_decrypt`, …) — every `.c` includes it. The cipher-agnostic core (`engine`,
  `scoring`, `trans_common`) and each per-cipher-type solver also get a *thin* `.h`
  exposing only that module's public API (`solve_<type>()`, the engine/scoring entry
  points); a `.c` includes the module headers it calls into. Put new solver prototypes
  in the module header, new shared structs/constants/primitive prototypes in `colossus.h`.
  (The already-split primitive files — `vigenere.c`, `beaufort.c`, … — keep their
  prototypes in `colossus.h` rather than carrying their own headers.)
- **`rng_state`** is a global in `utils.c`; the RNG (`fast_rand`, `frand`, `rand_int`,
  `rand_bounded`) is `static inline` in the header, seeded once in `main`. The
  `srand()` call in `main` is dead code — nothing uses libc `rand()`.
- **Stack-heavy.** `solve_cipher` and the hill climber declare several
  `MAX_CIPHER_LENGTH` (10000) int arrays on the stack, and there is **no bounds check**
  after `fscanf("%s", ...)`, so inputs must stay under the limit.
- **`-variant`** swaps decryption for encryption in the Quagmire/Vigenère math
  (reciprocal tableau). **`-samekey`** ties keyword and cycleword together.
- **Transposition is a post-decrypt stage**: `-transperoffset <offset> <period>` or
  `-transmatrix <w1> <w2> <cw|ccw>`. Crib positions are un-mapped back through it via
  `map_crib_to_cipher_pos` so cribs still line up.

## Fixed issues (end-to-end tests in `ciphers/tests/bugfixes/`)

- The trailing partial-crib match line (the `_`/digit/`*` row) used to index the
  packed `crib_indices` array positionally (testing an uninitialized `-1`), printing
  garbage at every position. Now indexed by cipher position via `cribtext_str`.
  Test: `bug1_partial_crib.sh`.
- The `-transmatrix` `>>>` summary (no-dictionary branch) printed period/offset
  instead of `w1`/`w2`/`clockwise`. Fixed to match the dictionary branch.
  Test: `bug2_transmatrix_summary.sh`.
- `load_ngrams` looped on `while(!feof(fp))`, re-reading the final line and
  mis-assigning a stale `freq` on any trailing/malformed line. Now loops on
  `fscanf(...) == 2`. Test: `bug3_ngram_load.sh`.
- A cipher/ngram/crib path longer than `MAX_FILENAME_LEN` overflowed the fixed
  `char[]` in `ColossusConfig` (`main()` `strcpy`s the CLI arg in unbounded),
  corrupting the struct and crashing with SIGILL — any absolute path past the old
  100-byte limit triggered it. `MAX_FILENAME_LEN` raised to 4096.
  Test: `bug4_long_path.sh`.
- `int_pow` did a final `base *= base` after the result was already accumulated;
  `int_pow(26, ngram_size)` overflowed signed `int` (e.g. `26^4` squared). The UB
  was benign at `-O0` but `-O3` could exploit it. Now skips the unused final squaring.

## Working agreements

- Match the existing style: 4-space indent, `snake_case`, integer-index text arrays,
  explicit per-cipher-type `switch`/`if` ladders. The code favors explicitness over
  abstraction — don't refactor the cipher-type dispatch into clever generic code.
- The binary `colossus`, `*.o`, and `.DS_Store` are git-ignored.
- Don't commit or push unless asked.
