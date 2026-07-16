# CLAUDE.md

Guidance for working in this repository.

> This file was condensed from a much longer version. The full per-cipher design
> rationale now lives where it belongs: in each solver module's header comment, in the
> unit tests, and in the auto-memory (`memory/MEMORY.md` indexes a note for most notable
> ciphers). Reach for those when you need the deep story on one type.

## Scope

This directory is the **entire project** and the git root. It tracks
`https://github.com/stblake/colossus` (branch `main`). Everything outside this
directory is out of scope — git can't see it, and neither should you. (The parent
folder holds unrelated experiment runs, logs, and candidate dumps; ignore it.)

## What this is

Colossus is a polyalphabetic substitution cipher solver in C by Sam Blake (started
14 July 2023). It began as a Vigenère-family solver (Vigenère, Gronsfeld, Beaufort,
Porta, Quagmire I–IV, Autokey, Progressive Key, Interrupted Key, Condi) optionally
composed with a transposition stage, and has grown to cover most ACA cipher types
(polygraphic squares, fractionation, Morse, transposition variants — see the type
reference below). The core engine is a **stochastic, slippery, shotgun-restarted hill
climber with backtracking** (plus annealing and PSO). Cipher conventions follow the
American Cryptogram Association (https://www.cryptogram.org/resource-area/cipher-types/).
It exists to crack the Kryptos sculpture's K1–K4. See `README.md` for the author's writeup.

## Layout

Sources are grouped by cipher class under `src/<class>/`; everything else (build,
tests, tools, data, ciphers) is at the repo root. All local `#include`s are **flat**
(`#include "foo.h"`, no dir prefix); the makefile's `INCLUDES` var supplies one `-I`
per `src/` subdir, so a header is found regardless of which subdir it's in. Add a new
`src/` subdir → add it to `INCLUDES`.

```
src/core/         # cipher-agnostic engine + shared infrastructure
  colossus.c        # main(): arg parsing, init_config(), solve_cipher() dispatcher
  colossus.h        # shared CORE header: config/ctx/model structs, constants, cipher-type
                    #   codes, globals, inline RNG, cipher-PRIMITIVE prototypes
  engine.c/.h       # search engine: run_solver(), run_one_config(), make_solver_ctx(),
                    #   SearchDefaults registry + apply_cipher_defaults(); anneal/shotgun/pso
  scoring.c/.h      # state_score / ngram_score / crib_score, load_ngrams, keyword RNG
  parse.c           # parse_cipher_type(): string/int aliases -> cipher-type code
  perioc.c          # estimate_cycleword_lengths(): IoC period estimation
  optimal_cycleword.c  # derive_optimal_cycleword(): deterministic per-column frequency attack
  dict.c            # dictionary load + word-finding
  utils.c           # ord/print, decode_cipher/print_cipher (symbol I/O), IoC, chi-squared

src/polyalphabetic/  # Vigenère family (mostly inside POLYALPHA_MODEL)
  polyalpha_solver.c/.h   # POLYALPHA_MODEL (vig/quag/beau/porta/autokey/gronsfeld) + solve_polyalpha()
  vigenere.c gronsfeld.c beaufort.c porta.c quagmire.c autokey.c   # per-cipher primitives
  gromark.c gromark_solver.c/.h        # Gromark + Periodic Gromark (own CipherModels)
  nicodemus.c nicodemus_solver.c/.h    # Nicodemus (substitution + per-block columnar)
  progkey.c progkey_solver.c/.h        # Progressive Key (3 base types)
  intkey.c intkey_solver.c/.h          # Interrupted Key (3 base types)
  condi.c condi_solver.c/.h            # Condi (plaintext-feedback substitution)

src/transposition/   # pure-transposition solvers + shared helpers
  trans_common.c/.h    # report_transposition(), TransKeyOps, perm helpers, held_karp_best_path(),
                       #   seam_best_row_order() (exact best within-column track order L), trans_word_set()
  transpositions.c     # transperoffset() (decimation), transmatrix() (K3-style double rotation),
                       #   decrypt_columnar(), decrypt_tile(), period_column_transform()
  transmatrix_solver / permutation_solver / columnar_solver
  railfence / route / amsco / myszkowski / redefence / cadenus / nihilist / swagman / grille _solver
  columnar_track_solver.c/.h        # transcol-L: columnar + within-column row perm L; -cribanchored matcher
  route_chain_solver.c/.h           # transroutecol: fixed read-route + searched column key
  tile_solver.c/.h                  # transtile: sub-grid h×w tile transposition
  period_column_solver.c/.h         # period-column: DETERMINISTIC EXHAUSTIVE depth<=2
  period_column_space_solver.c/.h   # period-column-space: space-robust indel-repair variant
  double_transposition_solver.c/.h  # transcol2-dc: double columnar, divide & conquer (IDP, Lasry 2014)

src/polygraphic/     # square/cube/matrix ciphers — each: primitive + a CipherModel solver
  playfair / bifid / trifid / hill / phillips / twosquare / foursquare / adfgvx
  nihilist_sub / bazeries / portax / slidefair / seriated_playfair / digrafid
  cm_bifid / trisquare / fracmorse / pollux / morbit / straddling_checkerboard

src/substitution/    # monoalphabetic / homophonic
  indep_solver / homophonic_solver / ragbaby (+ _solver)

makefile  README.md  LICENSE  example.sh
cipher.txt  crib.txt              # sample ciphertext + crib
tools/<type>_gen.c                # standalone per-type test-data generators (make <type>_gen)
english_quadgrams.txt             # default n-gram table; english_quintgrams.txt optional (with -logprob)
OxfordEnglishWords.txt            # default dictionary (auto-loaded if present in cwd)
ciphers/kryptos/                  # K1–K4 ciphertexts + run scripts
ciphers/tests/                    # per-cipher end-to-end cases + run_tests.sh
```

## Build

```bash
make            # builds ./colossus
make clean
```

- Build with **Homebrew gcc-16**, not Apple clang's `gcc` shim (see the build memory).
- The `CC` line does **not** include `-lm` — links on macOS (clang folds libm into
  libc) but **needs `-lm` on Linux**.
- `make` also runs `cp colossus ..` (and `../quagmire`), copying the binary outside
  this dir — predates the repo isolation; the in-tree `./colossus` is what matters.
- TU lists live in makefile vars, each prefixed with a class dir var (`$(CORE)`,
  `$(POLY)`, `$(TRANS)`, `$(GRAPH)`, `$(SUBST)`): `PRIMITIVES`, `SOLVERS`, and
  `SOLVER_SRC = $(PRIMITIVES) $(SOLVERS) $(CORE)/colossus.c`. Add a new solver module
  to `SOLVERS` with its `src/<class>/` prefix.

## Test

- **`make test`** — framework-free unit tests of the cipher **primitives**
  (`tests/test_<type>.c`): each pins an ACA/Wikipedia worked-example known-answer vector
  and does encrypt/decrypt round-trips (+ structural invariants) over random keys ×
  lengths × periods incl. ragged/edge cases.
- **`make testopt`** — additionally runs the in-process **solver** regressions
  (`tests/test_<type>_solver.c`): validate the `SearchDefaults` registry entry, measure
  a capability floor + length cliff, sweep keywords/periods (blind-selection asserted
  where applicable), and calibrate `-method anneal/shotgun/pso`. These pin each solver's
  real capability and tuning.
- **`ciphers/tests/run_tests.sh`** — the **accuracy regression suite**: a manifest of
  end-to-end cases (one per cipher family) each solved to ~100% with a **fixed `-seed`**
  + quadgrams. It runs the solver, pulls the recovered plaintext from the `>>>` CSV
  line, char-compares to a sibling `<name>.solution`, prints per-test accuracy + time,
  and exits non-zero if any drops below threshold (default 99%). Fixed seed ⇒ a
  bit-identical refactor keeps every score at 100%. Tiers: `--fast` (~64s, iterate with
  this), `--slow` (heavier square/fractionation solves), no flag runs both. Add a case:
  append a `tier|name|type|cipher|args` line, then `./run_tests.sh --generate <name>`
  once the recovered text is verified.

**When you add or change a solver:** update its module header comment, add/extend its
`test_<type>.c` + `test_<type>_solver.c`, add a `run_tests.sh` case if it can reach
~99%, and record the non-obvious findings in a memory note.

## Run

Run from this directory (the binary loads its n-gram table, dictionary, and ciphertext
from cwd).

```bash
./example.sh
# minimally:
./colossus -type q3 -cipher cipher.txt -ngramsize 4 -ngramfile english_quadgrams.txt
```

Required: `-type`, a source (`-cipher <file>` or `-batch <file>`), `-ngramsize`,
`-ngramfile`. Everything else defaults (see `init_config`). Output is a human block then
a `>>> ...` one-line CSV summary for batch grep/sort. `-type` accepts an alias or integer
code (full list in `parse.c`, codes in `colossus.h`). By default only the **first line**
of `-cipher` is read; `-multiline` reads the whole file (dropping newlines) so a
multi-line ciphertext is concatenated.

### Type reference (code — aliases — essence)

Polyalphabetic (in `POLYALPHA_MODEL`, share the IoC/optimal-cycleword pipeline):
- `0` vig · `1..4` q1–q4 (Quagmire) · `5` beau · `6` porta · `7` auto ·
  `8..11` auto1–4 · autobeau · autoporta · `34` gronsfeld/gron (Vigenère, digits 0–9).

Own CipherModels, polyalpha-adjacent:
- `49` gromark/gm, `50` gromark-periodic/pgromark — keyed alphabet + chain-addition running key.
- `51/52/53` nicodemus[-variant/-beaufort] — substitution + per-block columnar; sweeps (P, H).
- `56/57/58` progkey[-var/-beau] — periodic key + per-group drift; period × progression enumerated.
- `66/67/68` intkey[-var/-beau] — periodic keyword reset at break points; period swept + strategy enumerated.
- `69` condi — plaintext-feedback substitution over keyed σ.

Transposition (isolated by an early branch in `solve_cipher`, optimization not keyword-search):
- `14` transmatrix · `15` transperoffset · `16` transposition · `17` transcol · `18` transcol2.
- `41` transcol-l/coltrack · `42` transroutecol/routecol · `43` transtile/tile.
- `71` period-column/pcol · `72` period-column-space/pcolsp · `73` transcol2-dc/dcol.
- Plus railfence/route/amsco/myszkowski/redefence/cadenus/nihilist/swagman/grille solvers.

Polygraphic squares/cubes/matrix:
- `30` playfair/pf · `31` bifid/bf · `32` trifid/tf · `33` hill.
- `35/36/37` phillips[-c/-rc] · `38/39` twosquare[-v]/ts[v] · `40` foursquare/fs.
- `44` adfgx · `45` adfgvx/adfg · `46/47/48` nihilist-sub[-nc/-m100]/nihsub.
- `54` bazeries/baz · `55` portax/ptx · `59/60/61` slidefair[-var/-beau]/sf.
- `62` seriated-playfair/spf · `63` digrafid/df · `64` cm-bifid/cmb · `65` trisquare/3sq.
- `82` checkerboard/checker/cb (keyed 5x5 square, 25-letter J->I; plaintext letter -> (row label,
  col label) digraph). Case auto-detected PER AXIS from the ciphertext (an axis with >5 distinct
  labels is complex). Label ORDER is not identifiable (absorbed by a row/col permutation of the
  square, like nihilist-sub). SIMPLE (1 label/axis) ⇒ a free 25-code→25-letter bijection = an
  Aristocrat over the merged codes, on the homophonic incremental fast path. COMPLEX (2 labels/axis,
  homophonic) ⇒ a square-INDEPENDENT per-axis χ² homogeneity pre-pass ranks the label PAIRINGS
  (top-K per axis crossed into engine configs). Simple recovers ~100% from ~130 letters; complex
  needs ~400-600+ (below the ACA 60-90 range — see notable findings). -logprob.

Morse / checkerboard (digit-stream input parsed from `ciphertext_str`):
- `70` fractionated-morse/fm · `74` pollux/pol · `75` morbit/mor · `76` straddling/sc ·
  `78` monome-dinome/md (3x8 box, 24-letter J->I/Z->Y; needs quintgrams + dict — config
  selection is by dictionary coverage, n-gram alone is gamed cross-config).
- `81` tridigital/td (keyed 3x10 block, full 26-letter alphabet, one digit per letter + a
  word-separator digit). AMBIGUOUS 3-to-1 decode (unique in this family): the key is a
  partition of 26 letters into 9 column-groups, the plaintext is chosen per position by an
  inner beam-Viterbi (spaces transparent → context carries across words); separator picked by
  word-length fit, cross-config winner by WHOLE-WORD coverage. Dense polyphonic ⇒ partial,
  high-variance recovery below the 99% floor at all practical lengths; needs quintgrams + dict.

Substitution:
- `28` indep · `29` homophonic · `77` ragbaby/rag · `79` aristocrat/arist · `80` patristocrat/patri
  (one solver core: free 26-perm climbed by n-gram with the homophonic incremental fast path;
  word divisions preserved for the Aristocrat's spaced report, dropped/5-grouped for the
  Patristocrat; -logprob).

### Key global flags

- `-logprob` (a.k.a. `-azdecrypt`): AZDecrypt-style log10 n-gram fitness with an
  unseen-n-gram floor penalty, vs the default reward-only `log(1+count)` (unseen → 0).
  **Effectively required** for the square/fractionation types; pairs well with
  quintgrams (`-ngramsize 5 -ngramfile english_quintgrams.txt`). Default off ⇒ unchanged.
- `-reversengrams` (`-revngrams`): symmetrize the table so each n-gram and its reversed
  twin share the `max` weight — reads reversed-word text like clean English (for the W168
  alternate-word-reversal hypothesis). Roughly doubles the acceptable solution set.
  Default off ⇒ bit-identical.
- `-cribdrag WORD` (or `WORDA|WORDB`): position-free crib. Each word is slid across the
  decrypt and its best-offset partial match rewarded (`-weightcribdrag`, default 36);
  pipe = AND (mean over words). A global toggle in `state_score`, so it works for
  **every** solver. Steers score-driven keys most under `-stochasticcycle`. Default off ⇒
  bit-identical.
- `-crib`: fixed crib pinned to absolute cipher positions (blends in `state_score`).
- `-method shotgun|anneal|pso`: override the model's default search shape (below).
- `-nthreads N` (default 1): parallelize the restart loop across N pthreads (splits
  `-nrestarts`, ~N× faster). N=1 is the original sequential path, **bit-identical**; N>1
  deterministic per (seed, N) modulo tie-break order. Forced to N=1 for:
  deterministic-exhaustive solvers, standalone transposition climbers, the homophonic
  incremental fast-path.
- `-optimalcycle` (default) / `-stochasticcycle`: derive the cycleword by column
  monograms vs perturb it randomly.
- `-variant`: swap decrypt↔encrypt in the Quagmire/Vigenère math (reciprocal tableau).
  `-samekey`: tie keyword and cycleword together.
- `-multiline`, `-delimiter <char>` (tokenized symbol I/O — see below).
- Post-decrypt transposition stage: `-transperoffset <offset> <period>` /
  `-transmatrix <w1> <w2> <cw|ccw>` (distinct from the `-type` transposition solvers;
  cribs are un-mapped back through it via `map_crib_to_cipher_pos`).
- Sweep/estimator/search knobs: `-period`, `-cyclewordlen`, `-mincols`/`-maxcols`,
  `-maxperiod`, `-nperiods`, `-blockheight`/`-maxblockheight`, `-depth` (period-column),
  `-readdir tb|bt|both`, `-readrowdir`, `-nprimers`, `-nrestarts`/`-nhillclimbs`,
  `-inittemp`/`-mintemp`, `-nparticles`/`-inertia`/`-cognitive`/`-social`/`-refine` (PSO).
- Cipher-specific: `-progression`, `-intscheme ct|pt|breaks|joint`, `-breaks <file>`,
  `-interruptor <A-Z>`, `-startkey`, `-tile h w`, `-maxgaps`/`-maxdels`, `-cribanchored`,
  `-weightword`, `-weightmono`, `-weightstructure`.

**Tokenized symbol I/O.** `decode_cipher()` (utils.c) decodes ciphertext. For every type
except homophonic-with-no-`-delimiter` this is byte-identical to the historical per-char /
0..25 encoding (regression suite stays bit-identical). For `homophonic` (or any type with
`-delimiter <char>`) it tokenizes into a `SymbolTable` and emits one symbol id per
position, so a ciphertext alphabet larger than A..Z works. Default delimiter: auto (comma
if the homophonic input has one, else per-char).

## Cross-cutting design patterns

The recurring ideas behind the per-cipher solvers — apply these when adding a new type:

- **Optimization-only engine.** The whole key is optimized; we don't add exhaustive
  drivers even when the keyspace is a number/string (Bazeries climbs N's digits). The
  exception is genuine **needles** (below).
- **Deterministic-exhaustive for needles.** When one key change re-parses the whole
  decrypt (no gradient/basin) AND the keyspace is small, enumerate instead of climbing:
  Period column (depth ≤ 2), Pollux (3¹⁰), Morbit (9!). A stochastic climber flails where
  enumeration is certain. These take N=1 regardless of `-nthreads`.
- **Decoupling rewards.** For a coupled (square + key) search, find a statistic that
  depends on only ONE half and fold it into `score_adjust` to give that half a gradient
  flat in the other: ADFGVX's structural IoC (column order), Nihilist-Sub's validity
  (additive key), Bazeries' monogram fit (the square), the Morse types' validity reward.
- **Per-column monogram warm start.** When each column/position is enciphered by one key
  letter independently (Portax, Slidefair, Progressive Key after de-progressing,
  Interrupted-Key ct, Nicodemus after de-transposing), derive each column's shift by
  monogram fit to warm-start the seed; the n-gram anneal only corrects a few columns.
  Descends from `derive_optimal_cycleword`. These recover from short text, no `-logprob`.
- **Keyed-alphabet search, not free permutation.** ACA keyed alphabets are keyword +
  ascending tail. Search that structure (`*_move_seq` family: fracmorse/digrafid/ragbaby)
  rather than a free 26!/54-cell permutation — it tracks the keyword and drops the blind
  cliff dramatically (Digrafid ~700 → ~300 letters).
- **Joint multi-square anneal.** When no decoupling reward exists (both/all squares in the
  n-gram fitness), anneal all squares packed back-to-back, perturbing one per move:
  Four-Square, CM-Bifid, Tri-Square, Straddling.
- **Period: sweep vs estimate.** IoC period estimation works for stationary periodic
  ciphers (Bifid/Trifid columnar-IoC top-K annealed). It **fails** through digraphic
  pairing, transposition, or key drift — those **sweep** P (one engine config per P) and
  let the n-gram score pick (a wrong P decrypts to gibberish).
- **Length change** (fractionation/Morse: N pt ↔ C ct): pass the plaintext/scoring length
  to `make_solver_ctx` and either tile the variable-length decode to a fixed length
  (fracmorse) or let the mean-n-gram be length-fair (Pollux/Morbit).

## Notable per-type findings & limitations

Documented structural facts (asserted or characterized in the solver tests), not solver bugs:

- **Condi** — the plaintext feedback makes the true σ an **isolated needle** (one swap
  cascades the whole downstream decrypt): no local search cracks it blind at any budget.
  The untapped tractable attack is crib-anchored constraint solving of σ.
- **CM Bifid** — **even periods are degenerate ciphertext-only** (rows/cols never share an
  output pair → transpose-like square ambiguity, no budget escapes); odd periods recover
  from ~480 letters.
- **Straddling Checkerboard** — letters recover ~100% from ~100–150 chars, but
  **numeric/figure-shift is a documented limitation**. Solve the FREE code→cell bijection,
  not arrangement+labels (redundant, stalls). No cheap statistic ranks the 45 tokenization
  configs → an SA mini-solve pre-pass, keep top 12, warm-start each.
- **Playfair / Seriated Playfair** — rare-letter-X ambiguity pins some grids at ~92%.
  Square grids recover only up to a cyclic row/col rotation (plaintext is unique).
- **Nicodemus / Slidefair / Progressive Key** — Vigenère and Variant are not separately
  identifiable (a free derived shift absorbs the sign); only Beaufort is distinct.
- **Interrupted Key** — ct-interruptor is the reliable blind workhorse (decouples like
  Vigenère); pt-interruptor is fragile (causal reset, rugged basins); breaks/joint for
  random breaks.
- **Tri-Square** — easier than Four-Square despite 75 cells: the polyphonic c0/c2 letters
  must randomize on encode (a canonical choice starves the gradient).
- **Checkerboard** — the SIMPLE case reduces exactly to an Aristocrat over 25 merged codes
  (recovers ~100% from ~130 letters). The **COMPLEX case sits BELOW the ACA 60–90 range**: the
  per-axis pairing statistic has an O(N) bias favouring wrong pairings against an O(N²) signal, so
  the true pairing does not rank first until ~400–600 plaintext letters (the calibrated rank curve
  in `test_checkerboard_solver.c` — r323 at N=90 → r0 at N=600). Documented limitation, like
  Tridigital / CM-Bifid even periods — not a solver bug. The label KEYWORDS (BLACK/WHITE/…) are
  unrecoverable ciphertext-only (label order folds into the square). The ACA square is
  spiral-routed (keyword+tail read clockwise); the route matters only to the generator/tests, not
  the solver (which searches the composite code→letter map).

## SearchDefaults (per-type schedules)

`init_config()` globals suit the polyalphabetic/transposition score scale. A type whose
score lives on a different scale gets a tuned profile in the compiled-in registry
(`g_search_defaults[]` in `colossus.c`, keyed by cipher type), carrying anneal (`a_*`),
shotgun (`s_*`), and PSO (`p_*`) knobs. `main()` overlays the matching profile before the
arg loop, so precedence is **globals < registry < explicit CLI flags**. Types with no
entry keep the global defaults bit-for-bit (regression suite unaffected). This moves magic
per-type budgets out of run scripts into the binary; add entries incrementally. Validated
in `tests/test_playfair_solver.c`. (For the exact budget of a given type, read the
registry — don't hardcode it here.)

## Optimisation methods (`-method`, cipher-agnostic)

All three run over the *same* `run_solver`/`run_one_config` skeleton via the model hooks —
none know the cipher representation:
- **Shotgun** (`SHAPE_SHOTGUN`): greedy uphill + flat `slip_probability` accept-worse;
  escape via restarts + backtracking.
- **Anneal** (`SHAPE_ANNEAL`): greedy uphill + Metropolis `exp(Δ/temp)` on a geometric
  `inittemp → mintemp` schedule. Each model declares its default shape.
- **PSO** (`SHAPE_PSO`, only via `-method pso`): memetic discrete swap-sequence swarm in
  `run_one_config_pso`. A particle's position *is* a `SolverState`; "pull toward
  pbest/gbest" applies the model's own `perturb()` and keeps moves that reduce a generic
  Hamming distance (`state_distance`) — so a permutation stays a permutation, etc., with no
  per-cipher code. Works on every type; whether it beats annealing is a tuning question.
  Gated behind `-method pso`, so `METHOD_DEFAULT` and the regression suite stay byte-identical.

**Thread-safety** (for `-nthreads`): every hook-reachable static scratch/lazy cache written
during the search is `_Thread_local` (e.g. `bifid.c`'s `g_bifid_stream`, `phillips.c`'s
derived-square scratch, running-key chains, `scoring.c` scratch); setup-phase statics
written once on the main thread and only read in the search stay shared. Verified race-free
under ThreadSanitizer.

## How the solver works (mental model)

`solve_cipher()` (in `colossus.c`) dispatches: transposition / deterministic-exhaustive /
Morse-digit / space-significant types branch out early; the rest run the periodic pipeline:

1. **Period estimation** — `estimate_cycleword_lengths` (`perioc.c`) picks candidate
   lengths by columnar IoC Z-scores. For autokey / transposition-composed ciphers IoC is
   useless, so lengths `1..max_cycleword_len` are brute-forced.
2. **Shotgun loop** — nested loops over `(cycleword_len, pt_keyword_len, ct_keyword_len)`
   with per-type validity constraints (the dense `if (...) continue;` blocks; e.g.
   Vigenère/Beaufort/Porta force straight alphabets → length 1; Q3/A3 force `j==k`).
3. **`shotgun_hill_climber()`** — random restarts, per-iteration keyword perturbation,
   optional slip, backtracking. Cycleword strategy: `-optimalcycle` (default; derive each
   column's key deterministically by monogram fit) or `-stochasticcycle` (perturb it too).
4. **Scoring** (`state_score`) — n-gram log-prob backbone + optional `crib_score`;
   `weight_ioc`/`weight_entropy` default 0. Table modes via `load_ngrams` (reward-only
   default vs `-logprob`).
5. **Reporting** — re-decrypt best state, apply any transposition, count dictionary words, print.

Text is carried internally as **0–25 integer index arrays** (`ord()` in, `+ 'A'` out).
A "keyword" is a 26-entry keyed-alphabet permutation; a "cycleword" is the periodic key.
`MAX_ALPHABET_SIZE` is 36 (largest runtime alphabet: ADFGVX's A..Z+0..9); `ALPHABET_SIZE`
(26) stays the hardcoded mod base of the polyalphabetic primitives. Some types force a
different runtime alphabet before `load_ngrams` via an `init_alphabet*` call: Trifid 27
(A..Z+`+`), Digrafid 27 (A..Z+`#`), Playfair/Bifid/etc. 25 (J→I), Ragbaby 24 (I/J, W/X
paired).

## Conventions & gotchas

- **Header split.** `colossus.h` is the shared core (config/ctx/model structs, constants,
  cipher-type codes, globals, inline RNG, primitive prototypes) — every `.c` includes it.
  The cipher-agnostic core and each per-type solver also get a thin `.h` exposing only
  their public API. New solver prototypes → the module header; new shared
  structs/constants/primitive prototypes → `colossus.h`. Already-split primitive files
  (`vigenere.c`, …) keep their prototypes in `colossus.h`.
- **`rng_state`** is `_Thread_local` in `utils.c`; the RNG (`fast_rand`, `frand`,
  `rand_int`, `rand_bounded`) is `static inline` in the header, seeded once in `main`
  (per-worker via `rng_seed_thread`). The `srand()` in `main` is dead code.
- **Stack-heavy.** `solve_cipher` / the hill climber declare several `MAX_CIPHER_LENGTH`
  (10000) int arrays on the stack with **no bounds check** after `fscanf("%s", ...)` —
  inputs must stay under the limit.
- **Space-significant types** (period-column, period-column-space, transcol2-dc, ragbaby,
  and the digit-stream Morse/checkerboard types) carry spaces/punctuation as real grid
  cells or parse `ciphertext_str` directly; `main()` does not trim trailing whitespace for
  the pure-transposition types (a trailing space is a real cell).
- **Test cipher files** (`ciphers/tests/*.txt`) have trailing metadata after the cipher —
  readers must stop at the first newline.

## Fixed issues (regression tests in `ciphers/tests/bugfixes/`)

- Partial-crib line indexed the packed `crib_indices` array positionally (garbage); now by
  cipher position via `cribtext_str`. (`bug1_partial_crib.sh`)
- `-transmatrix` `>>>` summary (no-dict branch) printed period/offset instead of
  w1/w2/clockwise. (`bug2_transmatrix_summary.sh`)
- `load_ngrams` looped on `while(!feof(fp))`, re-reading/mis-assigning the last line; now
  loops on `fscanf(...) == 2`. (`bug3_ngram_load.sh`)
- A path longer than `MAX_FILENAME_LEN` overflowed the fixed `char[]` in `ColossusConfig`
  (unbounded `strcpy`), SIGILL; limit raised to 4096. (`bug4_long_path.sh`)
- `int_pow` did a final `base *= base` after accumulating; `int_pow(26,4)` overflowed
  signed int (benign at -O0, exploitable at -O3). Now skips the unused final squaring.

## Working agreements

- Match the existing style: 4-space indent, `snake_case`, integer-index text arrays,
  explicit per-cipher-type `switch`/`if` ladders. The code favors explicitness over
  abstraction — don't refactor the cipher-type dispatch into clever generic code.
- New optimisation methods go in the **engine** (work for all types via the existing
  hooks), never as per-cipher hooks.
- Refactors must stay bit-identical at fixed seed: `run_tests.sh --fast` should keep every
  score at 100%.
- The binary `colossus`, `*.o`, and `.DS_Store` are git-ignored.
- Don't commit or push unless asked.

<!-- code-review-graph MCP tools -->
## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

### When to use graph tools FIRST

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with callers_of/callees_of/imports_of/tests_for
- **Architecture questions**: `get_architecture_overview` + `list_communities`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.

### Key Tools

| Tool | Use when |
| ------ | ---------- |
| `detect_changes` | Reviewing code changes — gives risk-scored analysis |
| `get_review_context` | Need source snippets for review — token-efficient |
| `get_impact_radius` | Understanding blast radius of a change |
| `get_affected_flows` | Finding which execution paths are impacted |
| `query_graph` | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes` | Finding functions/classes by name or keyword |
| `get_architecture_overview` | Understanding high-level codebase structure |
| `refactor_tool` | Planning renames, finding dead code |

### Workflow

1. **The graph is built/refreshed at the start of every session**, automatically: the
   `SessionStart` hook in `.claude/settings.json` runs `code-review-graph update` (falling
   back to a full `build` if no graph exists yet, e.g. a fresh clone) and then prints
   `status` as context. It takes ~1s warm. You do **not** need to build it by hand — but if
   the graph ever looks stale or the status line is missing, run `code-review-graph update
   --repo .` (or `build --repo .` to re-parse everything from scratch).
2. It also auto-updates after each Edit/Write/Bash (the `PostToolUse` hook, `--skip-flows`).
3. Use `detect_changes` for code review.
4. Use `get_affected_flows` to understand impact.
5. Use `query_graph` pattern="tests_for" to check coverage.

### Scope in this repo

- Covers **C and bash**: ~308 files / ~1900 nodes / ~14k edges, spanning `src/`, `tests/`,
  and `tools/`. Call edges resolve across the flat-include layout (e.g. `ngram_score` in
  `src/core/scoring.c` correctly lists its callers in the per-cipher solvers).
- The graph DB lives in `.code-review-graph/` and is git-ignored.
- The MCP server is declared in `.mcp.json` (`uvx code-review-graph serve`). The
  `mcp__code-review-graph__*` tools only appear **after a Claude Code restart** following
  that file being added.
