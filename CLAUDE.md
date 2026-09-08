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
  cm_bifid / twin_bifid / twin_trifid / trisquare / fracmorse / pollux / morbit / straddling_checkerboard

src/substitution/    # monoalphabetic / homophonic
  indep_solver / homophonic_solver / ragbaby (+ _solver)

src/machine/         # rotor machines
  enigma.c/.h            # the Enigma machine primitive (rotors I-VIII, Greek Beta/Gamma,
                         #   reflectors B/C + thin; stepping w/ double-step; self-reciprocal)
  enigma_solver.c/.h     # ENIGMA_MODEL + solve_enigma(): ciphertext-only IoC/ring/plugboard
                         #   attack (Gillogly), known-key decrypt, dispatch to the Bombe
  enigma_bombe.c         # Turing-Welchman menu + diagonal-board crib attack (solve_enigma_bombe)

makefile  README.md  LICENSE  example.sh
cipher.txt  crib.txt              # sample ciphertext + crib
tools/<type>_gen.c                # standalone per-type test-data generators (make <type>_gen)
tools/build_corpus_ngrams.py      # build n-gram tables (orders 1-8) from raw corpus text (English format)
tools/compress_ngrams.py          # text n-gram table -> dense 8-bit .ngbin (stdlib-only): the 26-letter
                                  #   letter tables, or `--spaces` for the 27-symbol {A..Z,' '} tables; for
                                  #   reduced letter alphabets use `-writengrambin` — see the .ngbin note below
tools/build_corpus_spaces_ngrams.py # SPACE-INCLUSIVE sibling: keeps word spaces (window over {A..Z,' '}),
                                  #   builds <lang>_spaces_*grams.txt for the -spaces readability pass
tools/gutenberg_fetch.py          # enumerate + download Project Gutenberg books by language (gutendex)
tools/{english,french,german,italian,spanish,latin,danish,portuguese}_corpus/  # each: list.tsv + MANIFEST (Gutenberg books; books/ git-ignored, re-downloadable)
ngram_data/english/               # mono/bi/tri/six = corpus (798 books); quad+quint = EXTERNAL web-scale
                                  #   (quadgrams = tracked default; quad/quint kept external for test calibration)
ngram_data/french/                # french_{mono,bi,tri,quad,quint,six}grams.txt (built from ~448M letters, 976 books)
ngram_data/german/                # german_{mono,bi,tri,quad,quint,six}grams.txt (built from ~650M letters, 2365 books)
ngram_data/italian/               # italian_{mono,bi,tri,quad,quint,six}grams.txt (built from ~344M letters, 994 books)
ngram_data/spanish/               # spanish_{mono,bi,tri,quad,quint,six}grams.txt (built from ~283M letters, 883 books)
ngram_data/latin/                 # latin_{mono,bi,tri,quad,quint,six}grams.txt (built from ~24M letters, 88 books;
                                  #   small Gutenberg holding + a 2-stage English-apparatus filter — see MANIFEST)
ngram_data/danish/                # danish_{mono,bi,tri,quad,quint,six}grams.txt (built from ~21M letters, 83 books;
                                  #   fold() gained ø->o, since ø has no NFKD decomposition — see MANIFEST)
ngram_data/portuguese/            # portuguese_{mono,bi,tri,quad,quint,six}grams.txt (built from ~95M letters, 642 books;
                                  #   all diacritics á/ã/ç/õ/… fold via NFKD, NO code change — see MANIFEST)
ngram_data/dutch/                 # dutch_{mono,bi,tri,quad,quint,six}grams.txt (built from ~296M letters, 1096 books;
                                  #   all diacritics + the IJ ligature ĳ/Ĳ fold via NFKD, NO code change; 2 EN-NL
                                  #   dictionaries dropped by content not density — see MANIFEST)
                                  # Each lang dir ALSO holds <lang>_spaces_{bi,tri,quad,quint,six}grams.txt:
                                  #   SPACE-INCLUSIVE char n-grams over {A..Z,' '} (window may embed a literal
                                  #   space; count split on the RIGHTMOST space) feeding -spaces (src/core/spaces.c).
                                  #   All corpus-derived orders 2-6 (English spaces tables are corpus, unlike its
                                  #   external letter-only quad/quint). Built by tools/build_corpus_spaces_ngrams.py.
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
./colossus -type q3 -cipher cipher.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt
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
- `90` running-key/runningkey/rk — Vigenère-family with a running-TEXT key (key length ==
  message length, no period). ACA self-keyed by default (the plaintext's first half is its own
  running key, so the whole 2N passage is recovered) — `-indepkey` for an independent-key text,
  `-runningkeyfile <file>` for a KNOWN key (deterministic decrypt; the drag-K1/K2/K3 capability).
  Four families (Vigenère/Beaufort/Variant/Porta) swept, n-gram picks (`-variant`/`-beaufort`
  pin). Blind scores BOTH streams as English (beam warm start + anneal over the key stream);
  cribs are hard-ANCHORED (a known plaintext letter forces the key letter). BLIND recovery is a
  documented limitation (running key is gamed without a key or crib — see notable findings);
  known-key is exact and crib-anchoring recovers a bounded free gap. -logprob (+ quintgrams).

Transposition (isolated by an early branch in `solve_cipher`, optimization not keyword-search):
- `14` transmatrix · `15` transperoffset · `16` transposition · `17` transcol · `18` transcol2.
- `41` transcol-l/coltrack · `42` transroutecol/routecol · `43` transtile/tile.
- `71` period-column/pcol · `72` period-column-space/pcolsp · `73` transcol2-dc/dcol.
- `83` sequence-transposition/seqtrans/st (ACA): a Gromark chain-addition digit sequence (5-digit
  primer) buckets each plaintext letter into 1 of 10 columns; a 10-letter keyword ranks the columns
  into a read-out order. The unknown is the 10-bucket read-order permutation (a small transposition
  climb, CipherModel/SHAPE_ANNEAL so all three -method schedules apply). The primer is transmitted in
  ACA (`-primer 69315`) ⇒ a fast, reliable pure read-order search (~100% from ~60 letters, -logprob).
  Without `-primer` a Gromark-style pre-pass ranks the 10^5 primer space, but that is a documented
  limitation: 10^5 primers each game the n-gram, so the true primer rarely wins the global max even
  at 300 letters. Reuses `gromark_chain_key`. -logprob effectively required (the interleaving gams
  reward-only quadgrams).
- Plus railfence/route/amsco/myszkowski/redefence/cadenus/nihilist/swagman/grille solvers.

Polygraphic squares/cubes/matrix:
- `30` playfair/pf · `31` bifid/bf · `32` trifid/tf · `33` hill.
- `35/36/37` phillips[-c/-rc] · `38/39` twosquare[-v]/ts[v] · `40` foursquare/fs.
- `44` adfgx · `45` adfgvx/adfg · `46/47/48` nihilist-sub[-nc/-m100]/nihsub.
- `54` bazeries/baz · `55` portax/ptx · `59/60/61` slidefair[-var/-beau]/sf.
- `62` seriated-playfair/spf · `63` digrafid/df · `64` cm-bifid/cmb · `65` trisquare/3sq.
- `93` twin-bifid/tbf · `94` twin-trifid/ttf (ACA): TWO Bifid/Trifid messages sharing ONE keyed
  square/cube at DIFFERENT periods (the plaintexts share a common phrase; ACA 100-150 letters
  EACH). The shared key is the crack: anneal a SINGLE square/cube (Bifid/Trifid's own state +
  move set, `score_adjust=0`) while n-gram-scoring the CONCATENATED decrypt of both messages,
  so every key move is judged against ~2x the text (the Tri-Square "scoring length != raw
  length" wiring: both ciphertexts held in the scratch, engine scoring length n1+n2, the
  decrypt hook emits both plaintexts back-to-back). The SECOND ciphertext is supplied with
  `-cipher2 <file>` (config-carried `twincipher_str`, decoded by the same `decode_cipher` as the
  primary; tests inject it in-process). Two INDEPENDENT periods: `bifid/trifid_estimate_periods`
  ranks each message's period, the solver anneals the CROSS PRODUCT of the two top-K lists;
  `-period` pins message 1, `-period2` pins message 2. Twin Bifid recovers the shared square from
  ~110 letters each (well below a lone Bifid's ~500 floor -- the twin advantage) with quadgrams;
  Twin Trifid's 27-cell cube needs QUINTGRAMS and ~210 each. Blind two-period recovery is
  estimator-limited at ACA lengths (columnar IoC is unreliable on short Bifid/Trifid text, ~400+
  needed -- the same Bifid-family limitation) so PIN the periods, which the ACA con lets you find.
  Included in `-type all` only when `-cipher2` is present. `-logprob` (Trifid + quintgrams).
- `82` checkerboard/checker/cb (keyed 5x5 square, 25-letter J->I; plaintext letter -> (row label,
  col label) digraph). Case auto-detected PER AXIS from the ciphertext (an axis with >5 distinct
  labels is complex). Label ORDER is not identifiable (absorbed by a row/col permutation of the
  square, like nihilist-sub). SIMPLE (1 label/axis) ⇒ a free 25-code→25-letter bijection = an
  Aristocrat over the merged codes, on the homophonic incremental fast path. COMPLEX (2 labels/axis,
  homophonic) ⇒ a square-INDEPENDENT per-axis χ² homogeneity pre-pass ranks the label PAIRINGS
  (top-K per axis crossed into engine configs). Simple recovers ~100% from ~130 letters; complex
  needs ~400-600+ (below the ACA 60-90 range — see notable findings). -logprob.
- `84` grandpre/gp (ACA; digit-stream CT parsed from `ciphertext_str`): an N×N word square (N
  6..10, 8×8 standard; rows are words, first column a word, all 26 letters present). Each letter →
  a 2-digit (row,col) code of ANY cell holding it. DECODE is unique (code→letter), ENCODE is
  many-to-one ⇒ a HOMOPHONIC substitution over ≤N² numeric codes → 26 letters. Solver interns the
  codes and reuses the homophonic incremental fast path (the Checkerboard-SIMPLE twin); labels are
  the grid indices so the actual square is recovered directly. Recovers ~93% @285 / ~99% modulo
  rare-letter homophones (Q/X/Z/J) at longer lengths — characterized, below the 99% suite bar.
  -logprob (+ quintgrams).
- `85` syllabary/syll/sy (ACA; digit-stream CT): a 10×10 square of a FIXED 100-element syllabary
  alphabet (26 letters, digits 1-9, a null, 64 syllables like THE/ING/RED). Each plaintext element
  → a 2-digit (row-label,col-label) code; variant spellings (isologs) flatten frequency. DECODE is
  a unique bijection code→token (1-3 letters, LENGTH-CHANGING). Solver climbs the composite 100-token
  code→token permutation (SHAPE_ANNEAL, two-code swap) with the decode TILED into a fixed scoring
  buffer for a length-fair mean n-gram (Fractionated-Morse pattern); the label order folds into the
  map so this one search subsumes all three Known/Unknown variants. BLIND recovery is a documented
  limitation (see notable findings): the syllable tokens are n-gram magnets, so the search games the
  fitness at ACA lengths. -logprob (+ quintgrams).

Layered / concatenated (the Paradigm challenge, GitHub issue #5): an outer substitution
composed over an inner transposition/substitution, all over the KRYPTOS keyed alphabet. A
Quagmire III over a fixed keyed alphabet is a Vigenère in keyed-index space, so stacked Q
layers SUM to one Quagmire of period lcm — give the component periods with `-cyclewordlens
a,b` (searched over their components, ~N/len samples each, not ~N/lcm). Pin the alphabet with
`-plaintextkeyword KRYPTOS`. -logprob (+ quintgrams).
- `88` quagtrans/qtrans: outer Quagmire III o inner columnar transposition, `-depth {0,1,2}`
  (0 = plain Quagmire / PK3; 1 = single columnar / PK4; 2 = double columnar / PK6). Strips the
  outer Quagmire by the MONOGRAM statistic (transposition-invariant), solves the inner columnar
  by n-gram (depth 1: beam over the column order, K≤8 exhaustive, `-mincols/-maxcols`; depth 2:
  reuses `dct_solve_core`), then n-gram coordinate-ascends the cycleword COMPONENTS through the
  recovered transposition. Solves the real PK3/PK4/PK6 to ~100% in seconds.
- `89` hillquag/hq: outer Hill(k×k) o inner Quagmire III, `-period k` (default 3), `-cyclewordlen
  P` (inner Quag period; required — the Hill hides it). Recovers the Hill DECRYPTION matrix by a
  Quag-cycleword-INDEPENDENT statistic — when k|P each mod-P column of D·CT is fed by ONE matrix
  row, so rows are scored by the (cosine-normalised) monogram fit of their columns; anchor k−1
  rows and search the last exhaustively — then strips the Hill and solves the Quagmire. Tries the
  matrix in BOTH the plain A..Z and KRYPTOS-keyed index spaces (the real PK7 is KRYPTOS-keyed).
  Blind Hill is length-limited (floor ~280 chars); the real 279-char PK7 solves. -logprob.

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
- `92` compressocrat/compress/comp (ACA): the fractionation TWIN of Fractionated Morse. A FIXED
  prefix-free {1,2,3} Huffman code maps each letter to a variable-length code (E=31, T=12, …,
  Z=321113); concatenate, pad with `1` to a multiple of 3, group into trigraphs, map each to a
  ciphertext LETTER via a keyed 26-alphabet — the 26 trigraphs of {1,2,3}³ EXCLUDING `333` (never
  occurs; the fracmorse-`xxx` analogue). Ciphertext is normal A–Z and COMPRESSES (shorter than the
  plaintext). Solver = fracmorse verbatim: keyed-alphabet anneal (`*_move_seq`), tile the length-
  changing greedy-prefix-parse decode to C, fold `nv/nt` validity into `score_adjust`. KEY LIMITATION:
  the compression shortens the scoring signal, so short-length fractionation GAMING is worse than
  fracmorse — the ACA 110-150 range (and some keyed alphabets to ~250) is gamed (the true key is not
  the global n-gram max; more restarts/quintgrams don't help), reliable only from ~300; even @300 the
  anneal can stick in a gaming local optimum for an unlucky seed (truth IS the max), so the solver
  test asserts best-of-N seeds (grandpré pattern). -logprob (+ quadgrams; quints don't help).

Biliteral (concealment in cover text):
- `91` baconian/bacon/bac (ACA): each plaintext letter → a 5-symbol a/b group via the FIXED
  24-letter table (I=J, U=V; the a=0/b=1 code is the integer 0..23, and the 8 patterns 24..31 are
  IMPOSSIBLE ⇒ a strong biliteral-VALIDITY signal). The ciphertext is normal English COVER TEXT; a
  hidden CLASSIFIER (a 26-letter a/b labelling) says which cover letter/word stands for a vs b. The
  decode table is fixed, so the classifier is the ONLY unknown. Solver = CipherModel/SHAPE_ANNEAL
  over that labelling: the canonical rules (A-M/N-Z, vowel/consonant × polarity) are single SWEEP
  cells — so a canonical ACA Baconian decodes EXACTLY (~100%) — plus a free-label anneal per grouping
  mode (`-baconmode letter|word|auto`, auto sweeps both; a flip re-parses every group holding that
  letter, guided by the validity reward). Length-changing decode CYCLICALLY TILED to a capped fixed
  length + validity reward folded into score_adjust (the Fractionated-Morse pattern). BLIND
  non-canonical recovery at the ACA ≤25-letter maximum is a documented gaming limitation (n-gram is
  weak there; characterized in the solver test, not in run_tests). -logprob (+ quintgrams).

Rotor machine (own solver, branches early in `solve_cipher`; not the periodic pipeline):
- `95` enigma/enig: the German Enigma. Services Enigma I + naval M3/M4 — rotors I-V (single
  notch), naval VI-VIII (two notches), Greek 4th wheel Beta/Gamma (`-greek`, `-model m4`),
  reflectors UKW-B/C and thin B/C (`-reflector`), plugboard up to 10 pairs, ETW=identity. The
  machine (`enigma.c`) is self-reciprocal with the double-step anomaly. TWO attacks:
  (1) **ciphertext-only** (Gillogly/Ostwald-Weierud/Pound), the default: a two-pass IoC search
  (cheap rings-AAA order ranking, then the full fast-ring position search on the top orders —
  fast-ring searched because a middle-rotor turnover inside the message depends on it),
  middle-ring refinement, then a greedy plugboard warm start + engine anneal scored by n-grams
  (+ cribs). Threaded over the wheel-order search (`-nthreads`); blind 60-order is the ~60×26⁴
  "Bombe workload". Solves Gillogly's 647-letter example blind (`ciphers/tests/enigma_gillogly`).
  Length/plug-limited & probabilistic per key (Gillogly's documented weakness) — reliable for a
  few plugs / longer text, marginal short with many plugs. (2) **Turing-Welchman Bombe** (`-bombe`
  + a crib): a menu + diagonal-board (involution) constraint search recovers the rotor config
  (order + fast ring + start) from a crib in seconds/threaded, then the shared ring/plugboard
  completion finishes it — exact recovery even where the IoC attack is marginal. Pins:
  `-rotors II,I,III` (comma or space list), `-ring`/`-startpos` (letters A..Z or 1-based numbers),
  `-plugboard "EZ RW …"`; all four pinned ⇒ a deterministic known-key decrypt. `-ntopk`,
  `-maxplugs`. Blind M4 ciphertext-only is impractical (26⁴×orders) — pin `-rotors` or use `-bombe`.
  -logprob recommended. KATs pin Gillogly + the Ostwald-Weierud B432 vector. `-enigmaadaptive`
  (Ostwald-Weierud short-message selection) reranks the top rotor configs by a plugboard-COMPLETED
  n-gram before the climb (the default empty-plugboard-IoC ranking drops the true config on
  short/many-plug messages), plus an E-Stecker partial exhaustion below 300 letters — lifts the
  few-plug short floor substantially and roughly doubles many-plug (≥6) short recovery (though full
  ≥6-plug short solves stay rare, near the fundamental floor). Length-guarded so long messages stay
  fast; default off ⇒ bit-identical. See [[enigma-adaptive-ranking]].
- `96` chaocipher/chao: John F. Byrne's Chaocipher (1918). Two 26-letter alphabets — LEFT
  (ciphertext) and RIGHT (plaintext) — each PERMUTED after every enciphered letter, so the
  substitution drifts per position; the KEY is the pair of STARTING alphabets (`chaocipher.c`,
  algorithm per Rubin's "Chaocipher Revealed"; zenith=0/nadir=13, left block [1..13] and right
  block [2..13] each rotate by one after an extra right-disk shift). Rotating BOTH alphabets by
  the same offset is an equivalent key (keyspace 26!·25!; the report canonicalises so LEFT starts
  with A). KEY FINDING: Chaocipher is a NEEDLE for local search — a single starting-alphabet swap
  cascades the whole downstream decrypt, so blind n-gram annealing does NOT reliably recover the
  key at practical lengths (~random; the apparent "solves" in early testing were a plant/solve
  seed collision starting AT the key), like Condi. So two modes: (1) **blind ciphertext-only** —
  a SHAPE_ANNEAL n-gram climb over the two starting alphabets, shipped for completeness but a
  DOCUMENTED LIMITATION (characterised in the solver test, not run_tests). (2) **known-plaintext**
  (a contiguous crib prefix via `-crib`) — a DETERMINISTIC BACKTRACKING RECONSTRUCTION (the method
  that solved the real exhibits; the repo's "deterministic/constructive for needles" rule): walk
  the known pt+ct streams, pinning the shared position of each pt/ct pair into the starting
  alphabets (with an origin-index map back to the start frame) and pruning on contradictions.
  Recovers the key from a crib prefix and decrypts the whole message (full crib ⇒ 100%; a partial
  prefix recovers the tail down to rare unexercised-letter cells). Typical solves <0.1s, worst
  ~10s (key-dependent early branching, node-capped). `run_tests` case `chaocipher_kpa` is a
  full-crib solve. See [[chaocipher-cipher]].

Substitution:
- `28` indep · `29` homophonic · `77` ragbaby/rag · `79` aristocrat/arist · `80` patristocrat/patri
  (one solver core: free 26-perm climbed by n-gram with the homophonic incremental fast path;
  word divisions preserved for the Aristocrat's spaced report, dropped/5-grouped for the
  Patristocrat; -logprob).
- `86` keyphrase/key-phrase/kp (ACA): a 26-letter key phrase IS the cipher alphabet matched to
  straight a..z; the phrase repeats letters so decode is many-to-one (ambiguous) → a partition of
  a..z among the observed ct letters + an inner beam-Viterbi. · `87` affine/af: monoalphabetic
  CT = (a·PT + b) mod 26, gcd(a,26)=1; deterministic-exhaustive 12 multipliers × 26 shifts = 312 keys.

### Key global flags

- `-logprob` (a.k.a. `-azdecrypt`): AZDecrypt-style log10 n-gram fitness with an
  unseen-n-gram floor penalty, vs the default reward-only `log(1+count)` (unseen → 0).
  **Effectively required** for the square/fractionation types; pairs well with
  quintgrams (`-ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt`). Default off ⇒ unchanged.
- `-reversengrams` (`-revngrams`): symmetrize the table so each n-gram and its reversed
  twin share the `max` weight — reads reversed-word text like clean English (for the W168
  alternate-word-reversal hypothesis). Roughly doubles the acceptable solution set.
  Default off ⇒ bit-identical.
- **Compressed `.ngbin` n-gram tables.** A `-ngramfile` whose first 8 bytes are the magic
  `COLNGBIN` is a **dense 8-bit** table (64-byte header + `g_alpha^n` bytes; byte *i* = the
  8-bit-quantized log10 weight of the n-gram at big-endian index *i*, exactly
  `ngram_index_str`'s packing). `load_ngrams` **mmaps** it and `ngram_score` scores via a
  256-entry LUT (`g_ngram_lut[g_ngram_u8[idx]]`) — 1 byte/entry vs 4 (≈4× less RAM), and a
  single mmap vs parsing millions of text lines (quintgram load+RSS ≈ 3×/7× cheaper). It
  **implies `-logprob`** (the only mode stored) and is **alphabet-specific** (a 25-letter
  J→I Bifid table is a different image from the 26-letter one; the reader validates
  `alphabet_size == g_alpha`). Dense ⇒ practical only through order 6 (26⁶ = 309 MB;
  26⁷ overflows the `int` walk index). Build one with `-writengrambin <file>` (dumps the
  exact per-`-type` table then exits — correct for ANY alphabet) or, for plain 26-letter
  tables, `tools/compress_ngrams.py`. `-reversengrams` is rejected with a `.ngbin` (the
  mmap is read-only). The g_ngram_u8==NULL float path is byte-for-byte unchanged, so the
  regression suite stays bit-identical. The **`-spaces` tables** (`-spacesngramfile`, the
  27-symbol {A..Z,' '} alphabet, `load_spaces_ngrams` in `src/core/spaces.c`) accept the
  **same** `.ngbin` format (header `alphabet_size == 27`); build them with
  `compress_ngrams.py --spaces` (the spaces alphabet is fixed, so no per-type variant). The
  `-spaces` pass is a one-shot report Viterbi, so it dequantizes through a per-table LUT.
- `-cribdrag WORD` (or `WORDA|WORDB`): position-free crib. Each word is slid across the
  decrypt and its best-offset partial match rewarded (`-weightcribdrag`, default 36);
  pipe = AND (mean over words). A global toggle in `state_score`, so it works for
  **every** solver. Steers score-driven keys most under `-stochasticcycle`. Default off ⇒
  bit-identical.
- `-crib`: fixed crib pinned to absolute cipher positions (blends in `state_score`).
- `-check-solution-file <file>`: known-plaintext solution, pre-loaded once (uppercased,
  whitespace stripped) into `g_check_solution`. `print_solution_check()` (utils.c) is
  called right alongside every solver's `print_text(decrypted, len)` — both the
  `-verbose` best-improvement dialog and the final report — diffing the candidate
  letter-by-letter (matches print as-is, mismatches as `.`) plus a `NN.NN% correct`
  line. Purely a reporting aid (never touches scoring). Default off ⇒ bit-identical.
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
- `-cyclewordlens a,b[,...]` (layered `quagtrans`): component periods of a composed
  multi-Quagmire; the effective period is their lcm and the cycleword is searched over the
  components (far more samples/parameter than the collapsed period).
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
- **Twin Bifid / Twin Trifid** — the shared key doubles the n-gram signal, so the square/cube
  recovers from far SHORTER text than a lone message (Twin Bifid ~110 letters each vs a lone
  Bifid's ~500; Twin Trifid ~210 each vs a lone Trifid's ~500+). Two independent limitations,
  both asserted/characterized in the solver tests, not solver bugs: (1) **blind PERIOD
  estimation is the bottleneck**, not square recovery — the columnar-IoC estimator is unreliable
  on short Bifid/Trifid ciphertext (the true period misses the top-K until ~400 letters), so a
  fully-blind ACA-length solve is estimator-limited; PIN the periods (`-period`/`-period2`), and
  a blind sweep of the cross-product of the two top-K period lists otherwise multiplies the
  budget. (2) **Twin Trifid effectively needs QUINTGRAMS** (the 27-cell cube's rugged landscape
  is under-signalled by quadgrams at ACA lengths) and is seed-fragile near its ~210 floor, so its
  registry budget is large (`12×300000`) and its run_tests case pins the periods + uses
  quintgrams. Second ciphertext is `-cipher2`; the recovered plaintext (and `.solution`) is the
  two decrypts CONCATENATED. Cribs are not wired (they do not map cleanly onto the concatenated
  two-message plaintext).
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
- **Grandpré** — cryptanalytically a homophonic substitution over ≤N² numeric codes → 26 letters,
  so it reuses the homophonic incremental fast path verbatim. Strong solver (calibrated curve in
  `test_grandpre_solver.c`: ~0.36 @150 → ~0.84 @200 → ~0.93 @285), but plateaus ~98.7% at longer
  lengths — the residual errors are the RARE-LETTER homophones (Q/X/Z/J appear too few times to pin
  their cells), so it stays just under the 99% suite bar. Omitted from `run_tests.sh` and
  characterized via the solver test, like Tridigital. Labels are the grid indices (not scrambled),
  so the actual square is recovered directly — no label ambiguity, unlike Checkerboard/Syllabary.
- **Syllabary** — a substitution over 100 codes → 100 KNOWN 1-3 letter tokens, with a
  length-changing decode tiled into a fixed scoring buffer (Fractionated-Morse pattern). BLIND
  recovery is a documented limitation MORE severe than Checkerboard-complex: the 100-token set
  includes common syllables (THE, AND, ING, RE, …) that act as **n-gram magnets**, so the search
  assembles fluent-but-wrong English and the true map never wins the global n-gram max — it games
  even a 22-symbol (single-letter-token) instance (`test_syllabary_solver.c` prints ~0.06-0.12
  true-letter recovery at ACA lengths; decode correctness itself is pinned by the four ACA isolog
  KAT vectors in `test_syllabary.c`). Whole-word coverage does not discriminate either (the gamed
  output is full of real words). The 100-token composite-map search subsumes all three ACA
  Known/Unknown Coordinates × Keysquare variants; the label order and canonical unmixed token order
  fold into the map and are not identifiable ciphertext-only. The tractable attack (untapped) is the
  KNOWN-KEYSQUARE variant (search only the 10!×10! label perms).
- **Running Key** — the classic "two English streams added together". Fixing the key stream K
  determines the plaintext P = decode(CT, K), so the objective scores BOTH streams jointly (raw
  n-gram of K folded into `*score_adjust` at the SAME scale as `state_score`'s raw n-gram of P —
  weighting K by `weight_ngram` swamps P ~12:1 and induces gaming). BLIND recovery is a documented
  limitation: the joint mean-n-gram does not uniquely pin a solution, so a fluent-but-wrong (K,P)
  pair out-scores the truth (~25-35% recovery, NOT length-improving — running key's classical
  security). The RELIABLE modes are KNOWN-KEY (`-runningkeyfile`, exact — the book-cipher attack /
  the drag-K1/K2/K3 capability) and crib-ANCHORING (a crib letter forces the key letter via
  `rk_key_from_pt`, so pinned positions are held fixed and only the free positions searched —
  recovers a bounded free gap: gap ≤ 8 exact, degrading monotonically with gap size). Asserted
  vs characterized in `test_running_key_solver.c`; `run_tests.sh` case `running_key_known` is the
  known-key mode. Porta maps between alphabet halves, so a crib is only consistent with the
  ciphertext in the opposite half (its key is /2-folded).
- **Enigma** — the ciphertext-only IoC attack is **length- and plug-limited and probabilistic per
  key** (Gillogly's documented weakness): reliable for a few plugs / longer text, marginal for
  short text with many plugs; a near-solution (all-but-one stecker) is common and refinable. The
  correctness-critical subtlety is the **fast-rotor RING** — a coarse rings-AAA position sweep
  garbles right after the first middle-rotor turnover, so the true position is not reliably rank-1;
  the solver must search the fast ring (the 60×26⁴ locations) for the decrypt to stay clean to the
  much-later middle turnover. Done as a two-pass hybrid (cheap AAA order-rank → full fast-ring
  search on the top orders) so blind stays ~30s threaded instead of minutes. The **left (slow)
  rotor ring is unidentifiable** (it never steps) — recovered up to the pos−ring offset, so
  reported rings/pos may differ from the true key while decrypting identically. Ring settings are
  required for the exact plaintext (they set the turnover timing — a rings-AAA-only solve is
  correct only up to the first turnover). **Blind M4 is impractical** (26⁴×orders); M4 is served by
  the machine, the Bombe, and known-key decrypt. The **Bombe (`-bombe` + crib)** is the reliable
  workhorse where IoC is marginal: it recovers the rotor config + plugboard exactly from a crib,
  and is the tool for short messages. Solver test asserts the Bombe + a characterised ciphertext-
  only length/plug curve; `run_tests.sh` case `enigma_gillogly` is the 647-letter Gillogly example
  (100%, wheel order pinned for speed).

## SearchDefaults (per-type schedules)

`init_config()` globals suit the polyalphabetic/transposition score scale. A type whose
score lives on a different scale gets a tuned profile in the compiled-in registry
(`g_search_defaults[]` in `engine.c`, keyed by cipher type), carrying anneal (`a_*`),
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
