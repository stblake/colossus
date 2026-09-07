#!/usr/bin/env python3
"""Build dense letter n-gram tables (colossus native format) from raw corpus text.

Output format is IDENTICAL to the shipped English tables (ngram_data/english/english_quadgrams.txt,
ngram_data/english/english_quintgrams.txt): one line per observed n-gram,

    NGRAM<space>RAWCOUNT

space-separated, RAW occurrence counts (not rescaled), sorted by count descending.
colossus' load_ngrams() parses this with fscanf("%s\\t%d", ...) -- the "\\t" matches
any whitespace run, so a single space works exactly like the historical tab. Under
-logprob each cell becomes log10(count / total_count), which is scale-invariant, so
raw counts and any global rescaling score identically; raw is chosen here purely to
match the English tables byte-for-byte in structure.

The TAIL IS NEVER PRUNED: every n-gram observed even once is written. Table size is
therefore governed by corpus size, not by the theoretical 26**N space.

Corpus handling (tuned for cipher n-grams, matching a spaceless transposition
plaintext): Project Gutenberg boilerplate is stripped, accents are NFKD-folded to
A-Z (oe/ae ligatures expanded), everything non-letter is dropped, and the window
slides over the CONTINUOUS letter stream (spaces removed => cross-word n-grams are
included). Feed only single-language text -- mixing languages pollutes the table
(e.g. a Spanish Don Quijote leaking into a French corpus).

Usage:
    python3 build_corpus_ngrams.py <files> <sizes> <prefix> [outdir]

    <files>   comma-separated file paths or globs   (e.g. "french_corpus/*.txt")
    <sizes>   comma-separated n-gram orders          (e.g. "3,4,5,6")
    <prefix>  output basename prefix                 (e.g. "frcorp")
    [outdir]  output directory (default: cwd)

Example (rebuild the shipped French tables from the 12-work corpus):
    python3 tools/build_corpus_ngrams.py "tools/french_corpus/*.txt" 3,4,5,6 frcorp ngram_data/french

Names: 3->trigrams, 4->quadgrams, 5->quintgrams, 6->sixgrams (extend NAMES below).
Note: colossus reads counts as int (%d); n=6 works (MAX_NGRAM_SIZE=8) but the runtime
scoring array is g_alpha**6 * 4 bytes (~1.2 GB at g_alpha=26).
"""
import sys, glob, unicodedata, collections, re

NAMES = {1: 'monograms', 2: 'bigrams', 3: 'trigrams', 4: 'quadgrams',
         5: 'quintgrams', 6: 'sixgrams', 7: 'septgrams', 8: 'octgrams'}


def strip_gutenberg(t):
    """Drop Project Gutenberg header/footer boilerplate if the markers are present."""
    m1 = re.search(r'\*\*\* ?START OF.*?\*\*\*', t, re.S)
    if m1:
        t = t[m1.end():]
    if '*** END OF' in t:
        t = t[:t.rfind('*** END OF')]
    return t


def fold(t):
    """Uppercase, expand ligatures, strip accents to A-Z, map all else to a space."""
    t = t.replace('œ', 'oe').replace('Œ', 'OE')   # oe / OE ligatures
    t = t.replace('æ', 'ae').replace('Æ', 'AE')   # ae / AE ligatures
    t = t.replace('ø', 'o').replace('Ø', 'O')     # Danish/Norwegian o-slash (no NFKD decomp)
    t = unicodedata.normalize('NFKD', t)
    t = ''.join(c for c in t if not unicodedata.combining(c))
    return ''.join(c if 'A' <= c <= 'Z' else ' ' for c in t.upper())


def language_of(t):
    """Best-effort Gutenberg 'Language:' header, for a contamination sanity check."""
    m = re.search(r'^Language:[ \t]*(.+)$', t, re.M)
    return m.group(1).strip() if m else '?'


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    patterns = sys.argv[1].split(',')
    sizes = [int(x) for x in sys.argv[2].split(',')]
    prefix = sys.argv[3]
    outdir = sys.argv[4].rstrip('/') if len(sys.argv) > 4 else '.'

    stream = []
    total_letters = 0
    files = sorted({fn for pat in patterns for fn in glob.glob(pat)})
    if not files:
        sys.exit(f"no files matched: {sys.argv[1]}")
    for fn in files:
        raw = open(fn, encoding='utf-8', errors='ignore').read()
        s = fold(strip_gutenberg(raw)).replace(' ', '')
        stream.append(s)
        total_letters += len(s)
        print(f"  {fn}: {len(s):>9} letters  (Language: {language_of(raw)})")
    S = ''.join(stream)
    print(f"corpus letters: {total_letters}  ({len(files)} files)")

    for N in sizes:
        counts = collections.Counter(S[i:i + N] for i in range(len(S) - N + 1))
        out = f"{outdir}/{prefix}_{NAMES[N]}.txt"
        with open(out, 'w') as f:
            for ng, c in counts.most_common():          # descending, like English tables
                f.write(f"{ng} {c}\n")                   # space + RAW count, tail NOT pruned
        mx = max(counts.values())
        print(f"{out}: {len(counts)} {N}-grams (of {26**N} possible), "
              f"max_raw={mx}, min_raw=1 (full tail)")


if __name__ == '__main__':
    main()
