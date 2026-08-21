#!/usr/bin/env python3
"""Build SPACE-INCLUSIVE character n-gram tables (colossus -spaces format) from raw corpus text.

Sibling of build_corpus_ngrams.py. That script slides its window over a CONTINUOUS letter
stream (spaces removed, cross-word n-grams included) and feeds load_ngrams(). THIS script
KEEPS the word spaces: the window slides over a {A..Z, ' '} stream, so a window may contain
embedded single spaces, and the resulting table feeds the -spaces readability pass
(src/core/spaces.c load_spaces_ngrams()), NOT load_ngrams().

Output format is exactly what load_spaces_ngrams() parses -- one line per observed window,

    WINDOW<space>RAWCOUNT

where WINDOW is `order` characters of A-Z and/or embedded LITERAL single spaces, and the
count is separated by the RIGHTMOST space (so a window that itself ends in a space, e.g.
"THE ", prints as "THE  <count>" -- two spaces -- and still parses correctly). Examples
(order 5): "AAAAA 19157" (no space) and "AAA A 1264" (window "AAA" + space + "A"). The space
MUST be a literal space: load_spaces_ngrams() maps ' '->symbol 26 and A-Z->0..25, and skips
any line containing another character, so an underscore (or any other sentinel) would yield
an empty table.

Counts are RAW occurrence counts, sorted descending, tail NEVER pruned -- same convention as
the letter-only tables. Under the loader's log10(count/total) conversion this is
scale-invariant, so raw counts score identically to any rescaling; raw is chosen to match the
existing tables structurally.

Corpus handling matches build_corpus_ngrams.py's fold() (Gutenberg boilerplate stripped,
accents NFKD-folded to A-Z, oe/ae ligatures expanded), with ONE difference: instead of
removing every space, runs of spaces (fold() emits one space per non-letter char) are
COLLAPSED to a single space and the ends stripped, and books are joined with a single space
(a real word boundary, and it blocks false cross-book letter adjacency). Consecutive spaces
therefore never occur in a window -- matching the -spaces Viterbi's "no double spaces" rule.
Feed only single-language text.

Usage:
    python3 build_corpus_spaces_ngrams.py <files> <sizes> <prefix> [outdir]

    <files>   comma-separated file paths or globs   (e.g. "tools/english_corpus/books/*.txt")
    <sizes>   comma-separated n-gram orders          (e.g. "2,3,4,5,6")
    <prefix>  output basename prefix                 (e.g. "english_spaces")
    [outdir]  output directory (default: cwd)

Example (build the shipped English space-inclusive tables, orders 2-6):
    python3 tools/build_corpus_spaces_ngrams.py "tools/english_corpus/books/*.txt" \\
            2,3,4,5,6 english_spaces ngram_data/english

Names: 2->bigrams, 3->trigrams, 4->quadgrams, 5->quintgrams, 6->sixgrams (see NAMES).
Cap orders at 6: load_spaces_ngrams() rejects a dense 27^order table above ~6e8 cells
(27^6 ~ 3.87e8 is the max feasible; 27^7 ~ 1.05e10 is rejected).
"""
import sys, os, glob, unicodedata, collections, re

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
    os.makedirs(outdir, exist_ok=True)

    stream = []
    total_chars = 0
    files = sorted({fn for pat in patterns for fn in glob.glob(pat)})
    if not files:
        sys.exit(f"no files matched: {sys.argv[1]}")
    for fn in files:
        raw = open(fn, encoding='utf-8', errors='ignore').read()
        # KEEP spaces: collapse fold()'s per-non-letter spaces to one, strip the ends.
        s = re.sub(r' +', ' ', fold(strip_gutenberg(raw))).strip()
        stream.append(s)
        total_chars += len(s)
        print(f"  {fn}: {len(s):>9} chars  (Language: {language_of(raw)})")
    # Join books with a single space -- a real word boundary, no false cross-book adjacency.
    S = ' '.join(stream)
    print(f"corpus chars (incl. spaces): {len(S)}  ({len(files)} files)")

    for N in sizes:
        counts = collections.Counter(S[i:i + N] for i in range(len(S) - N + 1))
        out = f"{outdir}/{prefix}_{NAMES[N]}.txt"
        with open(out, 'w') as f:
            for ng, c in counts.most_common():          # descending, like the letter tables
                f.write(f"{ng} {c}\n")                   # window + space + RAW count (rightmost
                                                          # space splits the count on load)
        mx = max(counts.values())
        print(f"{out}: {len(counts)} {N}-grams (of {27**N} possible over 27 symbols), "
              f"max_raw={mx}, min_raw=1 (full tail)")


if __name__ == '__main__':
    main()
