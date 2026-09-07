#!/usr/bin/env python3
"""Compress a plain-text n-gram table into a dense 8-bit .ngbin (colossus native).

The .ngbin format is a 64-byte little-endian header followed by ALPHA**order bytes:
byte i is the 8-bit quantization of the log10 weight of the n-gram whose BIG-ENDIAN
base-ALPHA index is i (AAA=0, AAB=1, ... -- exactly colossus' index packing). colossus
mmaps the payload and dequantizes each byte through a 256-entry LUT
(w = w_floor + byte*w_scale), so scoring is a single O(1) gather at 1 byte/entry instead
of 4 (4x less RAM; a big table loads by one mmap instead of parsing millions of text
lines). See src/core/scoring.c (load_ngrams / write_ngram_bin) and src/core/spaces.c
(load_spaces_ngrams) for the two readers.

WEIGHTING is the AZDecrypt / -logprob scheme: w = log10(count/total) for seen n-grams,
and the floor w = log10(0.01/total) (byte 0) for every unseen one. Log compresses the
dynamic range, which is what makes 8-bit faithful. Pure stdlib -- no numpy.

TWO ALPHABETS:
  (default)  the plain 26-letter A..Z letter tables (ngram_data/<lang>/<lang>_*grams.txt),
             read by load_ngrams. A cipher on a REDUCED alphabet (25-letter J->I squares,
             27/36-symbol types) needs an alphabet-SPECIFIC table -- build those with
             colossus itself (`-writengrambin`), which reuses the exact runtime fold rules.
  --spaces   the 27-symbol {A..Z, ' '} space-inclusive tables
             (ngram_data/<lang>/<lang>_spaces_*grams.txt), read by load_spaces_ngrams for
             the -spaces readability pass. The space is symbol 26; each line's frequency is
             split on the RIGHTMOST space (the window may embed literal spaces). This
             alphabet is fixed (no per-type variants), so this tool covers it exactly.

DENSE => practical only through order 6 (letter 26**6 = 309 MB; spaces 27**6 = 387 MB).

Usage:
    python3 tools/compress_ngrams.py [--spaces] <in.txt> <order> <out.ngbin>

Examples:
    python3 tools/compress_ngrams.py ngram_data/english/english_quintgrams.txt 5 \\
        ngram_data/english/english_quintgrams.ngbin
    python3 tools/compress_ngrams.py --spaces \\
        ngram_data/english/english_spaces_quintgrams.txt 5 \\
        ngram_data/english/english_spaces_quintgrams.ngbin
"""
import sys
import math
import struct

MAGIC = b"COLNGBIN"          # 8 bytes
VERSION = 1
MODE_LOGPROB = 0
HEADER_FMT = "<8sHBBB3sQddd16s"   # magic, version, order, alpha, mode, rsv, n, floor, scale, total, rsv

assert struct.calcsize(HEADER_FMT) == 64, "header must be 64 bytes"


def index_of(window, order, alpha, spaces):
    """Big-endian base-`alpha` index of a window (bytes), or None if it has a symbol
    outside the alphabet -- mirrors ngram_index_str / sp_char_to_sym."""
    idx = 0
    for by in window:
        if spaces and by == 0x20:          # literal space -> symbol 26
            v = 26
        else:
            v = (by & 0xDF) - 65           # uppercase ASCII letter -> 0..25
            if v < 0 or v > 25:
                return None
        idx = idx * alpha + v
    return idx


def main():
    spaces = "--spaces" in sys.argv[1:]
    pos = [a for a in sys.argv[1:] if a != "--spaces"]
    if len(pos) != 3:
        sys.exit(__doc__)
    infile, order, outfile = pos[0], int(pos[1]), pos[2]
    if not 1 <= order <= 6:
        sys.exit("order must be 1..6 (dense 8-bit is impractical above 6)")

    alpha = 27 if spaces else 26
    n_entries = alpha ** order

    # One pass: map each observed window to its base-`alpha` index (skipping any with an
    # out-of-alphabet symbol or the wrong length, exactly like the C loaders). In --spaces
    # mode the frequency is split on the RIGHTMOST space (windows embed literal spaces);
    # otherwise a plain whitespace split. Keep only OBSERVED (index, count) pairs so
    # nothing scales with ALPHA**order except the final byte buffer.
    idxs, cnts = [], []
    total = 0
    kept = skipped = 0
    with open(infile, "rb") as f:
        for line in f:
            line = line.rstrip(b"\r\n")
            if spaces:
                sp = line.rfind(b" ")
                if sp <= 0:
                    continue
                window, freq_b = line[:sp], line[sp + 1:]
            else:
                parts = line.split()
                if len(parts) != 2:
                    continue
                window, freq_b = parts[0], parts[1]
            if len(window) != order or not freq_b.isdigit():
                skipped += 1
                continue
            idx = index_of(window, order, alpha, spaces)
            if idx is None:
                skipped += 1
                continue
            c = int(freq_b)
            idxs.append(idx)
            cnts.append(c)
            total += c
            kept += 1

    if total <= 0:
        total = 1
    floor = math.log10(0.01 / total)
    wmax = math.log10(max(cnts) / total) if cnts else floor + 1.0
    if wmax <= floor:
        wmax = floor + 1.0                    # degenerate guard (empty table)
    scale = (wmax - floor) / 255.0
    inv = 1.0 / scale

    # bytearray defaults to 0 == floor == unseen; only the observed cells get quantized.
    buf = bytearray(n_entries)
    log10 = math.log10
    max_err = 0.0
    for idx, c in zip(idxs, cnts):
        w = log10(c / total)
        b = int(round((w - floor) * inv))
        if b < 0:
            b = 0
        elif b > 255:
            b = 255
        buf[idx] = b
        err = abs((floor + b * scale) - w)
        if err > max_err:
            max_err = err

    header = struct.pack(HEADER_FMT, MAGIC, VERSION, order, alpha, MODE_LOGPROB,
                         b"\x00" * 3, n_entries, floor, scale, float(total), b"\x00" * 16)
    with open(outfile, "wb") as fo:
        fo.write(header)
        fo.write(buf)

    size_mb = (64 + n_entries) / (1024 * 1024)
    kind = "spaces (27-symbol)" if spaces else "letter (26)"
    print(f"wrote {outfile}: {kind}, order {order}, {n_entries} entries ({size_mb:.2f} MB)")
    print(f"  kept {kept} n-grams ({skipped} skipped), total count {total}")
    print(f"  floor={floor:.4f} scale={scale:.6g} max_quant_err={max_err:.4f}")


if __name__ == "__main__":
    main()
