#!/usr/bin/env python3
"""Enumerate + download Project Gutenberg books for one language into a corpus dir.

Uses the gutendex API (https://gutendex.com) to list every book for a language,
picks each book's plain-text UTF-8 format, and downloads it politely (small delay,
skip-if-present, skip failures). Feeds tools/build_corpus_ngrams.py.

Usage:
    python3 gutenberg_fetch.py enumerate <lang> <out_list> [max_books]         # write id<TAB>url<TAB>title (top-N by popularity)
    python3 gutenberg_fetch.py download  <out_list> <corpus_dir> [delay_s] [max_books]

Example (German):
    python3 tools/gutenberg_fetch.py enumerate de tools/german_corpus/list.tsv
    python3 tools/gutenberg_fetch.py download tools/german_corpus/list.tsv \\
            tools/german_corpus/books 0.5
"""
import sys, os, time, json, urllib.request, urllib.error

UA = "colossus-ngram-builder/1.0 (research; polite)"


def _get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout).read()


def _get_retry(url, tries=4):
    for a in range(tries):
        try:
            return _get(url)
        except Exception as e:                       # transient throttle/timeout -> back off
            if a == tries - 1:
                raise
            time.sleep(1.5 * (a + 1))


def enumerate_lang(lang, out_list, max_books=0):
    # Write rows incrementally so partial progress survives an interruption.
    # gutendex is popularity-ordered, so max_books stops early with the top-N (essential
    # for huge languages like English ~70k books when only the top ~1000 are wanted).
    max_books = int(max_books)
    url = f"https://gutendex.com/books/?languages={lang}"
    n, page = 0, 0
    with open(out_list, "w") as f:
        while url:
            page += 1
            d = json.loads(_get_retry(url))
            for b in d["results"]:
                fmts = b.get("formats", {})
                txt = fmts.get("text/plain; charset=utf-8")
                if not txt:
                    for k, v in fmts.items():
                        if k.startswith("text/plain") and not v.endswith(".zip"):
                            txt = v
                            break
                if txt and not txt.endswith(".zip"):
                    title = (b.get("title") or "").replace("\t", " ").replace("\n", " ")
                    f.write(f"{b['id']}\t{txt}\t{title[:80]}\n")
                    n += 1
            f.flush()
            print(f"  page {page}: {n} books so far", flush=True)
            if max_books and n >= max_books:
                break
            url = d.get("next")
    print(f"enumerated {n} {lang} books -> {out_list}")


def download(out_list, corpus_dir, delay=0.5, max_books=0):
    os.makedirs(corpus_dir, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(out_list) if l.strip()]
    if max_books:                          # cap to the first N (list is popularity-ordered)
        rows = rows[:int(max_books)]
    ok = skip = fail = 0
    for n, parts in enumerate(rows, 1):
        bid, url = parts[0], parts[1]
        dest = os.path.join(corpus_dir, f"g_{bid}.txt")
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            skip += 1
            continue
        try:
            data = _get(url)
            with open(dest, "wb") as f:
                f.write(data)
            ok += 1
        except Exception as e:                      # network/404/etc -> skip, keep going
            fail += 1
            print(f"  FAIL g_{bid}: {e}", flush=True)
        if n % 50 == 0:
            print(f"  {n}/{len(rows)}  ok={ok} skip={skip} fail={fail}", flush=True)
        time.sleep(float(delay))
    print(f"DONE: ok={ok} skip={skip} fail={fail} of {len(rows)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "enumerate":
        enumerate_lang(sys.argv[2], sys.argv[3],
                       sys.argv[4] if len(sys.argv) > 4 else 0)
    elif cmd == "download":
        download(sys.argv[2], sys.argv[3],
                 sys.argv[4] if len(sys.argv) > 4 else 0.5,
                 sys.argv[5] if len(sys.argv) > 5 else 0)
    else:
        sys.exit(__doc__)
