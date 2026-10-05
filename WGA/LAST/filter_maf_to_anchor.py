#!/usr/bin/env python3
"""
Keep only MAF blocks whose anchor (first 's' row) scaffold is in the filtered anchor set.

The lastdb index was built before the >=1 Mb anchor filter existed, and prepare_anchor.sh's
skip-if-exists guard then kept that stale index -- so every chunk was aligned against all 73,514
Panda sequences rather than the 24 we kept. chain/net then dies with e.g.
"NC_009492.1 is not in Panda.2bit" (that one is the 16.8 kb mitochondrion).

Rebuilding the index would mean redoing all 303 chunk alignments. Those extra blocks are precisely
what the >=1 Mb filter was supposed to remove, so dropping them here reaches the same end state
without repeating days of compute.

Usage: filter_maf_to_anchor.py --sizes Panda.sizes < in.maf > out.maf
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", required=True, help="anchor .sizes (or .fai); first column is kept")
    args = ap.parse_args()

    keep = set()
    with open(args.sizes) as fh:
        for line in fh:
            f = line.split()
            if f:
                keep.add(f[0])

    out = sys.stdout
    block, s_rows, anchor_ok = [], 0, False
    kept = dropped = 0

    def flush():
        nonlocal block, s_rows, anchor_ok, kept, dropped
        if block:
            if anchor_ok and s_rows >= 2:
                out.write("".join(block))
                out.write("\n")
                kept += 1
            else:
                dropped += 1
        block, s_rows, anchor_ok = [], 0, False

    for line in sys.stdin:
        if line.startswith("#"):
            if not block:
                out.write(line)
            continue
        if line.startswith("a"):
            flush()
            block = [line]
            continue
        if line.startswith("s") and block:
            s_rows += 1
            if s_rows == 1:
                f = line.split()
                anchor_ok = len(f) >= 2 and f[1] in keep
            block.append(line)
            continue
        if not line.strip():
            flush()
            continue
        if block:
            block.append(line)
    flush()

    print(f"  [filter] kept {kept:,} blocks, dropped {dropped:,} off-anchor blocks",
          file=sys.stderr)


if __name__ == "__main__":
    main()
