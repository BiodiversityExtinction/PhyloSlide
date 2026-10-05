#!/usr/bin/env python3
"""
Project a pairwise syntenic MAF onto anchor coordinates, producing a per-species pseudo-genome
FASTA the same length as the anchor.

This is the "option b" route for the WGA arm: instead of extracting windows out of a multi-way
alignment, we render each species as an anchor-coordinate FASTA, so the WGA arm can be fed to
PhyloSlide through the exact same code path as the reference-mapped arms. The two then differ in
exactly one respect -- whether a species' sequence came from genome alignment or from read
mapping -- with identical windowing, filtering, tree inference and concordance code downstream.

Note this makes MULTIZ unnecessary. A multiple alignment projected back onto reference
coordinates is just the union of the pairwise projections, for every column the reference is
present in; the multi-way merge only adds columns where the reference has a gap, and those have
no anchor coordinate to live at. Skipping it also avoids roast's guide tree.

Deliberate choices, both to match what read mapping does:
  - query gaps (deletions relative to the anchor) are written as N, exactly as angsd -dofasta
    leaves uncovered sites as N. This conflates deletion with missing data, but the mapped arms
    conflate them the same way, which is the point.
  - anchor-gap columns (insertions in the query) are dropped; they have no anchor coordinate.
    Read mapping discards these too.

Usage: project_maf_to_anchor.py --maf X.syntenic.maf --anchor-fai Panda.filtered.fa.fai \\
                                --anchor-prefix Panda. --query-name Brown --out Brown.fa
"""
import argparse
import gzip
import sys


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maf", required=True)
    ap.add_argument("--anchor-fai", required=True)
    ap.add_argument("--anchor-prefix", default="")
    ap.add_argument("--query-name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--line-width", type=int, default=60)
    args = ap.parse_args()

    # anchor scaffolds, in .fai order, pre-filled with N
    order, seqs = [], {}
    with open(args.anchor_fai) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 2:
                order.append(p[0])
                seqs[p[0]] = bytearray(b"N" * int(p[1]))

    blocks = filled = conflicts = skipped = 0
    with open_maybe_gz(args.maf) as fh:
        rows = []
        for line in fh:
            if line.startswith("a"):
                rows = []
            elif line.startswith("s"):
                rows.append(line.split())
                if len(rows) == 2:
                    blocks += 1
                    t, q = rows[0], rows[1]
                    # s src start size strand srcSize text
                    tname = t[1]
                    if args.anchor_prefix and tname.startswith(args.anchor_prefix):
                        tname = tname[len(args.anchor_prefix):]
                    if tname not in seqs:
                        skipped += 1
                        rows = []
                        continue
                    if t[4] != "+":
                        # axtToMaf emits the target on +; a - target would need the whole block
                        # reverse-complemented, so refuse rather than silently mis-place bases.
                        sys.exit(f"Unexpected '-' strand on anchor row for {tname}")
                    tpos = int(t[2])
                    ttext, qtext = t[6], q[6]
                    dest = seqs[tname]
                    for tc, qc in zip(ttext, qtext):
                        if tc == "-":
                            continue          # insertion in query: no anchor coordinate
                        if qc != "-":
                            b = ord(qc.upper())
                            if dest[tpos] != 78:   # already written by another chain
                                conflicts += 1
                            else:
                                dest[tpos] = b
                                filled += 1
                        tpos += 1
                    rows = []

    with open(args.out, "w") as out:
        for name in order:
            out.write(f">{name}\n")
            s = seqs[name]
            for i in range(0, len(s), args.line_width):
                out.write(s[i:i + args.line_width].decode() + "\n")

    total = sum(len(v) for v in seqs.values())
    print(f"[{args.query_name}] blocks={blocks} filled={filled:,} "
          f"({100.0*filled/total:.1f}% of anchor) conflicts={conflicts} "
          f"skipped_blocks_off_anchor={skipped}", file=sys.stderr)


if __name__ == "__main__":
    main()
