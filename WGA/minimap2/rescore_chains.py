#!/usr/bin/env python3
"""
Rewrite chain scores (and renumber ids) in a chain file produced by paf2chain.

paf2chain copies the PAF mapq into the chain score field, so every chain comes out with score 255.
chainNet defaults to -minScore=2000 and therefore discards all of them, producing an empty net and
ultimately an empty MAF -- silently, with a zero exit status.

Lowering -minScore would fix the emptiness but not the underlying problem: netting ranks competing
chains by score, so leaving every chain on 255 makes that ranking arbitrary. Instead we score each
chain by its number of aligned bases (the sum of the block sizes in the chain body), which is a
sensible proxy for alignment quality and is on a comparable scale to what axtChain computes for
the LAST arm.

Ids are renumbered from 1 at the same time: paf2chain numbers from 0, and chainSort then collides
that chain onto id 1, which netToAxt rejects as a duplicate.

Usage: rescore_chains.py < in.chain > out.chain
"""
import sys


def main():
    header = None
    body = []
    total = 0
    cid = 0
    out = sys.stdout

    def flush():
        nonlocal header, body, total, cid
        if header is None:
            return
        cid += 1
        f = header.split("\t")
        if len(f) >= 13:
            f[1] = str(total)   # score <- aligned bases
            f[12] = str(cid)    # id   <- renumber from 1
        out.write("\t".join(f) + "\n")
        for b in body:
            out.write(b + "\n")
        out.write("\n")
        header, body, total = None, [], 0

    for line in sys.stdin:
        line = line.rstrip("\n")
        if line.startswith("chain\t") or line.startswith("chain "):
            flush()
            header = line.replace(" ", "\t") if not line.startswith("chain\t") else line
        elif not line.strip():
            continue
        elif header is not None:
            body.append(line)
            try:
                total += int(line.split("\t")[0].split()[0])
            except (ValueError, IndexError):
                pass
    flush()


if __name__ == "__main__":
    main()
