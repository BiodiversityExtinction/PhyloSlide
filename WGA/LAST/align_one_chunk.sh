#!/bin/bash
# Worker: LAST-aligns ONE chunk of a query genome (a chromosome, or a bin of smaller scaffolds)
# against the anchor index. One SLURM array task per chunk.
#
# Chunking exists because memory tracks query scaffold size: a measured 126 Mb scaffold peaked at
# 65 GB at -P1, and whole-genome jobs OOM-killed repeatedly even at 150 GB. Per-chromosome tasks
# bound that structurally instead of relying on tuning.
#
# -P 1 deliberately: measured -P8 vs -P1 on the same scaffold was 475 GB / 7h54m versus
# 65 GB / 7h41m -- eight threads cost ~7x the memory for a 3% slowdown.
#
# Usage: align_one_chunk.sh <anchor_db> <query_fa> <scaffold_list> <out_chunk_maf>
set -euo pipefail

anchor_db="$1"
query_fa="$2"
scaffold_list="$3"
out_maf="$4"

SAMTOOLS=/home/ctools/samtools-1.13/samtools

if [ -s "$out_maf" ]; then
    echo "[skip] chunk already aligned: $(basename "$out_maf")"
    exit 0
fi

work="${out_maf}.tmp.$$"
mkdir -p "$work"
trap 'rm -rf "$work"' EXIT

n_scaf=$(wc -l < "$scaffold_list")
echo "[chunk] $(basename "$out_maf"): $n_scaf scaffold(s) from $(basename "$query_fa")"
$SAMTOOLS faidx -r "$scaffold_list" "$query_fa" > "$work/chunk.fa"
awk '/^>/{next}{t+=length($0)} END{printf "  chunk size: %.1f Mb\n", t/1e6}' "$work/chunk.fa"

lastal -P1 -m 50 -E 0.05 -u 2 "$anchor_db" "$work/chunk.fa" \
    | last-split > "$work/chunk.maf"

mv "$work/chunk.maf" "$out_maf"
echo "[done] $(basename "$out_maf"): $(du -h "$out_maf" | cut -f1)"
