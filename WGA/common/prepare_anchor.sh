#!/bin/bash
# Worker: builds the filtered WGA anchor (scaffolds >= MIN_SCAF) and its 2bit/sizes/lastdb index.
#
# Panda carries 73,514 sequences but 94.2% of its length sits in just 24 scaffolds >= 1 Mb.
# Dropping the rest shrinks the index, removes junk target sequence, and keeps the downstream
# MULTIZ and window steps from having to deal with tens of thousands of unplaced contigs --
# windows on those were never usable anyway.
#
# lastdb -c is what makes LAST act on the RepeatMasker soft-masking already present in the FASTA
# (17-58% lowercase across these genomes); without it the masking is recorded but never applied,
# which is what produced TB-scale MAF output on the first attempt.
#
# Usage: prepare_anchor.sh <anchor_fa> <min_scaf_bp> <out_dir> <anchor_name> <threads>
set -euo pipefail

anchor_fa="$1"
min_scaf="$2"
out_dir="$3"
anchor_name="$4"
threads="${5:-1}"

SAMTOOLS=/home/ctools/samtools-1.13/samtools

filtered_fa="$out_dir/${anchor_name}.filtered.fa"
anchor_2bit="$out_dir/${anchor_name}.2bit"
anchor_sizes="$out_dir/${anchor_name}.sizes"
anchor_db="$out_dir/${anchor_name}.lastdb"

mkdir -p "$out_dir"

if [ ! -s "$filtered_fa" ]; then
    echo "[filter] $anchor_name: keeping scaffolds >= ${min_scaf} bp"
    awk -v m="$min_scaf" '$2 >= m {print $1}' "${anchor_fa}.fai" > "$out_dir/${anchor_name}.keep.txt"
    n_keep=$(wc -l < "$out_dir/${anchor_name}.keep.txt")
    n_all=$(wc -l < "${anchor_fa}.fai")
    echo "  keeping $n_keep of $n_all sequences"
    $SAMTOOLS faidx -r "$out_dir/${anchor_name}.keep.txt" "$anchor_fa" > "$filtered_fa"
    $SAMTOOLS faidx "$filtered_fa"
    awk '{t+=$2} END{printf "  filtered anchor: %.2f Gb\n", t/1e9}' "${filtered_fa}.fai"
fi

# Rebuild whenever the 2bit is older than the filtered FASTA: pre-filter copies of these files
# survived an earlier cleanup and a plain existence check silently kept them, leaving the chain/net
# stage with a 73,514-sequence sizes file while the index held 24.
if [ ! -s "$anchor_2bit" ] || [ "$filtered_fa" -nt "$anchor_2bit" ]; then
    echo "[anchor] (re)building 2bit/sizes from the filtered FASTA"
    faToTwoBit "$filtered_fa" "$anchor_2bit"
    twoBitInfo "$anchor_2bit" "$anchor_sizes"
fi

if [ -s "${anchor_db}.prj" ] && [ ! "$filtered_fa" -nt "${anchor_db}.prj" ]; then
    echo "[skip] lastdb index already present and newer than the filtered anchor: $anchor_db"
else
    echo "[index] lastdb -P$threads -c on the filtered anchor"
    lastdb -P"$threads" -c "$anchor_db" "$filtered_fa"
fi

echo "[done] anchor ready: $anchor_db"
