#!/bin/bash
# Worker: concatenates one species' per-chunk MAFs, then chains, nets, and filters to syntenic
# blocks -- producing the single <species>.syntenic.maf that MULTIZ consumes. One task per species.
#
# Chunking the query is invisible downstream: every chunk was aligned against the same anchor
# index, so concatenating them reconstitutes the full pairwise alignment before chain/net.
#
# Usage: finalize_species_wga.sh <anchor_name> <anchor_2bit> <anchor_sizes> <codename> \
#                                <query_fa> <chunk_dir> <out_maf>
set -euo pipefail

anchor_name="$1"
anchor_2bit="$2"
anchor_sizes="$3"
codename="$4"
query_fa="$5"
chunk_dir="$6"
out_maf="$7"

if [ -s "$out_maf" ]; then
    echo "[skip] $codename already finalized: $out_maf"
    exit 0
fi

work="$(dirname "$out_maf")/${codename}.finalize"
rm -rf "$work"; mkdir -p "$work"

chunks=$(find "$chunk_dir" -name "*.maf" -size +0 | sort)
n=$(echo "$chunks" | grep -c . || true)
[ "$n" -gt 0 ] || { echo "[error] $codename: no chunk MAFs in $chunk_dir" >&2; exit 1; }
echo "[finalize] $codename: merging $n chunk MAF(s)"

# Keep the header from the first chunk only; strip it from the rest.
first=1
: > "$work/all.maf"
while read -r m; do
    [ -z "$m" ] && continue
    if [ "$first" -eq 1 ]; then
        cat "$m" >> "$work/all.maf"
        first=0
    else
        grep -v "^#" "$m" >> "$work/all.maf" || true
    fi
done <<< "$chunks"
echo "  merged MAF: $(du -h "$work/all.maf" | cut -f1)"

# Drop blocks whose anchor scaffold is not in the filtered anchor (the lastdb index predates the
# >=1 Mb filter, so chunks also aligned to small scaffolds like the mitochondrion, which chain/net
# rejects: "NC_009492.1 is not in Panda.2bit").
python3 "$(dirname "${BASH_SOURCE[0]}")/filter_maf_to_anchor.py" --sizes "$anchor_sizes" \
    < "$work/all.maf" > "$work/all.filtered.maf"
mv "$work/all.filtered.maf" "$work/all.maf"
echo "  after anchor filter: $(du -h "$work/all.maf" | cut -f1)"

query_2bit="$work/${codename}.2bit"
query_sizes="$work/${codename}.sizes"
faToTwoBit "$query_fa" "$query_2bit"
twoBitInfo "$query_2bit" "$query_sizes"

maf-convert axt "$work/all.maf" > "$work/all.axt"
axtChain -linearGap=medium "$work/all.axt" "$anchor_2bit" "$query_2bit" "$work/all.chain"
chainSort "$work/all.chain" "$work/all.sorted.chain"
chainNet "$work/all.sorted.chain" "$anchor_sizes" "$query_sizes" \
    "$work/target.net" "$work/query.net"
netSyntenic "$work/target.net" "$work/syntenic.net"
netToAxt "$work/syntenic.net" "$work/all.sorted.chain" "$anchor_2bit" "$query_2bit" \
    "$work/syntenic.axt"
axtToMaf "$work/syntenic.axt" "$anchor_sizes" "$query_sizes" "$out_maf" \
    -tPrefix="${anchor_name}." -qPrefix="${codename}."

rm -rf "$work"
echo "[done] $codename -> $out_maf ($(du -h "$out_maf" | cut -f1))"
