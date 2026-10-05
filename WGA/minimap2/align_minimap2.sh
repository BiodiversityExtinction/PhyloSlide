#!/bin/bash
# Worker: aligns one assembly to the anchor with minimap2, then runs the SAME chain/net/syntenic
# tail as the LAST arm, so the two are directly comparable downstream.
#
# This is a cross-check arm, not the primary: LAST is the standard for this and matches the
# rhino-paper method. minimap2 is ~50x faster but less sensitive in divergent and repetitive
# regions, so agreement between the two is evidence the WGA topology is not an artefact of the
# aligner. asm10 suits the ~4-6% divergence between Panda and the other ursids.
#
# No chunking needed: minimap2 handles whole mammalian genomes in hours at modest memory.
#
# Usage: align_minimap2.sh <anchor_name> <anchor_fa> <anchor_2bit> <anchor_sizes> \
#                          <codename> <query_fa> <out_maf> <threads>
set -euo pipefail

anchor_name="$1"
anchor_fa="$2"
anchor_2bit="$3"
anchor_sizes="$4"
codename="$5"
query_fa="$6"
out_maf="$7"
threads="${8:-8}"

if [ -s "$out_maf" ]; then
    echo "[skip] $codename already aligned (minimap2): $out_maf"
    exit 0
fi

work="$(dirname "$out_maf")/${codename}.mm2work"
mkdir -p "$work"

# Reuse a PAF left by a previous attempt: the alignment itself takes ~30 min per species and is
# unaffected by the downstream chain/net problems that have caused the retries here.
if [ -s "$work/${codename}.paf" ]; then
    echo "[minimap2] $codename: reusing existing PAF ($(du -h "$work/${codename}.paf" | cut -f1))"
else
    echo "[minimap2] $anchor_name vs $codename (asm10, threads=$threads)"
    minimap2 -cx asm10 -t "$threads" --cs "$anchor_fa" "$query_fa" > "$work/${codename}.paf.tmp"
    mv "$work/${codename}.paf.tmp" "$work/${codename}.paf"
    echo "  PAF: $(du -h "$work/${codename}.paf" | cut -f1)"
fi

# paf2chain copies the PAF mapq into the chain score, so every chain arrives with score 255 --
# and chainNet defaults to -minScore=2000, so it silently discarded all 75k chains and produced an
# empty net, an empty MAF, and a zero exit status. rescore_chains.py scores each chain by its
# aligned bases (restoring a meaningful ranking for netting, rather than just lowering the
# threshold) and renumbers ids from 1 (paf2chain starts at 0, which chainSort collides onto 1).
paf2chain -i "$work/${codename}.paf" \
  | python3 "$(dirname "${BASH_SOURCE[0]}")/rescore_chains.py" \
  > "$work/${codename}.chain"
n_chain=$(awk 'BEGIN{FS="\t"} /^chain\t/{n++} END{print n+0}' "$work/${codename}.chain")
n_ids=$(awk 'BEGIN{FS="\t"} /^chain\t/{ids[$NF]++} END{print length(ids)}' "$work/${codename}.chain")
echo "  chains: $n_chain (distinct ids: $n_ids)"
[ "$n_chain" = "$n_ids" ] || { echo "[error] $codename: chain ids not unique" >&2; exit 1; }

query_2bit="$work/${codename}.2bit"
query_sizes="$work/${codename}.sizes"
faToTwoBit "$query_fa" "$query_2bit"
twoBitInfo "$query_2bit" "$query_sizes"

chainSort "$work/${codename}.chain" "$work/${codename}.sorted.chain"
echo "  [size] sorted.chain: $(wc -l < "$work/${codename}.sorted.chain" 2>/dev/null || echo NA) lines"
chainNet "$work/${codename}.sorted.chain" "$anchor_sizes" "$query_sizes" \
    "$work/target.net" "$work/query.net"
echo "  [size] target.net: $(wc -l < "$work/target.net") lines"
netSyntenic "$work/target.net" "$work/syntenic.net"
echo "  [size] syntenic.net: $(wc -l < "$work/syntenic.net" 2>/dev/null || echo NA) lines"
netToAxt "$work/syntenic.net" "$work/${codename}.sorted.chain" "$anchor_2bit" "$query_2bit" \
    "$work/syntenic.axt"
echo "  [size] syntenic.axt: $(wc -l < "$work/syntenic.axt" 2>/dev/null || echo NA) lines"
axtToMaf "$work/syntenic.axt" "$anchor_sizes" "$query_sizes" "$out_maf" \
    -tPrefix="${anchor_name}." -qPrefix="${codename}."

n_blocks=$(grep -c "^a " "$out_maf" 2>/dev/null || echo 0)
[ "$n_blocks" -gt 0 ] || { echo "[error] $codename: MAF has no alignment blocks" >&2; exit 1; }
echo "  MAF blocks: $n_blocks"
rm -rf "$work"
echo "[done] $codename -> $out_maf ($(du -h "$out_maf" | cut -f1))"
