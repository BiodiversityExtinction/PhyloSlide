#!/bin/bash
# Submits the minimap2 cross-check arm: every assembly aligned to the same filtered Panda anchor,
# through the same chain/net/syntenic tail as the LAST arm, into results/wga/pairwise_minimap2/.
#
# Purpose is methodological robustness, not replacement -- if the WGA topology is the same under
# LAST and under a completely different aligner, that is a strong statement for the paper. LAST
# stays the primary (it is the standard for this, and matches the rhino-paper method).
#
# Depends on LAST/align_to_anchor.sh stage 1 having built the filtered anchor + 2bit/sizes; pass
# --after <jobid> to chain onto it, or run once that has completed.
#
# Usage: ./align_to_anchor.sh --projdir <analysis dir> [--dry-run] [--after <jobid>] [--concurrent N]
set -euo pipefail

# Archived from the reference-bias benchmark. PROJDIR is the ANALYSIS directory (the one holding
# data/ and results/), which is not this repository: set it with --projdir or PHYLOSLIDE_PROJDIR.
PROJDIR="${PHYLOSLIDE_PROJDIR:-}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANCHOR="Panda"
THREADS=16
MEM=64G
CONCURRENT=4
EXCLUDE="node01,node02,node03,node04,node05,node06,node08,node09,node12,node13"
CONDA_SH="/home/people/micwe/Software/miniconda3/etc/profile.d/conda.sh"
DRY_RUN=""
AFTER=""

while [ $# -gt 0 ]; do
    case "$1" in
        --projdir)  PROJDIR="$2"; shift 2 ;;
        --projdir) PROJDIR="$2"; shift 2 ;;
        --dry-run)    DRY_RUN="yes"; shift ;;
        --after)      AFTER="$2"; shift 2 ;;
        --concurrent) CONCURRENT="$2"; shift 2 ;;
        *) echo "[error] unknown argument: $1" >&2; exit 1 ;;
    esac
done
[ -n "$PROJDIR" ] || { echo "[error] set --projdir or PHYLOSLIDE_PROJDIR (the analysis directory holding data/ and results/)" >&2; exit 1; }

GENOMES="$PROJDIR/data/genomes"
LASTDIR="$PROJDIR/results/wga/pairwise"
OUTDIR="$PROJDIR/results/wga/pairwise_minimap2"
META="$PROJDIR/data/metadata/genomes.tsv"

RUNDIR="$PROJDIR/logs/wga_minimap2/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTDIR" "$RUNDIR"
cp "$SCRIPTDIR/align_minimap2.sh" "$SCRIPTDIR/rescore_chains.py" "$RUNDIR/"
# rescore_chains.py MUST be snapshotted too: the worker resolves it relative to itself, so
# copying only the worker left it looking in RUNDIR for a file that was never there.

# built by 03's stage 1 -- the same filtered anchor both arms align against
ANCHOR_FA="$LASTDIR/${ANCHOR}.filtered.fa"
ANCHOR_2BIT="$LASTDIR/${ANCHOR}.2bit"
ANCHOR_SIZES="$LASTDIR/${ANCHOR}.sizes"

joblist="$RUNDIR/joblist.txt"
: > "$joblist"
tail -n +2 "$META" | while IFS=$'\t' read -r codename accession common sci source role local_path; do
    [ -z "$codename" ] && continue
    [ "$codename" = "$ANCHOR" ] && continue
    if [ -s "$OUTDIR/${codename}.syntenic.maf" ]; then
        echo "[done]  $codename (minimap2 already aligned)"
        continue
    fi
    echo "$codename" >> "$joblist"
    echo "[queue] $codename vs $ANCHOR (minimap2 asm10)"
done

n=$(wc -l < "$joblist")
[ "$n" -gt 0 ] || { echo "Nothing to align."; exit 0; }

array="$RUNDIR/minimap2_array.sbatch.sh"
cat > "$array" <<EOF
#!/bin/bash
#SBATCH --job-name=mm2_wga
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$MEM
#SBATCH --array=1-${n}%${CONCURRENT}
#SBATCH --exclude=$EXCLUDE
#SBATCH --output=$RUNDIR/mm2_%a.%j.out
source "$CONDA_SH"
conda activate phyloslide_wga
set -euo pipefail

codename=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$joblist")
bash "$RUNDIR/align_minimap2.sh" "$ANCHOR" "$ANCHOR_FA" "$ANCHOR_2BIT" "$ANCHOR_SIZES" \\
    "\$codename" "$GENOMES/\${codename}.fa" "$OUTDIR/\${codename}.syntenic.maf" $THREADS
EOF

echo
echo "[submit] minimap2 array: $n species, max $CONCURRENT concurrent, $THREADS threads, $MEM"
if [ -z "$DRY_RUN" ]; then
    args=(); [ -n "$AFTER" ] && args=(--dependency=afterok:"$AFTER")
    sbatch --parsable "${args[@]}" "$array"
fi
echo "Run dir: $RUNDIR"
