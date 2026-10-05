#!/bin/bash
# Projects each pairwise syntenic MAF onto Panda anchor coordinates, giving one pseudo-genome
# FASTA per species in results/wga/projected/. Those are then fed to PhyloSlide exactly like the
# reference-mapped consensus FASTAs, so the WGA arm and the Panda reference arm run through an
# identical code path on identical window coordinates and differ only in how each species'
# sequence was obtained (genome alignment vs read mapping).
#
# The anchor itself needs no projection -- it is its own coordinate system -- so Panda.filtered.fa
# is simply linked in as the Panda taxon.
#
# MULTIZ is deliberately not used; see common/project_maf_to_anchor.py for why the multi-way
# merge adds nothing once you project back to anchor coordinates.
#
# Usage: ./project_to_anchor.sh --projdir <analysis dir> [--dry-run] [--arm last|minimap2] [--after <jobid>]
set -euo pipefail

# Archived from the reference-bias benchmark. PROJDIR is the ANALYSIS directory (the one holding
# data/ and results/), which is not this repository: set it with --projdir or PHYLOSLIDE_PROJDIR.
PROJDIR="${PHYLOSLIDE_PROJDIR:-}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANCHOR="Panda"
ARM="last"
DRY_RUN=""
AFTER=""

while [ $# -gt 0 ]; do
    case "$1" in
        --projdir)  PROJDIR="$2"; shift 2 ;;
        --projdir) PROJDIR="$2"; shift 2 ;;
        --dry-run) DRY_RUN="yes"; shift ;;
        --arm)     ARM="$2"; shift 2 ;;
        --after)   AFTER="$2"; shift 2 ;;
        *) echo "[error] unknown argument: $1" >&2; exit 1 ;;
    esac
done
[ -n "$PROJDIR" ] || { echo "[error] set --projdir or PHYLOSLIDE_PROJDIR (the analysis directory holding data/ and results/)" >&2; exit 1; }

case "$ARM" in
    last)     MAFDIR="$PROJDIR/results/wga/pairwise" ;;
    minimap2) MAFDIR="$PROJDIR/results/wga/pairwise_minimap2" ;;
    *) echo "[error] --arm must be last or minimap2" >&2; exit 1 ;;
esac
OUTDIR="$PROJDIR/results/wga/projected_${ARM}"
ANCHOR_FA="$PROJDIR/results/wga/pairwise/${ANCHOR}.filtered.fa"
EXCLUDE="compute04,node01,node02,node03,node04,node05,node06,node08,node09,node12,node13"
CONDA_SH="/home/people/micwe/Software/miniconda3/etc/profile.d/conda.sh"

RUNDIR="$PROJDIR/logs/wga_project/run_${ARM}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTDIR" "$RUNDIR"
cp "$SCRIPTDIR/project_maf_to_anchor.py" "$RUNDIR/"

joblist="$RUNDIR/joblist.tsv"
: > "$joblist"
for maf in "$MAFDIR"/*.syntenic.maf; do
    [ -e "$maf" ] || continue
    codename=$(basename "$maf" .syntenic.maf)
    out="$OUTDIR/${codename}.fa"
    if [ -s "$out" ]; then echo "[done]  $codename"; continue; fi
    printf '%s\t%s\t%s\n' "$codename" "$maf" "$out" >> "$joblist"
    echo "[queue] $codename"
done

n=$(wc -l < "$joblist")
if [ "$n" -eq 0 ]; then
    echo "Nothing to project (no syntenic MAFs yet, or all done)."
    exit 0
fi

array="$RUNDIR/project_array.sbatch.sh"
# $n must already be set here: it becomes the --array bound below. Without that directive SLURM
# runs the script once with SLURM_ARRAY_TASK_ID unset, which under `set -u` dies immediately.
[ "$n" -ge 1 ] || { echo "[error] job count not set before writing the array script" >&2; exit 1; }
cat > "$array" <<EOF
#!/bin/bash
#SBATCH --job-name=wga_proj_${ARM}
#SBATCH --array=1-${n}%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --exclude=$EXCLUDE
#SBATCH --output=$RUNDIR/proj_%a.%j.out
source "$CONDA_SH"
conda activate phyloslide
set -euo pipefail

line=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$joblist")
IFS=\$'\t' read -r codename maf out <<< "\$line"

python3 "$RUNDIR/project_maf_to_anchor.py" \\
    --maf "\$maf" --anchor-fai "${ANCHOR_FA}.fai" --anchor-prefix "${ANCHOR}." \\
    --query-name "\$codename" --out "\$out"
/home/people/micwe/Software/miniconda3/envs/phyloslide/bin/samtools faidx "\$out"
EOF

echo
echo "[submit] projection array (${ARM}): $n species, 16G each"
if [ -z "$DRY_RUN" ]; then
    args=(); [ -n "$AFTER" ] && args=(--dependency=afterok:"$AFTER")
    sbatch --parsable "${args[@]}" "$array"
fi

# the anchor is its own projection
if [ ! -e "$OUTDIR/${ANCHOR}.fa" ] && [ -s "$ANCHOR_FA" ]; then
    ln -s "$ANCHOR_FA" "$OUTDIR/${ANCHOR}.fa"
    ln -s "${ANCHOR_FA}.fai" "$OUTDIR/${ANCHOR}.fa.fai"
    echo "[link]   $ANCHOR (anchor needs no projection)"
fi
echo "Run dir: $RUNDIR"
