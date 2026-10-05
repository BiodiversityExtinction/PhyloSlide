#!/bin/bash
# Submits the LAST whole-genome-alignment arm to SLURM as three dependent stages:
#
#   1. prepare the anchor  : filter Panda to scaffolds >= 1 Mb, build 2bit/sizes + lastdb -c
#   2. align chunks (array): one task per chromosome (>= 10 Mb), plus binned tasks for the
#                            1-10 Mb remainder, each aligned against the anchor index at -P 1
#   3. finalize (array)    : per species, concatenate chunk MAFs then chain/net/syntenic-filter
#                            into results/wga/pairwise/<species>.syntenic.maf  (MULTIZ input)
#
# Why chunked: memory tracks query scaffold size (a measured 126 Mb scaffold peaked at 65 GB at
# -P 1), so whole-genome jobs OOM-killed repeatedly even at 150 GB. Chunking bounds memory
# structurally and turns ~1000 sequential CPU-hours into a few hundred parallel tasks.
#
# Why >= 1 Mb: keeps 92-99% of every genome while dropping ~95% of the sequence count. A 10 Mb
# cutoff was considered and rejected -- AmBlack is fragmented enough that it would lose 28%.
#
# Worker scripts and job lists are SNAPSHOTTED into logs/wga/run_<timestamp>/ at submit time so
# later edits to these scripts cannot yank the file out from under an in-flight job.
#
# Usage: ./align_to_anchor.sh --projdir <analysis dir> [--dry-run] [--concurrent-big N] [--concurrent-small N]
set -euo pipefail

# Archived from the reference-bias benchmark. PROJDIR is the ANALYSIS directory (the one holding
# data/ and results/), which is not this repository: set it with --projdir or PHYLOSLIDE_PROJDIR.
PROJDIR="${PHYLOSLIDE_PROJDIR:-}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANCHOR="Panda"
MIN_SCAF=1000000        # keep scaffolds >= 1 Mb
BIG_SCAF=10000000       # >= 10 Mb gets its own task (these are the chromosomes)
BIN_TARGET=50000000     # bin the 1-10 Mb remainder into ~50 Mb tasks

BIG_MEM=250G            # 120G OOM-killed 50 of 303 chunks (90-167 Mb chromosomes). Both
                        # lastal AND last-split are heavy here -- last-split holds every
                        # alignment for a query sequence and reported "out of memory" itself.
SMALL_MEM=48G
CONCURRENT_BIG=8
CONCURRENT_SMALL=8
INDEX_THREADS=8

# AVX2 is required by the bioconda-built LAST/UCSC-kent binaries; these nodes only have AVX and
# SIGILL on them. node08 also excluded: 124G RAM, under the big-task request.
# Big chunks now need 250G, so restrict to nodes that actually have it (compute02/03/05/06
# 0.5-1TB, node07 480G); node08 (124G) and node11 (250G) cannot host the request.
# 500G for the two largest chromosomes (124/167 Mb) where last-split exceeded 250G;
# only the 1TB compute nodes can host that (compute04 is reserved for the GPU seminar).
# 250G suits ~25 of 27 chromosome chunks; the two or three largest (>150 Mb) exceed it and
# get rerun at 500G afterwards. Requesting 500G up front confined everything to the three
# 1TB nodes and dropped concurrency to 2.
EXCLUDE="compute04,node01,node02,node03,node04,node05,node06,node08,node09,node10,node11,node12,node13"
CONDA_SH="/home/people/micwe/Software/miniconda3/etc/profile.d/conda.sh"
DRY_RUN=""

while [ $# -gt 0 ]; do
    case "$1" in
        --projdir)  PROJDIR="$2"; shift 2 ;;
        --projdir) PROJDIR="$2"; shift 2 ;;
        --dry-run)          DRY_RUN="yes"; shift ;;
        --concurrent-big)   CONCURRENT_BIG="$2"; shift 2 ;;
        --concurrent-small) CONCURRENT_SMALL="$2"; shift 2 ;;
        *) echo "[error] unknown argument: $1" >&2; exit 1 ;;
    esac
done
[ -n "$PROJDIR" ] || { echo "[error] set --projdir or PHYLOSLIDE_PROJDIR (the analysis directory holding data/ and results/)" >&2; exit 1; }

GENOMES="$PROJDIR/data/genomes"
OUTDIR="$PROJDIR/results/wga/pairwise"
CHUNKROOT="$PROJDIR/results/wga/chunks"
META="$PROJDIR/data/metadata/genomes.tsv"

RUNDIR="$PROJDIR/logs/wga/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTDIR" "$CHUNKROOT" "$RUNDIR"
cp "$SCRIPTDIR/../common/prepare_anchor.sh" \
   "$SCRIPTDIR/align_one_chunk.sh" \
   "$SCRIPTDIR/finalize_species_wga.sh" \
   "$SCRIPTDIR/filter_maf_to_anchor.py" "$RUNDIR/"

ANCHOR_DB="$OUTDIR/${ANCHOR}.lastdb"
ANCHOR_2BIT="$OUTDIR/${ANCHOR}.2bit"
ANCHOR_SIZES="$OUTDIR/${ANCHOR}.sizes"

# ---------- stage 1: anchor ----------
prep="$RUNDIR/prepare_anchor.sbatch.sh"
cat > "$prep" <<EOF
#!/bin/bash
#SBATCH --job-name=wga_anchor
#SBATCH --cpus-per-task=$INDEX_THREADS
#SBATCH --mem=64G
#SBATCH --exclude=$EXCLUDE
#SBATCH --output=$RUNDIR/prepare_anchor.%j.out
source "$CONDA_SH"
conda activate phyloslide_wga
set -euo pipefail

bash "$RUNDIR/prepare_anchor.sh" "$GENOMES/${ANCHOR}.fa" $MIN_SCAF "$OUTDIR" "$ANCHOR" $INDEX_THREADS
EOF

echo "[stage 1] anchor prep (filter >= $((MIN_SCAF/1000000)) Mb + lastdb -c)"
prep_id=""
if [ -z "$DRY_RUN" ]; then
    prep_id=$(sbatch --parsable "$prep")
    echo "          -> job $prep_id"
fi

# ---------- build chunk lists from the existing .fai (no need to wait on stage 1) ----------
joblist_big="$RUNDIR/joblist_big.tsv"
joblist_small="$RUNDIR/joblist_small.tsv"
: > "$joblist_big"; : > "$joblist_small"

species_list="$RUNDIR/species.txt"
: > "$species_list"

tail -n +2 "$META" | while IFS=$'\t' read -r codename accession common sci source role local_path; do
    [ -z "$codename" ] && continue
    [ "$codename" = "$ANCHOR" ] && continue
    [ -s "$OUTDIR/${codename}.syntenic.maf" ] && { echo "[done]  $codename already aligned"; continue; }
    echo "$codename" >> "$species_list"

    fai="$GENOMES/${codename}.fa.fai"
    cdir="$CHUNKROOT/$codename"
    mkdir -p "$cdir"

    python3 - "$fai" "$cdir" "$MIN_SCAF" "$BIG_SCAF" "$BIN_TARGET" "$codename" \
              "$GENOMES/${codename}.fa" "$joblist_big" "$joblist_small" <<'PY'
import sys, os
fai, cdir, min_scaf, big_scaf, bin_target, codename, query_fa, jb_big, jb_small = sys.argv[1:]
min_scaf, big_scaf, bin_target = int(min_scaf), int(big_scaf), int(bin_target)

scaf = []
with open(fai) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) >= 2 and int(p[1]) >= min_scaf:
            scaf.append((p[0], int(p[1])))
scaf.sort(key=lambda x: -x[1])

big = [s for s in scaf if s[1] >= big_scaf]
small = [s for s in scaf if s[1] < big_scaf]

def write(chunk_id, names, joblist):
    path = os.path.join(cdir, f"{chunk_id}.scaffolds")
    with open(path, "w") as fh:
        fh.write("\n".join(names) + "\n")
    out_maf = os.path.join(cdir, f"{chunk_id}.maf")
    # skip chunks already aligned, so a rerun after failures queues only what is missing
    if os.path.exists(out_maf) and os.path.getsize(out_maf) > 0:
        return
    with open(joblist, "a") as fh:
        fh.write(f"{codename}\t{path}\t{query_fa}\t{out_maf}\n")

for i, (name, _) in enumerate(big):
    write(f"big{i:04d}", [name], jb_big)

bin_names, bin_bp, idx = [], 0, 0
for name, ln in small:
    bin_names.append(name); bin_bp += ln
    if bin_bp >= bin_target:
        write(f"small{idx:04d}", bin_names, jb_small); idx += 1
        bin_names, bin_bp = [], 0
if bin_names:
    write(f"small{idx:04d}", bin_names, jb_small)

total = sum(s[1] for s in scaf)
print(f"[queue] {codename}: {len(big)} chromosome task(s) + {idx + (1 if bin_names else 0)} "
      f"binned task(s), {total/1e9:.2f} Gb retained of {len(scaf)} scaffolds >= {min_scaf//10**6} Mb")
PY
done

n_big=$(wc -l < "$joblist_big")
n_small=$(wc -l < "$joblist_small")
n_species=$(wc -l < "$species_list")
echo
echo "[stage 2] chunk alignment: $n_big chromosome task(s) @ $BIG_MEM, $n_small binned task(s) @ $SMALL_MEM"

# Do NOT exit here when there is nothing left to align: once every chunk is done, finalize is
# precisely the stage that still needs submitting.
if [ "$n_big" -eq 0 ] && [ "$n_small" -eq 0 ]; then
    echo "  all chunks already aligned; going straight to finalize"
fi

submit_chunk_array() {
    local name="$1" list="$2" n="$3" mem="$4" conc="$5" dep="$6"
    local script="$RUNDIR/${name}_array.sbatch.sh"
    cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=wga_$name
#SBATCH --cpus-per-task=1
#SBATCH --mem=$mem
#SBATCH --array=1-${n}%${conc}
#SBATCH --exclude=$EXCLUDE
#SBATCH --output=$RUNDIR/${name}_%a.%j.out
source "$CONDA_SH"
conda activate phyloslide_wga
set -euo pipefail

line=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$list")
IFS=\$'\t' read -r codename scaffold_list query_fa out_maf <<< "\$line"

bash "$RUNDIR/align_one_chunk.sh" "$ANCHOR_DB" "\$query_fa" "\$scaffold_list" "\$out_maf"
EOF
    if [ -z "$DRY_RUN" ]; then
        local args=()
        [ -n "$dep" ] && args=(--dependency=afterok:"$dep")
        sbatch --parsable "${args[@]}" "$script"
    fi
}

big_id=""; small_id=""
[ "$n_big" -gt 0 ]   && big_id=$(submit_chunk_array big "$joblist_big" "$n_big" "$BIG_MEM" "$CONCURRENT_BIG" "$prep_id")
[ "$n_small" -gt 0 ] && small_id=$(submit_chunk_array small "$joblist_small" "$n_small" "$SMALL_MEM" "$CONCURRENT_SMALL" "$prep_id")
[ -n "$big_id" ]   && echo "          -> big array   job $big_id"
[ -n "$small_id" ] && echo "          -> small array job $small_id"

# ---------- stage 3: per-species finalize ----------
echo "[stage 3] finalize: $n_species species (concat chunks + chain/net)"
fin="$RUNDIR/finalize_array.sbatch.sh"
cat > "$fin" <<EOF
#!/bin/bash
#SBATCH --job-name=wga_final
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=1-${n_species}%4
#SBATCH --exclude=$EXCLUDE
#SBATCH --output=$RUNDIR/finalize_%a.%j.out
source "$CONDA_SH"
conda activate phyloslide_wga
set -euo pipefail

codename=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$species_list")
bash "$RUNDIR/finalize_species_wga.sh" "$ANCHOR" "$ANCHOR_2BIT" "$ANCHOR_SIZES" \\
    "\$codename" "$GENOMES/\${codename}.fa" "$CHUNKROOT/\$codename" "$OUTDIR/\${codename}.syntenic.maf"
EOF

if [ -z "$DRY_RUN" ]; then
    deps=""
    [ -n "$big_id" ]   && deps="afterok:$big_id"
    [ -n "$small_id" ] && deps="${deps:+$deps,}afterok:$small_id"
    args=(); [ -n "$deps" ] && args=(--dependency="$deps")
    fin_id=$(sbatch --parsable "${args[@]}" "$fin")
    echo "          -> job $fin_id"
fi
echo
echo "Run dir: $RUNDIR"
