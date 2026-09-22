# PhyloSlide

PhyloSlide is a sliding-window phylogenomics pipeline for extracting genomic windows from multiple individuals and (optionally) inferring per-window phylogenetic trees in parallel across windows. It includes filters for missing data and low-information windows, optional aDNA-oriented masking (transversions-only), concordance factor analysis, and topology-based window selection for downstream dating analyses.

---

## What PhyloSlide Does

### Extraction Stage (always runs)

- Reads a sample table (`CODENAME<TAB>genome.fasta`)
- Uses a regions list in samtools faidx coordinate format (`chr:start-end`, 1-based inclusive)
- For each sample and region:
  - extracts the window using `samtools faidx`
  - pads with Ns if the region is missing
  - writes window FASTAs with header format: `>CODENAME`
- Writes a per-sample concatenated sequence across all regions:
  - `OUT/<CODENAME>/<CODENAME>_concat.fasta`

### Optional Tree Pipeline (`--runtrees`)

- Filters windows by missing data (`--maxN`)
- Builds per-window multi-FASTA alignments
- Optional conservative aDNA mode (`--transversions`)
  - masks columns containing A/G or C/T variation to N
  - ensures all downstream alignments are built consistently from masked data
- Filters windows by minimum parsimony-informative sites (`--minpi`)
- Runs one IQ-TREE job per window in parallel (`--jobs`)
- Builds a reference tree:
  - concatenated IQ-TREE (`--ref concat`)
  - or ASTRAL species tree (`--ref astral`)
- Computes gCF and sCF using IQ-TREE2

### Optional Topology Filtering (`--topofilter`)

Builds a dating supermatrix using only windows that:
- passed missingness and minPI filters
- are within `--maxrf` Robinson-Foulds distance of the reference topology
- are clock-like: root-to-tip coefficient of variation ≤ `--maxcov`
- have **mean** internal bootstrap ≥ `--minbs`

Columns that are `N` in every taxon are dropped from the dating supermatrix
(they carry no information and upset PAML/baseml). The per-window filters are
per-sample and never inspect individual alignment columns, so this is the only
place columns are removed.

---

## Installation

### Recommended: Conda / Mamba Environment

Create a file named `environment.yml`:

```yaml
name: phyloslide
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - biopython>=1.80
  - samtools
  - seqtk
  - bedtools
  - iqtree
  - openjdk
  - astral-tree
```

Create and activate the environment:

```bash
mamba env create -f environment.yml
mamba activate phyloslide
```

If using conda instead of mamba:

```bash
conda env create -f environment.yml
conda activate phyloslide
```

Notes:
- `bedtools` is required only for `--makewindows`
- `astral-tree` and `openjdk` are required only for `--ref astral`
- `biopython` is required for `--ref astral` and `--topofilter`

---

## Input Format

### Sample Table (`--input`)

Tab-separated file:

```
American_black    /path/to/American_black.fa
Asian_black       /path/to/Asian_black.fa
Brown             /path/to/Brown.fa
Polar             /path/to/Polar.fa
```

- No spaces in codenames
- Codenames must be unique (duplicates are rejected)
- FASTA files must be indexed (`samtools faidx`) — PhyloSlide will create `.fai` if missing

### Regions File (`--regions`)

Each line must be:

```
chr:start-end
```

Coordinates must be 1-based inclusive.

Example:

```
chr1:1-20000
chr1:1000001-1020000
```

---

## Generating Windows Internally (`--makewindows`)

Instead of manually creating a regions file, PhyloSlide can generate one using `bedtools makewindows`.

Example: 20 kb windows sliding every 1 Mb:

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000
```

Restrict to long scaffolds:

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000 \
  --min_scaffold_len 14000000
```

Restrict to specific chromosomes:

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows \
  --refgenome ref.fa \
  --window 20000 \
  --step 1000000 \
  --chroms chroms.txt
```

`chroms.txt` should contain one contig name per line.

---

## Running the Tree Pipeline

Example full workflow:

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows --refgenome ref.fa --window 20000 --step 1000000 \
  --runtrees \
  --jobs 32 \
  --iqtree_threads_per_job 1
```

Parallel recommendations:
- Prefer `--iqtree_threads_per_job 1`
- Set `--jobs` equal to available cores (or slightly below)
- Reduce `--jobs` if running on slow network storage

---

## aDNA Conservative Mode (`--transversions`)

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --regions regions.txt \
  --outdir OUT \
  --runtrees \
  --transversions
```

Behavior:
- Masks entire columns that contain both A and G OR both C and T
- All downstream alignments remain consistent with TV-only filtering

---

## Filters

### Missing Data Filter (`--maxN`)
- Default: 0.5
- Drops window if ANY sample has fraction of N > maxN

### Minimum Parsimony-Informative Sites (`--minpi`)
- Default:
  - 20 (full data)
  - 10 (with `--transversions`)
- Applied after TV masking (if enabled)

---

## Reference Trees

### Concatenated (default)

```
--ref concat
```

Runs IQ-TREE on:
- `All_concat.filtered.fasta`
- or `All_concat.filtered.tv.fasta`

### ASTRAL

```
--ref astral --astral /path/to/astral.jar
```

- Collapses branches with support < `--collapse`
- Requires Java
- Requires biopython

---

## Concordance Factors

When `--runtrees` is enabled:
- gCF and sCF are computed using IQ-TREE2
- Outputs in `OUT/trees/concordance/`

---

## Topology Filtering

```bash
python3 PhyloSlide.py \
  --input samples.tsv \
  --regions regions.txt \
  --outdir OUT \
  --runtrees \
  --topofilter \
  --minbs 90 \
  --topomode exact
```

Options:

- `--maxrf` default **2** — maximum Robinson-Foulds distance between a window
  tree and the reference tree.
  - `0` = **exact match**: every bipartition identical to the reference.
  - `2` = one bipartition may differ, i.e. one NNI move from the reference.
  - RF counts differing bipartitions in *both* directions, so for two fully
    resolved trees over the same taxa it is always **even**. `--maxrf 1` behaves
    identically to `--maxrf 0`, and `--maxrf 3` identically to `--maxrf 2`.

  Exact matching is strict: on an 18-taxon dataset only ~13% of windows matched
  exactly, while ~42% were within one NNI move.

- `--maxcov` default **0.1** — maximum coefficient of variation in root-to-tip
  length ("non-clocklikeness") of the window tree, after midpoint rooting.
  A strict molecular clock gives 0; larger values mean more rate variation among
  lineages. Windows with erratic rates distort the relaxed-clock model used for
  dating. Set to a large number (e.g. `--maxcov 999`) to disable.

- `--minbs` default 90 — applied to the **mean** internal bootstrap of the
  window tree, averaged over all internal branches.
  (Before 2026-09 this was the *minimum* internal bootstrap. That required
  every one of the N-3 internal nodes to clear the threshold, which gets
  steadily harsher as taxa are added and discarded most windows on datasets
  with more than ~10 taxa.)

- `--topomode` — **deprecated**, kept for backwards compatibility.
  `--topomode exact` is equivalent to `--maxrf 0`; `--topomode compatible`
  ignores `--maxrf` and tests whether the window splits are a subset of the
  reference splits. Prefer `--maxrf`.

Defaults follow published practice for window-based dating (RF ≤ 2 and
root-to-tip CoV < 0.1).

Outputs (in `filtering/` and `Combined/`, tagged with the filter settings):
- `regions.topomatch.rf{R}.cov{C}.minbs{B}.{full|tv}.txt` — windows kept
- `regions.topofail.rf{R}.cov{C}.minbs{B}.{full|tv}.txt` — windows dropped, with reason
- `window_tree_stats.{full|tv}.tsv` — per-window `rf`, `rtt_cov`, `mean_bs` for
  **every** window, so thresholds can be re-chosen without rerunning anything
- `All_concat.topomatch.rf{R}.cov{C}.minbs{B}.{full|tv}.fasta` — dating supermatrix
- `All_concat.topomatch.rf{R}.cov{C}.minbs{B}.{full|tv}.phy` — same, PHYLIP
  (with `--dating_phylip`)
- `All_concat.topomatch.rf{R}.cov{C}.minbs{B}.{full|tv}.chrom{N}.phy` — PHYLIP with
  one block per chromosome (with `--dating_partition chrom`)

### PAML / MCMCtree output

```
--dating_phylip                  also write sequential PHYLIP
--dating_partition {none,chrom}  PHYLIP layout (default: none)
```

`--dating_phylip` writes the supermatrix in the sequential PHYLIP format read by
PAML's `baseml` and `mcmctree`, so no external conversion step is needed.

`--dating_partition` controls the layout and is implied by `chrom`:

- `none` — one alignment block. Use `ndata = 1` in the MCMCtree control file.
- `chrom` — one block per chromosome/scaffold, ordered by first appearance in the
  regions file. The number of blocks is written to the log; use it as `ndata`.

Partitioning is **opt-in on purpose**. Per-partition rate parameters noticeably
narrow the posterior on node ages, but that narrowing comes from the model rather
than from additional data, and the effect is largest on exactly the nodes that
carry no calibration. Prefer `none` unless the unpartitioned uncertainty is
genuinely unusable, and say which you used in your methods.

MCMCtree then runs in two passes (see the PAML documentation):

1. `usedata = 3` — `baseml` estimates branch lengths plus the gradient and
   Hessian, writing `out.BV`. Cost scales with the number of partitions.
2. `mv out.BV in.BV`, set `usedata = 2` — approximate-likelihood MCMC.

`mcmctree` invokes `baseml` through a shell call, so **`baseml` must be on
`PATH`**. If it is not, every locus fails with `file rst2 not found!`, `mcmctree`
still exits 0, and `out.BV` is written but invalid. Check the log before step 2.

---

## Output Structure

```
OUT/
  phyloslide.log
  filtering/
  Combined_windows/
  Combined/
  trees/
  <CODENAME>/
```

Key files:
- `regions.kept.final.txt`
- `All_concat.filtered.fasta`
- `all_window_trees.trs`
- Reference tree
- Concordance outputs

---

## Reruns, Resume, and Restart

PhyloSlide supports resumable runs and now includes parameter-safety checks.

- Default rerun behavior:
  - Reuses existing non-empty intermediates when possible.
  - Writes and checks `OUT/phyloslide.run_manifest.json`.
  - If parameters differ from a previous run in the same `--outdir`, execution stops to avoid mixed outputs.

- `--force_rebuild`:
  - Recomputes and overwrites intermediates in-place, even if files already exist.
  - Allows rerun despite manifest parameter differences.

- `--restart`:
  - Cleans prior PhyloSlide outputs in `--outdir` and reruns from scratch.
  - Implies rebuild behavior for that run.

Recommended:
- Use a new `--outdir` for each parameter set, OR
- Use `--restart` / `--force_rebuild` explicitly when changing parameters.

---

## Troubleshooting

If slow:
- Reduce `--jobs`
- Run on local scratch disk

If no windows left:
- Increase `--maxN`
- Decrease `--minpi`
- Increase window size
- Reduce `--minbs`

If IQ-TREE fails:
- Check `trees/window_logs/`
- Check `failed_windows.txt`

Windows filenames contain colons (`chr:start-end`).
Use Linux/macOS; Windows may fail due to colon in filenames.

---

## Recommended Student Workflow

1. Test on ~50 windows first.
2. Confirm filtering behavior.
3. Scale up gradually.
4. Record:
   - window parameters
   - filtering parameters
   - model and bootstrap settings
   - reference method used

---

## Dependencies Summary

Always required:
- samtools
- seqtk

If using:
- `--makewindows` → bedtools
- `--runtrees` → iqtree
- `--ref astral` → java + astral-tree
- `--topofilter` or `--ref astral` → biopython

---

## Citation

Please cite:
- samtools
- bedtools
- IQ-TREE2
- ASTRAL

PhyloSlide orchestrates these tools but does not replace them.
