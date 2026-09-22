# PhyloSlide

PhyloSlide is a sliding-window phylogenomics pipeline for extracting genomic
windows from multiple individuals and (optionally) inferring per-window
phylogenetic trees in parallel across windows. It includes filters for missing
data and low-information windows, optional aDNA-oriented masking
(transversions-only), concordance factor analysis, and topology- and
clock-based window selection for downstream dating analyses.

It can also prepare a complete MCMCtree starter kit — alignment, rooted tree,
node key and control files — so the only thing left to supply is your own fossil
calibrations.

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

Optionally it also writes:
- PHYLIP for PAML/MCMCtree (`--dating_phylip`), partitioned per chromosome if
  wanted (`--dating_partition chrom`)
- an MCMCtree starter kit (`--dating_template`): node key, rooted template tree,
  a calibration file to fill in, and a control file with `ndata` pre-filled

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

Other `--makewindows` options:

| Option | Default | Meaning |
|---|---|---|
| `--fai` | — | use an existing `.fai` index instead of `--refgenome` |
| `--regions_out` | `OUT/regions.generated.txt` | where to write the generated regions file |

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

### Windows to dated tree, end to end

```bash
# 1. windows, window trees, species tree, concordance factors,
#    filtered dating alignment, and an MCMCtree starter kit
python3 PhyloSlide.py \
  --input samples.tsv \
  --outdir OUT \
  --makewindows --refgenome ref.fa --window 50000 --step 1000000 \
  --runtrees \
  --ref astral --astral /path/to/astral.jar \
  --topofilter --maxrf 2 --maxcov 0.1 --minbs 90 \
  --dating_phylip \
  --dating_template --dating_outgroup Bos_taurus \
  --jobs 32 --iqtree_threads_per_job 1

# 2. read OUT/Combined/dating_template/node_key.txt, then add your fossils to
#    OUT/Combined/dating_template/calibrations.txt, e.g.
#       N1   B(0.172, 0.195, 1e-300, 0.025)

# 3. turn that into MCMCtree inputs
python3 prepare_mcmctree.py OUT/Combined/dating_template

# 4. run MCMCtree yourself (baseml must be on PATH)
cd OUT/Combined/dating_template
cp mcmctree_step1_hessian.ctl mcmctree.ctl && mcmctree mcmctree.ctl
grep -c 'rst2 not found' *.log          # must be 0
mv out.BV in.BV
cp mcmctree_step2_mcmc.ctl mcmctree.ctl && mcmctree mcmctree.ctl
```

PhyloSlide stops at step 3. Steps 2 and 4 are yours: the calibrations because
they are specific to your fossils, and the MCMCtree runs because their cost and
chain length depend on your data.

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

## IQ-TREE Parameters

| Option | Default | Applies to |
|---|---|---|
| `--model` | `GTR+R6` | window trees and the concatenated reference tree |
| `--bootstrap` | `1000` | IQ-TREE ultrafast bootstrap replicates (`-bb`) |
| `--scf` | `100` | site concordance factor replicates (`--scf`) |

`--model` is passed straight to IQ-TREE, so any model string IQ-TREE accepts will
work. Note that window trees and the concatenated reference use the same model.

---

## Disk and File Management

By default PhyloSlide deletes each per-sample window FASTA
(`OUT/<CODENAME>/<region>.fasta`) once the corresponding multi-FASTA alignment in
`Combined_windows/` has been built, since they are redundant at that point.

| Option | Effect |
|---|---|
| `--no_delete_windows` | keep the per-sample window FASTAs |
| `--archive` | on success, tar.gz `Combined_windows/`, `trees/window_trees/` and `trees/window_logs/` into `OUT/phyloslide_intermediates.tar.gz` and delete the originals |

`--no_delete_windows` is worth using if you expect to rerun with different
`--maxN`: the missingness filter reads those per-sample files, so without them a
rerun has to re-extract from the genomes. `--archive` keeps `Combined/`,
`filtering/`, `trees/reference/`, `trees/concordance/` and
`trees/all_window_trees.trs` unarchived.

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
    resolved trees over the same taxa it is always **even**. Odd values are
    rejected at startup, since they would silently behave like the next value
    down.

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

### MCMCtree starter kit

```
--dating_template --dating_outgroup CODENAME
```

Fossil calibrations are the one part of dating PhyloSlide cannot supply: they
depend on your taxa, your fossils, and a palaeontological judgement about which
node each fossil diagnoses. What PhyloSlide *can* do is remove the fiddly
mechanics around them. `--dating_template` writes:

```
Combined/dating_template/
    node_key.txt        every internal node, with its taxa spelled out
    tree.template.nwk   rooted tree with @N1@.. placeholders
    calibrations.txt    for you to fill in
    mcmctree.ctl        control file with ndata already matching the PHYLIP
```

Open `node_key.txt`, find the node you have a fossil for, and add one line per
calibration to `calibrations.txt`:

```
N1    B(0.172, 0.195, 1e-300, 0.025)
N13   L(0.03)
```

Then run the companion script, which substitutes them into the tree, strips the
unused placeholders, and writes both control files:

```
python3 prepare_mcmctree.py Combined/dating_template
```

It echoes back the taxa in every clade you calibrated, because **node numbers are
specific to one tree** and mean something different in any other analysis. Check
that echo before running anything.

Calibration syntax (times in units of **100 Myr**, so 17.2 Ma is `0.172`):

| Form | Meaning |
|---|---|
| `L(lo)` | minimum age only — the usual choice for a fossil |
| `U(hi)` | maximum age only |
| `B(lo, hi)` | bounded both sides, soft tails |
| `B(lo, hi, pL, pU)` | bounded, explicit tail probabilities; `pL=1e-300` = hard minimum |
| `G(alpha, beta)` | gamma prior |

You need at least one calibration, and the root must be constrained either by a
calibration or by `RootAge` in the control file. Put each fossil on the node it
actually diagnoses: a **stem** fossil of a group dates that group's split from
its sister lineage, not the group's crown node. Getting that wrong shifts every
date in the tree.

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
  phyloslide.run_manifest.json
  <CODENAME>/                      per-sample window FASTAs + <CODENAME>_concat.fasta
  Combined_windows/                <region>.fa  per-window multi-FASTA alignments
  filtering/
    regions.kept.maxN.txt
    regions.dropped.maxN.txt
    regions.kept.final.txt                     after --minpi
    regions.dropped.minpi.txt
    window_tree_stats.{full|tv}.tsv            rf, rtt_cov, mean_bs for EVERY window
    regions.topomatch.<tag>.txt                windows kept for dating
    regions.topofail.<tag>.txt                 windows dropped, with the reason
  Combined/
    All_concat.filtered.fasta                  all windows passing maxN + minpi
    All_concat.topomatch.<tag>.fasta           dating supermatrix
    All_concat.topomatch.<tag>.phy             --dating_phylip
    All_concat.topomatch.<tag>.chrom<N>.phy    --dating_partition chrom
    dating_template/                           --dating_template
      node_key.txt
      tree.template.nwk
      calibrations.txt
      mcmctree.ctl
  trees/
    window_trees/                  per-window IQ-TREE output
    window_logs/
    all_window_trees.trs           gene trees, in kept order
    reference/                     concatenated or ASTRAL reference tree
    concordance/                   gCF / sCF / discordance factors
```

`<tag>` encodes the filter settings, e.g. `rf2.cov0.1.minbs90.full`, so runs with
different thresholds do not overwrite each other.

`window_tree_stats.*.tsv` is worth knowing about: it records `rf`, `rtt_cov` and
`mean_bs` for **every** window that reached the topology filter, not just the ones
kept. Thresholds can therefore be re-chosen from the table without rerunning the
pipeline.

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

If no windows left after `--maxN` / `--minpi`:
- Increase `--maxN`
- Decrease `--minpi`
- Increase window size

If few or no windows survive `--topofilter`:
- Check `filtering/regions.topofail.<tag>.txt` — it gives the reason per window
- Consult `filtering/window_tree_stats.<suffix>.tsv` and pick thresholds from the
  actual distribution rather than guessing
- `--maxrf 0` (exact topology match) is strict; `--maxrf 2` typically keeps
  several times more windows
- `--maxcov` is often the binding filter. If your dataset contains a long-branch
  or high-missingness taxon, the whole root-to-tip CoV distribution shifts up
- `--minbs` applies to the **mean** internal bootstrap; windows matching the
  reference topology nearly always pass it, so it is rarely the culprit

`--maxrf must be an even number`:
- Robinson-Foulds distance counts differing bipartitions in both directions, so
  it is always even between two fully resolved trees. Use `0` or `2`.

MCMCtree reports `file rst2 not found!` for every locus:
- `mcmctree` calls `baseml` through the shell, so **`baseml` must be on `PATH`**.
  When it is not, `mcmctree` still exits 0 and writes an `out.BV` that is
  silently invalid. Always check the log before moving to `usedata=2`.

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
   - filtering parameters (`--maxN`, `--minpi`)
   - dating filters (`--maxrf`, `--maxcov`, `--minbs`)
   - model and bootstrap settings
   - reference method used (`--ref concat` or `--ref astral`)
   - whether the dating alignment was partitioned, and the calibrations used

---

## Dependencies Summary

Always required:
- samtools
- seqtk

If using:
- `--makewindows` → bedtools
- `--runtrees` → iqtree
- `--ref astral` → java + astral-tree
- `--topofilter`, `--dating_template` or `--ref astral` → biopython

`prepare_mcmctree.py` uses only the Python standard library.

Running the dating itself needs **PAML** (`mcmctree` and `baseml`), which
PhyloSlide does not call — it only prepares the inputs.

---

## Citation

Please cite:
- samtools
- bedtools
- IQ-TREE2
- ASTRAL
- PAML / MCMCtree, if you use the dating outputs

PhyloSlide orchestrates these tools but does not replace them.
