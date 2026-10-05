# Whole-genome alignment pipelines

Scripts used to build the whole-genome-alignment arms of the reference-bias benchmark. Each
assembly is aligned to a single anchor genome, reduced to one-to-one syntenic blocks, and then
projected back onto anchor coordinates to give one pseudo-genome FASTA per taxon. Those
pseudo-genomes are drop-in replacements for reference-mapped consensus sequences, so the
alignment and mapping arms can be pushed through an identical downstream analysis on identical
window coordinates.

Two aligners are implemented independently so that any disagreement with the mapped arms can be
attributed to the mapping step rather than to one aligner's behaviour.

```
WGA/
  common/     anchor preparation and projection, shared by both arms
  LAST/       arm A -- lastal + last-split
  minimap2/   arm B -- minimap2 asm10 + paf2chain
```

## Running

Run order, for each arm:

```bash
# 1. align every assembly to the anchor (submits a SLURM array)
LAST/align_to_anchor.sh      --projdir /path/to/analysis
minimap2/align_to_anchor.sh  --projdir /path/to/analysis

# 2. project the syntenic MAFs onto anchor coordinates
common/project_to_anchor.sh  --projdir /path/to/analysis --arm last
common/project_to_anchor.sh  --projdir /path/to/analysis --arm minimap2
```

`--projdir` is the analysis directory holding `data/` and `results/`, not this repository; it can
also be given as `PHYLOSLIDE_PROJDIR`. Every driver takes `--dry-run`, which prints the work it
would submit without submitting it. `LAST/align_to_anchor.sh` also prepares the anchor on first
run (`common/prepare_anchor.sh`).

Inputs expected under `--projdir`:

| path | contents |
| --- | --- |
| `data/genomes/<codename>.fa` | one soft-masked FASTA per taxon, indexed |
| `data/metadata/genomes.tsv` | `codename`, `role` (the anchor is marked in `role`) |

Outputs: `results/wga/pairwise/` and `results/wga/pairwise_minimap2/` (one `<codename>.syntenic.maf`
per taxon), then `results/wga/projected_last/` and `results/wga/projected_minimap2/` (one
`<codename>.fa` per taxon, all in anchor coordinates and of identical length).

## How it works

**Anchor.** `common/prepare_anchor.sh` filters the anchor assembly to scaffolds above a length
cutoff and builds the 2bit, sizes and `lastdb` indexes. `lastdb -c` is what makes LAST act on the
soft-masking already present in the FASTA; without it the masking is read but never applied, and
the alignment becomes intractable.

**Arm A, LAST.** `lastal -P1 -m 50 -E 0.05 -u 2 | last-split`. Queries are split per chromosome
(large scaffolds individually, small ones binned) and aligned independently before merging, which
bounds peak memory. Single-threaded alignment is deliberate: `-P8` raised peak memory from 65 GB
to 475 GB with no reduction in wall time.

**Arm B, minimap2.** `minimap2 -cx asm10 --cs | paf2chain`, then `rescore_chains.py`. The rescoring
is required, not optional. `paf2chain` copies the PAF mapping quality into the chain score, so
every chain scores 255 and all are discarded by `chainNet`, whose default `-minScore` is 2000.
The result is an empty net and a near-empty MAF **produced with exit status 0** — a silent failure
that yields plausible-looking output. `rescore_chains.py` rescores each chain by its total aligned
bases and renumbers the chain IDs.

**Both arms** then run an identical UCSC chain/net pipeline to extract one-to-one syntenic blocks:
`axtChain -linearGap=medium` → `chainSort` → `chainNet` → `netSyntenic` → `netToAxt` → `axtToMaf`.

**Projection.** `common/project_maf_to_anchor.py` writes each taxon into anchor coordinates:
positions deleted relative to the anchor become `N` (matching the convention ANGSD uses for
uncovered positions), insertions relative to the anchor are dropped, and anchor positions covered
by no alignment block become `N`. This replaces the multiple-alignment stage (MULTIZ) used in
comparable pipelines: once every taxon is in anchor coordinates the taxa are already mutually
aligned column for column, and a progressive multi-way merge adds nothing.

## Dependencies

| tool | version used |
| --- | --- |
| LAST | 1654 |
| minimap2 | 2.31 |
| paf2chain | 0.1.1 |
| UCSC genome browser tools | 482 |
| samtools | 1.13 |
| Python | 3 (standard library only) |

## Portability

These were written for one SLURM cluster and are archived as run, not generalised. Adapting them
elsewhere means reviewing at least the following, all of which are near the top of each driver:

- `#SBATCH` partitions, memory and `--exclude` node lists. The exclusions are not arbitrary:
  bioconda builds of LAST and the UCSC tools use AVX2, and nodes without it fail with
  `Illegal instruction` rather than a clean error.
- Absolute paths to conda and to tools (`samtools` in particular is pinned by path).
- Memory ceilings. `last-split` is the peak consumer; chunks above roughly 150 Mb needed 500 GB.
- Workers are snapshotted into a timestamped run directory at submit time, so that editing a
  script cannot disturb a job already executing it over NFS.
