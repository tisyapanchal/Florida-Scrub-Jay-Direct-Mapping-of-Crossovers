# FSJ Direct Recombination Mapping — Snakemake Pipeline

This is a Snakemake port of [`notebooks/scrubjay.ipynb`](../notebooks/scrubjay.ipynb),
generalized to run on any Florida Scrub Jay (FSJ) chromosome instead of only
chromosome 1. Every stage takes a chromosome label as a wildcard, and every
CSV/plot it produces is written under `results/chr<chrom>/` with `_chr<chrom>`
in the filename, so outputs from different chromosomes never collide or get
mixed up.

## What the pipeline does

Given a phased, pedigree-linked tree sequence for one chromosome (from
`shapeit5`/`tsinfer`-style phasing over the FSJ pedigree) plus the pedigree
itself, the pipeline:

1. **Builds a null model.** Constructs an `msprime` `PedigreeBuilder` tree
   sequence from the FSJ pedigree (`fsj-anon.ped`) at this chromosome's
   length, then simulates 100 replicate meioses over that fixed pedigree
   under a plain recombination process (no phasing, no biology beyond a flat
   recombination rate). This is the "what would breakpoint counts look like
   with no real signal, just this pedigree and chromosome length" baseline.
2. **Counts breakpoints in the real data.** For every individual, counts how
   many times each of their two haplotypes switches tree-parent across the
   chromosome — that's a recombination breakpoint.
3. **Cleans up phasing errors.** Raw breakpoint calls in the real data are
   inflated by phasing errors, not just true crossovers. Three filters strip
   these out:
   - **full-sibling recurrence** — a breakpoint on the same parental
     haplotype shared by ≥2 full siblings within a small window is far more
     likely a shared phasing artifact than 2+ independent crossovers at the
     same spot.
   - **opposite-direction switches** — two siblings switching parental
     haplotype in *opposite* directions at (nearly) the same position.
   - **half-sibling recurrence** — the same recurrence idea, but across all
     offspring of a single parent (catches what full-sib grouping misses).
4. **Filters by spacing.** After the sibling filters, breakpoints that still
   sit implausibly close to their nearest neighbour on the same haplotype are
   dropped, using a cutoff found automatically from that chromosome's own
   bimodal nearest-neighbour-distance distribution (the valley/antimode
   between "true spacing" and "clustered artifacts").
5. **Builds recombination maps** (direct mapping from the filtered real
   breakpoints; comparison against the Romero et al. (2024) linkage
   map for the same chromosome) and checks how they relate to gene
   annotation (genic/intergenic/exon/intron enrichment, distance to
   nearest TSS/TES vs. the null replicates) and to full-sib family size.

## Layout

```
workflow/
  Snakefile              # rule definitions, wired chrom -> chrom
  config.yaml             # chromosomes to run, input paths, filter parameters
  scripts/
    common.py             # shared logic ported from the notebook
    plot_ibd_tracts.py    # standalone IBD-tract plotter (unchanged from notebooks/)
    build_pedigree.py
    simulate_null.py
    null_summary.py
    real_breakpoints.py
    sibling_filters.py
    distance_filter.py
    density_plots.py
    marey_maps.py
    gene_annotation_overlap.py
    tss_tes_distance.py
    family_size_maps.py
    ibd_plots.py
results/
  chr<chrom>/              # one directory per chromosome, created on first run
```

## Running it

The pipeline needs [Snakemake](https://snakemake.readthedocs.io/) and the same
packages the notebook uses (`msprime`, `tskit`, `pandas`, `numpy`,
`matplotlib`, `networkx`, `pyyaml`) — the `scrubjay2` conda environment
already has all of these plus Snakemake installed.

```bash
conda activate scrubjay2

# from the repo root (or anywhere — the Snakefile pins itself to the repo root):
snakemake -s workflow/Snakefile --cores 4          # run everything in config.yaml
snakemake -s workflow/Snakefile --cores 4 -n        # dry run — see the plan without executing
snakemake -s workflow/Snakefile --cores 4 -- results/chr1/all_filters_breakpoints_chr1.csv
                                                     # build just one target (and its dependencies)
```

### Adding the other 35 chromosomes

Only chromosome 1's tree sequence
(`notebooks/1.shapeit5.hmm.v3.pedigree.trees`) is present in this repo. To
run another chromosome:

1. Place its phased tree sequence where `real_trees_pattern` in
   `config.yaml` expects it — by default
   `notebooks/<chrom>.shapeit5.hmm.v3.pedigree.trees` — using the same
   chromosome label used in the linkage-map `Scaffold` column and the GFF's
   `Chr<chrom>` naming (e.g. `"1"`, `"1A"`, `"4A"`, `"Z"`, `"10"`, ... — this
   is the FSJ karyotype labeling, not necessarily sequential integers).
2. Add that label to the `chromosomes:` list in `config.yaml`.
3. Re-run `snakemake`. Only the new chromosome's jobs will run — chromosome
   1's outputs already exist and won't be recomputed.

You can run several chromosomes in one invocation by listing them all in
`chromosomes:` and giving Snakemake enough `--cores` to parallelize across
them (each chromosome's rules are otherwise fully independent).

### IBD diagnostic plots

`ibd_plots` (IBD-tract figures for one focal individual, at three filtering
stages) is *not* run for all 2,602 individuals by default — that's a
per-individual diagnostic, not a batch analysis. It only runs for the IDs
listed in `config["ibd_plot_individuals"]` (default: `FSJ3046`, matching the
notebook's example). Add more IDs there, or point directly at
`results/chr<chrom>/ibd_plots/<id>_all_filters_ibd_chr<chrom>.png` to
generate one on demand.

## Pipeline stages and outputs

Every output below is written to `results/chr<chrom>/` and every filename
(and every CSV's `chrom` column) is tagged with the chromosome it came from.

| Rule | Outputs | What's in them |
|---|---|---|
| `build_pedigree` | `final_null_pedigree_chr<c>.trees` | msprime pedigree, finalised at this chromosome's length |
| `simulate_null` | `finalized_null_breakpoints_100sims_chr<c>.csv`, `null_trees/null_seed<1‑100>.trees` | per-individual breakpoint counts for each of the 100 null replicates, and each replicate's tree sequence |
| `null_summary` | `per_individual_variance_chr<c>.csv` | mean/var/std/min/max breakpoint count per individual, across null replicates |
| `real_breakpoints` | `real_id_table_chr<c>.csv`, `real_breakpoints_table_chr<c>.csv`, `qc_summary_chr<c>.txt` | pedigree-ID ↔ tree-sequence-ID lookup; raw (unfiltered) per-individual breakpoint counts; REAL vs. NULL crossovers-per-gamete sanity check |
| `sibling_filters` | `breakpoints_sibling_recurrence_chr<c>.csv`, `sibship_sharing_summary_chr<c>.csv`, `opposite_direction_switches_chr<c>.csv`, `sibling_filtered_bps_chr<c>.csv` | full-sibling recurrence + opposite-direction-switch + half-sibling recurrence filters applied; `sibling_filtered_bps` is the survivor set used downstream |
| `distance_filter` | `distances_btwn_sibling_filtered_breakpoints_chr<c>.csv`, `distance_btwn_breakpoint_distribution_chr<c>.png`, `nearest_neighbor_cutoff_chr<c>.txt`, `all_filters_breakpoints_chr<c>.csv`, `filtered_bps_chr<c>.csv` | nearest-neighbour spacing filter (cutoff auto-detected per chromosome); `all_filters_breakpoints` is the human-readable (Mb-rounded) final table, `filtered_bps` is its full-bp-precision companion used by every later stage |
| `density_plots` | `mat_vs_pat_chr<c>.png`, `mat_vs_pat_filtered_chr<c>.png` | crossover density along the chromosome, maternal vs. paternal, before/after the spacing filter |
| `marey_maps` | `dm_chr<c>.png`, `lm_map_chr<c>.png`, `dm_vs_lm_chr<c>.png`, `direct_map_cM_summary_chr<c>.csv` | direct-mapping Marey map (female/male), the Romero et al. (2024) linkage map for this chromosome, and their overlay |
| `gene_annotation_overlap` | `feature_enrichment_chr<c>.csv`, `breakpoints_by_region_chr<c>.png` | exon/CDS/gene enrichment of filtered breakpoints relative to genome proportion; weighted gene/intergenic/exon/intron proportions |
| `tss_tes_distance` | `bp_dist_to_tss_chr<c>.png`, `bp_dist_to_tes_chr<c>.png`, `tss_tes_distance_summary_chr<c>.csv` | real (pooled) median distance to nearest TSS/TES vs. the per-replicate null distribution |
| `family_size_maps` | `marey_by_family_size_chr<c>.png`, `breakpoints_{lt5,5to10,11to15,16to20,gt20}_chr<c>.csv` | recombination maps and breakpoint tables split by full-sib family size |
| `ibd_plots` (opt-in) | `ibd_plots/<id>_{kept,sibling_filtered,all_filters}_ibd_chr<c>.png` | IBD-tract diagnostic plots for one individual at three filtering stages |

## Notes on changes from the notebook

The notebook was written and tuned against chromosome 1 specifically, so a
few things had to be made chromosome-agnostic (or fixed) to generalize
correctly rather than silently reusing chromosome-1-specific numbers on other
chromosomes:

- **Chromosome length** is always read from the tree sequence
  (`ts.sequence_length`) rather than the hardcoded `120748055` (chr1's
  length) used in several notebook cells.
- **The nearest-neighbour spacing cutoff** (chromosome 1 resolved to
  ~258 kb) is recomputed per chromosome from that chromosome's own
  breakpoint-spacing distribution, instead of reusing chromosome 1's value.
  The configured `default_nearest_neighbor_cutoff_bp` is only a fallback for
  the rare case a chromosome has no clear bimodal antimode to detect.
- **Null replicate tree sequences are each kept** (`null_trees/null_seed<N>.trees`
  per seed) rather than all 100 seeds overwriting a single `output.trees`
  path, since the TSS/TES-distance comparison needs every replicate
  independently.
- **`per_individual_variance`** now reads the null-breakpoints CSV this
  pipeline just generated for the chromosome in question (the notebook cell
  read a stale filename left over from an earlier iteration of the
  notebook).
- Every CSV gains a leading `chrom` column, in addition to the `_chr<c>`
  filename suffix, so per-chromosome files can be concatenated later without
  losing track of which chromosome each row came from.
