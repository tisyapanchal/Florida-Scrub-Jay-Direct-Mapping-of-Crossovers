"""Rule: tss_tes_distance — is the real breakpoint set closer to gene
transcription start/end sites than expected under the null model? Computes
the real (pooled) median distance to nearest TSS/TES and compares it against
the per-replicate median from each of the null simulations."""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

from common import dist_to, informative_nodes, load_genes, pooled_breakpoints

chrom = snakemake.wildcards.chrom
gff_chrom_name = snakemake.params.gff_chrom_name

genes = load_genes(snakemake.input.gff, gff_chrom_name)
tss = np.where(genes.strand == "+", genes.start.values, genes.end.values)
tes = np.where(genes.strand == "+", genes.end.values, genes.start.values)
tss_sorted = np.sort(tss)
tes_sorted = np.sort(tes)

filtered_table = pd.read_csv(snakemake.input.filtered)
real_bp = np.sort(filtered_table.position_bp.to_numpy())
real_median = {site: np.median(dist_to(real_bp, sites))
                for site, sites in [("TSS", tss_sorted), ("TES", tes_sorted)]}

null_files = sorted(glob.glob(os.path.join(snakemake.input.null_trees_dir, "null_seed*.trees")))
print(f"[chr{chrom}] found {len(null_files)} null replicates")

null_medians = {"TSS": [], "TES": []}
for path in null_files:
    seed_ts = tskit.load(path)
    seed_bp = pooled_breakpoints(informative_nodes(seed_ts))
    if seed_bp.size == 0:
        continue
    null_medians["TSS"].append(np.median(dist_to(seed_bp, tss_sorted)))
    null_medians["TES"].append(np.median(dist_to(seed_bp, tes_sorted)))
null_medians = {site: np.asarray(v) for site, v in null_medians.items()}

summary_rows = []
for site in ["TSS", "TES"]:
    n = len(null_medians[site])
    below = int((null_medians[site] >= real_median[site]).sum()) if n else None
    print(f"[chr{chrom}] real -> {site}: {real_median[site]:,.0f} bp")
    summary_rows.append({
        "chrom": chrom, "site": site,
        "real_median_bp": real_median[site],
        "n_null_replicates": n,
        "null_mean_median_bp": null_medians[site].mean() if n else float("nan"),
        "n_null_more_extreme_or_equal": below,
    })
pd.DataFrame(summary_rows).to_csv(snakemake.output.summary_csv, index=False)

plt.rcParams["font.size"] = 20
GREY, GREEN = "#5B595D", "#5903b5"
X_STEP, N_BINS, Y_STEP = 3000, 20, 5

outputs = {"TSS": snakemake.output.tss_png, "TES": snakemake.output.tes_png}
for site in ["TSS", "TES"]:
    vals = null_medians[site]
    real_val = real_median[site]
    out_path = outputs[site]

    if len(vals) == 0:
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.text(0.5, 0.5, "no null replicates available", ha="center", va="center")
        plt.savefig(out_path, dpi=150)
        plt.close(fig)
        continue

    n = len(vals)
    x_min = np.floor(min(vals.min(), real_val) / X_STEP) * X_STEP
    x_max = np.ceil(max(vals.max(), real_val) / X_STEP) * X_STEP
    if x_max == x_min:
        x_max += X_STEP
    x_ticks = np.arange(x_min, x_max + 1, X_STEP)
    bin_edges = np.linspace(x_min, x_max, N_BINS + 1)

    fig, ax = plt.subplots(figsize=(9, 7))
    counts, bin_edges, patches = ax.hist(
        vals, bins=bin_edges, color=GREY, alpha=0.8, edgecolor="#000000", linewidth=0.4,
        label=f"null (n={n}, mean {vals.mean():,.0f} bp)",
    )
    y_max = max(Y_STEP, np.ceil(counts.max() / Y_STEP) * Y_STEP)
    ax.vlines(real_val, ymin=0, ymax=y_max, color=GREEN, lw=3)
    ax.text(real_val, y_max + 0.18, f"Real Median: {real_val / 1e3:,.0f} kb", ha="center", va="bottom", fontsize=18)
    ax.set_xlim(x_ticks[0], x_ticks[-1])
    ax.set_xticks(x_ticks)
    ax.set_xlabel(f"Median Distance to Nearest {site} (bp)", fontsize=22)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Frequency", fontsize=22)
    ax.set_title(f"chr{chrom}")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
