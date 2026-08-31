"""Rule: distance_filter — final filtering stage: drop breakpoints that sit
implausibly close to their nearest neighbour on the same haplotype.

The cutoff is *not* hardcoded — it's the antimode (valley) of this
chromosome's own bimodal nearest-neighbour-distance distribution, found the
same way the original notebook found chromosome 1's ~258 kb cutoff. Each
chromosome therefore gets its own data-driven cutoff instead of reusing
chromosome 1's value, which is written out for inspection alongside the plot.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit
from matplotlib.ticker import FuncFormatter

from common import breakpoint_distance_table, find_antimode_cutoff, nearest_neighbor_filter

chrom = snakemake.wildcards.chrom
default_min_dist = float(snakemake.params.default_min_dist)

real = tskit.load(snakemake.input.real_trees)
chrom_length = int(real.sequence_length)

final_bps_table = pd.read_csv(snakemake.input.final_bps)
final_bps = {node: g.position_bp.to_numpy() for node, g in final_bps_table.groupby("node")}
nodes_meta = {node: (g.fsj_id.iloc[0],) for node, g in final_bps_table.groupby("node")}

bp_table = breakpoint_distance_table(final_bps, nodes_meta, chrom_length=chrom_length)
bp_table.insert(0, "chrom", chrom)
bp_table.to_csv(snakemake.output.dist_table, index=False)
print(f"[chr{chrom}] {len(bp_table)} breakpoints across {bp_table.fsj_id.nunique()} individuals")

cutoff = find_antimode_cutoff(bp_table["nearest_bp"].to_numpy(), hi_bp=0.10 * chrom_length)
used_default = cutoff is None
if used_default:
    cutoff = max(1e4, 0.002 * chrom_length)
else:
    print(f"[chr{chrom}] antimode threshold: {cutoff / 1e3:.1f} kb")

with open(snakemake.output.cutoff_txt, "w") as f:
    f.write(f"chrom: {chrom}\n")
    f.write(f"nearest_neighbor_cutoff_bp: {cutoff:.0f}\n")
    f.write(f"source: {'configured default_min_dist (no clear antimode)' if used_default else 'data-driven antimode'}\n")

# --- histogram of nearest-neighbour distances with the cutoff marked -------
d = bp_table["nearest_bp"].to_numpy()
d = d[d > 0]

if len(d) > 1:
    plt.rcParams["font.size"] = 18
    fig, ax = plt.subplots(figsize=(15, 5))
    bins = np.logspace(np.log10(d.min()), np.log10(d.max()), 286)
    bins = np.unique(np.append(bins, cutoff))

    counts, edges, patches = ax.hist(d, bins=bins)
    for patch, left in zip(patches, edges[:-1]):
        if left < cutoff:
            patch.set_facecolor("#181717"); patch.set_alpha(1.0)
        else:
            patch.set_facecolor("#5903b5"); patch.set_alpha(0.9)
    ax.axvline(cutoff, color="black", ls="--", lw=3)
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1000:g}" if x >= 1000 else f"{x:g} bp"))
    ax.set_xlabel("Distance to the Nearest Breakpoint (kb)", fontsize=25)
    ax.set_ylabel("Number of Breakpoints", fontsize=25)
    ax.set_title(f"chr{chrom}")
    plt.tight_layout()
    plt.savefig(snakemake.output.hist_png, dpi=150)
    plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.text(0.5, 0.5, "not enough breakpoints to plot", ha="center", va="center")
    plt.savefig(snakemake.output.hist_png, dpi=150)
    plt.close(fig)

# --- apply the filter and write the fully-filtered breakpoint set ----------
filtered = nearest_neighbor_filter(final_bps, min_dist=cutoff, chrom_length=chrom_length)
total = sum(len(b) for b in filtered.values())
n_hap = len(filtered)
print(f"[chr{chrom}] surviving breakpoints: {total} ({total / n_hap:.2f}/hap)")

filtered_table = breakpoint_distance_table(filtered, nodes_meta, chrom_length=chrom_length)
filtered_table.insert(0, "chrom", chrom)
filtered_table.to_csv(snakemake.output.filtered, index=False)

# bp-precision companion to `filtered` (the report table above rounds to Mb
# for readability) — this is what every downstream stage actually consumes.
filtered_bps_table = pd.DataFrame(
    [{"fsj_id": nodes_meta[node][0], "node": int(node), "position_bp": int(p)}
     for node, bps in filtered.items() for p in bps]
).sort_values(["fsj_id", "node", "position_bp"]).reset_index(drop=True)
filtered_bps_table.insert(0, "chrom", chrom)
filtered_bps_table.to_csv(snakemake.output.filtered_bps, index=False)
