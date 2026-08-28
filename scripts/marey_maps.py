"""Rule: marey_maps — direct-mapping Marey maps (female vs male, from the
fully-filtered real breakpoints) and their comparison against the
Romero et al. (2024) linkage-disequilibrium map for this chromosome."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

from common import build_mareymap, informative_nodes, is_maternal

chrom = snakemake.wildcards.chrom

real = tskit.load(snakemake.input.real_trees)
chrom_len = int(real.sequence_length)

filtered_table = pd.read_csv(snakemake.input.filtered)
filtered = {int(node): g.position_bp.to_numpy()
            for node, g in filtered_table.groupby("node")}

# Every informative gamete counts toward the denominator, including the ones
# whose breakpoints were all filtered out — those rows are absent from the CSV.
nodes = informative_nodes(real)
mat = {int(n): filtered.get(int(n), np.empty(0)) for n in nodes if is_maternal(real, n) is True}
pat = {int(n): filtered.get(int(n), np.empty(0)) for n in nodes if is_maternal(real, n) is False}

edges_f, cM_f = build_mareymap(mat, chrom_len)
edges_m, cM_m = build_mareymap(pat, chrom_len)

# --- direct-mapping plot ----------------------------------------------------
plt.rcParams.update(plt.rcParamsDefault)
fig, ax = plt.subplots(figsize=(11, 6))
ax.plot(edges_f / 1e6, cM_f, "-", color="#ff5500", lw=3, label=f"female ({cM_f[-1]:.0f} cM)")
ax.plot(edges_m / 1e6, cM_m, "-", color="#0055ff", lw=3, label=f"male ({cM_m[-1]:.0f} cM)")
ax.set_xlabel("physical position (Mb)")
ax.set_ylabel("cumulative genetic position (cM)")
ax.set_title(f"Direct Mapping chr{chrom} Female vs Male Recombination Map")
ax.legend()
plt.tight_layout()
plt.savefig(snakemake.output.dm_png, dpi=150)
plt.close(fig)

# --- linkage map for this chromosome ---------------------------------------
lm = pd.read_csv(snakemake.input.linkage_map, sep="\t")
scaffold_stripped = (
    lm.Scaffold.astype(str).str.strip().str.replace(r"^(chr|Chr|scaffold_?)", "", regex=True)
)
lm1 = lm[scaffold_stripped == str(chrom)].copy()
for c in ["bpPosition", "cMPosition.Male", "cMPosition.Female"]:
    lm1[c] = pd.to_numeric(lm1[c], errors="coerce")
lm1 = lm1.dropna(subset=["bpPosition"]).sort_values("bpPosition").reset_index(drop=True)
print(f"[chr{chrom}] linkage map: {len(lm1)} markers")

if len(lm1):
    map_female_cM = lm1["cMPosition.Female"].max()
    map_male_cM = lm1["cMPosition.Male"].max()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(lm1.bpPosition / 1e6, lm1["cMPosition.Male"], color="#0055ff", lw=3, alpha=0.8, label="Male Recombination Rate")
    ax.plot(lm1.bpPosition / 1e6, lm1["cMPosition.Female"], color="#ff5500", lw=3, alpha=0.8, label="Female Recombination Rate")
    ax.set_xlabel("Physical Chromosomal Position (Mb)", fontsize=14)
    ax.set_ylabel("Cumulative Genetic Position (cM)", fontsize=14)
    ax.set_title(f"chr{chrom}")
    ax.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(snakemake.output.lm_png, dpi=150)
    plt.close(fig)

    # --- direct mapping vs linkage map, combined ----------------------------
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.plot(edges_f / 1e6, cM_f, "-", color="#ef8636", lw=5, label=f"Female Direct Mapping ({cM_f[-1]:.0f} cM)")
    ax.plot(edges_m / 1e6, cM_m, "-", color="#3b75af", lw=5, label=f"Male Direct Mapping ({cM_m[-1]:.0f} cM)")
    ax.plot(lm1.bpPosition / 1e6, lm1["cMPosition.Male"], color="#ef8636", ls="--", lw=4,
            label=f"Male Linkage Map ({map_male_cM:.1f} cM)")
    ax.plot(lm1.bpPosition / 1e6, lm1["cMPosition.Female"], color="#3b75af", ls="--", lw=4,
            label=f"Female Linkage Map ({map_female_cM:.1f} cM)")
    ax.set_xlabel("Chromosome Length (Mb)", fontsize=30)
    ax.set_ylabel("Cumulative Genetic Length (cM)", fontsize=30)
    ax.set_title(f"chr{chrom}", fontsize=30)
    ax.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(snakemake.output.combined_png, dpi=200)
    plt.close(fig)
else:
    map_female_cM = map_male_cM = float("nan")
    for out in (snakemake.output.lm_png, snakemake.output.combined_png):
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, f"no linkage-map markers found for chrom {chrom!r}", ha="center", va="center")
        plt.savefig(out, dpi=150)
        plt.close(fig)

summary = pd.DataFrame([{
    "chrom": chrom,
    "direct_map_female_cM": cM_f[-1],
    "direct_map_male_cM": cM_m[-1],
    "linkage_map_female_cM": map_female_cM,
    "linkage_map_male_cM": map_male_cM,
    "linkage_map_n_markers": len(lm1),
}])
summary.to_csv(snakemake.output.cM_summary, index=False)
