"""Rule: gene_annotation_overlap — how filtered real breakpoints distribute
across genic/intergenic/exon/intron/CDS regions, relative to what's expected
by chance given how much of the chromosome each region covers."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

from common import frac_in, load_features, merge_intervals

chrom = snakemake.wildcards.chrom
gff_chrom_name = snakemake.params.gff_chrom_name
gff_path = snakemake.input.gff

real = tskit.load(snakemake.input.real_trees)
chrom_len = int(real.sequence_length)

filtered_table = pd.read_csv(snakemake.input.filtered)
filt_pos = np.sort(filtered_table.position_bp.to_numpy())

feats = load_features(gff_path, gff_chrom_name, ["exon", "CDS", "gene", "intron"])
merged_feats = {ft: merge_intervals(iv) for ft, iv in feats.items()}

# --- feature enrichment (exon / CDS / gene) --------------------------------
rows = []
for name in ["exon", "CDS", "gene"]:
    merged, cov = merged_feats[name]
    hit, obs = frac_in(filt_pos, merged)
    exp = cov / chrom_len
    enr = obs / exp if exp else float("nan")
    rows.append({
        "chrom": chrom, "feature": name, "n_breakpoints_in_feature": hit,
        "observed_frac": obs, "genome_frac": exp, "enrichment": enr,
    })
    print(f"[chr{chrom}] {name:16s}: {hit:5d} bp in feature | obs {obs:.1%} | genome {exp:.1%} | enrichment {enr:.2f}x")
pd.DataFrame(rows).to_csv(snakemake.output.enrichment, index=False)

# --- weighted proportion of breakpoints per region --------------------------
n_real = len(filt_pos)
regions = ["gene", "intergenic", "exon", "intron"]
weighted = []
for ft in regions:
    if ft == "intergenic":
        gene_merged, gene_cov = merged_feats["gene"]
        _, obs_frac = frac_in(filt_pos, gene_merged)
        obs_frac = 1 - obs_frac
        genome_frac = 1 - gene_cov / chrom_len
    else:
        merged, cov = merged_feats[ft]
        _, obs_frac = frac_in(filt_pos, merged)
        genome_frac = cov / chrom_len
    observed = obs_frac * n_real
    weighted.append(observed / genome_frac if genome_frac else float("nan"))

x = np.arange(len(regions))
PINK, PURPLE = "#e377c2", "#9467bd"
bar_colors = [PURPLE if r in ("intergenic", "intron") else PINK for r in regions]

fig, ax = plt.subplots(figsize=(18, 8))
ax.bar(x, weighted, 0.40, color=bar_colors, edgecolor="black", linewidth=0.5)
ax.set_xticks(x); ax.set_xticklabels(regions, fontsize=24)
ax.set_ylabel("Crossover Count weighted by genome proportion", fontsize=24)
ax.set_title(f"chr{chrom}")
finite_weighted = [w for w in weighted if np.isfinite(w)]
for i, wt in enumerate(weighted):
    if np.isfinite(wt):
        ax.text(i, wt + max(finite_weighted, default=1) * 0.01, f"{wt:.0f}", ha="center", fontsize=20, color="#000000")
plt.tight_layout()
plt.savefig(snakemake.output.region_png, dpi=150)
plt.close(fig)
