"""Rule: density_plots — crossover density along the chromosome, split by
parent of origin (maternal vs paternal), for both the sibling-filtered and
the fully-filtered (sibling + nearest-neighbour) breakpoint sets."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

from common import is_maternal

chrom = snakemake.wildcards.chrom
real = tskit.load(snakemake.input.real_trees)
chrom_len_mb = real.sequence_length / 1e6


def load_bps_by_node(path):
    df = pd.read_csv(path)
    return {node: g.position_bp.to_numpy() for node, g in df.groupby("node")}


def split_by_sex(bps_by_node):
    mat, pat, n_mat, n_pat = [], [], 0, 0
    for node, bps in bps_by_node.items():
        if len(bps) == 0:
            continue
        if is_maternal(real, node):
            mat.extend(bps); n_mat += 1
        else:
            pat.extend(bps); n_pat += 1
    return np.array(mat), np.array(pat), n_mat, n_pat


def plot_density(bps_by_node, title, out_path):
    mat, pat, n_mat, n_pat = split_by_sex(bps_by_node)
    print(f"[chr{chrom}] {title}: maternal {len(mat)} on {n_mat} hap, paternal {len(pat)} on {n_pat} hap")

    nbins = 120
    bins = np.linspace(0, chrom_len_mb + 1, nbins + 1)
    h_mat, _ = np.histogram(mat / 1e6, bins=bins, density=True) if len(mat) else (np.zeros(nbins), None)
    h_pat, _ = np.histogram(pat / 1e6, bins=bins, density=True) if len(pat) else (np.zeros(nbins), None)
    centers = (bins[:-1] + bins[1:]) / 2

    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(centers, h_mat, "-o", color="#e34948", alpha=0.7, ms=4, label=f"maternal (n={n_mat})")
    ax.plot(centers, h_pat, "-s", color="#2a78d6", alpha=0.7, ms=4, label=f"paternal (n={n_pat})")
    ax.set_xlabel("chromosome position (Mb)")
    ax.set_ylabel("density")
    ax.set_title(f"{title} — chr{chrom}")
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


plot_density(
    load_bps_by_node(snakemake.input.final_bps),
    "Recombination: maternal vs paternal (sibling-filtered)",
    snakemake.output.raw_png,
)
plot_density(
    load_bps_by_node(snakemake.input.filtered),
    "Recombination: maternal vs paternal (fully filtered)",
    snakemake.output.filt_png,
)
