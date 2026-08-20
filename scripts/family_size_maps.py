"""Rule: family_size_maps — recombination maps and breakpoint tables split by
full-sib family size (<5, 5-10, 11-15, 16-20, >20 siblings), to check whether
larger families behave differently (e.g. more power to detect rare events,
or family-specific artifacts)."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

from common import breakpoint_distance_table, get_id, read_sibships

chrom = snakemake.wildcards.chrom

real = tskit.load(snakemake.input.real_trees)
chrom_len = int(real.sequence_length)
_, sibgroups = read_sibships(snakemake.input.ped)
id_to_ind = {get_id(ind): ind for ind in real.individuals()}

filtered_table = pd.read_csv(snakemake.input.filtered)
filtered = {node: g.position_bp.to_numpy() for node, g in filtered_table.groupby("node")}
nodes_meta = {node: (g.fsj_id.iloc[0],) for node, g in filtered_table.groupby("node")}

# mutually exclusive bins, first-match-wins so shared boundaries (5, 10, 15, 20)
# fall into the lower bucket only, never double-counted
size_bins = {
    "lt5": (0, 4),
    "5to10": (5, 10),
    "11to15": (11, 15),
    "16to20": (16, 20),
    "gt20": (21, np.inf),
}

bps_by_size = {label: {} for label in size_bins}
for (father, mother), offspring in sibgroups.items():
    n_sibs = len(offspring)
    label = next(l for l, (lo, hi) in size_bins.items() if lo <= n_sibs <= hi)
    for fsj in offspring:
        if fsj not in id_to_ind:
            continue
        for node in id_to_ind[fsj].nodes:
            if node in filtered:
                bps_by_size[label][node] = filtered[node]

for label, d in bps_by_size.items():
    n_bp = sum(len(b) for b in d.values())
    print(f"[chr{chrom}] {label:8s}: {len(d):4d} haplotypes, {n_bp:6d} breakpoints")

bin_output_map = {
    "lt5": snakemake.output.lt5,
    "5to10": snakemake.output.bin_5to10,
    "11to15": snakemake.output.bin_11to15,
    "16to20": snakemake.output.bin_16to20,
    "gt20": snakemake.output.gt20,
}
for label, d in bps_by_size.items():
    bin_table = breakpoint_distance_table(d, nodes_meta, chrom_length=chrom_len)
    bin_table.insert(0, "chrom", chrom)
    bin_table.to_csv(bin_output_map[label], index=False)

# --- averaged Marey map per family-size bucket ------------------------------
N_BINS = 200
EDGES = np.linspace(0, chrom_len, N_BINS)
CENTERS = (EDGES[:-1] + EDGES[1:]) / 2


def family_cumulative_marey(nodes):
    nodes = [n for n in nodes if n in filtered]
    n_gam = len(nodes)
    if n_gam == 0:
        return None
    allpos = np.sort(np.concatenate([filtered[n] for n in nodes]).astype(int))
    counts, _ = np.histogram(allpos, bins=EDGES)
    return np.cumsum(counts / n_gam) * 100


def family_nodes(offspring):
    nodes = []
    for fsj in offspring:
        if fsj not in id_to_ind:
            continue
        nodes.extend(n for n in id_to_ind[fsj].nodes if n in filtered)
    return nodes


bucket_family_curves = {label: [] for label in size_bins}
for (father, mother), offspring in sibgroups.items():
    n_sibs = len(offspring)
    for label, (lo, hi) in size_bins.items():
        if lo <= n_sibs <= hi:
            curve = family_cumulative_marey(family_nodes(offspring))
            if curve is not None:
                bucket_family_curves[label].append(curve)
            break

fig, ax = plt.subplots(figsize=(13, 7))
cmap = plt.cm.viridis(np.linspace(0, 1, len(size_bins)))
colors = dict(zip(size_bins, cmap))
labels = {
    "lt5": "< 5 siblings", "5to10": "5-10 siblings", "11to15": "11-15 siblings",
    "16to20": "16-20 siblings", "gt20": "> 20 siblings",
}
for label, color in colors.items():
    curves = np.array(bucket_family_curves[label])
    if len(curves) == 0:
        continue
    ax.plot(CENTERS / 1e6, curves.mean(axis=0), "-", color=color, lw=2.5, label=labels[label])

ax.set_xlabel("Chromosome Length (Mb)")
ax.set_ylabel("Cumulative Genetic Length (cM)")
ax.set_title(f"chr{chrom}")
ax.legend(fontsize=14)
plt.tight_layout()
plt.savefig(snakemake.output.marey_png, dpi=150)
plt.close(fig)
