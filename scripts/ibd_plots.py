"""Rule: ibd_plots — per-individual IBD-tract diagnostic plots for a focal
individual, at three filtering stages: unfiltered ("kept"), after the
full-sibling filters, and after every filter (sibling + half-sib + nearest
neighbour). Not run for all 2,602 individuals by default — opt in specific
IDs via config["ibd_plot_individuals"] (see Snakefile / config.yaml)."""

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import tskit
from plot_ibd_tracts import plot_ibd

chrom = snakemake.wildcards.chrom
focal_id = snakemake.wildcards.individual

real = tskit.load(snakemake.input.real_trees)

rec = pd.read_csv(snakemake.input.rec)
opp = pd.read_csv(snakemake.input.opp)
final_bps_table = pd.read_csv(snakemake.input.final_bps)

# unfiltered
fig = plot_ibd(real, focal_id, chrom=chrom, title=f"{focal_id} — unfiltered  chr{chrom}")
fig.savefig(snakemake.output.kept_png, dpi=150, bbox_inches="tight")

# excluded by the full-sibling filters (shared recurrence + opposite-direction)
shared = set(rec[(rec.fsj_id == focal_id) & (rec.n_sibs_matching >= 1)]["position_bp"])
opp_excl = set(opp[opp.sib_a == focal_id]["pos_a"]) | set(opp[opp.sib_b == focal_id]["pos_b"])
sib_excl = shared | opp_excl
fig = plot_ibd(real, focal_id, chrom=chrom, exclude_positions=sib_excl, exclude_tol=1.0,
               title=f"{focal_id} — sibling-filtered  chr{chrom}")
fig.savefig(snakemake.output.sibling_filtered_png, dpi=150, bbox_inches="tight")

# excluded by every filter: anything not present in the final survivor set
survivor_pos = set(final_bps_table[final_bps_table.fsj_id == focal_id]["position_bp"])
from common import informative_nodes

nodes = informative_nodes(real)
all_raw_pos = set()
for node, (fsj, bps) in nodes.items():
    if fsj == focal_id:
        all_raw_pos.update(int(p) for p in bps)
all_excl = all_raw_pos - survivor_pos
fig = plot_ibd(real, focal_id, chrom=chrom, exclude_positions=all_excl, exclude_tol=1.0,
               title=f"{focal_id} — all filters  chr{chrom}")
fig.savefig(snakemake.output.all_filters_png, dpi=150, bbox_inches="tight")
