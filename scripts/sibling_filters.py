"""Rule: sibling_filters — clean phasing errors out of the real breakpoint set
using two complementary sibling-recurrence filters:

  1. full-sibling recurrence: a breakpoint shared (same parental haplotype,
     within `sibling_window` bp) by >=2 full siblings is flagged as a likely
     shared phasing error rather than an independent crossover.
  2. opposite-direction switches: a breakpoint where two full siblings switch
     parental haplotype in opposite directions at (nearly) the same position.
  3. half-sibling recurrence: a breakpoint shared by >=2 offspring of the same
     single parent (catches recurrence full-sib filtering misses).

The survivors of all three filters are the "sibling_filtered_bps" set used by
every downstream analysis stage.
"""

import numpy as np
import pandas as pd
import tskit
from common import (
    all_parental_switches,
    get_id,
    informative_nodes,
    opposite_direction_flags,
    per_parent_recurrence,
    read_sibships,
    sibling_recurrence_pairwise,
)

chrom = snakemake.wildcards.chrom
sibling_window = float(snakemake.params.sibling_window)
half_sib_window = float(snakemake.params.half_sib_window)
half_sib_min_share = int(snakemake.params.half_sib_min_share)

real = tskit.load(snakemake.input.real_trees)
parents, sibgroups = read_sibships(snakemake.input.ped)
print(f"[chr{chrom}] {len(sibgroups)} sibships with >=2 sibs")

nodes_real = informative_nodes(real)

# --- full-sibling recurrence filter -----------------------------------------
rec = sibling_recurrence_pairwise(nodes_real, sibgroups, real, window=sibling_window)
rec.insert(0, "chrom", chrom)
rec.to_csv(snakemake.output.rec, index=False)

in_sib = rec.sibship.notna()
sib_summary = (
    rec[in_sib]
    .groupby("sibship")
    .agg(n_sibs=("fsj_id", "nunique"),
         total_bp=("position_bp", "size"),
         shared_bp=("n_sibs_matching", lambda s: (s >= 1).sum()),
         high_conf_shared=("n_sibs_matching", lambda s: (s >= 2).sum()))
    .reset_index()
)
sib_summary["parents"] = sib_summary["sibship"]
sib_summary["siblings"] = sib_summary["sibship"].map(
    lambda s: ", ".join(sibgroups.get(tuple(s.split("x")), [])))
sib_summary["frac_shared"] = (sib_summary.shared_bp / sib_summary.total_bp).round(3)
sib_summary = sib_summary.sort_values("shared_bp", ascending=False)
sib_summary = sib_summary[
    ["parents", "siblings", "n_sibs", "total_bp", "shared_bp", "high_conf_shared", "frac_shared"]
]
sib_summary.insert(0, "chrom", chrom)
sib_summary.to_csv(snakemake.output.sib_summary, index=False)
print(f"[chr{chrom}] {len(sib_summary)} sibships | {sib_summary.shared_bp.sum()} shared breakpoints")

# --- opposite-direction-switch filter ---------------------------------------
id_to_ind = {get_id(ind): ind for ind in real.individuals()}
switches, off_to_sib = all_parental_switches(real, sibgroups, parents, id_to_ind)
flags = opposite_direction_flags(switches, off_to_sib, sibgroups, id_to_ind, window=sibling_window)
flags.insert(0, "chrom", chrom)
flags.to_csv(snakemake.output.opp, index=False)
print(f"[chr{chrom}] {len(flags)} opposite-direction coincidences")

from collections import defaultdict

opp_by_ind = defaultdict(set)
for _, r in flags.iterrows():
    opp_by_ind[r.sib_a].add(int(r.pos_a))
    opp_by_ind[r.sib_b].add(int(r.pos_b))

shared_by_ind = defaultdict(set)
for fsj_id, sub in rec[rec.n_sibs_matching >= 1].groupby("fsj_id"):
    shared_by_ind[fsj_id] = set(sub.position_bp)

survivor_bps = {}
for node, (fsj, bps) in nodes_real.items():
    drop = shared_by_ind.get(fsj, set()) | opp_by_ind.get(fsj, set())
    kept = np.array([p for p in bps if int(p) not in drop])
    survivor_bps[node] = (fsj, kept)

total = sum(len(b) for _, b in survivor_bps.values())
print(f"[chr{chrom}] survivor_bps (post full-sib filters): {total} breakpoints")

# --- half-sibling (per-parent) recurrence filter ----------------------------
flagged = per_parent_recurrence(survivor_bps, parents, window=half_sib_window, min_share=half_sib_min_share)

rows = []
for node, (fsj, bps) in survivor_bps.items():
    drop = flagged.get(fsj, set())
    for p in bps:
        if p not in drop:
            rows.append({"fsj_id": fsj, "node": int(node), "position_bp": int(p),
                         "position_Mb": round(p / 1e6, 3)})
survivors_df = pd.DataFrame(rows, columns=["fsj_id", "node", "position_bp", "position_Mb"])
if len(survivors_df):
    survivors_df = survivors_df.sort_values(["fsj_id", "position_bp"]).reset_index(drop=True)
print(f"[chr{chrom}] {len(survivors_df)} breakpoints survive the per-parent (half-sib) filter")

# --- final_bps: the fully sibling+half-sib filtered set, used everywhere downstream
final_bps = {
    node: np.sort(g.position_bp.to_numpy())
    for node, g in survivors_df.groupby("node")
}
for node in survivor_bps:
    final_bps.setdefault(node, np.array([], dtype=int))

final_bps_table = pd.DataFrame(
    [{"fsj_id": nodes_real[node][0], "node": int(node), "position_bp": int(p)}
     for node, bps in final_bps.items() for p in bps]
).sort_values(["fsj_id", "node", "position_bp"]).reset_index(drop=True)
final_bps_table.insert(0, "chrom", chrom)
final_bps_table.to_csv(snakemake.output.final_bps, index=False)

total = len(final_bps_table)
n_hap = len(final_bps)
print(f"[chr{chrom}] final_bps: {total} breakpoints across {n_hap} haplotypes "
      f"({total / n_hap:.2f}/hap)")
