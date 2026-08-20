"""Rule: simulate_null — run N fixed-pedigree msprime simulations (one per
random seed) against the chromosome-specific pedigree, and tabulate
per-individual breakpoint counts for each.

Each seed's tree sequence is dumped to its own file under the null_trees
output directory (null_seed<seed>.trees) so the tss_tes_distance stage can
treat every replicate independently, rather than only keeping the last seed
(as a single fixed "output.trees" path would silently do).
"""

import os

import msprime
import pandas as pd
import tskit
from common import build_breakpoints_table_fast, build_id_table

chrom = snakemake.wildcards.chrom
recombination_rate = float(snakemake.params.recombination_rate)
seeds = range(1, int(snakemake.params.n_sims) + 1)

pedigree_tables = tskit.load(snakemake.input.pedigree).tables
sequence_length = pedigree_tables.sequence_length

null_trees_dir = snakemake.output.null_trees_dir
os.makedirs(null_trees_dir, exist_ok=True)

breakpoints_out_path = snakemake.output.breakpoints
write_header = True

for seed in seeds:
    null_ts = msprime.sim_ancestry(
        initial_state=pedigree_tables,
        model="fixed_pedigree",
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=seed,
        additional_nodes=(
            msprime.NodeType.RECOMBINANT
            | msprime.NodeType.PASS_THROUGH
            | msprime.NodeType.COMMON_ANCESTOR
        ),
        coalescing_segments_only=False,
    )
    null_ts.dump(os.path.join(null_trees_dir, f"null_seed{seed}.trees"))

    _, ped_to_ts, _ = build_id_table(null_ts)
    bp_df = build_breakpoints_table_fast(null_ts, ped_to_ts)
    bp_df["seed"] = seed
    bp_df.insert(0, "chrom", chrom)
    bp_df.to_csv(breakpoints_out_path, mode="a", header=write_header, index=False)
    write_header = False

    print(f"[chr{chrom}] seed {seed} done")

print(f"Done. Results saved to {breakpoints_out_path}")
