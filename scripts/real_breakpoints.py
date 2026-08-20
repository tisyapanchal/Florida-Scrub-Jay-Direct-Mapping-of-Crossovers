"""Rule: real_breakpoints — per-individual breakpoint counts from the real
(phased) tree sequence, plus a REAL-vs-NULL crossovers-per-gamete QC summary.
"""

import glob
import os

import tskit
from common import build_breakpoints_table_fast, build_id_table, co_per_gamete

chrom = snakemake.wildcards.chrom

real = tskit.load(snakemake.input.real_trees)

id_table, ped_to_ts, ts_to_ped = build_id_table(real)
id_table.insert(0, "chrom", chrom)
id_table.to_csv(snakemake.output.id_table, index=False)
print(f"[chr{chrom}] {len(id_table)} real individuals mapped")

breakpoints_df = build_breakpoints_table_fast(real, ped_to_ts)
breakpoints_df.insert(0, "chrom", chrom)
breakpoints_df.to_csv(snakemake.output.bp_table, index=False)

# QC: crossovers/gamete, real vs the null replicates simulated for this chromosome
co_real, _ = co_per_gamete(real)
null_files = sorted(glob.glob(os.path.join(snakemake.input.null_trees_dir, "null_seed*.trees")))
co_null_means = []
for path in null_files:
    null_ts = tskit.load(path)
    co_null, _ = co_per_gamete(null_ts)
    co_null_means.append(co_null.mean())

with open(snakemake.output.qc, "w") as f:
    f.write(f"chrom: {chrom}\n")
    f.write(f"REAL: {real.num_individuals} ind | {real.sequence_length:.0f} bp | {real.num_trees} trees\n")
    f.write(f"REAL CO/gamete: {co_real.mean():.3f}\n")
    if co_null_means:
        f.write(
            f"NULL CO/gamete: {sum(co_null_means) / len(co_null_means):.3f} "
            f"(mean over {len(co_null_means)} replicates)\n"
        )
    else:
        f.write("NULL CO/gamete: no null replicates found\n")

print(open(snakemake.output.qc).read())
