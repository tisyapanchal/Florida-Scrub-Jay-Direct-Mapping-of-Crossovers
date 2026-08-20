"""Rule: null_summary — per-individual variance of breakpoint counts across
the N null replicates (sanity check that the null model behaves consistently)."""

import pandas as pd

chrom = snakemake.wildcards.chrom

null_df = pd.read_csv(snakemake.input.breakpoints)

print(null_df["total_breakpoints"].describe())
print("Variance:", null_df["total_breakpoints"].var())

per_individual_variance = (
    null_df.groupby("ped_id")["total_breakpoints"]
    .agg(["mean", "var", "std", "min", "max"])
    .reset_index()
)
per_individual_variance.insert(0, "chrom", chrom)
per_individual_variance.to_csv(snakemake.output[0], index=False)
