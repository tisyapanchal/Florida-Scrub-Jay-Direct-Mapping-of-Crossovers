"""Rule: build_pedigree — build the msprime PedigreeBuilder tree sequence for
one chromosome's sequence length, from the shared .ped file."""

import tskit
from common import build_and_save_pedigree

real = tskit.load(snakemake.input.real_trees)
build_and_save_pedigree(
    ped_path=snakemake.input.ped,
    pedigree_out_path=snakemake.output.pedigree,
    sequence_length=real.sequence_length,
)
