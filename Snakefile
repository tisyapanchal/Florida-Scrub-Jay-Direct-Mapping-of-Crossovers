"""FSJ direct recombination mapping pipeline.

Snakemake port of notebooks/scrubjay.ipynb, parameterized over chromosome so
the same pipeline runs on any chromosome listed in config["chromosomes"]
instead of being hardcoded to chromosome 1. See workflow/README.md for a full
description of each stage and its outputs.

Run from anywhere; `workdir` below always pins execution to the repo root so
paths in config.yaml can stay simple ("notebooks/...", "results/...").
"""

import os

workdir: os.path.normpath(os.path.join(os.path.dirname(workflow.snakefile), ".."))

configfile: "workflow/config.yaml"

CHROMS = [str(c) for c in config["chromosomes"]]
RESULTS = config["results_dir"]


def gff_chrom_name(chrom):
    return config["gff_chrom_pattern"].format(chrom=chrom)


rule all:
    input:
        expand(RESULTS + "/chr{chrom}/per_individual_variance_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/qc_summary_chr{chrom}.txt", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/sibship_sharing_summary_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/all_filters_breakpoints_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/mat_vs_pat_filtered_chr{chrom}.png", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/dm_vs_lm_chr{chrom}.png", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/direct_map_cM_summary_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/feature_enrichment_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/breakpoints_by_region_chr{chrom}.png", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/tss_tes_distance_summary_chr{chrom}.csv", chrom=CHROMS),
        expand(RESULTS + "/chr{chrom}/marey_by_family_size_chr{chrom}.png", chrom=CHROMS),
        expand(
            RESULTS + "/chr{chrom}/breakpoints_{bin}_chr{chrom}.csv",
            chrom=CHROMS, bin=["lt5", "5to10", "11to15", "16to20", "gt20"],
        ),
        expand(
            RESULTS + "/chr{chrom}/ibd_plots/{individual}_all_filters_ibd_chr{chrom}.png",
            chrom=CHROMS, individual=config["ibd_plot_individuals"],
        ) if config["ibd_plot_individuals"] else [],


rule build_pedigree:
    input:
        real_trees=config["real_trees_pattern"],
        ped=config["ped_path"],
    output:
        pedigree=RESULTS + "/chr{chrom}/final_null_pedigree_chr{chrom}.trees",
    script:
        "scripts/build_pedigree.py"


rule simulate_null:
    input:
        pedigree=rules.build_pedigree.output.pedigree,
    params:
        recombination_rate=config["null_recombination_rate"],
        n_sims=config["n_null_simulations"],
    output:
        breakpoints=RESULTS + "/chr{chrom}/finalized_null_breakpoints_100sims_chr{chrom}.csv",
        null_trees_dir=directory(RESULTS + "/chr{chrom}/null_trees"),
    script:
        "scripts/simulate_null.py"


rule null_summary:
    input:
        breakpoints=rules.simulate_null.output.breakpoints,
    output:
        RESULTS + "/chr{chrom}/per_individual_variance_chr{chrom}.csv",
    script:
        "scripts/null_summary.py"


rule real_breakpoints:
    input:
        real_trees=config["real_trees_pattern"],
        null_trees_dir=rules.simulate_null.output.null_trees_dir,
    output:
        id_table=RESULTS + "/chr{chrom}/real_id_table_chr{chrom}.csv",
        bp_table=RESULTS + "/chr{chrom}/real_breakpoints_table_chr{chrom}.csv",
        qc=RESULTS + "/chr{chrom}/qc_summary_chr{chrom}.txt",
    script:
        "scripts/real_breakpoints.py"


rule sibling_filters:
    input:
        real_trees=config["real_trees_pattern"],
        ped=config["ped_path"],
    params:
        sibling_window=config["sibling_window_bp"],
        half_sib_window=config["half_sibling_window_bp"],
        half_sib_min_share=config["half_sibling_min_share"],
    output:
        rec=RESULTS + "/chr{chrom}/breakpoints_sibling_recurrence_chr{chrom}.csv",
        sib_summary=RESULTS + "/chr{chrom}/sibship_sharing_summary_chr{chrom}.csv",
        opp=RESULTS + "/chr{chrom}/opposite_direction_switches_chr{chrom}.csv",
        final_bps=RESULTS + "/chr{chrom}/sibling_filtered_bps_chr{chrom}.csv",
    script:
        "scripts/sibling_filters.py"


rule distance_filter:
    input:
        real_trees=config["real_trees_pattern"],
        final_bps=rules.sibling_filters.output.final_bps,
    params:
        default_min_dist=config["default_nearest_neighbor_cutoff_bp"],
    output:
        dist_table=RESULTS + "/chr{chrom}/distances_btwn_sibling_filtered_breakpoints_chr{chrom}.csv",
        hist_png=RESULTS + "/chr{chrom}/distance_btwn_breakpoint_distribution_chr{chrom}.png",
        filtered=RESULTS + "/chr{chrom}/all_filters_breakpoints_chr{chrom}.csv",
        filtered_bps=RESULTS + "/chr{chrom}/filtered_bps_chr{chrom}.csv",
        cutoff_txt=RESULTS + "/chr{chrom}/nearest_neighbor_cutoff_chr{chrom}.txt",
    script:
        "scripts/distance_filter.py"


rule density_plots:
    input:
        real_trees=config["real_trees_pattern"],
        final_bps=rules.sibling_filters.output.final_bps,
        filtered=rules.distance_filter.output.filtered_bps,
    output:
        raw_png=RESULTS + "/chr{chrom}/mat_vs_pat_chr{chrom}.png",
        filt_png=RESULTS + "/chr{chrom}/mat_vs_pat_filtered_chr{chrom}.png",
    script:
        "scripts/density_plots.py"


rule marey_maps:
    input:
        real_trees=config["real_trees_pattern"],
        filtered=rules.distance_filter.output.filtered_bps,
        linkage_map=config["linkage_map_path"],
    output:
        dm_png=RESULTS + "/chr{chrom}/dm_chr{chrom}.png",
        lm_png=RESULTS + "/chr{chrom}/lm_map_chr{chrom}.png",
        combined_png=RESULTS + "/chr{chrom}/dm_vs_lm_chr{chrom}.png",
        cM_summary=RESULTS + "/chr{chrom}/direct_map_cM_summary_chr{chrom}.csv",
    script:
        "scripts/marey_maps.py"


rule gene_annotation_overlap:
    input:
        real_trees=config["real_trees_pattern"],
        filtered=rules.distance_filter.output.filtered_bps,
        gff=config["gff_path"],
    params:
        gff_chrom_name=lambda wc: gff_chrom_name(wc.chrom),
    output:
        enrichment=RESULTS + "/chr{chrom}/feature_enrichment_chr{chrom}.csv",
        region_png=RESULTS + "/chr{chrom}/breakpoints_by_region_chr{chrom}.png",
    script:
        "scripts/gene_annotation_overlap.py"


rule tss_tes_distance:
    input:
        filtered=rules.distance_filter.output.filtered_bps,
        gff=config["gff_path"],
        null_trees_dir=rules.simulate_null.output.null_trees_dir,
    params:
        gff_chrom_name=lambda wc: gff_chrom_name(wc.chrom),
    output:
        tss_png=RESULTS + "/chr{chrom}/bp_dist_to_tss_chr{chrom}.png",
        tes_png=RESULTS + "/chr{chrom}/bp_dist_to_tes_chr{chrom}.png",
        summary_csv=RESULTS + "/chr{chrom}/tss_tes_distance_summary_chr{chrom}.csv",
    script:
        "scripts/tss_tes_distance.py"


rule family_size_maps:
    input:
        real_trees=config["real_trees_pattern"],
        ped=config["ped_path"],
        filtered=rules.distance_filter.output.filtered_bps,
    output:
        marey_png=RESULTS + "/chr{chrom}/marey_by_family_size_chr{chrom}.png",
        lt5=RESULTS + "/chr{chrom}/breakpoints_lt5_chr{chrom}.csv",
        bin_5to10=RESULTS + "/chr{chrom}/breakpoints_5to10_chr{chrom}.csv",
        bin_11to15=RESULTS + "/chr{chrom}/breakpoints_11to15_chr{chrom}.csv",
        bin_16to20=RESULTS + "/chr{chrom}/breakpoints_16to20_chr{chrom}.csv",
        gt20=RESULTS + "/chr{chrom}/breakpoints_gt20_chr{chrom}.csv",
    script:
        "scripts/family_size_maps.py"


rule ibd_plots:
    input:
        real_trees=config["real_trees_pattern"],
        rec=rules.sibling_filters.output.rec,
        opp=rules.sibling_filters.output.opp,
        final_bps=rules.sibling_filters.output.final_bps,
    output:
        kept_png=RESULTS + "/chr{chrom}/ibd_plots/{individual}_kept_ibd_chr{chrom}.png",
        sibling_filtered_png=RESULTS + "/chr{chrom}/ibd_plots/{individual}_sibling_filtered_ibd_chr{chrom}.png",
        all_filters_png=RESULTS + "/chr{chrom}/ibd_plots/{individual}_all_filters_ibd_chr{chrom}.png",
    script:
        "scripts/ibd_plots.py"
