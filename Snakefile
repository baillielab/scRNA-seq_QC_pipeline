# ----------------------------
# Snakefile
# ----------------------------

import os
from pathlib import Path

configfile: "config.yaml"

BASE_DIR = Path(config["input"]["base_dir"])
DATA_DIR = BASE_DIR / config["input"]["data_dir"]
INPUT_TYPE = config["input"]["type"]

# List of samples
SAMPLES = [
    p.name 
    for p in DATA_DIR.iterdir() 
    if (
        p.is_dir() and 
        p.name != "runs"
        and not p.name.startswith(".")
    )
]

# ----------------------------
# Rule all
# ----------------------------

rule all:
    input:
        # Mapping
        str(BASE_DIR / "outputs" / "QC" / "mapping_stats" / "plate_mapping_stats.png"),

        # Initial QC
        str(BASE_DIR / "outputs" / "QC" / "knee_inclusion.png"),

        # Seurat steps
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered_doublets_dead.rds"),

        # Summarise
        str(BASE_DIR / "results" / "Plate" / "matrix.mtx.gz")

# ----------------------------------------------------------------------------------------------------------
#            SPLITPIPE
# ----------------------------------------------------------------------------------------------------------
if INPUT_TYPE == "splitpipe":

    # ----------------------------
    # S01 Mapping
    # ----------------------------
    rule mapping_SP:
        input:
            str(DATA_DIR / "sample-list.txt"),
            expand(str(DATA_DIR / "{sample}" / "report" / "analysis_summary.csv"), sample=SAMPLES),
            config="config.yaml"
        output:
            str(BASE_DIR / "outputs" / "QC" / "mapping_stats" / "plate_mapping_stats.png")
        shell:
            "Rscript scripts/S01_mapping.R {input.config}"


    # ----------------------------
    # S02 Initial_QC
    # ----------------------------
    rule initial_QC_SP:
        input:
            expand(str(DATA_DIR / "{sample}" / "DGE_unfiltered" / "cell_metadata.csv"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "DGE_unfiltered" / "count_matrix.mtx"), sample=SAMPLES),
            str(DATA_DIR / "sample-list.txt"),
            config="config.yaml"
        output:
            inclusion=str(BASE_DIR / "outputs" / "QC" / "knee_inclusion.png"),
            filter=expand(str(BASE_DIR / "outputs" / "QC" / "samples" / "{sample}" / "barcode_keep_filter.csv"), sample=SAMPLES)
        shell:
            "Rscript scripts/S02_initial_QC.R {input.config}"


    # ----------------------------
    # S03 Create Seurat
    # ----------------------------
    rule create_seurat_SP:
        input:
            expand(str(DATA_DIR / "{sample}" / "DGE_unfiltered" / "cell_metadata.csv"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "DGE_unfiltered" / "all_genes.csv"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "DGE_unfiltered" / "count_matrix.mtx"), sample=SAMPLES),
            str(DATA_DIR / "sample-list.txt"),
            expand(str(BASE_DIR / "outputs" / "QC" / "samples" / "{sample}" / "barcode_keep_filter.csv"), sample=SAMPLES),
            config="config.yaml"
        output:
            seurat_list=str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list.rds")
        shell:
            "Rscript scripts/S03_create_seurat.R {input.config}"

# ----------------------------------------------------------------------------------------------------------
#            ARCAS
# ----------------------------------------------------------------------------------------------------------
if INPUT_TYPE == "arcas":

    # ----------------------------
    # A00 Bustools
    # ----------------------------
    rule bustools_A:
        input:
            str(DATA_DIR / "{sample}" / "kallisto.sorted_bus")
        output:
            str(DATA_DIR / "{sample}" / "sorted_bus.txt")
        shell:
            "bustools text -o {output} {input}"

    # ----------------------------
    # A01 Mapping
    # ----------------------------
    RUN_SUMMARIES = list((DATA_DIR / "runs").glob("*.extract.summary.json"))

    rule mapping_A:
        input:
            str(DATA_DIR / "run_sheet.tsv"),
            str(DATA_DIR / "samples.tsv"),
            RUN_SUMMARIES,
            expand(str(DATA_DIR / "{sample}" / "kallisto.run_info.json"), sample=SAMPLES),
            config="config.yaml"
        output:
            str(BASE_DIR / "outputs" / "QC" / "mapping_stats" / "plate_mapping_stats.png")
        shell:
            "Rscript scripts/A01_mapping.R {input.config}"


    # ----------------------------
    # A02 Initial_QC
    # ----------------------------
    rule initial_A:
        input:
            expand(str(DATA_DIR / "{sample}" / "sorted_bus.txt"), sample=SAMPLES),
            str(DATA_DIR / "samples.tsv"),
            config="config.yaml"
        output:
            inclusion=str(BASE_DIR / "outputs" / "QC" / "knee_inclusion.png"),
            filter=expand(str(BASE_DIR / "outputs" / "QC" / "samples" / "{sample}" / "barcode_keep_filter.csv"), sample=SAMPLES)
        shell:
            "Rscript scripts/A02_initial_QC.R {input.config}"


    # ----------------------------
    # A03.0 Combine Matrices
    # ----------------------------
    rule combine_matrices_A:
        input:
            expand(str(DATA_DIR / "{sample}" / "count.genes.names.txt"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "sorted_bus.txt"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "count.genes.mature.mtx"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "count.genes.ambiguous.mtx"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "count.genes.nascent.mtx"), sample=SAMPLES),
            str(DATA_DIR / "samples.tsv"),
            config="config.yaml"
        output:
            seurat_list=expand(str(DATA_DIR / "{sample}" / "combined_matrix.rds"), sample=SAMPLES)
        shell:
            "Rscript scripts/A03.0_combine_matrices.R {input.config}"

    # ----------------------------
    # A03 Create Seurat
    # ----------------------------
    rule create_seurat_A:
        input:
            expand(str(DATA_DIR / "{sample}" / "count.genes.names.txt"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "sorted_bus.txt"), sample=SAMPLES),
            expand(str(DATA_DIR / "{sample}" / "combined_matrix.rds"), sample=SAMPLES),
            str(DATA_DIR / "samples.tsv"),
            expand(str(BASE_DIR / "outputs" / "QC" / "samples" / "{sample}" / "barcode_keep_filter.csv"), sample=SAMPLES),
            config="config.yaml"
        output:
            seurat_list=str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list.rds")
        shell:
            "Rscript scripts/A03_create_seurat.R {input.config}"


# ----------------------------------------------------------------------------------------------------------
#            CONVERGE
# ----------------------------------------------------------------------------------------------------------

# ----------------------------
# 04 Filter Seurat
# ----------------------------
rule filter_seurat:
    input:
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list.rds"),
        config="config.yaml"
    output:
        temp(expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering_S04.csv"), sample=SAMPLES)),
        all_filter=temp(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering_S04.csv")),
        filtering_criteria=temp(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria_S04.csv")),
        filter_list=str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered.rds")
    shell:
        "Rscript scripts/04_filter_seurat.R {input.config}"

# ----------------------------
# 05 Filter Doublets
# ----------------------------
rule doublets:
    input:
        expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering_S04.csv"), sample=SAMPLES),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering_S04.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria_S04.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered.rds"),
        config="config.yaml"
    output:
        temp(expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering_S05.csv"), sample=SAMPLES)),
        temp(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering_S05.csv")),
        temp(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria_S05.csv")),
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered_doublets.rds")
    shell:
        "Rscript scripts/05_doublets.R {input.config}"

# ----------------------------
# 06 Filter Death
# ----------------------------
rule death:
    input:
        expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering_S05.csv"), sample=SAMPLES),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering_S05.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria_S05.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered_doublets.rds"),
        str(BASE_DIR / "t2g.csv"),
        config="config.yaml"
    output:
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered_doublets_dead.rds"),
        expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering.csv"), sample=SAMPLES)
    shell:
        "Rscript scripts/06_cell_death.R {input.config}"

# ----------------------------
# S07 Summarise
# ----------------------------
rule summarise:
    input:
        expand(str(BASE_DIR / "outputs" / "Seurat" / "QC" / "samples" / "{sample}_cell_filtering.csv"), sample=SAMPLES),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "cell_filtering.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "QC" / "filtering_criteria.csv"),
        str(BASE_DIR / "outputs" / "Seurat" / "intermediate" / "seurat_list_filtered_doublets_dead.rds"),
        config="config.yaml"
    output:
        plate_counts=str(BASE_DIR / "results" / "Plate" / "matrix.mtx.gz")
    shell:
        "Rscript scripts/07_summarise_QC.R {input.config}"

