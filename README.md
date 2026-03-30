# snRNA-seq QC Pipeline (Snakemake)

This repository contains the Snakemake workflow and scripts for running the ODAP single-cell / single-nucleus RNA-seq QC pipeline.

## Overview

This pipeline performs quality control on data processed with:

* Parse Splitpipe
* Arcas

It removes low-quality cells, doublets, and dead cells, and produces cleaned count matrices and annotations for downstream analysis.

The workflow is implemented in **Snakemake** for reproducibility and scalability on HPC systems.

## Pipeline Summary

The pipeline includes:

* Run-level statistics (reads, barcoding, mapping)
* Initial QC and knee-point filtering
* Seurat object creation
* Filtering based on:

  * number of detected genes (`nFeature_RNA`)
  * mitochondrial content
* Doublet detection (DoubletFinder)
* Cell death detection (marker-based)
* Per-sample and merged outputs

## Input Support

Two input types are supported:

* **Splitpipe**
* **Arcas** (includes preprocessing step)

Both converge to a shared downstream QC and analysis workflow.

## Usage

This repository is designed to be run on an HPC environment using Slurm.

Basic steps:


Edit configuration:

```bash
nano config.yaml
```

Run pipeline:

```bash
sbatch run_pipeline.slurm
```

## Configuration

Key parameters are defined in `config.yaml`, including:

* input type (`splitpipe` or `arcas`)
* data directory (containing arcas / splitpipe outputs)
* QC filtering thresholds

Filtering is adaptive and designed to handle both scRNA-seq and **snRNA-seq (low UMI)** datasets.

## Outputs

* Intermediate QC metrics and Seurat objects (`outputs/`)
* Final cleaned count matrices and annotations (`results/`)
* Per-sample and merged outputs
* Summary reports and figures

## Documentation

Full user guide (including setup, data structure, and parameter details) is available in the ODAP documentation:
[Full documentation and usage instructions](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/Analysis-in-ODAP/scRNA-Seq-QC-Pipeline)

## Notes

* Designed for automated, multi-sample QC workflows
* Supports heterogeneous datasets and sequencing depths
* Includes safeguards to prevent over-filtering in shallow snRNA-seq data

## Requirements

These are all available in ODAP containers.

* Snakemake
* Python
* bustools
* Slurm (for HPC execution)
* R
  * yaml
  * jsonlite
  * ggplot2
  * dplyr
  * tidyr
  * data.table
  * Matrix
  * Seurat
  * DoubletFinder

## Author

**Kathryn Campbell**
on behalf of The ODAP Team

30th March 2026

