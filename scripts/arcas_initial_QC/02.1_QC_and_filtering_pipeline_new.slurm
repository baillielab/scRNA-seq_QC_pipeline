#!/bin/bash                                                                     

# === SLURM SETTINGS ===
#SBATCH --job-name=QC_and_filtering4                                                      
#SBATCH --time=72:00:00                                                         
#SBATCH --mem=64G                                                               
#SBATCH --output=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.out           
#SBATCH --error=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.err            
#SBATCH --ntasks=1                                                              
#SBATCH --cpus-per-task=16  

# =========================


# === USER SETTINGS ===

# Set base paths                                                                

# File path to where analysis will take place, same as datadir in config file
BASE=/mnt/odap-beegfs/a008/pilot_pilot_new



# =================================================

# Set base paths                                                                
CONT=/mnt/odap-beegfs/software/singularity-images

#Set specific paths                                                             
CONT_ARCAS=${CONT}/arcas-container-v2.sif
CONT_R=${CONT}/mmc-rstudio-container-v3.sif
SNAKEFILE=${BASE}/scripts/02_QC_and_filtering/02.1.1_Snakefile_bustools

# Set config arguments
SUBDIR=MMC
IN_DIR=outputs/arcas_pipeline/* 
OUTDIR=outputs/QC_and_filtering_new

# SAMPLES
SAMPLE_SHEET=${BASE}/data/samples.tsv
RUN_SHEET=${BASE}/data/run_sheet.tsv


# === QC AND FILTERING ===

# Create t2g
#@@@singularity exec --bind ${BASE},${CONT} ${CONT_R} /usr/bin/Rscript ${BASE}/scripts/02_QC_and_filtering/02.1.0_t2g_generate.R "${BASE}"
 
# Convert kallisto.sorted_bus to sorted_bus.txt - readable text file 
#@@@singularity exec --bind ${BASE},${CONT} ${CONT_ARCAS} snakemake -s ${SNAKEFILE} --config base="${BASE}" in_dir="${IN_DIR}" subdir="${SUBDIR}" --cores 32 --use-conda --use-singularity


# Extract number (and percentage) of total reads, barcoded reads, and mapped reads 
#@@@mkdir -p "${BASE}/outputs/QC_and_filtering_new"
#@@@singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02_QC_and_filtering/02.1.2_mapping_stats_command.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}" "${SAMPLE_SHEET}" "${RUN_SHEET}"

# Extract sequencing and run QC metrics from sorted_bus.txt files
#@@@singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02_QC_and_filtering/02.1.3_busfile_QC_command.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}" "${SAMPLE_SHEET}"

singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02_QC_and_filtering/new_knee.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}" "${SAMPLE_SHEET}"


# Create seurat objects for each sample from count matricies
# QC the seurat objects and create figures, filter for nFeature_RNA, n_Count_RNA and percent.mt
mkdir -p "${BASE}/outputs/seurat_processing_new/intermediate"
mkdir -p "${BASE}/outputs/seurat_processing_new/figures/qc"
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02_QC_and_filtering/02.1.4_seurat_QC_command_new.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}" "${SAMPLE_SHEET}"

# Using filtered seurat data, calculate number of usable cells remaining and reads per cell
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02_QC_and_filtering/02.1.5_cells_and_reads_command_new.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}" "${SAMPLE_SHEET}"

