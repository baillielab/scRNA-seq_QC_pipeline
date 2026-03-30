#!/bin/bash                                                                     

# === SLURM SETTINGS ===
#SBATCH --job-name=QC_and_filtering_2                                                      
#SBATCH --time=72:00:00                                                         
#SBATCH --mem=64G                                                               
#SBATCH --output=/mnt/odap-beegfs/a008/Kat/logs/%x/%x_%j.out           
#SBATCH --error=/mnt/odap-beegfs/a008/Kat/logs/%x/%x_%j.err            
#SBATCH --ntasks=1                                                              
#SBATCH --cpus-per-task=16  

# =========================


# === USER SETTINGS ===

# Set base paths                                                                

# File path to where analysis will take place, same as datadir in config file
ROOT=/mnt/odap-beegfs/a008
BASE=${ROOT}/Kat/plate_2
PLATE=2

# =================================================

# Set base paths                                                                
CONT=/mnt/odap-beegfs/software/singularity-images

#Set specific paths                                                             
CONT_R=${CONT}/mmc-rstudio-container-v7.sif

# Set SP output path
SPLITPIPE=/mnt/odap-beegfs/a008/exploratory/plate2/outputs/01_splitpipe/splitpipe_outputs_new/combined_new

# === QC AND FILTERING ===

# Extract number (and percentage) of total reads, barcoded reads, and mapped reads 
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.1_mapping_stats.R "${BASE}" "${SPLITPIPE}" "${BASE}/outputs" "${PLATE}"

# Extract sequencing and run QC metrics
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.2_initial_QC.R "${BASE}" "${SPLITPIPE}" "${BASE}/outputs/QC" "${PLATE}"

# Create seurat objects for each sample from count matricies
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.3_create_seurat.R "${BASE}" "${SPLITPIPE}" "${BASE}/outputs" "${PLATE}"

# Filter Seurat objects based on nFeature_RNA and percent.mt
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.4_filter_seurat.R "${BASE}" "${BASE}/outputs/Seurat/intermediate" "${BASE}/outputs/Seurat" "${PLATE}"

# Filter Seurat objects based on doublets
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.5_doublets.R "${BASE}" "${BASE}/outputs/Seurat/intermediate" "${BASE}/outputs/Seurat" "${PLATE}"

# Filter Seurat objects based on cell death
#singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.6_cell_death.R "${BASE}" "${BASE}/outputs/Seurat/intermediate" "${BASE}/outputs/Seurat" "${PLATE}"

# Summarise QC, make report figures and create merged object
singularity exec --bind ${ROOT},${CONT} ${CONT_R} Rscript ${BASE}/scripts/02-SP_QC_and_filtering/02.7_summarise_QC.R "${BASE}" "${BASE}/outputs/Seurat" "${BASE}/results" "${PLATE}"



