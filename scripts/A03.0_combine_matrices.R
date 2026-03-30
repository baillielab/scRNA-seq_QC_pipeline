rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)

# Ingest config
config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

input_dir <- file.path(base, config$input$data_dir)
if (is.null(input_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate

out_dirs<-file.path(base, "outputs")


# Sample names - import from run and sample sheet
sample_names<-read.table(file.path(input_dir, "samples.tsv"))

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]


  
for (i in 1:length(sample_names)) {
    # for each sample do:
    sample<-sample_names[i]
    
    print(paste("Processing ", sample, " ; ", i, " of ", length(sample_names), sep = ""))
    
    # Read in each count matrix for the sample
    mature<-readMM(file.path(input_dir, sample, "count.genes.mature.mtx"))
    nascent<-readMM(file.path(input_dir, sample, "count.genes.nascent.mtx"))
    ambiguous<-readMM(file.path(input_dir, sample, "count.genes.ambiguous.mtx"))
    
    # Extract gene names
    genes<-read.delim(file.path(input_dir, sample, "count.genes.names.txt"), header = FALSE)
    genes$name<-genes$V1
    
    # Read in busfile to extract barcodes
    bus<-fread(file.path(input_dir, sample, "sorted_bus.txt"))
    colnames(bus)<-c("barcode", "UMI", "ec", "count")
    
    barcodes<-unique(bus$barcode)
    
    # Name count matrix rows with barcodes
    rownames(mature)<-barcodes
    rownames(nascent)<-barcodes
    rownames(ambiguous)<-barcodes
    
    # Name count matrix columns with genes
    colnames(mature)<-genes$name
    colnames(nascent)<-genes$name
    colnames(ambiguous)<-genes$name
    
    # Add counts across all matricies
    counts_total <- mature + nascent + ambiguous
    umis_per_cell <- rowSums(counts_total)
    
    
    # Save the combined (unfiltered) matrix for log and the filtered matrix to use in subsequent steps
    saveRDS(counts_total, file.path(input_dir, sample, "combined_matrix.rds"))
    
}
