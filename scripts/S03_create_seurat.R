rm(list=ls())

print("Starting 02.3_create_seurat.R")

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(data.table)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)

# Ingest config
config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

in_dir <- file.path(base, config$input$data_dir)
if (is.null(in_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate

out_dir<-file.path(base, "outputs")

sample_names<-fread(file.path(in_dir, "sample-list.txt"), header = F)
sample_names<-sample_names$V1


### FUNCTIONS ###
# Takes the count matrix count.genes.mature.mtx, and turns it into an annotated count matrix using gene ids from count.genes.names.txt and gene names from gene_names input table as
# row (gene) names and extracts barcodes from busfile to use as column (barcode) names
make_annotated_matrix<-function(sample) {
  # Read in count matrix (mature), annotate gene names
  counts<-readMM(file.path(in_dir, sample, "DGE_unfiltered", "count_matrix.mtx"))
  counts<-t(counts)
  
  genes<-read.csv(file.path(in_dir, sample, "DGE_unfiltered", "all_genes.csv"))

  rownames(counts)<-make.unique(genes$gene_name)
  
  # Read in cell_metadata to extract barcodes
  bus<-read.csv(file.path(in_dir, sample, "DGE_unfiltered", "cell_metadata.csv"))

  barcodes<-unique(bus$bc_wells)
  
  if(length(barcodes) != ncol(counts)) {
    print("Number of barcodes does not match number of columns in count matrix!")
  } else {
    colnames(counts)<-barcodes
  }
  
  return(counts)
}

# Takes the count matrix and creates a seurat object, calculates % mitochondrial genes, outputs summary stats (to screen - saved in next function), creates and saves violin plots of
# nFeature_RNA, nCount_RNA and percent.mt to use in QC filtering, also creates and saves counts vs genes plot. Saves in Seurat directory in sample directory (must exist)
# Returns Seurat object. Flowcell argument for unique labelling of samples.
multisample_seurat_QC<-function(count_matrix, sample) {

  unique_name<-(sample)
  
  seurat_obj<-CreateSeuratObject(counts = count_matrix, project = sample)
  
  seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  summary_stats <- data.frame(
    sample = sample,
    n_cells = ncol(seurat_obj),
    median_umis = median(seurat_obj$nCount_RNA),
    median_genes = median(seurat_obj$nFeature_RNA),
    median_percent_mt = median(seurat_obj$percent.mt, na.rm = T)
  )
  
  print(summary_stats)
  
  df<-seurat_obj@meta.data %>%
    select(nFeature_RNA, nCount_RNA, percent.mt) %>%
    mutate(cell = rownames(seurat_obj@meta.data)) %>%
    pivot_longer(cols = -cell, names_to = "feature", values_to = "value") %>%
    mutate(dummy = 1)
  
  p1<-ggplot(df, aes(x = factor(dummy), y = value)) +
    geom_violin(fill = "skyblue", alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.1, alpha = 0.3) +
    facet_wrap(~feature, scales = "free_y") +
    theme_minimal() +
    ylab("Counts / Percent") +
    xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.length.x = element_blank())
  ggsave(file.path(out_dir, "QC", "samples", sample, "seurat_violins.png"), p1, width = 10, height = 6)
  
  p2<-FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(paste("Counts vs Genes for", sample))
  ggsave(file.path(out_dir, "QC","samples", sample, "counts_vs_genes.png"), p2, width = 10, height = 6)
  
  return(seurat_obj)
  
}

### EXTRACT DATA ###
# Loop over all samples across both flowcells and create annotated count matrices, then create seurat objects from these (with samples differentiated by flowcell) and combine all in a list
seurat_list<-list()


for (i in 1:length(sample_names)) {
    sample<-sample_names[i]
    print(paste("Processing sample", sample, "number", i, "of", length(sample_names)))
    
    # Make sample specific output directory
    dir.create(file.path(out_dir, "Seurat"))
    dir.create(file.path(out_dir, "Seurat", "Samples"))
    
    # Run the functions made above
    counts<-make_annotated_matrix(sample)
    
    filtered_barcodes<-read.csv(file.path(out_dir, "QC", "samples", sample, "barcode_keep_filter.csv"))
    counts<-counts[,colnames(counts) %in% filtered_barcodes$barcode]
    
    seurat_obj<-multisample_seurat_QC(counts, sample)
    
    # Adjust for same samples on different flowcells - need unique names for list
    unique_name<-paste(sample)
    
    # Add flowcell column to seurat object so they can be differentiated later
    seurat_obj$plate<-plate
    
    # Save sample specific seurat object
    saveRDS(seurat_obj, file.path(out_dir, "Seurat", "Samples", paste0(sample, ".rds")))
    
    # Add sample specif seurat object to list of all seurat objects
    seurat_list[[unique_name]]<-seurat_obj
  } 

dir.create(file.path(out_dir, "Seurat", "intermediate"))
saveRDS(seurat_list, file.path(out_dir, "Seurat", "intermediate", "seurat_list.rds"))

# Create table of QC stats across all seurat objects in list
qc_stats <- lapply(names(seurat_list), function(sample) {
  data.frame(
    sample = sample,
    nCount_RNA = seurat_list[[unique_name]]$nCount_RNA,
    nFeature_RNA = seurat_list[[unique_name]]$nFeature_RNA,
    percent.mt = seurat_list[[unique_name]]$percent.mt
  )
}) %>% bind_rows()

dir.create(file.path(out_dir, "QC"))
write.csv(qc_stats, file.path(out_dir, "QC", "qc_stats.csv"), row.names = F)

### ANALYSE DATA ###

# View QC metric plots to identify suitable cutoffs for filtering 
count<-ggplot(qc_stats, aes(x=nCount_RNA, fill = sample)) + 
  geom_density(alpha=0.3) + 
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); count
ggsave(file.path(out_dir, "QC", "nCount_RNA.png"), count, bg="white", width = 10, height = 6)

feature<-ggplot(qc_stats, aes(x=nFeature_RNA, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); feature
ggsave(file.path(out_dir, "QC", "nFeature_RNA.png"), count, bg="white", width = 10, height = 6)

mt<-ggplot(qc_stats, aes(x=percent.mt, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); mt
ggsave(file.path(out_dir, "QC", "percent_mt.png"), count, bg="white", width = 10, height = 6)


