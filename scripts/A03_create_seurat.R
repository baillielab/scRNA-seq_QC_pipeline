rm(list=ls())

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

input_dir <- file.path(base, config$input$data_dir)
if (is.null(input_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate

out_dirs<-file.path(base, "outputs")


# Sample names - import from run and sample sheet
sample_names<-read.table(file.path(input_dir, "samples.tsv"))
gene_names<-read.csv(paste(base, "/t2g.csv", sep = ""))

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]


# ALL count.genes.names.txt FILES IDENTICAL - ONLY DO ONCE
# Read in gene ids to add to matrix
genes<-read.delim(file.path(input_dir, sample_names[1], "count.genes.names.txt"), header = FALSE)
genes$name<-genes$V1

MT<-gene_names[grep("^MT-", gene_names$gene_name),]
# To use gene names instead of IDs:

# for (i in 1:length(genes$V1)) {
#   if (length(which(gene_names$gene == genes$V1[i])) != 0) {
#     genes$name[i]<-unique(gene_names$gene_name[which(gene_names$gene == genes$V1[i])])
#   }
# }


### FUNCTIONS ###
# Takes the count matrix count.genes.mature.mtx, and turns it into an annotated count matrix using gene ids from count.genes.names.txt and gene names from gene_names input table as
# row (gene) names and extracts barcodes from busfile to use as column (barcode) names
# NOT NEEDED ANY MORE - REPLACED WITH NEW_KNEE; CODE LEFT IN BUT FUNCTION NEVER CALLED
make_annotated_matrix<-function(sample) {
  # Read in count matrix (mature), annotate gene names
  counts<-readMM(file.path(input_dir, sample, "count.genes.mature.mtx"))
  counts<-t(counts)
  
  rownames(counts)<-make.unique(genes$name)

  # Read in busfile to extract barcodes
  bus<-fread(file.path(input_dir, sample, "sorted_bus.txt"))
  colnames(bus)<-c("barcode", "UMI", "ec", "count")
  
  barcodes<-unique(bus$barcode)
  
  if(length(barcodes) != ncol(counts)) {
    print("Number of barcodes does not match number of columns in count matrix!")
  } else {
    colnames(counts)<-barcodes
  }
  
  barcode_file<-file.path(input_dir, sample, "QC", "barcode_filter.csv")
  
  if(file.exists(barcode_file)) {
    message("Using filtered barcodes for ", sample)
    
    filtered_barcodes<-fread(barcode_file, header = FALSE)$V1
    counts<-counts[,colnames(counts) %in% filtered_barcodes]
  }
  return(counts)
}

# Takes the count matrix and creates a seurat object, calculates % mitochondrial genes, outputs summary stats (to screen - saved in next function), creates and saves violin plots of
# nFeature_RNA, nCount_RNA and percent.mt to use in QC filtering, also creates and saves counts vs genes plot. Saves in Seurat directory in sample directory (must exist)
# Returns Seurat object. Flowcell argument for unique labelling of samples.
multisample_seurat_QC<-function(count_matrix, sample, output1, output2) {

  unique_name<-sample
  
  seurat_obj<-CreateSeuratObject(counts = count_matrix, project = sample)
  
  seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj, features = MT$gene)
  
  summary_stats <- data.frame(
    sample = sample,
    n_cells = ncol(seurat_obj),
    median_umis = median(seurat_obj$nCount_RNA),
    median_genes = median(seurat_obj$nFeature_RNA),
    median_percent_mt = median(seurat_obj$percent.mt, na.rm = T)
  )
  
  if (is.na(sum(seurat_obj$percent.mt))) {
    p1<- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3) + ggtitle(paste("QC for", sample))
  } else {
    p1<- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + ggtitle(paste("QC for", sample))
  }
  ggsave(file.path(out_dirs, "QC", "samples", sample, output1), width = 10, height = 6)
  
  p2<-FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(paste("Counts vs Genes for", sample))
  ggsave(file.path(out_dirs, "QC", "samples", sample, output2), width = 10, height = 6)
  
  return(seurat_obj)

}



### EXTRACT DATA ###
# Loop over all samples across both flowcells and create annotated count matrices, then create seurat objects from these (with samples differentiated by flowcell) and combine all in a list
seurat_list<-list()

dir.create(file.path(out_dirs, "Seurat", "Samples"), recursive = T)
    
for (i in 1:length(sample_names)) {
      sample<-sample_names[i]
      print(paste("Processing sample", sample, "number", i, "of", length(sample_names)))
      
      dir.create(file.path(base, "outputs", "QC", "samples", sample), recursive = T)
      # Run the functions made above
      #counts<-make_annotated_matrix(sample)
      counts<-readRDS(file.path(input_dir, sample, "combined_matrix.rds"))
      counts<-t(counts)
      
      filtered_barcodes<-read.csv(file.path(out_dirs, "QC", "samples", sample, "barcode_keep_filter.csv"))
      counts<-counts[,colnames(counts) %in% filtered_barcodes$barcode]
      
      
      seurat_obj<-multisample_seurat_QC(counts, sample, "seurat_violns.png", "counts_vs_genes.png")
      
      # Adjust for same samples on different flowcells - need unique names for list
      unique_name<-sample
      
      # Add flowcell column to seurat object so they can be differentiated later
      seurat_obj$plate<-plate
      
      # Save sample specific seurat object
      saveRDS(seurat_obj, file.path(out_dirs, "Seurat", "Samples", paste0(sample, ".rds")))
      
      # Add sample specif seurat object to list of all seurat objects
      seurat_list[[unique_name]]<-seurat_obj
} 
dir.create(file.path(out_dirs, "Seurat", "intermediate"), recursive = T)
saveRDS(seurat_list, file.path(out_dirs, "Seurat", "intermediate", "seurat_list.rds"))

# Create table of QC stats across all seurat objects in list
qc_stats <- lapply(names(seurat_list), function(sample) {
  data.frame(
    sample = sample,
    nCount_RNA = seurat_list[[unique_name]]$nCount_RNA,
    nFeature_RNA = seurat_list[[unique_name]]$nFeature_RNA,
    percent.mt = seurat_list[[unique_name]]$percent.mt
  )
}) %>% bind_rows()



### ANALYSE DATA ###

# View QC metric plots to identify suitable cutoffs for filtering 
count<-ggplot(qc_stats, aes(x=nCount_RNA, fill = sample)) + 
  geom_density(alpha=0.3) + 
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); count
ggsave(file.path(out_dirs, "QC", "nCount_RNA.png"), count, bg="white", width = 10, height = 6)

feature<-ggplot(qc_stats, aes(x=nFeature_RNA, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); feature
ggsave(file.path(out_dirs, "QC", "nFeature_RNA.png"), count, bg="white", width = 10, height = 6)

mt<-ggplot(qc_stats, aes(x=percent.mt, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); mt
ggsave(file.path(out_dirs, "QC", "percent_mt.png"), count, bg="white", width = 10, height = 6)
