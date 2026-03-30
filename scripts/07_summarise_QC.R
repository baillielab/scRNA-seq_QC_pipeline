rm(list=ls())

print("Starting 02.7_summarise_QC.R")


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

config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

plate<-config$input$plate

input_dir<-file.path(base, "outputs", "Seurat")

out_dir<-file.path(base, "results")

seurat_list<-readRDS(file.path(input_dir, "intermediate", "seurat_list_filtered_doublets_dead.rds"))
filtering_criteria<-read.csv(file.path(input_dir, "QC", "filtering_criteria.csv"))
cell_filtering<-read.csv(file.path(input_dir, "QC", "cell_filtering.csv"))

postQC_summary<-data.frame(
  sample = NA,
  plate = NA,
  cells = NA,
  median_genes_per_cell = NA,
  median_UMIs_per_cell = NA,
  median_seq_sat = NA,
  median_mt_pct = NA
)

for (i in 1:length(seurat_list)) {
  seu<-seurat_list[[i]]
  
  sample<-names(seurat_list)[i]
  
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(seurat_list)))
  
  meta<-seu@meta.data
  meta$orig.ident<-rownames(meta)
  
  postQC_sample<-data.frame(
    sample = sample,
    plate = plate,
    cells = nrow(meta),
    median_genes_per_cell = median(meta$nFeature_RNA),
    median_UMIs_per_cell = median(meta$nCount_RNA),
    median_seq_sat = median(meta$nFeature_RNA / meta$nCount_RNA),
    median_mt_pct = median(meta$percent.mt)
  )
  
  meta$tissue<-NA
  
  if (length(grep("SKN", sample)) != 0) {
    meta$tissue<-"Skin"
  } else if (length(grep("BRN", sample)) != 0) {
    meta$tissue<-"Brain"
  } else if (length(grep("LGS", sample)) != 0) {
    meta$tissue<-"Lung"
  } else if (length(grep("BLD", sample)) != 0) {
    meta$tissue<-"Blood"
  } else {
    meta$tissue<-"Unknown"
  }
  
  
  seurat_list[[i]]@meta.data<-meta
  
  meta$donor_id<-strsplit(sample, "_")[[1]][1]
  meta$condition<-strsplit(sample, "_")[[1]][2]
  meta$condition_extra<-gsub(".*_", "", sample)
  
  counts<-GetAssayData(seu, assay = "RNA", layer = "counts")
  
  dir.create(file.path(out_dir, "Samples", sample), recursive = T)
  writeMM(counts, file.path(out_dir, "Samples", sample, "matrix.mtx"))
  
  write.table(
    colnames(counts),
    file = file.path(out_dir, "Samples", sample, "barcodes.tsv"),
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  
  write.table(
    rownames(counts),
    file = file.path(out_dir, "Samples", sample, "genes.tsv"),
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  
  write.table(
    meta,
    file = file.path(out_dir, "Samples", sample, "metadata.tsv"),
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  
  system(paste("gzip", file.path(out_dir, "Samples", sample, "matrix.mtx")))
  system(paste("gzip", file.path(out_dir, "Samples", sample, "barcodes.tsv")))
  system(paste("gzip", file.path(out_dir, "Samples", sample, "genes.tsv")))
  
  postQC_summary<-rbind(postQC_summary, postQC_sample)
}

postQC_summary<-postQC_summary[-1,]

cells<-ggplot(filtering_criteria, aes(x = sample, y = cells_post_death)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0, 5000)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Plate", plate, "post QC cell counts"), 
       y = "Number of cells",
       x = NULL
  ); cells
ggsave(file.path(out_dir, "postQC_cell_counts.png"), cells, bg = "white", width = 10, height=6)

genes<-ggplot(postQC_summary, aes(x = sample, y = median_genes_per_cell)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0, 4000)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Plate", plate, "post QC median genes per cell"), 
       y = "Number of genes",
       x = NULL
  ); genes
ggsave(file.path(out_dir, "postQC_gene_counts.png"), genes, bg = "white", width = 10, height=6)


doublets<-ggplot(filtering_criteria, aes(x = sample, y = doublet_pct)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0, 50)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Plate", plate, "% doublets"), 
       y = "% doublets",
       x = NULL
  ); doublets
ggsave(file.path(out_dir, "pct_doublets.png"), doublets, bg = "white", width = 10, height=6)


death<-ggplot(filtering_criteria, aes(x = sample, y = dead_pct)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0, 10)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = paste("Plate", plate, "% dead"), 
       y = "% dead",
       x = NULL
  ); death
ggsave(file.path(out_dir, "pct_dead.png"), death, bg = "white", width = 10, height=6)

write.csv(filtering_criteria, file.path(out_dir, "filtering_criteria.csv"), row.names = F)
write.csv(cell_filtering, file.path(out_dir, "cell_filtering.csv"), row.names = F)

print("Merging...")

counts_list<-lapply(seurat_list, function(x) GetAssayData(x, assay = "RNA", layer = "counts"))
meta_list<-lapply(seurat_list, function(x) x@meta.data)
meta_list<-lapply(meta_list, function(df) {
  pANN_col<-grep("^pANN", colnames(df), value = T)
  colnames(df)[which(colnames(df) == pANN_col)]<-"pANN"
  
  DF_col<-grep("^DF", colnames(df), value = T)
  colnames(df)[which(colnames(df) == DF_col)]<-"DF.classifications"
  
  return(df)
})

for (i in seq_along(counts_list)) {
  sample_id<-names(seurat_list)[i]
  colnames(counts_list[[i]])<-paste0(sample_id, "_", colnames(counts_list[[i]]))
  rownames(meta_list[[i]])<-colnames(counts_list[[i]])
  meta_list[[i]]$sample_id<-sample_id
}

combined_counts<-do.call(cbind, counts_list)

combined_meta<-do.call(rbind, meta_list)

dir.create(file.path(out_dir, "Plate"))

writeMM(combined_counts, file.path(out_dir, "Plate", "matrix.mtx"))

write.table(
  colnames(combined_counts),
  file = file.path(out_dir, "Plate", "barcodes.tsv"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)

write.table(
  rownames(combined_counts),
  file = file.path(out_dir, "Plate", "genes.tsv"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)

write.table(
  combined_meta,
  file = file.path(out_dir, "Plate", "metadata.tsv"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

system(paste("gzip", file.path(out_dir, "Plate", "matrix.mtx")))
system(paste("gzip", file.path(out_dir, "Plate", "barcodes.tsv")))
system(paste("gzip", file.path(out_dir, "Plate", "genes.tsv")))

print("Complete!")
