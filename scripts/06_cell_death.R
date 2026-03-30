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

config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

plate<-config$input$plate

q_dead<-config$filtering$quantile_dead

out_dir<-file.path(base, "outputs", "Seurat")

seurat_list<-readRDS(file.path(base, "outputs", "Seurat", "intermediate", "seurat_list_filtered_doublets.rds"))

death_markers<-config$filtering$markers

filtering_criteria<-read.csv(file.path(base, "outputs", "Seurat", "QC", "filtering_criteria_S05.csv"))
filtering_criteria$celldeath_modulescore<-NA
filtering_criteria$dead_pct<-NA
filtering_criteria$cells_post_death<-NA
cell_filtering<-read.csv(file.path(base, "outputs", "Seurat", "QC", "cell_filtering_S05.csv"))

if (length(grep("ENSG000", rownames(seurat_list[[1]]))) != 0) {
  gene_names<-read.csv(file.path(base, "t2g.csv"))
  death_markers<-c(death_markers, gene_names$gene[which(gene_names$gene_name %in% death_markers)])
}


seurat_list_filtered<-seurat_list
for (i in 1:length(seurat_list)) {
  seu<-seurat_list[[i]]
  
  
  sample<-names(seurat_list)[i]
  
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(seurat_list)))
  
  cell_sample<-read.csv(file.path(base, "outputs", "Seurat", "QC", "samples", paste0(sample, "_cell_filtering_S05.csv")))
  
  death_markers<-death_markers[which(death_markers %in% rownames(seu))]
  
  if(length(death_markers) != 0) {
    if (length(which(cell_sample$status == "included")) > 100) {
      seu<-AddModuleScore(seu, features = list(death_markers), name = "CellDeath")
      meta<-seu@meta.data
      meta$orig.ident<-rownames(meta)
      
      q<-quantile(meta$CellDeath1, q_dead)
      
      x<-meta$orig.ident[which(meta$CellDeath1 > q)]
      
      cell_sample$status[which(cell_sample$cell %in% x)]<-"excluded"
      cell_sample$reason[which(cell_sample$cell %in% x)]<-"dead"
      
      cell_filtering$status[which(cell_filtering$cell %in% x)]<-"excluded"
      cell_filtering$reason[which(cell_filtering$cell %in% x)]<-"dead"
      
      check<-cell_sample$cell[which(cell_sample$status == "included")]
      
      
      singlet_idx<-which(rownames(seu@meta.data) %in% check)
      
      seu<-seu[,singlet_idx]
      
      if (length(rownames(seu@meta.data)) != length(which(rownames(seu@meta.data) %in% check))) {
        stop("Cells in filtered object do not match included cell list")
      }  
      
      filtering_criteria$celldeath_modulescore[which(filtering_criteria$sample == sample)]<-q
      filtering_criteria$dead_pct[which(filtering_criteria$sample == sample)]<-round(100*length(x)/nrow(meta),2)
      filtering_criteria$cells_post_death[which(filtering_criteria$sample == sample)]<-nrow(seu@meta.data)
      
      seurat_list_filtered[[sample]]<-seu
      
    } else {
      x<-cell_sample$cell[which(cell_sample$status == "included")]
      
      cell_sample$status[which(cell_sample$cell %in% x)]<-"excluded"
      cell_sample$reason[which(cell_sample$cell %in% x)]<-"sample_fail"
      
      cell_filtering$status[which(cell_filtering$cell %in% x)]<-"excluded"
      cell_filtering$reason[which(cell_filtering$cell %in% x)]<-"sample_fail"
      
      n<-which(names(seurat_list_filtered) == sample)
      
      seurat_list_filtered<-seurat_list_filtered[-(n)]
    }
    
    write.csv(cell_sample, file.path(out_dir, "QC", "samples", paste0(sample, "_cell_filtering.csv")), row.names = F)
    
  } else {
    print(paste("No death markers found in sample", sample))
    
    seurat_list_filtered[[sample]]<-seu
    write.csv(cell_sample, file.path(out_dir, "QC", "samples", paste0(sample, "_cell_filtering.csv")), row.names = F)
    
  }
  

  
}


write.csv(filtering_criteria, file.path(out_dir, "QC", "filtering_criteria.csv"), row.names = F)
write.csv(cell_filtering, file.path(out_dir, "QC", "cell_filtering.csv"), row.names = F)

saveRDS(seurat_list_filtered, file.path(out_dir, "intermediate", "seurat_list_filtered_doublets_dead.rds"))

