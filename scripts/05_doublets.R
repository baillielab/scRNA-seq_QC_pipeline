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
library(DoubletFinder)

config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

plate<-config$input$plate

method<-config$doubletFinder$findVariableFeatures$selectionMethod
features<-config$doubletFinder$findVariableFeatures$nFeatures
npcs<-config$doubletFinder$runPCA$npcs
prop_exp<-config$doubletFinder$prop_exp
pN<-config$doubletFinder$pN



out_dir<-file.path(base, "outputs", "Seurat")

seurat_list<-readRDS(file.path(base, "outputs", "Seurat", "intermediate", "seurat_list_filtered.rds"))

filtering_criteria<-read.csv(file.path(base, "outputs", "Seurat", "QC", "filtering_criteria_S04.csv"))
filtering_criteria$nExp<-NA
filtering_criteria$pK<-NA
filtering_criteria$doublet_pct<-NA
filtering_criteria$cells_post_doublet<-NA

cell_filtering<-read.csv(file.path(base, "outputs", "Seurat", "QC", "cell_filtering_S04.csv"))

seurat_list_filtered<-seurat_list

for (i in 1:length(seurat_list)) {
  seu<-seurat_list[[i]]
  sample<-names(seurat_list)[i]
  
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(seurat_list)))
  
  cell_sample<-read.csv(file.path(base, "outputs", "Seurat", "QC", "samples", paste0(sample, "_cell_filtering_S04.csv")))
  
  if (length(which(cell_sample$status == "included")) > 100) {
    seu<-NormalizeData(seu)
    seu<-FindVariableFeatures(seu, selection.method = method, nfeatures = features)
    seu<-ScaleData(seu)
    seu<-RunPCA(seu, npcs = npcs, verbose = F)
    
    nExp<-round(prop_exp * ncol(seu))
    ps<-suppressMessages(
      paramSweep(seu, PCs = 1:npcs)
    )
    ss<-summarizeSweep(ps)
    pK<-ss$pK[which.max(ss$BCreal)][1]
    pK<-as.vector(pK)
    pK<-as.numeric(pK)
    
    seu<-doubletFinder(seu, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp)
    
    meta<-seu@meta.data
    meta$orig.ident<-rownames(meta)
    
    doublet_col<-grep("DF.classifications", colnames(meta), value = T)
    dc<-which(colnames(meta) == doublet_col)
    
    x<-meta$orig.ident[which(meta[,dc] == "Doublet")]
    
    doublet_pct<-round(100*length(x)/nrow(meta),2)
    cell_sample$status[which(cell_sample$cell %in% x)]<-"excluded"
    cell_sample$reason[which(cell_sample$cell %in% x)]<-"doublet"
    
    cell_filtering$status[which(cell_filtering$cell %in% x)]<-"excluded"
    cell_filtering$reason[which(cell_filtering$cell %in% x)]<-"doublet"
    
    check<-cell_sample$cell[which(cell_sample$status == "included")]
    
    
    singlet_idx<-which(rownames(seu@meta.data) %in% check)
    
    seu<-seu[,singlet_idx]
    
    if (length(rownames(seu@meta.data)) != length(which(rownames(seu@meta.data) %in% check))) {
      stop("Cells in filtered object do not match included cell list")
    }   
    
    filtering_criteria$nExp[which(filtering_criteria$sample == sample)]<-nExp
    filtering_criteria$pK[which(filtering_criteria$sample == sample)]<-pK
    filtering_criteria$doublet_pct[which(filtering_criteria$sample == sample)]<-doublet_pct
    filtering_criteria$cells_post_doublet[which(filtering_criteria$sample == sample)]<-nrow(seu@meta.data)
    
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
  
  write.csv(cell_sample, file.path(out_dir, "QC", "samples", paste0(sample, "_cell_filtering_S05.csv")), row.names = F)
  
}

write.csv(filtering_criteria, file.path(out_dir, "QC", "filtering_criteria_S05.csv"), row.names = F)
write.csv(cell_filtering, file.path(out_dir, "QC", "cell_filtering_S05.csv"), row.names = F)

saveRDS(seurat_list_filtered, file.path(out_dir, "intermediate", "seurat_list_filtered_doublets.rds"))




