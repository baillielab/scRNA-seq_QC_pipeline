rm(list=ls())

print("Starting 02.4_filter_seurat.R")

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

plate<-config$input$plate

q_feature<-config$filtering$quantile_feature
feature_lim<-config$filtering$feature_lim
feature_max<-config$filtering$feature_lim_max

q_mt<-config$filtering$quantile_mt
mt_lim<-config$filtering$mt_lim

out_dir<-file.path(base, "outputs", "Seurat")
seurat_list<-readRDS(file.path(base, "outputs", "Seurat", "intermediate", "seurat_list.rds"))

dir.create(file.path(out_dir, "QC", "samples"), recursive = T)

filtering_criteria<-data.frame(
  sample = NA,
  total_cells = NA,
  q95_nFeature = NA,
  used_nFeature = NA,
  cells_post_nFeature = NA,
  q95_mt = NA,
  used_mt = NA,
  cells_post_mt = NA
)

cell_filtering<-data.frame(
  sample = NA,
  cell = NA,
  status = NA,
  reason = NA
)

seurat_list_filtered<-seurat_list

for(i in 1:length(seurat_list)) {
  
  seu<-seurat_list[[i]]
  
  sample<-names(seurat_list)[i]
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(seurat_list)))
  meta<-seu@meta.data
  
  meta$orig.ident<-rownames(meta)
  
  sample_cell<-data.frame(
    sample = sample,
    cell = meta$orig.ident,
    status = "included",
    reason = NA
  )
  
  total_cells<-nrow(meta)
  q95_nFeature<-quantile(meta$nFeature_RNA, q_feature)
  used_nFeature<-max(q95_nFeature, feature_lim)
  used_nFeature<-min(used_nFeature, feature_max)
  
  x<-meta$orig.ident[which(meta$nFeature_RNA < used_nFeature)]
  
  sample_cell$status[which(sample_cell$cell %in% x)]<-"excluded"
  sample_cell$reason[which(sample_cell$cell %in% x)]<-"nFeature_RNA"
  
  seu<-subset(seu, 
              subset = nFeature_RNA > used_nFeature)
  
  meta<-seu@meta.data
  
  meta$orig.ident<-rownames(meta)
  
  cells_post_nFeature<-nrow(meta)
  q95_mt<-quantile(meta$percent.mt, q_mt)
  used_mt<-min(q95_mt, mt_lim)
  
  x<-meta$orig.ident[which(meta$percent.mt > used_mt)]
  
  sample_cell$status[which(sample_cell$cell %in% x)]<-"excluded"
  sample_cell$reason[which(sample_cell$cell %in% x)]<-"percent.mt"
  
  if(length(x) != 0) {
    seu<-subset(seu, 
                subset = percent.mt < used_mt)
  }

  
  check<-sample_cell$cell[which(sample_cell$status == "included")]
  
  if (length(rownames(seu@meta.data)) != length(which(rownames(seu@meta.data) %in% check))) {
    stop("Cells in filtered object do not match included cell list")
  }    
  
  filtering_sample<-data.frame(
    sample = sample,
    total_cells = total_cells,
    q95_nFeature = q95_nFeature,
    used_nFeature = used_nFeature,
    cells_post_nFeature = cells_post_nFeature,
    q95_mt = q95_mt,
    used_mt = used_mt,
    cells_post_mt = nrow(seu@meta.data)
  )
  
  write.csv(sample_cell, file.path(out_dir, "QC", "samples", paste0(sample, "_cell_filtering_S04.csv")), row.names = F)
  
  filtering_criteria<-rbind(filtering_criteria, filtering_sample)
  cell_filtering<-rbind(cell_filtering, sample_cell)
  
  seurat_list_filtered[[sample]]<-seu
  
}

filtering_criteria<-filtering_criteria[-1,]
cell_filtering<-cell_filtering[-1,]

write.csv(filtering_criteria, file.path(out_dir, "QC", "filtering_criteria_S04.csv"), row.names = F)
write.csv(cell_filtering, file.path(out_dir, "QC", "cell_filtering_S04.csv"), row.names = F)

saveRDS(seurat_list_filtered, file.path(out_dir, "intermediate", "seurat_list_filtered.rds"))



