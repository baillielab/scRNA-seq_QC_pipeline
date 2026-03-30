rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)

# Ingest config
config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

input_dir <- file.path(base, config$input$data_dir)
if (is.null(input_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate

out_dirs<-file.path(base, "outputs", "QC")
dir.create(file.path(out_dirs, "mapping_stats"), recursive = T)

# fastq and sample names - import from run and sample sheet
fastqs<-read.delim(file.path(input_dir, "run_sheet.tsv"))
sample_names<-read.table(file.path(input_dir, "samples.tsv"))

# Fix and standardise fastq names 
for (i in 1:length(fastqs$run_id)) {
  fastqs$run_id[i]<-strsplit(fastqs$run_id[i], "[.]")[[1]][1]
}

fastqs<-fastqs$run_id

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]


### EXTRACT DATA ###
# Set up empty dataframes to populate
mapping_stats<-data.frame(plate = plate, 
                          n_processed = NA, 
                          n_assigned = NA, 
                          n_pseudoaligned = NA,
                          percent_barcoded = NA,
                          percent_all_mapped = NA,
                          percent_bc_mapped = NA
)

sample_stats<-data.frame(sample = NA,
                         plate = plate,
                         assigned = NA,
                         pseudoaligned = NA,
                         percent_bc_mapped = NA)

processed<-0
assigned<-0
pseudoaligned<-0
  
for (i in 1:length(fastqs)) {
    # Process each fastq 
    fastq<-fastqs[i]
    
    print(paste("Processing sample", fastq, "number", i, "of", length(fastqs)))
    
    # Read in json file with summaries of reads for each fastq
    json<-read_json(file.path(input_dir, "runs", paste0(fastq, ".extract.summary.json")))
    
    # Extract number of total reads in fastq and number of reads assigned to a barcode
    # Add to running total of processed and assigned
    processed<-processed+json[["n_processed"]]
    assigned<-assigned+json[["n_assigned"]]
    
  }
  
  for (j in 1:length(sample_names)) {
    # Process each sample
    sample<-sample_names[j]
    
    print(paste("Processing sample", sample, "number", j, "of", length(sample_names)))
    
    # Read in json with summary of demultiplexed reads per sample
    json<-read_json(file.path(input_dir, sample, "kallisto.run_info.json"))
    
    # Extract number of pseudoaligned (mapped) reads for each sample
    # Add to running total of pseudoaligned
    pseudoaligned<-pseudoaligned+json[["n_pseudoaligned"]]
    
    # Make a dataframe with number of processed (here meaning barcoded - this is the starting point for this json), pseudoaligned
    # and percentage of pseudoaligned (mapped/barcoded) reads per sample
    sample_add<-data.frame(sample = sample,
                           plate = plate,
                           assigned = json[["n_processed"]],
                           pseudoaligned = json[["n_pseudoaligned"]],
                           percent_bc_mapped = json[["p_pseudoaligned"]])
    
    # Add dataframe row of sample to overall sample_stats dataframe of all samples
    sample_stats<-rbind(sample_stats, sample_add)
    
}
  
  # Make a dataframe with total number of reads, number of assigned (barcoded) reads and number of pseudoaligned (mapped) reads for flowcell
  # Calculate percentage barcoded, percentage mapped total (mapped/total reads) and percentage mapped from barcoded reads (mapped/barcoded)
mapping_stats<-data.frame(plate = plate, 
                          n_processed = processed, 
                          n_assigned = assigned, 
                          n_pseudoaligned = pseudoaligned,
                          percent_barcoded = round(assigned/processed * 100, 2),
                          percent_all_mapped = round(pseudoaligned/processed * 100, 2),
                          percent_bc_mapped = round(pseudoaligned/assigned * 100, 2))
  



# Remove empty rows - used to create dataframes
sample_stats<-sample_stats[-1,]

# Save output tables
write.csv(mapping_stats, file.path(out_dirs, "mapping_stats", "plate_mapping_stats.csv"), row.names = F)
write.csv(sample_stats, file.path(out_dirs, "mapping_stats", "sample_mapping_stats.csv"), row.names = F)
 

### ANALYSE DATA ###

# Create bar plot showing percentage of barcoded reads mapped per sample
# Compare across flowcells
# Save figure
sample_mapping<-ggplot(sample_stats, aes(x=sample, y=percent_bc_mapped)) +
  geom_col(position = "dodge") +
  labs(title = "% of barcoded reads mapped", y = "%", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); sample_mapping
ggsave(file.path(out_dirs, "mapping_stats", "sample_mapping_stats.png"), sample_mapping, bg = "white", width = 10, height=6)

# Mutate mapping_stats dataframe to prepare for stacked barplot
# For each flowcell, calculate number of reads mapped (a), number of reads barcoded but not mapped (b) and number of reads neither barcoded or mapped (c)
df_long <- mapping_stats %>%
  mutate(
    a = n_pseudoaligned,
    b = n_assigned - n_pseudoaligned,
    c = n_processed - n_assigned
  ) %>%
  select(plate, a, b, c) %>%
  pivot_longer(cols = c(a, b, c), names_to = "status", values_to = "count")

# Dynamic label creation for text above bars to show percentages of mapped (for all and for barcoded) and barcoded reads
df_labels <- df_long %>%
  group_by(plate) %>%
  summarise(total=sum(count), .groups = "drop") %>%
  left_join(mapping_stats, by = "plate") %>%
  mutate(label1 = paste0(percent_all_mapped, "% of all reads mapped"),
         label2 = paste0(percent_bc_mapped, "% of barcoded reads mapped"),
         label3 = paste0(percent_barcoded, "% of all reads barcoded"))

# Create stacked barplot to show raw numbers of reads mapped, barcoded but not mapped, and not mapped or barcoded with dynamic % labels
# Save figure
flowcell_mapping<-ggplot(df_long, aes(x=factor(plate), y=count, fill=status)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("darkgreen","orange","darkred"), labels = c("Mapped", "Barcoded, not mapped", "Not mapped or barcoded")) +
  geom_text(data = df_labels, aes(x = plate, y = total, label = label1), vjust = -2.6, inherit.aes = FALSE) +
  geom_text(data = df_labels, aes(x = plate, y = total, label = label2), vjust = -1.55, inherit.aes = FALSE) +
  geom_text(data = df_labels, aes(x = plate, y = total, label = label3), vjust = -0.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))+
  theme_minimal() +
  labs(title = "Number of reads barcoded and mapped", y = "Number of reads", x = "Plate"); flowcell_mapping
ggsave(file.path(out_dirs, "mapping_stats", "plate_mapping_stats.png"), flowcell_mapping, bg = "white", width = 10, height=6)



