rm(list=ls())

print("Starting 02.1_mapping_stats.R")

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(yaml)

# Ingest config
config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

input_dir <- file.path(base, config$input$data_dir)
if (is.null(input_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate



out_dirs<-file.path(base, "outputs")
dir.create(out_dirs, recursive = T)

samples<-fread(file.path(input_dir, "sample-list.txt"), header = F)
samples<-samples$V1

# Set up empty dataframes to populate
sample_stats<-data.frame(sample = NA,
                         plate = NA,
                         processed = NA,
                         assigned = NA,
                         pseudoaligned = NA,
                         percent_bc_mapped = NA)

processed_sum<-0
assigned_sum<-0
pseudoaligned_sum<-0

for (i in 1:length(samples)) {
  sample<-samples[i]
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(samples)))
  
  stats<-read.csv(file.path(input_dir, sample, "report", "analysis_summary.csv"))
  
  p<-which(stats$statistic == "number_of_reads")
  b<-which(stats$statistic == "valid_barcode_fraction")
  a<-which(stats$statistic == "transcriptome_map_fraction")
  
  
  processed<-stats$combined[p]
  assigned<-processed * stats$combined[b]
  pseudoaligned<-assigned * stats$combined[a]
  
  processed_sum<-processed_sum+processed
  assigned_sum<-assigned_sum+assigned
  pseudoaligned_sum<-pseudoaligned_sum+pseudoaligned
  
  sample_add<-data.frame(
    sample = sample,
    plate = plate,
    processed = processed,
    assigned = assigned,
    pseudoaligned = pseudoaligned,
    percent_bc_mapped = 100 * (pseudoaligned/assigned)
  )
  
  sample_stats<-rbind(sample_stats, sample_add)
  
}

sample_stats<-sample_stats[-1,]

mapping_stats<-data.frame(plate = plate, 
                          n_processed = processed_sum, 
                          n_assigned = assigned_sum, 
                          n_pseudoaligned = pseudoaligned_sum,
                          percent_barcoded = round(100 * (assigned_sum/processed_sum),2),
                          percent_all_mapped = round(100 * (pseudoaligned_sum/processed_sum),2),
                          percent_bc_mapped = round(100 * (pseudoaligned_sum/assigned_sum),2)
)

# Save output tables
summary_out<-file.path(out_dirs, "QC", "mapping_stats")
dir.create(summary_out, recursive = T)

write.csv(mapping_stats, file.path(summary_out, "plate_mapping_stats.csv"), row.names = F)
write.csv(sample_stats, file.path(summary_out, "sample_mapping_stats.csv"), row.names = F)

### ANALYSE DATA ###

# Create bar plot showing percentage of barcoded reads mapped per sample
# Compare across flowcells
# Save figure
sample_mapping<-ggplot(sample_stats, aes(x=sample, y=percent_bc_mapped)) +
  geom_col(position = "dodge") +
  labs(title = "% of barcoded reads mapped", y = "%", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); sample_mapping
ggsave(file.path(summary_out, "sample_mapping_stats.png"), sample_mapping, bg = "white", width = 10, height=6)

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
  geom_text(data = df_labels, aes(x = 1, y = total, label = label1), vjust = -2.6, inherit.aes = FALSE) +
  geom_text(data = df_labels, aes(x = 1, y = total, label = label2), vjust = -1.55, inherit.aes = FALSE) +
  geom_text(data = df_labels, aes(x = 1, y = total, label = label3), vjust = -0.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))+
  theme_minimal() +
  labs(title = "Number of reads barcoded and mapped", y = "Number of reads", x = "Plate"); flowcell_mapping
ggsave(file.path(summary_out, "plate_mapping_stats.png"), flowcell_mapping, bg = "white", width = 10, height=6)

