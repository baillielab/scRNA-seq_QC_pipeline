rm(list=ls())

print("Starting 02.2_initial_QC.R")

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(yaml)

# Ingest config
config <- yaml::read_yaml(args[1])

base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

in_dir <- file.path(base, config$input$data_dir)
if (is.null(in_dir)) stop("config$input$data_dir is required")

plate<-config$input$plate

out_dir<-file.path(base, "outputs", "QC")

sample_out<-file.path(out_dir, "samples")
dir.create(sample_out)

samples<-fread(file.path(in_dir, "sample-list.txt"), header = F)
samples<-samples$V1

# Create function to detect inflection point (knee)
get_knee<-function(x,y) {
  # Define first and last points on curve to create line between
  p1<-c(x[1], y[1])
  p2<-c(x[length(x)], y[length(y)])
  # Vector between points 1 and 2, normalised to project other points on to
  line_vec<-p2-p1
  line_vec<-line_vec / sqrt(sum(line_vec^2))
  
  distances <- sapply(1:length(x), function(i) {
    # For every point on the curve
    p <- c(x[i], y[i])
    # Get the vector from p1 to the current point
    vec <- p-p1
    # Project the vector onto the line = scalar projection/dot product and get the point on the line where the vector lands
    proj_len<-sum(vec*line_vec)
    proj_point<-p1 + proj_len * line_vec
    # Calculate Euclidean (perpendicular) distance distance between p and projection
    dist<-sqrt(sum((p - proj_point)^2))
    return(dist)
  })
  # Find the index of the point with the greatest distance from the line = knee point 
  return(which.max(distances))
}



qc_metrics<-data.frame(Sample=NA, 
                       Plate = NA,
                       Total_barcodes = NA,
                       Mean_UMIs_per_barcode = NA, 
                       Median_UMIs_per_barcode= NA, 
                       Total_raw_BUS_records = NA, 
                       Total_unique_UMI_barcode_pairs = NA,
                       Knee_rank = NA)

for (i in 1:length(samples)) {
  sample<-samples[i]
  
  # Identify sample and create QC output directory
  print(paste0("Processing sample ", sample, ". Number ", i, " of ", length(samples)))
  dir.create(file.path(sample_out, sample))
  
  barcodes<-read.csv(file.path(in_dir, sample, "DGE_unfiltered", "cell_metadata.csv"))
  matrix<-readMM(file.path(in_dir, sample, "DGE_unfiltered", "count_matrix.mtx"))
  
  # Create summary table of number of number of UMIs and total numbers of each barcode 
  # Rank from most to least
  barcode_summary<-tibble(
    barcode = barcodes$bc_wells,
    n_umis = rowSums(matrix),
    n_genes = rowSums(matrix > 0)
  ) %>%
    arrange(desc(n_umis)) %>%
    mutate(rank = row_number())
  
  # For each sample, save this summary
  write.csv(barcode_summary, file.path(sample_out, sample, "barcode_summary.csv"), row.names = F)
  
  # Prep log10 of x and y variables (rank and number of UMIs) to put in to function to detect inflection point (knee)
  x<-log10(barcode_summary$rank)
  y<-log10(barcode_summary$n_umis + 1)
  
  # Apply function to x and y generated in previous step
  knee_index<-get_knee(x,y)
  
  # Find corresponding rank and n_umis for knee index
  knee_rank<-barcode_summary$rank[knee_index]
  knee_umis<-barcode_summary$n_umis[knee_index]
  
  
  knee_rank2<-barcode_summary$rank[knee_index*5]
  knee_umis2<-barcode_summary$n_umis[knee_index*5]
  
  knee_filter<-knee_rank*5
  
  keep<-data.frame(barcode =
                     barcode_summary$barcode[which(barcode_summary$rank < knee_filter)])
  
  write.csv(keep, file.path(out_dir, "samples", sample, "barcode_keep_filter.csv"), row.names = F)
  
  
  # Generate knee plot with knee rank annotated
  # Save plot
  knee<-ggplot(barcode_summary, aes(x=rank, y=n_umis)) +
    geom_line() +
    scale_x_log10() + scale_y_log10() +
    geom_point(data = barcode_summary[knee_index, ], aes(x=rank, y=n_umis), color = "deepskyblue", size = 2) +
    geom_point(data = barcode_summary[knee_index*5, ], aes(x=rank, y=n_umis), color = "blue", size = 2) +
    labs(title = "Barcode Rank Plot (log-log)", x = "Rank", y = "UMI count") +
    annotate("text", x = knee_rank * 1.2, y = knee_umis, label = paste("Knee @ rank", knee_rank), color = "deepskyblue") +
    annotate("text", x = knee_rank2 / 2, y = knee_umis2, label = paste("Knee filter @ rank", knee_rank2), color = "blue")
  ggsave(file.path(sample_out, sample, "knee_plot.png"), knee, width = 10, height = 6)
  
  # Create histogram of UMI counts per barcode to check rough normal distribution
  # Save plot
  UMI<-ggplot(barcode_summary, aes(x=n_umis)) +
    geom_histogram(bins=50, fill = "steelblue") +
    scale_x_log10() +
    labs(title = "Distribution of UMI counts per barcode", x = "UMIs", y = "Frequency")
  ggsave(file.path(sample_out, sample, "UMI_per_barcode.png"), UMI, width = 10, height = 6)
  
  # Create a dataframe with QC summary stats for sample, including total barcodes (number of cells), mean UMIs (unique molecules per cell),
  # total records (reads per sample), unique UMIs (unique reads - without PCR duplicates) and knee rank (number of barcodes before inflection)
  sample_QC<-data.frame(Sample=sample, 
                        Plate = plate,
                        Total_barcodes = nrow(barcodes),
                        Mean_UMIs_per_barcode = mean(rowSums(matrix)), 
                        Median_UMIs_per_barcode= median(rowSums(matrix)), 
                        Total_raw_BUS_records = sum(matrix), 
                        Total_unique_UMI_barcode_pairs = length(matrix@x),
                        Knee_rank = knee_rank
  )
  
  # Add dataframe row for sample to overall qc_metrics dataframe for all samples
  qc_metrics<-rbind(qc_metrics, sample_QC)
}

# Remove empty row used to create dataframe
qc_metrics<-qc_metrics[-1,]

# Calculate reads per cell and library complexity
qc_metrics<-qc_metrics %>%
  mutate(reads_per_barcode = Total_raw_BUS_records / Total_barcodes) %>%
  mutate(library_complexity = Total_unique_UMI_barcode_pairs / Total_raw_BUS_records)

# Save qc_metrics table with info for all samples
write.csv(qc_metrics, file.path(out_dir, "qc_metrics.csv"), row.names = F)

### ANALYSE DATA ###

# Bar plot of the number of complexity per sample, compared across flowcells
# Save the plot
complexity<-ggplot(qc_metrics, aes(x=Sample, y=library_complexity)) +
  geom_col(position = "dodge")+
  labs(title = "Library complexity", y = "Proportion of unique molecules", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); complexity
ggsave(file.path(out_dir, "library_complexity.png"), complexity, bg = "white", width = 10, height = 6)
  
# Bar plot of mean UMIs (unique molecules) per barcode per sample compared across flowcells 
# Save the plot
UMIs<-ggplot(qc_metrics, aes(x=Sample, y=Mean_UMIs_per_barcode)) +
  geom_col(position = "dodge")+
  labs(title = "Mean UMIs per barcode per flowcell", y = "Mean UMIs", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); UMIs
ggsave(file.path(out_dir, "mean_UMIs.png"), UMIs, bg = "white", width = 10, height = 6)
  
# Bar plot of knee ranks (number of barcodes above inflection point) per sample, compared across flowcells
# Save the plot
knee<-ggplot(qc_metrics, aes(x=Sample, y=Knee_rank)) +
  geom_col(position = "dodge")+
  labs(title = "Knee ranks per flowcell", y = "Knee rank", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); knee
ggsave(file.path(out_dir, "knee_ranks.png"), knee, bg = "white", width = 10, height = 6)
  
# Mutate qc_metrics to allow for stacked barplot
# For each sample, identify number of barcodes above knee rank and number below
df_long <- qc_metrics %>%
  mutate(
    included = Knee_rank,
    excluded = Total_barcodes - Knee_rank
  ) %>%
  select(Sample, Plate, included, excluded) %>%
  pivot_longer(cols = c(included, excluded), names_to = "status", values_to = "count")

# Stacked bar plot of number of barcodes above and below knee rank for each sample, split by flowcell
# Save the plot
inclusion<-ggplot(df_long, aes(x=Sample, y=count, fill=status)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(included = "steelblue", excluded = "grey70")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Number of barcodes included and excluded using knee rank as filtering cut off", y = "Number of barcodes"); inclusion
ggsave(file.path(out_dir, "knee_inclusion.png"), inclusion,  bg="white", width = 10, height = 6)

