rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(data.table)
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

out_dirs<-file.path(base, "outputs")

# sample names - import from sample sheet
sample_names<-read.table(file.path(input_dir, "samples.tsv"))

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]

### EXTRACT DATA ###
# Set up empty dataframe to populate
qc_metrics<-data.frame(Sample=NA, 
                       Plate = plate, 
                       Total_barcodes = NA, 
                       Mean_UMIs_per_barcode = NA, 
                       Median_UMIs_per_barcode = NA, 
                       Total_raw_BUS_records = NA, 
                       Total_unique_UMI_barcode_pairs = NA, 
                       Knee_rank = NA)


  
for (i in 1:length(sample_names)) {
    # Process each sample
    sample<-sample_names[i]
    print(paste("Processing sample", sample, "number", i, "of", length(sample_names)))
    
    # Create sample specific QC directory to store metrics in
    dir.create(file.path(out_dirs, "QC", "samples", sample))
    
    # Import the file and set column names -SWAP UMI/BARCODE FOR FR RUN
    bus_text=file.path(input_dir, sample, "sorted_bus.txt")
    bus<-fread(bus_text, col.names = c("barcode", "UMI", "ec","count"))
    
    # Create summary table of number of number of UMIs and total numbers of each barcode 
    # Rank from most to least
    barcode_summary <- bus %>%
      group_by(barcode) %>%
      summarise(n_umis = n(),
                total_count = sum(count)) %>%
      arrange(desc(n_umis)) %>%
      mutate(rank = row_number())
    
    # For each sample, save this summary
    write.table(barcode_summary, file.path(out_dirs, "QC", "samples", sample, "barcode_summary.tsv"), sep = "\t", row.names = F, quote = F)
    
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
    
    write.csv(keep, file.path(out_dirs, "QC", "samples", sample, "barcode_keep_filter.csv"), row.names = F)
    
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
    ggsave(file.path(out_dirs, "QC", "samples", sample, "knee_plot.png"), knee, width = 10, height = 6)

    # Create histogram of UMI counts per barcode to check rough normal distribution
    # Save plot
    UMI<-ggplot(barcode_summary, aes(x=n_umis)) +
      geom_histogram(bins=50, fill = "steelblue") +
      scale_x_log10() +
      labs(title = "Distribution of UMI counts per barcode", x = "UMIs", y = "Frequency")
    ggsave(file.path(out_dirs, "QC", "samples", sample, "UMI_per_barcode.png"), UMI, width = 10, height = 6)
    
    # Create a dataframe with QC summary stats for sample, including total barcodes (number of cells), mean UMIs (unique molecules per cell),
    # total records (reads per sample), unique UMIs (unique reads - without PCR duplicates) and knee rank (number of barcodes before inflection)
    qc_add<-data.frame(Sample = sample,
                       Plate = plate,
                       Total_barcodes = nrow(barcode_summary), 
                       Mean_UMIs_per_barcode = mean(barcode_summary$n_umis), 
                       Median_UMIs_per_barcode =  median(barcode_summary$n_umis),
                       Total_raw_BUS_records = sum(bus$count),
                       Total_unique_UMI_barcode_pairs = nrow(bus),
                       Knee_rank = knee_rank)
    
    # Add dataframe row for sample to overall qc_metrics dataframe for all samples
    qc_metrics<-rbind(qc_metrics, qc_add)
    
}

# Remove empty row used to create dataframe
qc_metrics<-qc_metrics[-1,]

# Calculate reads per cell and library complexity
qc_metrics<-qc_metrics %>%
  mutate(reads_per_barcode = Total_raw_BUS_records / Total_barcodes) %>%
  mutate(library_complexity = Total_unique_UMI_barcode_pairs / Total_raw_BUS_records)

# Save qc_metrics table with info for all samples
write.csv(qc_metrics, file.path(out_dirs, "QC", "qc_metrics.csv"), row.names = F)




### ANALYSE DATA ###

# Bar plot of the number of complexity per sample, compared across flowcells
# Save the plot
complexity<-ggplot(qc_metrics, aes(x=Sample, y=library_complexity)) +
  geom_col(position = "dodge")+
  labs(title = "Library complexity", y = "Proportion of unique molecules", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); complexity
ggsave(file.path(out_dirs, "QC", "library_complexity.png"), complexity, bg="white", width=10, height=6)

# Bar plot of mean UMIs (unique molecules) per barcode per sample compared across flowcells 
# Save the plot
UMIs<-ggplot(qc_metrics, aes(x=Sample, y=Mean_UMIs_per_barcode)) +
  geom_col(position = "dodge")+
  labs(title = "Mean UMIs per barcode per flowcell", y = "Mean UMIs", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); UMIs
ggsave(file.path(out_dirs, "QC", "mean_UMIs.png"), UMIs, bg="white", width=10, height=6)

# Bar plot of knee ranks (number of barcodes above inflection point) per sample, compared across flowcells
# Save the plot
knee<-ggplot(qc_metrics, aes(x=Sample, y=Knee_rank)) +
  geom_col(position = "dodge")+
  labs(title = "Knee ranks per flowcell", y = "Knee rank", x = "Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); knee
ggsave(file.path(out_dirs, "QC", "knee_ranks.png"), knee, bg="white", width=10, height=6)

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
ggsave(file.path(out_dirs, "QC", "knee_inclusion.png"), inclusion, bg="white", width=10, height=6)
