rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

p1_fc1<-read.csv("/odp-beegfs/a008/02_QC/plate1/FC1/results/filtering_criteria.csv")
p1_fc2<-read.csv("/odp-beegfs/a008/02_QC/plate1/FC2/results/filtering_criteria.csv")
p2<-read.csv("/odp-beegfs/a008/02_QC/plate2/results/filtering_criteria.csv")
p3<-read.csv("/odp-beegfs/a008/02_QC/plate3/results/filtering_criteria.csv")
p4<-read.csv("/odp-beegfs/a008/02_QC/plate4/results/filtering_criteria.csv")
p5<-read.csv("/odp-beegfs/a008/02_QC/plate5/results/filtering_criteria.csv")
p6<-read.csv("/odp-beegfs/a008/02_QC/plate6/results/filtering_criteria.csv")


all_samples<-rbind(p1_fc1, p1_fc2, p2, p3, p4, p5,p6)
rm(p1_fc1, p1_fc2, p2, p3, p4, p5,p6)

staged<-data.frame(sample = all_samples$sample, start = all_samples$total_cells, 
                   gene_number_filtering = all_samples$cells_post_nFeature, mitochondrial_filtering = all_samples$cells_post_mt,
                   doublet_filtering = all_samples$cells_post_doublet, death_filtering = all_samples$cells_post_death, tissue = NA)

staged$tissue[grep("BLD", staged$sample)]<-"Blood"
staged$tissue[grep("BRN", staged$sample)]<-"Brain"
staged$tissue[grep("SKN", staged$sample)]<-"Skin"
staged$tissue[grep("LGS", staged$sample)]<-"Lung"
staged$tissue[grep("EVL", staged$sample)]<-"EVLP"
staged$tissue[grep("TSS", staged$sample)]<-"TSW_Skin"
staged$tissue[grep("TSB", staged$sample)]<-"TSW_Blood"
staged$tissue[grep("UCB", staged$sample)]<-"Umbilical_blood"


stage_order<-c(
  "start",
  "gene_number_filtering",
  "mitochondrial_filtering",
  "doublet_filtering",
  "death_filtering"
)

qc_long <- staged %>%
  mutate(sample = sample) %>%
  pivot_longer(
    cols = all_of(stage_order),
    names_to = "stage",
    values_to = "n_cells"
  ) %>%
  mutate(
    stage = factor(stage, levels = stage_order)
  )

qc_summary <- qc_long %>%
  group_by(tissue, stage) %>%
  summarise(
    total_cells = sum(n_cells, na.rm = T),
    median_cells = median(n_cells, na.rm = T),
    q25 = quantile(n_cells, 0.25, na.rm = T),
    q75 = quantile(n_cells, 0.75, na.rm = T),
    .groups = "drop"
  ) %>%
  group_by(tissue) %>%
  mutate(
    pct_start = 100 * (total_cells / first(total_cells))
  ) %>%
  ungroup()

global_ymin <- min(qc_long$n_cells[qc_long$n_cells > 0], na.rm = T)
global_ymax <- max(
  qc_long$n_cells,
  na.rm = T
)

make_qc_panel <- function(tissue_name, tissue_colour) {
  df <- qc_long %>%
    filter(tissue == tissue_name)
  
  sm <- qc_summary %>%
    filter(tissue == tissue_name)
  
  log_low <- log10(global_ymin)
  log_high <- log10(global_ymax)
  
  total_low <- log_high + 0.25
  total_high <- log_high + 1.00
  
  sm <- sm %>%
    mutate(
      total_y = 10^rescale(
        pct_start,
        from = c(0, 100),
        to = c(
          log_high + 0.15,
          log_high + 1
        )
      )
    )
  
  p_main <- ggplot(df, aes(x=stage, y = n_cells)) +
    
    geom_violin(
      fill = tissue_colour,
      colour = "black",
      width = 0.72,
      trim = T,
      alpha = 0.75
    ) +
    
    geom_errorbar(
      data = sm,
      aes(
        x = stage,
        ymin = q25,
        ymax = q75
      ),
      inherit.aes = F,
      width = 0.08,
      linewidth = 2.5,
      colour = "white"
    ) +
    
    geom_point(
      data = sm,
      aes(
        x = stage,
        y = median_cells
      ),
      inherit.aes = F,
      size = 2.5,
      colour = "black"
    ) +
  
  geom_line(
    data = sm,
    aes(x = stage,
        y = total_y,
        group = 1),
    inherit.aes = F,
    linewidth = 0.9,
    colour = "black"
  ) +
  
  geom_point(
    data = sm,
    aes(
      x = stage,
      y = total_y
    ),
    inherit.aes = F,
    size = 2.5,
    colour = "black"
  ) +
  
  geom_text(
    data = sm,
    aes(
      x = stage,
      y = total_y,
      label = comma(total_cells)
    ),
    inherit.aes = F,
    vjust = -0.7,
    size = 3
    ) +
  
  scale_y_log10(
    labels = comma,
    limits = c(
      global_ymin,
      10^(log_high + 2.15)
    ),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    title = tissue_name,
    x = NULL,
    y = "Cells per sample"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =  element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    plot.margin = margin(5, 5, 0, 5)
  )
  
  table_df <- sm %>%
  select(stage, median_cells, pct_start) %>%
  pivot_longer(
    cols = c(median_cells, pct_start),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    label = case_when(
      metric == "median_cells" ~ comma(round(value)),
      metric == "pct_start" ~ paste0(round(value, 1), "%")
    ),
    metric = factor(
      metric,
      levels = c("median_cells", "pct_start"),
      labels = c(
        "Median cells\nper sample",
        "% of start\n(total cells)"
      )
    )
  )
  
  p_table <- ggplot(
  table_df,
  aes(x = stage, y = metric)
  ) + 
  
  geom_tile(
    fill = "white",
    colour = "grey60",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(label = label),
    size = 3.4
  ) +
  
  scale_x_discrete(drop = F) +

  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_void() +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1,
    colour = "black",
    size = 9
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 9
  ),
  plot.margin = margin(0,5,5,5)
  )

p_main / p_table +
  plot_layout(
    heights = c(4.5, 1.25)
  )
}

p_blood <- make_qc_panel("Blood", "#7FCDBB")
p_brain <- make_qc_panel("Brain", "#FC8D62")
p_lung <- make_qc_panel("Lung", "#8DA0CB")
p_skin <- make_qc_panel("Skin", "#E78AC3")
p_evlp <- make_qc_panel("EVLP", "#C6BBFA")
p_tss <- make_qc_panel("TSW_Skin", "#A8224C")
p_tsb <- make_qc_panel("TSW_Blood", "#749e60")

final_plot<-
  (p_blood  | p_brain  | p_skin  | p_lung) /
  (p_evlp   | p_tss    | p_tsb   | plot_spacer()) +
  plot_annotation(
    caption = "Violins show the distribution of cells per sample at each QC stage. Black line and labels show total cells across all samples within each tissue"
  ); final_plot

ggsave(
  "/odp-beegfs/a008/02_QC/all_QC_summary.png",
  final_plot,
  width = 20,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)


