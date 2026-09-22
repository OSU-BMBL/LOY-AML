source("analysis/human_bulk/config.R")

# R Script for Overlap Analysis: LOY AML DEGs vs. ZFY Knockdown DEGs
#
# Description:
# This script compares differentially expressed genes (DEGs) from a Loss-of-Y (LOY)
# vs. Control-Normal (CN) AML bulk RNA-seq experiment with a published list of genes
# affected by ZFY knockdown (San Roman et al., Cell Genomics, 2024).
#
# The primary analysis compares down-regulated genes in both datasets. An
# exploratory analysis for up-regulated genes is also included.
#
# Steps:
# 1. Load required libraries.
# 2. Define file paths and analysis parameters.
# 3. Load and process user's LOY vs. CN DEG data (up and down).
# 4. Load and process the ZFY knockdown data (up and down).
# 5. Perform overlap analysis for both down- and up-regulated gene sets.
# 6. Generate and save visualizations: two Venn diagrams and one volcano plot.
# 7. Save the lists of overlapping genes.
#

# 1. SETUP: Load Libraries

library(readxl)
library(dplyr)
library(ggplot2)
library(ggvenn)
library(ggrepel)
library(ggpubr)
# 2. SETUP: Define File Paths and Parameters

deg_file_path <- bulk_files$zfy_deg
zfy_paper_file_path <- bulk_files$zfy_reference
output_directory <- bulk_output("figures", "zfy_overlap")
table_directory <- bulk_output("tables", "zfy_overlap")

# Create the output directory if it doesn't exist
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# 3. LOAD & PROCESS: LOY vs. CN DEG Data

cat("--- Loading and processing the DEG table... ---\n")

# Read the DEG file
loy_degs_all <- read.csv(deg_file_path)
colnames(loy_degs_all)[1] <- "Gene"

# Filter for significantly DOWN-REGULATED genes in LOY samples
loy_degs_down <- loy_degs_all %>%
  filter(padj < 0.05, log2FoldChange < 0)
loy_down_genes <- loy_degs_down$Gene

# Filter for significantly UP-REGULATED genes in LOY samples
loy_degs_up <- loy_degs_all %>%
  filter(padj < 0.05, log2FoldChange > 0)
loy_up_genes <- loy_degs_up$Gene

cat(paste("Found", nrow(loy_degs_all), "total genes tested in the dataset.\n"))
cat(paste("Found", length(loy_down_genes), "genes significantly down-regulated in LOY AML.\n"))
cat(paste("Found", length(loy_up_genes), "genes significantly up-regulated in LOY AML.\n\n"))

# 4. LOAD & PROCESS: Published ZFY Knockdown Data

cat("--- Loading and processing ZFY knockdown data from San Roman et al. (2024)... ---\n")

# The most relevant data is the ZFY knockdown in XY cells.
# This corresponds to "Table S11E" in the file provided with the paper.
sheet_names <- excel_sheets(zfy_paper_file_path)
zfy_sheet_name <- "S11A"

if(!zfy_sheet_name %in% sheet_names){
  stop(paste("The specified sheet '", zfy_sheet_name, "' was not found. Please check the Excel file. Available sheets are:", paste(sheet_names, collapse=", ")))
}

# Read the specific sheet for ZFY knockdown
zfy_kd_all <- read_excel(zfy_paper_file_path, sheet = zfy_sheet_name, skip = 3)

# Filter for genes significantly DOWN-REGULATED upon ZFY knockdown
zfy_kd_down <- zfy_kd_all %>%
  filter(padj < 0.05, log2FoldChange < 0)
zfy_down_genes <- zfy_kd_down$Gene

# Filter for genes significantly UP-REGULATED upon ZFY knockdown
zfy_kd_up <- zfy_kd_all %>%
  filter(padj < 0.05, log2FoldChange > 0)
zfy_up_genes <- zfy_kd_up$Gene

cat(paste("Found", length(zfy_down_genes), "genes significantly down-regulated upon ZFY knockdown.\n"))
cat(paste("Found", length(zfy_up_genes), "genes significantly up-regulated upon ZFY knockdown.\n\n"))

# 5. OVERLAP ANALYSIS

cat("--- Performing overlap analysis... ---\n")

# a) Overlap for DOWN-REGULATED genes (Primary hypothesis)
down_gene_lists <- list(
  `Down in LOY AML` = loy_down_genes,
  `Down in ZFY Knockdown` = zfy_down_genes
)
overlapping_down_genes <- intersect(loy_down_genes, zfy_down_genes)
n_overlap_down <- length(overlapping_down_genes)
cat(paste("Found", n_overlap_down, "overlapping DOWN-regulated genes.\n"))

# b) Overlap for UP-REGULATED genes (Exploratory)
up_gene_lists <- list(
  `Up in LOY AML` = loy_up_genes,
  `Up in ZFY Knockdown` = zfy_up_genes
)
overlapping_up_genes <- intersect(loy_up_genes, zfy_up_genes)
n_overlap_up <- length(overlapping_up_genes)
cat(paste("Found", n_overlap_up, "overlapping UP-regulated genes.\n\n"))

# 6. VISUALIZATION
cat("--- Generating and saving plots... ---\n")

# a) Venn Diagram for DOWN-regulated genes
venn_plot_down <- ggvenn(
  down_gene_lists,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4, text_size = 3
) + labs(title = "Overlap of Down-regulated Genes")

ggsave(file.path(output_directory, "Venn_Diagram_DOWN_LOY_vs_ZFYkd.png"), plot = venn_plot_down, width = 6, height = 5, dpi = 300)
cat("Saved Venn Diagram for down-regulated genes.\n")

# b) Venn Diagram for UP-regulated genes
venn_plot_up <- ggvenn(
  up_gene_lists,
  fill_color = c("#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4, text_size = 3
) + labs(title = "Overlap of Up-regulated Genes")

ggsave(file.path(output_directory, "Venn_Diagram_UP_LOY_vs_ZFYkd.png"), plot = venn_plot_up, width = 6, height = 5, dpi = 300)
cat("Saved Venn Diagram for up-regulated genes.\n")

# c) Volcano Plot highlighting both overlaps
loy_degs_all <- loy_degs_all %>%
  mutate(Highlight = case_when(
    Gene %in% overlapping_down_genes ~ "Overlap (Down)",
    Gene %in% overlapping_up_genes   ~ "Overlap (Up)",
    padj < 0.05 & log2FoldChange < 0 ~ "Down in LOY AML",
    padj < 0.05 & log2FoldChange > 0 ~ "Up in LOY AML",
    TRUE ~ "Not Significant"
  ))

loy_degs_all$Highlight <- factor(loy_degs_all$Highlight, levels = c("Overlap (Down)", "Overlap (Up)", "Down in LOY AML", "Up in LOY AML", "Not Significant"))

volcano_plot <- ggplot(loy_degs_all, aes(x = log2FoldChange, y = -log10(pvalue), color = Highlight)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Overlap (Down)" = "red",
                                "Overlap (Up)" = "purple",
                                "Down in LOY AML" = "#0073C2FF",
                                "Up in LOY AML" = "orange",
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot: LOY-Male vs. CN-Male AML",
    subtitle = "Genes overlapping with ZFY knockdown data are highlighted",
    x = "Log2 Fold Change",
    y = "-Log10 P-value"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_text_repel(
    data = subset(loy_degs_all, Highlight %in% c("Overlap (Down)", "Overlap (Up)")) %>% arrange(padj) %>% head(20),
    aes(label = Gene),
    size = 3, box.padding = 0.5, point.padding = 0.5,
    segment.color = 'grey50', max.overlaps = 20
  )

ggsave(file.path(output_directory, "Volcano_Plot_highlighted.png"), plot = volcano_plot, width = 10, height = 8, dpi = 300)
cat("Saved Volcano Plot.\n\n")

# 7. SAVE OVERLAPPING GENE LISTS
cat("--- Saving lists of overlapping genes... ---\n")

# Save DOWN-regulated overlapping genes
overlapping_down_data <- loy_degs_all %>%
  filter(Gene %in% overlapping_down_genes) %>%
  arrange(padj)
write.csv(overlapping_down_data, file.path(table_directory, "overlapping_DOWNregulated_genes.csv"), row.names = FALSE)
cat(paste("Saved overlapping DOWN-regulated gene data.\n"))

# Save UP-regulated overlapping genes
overlapping_up_data <- loy_degs_all %>%
  filter(Gene %in% overlapping_up_genes) %>%
  arrange(padj)
write.csv(overlapping_up_data, file.path(table_directory, "overlapping_UPregulated_genes.csv"), row.names = FALSE)
cat(paste("Saved overlapping UP-regulated gene data.\n"))

cat("\n--- Analysis Complete! ---\n")

# 8. SCATTER PLOT VISUALIZATION: LOY vs. ZFY Fold Changes

cat("\n--- Generating scatter plot comparing Log2 Fold Changes... ---\n")

# a) Merge the two datasets to get fold changes for common genes
# We use inner_join to keep only genes present in both datasets.
# Suffixes are added to distinguish columns with the same name (e.g., log2FoldChange).
merged_degs <- inner_join(loy_degs_all, zfy_kd_all, by = "Gene", suffix = c("_LOY", "_ZFY"))

# b) Add the same highlight column for coloring the points
merged_degs <- merged_degs %>%
  mutate(Highlight = case_when(
    Gene %in% overlapping_down_genes ~ "Overlap (Down)",
    Gene %in% overlapping_up_genes   ~ "Overlap (Up)",
    TRUE                             ~ "Other Genes"
  )) %>%
  # Make it a factor to control the plotting order (highlights on top)
  mutate(Highlight = factor(Highlight, levels = c("Overlap (Up)", "Overlap (Down)", "Other Genes"))) %>%
  dplyr::filter(padj_LOY < 0.05 & padj_ZFY < 0.05 )

merged_degs2 <- merged_degs %>%
  dplyr::filter(padj_LOY < 0.05 & padj_ZFY < 0.05 & log2FoldChange_LOY > 0 & log2FoldChange_ZFY > 0)

cor_test_result <- cor.test(merged_degs$log2FoldChange_LOY, merged_degs$log2FoldChange_ZFY)
correlation_label <- sprintf("Pearson's R = %.2f\np-value < 2.2e-16", cor_test_result$estimate)

# d) Create the scatter plot
scatter_plot <- ggplot(merged_degs, aes(x = log2FoldChange_LOY, y = log2FoldChange_ZFY)) +
  # Draw points for all genes, with non-highlighted ones being more transparent
  geom_point(aes(color = Highlight), alpha = 0.7, size = 2) +
  # Define colors and labels for the legend
  scale_color_manual(
    name = "Gene Set",
    values = c(
      "Overlap (Down)" = "red",
      "Overlap (Up)" = "purple",
      "Other Genes" = "grey70"
    )
  ) +
  # Add quadrant lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Add labels for the highlighted overlapping genes
  geom_text_repel(
    data = subset(merged_degs, Highlight != "Other Genes"),
    aes(label = Gene),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = 20
  ) +
  # Add the correlation text to the plot
  annotate(
    "text", x = Inf, y = -Inf,
    label = correlation_label,
    hjust = 1.05, vjust = -0.5,
    size = 4, color = "black"
  ) +
  # Add titles and labels
  labs(
    title = "Correlation of Gene Expression Changes",
    subtitle = "LOY AML vs. ZFY Knockdown",
    x = "Log2 Fold Change (Down/Up in LOY AML)",
    y = "Log2 Fold Change (Down/Up in ZFY Knockdown)"
  ) +
  theme_classic2() + # A clean theme for scatter plots
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# e) Save the scatter plot
ggsave(
  file.path(output_directory, "Scatter_Plot_LOY_vs_ZFYkd.png"),
  plot = scatter_plot,
  width = 8, height = 8, dpi = 300
)

cat("Saved scatter plot comparing fold changes.\n")

cat("\n--- Analysis Fully Complete! ---\n")

# Quadrant schematic (no points)
library(ggplot2)
library(grid)  # for unit()

# Where to save
outdir <- output_directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Axis limits for the schematic (tweak if you want longer/shorter arrows)
xlim <- c(-3, 3)
ylim <- c(-3, 3)

quad_schematic <- ggplot() +
  # horizontal & vertical arrow axes
  annotate("segment", x = min(xlim), xend = max(xlim), y = 0, yend = 0,
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +
  annotate("segment", x = 0, xend = 0, y = min(ylim), yend = max(ylim),
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +
  # labels
  annotate("text", x = 0, y = max(ylim), label = "Up in ZFY KO", vjust = -0.6, size = 5.5) +
  annotate("text", x = 0, y = min(ylim), label = "Down in ZFY KO", vjust = 1.6, size = 5.5) +
  annotate("text", x = max(xlim), y = 0, label = "Up in LOY AML", hjust = -0.06, size = 5.5) +
  annotate("text", x = min(xlim), y = 0, label = "Down in LOY AML", hjust = 1.06, size = 5.5) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(10, 20, 10, 20),
    plot.title = element_text(hjust = 0, face = "bold", size = 13)
  )

# Save PDF (embed fonts nicely with cairo_pdf if available)
ggsave(file.path(outdir, "Quadrant_Schematic_LOY_ZFY.pdf"),
       quad_schematic, width = 4.2, height = 3.2, device = cairo_pdf)

# Optional PNG for quick preview
ggsave(file.path(outdir, "Quadrant_Schematic_LOY_ZFY.png"),
       quad_schematic, width = 4.2, height = 3.2, dpi = 300)

# Scatter with quadrant arrows/labels (builds on merged_degs)
# Assumes 'merged_degs' from the preceding analysis exists with columns:
# log2FoldChange_LOY, log2FoldChange_ZFY, and 'Highlight'

# Set symmetric limits around 0 based on data (with a little padding)
pad <- 0.2
lim_x <- max(abs(range(merged_degs$log2FoldChange_LOY, na.rm = TRUE))) + pad
lim_y <- max(abs(range(merged_degs$log2FoldChange_ZFY, na.rm = TRUE))) + pad
scatter_quad <- ggplot(merged_degs,
                       aes(x = log2FoldChange_LOY, y = log2FoldChange_ZFY)) +
  geom_point(aes(color = Highlight), alpha = 0.7, size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    name = NULL,
    values = c("Overlap (Down)" = "red",
               "Overlap (Up)"   = "purple",
               "Other Genes"    = "grey70")
  ) +
  # remove dashed quadrant lines
  #geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black") +

  # heavy arrow axes
  annotate("segment", x = -lim_x, xend =  lim_x, y = 0, yend = 0,
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +
  annotate("segment", x = 0, xend = 0, y = -lim_y, yend =  lim_y,
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +

  # labels
  annotate("text", x = 0, y =  lim_y, label = "Up in ZFY KO", vjust = -0.6, size = 5.0) +
  annotate("text", x = 0, y = -lim_y, label = "Down in ZFY KO", vjust = 1.6, size = 5.0) +
  annotate("text", x =  lim_x, y = 0, label = "Up in LOY AML", hjust = -0.06, size = 5.0) +
  annotate("text", x = -lim_x, y = 0, label = "Down in LOY AML", hjust = 1.06, size = 5.0) +

  coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y), expand = FALSE) +
  labs(
    title = "LOY AML vs ZFY Knockout",
    x = "Log2 Fold Change (LOY AML)",
    y = "Log2 Fold Change (ZFY KO)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(10, 20, 10, 20),
    # strip axis lines/ticks
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file.path(outdir, "Scatter_Quadrant_LOY_vs_ZFYKO.pdf"),
       scatter_quad, width = 5.0, height = 5.0, device = cairo_pdf)
ggsave(file.path(outdir, "Scatter_Quadrant_LOY_vs_ZFYKO.png"),
       scatter_quad, width = 5.0, height = 5.0, dpi = 300)

library(ggrepel)

scatter_quad_annot <- ggplot(merged_degs,
                             aes(x = log2FoldChange_LOY, y = log2FoldChange_ZFY)) +
  geom_point(aes(color = Highlight), alpha = 0.7, size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    name = NULL,
    values = c("Overlap (Down)" = "red",
               "Overlap (Up)"   = "purple",
               "Other Genes"    = "grey70")
  ) +
  # custom arrows only
  annotate("segment", x = -lim_x, xend =  lim_x, y = 0, yend = 0,
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +
  annotate("segment", x = 0, xend = 0, y = -lim_y, yend =  lim_y,
           size = 1.2, colour = "#053B4C",
           arrow = arrow(type = "closed", length = unit(10, "pt"))) +
  # labels for quadrants
  annotate("text", x = 0, y =  lim_y, label = "Up in ZFY KO", vjust = -0.6, size = 5.0) +
  annotate("text", x = 0, y = -lim_y, label = "Down in ZFY KO", vjust = 1.6, size = 5.0) +
  annotate("text", x =  lim_x, y = 0, label = "Up in LOY AML", hjust = -0.06, size = 5.0) +
  annotate("text", x = -lim_x, y = 0, label = "Down in LOY AML", hjust = 1.06, size = 5.0) +
  # gene labels for overlaps
  geom_text_repel(
    data = subset(merged_degs, Highlight %in% c("Overlap (Down)", "Overlap (Up)")),
    aes(label = Gene, color = Highlight),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf   # allow all labels
  ) +
  coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y), expand = FALSE) +
  labs(
    title = "LOY AML vs ZFY Knockout",
    x = "Log2 Fold Change (LOY AML)",
    y = "Log2 Fold Change (ZFY KO)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(10, 20, 10, 20),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

# save
ggsave(file.path(outdir, "Scatter_Quadrant_LOY_vs_ZFYKO_annotated.pdf"),
       scatter_quad_annot, width = 5, height = 5, device = cairo_pdf)

