source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "volcano")

# Packages
library(EnhancedVolcano)
library(tidyverse)
library(ggrepel)

# Function to create a volcano plot WITH gene annotations
plot_volcano_annotated <- function(deg_markers, genes_to_color_red, genes_to_color_blue, genes_to_color_grey, plot_title) {

  # Ensure the gene column is named 'gene'
  colnames(deg_markers)[1] <- "gene"

  # Check which of the genes you want to label are actually present in the data.
  genes_to_label <- c(genes_to_color_red, genes_to_color_blue, genes_to_color_grey)
  genes_not_found <- genes_to_label[!genes_to_label %in% deg_markers$gene]

  if (length(genes_not_found) > 0) {
    cat("\n--- NOTE for plot:", plot_title, "---\n")
    cat("The following genes to be labeled were NOT FOUND in the input dataset:\n")
    cat(paste(genes_not_found, collapse = ", "), "\n")
    cat("These genes will not be plotted or labeled. Check the input CSV file.\n---\n")
  }

  # Create a color vector for every gene
  keyvals <- deg_markers %>%
    dplyr::mutate(color = case_when(
      gene %in% genes_to_color_red   ~ "red",
      gene %in% genes_to_color_blue  ~ "blue",
      gene %in% genes_to_color_grey  ~ "grey30",
      padj < 0.05 & log2FoldChange > 1.0  ~ "lightcoral", # Other significant up-regulated
      padj < 0.05 & log2FoldChange < -1.0 ~ "lightblue",    # Other significant down-regulated
      TRUE                               ~ "grey85"
    )) %>%
    dplyr::pull(color)

  names(keyvals) <- deg_markers$gene

  # Generate the volcano plot with labels
  p <- EnhancedVolcano(
    deg_markers,
    lab = deg_markers$gene,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = genes_to_label, # This line ensures specified genes are labeled
    colCustom = keyvals,
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 3.0,
    labSize = 6.0,
    labFace = "bold",
    colAlpha = 0.8,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    colConnectors = 'grey50',
    max.overlaps = Inf,
    title = plot_title,
    subtitle = "",
    legendPosition = 'none'
  ) +
    # Set the y-axis title with proper formatting
    ylab(expression(-log[10]~italic(p)[adj])) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 22, hjust = 0.5),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      # Make y-axis title bold
      axis.title.y = element_text(face = "bold")
    )

  return(p)
}

# Function to create a volcano plot WITHOUT gene annotations
plot_volcano_unannotated <- function(deg_markers, genes_to_color_red, genes_to_color_blue, genes_to_color_grey, plot_title) {

  # Ensure the gene column is named 'gene'
  colnames(deg_markers)[1] <- "gene"

  # Create a color vector for every gene
  keyvals <- deg_markers %>%
    dplyr::mutate(color = case_when(
      gene %in% genes_to_color_red   ~ "red",
      gene %in% genes_to_color_blue  ~ "blue",
      gene %in% genes_to_color_grey  ~ "grey30",
      padj < 0.05 & log2FoldChange > 1.0  ~ "lightcoral", # Other significant up-regulated
      padj < 0.05 & log2FoldChange < -1.0 ~ "lightblue",    # Other significant down-regulated
      TRUE                               ~ "grey85"
    )) %>%
    dplyr::pull(color)

  names(keyvals) <- deg_markers$gene

  # Generate the volcano plot WITHOUT labels
  p <- EnhancedVolcano(
    deg_markers,
    lab = deg_markers$gene,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = c(),
    colCustom = keyvals,
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 3.0,
    labSize = 0.0,
    labFace = "bold",
    colAlpha = 0.8,
    drawConnectors = FALSE,
    widthConnectors = 0.75,
    colConnectors = 'grey50',
    max.overlaps = Inf,
    title = plot_title,
    subtitle = "",
    legendPosition = 'none'
  ) +
    # Set the y-axis title with proper formatting
    ylab(expression(-log[10]~italic(p)[adj])) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 22, hjust = 0.5),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      # Make y-axis title bold
      axis.title.y = element_text(face = "bold")
    )

  return(p)
}

# Define Gene Lists

genes_to_annotate_loy_cn_male <- c("MECOM", "SOX17", "OLIG2", "ID4")
genes_Y_chromosome <- c("DDX3Y", "EIF1AY", "KDM5D", "USP9Y", "ZFY")
genes_X_inactivation <- c("XIST", "TSIX", "CRLF2", "CD1C", "STC1", "LILRA4")

# Define File Paths

base_path <- bulk_files$volcano
output_path <- figure_dir
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Generate Plots for Each Comparison

### Comparison 1: LOY Male vs CN Male
cat("\nProcessing: LOY Male vs CN Male...\n")
comp_name_1 <- "LOY-Male_vs_CN-Male"
deg_results_1 <- read.csv(file.path(base_path, "loy_vs_cn_male.csv"))
cat("Loaded", nrow(deg_results_1), "total genes for", comp_name_1, "\n")

# Annotated plot
annotated_plot_1 <- plot_volcano_annotated(
  deg_markers = deg_results_1,
  genes_to_color_red = genes_to_annotate_loy_cn_male,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "LOY Male vs CN Male"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_1, "_annotated_volcano.pdf")),
  plot = annotated_plot_1,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

# Unannotated plot
unannotated_plot_1 <- plot_volcano_unannotated(
  deg_markers = deg_results_1,
  genes_to_color_red = genes_to_annotate_loy_cn_male,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "LOY Male vs CN Male"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_1, "_unannotated_volcano.pdf")),
  plot = unannotated_plot_1,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

### Comparison 2: LOY Male vs CN Female
cat("\nProcessing: LOY Male vs CN Female...\n")
comp_name_2 <- "LOY-Male_vs_CN-Female"
deg_results_2 <- read.csv(file.path(base_path, "loy_vs_cn_female.csv"))
cat("Loaded", nrow(deg_results_2), "total genes for", comp_name_2, "\n")

# Annotated plot
annotated_plot_2 <- plot_volcano_annotated(
  deg_markers = deg_results_2,
  genes_to_color_red = genes_Y_chromosome,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "LOY Male vs CN Female"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_2, "_annotated_volcano.pdf")),
  plot = annotated_plot_2,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

# Unannotated plot
unannotated_plot_2 <- plot_volcano_unannotated(
  deg_markers = deg_results_2,
  genes_to_color_red = genes_Y_chromosome,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "LOY Male vs CN Female"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_2, "_unannotated_volcano.pdf")),
  plot = unannotated_plot_2,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

### Comparison 3: CN Male vs CN Female
cat("\nProcessing: CN Male vs CN Female...\n")
comp_name_3 <- "CN-Male_vs_CN-Female"
deg_results_3 <- read.csv(file.path(base_path, "cn_male_vs_cn_female.csv"))
cat("Loaded", nrow(deg_results_3), "total genes for", comp_name_3, "\n")

# Annotated plot
annotated_plot_3 <- plot_volcano_annotated(
  deg_markers = deg_results_3,
  genes_to_color_red = genes_Y_chromosome,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "CN Male vs CN Female"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_3, "_annotated_volcano.pdf")),
  plot = annotated_plot_3,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

# Unannotated plot
unannotated_plot_3 <- plot_volcano_unannotated(
  deg_markers = deg_results_3,
  genes_to_color_red = genes_Y_chromosome,
  genes_to_color_blue = genes_X_inactivation,
  genes_to_color_grey = c(),
  plot_title = "CN Male vs CN Female"
)
ggsave(
  filename = file.path(output_path, paste0(comp_name_3, "_unannotated_volcano.pdf")),
  plot = unannotated_plot_3,
  width = 6, height = 6, units = "in", device = cairo_pdf
)

cat("\nScript finished. Both annotated and unannotated plots have been saved as PDF files.\n")

