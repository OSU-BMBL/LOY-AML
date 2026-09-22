source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "zfy_pathways")

# Load Libraries

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(stringr)

# Define File Paths
# Path to the directory with the input data files
# Path where the final PDF plots will be saved
output_path <- figure_dir
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Define Helper Functions and Themes

# Function to add newlines to long pathway names for better plotting
add_newlines <- function(string, width = 30) {
  # First, remove prefixes like "GOBP_" for cleaner labels
  string <- gsub("GOBP_", "", string)
  # Then, replace underscores with spaces
  string <- gsub("_", " ", string)
  # Capitalize the first letter of each word
  string <- str_to_title(string)
  # Finally, wrap the text
  wrapped <- strwrap(string, width = width, simplify = FALSE)
  sapply(wrapped, paste, collapse = "\n")
}

# Define a theme for the dot plots suitable for manuscripts
manuscript_theme_for_dots <- theme_classic2(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Dot Plot Generation Function (for GSEA data)
# A reusable function to process data and create a dot plot
generate_gsea_dotplot <- function(filepath, plot_title, output_filename) {

  cat(sprintf("Processing GSEA file: %s\n", filepath))

  gsea_data_raw <- read_excel(filepath)

  gsea_data_processed <- gsea_data_raw %>%
    mutate(
      pathway_label = sapply(pathway, add_newlines),
      gene_count = sapply(str_split(leadingEdge, " "), length),
      log_padj = -log10(padj)
    ) %>%
    arrange(NES) %>%
    mutate(pathway_label = factor(pathway_label, levels = pathway_label))

  dot_plot <- ggplot(gsea_data_processed, aes(x = NES, y = pathway_label)) +
    geom_point(aes(size = gene_count, color = log_padj)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
    scale_size(range = c(4, 10), name = "Gene Count") +
    labs(
      title = plot_title,
      x = "Normalized Enrichment Score",
      y = NULL
    ) +
    manuscript_theme_for_dots

  ggsave(
    filename = file.path(output_path, output_filename),
    plot = dot_plot,
    width = 10,
    height = 18,
    device = cairo_pdf
  )

  cat(sprintf("GSEA plot saved to: %s\n", file.path(output_path, output_filename)))
}

# Dot Plot Generation Function (for Enrichr data)
# Read and plot Enrichr results
generate_enrichr_dotplot <- function(filepath, plot_title, output_filename) {

  cat(sprintf("Processing Enrichr file: %s\n", filepath))

  # Load the enrichment results from the tab-separated text file
  enrichr_data_raw <- read.delim(filepath, sep = "\t")

  # Process the data for plotting
  enrichr_data_processed <- enrichr_data_raw %>%
    # Filter for significant pathways if desired (e.g., Adjusted P-value < 0.05)
    filter(Adjusted.P.value < 0.05) %>%
    # Select the ten highest combined scores
    slice_max(order_by = Combined.Score, n = 10) %>%
    mutate(
      # Clean up pathway names and wrap them for plotting
      pathway_label = sapply(Term, add_newlines),
      # Calculate the number of genes from the 'Genes' column
      gene_count = sapply(str_split(Genes, ";"), length),
      # Calculate -log10(Adjusted.P.value) for color scaling
      log_padj = -log10(Adjusted.P.value)
    ) %>%
    # Order pathways by Odds.Ratio for plotting
    arrange(Odds.Ratio) %>%
    mutate(pathway_label = factor(pathway_label, levels = pathway_label))

  # Create the dot plot
  dot_plot <- ggplot(enrichr_data_processed, aes(x = Odds.Ratio, y = pathway_label)) +
    geom_point(aes(size = gene_count, color = log_padj)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
    scale_size(range = c(4, 10), name = "Gene Count") +
    labs(
      title = plot_title,
      x = "Odds Ratio", # X-axis is now Odds Ratio
      y = NULL
    ) +
    manuscript_theme_for_dots

  # Save the plot as a high-quality PDF
  ggsave(
    filename = file.path(output_path, output_filename),
    plot = dot_plot,
    width = 6,
    height = 5,
    device = cairo_pdf
  )

  cat(sprintf("Enrichr plot saved to: %s\n\n", file.path(output_path, output_filename)))
}

# Generate Plots for Both GSEA Files

# Developmental pathways
generate_gsea_dotplot(
  filepath = bulk_files$go_development,
  plot_title = "Developmental pathways",
  output_filename = "go_development.pdf"
)

# Immune pathways
generate_gsea_dotplot(
  filepath = bulk_files$go_immune,
  plot_title = "Immune pathways",
  output_filename = "go_immune.pdf"
)

# Generate Plot for ZFY Overlap Data

# Plot for the Reactome pathway analysis of overlapping up-regulated genes
generate_enrichr_dotplot(
  filepath = bulk_files$zfy_reactome,
  plot_title = "",
  output_filename = "zfy_overlap_reactome.pdf"
)

cat("Script finished successfully.\n")
