source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "go_enrichment")

# Load Libraries

library(tidyverse)
library(ggplot2)
library(readxl)

# Define File Paths
# Curated pathway workbook directory.
# Set the directory where the final PDF plots will be saved
output_directory <- figure_dir

# Create the output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Define Helper Functions and Themes

#' Clean and format pathway names for plotting.
#'
#' This function removes common prefixes, replaces underscores with spaces,
#' converts the string to sentence case, and then wraps long lines.
#' @param pathway_name The raw pathway name string.
#' @param wrap_width The character width at which to wrap the string.
#' @return A formatted, wrapped string suitable for plot labels.
format_pathway_name <- function(pathway_name, wrap_width = 50) {
  cleaned_name <- pathway_name %>%
    str_remove("^(GOBP|KEGG|REACTOME)_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()

  wrapped <- strwrap(cleaned_name, width = wrap_width, simplify = FALSE)
  sapply(wrapped, paste, collapse = "\n")
}

# Define a consistent theme for the bar plots
manuscript_theme <- theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

# Define a theme specifically for dot plots
manuscript_theme_for_dots <- theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, color = "black", hjust = 1),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

# GSEA Plotting Functions

#' Generate a GSEA bar plot from a single file.
#' @param filepath Path to the GSEA results Excel file.
#' @param plot_title The main title for the plot.
#' @param output_filename Name for the output PDF file.
generate_gsea_barplot <- function(filepath, plot_title, output_filename) {

  cat(paste0("\n--- Generating GSEA Bar Plot: ", plot_title, " ---\n"))

  gsea_data <- read_excel(filepath) %>%
    select(pathway, NES) %>%
    mutate(Direction = ifelse(NES > 0, "Up in LOY", "Down in LOY")) %>%
    mutate(pathway_fmt = map_chr(pathway, format_pathway_name)) %>%
    arrange(NES) %>%
    mutate(pathway_fmt = factor(pathway_fmt, levels = unique(.$pathway_fmt)))

  # Generate the plot
  gsea_barplot <- ggplot(gsea_data, aes(x = NES, y = pathway_fmt, fill = Direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    scale_fill_manual(name = "Direction", values = c("Up in LOY" = "#B2182B", "Down in LOY" = "#2166AC")) +
    labs(
      title = plot_title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    manuscript_theme

  # Save the plot
  plot_height <- 0.25 * nrow(gsea_data) + 2.0
  output_path <- file.path(output_directory, output_filename)

  ggsave(
    filename = output_path,
    plot = gsea_barplot,
    width = 8,
    height = max(plot_height, 5),
    device = cairo_pdf,
    limitsize = FALSE
  )
  cat(paste("  Saved GSEA bar plot to:", output_path, "\n"))
}

#' Generate a GSEA dot plot from a single file.
#' @param filepath Path to the GSEA results Excel file.
#' @param plot_title The main title for the plot.
#' @param output_filename Name for the output PDF file.
generate_gsea_dotplot <- function(filepath, plot_title, output_filename) {

  cat(paste0("\n--- Generating GSEA Dot Plot: ", plot_title, " ---\n"))

  gsea_data <- read_excel(filepath) %>%
    # Select the 'size' column as well
    select(pathway, padj, NES, leadingEdge, size) %>%
    mutate(
      pathway_label = map_chr(pathway, format_pathway_name),
      # Calculate gene count from leadingEdge
      gene_count = sapply(str_split(str_trim(leadingEdge), "\\s+"), length),
      # Calculate the gene ratio
      gene_ratio = gene_count / size,
      log_padj = -log10(padj),
      # Determine pathway color based on NES direction
      pathway_color = ifelse(NES > 0, "#B2182B", "#2166AC"),
      # Direction label
      Direction = ifelse(NES > 0, "Up in LOY", "Down in LOY")
    ) %>%
    arrange(NES) %>%
    mutate(pathway_label = factor(pathway_label, levels = unique(.$pathway_label))) %>%
    # Ensure Direction is a factor with the correct order for the legend
    mutate(Direction = factor(Direction, levels = c("Up in LOY", "Down in LOY")))

  # Generate the plot
  dot_plot <- ggplot(gsea_data, aes(x = NES, y = pathway_label)) +
    # Add 'shape = Direction' to aesthetics to create a basis for the legend
    # Change size aesthetic to use gene_ratio
    geom_point(aes(size = gene_ratio, color = log_padj, shape = Direction)) +

    # Adjusted P-value colour scale
    scale_color_viridis_c(option = "magma", name = bquote(~-log[10]~'(p.adj)')) +

    # Gene-ratio scale
    scale_size(range = c(4, 10), name = "Gene Ratio") +

    # Direction legend
    scale_shape_manual(name = "Direction", values = c("Up in LOY" = 19, "Down in LOY" = 19)) +

    # Override the color in the 'shape' legend to show direction colors
    guides(shape = guide_legend(override.aes = list(color = c("Up in LOY" = "#B2182B", "Down in LOY" = "#2166AC")))) +

    labs(
      title = plot_title,
      x = "Normalized Enrichment Score",
      y = NULL
    ) +
    manuscript_theme_for_dots +
    # Override the theme to color the y-axis labels based on direction
    theme(axis.text.y = element_text(color = gsea_data$pathway_color))

  # Save the plot
  plot_height <- 0.25 * nrow(gsea_data) + 2.0
  output_path <- file.path(output_directory, output_filename)

  ggsave(
    filename = output_path,
    plot = dot_plot,
    width = 8,
    height = max(plot_height, 6),
    device = cairo_pdf,
    limitsize = FALSE
  )
  cat(paste("  Saved GSEA dot plot to:", output_path, "\n"))
}

# Execute the GSEA Analysis Script

# Define the input file for the GSEA analysis
gsea_file_leukemic <- bulk_files$go_pathways

# Generate the GSEA plots
cat("\nStarting GSEA analysis for Leukemic Cells...\n")

# Generate the GSEA bar plot
generate_gsea_barplot(
  filepath = gsea_file_leukemic,
  plot_title = "GSEA: LOY vs CN Males",
  output_filename = "go_barplot.pdf"
)

# Generate the GSEA dot plot
generate_gsea_dotplot(
  filepath = gsea_file_leukemic,
  plot_title = "GSEA: LOY vs CN Males",
  output_filename = "go_dotplot.pdf"
)

cat("\nGSEA analysis script finished successfully.\n")

