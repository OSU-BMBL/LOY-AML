source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "hallmark_dotplot")

# Hallmark enrichment dot plot
#
# Description:
# This script reads GSEA results for the MSigDB Hallmark gene set collection,
# filters for a predefined list of hallmark pathways, and generates a
# single dot plot to visualize their results. The plot uses:
# - X-axis: Normalized Enrichment Score (NES)
# - Y-axis: Pathway Name
# - Point Size: Number of genes in the leading edge (Gene Count)
# - Point Color: Statistical significance (-log10 of the adjusted p-value)
#

# Load Libraries

library(tidyverse)
library(ggplot2)
library(readr)

# Define File Paths
# Path to the GSEA results CSV file
input_file_path <- bulk_files$hallmark

# Path where the final PDF plot will be saved
output_path <- figure_dir
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Define a Custom Theme
# Define a theme for the dot plot suitable for manuscripts
manuscript_theme_for_dots <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Load and Process Hallmark GSEA Data

cat("Processing Hallmark GSEA results for dot plot...\n")

# Load the GSEA results from the CSV file
gsea_results_raw <- read_csv(input_file_path)

# Define the specific Hallmark pathways to be plotted
hallmark_paths_to_plot <- c(
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

# Filter for the specified pathways and prepare data for plotting
plot_data <- gsea_results_raw %>%
  # Keep only the rows corresponding to the pathways of interest
  filter(pathway %in% hallmark_paths_to_plot) %>%
  # Create new columns for plotting and apply custom names
  mutate(
    # Create pretty pathway labels for the y-axis
    Term = case_when(
      pathway == "HALLMARK_HEDGEHOG_SIGNALING"        ~ "Hedgehog Signaling",
      pathway == "HALLMARK_TGF_BETA_SIGNALING"        ~ "TGF-beta Signaling",
      pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ "Epithelial Mesenchymal Transition",
      pathway == "HALLMARK_HYPOXIA"                   ~ "Hypoxia",
      pathway == "HALLMARK_KRAS_SIGNALING_UP"         ~ "KRAS Signaling Up",
      pathway == "HALLMARK_MTORC1_SIGNALING"          ~ "mTORC1 Signaling",
      pathway == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"~ "Wnt/Beta-Catenin Signaling",
      TRUE ~ pathway # Fallback to the original name if no match
    ),
    # Calculate the number of genes in the leading edge.
    # GSEA 'leadingEdge' genes can be separated by spaces, slashes, or other delimiters.
    # Adjust str_split(" ") if the file uses a different separator (e.g., "/").
    gene_count = sapply(str_split(leadingEdge, " "), length),
    # Calculate -log10(padj) for the color scale.
    # Add a small constant to avoid log10(0) if padj is 0.
    log_padj = -log10(padj + 1e-10)
  ) %>%
  arrange(NES) %>%
  # Convert 'Term' to a factor with levels set by the new order.
  # This locks in the sorting for ggplot.
  mutate(Term = factor(Term, levels = Term))

cat(paste("Filtered for", nrow(plot_data), "specified hallmark pathways for plotting.\n"))

# Create and Save the Dot Plot

# Generate the dot plot
hallmark_dotplot <- ggplot(plot_data, aes(x = NES, y = Term)) +
  # Add points with size mapped to gene count and color to significance
  geom_point(aes(size = gene_count, color = log_padj)) +
  # Use a gradient color scale for the -log10(padj) values
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "-log10(p.adj)") +
  # Define the range of point sizes
  scale_size(range = c(4, 10), name = "Gene Count") +
  # Add a vertical dashed line at x=0 for reference
  geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
  # Add titles and labels
  labs(
    title = "",
    x = "Normalized Enrichment Score (NES)",
    y = NULL # Remove y-axis title
  ) +
  # Apply the custom manuscript theme
  manuscript_theme_for_dots

# Save the plot as a high-quality PDF
output_filename <- file.path(output_path, "hallmark_dotplot.pdf")
ggsave(
  filename = output_filename,
  plot = hallmark_dotplot,
  width = 7,  # Adjusted width to accommodate the legend
  height = 4.5,
  device = cairo_pdf
)

cat(paste("Hallmark pathway dot plot saved to:", output_filename, "\n"))
cat("Script finished successfully.\n")
