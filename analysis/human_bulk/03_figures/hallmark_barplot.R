source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "hallmark_barplot")

# R Script for Visualizing Specific Hallmark GSEA Results
#
# Description:
# This script reads GSEA results for the MSigDB Hallmark gene set collection,
# filters for a predefined list of hallmark pathways, and generates a
# single bar plot to visualize their Normalized Enrichment Scores (NES)
# using manually specified pathway names.
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
# Define a theme for the bar plot suitable for manuscripts
manuscript_theme_for_gsea <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Load and Process Hallmark GSEA Data

cat("Processing Hallmark GSEA results...\n")

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

# Filter for the specified pathways
selected_pathways <- gsea_results_raw %>%
  filter(pathway %in% hallmark_paths_to_plot)

# Prepare Data for Plotting with Hard-coded Names
plot_data <- selected_pathways %>%
  # Create a 'Direction' column for coloring the bars
  mutate(Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>%
  # Hard-code the desired pathway names for the plot axis
  mutate(Term = case_when(
    pathway == "HALLMARK_HEDGEHOG_SIGNALING"        ~ "Hedgehog Signaling",
    pathway == "HALLMARK_TGF_BETA_SIGNALING"        ~ "TGF-beta Signaling",
    pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ "Epithelial Mesenchymal Transition",
    pathway == "HALLMARK_HYPOXIA"                   ~ "Hypoxia",
    pathway == "HALLMARK_KRAS_SIGNALING_UP"         ~ "KRAS Signaling Up",
    pathway == "HALLMARK_MTORC1_SIGNALING"          ~ "mTORC1 Signaling",
    pathway == "HALLMARK_WNT_BETA_CATENIN_SIGNALING" ~ "Wnt/Beta-Catenin Signaling",
    TRUE ~ pathway # Fallback to the original name if no match is found
  )) %>%
  # This ensures the pathways are plotted in the correct order on the y-axis.
  arrange(NES) %>%
  mutate(Term = factor(Term, levels = Term))

cat(paste("Filtered for", nrow(plot_data), "specified hallmark pathways for plotting.\n"))

# Create and Save the Bar Plot

# Generate the bar plot
hallmark_barplot <- ggplot(plot_data, aes(x = NES, y = Term, fill = Direction)) +
  geom_col() + # Use geom_col for pre-calculated values
  # Use a diverging color scheme for up vs. down
  scale_fill_manual(values = c("Up-regulated" = "#B2182B", "Down-regulated" = "#2166AC")) +
  # Add a vertical line at x=0 for reference
  geom_vline(xintercept = 0, color = "grey20") +
  # Add titles and labels
  labs(
    title = "",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  # Apply the custom manuscript theme
  manuscript_theme_for_gsea

# Save the plot as a high-quality PDF
output_filename <- file.path(output_path, "hallmark_barplot.pdf")
ggsave(
  filename = output_filename,
  plot = hallmark_barplot,
  width = 5,  # Adjusted width for potentially longer labels
  height = 4,
  device = cairo_pdf
)

cat(paste("Hallmark pathway bar plot saved to:", output_filename, "\n"))
cat("Script finished successfully.\n")
