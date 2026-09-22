source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "cell_fractions")

# 1. Setup: Load Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)

# 2. Load and Prepare Data

# File paths to the input data
path_to_ciber <- bulk_files$fractions
path_to_meta <- bulk_files$metadata

# Load metadata and create group names by combining Cyto.Group and Sex
meta_df <- read.csv(path_to_meta, stringsAsFactors = F, row.names = 1)
meta_df$Group <- paste0(meta_df$Cyto.Group, "-", meta_df$Sex)

# Load CIBERSORTx results
ciber_res <- read.csv(path_to_ciber)
colnames(ciber_res)[1] <- "sample_id"

# 3. Merge and Tidy Data for Plotting

# Join CIBERSORTx results with metadata
plot_data_wide <- meta_df %>%
  left_join(ciber_res, by = "sample_id")

# Convert the data to a long format
plot_data_long <- plot_data_wide %>%
  pivot_longer(
    cols = Prog.like:GMP,
    names_to = "Cell_Type",
    values_to = "Fraction"
  ) %>%
  # Set the desired order for the groups on the plot's x-axis
  mutate(Group = factor(Group, levels = c("CN_AML-Male", "LOY-Male", "CN_AML-Female")))

# 4. Generate and Save Individual Plots in a Loop

# Get a unique list of all cell types to loop through
all_cell_types <- unique(plot_data_long$Cell_Type)

# Define the pairwise comparisons
my_comparisons <- list(
  c("CN_AML-Male", "LOY-Male"),
  c("CN_AML-Male", "CN_AML-Female")
)

# Define the color mapping
# Blue: CN_AML-Male
# RED: LOY-Male
# Green: CN_AML-Female
custom_colors <- c(
  "CN_AML-Male" = "#0072B2",   # Blue
  "LOY-Male" = "#D55E00",      # Red/Orange
  "CN_AML-Female" = "#009E73"  # Green
)

# Create a directory for box plots
output_dir <- figure_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through each cell type, create a plot, and save it
for (current_cell_type in all_cell_types) {

  # Filter data for the current cell type
  data_for_plot <- plot_data_long %>%
    filter(Cell_Type == current_cell_type)

  p <- ggboxplot(
    data_for_plot,
    x = "Group",
    y = "Fraction",
    fill = "Group",
    add = c("jitter"),
    add.params = list(fill = "white")
  ) +
    # Apply the colors
    scale_fill_manual(values = custom_colors) +

    # This will display the actual p-value instead of stars
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.format"
    ) +
    labs(
      title = current_cell_type,
      x = "Sample Group",
      y = "CIBERSORTx Estimated Fraction"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Create a clean filename
  clean_filename <- str_replace_all(current_cell_type, "[^a-zA-Z0-9_]", "_")

  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0("CIBERSORTx_Boxplot_", clean_filename, ".pdf")),
    plot = p,
    width = 8,
    height = 7
  )
}
