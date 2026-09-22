source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "deconvolution_correlations")

# 1. Setup: Load Libraries
library(tidyverse)
library(corrr)
library(ggcorrplot)
library(patchwork)

# 2. Load and Prepare Data
# Define file paths
path_to_ciber <- bulk_files$fractions
path_to_meta <- bulk_files$metadata

# Load metadata and create group names
meta_df <- read.csv(path_to_meta, stringsAsFactors = F, row.names = 1) %>%
  mutate(Group = paste0(Cyto.Group, "-", Sex)) %>%
  filter(Group %in% c("CN_AML-Male", "LOY-Male", "CN_AML-Female"))

# Load CIBERSORTx results
ciber_res <- read.csv(path_to_ciber)
colnames(ciber_res)[1] <- "sample_id"

# Merge T and CTL columns into "T cells"
t_cell_cols <- c("T", "CTL")
# Check if both columns exist before attempting to merge
if (all(t_cell_cols %in% colnames(ciber_res))) {
  ciber_res <- ciber_res %>%
    # Create the new "T cells" column by summing the old ones
    mutate(`T cells` = T + CTL) %>%
    # Remove the original 'T' and 'CTL' columns
    select(-all_of(t_cell_cols))

  cat("Successfully merged 'T' and 'CTL' columns into 'T cells'.\n")
} else {
  cat("Warning: 'T' and/or 'CTL' columns not found. Skipping merge step.\n")
}

# Explicitly define cell type columns after potential merge
all_ciber_cols <- colnames(ciber_res)
cell_type_cols <- setdiff(all_ciber_cols, c("sample_id", "P.value", "Correlation", "RMSE"))

# Join CIBERSORTx results with metadata
full_data <- meta_df %>%
  left_join(ciber_res, by = "sample_id")

# 3. Determine the Master Order from the LOY Group
# Isolate the LOY-Male data
loy_df <- full_data %>%
  filter(Group == "LOY-Male") %>%
  select(all_of(cell_type_cols)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0)) # Remove zero-variance columns

# Calculate correlation and clustering for the LOY group
loy_corr_matrix <- cor(loy_df, method = "pearson", use = "pairwise.complete.obs")
loy_hclust <- hclust(dist(loy_corr_matrix), method = "complete")
master_cell_order <- rownames(loy_corr_matrix)[loy_hclust$order]

# 4. Function to Create Heatmaps with a Fixed Order
create_ordered_corr_plot <- function(df, group_name, cell_order) {
  # Select cell types and remove zero-variance columns for this specific group
  cell_type_df <- df %>%
    select(all_of(cell_type_cols)) %>%
    select(where(~ sd(.x, na.rm = TRUE) > 0))

  # Check if there are at least 2 columns left to correlate
  if (ncol(cell_type_df) < 2) {
    return(
      ggplot() +
        annotate("text", x = 1, y = 1, size = 4, label = "Not enough varying cell types\nto create a correlation matrix.") +
        ggtitle(group_name) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  }

  # Calculate correlation matrix
  corr_data <- cor(cell_type_df, method = "pearson", use = "pairwise.complete.obs")

  # Use the master order, but only for cells present in this group's matrix
  valid_order <- cell_order[cell_order %in% rownames(corr_data)]
  corr_data_ordered <- corr_data[valid_order, valid_order]

  # Create the heatmap plot with clustering turned OFF
  ggcorrplot(
    corr_data_ordered,
    hc.order = FALSE, # Use the provided order
    type = "lower",
    lab = TRUE,
    lab_size = 2.5,
    method = "circle",
    colors = c("#6D9EC1", "white", "#E46726"),
    title = group_name,
    ggtheme = theme_minimal()
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# 5. Generate Plots for Each Group Using the Master Order
# Split data by group
list_of_dfs <- full_data %>%
  group_split(Group)

if (length(list_of_dfs) > 0) {
  names(list_of_dfs) <- sapply(list_of_dfs, function(df) unique(df$Group))
}

# Create plots, passing the master_cell_order to each
plot_loy_male <- if ("LOY-Male" %in% names(list_of_dfs)) create_ordered_corr_plot(list_of_dfs[["LOY-Male"]], "LOY Male", master_cell_order) else NULL
plot_cn_male <- if ("CN_AML-Male" %in% names(list_of_dfs)) create_ordered_corr_plot(list_of_dfs[["CN_AML-Male"]], "CN-AML Male", master_cell_order) else NULL
plot_cn_female <- if ("CN_AML-Female" %in% names(list_of_dfs)) create_ordered_corr_plot(list_of_dfs[["CN_AML-Female"]], "CN-AML Female", master_cell_order) else NULL

# 6. Combine and Save the Final Figure
# Combine only the non-null plots
plot_list <- compact(list(plot_loy_male, plot_cn_male, plot_cn_female))

if (length(plot_list) > 0) {
  final_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom')

  # Save as PDF
  ggsave(
    file.path(figure_dir, "cell_fraction_correlations.pdf"),
    plot = final_plot,
    width = 7 * length(plot_list),
    height = 7,
    device = cairo_pdf # Using cairo_pdf for better font embedding
  )

  message("Saved correlation heatmap to ", figure_dir)
} else {
  cat("Script finished. No plots were generated as no data was available for the specified groups.\n")
}
