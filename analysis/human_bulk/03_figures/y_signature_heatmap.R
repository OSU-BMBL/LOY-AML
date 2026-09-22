source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "y_signature_heatmap")

# 1. Load Libraries
library(tidyverse)
library(dplyr)
library(edgeR)
library(pheatmap)
library(GSVA)
library(grid)
library(RColorBrewer)

# 2. Load Data
path_to_counts <- bulk_files$counts
path_to_meta <- bulk_files$metadata
path_to_ciber <- bulk_files$fractions
path_to_loy_genes <- bulk_files$signature

counts_df <- read.csv(path_to_counts, row.names = 1)
meta_df <- read.csv(path_to_meta, stringsAsFactors = FALSE)
ciber_res <- read.csv(path_to_ciber)

#colnames(meta_df)[1] <- "sample_id"
colnames(ciber_res)[1] <- "sample_id"
loy_genes <- readLines(path_to_loy_genes)

# 3. Filter Metadata for LOY and CN Males
meta_filtered <- meta_df %>%
  dplyr::mutate(Group = dplyr::case_when(
    Cyto.Group == "CN_AML" & Sex == "Male"   ~ "CN Males",
    Cyto.Group == "LOY"    & Sex == "Male"   ~ "LOY Males",
    TRUE                                     ~ "Other"
  )) %>%
  dplyr::filter(Group %in% c("LOY Males", "CN Males"))

combined_samples <- meta_filtered$sample_id

# 4. Prepare Expression Matrix for the Combined Cohort
counts_filtered <- counts_df[, combined_samples]
counts_filtered[counts_filtered < 0] <- 0

dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# 5. Calculate and Scale LOY Signature Score on the Combined Cohort
y_signature_list <- list("Y_Signature" = loy_genes)

gsvapar <- ssgseaParam(expr = as.matrix(logcpm), geneSets = y_signature_list)
gsva_scores <- gsva(gsvapar, verbose = FALSE)

loy_score_gsva <- gsva_scores[1, ]
loy_score_scaled <- scale(loy_score_gsva)

cat("Y-Signature scores calculated and scaled using LOY and CN Male samples.\n")

# 6. Match CIBERSORTx Fractions
frac_mat <- ciber_res %>%
  dplyr::filter(sample_id %in% combined_samples) %>%
  tibble::column_to_rownames("sample_id") %>%
  dplyr::select(-c(P.value, Correlation, RMSE))

frac_mat <- frac_mat[names(loy_score_gsva), ]

# 7. Compute Estimated LOY per Cell Type
loy_score_celltype_full <- frac_mat * as.vector(loy_score_scaled)

# 8. Prepare Final Data for Combined Heatmap
annotation_df <- meta_filtered %>%
  dplyr::select(sample_id, Group, Age, Race, Batch) %>%
  tibble::column_to_rownames("sample_id")

annotation_df <- annotation_df[rownames(loy_score_celltype_full), ]

# 9. Heatmap Visualization
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- grid::textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"),
                        vjust = 0.5, hjust = 1, rot = 45, gp = grid::gpar(...))
  return(res)
}
assignInNamespace("draw_colnames", draw_colnames_45, ns = asNamespace("pheatmap"))

ann_colors <- list(
  Group = c("LOY Males" = "#D55E00", "CN Males" = "#0072B2")
)

# CHANGE: Use the svg() device for saving
# Step 1: Open the SVG file device
svg("Combined_LOY_and_CN_Heatmap.svg", width = 10, height = 24)

# Step 2: Create the heatmap. It will be drawn to the open SVG file.
# Note that the 'filename' argument has been removed.
pheatmap::pheatmap(
  loy_score_celltype_full,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Estimated LOY Signature per Cell Type (LOY vs. CN Males)",
  fontsize_row = 8
)

# Step 3: Close the device to finalize and save the file.
dev.off()

cat("Combined heatmap saved as 'Combined_LOY_and_CN_Heatmap.svg'.\n")

# Save to PDF
cat("Saving heatmap to PDF file...\n")
pdf(file.path(figure_dir, "y_signature_sample_heatmap.pdf"), width = 10, height = 24)
pheatmap::pheatmap(
  loy_score_celltype_full,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Estimated LOY Signature per Cell Type (LOY vs. CN Males)",
  fontsize_row = 8
)
dev.off() # Close the PDF device
cat("Combined heatmap saved as 'y_signature_sample_heatmap.pdf'.\n")

cat("--- Creating 3rd Heatmap: Summarized View ---\n")

# 10. Summarize Data by Group
# Combine the scores and the group annotation
data_to_summarize <- as.data.frame(loy_score_celltype_full)
data_to_summarize$Group <- annotation_df$Group

# Calculate the mean score for each cell type, grouped by "Group"
summary_mat <- data_to_summarize %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  tibble::column_to_rownames("Group") # pheatmap needs row names

# Ensure the order is consistent for comparison
summary_mat <- summary_mat[c("LOY Males", "CN Males"), ]

# 11. Generate the Summarized Heatmap
# This heatmap will have only two rows: one for LOY Males (average)

# Open PDF device
pdf(file.path(figure_dir, "y_signature_group_heatmap.pdf"), width = 10, height = 4)

pheatmap::pheatmap(
  summary_mat,
  cluster_rows = FALSE, # We only have two rows, no need to cluster
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Mean Estimated LOY Signature per Cell Type\n(LOY vs. CN Males)",
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45,
  display_numbers = TRUE, # Display the average values directly on the heatmap
  number_format = "%.2f"   # Format numbers to two decimal places
)

# Close the device
dev.off()
cat("Summarized heatmap saved as 'y_signature_group_heatmap.pdf'.\n")

# A boxplot is a great way to see the distribution that the mean summarizes.
# It shows variability, outliers, and statistical difference more clearly.

cat("--- Creating Bonus Visualization: Boxplots ---\n")

# Convert the summarized data back to a "long" format suitable for ggplot
long_df <- data_to_summarize %>%
  tidyr::pivot_longer(
    cols = -Group,
    names_to = "Cell_Type",
    values_to = "Estimated_LOY_Score"
  )

# Create the plot
boxplot_comparison <- ggplot(long_df, aes(x = Cell_Type, y = Estimated_LOY_Score, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + # Hiding outliers for clarity, geom_jitter can show them
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.3, size=1.5) +
  scale_fill_manual(values = c("LOY Males" = "#D55E00", "CN Males" = "#0072B2")) +
  theme_bw() +
  labs(
    title = "Distribution of Estimated LOY Signature by Cell Type",
    x = "Cell Type",
    y = "Estimated LOY Score (Scaled)",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  )

# Save the boxplot
ggsave(file.path(figure_dir, "y_signature_boxplots.pdf"), plot = boxplot_comparison, width = 12, height = 7)
cat("Comparison boxplots saved as 'y_signature_boxplots.pdf'.\n")

# We will use 'loy_score_celltype_full', 'annotation_df', and 'summary_mat'

cat("--- Creating 4th Heatmap: Individual LOY vs. Pooled CN ---\n")

# 13. Prepare Data for the Hybrid Heatmap

# First, get the matrix for only the individual LOY samples
loy_individual_mat <- loy_score_celltype_full[annotation_df$Group == "LOY Males", ]

# Second, get the single row for the pooled (mean) CN samples from our previous summary
# We use drop = FALSE to ensure it stays as a matrix with 1 row
cn_pooled_row <- summary_mat["CN Males", , drop = FALSE]

# To make the plot clear, let's give this pooled row a descriptive name
rownames(cn_pooled_row) <- "CN Males (Pooled Mean)"

# Now, combine the individual LOY samples and the single pooled CN row
hybrid_mat <- rbind(loy_individual_mat, cn_pooled_row)

# 14. Create a new Annotation for this Hybrid Heatmap
# We need to distinguish the individual samples from the pooled reference row.

# Create a data frame for the annotation
hybrid_annotation_df <- data.frame(
  SampleType = c(
    rep("Individual LOY", nrow(loy_individual_mat)),
    "Pooled CN"
  )
)
# The row names must match the matrix we just created
rownames(hybrid_annotation_df) <- rownames(hybrid_mat)

# Define colors for our new annotation
hybrid_ann_colors <- list(
  SampleType = c("Individual LOY" = "#D55E00", "Pooled CN" = "black")
)

# 15. Generate and Save the Hybrid Heatmap

# Open PDF device for saving
pdf(file.path(figure_dir, "y_signature_pooled_cn_heatmap.pdf"), width = 10, height = 9)

pheatmap::pheatmap(
  hybrid_mat,
  cluster_rows = TRUE,  # Clustering rows is very useful here!
  cluster_cols = TRUE,
  annotation_row = hybrid_annotation_df,
  annotation_colors = hybrid_ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Individual LOY Samples vs. Pooled CN-Male Baseline",
  fontsize_row = 14,
  fontsize_col = 14,
  angle_col = 45
  # We remove the cell numbers here as they would only show for the one pooled row
  # and make the plot look cluttered.
  # display_numbers = FALSE
)

# Close the device
dev.off()

cat("Hybrid heatmap saved as 'y_signature_pooled_cn_heatmap.pdf'.\n")

# We will use 'hybrid_mat', 'frac_mat', and 'meta_filtered' from the previous steps.

cat("--- Creating 5th Heatmap: Merging T and CTL columns ---\n")

# 16. Prepare data by merging T and CTL columns

# Define the columns to merge
t_cell_cols <- c("T", "CTL")

# Get the original cell fractions for all relevant samples (LOY and CN Males)
original_fractions <- ciber_res %>%
  dplyr::filter(sample_id %in% rownames(loy_score_celltype_full)) %>%
  tibble::column_to_rownames("sample_id")

# First, calculate the average fractions for T and CTL across all CN Males.
# This is needed for weighting the "CN Males (Pooled Mean)" row.
cn_male_samples <- meta_filtered %>%
  filter(Group == "CN Males") %>%
  pull(sample_id)
avg_cn_fractions <- colMeans(original_fractions[cn_male_samples, t_cell_cols])

# Initialize an empty vector to store the new, combined scores
new_combined_col <- numeric(nrow(hybrid_mat))
names(new_combined_col) <- rownames(hybrid_mat)

# Loop through each row of the hybrid matrix (each sample)
for (sample_name in rownames(hybrid_mat)) {

  # Get the scores for T and CTL for the current sample
  scores_to_combine <- hybrid_mat[sample_name, t_cell_cols]

  # Check if it's the special pooled row or an individual sample
  if (sample_name == "CN Males (Pooled Mean)") {
    # For the pooled row, use the average CN fractions as weights
    fractions_to_use <- avg_cn_fractions
  } else {
    # For individual LOY samples, use their own specific fractions
    fractions_to_use <- original_fractions[sample_name, t_cell_cols]
  }

  # Calculate the weighted average: sum(score * fraction) / sum(fraction)
  total_fraction <- sum(fractions_to_use)

  if (total_fraction > 0) {
    weighted_avg <- sum(scores_to_combine * fractions_to_use) / total_fraction
  } else {
    # If both fractions are 0, the combined score is 0
    weighted_avg <- 0
  }

  # Store the result
  new_combined_col[sample_name] <- weighted_avg
}

# Create the new matrix for the heatmap
# 1. Remove the old 'T' and 'CTL' columns
mat_without_t_cells <- hybrid_mat[, !colnames(hybrid_mat) %in% t_cell_cols]

# 2. Add the new combined column
hybrid_mat_merged <- cbind(mat_without_t_cells, "T cells" = new_combined_col)

cat("Data prepared with T and CTL columns merged.\n")

# 17. Generate and Save the Merged Heatmap

# Define the custom order for the individual samples
custom_sample_order <- c(
  "C_12_4569", "C_99_1502", "C_12_6815", "C_98_0327", "XAKMT02",
  "XAKMT80", "C_98_0705", "XAKMT25", "C_06_0996", "XAKMT74",
  "TDP93TDPC014", "TDPC050564", "TDPPS920020", "PS89_0219", "XAKMT16",
  "C14_93_C_014", "XAKMT35", "PS88_0279", "TDP93TDPC134", "XAKMT21",
  "C_09_4183", "C_08_1102", "TDPC081084"
)

# Combine the custom sample order with the pooled mean row at the end
final_row_order <- c(custom_sample_order, "CN Males (Pooled Mean)")

# Reorder the matrix and annotation data frame based on the final order
hybrid_mat_ordered <- hybrid_mat_merged[final_row_order, ]
hybrid_annotation_ordered <- hybrid_annotation_df[final_row_order, , drop = FALSE]

# Open PDF device for saving
pdf(file.path(figure_dir, "y_signature_pooled_cn_merged_t_heatmap.pdf"), width = 10, height = 9)

pheatmap::pheatmap(
  hybrid_mat_ordered,        # Use the reordered matrix
  cluster_rows = FALSE,      # Disable row clustering
  cluster_cols = TRUE,
  annotation_row = hybrid_annotation_ordered, # Use the reordered annotation
  annotation_colors = hybrid_ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Individual LOY vs. Pooled CN-Male",
  fontsize_row = 8,         # Adjusted font size for readability
  fontsize_col = 12,
  angle_col = 45
)

# Close the device
dev.off()

cat("Merged T/CTL heatmap with custom order saved as 'y_signature_pooled_cn_merged_t_heatmap.pdf'.\n")
