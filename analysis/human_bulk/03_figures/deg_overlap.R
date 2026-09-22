source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "deg_overlap")

# Load necessary libraries
library(tidyverse)
library(VennDiagram)
library(grid) # Needed for drawing the plot to the PDF device

# Step 1: Load the DEG data for the three comparisons
deg_loy_vs_cn_male_file <- file.path(bulk_files$deg_overlap, "loy_vs_cn_male.csv")
deg_loy_vs_cn_female_file <- file.path(bulk_files$deg_overlap, "loy_vs_cn_female.csv")
deg_male_vs_female_file <- file.path(bulk_files$deg_overlap, "cn_male_vs_cn_female.csv")

deg1 <- read.csv(deg_loy_vs_cn_male_file)
deg2 <- read.csv(deg_loy_vs_cn_female_file)
deg3 <- read.csv(deg_male_vs_female_file)

# Step 2: Define a function to extract DEGs by direction
get_degs_by_direction <- function(df, direction) {
  df_filtered <- df %>% filter(padj < 0.05)
  if (direction == "up") {
    df_filtered <- df_filtered %>% filter(log2FoldChange > 0)
  } else if (direction == "down") {
    df_filtered <- df_filtered %>% filter(log2FoldChange < 0)
  }
  df_filtered %>% pull(X)
}

# Define the color palette
color_palette <- c("#8A2BE2", "#FF8C00", "#008080")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Part 1: Venn Diagram for UP-REGULATED Genes

cat("============================================\n")
cat("ANALYSIS FOR UP-REGULATED GENES\n")
cat("============================================\n\n")

up_genes1 <- get_degs_by_direction(deg1, "up")
up_genes2 <- get_degs_by_direction(deg2, "up")
up_genes3 <- get_degs_by_direction(deg3, "up")

# DIAGNOSTIC STEP
# Check how many genes were found for each list
cat("Number of UP-regulated genes in 'LOY Male vs CN Male':", length(up_genes1), "\n")
cat("Number of UP-regulated genes in 'LOY Male vs CN Female':", length(up_genes2), "\n")
cat("Number of UP-regulated genes in 'CN Male vs CN Female':", length(up_genes3), "\n\n")

up_gene_list_for_venn <- list(
  `LOY Male vs CN Male` = up_genes1,
  `LOY Male vs CN Female` = up_genes2,
  `CN Male vs CN Female` = up_genes3
)

# Generate the Venn diagram object (without saving to file yet)
venn_object_up <- venn.diagram(
  x = up_gene_list_for_venn,
  filename = NULL, # Set filename to NULL to return a plot object
  category.names = c("LOY Male vs CN Male", "LOY Male vs CN Female", "CN Male vs CN Female"),
  cat.cex = 1.1, cat.fontface = "bold", cat.col = color_palette,
  fill = color_palette, alpha = 0.5, cex = 1.2, fontface = "plain",
  margin = 0.2,
  cat.dist = c(0.06, 0.06, 0.06)
)

# Save to PDF
pdf(file.path(figure_dir, "deg_overlap_upregulated.pdf"), height = 7, width = 7)
grid.draw(venn_object_up)
dev.off()

# Part 2: Venn Diagram for DOWN-REGULATED Genes

cat("\n============================================\n")
cat("ANALYSIS FOR DOWN-REGULATED GENES\n")
cat("============================================\n\n")

down_genes1 <- get_degs_by_direction(deg1, "down")
down_genes2 <- get_degs_by_direction(deg2, "down")
down_genes3 <- get_degs_by_direction(deg3, "down")

# DIAGNOSTIC STEP
cat("Number of DOWN-regulated genes in 'LOY Male vs CN Male':", length(down_genes1), "\n")
cat("Number of DOWN-regulated genes in 'LOY Male vs CN Female':", length(down_genes2), "\n")
cat("Number of DOWN-regulated genes in 'CN Male vs CN Female':", length(down_genes3), "\n\n")

down_gene_list_for_venn <- list(
  `LOY Male vs CN Male` = down_genes1,
  `LOY Male vs CN Female` = down_genes2,
  `CN Male vs CN Female` = down_genes3
)

# Generate the Venn diagram object
venn_object_down <- venn.diagram(
  x = down_gene_list_for_venn,
  filename = NULL, # Set filename to NULL
  category.names = c("LOY Male vs CN Male", "LOY Male vs CN Female", "CN Male vs CN Female"),
  cat.cex = 1.1, cat.fontface = "bold", cat.col = color_palette,
  fill = color_palette, alpha = 0.5, cex = 1.2, fontface = "plain",
  margin = 0.2,
  cat.dist = c(0.06, 0.06, 0.06)
)

# Save to PDF
pdf(file.path(figure_dir, "deg_overlap_downregulated.pdf"), height = 7, width = 7)
grid.draw(venn_object_down)
dev.off()

cat("Script finished. Venn diagrams have been saved as PDF files.\n")
