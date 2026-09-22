source("analysis/human_bulk/config.R")
figure_dir <- bulk_output("figures", "y_gene_expression")

# Y-chromosome expression

# 1. Packages
message("--- Loading required libraries... ---")
need_pkgs <- c(
  "tidyverse","edgeR","ggpubr","ggsci","pheatmap",
  "RColorBrewer","biomaRt","broom"
)
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install the documented dependencies before running: ", paste(missing_pkgs, collapse = ", "))
}
library(tidyverse)
library(edgeR)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(broom)
message("Libraries loaded successfully.")

# 2. Parameters and inputs

# Input paths
path_to_counts <- bulk_files$counts
path_to_meta   <- bulk_files$metadata
path_to_ciber  <- bulk_files$fractions

# Sentinel markers (robust single-copy Y genes)
sentinel_y <- c("EIF1AY","DDX3Y","RPS4Y1","KDM5D","UTY","ZFY")

# HSC-like columns in CIBERSORTx
hsc_columns <- c("HSC","HSC.like","Prog.like")

# 3. Load data
message("--- Loading data... ---")
counts_df <- read.csv(path_to_counts, row.names = 1, check.names = FALSE)
meta_df   <- read.csv(path_to_meta, check.names = FALSE)
ciber_res <- read.csv(path_to_ciber, check.names = FALSE)

# 2a. Standardize sample_id in meta
if (!"sample_id" %in% colnames(meta_df)) {
  # Try a few common fallbacks; otherwise stop with a clear message
  guess_cols <- c("Sample","sample","SAMPLE","sample_name","SampleID","sampleID","Sample_ID")
  hit <- guess_cols[guess_cols %in% colnames(meta_df)]
  if (length(hit) >= 1) {
    meta_df <- meta_df %>% rename(sample_id = !!hit[1])
    message(sprintf("'sample_id' not found; using column '%s' as sample_id.", hit[1]))
  } else {
    stop("meta_df needs a 'sample_id' column (or one of: Sample/sample/sample_name/SampleID/sampleID/Sample_ID).")
  }
}
if (!"sample_id" %in% colnames(ciber_res)) {
  # Use the first column as the sample identifier.
  colnames(ciber_res)[1] <- "sample_id"
  message("Renamed first column of CIBERSORTx results to 'sample_id'.")
}

# 3. Define Groups & Filter Samples
# Expect meta_df to have 'Cyto.Group' and 'Sex'
if (!all(c("Cyto.Group","Sex") %in% colnames(meta_df))) {
  stop("meta_df must contain 'Cyto.Group' and 'Sex' columns.")
}
# 1) Fix empty/NA/duplicate column names
nms <- names(meta_df)
bad <- is.na(nms) | trimws(nms) == ""
if (any(bad)) {
  nms[bad] <- paste0("X", seq_len(sum(bad)))  # X1, X2, ...
}
names(meta_df) <- make.unique(nms)

# (Optional) normalize whitespace and dots in names
names(meta_df) <- gsub("\\s+", "_", names(meta_df))
names(meta_df) <- gsub("\\.+", ".", names(meta_df))

pick_first <- function(cands) { hit <- intersect(cands, names(meta_df)); if (length(hit)) hit[1] else NA_character_ }

cyto_col <- pick_first(c("Cyto.Group","Cyto_Group","CytoGroup","Cytogenetics","Cyto.Group.","Cyto Group"))
sex_col  <- pick_first(c("Sex","SEX","sex","Gender","GENDER"))

if (is.na(cyto_col) || is.na(sex_col)) {
  stop(sprintf("Couldn't find Cyto.Group/Sex. Columns available: %s", paste(names(meta_df), collapse=", ")))
}

meta_df <- meta_df %>%
  dplyr::rename(Cyto.Group = !!cyto_col, Sex = !!sex_col)

meta_df <- meta_df %>%
  dplyr::mutate(
    Group = dplyr::case_when(
      Cyto.Group == "CN_AML" & Sex == "Male"   ~ "CN Males",
      Cyto.Group == "LOY"    & Sex == "Male"   ~ "LOY Males",
      Cyto.Group == "CN_AML" & Sex == "Female" ~ "CN Females",
      TRUE                                     ~ "Other"
    )
  ) %>%
  dplyr::filter(Group %in% c("CN Males","LOY Males","CN Females"))
table(meta_df$Cyto.Group, meta_df$Sex)
# Keep only samples present in counts and meta
keep <- intersect(colnames(counts_df), meta_df$sample_id)
if (length(keep) == 0) stop("No overlapping samples between counts_df columns and meta_df$sample_id.")
meta_df <- meta_df %>% filter(sample_id %in% keep)

# Reorder counts to match meta_df
counts_df <- counts_df[, meta_df$sample_id, drop = FALSE]

# Align CIBERSORTx to same order (and filter)
ciber_res <- ciber_res %>%
  filter(sample_id %in% meta_df$sample_id) %>%
  arrange(factor(sample_id, levels = meta_df$sample_id))

stopifnot(identical(meta_df$sample_id, colnames(counts_df)))
stopifnot(identical(ciber_res$sample_id, meta_df$sample_id))

# Sanity: no negative counts
counts_df[counts_df < 0] <- 0

# Y-chromosome genes from Ensembl
message("Fetching Y-chromosome annotations from Ensembl.")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
y_genes_df <- getBM(
  attributes = c("hgnc_symbol","chromosome_name","gene_biotype",
                 "start_position","end_position"),
  filters = "chromosome_name", values = "Y", mart = ensembl
) %>%
  distinct()
y_genes_all <- sort(unique(y_genes_df$hgnc_symbol))

# 5. HSC Fractions from CIBERSORTx
hsc_in <- intersect(hsc_columns, colnames(ciber_res))
if (length(hsc_in) == 0) {
  stop(sprintf("None of the HSC columns found in CIBERSORTx file. Looked for: %s", paste(hsc_columns, collapse=", ")))
}
hsc_fractions <- ciber_res %>%
  dplyr::select(sample_id, all_of(hsc_in)) %>%
  mutate(HSC_Total_Fraction = rowSums(across(all_of(hsc_in)))) %>%
  dplyr::select(sample_id, HSC_Total_Fraction)

# 6. edgeR Normalization -> Linear CPM
message("--- Normalizing counts with TMM and computing linear CPM... ---")
dge <- DGEList(counts = counts_df)
dge <- calcNormFactors(dge, method = "TMM")
linear_cpm <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)

# 7. Prepare Y-Gene Matrix & Group Labels
y_genes_in_data <- intersect(y_genes_all, rownames(linear_cpm))
if (length(y_genes_in_data) == 0) stop("No Y genes found in the count matrix after filtering.")

# Merge metadata and HSC fractions for modeling/plotting
meta_use <- meta_df %>%
  left_join(hsc_fractions, by = "sample_id") %>%
  mutate(Group = factor(Group, levels = c("CN Males","LOY Males","CN Females")))

# 8. Primary Figure: Boxplot of Sentinel Y Genes (bulk CPM, log2)
message("--- Building boxplot for sentinel Y genes... ---")
sentinel_in_data <- intersect(sentinel_y, rownames(linear_cpm))
if (length(sentinel_in_data) < 2) {
  warning("Few sentinel Y genes found in data; boxplot will include whatever is present.")
}
plot_df <- as.data.frame(linear_cpm[sentinel_in_data, , drop = FALSE]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "CPM") %>%
  left_join(meta_use, by = "sample_id") %>%
  mutate(Log2CPM = log2(CPM + 1),
         Gene = factor(Gene, levels = sentinel_in_data))

# Quick pairwise comparison we care about (LOY vs CN male) per gene
comparisons <- list(c("CN Males","LOY Males"))

p_box <- ggplot(plot_df, aes(x = Gene, y = Log2CPM, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.9) +
  geom_point(aes(group = Group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
             alpha = 0.5, size = 1) +
  stat_compare_means(
    comparisons = comparisons, method = "wilcox.test",
    label = "p.signif", hide.ns = TRUE, size = 4,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = c("CN Males"="#2E86C1","LOY Males"="#D35400","CN Females"="#239B56")) +
  labs(
    title = "Sentinel Y-Chromosome Genes",
    subtitle = "log2(CPM+1); edgeR TMM-normalized",
    x = "Gene", y = "Expression (log2 CPM + 1)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(file.path(figure_dir, "sentinel_y_genes.pdf"), p_box, width = 16, height = 7)

# 9. Heatmap: All Y Genes, Group Means (average on linear scale, then log)
message("--- Building heatmap for all Y genes (group means on linear scale)... ---")
cn_male_samples <- meta_use %>% filter(Group == "CN Males")  %>% pull(sample_id)
loy_male_samples <- meta_use %>% filter(Group == "LOY Males") %>% pull(sample_id)
cn_female_samples <- meta_use %>% filter(Group == "CN Females") %>% pull(sample_id)

Y_lin <- linear_cpm[y_genes_in_data, , drop = FALSE]
cn_male_avg_lin <- if (length(cn_male_samples)) rowMeans(Y_lin[, cn_male_samples, drop=FALSE]) else rep(NA_real_, length(y_genes_in_data))
loy_male_avg_lin <- if (length(loy_male_samples)) rowMeans(Y_lin[, loy_male_samples, drop=FALSE]) else rep(NA_real_, length(y_genes_in_data))
cn_female_avg_lin <- if (length(cn_female_samples)) rowMeans(Y_lin[, cn_female_samples, drop=FALSE]) else rep(NA_real_, length(y_genes_in_data))

summary_heatmap_matrix <- cbind(
  "CN Males"  = log2(cn_male_avg_lin + 1),
  "LOY Males" = log2(loy_male_avg_lin + 1),
  "CN Females"= log2(cn_female_avg_lin + 1)
)
rownames(summary_heatmap_matrix) <- y_genes_in_data

# Order genes by CN male mean for readability
ord <- order(summary_heatmap_matrix[, "CN Males"], decreasing = TRUE)
summary_heatmap_matrix <- summary_heatmap_matrix[ord, , drop = FALSE]

pheatmap(
  summary_heatmap_matrix,
  main = "Y-Chromosome Genes (Group Means; log2(CPM+1))",
  fontsize_row = 6,
  fontsize_col = 11,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  scale = "row",
  filename = file.path(figure_dir, "y_genes_group_means.pdf"),
  width = 7, height = max(8, min(36, 0.22 * nrow(summary_heatmap_matrix)))
)

# Annotate MSY vs PAR for the heatmap rows

# We already have y_genes_df with coordinates from biomaRt
# Make a lookup for the genes actually in the heatmap
y_list <- rownames(summary_heatmap_matrix)

y_info <- y_genes_df %>%
  dplyr::filter(hgnc_symbol %in% y_list) %>%
  dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

# GRCh38 PAR coordinates on Y
# PAR1: 10,000–2,781,479;  PAR2: 56,887,902–57,217,415
is_PAR <- with(y_info,
               (start_position <= 2781479 & end_position >= 10000) |
                 (start_position <= 57217415 & end_position >= 56887902)
)

y_info$Type <- ifelse(is_PAR, "PAR", "MSY")

# Build an ordered vector aligned to heatmap rows
y_type <- setNames(y_info$Type, y_info$hgnc_symbol)
y_type <- y_type[y_list]
# Any genes not found in y_info (rare) → default to MSY
y_type[is.na(y_type)] <- "MSY"

ann_row <- data.frame(Type = factor(y_type, levels = c("MSY","PAR")))
rownames(ann_row) <- y_list

# Optional: define colors for the annotation
ann_colors <- list(Type = c(MSY = "#1f77b4", PAR = "#ff7f0e"))

pdf(file.path(figure_dir, "y_genes_annotated_heatmap.pdf"), width = 5, height = 8)

pheatmap(
  summary_heatmap_matrix,
  main = "Y-Chromosome Genes (Group Means; log2(CPM+1))",
  fontsize_row = 11,
  fontsize_col = 11,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  scale = "row",
  angle_col = "45",
  width = 7, height = max(8, min(36, 0.22 * nrow(summary_heatmap_matrix)))
)
dev.off()

library(writexl)
library(tibble)
save_df <- as.data.frame(summary_heatmap_matrix) %>% rownames_to_column("Gene")

write_xlsx(
  save_df,
  path = file.path(figure_dir, "y_genes_group_means.xlsx")
)

# (Optional) save the MSY/PAR table you used
write.csv(
  data.frame(Gene = y_list, Type = as.character(ann_row$Type)),
  file.path(figure_dir, "y_genes_msy_par_annotation.csv"),
  row.names = FALSE
)

# 8b. Boxplot for ALL Y genes used (BH-adjusted Wilcoxon per gene)

# Use the same ordered genes as the heatmap input, if available
genes_for_boxplot <- if (exists("summary_heatmap_matrix")) {
  rownames(summary_heatmap_matrix)
} else {
  y_genes_in_data
}

# Male groups only (same as the sentinel plot)
meta_males <- meta_use %>% dplyr::filter(Group %in% c("CN Males","LOY Males"))

plot_all_df <- as.data.frame(linear_cpm[genes_for_boxplot, meta_males$sample_id, drop = FALSE]) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "sample_id", values_to = "CPM") %>%
  dplyr::left_join(meta_males, by = "sample_id") %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = genes_for_boxplot),
    Group = factor(Group, levels = c("CN Males","LOY Males")),
    Log2CPM = log2(CPM + 1)
  )

# Per-gene Wilcoxon LOY vs CN and BH adjust
stats_all <- plot_all_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    p = tryCatch(
      wilcox.test(Log2CPM[Group=="CN Males"], Log2CPM[Group=="LOY Males"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH"),
    p_label = dplyr::case_when(
      is.na(p_adj)   ~ "NA",
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Y positions for labels (a bit above the per-gene max)
ypos_all <- plot_all_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(y = max(Log2CPM, na.rm = TRUE) + 0.12, .groups = "drop")

label_all <- dplyr::left_join(stats_all, ypos_all, by = "Gene")
plot_all_df$Type <- y_type[as.character(plot_all_df$Gene)]

p_box_all <- ggplot(plot_all_df, aes(x = Gene, y = Log2CPM, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.9) +
  geom_point(
    aes(color = Type),
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.25),
    alpha = 0.75, size = 0.7
  ) +
  geom_text(data = label_all, aes(Gene, y, label = p_label),
            inherit.aes = FALSE, vjust = 0, size = 3.1) +
  scale_fill_manual(values = c("CN Males"="#2E86C1","LOY Males"="#D35400")) +
  scale_color_manual(values = c(MSY = "#1f77b4", PAR = "#ff7f0e"), name = "Type") +
  labs(
    title = "All Y Genes • CN vs LOY (males)",
    subtitle = "Points colored by Type (MSY/PAR) • Wilcoxon BH • edgeR TMM → log2(CPM+1)",
    x = "Gene", y = "Expression (log2 CPM + 1)"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p_box_all)
ggsave(file.path(figure_dir, "y_genes_boxplots.pdf"), p_box_all,
       width = max(16, 0.35 * length(genes_for_boxplot)), height = 7, limitsize = FALSE)
