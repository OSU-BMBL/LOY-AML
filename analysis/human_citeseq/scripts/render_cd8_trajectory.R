#!/usr/bin/env Rscript
# CD8 trajectory on the RPCA-integrated embedding with Monocle3.
# Run from the repository root after generating the cd8_states target.

source("config/paths.R")
paths <- loy_paths("single_cell")
args <- getOption("loy.script_args", commandArgs(trailingOnly = TRUE))
if (length(args) > 2L) {
  stop("Usage: Rscript scripts/run.R cd8-trajectory ",
       "[targets_store] [output_directory]")
}
targets_store <- if (length(args) >= 1L) args[[1L]] else
  file.path(paths$results, "_targets")
BASE <- if (length(args) >= 2L) args[[2L]] else
  file.path(paths$results, "figures")
if (!dir.exists(targets_store)) {
  stop("Target store not found: ", targets_store,
       ". Generate the cd8_states target before rendering the trajectory.")
}
source("analysis/human_citeseq/R/plot_style.R")
source("analysis/human_citeseq/R/plot_helpers.R")
suppressPackageStartupMessages({
  library(targets)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(SeuratWrappers)
  library(dplyr)
  library(igraph)
})
options(future.globals.maxSize = 4e9)
SEED <- 42
dir.create(file.path(BASE, "fig5"), showWarnings = FALSE, recursive = TRUE)

# RPCA integration ---------------------------------------------------------
cd8 <- tar_read(cd8_states, store = targets_store)$cd8_obj
cd8 <- subset(cd8, subset = group %in% c("LOY", "CN-Male"))
cd8$group <- factor(as.character(cd8$group), levels = c("LOY", "CN-Male"))
cd8$assigned_state <- factor(as.character(cd8$assigned_state),
                             levels = c("Naive", "Effector", "Dysfunctional"))
cd8$study <- ifelse(grepl("^Pool", cd8$orig.ident.x), "Velegraki", "GSE185381")
DefaultAssay(cd8) <- "RNA"

lst <- SplitObject(cd8, split.by = "study")
lst <- lapply(lst, function(o) FindVariableFeatures(NormalizeData(o, verbose = FALSE),
                                                    nfeatures = 2000, verbose = FALSE))
feats <- SelectIntegrationFeatures(lst, nfeatures = 2000)
lst <- lapply(lst, function(o) RunPCA(ScaleData(o, features = feats, verbose = FALSE),
                                      features = feats, npcs = 20, verbose = FALSE))
anchors <- FindIntegrationAnchors(lst, reduction = "rpca",
                                  anchor.features = feats, dims = 1:20)
intg <- IntegrateData(anchors, dims = 1:20)
DefaultAssay(intg) <- "integrated"
intg <- ScaleData(intg, verbose = FALSE)
intg <- RunPCA(intg, npcs = 20, verbose = FALSE)
set.seed(SEED)
# Construct the UMAP used for trajectory inference.
intg <- RunUMAP(intg, reduction = "pca", dims = 1:20, verbose = FALSE,
                n.neighbors = 10, min.dist = 0.05, seed.use = SEED)
intg$assigned_state <- factor(as.character(intg$assigned_state),
                              levels = c("Naive", "Effector", "Dysfunctional"))
intg$group <- factor(as.character(intg$group), levels = c("LOY", "CN-Male"))

# Monocle3 on the RPCA UMAP ------------------------------------------------
cds <- SeuratWrappers::as.cell_data_set(intg)
set.seed(SEED)
cds <- cluster_cells(cds, reduction_method = "UMAP", random_seed = SEED)
set.seed(SEED)
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE,
                   learn_graph_control = list(minimal_branch_len = 10))

# Root at the principal-graph vertex with the most Naive cells among
# vertices containing at least 20 cells.
cvm <- as.data.frame(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
colnames(cvm)[1] <- "vertex"
cvm$cell <- rownames(cvm)
cvm$state <- as.character(intg$assigned_state[cvm$cell])
vst <- cvm %>% group_by(vertex) %>%
  summarise(n = n(), n_naive = sum(state == "Naive", na.rm = TRUE),
            pct_naive = n_naive / n(), .groups = "drop") %>%
  arrange(desc(n_naive), desc(pct_naive))
root_vertex <- vst %>% filter(n >= 20) %>% slice(1) %>% pull(vertex)
root_node <- paste0("Y_", root_vertex)
cat("Root node (most-Naive vertex):", root_node, "\n")
cds <- order_cells(cds, root_pr_nodes = root_node)
intg$pseudotime <- pseudotime(cds)

# Pseudotime summaries -----------------------------------------------------
ps <- intg@meta.data %>% filter(is.finite(pseudotime)) %>%
  group_by(assigned_state) %>%
  summarise(median_pt = round(median(pseudotime), 2),
            mean_pt = round(mean(pseudotime), 2), n = n(), .groups = "drop")
cat("\nPseudotime by state\n")
print(ps)
gp <- intg@meta.data %>% filter(is.finite(pseudotime)) %>%
  group_by(group) %>% summarise(median_pt = round(median(pseudotime), 2), n = n(), .groups = "drop")
cat("\nPseudotime by group\n")
print(gp)

# Figure panels -----------------------------------------------------------
base_theme <- theme_nature_contract() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), axis.title = element_text(size = 5.5),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5.5), plot.title = element_blank())

p_pt <- plot_cells(cds, color_cells_by = "pseudotime",
                   label_cell_groups = FALSE, label_leaves = FALSE,
                   label_branch_points = FALSE, label_roots = TRUE,
                   cell_size = 0.45, trajectory_graph_segment_size = 0.6,
                   graph_label_size = 2) +
  labs(x = "UMAP 1", y = "UMAP 2") + base_theme

p_state <- plot_cells(cds, color_cells_by = "assigned_state",
                      label_cell_groups = FALSE, label_leaves = FALSE,
                      label_branch_points = FALSE, label_roots = TRUE,
                      cell_size = 0.45, trajectory_graph_segment_size = 0.6) +
  scale_color_manual(values = COLORS$cd8_state, name = "CD8 state") +
  labs(x = "UMAP 1", y = "UMAP 2") + base_theme

p_group <- plot_cells(cds, color_cells_by = "group",
                      label_cell_groups = FALSE, label_leaves = FALSE,
                      label_branch_points = FALSE, label_roots = TRUE,
                      cell_size = 0.45, trajectory_graph_segment_size = 0.6) +
  scale_color_manual(values = COLORS$group[c("LOY", "CN-Male")],
                     breaks = c("LOY", "CN-Male"),
                     labels = parse(text = GROUP_LABELS[c("LOY", "CN-Male")]),
                     name = "Group") +
  labs(x = "UMAP 1", y = "UMAP 2") + base_theme

save_pub_r(p_pt,    file.path(BASE, "fig5/F5A2_traj_pseudotime_rpca"), width_mm = 80, height_mm = 70)
save_pub_r(p_state, file.path(BASE, "fig5/F5A2_traj_state_rpca"),      width_mm = 80, height_mm = 70)
save_pub_r(p_group, file.path(BASE, "fig5/F5A2_traj_group_rpca"),      width_mm = 80, height_mm = 70)
ggsave(file.path(BASE, "fig5/F5A2_traj_rpca_preview.png"),
       p_state + p_group + p_pt, width = 15, height = 5.2, dpi = 170, bg = "white")

# Pseudotime density: LOY versus CN-Male ------------------------------------
grp_lvls <- c("CN-Male", "LOY")
dens_df <- intg@meta.data %>%
  filter(is.finite(pseudotime), group %in% grp_lvls)
dens_df$group <- factor(as.character(dens_df$group), levels = grp_lvls)
med_df <- dens_df %>% group_by(group) %>%
  summarise(median_pt = median(pseudotime), .groups = "drop")
# Wilcoxon rank-sum on pseudotime, M^LOY vs M^CN (shift in distribution)
wt   <- wilcox.test(pseudotime ~ group, data = dens_df)
plab <- if (wt$p.value < 2.2e-16) "Wilcoxon p < 2.2e-16" else
  sprintf("Wilcoxon p = %s", formatC(wt$p.value, format = "e", digits = 1))
p_dens <- ggplot(dens_df, aes(x = pseudotime, fill = group, color = group)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(data = med_df, aes(xintercept = median_pt, color = group),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  annotate("text", x = Inf, y = Inf, label = plab, hjust = 1.08, vjust = 1.6,
           size = 1.9, colour = "grey20") +
  scale_fill_manual(values = COLORS$group[grp_lvls], breaks = grp_lvls,
                    labels = parse(text = GROUP_LABELS[grp_lvls])) +
  scale_color_manual(values = COLORS$group[grp_lvls], breaks = grp_lvls,
                     labels = parse(text = GROUP_LABELS[grp_lvls])) +
  labs(x = "Pseudotime", y = "Density", title = NULL) +
  theme_nature_contract() +
  theme(legend.position = "top", legend.title = element_blank())
save_pub_r(p_dens, file.path(BASE, "fig5/F5B_pseudotime_density_rpca"),
           width_mm = 70, height_mm = 55)

cat("\nWrote fig5/F5A2_traj_{pseudotime,state}_rpca + F5B_pseudotime_density_rpca (+ preview)\n")
