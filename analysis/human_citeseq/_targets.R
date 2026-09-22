# Human CITE-seq analysis and manuscript figure workflow.
# Run from the repository root using scripts/run.R.
library(targets)
source("config/paths.R")

tar_option_set(
  packages = c("Seurat", "qs", "dplyr", "ggplot2", "tidyr", "tibble",
               "DESeq2", "fgsea", "msigdbr", "presto", "Matrix"),
  format = "qs", memory = "transient", garbage_collection = TRUE,
  error = "continue"
)
tar_source("analysis/human_citeseq/R/")

list(
  # Input and annotation ----------------------------------------------------
  tar_target(seurat_input_file, file.path(loy_paths("single_cell")$results, "objects", "cohort.qs"), format = "file"),
  tar_target(seurat_filtered, load_filtered_seurat(seurat_input_file)),
  tar_target(nk_annotated, annotate_nk_lymphoid(seurat_filtered)),
  tar_target(lineage_annotation, annotate_lineage(seurat_filtered, nk_annotated)),
  tar_target(cd8_states, classify_cd8_3state(seurat_filtered)),
  tar_target(cd8_states_clean, classify_cd8_3state_clean(seurat_filtered, lineage_annotation)),

  # Differential expression and enrichment ---------------------------------
  tar_target(degs_by_celltype, run_degs_by_celltype(seurat_filtered)),
  tar_target(degs_by_cd8_state_clean, run_degs_by_cd8_state(cd8_states_clean)),
  tar_target(gene_dosage_results, run_gene_dosage_analysis(seurat_filtered)),
  tar_target(deg_pathway_enrichment, run_deg_pathway_enrichment(degs_by_celltype)),
  tar_target(deg_pathway_enrichment_by_cd8_state_clean,
             run_deg_pathway_enrichment_by_cd8_state(degs_by_cd8_state_clean)),
  # Blast ORA uses the cell-level Wilcoxon table.
  tar_target(blast_ora_f1d,
    run_blast_ora(degs_by_celltype$blast$wilcox,
      padj_cut = 0.05, lfc_cut = 0.5, simplify_cutoff = 0.7),
    packages = c(tar_option_get("packages"), "clusterProfiler", "org.Hs.eg.db")),
  tar_target(blast_ora_f1e,
    run_blast_ora(degs_by_celltype$blast$wilcox,
      padj_cut = 0.05, lfc_cut = 0.5, pvalueCutoff = 1, qvalueCutoff = 1,
      simplify_cutoff = NULL, minGSSize = 3, maxGSSize = 2000),
    packages = c(tar_option_get("packages"), "clusterProfiler", "org.Hs.eg.db")),

  # Module scores and external regulon results -------------------------------
  tar_target(blast_cytokines, check_blast_cytokines(degs_by_celltype, seurat_filtered)),
  tar_target(blast_cue_scores,
    compute_blast_cue_score(seurat_filtered, cd8_states, blast_cytokines)),
  tar_target(scenic_cd8_by_state_dir, file.path(loy_paths("single_cell")$data, "regulons", "cd8_states"), format = "file"),
  tar_target(scenic_cd8_by_state_imported,
    import_scenic_cd8_by_state_results(scenic_cd8_by_state_dir)),
  tar_target(scenic_blast_diff_file,
    file.path(loy_paths("single_cell")$data, "regulons", "blasts",
              "Leukemia_LOY_vs_CN_diff_regulons.csv"), format = "file"),

  # Figure outputs ---------------------------------------------------------
  tar_target(fig_final_F1A_umap_celltype_file,
    fig_final_F1A_umap_celltype(seurat_filtered), format = "file"),
  tar_target(fig_final_F1A_umap_group_file,
    fig_final_F1A_umap_group(seurat_filtered), format = "file"),
  tar_target(fig_final_F1B_blast_volcano_file,
    fig_final_F1B_blast_volcano(degs_by_celltype), format = "file"),
  tar_target(fig_final_F1C_genepair_leukemia_file,
    fig_final_F1C_genepair_leukemia(gene_dosage_results, degs_by_celltype), format = "file"),
  tar_target(fig_final_F1D_pathway_butterfly_file,
    fig_final_F1D_pathway_butterfly(blast_ora_f1d), format = "file"),
  tar_target(fig_final_F1E_pathway_origterms_file,
    fig_final_F1E_pathway_origterms(blast_ora_f1e), format = "file"),

  # Y-linked expression and paired-gene dosage.
  tar_target(fig_final_F2cd_paired_YX_file,
    fig_final_F2cd_paired_YX(gene_dosage_results, degs_by_celltype), format = "file"),
  tar_target(fig_final_F2CD_dosage_compensation_file,
    fig_final_F2CD_dosage_compensation(gene_dosage_results), format = "file"),
  tar_target(fig_final_F3a_ysig_umap_file,
    fig_final_F3a_ysig_umap(seurat_filtered, lineage_annotation), format = "file"),
  tar_target(fig_final_F3b_ygene_dotplot_file,
    fig_final_F3b_ygene_dotplot(seurat_filtered, lineage_annotation), format = "file"),
  tar_target(fig_final_F3d_ysig_persample_violin_file,
    fig_final_F3d_ysig_persample_violin(seurat_filtered, lineage_annotation), format = "file"),

  # Blast pathways and antigen-presentation scores.
  tar_target(fig_final_F4B_blast_gsea_hallmark_file,
    fig_final_F4B_blast_gsea_hallmark(deg_pathway_enrichment), format = "file"),
  tar_target(fig_final_F4D_blast_mhc_ap_file,
    fig_final_F4D_blast_mhc_ap(seurat_filtered, degs_by_celltype), format = "file"),
  tar_target(fig_final_F4E_nlrc5_ciita_file,
    fig_final_F4E_nlrc5_ciita(seurat_filtered, degs_by_celltype), format = "file"),
  tar_target(fig_final_F4F_blast_apm_score_file,
    fig_final_F4F_blast_apm_score(blast_cue_scores), format = "file"),
  tar_target(fig_final_F4opt_scenic_blast_tf_file,
    fig_final_F4opt_scenic_blast_tf(scenic_blast_diff_file), format = "file"),

  # Clean CD8 states, differential expression and BACH2.
  tar_target(fig_final_F5C1_state_proportions_file,
    fig_final_F5C1_state_proportions(cd8_states_clean), format = "file"),
  tar_target(fig_final_F5b_naive_volcano_file,
    fig_final_F5b_naive_volcano(degs_by_cd8_state_clean), format = "file"),
  tar_target(fig_final_F5g_naive_pathway_file,
    fig_final_F5g_naive_pathway(deg_pathway_enrichment_by_cd8_state_clean), format = "file"),
  tar_target(fig_final_F5g2_naive_pathway_gobp_file,
    fig_final_F5g2_naive_pathway_gobp(deg_pathway_enrichment_by_cd8_state_clean), format = "file"),
  tar_target(fig_final_F5d_naive_scenic_file,
    fig_final_F5d_naive_scenic(scenic_cd8_by_state_imported), format = "file"),
  tar_target(fig_final_F5E1_bach2_violin_celltype_file,
    fig_final_F5E1_bach2_violin_celltype(seurat_filtered, lineage_annotation), format = "file"),
  tar_target(fig_final_F5E2_bach2_violin_state_file,
    fig_final_F5E2_bach2_violin_state(cd8_states_clean, degs_by_cd8_state_clean), format = "file"),
  tar_target(fig_final_F5E3_bach2_regulon_file,
    fig_final_F5E3_bach2_regulon(scenic_cd8_by_state_imported, cd8_states_clean), format = "file")
)
