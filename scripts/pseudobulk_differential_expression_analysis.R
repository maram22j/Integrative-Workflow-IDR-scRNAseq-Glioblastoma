# --- 0. Packages ---------------------------------------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(limma)
library(ggplot2)
library(ggrepel)
library(reshape2)

# ==========================================================================
# 1. Data preparation
# ==========================================================================

seurat_obj <- readRDS("object_malignant.rds")  

stopifnot(all(c("tumour.name", "adult.pediatric") %in% colnames(seurat_obj@meta.data)))

# Vérification rapide : combien de patients par groupe ?
table(unique(seurat_obj@meta.data[, c("tumour.name", "adult.pediatric")])$age_group)

# ==========================================================================
# 2.Expression matrix expression
# ==========================================================================
expr_mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
meta <- seurat_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, tumour.name, adult.pediatric)

# ==========================================================================
# 3. Aggregation by patients
# ==========================================================================

patients <- unique(meta$tumour.name)

pseudobulk_mat <- sapply(patients, function(p) {
  cells_p <- meta$cell_id[meta$tumour.name == p]
  Matrix::rowMeans(expr_mat[, cells_p, drop = FALSE])
})
colnames(pseudobulk_mat) <- patients

# Desing table by age_group
design_meta <- meta %>%
  distinct(tumour.name, adult.pediatric) %>%
  column_to_rownames("tumour.name")
design_meta <- design_meta[colnames(pseudobulk_mat), , drop = FALSE]

cat("Nombre de patients adultes :", sum(design_meta$adult.pediatric == "adult"), "\n")
cat("Nombre de patients pédiatriques :", sum(design_meta$adult.pediatric == "pediatric"), "\n")

# --- Filtering of lowly/not expressed genes (before the test, not after) ------
# We keep genes that are expressed (log-TPM > 0) in at least 20% of patients
keep_genes <- rowSums(pseudobulk_mat > 0) >= ceiling(0.2 * ncol(pseudobulk_mat))
pseudobulk_mat_filt <- pseudobulk_mat[keep_genes, ]
cat("Gènes conservés après filtrage :", nrow(pseudobulk_mat_filt),
    "/", nrow(pseudobulk_mat), "\n")

# ==========================================================================
# 4. DE analysis using Limma
# ==========================================================================

group <- factor(design_meta$adult.pediatric, levels = c("pediatric", "adult"))
design <- model.matrix(~ group)

fit <- lmFit(pseudobulk_mat_filt, design)
fit <- eBayes(fit, trend = TRUE,robust = TRUE)  

res_limma <- topTable(fit, coef = "groupadult", number = Inf) %>%
  rownames_to_column("gene") %>%
  rename(logFC = logFC, pval = P.Value, padj = adj.P.Val)

# ==========================================================================
# DE test option B (alternative) : Wilcoxon patient-level
# ==========================================================================

run_wilcoxon <- function(mat, group) {
  pvals <- apply(mat, 1, function(x) {
    tryCatch(wilcox.test(x[group == "adult"], x[group == "pediatric"])$p.value,
             error = function(e) NA)
  })
  logfc <- apply(mat, 1, function(x) {
    mean(x[group == "adult"]) - mean(x[group == "pediatric"])
  })
  data.frame(gene = rownames(mat), logFC = logfc, pval = pvals) %>%
    mutate(padj = p.adjust(pval, method = "BH"))
}

res_wilcox <- run_wilcoxon(pseudobulk_mat_filt, design_meta$adult.pediatric)

# ==========================================================================
# 5.DEG filtration (padj < 0.05, |logFC| >= 0.5)
# ==========================================================================
PADJ_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

degs_limma <- res_limma %>%
  filter(padj < PADJ_CUTOFF, abs(logFC) >= LOGFC_CUTOFF) %>%
  arrange(padj)

degs_wilcox <- res_wilcox %>%
  filter(padj < PADJ_CUTOFF, abs(logFC) >= LOGFC_CUTOFF) %>%
  arrange(padj)

cat("DEGs (limma)   :", nrow(degs_limma), "\n")
cat("DEGs (wilcoxon):", nrow(degs_wilcox), "\n")

overlap <- intersect(degs_limma$gene, degs_wilcox$gene)
cat("Recouvrement limma / wilcoxon :", length(overlap),
    "gènes (", round(100 * length(overlap) / nrow(degs_limma), 1), "% des DEGs limma)\n")

# ==========================================================================
# 6. VOLCANO PLOT
# ==========================================================================
res_limma <- res_limma %>%
  mutate(sig = case_when(
    padj < PADJ_CUTOFF & logFC >=  LOGFC_CUTOFF ~ "Up (adult)",
    padj < PADJ_CUTOFF & logFC <= -LOGFC_CUTOFF ~ "Up (pediatric)",
    TRUE ~ "NS"
  ))

volcano <- ggplot(res_limma, aes(x = logFC, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Up (adult)" = "firebrick",
                                 "Up (pediatric)" = "steelblue",
                                 "NS" = "grey70")) +
  geom_text_repel(data = degs_limma %>% slice_min(padj, n = 15),
                   aes(label = gene), size = 3, color = "black", max.overlaps = 20) +
  labs(title = "Volcano plot — pseudobulk par patient (limma)",
       subtitle = paste0("padj < ", PADJ_CUTOFF, ", |log2FC| >= ", LOGFC_CUTOFF),
       x = "log2 Fold Change (adulte vs pédiatrique)",
       y = "-log10(padj)") +
  theme_minimal()

print(volcano)
ggsave("DE analysis/volcano_plot_DE_with_logFC_1.png", width = 8, height = 4)

# ==========================================================================
# 7. Results
# ==========================================================================
write.csv(res_limma, "DE analysis/DE_pseudobulk_patient_limma_full.csv", row.names = FALSE)
write.csv(degs_limma, "DE analysis/DE_pseudobulk_patient_limma_DEGs.csv", row.names = FALSE)
write.csv(res_wilcox, "DE analysis/DE_pseudobulk_patient_wilcoxon_full.csv", row.names = FALSE)
write.csv(degs_wilcox, "DE analysis/DE_pseudobulk_patient_wilcoxon_DEGs.csv", row.names = FALSE)


cat("\nPipeline terminé. Fichiers exportés :\n",
    "- DE_pseudobulk_patient_limma_full.csv\n",
    "- DE_pseudobulk_patient_limma_DEGs.csv\n",
    "- DE_pseudobulk_patient_wilcoxon_full.csv\n",
    "- DE_pseudobulk_patient_wilcoxon_DEGs.csv\n",
    "- volcano_pseudobulk_patient.png\n")
# ==========================================================================
# 8. PCA 
# ==========================================================================
deg_genes <- degs_limma$gene
pca_mat <- pseudobulk_mat_filt[deg_genes, ]
pca_input <- t(pca_mat)
# PCA
pca_res <- prcomp(pca_input,
                  center = TRUE,
                  scale. = TRUE)
# Dataframe PCA
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  patient = rownames(pca_res$x),
  group = design_meta$adult.pediatric
)
# Variance expliquée
percent_var <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pca_plot <- ggplot(pca_df,
                   aes(x = PC1,
                       y = PC2,
                       color = group,
                       label = patient)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3,
                  max.overlaps = 30) +
  labs(
    title = "PCA based on limma Differentially Expressed Genes",
    x = paste0("PC1 (", round(percent_var[1],1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2],1), "% variance)")
  ) +
  theme_classic()

print(pca_plot)
ggsave("DE analysis/PCA_plot_of_DE_using_limma_with_logFC_1.png", width = 8, height = 4)
# ==========================================================================
# Boxplots and Density plots
# ==========================================================================
expr_deg <- pseudobulk_mat_filt[deg_genes, ]

# Format long
expr_long <- melt(as.data.frame(expr_deg),
                  variable.name = "Patient",
                  value.name = "Expression")

expr_long$Group <- design_meta$adult.pediatric[
  match(expr_long$Patient, rownames(design_meta))
]
# Boxplot
# ==========================================================================

boxplot_deg <- ggplot(expr_long,
                      aes(x = Patient,
                          y = Expression,
                          fill = Group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Expression distribution of limma DEGs",
    x = "Patients",
    y = "Log-normalized expression"
  )

print(boxplot_deg)

ggsave("DE analysis/Boxplot_limma_DEGs_with_logFC_1.png",boxplot_deg,width = 10,height = 5)

# Density plot
# ==========================================================================

density_deg <- ggplot(expr_long,
                      aes(x = Expression,
                          colour = Patient,
                          group = Patient)) +
  geom_density(linewidth = 0.7) +
  theme_classic() +
  labs(
    title = "Density plot of limma DEGs",
    x = "Log-normalized expression",
    y = "Density"
  ) +
  theme(legend.position = "none")

print(density_deg)

ggsave("DE analysis/Density_limma_DEGs_with_logFC_1.png",density_deg,width = 8,height = 5)

# ==========================================================================
# 9. HEATMAP — Top DEGs
# ==========================================================================

library(pheatmap)

# Select top 50 DEGs based on adjusted p-value
top_genes <- degs_limma %>%
  arrange(padj) %>%
  slice_head(n = 70) %>%
  pull(gene)

# Expression matrix
heatmap_mat <- pseudobulk_mat_filt[top_genes, , drop = FALSE]

# --------------------------------------------------------------------------
# Row-wise Z-score normalization
# --------------------------------------------------------------------------

heatmap_mat_z <- t(scale(t(heatmap_mat)))

# Replace NA values if a gene has zero variance
heatmap_mat_z[is.na(heatmap_mat_z)] <- 0

# --------------------------------------------------------------------------
# Patient annotation
# --------------------------------------------------------------------------

annotation_col <- data.frame(
  Group = design_meta$adult.pediatric
)

rownames(annotation_col) <- rownames(design_meta)

# Make sure annotation order matches heatmap columns
annotation_col <- annotation_col[colnames(heatmap_mat_z), , drop = FALSE]

# --------------------------------------------------------------------------
# Heatmap
# --------------------------------------------------------------------------

pheatmap(
  heatmap_mat_z,
  scale = "none",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  main = "Differentially Expressed Genes"
)

# --------------------------------------------------------------------------
# Save heatmap
# --------------------------------------------------------------------------

png(
  "DE analysis/Heatmap_limma_DEGs.png",
  width = 220,
  height = 180,
  res = 25
)

pheatmap(
  heatmap_mat_z,
  scale = "none",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  main = "Differentially Expressed Genes"
)

dev.off()
save.image("DE analysis/DE_analysis.RData")
