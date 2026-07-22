#loading the necessary packages: 
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(data.table)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
#loading expression matrix 
data<-read.delim("C:/Users/maram/Desktop/intital 1/Internships/Centre Biotechnologique de Sfax/codes/single cell/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",header=T)
dim(data)
head(data)
#load metadata--------------------
metadata<-read.csv("C:/Users/maram/Desktop/intital 1/Internships/Centre Biotechnologique de Sfax/codes/single cell/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.csv",header = T,sep =";" )
#convert rownames to sample name
rownames(metadata)<-metadata$Sample.name
rownames(metadata) <- gsub("-", ".", c(rownames(metadata)))
head(rownames(metadata))
all(colnames(data) %in% rownames(metadata))
head(colnames(data))
missing_cells <- setdiff(colnames(data), rownames(metadata))
head(missing_cells)
rownames(data) <- data$GENE
data$GENE <- NULL
meta_ss2 <- metadata %>%
  dplyr::filter(processed.data.file == "Smartseq2_GBM_IDHwt_processed_TPM.tsv") %>%
  dplyr::select(Sample.name, tumour.name, adult.pediatric) %>%
  as.data.frame()
# Keep metadata in the same order as the expression matrix columns
meta_ss2<-meta_ss2[colnames(data),]
#creating a seurat object:
object<-CreateSeuratObject(counts = data,project = "Smart-seq",min.cells = 3,min.features = 0)
object<-SetAssayData(object,layer = "data",new.data = data)
object <- AddMetaData(object, metadata = meta_ss2[colnames(object), c("tumour.name", "adult.pediatric")])
# QC ----------------------------
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ggsave("qc_violin_inspect_only.png", width = 8, height = 4)
#plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
ggsave("scatter_plot_of_QC.png", width = 8, height = 4)
summary(object$nFeature_RNA)
summary(object$nCount_RNA)
#quantile(object$nFeature_RNA,probs = c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))
#quantile(object$nCount_RNA,probs = c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))
object <- subset(
  object,
  subset = nFeature_RNA > 1000 &
    nFeature_RNA < 10000
)
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ggsave("Volcano_plot_after_filtration.png", width = 8, height = 4)
#Cell filtration:---------------------------------
#    Exact marker sets used in the paper's STAR Methods
# -------------------------------------------------------------
macrophage_markers     <- c("CD14","AIF1","FCER1G","FCGR3A","TYROBP","CSF1R")
tcell_markers           <- c("CD2","CD3D","CD3E","CD3G")
oligodendrocyte_markers <- c("MBP","TF","PLP1","MAG","MOG","CLDN11")

object <- AddModuleScore(object, features = list(macrophage_markers),
                      name = "MacrophageScore", ctrl = 50)
object <- AddModuleScore(object, features = list(tcell_markers),
                      name = "TcellScore", ctrl = 50)
object <- AddModuleScore(object, features = list(oligodendrocyte_markers),
                      name = "OligoScore", ctrl = 50)
score_cut <- function(x) mean(x) + 2 * sd(x)

object$is_macrophage <- object$MacrophageScore1 > score_cut(object$MacrophageScore1)
object$is_tcell       <- object$TcellScore1       > score_cut(object$TcellScore1)
object$is_oligo        <- object$OligoScore1       > score_cut(object$OligoScore1)

object$marker_nonmalignant <- object$is_macrophage | object$is_tcell | object$is_oligo
edb <- EnsDb.Hsapiens.v86

genes_in_object <- rownames(object)

gene_pos <- ensembldb::genes(
  edb,
  filter = ~ gene_biotype == "protein_coding",
  columns = c("gene_name")
) %>%
  as.data.frame() %>%
  dplyr::mutate(seq_name = as.character(seqnames)) %>%
  dplyr::filter(
    gene_name %in% genes_in_object,
    seq_name %in% as.character(c(1:22, "X"))
  ) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::arrange(seq_name, start)
avg_expr <- Matrix::rowMeans(GetAssayData(object, layer = "data"))
expressed_genes <- names(avg_expr[avg_expr > 0.1])   # tune to your data's scale
gene_pos <- gene_pos %>% dplyr::filter(gene_name %in% expressed_genes)

window_size <- 100
sliding_window_cna <- function(mat_sparse, gene_order, window = 100) {
  # mat_sparse: genes x cells sparse matrix, already centered per gene
  # gene_order: character vector of genes in chromosomal order
  mat_ord <- mat_sparse[gene_order, , drop = FALSE]
  n <- nrow(mat_ord)
  # rolling mean via cumsum on a sparse-friendly dense-per-chunk basis
  cs <- apply(as.matrix(mat_ord), 2, function(col) {
    cx <- cumsum(c(0, col))
    (cx[(window+1):(n+1)] - cx[1:(n-window+1)]) / window
  })
  cs
}
run_cna_for_tumor <- function(object, tumor_id, gene_pos, window = 100) {
  cells <- colnames(object)[object$orig.ident == tumor_id]
  sub <- subset(object, cells = cells)
  expr <- GetAssayData(sub, layer = "data")[gene_pos$gene_name, , drop = FALSE]
  gene_means <- Matrix::rowMeans(expr)
  expr_centered <- expr - gene_means
  
  cna_list <- list()
  for (chr in unique(gene_pos$seq_name)) {
    genes_chr <- gene_pos$gene_name[gene_pos$seq_name == chr]
    if (length(genes_chr) < window) next
    cna_list[[chr]] <- sliding_window_cna(expr_centered, genes_chr, window)
  }
  cna_all <- do.call(rbind, cna_list)
  
  cna_signal <- Matrix::colMeans(cna_all^2)
  candidate_malignant <- !sub$marker_nonmalignant
  tumor_avg_profile <- Matrix::rowMeans(cna_all[, candidate_malignant, drop = FALSE])
  cna_correlation <- apply(cna_all, 2, function(col) cor(col, tumor_avg_profile))
  
  data.frame(cell = colnames(sub), cna_signal, cna_correlation)
}
tumors <- unique(object$orig.ident)
cna_results <- lapply(tumors, function(t) {
  message("Processing tumor: ", t)
  res <- run_cna_for_tumor(object, t, gene_pos, window_size)
  gc()   # free memory between tumors -- important on 8GB
  res
}) %>% bind_rows()
threshold <- quantile(object$cna_correlation[object$marker_nonmalignant == TRUE], 0.95)
cna_results$cna_malignant <- cna_results$cna_signal > 0.02 &
  cna_results$cna_correlation > threshold
object$cna_signal      <- cna_results$cna_signal[match(colnames(object), cna_results$cell)]
object$cna_correlation <- cna_results$cna_correlation[match(colnames(object), cna_results$cell)]
object$cna_malignant    <- cna_results$cna_malignant[match(colnames(object), cna_results$cell)]
object$cna_signal <- cna_results$cna_signal[match(colnames(object), cna_results$cell)]
# Flag clusters that are >80% marker-positive non-malignant as non-malignant
cluster_nonmal_frac <- object@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(frac_nonmal = mean(marker_nonmalignant))

nonmal_clusters <- cluster_nonmal_frac$seurat_clusters[cluster_nonmal_frac$frac_nonmal > 0.8]
object$cluster_nonmalignant <- object$seurat_clusters %in% nonmal_clusters
# 4. INTEGRATED CALL (all three must agree, as in the paper)
# -------------------------------------------------------------
object$malignant <- object$cna_malignant &
  !object$marker_nonmalignant

object_malignant <- subset(object, subset = malignant == TRUE)

table(object$malignant)
saveRDS(object_malignant, "malignant_cells_only.rds")

#Variable features+scaling:-------------------------
object_malignant <- FindVariableFeatures(object_malignant,selection.method = "vst",nfeatures=2000)
top10<-head(VariableFeatures(object_malignant),10)
print(top10)
plot3<-VariableFeaturePlot(object_malignant)
plot4<-LabelPoints(plot=plot3,points=top10,repel = T,max.overlaps=20)
print(plot3)
print(plot4)
ggsave("variable_features_annotated.png", width = 8, height = 4)
#scaling:
allgenes<-rownames(object_malignant)
object_malignant <- ScaleData(object_malignant, features = allgenes)
zero_var_genes <- setdiff(rownames(object_malignant), rownames(LayerData(object_malignant, assay = "RNA", layer = "scale.data")))
object_malignant <- subset(object_malignant, features = setdiff(rownames(object_malignant), zero_var_genes))
#PCA+clustering+ UMAP+tSNE-----------
#PCA
object_malignant <- RunPCA(object_malignant, features = VariableFeatures(object = object_malignant))
ElbowPlot(object_malignant, ndims = 30)
ggsave("elbow_plot.png", width = 6, height = 4)
#DimPlot(object, reduction = "pca") + NoLegend()
n_pcs <- 10
object_malignant <- FindNeighbors(object_malignant, dims = 1:n_pcs)
object_malignant <- FindClusters(object_malignant, resolution = 0.5)
# PCA Scatter Plots
##############################

plots <- list(
  
  DimPlot(
    object_malignant,
    reduction="pca",
    group.by="seurat_clusters",
    label=TRUE
  ) + ggtitle("PCA - Clusters"),
  
  DimPlot(
    object_malignant,
    reduction="pca",
    group.by="tumour.name"
  ) + ggtitle("PCA - Tumor"),
  
  DimPlot(
    object_malignant,
    reduction="pca",
    group.by="adult.pediatric"
  ) + ggtitle("PCA - Age"),
  
  DimPlot(
    object_malignant,
    reduction="pca",
    dims=c(2,3),
    group.by="seurat_clusters",
    label=TRUE
  ) + ggtitle("PC2 vs PC3"),
  
  DimPlot(
    object_malignant,
    reduction="pca",
    dims=c(3,4),
    group.by="seurat_clusters",
    label=TRUE
  ) + ggtitle("PC3 vs PC4")
  
)

names(plots) <- c(
  "02_PCA_Clusters",
  "03_PCA_Tumor",
  "04_PCA_Age",
  "05_PC2_PC3",
  "06_PC3_PC4"
)

for(i in seq_along(plots)){
  ggsave(
    paste0("PCA_Results/",names(plots)[i],".png"),
    plots[[i]],
    width=7,
    height=6,
    dpi=300
  )
}

##############################
# Gene Loading Plots
##############################

pdf("PCA_Results/Loadings_PC1_PC5.pdf", width = 12, height = 10)
VizDimLoadings(object_malignant, dims = 1:5, reduction = "pca",nfeatures = 10)
dev.off()

pdf("PCA_Results/Loadings_PC6_PC10.pdf", width = 12, height = 10)
VizDimLoadings(object_malignant, dims = 6:10, reduction = "pca",nfeatures=10)
dev.off()
pdf("PCA_Results/Loadings_PC1_PC10.pdf", width = 12, height = 10)
VizDimLoadings(object_malignant, dims = 1:10, reduction = "pca",nfeatures=10)
dev.off()
# PCA Heatmap
##############################

pdf(
  "PCA_Results/08_PCA_Heatmap_PC1_PC10.pdf",
  width=12,
  height=12
)

DimHeatmap(
  object_malignant,
  dims=1:10,
  cells=500,
  balanced=TRUE
)

dev.off()
# Variance Explained
##############################

stdev <- object_malignant[["pca"]]@stdev

variance <- data.frame(
  
  PC = paste0("PC",1:length(stdev)),
  Standard_Deviation = stdev,
  Variance_Explained = (stdev^2)/sum(stdev^2)*100,
  Cumulative = cumsum((stdev^2)/sum(stdev^2)*100)
  
)

write.csv(
  variance,
  "PCA_Results/PCA_Variance_Explained.csv",
  row.names=FALSE
)
# Variance Plot
##############################

p <- ggplot(
  variance[1:30,],
  aes(
    x=1:30,
    y=Variance_Explained
  ))+
  
  geom_point(size=3)+
  geom_line()+
  theme_classic(base_size=14)+
  xlab("Principal Component")+
  ylab("% Variance Explained")

ggsave(
  "PCA_Results/09_Variance_Explained.png",
  p,
  width=7,
  height=5,
  dpi=300
)
# PCA Coordinates
##############################

coords <- Embeddings(
  object_malignant,
  reduction="pca"
)

write.csv(
  coords,
  "PCA_Results/PCA_Cell_Coordinates.csv"
)
# PCA Loadings
##############################

loadings <- Loadings(
  object_malignant[["pca"]]
)

write.csv(
  loadings,
  "PCA_Results/PCA_Gene_Loadings.csv"
)
# Top Genes per PC
##############################

top.genes <- data.frame()

for(i in 1:10){
  
  L <- loadings[,i]
  
  positive <- sort(L,decreasing=TRUE)[1:30]
  negative <- sort(L)[1:30]
  
  tmp <- data.frame(
    
    PC = paste0("PC",i),
    
    Positive_Gene = names(positive),
    Positive_Loading = positive,
    
    Negative_Gene = names(negative),
    Negative_Loading = negative
    
  )
  
  top.genes <- rbind(top.genes,tmp)
  
}

write.csv(
  top.genes,
  "PCA_Results/Top_Genes_Per_PC.csv",
  row.names=FALSE
)

#Cluster marker genes-----------------------------
all_markers <- FindAllMarkers(object_malignant, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#saving the markers
saveRDS(all_markers, "all_markers.rds")
all_markers <- readRDS("all_markers.rds")
top50 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50)
write.csv(top50, "top50_markers.csv", row.names = FALSE)
#cell annotation
annotation<- as.SingleCellExperiment(object_malignant)
ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(
  test = annotation,
  ref = ref,
  labels = ref$label.main
)
object_malignant$SingleR.labels <- pred$labels
object_malignant <- RunUMAP(object_malignant, dims = 1:n_pcs)
DimPlot(
  object_malignant,
  reduction = "umap",
  group.by = "SingleR.labels",
  label = TRUE,
  repel = TRUE,
)
ggsave("umap_by_cell_annotation.png", width = 6, height = 5)
cell_counts <- table(object_malignant$SingleR.labels)
print(cell_counts)
cell_percentages <- prop.table(cell_counts) * 100
print(round(cell_percentages, 2))
#save UMAP

DimPlot(object_malignant, reduction = "umap")
DimPlot(object_malignant, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")
ggsave("umap_by_clusters_annotation.png", width = 6, height = 5)
DimPlot(object_malignant, reduction = "umap", group.by = "tumour.name") +
  ggtitle("By tumor of origin")
ggsave("umap_by_tumour_annotation.png", width = 6, height = 5)
DimPlot(object_malignant, reduction = "umap", group.by = "adult.pediatric") +
  ggtitle("Adult vs Pediatric")
ggsave("umap_by_age_annotation.png", width = 6, height = 5)
#tSNE plot
object_malignant<-RunTSNE(object_malignant,dims = 1:n_pcs)
DimPlot(object_malignant, reduction="tsne",group.by="SingleR.labels")+
  ggtitle("Grouped by annotated cells")
ggsave("tsne_by_cell_clustering.png", width = 6, height = 5)
DimPlot(object_malignant,reduction = "tsne",group.by = "tumour.name")+
  ggtitle("Grouped by tumor")
ggsave("tsne_by_classified_by_tumour.png", width = 6, height = 5)
DimPlot(object_malignant,reduction = "tsne",group.by = "seurat_clusters",label = T)+
  ggtitle("Grouped by clusters")
ggsave("tsne_classified_by_clusters.png", width = 6, height = 5)
 
#differential expression analysis for age_group as a feature:
Idents(object_malignant)<-"adult.pediatric"
if (length(unique(object_malignant$adult.pediatric)) == 2) {
  de_age <- FindMarkers(object_malignant, ident.1 = "adult", ident.2 = "pediatric")
  write.csv(de_age, "adult_vs_pediatric_DE.csv")
  cat("\nTop adult-vs-pediatric DE genes:\n")
  print(head(de_age[order(de_age$p_val_adj), ], 20))
} else {
  cat("\nOnly one age group present in this subset - skipping adult vs pediatric DE.\n")
}
# Restore cluster identities for any further downstream work
Idents(object_malignant) <- "seurat_clusters"
## 10. Save processed object for downstream IDR-integration work --------------
saveRDS(object_malignant, file = "GSE131928_SS2_seurat_processed.rds")
object_malignant <- readRDS("GSE131928_SS2_seurat_processed.rds")
expr_normalized <- GetAssayData(object_malignant, layer = "data")
write.csv(as.matrix(expr_normalized), "GSE131928_SS2_expression_normalized.csv")
cell_metadata <- object_malignant@meta.data
write.csv(cell_metadata, "GSE131928_SS2_cell_metadata.csv")
#saving the code session:
save.image("GBM_workspace.RData")
