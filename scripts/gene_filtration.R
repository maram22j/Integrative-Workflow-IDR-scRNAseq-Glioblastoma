library(dplyr)
library(EnhancedVolcano)
data<-read.csv("DE analysis/DE_pseudobulk_patient_limma_DEGs.csv")
PADJ_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5
DE_adult<-data%>%
  filter(padj < PADJ_CUTOFF & logFC >=  LOGFC_CUTOFF)
DE_pediatric<-data%>%
  filter(padj < PADJ_CUTOFF & logFC <=  -LOGFC_CUTOFF)
write.csv(DE_adult,"gene filtration/DE_adult.csv",row.names = T)
write.csv(DE_pediatric,"gene filtration/DE_pediatric.csv",row.names = T)
save.image("gene filtration/Gene_filtration.RData")
