library(Seurat)
library(EnhancedVolcano)
# loading data:
de_age <- read.csv("DE version 1/adult_vs_pediatric_DE.csv", row.names = 1)
de_age$gene <- rownames(de_age)

# Detection gap: how "exclusive" is expression to one group
de_age$pct_diff <- de_age$pct.1 - de_age$pct.2   # positive = more exclusive to adult

# Specificity ratio (handles cases where pct.2 is near zero better than subtraction alone)
de_age$pct_ratio <- (de_age$pct.1 + 0.01) / (de_age$pct.2 + 0.01)
#Remove technical artifacts

artifact_pattern <- "^RP[SL][0-9]|^MT-|^MALAT1|^XIST|^HB[ABDEGMQZ][0-9]?$"
de_age$is_artifact <- grepl(artifact_pattern, de_age$gene)
de_clean <- de_age[!de_age$is_artifact, ]

# Genes specifically ON in adult, largely OFF in pediatric
adult_specific <- de_clean[
  de_clean$p_val_adj < 0.05 &
  de_clean$avg_log2FC >= 0.25 &          # meaningfully higher in adult
  de_clean$pct.1 >= 0.10 &              # broadly expressed in adult (majority of cells)
  de_clean$pct.2 >=0.15 &             # largely absent in pediatric
  de_clean$pct_diff >= 0.10,            # substantial detection gap
]
adult_specific <- adult_specific[order(-adult_specific$pct_diff), ]
write.csv(adult_specific,"filtred_genes_adult.csv")

# Genes specifically ON in pediatric, largely OFF in adult
pediatric_specific <- de_clean[
  de_clean$p_val_adj < 0.05 &
  de_clean$avg_log2FC <= -0.25 &
  de_clean$pct.2 >= 0.10 &
  de_clean$pct.1 <= 0.10 &
  de_clean$pct_diff <= -0.10,
]
pediatric_specific <- pediatric_specific[order(pediatric_specific$pct_diff), ]
write.csv(pediatric_specific,"filtred_genes_pediatric.csv")

cat("Adult-specific genes:", nrow(adult_specific), "\n")
cat("Pediatric-specific genes:", nrow(pediatric_specific), "\n")
#volcano plot
keyvals<-rep("grey75",nrow(de_clean))
names(keyvals)<-rep("Not significant",nrow(de_clean))
keyvals[de_clean$p_val_adj<0.05 &de_clean$avg_log2FC>=0.25]<-"firebrick"
names(keyvals)[de_clean$p_val_adj<0.05 &de_clean$avg_log2FC>=0.25]<-"Up in Adult"
keyvals[de_clean$p_val_adj<0.05 &de_clean$avg_log2FC<=-0.25]<-"steelblue"
names(keyvals)[de_clean$p_val_adj<0.05 &de_clean$avg_log2FC<=-0.25]<-"Up in Pediatric"

EnhancedVolcano(de_clean,
                lab = de_clean$gene,
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Adult vs Pediatric GBM",
                pointSize = 1.5,
                labSize = 3.0,
                colCustom = keyvals,
                colAlpha = 0.7,
                legendPosition = "right",
                legendLabSize = 10)
ggsave("volcano_enhanced_avec__filtrage_des_cellules.png", width = 8, height = 7, dpi = 300)
save.image("GBM_DE_analysis.RData")
