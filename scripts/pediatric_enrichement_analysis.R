library(dplyr)
library(ggplot2)
library(HGNChelper)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyr)
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)
library(ReactomePA)
library(ggtree)
library(pathview)
library(httr)
library(STRINGdb)

#Read pediatric expression data

pediatric_expression_data <- read.csv("gene filtration/DE_pediatric.csv",header = TRUE,sep = ",")

colnames(pediatric_expression_data)[colnames(pediatric_expression_data) == "gene"] <- "Gene"

pediatric_expression_data$X <- NULL

#Read pediatric IDR data

IDR_data_pediatric <- read.csv("IDR conversion/adult/AIUpred_pediatric_results.csv", sep = ",", header = TRUE)

colnames(IDR_data_pediatric)[colnames(IDR_data_pediatric) == "gene"] <- "Gene"

#Find genes missing from IDR dataset

missing_IDR_genes_pediatric <- pediatric_expression_data %>%
  anti_join(IDR_data_pediatric, by = "Gene") %>%
  pull(Gene)
# Check possible gene-symbol corrections
check_pediatric <- checkGeneSymbols(missing_IDR_genes_pediatric)

gene_corrections_pediatric <- check_pediatric %>%
  filter(!is.na(Suggested.Symbol)) %>%
  dplyr::select(x, Suggested.Symbol)
# Correct symbols in expression dataset

pediatric_expression_data <- pediatric_expression_data %>%
  left_join(
    gene_corrections_pediatric,
    by = c("Gene" = "x")
  ) %>%
  mutate(
    Gene = ifelse(
      !is.na(Suggested.Symbol),
      Suggested.Symbol,
      Gene
    )
  ) %>%
  dplyr::select(-Suggested.Symbol)

# 6. FINAL JOIN

IDR_data_pediatric <- pediatric_expression_data %>%
  left_join(
    IDR_data_pediatric,
    by = "Gene"
  )

# Check remaining missing IDR values
sum(is.na(IDR_data_pediatric$disorder_percentage)) 
# Gene ontology -----------------------------------------------------------------------------
genes_IDR_pediatric <- IDR_data_pediatric %>%
  pull(Gene)
gene_entrez_pediatric <- bitr(
  genes_IDR_pediatric,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Biological process
ego_BP_pediatric<- enrichGO(
  gene = gene_entrez_pediatric$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.35,
  qvalueCutoff = 0.35,
  readable = TRUE
)
BP<-as.data.frame(ego_BP_pediatric)
#Molecular Function
ego_MF_pediatric<- enrichGO(
  gene = gene_entrez_pediatric$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.3,
  qvalueCutoff = 0.3,
  readable = TRUE
)
MF<-as.data.frame(ego_MF_pediatric)
#Cellular component
ego_CC_pediatric<-enrichGO(
  gene = gene_entrez_pediatric$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.3,
  qvalueCutoff = 0.3,
  readable = TRUE
)
CC<-as.data.frame(ego_CC_pediatric)
extract_GO_by_gene <- function(ego){
  
  if (nrow(as.data.frame(ego)) == 0) {
    return(data.frame(
      Gene = character(),
      GO_terms = character()
    ))
  }
  
  as.data.frame(ego) %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(
      Gene = geneID,
      GO_terms = Description
    ) %>%
    group_by(Gene) %>%
    summarise(
      GO_terms = paste(unique(GO_terms), collapse = "; "),
      .groups = "drop"
    )
}
BP_genes <- extract_GO_by_gene(ego_BP_pediatric)
colnames(BP_genes)[2] <- "BP"

MF_genes <- extract_GO_by_gene(ego_MF_pediatric)
colnames(MF_genes)[2] <- "MF"

CC_genes <- extract_GO_by_gene(ego_CC_pediatric)
colnames(CC_genes)[2] <- "CC"
# Joining GO with the expression dataset
 IDR_pediatric_GO<- IDR_data_pediatric %>%
  left_join(BP_genes, by = "Gene") %>%
  left_join(MF_genes, by = "Gene") %>%
  left_join(CC_genes, by = "Gene")
colnames(IDR_pediatric_GO)[colnames(IDR_pediatric_GO)=="BP"]<-"Biological Process"
colnames(IDR_pediatric_GO)[colnames(IDR_pediatric_GO)=="MF"]<-"Molecular Function"
colnames(IDR_pediatric_GO)[colnames(IDR_pediatric_GO)=="CC"]<-"Cellular Component"
write.csv(IDR_pediatric_GO,file = "enrichement analysis/pediatric/pediatric_IDR_GO.csv",row.names = F)

#plots --------------
# Biological process
dotplot(
  ego_BP_pediatric,
  showCategory = 18,
  title= "Gene Ontology of Pediatric for Biological Process",
  orderBy = "pvalue",
  font.size = 10
)+
  theme_bw() +
  theme(
  axis.text.y = element_text(
    size = 9
  ),
  axis.text.x = element_text(
  size = 9
))
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_biological_process.png",width = 10, height = 10)
barplot(
  ego_BP_pediatric,
  showCategory = 20,
  title = "Gene Ontology for Pediatric for Biological Process",
  font.size=10
)+
  theme_bw() +
  theme(
    axis.text.y = element_text(
      size = 9
    ),
    axis.text.x = element_text(
      size = 9
    ))
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_barplot_biological_process.png",width = 10, height = 10)
#clustering plot
ego_BP_simplified <- clusterProfiler::simplify(
  ego_BP_pediatric,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_BP_pediatric_top <- pairwise_termsim(
  ego_BP_simplified)
p <- treeplot(
  ego_BP_pediatric_top,
  cluster_method = "average",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize_cladelab = 4
) +
  hexpand(0.3)

print(p)
ggsave(filename = "enrichement analysis/pediatric/treeplot_BP_pediatric.png",plot = p,width = 15,height = 15,dpi = 300)
#Cellular Component
dotplot(
  ego_CC_pediatric,
  showCategory = 18,
  title="Gene Ontology for Pediatric for Cellular Component",
  font.size=10
)+
  theme_bw() +
  theme(
    axis.text.y = element_text(
      size = 9
    ),
    axis.text.x = element_text(
      size = 9
    ))
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_cellular_component.png",width = 10, height = 10)
barplot(
  ego_CC_pediatric,
  showCategory = 20,
  title = "Gene Ontology for Pediatric for Cellular Component"
)
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_barplot_cellular_component.png",width = 10, height = 10)
#clustering plot
ego_CC_simplified <- clusterProfiler::simplify(
  ego_CC_pediatric,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_CC_pediatric_top <- pairwise_termsim(
  ego_CC_simplified)
c <- treeplot(
  ego_CC_pediatric_top,
  cluster_method = "average",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize_cladelab = 4,
) +
  hexpand(0.3)
print(p)
ggsave(filename = "enrichement analysis/pediatric/treeplot_CC_pediatric.png",plot = p,width = 15,height = 15,dpi = 300)
# Molecular Function
dotplot(
  ego_MF_pediatric,
  showCategory = 18,
  title="Gene Ontology for Pediatric for Molecular Function",
  orderBy="pvalue"
)
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_molecular_function.png",width = 10, height = 10)
barplot(
  ego_MF_pediatric,
  showCategory = 20,
  title = "Gene Ontology for pediatric for Molecular Function"
)
ggsave("enrichement analysis/pediatric/GO_enrichement_pediatric_barplot_Molecular_function.png",width = 10, height = 10)
#clustering plot
ego_MF_simplified <- clusterProfiler::simplify(
  ego_MF_pediatric,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_MF_pediatric_top <- pairwise_termsim(
  ego_MF_simplified)
M <- treeplot(
  ego_MF_pediatric_top,
  cluster_method = "average",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize_cladelab = 4,
) +
  hexpand(0.3)
print(M)
ggsave(filename = "enrichement analysis/pediatric/treeplot_MF_pediatric.png",plot = p,width = 15,height = 15,dpi = 300)

#KEGG enrichement----------------------------------------------------------------------------------------------
ekegg_pediatric <- enrichKEGG(
  gene = gene_entrez_pediatric$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.3,
  qvalueCutoff = 0.3
)
kk<-as.data.frame(ekegg_pediatric)
# Convert Entrez IDs to gene symbols
ekegg_pediatric_readable <- setReadable(
  ekegg_pediatric,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)
# Extract KEGG pathways associated with each gene
KEGG_genes <- as.data.frame(ekegg_pediatric_readable) %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(
    Gene = geneID,
    KEGG = Description
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    KEGG = paste(unique(KEGG), collapse = "; "),
    .groups = "drop"
  )
IDR_pediatric_GO <- IDR_pediatric_GO %>%
  left_join(KEGG_genes, by = "Gene")
write.csv(IDR_pediatric_GO,file = "enrichement analysis/pediatric/Pediatric_IDR_KEGG.csv",row.names = F)
#plotting
dotplot(
  ekegg_pediatric,
  showCategory = 18,
  title="KEGG for Pediatric",
  orderBy="pvalue"
)
ggsave("enrichement analysis/pediatric/KEGG_pediatric.png",width = 10, height = 10)
barplot(
  ekegg_pediatric,
  showCategory = 20,
    title = "KEGG for pediatric"
)
ggsave("enrichement analysis/pediatric/KEGG_barplot.png",width = 10, height = 10)
# Reactome enrichement analysis-------------------------------------------------------------------------------
reactome_pediatric <- enrichPathway(
  gene = gene_entrez_pediatric$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.3,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.3,
  readable = TRUE
)
r<-as.data.frame(reactome_pediatric)
Reactome_genes <- as.data.frame(reactome_pediatric) %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(
    Gene = geneID,
    Reactome = Description
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Reactome = paste(unique(Reactome), collapse = "; "),
    .groups = "drop"
  )
IDR_pediatric_GO <- IDR_pediatric_GO %>%
  left_join(Reactome_genes, by = "Gene")
write.csv(IDR_pediatric_GO, "enrichement analysis/pediatric/pediatric_IDR_GO_Reactome.csv", row.names = FALSE)
#plots
dotplot(
  reactome_pediatric,
  showCategory = 20,
  title = "Reactome Pathway Enrichment - pediatric"
)
ggsave("enrichement analysis/pediatric/Reactome_enrichment_pediatric_dotplot.png", width = 10, height = 10)
barplot(
  reactome_pediatric,
  showCategory = 20,
  title = "Reactome Pathway Enrichment - Pediatric"
)
ggsave("enrichement analysis/pediatric/Reactome_enrichment_pediatric_barplot.png", width = 10, height = 10)
save.image("enrichement analysis/pediatric/pediatric_enrichement.RData")
#STRING analysis--------------------------------------------------------------------------------------------
# STRING API: Download individual high-resolution networks
# for every adult protein
# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "12",
  species = 9606,
  score_threshold = 200
)

# Your adult genes
genes_IDR_pediatric <- unique(pediatric_expression_data$Gene)

# Map gene symbols to STRING IDs
pediatric_string_mapping <- string_db$map(
  data.frame(Gene = genes_IDR_pediatric),
  "Gene",
  removeUnmappedRows = TRUE
)
# Output directory
string_network_dir <- "enrichement analysis/pediatric/STRING_networks"

dir.create(
  string_network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# STRING version 12
string_api_url <- "https://version-12-0.string-db.org/api/highres_image/network"

# Keep track of successful and failed downloads
download_results <- data.frame(
  Gene = character(),
  STRING_ID = character(),
  Status = character(),
  File = character(),
  stringsAsFactors = FALSE
)

# Loop through every mapped adult protein

for (i in seq_len(nrow(pediatric_string_mapping))) {
  
  gene <- pediatric_string_mapping$Gene[i]
  string_id <- pediatric_string_mapping$STRING_id[i]
  
  cat(
    "\nDownloading:", gene,
    "| STRING ID:", string_id,
    "|", i, "of", nrow(pediatric_string_mapping), "\n"
  )
  
  # Output filename
  file_name <- file.path(
    string_network_dir,
    paste0(gene, "_STRING_network_highres.png")
  )
  
  # Parameters
  params <- list(
    identifiers = string_id,
    species = 9606,
    required_score = 200,
    network_flavor = "confidence",
    network_type = "functional",
    show_query_node_labels = 1,
    caller_identity = "pediatric_IDR_glioblastoma_analysis"
  )
  
  # Query STRING
  response <- tryCatch(
    
    POST(
      url = string_api_url,
      body = params,
      encode = "form"
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  # Check response
  if (!is.null(response) && status_code(response) == 200) {
    
    # Save PNG
    writeBin(
      content(response, as = "raw"),
      file_name
    )
    
    cat("  ✓ Saved:", file_name, "\n")
    
    download_results <- rbind(
      download_results,
      data.frame(
        Gene = gene,
        STRING_ID = string_id,
        Status = "Success",
        File = file_name,
        stringsAsFactors = FALSE
      )
    )
    
  } else {
    
    cat("  ✗ Failed:", gene, "\n")
    
    download_results <- rbind(
      download_results,
      data.frame(
        Gene = gene,
        STRING_ID = string_id,
        Status = "Failed",
        File = file_name,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # STRING requests should be separated by approximately 1 second
  Sys.sleep(1)
}
# Save download log
write.csv(
  download_results,
  file.path(
    "enrichement analysis/pediatric/STRING_networks",
    "STRING_network_download_log.csv"
  ),
  row.names = FALSE
)

# Summary
table(download_results$Status)


# ============================================================
# STRING interaction networks and degree for ALL adult genes
# ============================================================

string_api_url <- 
  "https://version-12-0.string-db.org/api/tsv/network"

# Results table
STRING_metrics_pediatric <- data.frame(
  Gene = character(),
  STRING_ID = character(),
  Nodes = numeric(),
  Edges = numeric(),
  Degree = numeric(),
  Clustering_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Store individual networks
STRING_networks_pediatric <- list()


# Loop through all adult genes


for (i in seq_len(nrow(pediatric_string_mapping))) {
  
  gene <- pediatric_string_mapping$Gene[i]
  string_id <- pediatric_string_mapping$STRING_id[i]
  
  cat(
    "\nAnalyzing:", gene,
    "|", i, "of", nrow(pediatric_string_mapping), "\n"
  )
  
  response <- tryCatch(
    
    POST(
      string_api_url,
      body = list(
        identifiers = string_id,
        species = 9606,
        
        # STRING confidence threshold
        required_score = 200,
        
        # IMPORTANT:
        # request a large neighborhood
        add_white_nodes = 1000,
        
        network_type = "functional"
      ),
      encode = "form"
    ),
    
    error = function(e) {
      message("ERROR: ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(response) && status_code(response) == 200) {
    
    txt <- content(
      response,
      as = "text",
      encoding = "UTF-8"
    )
    
    network_data <- read.delim(
      text = txt,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    
    if (nrow(network_data) > 0) {
      
      
      # Save individual STRING network
      
      
      STRING_networks_pediatric[[gene]] <- network_data
      
      
      # Build graph
      
      
      edges <- network_data %>%
        select(
          from = preferredName_A,
          to = preferredName_B
        ) %>%
        distinct()
      
      g <- graph_from_data_frame(
        edges,
        directed = FALSE
      )
      
      
      # Number of nodes
      
      
      number_nodes <- vcount(g)
      
      
      # Number of edges
      
      
      number_edges <- ecount(g)
      
      
      # Degree of the QUERY protein
      
      
      degree_all <- degree(
        g,
        mode = "all"
      )
      
      query_degree <- if (
        gene %in% names(degree_all)
      ) {
        degree_all[gene]
      } else {
        NA
      }
      
      
      # Clustering coefficient of QUERY protein
      
      
      clustering_all <- transitivity(
        g,
        type = "local",
        isolates = "zero"
      )
      
      query_clustering <- if (
        gene %in% names(clustering_all)
      ) {
        clustering_all[gene]
      } else {
        NA
      }
      
      
      # Save metrics
      
      
      STRING_metrics_pediatric <- rbind(
        STRING_metrics_pediatric,
        data.frame(
          Gene = gene,
          STRING_ID = string_id,
          Nodes = number_nodes,
          Edges = number_edges,
          Degree = as.numeric(query_degree),
          Clustering_coefficient = as.numeric(query_clustering),
          stringsAsFactors = FALSE
        )
      )
      
      cat(
        "   Nodes:", number_nodes,
        "| Edges:", number_edges,
        "| Degree:", query_degree,
        "| Clustering:", query_clustering,
        "\n"
      )
      
    } else {
      
      cat("   No interactions returned\n")
    }
    
  } else {
    
    cat("   STRING request failed\n")
  }
  
  Sys.sleep(1)
}


# Save results
write.csv(
  STRING_metrics_pediatric,
  "enrichement analysis/pediatric/STRING_networks/STRING_network_metrics_adult.csv",
  row.names = FALSE
)

# View
STRING_metrics_pediatric
#merging with the original dataset
IDR_pediatric_GO<-IDR_pediatric_GO%>%
  left_join(STRING_metrics_pediatric,
            by="Gene")

write.csv(IDR_pediatric_GO,"enrichement analysis/pediatric/STRING_analysis.csv")
save.image("enrichement analysis/pediatric/pediatric_enrichement.RData")
##hUB genes-------------------------------------------------------------------------------------------------
#Method A
mean_degree <- mean(IDR_pediatric_GO$Degree)
sd_degree   <- sd(IDR_pediatric_GO$Degree)
statistical_cutoff <- mean_degree + (2 * sd_degree)

hubs_statistical <- IDR_pediatric_GO %>%
  filter(Degree >= statistical_cutoff) %>%
  arrange(desc(Degree))
#Method B
IDR_pediatric_GO1<-IDR_pediatric_GO%>%
  filter(!is.na(Degree))
percentile_cutoff <- quantile(IDR_pediatric_GO1$Degree, 0.90)

hubs_percentile <- IDR_pediatric_GO1 %>%
  filter(Degree >= percentile_cutoff) %>%
  arrange(desc(Degree))
save.image("enrichement analysis/pediatric/pediatric_enrichement.RData")
