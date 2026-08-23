library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ReactomePA)
library(ggtree)
library(pathview)
library(httr)
library(STRINGdb)
#adult joining tables
IDR_adult<-read.csv("IDR conversion/adult/AIUpred_adult_results.csv",header = T,sep = ",")
colnames(IDR_adult)[colnames(IDR_adult)=="gene"]<-"Gene"
adult_expression_data<-read.csv("gene filtration/DE_adult.csv",header=T,sep = ",")
adult_expression_data$X<-NULL
colnames(adult_expression_data)[colnames(adult_expression_data)=="gene"]<-"Gene"
adult_expression_data<-adult_expression_data%>%
  left_join(IDR_adult,by="Gene")
genes_IDR_adult <- adult_expression_data %>%
  pull(Gene)
writeLines(unique(genes_IDR_adult), "enrichement analysis/adult/adult_genes_PANTHER.txt")
# Gene ontology ------------------------
genes_IDR_adult <- adult_expression_data %>%
  pull(Gene)
gene_entrez_adult <- bitr(
  genes_IDR_adult,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Biological process
ego_BP_adult <- enrichGO(
  gene = gene_entrez_adult$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.06 ,
  qvalueCutoff = 0.06,
  readable = TRUE
)
BP<-as.data.frame(ego_BP_adult)
#Molecular Function
ego_MF_adult<- enrichGO(
  gene = gene_entrez_adult$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.06,
  qvalueCutoff = 0.06,
  readable = TRUE
)
MF<-as.data.frame(ego_MF_adult)
#Cellular component
ego_CC_adult<-enrichGO(
  gene = gene_entrez_adult$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
CC<-as.data.frame(ego_CC_adult)
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
BP_genes <- extract_GO_by_gene(ego_BP_adult)
colnames(BP_genes)[2] <- "BP"
MF_genes <- extract_GO_by_gene(ego_MF_adult)
colnames(MF_genes)[2] <- "MF"
CC_genes <- extract_GO_by_gene(ego_CC_adult)
colnames(CC_genes)[2] <- "CC"
# Joining GO with the expression dataset
adult_IDR_GO <- adult_expression_data %>%
  left_join(BP_genes, by = "Gene") %>%
  left_join(MF_genes, by = "Gene") %>%
  left_join(CC_genes, by = "Gene")
colnames(adult_IDR_GO)[colnames(adult_IDR_GO)=="BP"]<-"Biological Process"
colnames(adult_IDR_GO)[colnames(adult_IDR_GO)=="MF"]<-"Molecular Function"
colnames(adult_IDR_GO)[colnames(adult_IDR_GO)=="CC"]<-"Cellular Component"
write.csv(adult_IDR_GO,file = "enrichement analysis/adult/adult_IDR_GO.csv",row.names = F)

#plots --------------
# Biological process
dotplot(
  ego_BP_adult,
  showCategory = 20,
  title= "Gene Ontology for Adult for Biological Process"
)
ggsave("enrichement analysis/adult/GO_enrichement_adult_biological_process.png",width = 10, height = 10)
barplot(
  ego_BP_adult,
  showCategory = 20,
  title = "Gene Ontology for Adult for Biological Process"
)
ggsave("enrichement analysis/adult/GO_enrichement_adult_barplot_biological_process.png",width = 10, height = 10)
#Clustering plots 
# Calculate similarity
ego_BP_adult_sim <- pairwise_termsim(
  ego_BP_adult)

# Network
emapplot(
  ego_BP_adult_sim,
  showCategory = 30,
  layout = "kk",
  cex_category = 0.5,
  cex_line = 0.5
)
ggsave(filename = "enrichement analysis/adult/emapplot_BP_adult.png",width = 15,height = 15,dpi = 300)

ego_BP_simplified <- clusterProfiler::simplify(
  ego_BP_adult,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_BP_adult_top <- pairwise_termsim(
  ego_BP_simplified)
p <- treeplot(
  ego_BP_adult_top,
  cluster_method = "average",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize_cladelab = 4
) +
  hexpand(0.3)

print(p)
ggsave(filename = "enrichement analysis/adult/treeplot_BP_adult.png",plot = p,width = 15,height = 15,dpi = 300)
#Cellular Component
dotplot(
  ego_CC_adult,
  showCategory = 18,
  title="Gene Ontology for Adult for Cellular Component"
)
ggsave("enrichement analysis/adult/GO_enrichement_adult_cellular_component.png",width = 10, height = 10)
barplot(
  ego_CC_adult,
  showCategory = 20,
  title = "Gene Ontology for Adult for Cellular Component"
)
ggsave("enrichement analysis/adult/GO_enrichement_adult_barplot_cellular_component.png",width = 10, height = 10)
#Clustering plots 
# Calculate similarity
ego_CC_adult_sim <- pairwise_termsim(
  ego_CC_adult)

# Network
emapplot(
  ego_CC_adult_sim,
  showCategory = 30,
  layout = "kk",
  cex_category = 0.5,
  cex_line = 0.5
)
ggsave(filename = "enrichement analysis/adult/emapplot_CC_adult.png",width = 15,height = 15,dpi = 300)
#treeplot
ego_CC_simplified <- clusterProfiler::simplify(
  ego_CC_adult,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_CC_adult_top <- pairwise_termsim(
  ego_CC_simplified)
p <- treeplot(
  ego_CC_adult_top,
  cluster_method = "average",
  cladelab_offset = 8,
  tiplab_offset = 0.3,
  fontsize_cladelab = 4
) +
  hexpand(0.3)

print(p)
ggsave(filename = "enrichement analysis/adult/treeplot_CC_adult.png",plot = p,width = 15,height = 15,dpi = 300)
# Molecular Function
"""dotplot(
  ego_MF_adult,
  showCategory = 18,
  title="Gene Ontology for Adult for Molecular Function"
)
ggsave("satistical_test/Adult/GO_enrichement_adult_molecular_function.png",width = 10, height = 10)
barplot(
  ego_MF_adult,
  showCategory = 20,
  title = "Gene Ontology for Adult for Molecular Function"
)
ggsave("satistical_test/Adult/GO_enrichement_adult_barplot_Molecular_function.png",width = 10, height = 10)"""
#KEGG enrichement----------------------------------------------------------------------------------------------
ekegg_adult <- enrichKEGG(
  gene = gene_entrez_adult$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.06,
  qvalueCutoff = 0.06
)
kk<-as.data.frame(ekegg_adult)
# Convert Entrez IDs to gene symbols
ekegg_adult_readable <- setReadable(
  ekegg_adult,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)
# Extract KEGG pathways associated with each gene
KEGG_genes <- as.data.frame(ekegg_adult_readable) %>%
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
#joining the KEGG coloum with original dataset
adult_IDR_GO <- adult_IDR_GO %>%
  left_join(KEGG_genes, by = "Gene")
write.csv(adult_IDR_GO,file = "enrichement analysis/adult/adult_IDR_KEGG.csv",row.names = F)
#plotting pathview of KEGG
gene_fc <- adult_IDR_GO$logFC
names(gene_fc) <- gene_entrez_adult$ENTREZID
pathview(
  gene.data = gene_fc,
  pathway.id = "hsa04146",
  species = "hsa",
  gene.idtype = "entrez",
  out.suffix = "adult" 
)
# Reactome enrichement analysis----------------------------------------------------------------------------
reactome_adult <- enrichPathway(
  gene = gene_entrez_adult$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.06,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.06,
  readable = TRUE
)
reactome<-as.data.frame(reactome_adult)
Reactome_genes <- as.data.frame(reactome_adult) %>%
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
adult_IDR_GO <- adult_IDR_GO %>%
  left_join(Reactome_genes, by = "Gene")
#plots
dotplot(
  reactome_adult,
  showCategory = 20,
  title = "Reactome Pathway Enrichment - Adult"
)
ggsave("enrichement analysis/adult/Reactome_enrichment_adult_dotplot.png", width = 10, height = 10)
barplot(
  reactome_adult,
  showCategory = 20,
  title = "Reactome Pathway Enrichment - Adult"
)
ggsave("enrichement analysis/adult/Reactome_enrichment_adult_barplot.png", width = 10, height = 10)

write.csv( adult_IDR_GO, "enrichement analysis/adult/adult_IDR_GO_Reactome.csv", row.names = FALSE)

#STRING analysis--------------------------------------------------------------------------------------------
# STRING API: Download individual high-resolution networks
# for every adult protein
# Output directory
# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "12",
  species = 9606,
  score_threshold = 200
)

# Your adult genes
genes_IDR_adult <- unique(adult_expression_data$Gene)

# Map gene symbols to STRING IDs
adult_string_mapping <- string_db$map(
  data.frame(Gene = genes_IDR_adult),
  "Gene",
  removeUnmappedRows = TRUE
)
string_network_dir <- "enrichement analysis/adult/STRING_networks"

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

for (i in seq_len(nrow(adult_string_mapping))) {
  
  gene <- adult_string_mapping$Gene[i]
  string_id <- adult_string_mapping$STRING_id[i]
  
  cat(
    "\nDownloading:", gene,
    "| STRING ID:", string_id,
    "|", i, "of", nrow(adult_string_mapping), "\n"
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
    caller_identity = "adult_IDR_glioblastoma_analysis"
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
      httr::content(response, as = "raw"),
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
    "enrichement analysis/adult",
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
STRING_metrics_adult <- data.frame(
  Gene = character(),
  STRING_ID = character(),
  Nodes = numeric(),
  Edges = numeric(),
  Degree = numeric(),
  Clustering_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Store individual networks
STRING_networks_adult <- list()


# Loop through all adult genes


for (i in seq_len(nrow(adult_string_mapping))) {
  
  gene <- adult_string_mapping$Gene[i]
  string_id <- adult_string_mapping$STRING_id[i]
  
  cat(
    "\nAnalyzing:", gene,
    "|", i, "of", nrow(adult_string_mapping), "\n"
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
    
    txt <- httr::content(
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
      
      
      STRING_networks_adult[[gene]] <- network_data
      
      
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
      
      
      STRING_metrics_adult <- rbind(
        STRING_metrics_adult,
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
  STRING_metrics_adult,
  "enrichement analysis/adult/STRING_network_metrics_adult.csv",
  row.names = FALSE
)

# View
STRING_metrics_adult
#merging with the original dataset
adult_IDR_GO<-adult_IDR_GO%>%
  left_join(STRING_metrics_adult,
            by="Gene")

write.csv(adult_IDR_GO,"enrichement analysis/adult/STRING_analysis.csv")
#HUB genes-------------------------------------------------------------------------------------------------
IDR_adult_GO1<-adult_IDR_GO%>%
  filter(!is.na(Degree))
percentile_cutoff <- quantile(IDR_adult_GO1$Degree, 0.90)

hubs_percentile <- IDR_adult_GO1 %>%
  filter(Degree >= percentile_cutoff) %>%
  arrange(desc(Degree))

save.image("enrichement analysis/adult/Adult_enrichement.RData")


































