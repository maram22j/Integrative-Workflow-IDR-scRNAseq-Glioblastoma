library(HGNChelper)
library(Biostrings)
library(UniProt.ws)
library(biomaRt)
library(protti)
library(dplyr)
library(IRanges)
#loading dataset
data<-read.csv("gene filtration/filtred_genes_adult.csv",header = T,sep = ",",row.names=1)
write.table(data$gene, file = "IDR conversion/DE_genes_in_adult_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_adult<-readLines("DE_genes_in_adult_v2.txt")
#Selecting an Ensembl BioMart database and dataset
mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)
# Map HGNC gene symbols to Ensembl and UniProt protein identifiers
conversion <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "ensembl_peptide_id",
    "uniprotswissprot"
  ),
  filters = "hgnc_symbol",
  values = genes_adult,
  mart = mart
)
head(conversion)
#renaming hgnc_symbol and uniprotswissprot
conversion <- dplyr::rename(
  conversion,
  gene = hgnc_symbol,
  UniProt_ID = uniprotswissprot
)
#filtering genes that do not have uniprot id 
conversion_clean <- conversion %>%
  dplyr::filter(UniProt_ID != "") %>%
  dplyr::distinct(gene, .keep_all = TRUE)
#seeing the missing genes
missing_genes <- setdiff(
  genes_adult,
  conversion_clean$gene
)
missing_genes
write.table(missing_genes,"IDR conversion/missing_genes_in_adult.txt",quote = FALSE,row.names = FALSE, col.names = FALSE)
# Check for outdated HGNC symbols
check <- checkGeneSymbols(missing_genes)
# Keep only corrected symbols
updated_symbols <- check %>%
  dplyr::filter(!is.na(Suggested.Symbol))
# Retrieve UniProt IDs for corrected symbols
conversion_updated <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "ensembl_peptide_id",
    "uniprotswissprot"
  ),
  filters = "hgnc_symbol",
  values = updated_symbols$Suggested.Symbol,
  mart = mart
)
conversion_updated <- conversion_updated %>%
  dplyr::rename(
    gene = hgnc_symbol,
    UniProt_ID = uniprotswissprot
  ) %>%
  dplyr::filter(UniProt_ID != "") %>%
  dplyr::distinct(gene, .keep_all = TRUE)
# Add the recovered genes to the main table
conversion_clean <- dplyr::bind_rows(
  conversion_clean,
  conversion_updated
) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
# Connect to the UniProt database for Homo sapiens (NCBI taxonomy ID: 9606) to retrieve human protein annotation information
up <- UniProt.ws(taxId = 9606)
#Extracting protein length
protein_length <- AnnotationDbi::select(
  up,
  keys = conversion_clean$UniProt_ID,
  columns = c("accession", "length"),
  keytype = "UniProtKB"
)
head(protein_length)
#renaming the coloumn entry to accession 
protein_length <- dplyr::rename(protein_length,
                                accession = Entry,
                                length = Length
)
protein_length <- protein_length %>%
  dplyr::select(-From)
protein_length <- protein_length %>%
  mutate(
    length = as.numeric(length))
#Joining the conversion_clean with final_IDR 
conversion_clean <- conversion_clean %>%
  left_join(
    protein_length,
    by = c("UniProt_ID" = "accession")
  )

#saving genes converted to uniprot_id
write.csv(conversion_clean,file = "IDR conversion/conversion_gene_adult_to_uniprot_id.csv",row.names = FALSE)
#creating a new folder for storing fasta files of related genes
dir.create("FASTA_sequences_for_adult")
# Function to retrieve UniProt FASTA sequences

get_uniprot_fasta <- function(uniprot_id){
  
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/",
    uniprot_id,
    ".fasta"
  )
  
  response <- GET(url)
  
  if(response$status_code == 200){
    
    fasta_text <- content(response, "text")
    
    return(fasta_text)
    
  } else {
    
    return(NULL)
    
  }
}
# Download and save UniProt FASTA sequences
for(id in conversion_clean$UniProt_ID){
  
  fasta <- get_uniprot_fasta(id)
  
  if(!is.null(fasta)){
    
    file_name <- paste0(
      "FASTA_sequences_for_adult/",
      id,
      ".fasta"
    )
    
    writeLines(
      fasta,
      file_name
    )
  }
}
length(list.files("FASTA_sequences_for_adult",pattern = "\\.fasta$"))
save.image("IDR conversion/IDR_protein_for_adult_V2.RData")
