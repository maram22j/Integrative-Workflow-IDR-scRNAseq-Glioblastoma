library(HGNChelper)
library(Biostrings)
library(UniProt.ws)
library(biomaRt)
library(protti)
library(dplyr)
library(httr)
system("aiupred -h")
#loading dataset
data<-read.csv("gene filtration/DE_adult.csv",header = T,sep = ",",row.names=1)
write.table(data$gene, file = "IDR conversion/adult/DE_genes_in_adult_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_adult<-readLines("IDR conversion/adult/DE_genes_in_adult_v2.txt")
#Selecting an Ensembl BioMart database and dataset
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
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
write.csv(conversion_clean,file = "IDR conversion/adult/conversion_gene_adult_to_uniprot_id.csv",row.names = FALSE)
#creating a new folder for storing fasta files of related genes
dir.create("IDR conversion/adult/FASTA_sequences_for_adult")
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
      "IDR conversion/adult/FASTA_sequences_for_adult/",
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
save.image("IDR conversion/adult/IDR_protein_for_adult_V2.RData")
# creating a new folder
dir.create("IDR conversion/adult/AIUPred_results")
# listing fasta files
fasta_files <- list.files(
  "IDR conversion/adult/FASTA_sequences_for_adult",
  pattern = "\\.fasta$",
  full.names = TRUE
)
#Run AIUPred on all FASTA files
#Disorder prediction
for(f in fasta_files){
  
  protein_id <- tools::file_path_sans_ext(
    basename(f)
  )
  
  output_file <- file.path(
    "IDR conversion/adult/AIUPred_results",
    paste0(protein_id, ".res")
  )
  
  
  command <- paste(
    "aiupred",
    "-i",
    shQuote(f),
    "-o",
    shQuote(output_file)
  )
  
  
  system(command)
  
}
#Binding prediction
dir.create("IDR conversion/adult/AIUPred_results/binding_prediction")
for(f in fasta_files){
  
  id <- tools::file_path_sans_ext(basename(f))
  
  system(
    paste(
      "aiupred",
      "-i", shQuote(f),
      "-b",
      "-o",
      shQuote(
        paste0(
          "IDR conversion/adult/AIUPred_results/binding_prediction/",
          id,
          "_binding.res"
        )
      )
    )
  )
}
dir.create("IDR conversion/adult/AIUPred_results/linker_protein")
#linker protein
for(f in fasta_files){
  
  id <- tools::file_path_sans_ext(basename(f))
  
  system(
    paste(
      "aiupred",
      "-i", shQuote(f),
      "-l",
      "-o",
      shQuote(
        paste0(
          "IDR conversion/adult/AIUPred_results/linker_protein/",
          id,
          "_linker.res"
        )
      )
    )
  )
}
# redox state
dir.create("IDR conversion/adult/AIUPred_results/redox_state")
for (f in fasta_files) {
  
  id <- tools::file_path_sans_ext(basename(f))
  
  output_file <- file.path(
    "IDR conversion/adult/AIUPred_results/redox_state/",
    paste0(id, "_redox.res")
  )
  
  system(
    paste(
      "aiupred",
      "-i", shQuote(f),
      "-r",
      "-o", shQuote(output_file)
    )
  )
}
# Convert residue predictions into protein-level scores
calculate_AIUPred_summary <- function(file){
  
  data <- read.table(
    file,
    header = TRUE,
    comment.char = "#"
  )
  
  score <- data[[3]]
  
  data.frame(
    
    UniProt_ID =
      tools::file_path_sans_ext(
        basename(file)
      ),
    
    mean_score =
      mean(score, na.rm = TRUE),
    
    disordered_residues =
      sum(score > 0.5, na.rm = TRUE),
    
    disorder_percentage =
      100 * mean(score > 0.5, na.rm = TRUE)
    
  )
}

#apply to all proteins
disorder_files <- list.files(
  "IDR conversion/adult/AIUPred_results",
  pattern="\\.res$",
  full.names=TRUE
)


AIUPred_summary <- do.call(
  rbind,
  lapply(
    disorder_files,
    calculate_AIUPred_summary
  )
)
final_AIUPred <- conversion_clean %>%
  left_join(
    AIUPred_summary,
    by="UniProt_ID"
  )
write.csv(final_AIUPred,file = "IDR conversion/adult/AIUpred_adult_results.csv",row.names = FALSE)
save.image("IDR conversion/adult/IDR_conversion_adult.RData")
 