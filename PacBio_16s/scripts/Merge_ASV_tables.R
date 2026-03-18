if(!require(dada2)){
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}
library("tidyverse")
library("data.table")
library(Biostrings)


# --------------------------
# Set Variables
# --------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6){
  stop(" Usage: Merge_asv_tables <asv_table1> <asv_table2> <asv_table3> <tax_db_genus> <tax_db_species> <output_dir>", call.=FALSE)
} else {
  ASV_samples_table_noChim1_p <- args[1] # folder with all pre-processed reads
  ASV_samples_table_noChim2_p <- args[2] # folder with all pre-processed reads
  ASV_samples_table_noChim3_p <- args[3] # folder with all pre-processed reads
  Tax_db_toGenus <- args[4] # path to fasta file with species-level taxonomy for assignSpecies
  Tax_db_toSpecies <- args[5] # path to fasta file with genus-level taxonomy for assign
  output_dir <- args[6] # output directory
}


#outpout directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --------------------------
# Load Data
# --------------------------
# open ASV table
ASV_samples_table_noChim1<- readRDS(ASV_samples_table_noChim1_p) 
ASV_samples_table_noChim2<- readRDS(ASV_samples_table_noChim2_p) 
ASV_samples_table_noChim3<- readRDS(ASV_samples_table_noChim3_p) 
# make namees of table 3 unique
rownames(ASV_samples_table_noChim3) <- paste0(rownames(ASV_samples_table_noChim3), "-Bicom6")
# ------------------------------------------------------------------

cat("=== Merging DADA2 ASV tables from two datasets ===\n")

# 1. Report ASV counts before merging
n_asv1 <- ncol(ASV_samples_table_noChim1)
n_asv2 <- ncol(ASV_samples_table_noChim2)
n_asv3 <- ncol(ASV_samples_table_noChim3)

cat("ASVs in table 1:", n_asv1, "\n")
cat("ASVs in table 2:", n_asv2, "\n")
cat("ASVs in table 3:", n_asv3, "\n")


# 2. Merge the ASV tables by sequence
seqtab_merged <- mergeSequenceTables(ASV_samples_table_noChim1,
                                     ASV_samples_table_noChim2,
                                     ASV_samples_table_noChim3)

n_asv_merged <- ncol(seqtab_merged)
cat("Unique ASVs after merging:", n_asv_merged, "\n")


# 3. Remove chimeras again globally
seqtab_nochim <-  removeBimeraDenovo(seqtab_merged, method="pool", multithread=TRUE)

n_asv_final <- ncol(seqtab_nochim)
cat("ASVs after second chimera removal:", n_asv_final, "\n")
cat("Chimeras removed in second pass:", n_asv_merged - n_asv_final, "\n")


# 4. Assign taxonomy
ASV_taxonomy <- assignTaxonomy(seqs = seqtab_nochim,
                               refFasta = Tax_db_toGenus)
ASV_taxonomy<-ASV_taxonomy[, colnames(ASV_taxonomy) != "Species"]

ASV_taxonomy_sp<- addSpecies(ASV_taxonomy,
                             refFasta = Tax_db_toSpecies )


# 5. Save final Tables
saveRDS(seqtab_nochim, file = file.path(output_dir, "ASV_table_merged.rds"))
saveRDS(ASV_taxonomy, file.path(output_dir,"ASV_Taxonomy.RDS"))
saveRDS(ASV_taxonomy_sp, file.path(output_dir,"ASV_Taxonomy_sp.RDS"))
