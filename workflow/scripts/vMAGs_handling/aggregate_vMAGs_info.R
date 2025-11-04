library(optparse)
library(data.table)
library(tidyverse)

option_list <- list(
    make_option(c("-c", "--checkv"), type = "character",
                help = "checkv folder", metavar = "character"),
    make_option(c("-l", "--lifestyle"), type = "character",
                help = "lifestlye folder", metavar = "character"),
    make_option(c("-t", "--taxonomy"), type = "character",
                help = "taxonomy folder", metavar = "character"),
    make_option(c("-p", "--PhageHost"), type = "character",
                help = "phage-host linkage", metavar = "character"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output file for parsed genomic information", metavar = "character")
)

parser <-OptionParser(option_list=option_list)
args <- parse_args(parser)

cat("Reading input files...\n")
# Read CheckV output
checkv <- read.table(file.path(args$checkv, "quality_summary.tsv"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

lifeyle<- read.table(file.path(args$lifestyle, "all_viral_contigs.fasta.bacphlip"),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

taxonomy <- read.table(file.path(args$taxonomy, "final_prediction/phagcn_prediction.tsv"), 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

phage_host <- read.table(args$PhageHost,
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
              rename(contig_id=virus)

# filter tables
cat("Processing CheckV information...\n")
checkv_filt <- checkv %>%
  select(contig_id, contig_length, checkv_quality)

cat("Processing lifestyle information...\n")
colnames(lifeyle)[1] <- "contig_id"
lifestyle_filt <- lifeyle %>%
    mutate(Lifestyle = ifelse(Virulent >= Temperate, "Lytic", "Temperate"),
           Lifestyle=ifelse(grepl("ESL", contig_id), "Temperate", Lifestyle)) %>%
    select(contig_id, Lifestyle)


format_lineage <- function(lineage) {
  lineage_split <- unlist(strsplit(lineage, ";"))
  
  lineage_list <- lapply(lineage_split, function(x) {
    parts <- unlist(strsplit(x, ":"))
    if(length(parts) == 2) {
      setNames(parts[2], parts[1])
    } else {
      NULL
    }
  })
  
  # Combine into a single named vector
  lineage_vec <- unlist(lineage_list, use.names = TRUE)
  
  # Keep only the first occurrence of each rank
  lineage_vec <- lineage_vec[!duplicated(names(lineage_vec))]
  
  return(lineage_vec)
}

cat("Processing taxonomy information...\n")
taxonomy_filt <- taxonomy %>%
  select(Accession, Lineage) %>%
  mutate(Lineage = lapply(Lineage, format_lineage)) %>%
  tidyr::unnest_wider(Lineage) %>%
  rename(contig_id = Accession) %>%
  relocate(order, family, subfamily, .before = genus)

# Merge tables
cat("Merging information...\n")
merged <- checkv_filt %>%
    left_join(lifestyle_filt, by = "contig_id") %>%
    left_join(taxonomy_filt, by = "contig_id") %>%
    left_join(phage_host, by = "contig_id") %>%
    rename(genome=contig_id)
    
# Write output
print("Writing output...")
write.table(merged, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)

print("Done.")