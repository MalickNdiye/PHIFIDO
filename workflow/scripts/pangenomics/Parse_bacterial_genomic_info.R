library(optparse)
library(data.table)
library(tidyverse)

option_list <- list(
    make_option(c("-g", "--genomic_info"), type = "character",
                help = "Genomic information file", metavar = "character"),
    make_option(c("-c", "--checkm"), type = "character",
                help = "CheckM file", metavar = "character"),
    make_option(c("-d", "--gtdb"), type = "character",
                help = "GTDB classification file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output file for parsed genomic information", metavar = "character")
)

parser <-OptionParser(option_list=option_list)
args <- parse_args(parser)

# read input file
print("Reading input files...")
print(args$genomic_info)
print(args$checkm)
print(args$gtdb)
gen_info <- fread(args$genomic_info)
checkm <- fread(args$checkm)
gtdb<- fread(args$gtdb)

print("Input files read successfully.")

# function to parse GTDB classification
print("Parsing GTDB classification...")
# this function will split the classification string into a list and clean it up
parse_gtdb_classification <- function(classification) {
    class_list <- strsplit(classification, ";")[[1]]
    class_list <- trimws(class_list)
    class_list <- gsub("d__|p__|c__|o__|f__|g__|s__", "", class_list)
    names(class_list) <- c("domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


    return(class_list)
}

gtdb_parsed <- gtdb %>%
    rowwise()%>%
    mutate(Classification = list(sapply(classification, parse_gtdb_classification))) %>%
    unnest_wider(Classification) %>%
    mutate(Genus= gsub("Bombiscardovia", "Bifidobacterium", Genus),
           Species= gsub("Bombiscardovia", "Bifidobacterium", Species),
           Species= gsub("polysaccharolytica", "polysaccharolyticum", Species)) 

print(colnames(gtdb_parsed))


# select columns of interest
print("Selecting relevant columns...")
checkm_selected <- checkm %>%
    select("Bin Id", "Completeness", "Contamination", "Genome size (bp)", "# contigs", "N50 (contigs)") %>%
    rename("Genome" = "Bin Id", "Genome_Size_bp"= "Genome size (bp)", "Nr_contigs" = "# contigs", "N50_contigs" = "N50 (contigs)")

gen_info_selected <- gen_info %>%
    select("genome", "strain", "Isolation_source", "BioSample") %>%
    rename("Genome"="genome", "Strain"="strain") 

gtdb_selected <- gtdb_parsed %>%
    select("user_genome", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>%
    rename("Genome" = "user_genome")

# merge data frames
print("Merging data frames...")
merged_data <- checkm_selected %>%
    left_join(., gen_info_selected, by = "Genome") %>%
    left_join(., gtdb_selected, by = "Genome") %>%
    mutate(Isolation_source = ifelse(is.na(Isolation_source), "Apis Mellifera", Isolation_source),
            Isolation_source = ifelse(Species=="Bifidobacterium breve", "Unkown (Type Strain)", Isolation_source),
           BioSample = ifelse(is.na(BioSample), "TBD", BioSample),
           Strain = ifelse(is.na(Strain), Genome, Strain)) %>%
    select(Genome, Strain, Phylum, Class, Order, Family, Genus, Species, Isolation_source, Completeness, Contamination, Genome_Size_bp, Nr_contigs, N50_contigs, BioSample)


# write output to file
write.table(merged_data, args$output, sep = ",", quote = FALSE, row.names = FALSE)
message("Parsed bacterial genomic information saved to ", args$output)
