library(data.table)
library(tidyverse)
library(optparse)

# this is the command: Rscript scripts/phage_host_link/parse_phage_host_links.R -b {input.bacteria_metadata} -s {input.spacers_metadata} -a {input.blastout} -o {output}
# Define command line arguments
option_list <- list(
    make_option(c("-b", "--bacteria_metadata"), type = "character",
                help = "Path to bacteria metadata file", metavar = "character"),
    make_option(c("-s", "--spacers_metadata"), type = "character",
                help = "Path to spacers metadata file", metavar = "character"),
    make_option(c("-a", "--blastout"), type = "character",
                help = "Path to BLAST output file", metavar = "character"), 
    make_option(c("-o", "--output"), type = "character",
                help = "Output file for parsed phage-host links", metavar = "character")    
)   

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

cat("Reading input files...\n")
# Read input files
bacteria_metadata <- fread(args$bacteria_metadata, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
spacers_metadata <- fread(args$spacers_metadata, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)     
blastout <- fread(args$blastout, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)

# use bacteria metadata to get named vector of genome to genus
bacteria_metadata_vec <- setNames(bacteria_metadata$genus, bacteria_metadata$genome)

# use spacers metadata to get named vector of spacer_id to length
spacers_metadata<- spacers_metadata %>%
    mutate(spacer_length = nchar(SpacerSeq),
            spacer_id=paste(Strain, Orientation, ShortID, sep="_"),
            genus=bacteria_metadata_vec[Strain])
spacer_length_vec <- setNames(spacers_metadata$spacer_length, spacers_metadata$spacer_id)
spacer_genus_vec <- setNames(spacers_metadata$genus, spacers_metadata$spacer_id)

# add spacer length to blastout
blastout <- blastout %>%
    mutate(spacer_length=spacer_length_vec[SPACER_ID],
           genus=spacer_genus_vec[SPACER_ID],
           true_mm=spacer_length - alignement_length + mismatch)

# filter blastout for true mismatches <=2
cat("Total number of phage-host links found before filtering:", nrow(blastout), "\n")
blastout_filtered <- blastout %>%
    filter(true_mm <= 2) 

# summarize results in table
cat("Number of unique phage-host links found:", length(unique(blastout_filtered$Query)), "\n")
phage_host_links <- blastout_filtered %>%
    group_by(Query) %>%
    summarize(host_genus=paste(unique(genus), collapse=";")) %>%
    ungroup() %>%
    rename(virus=Query)

# write output
fwrite(phage_host_links, args$output, sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = TRUE)
cat("Parsed phage-host links written to", args$output, "\n")


