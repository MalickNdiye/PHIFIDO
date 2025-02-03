library(argparse)
library(tidyverse)
library(data.table)

# Parse arguments "-d, --drep", "-v, --votu", "-a, --ANI" "-o, --output"
parser <- ArgumentParser()
parser$add_argument("-d", "--drep", help = "dRep output directory")
parser$add_argument("-v", "--votu", help = "vOTU target vOTU")
parser$add_argument("-a", "--ANI", help = "ANI table")
parser$add_argument("-o", "--output", help = "output directory")
args <- parser$parse_args()

# gbk path
gbk_path <- file.path("../scratch_link/annotations_vOTU", args$votu, "gbk")

# Load drep table
drep <- read_tsv(args$drep)
drep_filtered <- drep %>% filter(vOTU == args$votu)

genomes<- unique(drep_filtered$genome)
genome_lengths<- unique(drep_filtered$length)

genome2vOTU<- split(drep$vOTU, drep$genome)


# Load ANI table
ANI <- fread(args$ANI, header = F)
colnames(ANI) <- c("genome1", "genome2", "ANI", "alignmentd_fragments", "total_fragments")
ANI$AF=ANI$alignmentd_fragments/ANI$total_fragments
ANI$ANI=ANI$ANI*ANI$AF

# genome1 and genome2 are paths to the genomes, we only need the genome names without the path and extension
ANI$genome1 <- basename(ANI$genome1)
ANI$genome2 <- basename(ANI$genome2)
# remove extension
ANI$genome1 <- gsub("\\.fasta", "", ANI$genome1)
ANI$genome2 <- gsub("\\.fasta", "", ANI$genome2)
# Add vOTU to ANI table
ANI$vOTU1 <- genome2vOTU[ANI$genome1]
ANI$vOTU2 <- genome2vOTU[ANI$genome2]

head(ANI)


# transform ANI table to matrix
ANI_matrtix <- ANI %>% 
    filter(vOTU1 == args$votu | vOTU2 == args$votu) %>%
    select(genome1, genome2, ANI) %>%
    pivot_wider(names_from = genome2, values_from = ANI, values_fill=0) %>%
    column_to_rownames("genome1") %>%
    as.matrix()

# The diagonal of the matrix should be 100% (ANI of a genome with itself)
diag(ANI_matrtix) <- 100

# order ANI matrix as genome
ANI_matrtix <- ANI_matrtix[genomes, genomes]

# Function to order genomes based on ANI and length
# start taking longest genome
# then take the genome with the highest ANI to the longest genome (if tie, take the longest genome)
# then take the genome with the highest ANI to the second genome (if tie, take the longest genome)
# and so on
order_genomes_by_ani <- function(genomes, genome_lengths, ANI_matrtix){
    ordered_genomes <- c()
    while(length(ordered_genomes) < length(genomes)){
        if(length(ordered_genomes) == 0){
            # take the longest genome
            ordered_genomes <- c(ordered_genomes, genomes[which.max(genome_lengths)])
        }else{
            # take the genome with the highest ANI to the last genome in the ordered list
            last_genome <- ordered_genomes[length(ordered_genomes)]
            ANI_to_last_genome <- ANI_matrtix[last_genome, genomes]
            ANI_to_last_genome[ordered_genomes] <- 0
            next_genome <- genomes[which.max(ANI_to_last_genome)]
            ordered_genomes <- c(ordered_genomes, next_genome)
        }
    }
    return(ordered_genomes)
}


# Get the ordered genomes based on ANI
ordered_genomes <- order_genomes_by_ani(genomes, genome_lengths, ANI_matrtix)
ordered_genomes<- unique(ordered_genomes)
# add gbk path to the ordered genomes
ordered_genomes <- paste0(gbk_path, "/", ordered_genomes, ".gbk")
print("Ordered Genomes:")
print(ordered_genomes)

# Save ordered genomes in a txt file. all genomes on a single line separated by a space. the file will have only one line
write.table(t(ordered_genomes), file = args$output, row.names = F, col.names = F, quote = F, sep = " ")
