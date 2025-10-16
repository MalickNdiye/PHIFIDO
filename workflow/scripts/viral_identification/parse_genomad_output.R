library(data.table)
library(optparse)
library(tidyverse)
library(Biostrings)

# Aguments parser
option_list <- list ( make_option (c("-d","--dir"), type="character",
                                   help="genomad directory", metavar="character"),
                      make_option (c("-o","--outtab"), type="character",
                                   help="genomad output table", metavar="character"),
                        make_option (c("-f","--outfasta"), type="character",
                                      help="genomad output fasta", metavar="character")
)

parser <-OptionParser(option_list=option_list)
args <- parse_args(parser)

# setup, input directory
dir<- args$dir
print(dir)
sample<- unlist(strsplit(basename(args$dir), split="_"))[2]

# set path to relevant files
viruses_info<- file.path(dir, paste0(sample, "_summary"), paste0(sample, "_virus_summary.tsv"))
viruses_fasta<- file.path(dir, paste0(sample, "_summary"), paste0(sample, "_virus.fna"))

# open tablee
if (file.exists(viruses_info)){
    virus_df<- fread(viruses_info) 
} else {
    print(paste0("No virus summary file found for ", sample))
    quit(save = "no", status = 1)
}

virus_df_formatted<- virus_df %>%
    mutate(sample= sample,
           scaffold=seq_name,
           new_scaffold_name=paste(sample, "_viral_contig_", row_number(), sep=""),
           type= "lysogenic",
           vibrant_quality= "GeNomad prophage",
           total_genes=n_genes,
           multiplicity=NA) %>%
        relocate(sample, scaffold, new_scaffold_name, type, vibrant_quality, total_genes, length, multiplicity) %>%
        select(sample, scaffold, new_scaffold_name, type, vibrant_quality, total_genes, length, multiplicity)

print(head(virus_df_formatted))

# write table
print(paste0("Writing table file: ", viruses_info,  " to ", args$outtab, "with ", nrow(virus_df_formatted), " rows"))
fwrite(virus_df_formatted, args$outtab, quote=F, row.names = F, sep="\t" )

# Open fasta file
open_fasta<- readDNAStringSet(viruses_fasta)

# use table to rename sequences
names(open_fasta)<- virus_df_formatted$new_scaffold_name[match(names(open_fasta), virus_df_formatted$scaffold)]

# write fasta file
print(paste0("Writing fasta file: ", viruses_fasta,  " to ", args$outfasta, "with ", length(open_fasta), " sequences with changed names"))
writeXStringSet(open_fasta, args$outfasta)