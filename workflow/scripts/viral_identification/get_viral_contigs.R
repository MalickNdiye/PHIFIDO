library(data.table)
library(optparse)
library(tidyverse)
library(Biostrings)

# Aguments parser
option_list <- list (make_option (c("-i","--infasta"), type="character",
                                   help="assembly file", metavar="character"),
                        make_option (c("-t","--intab"), type="character",
                                   help="viral contigs table", metavar="character"),            
                        make_option (c("-o","--outfasta"), type="character",
                                   help="filtered fasta file", metavar="character")
)

parser <-OptionParser(option_list=option_list)
args <- parse_args(parser)

assembly_df<- fread(args$intab)
new_names<- unlist(split(assembly_df$new_scaffold_name, as.character(assembly_df$scaffold)))

assembly<-readDNAStringSet(args$infasta)
names(assembly)<- new_names[names(assembly)]

to_keep<- names(assembly)[!is.na(names(assembly))]

viral_contigs<- assembly[to_keep]

tryCatch({
    writeXStringSet(viral_contigs, args$outfasta)
}, error = function(e) {
  # Create an empty file if an error occurs
  file.create(args$outfasta)
  message("No viral Contigs. Created an empty file.")
})
