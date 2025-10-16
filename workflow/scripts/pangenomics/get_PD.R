suppressPackageStartupMessages(library(argparse))


# Packages
library("ape")
library(tidyverse)

if (!require("ggtree", character.only = TRUE)) {
  BiocManager::install("ggtree", ask = FALSE)
  require("ggtree", character.only = TRUE)
}

if (!require("picante", character.only = TRUE)) {
  install.packages("picante", repos = "https://cloud.r-project.org")
  require("picante", character.only = TRUE)
}



source("scripts/General_functions.R")

# Create parser
parser <- ArgumentParser(description = "Calculate phylogenetic diversity (PD) from a tree")

# Add arguments
parser$add_argument(
  "-i", "--input",
  required = TRUE,
  help = "Input tree file (Newick format)"
)

parser$add_argument(
  "-o", "--output",
  required = TRUE,
  help = "Output file for PD results"
)

# Parse arguments
args <- parser$parse_args()

# Access them like this:
input_tree <- args$input
output_file <- args$output

# Read in the tree
tree<- read.tree(input_tree) %>%
  as_tibble() %>%
  mutate(label=gsub("_genes", "", label)) %>%
  as.phylo()


# Function to compute Faith's PD between two strains
faith_pd_pair <- function(tree, strain1, strain2) {
  # check strains exist in tree
  if (!(strain1 %in% tree$tip.label) | !(strain2 %in% tree$tip.label)) {
    print(strain1)
    print(strain2)
    stop("One or both strains not found in tree")
  }
  
  # make a "community matrix": rows = communities, cols = species
  comm <- matrix(0, nrow = 1, ncol = length(tree$tip.label),
                 dimnames = list("pair", tree$tip.label))
  comm[1, strain1] <- 1
  comm[1, strain2] <- 1
  
  # compute Faith's PD
  pd_val <- pd(comm, tree, include.root = F)
  
  return(pd_val$PD)
}

strain_combos <- t(combn(as.character(Bifido_strain_order), 2))

pairwise_pd_df <- data.frame(
  strain1 = strain_combos[,1],
  strain2 = strain_combos[,2],
  PD = apply(strain_combos, 1, function(x) 
    faith_pd_pair(tree, x[1], x[2]))
)

# Write output
write.table(pairwise_pd_df, file = output_file, sep = ",", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)