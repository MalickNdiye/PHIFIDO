# Load libraries
library(data.table)
library(tidyverse)
library(ape)
library(optparse)

# set options
options(stringsAsFactors = FALSE)
option_list <- list(
  make_option(c("-m", "--matrix"), type = "character", default = NULL,
              help = "Path to VIRIDIC similarity matrix file", metavar = "character"),
    make_option(c("-v", "--vircom"), type = "character", default = NULL,
              help = "Path to AmpliPhage community vOTU file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "../output/",
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

output_dir <- opt$output
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# 1. Open similarity matrix
# -----------------------------
message("Loading VIRIDIC similarity matrix...")
sim_mat_p <- opt$matrix
sim_mat <- fread(sim_mat_p) %>%
  column_to_rownames("genome") %>%
  as.matrix()

# -----------------------------
# 2. Load AmpliPhage community data
# -----------------------------
message("Loading AmpliPhage community data...")
ampliphage_vircom <- fread(opt$vircom) %>%
  filter(
    Exp == "AmpliPhage",
    Cocktail != "none",
    !grepl("Pc", sample)
  )

vOTUs <- unique(ampliphage_vircom$vOTU)
target_genomes <- unique(ampliphage_vircom$genome)

gen_to_vOTU <- ampliphage_vircom %>%
  select(genome, vOTU) %>%
  distinct() %>%
  deframe()

# -----------------------------
# 3. Filter similarity matrix
# -----------------------------$
message("Filtering similarity matrix for target genomes...")
sim_mat_filt <- sim_mat[target_genomes, target_genomes]
colnames(sim_mat_filt) <- gen_to_vOTU[colnames(sim_mat_filt)]
rownames(sim_mat_filt) <- gen_to_vOTU[rownames(sim_mat_filt)]

# -----------------------------
# 4. Transform similarity into dissimilarity matrix, aggregate vOTUs
# -----------------------------
if (max(sim_mat_filt, na.rm = TRUE) > 1) sim_mat_filt <- sim_mat_filt / 100
diag(sim_mat_filt) <- 1
Dmat <- 1 - sim_mat_filt   # Distance: 0 = identical, 1 = fully different

long_df <- data.frame(
  vOTU_1 = rep(rownames(Dmat), times = ncol(Dmat)),
  vOTU_2 = rep(colnames(Dmat), each = nrow(Dmat)),
  distance = as.vector(Dmat)
)
Agg_long_df <- aggregate(distance ~ vOTU_1 + vOTU_2, data = long_df, FUN = mean)
Agg_Dmat <- xtabs(distance ~ vOTU_1 + vOTU_2, data = Agg_long_df)

# save aggregated distance matrix in output directory
write.csv(Agg_Dmat, file = file.path(output_dir, "vOTU_distance_matrix.csv"), row.names = TRUE)

D <- as.dist(Agg_Dmat)

# -----------------------------
# 6. Perform UPGMA clustering
# -----------------------------
hc_upgma <- hclust(D, method = "average")

# Convert to phylo object and ladderize
tree_upgma <- as.phylo(hc_upgma)
tree_upgma <- ladderize(tree_upgma)

# -----------------------------
# 8. Save tree in Newick format
# -----------------------------
write.tree(tree_upgma, file = file.path(output_dir, "vOTU_upgma_tree.nwk", sep=""))
