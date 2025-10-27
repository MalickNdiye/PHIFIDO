#!/usr/bin/env Rscript
library("optparse")
library(data.table)
library(dplyr)
library(igraph)

# -----------------------------
# Options
# -----------------------------
option_list = list(
  make_option(c("-i", "--vOTU_table"), type="character", default=NULL,
              help="Parsed vOTU table", metavar="character"),
  make_option(c("-p", "--problematic"), type="character", default=NULL,
              help="Problematic vOTU pairs", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="Output updated vOTU table", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# -----------------------------
# 1. Read vOTU table
# -----------------------------
cat("--- Reading vOTU table\n")
vOTU_table <- fread(opt$vOTU_table, data.table = FALSE)

# Ensure qual_score numeric
if(!"qual_score" %in% colnames(vOTU_table)){
  qual_score_map <- c("Low-quality"=1, "Medium-quality"=2, "High-quality"=3, "Complete"=4)
  vOTU_table$qual_score <- qual_score_map[vOTU_table$checkv_quality]
  vOTU_table$qual_score[is.na(vOTU_table$qual_score)] <- 0
}

# -----------------------------
# 2. Read problematic genome pairs
# -----------------------------
cat("--- Reading problematic vOTU pairs\n")
problem_pairs <- fread(opt$problematic, data.table = FALSE)

# -----------------------------
# 3. Build similarity graph
# -----------------------------
cat("--- Building similarity graph\n")
edges <- unique(problem_pairs %>% select(genome1, genome2))
g <- graph_from_data_frame(edges, directed = FALSE)

cat("--- Finding connected components\n")
components <- components(g)
membership <- components$membership

# Assign a new vOTU ID for each component
new_vOTU_map <- tibble(
  genome = names(membership),
  new_vOTU = paste0("vOTU_", membership)
)

# -----------------------------
# 4. Merge new vOTU IDs into original table
# -----------------------------
cat("--- Merging with original vOTU table\n")
vOTU_updated <- vOTU_table %>%
  left_join(new_vOTU_map, by="genome") %>%
  mutate(vOTU = ifelse(!is.na(new_vOTU), new_vOTU, vOTU)) %>%
  select(-new_vOTU)

# -----------------------------
# 5. Choose representatives per vOTU
# -----------------------------
cat("--- Choosing representatives for each vOTU\n")
vOTU_updated <- vOTU_updated %>%
  group_by(vOTU) %>%
  arrange(desc(qual_score), desc(length), .by_group=TRUE) %>%
  mutate(representative = row_number() == 1) %>%
  ungroup()

# -----------------------------
# 6. Select and reorder columns
# -----------------------------
vOTU_updated <- vOTU_updated %>%
  select(genome, length, checkv_quality, primary_cluster, secondary_cluster, vOTU, representative)

# -----------------------------
# 7. Write updated table
# -----------------------------
cat("--- Writing updated vOTU table\n")
fwrite(vOTU_updated, opt$outfile, sep="\t")
cat("--- Done! Updated vOTU table written to", opt$outfile, "\n")
