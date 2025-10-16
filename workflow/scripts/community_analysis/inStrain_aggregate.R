library(argparse)
library(dplyr)
library(tidyverse)

parser <- ArgumentParser(description = "Aggregate inStrain results with dRep and metadata")

parser$add_argument("-i", "--input.dir",
                    nargs = "+",      # accept multiple values (space-separated list)
                    required = TRUE,
                    help = "List of inStrain output directories")
parser$add_argument("-d", "--drep",
                    required = TRUE,
                    help = "dRep clustering file (e.g., genome_info.tsv)")
parser$add_argument("-m", "--metadata",
                    required = TRUE,
                    help = "Sample metadata file")
parser$add_argument("-g", "--gen_func",
                    required = TRUE,
                    help = "script with general functions used in the script")
parser$add_argument("-o", "--output",
                    required = TRUE,
                    help = "Output table path")


# Parse arguments
args <- parser$parse_args()
source(args$gen_func)

# create output directory if it does not exist
if (!dir.exists(args$output)) {
    dir.create(args$output, recursive = TRUE)
}

# recover the files "output/genome_info.tsv" in each input directory
print("recovering genome_info.tsv files from input directories")
input_dirs <- args$input.dir
print(paste("Input directories:", paste(input_dirs, collapse = ",")))

input_files <- lapply(input_dirs, function(dir) {
    sample<- strsplit(basename(dir), split = "_")[[1]][1]
    file.path(dir, "output", paste0(sample, "_profile_genome_info.tsv"))
})

find_vOTU_merges_by_cooccurrence <- function(filtered_data,
                                             cor_threshold = 0.9,
                                             dist_threshold_factor = 1,
                                             min_shared = 3,
                                             cooccur_threshold = 0.7) {
  message("Detecting co-occurring genomes to merge vOTUs...")
  
  # ---- Prepare coverage matrix ----
  cov_mat <- filtered_data %>%
    select(sample, genome, coverage) %>%
    pivot_wider(names_from = genome, values_from = coverage, values_fill = 0) %>%
    column_to_rownames("sample")
  
  genomes <- colnames(cov_mat)
  if (length(genomes) < 2) {
    message("Not enough genomes to compare — skipping merging.")
    return(data.frame(genome = character(), old_vOTU = character(), new_vOTU = character()))
  }
  
  # ---- Compute pairwise metrics ----
  cor_mat <- cor(cov_mat, method = "pearson", use = "pairwise.complete.obs")
  dist_mat <- as.matrix(dist(t(cov_mat), method = "euclidean"))
  threshold_dist <- dist_threshold_factor * mean(dist_mat)
  
  shared_nonzero <- function(a, b) sum(a > 0 & b > 0)
  union_nonzero  <- function(a, b) sum(a > 0 | b > 0)
  
  shared_mat <- outer(genomes, genomes, Vectorize(function(g1, g2) shared_nonzero(cov_mat[, g1], cov_mat[, g2])))
  union_mat  <- outer(genomes, genomes, Vectorize(function(g1, g2) union_nonzero(cov_mat[, g1], cov_mat[, g2])))
  dimnames(shared_mat) <- list(genomes, genomes)
  dimnames(union_mat)  <- list(genomes, genomes)
  
  # ---- Identify valid pairs ----
  pairs <- expand.grid(genome1 = genomes, genome2 = genomes, stringsAsFactors = FALSE) %>%
    filter(genome1 < genome2) %>%
    mutate(
      cor = mapply(function(a, b) cor_mat[a, b], genome1, genome2),
      dist = mapply(function(a, b) dist_mat[a, b], genome1, genome2),
      shared = mapply(function(a, b) shared_mat[a, b], genome1, genome2),
      union  = mapply(function(a, b) union_mat[a, b], genome1, genome2),
      cooccur_ratio = shared / union
    ) %>%
    filter(
      cor > cor_threshold,
      dist < threshold_dist,
      shared >= min_shared,
      cooccur_ratio >= cooccur_threshold
    )
  
  if (nrow(pairs) == 0) {
    message("No genome pairs met co-occurrence criteria — no merges needed.")
    return(data.frame(genome = character(), old_vOTU = character(), new_vOTU = character()))
  }
  
  # ---- Find connected genome clusters (without igraph) ----
  genomes_all <- unique(c(pairs$genome1, pairs$genome2))
  
  adjacency <- lapply(genomes_all, function(g) {
    c(pairs$genome2[pairs$genome1 == g], pairs$genome1[pairs$genome2 == g])
  })
  names(adjacency) <- genomes_all
  
  visited <- setNames(rep(FALSE, length(genomes_all)), genomes_all)
  clusters <- list()
  cluster_id <- 0
  
  dfs <- function(node, cluster_idx) {
    visited[node] <<- TRUE
    clusters[[cluster_idx]] <<- c(clusters[[cluster_idx]], node)
    for (neighbor in adjacency[[node]]) {
      if (!visited[neighbor]) dfs(neighbor, cluster_idx)
    }
  }
  
  for (g in genomes_all) {
    if (!visited[g]) {
      cluster_id <- cluster_id + 1
      clusters[[cluster_id]] <- character()
      dfs(g, cluster_id)
    }
  }
  
  cluster_map <- do.call(rbind, lapply(seq_along(clusters), function(i) {
    data.frame(genome = clusters[[i]], cluster = i, stringsAsFactors = FALSE)
  }))
  
  # ---- Assign new vOTU per cluster ----
  vOTU_map <- filtered_data %>%
    select(genome, vOTU) %>%
    distinct()
  
  cluster_vOTU <- cluster_map %>%
    left_join(vOTU_map, by = "genome") %>%
    group_by(cluster) %>%
    mutate(new_vOTU = min(vOTU, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(vOTU != new_vOTU) %>%     # only genomes that change
    select(genome, old_vOTU = vOTU, new_vOTU)
  
  message("vOTU merge mapping completed.")
  return(cluster_vOTU)
}

comm_data<- do.call(bind_rows, lapply(input_files, function(file) {
    if (file.exists(file)) {
        # if a file that ends with ".done" (need to find the extension) is in the same directory as "file", return a empty dataframe with only the sample name in one row
        dir <- dirname(file)
        done_file <- list.files(dir, pattern = "\\.done$", full.names = TRUE)
        if (length(done_file) > 0) {
            sample_name <- strsplit(basename(file), split = "_")[[1]][1]
            return(data.frame(sample = sample_name))
        } else {
            sample_name <- strsplit(basename(file), split = "_")[[1]][1]
            df<- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
            mutate(sample = sample_name) %>%
            select(sample, everything())

            return(df)
        }
    } else {
        warning(sprintf("File %s does not exist.", file))
        NULL
    }
})) 


# Read dRep clustering file
drep_file <- args$drep
drep<- read.table(drep_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quality_threshold <- c("Complete",  "High-quality", "Medium-quality")

quality_threshold_quantitative <- c("Complete"=4,  "High-quality"=3, "Medium-quality"=2, "Low-quality"=1)
num_to_quality <- setNames(names(quality_threshold_quantitative),
                           as.character(unname(quality_threshold_quantitative)))
# VOTUs to keep based on quality
vOTU_qual_keep<- drep %>% filter(representative==T, checkv_quality %in% quality_threshold) %>%
  pull(vOTU)

# Read metadata file
metadata_file <- args$metadata
metadata <- read.table(metadata_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# function to format community data
format_community_data <- function(df, drep_data, mtdata) {
    qual_genome<- setNames(drep_data$checkv_quality, drep_data$genome)
    genome_vOTU<- setNames(drep_data$vOTU, drep_data$genome)

    df_metadata<- df %>%
    mutate(genome=gsub(".fasta", "", genome),
           Exp=get_mtadata(mtdata, "Exp", "SampleName", sample),
           Batch= get_mtadata(mtdata, "Batch", "SampleName", sample),
           Passage= get_mtadata(mtdata, "Passage", "SampleName", sample),
           Media= get_mtadata(mtdata, "Media", "SampleName", sample),
           Replicate= get_mtadata(mtdata, "Replicate", "SampleName", sample),
           Cocktail= get_mtadata(mtdata, "Cocktail", "SampleName", sample),
           checkv_quality= as.character(qual_genome[genome]),
           vOTU= genome_vOTU[genome]) %>%
    relocate(sample, genome, vOTU, checkv_quality, Media, Passage, Cocktail) %>%
    relocate(length, .after = genome) %>%
    arrange(Exp, sample, vOTU)

    return(df_metadata)
}

# Format the community data
print("Formatting community data")
formatted_data <- format_community_data(comm_data, drep, metadata)

print(head(formatted_data))

# Write the formatted data to the output file
print(paste("Writing formatted data to", args$output))
output_file <- file.path(args$output, "vircom_data.csv")
write.csv(formatted_data, output_file, row.names = FALSE, quote = FALSE)

# filter data
print("Filtering data for vOTUs with breadth >= 0.7 and quality")
filtered_data<- formatted_data %>%
    filter(breadth>=0.7, vOTU %in% vOTU_qual_keep) %>%
    group_by(sample) %>%
    mutate(rel_ab= coverage / sum(coverage, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Exp, sample, vOTU, rel_ab)

print(colnames(filtered_data))

merge_map <- find_vOTU_merges_by_cooccurrence(filtered_data)
merge_map_named_vector <- setNames(merge_map$new_vOTU, merge_map$old_vOTU)

filtered_data_updated <- filtered_data %>%
  mutate(vOTU = ifelse(
    vOTU %in% merge_map$old_vOTU,
    merge_map_named_vector[vOTU],  # replace with new vOTU
    vOTU                           # keep original if not merged
  ))
    
drep_updated <- drep %>%
  mutate(vOTU = ifelse(
    vOTU %in% merge_map$old_vOTU,
    merge_map_named_vector[vOTU],  # replace with new vOTU
    vOTU                           # keep original if not merged
  ))
  
print("Aggregating data by sample and vOTU")
filtered_data_vOTU<- filtered_data_updated %>%
    group_by(Exp, sample, Media,  Cocktail, Passage, Batch, Replicate, vOTU) %>%
    reframe(rel_ab= sum(rel_ab, na.rm = TRUE),
            nucl_diversity= mean(nucl_diversity, na.rm = TRUE),
            mean_breadth= mean(breadth, na.rm = TRUE),
            max_length= max(length, na.rm = TRUE),
            checkv_quality=max(quality_threshold_quantitative[checkv_quality]),
            checkv_quality= num_to_quality[as.character(checkv_quality)]) %>%
    ungroup() 
print(colnames(filtered_data_vOTU))

# Write the filtered data to the output file
print(paste("Writing filtered data to", args$output))
output_filtered_file <- file.path(args$output, "vircom_data_filtered.csv")
write.csv(filtered_data, output_filtered_file, row.names = FALSE, quote = FALSE)

output_filtered_file_vOTU <- file.path(args$output, "vircom_data_filtered_vOTU.csv")
write.csv(filtered_data_vOTU, output_filtered_file_vOTU, row.names = FALSE, quote = FALSE)

print("Writing updated dRep file to output directory")
output_drep_file <- file.path(args$output, "dRep_summary_average_updated.tsv")
write.table(drep_updated, output_drep_file, sep = "\t", row.names = FALSE, quote = FALSE)

print("writing vOTU merge map to output directory")
output_merge_map_file <- file.path(args$output, "vOTU_merge_map.tsv")
write.table(merge_map, output_merge_map_file, sep = "\t", row.names = FALSE, quote = FALSE)

print("inStrain aggregation completed successfully.")