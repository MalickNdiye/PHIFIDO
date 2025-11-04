#!/usr/bin/env Rscript

# ============================================================
# Aggregate inStrain Results with dRep and Metadata
# ============================================================
# This script:
#   1. Aggregates per-sample inStrain genome profiles
#   2. Integrates them with dRep clustering and metadata
#   3. Detects and merges co-occurring vOTUs based on coverage patterns
#   4. Produces summary tables and diagnostic plots
# ============================================================

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(vegan)
  library(purrr)
  library(tibble)
})

# ============================================================
# Argument Parsing
# ============================================================

parser <- ArgumentParser(description = "Aggregate inStrain results with dRep and metadata")

parser$add_argument("-i", "--input.dir", nargs = "+", required = TRUE,
                    help = "List of inStrain output directories")
parser$add_argument("-d", "--drep", required = TRUE,
                    help = "dRep clustering file (e.g., genome_info.tsv)")
parser$add_argument("-m", "--metadata", required = TRUE,
                    help = "Sample metadata file")
parser$add_argument("-g", "--gen_func", required = TRUE,
                    help = "Script with general helper functions used here")
parser$add_argument("-a", "--ani", required = TRUE,
                    help = "Pairwise ANI comparison file")
parser$add_argument("-o", "--output", required = TRUE,
                    help = "Output directory path")

args <- parser$parse_args()
source(args$gen_func)

# Create output directory
if (!dir.exists(args$output)) dir.create(args$output, recursive = TRUE)

# ============================================================
# Input Handling
# ============================================================

message("Recovering genome_info.tsv files from input directories...")

input_dirs <- args$input.dir
input_files <- lapply(input_dirs, function(dir) {
  sample <- strsplit(basename(dir), "_")[[1]][1]
  file.path(dir, "output", paste0(sample, "_profile_genome_info.tsv"))
})

# ============================================================
# Helper Functions
# ============================================================

# -----------------------------
# Coverage Matrix Preparation
# -----------------------------
prepare_cov_matrix <- function(data) {
  data %>%
    select(sample, genome, rel_ab) %>%
    rename(coverage = rel_ab) %>%
    ungroup() %>%
    pivot_wider(names_from = genome, values_from = coverage, values_fill = 0) %>%
    column_to_rownames("sample")
}

# -----------------------------
# Pairwise Stats (correlation, distance, co-occurrence)
# -----------------------------
compute_pairwise_stats <- function(cov_mat) {
  genomes <- colnames(cov_mat)
  if (length(genomes) < 2) return(NULL)

  cor_mat  <- cor(cov_mat, method = "spearman", use = "pairwise.complete.obs")
  dist_mat <- as.matrix(vegdist(t(cov_mat), method = "bray"))

  shared_nonzero <- function(a, b) sum(a > 0 & b > 0)
  union_nonzero  <- function(a, b) sum(a > 0 | b > 0)

  shared_mat <- outer(genomes, genomes, Vectorize(function(g1, g2)
    shared_nonzero(cov_mat[, g1], cov_mat[, g2])))
  union_mat  <- outer(genomes, genomes, Vectorize(function(g1, g2)
    union_nonzero(cov_mat[, g1], cov_mat[, g2])))

  dimnames(shared_mat) <- list(genomes, genomes)
  dimnames(union_mat)  <- list(genomes, genomes)

  expand.grid(genome1 = genomes, genome2 = genomes, stringsAsFactors = FALSE) %>%
    filter(genome1 < genome2) %>%
    mutate(
      cor = mapply(function(a, b) cor_mat[a, b], genome1, genome2),
      dist = mapply(function(a, b) dist_mat[a, b], genome1, genome2),
      shared = mapply(function(a, b) shared_mat[a, b], genome1, genome2),
      union  = mapply(function(a, b) union_mat[a, b], genome1, genome2),
      cooccur_ratio = shared / union
    )
}

# -----------------------------
# Connected Component Clustering (DFS)
# -----------------------------
find_clusters <- function(pairs) {
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

  do.call(rbind, lapply(seq_along(clusters), function(i) {
    data.frame(genome = clusters[[i]], cluster = i, stringsAsFactors = FALSE)
  }))
}

# -----------------------------
# Map clusters to new vOTUs
# -----------------------------
map_clusters_to_vOTUs <- function(cluster_map, filtered_data) {
  vOTU_map <- filtered_data %>%
    select(genome, vOTU) %>%
    distinct()

  cluster_map %>%
    left_join(vOTU_map, by = "genome") %>%
    group_by(cluster) %>%
    mutate(new_vOTU = min(vOTU, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(vOTU != new_vOTU) %>%
    select(genome, old_vOTU = vOTU, new_vOTU)
}

# -----------------------------
# Diagnostic Histograms
# -----------------------------
plot_metric_histograms <- function(pairs, cor_threshold, dist_threshold_factor,
                                   cooccur_threshold, outdir) {
                                    
  plotdir<- file.path(outdir, "plots")
  if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)
  mean_dist <- mean(pairs$dist, na.rm = TRUE)
  dist_thresh <- dist_threshold_factor * mean_dist

  # Plot helper
  make_plot <- function(x, vline, title, xlabel) {
    ggplot(pairs, aes(x = {{ x }})) +
      geom_histogram(bins = 50, fill = "grey70", color = "black") +
      geom_vline(xintercept = vline, color = "red", lwd = 1) +
      labs(title = title, x = xlabel, y = "Count")
  }

  # make scatter plot with x=dist, y=cor, color=cooccur_ratio
  # color above threshold in red
  p_scatter <- ggplot(pairs, aes(x = dist, y = cor, color = cooccur_ratio, size = pmin(union, 5))) +
    geom_point(alpha=0.7) +
    scale_color_gradient(
        low = "blue", 
        high = "green",
        limits = c(0, cooccur_threshold), # Gradient ends at the threshold
        oob = scales::censor,            # Censor values above it
        na.value = "red"                 # Color the censored values red
    ) +
    geom_hline(yintercept = cor_threshold, color = "red", lwd = 1) +
    geom_vline(xintercept = dist_thresh, color = "red", lwd = 1) +
    labs(title = "Pairwise correlation vs distance",
         x = "Bray–Curtis distance",
         y = "Spearman correlation",
         color = "Co-occurrence ratio") +
    theme_minimal()

  ggsave(file.path(plotdir, "correlation_histogram.png"),
         make_plot(pairs$cor, cor_threshold, "Pairwise correlations", "Spearman correlation"))
  ggsave(file.path(plotdir, "distance_histogram.png"),
         make_plot(pairs$dist, dist_thresh, "Pairwise distances", "Bray–Curtis distance"))
  ggsave(file.path(plotdir, "cooccurrence_histogram.png"),
         make_plot(pairs$cooccur_ratio, cooccur_threshold, "Pairwise co-occurrence ratios", "Shared / Union"))
  ggsave(file.path(plotdir, "correlation_distance_scatter.png"), p_scatter) 
}

# ============================================================
# Main Function: vOTU Merge Detection by Co-occurrence
# ============================================================

find_vOTU_merges_by_cooccurrence <- function(filtered_data,
                                             cor_threshold = 0.95,
                                             dist_threshold_factor = 0.2,
                                             min_shared = 5,
                                             cooccur_threshold = 0.9,
                                             outdir) {
  message("Detecting co-occurring genomes to merge vOTUs...")

  cov_mat <- prepare_cov_matrix(filtered_data)
  pairs <- compute_pairwise_stats(cov_mat)

  if (is.null(pairs) || nrow(pairs) == 0) {
    message("Not enough genomes to compare — skipping merging.")
    return(data.frame(genome = character(), old_vOTU = character(), new_vOTU = character()))
  }

  plot_metric_histograms(pairs, cor_threshold, dist_threshold_factor, cooccur_threshold, outdir)
  threshold_dist <- dist_threshold_factor * mean(pairs$dist, na.rm = TRUE)

  pairs_filtered <- pairs %>%
    filter(
      cor > cor_threshold,
      dist < threshold_dist,
      shared >= min_shared,
      cooccur_ratio >= cooccur_threshold
    )

  if (nrow(pairs_filtered) == 0) {
    message("No genome pairs met co-occurrence criteria — no merges needed.")
    return(data.frame(genome = character(), old_vOTU = character(), new_vOTU = character()))
  }

  cluster_map <- find_clusters(pairs_filtered)
  cluster_vOTU <- map_clusters_to_vOTUs(cluster_map, filtered_data)

  message("vOTU merge mapping completed. Diagnostic plots saved in: ", outdir)
  return(cluster_vOTU)
}

# ============================================================
# Plot vOTU Coverage
# ============================================================

plot_vOTU_coverage <- function(filtered_data, merge_map, outdir) {
  message("Plotting coverage across samples for merged vOTUs...")
  plotdir <- file.path(outdir, "plots")
  if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

  vOTU_to_genome <- filtered_data %>% select(genome, vOTU) %>% distinct()

  cluster_genomes <- merge_map %>%
    left_join(vOTU_to_genome, by = c("new_vOTU" = "vOTU")) %>%
    group_by(new_vOTU) %>%
    summarise(
      genomes = list(unique(c(genome.x, genome.y))),
      representative_genome = unique(genome.y),
      .groups = "drop"
    )

  pdf(file.path(plotdir, "merged_vOTU_coverage_plots.pdf"), width = 9, height = 5)

  for (i in seq_len(nrow(cluster_genomes))) {
    new_v <- cluster_genomes$new_vOTU[i]
    genomes_in_cluster <- unique(na.omit(cluster_genomes$genomes[[i]]))
    representative <- cluster_genomes$representative_genome[i]

    df_cluster <- filtered_data %>% filter(genome %in% genomes_in_cluster)
    if (nrow(df_cluster) == 0) next

    df_cluster <- df_cluster %>%
      mutate(
        is_representative = genome == representative,
        sample = factor(sample, levels = unique(sample))
      )

    genome_colors <- setNames(
      scales::hue_pal()(length(unique(df_cluster$genome))),
      unique(df_cluster$genome)
    )

    p <- ggplot(df_cluster, aes(x = sample, y = rel_ab, group = genome, color = genome)) +
      geom_line(aes(size = is_representative)) +
      geom_point(aes(size = is_representative)) +
      scale_color_manual(values = genome_colors) +
      scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.8), guide = "none") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(paste0("Relative Abundance across samples — merged vOTU: ", new_v,
                     "\nRepresentative genome: ", representative)) +
      xlab("Sample") +
      ylab("Relative Abundance")

    print(p)
  }

  dev.off()
  message("Saved merged vOTU coverage plots to: ",
          file.path(plotdir, "merged_vOTU_coverage_plots.pdf"))
}

# ============================================================
# Main Workflow
# ============================================================

# ---- Read Input Data ----
comm_data <- do.call(bind_rows, lapply(input_files, function(file) {
  if (file.exists(file)) {
    dir <- dirname(file)
    done_file <- list.files(dir, pattern = "\\.done$", full.names = TRUE)
    sample_name <- strsplit(basename(file), "_")[[1]][1]
    if (length(done_file) > 0) {
      return(data.frame(sample = sample_name))
    } else {
      read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
        mutate(sample = sample_name) %>%
        select(sample, everything())
    }
  } else {
    warning(sprintf("File %s does not exist.", file))
    NULL
  }
}))

# ---- Load dRep & Metadata ----
drep <- read.table(args$drep, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- read.table(args$metadata, header = TRUE, sep = ",", stringsAsFactors = FALSE)

quality_threshold <- c("Complete", "High-quality", "Medium-quality")
quality_threshold_quantitative <- c("Complete" = 4, "High-quality" = 3, "Medium-quality" = 2, "Low-quality" = 1)
num_to_quality <- setNames(names(quality_threshold_quantitative), as.character(unname(quality_threshold_quantitative)))

vOTU_qual_keep <- drep %>%
  filter(representative == TRUE, checkv_quality %in% quality_threshold) %>%
  pull(vOTU)

rep_genomes_vector <- setNames(drep$genome[drep$representative], drep$vOTU[drep$representative])

# ---- Format Community Data ----
format_community_data <- function(df, drep_data, mtdata) {
  qual_genome <- setNames(drep_data$checkv_quality, drep_data$genome)
  genome_vOTU <- setNames(drep_data$vOTU, drep_data$genome)

  df %>%
    mutate(
      genome = gsub(".fasta", "", genome),
      Exp = get_mtadata(mtdata, "Exp", "SampleName", sample),
      Batch = get_mtadata(mtdata, "Batch", "SampleName", sample),
      Passage = get_mtadata(mtdata, "Passage", "SampleName", sample),
      Media = get_mtadata(mtdata, "Media", "SampleName", sample),
      Replicate = get_mtadata(mtdata, "Replicate", "SampleName", sample),
      Cocktail = get_mtadata(mtdata, "Cocktail", "SampleName", sample),
      checkv_quality = as.character(qual_genome[genome]),
      vOTU = genome_vOTU[genome]
    ) %>%
    relocate(sample, genome, vOTU, checkv_quality, Media, Passage, Cocktail) %>%
    relocate(length, .after = genome) %>%
    arrange(Exp, sample, vOTU)
}

formatted_data <- format_community_data(comm_data, drep, metadata)
write.csv(formatted_data, file.path(args$output, "vircom_data.csv"), row.names = FALSE, quote = FALSE)

# ---- Filter and Merge ----
filtered_data <- formatted_data %>%
  filter(round(breadth, 1) >= 0.7, vOTU %in% vOTU_qual_keep) %>%
  group_by(sample) %>%
  mutate(rel_ab = coverage / sum(coverage, na.rm = TRUE)) %>%
  ungroup()

merge_map <- find_vOTU_merges_by_cooccurrence(filtered_data, outdir = args$output)

merge_map$merged_to <- rep_genomes_vector[merge_map$new_vOTU]
merge_map_named_vector <- setNames(merge_map$new_vOTU, merge_map$old_vOTU)


# add ANI info to merge_map
ani.data<- read.table(args$ani, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(ani.data)<- c("genome1", "genome2", "ANI", "Aligned_fragments", "total_fragments") 

ani.data<- ani.data %>%
  mutate(genome1 = gsub(".fasta", "", genome1),
         genome2 = gsub(".fasta", "", genome2),
         genome1 = gsub("../results/assembly/viral/single_genomes/", "", genome1),
         genome2 = gsub("../results/assembly/viral/single_genomes/", "", genome2),
         AF=Aligned_fragments/total_fragments)%>% 
  select(genome1, genome2, ANI, AF) 

merge_map_ani <- merge_map %>%
  left_join(ani.data, by = c("genome" = "genome1", "merged_to" = "genome2")) %>%
  rename(ANI_fw = ANI, AF_fw = AF) %>%
  left_join(ani.data, by = c("genome" = "genome2", "merged_to" = "genome1")) %>%
  rename(ANI_rev = ANI, AF_rev = AF)  

plot_vOTU_coverage(filtered_data, merge_map, outdir = args$output)

# ---- Update vOTU Assignments ----
filtered_data_updated <- filtered_data %>%
  mutate(vOTU = ifelse(vOTU %in% merge_map$old_vOTU, merge_map_named_vector[vOTU], vOTU))

drep_updated <- drep %>%
  mutate(vOTU = ifelse(vOTU %in% merge_map$old_vOTU, merge_map_named_vector[vOTU], vOTU))

# ---- Aggregate and Write Outputs ----
filtered_data_vOTU <- filtered_data_updated %>%
  group_by(Exp, sample, Media, Cocktail, Passage, Batch, Replicate, vOTU) %>%
  reframe(
    rel_ab = sum(rel_ab, na.rm = TRUE),
    nucl_diversity = mean(nucl_diversity, na.rm = TRUE),
    mean_breadth = mean(breadth, na.rm = TRUE),
    max_length = max(length, na.rm = TRUE),
    checkv_quality = max(quality_threshold_quantitative[checkv_quality]),
    checkv_quality = num_to_quality[as.character(checkv_quality)]
  ) %>%
  ungroup()

write.csv(filtered_data, file.path(args$output, "vircom_data_filtered.csv"), row.names = FALSE, quote = FALSE)
write.csv(filtered_data_updated, file.path(args$output, "vircom_data_filtered_updated.csv"), row.names = FALSE, quote = FALSE)
write.csv(filtered_data_vOTU, file.path(args$output, "vircom_data_filtered_vOTU.csv"), row.names = FALSE, quote = FALSE)
write.table(drep_updated, file.path(args$output, "dRep_summary_average_updated.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(merge_map_ani, file.path(args$output, "vOTU_merge_map.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("inStrain aggregation completed successfully.")
