#!/usr/bin/env Rscript

for (pkg in c("optparse", "extrafont")) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://stat.ethz.ch/CRAN/")
    library(pkg, character.only = TRUE)
  }
}
library(vegan)
library(tidyverse)
library(data.table)
library(ggplot2)

if (!require(phyloseq)) {
  if (!require(BiocManager)) install.packages("BiocManager", repos = "https://stat.ethz.ch/CRAN/")
  source("http://bioconductor.org/biocLite.R")
  biocLite("phyloseq")
}

FONT_FAMILY <- tryCatch({
  extrafont::loadfonts(quiet = TRUE)
  if ("Arial" %in% extrafont::fonts()) "Arial" else "sans"
}, error = function(e) "sans")

# --------------------------
# 1. Command-Line Arguments
# --------------------------
option_list <- list(
  make_option(c("-a", "--asv_table"),      type = "character", default = "processed_data/ASV_table_merged_new.rds"),
  make_option(c("-t", "--taxonomy"),       type = "character", default = "processed_data/ASV_Taxonomy_sp_new.RDS"),
  make_option(c("-m", "--metadata"),       type = "character", default = "../data/sample_metadata.csv"),
  make_option(c("-e", "--experiment"),     type = "character", default = "Bicom6"),
  make_option(c("-o", "--output_dir"),     type = "character", default = "1_Output/"),
  make_option(c("-p", "--prev_thresh"),    type = "numeric",   default = 0.1),
  make_option(c("-r", "--relabund_thresh"),type = "numeric",   default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

ASV_samples_table_noChim_p <- opt$asv_table
ASV_taxonomy_sp_p           <- opt$taxonomy
metadata_p                  <- opt$metadata
experiment                  <- opt$experiment
output_dir                  <- opt$output_dir
prev_thresh                 <- opt$prev_thresh
relabund_thresh             <- opt$relabund_thresh

tables_dir  <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

bicom_compositions <- list(
  Bicom1 = c("ESL0820", "ESL0822", "ESL0198", "ESL0199", "ESL0200", "ESL0819"),
  Bicom2 = c("ESL0820", "ESL0822", "ESL0198", "ESL0824", "ESL0827", "ESL1060"),
  Bicom3 = c("ESL0199", "ESL0200", "ESL0819", "ESL0824", "ESL0827", "ESL1069"),
  Bicom4 = c("ESL0820", "ESL0199", "ESL0819", "ESL1060", "ESL1069", "ESL1077"),
  Bicom5 = c("ESL0822", "ESL0200", "ESL0824", "ESL1060", "ESL1069", "ESL1077"),
  Bicom6 = c("ESL0198", "ESL0827", "ESL1060", "ESL1069", "ESL1077", "ESL0819")
)

# Figure dimensions: A4 portrait (210 mm wide), max 1/3 page height (99 mm)
FIG_W <- 210
FIG_H <- 99

base_fig_theme <- theme_minimal() +
  theme(
    text         = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    axis.text    = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    axis.title   = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    plot.title   = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    legend.text  = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    legend.title = element_text(family = FONT_FAMILY, size = 12, face = "bold"),
    strip.text   = element_text(family = FONT_FAMILY, size = 12, face = "bold")
  )

# --------------------------
# 2. Build Phyloseq
# --------------------------
ASV_samples_table_noChim <- readRDS(ASV_samples_table_noChim_p)
rownames(ASV_samples_table_noChim) <- gsub(".hifi_reads", "", rownames(ASV_samples_table_noChim))
ASV_taxonomy_sp <- readRDS(ASV_taxonomy_sp_p)

metadata <- read.csv(metadata_p, sep = ",")
samples_to_remove <- metadata$SampleID[metadata$Exp != experiment]
metadata <- metadata %>% column_to_rownames("SampleID")

ps <- phyloseq(otu_table(ASV_samples_table_noChim, taxa_are_rows = FALSE),
               sample_data(metadata), tax_table(ASV_taxonomy_sp))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps <- prune_samples(!sample_names(ps) %in% samples_to_remove, ps)

tmp <- as.data.frame(sample_data(ps))
tmp$original_lib_size <- sample_sums(ps)
sample_data(ps) <- sample_data(tmp)

# --------------------------
# 3. Prevalence & Abundance Filtering
# --------------------------
asv_stats <- psmelt(ps) %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  dplyr::rename(ASV = OTU) %>%
  filter(!grepl("Blank", Sample)) %>%
  group_by(ASV) %>%
  summarise(
    genus    = unique(Genus),
    species  = unique(Species),
    ASV_type = case_when(
      grepl("ESL", species)                        ~ "BiCom",
      genus == "Bifidobacterium" & is.na(species)  ~ "Unassigned Bifidobacterium ASV",
      is.na(genus)                                 ~ "Genus Unassigned",
      TRUE                                         ~ "Other Genus"
    ),
    prevalence  = mean(Abundance > 0),
    max_rel     = max(rel_abund),
    median_rel  = median(rel_abund),
    .groups = "drop"
  ) %>%
  filter(!is.na(median_rel)) %>%
  arrange(median_rel)

plt <- asv_stats %>%
  mutate(median_rel_bin = cut(median_rel, breaks = c(0, 0.001, 0.01, 0.1, 1),
                              include.lowest = TRUE, dig.lab = 3)) %>%
  ggplot(aes(x = prevalence, y = max_rel, size = median_rel_bin, col = ASV_type, label = ASV)) +
  geom_text(data = subset(asv_stats, prevalence > 0.05), hjust = 1.5, size = 3) +
  geom_point() +
  geom_hline(yintercept = relabund_thresh, linetype = "dashed", col = "darkred") +
  geom_vline(xintercept = prev_thresh,     linetype = "dashed", col = "darkred") +
  labs(x = "Prevalence", y = "Maximum relative abundance", title = experiment,
       size = "Median rel abundance", col = "ASV origin") +
  scale_size_manual(values = c(2, 3, 4, 5, 6), drop = FALSE) +
  base_fig_theme

keep_asvs  <- asv_stats %>% filter(!(prevalence < prev_thresh & max_rel < relabund_thresh)) %>% pull(ASV)
ps_filtered <- prune_taxa(keep_asvs, ps)

# --------------------------
# 4. Decontamination Functions (library size only)
# --------------------------
prepare_metadata <- function(ps, group_col) {
  sd     <- sample_data(ps)
  sd_df  <- as.data.frame(sd)
  sd$Group <- apply(sd_df[, group_col, drop = FALSE], 1, paste, collapse = "_")
  sample_data(ps) <- sd
  return(ps)
}

# Returns Spearman correlation with library size for every species group present
find_via_correlation <- function(seqtab_sub, lib_size_vals) {
  seqtab_rel <- sweep(seqtab_sub, 1, pmax(rowSums(seqtab_sub), 1), "/")
  sp_cors <- apply(seqtab_rel, 2, function(rel_abund) {
    valid_idx <- which(rel_abund > 0 & !is.na(lib_size_vals) & !is.na(rel_abund))
    if (length(valid_idx) < 3 || sd(rel_abund[valid_idx], na.rm = TRUE) == 0) return(NA)
    lib_log <- log10(pmax(lib_size_vals[valid_idx], 1e-6))
    rel_log <- log10(rel_abund[valid_idx])
    if (sd(lib_log, na.rm = TRUE) == 0) return(NA)
    cor(rel_log, lib_log, method = "spearman", use = "complete.obs")
  })
  data.frame(Species_Group = names(sp_cors), cor = unname(sp_cors), row.names = NULL)
}

# Returns prevalence for every ASV present
find_via_prevalence <- function(seqtab_sub) {
  prev_calc <- colSums(seqtab_sub > 0) / nrow(seqtab_sub)
  data.frame(ASV = names(prev_calc), prevalence = unname(prev_calc), row.names = NULL)
}

identify_group_contaminants <- function(ps, cor_threshold, prev_threshold = 0.1) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("ASV") %>%
    mutate(Species_Group = case_when(
      !is.na(Species) & Species != "" & Species != "NA" & Species != "unassigned" ~ paste(Genus, Species, sep = "_"),
      TRUE ~ paste("Unclassified_ASV", ASV, sep = "_")
    ))

  agg_otu <- t(rowsum(t(otu_mat), group = tax_df$Species_Group[match(colnames(otu_mat), tax_df$ASV)]))
  sd_df   <- as.data.frame(sample_data(ps))

  list_of_contaminants <- lapply(unique(as.character(sample_data(ps)$Group)), function(grp) {
    if (is.na(grp)) return(NULL)
    samples_to_keep <- sample_names(ps)[sample_data(ps)$Group == grp]
    if (length(samples_to_keep) < 3) return(NULL)

    seqtab_sub_asv <- otu_mat[samples_to_keep, colSums(otu_mat[samples_to_keep, , drop = FALSE]) > 0, drop = FALSE]
    seqtab_sub_sp  <- agg_otu[samples_to_keep, colSums(agg_otu[samples_to_keep, , drop = FALSE]) > 0, drop = FALSE]
    if (ncol(seqtab_sub_asv) == 0 || ncol(seqtab_sub_sp) == 0) return(NULL)

    lib_size_vals <- as.numeric(unlist(sd_df$original_lib_size[match(samples_to_keep, rownames(sd_df))]))

    # Compute values for ALL ASVs / species groups
    all_cors  <- find_via_correlation(seqtab_sub_sp, lib_size_vals)
    all_prevs <- find_via_prevalence(seqtab_sub_asv)

    # Determine which are flagged
    contam_sp_groups <- all_cors$Species_Group[!is.na(all_cors$cor) & all_cors$cor < cor_threshold]
    contam_asv_prev  <- all_prevs$ASV[all_prevs$prevalence < prev_threshold]

    flagged_sp <- unique(c(contam_sp_groups, tax_df$Species_Group[tax_df$ASV %in% contam_asv_prev]))
    if (length(flagged_sp) == 0) return(NULL)

    # Build result: cor and prevalence values included for all flagged ASVs regardless of method
    tax_df %>%
      filter(Species_Group %in% flagged_sp) %>%
      select(ASV, Species_Group, Genus, Species) %>%
      mutate(Condition = grp) %>%
      left_join(all_cors,  by = "Species_Group") %>%
      left_join(all_prevs, by = "ASV") %>%
      mutate(Method = sub(";+$", "", paste0(
        ifelse(Species_Group %in% contam_sp_groups, "correlation;", ""),
        ifelse(ASV %in% contam_asv_prev,            "prevalence",   "")
      )))
  })

  bind_rows(list_of_contaminants) %>%
    select(Condition, ASV, Species_Group, Genus, Species, cor, prevalence, Method)
}

remove_targeted_contaminants <- function(ps, contam_df) {
  otu_mat      <- as(otu_table(ps), "matrix")
  is_taxa_rows <- taxa_are_rows(ps)
  if (is_taxa_rows) otu_mat <- t(otu_mat)

  for (i in seq_len(nrow(contam_df))) {
    target_samples <- intersect(sample_names(ps)[sample_data(ps)$Group == contam_df$Condition[i]], rownames(otu_mat))
    if (contam_df$ASV[i] %in% colnames(otu_mat) && length(target_samples) > 0)
      otu_mat[target_samples, contam_df$ASV[i]] <- 0
  }
  otu_table(ps) <- if (is_taxa_rows) otu_table(t(otu_mat), taxa_are_rows = TRUE) else otu_table(otu_mat, taxa_are_rows = FALSE)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

clean_microbiome_data <- function(ps, group_col = "Replicate", cor_threshold = -0.5, prev_threshold = 0.1,
                                  expected_comp = NULL, experiment_name = NULL) {
  cat("1. Formatting metadata...\n")
  ps <- prepare_metadata(ps, group_col)

  cat("2. Identifying contaminants using library size...\n")
  contam_df <- identify_group_contaminants(ps, cor_threshold, prev_threshold)

  if (!is.null(contam_df) && nrow(contam_df) > 0) {

    otu_mat       <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    rel_abund_mat <- sweep(otu_mat, 1, pmax(rowSums(otu_mat), 1), "/")
    samp_condition <- as.character(sample_data(ps)$Group)
    samp_names     <- rownames(otu_mat)

    contam_df <- contam_df %>%
      rowwise() %>%
      mutate(
        mean_rel_abund = {
          vs <- samp_names[samp_condition == Condition]
          if (length(vs) > 0 && ASV %in% colnames(rel_abund_mat)) mean(rel_abund_mat[vs, ASV], na.rm = TRUE) else NA_real_
        },
        median_rel_abund = {
          vs <- samp_names[samp_condition == Condition]
          if (length(vs) > 0 && ASV %in% colnames(rel_abund_mat)) median(rel_abund_mat[vs, ASV], na.rm = TRUE) else NA_real_
        }
      ) %>%
      ungroup()

    if (!is.null(expected_comp) && identical(experiment_name, "Bicom6")) {
      contam_df <- contam_df %>% rowwise() %>%
        mutate(FP = {
          pieces  <- strsplit(as.character(Condition), "_")[[1]]
          matched <- intersect(pieces, names(expected_comp))
          if (length(matched) > 0) Species %in% expected_comp[[matched[1]]] else FALSE
        }) %>% ungroup()
    } else {
      contam_df <- contam_df %>% mutate(FP = FALSE)
    }

    true_contam <- contam_df %>% filter(!FP)

    cat(sprintf("\n--- Decontamination Summary ---\nTotal Flagged: %d\nProtected: %d\nRemoved: %d\n\n",
                nrow(contam_df), nrow(contam_df) - nrow(true_contam), nrow(true_contam)))

    ps_clean <- if (nrow(true_contam) > 0) remove_targeted_contaminants(ps, true_contam) else ps
    return(list(clean_phyloseq = ps_clean, contaminants_removed = contam_df))
  }

  cat("\nNo contaminants found. Pipeline complete.\n")
  return(list(clean_phyloseq = ps, contaminants_removed = NULL))
}

plot_bicom_correlation_spectrum <- function(ps, expected_comp) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_rel <- sweep(otu_mat, 1, pmax(rowSums(otu_mat), 1), "/")

  sd_df  <- as.data.frame(sample_data(ps))
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE) %>% rownames_to_column("ASV")

  all_cor_data <- do.call(rbind, lapply(unique(as.character(sd_df$Group)), function(grp) {
    samps <- rownames(sd_df)[sd_df$Group == grp]
    if (length(samps) < 3) return(NULL)
    pieces    <- strsplit(grp, "_")[[1]]
    bicom_key <- intersect(pieces, names(expected_comp))
    if (length(bicom_key) == 0) return(NULL)

    allowed_species <- expected_comp[[bicom_key[1]]]
    lib_vals  <- as.numeric(unlist(sd_df[samps, "original_lib_size"]))
    lib_log   <- log10(pmax(lib_vals, 1e-6))
    sub_rel   <- otu_rel[samps, , drop = FALSE]
    present_asvs <- colnames(sub_rel)[colSums(sub_rel) > 0]

    do.call(rbind, lapply(present_asvs, function(asv_id) {
      x       <- as.numeric(sub_rel[, asv_id])
      v_idx   <- which(x > 0)
      sp_name <- as.character(tax_df$Species[tax_df$ASV == asv_id])[1]
      rho     <- NA
      if (length(v_idx) >= 3 && sd(x[v_idx]) > 0 && sd(lib_log[v_idx]) > 0)
        rho <- cor(log10(x[v_idx]), lib_log[v_idx], method = "spearman")
      data.frame(ASV = asv_id, Species = sp_name, Condition = grp, Correlation = rho,
                 Type = ifelse(sp_name %in% allowed_species, "Expected Member", "Unexpected (Contaminant)"),
                 stringsAsFactors = FALSE)
    }))
  }))

  plot_df <- all_cor_data %>% filter(!is.na(Correlation))

  plt <- ggplot(plot_df, aes(x = Correlation, fill = Type)) +
    geom_histogram(binwidth = 0.05, alpha = 0.7, color = "white", position = "identity") +
    scale_fill_manual(values = c("Expected Member" = "#4DAF4A", "Unexpected (Contaminant)" = "#E41A1C")) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    labs(x = "Spearman Rho", y = "ASV Count", fill = "Biological Status") +
    base_fig_theme + theme(legend.position = "bottom")

  return(list(plot = plt, data = plot_df))
}

# --------------------------
# 5. Execute Pipeline
# --------------------------
if (experiment == "Bicom6") {
  decontam_results <- clean_microbiome_data(ps_filtered, group_col = "Replicate",
                                            cor_threshold = -0.5, prev_threshold = prev_thresh,
                                            expected_comp = bicom_compositions, experiment_name = experiment)
  contam_df <- decontam_results$contaminants_removed
  ps_clean  <- decontam_results$clean_phyloseq

  ps_nini     <- prepare_metadata(ps_filtered, group_col = "Replicate")
  cor_results <- plot_bicom_correlation_spectrum(ps_nini, bicom_compositions)
  cor_hst     <- cor_results$plot
} else {
  ps_clean  <- ps_filtered
  contam_df <- data.frame()
  cor_hst   <- NULL
}

# --------------------------
# 6. Format Community Tables
# --------------------------
tax_tab <- as.data.frame(tax_table(ps_clean))
sp_vec  <- setNames(tax_tab$Species, rownames(tax_tab))
gn_vec  <- setNames(tax_tab$Genus,   rownames(tax_tab))

mtdata_final <- data.frame(sample_data(ps_clean)) %>% rownames_to_column("sample")

build_comm_tab <- function(physeq) {
  as.data.frame(otu_table(physeq)) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "ASV", values_to = "count") %>%
    mutate(species = replace_na(sp_vec[ASV], "unassigned"), genus = gn_vec[ASV]) %>%
    left_join(mtdata_final, by = "sample") %>%
    group_by(sample) %>%
    mutate(total_reads = sum(count), rel_ab = count / total_reads) %>%
    ungroup()
}

comm_tab_all_asv_dirty <- build_comm_tab(ps_filtered)
comm_tab_all_asv       <- build_comm_tab(ps_clean)

comm_tab_all_species <- comm_tab_all_asv %>%
  select(-ASV) %>%
  group_by(across(-c(count, rel_ab))) %>%
  summarise(count = sum(count, na.rm = TRUE), rel_ab = sum(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  relocate(sample, SampleName, Passage, Media, Cocktail, Replicate, genus, species, count, rel_ab)

# --------------------------
# 7. Save Outputs
# --------------------------
saveRDS(ps,       file.path(tables_dir, paste0(experiment, "_phyloseq_object.rds")))
saveRDS(ps_clean, file.path(tables_dir, paste0(experiment, "_phyloseq_object_clean.rds")))

write.csv(asv_stats, file.path(tables_dir, paste0(experiment, "_ASV_stats.csv")), row.names = FALSE)
if (nrow(contam_df) > 0) write.csv(contam_df, file.path(tables_dir, paste0(experiment, "_Contaminants.csv")), row.names = FALSE)

write.csv(comm_tab_all_asv,       file.path(tables_dir, paste0(experiment, "_ASV_community.csv")),       row.names = FALSE)
write.csv(comm_tab_all_asv_dirty, file.path(tables_dir, paste0(experiment, "_ASV_community_dirty.csv")), row.names = FALSE)
write.csv(comm_tab_all_species,   file.path(tables_dir, paste0(experiment, "_species_community.csv")),   row.names = FALSE)

ggsave(file.path(figures_dir, paste0(experiment, "_ASV_filtering_plot.pdf")),
       plt, width = FIG_W, height = FIG_H, units = "mm", dpi = 300, device = cairo_pdf)
if (!is.null(cor_hst))
  ggsave(file.path(figures_dir, paste0(experiment, "_ASV_corr_histogram_plot.pdf")),
         cor_hst, width = FIG_W, height = FIG_H, units = "mm", dpi = 300, device = cairo_pdf)