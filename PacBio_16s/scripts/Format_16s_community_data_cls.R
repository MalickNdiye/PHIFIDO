#!/usr/bin/env Rscript

if(!require(optparse)){
  install.packages(pkgs = 'optparse', repos = 'https://stat.ethz.ch/CRAN/')
  library(optparse)
}

library(vegan)
library(tidyverse)
library(data.table)
library(ggplot2)

if(!require(phyloseq)){
  if(!require(BiocManager)){
    install.packages(pkgs = 'BiocManager', repos = 'https://stat.ethz.ch/CRAN/')
  }
    source("http://bioconductor.org/biocLite.R")
    biocLite("phyloseq")
}

# --------------------------
# 1. Command-Line Arguments
# --------------------------
option_list <- list(
  make_option(c("-a", "--asv_table"), type="character", default="processed_data/ASV_table_merged_new.rds", 
              help="Path to ASV table RDS [default= %default]"),
  make_option(c("-t", "--taxonomy"), type="character", default="processed_data/ASV_Taxonomy_sp_new.RDS", 
              help="Path to Taxonomy RDS [default= %default]"),
  make_option(c("-m", "--metadata"), type="character", default="../data/sample_metadata.csv", 
              help="Path to sample metadata CSV [default= %default]"),
  make_option(c("-q", "--qpcr"), type="character", default="NULL", 
              help="Path to qPCR data CSV. Use 'NULL' if no data [default= %default]"),
  make_option(c("-e", "--experiment"), type="character", default="Bicom6", 
              help="Experiment name for filtering and file prefixing [default= %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="1_Output/", 
              help="Base output directory [default= %default]"),
  make_option(c("-p", "--prev_thresh"), type="numeric", default=0.1, 
              help="Prevalence threshold for filtering [default= %default]"),
  make_option(c("-r", "--relabund_thresh"), type="numeric", default=0.1, 
              help="Relative abundance threshold for filtering [default= %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Map arguments to variables
ASV_samples_table_noChim_p <- opt$asv_table
ASV_taxonomy_sp_p <- opt$taxonomy
metadata_p <- opt$metadata
qpcr_p <- if(opt$qpcr == "NULL" || opt$qpcr == "") NULL else opt$qpcr
experiment <- opt$experiment
output_dir <- opt$output_dir
prev_thresh <- opt$prev_thresh     
relabund_thresh <- opt$relabund_thresh 

# Directory Setup
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Expected Bicom members 
bicom_compositions <- list(
  Bicom1 = c("ESL0820", "ESL0822", "ESL0198", "ESL0199", "ESL0200", "ESL0819"),
  Bicom2 = c("ESL0820", "ESL0822", "ESL0198", "ESL0824", "ESL0827", "ESL1060"),
  Bicom3 = c("ESL0199", "ESL0200", "ESL0819", "ESL0824", "ESL0827", "ESL1069"),
  Bicom4 = c("ESL0820", "ESL0199", "ESL0819", "ESL1060", "ESL1069", "ESL1077"),
  Bicom5 = c("ESL0822", "ESL0200", "ESL0824", "ESL1060", "ESL1069", "ESL1077"),
  Bicom6 = c("ESL0198", "ESL0827", "ESL1060", "ESL1069", "ESL1077", "ESL0819")
)

# --------------------------
# 2. Build Phyloseq
# --------------------------
ASV_samples_table_noChim <- readRDS(ASV_samples_table_noChim_p) 
rownames(ASV_samples_table_noChim) <- gsub(".hifi_reads", "", rownames(ASV_samples_table_noChim))
ASV_taxonomy_sp <- readRDS(ASV_taxonomy_sp_p)

metadata <- read.csv(metadata_p, sep = ",") 
samples_to_remove <- metadata$SampleID[metadata$Exp != experiment]


qpcr <- fread(qpcr_p) %>% select(SampleID, copy_number)
metadata <- left_join(metadata, qpcr, by = "SampleID") %>% column_to_rownames("SampleID")

ps <- phyloseq(otu_table(ASV_samples_table_noChim, taxa_are_rows = FALSE), 
               sample_data(metadata), tax_table(ASV_taxonomy_sp))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps <- prune_samples(!sample_names(ps) %in% samples_to_remove, ps)

tmp <- as.data.frame(sample_data(ps))
tmp$original_lib_size <- sample_sums(ps)
tmp$copy_number <- if(!is.null(qpcr_p)) as.numeric(tmp$copy_number) else NA
tmp$limit_of_detection <- if(!is.null(qpcr_p)) tmp$copy_number / tmp$original_lib_size else NA
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
    genus = unique(Genus),
    species = unique(Species),
    ASV_type = case_when(
      grepl("ESL", species) ~ "BiCom",
      genus == "Bifidobacterium" & is.na(species) ~ "Unassigned Bifidobacterium ASV",
      is.na(genus) ~ "Genus Unassigned",
      TRUE ~ "Other Genus"
    ),
    prevalence = mean(Abundance > 0),
    max_rel = max(rel_abund),
    median_rel = median(rel_abund),
    .groups = "drop"
  ) %>%
  filter(!is.na(median_rel)) %>%
  arrange(median_rel) 

plt <- asv_stats %>%
  mutate(median_rel_bin = cut(median_rel, breaks = c(0, 0.001, 0.01, 0.1, 1), include.lowest = TRUE, dig.lab = 3)) %>%
  ggplot(aes(x = prevalence, y = max_rel, size = median_rel_bin, col = ASV_type, label = ASV)) +
  geom_text(data = subset(asv_stats, prevalence > 0.05), hjust = 1.5, size = 3) +
  geom_point() +
  geom_hline(yintercept = relabund_thresh, linetype = "dashed", col = "darkred") +
  geom_vline(xintercept = prev_thresh, linetype = "dashed", col = "darkred") +
  labs(x = "Prevalence", y = "Maximum relative abundance (%)", title = experiment, size = "Median rel abundance (%)", col = "ASV origin") +
  scale_size_manual(values = c(2,3,4,5,6), drop = FALSE) +
  theme_minimal()

keep_asvs <- asv_stats %>% filter(!(prevalence < prev_thresh & max_rel < relabund_thresh)) %>% pull(ASV)
ps_filtered <- prune_taxa(keep_asvs, ps)

# --------------------------
# 4. Decontamination Core Functions
# --------------------------
prepare_metadata <- function(ps, group_col, conc_col) {
  sd <- sample_data(ps)
  if (conc_col %in% sample_variables(ps) && conc_col != "original_lib_size") {
    sd[[conc_col]] <- as.numeric(sd[[conc_col]])
    sd[[conc_col]][is.na(sd[[conc_col]]) | sd[[conc_col]] <= 0] <- 1e-6
  }
  sd_df <- as.data.frame(sd)
  sd$Group <- apply(sd_df[, group_col, drop = FALSE], 1, paste, collapse = "_")
  sample_data(ps) <- sd
  return(ps)
}

find_via_correlation <- function(seqtab_sub, conc_var_vals, cor_threshold = -0.5) {
  seqtab_rel <- sweep(seqtab_sub, 1, pmax(rowSums(seqtab_sub), 1), "/")
  sp_cors <- apply(seqtab_rel, 2, function(rel_abund) {
    valid_idx <- which(rel_abund > 0 & !is.na(conc_var_vals) & !is.na(rel_abund))
    if (length(valid_idx) < 3 || sd(rel_abund[valid_idx], na.rm = TRUE) == 0) return(NA)
    conc_log <- log10(pmax(conc_var_vals[valid_idx], 1e-6))
    rel_log <- log10(rel_abund[valid_idx])
    if (sd(conc_log, na.rm = TRUE) == 0) return(NA)
    cor(rel_log, conc_log, method = "spearman", use = "complete.obs")
  })
  
  contam_sp <- names(sp_cors)[which(sp_cors < cor_threshold)]
  if (length(contam_sp) == 0) return(data.frame(Species_Group = character(), cor = numeric()))
  return(data.frame(Species_Group = contam_sp, cor = sp_cors[contam_sp], row.names = NULL))
}

find_via_prevalence <- function(seqtab_sub, prev_threshold = 0.1) {
  prev_calc <- colSums(seqtab_sub > 0) / nrow(seqtab_sub)
  contam_asvs <- names(prev_calc)[prev_calc < prev_threshold]
  if (length(contam_asvs) == 0) return(data.frame(ASV = character(), prevalence = numeric()))
  return(data.frame(ASV = contam_asvs, prevalence = prev_calc[contam_asvs], row.names = NULL))
}

identify_group_contaminants <- function(ps, conc_var_name, cor_threshold, prev_threshold = 0.1) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("ASV") %>%
    mutate(Species_Group = case_when(
      !is.na(Species) & Species != "" & Species != "NA" & Species != "unassigned" ~ paste(Genus, Species, sep = "_"),
      TRUE ~ paste("Unclassified_ASV", ASV, sep = "_") 
    ))
  
  agg_otu <- t(rowsum(t(otu_mat), group = tax_df$Species_Group[match(colnames(otu_mat), tax_df$ASV)]))
  
  list_of_contaminants <- lapply(unique(sample_data(ps)$Group), function(grp) {
    if (is.na(grp)) return(NULL)
    samples_to_keep <- sample_names(ps)[sample_data(ps)$Group == grp]
    if (length(samples_to_keep) < 3) return(NULL)
    
    seqtab_sub_asv <- otu_mat[samples_to_keep, colSums(otu_mat[samples_to_keep, , drop = FALSE]) > 0, drop = FALSE]
    seqtab_sub_sp <- agg_otu[samples_to_keep, colSums(agg_otu[samples_to_keep, , drop = FALSE]) > 0, drop = FALSE]
    if (ncol(seqtab_sub_asv) == 0 || ncol(seqtab_sub_sp) == 0) return(NULL)
    
    conc_vals <- sample_data(ps)[[conc_var_name]][match(rownames(seqtab_sub_asv), sample_names(ps))]
    
    res_cor <- find_via_correlation(seqtab_sub_sp, conc_vals, cor_threshold)
    res_prev <- find_via_prevalence(seqtab_sub_asv, prev_threshold)
    if(nrow(res_cor) == 0) res_cor <- data.frame(Species_Group = character(), cor = numeric())
    if(nrow(res_prev) == 0) res_prev <- data.frame(ASV = character(), prevalence = numeric())
    
    flagged_sp <- unique(c(res_cor$Species_Group, tax_df$Species_Group[tax_df$ASV %in% res_prev$ASV]))
    if (length(flagged_sp) == 0) return(NULL)
    
    tax_df %>% filter(Species_Group %in% flagged_sp) %>% select(ASV, Species_Group, Genus, Species) %>%
      mutate(Condition = grp) %>%
      left_join(res_cor, by = "Species_Group") %>% 
      left_join(res_prev, by = "ASV") %>%
      mutate(Method = sub(";+$", "", paste0(
        ifelse(Species_Group %in% res_cor$Species_Group, "correlation;", ""),
        ifelse(ASV %in% res_prev$ASV, "prevalence", "")
      )))
  })
  bind_rows(list_of_contaminants) %>% select(Condition, ASV, Species_Group, Genus, Species, cor, prevalence, Method)
}

remove_targeted_contaminants <- function(ps, contam_df) {
  otu_mat <- as(otu_table(ps), "matrix")
  is_taxa_rows <- taxa_are_rows(ps)
  if (is_taxa_rows) otu_mat <- t(otu_mat) 
  
  for (i in seq_len(nrow(contam_df))) {
    target_samples <- intersect(sample_names(ps)[sample_data(ps)$Group == contam_df$Condition[i]], rownames(otu_mat))
    if (contam_df$ASV[i] %in% colnames(otu_mat) && length(target_samples) > 0) {
      otu_mat[target_samples, contam_df$ASV[i]] <- 0
    }
  }
  otu_table(ps) <- if(is_taxa_rows) otu_table(t(otu_mat), taxa_are_rows = TRUE) else otu_table(otu_mat, taxa_are_rows = FALSE)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

clean_microbiome_data <- function(ps, group_col = "Replicate", conc_var = "original_lib_size", 
                                  cor_threshold = -0.5, prev_threshold = 0.1, 
                                  expected_comp = NULL, experiment_name = NULL) { 
  
  cat("1. Formatting metadata...\n")
  ps <- prepare_metadata(ps, group_col, conc_var)
  
  cat(sprintf("2. Identifying contaminants using %s...\n", conc_var))
  contam_df <- identify_group_contaminants(ps, conc_var, cor_threshold, prev_threshold)
  
  if (!is.null(contam_df) && nrow(contam_df) > 0) {
    
    # --- NEW: Calculate Mean and Median Relative Abundance ---
    # 1. Extract OTU matrix and convert to relative abundance
    otu_mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    rel_abund_mat <- sweep(otu_mat, 1, pmax(rowSums(otu_mat), 1), "/")
    
    # 2. Get sample mapping
    samp_condition <- as.character(sample_data(ps)$Group)
    samp_names <- rownames(otu_mat)
    
    # 3. Append the stats to the dataframe
    contam_df <- contam_df %>% 
      rowwise() %>%
      mutate(
        mean_rel_abund = {
          valid_samps <- samp_names[samp_condition == Condition]
          if (length(valid_samps) > 0 && ASV %in% colnames(rel_abund_mat)) {
            mean(rel_abund_mat[valid_samps, ASV], na.rm = TRUE)
          } else { NA_real_ }
        },
        median_rel_abund = {
          valid_samps <- samp_names[samp_condition == Condition]
          if (length(valid_samps) > 0 && ASV %in% colnames(rel_abund_mat)) {
            median(rel_abund_mat[valid_samps, ASV], na.rm = TRUE)
          } else { NA_real_ }
        }
      ) %>%
      ungroup()
    # ---------------------------------------------------------
    
    # 3. Apply False Positive (FP) logic ONLY if experiment is "Bicom6" AND a list is provided
    if (!is.null(expected_comp) && identical(experiment_name, "Bicom6")) {
      contam_df <- contam_df %>% rowwise() %>%
        mutate(FP = {
          pieces <- strsplit(as.character(Condition), "_")[[1]]
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

plot_bicom_correlation_spectrum <- function(ps, expected_comp, conc_var = "original_lib_size") {
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_rel <- sweep(otu_mat, 1, pmax(rowSums(otu_mat), 1), "/")
  
  sd_df <- as.data.frame(sample_data(ps))
  tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("ASV")
  all_groups <- unique(as.character(sd_df$Group))
  
  all_cor_data <- lapply(all_groups, function(grp) {
    samps <- rownames(sd_df)[sd_df$Group == grp]
    if(length(samps) < 3) return(NULL)
    
    pieces <- strsplit(grp, "_")[[1]]
    bicom_key <- intersect(pieces, names(expected_comp))
    if(length(bicom_key) == 0) return(NULL) 
    
    allowed_species <- expected_comp[[bicom_key[1]]]
    conc_vals <- as.numeric(sd_df[samps, ][[conc_var]])
    conc_log <- log10(pmax(conc_vals, 1e-6))
    sub_rel <- otu_rel[samps, , drop = FALSE]
    present_asvs <- colnames(sub_rel)[colSums(sub_rel) > 0]
    
    lapply(present_asvs, function(asv_id) {
      x <- sub_rel[, asv_id]
      v_idx <- which(x > 0)
      sp_name <- tax_df$Species[tax_df$ASV == asv_id]
      
      rho <- NA
      if(length(v_idx) >= 3 && sd(x[v_idx]) > 0 && sd(conc_log[v_idx]) > 0) {
        rho <- cor(log10(x[v_idx]), conc_log[v_idx], method = "spearman")
      }
      data.frame(ASV = asv_id, Species = sp_name, Condition = grp, Correlation = rho, 
                 Type = ifelse(sp_name %in% allowed_species, "Expected Member", "Unexpected (Contaminant)"))
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  plot_df <- all_cor_data %>% filter(!is.na(Correlation))
  
  plt <- ggplot(plot_df, aes(x = Correlation, fill = Type)) +
    geom_histogram(binwidth = 0.05, alpha = 0.7, color = "white", position = "identity") +
    scale_fill_manual(values = c("Expected Member" = "#4DAF4A", "Unexpected (Contaminant)" = "#E41A1C")) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    labs(x = "Spearman Rho", y = "ASV Count", fill = "Biological Status") +
    theme_minimal() + theme(legend.position = "bottom")
  
  return(list(plot = plt, data = plot_df))
}

# --------------------------
# 5. Execute Pipeline
# --------------------------
if (experiment == "Bicom6") {
  decontam_results <- clean_microbiome_data(ps_filtered, group_col = "Replicate", conc_var = "original_lib_size", 
                                            expected_comp = bicom_compositions, experiment_name = experiment, 
                                            cor_threshold = -0.5, prev_threshold = prev_thresh) 
  contam_df <- decontam_results$contaminants_removed 
  ps_clean <- decontam_results$clean_phyloseq
  
  ps_nini <- prepare_metadata(ps_filtered, group_col = "Replicate", conc_col = "original_lib_size")
  cor_results <- plot_bicom_correlation_spectrum(ps_nini, bicom_compositions, conc_var = "original_lib_size")
  cor_hst <- cor_results$plot
} else {
  ps_clean <- ps_filtered
  contam_df <- data.frame() 
  cor_hst<- NULL
}

# --------------------------
# 6. Format Community Tables
# --------------------------
tax_tab <- as.data.frame(tax_table(ps_clean))
sp_vec <- setNames(tax_tab$Species, rownames(tax_tab))
gn_vec <- setNames(tax_tab$Genus, rownames(tax_tab))

mtdata_final <- as.data.frame(sample_data(ps_clean)) %>% rownames_to_column("sample")

build_comm_tab <- function(physeq) {
  as.data.frame(otu_table(physeq)) %>% 
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "ASV", values_to = "count") %>%
    mutate(species = replace_na(sp_vec[ASV], "unassigned"), genus = gn_vec[ASV]) %>%
    left_join(mtdata_final, by = "sample") %>%
    group_by(sample) %>%
    mutate(total_reads = sum(count), limit_of_detection = copy_number / total_reads, 
           rel_ab = count / total_reads, abs_ab = rel_ab * copy_number) %>% 
    ungroup() 
}

comm_tab_all_asv_dirty <- build_comm_tab(ps_filtered)
comm_tab_all_asv <- build_comm_tab(ps_clean)

comm_tab_all_species <- comm_tab_all_asv %>%
  select(-ASV) %>%
  group_by(across(-c(count, rel_ab, abs_ab))) %>%
  summarise(count = sum(count, na.rm = TRUE), rel_ab = sum(rel_ab, na.rm = TRUE), 
            abs_ab = sum(abs_ab, na.rm = TRUE), .groups = "drop") %>%
  relocate(sample, SampleName, Passage, Media, Cocktail, Replicate, genus, species, count, rel_ab, abs_ab)

# --------------------------
# 7. Save Outputs
# --------------------------
saveRDS(ps, file.path(tables_dir, paste0(experiment, "_phyloseq_object.rds")))
saveRDS(ps_clean, file.path(tables_dir, paste0(experiment, "_phyloseq_object_clean.rds")))

write.csv(asv_stats, file.path(tables_dir, paste0(experiment, "_ASV_stats.csv")), row.names = FALSE)
if(nrow(contam_df) > 0) write.csv(contam_df, file.path(tables_dir, paste0(experiment, "_Contaminants.csv")), row.names = FALSE)

write.csv(comm_tab_all_asv, file.path(tables_dir, paste0(experiment, "_ASV_community.csv")), row.names = FALSE)
write.csv(comm_tab_all_asv_dirty, file.path(tables_dir, paste0(experiment, "_ASV_community_dirty.csv")), row.names = FALSE)
write.csv(comm_tab_all_species, file.path(tables_dir, paste0(experiment, "_species_community.csv")), row.names = FALSE)

ggsave(file.path(figures_dir, paste0(experiment, "_ASV_filtering_plot.pdf")), plt, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(file.path(figures_dir, paste0(experiment, "_ASV_corr_histogram_plot.pdf")), cor_hst, width = 297, height = 210, units = "mm", dpi = 300)