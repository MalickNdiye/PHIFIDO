# Load required libraries

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(stringr)){
  install.packages(pkgs = 'stringr', repos = 'https://stat.ethz.ch/CRAN/')
  library(stringr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(dada2)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}

if(!require(iNEXT)){
  install.packages(pkgs = 'iNEXT', repos = 'https://stat.ethz.ch/CRAN/')
  library(iNEXT)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10){
  stop(" Usage: 04_denoising.R <processed_reads_dir> <read_count_table> <max_reads_derep> <error_model> <max_bases_errormodel> <priors_db> <pool_T_F> <max_reads_raref> <denoise_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.reads <- args[1] # folder with all pre-processed reads
  input.readcounts <- args[2] # table reporting the number of reads at each step
  maxReads <- args[3] # max number of reads to load at once for dereplication
  errModel <- args[4] # dada2-provided function to estimate the error model
  maxBases <- args[5] # max number of bases to use for error model inference
  db2 <- args[6] # database of ASVs expected in the samples
  pool <- args[7] # "T" or "pseudo" or "F", whether to pool samples for ASV inference
  maxraref <- args[8] # maximum number of reads to extrapolate rarefaction curves
  out.denois <- args[9] # folder to write denoising results
  out.plots <- args[10] # folder to write plots
}

# root <- file.path("/Volumes", "D2c", "mgarcia", "20240708_mgarcia_syncom_assembly", "pacbio_analysis", "run1_bees")
# input.reads <- file.path(root, "results", "preprocessing", "trimmed_filtered_reads")
# input.readcounts <- file.path(root, "results", "preprocessing", "read_count_before_after.tsv")
# maxReads <- 1E6
# errModel <- "binnedQualErrfun"
# maxBases <- 1E10
# db2 <- file.path(root, "data", "databases", "amplicon_based_db/syncom_custom_db_addSpecies.fa")
# pool <- "F"
# maxraref <- 8000
# out.denois <- file.path(root, "results", "denoising")
# out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out.denois, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

maxReads <- as.numeric(maxReads)
maxBases <- as.numeric(maxBases)
maxraref <- as.numeric(maxraref)

filtered_trimmed_reads_paths <- list.files(input.reads, full.names=TRUE)
reads_df <- read.table(input.readcounts, sep = "\t", header = T)

### Dereplicate sequences ###

# if its exists, load pre-existing dereplication object
if (file.exists(file.path(out.denois, "dereps.RDS"))){
  cat("Loading pre-existing dereplication object\n")
  dereps <- readRDS(file.path(out.denois, "dereps.RDS"))
} else {
  cat("Dereplicating sequences\n")
  dereps <- derepFastq(filtered_trimmed_reads_paths, verbose=TRUE, n = maxReads)
  saveRDS(dereps, file.path(out.denois, "dereps.RDS"))
  cat("Dereplication done\n")
}


### Build error model ###

# utility
checkConvergence <- function(dadaO) {
  sapply(dadaO$err_in, function(x) sum(abs(dadaO$err_out-x)))
}

# sample quality scores from 5 random samples to check whether they are binned or not
quals <- c()
for (file in sample(filtered_trimmed_reads_paths, 5, replace = F)){
  quals <- unique(append(quals, unique(plotQualityProfile(file)$data$Score)))
}

cat("Quality scores detected in the trimmed and filtered reads (five randomly chosen samples):\n")
cat(sort(quals))
cat("\n")


# if its exists, load pre-existing error model
if (file.exists(file.path(out.denois, "dada2_error_model.RDS"))){
  cat("Loading pre-existing error model\n")
  error_model <- readRDS(file.path(out.denois, "dada2_error_model.RDS"))
} else {
  cat("Building error model\n")
  
  set.seed(42)
  
  if (errModel == "binnedQualErrfun"){
    cat("Building error model with bins [3, 10, 17, 22, 27, 35, 40]\n")
    binnedQs <- c(3, 10, 17, 22, 27, 35, 40)
    binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)
    error_model <- learnErrors(
      dereps,
      errorEstimationFunction = binnedQualErrfun,
      nbases = maxBases,
      randomize = T,
      BAND_SIZE = 32,
      multithread = T,
      verbose = T
    )
  } else if (errModel == "PacBioErrfun") {
    cat("Building error model with standard PacBio model\n")
    error_model <- learnErrors(
      dereps,
      errorEstimationFunction = PacBioErrfun,
      nbases = maxBases,
      randomize = T,
      BAND_SIZE = 32,
      multithread = T,
      verbose = T
    )
  }
  saveRDS(error_model, file.path(out.denois, "dada2_error_model.RDS"))
  cat("Error model built\n")
}

cat("Error model residuals:\n")
print(checkConvergence(error_model))

err.plot <- plotErrors(error_model, nominalQ=TRUE)
ggsave(file.path(out.plots, paste0("04_error_plot_", errModel, ".pdf")), err.plot, device="pdf", width = 10, height = 8)

### ASV inference ###

cat("Denoising into ASVs\n")

# options:
# pooling T/pseudo/F
# priors yes/no

if (db2 == ""){
  cat("No priors given\n")
  if (pool=="T"){
    dds <- dada(dereps, err=error_model, pool=TRUE, multithread=TRUE, verbose = T)
    } else if (pool=="pseudo") {
    dds <- dada(dereps, err=error_model, pool="pseudo", multithread=TRUE, verbose = T)
    } else if (pool=="F") {
    dds <- dada(dereps, err=error_model, pool=FALSE, multithread=TRUE, verbose = T)
    } else {
      stop("Incorrect pooling option provided\n")
    }
  } else {
    cat(paste0("Using ", db2, " as priors\n"))
    my_priors <- getSequences(db2)
    if (pool=="T"){
      dds <- dada(dereps, err=error_model, pool=TRUE, priors=my_priors, multithread=TRUE, verbose = T)
    } else if (pool=="pseudo") {
      dds <- dada(dereps, err=error_model, pool="pseudo", priors=my_priors, multithread=TRUE, verbose = T)
    } else if (pool=="F") {
      dds <- dada(dereps, err=error_model, pool=FALSE, priors=my_priors, multithread=TRUE, verbose = T)
    } else {
      stop("Incorrect pooling option provided\n")
    }
  }

saveRDS(dds, file.path(out.denois, "denoised_seqs.rds"))

cat("Denoising done\n")


### Generate ASV table ###

cat("Generating ASV table\n")

ASV_samples_table <- makeSequenceTable(dds)

cat(paste0("Found ", ncol(ASV_samples_table), " ASVs across the ", nrow(ASV_samples_table), " samples", "\n"))

### Remove chimera ###

cat("Removing chimera\n")

if (pool %in% c("pseudo", "F")){
  ASV_samples_table_noChim <- removeBimeraDenovo(ASV_samples_table, verbose = T, multithread = T)
} else {
  cat("Removing chimera in pooled mode")
  ASV_samples_table_noChim <- removeBimeraDenovo(ASV_samples_table, method = "pooled", verbose = T, multithread = T)
}

#clean sample names
names <- basename(rownames(ASV_samples_table_noChim))
names <- gsub(".fastq.gz", "", names)
rownames(ASV_samples_table_noChim) <- names

saveRDS(ASV_samples_table_noChim, file.path(out.denois, "ASV_samples_table_noChim.rds"))

cat(paste0("Number of reads removed: ", sum(ASV_samples_table)-sum(ASV_samples_table_noChim), "\n"))
cat(paste0("Proportion of reads removed: ", round(1-sum(ASV_samples_table_noChim)/sum(ASV_samples_table),5), "\n"))

cat("Chimera removal done\n")


# ### Rarefaction curves ###

# cat("Generating rarefaction curves\n")

# # remove zero-abundance samples
# totab <- apply(ASV_samples_table_noChim, 1, sum)
# zeroes <- names(which(totab == 0))

# inext_input <- t(ASV_samples_table_noChim[setdiff(names(totab), zeroes),])

# # using iNEXT so I can also estimate the sampling coverage
# dt <- iNEXT(
#   inext_input, # samples must be as columns
#   q = c(0,1),
#   datatype = "abundance",
#   endpoint = maxraref,
#   knots = 50,
#   se = TRUE,
#   conf = 0.95,
#   nboot = 20
# )

# inextqd <- dt$iNextEst$size_based %>%
#   dplyr::rename(SampleID = Assemblage)

# qd.plot <- ggplot(
#   inextqd[inextqd$Method != "Extrapolation", ],
#   aes(
#     x = m,
#     y = qD,
#     group = SampleID
#   )
# ) +
#   geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
#   geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
#   scale_color_manual(name = "", values = c(low = "#669bbc", high = "#e76f51"), labels = c(low = "Lowest depth", high = "Highest depth")) +
#   geom_line(alpha = 0.4) +
#   geom_ribbon(
#     aes(
#       ymin = qD.LCL,
#       ymax = qD.UCL
#     ), alpha = 0.1
#   ) +
#   geom_label(
#     data = inextqd[inextqd$Method == "Observed", ],
#     aes(
#       x = m,
#       y = qD,
#       label = SampleID
#     ), size = 4, nudge_x = 70
#   ) +
#   theme_bw() +
#   labs(
#     x = "# of reads",
#     y = "# of ASVs"
#   ) +
#   theme(
#     legend.position = "inside",
#     legend.position.inside = c(0.1,0.9),
#     legend.background = element_rect(fill=alpha('white', 0.4))
#   ) +
#   facet_wrap( ~ Order.q, scales = "free_y")

# ggsave(file.path(out.plots, "04_rarefaction_curves_ASVs.pdf"), qd.plot, device="pdf", width = 12, height = 8)

# write.table(
#   inextqd,
#   file = file.path(out.denois, "inext_data.tsv"),
#   sep = "\t",
#   quote = F,
#   row.names = F,
#   col.names = T
# )

### Compute general statistics ###

cat("Plotting global statistics\n")

parts <- str_count(pattern = "/", rownames(ASV_samples_table_noChim)[[1]])

track <- tibble(sample = rownames(ASV_samples_table_noChim),
                Denoising = rowSums(ASV_samples_table), 
                `Chimera removal` = rowSums(ASV_samples_table_noChim)) %>% 
  mutate(basename = paste(sample, ".fastq.gz", sep = "")) %>%
  dplyr::select(-"sample") %>%
  pivot_longer(c("Denoising", "Chimera removal"), names_to = "stage", values_to = "n_reads")

reads_df2 <- rbind(reads_df, track) %>% arrange(basename) %>% 
  separate(basename, into = c("sample", NA), sep = ".fastq.gz", remove = F)

reads_df2$stage <- factor(reads_df2$stage, levels = c("Input", "Primers removed", "Post processing", "Denoising", "Chimera removal"), ordered = T)

track_medians <- reads_df2 %>% 
  group_by(stage) %>% 
  summarize(median = median(n_reads))

# plots
reads.plot1 <- ggplot(
  reads_df2,
  aes(
    y=n_reads,
    x=stage,
    group=basename
  )
) +
  geom_line() +
  geom_label(
    data = reads_df2[reads_df2$stage == "Chimera removal", ],
    aes(
      y = n_reads,
      x = stage,
      label = sample
    ), size = 3, nudge_x = 0.1
  ) +
  geom_point(
    data = track_medians,
    aes(
      y=median,
      x=stage,
      color="median"
    ), inherit.aes = F, size = 3
  ) +
  scale_color_manual(values = c(median = "red"), label = c(median = "Median"), name = "") +
  theme_bw() +
  labs(
    x = "Step",
    y = "# of sequences"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.1,0.1),
    legend.background = element_rect(fill=alpha('white', 0.4))
  )

reads.plot2 <- ggplot(
  reads_df2,
  aes(
    x=n_reads
  )
) +
  geom_histogram(bins = 20, boundary = 0, closed = "left") +
  theme_bw() +
  labs(
    x = "# of sequences",
    y = "# of samples"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = length(unique(reads_df2$stage)), scales = "free_y")

ggsave(file.path(out.plots, "04_read_number_lines.pdf"), reads.plot1, device="pdf", width = 10, height = 8)
ggsave(file.path(out.plots, "04_read_number_hist.pdf"), reads.plot2, device="pdf", width = 4, height = 8)

write.table(
  reads_df2,
  file = file.path(out.denois, "read_counts_steps.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cat("Finished plotting global statistics\n")


cat("Denoising done\n")