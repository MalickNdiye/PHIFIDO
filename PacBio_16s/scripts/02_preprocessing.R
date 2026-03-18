# Load required libraries

if(!require(foreach)){
  install.packages(pkgs = 'foreach', repos = 'https://stat.ethz.ch/CRAN/')
  library(foreach)
}

if(!require(doParallel)){
  install.packages(pkgs = 'doParallel', repos = 'https://stat.ethz.ch/CRAN/')
  library(doParallel)
}

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(gridExtra)){
  install.packages(pkgs = 'gridExtra', repos = 'https://stat.ethz.ch/CRAN/')
  library(gridExtra)
}

# if(!require(data.table)){
#   install.packages(pkgs = 'data.table', repos = 'https://stat.ethz.ch/CRAN/')
#   library(data.table)
# }

if(!require(dada2)){
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8){
  stop(" Usage: 02_preprocessing.R <raw_reads_dir> <fwd_primer_sequence> <rev_primer_sequence> <min_read_length> <max_read_length> <maxEE> <preproc_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.raw <- args[1] # folder with all raw read files (if you have pre-rarefied the reads, give the folder with the pre-rarefied fastq files)
  fwd.primer <- args[2] # forward primer sequence
  rev.primer <- args[3] # reverse primer sequence
  minLen <- args[4] # lower bound for read length, for filterAndTrim
  maxLen <- args[5] # upper bound for read length, for filterAndTrim
  maxEE <- args[6] # max number of errors per read
  out.preproc <- args[7] # folder to write pre-processing results
  out.plots <- args[8] # folder to write plots (can be the same folder throughout the pipeline)
}

### Create outdirs ###

print("Creating directories")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out.preproc, "primerfree_reads"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out.preproc, "trimmed_filtered_reads"), recursive = TRUE, showWarnings = FALSE)

### Inputs ###

minLen <- as.numeric(minLen)
maxLen <- as.numeric(maxLen)

raw_reads_paths <- list.files(input.raw, full.names=TRUE)



### Set up parallel backend ###

numCores <- detectCores()
print(paste0("Using ", numCores, " cores for parallel processing"))
cl <- makeCluster(numCores)
registerDoParallel(cl)



### Orient sequences and trim primers ###

print("Orienting sequences and trimming primers")

# Define paths for primers-free reads
trimmed_reads_paths <- file.path(out.preproc, "primerfree_reads", basename(raw_reads_paths))

print(paste0("Primer-free reads written to: ", file.path(out.preproc, "primerfree_reads")))


primer_removal_summary <- foreach(i = seq_along(raw_reads_paths), .packages = c("dada2"), .combine = 'rbind') %dopar% {
  print(paste0("Processing ", basename(raw_reads_paths[i])))
  res <- removePrimers(fn = raw_reads_paths[i],
                       fout = trimmed_reads_paths[i], 
                       primer.fwd = fwd.primer,
                       primer.rev = rc(rev.primer),
                       orient=TRUE,
                       verbose=TRUE) # cannot multithread with dada2
}

print(paste0("Mean proportion of reads removed: ", round(mean(primer_removal_summary[,"reads.out"]/primer_removal_summary[,"reads.in"]),2)))

saveRDS(primer_removal_summary, file.path(out.preproc, "primer_removal_summary.rds"))
# primer_removal_summary <- readRDS(file.path(out.preproc, "primer_removal_summary.rds"))

primers_summary <- as.data.frame(primer_removal_summary)
primers_summary <- primers_summary %>%
  mutate(basename = rownames(primers_summary)) %>%
  pivot_longer(c("reads.in", "reads.out"), values_to = "n_reads", names_to = "stage") %>%
  mutate(stage = ifelse(stage == "reads.in", "Input", "Primers removed"))

print("Orientation and primer removal done")



### Trim and filter reads based on quality and length ###

print("Trimming and filtering reads")

# Define paths for trimmed and filtered reads
filtered_trimmed_reads_paths <- file.path(out.preproc, "trimmed_filtered_reads", basename(trimmed_reads_paths))

track_filtering <- foreach(i = seq_along(trimmed_reads_paths), .packages = c("dada2"), .combine = 'rbind') %dopar% {
  filterAndTrim(fwd = trimmed_reads_paths[i],
                filt = filtered_trimmed_reads_paths[i],
                minLen=minLen,
                maxLen=maxLen,
                maxN=0, # no ambiguous bases
                rm.phix=FALSE,
                maxEE=3, # max two expected errors per sequence, calculated from the quality score
                multithread = TRUE,
                verbose = T) 
}

print(paste0("Mean proportion of reads removed: ", round(mean(track_filtering[,"reads.out"]/track_filtering[,"reads.in"]),2)))

filtering_summary <- as.data.frame(track_filtering)
filtering_summary <- filtering_summary %>%
  mutate(basename = rownames(filtering_summary)) %>%
  pivot_longer(c("reads.in", "reads.out"), values_to = "n_reads", names_to = "stage") %>%
  mutate(stage = ifelse(stage == "reads.in", "Primers removed", "Post processing"))

print("Trimming and filtering reads done")

### Summary stats on number of reads and read length ###

print("Computing read statistics")

# Read length

compute_all_readL <- function(paths, stage){
  list <- foreach(i = seq_along(paths), .packages = c("dada2"), .inorder = FALSE) %dopar% {
    my_fastq_seq <- getSequences(paths[i])
    lens <- unname(sapply(my_fastq_seq, nchar))
    
    data.frame(
      basename = rep(basename(paths[i]), length(lens)),
      read = seq_along(lens),
      length = lens,
      stage = stage,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, list)
  return(df)
}

lens_raw_df <- compute_all_readL(raw_reads_paths, "Input")
lens_bf_df <- compute_all_readL(trimmed_reads_paths, "Primers removed")
lens_af_df <- compute_all_readL(filtered_trimmed_reads_paths, "Post processing")

lens_df <- rbind(lens_raw_df, lens_bf_df, lens_af_df) %>% arrange(basename) # combine before and after data
lens_df$stage <- factor(lens_df$stage, levels = c("Input", "Primers removed", "Post processing"), ordered = T)

# plot read length before and after pre-processing
length.plot <- ggplot(
  lens_df,
  aes(
    x = length
  )
) +
  geom_histogram(bins = 50, boundary = 0, closed = "left") +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Length (bp)",
    y = "# of sequences"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = 3)


# Read count

reads_df <- rbind(primers_summary, filtering_summary[filtering_summary$stage == "Post processing", ]) %>% arrange(basename)
reads_df$stage <- factor(reads_df$stage, levels = c("Input", "Primers removed", "Post processing"), ordered = T)

ns.plot <- ggplot(
  reads_df,
  aes(
    x = n_reads
  )
) +
  geom_histogram(bins = 20, boundary = 0, closed = "left") +
  theme_bw()+
  labs(
    x = "# of sequences",
    y = "# of samples"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = 3)

qc_afbf <- grid.arrange(length.plot, ns.plot, ncol=2)

print("Saving plots and tables")

# save plot
ggsave(file.path(out.plots, "02_read_stats_before_after.pdf"), qc_afbf, device="pdf", width = 8, height = 8)

# save stats
write.table(lens_df, file.path(out.preproc, "read_length_before_after.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(reads_df, file.path(out.preproc, "read_count_before_after.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


print("All pre-processing done")

### Stop cluster ###
stopCluster(cl)


