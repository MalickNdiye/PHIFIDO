library(tidyverse)
library(data.table)

if(!require(gtools)){
  install.packages(pkgs = 'gtools', repos = 'https://stat.ethz.ch/CRAN/')
  library(gtools)
}


################################################################################
# Set Variables
################################################################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5){
  stop(" Usage: format_qPCR_data.R <metadata> <std_curve> <pcr_data> <Experiment> <output_dir>", call.=FALSE)
} else {
  metadata_p <- args[1] # path to metadata
  std_curve_p <- args[2] # path to qPCR standard curve data
  pcr_data <- args[3] # path to raw qPCR data
  Experiment <- args[4] # name of experiment (e.g. "Bicom6")
  output_dir <- args[5] # output directory
}


if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

print("input parameters:")
print(paste("Metadata path:", metadata_p))
print(paste("Standard curve path:", std_curve_p))
print(paste("qPCR data path:", pcr_data))
print(paste("Experiment name:", Experiment))
print(paste("Output directory:", output_dir))

################################################################################
# Open Metadata
################################################################################
print(paste("=== Loading metadata from", metadata_p, "==="))
metadata<- fread(metadata_p)


################################################################################
# qPCR standard Curve
################################################################################
print(paste("=== Loading qPCR standard curve data from", std_curve_p, "==="))
std_curve<- fread(std_curve_p) %>%
  mutate(log10_copy_number=log10(copy_number))

lm_std <- lm(Cq~log10_copy_number, data=std_curve)
summary(lm_std)

int<- coef(lm_std)[1]
slope <- coef(lm_std)[2]

# Calculate efficiency
efficiency <- 10^(-1/slope) - 1
efficiency_percent <- efficiency * 100

print(paste("Standard Curve Equation: Cq =", round(slope, 4), "* log10(copy_number) +", round(int, 4)))
print(paste("qPCR Efficiency:", round(efficiency_percent, 2), "%"))


################################################################################
# Clean & Format qPCR data
print(paste("=== Loading and formatting qPCR data from", pcr_data, "for experiment", Experiment, "==="))
Absolute_quant_raw<- fread(pcr_data) %>%
  filter(Exp==Experiment)

if(nrow(Absolute_quant_raw) == 0){
  print(paste("No qPCR data found in", pcr_data, "for experiment", Experiment, ". Saving empty data frame."))
  empty_df <- data.frame(SampleID=character(), Passage=character(), Cocktail=character(), Media=character(), Replicate=integer(), copy_number=numeric(), suspicious=logical())
  write.csv(empty_df, file.path(output_dir, paste0(Experiment, "_Absolute_quant_final.csv")), row.names = F, quote = F)
  quit(save = "no")
}


print(paste("Total qPCR measurements for", Experiment, "before cleaning:", nrow(Absolute_quant_raw)))
Absolute_quant_raw2<- Absolute_quant_raw %>%
  group_by(SampleID) %>%
  #filter(!is.na(Cq), SampleID!="H2O", SampleID!="Blank", !(SampleID=="B4P4" & Plate==1)) %>%
  #slice_max(Plate, n = 1) %>%
  mutate(replicate= row_number(),
         Cq = as.numeric(Cq),
         SampleID = ifelse(
           grepl("^[A-Z][0-9][A-Z][0-9]+$", SampleID),  # check the pattern
           sub("^(.{2})(.*)$", "\\1-\\2", SampleID),    # insert dash after 2 chars
           SampleID
         ),
         n=n()) %>%
  arrange(SampleID) %>%
  ungroup()

if (Experiment == "Bicom6"){
  Absolute_quant_raw2$SampleID<- paste(Absolute_quant_raw2$SampleID, "Bicom6", sep="-")
}

print(paste("Total qPCR measurements for", Experiment, "after initial formatting:", nrow(Absolute_quant_raw2)))
Absolute_quant_clean<- Absolute_quant_raw2 %>%
  group_by(SampleID) %>%
  mutate(.orig_row = row_number()) %>%
  group_modify(~{
    x <- .x
    n <- nrow(x)
    cq <- x$Cq
    
    # handle 1-replicate case: keep it but mark suspicious (no replicate to compare)
    if (n == 1) {
      x$suspicious <- TRUE
      return(x)
    }
    
    # pairwise absolute differences matrix
    D <- abs(outer(cq, cq, "-"))
    diag(D) <- Inf   # ignore diagonal when searching for min / close pairs
    
    # find indices that have at least one neighbor within 0.5 Cq
    has_close <- apply(D, 1, function(r) any(r <= 0.5, na.rm = TRUE))
    
    if (any(has_close)) {
      # keep all replicates that have >=1 neighbor within 0.5
      keep_idx <- which(has_close)
      out <- x[keep_idx, , drop = FALSE]
      out$suspicious <- FALSE
      return(out)
    } else {
      # no pair within 0.5 -> pick the closest pair (keep exactly two) and mark suspicious TRUE
      # search only upper triangle for the min to avoid duplicate symmetric hits
      D_ut <- D
      D_ut[lower.tri(D_ut, diag = TRUE)] <- Inf
      min_val <- min(D_ut, na.rm = TRUE)
      pos <- which(D_ut == min_val, arr.ind = TRUE)[1, ]   # first occurrence
      keep_idx <- sort(unique(c(pos[1], pos[2])))
      out <- x[keep_idx, , drop = FALSE]
      out$suspicious <- TRUE
      return(out)
    }
  }) %>%
  ungroup() %>%
  # restore original within-sample ordering
  arrange(SampleID, .orig_row) %>%
  select(-.orig_row)

print(paste("Total qPCR measurements for", Experiment, "after cleaning:", nrow(Absolute_quant_clean)))
Absolute_quant_final <- Absolute_quant_clean %>%
  mutate(log10_copies=(Cq - int) / slope,
         copy_number_2=(10^log10_copies)*10,
         copy_number=copy_number_2/2
         ) %>%
  left_join(., metadata, by=c("SampleID", "Exp")) %>%
  group_by(SampleID, Passage, Cocktail, Media, Replicate)%>%
  reframe(copy_number=mean(copy_number),
          suspicious=any(suspicious)) %>%
  mutate(Passage=factor(Passage, levels=mixedsort(unique(Passage)))) %>%
  arrange(SampleID, Replicate)


################################################################################
# Save PCR Data
################################################################################
print(paste("=== Saving formatted qPCR data for", Experiment, "to", output_dir, "==="))
write.csv(Absolute_quant_raw,   file.path(output_dir, paste0(Experiment, "_Absolute_quant_raw.csv")),
          row.names = F, quote = F)
write.csv(Absolute_quant_clean, file.path(output_dir, paste0(Experiment, "_Absolute_quant_clean.csv")),
          row.names = F, quote = F)
write.csv(Absolute_quant_final, file.path(output_dir, paste0(Experiment, "_Absolute_quant_final.csv")),
          row.names = F, quote = F)

print(paste("=== Finished formatting qPCR data for", Experiment, "==="))