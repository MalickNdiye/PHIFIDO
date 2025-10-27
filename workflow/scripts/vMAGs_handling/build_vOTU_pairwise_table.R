#!/usr/bin/env Rscript
library("optparse")
library(data.table)
library(tidyverse)

option_list = list(
  make_option(c("-i", "--drep_data"), type="character", default=NULL,
              help="dRep directory", metavar="character"),
  make_option(c("-s", "--vOTU_summary"), type="character", default=NULL,
              help="Output of parse_drep_viruses_average R script", metavar="character"),
  make_option(c("-o", "--outfile_all"), type="character", default=NULL,
              help="All-vs-all output table", metavar="character"),
  make_option(c("-p", "--outfile_problem"), type="character", default=NULL,
              help="Problematic vOTU members", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("--- Reading vOTU summary\n")
vOTU_data <- fread(opt$vOTU_summary) %>% select(genome, vOTU, checkv_quality)

cat("--- Reading dRep pairwise data\n")

