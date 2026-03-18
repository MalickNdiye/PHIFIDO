#!/bin/bash

#SBATCH --account pengel_beemicrophage
#SBATCH --job-name PacBio_16s_preprocessing
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 500000
#SBATCH --partition cpu
#SBATCH --time 02:00:00
#SBATCH --error %x.%j_Bicom6.log
#SBATCH --output %x.%j_Bicom6.out

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

source ~/.bashrc # Load conda
conda activate dada2_env # Activate Conda env

# Variables: modify these paths to your own
root=$1
script=$2

#raw_fastq_dir="$root"/data/raw_reads # if you did not pre-rarefy
raw_fastq_dir=$3 

#output directory
outdir=$4

rm -rf "$outdir" # Remove the output directory if it exists
mkdir -p "$outdir"



fwd_primer=AGRGTTYGATYMTGGCTCAG
rev_primer=RGYTACCTTGTTACGACTT
minLen=1400
maxLen=1600
maxEE=3 # use 2 for normal PacBio, 3 for Kinnex
out_preproc="$outdir"/preprocessing
out_plots="$outdir"/plots

# Execute the R script

echo "Parameters:"
echo input.raw: "$raw_fastq_dir"
echo fwd.primer: "$fwd_primer"
echo rev.primer: "$rev_primer"
echo minLen: "$minLen"
echo maxLen: "$maxLen"
echo maxEE: "$maxEE"
echo out.preproc: "$out_preproc"
echo out.plots: "$out_plots"

echo "Refer to the Rscript for information on the parameters"

Rscript --vanilla "$script" \
    "$raw_fastq_dir" \
    "$fwd_primer" \
    "$rev_primer" \
    "$minLen" \
    "$maxLen" \
    "$maxEE" \
    "$out_preproc" \
    "$out_plots"

echo -e "$(date)"