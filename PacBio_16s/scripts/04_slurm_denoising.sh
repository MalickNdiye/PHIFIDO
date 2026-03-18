#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 500000
#SBATCH --partition cpu
#SBATCH --time 5:00:00
#SBATCH --error 2_Bicom8_denoising.err
#SBATCH --output 2_Bicom8_denoising.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

source ~/.bashrc # Load conda
conda activate R # Activate Conda env

# Variables to modify
root=$1
script=$2
outdir=$3
processed_fastq_dir=$4
readcounts=$5
errModel=binnedQualErrfun # use 'binnedQualErrfun' if you have binned quality score, or else 'PacBioErrfun'
db2="$root"/MyBifido_toSpecies.fasta # give dada a set of expected ASVs, or set to ""
pool=T # "T" or "pseudo" or "F", whether to pool samples for ASV inference
maxraref=50000 # use the multiqc output to set this value close to the highest number of reads in your samples

# do not modify below this line, unless you know what you are doing
maxReads=1000000 # reduce if memory issues arise
maxBases=10000000000 # 1E10 strongly recommended
out_denois="$3"/denois_outputs
out_plots="$3"/denois_plots

#rm -rf "$out_denois" # Remove the output directory if it exists
#rm -rf "$out_plots" # Remove the output directory if it exists


# Execute the R script

echo "Parameters:"
echo input.reads: "$processed_fastq_dir"
echo input.readcounts: "$readcounts"
echo maxReads: "$maxReads"
echo errModel: "$errModel"
echo maxBases: "$maxBases"
echo db2: "$db2"
echo pool: "$pool"
echo maxraref: "$maxraref"
echo out.denois: "$out_denois"
echo out.plots: "$out_plots"

Rscript --vanilla "$script" \
    "$processed_fastq_dir" \
    "$readcounts" \
    "$maxReads" \
    "$errModel" \
    "$maxBases" \
    "$db2" \
    "$pool" \
    "$maxraref" \
    "$out_denois" \
    "$out_plots"

echo -e "$(date)"