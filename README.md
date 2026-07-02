# PHIFIDO
*Author: Malick Ndiaye*

# Introduction
As top predators of the microbial world, phages are supposed to impact microbial diversity. The “kill the winner” hypothesis suggests that phage predation modulates bacterial competition through competitive release, promoting coexistence. In this scenario, higher phage diversity can maintain higher bacterial diversity. However, most evidence supporting this idea comes from ecological surveys documenting oscillations in phage and bacterial densities, without establishing causality. While such patterns could arise from top-down effects of phage predation, they might also result from bottom-up dynamics, where favorable bacterial growth conditions indirectly boost phage diversity. To establish the direction of causality between viral and bacterial diversity, we focused on the simple and tractable gut microbiota of honey bees, combining ecological surveys with experiments. First, paired viral and bacterial metagenomes from 49 individual bees revealed a microbiome-wide interaction network formed of distinct phage-bacteria interaction modules. Within these modules, we found strong correlations between phage and bacterial composition and diversity. Strikingly, these correlations were stronger when considering the bacterial strain-level composition, underscoring the importance of phage strain-specificity in natural systems. Building on these findings, we aim to experimentally test whether diversity correlations arise from top-down phage predation or bottom-up resource effects. We assembled a synthetic community of Bifidobacteria strains isolated from bee guts and isolated phages targeting all strains. Consistent with our metagenomic observations, these phages displayed broad, yet strain-specific host ranges, infecting multiple species but only particular strains within them. Next, we aim to cultivate the synthetic bacterial community with and without phages, in environments varying in nutrient diversity. This full factorial design will enable us to determine whether phages promote bacterial diversity and under which environmental context. Our preliminary results suggest that top-down competitive release is most pronounced in nutrient-poor conditions, where resource competition is intense. Conversely, in nutrient-diverse environments, phages exert a relatively limited impact on bacterial hosts. Together, our findings suggest that complex trophic interactions shape biodiversity.

# Pipeline Overview
Snakemake pipeline supporting the manuscript **"Bacterial Niche Overlap Determines the Ecological Impact of Phage Predation"** (Ndiaye, Sbaghdi, et al., Engel lab, University of Lausanne).

The pipeline processes shotgun metagenomic and PacBio full-length 16S amplicon sequencing data from a passaging experiment in which synthetic *Bifidobacterium* communities (isolated from the honeybee gut) were grown with and without phages, under two resource conditions (sucrose water, low complexity; pollen extract, high complexity). It produces the bacterial community composition, viral (vMAG) catalog, phage-host links, and strain/vOTU-level population genomics data underlying the paper's figures.

## Pipeline structure

The workflow is organized as a `Snakefile` that includes several rule modules, each covering one analysis stage:

| File | Stage |
|---|---|
| `Snakefile` | Entry point: config loading, shared helper functions, read QC/trimming/host-filtering, viral assembly (MEGAHIT) and identification (VIBRANT, GeNomad), contig aggregation/splitting, dereplication (dRep) |
| `Host_analysis.smk` | Bacterial (*Bifidobacterium*) genome QC (CheckM), taxonomy (GTDB-Tk), functional annotation (DRAM), pangenome/phylogeny (OrthoFinder, IQ-TREE, phylogenetic diversity) |
| `vMAGs_info.smk` | Viral genome QC (CheckV), lifestyle prediction (Bacphlip), taxonomy (PhaBOX), phage-host prediction via CRISPR spacer matching, and aggregation/quality filtering of the viral catalog |
| `pangenome_phages.smk` | Phage gene annotation (Pharokka) and genome-genome similarity clustering (VIRIDIC) used to define viral OTUs (vOTUs) |
| `instrain.smk` | Read mapping (Bowtie2) and strain/population-level profiling (inStrain) of representative vMAGs across all metagenomic samples |
| `16s_workflow.smk` | PacBio full-length 16S rRNA amplicon processing (DADA2 denoising, taxonomy, qPCR-based absolute abundance) used to track total community composition over the 10 passages |
|||

Rules are annotated in place with short docstrings describing what each step does and, where relevant, how it connects to neighboring steps.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) with conda/mamba integration (each rule specifies its own `conda:` environment file under `envs/`)
- Access to a SLURM cluster with the `pengel_beemicrophage` account, or equivalent — resource directives (`account`, `mem_mb`, `runtime`) are set for this cluster and will need adjusting for other systems
- Reference databases configured in `../config/config.yaml` (Kraken2, CheckM, GTDB-Tk, CheckV, PhaBOX, GeNomad, DRAM, CRISPR spacers database) and a few tool databases hardcoded as absolute paths in individual rules (e.g. Pharokka DB in `pangenome_phages.smk`) .
- A `config.yaml` defining, at minimum: `samples` (per-experiment R1/R2 paths), `BiCom` (bacterial reference genomes), `sample_metadata`, and the various database paths referenced above

## Running the pipeline

From the `workflow/` directory (where the `Snakefile` lives):

```bash
snakemake --use-conda --profile <cluster_profile> -j <n_jobs>
```

`rule all` defines the default targets: the viral lifestyle/taxonomy tables, the aggregated inStrain profiles, the Kraken2 QC summary, and the mapping statistics. The 16S workflow is not currently listed in `rule all` and must be requested explicitly by target path if needed.

## Key outputs

- `../results/data_validation/` — read QC, trimming, and host-filtering statistics
- `../results/assembly/viral/all_viral_contigs.fasta` / `all_HQ_viral_contigs.fasta` — pooled and quality-filtered viral contig catalog
- `../results/vMAGs/vMAGs_summary.tsv` — master table of viral genome quality, lifestyle, taxonomy, and predicted hosts
- `../results/vMAGs/dereplication/` — dRep clustering results and representative genome lists
- `../results/pangenomics/viruses/viridic_parsed/` — vOTU clustering combined with community composition
- `../results/inStrain/aggregated_data_sing/` — per-sample, per-genome population genomics profiles
- `../PacBio_16s/results/community_data_formatted/` — per-experiment 16S community composition tables (ASV abundance + qPCR-scaled absolute abundance)
