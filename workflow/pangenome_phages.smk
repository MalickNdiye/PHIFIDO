"""
Snakemake Pipeline for Phage Genome Annotation & Pangenome/Clustering Analysis
===============================================================================
This module operates on the dereplicated, high-quality (HQ) viral contigs
produced upstream (see vMAGs_info.smk) and:
1. Annotates phage genes/functions using Pharokka.
2. Computes pairwise intergenomic similarities between all HQ phage genomes
   using VIRIDIC, used to define viral clusters (e.g. vOTUs) for the
   pangenome/taxonomy analyses reported in the manuscript.
"""

# Annotate phage genomes
checkpoint run_pharokka:
    """
    Checkpoint: Annotate all HQ viral genomes (genes, functions, taxonomy)
    using Pharokka, splitting output per-genome ("--split") for later use
    as gene lists in inStrain (see instrain.smk).
    """
    input:
        assembly = "../results/assembly/viral/all_HQ_viral_contigs.fasta"
    output:
        dir=directory("../results/vMAGs/annotations/pharokka/all_viruses")
    params:
        db="/work/FAC/FBM/DMF/pengel/general_data/mndiaye1/databases/pharokka_db"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "1h"
    threads:16
    conda:
        "envs/pharokka.yaml"
    log:
        "logs/vMAGs/pharokka_allgenomes.log"
    shell:
        "pharokka.py -i {input} -o {output} -d {params.db} -t {threads} --meta --split"



###############################  Viridic #############################################################
rule viridic_all_genomes:
    """
    Compute pairwise intergenomic similarities among all HQ phage genomes
    using the VIRIDIC singularity container, used downstream to define
    viral clusters/OTUs (vOTUs).
    """
    input:
        "../results/assembly/viral/all_HQ_viral_contigs.fasta"
    output:
        directory("../results/pangenomics/viruses/viridic")
    threads: 20
    log:
        "logs/pangenomics/phages/viridic_allGenomes.log"
    conda:
        "envs/viridic.yaml"
    params:
        viridic_sing="../resources/databases/containers/viridic_v1.1",
        # NOTE: absolute path required by VIRIDIC's singularity bind mounts;
        # update this if the pipeline is relocated to a different directory.
        abs_path="/work/FAC/FBM/DMF/pengel/general_data/mndiaye1/20241210_PHIFIDO_AmpliPhage_pipeline/workflow"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "12h"
    shell:
        "mkdir -p {output}; "
        "cd {params.viridic_sing}; "
        "in=$(basename {input}); "
        "indir=$(dirname {input}); "
        "singularity run -B \"{params.abs_path}/${{indir}}:/viridic/viridic_scripts/in\" -B \"{params.abs_path}/{output}:/viridic/viridic_scripts/out\" viridic_singularity_v1.1.simg projdir=/viridic/viridic_scripts/out in=/viridic/viridic_scripts/in/${{in}} ncor={threads}"


rule parse_viridic:
    """
    Parse VIRIDIC's intergenomic similarity matrix and combine it with
    inStrain-derived viral community composition data to produce the
    final viral clustering/pangenome summary table.
    """
    input:
        viridic="../results/pangenomics/viruses/viridic/",
        vircom="../results/inStrain/aggregated_data_sing"
    output:
        directory("../results/pangenomics/viruses/viridic_parsed")
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "10m"
    conda:
        "envs/base_R_env.yaml"
    params:
        script="scripts/pangenomics/parse_viridic_results.R"
    log:
        "logs/pangenomics/phages/parse_viridic.log"
    shell:
        "Rscript {params.script} -m {input.viridic}/04_VIRIDIC_out/sim_MA_genCol.csv -v {input.vircom}/vircom_data_filtered_updated.csv -o {output}"
