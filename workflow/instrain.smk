"""
Snakemake Pipeline for Population Genomics via inStrain
========================================================
This module quantifies within-population genetic diversity of the
dereplicated viral representative genomes (vMAGs) across all metagenomic
samples. It:
1. Builds the reference database (concatenated representative vMAG FASTA,
   scaffold-to-bin file, and gene list) needed by inStrain.
2. Maps trimmed reads from each sample to this reference (Bowtie2).
3. Runs inStrain profile per sample to compute coverage, breadth, and
   nucleotide diversity metrics.
4. Aggregates per-sample profiles into a single dataset for downstream
   community/population analyses (e.g. input to VIRIDIC parsing in
   pangenome_phages.smk).
"""

import os


def get_files_commas(path, sep=",", remove_hidden=True):
    """
    Returns a comma-separated string of all (non-hidden) file paths in
    a directory. Utility used when a shell command needs a single
    delimited string rather than a list of paths.
    """
    file_l=os.listdir(path)

    if remove_hidden:
        file_l=[f for f in file_l if not f.startswith(".")]

    file_l2=[]
    for f in file_l:
        file_l2.append(os.path.join(path,f))
    out=sep.join(file_l2)
    return(out)

############################################# InStrain SetUp ##########################################################
rule get_stb_viral:
    """
    Concatenate all dRep-representative vMAG genomes into a single FASTA
    reference and generate the scaffold-to-bin (.stb) file inStrain needs
    to map contigs back to their genome of origin.
    """
    input:
        refs=get_avg_representative_fasta
    output:
        concat="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        stb="../results/inStrain/all_drep_average_vMAG_representatives.stb"
    log:
        "logs/instrain/concat_refs.log"
    threads: 2
    conda:
        "envs/drep_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "30m"
    shell:
        "cat {input.refs} >> {output.concat}; "
        "parse_stb.py --reverse -f {input.refs}  -o {output.stb}"


def get_pharokka_gbks(wildcards):
    """
    Resolve the Pharokka per-genome GenBank (.gbk) file paths for each
    dRep-representative vMAG, once both the run_pharokka and
    parse_dRep_viruses_average checkpoints have been evaluated.
    """
    checkpoint_output = checkpoints.run_pharokka.get(**wildcards).output.dir
    gbk_dir = os.path.join(checkpoint_output, "single_gbks")

    # Get list of representative genomes (depends on previous checkpoint)
    rep_genomes = get_representative_genomes_average(wildcards)

    # Build full paths
    gbk_paths = [os.path.join(gbk_dir, f"{g}.gbk") for g in rep_genomes]
    return gbk_paths


rule generate_genelist_viral:
    """
    Build the combined gene-list (.gbk) file for all representative vMAGs
    from their Pharokka annotations, renaming "locus_tag" to "gene" so the
    file is compatible with inStrain's expected gene-calling format.
    """
    input:
        single_genomes="../results/assembly/viral/single_genomes/",
        ref=get_pharokka_gbks
    output:
        "../results/inStrain/all_drep_average_vMAG_representatives_cds.gbk"
    log:
        "logs/instrain/generate_gene_list_viruses.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "10m"
    shell:
        "cat {input.ref} > {output}; "
        # replace "locus_tag" with "gene" in the gene list
        "sed -i 's/locus_tag/gene/g' {output}; "

############################################# Mapping #################################################################

rule build_ref_index:
    """
    Build a Bowtie2 index from the concatenated representative vMAG FASTA,
    used to map each sample's reads against the dereplicated viral database.
    """
    input:
        ref="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
    output:
        index=directory("../results/inStrain/drep_average_vMAG_representatives_index/")
    conda:
        "envs/map_env.yaml"
    threads: 8
    log:
        "logs/instrain/mapping/build_bowtie_index.log"
    params:
        basename="drep_average_vMAG_representatives_index",
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "1h"
    shell:
        "mkdir -p {output.index}; "
        "bowtie2-build {input.ref} {output.index}/{params.basename} --threads {threads}"

# Map reads of viral fraction to dereplicated phage genome database and reads of bacterial fraction to dereplicated bacterial genome database
rule MapReads:
    """
    Map each sample's trimmed paired-end reads to the dereplicated vMAG
    reference index using Bowtie2.
    """
    input:
        index="../results/inStrain/drep_average_vMAG_representatives_index/",
        R1="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2="../data/trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        sam=temp("../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.sam"),
    conda:
        "envs/map_env.yaml"
    threads: 8
    params:
        basename="drep_average_vMAG_representatives_index"
    log:
        "logs/instrain/mapping/{sample}_bowtie_mapping.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 200000,
        runtime= "2h"
    shell:
        "bowtie2 -x {input.index}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"

# generate bam file and count the number of reads mapped to each library
rule generate_bam:
    """
    Convert SAM to BAM and record the total number of reads in each
    sample's mapping library.
    """
    input:
        sam="../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.sam"
    output:
        bam="../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.bam",
        cf=temp("../results/inStrain/mapping/{sample}_library_count.tsv")
    conda:
        "envs/map_env.yaml"
    threads: 2
    log:
        "logs/instrain/mapping/{sample}_generateBAM.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "1h"
    shell:
        "samtools view -S -b {input.sam} > {output.bam}; "
        "touch {output.cf}; "
        "s=$(basename {output.bam}); "
        "lib=${{s%.*}}; "
        "count=$(samtools view -c {output.bam}); "
        'echo -e "${{lib}}\t${{count}}" >> {output.cf}'

rule mapping_stats:
    """
    Compute per-sample mapping statistics: total reads vs. reads mapped
    to the viral reference database.
    """
    input:
        bam="../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.bam"
    output:
        temp("../results/inStrain/mapping/mapping_stats_{sample}.tsv")
    conda:
        "envs/map_env.yaml"
    threads: 1
    log:
        "logs/instrain/mapping/mapping_stats_{sample}.log"
    shell:
        "total=$(samtools view -c {input.bam}); "
        "mapped=$(samtools view -c -F 4 {input.bam}); "
        "echo -e '{wildcards.sample}\t'${{mapped}}'\t'${{total}} > {output}"

rule concat_mapping_stats:
    """
    Concatenate per-sample mapping statistics into a single summary table
    (used in the pipeline's final "rule all" target).
    """
    input:
        expand("../results/inStrain/mapping/mapping_stats_{sample}.tsv", sample=[s for grp in config["samples"].values() for s in grp.keys()])
    output:
        "../results/inStrain/mapping/all_mapping_stats.tsv"
    threads: 1
    log:
        "logs/instrain/mapping/concat_mapping_stats.log"
    shell:
        # echo header first and then append all files together
        "echo -e 'Sample\tMapped_reads\tTotal_reads' > {output}; "
        "cat {input} >> {output}; "


############################################# InStrain Profile ##########################################################
rule instrain_profile:
    """
    Run inStrain profile on each sample's BAM against the dereplicated
    vMAG reference to compute per-genome and per-gene coverage, breadth,
    and nucleotide diversity metrics. Falls back to creating empty
    placeholder outputs (and a "*_nophage.done" flag) if a sample has no
    detectable phage reads, so the DAG can still complete.
    """
    input:
        bam="../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.bam",
        ref="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        genL="../results/inStrain/all_drep_average_vMAG_representatives_cds.gbk",
        stb="../results/inStrain/all_drep_average_vMAG_representatives.stb"
    output:
        dir=directory("../results/inStrain/profiles/{sample}_profile/")
    threads: 15
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/profiles/{sample}_profile.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "10h"
    params:
        min_ANI=0.92
    shell:
        "(inStrain profile {input.bam} {input.ref} -o {output.dir} --min_read_ani {params.min_ANI} -p {threads} -g {input.genL} -s {input.stb}) "
        "|| (mkdir -p {output.dir}/output "
        "&& touch {output.dir}/output/{wildcards.sample}_profile_genome_info.tsv "
        "&& touch {output.dir}/output/{wildcards.sample}_profile_gene_info.tsv "
        "&& touch {output.dir}/output/{wildcards.sample}_nophage.done)"


rule aggregate_inStrain_sing:
    """
    Aggregate all per-sample inStrain profiles into a single dataset,
    merging in sample metadata, dRep cluster assignments, and ANI
    comparisons (produces the input used by parse_viridic in
    pangenome_phages.smk).
    """
    input:
        dir=expand("../results/inStrain/profiles/{sample}_profile/", sample=[s for grp in config["samples"].values() for s in grp.keys()]),
        metadata="../data/metadata/sample_metadata.csv",
        drep="../results/vMAGs/dereplication/dRep_summary_single.tsv",
        ani="../results/vMAGs/vMAGs_ANI_comparison.txt"
    output:
        directory("../results/inStrain/aggregated_data_sing/")
    conda:
        "envs/base_R_env.yaml"
    params:
        gen_func="scripts/General_functions.R"
    threads: 1
    log:
        "logs/instrain/aggregate_profile.log"
    shell:
        "Rscript scripts/community_analysis/inStrain_aggregate.R -i {input.dir} -d {input.drep} -m {input.metadata} -g {params.gen_func} -a {input.ani} -o {output}"
