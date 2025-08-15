import os


def get_files_commas(path, sep=",", remove_hidden=True):
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
    input:
        refs="../results/vMAGs/dereplication/dRep_average/dereplicated_genomes/"
    output:
        concat="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        stb="../results/inStrain/all_drep_average_vMAG_representatives.stb"
    log:
        "logs/instrain/concat_refs.log"
    threads: 2
    conda:
        "envs/drep_env.yaml"
    params:
        refs=lambda wildcards, input: get_files_commas(input[0], sep=" ")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "30m"
    shell:
        "cat {input.refs}/*.f* >> {output.concat}; "
        "parse_stb.py --reverse -f {params.refs}  -o {output.stb}"

rule generate_genelist_viral:
    input:
        ref=lambda wildcards: [
            os.path.join(
                "../results/vMAGs/annotations/pharokka/all_viruses/single_fastas",
                f"{g}.fasta"
            )
            for g in get_representative_genomes_average(wildcards)
        ]
    output:
        "../results/inStrain/all_drep_average_vMAG_representatives_cds.fna"
    log:
        "logs/instrain/generate_gene_list_viruses.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "10m"
    shell:
        "cat {input.ref} > {output}"

############################################# Mapping #################################################################

rule build_ref_index:
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

############################################# InStrain Profile ##########################################################
rule instrain_profile:
    input:
        bam="../scratch_link/mapping/mapdata/{sample}_bowtie_mapping.bam",
        ref="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        genL="../results/inStrain/all_drep_average_vMAG_representatives_cds.fna",
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
        "(inStrain profile {input.bam} {input.ref} -o {output.dir} --min_read_ani {params.min_ANI} -p {threads} -g {input.genL} -s {input.stb})2> {log}"

rule aggregate_inStrain:
    input:
        dir=expand("../results/inStrain/profiles/{sample}_profile/", sample=config["samples"])
    output:
        tab="../results/inStrain/aggregate_profile.tsv"
    shell:
        "touch {output.tab}; "