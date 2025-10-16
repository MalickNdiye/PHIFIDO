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
        refs="../results/vMAGs/dereplication/dRep_average"
    output:
        concat="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        stb="../results/inStrain/all_drep_average_vMAG_representatives.stb"
    log:
        "logs/instrain/concat_refs.log"
    threads: 2
    conda:
        "envs/drep_env.yaml"
    params:
        refs=lambda wildcards, input: get_files_commas(input[0] + "/dereplicated_genomes", sep=" ")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "30m"
    shell:
        "cat {input.refs}/dereplicated_genomes/*.f* >> {output.concat}; "
        "parse_stb.py --reverse -f {params.refs}  -o {output.stb}"

def get_pharokka_gbks(wildcards):
    checkpoint_output = checkpoints.run_pharokka.get(**wildcards).output.dir
    gbk_dir = os.path.join(checkpoint_output, "single_gbks")

    # Get list of representative genomes (depends on previous checkpoint)
    rep_genomes = get_representative_genomes_average(wildcards)

    # Build full paths
    gbk_paths = [os.path.join(gbk_dir, f"{g}.gbk") for g in rep_genomes]
    return gbk_paths

rule generate_genelist_viral:
    input:
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

rule mapping_stats:
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

rule aggregate_inStrain:
    input:
        dir=expand("../results/inStrain/profiles/{sample}_profile/", sample=[s for grp in config["samples"].values() for s in grp.keys()]),
        metadata="../data/metadata/sample_metadata.csv",
        drep="../results/vMAGs/dereplication/dRep_summary_average.tsv"
    output:
        directory("../results/inStrain/aggregated_data/")
    conda:
        "envs/base_R_env.yaml"
    params:
        gen_func="scripts/General_functions.R"
    threads: 1
    log:
        "logs/instrain/aggregate_profile.log"
    shell:
        "Rscript scripts/community_analysis/inStrain_aggregate.R -i {input.dir} -d {input.drep} -m {input.metadata} -g {params.gen_func} -o {output}"


############################################# InStrain Compare #################################################################
rule instrain_compare:
    input:
        IS=expand(
            "../results/inStrain/profiles/{sample}_profile/",
            sample=[
                s
                for grp in config["samples"].values()
                for s in grp.keys()
                if all(excl not in s for excl in ["Blank", "none", "CTRL"])
            ]
        ),
        ref="../results/inStrain/all_drep_average_vMAG_representatives.fasta",
        stb="../results/inStrain/all_drep_average_vMAG_representatives.stb"
    output:
        directory("../results/inStrain/compare")
    threads: 25
    conda:
         "envs/inStrain.yaml"
    log:
        "logs/instrain/compare/compare.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 600000,
        runtime= "7h"
    shell:
        "inStrain compare -i {input.IS} -o {output} -p {threads} -s {input.stb} --database_mode --store_mismatch_locations"