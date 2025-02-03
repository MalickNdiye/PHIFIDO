rule fastani_phages:
    input:
        viruses="../results/assembly/viral/single_genomes_list.txt"
    output:
        "../results/vMAGs/vMAGs_ANI_comparison.txt"
    threads: 15
    log:
        "logs/vMAGs/fastani.log"
    conda:
        "envs/fastani.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "30m"
    shell:
        "fastANI --ql {input.viruses} --rl {input.viruses} -t {threads} --fragLen 100 -o {output}; "

rule get_vOTU_order:
    input:
        drep="../results/vMAGs/dereplication/dRep_summary_single.tsv",
        ANI="../results/vMAGs/vMAGs_ANI_comparison.txt"
    output:
        temp("../scratch_link/{vOTU}_align_order.txt")
    resources:
        account="pengel_beemicrophage",
        mem_mb= 5000,
        runtime = "20m"
    threads:1
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/vMAGs/order_{vOTU}.log"
    shell:
        "Rscript scripts/vMAGs_handling/get_vOTU_order.R -d {input.drep} -a {input.ANI} -v {wildcards.vOTU} -o {output}"

rule mafft_VOTU:
    input:
        annot="../scratch_link/annotations_vOTU/{vOTU}",
        order="../scratch_link/{vOTU}_align_order.txt"
    output:
        "../results/alignments_vOTUs/{vOTU}_mafft_aln.fasta"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "1h"
    threads:5
    params:
        tmp="{vOTU}_all_genomes.fasta"
    log:
        "logs/aligments/maftt_{vOTU}.log"
    conda:
        "envs/mafft.yaml"
    shell:
        "file_list=$(cat {input.order}); "
        # replace all the ".gbk" by ".fasta" in the file list
        "file_list=$(echo $file_list | sed 's/\.gbk/.fasta/g'); "
        # replace all the "gbk" by "fna" in the file list
        "file_list=$(echo $file_list | sed 's/gbk/fna/g'); "
        "cat ${{file_list}} > {params.tmp}; "
        "mafft --auto --thread {threads} --adjustdirection {params.tmp}> {output} || cat {params.tmp} > {output}; " 
        "rm {params.tmp}"


rule pgv_VOTU:
    input:
        annot="../scratch_link/annotations_vOTU/{vOTU}",
        order="../scratch_link/{vOTU}_align_order.txt"
    output:
        directory("../scratch_link/pgv_vOTUs/{vOTU}_pgv")
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "1h"
    threads:5
    log:
        "logs/aligments/pgv_{vOTU}.log"
    conda:
        "envs/pgv.yaml"
    shell:
        "file_list=$(cat {input.order}); "
        "pgv-mummer  --gbk_resources ${{file_list}} -t {threads} \
           -o {output} --seqtype nucleotide --curve \
           --feature_track_ratio 0.15 --fig_track_height 0.7 --feature_linewidth 0.5 --feature_plotstyle bigarrow \
           --normal_link_color chocolate --inverted_link_color limegreen || mkdir -p {output}"

def get_vOTU_names(wildcards):
    import pandas

    file=checkpoints.parse_dRep_viruses_single.get(**wildcards).output[0]
    df=pandas.read_csv(file, sep="\t")
    vOTU=df["vOTU"].tolist()
    return expand("../results/alignments_vOTUs/{vOTU}_mafft_aln.fasta", vOTU=vOTU) + \
           expand("../scratch_link/pgv_vOTUs/{vOTU}_pgv", vOTU=vOTU)

def get_tageted_vOTU_names(wildcards):
    import pandas

    file=checkpoints.parse_dRep_viruses_single.get(**wildcards).output[0]
    df=pandas.read_csv(file, sep="\t")

    # remove all rows where checkv_quality is "Low-quality"
    df_filt = df[df.checkv_quality != "Low-quality"]
    # remove all rows where genome contains the pattern "Pc"
    df_filt = df_filt[~df_filt.genome.str.contains("Pc")]


    vOTU=df_filt["vOTU"].tolist()
    return expand("../results/alignments_vOTUs/{vOTU}_mafft_aln.fasta", vOTU=vOTU) + \
           expand("../scratch_link/alignments_vOTUs/{vOTU}_pgv", vOTU=vOTU)

rule gather_vOTU_aln:
    input: 
        get_tageted_vOTU_names
    output:
        "../results/alignments_vOTUs/aln_vOTUs.done"
    shell:
        "touch {output}"

