def get_files_commas(path, sep=",", remove_hidden=True):
    file_l=os.listdir(path)

    if remove_hidden:
        file_l=[f for f in file_l if not f.startswith(".")]

    file_l2=[]
    for f in file_l:
        file_l2.append(os.path.join(path,f))
    out=sep.join(file_l2)
    return(out)

rule get_stb_viral:
    input:
        refs="../results/vMAGs/dereplication/dRep_average/dereplicated_genomes/"
    output:
        concat="../results/vMAGs/dereplication/dRep_average/all_vMAGs.fasta",
        stb="../results/vMAGs/dereplication/dRep_average/all_vMAGss.stb"
    log:
        "logs/ref_db/concat_refs.log"
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
        ref="../results/pangenomes/P/annotations_reference/genes_fna"
    output:
        "../results/inStrain/all_P_RefGenomes_genes.fna"
    log:
        "logs/instrain/generate_gene_list_P.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "10m"
    shell:
        "cat {input.ref}/* > {output}"

# run instrain profile
rule instrain_profile:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam",
        ref="../results/reference_db_filtered/{type}/all_{type}_RefGenomes.fasta",
        genL="../results/inStrain/all_{type}_RefGenomes_genes.fna",
        stb="../results/reference_db_filtered/{type}/all_{type}_RefGenomes.stb" 
    output:
        dir=directory("../results/inStrain/profile_{type}/{sample}{type}_profile/")
    threads: 32
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/{sample}{type}_profile.log"
    benchmark:
        "logs/instrain/{sample}{type}_profile.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "10h"
    params:
        min_ANI=0.92
    shell:
        "(inStrain profile {input.bam} {input.ref} -o {output.dir} --min_read_ani {params.min_ANI} -p {threads} -g {input.genL} -s {input.stb})2> {log}"
