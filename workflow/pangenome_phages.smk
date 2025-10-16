# Annotate phage genomes
checkpoint run_pharokka:
    input:
        assembly = "../results/assembly/viral/all_viral_contigs.fasta"
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


rule defense_finder_viruses: 
    input:
        annots="../results/vMAGs/annotations/pharokka/all_viruses"
    output:
        dir=directory("../results/pangenomics/viruses/defense_finder"),
        tmp=temp("../results/pangenomics/viruses/all_viruses.faa")
    threads: 8
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "1h"
    log: 
        "logs/defense_finder/viruses_defense_finder.log"
    conda:
        "envs/defense_finder.yaml"
    params:
        models=config["DF_MODELS"]
    shell:
        "cat  {input.annots}/single_faas/*.faa > {output.tmp}; "
        "defense-finder run -o {output.dir} -w {threads} --models-dir {params.models} {output.tmp}"

rule split_pharokka_vOTU:
    input:
        annot="../results/vMAGs/annotations/pharokka/all_viruses",
        drep="../results/vMAGs/dereplication/dRep_summary_average.tsv"
    output:
        directory("../scratch_link/annotations_vOTU/{vOTU}")
    resources:
        account="pengel_beemicrophage",
        mem_mb= 5000,
        runtime = "20m"
    threads:1
    log:
        "logs/vMAGs/split_pharokka_{vOTU}.log"
    shell:
        "mkdir -p {output}; "
        "python scripts/annotations/pharokka_splitter.py -i {input.annot} -d {input.drep} -v {wildcards.vOTU}  -o {output}"


################################ Vcontact2 #############################################################
rule gene_2_genome:
    input:
        all_vprot = "../results/vMAGs/annotations/pharokka/all_viruses"
    output:
        all_prot="../results/pangenomics/viruses/Vcontact2/all_viral_proteins.faa",
        gene_2_genome = "../results/pangenomics/viruses/Vcontact2/gene_to_genome_drep.csv"
    threads: 1
    conda:
        "envs/mags_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "10m"
    log:
        "logs/vcontact/gene_to_genome.log"
    shell:
        "cat {input.all_vprot}/single_faas/*.faa > {output.all_prot}; "
        "python scripts/Viral_classification/gene2genome.py -p {output.all_prot} -o {output.gene_2_genome} -s 'Prodigal-FAA'"

# Run vContact2
rule run_vcontact:
    input:
        all_vprot = "../results/pangenomics/viruses/Vcontact2/all_viral_proteins.faa",
        gene_2_genome = "../results/pangenomics/viruses/Vcontact2/gene_to_genome_drep.csv"
    output:
        directory("../results/pangenomics/viruses/Vcontact2/vCONTACT_results")
    threads: 48
    params:
        condaenv=config["CONDA_ENV_VCONTACT2"]
    log:
        "logs/vcontact/run_vcontact2.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 512000,
        runtime= "2d"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.condaenv}
            vcontact2 -t {threads} --raw-proteins {input.all_vprot} --rel-mode 'Diamond' --proteins-fp {input.gene_2_genome} --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.condaenv}/bin/cluster_one-1.0.jar --output-dir {output} -e 'cytoscape' -e 'csv''
        """


###############################  Phylogeny #############################################################
rule viridic_all_genomes:
    input:
        "../results/assembly/viral/all_viral_contigs.fasta"
    output:
        directory("../results/pangenomics/viruses/viridic")
    threads: 20
    log:
        "logs/pangenomics/phages/viridic_allGenomes.log"
    conda:
        "envs/viridic.yaml"
    params:
        viridic_sing="../resources/databases/containers/viridic_v1.1",
        abs_path="/work/FAC/FBM/DMF/pengel/general_data/mndiaye1/20241210_PHIFIDO_AmpliPhage_pipeline/workflow"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "12h"
    shell:
        "mkdir -p {output}; "
        "cd {params.viridic_sing}; "
        "in=$(basename {input}); "
        "indir=$(dirname {input}); "
        "singularity run -B \"{params.abs_path}/${{indir}}:/viridic/viridic_scripts/in\" -B \"{params.abs_path}/{output}:/viridic/viridic_scripts/out\" viridic_singularity_v1.1.simg projdir=/viridic/viridic_scripts/out in=/viridic/viridic_scripts/in/${{in}} ncor={threads}"


rule viral_phylogeny: 
    input:
        viruses="../results/assembly/viral/all_viral_contigs.fasta"
    output:
        directory("../results/vMAGs/vMAGs_phylogeny")
    threads: 20
    log:
        "logs/vMAGs/viral_phylogeny.log"
    conda:
        "envs/phabox.yaml"
    params:
        db=config["PHABOX_DB"]
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "10h"
    shell:
        "phabox2 --task tree --dbdir {params.db} \
        --len 5000 \
        --outpth  {output} \
        --contigs {input} \
        --threads {threads} --tree Y --msa Y"