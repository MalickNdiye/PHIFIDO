
# Annotate phage genomes
checkpoint run_pharokka:
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

rule gene_2_genome:
    input:
        all_vprot = "../results/vMAGs/annotations/pharokka/all_viruses",
        phoster_ref= "../data/references/PHOSTER_references/PHOSTER_vOTUs.faa"
    output:
        all_prot="../results/pangenomics/viruses/Vcontact2/all_viral_proteins.faa",
        all_prot_and_ref="../results/pangenomics/viruses/Vcontact2/all_viral_proteins_PlusPhoster.faa",
        gene_2_genome = "../results/pangenomics/viruses/Vcontact2/gene_to_genome_drep.csv"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "10m"
    conda:
        "envs/bacphlip.yaml"
    log:
        "logs/vcontact/gene_to_genome.log"
    shell:
        "cat {input.all_vprot}/single_faas/*.faa > {output.all_prot}; "
        "cat {output.all_prot} {input.phoster_ref} > {output.all_prot_and_ref}; "
        "python scripts/pangenomics/gene2genome.py -p {output.all_prot_and_ref} -o {output.gene_2_genome} -s 'Prodigal-FAA'"


rule defense_finder_viruses: 
    input:
        annots="../results/pangenomics/viruses/Vcontact2/all_viral_proteins_PlusPhoster.faa"
    output:
        dir=directory("../results/pangenomics/viruses/defense_finder")
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
        "defense-finder run -o {output.dir} -w {threads} --models-dir {params.models} {input.annots}"


################################ Vcontact2 #############################################################

# Run vContact2
rule run_vcontact:
    input:
        all_vprot = "../results/pangenomics/viruses/Vcontact2/all_viral_proteins_PlusPhoster.faa",
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
        runtime= "1d"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.condaenv}
            vcontact2 -t {threads} --raw-proteins {input.all_vprot} --rel-mode 'Diamond' --proteins-fp {input.gene_2_genome} --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.condaenv}/bin/cluster_one-1.0.jar --output-dir {output} -e 'cytoscape' -e 'csv''
        """

rule analyse_vcontact:
    input:
        vcontact="../results/pangenomics/viruses/Vcontact2/vCONTACT_results",
        comm="../results/vMAGs/dereplication/dRep_summary_average.tsv",
        gen_info="../results/vMAGs/vMAGs_summary.tsv"
    output:
        "../results/pangenomics/viruses/Vcontact2/vcontact_analysis/"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "10m"
    conda: "envs/base_R_env.yaml"
    params: 
        "scripts/pangenomics/analyze_vcontact.R"
    log:
        "logs/vcontact/analyze_vcontact.log"
    shell:
        "Rscript scripts/pangenomics/analyze_vcontact.R -i {input.vcontact}/c1.ntw -d {input.drep}/vircom_data_filtered_updated.csv -g {input.gen_info} -o {output}"


###############################  Phylogeny #############################################################
rule viridic_all_genomes:
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
    input:
        viridic="../results/pangenomics/viruses/viridic/",
        vircom="../results/inStrain/aggregated_data"
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