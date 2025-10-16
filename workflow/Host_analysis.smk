################################## Bacterial genome general Analysis ##################################


rule checkm_QC:
    input:
        assembly="../data/references/Bifido_genomes/"
    output:
        dir=directory("../results/pangenomics/bacteria/checkm/"),
        file="../results/pangenomics/bacteria/checkm/checkm_QC_stats.csv"
    log:
        "logs/Bacteria_genomes/checkm.log"
    threads: 8
    params:
        db=config["CHECKM_DB"]
    conda:
        "envs/checkm_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    shell:
        "export CHECKM_DATA_PATH={params.db}; "
        "checkm lineage_wf {input.assembly} {output.dir} -x fasta -t {threads}; "
        "checkm qa {output.dir}/lineage.ms {output.dir} -o 2 -f {output.file} --tab_table"

rule classify_gtdbtk:
    input:
        filtered_mags="../data/references/Bifido_genomes/"
    output:
        class_out=directory("../results/pangenomics/bacteria/gtdbtk_classification/")
    log:
        "logs/pangenomics/bacteria/gtdbtk_classification.log"
    threads: 8
    conda:
        "envs/gtdb-tk.yaml"
    params:
        db=config["GTDB"]
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 150000,
        runtime= "3h"
    shell:
        "export GTDBTK_DATA_PATH={params.db}; "
        "gtdbtk classify_wf --genome_dir {input.filtered_mags} --extension fasta --skip_ani_screen --out_dir {output.class_out} --cpus {threads}"


rule Parse_bacterial_genomic_info:
    input:
        genomic_info="../data/metadata/PHOSTER_bacteria_genomes_info.csv",
        checkm="../results/pangenomics/bacteria/checkm/checkm_QC_stats.csv",
        gtdb="../results/pangenomics/bacteria/gtdbtk_classification/"
    output:
        "../results/pangenomics/bacteria/bacterial_genomic_info.csv"
    log: 
        "logs/pangenomics/bacteria/Parse_bacterial_genomic_info.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "1h"
    conda: "envs/base_R_env.yaml"
    shell:
        "Rscript scripts/pangenomics/Parse_bacterial_genomic_info.R -g {input.genomic_info} -c {input.checkm} -d {input.gtdb}/gtdbtk.bac120.summary.tsv -o {output}"


################################ Annotation of Bacterial Genomes #####################################################
rule dram_annotate_genomes:
    input:
        dram_config =config["DRAM_CONFIG"],
        assembly =  "../data/references/Bifido_genomes/{bifido}.fasta"
    output:
        dram_annotations = directory("../results/pangenomics/bacteria/annotations/drammotate/{bifido}/annotations")
    threads: 8
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    log: "logs/annotation/drammotate_{bifido}.log"
    benchmark: "logs/annotation/drammotate_{bifido}.benchmark"
    conda: "envs/dram.yaml"
    shell:
        """
        rm -rf {output.dram_annotations}
        DRAM-setup.py import_config --config_loc {input.dram_config} 
        DRAM.py annotate -i '{input.assembly}' -o {output.dram_annotations} --threads {threads} 
        """

rule move_genes_sam:
    input:
        fasta="../results/pangenomics/bacteria/annotations/drammotate/{bifido}/annotations",
    output:
        faa="../results/pangenomics/bacteria/annotations/drammotate/cds/proteins/{bifido}_genes.faa",
        fna="../results/pangenomics/bacteria/annotations/drammotate/cds/nucleotides/{bifido}_genes.fasta",
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "1h"
    log: "logs/annotation/move_genes_{bifido}.log"
    shell:
        "cp {input.fasta}/genes.faa {output.faa}; "
        "cp {input.fasta}/genes.fna {output.fna}"

######################################### Pangenome Analysis & Phylogeny ##########################################################
rule run_orthofinder:
    input:
        faas=expand("../results/pangenomics/bacteria/annotations/drammotate/cds/proteins/{bifido}_genes.faa", bifido=config["all_bifidos"]),
    output:
        ortho_out=directory("../results/pangenomics/bacteria/Orthofinder/Bifidos")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/pangenomes/run_orthofinder.log"
    params:
        name="orthofinder_results",
        ulim=2000000,
    threads: 20
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "5h"
    shell:
        "faas=$(echo {input.faas} | tr ' ' '\n' | xargs -n 1 dirname | sort | uniq); "
        "orthofinder -f ${{faas}} -o {output.ortho_out} -n {params.name} -t {threads} -M msa"

# Trim MSA for position with more than 50% gaps
rule trim_msa:
    input:
        "../results/pangenomics/bacteria/Orthofinder/Bifidos"
    output:
        trimmed_msa="../results/pangenomics/bacteria/Orthofinder/Bifido_OG_msa_clean.fa"
    threads: 4
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "1h"
    log: "logs/pangenomes/trim_msa.log"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -in {input}/Results_orthofinder_results/MultipleSequenceAlignments/SpeciesTreeAlignment.fa -out {output.trimmed_msa} -clustal -gt 0.5"

rule make_phylogeny:
    input:
        "../results/pangenomics/bacteria/Orthofinder/Bifido_OG_msa_clean.fa"
    output:
        directory("../results/pangenomics/bacteria/phylogeny/species_tree/")
    threads: 25
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "5h"
    log: "logs/phylogeny/make_species_phylogeny.log"
    conda:
        "envs/iqtree.yaml"
    params:
        outgroup="Ga0098206_genes,Bombiscardovia_coagulans_DSM_22924_genes,Bifidobacterium_xylocopae_XV2_genes"
    shell:
        "mkdir -p {output}; "
        "iqtree3 -s {input} -m LG+F+I+G4 -bb 1000 -pre {output}/Bifido_species_tree -o {params.outgroup} -T {threads}"


rule get_PD:
    input:
        tree="../results/pangenomics/bacteria/phylogeny/species_tree/"
    output:
        "../results/pangenomics/bacteria/phylogeny/Bifido_species_tree_PD.csv"
    log: 
        "logs/phylogeny/get_PD.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "10m"
    conda: 
        "envs/base_R_env.yaml"
    threads: 1
    shell:
        "Rscript scripts/pangenomics/get_PD.R -i {input.tree}/Bifido_species_tree.treefile -o {output}"


################################# Defense Systems Identification ######################################################
rule defense_finder_bacteria:
    input:
        bifido="../results/pangenomics/bacteria/annotations/drammotate/cds/proteins/{bifido}_genes.faa"
    output:
        directory("../results/pangenomics/bacteria/defense_finder/{bifido}")
    threads: 8
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "1h"
    log: "logs/defense_finder/defense_finder_{bifido}.log"
    conda:
        "envs/defense_finder.yaml"
    params:
        models=config["DF_MODELS"]
    shell:
        "defense-finder run -o {output} -w {threads} --models-dir {params.models} {input.bifido}"

rule aggregate_DF_bacteria:
    input:
        expand("../results/pangenomics/bacteria/defense_finder/{bifido}", bifido=config["BiCom"])
    output:
        "../results/pangenomics/bacteria/defense_finder/defense_finder_systems.csv"
    log: "logs/defense_finder/aggregate_DF_results.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "1h"
    run:
        from pathlib import Path
        strains=list(config["BiCom"].keys())
        file_list = [str(file) for directory in input for file in Path(directory).glob("*_genes_defense_finder_systems.tsv")]
        df=concat_tables(file_list, skip=0, delimiter="\t", head=0, names=False, add={"genome": strains})
        df.to_csv(output[0], sep="\t", index=False)