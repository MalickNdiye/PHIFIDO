rule dram_annotate_genomes:
    input:
        dram_config =config["DRAM_CONFIG"],
        assembly =  expand("../data/references/Bifido_genomes/{bifido}.fasta", bifido=config["all_bifidos"])
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
        fna="../results/pangenomics/bacteria/annotations/drammotate/cds/nucleotides/{bifido}_genes.fna",
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    log: "logs/annotation/move_genes_{bifido}.log"
    shell:
        "cp {input.fasta}/genes.faa {output.faa}; "
        "cp {input.fasta}/genes.fna {output.fna}"

rule run_orthofinder:
    input:
        faas=expand("../results/pangenomics/bacteria/annotations/drammotate/cds/nucleotides/{bifido}_genes.fna", bifido=config["all_bifidos"]),
    output:
        ortho_out=directory("../results/pangenomics/bacteria/Orthofinder/Bifidos")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/pangenomes/run_orthofinder.log"
    params:
        name="orthofinder_results",
        ulim=2000000,
        direc="../results/pangenomics/bacteria/annotations/drammotate/cds/nucleotides_bifido/"
    threads: 30
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "12h"
    shell:
        "mkdir -p {params.direc}; "
        "cp {input.faas} {params.direc}; "
        "ulimit -n {params.ulim}; "
        "orthofinder -f {params.direc} -d -o {output.ortho_out} -n {params.name} -t {threads} -M msa; "
        "rm -rf {params.direc}"

rule make_phylogeny:
    input:
        "../results/pangenomics/bacteria/Orthofinder/Bifidos"
    output:
        directory("../results/pangenomics/bacteria/phylogeny/species_tree/")
    threads: 15
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "5h"
    log: "logs/phylogeny/make_species_phylogeny.log"
    conda:
        "envs/iqtree.yaml"
    params:
        outgroup="Ga0098206_genes"
    shell:
        "mkdir -p {output}; "
        "iqtree -s {input}/Results_orthofinder_results/MultipleSequenceAlignments/SpeciesTreeAlignment.fa -nt {threads} -m MFP -bb 1000 -pre {output}/species_tree -o {params.outgroup}"