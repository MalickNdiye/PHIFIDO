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

rule run_orthofinder:
    input:
        faas=expand("../results/pangenomics/bacteria/annotations/drammotate/cds/nucleotides/{bifido}_genes.fasta", bifido=config["all_bifidos"]),
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
        "orthofinder -f ${{faas}} -d -o {output.ortho_out} -n {params.name} -t {threads} -M msa"

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

rule run_defense_finder:
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

rule aggregate_DF_results:
    input:
        expand("../results/pangenomics/bacteria/defense_finder/{bifido}", bifido=config["BiCom"])
    output:
        "../results/pangenomics/bacteria/defense_finder/defense_finder_systems.tsv"
    log: "logs/defense_finder/aggregate_DF_results.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "1h"
    run:
        # join {bifido}_genes_defense_finder.tsv files to the path of each input
        strains=config["BiCom"]
        renamed_inputs = [f"{input}/{strains[i]}_genes_defense_finder_systems.tsv" for i in range(len(strains))] # TODO fix this
        print(renamed_inputs)
        df=concat_tables(renamed_inputs, skip=0, delimiter="\t", head=0, names=False, add={"genome": lambda wildcards: wildcards.bifido})
        df.to_csv(output[0], sep="\t", index=False)