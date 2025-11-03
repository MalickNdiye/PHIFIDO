##############################################################################################################
# vMAGs info gathering: CheckV, Bacphlip, Phabox2
##############################################################################################################

# this rule runs checkV on the binned viral contigs
rule run_checkv:
    input:
        assembly = "../results/assembly/viral/all_viral_contigs.fasta"
    output:
        dir=directory("../results/vMAGs/QC/Checkv")
    params:
        db=config["CHECKV_DB"]
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "1h"
    threads:10
    conda:
        "envs/checkv.yaml"
    log:
        "logs/vMAGs/checkV.log"
    shell:
        "checkv end_to_end {input.assembly} {output.dir} -t {threads} -d {params.db}"

rule lifestyle:
    input:
        viruses="../results/assembly/viral/all_viral_contigs.fasta"
    output:
        directory("../results/vMAGs/lifestyle")
    threads: 10
    log:
        "logs/vMAGs/lifestyle.log"
    conda:
        "envs/bacphlip.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    shell:
        "bacphlip -i {input.viruses} --multi_fasta -f; "
        "mkdir -p {output}; "
        "mv {input.viruses}.BACPHLIP_DIR {input.viruses}.bacphlip* {input.viruses}.hmmsearch* {output}"

rule taxonomy:
    input:
        viruses="../results/assembly/viral/all_viral_contigs.fasta"
    output:
        directory("../results/vMAGs/taxonomy")
    threads: 10
    log:
        "logs/vMAGs/taxonomy.log"
    conda:
        "envs/phabox.yaml"
    params:
        db=config["PHABOX_DB"]
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    shell:
        "phabox2 --task phagcn --dbdir {params.db} --contigs {input.viruses} --outpth {output} --threads {threads}"

###############################  Phage-Host prediction #############################################################

rule assign_host:
    input:
        phageDB="../results/assembly/viral/{sample}_viral_contigs.fasta",
        spacers_db=config["SPACERS_DB"]
    output:
        blastout=temp("../results/phage_host_link/spacers_{sample}_blastout.txt")
    log:
        "logs/phage_host_link/CRIPR_match_{sample}.log"
    threads: 10
    conda:
        "envs/CrisprOpenDB.yaml"
    params:
        DB="resources/default_DBs/CrisprOpenDB"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "1h"
    shell:
        "blastn -query {input.phageDB} -task blastn-short -db {input.spacers_db}/mySpacersDB -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -out {output.blastout} -num_threads {threads}"

rule aggregate_assign_host:
    input:
        blastout=expand("../results/phage_host_link/spacers_{sample}_blastout.txt", sample=[s for grp in config["samples"].values() for s in grp.keys()])
    output:
        all_blastout="../results/phage_host_link/all_spacers_blastout.txt"
    log:
        "logs/phage_host_link/aggregate_assign_host.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "10m"
    shell:
        "echo -e 'Query\tSPACER_ID\tidentity\talignement_length\tmismatch\tgap\tq_start\tq_end\ts_start\ts_end\te_value\tscore' > {output.all_blastout}; "
        "tail -n +2 -q {input.blastout} >> {output.all_blastout}; "

rule parse_phage_host_links:
    input:
        blastout="../results/phage_host_link/all_spacers_blastout.txt",
        bacteria_metadata=config["PHOSTER_BACTERIA_METADATA"],
        spacers_metadata=config["SPACERS_METADATA"]
    output:
        "../results/phage_host_link/spacers_phage_host_links_summary.tsv"
    log:
        "logs/phage_host_link/parse_phage_host_links.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "10m"
    conda:
        "envs/base_R_env.yaml"
    shell:
        "Rscript scripts/phage_host_link/parse_phage_host_links.R -b {input.bacteria_metadata} -s {input.spacers_metadata} -a {input.blastout} -o {output}"

####################################  Aggregate vMAGs info #############################################################
rule aggregate_vMAGs_info:
    input:
        lifestyle="../results/vMAGs/lifestyle",
        taxonomy="../results/vMAGs/taxonomy",
        checkv="../results/vMAGs/QC/Checkv",
        phage_host="../results/phage_host_link/spacers_phage_host_links_summary.tsv",
        drep="../results/vMAGs/dereplication/dRep_summary_average.tsv"
    output:
        "../results/vMAGs/vMAGs_summary.tsv"
    log:
        "logs/vMAGs/aggregate_vMAGs_info.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "10m"
    threads: 1
    conda:
        "envs/base_R_env.yaml"
    shell:
        "Rscript scripts/vMAGs_handling/aggregate_vMAGs_info.R -l {input.lifestyle} -t {input.taxonomy} -c {input.checkv} -p {input.phage_host} -d {input.drep} -o {output}"


# rule to filter good quality vMAGs based on checkV results
rule filter_good_vMAGs:
    input:
        vMAGs_info="../results/vMAGs/vMAGs_summary.tsv",
        assembly="../results/assembly/viral/all_viral_contigs.fasta"
    output:
        "../results/assembly/viral/all_HQ_viral_contigs.fasta"
    log:
        "logs/vMAGs/filter_good_vMAGs.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "10m"
    threads: 1
    run:
        import pandas as pd
        from Bio import SeqIO

        vMAGs_info = pd.read_csv(input.vMAGs_info, sep="\t")
        good_quality_contigs = vMAGs_info[vMAGs_info['checkv_quality'].isin(['Medium-quality', 'High-quality', 'Complete'])]['genome'].tolist()

        with open(output[0], "w") as out_f:
            for record in SeqIO.parse(input.assembly, "fasta"):
                if record.id in good_quality_contigs:
                    SeqIO.write(record, out_f, "fasta")