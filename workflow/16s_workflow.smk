
rule preoprocess_raw_reads:
    input:
        raw_reads = lambda wildcards: config["Raw_16s_reads"][wildcards.Experiment]
    output:
        results = directory("../PacBio_16s/results/prepocessing/{Experiment}")
    conda:
        "envs/r_momsane.yaml"
    log:
        "../PacBio_16s/logs/preprocessing_{Experiment}.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "2h"
    threads: 12
    params:
        F= config["primers"]["F"],
        R= config["primers"]["R"],
        minlen=1400,
        maxlen=1600,
        maxEE=3,
        out_preproc="../PacBio_16s/results/prepocessing/{Experiment}/reads",
        out_plots="../PacBio_16s/results/prepocessing/{Experiment}/plots"
    shell:
        """
        Rscript --vanilla ../PacBio_16s/scripts/02_preprocessing.R \
            {input.raw_reads} \
            {params.F} \
            {params.R} \
            {params.minlen} \
            {params.maxlen} \
            {params.maxEE} \
            {params.out_preproc} \
            {params.out_plots}
        """

rule denoising:
    input:
        preproc="../PacBio_16s/results/prepocessing/{Experiment}"
    output:
        results=directory("../PacBio_16s/results/denoising/{Experiment}")
    conda:
        "envs/r_momsane.yaml"
    log:
        "../PacBio_16s/logs/denoising_{Experiment}.log"
    threads: 12
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "5h"
    params:
        processed_fastq="../PacBio_16s/results/prepocessing/{Experiment}/reads/trimmed_filtered_reads",
        readcounts="../PacBio_16s/results/prepocessing/{Experiment}/reads/read_count_before_after.tsv",
        errModel=config["Error_function"],
        db=config["Bifido_16s_seqs"],
        maxraref=50000,
        maxReads=1000000,
        maxBases=10000000000,
        out_denois="../PacBio_16s/results/denoising/{Experiment}/denois_outputs",
        out_plots="../PacBio_16s/results/denoising/{Experiment}/denois_plots",
    shell:
        """
        Rscript --vanilla ../PacBio_16s/scripts/04_denoising.R \
            {params.processed_fastq} \
            {params.readcounts} \
            {params.maxReads} \
            {params.errModel} \
            {params.maxBases} \
            {params.db} \
            F \
            {params.maxraref} \
            {params.out_denois} \
            {params.out_plots}
        """

rule Merge_ASV_tabs:
    input:
        expand("../PacBio_16s/results/denoising/{Experiment}", Experiment=config["Raw_16s_reads"].keys())
    output:
        directory("../PacBio_16s/results/merged_ASV_tables")
    conda:
        "envs/r_momsane.yaml"
    log:
        "../PacBio_16s/logs/merge_ASV_tables.log"
    threads: 4
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "1h"
    params:
        input_tables=expand("../PacBio_16s/results/denoising/{Experiment}/denois_outputs/ASV_samples_table_noChim.rds", Experiment=config["Raw_16s_reads"]),
        Tax_db_toGenus=config["Tax_db_toGenus"],
        Tax_db_toSpecies=config["Tax_db_toSpecies"]
    shell:
        "Rscript --vanilla ../PacBio_16s/scripts/Merge_ASV_tables.R {params.input_tables} \
         {params.Tax_db_toGenus} \
         {params.Tax_db_toSpecies} \
         {output}"

rule format_16s_qpcr:
    input:
        qPCR_data=config["qPCR_data"],
        std_curves=config["std_curves"],
        metdata="../data/metadata/sample_metadata.csv"
    output:
        directory("../PacBio_16s/results/qPCR_formatted_{Experiment}")
    conda:
        "envs/base_R_env.yaml"
    log:
        "../PacBio_16s/logs/format_qPCR_data_{Experiment}.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "20m"
    params:
        qPCR_data=config["qPCR_data"]
    shell:
        "Rscript --vanilla ../PacBio_16s/scripts/format_qPCR_data.R {input.metdata} \
        {input.std_curves} \
        {input.qPCR_data} \
        {wildcards.Experiment} \
        {output}"
    
rule format_community_data:
    input:
        ASV_tables="../PacBio_16s/results/merged_ASV_tables",
        metadata="../data/metadata/sample_metadata.csv",
        qPCR_data="../PacBio_16s/results/qPCR_formatted_{Experiment}",
    output:
        directory("../PacBio_16s/results/community_data_formatted/{Experiment}")
    conda:
        "envs/r_momsane.yaml"
    log:
        "../PacBio_16s/logs/format_community_data_{Experiment}.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "20m"
    params:
        ASV_table="../PacBio_16s/results/merged_ASV_tables/ASV_table_merged.rds",
        taxonomy="../PacBio_16s/results/merged_ASV_tables/ASV_Taxonomy_sp.RDS",
        metadata="../data/metadata/sample_metadata.csv"
    shell:
        "Rscript --vanilla ../PacBio_16s/scripts/Format_16s_community_data_final.R -a {params.ASV_table} \
         -t {params.taxonomy} \
         -m {params.metadata} \
         -e {wildcards.Experiment} \
         -o {output}"