# Include this file from dataset/Isoform/Snakefile:
# include: "../../local/rules/mad_ltmm.smk"
# Put the accompanying R script in ../../local/src/.
# Dry run: snakemake -np --cores 1 statistics_v2/mad_ltmm/run1/isoform_MAD_diff.rnk
# Choose a NEW run label for later reruns; the script refuses nonempty output folders.
rule plot_top_mad_transcripts_ltmm:
    input:
        genes="filtering/genes/GEP.count.exp_filter.ltmm.gz",
        transcripts="filtering/transcripts/GEP.count.exp_filter.ltmm.gz",
        mapping="statistics_v2/genes_20260911_104750/tx2gene_used.tsv",
        gtf="/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/46/primary_assembly.annotation.gtf",
        metadata="metadata.txt",
        raw_transcripts="filtering/transcripts/GEP.count.gz",
        rscript="../../local/src/top_mad_transcripts_ltmm.R"
    output:
        rnk="statistics_v2/mad_ltmm/{run}/isoform_MAD_diff.rnk",
        report="statistics_v2/mad_ltmm/{run}/analysis_report.txt",
        stages="statistics_v2/mad_ltmm/{run}/stage_counts.tsv",
        status="statistics_v2/mad_ltmm/{run}/STATUS.txt",
        genes="statistics_v2/mad_ltmm/{run}/gene_MAD.tsv",
        transcripts="statistics_v2/mad_ltmm/{run}/transcript_MAD.tsv",
        table="statistics_v2/mad_ltmm/{run}/transcript_gene_MAD_audit.tsv",
        all_mad="statistics_v2/mad_ltmm/{run}/all_MAD_diff.tsv",
        gene_scores="statistics_v2/mad_ltmm/{run}/protein_coding_gene_scores.tsv",
        input_list="statistics_v2/mad_ltmm/{run}/EnrichR_input_list.txt",
        background_list="statistics_v2/mad_ltmm/{run}/EnrichR_background_list.txt"
    params:
        outdir=lambda wildcards: "statistics_v2/mad_ltmm/" + wildcards.run
    threads: 1
    resources:
        mem_mb=8000
    shell:
        """
        Rscript --vanilla {input.rscript:q} {params.outdir:q} \
            {input.genes:q} {input.transcripts:q} {input.mapping:q} \
            {input.gtf:q} {input.metadata:q} {input.raw_transcripts:q}
        """


rule build_gene_level_mad_table_run27:
    input:
        gene_mad="statistics_v2/mad_ltmm/run27/gene_MAD.tsv",
        transcript_audit="statistics_v2/mad_ltmm/run27/transcript_gene_MAD_audit.tsv",
        raw_genes="filtering/genes/GEP.count.gz",
        metadata="metadata.txt",
        gtf="/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/46/primary_assembly.annotation.gtf",
        rscript="../../local/src/gene_level_mad_table.R"
    output:
        "statistics_v2/mad_ltmm/run27/gene_transcript_MAD_audit.tsv"
    params:
        min_samples=27
    threads: 1
    resources:
        mem_mb=8000
    shell:
        """
        Rscript --vanilla {input.rscript:q} {output:q} \
            {input.gene_mad:q} {input.transcript_audit:q} \
            {input.raw_genes:q} {input.metadata:q} {input.gtf:q} \
            {params.min_samples:q}
        """
