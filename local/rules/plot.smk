rule plot_read_distribution:
    input:
        all_gene_body_cov=expand(
            "rseqc/{sample}.read_distribution.txt",
            sample=SAMPLES
        ),
        rscript="../../local/src/plotReadDistributionSummary.R"
    output:
        "rseqc/reads_distribution_summary.pdf",
    shell:"""
        {CONDA_ACTIVATE} rstudio_Rv4.0.3; 
        Rscript --vanilla {input.rscript} rseqc
        """


rule generate_gene_counts:
    input:
        h5_files=expand("kallisto_PsiCLASS_idx_combined/{sample}/abundance.h5", sample=SAMPLES),
    output:
        counts="tximport_genes_counts.tsv",
        tpm="tximport_genes_tpm.tsv",
        gene_map="gene_id_to_symbol.tsv"
    shell:
        """
        Rscript --vanilla ../../local/src/run_tximport_gene.R
        """


rule plot_top_mad_transcripts:
    input:
        genes_tmm="tximport_genes_counts.tmm.tsv",
        tx_tmm="transcripts_counts.tmm.tsv",
        anno="gene_id_to_symbol.tsv",
        rscript="../../local/src/top_mad_transcripts.R"
    output:
        rnk="isoform_switching.rnk",
        top10_tx="plots/top10_highest_transcript_MAD.csv",
        top1000_switch="plots/top1000_isoform_switching.csv",
        all_mad="plots/all_MAD_diff.csv",
        enrichr_in="plots/EnrichR_input_list.txt",
        enrichr_bg="plots/EnrichR_background_list.txt"
    shell:
        """
        mkdir -p plots
        Rscript --vanilla {input.rscript}
        """