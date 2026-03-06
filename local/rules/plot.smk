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