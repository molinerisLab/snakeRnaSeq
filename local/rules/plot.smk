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

########################
### Rules for krona  ###
########################


rule krona_txt:
    input:
        "breports_filtered/{sample}.breport",
    output:
        "b_krona_txt/{sample}.b.krona.txt",
    shell:
        """
        kreport2krona.py -r {input} -o {output} --no-intermediate-ranks
    """


rule krona_html:
    input:
        "b_krona_txt/{sample}.b.krona.txt",
    output:
        "krona_html/{sample}.krona.html",
    shell:
        """
        ktImportText {input} -o {output}
        """

rule common_taxa_in_samples:
    input:
        expand("boutputs_filtered/{sample}.braken", sample=SAMPLES),
    output:
        "results/common_taxa.csv",
        directory("results"),
    params:
        samples=SAMPLES,
    shell:
        """
        mkdir -p results 
        Rscript ../../local/src/common_species.R
        """


rule reads_human_contam_classified:
    output:
        "classified_vs_human_contaminant_barplot_normalized.png",
        "classified_vs_human_contaminant_barplot_percentage.png",
    shell:
        """
        Rscript /home/molinerislab/NeriMetagenome/workflow/src/plot_reads.R
        """

rule extract_genes_heatmap:
    input:
        "top1000_genes.rds" # Placeholder
    output:
        "top1000_genes.txt"
    shell:
        """
        Rscript ../../local/src/extract_genes_from_heatmap_rds.R {input} {output}
        """


rule extract_fasta:
    input:
        rds="top1000_genes.rds", # Placeholder
        fasta="transcriptome.fa" # Placeholder
    output:
        "extracted_transcripts.fasta"
    shell:
        """
        Rscript ../../local/src/extract_fasta_from_rds.R {input.rds} {output} {input.fasta}
        """

rule generate_depth_stats:
    input:
        kreports=expand("kreports/{sample}.k2report", sample=SAMPLES)
    output:
        seq="plots/sequencing_depth.csv",
        unmap="plots/unmapped_depth.csv"
    shell:
        """
        mkdir -p plots
        echo "name,sequencing_depth" > {output.seq}
        echo "name,unmapped_depth" > {output.unmap}
        for f in {input.kreports}; do
            name=$(basename $f .k2report)
            u=$(awk '$6 == "U" {{print $2}}' $f)
            r=$(awk '$6 == "R" {{print $2}}' $f)
            if [ -z "$u" ]; then u=0; fi
            if [ -z "$r" ]; then r=0; fi
            total=$((u+r))
            echo "${{name}},${{total}}" >> {output.seq}
            echo "${{name}},${{u}}" >> {output.unmap}
        done
        """
rule merge_kraken_stats:
    input:
        seq="plots/sequencing_depth.csv",
        unmap="plots/unmapped_depth.csv",
        kreports=expand("kreports/{sample}.k2report", sample=SAMPLES)
    output:
        "plots/master_kraken_stats.csv"
    params:
        kreport_dir="kreports"
    shell:
        """
        python ../../local/src/merge_kraken_stats.py \
            --sequencing-csv {input.seq} \
            --unmapped-csv {input.unmap} \
            --kreport-dir {params.kreport_dir} \
            --output-csv {output}
        """

rule plot_depth:
    input:
        "plots/sequencing_depth.csv"
    output:
        "plots/sequencing_depth_overview.jpg",
        "plots/sequencing_depth_outliers_readable.jpg",
        "plots/sequencing_depth_density_16x9.jpg",
        "plots/sequencing_depth_violin_16x9.jpg"
    shell:
        """
        Rscript ../../local/src/plot_depth.R {input} plots/
        """
rule plot_stacked_depth:
    input:
        unmap="plots/unmapped_depth.csv",
        seq="plots/sequencing_depth.csv"
    output:
        "plots/unmapped_reads_overview.jpg",
        "plots/unmapped_reads_outliers_16x9.jpg",
        "plots/mapped_unmapped_percentage.jpg",
        "plots/top25_unmapped_percentage.jpg"
    shell:
        """
        Rscript ../../local/src/plot_stacked_depth.R {input.unmap} {input.seq} plots/
        """
rule plot_kraken_unmapped:
    input:
        "plots/master_kraken_stats.csv"
    output:
        "plots/stacked_kraken_unmapped_overview.jpg",
        "plots/stacked_kraken_unmapped_overview_linear.jpg",
        "plots/stacked_kraken_unmapped_outliers_16x9.jpg",
        "plots/stacked_kraken_unmapped_outliers_16x9_linear.jpg"
    shell:
        """
        Rscript ../../local/src/plot_kraken_unmapped.R {input} plots/
        """
rule plot_pct:
    input:
        "plots/master_kraken_stats.csv"
    output:
        "plots/stacked_unmapped_kraken_percent_with_9606.png"
    shell:
        """
        Rscript ../../local/src/plot_pct.R {input} plots/
        """
rule plot_bracken_overview:
    input:
        frac="bracken_merged_abbundances.frac.txt",
        num="bracken_merged_abbundances.num.txt"
    output:
        "plots/top10_species_overall_stacked.png"
    shell:
        """
        Rscript ../../local/src/plot.R {input.frac} {input.num} 10 0.9 plots/
        """
rule plot_bracken_detailed:
    input:
        "bracken_merged_abbundances.num.txt"
    output:
        "plots/heatmap_top30_species.png",
        "plots/staphylococcus_aureus_marginal_scatter.png"
    shell:
        """
        Rscript ../../local/src/bracken_plots.R {input} "Staphylococcus aureus" 30 plots/
        """
