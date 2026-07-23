# rule plot_read_distribution:
#     input:
#         all_gene_body_cov=expand(
#             "rseqc/{sample}.read_distribution.txt",
#             sample=SAMPLES
#         ),
#         rscript="../../local/src/plotReadDistributionSummary.R"
#     output:
#         "rseqc/reads_distribution_summary.pdf",
#     shell:"""
#         {CONDA_ACTIVATE} rstudio_Rv4.0.3; 
#         Rscript --vanilla {input.rscript} rseqc
#         """

####################################
### Rules for barplot generation ###
####################################

rule barplot:
   input:
      bkrona="b_krona_txt/{sample}.b.krona.txt"
   output:
      bkrona_phylum="barplot/{sample}.barplot.txt"
   shell:
     """
      mkdir -p barplot/	
      bawk '$1!=0 {{print $4,$1}}' {input.bkrona}  | sort |stat_base -g -t > {output.bkrona_phylum}
     """

rule order_level_barplot:
    input:
        barplot_dir = "barplot"
    output:
        pdf = "plots/order_level_barplot.pdf"
    params:
        script = "../../local/src/bar_plot.R",
        top_n = 20
    shell:
        """
        Rscript {params.script} \
            {input.barplot_dir} \
            {output.pdf} \
            {params.top_n}
        """

        
# rule common_taxa_in_samples:
#     params:
#         samples=SAMPLES
#     shell:
#         """
#         Rscript ../../local/src/bar_plot.R
#         """


# # rule reads_human_contam_classified:
# #     output:
# #             "Plots/classified_vs_human_contaminant_barplot_normalized.png", 
# #             "Plots/classified_vs_human_contaminant_barplot_percentage.png"       
# #     params:
# #         split_by=config["split_by"]
# #     shell: 
# #         """
# #         Rscript src/plot_reads.R
# #         """

rule top_taxa_per_sample:
    input:
        abundance="Reports/bracken_merged_abundances.num.txt.xlsx",
        reads="Reports/reads_summary.csv",
        metadata="metadata.txt",
        config="config.yaml"
    output:
        "Plots/barplot_CPM.png",
        "Plots/barplot_CPM_percentage.png"
    shell:
        """
        Rscript src/plot_sample_species.R
        """


# rule combined_analysis:
#     output:
#         [
#             "top10_species_per_sample.csv",
#             "top10_species_overall.csv",
#             "top10_species_per_sample_barplot.png", 
#             "plots/all_patients_top_species_faceted.png",
#             "plots/top10_species_reads_faceted.png"
#         ]
#     shell:
#         """
#         Rscript src/combined.R
#         """

#==================================================
# rule extract_reads_csv: 
#     input:
#         expand("fastq/{sample}_R1.fastq.gz", sample=SAMPLES)
#     output:
#         "reads_summary.csv"
#     script:
#         "../../local/src/read_number.py"



rule rel_abundance_barplot:
    input:
        relative="abundances.cleaned.relative_all.tsv",
        metadata="metadata_full.txt"
    output:
        barplot="plots/stacked_barplot.pdf"
    script:
        "../../local/src/rel_barplot1.R"


rule all_candidate_species_abundance:
    input:
        pdf=f"plots/{config['CONTRAST']}_candidate_species_abundance.pdf",
        selected=f"plots/{config['CONTRAST']}_selected_candidate_species.tsv",
        long=f"plots/{config['CONTRAST']}_candidate_species_abundance_long.tsv"


rule candidate_species_abundance:
    input:
        diff="DGE/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx",
        bracken="abundances.cleaned.relative.tsv",
        metadata="metadata.txt"
    output:
        pdf=f"plots/{config['CONTRAST']}_candidate_species_abundance.pdf",
        selected=f"plots/{config['CONTRAST']}_selected_candidate_species.tsv",
        long=f"plots/{config['CONTRAST']}_candidate_species_abundance_long.tsv"
    params:
        top_n=config.get("top_n_candidate_taxa"),
        padj=config.get("padj_cutoff"),
        group_col="condition"
    log:
        "logs/candidate_species_abundance.log"
    shell:
        """
        Rscript ../../local/src/plot_candidate_taxa_abundance.R \
            {input.diff:q} \
            {input.bracken:q} \
            {input.metadata:q} \
            {output.pdf:q} \
            {output.selected:q} \
            {output.long:q} \
            {params.top_n} \
            {params.padj} \
            {params.group_col:q} \
            > {log:q} 2>&1
        """


rule rel_abundance_heatmap:
    input:
        relative="abundances.cleaned.clr.tsv",
        metadata="metadata.txt"
    output:
        rel_heatmap="plots/relative_heatmap.pdf"
    script:
        "../../local/src/rel_heatmap.R"


rule clr_heatmap:
    input:
        clr="abundances.filtered.clr.tsv",
        metadata="metadata.txt"
    output:
        clr_heatmap="plots/clr_heatmap.pdf"
    script:
        "../../local/src/clr_heatmap.R"


rule distance_heatmap:
    input:
        clr="abundances.filtered.clr.tsv",
        metadata="metadata.txt"
    output:
        dist_heatmap="plots/distance_heatmap.pdf"
    script:
        "../../local/src/dist_heatmap.R"


rule clr_pca:
    input:
        clr="abundances.cleaned.clr_all.tsv",
        metadata="metadata_full.txt"
    output:
        clr_pca="plots/clr_pca.pdf"
    script:
        "../../local/src/clr_pca.R"



############################################
### Batch-aware summary plotting ###
############################################

rule prepare_inputs:
    input:
        abundance="bracken_merged_abundances.num.taxid_collapsed.cleaned.txt",
        metadata="metadata.txt",
        top_table="DGE/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx"
    output:
        filtered_abundance="abundance.noKIS.tsv",
        filtered_metadata="metadata.noKIS.tsv",
        top_table_tsv="top_table.tsv"
    params:
        script="../../local/src/prepare_inputs.R"
    shell:
        """
        mkdir -p plots/
        Rscript {params.script} \
            {input.abundance} \
            {input.metadata} \
            {input.top_table} \
            {output.filtered_abundance} \
            {output.filtered_metadata} \
            {output.top_table_tsv}
        """


rule pca:
    input:
        abundance="bracken_merged_abundances.num.taxid_collapsed.cleaned.txt",
        metadata="metadata.txt"
    output:
        raw_pca="plots/pca_clr_raw.pdf",
        batch_corrected_pca="plots/pca_clr_batch_corrected.pdf"
    params:
        script="../../local/src/pca.R"
    shell:
        """
        mkdir -p plots
        Rscript {params.script} \
            {input.abundance} \
            {input.metadata} \
            {output.raw_pca} \
            {output.batch_corrected_pca}
        """


rule volcano_ma:
    input:
        top_table="DGE/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx"
    output:
        volcano="plots/volcano_resistant_vs_control.pdf",
        ma="plots/ma_resistant_vs_control.pdf",
        pvalue_hist="plots/pvalue_histogram.pdf"
    params:
        script="../../local/src/volcano_ma.R"
    shell:
        """
        mkdir -p plots
        Rscript {params.script} \
            {input.top_table} \
            {output.volcano} \
            {output.ma} \
            {output.pvalue_hist}
        """


rule top_taxa_heatmap:
    input:
        abundance="bracken_merged_abundances.num.taxid_collapsed.cleaned.txt",
        metadata="metadata.txt",
        top_table="DGE/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx"
    output:
        heatmap="plots/top30_taxa_clr_heatmap.pdf"
    params:
        script="../../local/src/top_taxa_heatmap.R",
        top_n=30
    shell:
        """
        mkdir -p plots
        Rscript {params.script} \
            {input.abundance} \
            {input.metadata} \
            {input.top_table} \
            {output.heatmap} \
            {params.top_n}
        """


rule batch_boxplots:
    input:
        abundance="bracken_merged_abundances.num.taxid_collapsed.cleaned.txt",
        metadata="metadata.txt",
        top_table="DGE/edger.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx"
    output:
        boxplots="plots/top_taxa_batch_faceted_boxplots.pdf"
    params:
        script="../../local/src/batch_boxplots.R",
        top_n=8
    shell:
        """
        mkdir -p plots
        Rscript {params.script} \
            {input.abundance} \
            {input.metadata} \
            {input.top_table} \
            {output.boxplots} \
            {params.top_n}
        """


rule batch_direction_summary:
    input:
        abundance="abundance.noKIS.tsv",
        metadata="metadata.noKIS.tsv",
        top_table="top_table.tsv"
    output:
        summary="tables/top_taxa_batch_direction_summary.tsv",
        plot="plots/batch_logFC_consistency.pdf"
    params:
        script="../../local/src/batch_direction_summary.R",
        top_n=50
    shell:
        """
        mkdir -p tables plots
        Rscript {params.script} \
            {input.abundance} \
            {input.metadata} \
            {input.top_table} \
            {output.summary} \
            {output.plot} \
            {params.top_n}
        """


rule story13_summary_plots:
    input:
        "story13/plots/pca_clr_raw.pdf",
        "story13/plots/pca_clr_batch_corrected.pdf",
        "story13/plots/volcano_resistant_vs_control.pdf",
        "story13/plots/ma_resistant_vs_control.pdf",
        "story13/plots/pvalue_histogram.pdf",
        "story13/plots/top30_taxa_clr_heatmap.pdf",
        "story13/plots/top_taxa_batch_faceted_boxplots.pdf",
        "story13/tables/top_taxa_batch_direction_summary.tsv",
        "story13/plots/batch_logFC_consistency.pdf"

