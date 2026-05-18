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


rule explore_abundance:
    input:
        relative="abundances.filtered.relative.tsv",
        clr="abundances.filtered.clr.tsv",
        metadata="metadata_update.txt"
    output:
        outdir=directory("plots")
    script:
        "../../local/src/explore_abundance.R"


rule rel_abundance_barplot:
    input:
        relative="abundances.filtered.relative.tsv",
        metadata="metadata.txt"
    output:
        barplot="plots/stacked_barplot.pdf"
    script:
        "../../local/src/rel_barplot.R"


rule rel_abundance_heatmap:
    input:
        relative="abundances.filtered.relative.tsv",
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
        clr="abundances.filtered.clr.tsv",
        metadata="metadata.txt"
    output:
        clr_pca="plots/clr_pca.pdf"
    script:
        "../../local/src/clr_pca.R"

