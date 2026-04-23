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

# ####################################
# ### Rules for barplot generation ###
# ####################################
# rule barplot:
#    input:
#       bkrona="b_krona_txt/{sample}.b.krona.txt"
#    output:
#       bkrona_phylum="{sample}.barplot.txt"
#    shell:
#      """
#       mkdir -p barplot/	
#       bawk '$1!=0 {{print $4,$1}}' {input.bkrona}  | sort |stat_base -g -t > {output.bkrona_phylum}
#       mv {output.bkrona_phylum} barplot/{output.bkrona_phylum}
#      """

# rule common_taxa_in_samples:
#     params:
#         samples=config["samples"]
#     script:
#         """
#         Rscript src/common_species_samples.R
#         """


# rule reads_human_contam_classified:
#     output:
#             "Plots/classified_vs_human_contaminant_barplot_normalized.png", 
#             "Plots/classified_vs_human_contaminant_barplot_percentage.png"       
#     params:
#         split_by=config["split_by"]
#     shell: 
#         """
#         Rscript src/plot_reads.R
#         """

# rule top_taxa_per_sample:
#     shell:
#         """
#         Rscript src/plot_sample_species.R
#         """


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

rule explore_abundance:
    input:
        relative="abundances.filtered.relative.tsv",
        clr="abundances.filtered.clr.tsv",
        metadata="metadata.txt"
    output:
        outdir=directory("plots/exploration")
    script:
        "../../local/src/explore_abundance.R"

rule rel_abundance:
    input:
        relative="abundances.filtered.relative.tsv",
        metadata="metadata.txt"
    output:
        outdir=directory("plots")
    script:
        "../../local/src/abundance_barplot.R"

