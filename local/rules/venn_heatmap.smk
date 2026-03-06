# -------------------------
#  Venn comparison - genes

# tool_comparison.%.venn_in.gz: limma.toptable_clean_cut.%.gz edger.toptable_clean_cut.%.gz deseq2.toptable_clean_cut.%.gz
# 	(for i in limma edger deseq2; do bawk -v t=$$i '{print $$GeneID, t, "1"}' $$i.toptable_clean_cut.$*.gz; done) | tab2matrix -e 0 | gzip > $@
rule venn_prep_up:
    input:
        lambda wildcards: f"DGE/{config['DGE']['DGE_TOOL']}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
    output:
        temp("DGE/gene_comparison.{contrast}_up.venn_in.gz")
    shell:"""
        bawk -v contrast={wildcards.contrast} '{{if ($contrast == contrast) print $GeneID, $contrast, $significance}}' {input} \
        | bawk '{{up=($3 == -1 ? 0 : $3); print $0, up}}' | bawk '{{print $1, $2 "_up", $4}}' \
        | gzip > {output}
        """

rule venn_prep_down:
    input:
        lambda wildcards: f"DGE/{config['DGE']['DGE_TOOL']}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz"
    output:
        temp("DGE/gene_comparison.{contrast}_down.venn_in.gz")
    shell:"""
        bawk -v contrast={wildcards.contrast} '{{if ($contrast == contrast) print $GeneID, $contrast, $significance}}' {input} \
        | bawk '{{down=($3 == 1 ? 0 : $3); print $0, down}}' | bawk '{{print $1, $2 "_down", $4}}' \
        | gzip > {output}
        """


rule venn_input:
    input:
        expand("DGE/gene_comparison.{contrast}_{dge}.venn_in.gz", contrast = config["VENN_CONTRASTS"], dge = ["up", "down"])
    output:
        "DGE/gene_comparison.venn_in.gz"
    shell:"""
        zcat {input} | tab2matrix | gzip > {output}
        """

# rule venn_pdf:
#     input:
#         "DGE/gene_comparison.venn_in.gz"
#     output:
#         "DGE/gene_comparison.all_venn.pdf",
#     conda:
#         "../../local/env/bit_rnaseq_3_backup.yaml"
#     shell:"""
#         zcat {input} | venn -a > {output}
#         """

rule venn_genes_prep:
    input:
        "DGE/edger.toptable_clean.ALL_contrast.mark_seqc.gz"
    output:
        "DGE/edger.toptable_clean.ALL_contrast.mark_seqc.matrix.gz"
    shell:"""
        bawk '{{print $GeneID, $contrast, $significance}}' {input} | tab2matrix | gzip > {output}
    """


# ####################### HEATMAP CONTRASTI ######################
# ##################### heatmap_DGE_deseq2/ #####################

rule heatmap:
    input:
        gep = "GEP.count.exp_filter.ltmm.gz",
        metadata = "metadata_heatmap.txt"
    output:
        pdf = "heatmap.top100.pdf",
        rds = "heatmap-top100.Rds"
    params:
        sets = config["HEATMAP"]["COLUMNS"],
        n_top = config["HEATMAP"]["N_TOP"],
        color_set = config["HEATMAP"]["COLORS"],
        clusters = config["HEATMAP"]["CLUSTERS"]
    shell:"""
        heatmap {input.gep} -a {input.metadata} -C {params.sets} -c 12 -s -d correlation {params.color_set} --pdf -e 16 -w 12 {params.clusters} -N {params.n_top} -n -R {output.rds} {output.pdf}
    """



#/home/reference_data/bioinfotree/local/bin//heatmap

# heatmap_DGE_deseq2/deseq2_dds_rlog.txt: deseq2.RData
# 	r -e 'suppressMessages(library("DESeq2")); load("$<"); rld = rlog(dds, blind=FALSE); y = rld@assays$$data@listData; y = matrix(unlist(y),ncol=length(dimnames(dds@assays$$data@listData$$counts)[[2]])); colnames(y)=dimnames(dds@assays$$data@listData$$counts)[[2]]; rownames(y)=dimnames(dds@assays$$data@listData$$counts)[[1]]; write.table(y, "$@", sep="\t")'

# heatmap_DGE_deseq2/ctrl_rlog.txt: deseq2.toptable_clean_cut.contrast.CTRL_vs_onlyF.cut_seqc_both.header_added.xls heatmap_DGE_deseq2/deseq2_dds_rlog.txt
# 	r -e 'y = read.delim("$^2", header=TRUE); require(gdata); list = read.xls("$<"); list = as.character(list[,1]); index = sort(match(list, rownames(y))); ctrl = y[index,]; write.table(ctrl, "$@", sep="\t")'

# heatmap_DGE_deseq2/adenoca_rlog.txt: deseq2.toptable_clean_cut.contrast.adenoCA_vs_onlyF.cut_seqc_both.header_added.xls heatmap_DGE_deseq2/deseq2_dds_rlog.txt
# 	r -e 'y = read.delim("$^2", header=TRUE); require(gdata); list = read.xls("$<"); list = as.character(list[,1]); index = sort(match(list, rownames(y))); adenoca = y[index,]; write.table(adenoca, "$@", sep="\t")'

# heatmap_DGE_deseq2/igg4_rlog.txt: deseq2.toptable_clean_cut.contrast.IgG4_vs_onlyF.cut_seqc_both.header_added.xls heatmap_DGE_deseq2/deseq2_dds_rlog.txt
# 	r -e 'y = read.delim("$^2", header=TRUE); require(gdata); list = read.xls("$<"); list = as.character(list[,1]); index = sort(match(list, rownames(y))); igg4 = y[index,]; write.table(igg4, "$@", sep="\t")'

# heatmap_DGE_deseq2/metadata_to_heatmap.txt: metadata_clean.txt 
# 	cat $< | cut -f2- > $@

# heatmap_DGE_deseq2/%_rlog.png: heatmap_DGE_deseq2/%_rlog.txt heatmap_DGE_deseq2/metadata_to_heatmap.txt ../../local/src/heatmap.R
# 	$^3 -C coculture_method,condition_vs_onlyF,FBL_group,library $< $^2 $@

# heatmap_DGE_deseq2/%_rlog_dend.png: heatmap_DGE_deseq2/%_rlog.txt heatmap_DGE_deseq2/metadata_to_heatmap.txt ../../local/src/heatmap.R
# 	$^3 -C coculture_method,condition_vs_onlyF,FBL_group,library --cluster 12 $< $^2 $@

# heatmap_DGE_deseq2/%_rlog_dend_row.png: heatmap_DGE_deseq2/%_rlog.txt heatmap_DGE_deseq2/metadata_to_heatmap.txt ../../local/src/heatmap.R
# 	$^3 -C coculture_method,condition_vs_onlyF,FBL_group,library --cluster 1 $< $^2 $@







 