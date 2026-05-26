# DGE/edger.toptable_clean.ALL_contrast.mark_seqc.max_exp_in_condition.header_added.count.exp_filter.ltmm.gz
# edger.toptable_clean.ALL_contrast.mark_seqc.max_exp_in_condition.header_added.gz

rule ALL_DGE:
    input:
        expand(
            "{folder}/{tool}.toptable_clean.ALL_contrast.mark_seqc.header_added.xlsx",
            folder="DGE",
            tool=config["DGE"]["DGE_TOOL"]
        )

rule ALL_DGE_exp_added:
    input:
        expand(
            "{folder}/{tool}.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.xlsx",
            folder="DGE",
            tool=config["DGE"]["DGE_TOOL"]
        )


#TODO aggiungere counts_table2eset e append_each_row -> ora hanno env, ma in teoria non serve per forza
#Renamed bracken_merged_abundances.num.taxid_collapsed.cleaned.txt to GEP.count and gzipped
rule get_eset:
    input:
        gep = "GEP.count.gz",
        metadata = "metadata.txt"
    output:
        "{folder}/eset.rda"
    conda:
        "../../local/env/bit_rnaseq_3_backup.yaml"
    shell:"""
        counts_table2eset {input.gep} {input.metadata} > {output}
    """

rule get_rdata:
    input:
        "{folder}/eset.rda"
    output:
        "{folder}/" + config["DGE"]["DGE_TOOL"] + ".RData"
    params:
        dge_tool = config["DGE"]["DGE_TOOL"],
        min_cpm = config["DGE"]["EXPRESSED_GENES_MIN_CPM"],
        min_samples = config["DGE"]["MIN_NUM_OF_EXPRESSED_SAMPLE"],
        factors = config["DGE"]["LIMMA_FACTORS"],
        formula = "'" + config["DGE"]["LIMMA_DESIGN_FORMULA"] + "'",
        contrasts = " ".join(config["DGE"]["LIMMA_CONTRASTS"].values())
    conda:
        "../../local/env/bit_rnaseq_3_backup.yaml"
    shell:"""
        eset2toptable -t {params.dge_tool} -l {params.min_cpm} -n {params.min_samples} {params.factors} {input} {params.formula} {params.contrasts} > {output}
    """

#this code allows to call the target with a value used in LIMMA_CONTRASTS_NAMES and find data in $(DGE_TOOL).RData that are stored under a label given by LIMMA_CONTRASTS
#$(addprefix $(DGE_TOOL).top.ALL.contrast., $(addsuffix .gz, $(LIMMA_CONTRASTS_NAMES))): $(DGE_TOOL).top.ALL.contrast.%.gz: $(DGE_TOOL).RData

rule run_DGE:
    input:
        "{folder}/" + config["DGE"]["DGE_TOOL"] + ".RData"
    output:
        "{folder}/" + config["DGE"]["DGE_TOOL"] + ".toptable_clean.contrast_{contrast}.gz"
    params:
        contrast_value = lambda wildcards: config["DGE"]["LIMMA_CONTRASTS"][wildcards.contrast],
        dge_tool = config["DGE"]["DGE_TOOL"]
    conda:
        "../../local/env/bit_rnaseq_3_backup.yaml"
    shell:
        r"""
        C="{params.contrast_value}"
        echo "Using contrast: [$C]"

        if [[ "{params.dge_tool}" == "limma" ]]; then
            r -e "load('{input}'); write.table(file='{output}.tmp', top.list[['$C']][,c('logFC','P.Value','adj.P.Val')], sep='\t', row.names=TRUE, quote=FALSE)"
        elif [[ "{params.dge_tool}" == "edger" ]]; then
            r -e "load('{input}'); write.table(file='{output}.tmp', top.list[['$C']][,c('logFC','PValue','FDR')], sep='\t', row.names=TRUE, quote=FALSE)"
        elif [[ "{params.dge_tool}" == "deseq2" ]]; then
            r -e "load('{input}'); write.table(file='{output}.tmp', top.list[['$C']][,c('log2FoldChange','pvalue','padj')], sep='\t', row.names=TRUE, quote=FALSE)"
        else
            echo "ERROR: unknown DGE tool: {params.dge_tool}" >&2
            exit 1
        fi

        echo "Checking temporary output:"
        ls -lh {output}.tmp
        wc -l {output}.tmp
        sed -n '1,10p' {output}.tmp

        # Keep a human-readable debug version WITH header
        cp {output}.tmp {output}.debug.tsv

        # Remove header for downstream pipeline compatibility
        cat {output}.tmp | unhead | gzip > {output}

        # Remove temporary file
        rm {output}.tmp

        echo "Checking gzipped output:"
        gzip -cd {output} | wc -l
        gzip -cd {output} | sed -n '1,10p'
        """

rule all_contrasts:
    input:
        #expand("{path}.toptable_clean.contrast_{contrast}.gz", contrast=config["DGE"]["LIMMA_CONTRASTS_NAMES"])
        lambda wildcards: expand("{folder}/{path}.toptable_clean.contrast_{contrast}.gz", 
            folder=wildcards.folder,
            path=wildcards.path, 
            contrast=config["DGE"]["LIMMA_CONTRASTS_NAMES"])
    output:
        "{folder}/{path}.toptable_clean.ALL_contrast.gz"
    params:
        contrast=config["DGE"]["LIMMA_CONTRASTS_NAMES"]
    conda:
        "../../local/env/bit_rnaseq_3_backup.yaml"
    shell:"""
        for i in {params.contrast}; do
            zcat {wildcards.folder}/{wildcards.path}.toptable_clean.contrast_$i.gz | append_each_row -B $i;
        done | gzip > {output}
        """

rule mark_seqc:
    input:
        "{folder}/{path}.toptable_clean.ALL_contrast.gz"
    output:
        "{folder}/{path}.toptable_clean.ALL_contrast.mark_seqc.gz"
    params:
        logfc = config["MARK_SEQC"]["LOGFC"],
        pvalue = config["MARK_SEQC"]["PVALUE"],
        p_adj = config["MARK_SEQC"]["PVALUE_ADJ"]
    shell:"""
        bawk '{{M=0; \
            if($Pvalue_adj!="NA" && sqrt($logFC*$logFC)>{params.logfc} && $Pvalue+0<{params.pvalue} && $Pvalue_adj+0<{params.p_adj}) \
            {{ \
                if($logFC>0) {{M=1}}else{{M=-1}}\
            }} \
            print $0,M \
        }}' {input} | gzip > {output} 
    """
    #abs(a)>1 <==> a*a > 1

rule max_exp_in_cond:
    input:
        deg_out = "{folder}/{path}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz", 
        gep = "{folder}/GEP.count.exp_filter.ltmm.metadata.max_exp_in_condition.gz"
    output:
        "{folder}/{path}.toptable_clean.ALL_contrast.mark_seqc.max_exp_in_condition.header_added.gz"
    shell:"""
        bawk '$Pvalue_adj!="NA"' {input.deg_out} | translate -a -v <(zcat {input.gep}) 2 | gzip > {output}
    """

rule exp_in_cond:
    input:
        deg_out = "{folder}/{path}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz", 
        # gep = "{folder}/GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.gz"
        gep = "GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.gz"
    output:
        "{folder}/{path}.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.gz"
    shell:"""
        bawk '$Pvalue_adj!="NA"' {input.deg_out}| translate -a -v -e "NA" <(zcat {input.gep}) 2 | gzip > {output}
    """

rule filter_genes:
    input: 
        deg_out = "{path}.toptable_clean.ALL_contrast.mark_seqc.{expression}.header_added.gz", 
        gep = "GEP.{analysis}.gz"
    output: 
        "{path}.toptable_clean.ALL_contrast.mark_seqc.{expression}.header_added.{analysis}.gz"
    shell:"""
        zcat {input.deg_out} | translate -a -r <(zcat {input.gep} | sed s/Geneid/GeneID/) 2 | gzip > {output}
    """

rule filter_significant_genes:
    input:
        deg_out = "{path}.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.gz",
        significance = "{path}.toptable_clean.ALL_contrast.mark_seqc.gz"
    output:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.signif.gz"
    shell:"""
        zcat {input.deg_out} | filter_1col --header 1 2 <(bawk '$significance!=0 {{print $2}}' {input.significance}) | gzip > {output}
    """

rule DEG_count_matrix:
    input:
        "DGE/{path}.toptable_clean.ALL_contrast.mark_seqc.gz"
    output:
        "DGE/{path}.toptable_clean.ALL_contrast.mark_seqc.DEG_count_matrix"
    shell:"""
        bawk '{{print $contrast,$significance}}' {input} | symbol_count | tab2matrix -r contrast > {output}
    """
 
rule count_significance:
    input:
        "DGE/edger.toptable_clean.ALL_contrast.mark_seqc.gz"
    output:
        "DGE/edger.toptable_clean.ALL_contrast.mark_seqc.count"
    shell:
        """
        bawk '{{print $1, $6}}' {input} | symbol_count | tab2matrix -r contrast > {output}
        """

# -------------- #
# META           #
# -------------- #

#TODO move in correct .mk after https://www.pivotaltracker.com/story/show/188249716
"""
.META: *.toptable_clean.ALL_contrast.*mark_seqc.*gz 
	1	contrast
	2	GeneID
	3	logFC
	4	Pvalue
	5	Pvalue_adj
	6	significance

.META: *.toptable_clean.ALL_contrast.*gz 
	1	contrast
	2	GeneID
	3	logFC
	4	Pvalue
	5	Pvalue_adj

.META: *.toptable_clean*.gz
	1	GeneID
	2	logFC
	3	Pvalue
	4	Pvalue_adj

"""
