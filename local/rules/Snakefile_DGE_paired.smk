# --------------#
# Configuration #
# --------------#
import os

configfile: "../config.yaml"
if os.path.exists("Snakefile_versioned.sk"):
    include: "Snakefile_versioned.sk"


REFERENCE_ROOT=os.environ.get("REFERENCE_ROOT")

# GENERAL SETUP ---
CONDA_ROOT=config["CONDA_ROOT"]
CONDA_ACTIVATE="set +u; source %s/etc/profile.d/conda.sh ; conda activate ; conda activate" % config["CONDA_ROOT"]
RAW_DATA_DIR= ["."]

#TODO aggiungere counts_table2eset e append_each_row -> ora hanno env, ma in teoria non serve per forza

rule all:
    input:
        "edger.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.gz"

rule get_eset:
    input:
        gep = "../GEP.count.gz",
        metadata = "../metadata.txt"
    output:
        "eset.rda"
    conda:
        "../../../local/env/bit_rnaseq_3_backup.yaml"
    shell:"""
        counts_table2eset {input.gep} {input.metadata} > {output}
    """

rule get_rdata:
    input:
        "eset.rda"
    output:
        config["DGE"]["DGE_TOOL"] + ".RData"
    params:
        dge_tool = config["DGE"]["DGE_TOOL"],
        min_cpm = config["DGE"]["EXPRESSED_GENES_MIN_CPM"],
        min_samples = config["DGE"]["MIN_NUM_OF_EXPRESSED_SAMPLE"],
        factors = config["DGE"]["LIMMA_FACTORS"],
        formula = config["DGE"]["LIMMA_DESIGN_FORMULA"],
        contrasts = list(config["DGE"]["LIMMA_CONTRASTS"].values())
    conda:
        "../../../local/env/bit_rnaseq_3.yaml"
    shell:"""
        eset2toptable -t {params.dge_tool} -l {params.min_cpm} -n {params.min_samples} {params.factors} {input} {params.formula} {params.contrasts} > {output}
    """

rule test:
    output:
        "test.txt"
    params:
        contrasts = config["DGE"]["LIMMA_CONTRASTS"].values()
    shell:"""
        echo {params.contrasts} > {output}
    """


#this code allows to call the target with a value used in LIMMA_CONTRASTS_NAMES and find data in $(DGE_TOOL).RData that are stored under a label given by LIMMA_CONTRASTS
#$(addprefix $(DGE_TOOL).top.ALL.contrast., $(addsuffix .gz, $(LIMMA_CONTRASTS_NAMES))): $(DGE_TOOL).top.ALL.contrast.%.gz: $(DGE_TOOL).RData
rule run_DGE:
    input:
        config["DGE"]["DGE_TOOL"] + ".RData"
    output:
        config["DGE"]["DGE_TOOL"] + ".toptable_clean.contrast_{contrast}.gz"
    params:
        contrast_names = list(config["DGE"]["LIMMA_CONTRASTS"].keys()),
        contrasts = list(config["DGE"]["LIMMA_CONTRASTS"].values()),
        dge_tool = config["DGE"]["DGE_TOOL"] 
    conda:    
        "../../../local/env/bit_rnaseq_3.yaml"
    shell:"""
        C=$(\
            (\
                echo {params.contrast_names};\
                echo {params.contrasts}\
            ) | tr -d '"'\
            | perl -lane 'BEGIN{{$i=0}}; if($.==1){{for(@F){{if($_ ne "{wildcards.contrast}"){{$i++}}else{{last}}}}}} if($.==2){{print @F[$i]}}'\
        );\
        echo $C;\
        if [[ {params.dge_tool} == "limma" ]]; then\
            r -e "load('{input}');write.table(file='{output}.tmp',top.list[['$C']][,c('logFC','P.Value','adj.P.Val')], sep='\t', row.names=TRUE, quote=FALSE)"; \
        elif [[ {params.dge_tool} == "edger" ]]; then \
            r -e "load('{input}');write.table(file='{output}.tmp',top.list[['$C']][,c('logFC','PValue','FDR')], sep='\t', row.names=TRUE, quote=FALSE)"; \
        elif [[ {params.dge_tool} == "deseq2" ]]; then \
            r -e "load('{input}');write.table(file='{output}.tmp',top.list[['$C']][,c('log2FoldChange','pvalue','padj')], sep='\t', row.names=TRUE, quote=FALSE)"; \
        fi;
        cat {output}.tmp | unhead | gzip > {output};
        rm {output}.tmp
    """

rule all_contrasts:
    input:
        #expand("{path}.toptable_clean.contrast_{contrast}.gz", contrast=config["DGE"]["LIMMA_CONTRASTS_NAMES"])
        lambda wildcards: expand("{path}.toptable_clean.contrast_{contrast}.gz", 
            path=wildcards.path, 
            contrast=list(config["DGE"]["LIMMA_CONTRASTS"].keys()))
    output:
        "{path}.toptable_clean.ALL_contrast.gz"
    params:
        contrast=list(config["DGE"]["LIMMA_CONTRASTS"].keys())
    conda:
        "../../../local/env/bit_rnaseq_3.yaml"
    shell:"""
        for i in {params.contrast}; do
            zcat {wildcards.path}.toptable_clean.contrast_$i.gz | append_each_row -B $i;
        done | gzip > {output}
        """


# .META: *.toptable_clean*.gz
# 	1	GeneID
# 	2	logFC
# 	3	Pvalue
# 	4	Pvalue_adj
rule mark_seqc:
    input:
        "{path}.toptable_clean.ALL_contrast.gz"
    output:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.gz"
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
        deg_out = "{path}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz", 
        gep = "../GEP.count.exp_filter.ltmm.lfpkm.metadata.max_exp_in_condition.gz"
    output:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.max_exp_in_condition.header_added.gz"
    shell:"""
        bawk '$Pvalue_adj!="NA"' {input.deg_out} | translate -a -v <(zcat {input.gep}) 2 | gzip > {output}
    """

rule exp_in_cond:
    input:
        deg_out = "{path}.toptable_clean.ALL_contrast.mark_seqc.header_added.gz", 
        gep = "../GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.gz"
    output:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.exp_in_condition.header_added.gz"
    shell:"""
        bawk '$Pvalue_adj!="NA"' {input.deg_out}| translate -a <(zcat {input.gep}) 2 | gzip > {output}
    """

rule filter_genes:
    input: 
        deg_out = "{path}.toptable_clean.ALL_contrast.mark_seqc.{expression}.header_added.gz", 
        gep = "../GEP.{analysis}.gz"
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

rule count_significance:
    input:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.gz"
    output:
        "{path}.toptable_clean.ALL_contrast.mark_seqc.count"
    shell:
        """
        bawk '{{print $1, $6}}' {input} | symbol_count | tab2matrix -r contrast > {output}
        """

# -------------- #
# Generic rules  #
# -------------- #

rule tab2xlsx:
    input: 
        "{file}"
    output: 
        "{file}.xlsx"
    shell: 
        "cat < {input} | tab2xlsx > {output}"

rule gz2xlsx:
    input: 
        "{file}.gz"
    output: 
        "{file}.xlsx"
    shell: 
        "zcat < {input} | tab2xlsx > {output}"

rule add_header:
    input: 
        "{path}.gz"
    output: 
        "{path}.header_added.gz"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; zcat {input} ) | gzip > {output}"

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
