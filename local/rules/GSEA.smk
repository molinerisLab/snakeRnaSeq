import os

configfile:"config.yaml"

if os.path.exists("Snakefile_versioned.sk"):
    include: "Snakefile_versioned.sk"


# WARNING
# !!
# RNAseqRtools and bioinfotree tools must be installed in the base env
# !!

################
#
# Global Params
#

REFERENCE_ROOT=os.environ.get("REFERENCE_ROOT")
BIOINFO_REFERENCE_ROOT=REFERENCE_ROOT + "/bioinfotree/task"

#CONDA_ACTIVATE="set +u; source %sminiconda2/etc/profile.d/conda.sh ; conda activate ; conda activate" % REFERENCE_ROOT

MSIGDB_VERSION=config["MSIGDB_VERSION"]
MSIGDB_DIR=BIOINFO_REFERENCE_ROOT+"/msigdb/dataset/"+MSIGDB_VERSION

GMT=config["MSIGDB_TYPES"][config["MSIGDB_VERSION"]]

GSEA_DIR=os.environ.get("PWD")
GSEA_PUBLISH_USER=config["GSEA_PUBLISH_USER"]
GSEA_PUBLISH_HOST=config["GSEA_PUBLISH_HOST"]

RNK_FILES, = glob_wildcards("{rnk_file}.rnk")

if not RNK_FILES:
	GLOB = glob_wildcards("../../../{dgetool}.toptable_clean.contrast_{contrast}.gz")
	RNK_FILES = ["{d}.{c}.{r}".format(d=dge, c=cont, r=config["RNK_METRIC"]) for dge, cont in zip(GLOB.dgetool, GLOB.contrast)]
if not RNK_FILES:
	exit("Files .rnk not found and can not find suitable files in BiT RNAseq pipeline, see https://github.com/molinerisLab/SnakeGSEA/")

##############
#
# ALL targets
#

rule all:
	input:
		"multi_GSEA.link.header_added.gz",
		"multi_GSEA.link.header_added.xlsx",
		"multi_GSEA.html"
#		"gsea_publish_done"


################
#
# Generic rules
#

rule header_added_gz:
	input:
		"{file}.gz"
	output:
		"{file}.header_added.gz"
	shell:"""
		(bawk -M {input} | cut -f 2 | tr "\\n" "\\t" | perl -lne 's/\s+$//; print'; zcat {input} ) | gzip > {output}
	"""

rule tab_to_xlsx:
	input:
		"{file}.gz"
	output:
		"{file}.xlsx"
	shell: """
		zcat {input} | tab2xlsx > {output}
	"""

######################
#
# BiT RNAseq  specific rules
#

rule logFC_rnk:
	input:
		"../../../{dgetool}.toptable_clean.contrast_{contrast}.gz"
	output:
		"{dgetool}.{contrast}.logFC.rnk"
	shell: """
		bawk '$Pvalue_adj!="NA" {{print $GeneID, $logFC}}' {input} | sort -k2,2g > {output}
	"""

rule logFC_pval_rnk:
	input:
		"../../../{dgetool}.toptable_clean.contrast_{contrast}.gz"
	output:
		"{dgetool}.{contrast}.logFC_pval.rnk"
	shell: """
		bawk '$Pvalue_adj!="NA" {{print $GeneID, $logFC*(-log($Pvalue)/log(10))}}' {input} | sort -k2,2g > {output}
	"""

"""
.META: *.rnk
	1	GeneID
	2	value

"""

######################
#
# GESA specific rules
#

rule create_gmt_gz:
	input:
		MSIGDB_DIR+"/{gene_set}."+MSIGDB_VERSION+".symbols.gmt.gz"
	output:
		"{gene_set}.gmt"
	shell:"""
		gunzip < {input} > {output}
	"""

#rule create_gmt:
#	input:
#		MSIGDB_DIR+"/{gene_set}."+MSIGDB_VERSION+".symbols.gmt"
#	output:
#		"{gene_set}.gmt"
#	wildcard_constraints:
#	        gene_set="[^\/]+"
#	shell:"""
#		ln -sf {input} {output}
#	"""

rule run_fgsea:
	input:
		gmt_file="{gmt_file}.gmt",
		rnk_file="{rnk_file}.rnk"
	output:
		"fGSEA_Results/{rnk_file}.rnk.{gmt_file}.gmt.fGSEA.tab"
	threads: 24
	#conda:
	#	"../../local/env/fGSEA_Env.yaml"
	shell: """
		Rscript ../../local/bin/fgsea.CommandLine.R \
			-g {input.gmt_file} \
			-r {input.rnk_file} \
			-n 1000 \
			-b {config[MAX_GENES_IN_SET]} \
			-s {config[MIN_GENES_IN_SET]} \
			-d fGSEA_Results \
			-l Leading_Edge \
			--gseaParam {config[GSEA_WEIGHT]}\
			-j {threads}
	"""

rule multi_fGSEA:
	input:
		expand("fGSEA_Results/{rnk_file}.rnk.{gmt}.gmt.fGSEA.tab"
			, rnk_file=RNK_FILES, gmt=GMT)
	output:
		"multi_GSEA.gz"
	shell: """
		matrix_reduce -t 'fGSEA_Results/*.rnk.*.gmt.fGSEA.tab' \
		| tr ';' '\\t' | grep -v 'nMoreExtreme' | sort -k4,4g | gzip > {output}
	"""

"""
.META: multi_GSEA.gz
	1	contrast
	2	msigdb_type
	3	pathway
	4	pval
	5	padj
	6	ES
	7	NES
	8	nMoreExtreme
	9	size
	10	leadingEdge

"""

rule mdsum_file:
	input:
		expand("{rnk_file}.rnk"
			, rnk_file=RNK_FILES)
	output:
		"rnk_md5sum.txt"
	shell: """
		for i in $(ls *rnk); do \
			touch -a {output}; \
			a=$PWD'/'$i; \
			b=$(md5sum $a | tr " " "\\t" | cut -f 1); \
			echo | bawk -v b=$b -v i=$i '{{print "ln -sf "i" ash_rnk-"b}}' >> {output}; \
		done; \
		cat {output} | sh
	"""

rule multi_fGSEA_link:
	input:
		gsea="multi_GSEA.gz",
		ash_rnk="rnk_md5sum.txt"
	output:
		"multi_GSEA.link.gz"
	shell: """
		zcat multi_GSEA.gz \
		| translate -a -r <(tr " " "\\t" < {input.ash_rnk} | cut -f3- | sed "s/.rnk//") 1  \
		| bawk '{{HTTP="http://"; if(length($3)>100){{HTTP=""}} print $1~10, \
			"molinerislab.shinyapps.io/shinysea/?rnk=" $11 "&geneset=" $2 "_" $3 "&gseaParam=" {config[GSEA_WEIGHT]}, \
			"www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=" $3}}' \
		| gzip > {output}
	"""

"""
.META: multi_GSEA.link.gz
	1	contrast
	2	msigdb_type
	3	pathway
	4	pval
	5	padj
	6	ES
	7	NES
	8	nMoreExtreme
	9	size
	10	leadingEdge
	11	shinySea_link
	12	pathway_link

"""

rule fGSEA_Rmd_template:
	input:
		"../../local/src/fGSEA_Rmd_template.Rmd"
	output:
		"fGSEA_Rmd_template.Rmd"
	shell: """
		cp {input} {output}
	"""

rule multi_fGSEA_link_html:
	input:
		multi_GSEA="multi_GSEA.link.header_added.gz",
		rmd="fGSEA_Rmd_template.Rmd"
	output:
		"multi_GSEA.html"
	shell:"""
		Rscript -e 'library(rmarkdown); rmarkdown::render("{input.rmd}", output_file="{output}", quiet=TRUE)'
	"""

rule publish:
	input:
		xlsx="multi_GSEA.link.header_added.xlsx",
		html="multi_GSEA.html",
		ash_rnk="rnk_md5sum.txt"
	output:
		"gsea_publish_done"
	shell:"""
		GSEA_PUBLISH_DIR=$(pwd | sed 's/{config[GSEA_PUBLISH_REPLACE_FROM]}/{config[GSEA_PUBLISH_REPLACE_TO]}/');
		ssh {GSEA_PUBLISH_USER}@{GSEA_PUBLISH_HOST} mkdir -p $GSEA_PUBLISH_DIR; 
		rsync -avP {input.xlsx} {input.html} *.rnk ash_rnk-* {GSEA_PUBLISH_USER}@{GSEA_PUBLISH_HOST}:$GSEA_PUBLISH_DIR; 
		ssh {GSEA_PUBLISH_USER}@{GSEA_PUBLISH_HOST} ln -sf $GSEA_PUBLISH_DIR/ash_rnk-* /var/www/html/epigen/shinySea/hash; 
		echo {GSEA_PUBLISH_HOST}/$(sed 's|/var/www/html/||' <<<$GSEA_PUBLISH_DIR)/input.html | tee {output}
	"""

