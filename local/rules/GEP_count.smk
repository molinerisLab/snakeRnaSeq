# GEP.count.gz: $(FASTQ_FILTERING).featurecounts.ribo.ex.count.gz
# 	zgrep -v '^#' $< | cut -f 1,7- \
# 	| perl -pe 'if($$.==1){s|STAR/fastq/||g; s|.STAR[^\s]+.bam||g; s|_S\d+(\s)|\1|g}' \
# 	| gzip > $@


rule get_Gep:
    input: 
        featurecounts="featurecounts.ribo.ex.count.gz"
    output: 
        "GEP.count.gz"
    shell: """
        zgrep -v '^#' {input.featurecounts} | cut -f 1,7- \
        | perl -pe 'if($.==1){{s|star/||g;}}' \
        | gzip > {output}  
    """  


# GEP.count%exp_filter.gz: GEP.count%gz GEP.count.cpm.expressed_genes
# 	zcat $< | filter_1col --header 1 1 <(cut -f 1 $^2) | gzip > $@

rule Gep_filter:
    input:
        GEP_count="GEP.count{filter}gz",
        GEP_cpm="GEP.count.cpm.expressed_genes"
    output: 
        "GEP.count{filter}exp_filter.gz"
    shell: 
        "zcat {input.GEP_count} | filter_1col --header 1 1 <(cut -f 1 {input.GEP_cpm}) | gzip > {output}"

# GEP.count%.metadata.max_exp_in_condition.gz: GEP.count%.metadata.gz
# 	bawk 'NR>1 {print $$GeneID";"$$$(DGE_CONDITON),$$exp}' $< | bsort -k1,1 | stat_base -g -a | tr ";" "\t" | find_best 1 3 | gzip > $@

rule GEP_metadata:
    input: 
        "GEP.count{filter}.metadata.gz"
    output: 
        "GEP.count{filter}.metadata.max_exp_in_condition.gz"
    params: 
        DGE=config["DGE"]["DGE_CONDITION"]
    shell:""" 
        bawk 'NR>1 {{print $GeneID";"${params.DGE},$exp}}' {input} | bsort -k1,1 | stat_base -g -a | tr ";" "\t" | find_best 1 3 | gzip > {output}
    """

#tolto e non tradotto
# .META: GEP.count*metadata.max_exp_in_condition.gz
# 	1	GeneID
# 	2	best_condition
# 	3	exp

# GEP.count.cpm.expressed_genes: GEP.count.cpm.gz

# 	zcat $< | matrix2tab | bawk '$$3>$(EXPRESSED_GENES_MIN_CPM)  {print $$1}' | symbol_count |bawk '$$2>=$(MIN_NUM_OF_EXPRESSED_SAMPLE)'> $@

rule expressed_genes:
    input: 
        "GEP.count.cpm.gz"
    output: 
        "GEP.count.cpm.expressed_genes"
    params: 
        expressed=config["DGE"]["EXPRESSED_GENES_MIN_CPM"],
        min_sample=config["DGE"]["MIN_NUM_OF_EXPRESSED_SAMPLE"]
    shell: """
        zcat {input} | matrix2tab | bawk '$3>{params.expressed}  {{print $1}}' | symbol_count |bawk '$2>={params.min_sample}'> {output}
    """

# gene_len: $(FASTQ_FILTERING).featurecounts.ribo.ex.count.gz
# 	echo -e "Geneid\tlength" > $@
# 	zgrep -v '^#' $< | cut -f 1,6 | unhead -n 1 >> $@

rule gene_length:
    input: 
        featurecounts="featurecounts.ribo.ex.count.gz"
    output: 
        "gene_len"
    shell: """
        echo -e "Geneid\tlength" > {output}
        zgrep -v '^#' {input.featurecounts} | cut -f 1,6 | unhead -n 1 >> {output}
    """



# GEP.count.rpk.gz: GEP.count.gz gene_len
# 	zcat $< | translate -a $^2 1 \
# 	| perl -wlane 'BEGIN{$$,="\t"} {$$g=shift(@F); $$l=shift(@F); if($$.==1){print $$g,@F}else{@F=map {$$_=$$_/$$l} @F; print $$g,@F}}' \
# 	| gzip > $@

rule rpk:
    input: 
        GEP="GEP.count.gz",
        GenLen="gene_len"
    output: 
        "GEP.count.rpk.gz"
    shell: """
        zcat {input.GEP} | translate -a {input.GenLen} 1 \
        | perl -wlane 'BEGIN{$,="\t"} {$g=shift(@F); $l=shift(@F); if($.==1){print $g,@F}else{@F=map {$_=$_/$l} @F; print $g,@F}}' \
        | gzip > {output}
    """


# GEP.count%cpm.gz: GEP.count%gz
# 	r -e 'library(edgeR); x<-read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, normalized.lib.sizes=FALSE), "$@.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule cpm:
    input: 
        "GEP.count{filter}gz"
    output: 
        "GEP.count{filter}cpm.gz"
    shell: """
        zcat {input} | grep -v "Geneid"  > {input}.fixed;
        r -e 'library(edgeR); x<-read.table("{input}.fixed", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, normalized.lib.sizes=FALSE), "{output}.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {input}.fixed {output}.tmp
    """

# GEP.count.rpk.cpm.gz: GEP.count.rpk.gz
# 	echo -n "Geneid" | gzip >$@
# 	r -e 'library(edgeR); x<-read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, normalized.lib.sizes=FALSE), "$@.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule rpk_cpm:
    input: 
        "GEP.count.rpk.gz"
    output: 
        "GEP.count.rpk.cpm.gz"
    shell: """
        echo -n "Geneid" | gzip >{output}
        r -e 'library(edgeR); x<-read.table("{input}", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, normalized.lib.sizes=FALSE), "{output}.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """
# GEP.count%lcpm.gz: GEP.count%gz
# 	r -e 'library(edgeR); x<-read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, log=T, prior.count=$(LCPM_PRIOR_COUNT)), "$@.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule lcpm:
    input: 
        "GEP.count{filter}gz"
    output: 
        "GEP.count{filter}lcpm.gz"
    params: 
        LCPM=LCPM_PRIOR_COUNT
    shell: """
        r -e 'library(edgeR); x<-read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, log=T, prior.count=${params.LCPM}), "{output}.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """


# GEP.count.rpk.ltpm.gz: GEP.count.rpk.gz
# 	echo -n "Geneid" | gzip >$@
# 	r -e 'library(edgeR); x<-read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, log=T, prior.count=$(LTPM_PRIOR_COUNT)), "$@.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule rpk_ltpm:
    input:
        "GEP.count.rpk.gz"
    output:
        "GEP.count.rpk.ltpm.gz"
    params: 
        LTPM=LTPM_PRIOR_COUNT
    shell: """
        echo -n "Geneid" | gzip >{output}
        r -e 'library(edgeR); x<-read.table("{input}", header=T,check.names=FALSE,row.names=1); write.table(cpm(x, log=T, prior.count=${params.LTPM}), "{output}.tmp", append=T, sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp

    """

###ACHTUNG: doppio :, assumo che il vero input sia solo l'ultimo, non vedo più inputs nell'esecuzione
# GEP.count.tmm.gz GEP.count.exp_filter.tmm.gz: GEP.%.tmm.gz: GEP.%.gz
# 	r -e 'library(edgeR);\\
# x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
# y <- DGEList(counts=x);\\
# y <- calcNormFactors(y,method="TMM");\\
# write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
# write.table(cpm(y, normalized.lib.sizes=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule GEP_exp_filter:
    input:
        "GEP.{filter}.gz"
    output: 
        # "GEP.count.tmm.gz",
        # "GEP.count.exp_filter.tmm.gz",
        "GEP.{filter}.tmm.gz"
    shell: """
        r -e 'library(edgeR);\
        x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\
        y <- DGEList(counts=x);\
        y <- calcNormFactors(y,method="TMM");\
        write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\
        write.table(cpm(y, normalized.lib.sizes=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """


# GEP.count.ltmm.gz: GEP.count.gz
# 	r -e 'library(edgeR);\\
# x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
# y <- DGEList(counts=x);\\
# y <- calcNormFactors(y,method="TMM");\\
# write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
# write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp


rule GEP_count_ltmm:
    input: 
        "GEP.count.gz"
    output: 
        "GEP.count.ltmm.gz"
    shell: """
        r -e 'library(edgeR);\\
        x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\\
        y <- DGEList(counts=x);\\
        y <- calcNormFactors(y,method="TMM");\\
        write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
        write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

# GEP.count.exp_filter.ltmm.gz: GEP.count.exp_filter.gz
# 	r -e 'library(edgeR);\\
# x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
# y <- DGEList(counts=x);\\
# y <- calcNormFactors(y,method="TMM");\\
# write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
# write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule GEP_exp_filter_ltmm:
    input: 
        "GEP.count.exp_filter.gz"
    output: 
        "GEP.count.exp_filter.ltmm.gz"    
    shell: """
        r -e 'library(edgeR);\
        x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\
        y <- DGEList(counts=x);\
        y <- calcNormFactors(y,method="TMM");\
        write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\
        write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """


# GEP.count.exp_filter.luq.gz: GEP.count.exp_filter.gz
# 	r -e 'library(edgeR);\\
# x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
# y <- DGEList(counts=x);\\
# y <- calcNormFactors(y,method="upperquartile");\\
# write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
# write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule GEP_exp_filter_luq:
    input: 
        "GEP.count.exp_filter.gz"
    output: 
        "GEP.count.exp_filter.luq.gz"
    shell: """
        r -e 'library(edgeR);\\
        x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\\
        y <- DGEList(counts=x);\\
        y <- calcNormFactors(y,method="upperquartile");\\
        write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
        write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """



# GEP.count.exp_filter.lrle.gz: GEP.count.exp_filter.gz
# 	r -e 'library(edgeR);\\
# x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
# y <- DGEList(counts=x);\\
# y <- calcNormFactors(y,method="RLE");\\
# write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
# write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule GEP_exp_filter_lrle:
    input: 
        "GEP.count.exp_filter.gz"
    output: 
        "GEP.count.exp_filter.lrle.gz"
    shell: """
        r -e 'library(edgeR);\\
        x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\\
        y <- DGEList(counts=x);\\
        y <- calcNormFactors(y,method="RLE");\\
        write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
        write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

# GEP.count.exp_filter%vst.gz: GEP.count.exp_filter.gz metadata.txt
# 	r -e 'suppressMessages(library(DESeq2));\\
# x <- read.table("$<",  header=T,check.names=FALSE,row.names=1);\\
# z <- read.table("$^2", header=T,check.names=FALSE,row.names=1);\\
# dds<-DESeqDataSetFromMatrix(as.matrix(x),z,~1);\\
# vsd <- vst(dds, blind=TRUE);\\
# write.table(assay(vsd), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp

rule GEP_exp_filter_vst:
    input: 
        GEP="GEP.count.exp_filter.gz",
        metadata="metadata.txt"
    output: 
        "GEP.count.exp_filter{filter}vst.gz"
    shell: """
        r -e 'suppressMessages(library(DESeq2));\\
        x <- read.table("{input.GEP}",  header=T,check.names=FALSE,row.names=1);\\
        z <- read.table("{input.metadata}", header=T,check.names=FALSE,row.names=1);\\
        dds<-DESeqDataSetFromMatrix(as.matrix(x),z,~1);\\
        vsd <- vst(dds, blind=TRUE);\\
        write.table(assay(vsd), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

# GEP.count.exp_filter%rlog.gz: GEP.count.exp_filter.gz metadata.txt
# 	r -e 'suppressMessages(library(DESeq2));\\
# x <- read.table("$<",  header=T,check.names=FALSE,row.names=1);\\
# z <- read.table("$^2", header=T,check.names=FALSE,row.names=1);\\
# dds<-DESeqDataSetFromMatrix(as.matrix(x),z,~1);\\
# rld <- rlog(dds, blind=TRUE);\\
# write.table(assay(rld), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp


rule GEP_exp_filter_rlog:
    input: 
        GEP="GEP.count.exp_filter.gz",
        metadata="metadata.txt"
    output: 
        "GEP.count.exp_filter{filter}rlog.gz"
    shell: """
        r -e 'suppressMessages(library(DESeq2));\\
        x <- read.table("{input.GEP}",  header=T,check.names=FALSE,row.names=1);\\
        z <- read.table("{input.metadata}", header=T,check.names=FALSE,row.names=1);\\
        dds<-DESeqDataSetFromMatrix(as.matrix(x),z,~1);\\
        rld <- rlog(dds, blind=TRUE);\\
        write.table(assay(rld), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """



# GEP.count.tmm.gz.facotrs: GEP.count.tmm.gz
# 	@echo pass

rule GEP_factors:
    input:
        "GEP.count.tmm.gz"
    output:
        "GEP.count.tmm.gz.facotrs"
    shell:
        "echo pass"
#chatgpt dice che potrei anche scrivere shell:
#        "echo pass > /dev/null"
#per sopprimere l'output e non avere nulla che mi appaia nel terminale

###ACHTUNG anche qui c'è un doppio :. interpreto gli ultimi due come input dato che c'è un $^2 nella shell 
# GEP.count.cpm.fpkm.gz GEP.count.tmm.fpkm.gz: GEP.count.%.fpkm.gz: GEP.count.%.gz gene_len
# 	zcat $< | translate -a $^2 1 \
# 	| perl -wlane 'BEGIN{$$,="\t"} {$$g=shift(@F); $$l=shift(@F); if($$.==1){print $$g,@F}else{@F=map {$$_=($$_/$$l)*1000} @F; print $$g,@F}}' \
# 	| gzip > $@

rule GEP_count_fpkm:
    input: 
        GEP="GEP.count.{filter}.gz",
        GenLen="gene_len"
    output: 
        # "GEP.count.cpm.fpkm.gz",
        # "GEP.count.tmm.fpkm.gz" 
        "GEP.count.{filter}.fpkm.gz"
    shell: """
        zcat {input.GEP} | translate -a {input.GenLen} 1 \
        | perl -wlane 'BEGIN{$,="\t"} {$g=shift(@F); $l=shift(@F); if($.==1){print $g,@F}else{@F=map {$_=($_/$l)*1000} @F; print $g,@F}}' \
        | gzip > {output}
    """

###ACHTUNG di nuovo, essendoci due input, interpreto come input gli ultimi due 
# GEP.count.exp_filter.ltmm.lfpkm.gz GEP.count.exp_filter.lcpm.lfpkm.gz GEP.count.exp_filter.lcpm.lfpkm.gz: GEP.count.exp_filter.%.lfpkm.gz: GEP.count.exp_filter.%.gz gene_len
# 	zcat $< | translate -a $^2 1 \
# 	| perl -wlane 'BEGIN{$$,="\t"} {$$g=shift(@F); $$l=shift(@F); if($$.==1){print $$g,@F}else{@F=map {$$_=$$_ - log($$l/1000)/log(2)} @F; print $$g,@F}}' \
# 	| gzip > $@

rule GEP_count_lfpkm:
    input: 
        GEP1="GEP.count.exp_filter.{filter}.gz",
        GenLen="gene_len"
    output: 
    #    " GEP.count.exp_filter.ltmm.lfpkm.gz", 
    #    "GEP.count.exp_filter.lcpm.lfpkm.gz", 
    #    "GEP.count.exp_filter.lcpm.lfpkm.gz"
        " GEP.count.exp_filter.{filter}.lfpkm.gz"
    shell: """
        zcat {input.GEP1} | translate -a {input.GenLen} 1 \
        | perl -wlane 'BEGIN{$,="\t"} {$g=shift(@F); $l=shift(@F); if($.==1){print $g,@F}else{@F=map {$_=$_ - log($l/1000)/log(2)} @F; print $g,@F}}' \
        | gzip > {output}
    """
###ACHTUNG: questa era già commentata di suo. quindi non tradotta. 
#GEP.count.tmm.lfpkm.gz: GEP.count.tmm.fpkm.gz
#	perl -wlane 'BEGIN{$$,="\t"} {$$g=shift(@F); if($$.==1){print $$g,@F}else{@F=map {$$_=log($$_)/log(2)} @F; print $$g,@F}}' \
#	| gzip > $@

###ACHTUNG un solo input
# GEP.count.ltmm.log2r.gz GEP.count.lcpm.log2r.gz GEP.count.rpk.ltpm.log2r.gz: GEP.count.%.log2r.gz: GEP.count.%.gz
# 	r -e 'x <- read.table("$<", header=T,check.names=FALSE,row.names=1); write.table(x - rowMeans(x), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
# 	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
# 	rm $@.tmp


rule GEP_count_log2r:
    input: 
        "GEP.count.{filter}.gz"
    output: 
        # "GEP.count.ltmm.log2r.gz",
        # "GEP.count.lcpm.log2r.gz",
        # "GEP.count.rpk.ltpm.log2r.gz",
        "GEP.count.{filter}.log2r.gz"
    shell: """
        r -e 'x <- read.table("{input}", header=T,check.names=FALSE,row.names=1); write.table(x - rowMeans(x), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """  

# GEP.count%metadata.gz: GEP.count%gz metadata.txt
# 	(echo -e "GeneID\tsample\texp";zcat $< | matrix2tab ) | translate -a $^2 2 | gzip > $@

###non so perché sia rosso. idem per quello dopo 
rule GEP_count_metadata:
    input: 
        GEP="GEP.count{filter}gz",
        metadata="metadata.txt"
    output: 
        "GEP.count{filter}metadata.gz"
    shell:"""
        (echo -e "GeneID\tsample\texp"; zcat {input.GEP} | matrix2tab ) | translate -a {input.metadata} 2 | gzip > {output}
    """

###ACHTUNG un solo input
# GEP.count.whaterfall.gz GEP.count.frac.whaterfall.gz GEP.count.cpm.whaterfall.gz GEP.count.tmm.whaterfall.gz GEP.count.cpm.fpkm.whaterfall.gz: %.whaterfall.gz: %.gz
# 	zcat $< | matrix2tab | bsort -k2,2 -k3,3gr -S5% | repeat_group_pipe '\
# 		enumerate_rows -r \
# 		| enumerate_rows -n -r\
# 		| bawk '\''BEGIN{T=0} {T+=$$$$3; print $$$$0,T}'\'' \
# 	' 2 | gzip > $@

###qui ho lasciato gli $ come erano
rule GEP_whaterfall:
    input: 
        "{filter}.gz"
    output: 
        # "GEP.count.whaterfall.gz",
        # "GEP.count.frac.whaterfall.gz", 
        # "GEP.count.cpm.whaterfall.gz", 
        # "GEP.count.tmm.whaterfall.gz",
        # "GEP.count.cpm.fpkm.whaterfall.gz",
        "{filter}.whaterfall.gz"
    shell: """
        zcat {input} | matrix2tab | bsort -k2,2 -k3,3gr -S5% | repeat_group_pipe '\
        enumerate_rows -r \
        | enumerate_rows -n -r\
        | bawk '\''BEGIN{T=0} {T+=$$3; print $$0,T}'\'' \
        ' 2 | gzip > {output}
    """

# %.whaterfall.matrix.gz: %.whaterfall.gz
# 	bawk '{print $$order,$$sample,$$cumulat_count}' $< | tab2matrix -r count_captured_by_top_genes | gzip > $@

rule whaterfall_matrix:
    input: 
        "{filter}.whaterfall.gz"
    output: 
        "{filter}.whaterfall.matrix.gz"
    shell: 
        "bawk '{print $order,$sample,$cumulat_count}' {input} | tab2matrix -r count_captured_by_top_genes | gzip > {output}"

# GEP.count.frac.gz: GEP.count.gz
# 	zcat $< | normalize_columns -H -C | gzip > $@

rule GEP_count_frac:
    input: 
        "GEP.count.gz"
    output: 
        "GEP.count.frac.gz"
    shell: 
        "zcat {input} | normalize_columns -H -C | gzip > {output}"


# GEP.count.%.genome_plot.gz: GEP.count.%.gz /data/bioinfotree/task/annotations/dataset/gencode/hsapiens/29/basic.annotation.gene_positions
# 	zcat $< | matrix2tab | translate -a -d -j <(bawk '{print $$4,$$1,$$3-(($$3-$$2)/2)}' $^2 | uniq) 1 | bawk '{print $$4,$$2,$$3,$$5}' | uniq | bsort -k2,2V -k3,3n -S5% | gzip > $@

rule GEP_genome_plot:
    input: 
        GEP="GEP.count.{filter}.gz",
        GenePositions="/data/bioinfotree/task/annotations/dataset/gencode/hsapiens/29/basic.annotation.gene_positions"
    output: 
        "GEP.count.{filter}.genome_plot.gz"
    shell:
        "zcat {input.GEP} | matrix2tab | translate -a -d -j <(bawk '{print $4,$1,$3-(($3-$2)/2)}' {input.GenePositions} | uniq) 1 | bawk '{print $4,$2,$3,$5}' | uniq | bsort -k2,2V -k3,3n -S5% | gzip > {output}"

rule get_lfpkm:
    input:
        gep = "GEP.count.exp_filter.{normalization}.gz", 
        gene_len= "gene_len"
    output:
        "GEP.count.exp_filter.{normalization}.lfpkm.gz"
    shell:"""
        zcat {input.gep} | translate -a {input.gene_len} 1 \
	    | perl -wlane 'BEGIN{{$,="\t"}} {{$g=shift(@F); $l=shift(@F); if($.==1){{print $g,@F}}else{{@F=map {{$_=$_ - log($l/1000)/log(2)}} @F; print $g,@F}}}}' \
	    | gzip > {output}
    """

rule expression_add:
    input:
        "GEP.count.{normalization}.metadata.gz"
    output:
        "GEP.count.{normalization}.metadata.exp_genes_condition.gz"
    shell:"""
        bawk '{{print $condition ";" $GeneID, $exp}}' {input} | unhead | bsort -k1,1 | stat_base -g -a -m -l | tr ";" "\t" | gzip > {output}
    """

rule matrix_exp_condition:
    input:
        "GEP.{normalization}.metadata.exp_genes_condition.gz"
    output:
        "GEP.{normalization}.metadata.exp_genes_condition.matrix.gz"
    shell:"""
        bawk '{{print $gene,$condition,$avg_exp}}' {input} | tab2matrix -r GeneID | gzip > {output}
    """
