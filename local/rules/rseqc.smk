COUNT_REF_BED=GENCODE_DIR+"/"+config["GENCODE"]["ANNOTATION"]+".bed"
RSEQC_REF_BED=GENCODE_DIR+"/basic.annotation.bed"

# rseqc/%.infer_experiment.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam $(COUNT_REF_BED)
# 	mkdir -p $$(dirname $@)
# 	infer_experiment.py -r $^2 -i $< > $@
rule infer_experiment:
    input:
        bam = "star/{sample}.bam",
        ref_bed = COUNT_REF_BED
    output:
        "rseqc/{sample}.infer_experiment.txt"
    shell:
        "mkdir -p `dirname {output}`; "
        "infer_experiment.py -r {input.ref_bed} -i {input.bam} > {output}"


#Per ora ho fatto tutte le regole senza la directory star_{genome} perchè non avevamo ancora capito bene come integrarlo

# rseqc/%.bam_stat.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam
# 	mkdir -p $$(dirname $@)
# 	bam_stat.py -i $< > $@
rule get_bam_stat:
    input:
        "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/{sample}.bam_stat.txt"
    shell:
        """
        mkdir -p `dirname {output}`;
        bam_stat.py -i {input} > {output}
        """
#non sono sicurissima della prima linea di comando

# #Example 1
# #Fraction of reads explained by "1++,1--,2+-,2-+": 0.4992
# #Fraction of reads explained by "1+-,1-+,2++,2--": 0.5008
# #Fraction of reads explained by other combinations: 0.0000
# #Conclusion: We can infer that this is NOT a strand specific because 50% of reads can be explained by “1++,1–,2+-,2-+”, while the other 50% can be explained by “1+-,1-+,2++,2–”.
# #
# #Example 2:
# #Fraction of reads explained by "1+-,1-+,2++,2--": 0.0356
# #Fraction of reads explained by other combinations: 0.0000
# #Conclusion: We can infer that this is a strand-specific RNA-seq data. strandness of read1 is consistent with that of gene model, while strandness of read2 is opposite to the strand of reference gene model.
# #
# #Example 3:
# #Fraction of reads explained by "++,--": 0.9840
# #Fraction of reads explained by "+-,-+": 0.0160
# #Fraction of reads explained by other combinations: 0.0000
# #Conclusion: This is single-end, strand specific RNA-seq data. Strandness of reads are concordant with strandness of reference gene.

# rseqc/%.junctionSaturation_plot.r: STAR/%.STAR/Aligned.sortedByCoord.out.bam $(RSEQC_REF_BED)
# 	mkdir -p $$(dirname $@)
# 	touch $@
# 	junction_saturation.py -i $< -r $^2 -o rseqc/$*
rule get_junction_saturation:
    input:
        bam = "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam",
        ref_bed = COUNT_REF_BED
    output:
        "rseqc/{sample}.junctionSaturation_plot.r"
    shell:
        """
        mkdir -p `dirname {output}`;
        touch {output};
        junction_saturation.py -i {input.bam} -r {input.ref_bed} -o rseqc/{wildcards.sample}
        """

# rseqc/%.saturation.r: STAR/%.STAR/Aligned.sortedByCoord.out.bam $(RSEQC_REF_BED) rseqc/%.infer_experiment.txt
# 	mkdir -p $$(dirname $@)
# 	touch $@
# 	RPKM_saturation.py -r $^2 -i $< -o rseqc/$*\
# 		 $$(tr ":" "\t" < $^3 | bawk '$$2'  | bsort -k2,2n | tail -n 1 | perl -lne 'm/([^\s]+)\t/; print "--strand $1" if $$1!="determine"')
rule get_saturation:
    input:
        bam = "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam",
        ref_bed = COUNT_REF_BED,
        infer_exp = "rseqc/{sample}.infer_experiment.txt"
    output:
        "rseqc/{sample}.saturation.r"
    shell:
        """
        mkdir -p `dirname {output}`;
        touch {output};
        RPKM_saturation.py -r {input.ref_bed} -i {input.bam} -o rseqc/{wildcards.sample}\
            $(tr ":" "\t" < {input.infer_exp} | bawk '$2'  | bsort -k2,2n | tail -n 1 | perl -lne 'm/([^\s]+)\t/; print "--strand $1" if $1!="determine"')
        """
        
# rseqc/%.pos.DupRate.xls: STAR/%.STAR/Aligned.sortedByCoord.out.bam
# 	mkdir -p $$(dirname $@)
# 	read_duplication.py -i $< -o rseqc/$*
rule get_dup_rate:
    input:
        "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/{sample}.pos.DupRate.xls"
    shell:
        """
        mkdir -p `dirname {output}`;
        read_duplication.py -i {input} -o rseqc/{wildcards.sample}
        """

# rseqc/%.geneBodyCoverage.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam $(GENCODE_DIR)/rseqc.HouseKeepingGenes.bed.gz STAR/%.STAR/Aligned.sortedByCoord.out.bam.bai
# 	mkdir -p $$(dirname $@)
# 	docker run -u `id -u`:`id -g` --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0  bash -c "cd $(PWD); geneBody_coverage.py -i $< -r <(zcat $^2) -o rseqc/$*"
# 	sed -i 's|Aligned.sortedByCoord.out|$*|' $@
rule get_gene_coverage:
    input:
        bam = "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam",
        house_keepers = "{GENCODE_DIR}/rseqc.HouseKeepingGenes.bed.gz",
        bai = "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam.bai"
    output:
        "rseqc/{sample}.geneBodyCoverage.txt"
    params:
        docker_data_dir = config["DOCKER_DATA_DIR"],
        scratch_dir = config["TMPDIR"]
    shell:
        """
        mkdir -p `dirname {output}`;
        docker run -u `id -u`:`id -g` --rm -v {params.docker_data_dir}:{{params.docker_data_dir} -v {params.scratch_dir}:{params.scratch_dir} quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0  bash -c "cd $(PWD); geneBody_coverage.py -i {input.bam} -r <(zcat {input.house_keepers}) -o rseqc/{wildcards.sample}"
        sed -i 's|Aligned.sortedByCoord.out|{wildcards.sample}|' {output}
        """

# rseqc/%.inner_distance.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam
# 	mkdir -p $$(dirname $@)
# 	inner_distance.py -i $< -o rseqc/$* -r $(RSEQC_REF_BED)
rule get_inner_distance:
    input:
        "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/{sample}.inner_distance.txt"
    params:
        bed_ref = RSEQC_REF_BED
    shell:
        """
        mkdir -p `dirname {output}`;
        inner_distance.py -i {input} -o rseqc/{wildcards.sample} -r {params.bed_ref}
        """

# ALL.skewness: $(addprefix ./rseqc/$(FASTQ_FILTERING)/, $(addsuffix .geneBodyCoverage.txt,$(SAMPLES)))
# 	matrix_reduce './rseqc/$(FASTQ_FILTERING)/*.geneBodyCoverage.txt' \
# 	| fasta2tab | grep -v Percentile | cut -f 1,3- | tab2fasta | tr "\t" "\n" | fasta2tab | stat_base -o -g -k > $@
rule ALL_skewness:
    input:
        expand("./rseqc/{filter}/{samples}.geneBodyCoverage.txt", filter = config["FASTQ_FILTERING"], samples = FASTQ_SAMPLES)
    output:
        "ALL.skewness"
    shell:
        """
        matrix_reduce '{input}' \
        | fasta2tab | grep -v Percentile | cut -f 1,3- | tab2fasta | tr "\t" "\n" | fasta2tab | stat_base -o -g -k > {output}
        """

# rseqc/%.read_distribution.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam $(RSEQC_REF_BED)
# 	mkdir -p `dirname $@`
# 	read_distribution.py  -i $< -r $^2 > $@
rule get_read_distribution:
    input:
        bam = "STAR/{sample}.STAR/Aligned.sortedByCoord.out.bam",
        bed_ref = RSEQC_REF_BED
    output:
        "rseqc/{sample}.read_distribution.txt"
    shell:
        """
        mkdir -p `dirname {output}`;
        read_distribution.py  -i {input.bam} -r {input.bed_ref} > {output}
        """

# rseqc/$(FASTQ_FILTERING)/ALL.read_distribution.tagskb_matrix: $(addprefix rseqc/$(FASTQ_FILTERING)/,$(addsuffix .read_distribution.txt, $(SAMPLES))) 
# 	matrix_reduce 'rseqc/$(FASTQ_FILTERING)/*.read_distribution.txt' | fasta2tab | perl -ne 's/_S\d+(\s)/\1/; print if !m/===/ and !m/Group/ and !m/Total/' | perl -lpe 's/\s+/\t/g' | cut -f 1,2,5 | tab2matrix > $@
rule read_distribution_matrix:
    input:
        expand("rseqc/{filter}/{samples}.read_distribution.txt", filter = config['FASTQ_FILTERING'], samples = FASTQ_SAMPLES)
    output:
        "rseqc/{FASTQ_FILTERING}/ALL.read_distribution.tagskb_matrix"
    shell:
        """
        matrix_reduce '{input}' | fasta2tab | perl -ne 's/_S\d+(\s)/\1/; print if !m/===/ and !m/Group/ and !m/Total/' \
        | perl -lpe 's/\s+/\t/g' | cut -f 1,2,5 | tab2matrix > {output}
        """

# rseqc/$(FASTQ_FILTERING)/ALL.read_distribution.tagskb_tab_norm: $(addprefix rseqc/$(FASTQ_FILTERING)/,$(addsuffix .read_distribution.txt, $(SAMPLES)))
# 	matrix_reduce 'rseqc/$(FASTQ_FILTERING)/*.read_distribution.txt' | fasta2tab | perl -lne 'BEGIN{$$,="\t"} $$T=$$1 if m/Total Tags\s+(\d+)/; s/_S\d+(\s)/\1/; s/\s+/\t/g; @F=split("\t",$$_); print $$F[0],$$F[1],$$F[4],$$F[4]/$$T if !m/===/ and !m/Group/ and !m/Total/' > $@
rule norm_read_distribution_matrix:
    input:
        expand("rseqc/{filter}/{samples}.read_distribution.txt", filter = config['FASTQ_FILTERING'], samples = FASTQ_SAMPLES)
    output:
        "rseqc/{FASTQ_FILTERING}/ALL.read_distribution.tagskb_tab_norm"
    shell:
        """
        matrix_reduce '{input}' | fasta2tab | perl -lne 'BEGIN{{$,="\t"}} $T=$1 if m/Total Tags\s+(\d+)/; s/_S\d+(\s)/\1/; s/\\s+/\t/g; @F=split("\t",$_); print $F[0],$F[1],$F[4],$F[4]/$T if !m/===/ and !m/Group/ and !m/Total/' > {output}
        """
     
#Questo non sono bene sicura a cosa serva e come tradurlo 
# .META:	ALL.read_distribution.tagskb_tab_norm
# 	1	protocol	Access
# 	2	sample		1018
# 	3	seq_type	3UTR_Exons
# 	4	Tags_Kb 	89.48
# 	5	Tags_Kb_norm 	89.48



# #######################################
# #
# #	rseqc param
# #

# RSEQC_REF_BED=$(GENCODE_DIR)/basic.annotation.bed


# TRIM_GALORE_PARAM = --stringency=3

# #######################################

