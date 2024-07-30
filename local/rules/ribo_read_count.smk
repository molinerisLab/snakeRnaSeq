# ------------------------------- #
# Correct ribosomal read counting #
# ------------------------------- #

rule split_bam_ribo:
	input:
		bam="{path}.bam",
		bam_idx="{path}.bam.bai",
		ribosome_bed=GENCODE_DIR+"/primary_assembly.annotation.rRNA_complete.bed"
	output:
		ribo_ex="{path}.ribo.ex.bam",
		ribo_in="{path}.ribo.in.bam",
		ribo_log="{path}.summary"
	conda:
		"../../local/env/rseqc_v5.0.1.yaml"
	shell:
		"mkdir -p `dirname {output}`; "
		"split_bam.py -i {input.bam} -r {input.ribosome_bed} -o {wildcards.path}.ribo > {wildcards.path}.summary"

ruleorder: featurecounts > split_bam_ribo
rule featurecounts:
	input:
		bam="{path}.bam",
		annotation_gtf=GENCODE_ANNOTATION_GTF
	output:
		counts="{path}.bam.featurecounts.count",
		summary="{path}.bam.featurecounts.count.summary"
	conda:
		"../../local/env/subread_v2.0.3.yaml"
	params:
		cores=config["CORES"],
		tmpdir=config["TMPDIR"]
	shell:
		"featureCounts "
		"{input.bam} "
		"-o {output.counts} "
		"-a {input.annotation_gtf} "
		#"-t exon "
		#"-g gene_name "
		#"-s 2 "
		#"-p -C " ==> paired-end layout
		"--tmpDir {params.tmpdir} "
		"-T {params.cores} "

rule featurecounts_ribo_ex:
	input:
		file_all=expand("bam/{sample}.srt2.dedup.ribo.ex.bam.featurecounts.count", sample=FASTQ_SAMPLES),
		file_translate=expand("bam/{sample}.srt2.dedup.ribo.ex.bam.featurecounts.count", sample=FASTQ_SAMPLES[0])
	output:
		"fastq.featurecounts.ribo.ex.count.gz"
	#container:
	#	"../../local/images/bit.wip-rnaseq.0.8.sif"
	shell:
		"matrix_reduce 'bam/*.srt2.dedup.ribo.ex.bam.featurecounts.count' -l '{input.file_all}' "
		"| grep -v '^#|^Geneid' "
		"| fasta2tab "
		"| bawk '{{print $2,$1,$8}}' "
		"| tab2matrix -r Geneid "
		"| translate -a <(cut -f -6 {input.file_translate} | unhead) 1 "
		"| gzip > {output}"

rule featurecounts_ribo_ex_summary:
	input:
		expand("bam/{sample}.srt2.dedup.ribo.ex.bam.featurecounts.count.summary", sample=FASTQ_SAMPLES)
	output:
		"fastq.featurecounts.ribo.ex.count.gz.summary_matrix"
	#container:
	#	"../../local/share/images/bit.wip-rnaseq.0.8.sif"
	#conda:
	#	"../../local/bioinfotree.yaml"
	shell:
		"matrix_reduce -t 'bam/*.srt2.dedup.ribo.ex.bam.featurecounts.count.summary' "
		"| grep -v Status "
		"| tab2matrix -r Sample > {output}"

rule featurecounts_ribo_ex_summary_matrix:
	input:
		"fastq.featurecounts.ribo.ex.count.gz.summary_matrix"
	output:
		"fastq.featurecounts.ribo.ex.count.gz.summary_matrix.reduced"
#	conda:
#		"../../local/env/bioinfotree.yaml"
	shell:
		"matrix2tab {input} "
		"| bawk '$2==\"Unassigned_Ambiguity\" || $2==\"Assigned\" || $2==\"Unassigned_NoFeatures\"' "
		"| tab2matrix -r Sample > {output}"

rule read_count:
	input:
		"{path}.bam"
	output:
		"{path}.bam.read_count"
	shell:
		"samtools view {input} "
		"| cut -f -1 "
		"| count > {output}"

rule read_count_unmap:
	input:
		"{path}.bam"
	output:
		"{path}.unmap.bam.read_count"
	shell:
		"samtools view -f 4 {input} "
		"| cut -f -1 "
		"| count > {output}"

rule read_count_multimap:
	input:
		"{path}.bam"
	output:
		"{path}.multi_map.bam.read_count"
	shell:
		"samtools view -F 4 {input} "
		"| bawk '$1~/^@/ || $12!=\"NH:i:1\"' "
		"| cut -f -1 "
		"| count > {output}"

rule read_count_uniq:
	input:
		"{path}.bam"
	output:
		"{path}.uniq_map.bam.read_count"
	shell:
		"samtools view -F 260 {input} "
		"| bawk '$1~/^@/ || $12==\"NH:i:1\"' "
		"| cut -f -1 "
		"| count > {output}"

rule usable_reads:
	input:
		read_count="{path}.bam.read_count",
		read_count_unmap="{path}.unmap.bam.read_count",
		read_count_ribo_in="{path}.ribo.in.bam.read_count",
		read_count_ribo_ex_multimap="{path}.ribo.ex.multi_map.bam.read_count",
		read_count_ribo_ex_uniq="{path}.ribo.ex.uniq_map.bam.read_count"
	output:
		"{path}.usable_reads"
	shell:
		"cat {input} | transpose > {output}"

rule usable_reads_all:
	input:
		ribo_ex_matrix="fastq.featurecounts.ribo.ex.count.gz.summary_matrix.reduced",
		usable_reads=expand("bam/{sample}.srt2.dedup.usable_reads", sample=FASTQ_SAMPLES)
	output:
		"usable_reads.txt"
	shell:
		"matrix_reduce -t 'bam/*.srt2.dedup.usable_reads' "
		"| translate -a -r {input.ribo_ex_matrix} 1"
		"| bawk 'BEGIN{{print \"sample\",\"tot\",\"unmap\",\"ribo\",\"non_ribo_multi_map\",\"non_ribo_uniq_map\",\"Assigned\",\"Unassigned_Ambiguity\",\"Unassigned_NoFeatures\"}} {{print}}' > {output}"
