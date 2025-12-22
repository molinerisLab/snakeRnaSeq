##################
# Generic rules  #
##################

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

rule tab2xls:
    input: 
        "{file}"
    output: 
        "{file}.xls"
    shell: 
        "tab2xls < {input} > {output}"

# -------------- #
# Multiqc report #
# -------------- #

rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        html="fastqc/{fastq_filtering}/{sample}_{read}_fastqc.html",
        zip="fastqc/{fastq_filtering}/{sample}_{read}_fastqc.zip"
    threads: 2
    wrapper:
        "v3.3.6/bio/fastqc"

rule multiqc_fastq:
    input:
        fastqc_html=expand("fastqc/{fastq_filtering}/{s}_{p}_fastqc.html", s=SAMPLES, p=("R1","R2"), fastq_filtering=config["FASTQ_FILTERING"]),
        fastqc_zip=expand("fastqc/{fastq_filtering}/{s}_{p}_fastqc.zip",  s=SAMPLES, p=("R1","R2"), fastq_filtering=config["FASTQ_FILTERING"])
    output:
        "multiqc_report.html"
    params:
        data_dir = RAW_DATA_DIR
    shell:"""
        multiqc -f -n {output} {params.data_dir}
    """

rule multiqc_report_rseqc:
    input:
        expand("rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("rseqc/{sample}.infer_experiment.txt", sample=SAMPLES),
        expand("rseqc/{sample}.junctionSaturation_plot.r", sample=SAMPLES),
        expand("rseqc/{sample}.pos.DupRate.xls", sample=SAMPLES),
        expand("rseqc/{sample}.read_distribution.txt", sample=SAMPLES),
        expand("rseqc/{sample}.bam_stat.txt", sample=SAMPLES),
        expand("fastqc/{fastq_filtering}/{sample}_R1_fastqc.html", sample=SAMPLES, fastq_filtering=config["FASTQ_FILTERING"])
    output:
        "multiqc_report.rseqc.html"
    params:
        data_dir = RAW_DATA_DIR
    shell: """
        multiqc -f -n {output} {params.data_dir}
    """

 rule multiqc_alignment:
     input:
         fastqc_html=expand("fastq/{fastq}_fastqc.html", fastq=FASTQ_FILES),
         fastqc_zip=expand("fastq/{fastq}_fastqc.zip", fastq=FASTQ_FILES),
         bam=expand("bam/{sample}.bam", sample=SAMPLES),
         bai=expand("bam/{sample}.bam.bai", sample=SAMPLES)
     output:
         report="multiqc_report.alignment.html",
         star="multiqc_report.alignment_data/multiqc_star.txt"
     shell:
         "multiqc -f -n {output} ."

rule bam2cram:
    input: 
        "{file}.bam"
    output: 
        "{file}.cram"
    params:
        threads = config["CORES"],
        genome = GENCODE_GENOME_FASTA
    shell: 
        "samtools view -@ {params.threads} -T {params.genome} -C -o {output} {input}"

rule add_header:
    input: 
        "{path}.gz"
    output: 
        "{path}.header_added.gz"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; zcat {input} ) | gzip > {output}"

rule header_add:
    input: 
        "{file}"
    output: 
        "{file}.header_added"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; cat {input} ) > {output}"

rule get_bai:
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.bai"
    shell: 
        "samtools index {input}"

rule get_bam_id:
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.id"
    shell: 
        "samtools view {input} | cut -f 1 | bsort -S8% | uniq > {output}"

rule get_bw:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.bw"
    params:
        threads = config["CORES"],
        bin_size = "BIGWIG_BIN_SIZE"
    shell: 
        "bamCoverage --binSize={params.bin_size} -b {input.bam} -o {output} --numberOfProcessors={params.threads}"


if not 'BIGWIG_BIN_SIZE' in globals():
    BIGWIG_BIN_SIZE = 5

rule get_norm_bw:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.norm.bw"
    params:
        threads = config["CORES"],
        bin_size = BIGWIG_BIN_SIZE
    shell: 
        "bamCoverage --binSize={params.bin_size}  --normalizeUsing=CPM -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

rule get_bw_binsize:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.{bin_size}.bw"
    params:
        threads = config["CORES"],
        bin_size = BIGWIG_BIN_SIZE
    shell: 
        "bamCoverage --binSize={params.bin_size} -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

rule get_norm_bw_binsize:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.{bin_size}.norm.bw"
    params:
        threads = config["CORES"],
        bin_size = "BIGWIG_BIN_SIZE"
    shell: 
        "bamCoverage --binSize={params.bin_size}  --normalizeUsing=CPM -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

rule bam2bed:
    input: 
        "{file}.bam"
    output: 
        "{file}.bed"
    shell: 
        "bedtools bamtobed -splitD < {input} | bsort -k1,1V -k2,2n > {output}"

rule get_bedgraph_ranges:
    input: 
        "{file}.bedGraph"
    output: 
        "{file}.bedGraph.ranges"
    shell: 
        "bawk '{{print $1,$2; print $1,$3}}' {input} | stat_base -o -g -b | bawk {{print $1,0,$2}}' > {output}"

if not "GFF3_MIN_WINDOW_READS" in globals():
    GFF3_MIN_WINDOW_READS = 0

rule get_gff3:
    input: 
        bedGraph = "{file}.bedGraph",
        ranges = "{file}.bedGraph.ranges"
    output: 
        "{file}.bedGraph.gff3"
    params:
        min_window_reads = GFF3_MIN_WINDOW_READS
    shell: 
        """
        bawk '$4 > {params.min_window_reads}' {input.bedGraph} \
        | repeat_group_pipe '\
            grep -w $1 {input.ranges} | append_each_row -B "##sequence-region";\
               grep -w $1 {input.ranges} | bawk '\''{{print $1,".","chromosome",$2,$3,$4,".","+-",".","."}}'\'';\
               bawk '\''{{print $1,".","reads",$2,$3,$4,".","+-",".","."}}'\''\
        ' 1 > {output}
        """

if not "BEDGRAPH_FILTER" in globals():
    BEDGRAPH_FILTER = 100

rule get_filtered_bed:
    input: 
        "{file}.bedGraph"
    output: 
        "{file}.bedGraph.{BEDGRAPH_FILTER}.merged.bed"
    params:
        bedgraph_filter = BEDGRAPH_FILTER
    shell:
        """
        "bawk '$4 > {params.bedgraph_filter} {input} | union -s --allow-duplicates | bawk '{{print $1";"$2";"$3, $4}}' \
        | expandsets 2 | stat_base -o -g -b | tr ";" "\t" > {output}
        """

rule get_fa:
    input: 
        "{file}.fastq.gz"
    output: 
        "{file}.fa.gz"
    shell: 
        "zcat {input} | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip > {output}"

rule tmm:
    input:
        "{file}.gz"
    output: 
        "{file}.tmm.gz"
    shell: """
        r -e 'library(edgeR);\
            x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\
            y <- DGEList(counts=x);\
            y <- calcNormFactors(y,method="TMM");\
            write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\
            write.table(cpm(y, normalized.lib.sizes=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

rule ltmm:
    input: 
        "{file}.gz"
    output: 
        "{file}.ltmm.gz"
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