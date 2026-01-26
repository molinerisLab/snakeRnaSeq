# =============================================================================
# 1. MULTIQC REPORTING
# =============================================================================

rule fastqc:
    """Run FastQC on raw sequencing reads."""
    input:
        "fastq/{sample}_R1.fastq.gz"
    output:
        html="fastqc/{sample}_{read}_fastqc.html",
        zip="fastqc/{sample}_{read}_fastqc.zip"
    threads: 2
    wrapper:
        "v3.3.6/bio/fastqc"

rule multiqc_fastq:
    """Aggregate FastQC results into a single MultiQC report."""
    input:
        fastqc_html=expand("fastqc/{s}_{p}_fastqc.html", s=SAMPLES, p=("R1","R2")),
        fastqc_zip=expand("fastqc/{s}_{p}_fastqc.zip",  s=SAMPLES, p=("R1","R2"))
    output:
        "multiqc_report.html"
    params:
        data_dir = RAW_DATA_DIR
    shell:"""
        multiqc -f -n {output} {params.data_dir}
    """

rule multiqc_report_rseqc:
    """Aggregate RSeQC metrics and FastQC into a specialized QC report."""
    input:
        expand("rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("rseqc/{sample}.infer_experiment.txt", sample=SAMPLES),
        expand("rseqc/{sample}.junctionSaturation_plot.r", sample=SAMPLES),
        expand("rseqc/{sample}.pos.DupRate.xls", sample=SAMPLES),
        expand("rseqc/{sample}.read_distribution.txt", sample=SAMPLES),
        expand("rseqc/{sample}.bam_stat.txt", sample=SAMPLES),
        expand("fastqc/{sample}_R1_fastqc.html", sample=SAMPLES)
    output:
        "multiqc_report.rseqc.html"
    params:
        data_dir = RAW_DATA_DIR
    shell: """
        multiqc -f -n {output} {params.data_dir}
    """

rule multiqc_alignment:
    """Aggregate alignment stats (STAR, BAM) and FastQC."""
    input:
        fastqc_html=expand("fastqc/{sample}_{read}_fastqc.html",read=("R1", "R2"), sample=SAMPLES),
        fastqc_zip=expand("fastqc/{sample}_{read}_fastqc.zip", read=("R1", "R2"), sample=SAMPLES),
        # Inject the ALIGNER variable into the path
        bam=expand("{aligner}/{sample}.bam", aligner=config["aligner"], sample=SAMPLES),
        bai=expand("{aligner}/{sample}.bam.bai", aligner=config["aligner"], sample=SAMPLES)
    output:
        report="multiqc_report.alignment.html",
        star="multiqc_report.alignment_data/multiqc_star.txt"
    shell:
        "multiqc -f -n {output.report} ."
