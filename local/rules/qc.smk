# =============================================================================
# 1. MULTIQC REPORTING
# =============================================================================

QC_READS = ("R1", "R2") if config["LAYOUT"] == "PAIRED" else ("R1",)
QC_STAR_LOG = {
    "single_pass": "Results/star/{sample}/Log.final.out",
    "two_pass_manual": "Results/pass2/{sample}/Log.final.out",
    "two_pass_basic": "star_2pass/{sample}/Log.final.out",
}[config["STAR_MODE"]] if config["aligner"] == "star" else None
QC_TRIMMING_REPORTS = (
    expand("fastq/fastq_trimmed/fastp/{sample}.json", sample=SAMPLES)
    if config["TRIMMER"] == "fastp"
    else expand(
        "fastq/fastq_trimmed/trimgalore/{sample}_{read}.fastq.gz_trimming_report.txt",
        sample=SAMPLES, read=QC_READS,
    )
)

rule fastqc:
    """Run FastQC on raw sequencing reads."""
    input:
        "fastq/{sample}_{read}.fastq.gz"
    output:
        html="fastqc/{sample}_{read}_fastqc.html",
        zip="fastqc/{sample}_{read}_fastqc.zip"
    wildcard_constraints:
        read="R[12]"
    threads: 2
    log:
        "logs/fastqc/{sample}_{read}.log"
    wrapper:
        "v3.3.6/bio/fastqc"

rule multiqc_fastq:
    """Aggregate raw FastQC results and the selected trimmer's reports."""
    input:
        fastqc_html=expand("fastqc/{s}_{p}_fastqc.html", s=SAMPLES, p=QC_READS),
        fastqc_zip=expand("fastqc/{s}_{p}_fastqc.zip", s=SAMPLES, p=QC_READS),
        trimming_reports=QC_TRIMMING_REPORTS
    output:
        "multiqc_report.html"
    conda:
        "transcript_env.yaml"
    log:
        "logs/multiqc/fastq.log"
    shell:"""
        multiqc -f -n {output:q} {input.fastqc_zip:q} {input.trimming_reports:q} > {log:q} 2>&1
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
        expand("fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=QC_READS)
    output:
        "multiqc_report.rseqc.html"
    conda:
        "transcript_env.yaml"
    log:
        "logs/multiqc/rseqc.log"
    shell: """
        multiqc -f -n {output:q} {input:q} > {log:q} 2>&1
    """

rule samtools_flagstat:
    input:
        "Results/bam/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "qc/alignment/{sample}.flagstat.txt"
    conda:
        "transcript_env.yaml"
    log:
        "logs/samtools_flagstat/{sample}.log"
    shell:
        "samtools flagstat {input:q} > {output:q} 2> {log:q}"

rule multiqc_alignment:
    """Aggregate alignment stats (STAR, BAM) and FastQC."""
    input:
        fastqc_html=expand("fastqc/{sample}_{read}_fastqc.html", read=QC_READS, sample=SAMPLES),
        fastqc_zip=expand("fastqc/{sample}_{read}_fastqc.zip", read=QC_READS, sample=SAMPLES),
        flagstat=expand("qc/alignment/{sample}.flagstat.txt", sample=SAMPLES),
        star_logs=expand(QC_STAR_LOG, sample=SAMPLES) if QC_STAR_LOG else []
    output:
        report="multiqc_report.alignment.html"
    conda:
        "transcript_env.yaml"
    log:
        "logs/multiqc/alignment.log"
    shell:
        "multiqc -f -n {output.report:q} {input.fastqc_zip:q} {input.flagstat:q} {input.star_logs:q} > {log:q} 2>&1"
