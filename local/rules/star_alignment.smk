def choose_fastq_according_to_genome(wildcards, mate):
    sample = wildcards['sample']  # Definisce 'sample' usando 'wildcards'
    if wildcards['genome'] == "GRCh":
        return f"fastq_trimmed/{sample}_R{mate}.fastq.gz"
    elif wildcards['genome'] == "CHM":
        return f"star_GRCh/{sample}_unmapped_R{mate}.fastq.gz"
    else:
        raise Exception(f"Genome not valid: {wildcards['genome']}")

rule star_pe_multi:
    input:
        fq1=lambda wildcards: choose_fastq_according_to_genome(wildcards, 1),
        fq2=lambda wildcards: choose_fastq_according_to_genome(wildcards, 2),
        idx=lambda wildcards: config['STAR']['INDEX'][wildcards['genome']],
    output:
        aln="star_{genome}/{sample}.bam",
        log="star_{genome}/{sample}.Log.out",
        sj="star_{genome}/{sample}.SJ.out.tab",
        # Uncomment the next line if you want to handle unmapped reads
        # unmapped=["star_{genome}/unmapped/{sample}_R1.fastq.gz", "star_{genome}/unmapped/{sample}_R2.fastq.gz"],
        #unmapped read filtered after, sice by default STAR report as unmapped partially mapped (i.e. mapped only one mate of a paired end read)
    log:
        "star_{genome}/{sample}.log",
    params:
        extra=lambda wildcards: f"--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --chimOutType WithinBAM {config['STAR']['OPTIONS']}",
    threads: 16,
    wrapper:
        "v3.3.6/bio/star/align"

rule linl_unmapped:
    input: "star_CHM/{sample}_unmapped_R{mate}.fastq.gz"
    output: "fastq_unmapped/{sample}_R{mate}.fastq.gz"
    shell: "ln {input} {output}"

rule all_fastq_unmapped:
    input:
        expand("fastq_unmapped/{sample}_R1.fastq.gz", sample=config['samples']),
        expand("fastq_unmapped/{sample}_R2.fastq.gz", sample=config['samples'])

rule generate_unmapped_R1:
    input:
        "{sample}.bam"
    output:
        "{sample}_unmapped_R1.fastq.gz"
    shell:"""
        samtools view -f 76 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 76=4+8+64 = read unmapped AND mate unmapped AND first in pair, i.e., discard reads that are unmapped but that have mate mapped

rule generate_unmapped_R2:
    input:
        "{sample}.bam"
    output:
        "{sample}_unmapped_R2.fastq.gz"
    shell: """
        samtools view -f 140 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 140=4+8+128 = read unmapped AND mate unmapped AND second in pair, i.e., discard reads that are unmapped but that have mate mapped

