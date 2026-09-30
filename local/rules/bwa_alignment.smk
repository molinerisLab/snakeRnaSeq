if config["LAYOUT"] == "PAIRED":
    ruleorder: align_pe_bwa > align_se_bwa
elif config["LAYOUT"] == "SINGLE":
    ruleorder: align_se_bwa > align_pe_bwa


rule all_bwa:
    input:
        expand(
            "aligned_bwa_T30/{sample}/Aligned.sortedByCoord.out.bam",
            sample=BWA_SAMPLES
        )


rule align_pe_bwa:
    input:
        trimmed_fastq1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        trimmed_fastq2="fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        idx=multiext(
            config["bwa_db"],
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    output:
        bam="aligned_bwa_T30/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        reference=config["bwa_db"],
        min_score=config["BWA"]["MIN_ALIGNMENT_SCORE"],
    threads:
        config["CORES"]["bwa"]
    log:
        "aligned_bwa_T30/logs/{sample}_bwa.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p aligned_bwa_T30/{wildcards.sample}
        mkdir -p aligned_bwa_T30/logs
        exec 2> {log}

        bwa_tmpdir=$(mktemp -d)
        trap 'rm -rf -- "$bwa_tmpdir"' EXIT

        bwa-mem2 mem \
            -t {threads} \
            -T {params.min_score} \
            -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' \
            {params.reference} \
            {input.trimmed_fastq1} \
            {input.trimmed_fastq2} \
          | samtools view -b - \
          | sambamba sort \
                --tmpdir="$bwa_tmpdir" \
                -t {threads} \
                -o {output.bam} \
                /dev/stdin
        """

rule align_se_bwa:
    input:
        trimmed_fastq="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        idx=multiext(config["bwa_db"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam="aligned_bwa/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        reference=config["bwa_db"]
    threads: config["CORES"]["bwa"]
    conda:
        "transcript_env.yaml"
    log:
        "aligned_bwa/logs/{sample}_bwa.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p aligned_bwa/logs
        T=$(mktemp -d)
        trap 'rm -rf "$T"' EXIT

        bwa-mem2 mem -t {threads} \
            -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA" \
            {params.reference} {input.trimmed_fastq} 2> {log} \
          | samtools view -b - \
          | sambamba sort --tmpdir="$T" -t {threads} -o {output.bam} /dev/stdin
        """

CHM13_BWA_DB = config["CHM13_DB"]


rule all_bwa_chm13:
    input:
        expand(
            "aligned_bwa_CHM13_T30/{sample}/Aligned.sortedByCoord.out.bam",
            sample=BWA_SAMPLES,
        )


rule align_pe_bwa_chm13:
    input:
        fq1="aligned_bwa_T30/unmapped/{sample}_unmapped_R1.fastq.gz",
        fq2="aligned_bwa_T30/unmapped/{sample}_unmapped_R2.fastq.gz",
        idx=multiext(
            CHM13_BWA_DB,
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    output:
        bam="aligned_bwa_CHM13_T30/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        reference=CHM13_BWA_DB,
        min_score=config["BWA"]["MIN_ALIGNMENT_SCORE"],
    threads:
        config["CORES"]["bwa"]
    log:
        "aligned_bwa_CHM13_T30/logs/{sample}_bwa.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p aligned_bwa_CHM13_T30/{wildcards.sample}
        mkdir -p aligned_bwa_CHM13_T30/logs
        exec 2> {log}

        bwa_tmpdir=$(mktemp -d)
        trap 'rm -rf -- "$bwa_tmpdir"' EXIT

        bwa-mem2 mem \
            -t {threads} \
            -T {params.min_score} \
            -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' \
            {params.reference} \
            {input.fq1} \
            {input.fq2} \
          | samtools view -b - \
          | sambamba sort \
                --tmpdir="$bwa_tmpdir" \
                -t {threads} \
                -o {output.bam} \
                /dev/stdin

        samtools quickcheck -v {output.bam}
        """