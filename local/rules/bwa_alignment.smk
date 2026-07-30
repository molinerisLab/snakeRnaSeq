if config["LAYOUT"] == "PAIRED":
    ruleorder: align_pe_bwa > align_se_bwa
elif config["LAYOUT"] == "SINGLE":
    ruleorder: align_se_bwa > align_pe_bwa
DNA_SAMPLES = [s for s in SAMPLES if s.startswith("DNA")]


rule all_bwa:
    input:
        expand("aligned_bwa/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)


rule align_pe_bwa:
    input:
        trimmed_fastq1="Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
        trimmed_fastq2="Results/star/unmapped/{sample}_unmapped_R2.fastq.gz",
        idx=multiext(config["bwa_db"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam="aligned_bwa/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        reference=config["bwa_db"]          
    threads: config["CORES"]["bwa"]
    log:
        "aligned_bwa/logs/{sample}_bwa.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p aligned_bwa/logs
        T=$(mktemp -d)
        trap 'rm -rf "$T"' EXIT

        bwa-mem2 mem -t {threads} \
            {params.reference} {input.trimmed_fastq1} {input.trimmed_fastq2} 2> {log} \
          | samtools view -b - \
          | sambamba sort --tmpdir="$T" -t {threads} -o {output.bam} /dev/stdin
        """
#  -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA" \

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
