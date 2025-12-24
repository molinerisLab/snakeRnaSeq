# =============================================================================
# FASTP (QC AND TRIMMING)
# =============================================================================

rule fastp_se:
    """
    Run fastp on Single-End (SE) reads for quality control and adapter trimming.
    """
    input:
        "fastq/{sample}_R1.fastq.gz"
    output:
        trimmed="fastp/{sample}_R1.fastq.gz",
        #unpaired1="fastp/se/{sample}.u1.fastq",
        #merged="fastp/pe/{sample}.merged.fastq",
        #failed="fastp/pe/{sample}.failed.fastq",
        html="fastp/{sample}.html",
        json="fastp/{sample}.json"
    threads: 8
    log:
        "fastp/{sample}.log.txt"
    params:
        adapters_r1="--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
    shell: """
        fastp --thread {threads}  --html {output.html} \
        {params.adapters_r1} \
        --in1 {input} --out1 {output.trimmed} --json {output.json} \
        2> {log}
        """

rule fastp_pe:
    """
    Run fastp on Paired-End (PE) reads for quality control and adapter trimming.
    """
    input:
        sample=["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"]
    output:
        trimmed=["fastp/{sample}_R1.fastq.gz", "fastp/{sample}_R2.fastq.gz"],
        #unpaired1="fastp/pe/{sample}.u1.fastq",
        #unpaired2="fastp/pe/{sample}.u2.fastq",
        #merged="fastp/pe/{sample}.merged.fastq",
        #failed="fastp/pe/{sample}.failed.fastq",
        html="fastp/{sample}.html",
        json="fastp/{sample}.json"
    log:
        "fastp/{sample}.log"
    params:
        adapters="--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        #extra="--merge"
        extra=""
    threads: 8
    log:
        "fastp/{sample}.log.txt"
    shell: """
        fastp --thread {threads}  --html {output.html} \
        {params.adapters_r1} {params.adapters_r2} \
        --in1 {input.sample[0]} --in2 {input.sample[1]} --out1 {output.trimmed[0]} --out2 {output.trimmed[1]} --json {output.json} \
        2> {log}
        """


# =============================================================================
# TRIM GALORE (ALTERNATIVE TRIMMING)
# =============================================================================

rule trim_galore_se:
    """
    Run Trim Galore on Single-End reads (wrapper around Cutadapt and FastQC).
    """
    input:
        "fastq/{sample}_R1.fastq.gz"
    output:
        "trimgalore/{sample}_R1.fastq.gz"
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores=config["CORES"]
    shell:
        "mkdir -p `dirname {output}`; "
        "trim_galore -j {params.cores} "
        "-o trimgalore "
        "{params.trim_galore_params} "
        "{input}; "
        "mv trimgalore/{wildcards.sample}_R1_trimmed.fq.gz {output}"

rule trim_galore_pe:
    """
    Run Trim Galore on Paired-End reads in paired mode.
    """
    input:
        fastq_read1="fastq/{sample}_R1.fastq.gz",
        fastq_read2="fastq/{sample}_R2.fastq.gz"
    output:
        fastq_read1="trimgalore/{sample}_R1.fastq.gz",
        fastq_read2="trimgalore/{sample}_R2.fastq.gz"
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores=config["CORES"]
    shell:
        "mkdir -p `dirname {output}`; "
        "trim_galore -j {params.cores} "
        "-o trimgalore "
        "{params.trim_galore_params} --paired "
        "{input}; "
        "mv trimgalore/{wildcards.sample}_R1_val_1.fq.gz {output.fastq_read1}; "
        "mv trimgalore/{wildcards.sample}_R2_val_2.fq.gz {output.fastq_read2}"