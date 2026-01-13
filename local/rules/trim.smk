# =============================================================================
# FASTP (QC AND TRIMMING)
# =============================================================================

if config["TRIMMER"] == "fastp" and config["LAYOUT"] == "SINGLE":
    ruleorder: fastp_se > trim_galore_se > trim_galore_pe
elif config["TRIMMER"] == "fastp" and config["LAYOUT"] == "PAIRED":
    ruleorder: fastp_pe > trim_galore_pe > trim_galore_se
elif config["TRIMMER"] == "trim_galore" and config["LAYOUT"] == "SINGLE":
    ruleorder: trim_galore_se > fastp_se > fastp_pe
elif config["TRIMMER"] == "trim_galore" and config["LAYOUT"] == "PAIRED":
    ruleorder: trim_galore_pe > fastp_pe > fastp_se

rule fastp_se:
    input:
        "fastq/{sample}_R1.fastq.gz"
    output:
        trimmed="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        #unpaired1="fastq/fastq_trimmed/se/{sample}.u1.fastq",
        #merged="fastq/fastq_trimmed/pe/{sample}.merged.fastq",
        #failed="fastq/fastq_trimmed/pe/{sample}.failed.fastq",
        html="fastq/fastq_trimmed/{sample}.html",
        json="fastq/fastq_trimmed/{sample}.json"
    threads: 6
    log:
        "fastq/fastq_trimmed/{sample}.log.txt"
    params:
        #adapters_r1="--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" 
    conda: "transcript_env.yaml"
    wrapper:
        "v3.3.6/bio/fastp"

rule fastp_pe:
    input:
        sample=["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"]
    output:
        trimmed=["fastq/fastq_trimmed/{sample}_R1.fastq.gz", "fastq/fastq_trimmed/{sample}_R2.fastq.gz"],
        html="fastq/fastq_trimmed/{sample}.html",
        json="fastq/fastq_trimmed/{sample}.json"
    log:
        "fastq/fastq_trimmed/{sample}.log"
    conda: "transcript_env.yaml"
    params:
        extra=""
    threads: 6
    wrapper:
        "v3.3.6/bio/fastp"

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
        "fastq/fastq_trimmed/{sample}_R1.fastq.gz"
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores=config["CORES"]
    conda: "transcript_env.yaml"
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
        fastq_read1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        fastq_read2="fastq/fastq_trimmed/{sample}_R2.fastq.gz"
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores=config["CORES"]
    conda: "transcript_env.yaml"
    shell:
        "mkdir -p `dirname {output}`; "
        "trim_galore -j {params.cores} "
        "-o trimgalore "
        "{params.trim_galore_params} --paired "
        "{input}; "
        "mv trimgalore/{wildcards.sample}_R1_val_1.fq.gz {output.fastq_read1}; "
        "mv trimgalore/{wildcards.sample}_R2_val_2.fq.gz {output.fastq_read2}"