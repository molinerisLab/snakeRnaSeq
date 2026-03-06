# =============================================================================
# FASTP (QC AND TRIMMING)
# =============================================================================
from glob import glob

TRIMMER = config["TRIMMER"]
LAYOUT = config["LAYOUT"]


rule fastp_se:
    input:
        "fastq/{sample}_R1.fastq.gz"
    output:
        trimmed="fastq/fastq_trimmed/fastp/{sample}_R1.fastq.gz",
        html="fastq/fastq_trimmed/fastp/{sample}.html",
        json="fastq/fastq_trimmed/fastp/{sample}.json"
    threads: 6
    log:
        "fastq/fastq_trimmed/fastp/{sample}.log.txt"
    wrapper:
        "v3.3.6/bio/fastp"

rule fastp_pe:
    input:
        sample=["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"]
    output:
        trimmed=["fastq/fastq_trimmed/fastp/{sample}_R1.fastq.gz", "fastq/fastq_trimmed/fastp/{sample}_R2.fastq.gz"],
        html="fastq/fastq_trimmed/fastp/{sample}.html",
        json="fastq/fastq_trimmed/fastp/{sample}.json"
    log:
        "fastq/fastq_trimmed/fastp/{sample}.log"
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
        "fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz"
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
        fastq_read1="fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz",
        fastq_read2="fastq/fastq_trimmed/trimgalore/{sample}_R2.fastq.gz"
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


##################################
# Rule for linking fastq_trimmed #
##################################

if LAYOUT == "SINGLE":
    rule link_trimmed_se:
        input:
           r1=f"fastq/fastq_trimmed/{TRIMMER}/{{sample}}_R1.fastq.gz"
        output:
           tr1="fastq/fastq_trimmed/{sample}_R1.fastq.gz"
        wildcard_constraints:
            sample="[^/]+"
        shell:
            "ln -srf {input.r1} {output.tr1}"

elif LAYOUT == "PAIRED":
    rule link_trimmed_pe:
        input:
            r1=f"fastq/fastq_trimmed/{TRIMMER}/{{sample}}_R1.fastq.gz",
            r2=f"fastq/fastq_trimmed/{TRIMMER}/{{sample}}_R2.fastq.gz"
        output:
            tr1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
            tr2="fastq/fastq_trimmed/{sample}_R2.fastq.gz"
        wildcard_constraints:
            sample="[^/]+"
        shell:
            "ln -srf {input.r1} {output.tr1}; ln -srf {input.r2} {output.tr2}"