
# ========================
# CONFIGURATION VARIABLES
# ========================
wildcard_constraints:
    sample="[^/]+",


TRIMMER = config["TRIMMER"]
LAYOUT = config["LAYOUT"]


def _trimmer_out(sample, read, trimmer=TRIMMER):
    """Return the trimmer-specific output path for a given sample and read (R1/R2)."""
    return f"fastq/fastq_trimmed/{trimmer}/{sample}_{read}.fastq.gz"


# =============================================================================
# FASTP (QC AND TRIMMING)
# =============================================================================

rule fastp_se:
    input:
        sample=["fastq/{sample}_R1.fastq.gz"],
    output:
        trimmed="fastq/fastq_trimmed/fastp/{sample}_R1.fastq.gz",
        html="fastq/fastq_trimmed/fastp/{sample}.html",
        json="fastq/fastq_trimmed/fastp/{sample}.json",
    threads: config["CORES"]["fastp"]
    params:
        extra=config["FASTP"]["extra"]
    log:
        "fastq/fastq_trimmed/fastp/{sample}.log",
    conda:
        "transcript_env.yaml"
    wrapper:
        "v3.3.6/bio/fastp"


rule fastp_pe:
    input:
        sample=["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"],
    output:
        trimmed=[
            "fastq/fastq_trimmed/fastp/{sample}_R1.fastq.gz",
            "fastq/fastq_trimmed/fastp/{sample}_R2.fastq.gz",
        ],
        html="fastq/fastq_trimmed/fastp/{sample}.html",
        json="fastq/fastq_trimmed/fastp/{sample}.json",
    log:
        "fastq/fastq_trimmed/fastp/{sample}.log",
    params:
        extra=config["FASTP"]["extra"] 
    threads: config["CORES"]["fastp"]
    conda:
        "transcript_env.yaml"
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
        "fastq/{sample}_R1.fastq.gz",
    output:
        trimmed="fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz",
        report="fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz_trimming_report.txt",
    threads: config["CORES"]["trimgalore"]
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        # Trim Galore estimates up to 3 * -j + 3 CPUs; -j must remain at least 1.
        cores=lambda wildcards, threads: max(1, (threads - 3) // 3),
        outdir="fastq/fastq_trimmed/trimgalore",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/fastq_trimmed/trimgalore/{sample}_trim_galore.log",
    shell:
        """
        mkdir -p {params.outdir:q}
        
        trim_galore -j {params.cores} \
            --basename {wildcards.sample:q}_R1 \
            -o {params.outdir:q} \
            {params.trim_galore_params} \
            {input:q} > {log:q} 2>&1
            
        mv {params.outdir:q}/{wildcards.sample:q}_R1_trimmed.fq.gz {output.trimmed:q}
        """


rule trim_galore_pe:
    """
    Run Trim Galore on Paired-End reads in paired mode.
    """
    input:
        fastq_read1="fastq/{sample}_R1.fastq.gz",
        fastq_read2="fastq/{sample}_R2.fastq.gz",
    output:
        fastq_read1="fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz",
        fastq_read2="fastq/fastq_trimmed/trimgalore/{sample}_R2.fastq.gz",
        report1="fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz_trimming_report.txt",
        report2="fastq/fastq_trimmed/trimgalore/{sample}_R2.fastq.gz_trimming_report.txt",
    threads: config["CORES"]["trimgalore"]
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        # Trim Galore estimates up to 3 * -j + 3 CPUs; -j must remain at least 1.
        cores=lambda wildcards, threads: max(1, (threads - 3) // 3),
        outdir="fastq/fastq_trimmed/trimgalore",
    log:
        "fastq/fastq_trimmed/trimgalore/{sample}_trim_galore.log",
    conda:
        "transcript_env.yaml"
    shell:
        """
        mkdir -p {params.outdir:q}
        
        trim_galore -j {params.cores} \
            -o {params.outdir:q} \
            --basename {wildcards.sample:q} \
            {params.trim_galore_params} \
            --paired {input.fastq_read1:q} {input.fastq_read2:q} > {log:q} 2>&1
            
            mv {params.outdir:q}/{wildcards.sample:q}_val_1.fq.gz {output.fastq_read1:q}
            mv {params.outdir:q}/{wildcards.sample:q}_val_2.fq.gz {output.fastq_read2:q}

        """


##################################
# Rule for linking fastq_trimmed #
##################################

if LAYOUT == "SINGLE":

    rule link_trimmed_se:
        input:
            r1=_trimmer_out("{sample}", "R1"),
        output:
            r1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        log:
            "fastq/fastq_trimmed/{sample}_linking.log",
        shell:
            "ln -srf {input.r1:q} {output.r1:q} > {log:q} 2>&1"

elif LAYOUT == "PAIRED":

    rule link_trimmed_pe:
        input:
            r1=_trimmer_out("{sample}", "R1"),
            r2=_trimmer_out("{sample}", "R2"),
        output:
            r1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
            r2="fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        log:
            "fastq/fastq_trimmed/{sample}_linking.log",
        shell:
            """
            ln -srf {input.r1:q} {output.r1:q} > {log:q} 2>&1
            ln -srf {input.r2:q} {output.r2:q} >> {log:q} 2>&1
            """
