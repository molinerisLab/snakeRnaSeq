
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
#TODO: conda env per each rule, missing in fastp

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
        "fastq/fastq_trimmed/trimgalore/{sample}_R1.fastq.gz",
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores= config["CORES"]["trimgalore"],
        outdir="fastq/fastq_trimmed/trimgalore",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/fastq_trimmed/trimgalore/{sample}_trim_galore.log",
    shell:
        """
        mkdir -p {params.outdir}
        
        trim_galore -j {params.cores} \
            --basename {wildcards.sample}_R1 \
            -o {params.outdir} \
            {params.trim_galore_params} \
            {input} > {log} 2>&1
            
        mv {params.outdir}/{wildcards.sample}_R1_trimmed.fq.gz {output}
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
    params:
        trim_galore_params=config["TRIM_GALORE"]["PARAM"],
        cores=config["CORES"]["trimgalore"],
        outdir="fastq/fastq_trimmed/trimgalore",
    log:
        "fastq/fastq_trimmed/trimgalore/{sample}_trim_galore.log",
    conda:
        "transcript_env.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        trim_galore -j {params.cores} \
            -o {params.outdir} \
            --basename {wildcards.sample} \
            {params.trim_galore_params} \
            --paired {input.fastq_read1} {input.fastq_read2} > {log} 2>&1
            
            mv {params.outdir}/{wildcards.sample}_val_1.fq.gz {output.fastq_read1}
            mv {params.outdir}/{wildcards.sample}_val_2.fq.gz {output.fastq_read2}

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
        conda:
            "transcript_env.yaml"
        shell:
            "ln -srf {input.r1} {output.r1} > {log} 2>&1"

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
        conda:
            "transcript_env.yaml"
        shell:
            """
            ln -srf {input.r1} {output.r1} > {log} 2>&1
            ln -srf {input.r2} {output.r2} >> {log} 2>&1
            """
