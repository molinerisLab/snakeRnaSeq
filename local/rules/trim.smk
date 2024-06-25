# ----- #
# Fastp #
# ----- #

rule fastp_pe:
    input:
        sample=["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"]
    output:
        trimmed=["fastq_trimmed/{sample}_R1.fastq.gz", "fastq_trimmed/{sample}_R2.fastq.gz"],
        #unpaired1="trimmed/pe/{sample}.u1.fastq",
        #unpaired2="trimmed/pe/{sample}.u2.fastq",
        #merged="trimmed/pe/{sample}.merged.fastq",
        #failed="trimmed/pe/{sample}.failed.fastq",
        html="fastq_trimmed/{sample}.html",
        json="fastq_trimmed/{sample}.json"
    log:
        "fastq_trimmed/{sample}.log"
    params:
        adapters="--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        #extra="--merge"
        extra=""
    threads: 4
    wrapper:
        "v3.3.6/bio/fastp"


# ------------ #
# TRIM GALORE! #
# ------------ #

rule trim_galore_se:
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