if config["LAYOUT"] == "PAIRED":
    ruleorder: align_pe_bwa > align_se_bwa
elif config["LAYOUT"] == "SINGLE":
    ruleorder: align_se_bwa > align_pe_bwa


rule align_pe_bwa:
    input:
        trimmed_fastq1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        trimmed_fastq2="fastq/fastq_trimmed/{sample}_R2.fastq.gz"
    output:
        bam="aligned_bwa/{sample}.aligned.bam"
    params:
        reference="/home/molinerislab/NeriMetagenome/Reference_genome/bwa/GRCh38.p14.genome.fa", 
        threads=4  
    log:
        "aligned_bwa/logs/{sample}_bwa.log"
    shell:
        """
        T=$(mktemp -d) &&  \
        bwa mem -t {params.threads} {params.reference} {input.trimmed_fastq1} {input.trimmed_fastq2} 2> {log} | \
        samtools view -b | \
        sambamba sort --tmpdir=$T -t 4 --memory-limit 4GB -o {output.bam} /dev/stdin ; \
        rm -rf $T
        """


rule align_se_bwa:
    input:
        trimmed_fastq="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
    output:
        bam="aligned_bwa/{sample}.aligned.bam"
    params:
        reference="/home/molinerislab/NeriMetagenome/Reference_genome/bwa/GRCh38.p14.genome.fa", 
        threads=4  
    log:
        "aligned_bwa/logs/{sample}_bwa.log"
    shell:
        """
        T=$(mktemp -d) &&  \
        bwa mem -t {params.threads} {params.reference} {input.trimmed_fastq} 2> {log} | \
        samtools view -b | \
        sambamba sort --tmpdir=$T -t 4 --memory-limit 4GB -o {output.bam} /dev/stdin ; \
        rm -rf $T
        """

