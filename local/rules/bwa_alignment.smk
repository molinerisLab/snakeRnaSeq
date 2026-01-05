rule align_pe_bwa:
    input:
        trimmed_fastq1="{output_dir}/fastq_trimmed/{sample}_R1.fastq.gz",
        trimmed_fastq2="{output_dir}/fastq_trimmed/{sample}_R2.fastq.gz"
    output:
        bam="{output_dir}/aligned_bwa/{sample}.aligned.bam"
    params:
        reference="/home/molinerislab/NeriMetagenome/Reference_genome/bwa/GRCh38.p14.genome.fa", 
        threads=20  
    log:
        "{output_dir}/aligned_bwa/logs/{sample}_bwa.log"
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
        trimmed_fastq="{output_dir}/fastq_trimmed/{sample}.fastq.gz",
    output:
        bam="{output_dir}/aligned_bwa/{sample}.aligned.bam"
    params:
        reference="/home/molinerislab/NeriMetagenome/Reference_genome/bwa/GRCh38.p14.genome.fa", 
        threads=20  
    log:
        "{output_dir}/aligned_bwa/logs/{sample}_bwa.log"
    shell:
        """
        T=$(mktemp -d) &&  \
        bwa mem -t {params.threads} {params.reference} {input.trimmed_fastq} 2> {log} | \
        samtools view -b | \
        sambamba sort --tmpdir=$T -t 4 --memory-limit 4GB -o {output.bam} /dev/stdin ; \
        rm -rf $T
        """

