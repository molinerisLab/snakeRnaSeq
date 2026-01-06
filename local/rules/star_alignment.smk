## TODO: the genome load can be performed w/ the wrappers by using a dummy file 

def choose_fastq_according_to_genome(wildcards, mate):
    sample = wildcards['sample']  # Definisce 'sample' usando 'wildcards'
    if config['GENCODE']['ASSEMBLY'] == "GRCh":
        return f"fastq/{sample}_R{mate}_001.fastq.gz"
    # elif wildcards['genome'] == "CHM":
    #     return f"star_GRCh/{sample}_unmapped_R{mate}.fastq.gz"
    else:
        raise Exception(f"Genome not valid: {config['GENOME']}")


rule star_align_se:
    input:
        fq1 = "fastq/{sample}_R1.fastq.gz",
    output:
        aln="star/{sample}.bam",
        log="star/{sample}.Log.out",
        sj="star/{sample}.SJ.out.tab",
    threads: config["CORES"]
    params:
        genome_dir = STAR_GENOME_DIR,
        tmpdir="star/{sample}",
        out_sam_type = config["STAR"]["OUT_SAM_TYPE"],
        out_filter_multimap_nmax = config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        out_filter_multimap_score_range = config["STAR"]["MULTIMAP_SCORE_RANGE"],
        filter_mismatch = config["STAR"]["FILTER_MISMATCH"],
        additional_output = config["STAR"]["ADDITIONAL_OUTPUT"]
    shell:
        """
        mkdir -p {params.tmpdir};
        STAR \
            --genomeDir {params.genome_dir} \
            --genomeLoad LoadAndKeep \
            --runThreadN {threads} \
            --readFilesIn {input.fq1} \
            --readFilesCommand zcat \
            --outFileNamePrefix star/{wildcards.sample}. \
            --outTmpDir {params.tmpdir}/STARtmp \
            --outSAMtype {params.out_sam_type} \
            --limitBAMsortRAM 10000000000 \
            --outSAMunmapped Within \
            --outFilterMultimapNmax {params.out_filter_multimap_nmax} \
            --outFilterMultimapScoreRange {params.out_filter_multimap_score_range} \
            {params.filter_mismatch} \
            {params.additional_output}
        """

rule star_align_pe:
    input:
        fq1=lambda wildcards: f"{config['output_dir']}/fastq_trimmed/{wildcards.sample}_R1.fastq.gz",
        fq2=lambda wildcards: f"{config['output_dir']}/fastq_trimmed/{wildcards.sample}_R2.fastq.gz",
        idx=lambda wildcards: config['star_index'][wildcards['genome']],
    output:
        aln="{output_dir}/star_{genome}/{sample}.bam",
        log="{output_dir}/star_{genome}/{sample}.Log.out",
        sj="{output_dir}/star_{genome}/{sample}.SJ.out.tab",
        # Uncomment the next line if you want to handle unmapped reads
        # unmapped=["star_{genome}/unmapped/{sample}_R1.fastq.gz", "star_{genome}/unmapped/{sample}_R2.fastq.gz"],
        #unmapped read filtered after, sice by default STAR report as unmapped partially mapped (i.e. mapped only one mate of a paired end read)
        log_final="{output_dir}/star_{genome}/{sample}.Log.final.out"
    log:
        "{output_dir}/star_{genome}/{sample}.log",
    params:
        extra=lambda wildcards: f"--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --chimOutType WithinBAM {config['star_options']}",
    threads: 16,
    wrapper:
        "v3.3.6/bio/star/align"

#gestire single vs pair ends
rule star_pe_multi:
    input:
        fq1=config['FASTQ_FILTERING']+"/{sample}_R1.fastq.gz",
        fq2=config['FASTQ_FILTERING']+"/{sample}_R2.fastq.gz",
        idx=lambda wildcards: config['STAR']['INDEX']['GRCh'],
    output:
        aln="star/{sample}.bam",
        log="star/{sample}.Log.out",
        sj="star/{sample}.SJ.out.tab",
        # Uncomment the next line if you want to handle unmapped reads
        # unmapped=["star/unmapped/{sample}_R1.fastq.gz", "star/unmapped/{sample}_R2.fastq.gz"],
        #unmapped read filtered after, sice by default STAR report as unmapped partially mapped (i.e. mapped only one mate of a paired end read)
    log:
        "star/{sample}.log",
    params:
        extra=lambda wildcards: f"--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --chimOutType WithinBAM {config['STAR']['OPTIONS']}",
    threads: config["CORES"],
    wrapper:
        "v3.3.6/bio/star/align"

rule linl_unmapped:
    input: "star_CHM/{sample}_unmapped_R{mate}.fastq.gz"
    output: "fastq_unmapped/{sample}_R{mate}.fastq.gz"
    shell: "ln {input} {output}"


rule generate_unmapped_R1:
    input:
        "{sample}.bam"
    output:
        "{sample}_unmapped_R1.fastq.gz"
    shell:"""
        samtools view -f 76 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 76=4+8+64 = read unmapped AND mate unmapped AND first in pair, i.e., discard reads that are unmapped but that have mate mapped

rule generate_unmapped_R2:
    input:
        "{sample}.bam"
    output:
        "{sample}_unmapped_R2.fastq.gz"
    shell: """
        samtools view -f 140 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 140=4+8+128 = read unmapped AND mate unmapped AND second in pair, i.e., discard reads that are unmapped but that have mate mapped

# rule :
#     input: 
#     output: 
#     shell: 

