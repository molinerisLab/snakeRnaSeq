reference_fasta_path     = f"{config['REFERENCE_DIR']}/{config['GENOME_ASSEMBLY']}.genome.fa"
annotation_gtf_path      = f"{config['REFERENCE_DIR']}/primary_assembly.annotation.gtf"
transcriptome_fasta_path = f"{config['REFERENCE_DIR']}/gencode.v{config['GENCODE_RELEASE']}.transcripts.fa"


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
        idx = config['STAR']['INDEX']['GRCh'],
    output:
        aln="star/{sample}.bam",
        log="star/{sample}.Log.out",
        sj="star/{sample}.SJ.out.tab",
    threads: config["CORES"]
    params:
        genome_dir = input.idx,
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
        fq1 = "fastq/{sample}_R1.fastq.gz",
        fq2 = "fastq/{sample}_R2.fastq.gz",
        idx = config['STAR']['INDEX']['GRCh'],
    output:
        aln = "star/{sample}.bam",
        log = "star/{sample}.Log.out",
        sj ="star/{sample}.SJ.out.tab",
        # Uncomment the next line if you want to handle unmapped reads
        # unmapped=["star/unmapped/{sample}_R1.fastq.gz", "star/unmapped/{sample}_R2.fastq.gz"],
        #unmapped read filtered after, sice by default STAR report as unmapped partially mapped (i.e. mapped only one mate of a paired end read)
        log_final="star/{sample}.Log.final.out"
    log:
        "star/{sample}.log",
    params:
        extra=lambda wildcards: f"--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --chimOutType WithinBAM {config['STAR']['OPTIONS']}",
    threads: 16,
    wrapper:
        "v3.3.6/bio/star/align"

#rule star_pe_multi:
#    input:
#        fq1 = "fastq/{sample}_R1.fastq.gz",
#        fq2 = "fastq/{sample}_R2.fastq.gz",
#        idx = config['STAR']['INDEX']['GRCh'],
#    output:
#        aln = "star/{sample}.bam",
#        log = "star/{sample}.Log.out",
#        sj = "star/{sample}.SJ.out.tab",
        # Uncomment the next line if you want to handle unmapped reads
        # unmapped=["star/unmapped/{sample}_R1.fastq.gz", "star/unmapped/{sample}_R2.fastq.gz"],
        # unmapped read filtered after, sice by default STAR report as unmapped partially mapped (i.e. mapped only one mate of a paired end read)
#    log:
#        "star/{sample}.log",
#    params:
#        extra=lambda wildcards: f"--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --chimOutType WithinBAM {config['STAR']['OPTIONS']}",
#    threads: config["CORES"],
#    wrapper:
#        "v3.3.6/bio/star/align"

rule link_unmapped:
    input: "star/{sample}_unmapped_R{mate}.fastq.gz"
    output: "fastq_unmapped/{sample}_R{mate}.fastq.gz"
    shell: "ln {input} {output}"


rule generate_unmapped_R1:
    input:
        "star/{sample}.bam"
    output:
        "{sample}_unmapped_R1.fastq.gz"
    shell:"""
        samtools view -f 76 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 76=4+8+64 = read unmapped AND mate unmapped AND first in pair, i.e., discard reads that are unmapped but that have mate mapped

rule generate_unmapped_R2:
    input:
        "star/{sample}.bam"
    output:
        "{sample}_unmapped_R2.fastq.gz"
    shell: """
        samtools view -f 140 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 140=4+8+128 = read unmapped AND mate unmapped AND second in pair, i.e., discard reads that are unmapped but that have mate mapped



# ==============================================
# STAR DOUBLE PASS ALIGNMENT (if needed)
# ==============================================

rule star_align_first_pass:
    input:
        fq=lambda wc: (
            [
                f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
                f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz"
            ]
            if config["LAYOUT"] == "PAIRED"
            else [
                f"fastq/fastq_trimmed/{wc.sample}.fastq.gz"
            ]
        ),
        idx=lambda wc: config['STAR_GENOMEDIR']['GRCh']
    output:
        sj           = "Results/pass1_{genome}/{sample}/SJ.out.tab",
        log          = "Results/pass1_{genome}/{sample}/Log.out",
        log_final    = "Results/pass1_{genome}/{sample}/Log.final.out",
        log_progress = "Results/pass1_{genome}/{sample}/Log.progress.out"
    threads: 16
    params:
        tmpdir     = "Results/pass1_{genome}/{sample}",
        read_cmd   = config["STAR"]["readFilesCommand"],
        limitSjdb  = config["STAR"]["limitSjdbInsertNsj"],
        fq_join    = lambda wc, input: " ".join(input.fq)
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq_join} \
            --readFilesCommand {params.read_cmd} \
            --limitSjdbInsertNsj {params.limitSjdb} \
            --outFileNamePrefix {params.tmpdir}/ \
            --outSAMtype None
        """


        
rule merge_and_filter_sj:
    input:
        expand("Results/pass1_{{genome}}/{sample}/SJ.out.tab",
               sample= SAMPLES)
    output:
        "Results/pass1_{genome}/merged_filtered_SJ.out.tab"
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} \
          | awk '$5>=1 && $5<=6 && $6==0 && $7>2' \
          | sort -u \
          > {output}
        """



rule star_second_pass:
    input:
        fq=lambda wc: (
            [
                f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
                f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz"
            ]
            if config["LAYOUT"] == "PAIRED"
            else [
                f"fastq/fastq_trimmed/{wc.sample}.fastq.gz"
            ]
        ),
        idx=lambda wc: STAR_GENOMEDIR[wc.genome],
        sj=lambda wc: f"Results/pass1_{wc.genome}/merged_filtered_SJ.out.tab"
    output:
        bam         = "Results/pass2_{genome}/{sample}/Aligned.sortedByCoord.out.bam",
        gene_counts = "Results/pass2_{genome}/{sample}/ReadsPerGene.out.tab"
    threads: 16
    params:
        out_samtype = config["STAR"]["outSAMtype"],
        quant_mode  = config["STAR"]["quantMode"],
        sjdbOver    = config["STAR"]["sjdbOverhang"],
        read_cmd    = config["STAR"]["readFilesCommand"],
        limitSjdb   = config["STAR"]["limitSjdbInsertNsj"],
        fq_join     = lambda wc, input: " ".join(input.fq),
        tmpdir      = "Results/pass2_{genome}/{sample}"
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq_join} \
            --readFilesCommand {params.read_cmd} \
            --sjdbFileChrStartEnd {input.sj} \
            --sjdbOverhang {params.sjdbOver} \
            --limitSjdbInsertNsj {params.limitSjdb} \
            --outFileNamePrefix {params.tmpdir}/ \
            --outSAMtype {params.out_samtype} \
            --quantMode {params.quant_mode}
        """


rule star_twopass_basic_se:
    """
    STAR 2-pass Basic mode (single-end).
    STAR automatically:
    1. Maps reads (1st pass)
    2. Extracts junctions
    3. Updates splice junction DB
    4. Re-maps reads using updated junctions (2nd pass)
    """
    input:
        fq = lambda wc: f"fastq/fastq_trimmed/{wc.sample}.fastq.gz",
        idx = lambda wc: STAR_GENOMEDIR[GENOME_KEY]
    output:
        bam        = "star_2pass_{genome}/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass_{genome}/{sample}/SJ.out.tab",
        gene_counts= "star_2pass_{genome}/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass_{genome}/{sample}/Log.out",
        log_final  = "star_2pass_{genome}/{sample}/Log.final.out"
    threads: 16
    conda: "transcript_env.yaml"
    params:
        genome_dir   = lambda wc: STAR_GENOMEDIR[wc.genome_key],
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass_{genome}/{sample}/",
        gtf          = annotation_gtf_path,
        sjdbOverhang = config["STAR"]["sjdbOverhang"]
    shell:
        """
        mkdir -p star_2pass_{wildcards.genome}/{wildcards.sample}

        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts TranscriptomeSAM
        """


rule star_twopass_basic_pe:
    """
    STAR 2-pass Basic mode (paired-end).
    """
    input:
        fq1 = lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        fq2 = lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz",
        idx = lambda wc: STAR_GENOMEDIR[wc.genome]
    output:
        bam        = "star_2pass_{genome}/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass_{genome}/{sample}/SJ.out.tab",
        gene_counts= "star_2pass_{genome}/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass_{genome}/{sample}/Log.out",
        log_final  = "star_2pass_{genome}/{sample}/Log.final.out"
    threads: 16
    conda: "transcript_env.yaml"
    params:
        genome_dir   = lambda wc: STAR_GENOMEDIR[wc.genome],
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass_{genome}/{sample}/",
        twopass1readsN = config["STAR"].get("twopass1readsN", -1),
        gtf          = annotation_gtf_path,
        sjdbOverhang = config["STAR"]["sjdbOverhang"]
    shell:
        """
        mkdir -p star_2pass_{wildcards.genome}/{wildcards.sample}

        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --twopassMode Basic \
            --twopass1readsN {params.twopass1readsN} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts TranscriptomeSAM
        """
