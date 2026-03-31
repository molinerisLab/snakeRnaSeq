import os


#############
# Ruleorder #
#############
configfile: "/home/aleone/snakeRnaSeq/local/config/config_v1.yaml"


include: "/home/aleone/snakeRnaSeq/local/config/config_v1.sk"


if config["LAYOUT"] == "SINGLE":

    ruleorder: generate_unmapped_single > generate_unmapped_R1
    ruleorder: generate_unmapped_single > generate_unmapped_R2
    ruleorder: star_twopass_basic_se > star_twopass_basic_pe

elif config["LAYOUT"] == "PAIRED":

    ruleorder: generate_unmapped_R1 > generate_unmapped_single
    ruleorder: generate_unmapped_R2 > generate_unmapped_single
    ruleorder: star_twopass_basic_pe > star_twopass_basic_se


##############
# STAR RULES #
##############


rule star_align_se:
    input:
        fq1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        aln="Results/star/{sample}.bam",
        log="Results/star/{sample}.Log.out",
        sj="Results/star/{sample}.SJ.out.tab",
        unmapped=(
            "Results/star/unmapped/{sample}_unmapped_R1.fastq.gz"
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
        log_final="Results/star/{sample}.Log.final.out",
    log:
        "Results/star/{sample}.log",
    threads: 4
    conda:
        "transcript_env.yaml"
    params:
        multiscorerange=config["STAR"]["MULTIMAP_SCORE_RANGE"],
        outfiltermismatch=config["STAR"]["OUT_FILTER_MISMATCH_NMAX"],
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        outfiltermismatchnover=config["STAR"]["OUT_FILTER_MISMATCH_NOVER_LMAX"],
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        fq_join=lambda wc, input: " ".join(input.fq),
        tmpdir=lambda wc, output: os.path.dirname(output.aln),
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeLoad NoSharedMemory \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq_join} \
            --readFilesCommand {params.read_cmd} \
            --sjdbOverhang {params.sjdbOver} \
            --outFileNamePrefix {params.tmpdir}/ \
            --quantMode {params.quant_mode} \
            --outSAMtype {params.out_samtype} \
            --outSAMstrandField intronMotif \
            --outFilterMultimapScoreRange {params.multiscorerange} \
            --outFilterMismatchNmax {params.outfiltermismatch} \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --outFilterMismatchNoverLmax {params.outfiltermismatchnover} \
            --outSAMattributes All
        """


rule star_align_pe:
    input:
        fq1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        fq2="fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        aln="Results/star/{sample}.bam",
        log="Results/star/{sample}.Log.out",
        sj="Results/star/{sample}.SJ.out.tab",
        unmapped=(
            [
                "Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
                "Results/star/unmapped/{sample}_unmapped_R2.fastq.gz",
            ]
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
        log_final="Results/star/{sample}.Log.final.out",
    log:
        "Results/star/{sample}.log",
    threads: 4
    conda:
        "transcript_env.yaml"
    params:
        multiscorerange=config["STAR"]["MULTIMAP_SCORE_RANGE"],
        outfiltermismatch=config["STAR"]["OUT_FILTER_MISMATCH_NMAX"],
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        outfiltermismatchnover=config["STAR"]["OUT_FILTER_MISMATCH_NOVER_LMAX"],
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        tmpdir=lambda wc, output: os.path.dirname(output.aln),
        save_unmapped=config["STAR"]["SAVE_UNMAPPED"],
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeLoad NoSharedMemory \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --sjdbOverhang {params.sjdbOver} \
            --outFileNamePrefix {params.tmpdir}/ \
            --quantMode {params.quant_mode} \
            --outSAMtype {params.out_samtype} \
            --outSAMstrandField intronMotif \
            --outFilterMultimapScoreRange {params.multiscorerange} \
            --outFilterMismatchNmax {params.outfiltermismatch} \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --outFilterMismatchNoverLmax {params.outfiltermismatchnover} \
            --outSAMattributes All \
            $([ "{params.save_unmapped}" = "FASTQ" ] && echo "--outReadsUnmapped Fastx --outSAMunmapped Within" || echo "")
        """


rule link_unmapped:
    input:
        "Results/star/unmapped/{sample}_unmapped_R{mate}.fastq.gz",
    output:
        "fastq/unmapped/{sample}_unmapped_R{mate}.fastq.gz",
    log:
        "fastq/unmapped/{sample}_unmapped_R{mate}.log",
    conda:
        "transcript_env.yaml"
    shell:
        "ln {input} {output}"


rule generate_unmapped_single:
    input:
        lambda wildcards: (
            f"Results/star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        ),
    output:
        "fastq/unmapped/{sample}_unmapped.fastq.gz",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/unmapped/{sample}_unmapped.log",
    shell:
        """
        samtools view -f 4 {input} | awk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
        """


rule generate_unmapped_R1:
    input:
        lambda wildcards: (
            f"Results/star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        ),
    output:
        "fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/unmapped/{sample}_unmapped_R1.log",
    shell:
        """
        samtools view -f 76 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """


# 76=4+8+64 = read unmapped AND mate unmapped AND first in pair, i.e., discard reads that are unmapped but that have mate mapped


rule generate_unmapped_R2:
    input:
        lambda wildcards: (
            f"Results/star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        ),
    output:
        "fastq/unmapped/{sample}_unmapped_R2.fastq.gz",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/unmapped/{sample}_unmapped_R2.log",
    shell:
        """
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
                f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz",
            ]
            if config["LAYOUT"] == "PAIRED"
            else [f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz"]
        ),
        idx=STAR_GENOME_DIR,
    output:
        sj="Results/pass1/{sample}/{sample}_SJ.out.tab",
        log="Results/pass1/{sample}/{sample}_Log.out",
        log_final="Results/pass1/{sample}/{sample}_Log.final.out",
        log_progress="Results/pass1/{sample}/{sample}_Log.progress.out",
    threads: 16
    conda:
        "transcript_env.yaml"
    log:
        "Results/pass1/{sample}/{sample}_star_pass1.log",
    params:
        tmpdir=lambda wc, output: os.path.dirname(output.sj),
        read_cmd=config["STAR"]["readFilesCommand"],
        limitSjdb=config["STAR"]["limitSjdbInsertNsj"],
        fq_join=lambda wc, input: " ".join(input.fq),
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
            --
        """


rule merge_and_filter_sj:
    input:
        expand("Results/pass1/{sample}/{sample}_SJ.out.tab", sample=SAMPLES),
    output:
        "Results/pass1/merged_filtered_SJ.out.tab",
    conda:
        "transcript_env.yaml"
    log:
        "Results/pass1/merged_filtered_SJ.log",
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
                f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz",
            ]
            if config["LAYOUT"] == "PAIRED"
            else [f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz"]
        ),
        idx=STAR_GENOME_DIR,
        sj=lambda wc: f"Results/pass1/merged_filtered_SJ.out.tab",
    output:
        bam="Results/pass2/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        gene_counts="Results/pass2/{sample}/ReadsPerGene.out.tab",
    threads: 8
    conda:
        "transcript_env.yaml"
    log:
        "Results/pass2/{sample}/{sample}_star_pass2.log",
    params:
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        limitSjdb=config["STAR"]["limitSjdbInsertNsj"],
        fq_join=lambda wc, input: " ".join(input.fq),
        tmpdir=lambda wc, output: os.path.dirname(output.bam),
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
        fq=lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        bam="star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj="star_2pass/{sample}/SJ.out.tab",
        gene_counts="star_2pass/{sample}/ReadsPerGene.out.tab",
        log="star_2pass/{sample}/Log.out",
        log_final="star_2pass/{sample}/Log.final.out",
    threads: 8
    conda:
        "transcript_env.yaml"
    log:
        "star_2pass/{sample}/star_2pass_basic_se.log",
    params:
        read_cmd=config["STAR"]["readFilesCommand"],
        out_dir=lambda wc, output: os.path.dirname(output.bam),
        gtf=GENCODE_ANNOTATION_GTF,
        sjdbOverhang=config["STAR"]["sjdbOverhang"],
    shell:
        """
        mkdir -p {params.out_dir}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_dir}/ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts TranscriptomeSAM
        """


rule all_s2p_basic:
    input:
        expand("star_2pass/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),


rule star_twopass_basic_pe:
    """
    STAR 2-pass Basic mode (paired-end).
    """
    input:
        fq1=lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        fq2=lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R2.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        bam="star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj="star_2pass/{sample}/SJ.out.tab",
        gene_counts="star_2pass/{sample}/ReadsPerGene.out.tab",
        log="star_2pass/{sample}/Log.out",
        log_final="star_2pass/{sample}/Log.final.out",
    threads: 8
    conda:
        "transcript_env.yaml"
    log:
        "star_2pass/{sample}/star_2pass_basic_pe.log",
    params:
        read_cmd=config["STAR"]["readFilesCommand"],
        out_dir=lambda wc, output: os.path.dirname(output.bam),
        twopass1readsN=config["STAR"].get("twopass1readsN", -1),
        gtf=GENCODE_ANNOTATION_GTF,
        sjdbOverhang=config["STAR"]["sjdbOverhang"],
    shell:
        """
        mkdir -p {params.out_dir}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --twopassMode Basic \
            --twopass1readsN {params.twopass1readsN} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_dir}/ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts TranscriptomeSAM
        """
