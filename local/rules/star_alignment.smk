
if config["LAYOUT"] == "SINGLE":
     ruleorder: generate_unmapped_single > generate_unmapped_R1
     ruleorder: generate_unmapped_single > generate_unmapped_R2
elif config["LAYOUT"] == "PAIRED":
     ruleorder: generate_unmapped_R1 > generate_unmapped_single
     ruleorder: generate_unmapped_R2 > generate_unmapped_single

rule star_align_se:
    input:
        fq1 = "fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        idx = config['STAR']['INDEX'][config['GENCODE']['ASSEMBLY']]
    output:
        aln = "star/{sample}.bam",
        log = "star/{sample}.Log.out",
        sj  = "star/{sample}.SJ.out.tab",
        unmapped = (
            "star/unmapped/{sample}_unmapped_R1.fastq.gz"
            if config['STAR']['SAVE_UNMAPPED'] == "FASTQ" 
            else []
        ),
        log_final = "star/{sample}.Log.final.out"
    log:
        "star/{sample}.log"
    threads: config["CORES"]
    conda: "transcript_env.yaml"
    params:
        extra = lambda wildcards: (
            f"--outSAMtype {config['STAR']['OUT_SAM_TYPE']} "
            f"--limitBAMsortRAM 10000000000 "
            f"--genomeLoad LoadAndKeep "
            f"--chimOutType WithinBAM "   
            f"--outFilterMultimapNmax {config['STAR']['OUT_FILTER_MULTIMAP_NMAX']} "
            f"--outFilterMultimapScoreRange {config['STAR']['MULTIMAP_SCORE_RANGE']} "
            f"--outFilterMismatchNoverReadLmax {config['STAR']['OUT_FILTER_MISMATCH_NOVER_LMAX']} "
            f"--alignSJoverhangMin {config['STAR']['ALIGN_SJ_OVERHANG_MIN']} "
            f"--alignSJDBoverhangMin {config['STAR']['ALIGN_SJDB_OVERHANG_MIN']} "
            f"--alignIntronMin {config['STAR']['ALIGN_INTRON_MIN']} "
            f"--alignIntronMax {config['STAR']['ALIGN_INTRON_MAX']} "
            f"--alignMatesGapMax {config['STAR']['ALIGN_MATES_GAP_MAX']} "
            f"{config['STAR']['ADDITIONAL_OUTPUT']}"
            f"{'--outReadsUnmapped Fastx --outSAMunmapped None' if config['STAR']['SAVE_UNMAPPED'] == 'FASTQ' else '--outSAMunmapped Within'}"
        )
    wrapper:
        "v3.3.6/bio/star/align"

rule star_align_pe:
    input:
        fq1 = "fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        fq2 = "fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        idx = config['STAR']['INDEX'][config['GENCODE']['ASSEMBLY']]
    output:
        aln = "star/{sample}.bam",
        log = "star/{sample}.Log.out",
        sj  = "star/{sample}.SJ.out.tab",
        unmapped = (
            ["star/unmapped/{sample}_unmapped_R1.fastq.gz", 
             "star/unmapped/{sample}_unmapped_R2.fastq.gz"]
            if config['STAR']['SAVE_UNMAPPED'] == "FASTQ" 
            else []
        ),
        log_final = "star/{sample}.Log.final.out"
    log:
        "star/{sample}.log"
    conda: 
        "transcript_env.yaml"
    threads: 16
    params:
        extra = lambda wildcards: (
            f"--outSAMtype {config['STAR']['OUT_SAM_TYPE']} "
            f"--chimOutType WithinBAM "
            f"--outFilterMultimapNmax {config['STAR']['OUT_FILTER_MULTIMAP_NMAX']} "
            f"--outFilterMismatchNoverReadLmax {config['STAR']['OUT_FILTER_MISMATCH_NOVER_LMAX']} "
            f"--alignSJoverhangMin {config['STAR']['ALIGN_SJ_OVERHANG_MIN']} "
            f"--alignSJDBoverhangMin {config['STAR']['ALIGN_SJDB_OVERHANG_MIN']} "
            f"--alignIntronMin {config['STAR']['ALIGN_INTRON_MIN']} "
            f"--alignIntronMax {config['STAR']['ALIGN_INTRON_MAX']} "
            f"--alignMatesGapMax {config['STAR']['ALIGN_MATES_GAP_MAX']} "
            f"{config['STAR']['ADDITIONAL_OUTPUT']}"
            f"{'--outReadsUnmapped Fastx --outSAMunmapped None' if config['STAR']['SAVE_UNMAPPED'] == 'FASTQ' else '--outSAMunmapped Within'}"
        )
    wrapper:
        "v3.3.6/bio/star/align"


rule link_unmapped:
    input:
        "star/{sample}_unmapped_R{mate}.fastq.gz"
    output:
        "fastq/unmapped/{sample}_unmapped_R{mate}.fastq.gz"
    shell:
        "ln {input} {output}"

rule generate_unmapped_single:
    input:
        lambda wildcards: (
            f"star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        )
    output:
        "fastq/unmapped/{sample}_unmapped.fastq.gz"
    conda: "transcript_env.yaml"
    shell:
        """
        samtools view -f 4 {input} | awk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
        """


rule generate_unmapped_R1:
    input:
       lambda wildcards: (
            f"star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        )
    output:
        "fastq/unmapped/{sample}_unmapped_R1.fastq.gz"
    conda: "transcript_env.yaml"
    shell:"""
        samtools view -f 76 {input} | bawk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
    """
    # 76=4+8+64 = read unmapped AND mate unmapped AND first in pair, i.e., discard reads that are unmapped but that have mate mapped

rule generate_unmapped_R2:
    input:
        lambda wildcards: (
            f"star/{wildcards.sample}.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        )
    output:
        "fastq/unmapped/{sample}_unmapped_R2.fastq.gz"
    conda: "transcript_env.yaml"
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
                f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz"
            ]
        ),
        idx=lambda wc: config['STAR']['INDEX']['GRCh'],
    output:
        sj           = "Results/pass1/{sample}/SJ.out.tab",
        log          = "Results/pass1/{sample}/Log.out",
        log_final    = "Results/pass1/{sample}/Log.final.out",
        log_progress = "Results/pass1/{sample}/Log.progress.out"
    threads: 16
    conda: "transcript_env.yaml"
    params:
        tmpdir     = "Results/pass1/{sample}",
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
        expand("Results/pass1/{sample}/SJ.out.tab",
               sample= SAMPLES)
    output:
        "Results/pass1/merged_filtered_SJ.out.tab"
    conda: "transcript_env.yaml"
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
                f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz"
            ]
        ),
        idx=config['STAR']['INDEX']['GRCh'],
        sj=lambda wc: f"Results/pass1/merged_filtered_SJ.out.tab"
    output:
        bam         = "Results/pass2/{sample}/Aligned.sortedByCoord.out.bam",
        gene_counts = "Results/pass2/{sample}/ReadsPerGene.out.tab"
    threads: 8
    conda: "transcript_env.yaml"
    params:
        out_samtype = config["STAR"]["OUT_SAM_TYPE"],
        quant_mode  = config["STAR"]["quantMode"],
        sjdbOver    = config["STAR"]["sjdbOverhang"],
        read_cmd    = config["STAR"]["readFilesCommand"],
        limitSjdb   = config["STAR"]["limitSjdbInsertNsj"],
        fq_join     = lambda wc, input: " ".join(input.fq),
        tmpdir      = "Results/pass2/{sample}"
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
        fq = lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        idx = config['STAR']['INDEX']['GRCh'],
    output:
        bam        = "star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass/{sample}/SJ.out.tab",
        gene_counts= "star_2pass/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass/{sample}/Log.out",
        log_final  = "star_2pass/{sample}/Log.final.out"
    threads: 8
    conda: "transcript_env.yaml"
    params:
        genome_dir   = config['STAR']['INDEX']['GRCh'],
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass/{sample}/",
        gtf          = annotation_gtf_path,
        sjdbOverhang = config["STAR"]["sjdbOverhang"]
    shell:
        """
        mkdir -p star_2pass/{wildcards.sample}

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
        idx = config['STAR']['INDEX']['GRCh'],
    output:
        bam        = "star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass/{sample}/SJ.out.tab",
        gene_counts= "star_2pass/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass/{sample}/Log.out",
        log_final  = "star_2pass/{sample}/Log.final.out"
    threads: 8
    conda: "transcript_env.yaml"
    params:
        genome_dir   = config['STAR']['INDEX']['GRCh'],
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass/{sample}/",
        twopass1readsN = config["STAR"].get("twopass1readsN", -1),
        gtf          = annotation_gtf_path,
        sjdbOverhang = config["STAR"]["sjdbOverhang"]
    shell:
        """
        mkdir -p star_2pass/{wildcards.sample}

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
