
if config["LAYOUT"] == "SINGLE":
     ruleorder: generate_unmapped_single > generate_unmapped_R1
     ruleorder: generate_unmapped_single > generate_unmapped_R2
elif config["LAYOUT"] == "PAIRED":
     ruleorder: generate_unmapped_R1 > generate_unmapped_single
     ruleorder: generate_unmapped_R2 > generate_unmapped_single

rule star_align_se:
    input:
        fq1 = "fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        idx = STAR_GENOME_DIR
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
#            f"{'--outReadsUnmapped Fastx --outSAMunmapped None' if config['STAR']['SAVE_UNMAPPED'] == 'FASTQ' else '--outSAMunmapped Within'}"
        )
    wrapper:
        "v3.3.6/bio/star/align"

rule star_align_pe:
    input:
        fq1 = "fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        fq2 = "fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        idx = STAR_GENOME_DIR
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
#            f"{'--outReadsUnmapped Fastx --outSAMunmapped None' if config['STAR']['SAVE_UNMAPPED'] == 'FASTQ' else '--outSAMunmapped Within'}"
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
# STAR DOUBLE PASS ALIGNMENT USING ALL SJ OUT 
# ==============================================

rule load_star_genome:
    input:
        idx = STAR_GENOME_DIR
    output:
        flag = touch("Results/genome_loaded.flag")
    conda: "transcript_env.yaml"
    shell:
        """
        STAR --genomeDir {input.idx} --genomeLoad LoadAndExit
        """

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
        idx= STAR_GENOME_DIR,
        mem_flag = "Results/genome_loaded.flag"
    output:
        sj           = "Results/pass1/{sample}/SJ.out.tab",
        log          = "Results/pass1/{sample}/Log.out",
        log_final    = "Results/pass1/{sample}/Log.final.out",
        log_progress = "Results/pass1/{sample}/Log.progress.out",
    threads: 3
    conda: "transcript_env.yaml"
    params:
        tmpdir     = "Results/pass1/{sample}",
        read_cmd   = config["STAR"]["readFilesCommand"],
        fq_join    = lambda wc, input: " ".join(input.fq),
        sjdbOver    = config["STAR"]["sjdbOverhang"],
        multiscorerange = config["STAR"]["MULTIMAP_SCORE_RANGE"],
        outfiltermismatch = config["STAR"]["OUT_FILTER_MISMATCH_NMAX"],
        outfiltermultimapnmax = config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        outfiltermismatchnover = config["STAR"]["OUT_FILTER_MISMATCH_NOVER_LMAX"],
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeLoad LoadAndKeep \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq_join} \
            --readFilesCommand {params.read_cmd} \
            --outFileNamePrefix {params.tmpdir}/ \
            --sjdbOverhang {params.sjdbOver} \
            --outFilterMultimapScoreRange {params.multiscorerange} \
            --outFilterMismatchNoverLmax {params.outfiltermismatchnover} \
            --outFilterMismatchNmax {params.outfiltermismatch} \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --outSAMtype None
        """


        
rule filter_sj:
    input:
        expand("Results/pass1/{sample}/SJ.out.tab",
               sample= SAMPLES)
    output:
        sj_filtered = "Results/pass1/{sample}/SJ_filtered.out.tab"
    conda: "transcript_env.yaml"
    threads: 3
    shell:
        """
        cat {input} \
          | awk '$1!="HchrM" && $1!="MchrM" && $5>0 && $6==0 && $7>2' \
          | sort -u \
          > {output.sj_filtered}
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
        idx=STAR_GENOME_DIR,
        sj_files = expand("Results/pass1/{sample}/SJ_filtered.out.tab",
                          sample=SAMPLES)
    output:
        bam         = "Results/pass2/{sample}/Aligned.sortedByCoord.out.bam",
        gene_counts = "Results/pass2/{sample}/ReadsPerGene.out.tab"
    threads: 3
    conda: "transcript_env.yaml"
    params:
        multiscorerange = config["STAR"]["MULTIMAP_SCORE_RANGE"],
        outfiltermismatch = config["STAR"]["OUT_FILTER_MISMATCH_NMAX"],
        outfiltermultimapnmax = config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        outfiltermismatchnover = config["STAR"]["OUT_FILTER_MISMATCH_NOVER_LMAX"],
        out_samtype = config["STAR"]["OUT_SAM_TYPE"],
        quant_mode  = config["STAR"]["quantMode"],
        sjdbOver    = config["STAR"]["sjdbOverhang"],
        read_cmd    = config["STAR"]["readFilesCommand"],
        fq_join     = lambda wc, input: " ".join(input.fq),
        sj_join     = lambda wc, input: " ".join(input.sj_files),
        tmpdir      = "Results/pass2/{sample}"
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeLoad LoadAndKeep \
            --genomeDir {input.idx} \
            --readFilesIn {params.fq_join} \
            --readFilesCommand {params.read_cmd} \
            --sjdbFileChrStartEnd {params.sj_join} \
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

rule unload_star_genome:
    input:
        bams = expand("Results/pass2/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        flag = touch("Results/genome_unloaded.flag")
    conda: "transcript_env.yaml"
    shell:
        """
        STAR --genomeDir {STAR_GENOME_DIR} --genomeLoad Remove
        """





# ===============================================================================================
# STAR DOUBLE PASS ALIGNMENT USING A SINGLE SJ OUT EACH TIME (NOT RECOMMENDED FOR LARGE DATASETS)
# ===============================================================================================



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
        idx = STAR_GENOME_DIR,
    output:
        bam        = "star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass/{sample}/SJ.out.tab",
        gene_counts= "star_2pass/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass/{sample}/Log.out",
        log_final  = "star_2pass/{sample}/Log.final.out"
    threads: 6
    conda: "transcript_env.yaml"
    params:
        genome_dir   = STAR_GENOME_DIR,
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass/{sample}/",
        gtf          = GENCODE_ANNOTATION_GTF,
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
        idx = STAR_GENOME_DIR,
    output:
        bam        = "star_2pass/{sample}/Aligned.sortedByCoord.out.bam",
        sj         = "star_2pass/{sample}/SJ.out.tab",
        gene_counts= "star_2pass/{sample}/ReadsPerGene.out.tab",
        log        = "star_2pass/{sample}/Log.out",
        log_final  = "star_2pass/{sample}/Log.final.out"
    threads: 8
    conda: "transcript_env.yaml"
    params:
        genome_dir   = STAR_GENOME_DIR,
        read_cmd     = config["STAR"]["readFilesCommand"],
        out_prefix   = "star_2pass/{sample}/",
        twopass1readsN = config["STAR"].get("twopass1readsN", -1),
        gtf          = GENCODE_ANNOTATION_GTF,
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
