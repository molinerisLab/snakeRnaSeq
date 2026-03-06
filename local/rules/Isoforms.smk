#############################
### Salmon quantification ###
#############################

rule salmon_quant:
    input:
        fq = lambda wc: f"fastq/fastq_trimmed/{wc.sample}.fastq.gz",
        index = f"{GENCODE_DIR}/salmon_index"
    output:
        "salmon/{sample}/quant.sf"
    threads: 16
    conda: "transcript_env.yaml"
    shell:
        """
        mkdir -p salmon/{wildcards.sample}

        salmon quant \
            -i {input.index} \
            -l A \
            -r {input.fq} \
            -p {threads} \
            --validateMappings \
            -o salmon/{wildcards.sample}
        """


rule salmon_quant_bam:
    input:
        transcriptome = transcriptome_fasta_path,
        bam = lambda wc: f"star_2pass_{GENOME_KEY}/{wc.sample}/Aligned.toTranscriptome.out.bam"
    output:
        "salmon_bam/{sample}/quant.sf"
    threads: 16
    conda: "transcript_env.yaml"
    shell:
        """
        mkdir -p salmon_bam/{wildcards.sample}

        salmon quant \
            -t {input.transcriptome} \
            -l A \
            -a {input.bam} \
            -p {threads} \
            --gencode \
            -o salmon_bam/{wildcards.sample}
        """




rule bam_to_fastq:
    input:
        bam = lambda wc: f"star_2pass_{GENOME_KEY}/{wc.sample}/Aligned.sortedByCoord.out.bam"
    output:
        fq = "fastq_from_bam/{sample}.fastq.gz"
    threads: 8
    conda: "transcript_env.yaml"
    shell:
        """
        mkdir -p fastq_from_bam

        samtools fastq \
            -@ {threads} \
            -0 {output.fq} \
            {input.bam}
        """




###############################
### Kallisto quantification ###
###############################

rule kallisto_quant:
    input:
        fq = lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        index = f"{GENCODE_DIR}/kallisto_index/index_with_mask.idx" 
    output:
        "kallisto/{sample}/abundance.tsv"
    threads: 8
    conda: "transcript_env.yaml"
    shell:
        """
        mkdir -p kallisto/{wildcards.sample}

        kallisto quant \
            -i {input.index} \
            -o kallisto/{wildcards.sample} \
            -b 100 \
            -t {threads} \
            --single -l 200 -s 30 \
            {input.fq}
        """


rule merge_kallisto_transcripts:
    input:
        expand("kallisto/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        "transcripts_tpm.tsv",
        "transcripts_counts.tsv"
    params:
        names=",".join(SAMPLES),
        script="dataset/Isella/merge_kallisto.R",
        files=lambda wc, input: ",".join(map(str, input))
    shell:
        """
        Rscript {params.script} \
            --input "{params.files}" \
            --names "{params.names}" \
            --output .
        """


################################
### StringTie quantification ###
################################

rule stringtie_assemble:
    input:
        bam = lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
        gtf = GENCODE_ANNOTATION_GTF
    output:
        transcripts = "stringtie/{sample}/transcripts.gtf",
        gene_abund = "stringtie/{sample}/gene_abund.tab"
    threads: 4
    shell:
        """
        mkdir -p stringtie/{wildcards.sample}
        stringtie {input.bam} \
            -G {input.gtf} \
            -o {output.transcripts} \
            -A {output.gene_abund} \
            -p {threads}
        """


rule stringtie_merge:
    input:
        gtfs = expand("stringtie/{sample}/transcripts.gtf", sample=SAMPLES),
        gtf  = GENCODE_ANNOTATION_GTF
    output:
        merged = "stringtie/merged.gtf"
    params:
        mergelist = "stringtie/mergelist.txt"
    shell:
        """
        printf "%s\n" {input.gtfs} > {params.mergelist}

        stringtie --merge \
            -G {input.gtf} \
            -o {output.merged} \
            {params.mergelist}
        """


rule stringtie_quantify:
    input:
        bam    = lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
        merged = "stringtie/merged.gtf"
    output:
        quant = "stringtie/{sample}/abund_merged.tab", gtf = "stringtie/{sample}/quant_merged.gtf"
    threads: 4
    shell:
        """
        stringtie {input.bam} \
            -e -B \
            -G {input.merged} \
            -o {output.gtf} \
            -A {output.quant} \
            -p {threads}
        """

rule stringtie_prepDE:
    input:
#        gtfs = expand("stringtie/{sample}/quant_merged.gtf", sample=SAMPLES),
        gtfs = [ancient(f"stringtie/{sample}/quant_merged.gtf") for sample in SAMPLES]

    output:
        transcript_count_matrix = "transcript_count_matrix.csv"
    params:
        prepDEinput    = "stringtie/prepDE_input.txt"
    threads: 4
    shell:
        r"""
        > {params.prepDEinput}

        for gtf in {input.gtfs}; do
            sample_nam=$(basename $(dirname $gtf))
            dir=$(realpath $(dirname $gtf))
            echo -e "${{sample_nam}}\t${{gtf}}" >> {params.prepDEinput}
        done

        python ../../local/src/prepDE.py -i {params.prepDEinput}
        """

#################
### Rename BAM ###
#################



######################
### Spladder Rules ###
######################
rule change_bamname: 
    shell: 
        """
        ../../local/src/create_bamlist.sh Results/bam_pass2_renamed 
        """

rule spladder_build:
    input:
        bam = "Results/bam_pass2_renamed/{sample}.bam",
        bai = "Results/bam_pass2_renamed/{sample}.bam.bai",
        gtf = GENCODE_ANNOTATION_GTF
    output:
        "out_spladder/spladder/genes_graph_conf3.{sample}.pickle"
    threads: 3
    shell:
        """
        mkdir -p out_spladder
        spladder build \
            -o out_spladder \
            -a {input.gtf} \
            -b {input.bam} \
            --confidence 3 \
            --merge-strat single \
            --no-extract-ase \
	    --ignore-mismatches \
            --parallel {threads}
        """


        

    

#rule spladder_merge:
#    input: 
#        gtf = GENCODE_ANNOTATION_GTF
#        
#    output: 
#    shell: 
#        """
#        spladder build \ 
#            -o out_spladder \ 
#            -a {input.gtf}
#            -b alignment.txt


#rule spladder_quant:





###############################
### Sambamba sort and index BAM ###
###############################


rule sambamba_sort:
    input:
        "star/{sample}.bam"
    output:
        "sambamba/{sample}.sorted.bam"
    threads: 8
    shell:
        """
        sambamba sort \
            -o {output} \
            -t {threads} \
            {input}
        """


rule sambamba_index:
    input:
        "sambamba/{sample}.sorted.bam"
    output:
        "sambamba/{sample}.sorted.bam.bai"
    threads: 4
    shell:
        """
        sambamba index \
            -t {threads} \
            {input}
        """

###############
# Majiq Rules #
###############
rule all_Majiq_SJ:
    input:
        expand("Majiq/SJ/{sample}.sj", sample=SAMPLES)

rule majiq_zar_generate:
    input:
        gff3 = gff3
    output:
        "Majiq/annotations/sg.zarr"
    threads: 4
    shell:
        """
        mkdir -p Majiq/annotations
        export MAJIQ_LICENSE_FILE="/mnt/oncog/software/snakeRnaSeq/dataset/v1/profiles/majiq_license_academic_official.lic"
        majiq-v3 gff3 {input.gff3}  {output}
        """



rule majiq_SJ: 
    input: 
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
        zarr= "Majiq/annotations/sg.zarr"
    output: 
        sj= "Majiq/SJ/{sample}.sj"
    threads: 4
    shell:
        """
        mkdir --p Majiq/SJ
        majiq-v3 sj {input.bam} {input.zarr} --nthreads {threads} {output.sj}
        """
