#############################
### Salmon quantification ###
#############################
REFERENCE_DIR = config["REFERENCE_DIR"]

rule salmon_quant:
    input:
        fq = lambda wc: f"fastq/fastq_trimmed/{wc.sample}.fastq.gz",
        index = f"{REFERENCE_DIR}/salmon_index"
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
        index = f"{REFERENCE_DIR}/kallisto_index/index_with_mask.idx"
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
        bam = lambda wc: f"Results/pass2_{GENOME_KEY}/{wc.sample}/Aligned.sortedByCoord.out.bam",
        gtf = annotation_gtf_path
    output:
        "stringtie/{sample}/transcripts.gtf"
    threads: 8
    shell:
        """
        mkdir -p stringtie/{wildcards.sample}
        stringtie {input.bam} \
            -G {input.gtf} \
            -o {output} \
            -p {threads}
        """


rule stringtie_merge:
    input:
        gtfs = expand("stringtie/{sample}/transcripts.gtf", sample=SAMPLES),
        gtf  = annotation_gtf_path
    output:
        merged = "stringtie/merged/merged.gtf"
    params:
        mergelist = "stringtie/merged/mergelist.txt"
    shell:
        """
        mkdir -p stringtie/merged
        printf "%s\n" {input.gtfs} > {params.mergelist}

        stringtie --merge \
            -G {input.gtf} \
            -o {output.merged} \
            {params.mergelist}
        """


rule stringtie_quantify:
    input:
        bam    = lambda wc: f"Results/pass2_{GENOME_KEY}/{wc.sample}/Aligned.sortedByCoord.out.bam",
        merged = "stringtie/merged/merged.gtf"
    output:
        quant = "stringtie/{sample}/quant/abund.tab"
    params:
        outdir = "stringtie/{sample}/quant"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        stringtie {input.bam} \
            -e -B \
            -G {input.merged} \
            -o {params.outdir}/transcripts.gtf \
            -A {output.quant} \
            -p {threads}
        """





#################
### BAM index ###
#################
rule index_bam:
    input:
        "Results/pass2_{genome}/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "Results/pass2_{genome}/{sample}/Aligned.sortedByCoord.out.bam.bai"
    threads: 1
    shell:
        "samtools index {input}"


######################
### Spladder Rules ###
######################


rule spladder_build:
    input:
        bam = "Results/pass2_{genome}/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Results/pass2_{genome}/{sample}/Aligned.sortedByCoord.out.bam.bai",
        gtf = annotation_gtf_path
    output:
        "spladder/{genome}/{sample}/spladder/genes_graph_conf3.pickle"
    threads: 8
    shell:
        """
        mkdir -p spladder/{wildcards.genome}/{wildcards.sample}

        spladder build \
            -o spladder/{wildcards.genome}/{wildcards.sample} \
            -a {input.gtf} \
            -b {input.bam} \
            --confidence 3 \
            --ignore-mismatches \
            --parallel {threads}
        """



rule spladder_merge:
    input:
        graphs = expand(
            "spladder/{{genome}}/{sample}/spladder/genes_graph_conf3.pickle",
            sample=SAMPLES
        )
    output:
        "spladder/merged_{genome}/spladder/genes_graph_conf3.merge_graph.pickle"
    shell:
        """
        mkdir -p spladder/merged_{wildcards.genome}

        spladder merge \
            -o spladder/merged_{wildcards.genome} \
            -g {input.graphs}
        """


rule spladder_quant:
    input:
        graph = "spladder/merged_{genome}/spladder/genes_graph_conf3.merge_graph.pickle",
        bams  = expand(
            "Results/pass2/{{genome}}/{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        )
    output:
        "spladder/quantification_{genome}/spladder/genes_graph_conf3.quant.pickle"
    threads: 8
    shell:
        """
        mkdir -p spladder/quantification_{wildcards.genome}

        spladder quantify \
            -o spladder/quantification_{wildcards.genome} \
            -g {input.graph} \
            -b {input.bams} \
            --parallel {threads}
        """




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