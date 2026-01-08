#############################
### Salmon quantification ###
#############################

rule salmon_quant:
    input:
        fq = lambda wc: f"{FASTQ_DIR}/fastq_trimmed/{wc.sample}.fastq.gz",
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
        fq = lambda wc: f"{FASTQ_DIR}/fastq_trimmed/{wc.sample}.fastq.gz",
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
        expand("kallisto/{sample}/abundance.tsv", sample=config["samples"])
    output:
        "transcripts_tpm.tsv",
        "transcripts_counts.tsv"
    params:
        names=",".join(config["samples"]),
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
        gtfs = expand("stringtie/{sample}/transcripts.gtf", sample=config["samples"]),
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
            sample=config["samples"]
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
            sample=config["samples"]
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



#########################
### Ribosomal Removal ###
#########################

rule split_bam_ribo:
    input:
        bam = "{path}.bam",
        bam_idx = "{path}.bam.bai",
        ribosome_bed = f"{REFERENCE_DIR}/primary_assembly.annotation.rRNA_complete.bed"
    output:
        ribo_ex = "{path}.ribo.ex.bam",
        ribo_in = "{path}.ribo.in.bam",
        ribo_log = "{path}.summary"
    shell:
        """
        mkdir -p $(dirname {output.ribo_ex})

        split_bam.py \
            -i {input.bam} \
            -r {input.ribosome_bed} \
            -o {wildcards.path}.ribo \
            > {output.ribo_log}
        """


ruleorder: featurecounts > split_bam_ribo

rule featurecounts:
    input:
        bam = "{path}.bam",
        annotation_gtf = annotation_gtf_path
    output:
        counts  = "{path}.bam.featurecounts.count",
        summary = "{path}.bam.featurecounts.count.summary"
    params:
        cores  = config["CORES"],
        tmpdir = config["TMPDIR"],
        opts   = config["FEATURECOUNTS_PARAM"]
    shell:
        """
        featureCounts {input.bam} \
            -o {output.counts} \
            -a {input.annotation_gtf} \
            {params.opts} \
            --tmpDir {params.tmpdir} \
            -T {params.cores}
        """

rule featurecounts_ribo_ex:
    input:
        file_all = expand("star/{sample}.ribo.ex.bam.featurecounts.count", sample=config["samples"]),
        file_translate = lambda wc: f"star/{config['samples'][0]}.ribo.ex.bam.featurecounts.count"
    output:
        "featurecounts.ribo.ex.count.gz"
    shell:
        """
        matrix_reduce '*.ribo.ex.bam.featurecounts.count' -l '{input.file_all}' \
            | grep -v '^#' \
            | fasta2tab \
            | bawk '$2!="Geneid" {{print $2,$1,$8}}' \
            | tab2matrix -r Geneid \
            | translate -a <(cut -f -6 {input.file_translate} | unhead) 1 \
            | gzip > {output}
        """


rule featurecounts_ribo_ex_summary:
    input:
        expand("star/{sample}.ribo.ex.bam.featurecounts.count.summary", sample=config["samples"])
    output:
        "fastq.featurecounts.ribo.ex.count.gz.summary_matrix"
    shell:
        """
        matrix_reduce -t 'star/*.ribo.ex.bam.featurecounts.count.summary' \
            | grep -v Status \
            | tab2matrix -r Sample \
            > {output}
        """

rule featurecounts_ribo_ex_summary_matrix:
    input:
        "fastq.featurecounts.ribo.ex.count.gz.summary_matrix"
    output:
        "fastq.featurecounts.ribo.ex.count.gz.summary_matrix.reduced"
    shell:
        """
        matrix2tab {input} \
            | bawk '$2=="Unassigned_Ambiguity" || $2=="Assigned" || $2=="Unassigned_NoFeatures"' \
            | tab2matrix -r Sample \
            > {output}
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