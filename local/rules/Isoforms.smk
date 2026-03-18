#############################
### Salmon quantification ###
#############################
GENCODE_DIR = globals().get("GENCODE_DIR", config.get("GENCODE_DIR", "gencode"))
transcriptome_fasta_path = globals().get(
    "transcriptome_fasta_path",
    config.get("transcriptome_fasta_path", f"{GENCODE_DIR}/transcriptome.fa"),
)
SAMPLES = globals().get("SAMPLES", config.get("SAMPLES", []))

GENCODE_ANNOTATION_GTF = globals().get(
    "GENCODE_ANNOTATION_GTF",
    config.get("GENCODE_ANNOTATION_GTF", f"{GENCODE_DIR}/gencode.v46.annotation.gtf"),
)

gff3 = globals().get(
    "gff3",
    config.get("gff3", f"{GENCODE_DIR}/gencode.v46.annotation.gff3"),
)

ref_seq = globals().get(
    "REF_SEQ",
    config.get("REF_SEQ", f"{GENCODE_DIR}/ref_seq.fa"),
)

REF_SEQ = globals().get(
    "REF_SEQ",
    config.get("REF_SEQ", f"{GENCODE_DIR}/ref_seq.fa"),
)

STRANDED = globals().get("STRANDED", config.get("STRANDED", "reverse"))

rule salmon_quant:
    input:
        fq=lambda wc: f"fastq/fastq_trimmed/{wc.sample}.fastq.gz",
        index=f"{GENCODE_DIR}/salmon_index",
    output:
        "salmon/{sample}/quant.sf",
    threads: 16
    conda:
        "transcript_env.yaml"
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
        transcriptome=transcriptome_fasta_path,
        bam=lambda wc: f"star_2pass_{GENOME_KEY}/{wc.sample}/Aligned.toTranscriptome.out.bam",
    output:
        "salmon_bam/{sample}/quant.sf",
    threads: 16
    conda:
        "transcript_env.yaml"
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
        bam=lambda wc: f"star_2pass_{GENOME_KEY}/{wc.sample}/Aligned.sortedByCoord.out.bam",
    output:
        fq="fastq_from_bam/{sample}.fastq.gz",
    threads: 8
    conda:
        "transcript_env.yaml"
    shell:
        """
        mkdir -p fastq_from_bam

        samtools fastq \
            -@ {threads} \
            -0 {output.fq} \
            {input.bam}
        """


###############################
### rule bam_to_fastq_ prova ###
###############################
rule bam_to_fastq_prova:
    input:
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.ribo.ex.H.bam",
    output:
        fq="fastq_from_bam_prova/{sample}.fastq.gz",
    threads: 8
    conda:
        "transcript_env.yaml"
    shell:
        """
        mkdir -p fastq_from_bam_prova

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
        fq=lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        index=f"{GENCODE_DIR}/kallisto_index/index_with_mask.idx",
    output:
        "kallisto/{sample}/abundance.tsv",
    threads: 8
    conda:
        "transcript_env.yaml"
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
        expand("kallisto/{sample}/abundance.tsv", sample=SAMPLES),
    output:
        "transcripts_tpm.tsv",
        "transcripts_counts.tsv",
    params:
        names=",".join(SAMPLES),
        script="dataset/Isella/merge_kallisto.R",
        files=lambda wc, input: ",".join(map(str, input)),
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


rule stringtie_assemble_H:
    input:
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.ribo.ex.H.bam",
        gtf="splitted_human.gtf",
    output:
        transcripts="stringtie/human/{sample}/transcripts.gtf",
        gene_abund="stringtie/human/{sample}/gene_abund.tab",
    threads: 4
    shell:
        """
        mkdir -p stringtie/human/{wildcards.sample}
        stringtie {input.bam} \
            -G {input.gtf} \
            -o {output.transcripts} \
            -A {output.gene_abund} \
            -p {threads}
        """


rule stringtie_merge_H:
    input:
        gtfs=expand("stringtie/human/{sample}/transcripts.gtf", sample=SAMPLES),
        gtf="splitted_human.gtf",
    output:
        merged="stringtie/human/merged.gtf",
    params:
        mergelist="stringtie/human/mergelist.txt",
    shell:
        """
        printf "%s\n" {input.gtfs} > {params.mergelist}

        stringtie --merge \
            -G {input.gtf} \
            -o {output.merged} \
            {params.mergelist}
        """


rule stringtie_quantify_H:
    input:
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.ribo.ex.H.bam",
        merged="stringtie/human/merged.gtf",
    output:
        quant="stringtie/human/{sample}/abund_merged.tab",
        gtf="stringtie/human/{sample}/quant_merged.gtf",
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


rule stringtie_prepDE_H:
    input:
        gtfs=expand("stringtie/human/{sample}/quant_merged.gtf", sample=SAMPLES),
    output:
        transcript_count_matrix="stringtie/human/transcript_count_matrix.csv",
    params:
        prepDEinput="stringtie/human/prepDE_input.txt",
    threads: 4
    shell:
        r"""
        > {params.prepDEinput}

        for gtf in {input.gtfs}; do
            sample_nam=$(basename $(dirname $gtf))
            dir=$(realpath $(dirname $gtf))
            echo -e "${{sample_nam}}\t${{gtf}}" >> {params.prepDEinput}
        done

        python ../../local/src/prepDE.py \
            -i {params.prepDEinput} \
            -t stringtie/human/transcript_count_matrix.csv
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
        bam="Results/bam_pass2_renamed/{sample}.bam",
        bai="Results/bam_pass2_renamed/{sample}.bam.bai",
        gtf=GENCODE_ANNOTATION_GTF,
    output:
        "out_spladder/spladder/genes_graph_conf3.{sample}.pickle",
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


# rule spladder_merge:
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


# rule spladder_quant:


###############################
### Sambamba sort and index BAM ###
###############################


rule sambamba_sort:
    input:
        "star/{sample}.bam",
    output:
        "sambamba/{sample}.sorted.bam",
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
        "sambamba/{sample}.sorted.bam",
    output:
        "sambamba/{sample}.sorted.bam.bai",
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
        expand("Majiq/SJ/{sample}.sj", sample=SAMPLES),


rule majiq_zar_generate:
    input:
        gff3=gff3,
    output:
        "Majiq/annotations/sg.zarr",
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
        zarr="Majiq/annotations/sg.zarr",
    output:
        sj="Majiq/SJ/{sample}.sj",
    threads: 4
    shell:
        """
        export MAJIQ_LICENSE_FILE="/mnt/oncog/software/snakeRnaSeq/dataset/v1/profiles/majiq_license_academic_official.lic"
        mkdir --p Majiq/SJ
        majiq-v3 sj {input.bam} {input.zarr} --nthreads {threads} {output.sj}
        """


rule majiq_cohort:
    input:
        zarr="Majiq/annotations/sg.zarr",
        sj=expand("Majiq/SJ/{sample}.sj", sample=SAMPLES),
    output:
        "Majiq/cohort/sg.zarr",
    threads: 4
    shell:
        """
        export MAJIQ_LICENSE_FILE="/mnt/oncog/software/snakeRnaSeq/d..."
        mkdir -p Majiq/cohort
        majiq-v3 build {input.zarr} {output} --sj {input.sj} --nthreads {threads} --min-experiments 1 --mindenovo 3 --no-simplify
        """


rule all_majiq_cohort:
    input:
        expand("Majiq/Results/{sample}_cohort.psicov", sample=SAMPLES),


rule majiq_build_cohort:
    input:
        cohort="Majiq/cohort/sg.zarr",
        sj=expand("Majiq/SJ/{sample}.sj", sample=SAMPLES),
    output:
        "Majiq/Results/{sample}_cohort.psicov",
    threads: 4
    shell:
        """
        export MAJIQ_LICENSE_FILE="/mnt/oncog/software/snakeRnaSeq/dataset/v1/profiles/majiq_license_academic_official.lic"
        mkdir -p Majiq/Results
        majiq-v3 build {input.cohort} {output} {input.sj} 
        """


rule all_majiq_quantify:
    input:
        expand("Majiq/Results/tsv/{sample}_cohort.tsv", sample=SAMPLES),


rule majiq_quantify:
    input:
        psicov="Majiq/Results/{sample}_cohort.psicov",
        zarr="Majiq/cohort/sg.zarr",
    output:
        "Majiq/Results/tsv/{sample}_cohort.tsv",
    threads: 4
    shell:
        """
        export MAJIQ_LICENSE_FILE="/mnt/oncog/software/snakeRnaSeq/dataset/v1/profiles/majiq_license_academic_official.lic"
        mkdir -p Majiq/Results/tsv
        majiq-v3 quanitfy {input.zarr} --splicegraph {input.psicov} --output-tsv {output}
        """


#####################################
# Rules for PsiCLASS quantification #
#####################################


rule all_xs:
    input:
        expand("xs_bams/{sample}.xs.sorted.bam", sample=SAMPLES),


rule addXS_prepare:
    input:
        bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.bam",
        ref_seq=REF_SEQ,
    output:
        xs_bam="xs_bams/{sample}.out.ribo.ex.H.bam",
        bai="xs_bams/{sample}.out.ribo.ex.H.bai",
    threads: 4
    shell:
        """
        samtools view -h {input.bam} | \
        /home/aleone/psiclass/addXS {input.ref_seq} | \
        samtools sort -@ {threads} -o {output.xs_bam} -
        
        samtools index -@ {threads} {output.xs_bam}
        """


rule create_bamlist:
    input:
        bams=expand(
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.bam",
            sample=SAMPLES,
        ),
    output:
        txt="PsiCLASS/bamlist.txt",
    run:
        os.makedirs("PsiCLASS", exist_ok=True)
        with open(output.txt, "w") as f:
            for b in input.bams:
                f.write(os.path.abspath(b) + "\n")


rule psiclass_cohort:
    input:
        bamlist="PsiCLASS/bamlist.txt",
        bams=expand(
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.bam",
            sample=SAMPLES,
        ),
    output:
        vote="PsiCLASS_noTI/combined_vote.gtf",
    threads: 4
    params:
        outprefix="PsiCLASS_noTI/combined",
        vd=2.0,
        sa=1.0,
        c=0.05,
        stranded=STRANDED,
    shell:
        """
        /home/aleone/psiclass/psiclass \
            --lb {input.bamlist} \
            -p {threads} \
            -o {params.outprefix} \
            --stranded {params.stranded} \
            --vd {params.vd} \
            --sa {params.sa} \
            -c {params.c}
        """
