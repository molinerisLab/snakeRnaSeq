#############################
### Salmon quantification ###
#############################


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

################################
### StringTie quantification ###
################################


rule stringtie_assemble:
    input:
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
        gtf=GENCODE_ANNOTATION_GTF,
    output:
        transcripts="stringtie/{sample}/transcripts.gtf",
        gene_abund="stringtie/{sample}/gene_abund.tab",
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
        gtfs=expand("stringtie/{sample}/transcripts.gtf", sample=SAMPLES),
        gtf=GENCODE_ANNOTATION_GTF,
    output:
        merged="stringtie/merged.gtf",
    params:
        mergelist="stringtie/mergelist.txt",
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
        bam=lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
        merged="stringtie/merged.gtf",
    output:
        quant="stringtie/{sample}/abund_merged.tab",
        gtf="stringtie/{sample}/quant_merged.gtf",
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
        gtfs=[ancient(f"stringtie/{sample}/quant_merged.gtf") for sample in SAMPLES],
    output:
        transcript_count_matrix="transcript_count_matrix.csv",
    params:
        prepDEinput="stringtie/prepDE_input.txt",
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


###########
# PsiCLASS
###########


rule t_introns:
    input:
        gtf=GENCODE_ANNOTATION_GTF,
    output:
        txt="t_introns.txt",
    threads: 2
    run:
        import re
        from collections import defaultdict

        exons = defaultdict(list)
        with open(input.gtf, "r") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, src, feat, start, end, score, strand, frame, attrs = fields
                if feat != "exon":
                    continue
                m = re.search(r'transcript_id "([^"]+)"', attrs)
                if not m:
                    continue
                tid = m.group(1)
                exons[(chrom, strand, tid)].append((int(start), int(end)))

        introns = set()
        for (chrom, strand, tid), xs in exons.items():
            xs.sort()
            for (s1, e1), (s2, e2) in zip(xs, xs[1:]):
                istart = e1
                iend = s2
                if istart < iend:
                    introns.add((chrom, istart, iend, strand))

        with open(output.txt, "w") as w:
            for chrom, istart, iend, strand in sorted(introns):
                # Columns: 0=chr, 1=start, 2=end, 3=dummy, 4=strand
                w.write(f"{chrom}\t{istart}\t{iend}\t.\t{strand}\n")

        print(f"Wrote {len(introns)} introns -> {output.txt}")


rule all_xs:
    input:
        expand("xs_bams/{sample}.xs.sorted.bam", sample=SAMPLES),


rule addXS_prepare:
    input:
        bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.bam",
        ref_seq=GENCODE_GENOME_FASTA,
    output:
        xs_bam="xs_bams/{sample}.out.ribo.ex.bam",
        bai="xs_bams/{sample}.out.ribo.ex.bam.bai",
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
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.unique.bam",
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
        trusted="t_introns.txt" if config.get("USE_TRUSTED_INTRONS", False) else [],
        bams=expand(
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.unique.bam",
            sample=SAMPLES,
        ),
    output:
        vote="PsiCLASS/combined_vote.gtf",
    threads: 4
    conda:
        "../env/PsiCLASS.yaml"
    params:
        outprefix="PsiCLASS/combined",
        vd=2.0,
        sa=1.0,
        c=0.05,
        stranded=config["psitrand"],
        trusted_flag=lambda wc, input: f"-s {input.trusted}" if config.get("USE_TRUSTED_INTRONS", False) else "",
    shell:
        """
        psiclass \
            --lb {input.bamlist} \
            {params.trusted_flag} \
            -p {threads} \
            -o {params.outprefix} \
            --stranded {params.stranded} \
            --vd {params.vd} \
            --sa {params.sa} \
            -c {params.c}
        """




rule annotate_assembly:
    input:
        gtf="PsiCLASS/combined_vote.gtf",
        ref=GENCODE_ANNOTATION_GTF,
    output:
        annotated="Results/gffcmp/gffcmp.annotated.gtf",
        stats="Results/gffcmp/gffcmp.stats",
        tracking="Results/gffcmp/gffcmp.tracking",
        loci="Results/gffcmp/gffcmp.loci",
    params:
        prefix="gffcmp",
    shell:
        """
        mkdir -p Results/gffcmp
        gffcompare -r {input.ref} -o Results/gffcmp/{params.prefix} {input.gtf} 
        """


rule filter_kallisto_gtf:
    input:
        annotated_gtf="Results/gffcmp/gffcmp.annotated.gtf",
        original_gtf="PsiCLASS/combined_vote.gtf" 
    output:
        final_gtf="kallisto_output/cohort_kallisto_reference.gtf",
    params:
        min_sample_cnt=5,
    threads: 1
    run:
        import re
        import os

        os.makedirs(os.path.dirname(output.final_gtf), exist_ok=True)

        def get_attr(attrs, key):
            m = re.search(f'{key} "([^"]+)"', attrs)
            return m.group(1) if m else None

        # PASS 1: Map original transcript IDs to their sample counts using the PsiCLASS GTF
        sample_counts = {}
        with open(input.original_gtf, "r") as f_orig:
            for line in f_orig:
                if line.startswith("#") or "\ttranscript\t" not in line:
                    continue
                attrs = line.split("\t")[8]
                tid = get_attr(attrs, "transcript_id")
                cnt = get_attr(attrs, "sample_cnt")
                if tid and cnt:
                    sample_counts[tid] = int(cnt)

        trusted_transcripts = set()

        # PASS 2: Find non "=" transcripts in the annotated GTF that pass the threshold
        with open(input.annotated_gtf, "r") as f_in:
            for line in f_in:
                if line.startswith("#") or "\ttranscript\t" not in line:
                    continue
                attrs = line.split("\t")[8]
                code = get_attr(attrs, "class_code")
                
                # gffcompare moves the original PsiCLASS ID to 'oId'
                orig_id = get_attr(attrs, "oId") or get_attr(attrs, "transcript_id")
                curr_tid = get_attr(attrs, "transcript_id")
                if code != "=" and orig_id in sample_counts:
                    if sample_counts[orig_id] >= params.min_sample_cnt:
                        trusted_transcripts.add(curr_tid)

        with open(input.annotated_gtf, "r") as f_in, open(output.final_gtf, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    f_out.write(line)
                    continue
                
                curr_tid = get_attr(line.split("\t")[8], "transcript_id")
                if curr_tid and curr_tid in trusted_transcripts:
                    f_out.write(line)

        print(
            f"Kept {len(trusted_transcripts)} transcripts "
            f"(sample_cnt >= {params.min_sample_cnt})"
        )

rule extract_novel_transcripts:
    input:
        gtf="kallisto_output/cohort_kallisto_reference.gtf",
        genome=GENCODE_GENOME_FASTA  
    output:
        novel_fa="kallisto_output/cohort_transcriptome.fasta" 
    shell:
        """
        gffread -w {output.novel_fa} -g {input.genome} {input.gtf}
        """

rule build_kallisto_index_combined:
    input:
        novel_fa="kallisto_output/cohort_transcriptome.fasta",
        ref_tx_fa=transcriptome_fasta_path
    output:
        combined_fa="kallisto_output/combined_transcriptome.fa",
        idx="kallisto_output_combined/kallisto.idx"
    threads: 4
    shell:
        """
        mkdir -p kallisto_output_combined
        # 1. Concatenate the reference transcriptome with your novel transcripts
        cat {input.ref_tx_fa} {input.novel_fa} > {output.combined_fa}

        # 2. Build the index
        kallisto index -i {output.idx} {output.combined_fa}
        """


rule all_kaPSI:
    input:
        "transcripts_counts.tsv",
        "transcripts_tpm.tsv",



rule kallisto_quant_PsiCLASS:
    input:
        fq=lambda wc: f"fastq/fastq_trimmed/{wc.sample}_R1.fastq.gz",
        index="kallisto_output_combined/kallisto.idx",
    output:
            tsv="kallisto_PsiCLASS_idx_combined/{sample}/abundance.tsv",
            h5="kallisto_PsiCLASS_idx_combined/{sample}/abundance.h5",

    threads: 4
    shell:
        """
        mkdir -p kallisto_PsiCLASS_idx_combined/{wildcards.sample}

        kallisto quant \
            -i {input.index} \
            -o kallisto_PsiCLASS_idx_combined/{wildcards.sample} \
            -b 100 \
            -t {threads} \
            --single -l 200 -s 30 \
            {input.fq}
        """


rule merge_kallisto_transcripts_psiclass:
    input:
        expand("kallisto_PsiCLASS_idx_combined/{sample}/abundance.h5", sample=SAMPLES),
    output:
        counts="transcripts_counts.tsv",
        tpm="transcripts_tpm.tsv",
    params:
        names=",".join(SAMPLES),
        script="../../local/src/merge_kallisto.R",
        files=lambda wc, input: ",".join(map(str, input)),
    shell:
        """
        Rscript {params.script} \
            --input "{params.files}" \
            --names "{params.names}" \
            --output . \
            --verbose

        # Insert 'target_id' at the beginning of the first line for both output files
        sed -i '1s/^/target_id\t/' {output.counts}
        sed -i '1s/^/target_id\t/' {output.tpm}
        """



rule DREAMSEQ: 
    input: 
        "transcripts_counts.tsv"
    output:
        "DREAMSEQ/dreamseq_results.tsv"
    threads: 4
    shell: 
        """
        mkdir -p DREAMSEQ
        Rscript run_dreamseq.R --input {input} --output {output} --threads {threads}
        """


rule compare_filtered_to_gencode:
    input:
        query="kallisto_output/cohort_kallisto_reference.gtf",
        ref="Resources/gencode46_Hprefixed.annotation.gtf",
    output:
        annotated="Results/gffcmp_filtered/gffcmp_filtered.annotated.gtf",
        stats="Results/gffcmp_filtered/gffcmp_filtered.stats",
        tracking="Results/gffcmp_filtered/gffcmp_filtered.tracking",
        loci="Results/gffcmp_filtered/gffcmp_filtered.loci",
    params:
        prefix="gffcmp_filtered",
    shell:
        """
        mkdir -p Results/gffcmp_filtered
        gffcompare -r {input.ref} -o Results/gffcmp_filtered/{params.prefix} {input.query}
        """


rule filter_tmm:
    input:
        counts="{file}.tsv.gz"
    output:
        matrix="{file}.tmm.tsv.gz",
        factors="{file}.tmm.factors.tsv.gz"
    params:
        min_cpm=config["DGE"]["EXPRESSED_GENES_MIN_CPM"],
        min_samples= config["DGE"]["MIN_NUM_OF_EXPRESSED_SAMPLE"]
    shell:
        """
        Rscript ../../local/src/filter_tmm.R \
            {input.counts} \
            {output.matrix} \
            {output.factors} \
            {params.min_cpm} \
            {params.min_samples}
        """

rule plot_kallisto_compare:
    input: 
        "transcripts_counts.tmm.tsv.gz",
        "kallisto_only_transcripts_counts.tmm.tsv.gz"
    output: 
        "plots/kallisto_comparison_scatter.png",
        "plots/kallisto_comparison_scatter.pdf",
        
    shell: 
        "Rscript ../../local/src/compare_kallisto.R"



    
# Similar to above, it retains the header as well
rule transcripts_counts_not_in_kallisto:
    input:
        fp1="kallisto_only_transcripts_counts.tmm.tsv.gz",
        fp2="transcripts_counts.tmm.tsv.gz"
    output:
        fpFinal="transcripts_counts_exclude_kallisto_only.tmm.tsv"
    shell:
        """
        awk -F'\\t' 
            NR==FNR {{ id[$1]; next }}
            FNR==1  {{ print; next }}
            !($1 in id)
        ' <(zcat {input.fp1}) <(zcat {input.fp2}) > {output.fpFinal}
        """

rule plot_transcript_counts:
    input:
        tsv="{file}.tsv"
    output:
        plot="{file}_intersect_counts_plot.png"
    shell:
        """
        Rscript -e "
            library(ggplot2);
            library(data.table);
            
            # (salta la prima colonna degli ID se sono solo conte)
            df <- fread('{input.tsv}', header=FALSE);
            
            matrice_conte <- as.matrix(df[, -1, with=FALSE]);
            log_conte <- log2(matrice_conte + 1);
            
            png('{output.plot}', width=800, height=600);
            plot(density(log_conte), main='Distribuzione delle Conte (Log2)', 
                 xlab='Log2(Counts + 1)', col='blue', lwd=2);
            dev.off();
        "
        """

        

   

        