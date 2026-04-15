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
        import os
        import re
        from collections import defaultdict

        os.makedirs(PSI_DIR, exist_ok=True)
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
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.unique.bam",
            sample=SAMPLES,
        ),
    output:
        txt="PsiCLASS/bamlist.txt",
    run:
        os.makedirs("PsiCLASS", exist_ok=True)
        with open(output.txt, "w") as f:
            for b in input.bams:
                f.write(os.path.abspath(b) + "\n")


rule psiclass_cohort:  #todo, put the trusted introns as config, if i want it or not to be used
    input:
        bamlist="PsiCLASS/bamlist.txt",
        bams=expand(
            "Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.unique.bam",
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
    shell:
        """
        psiclass \
            --lb {input.bamlist} \
            -p {threads} \
            -o {params.outprefix} \
            --stranded {params.stranded} \
            --vd {params.vd} \
            --sa {params.sa} \
            -c {params.c}
        """


rule gtf_to_bed:
    input:
        gtf=GENCODE_ANNOTATION_GTF,
    output:
        bed="Resources/annotation/gencode.H.bed12",
    conda:
        "../env/PsiCLASS.yaml"
    log:
        "logs/gtf_to_bed.log",
    shell:
        """
    gtfToGenePred -genePredExt -geneNameAsName2 {input.gtf} /dev/stdout 2>{log} \
        | genePredToBed /dev/stdin /dev/stdout \
        | sed 's/^chr/Hchr/' \
        > {output.bed} 2>>{log}
        """


rule infer_experiment_psi:
    input:
        bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.unique.bam".format(
            sample=SAMPLES[0]
        ),
        bed="Resources/annotation/gencode.H.bed12",
    output:
        txt="QC/infer_experiment/{sample}.infer_experiment.txt".format(
            sample=SAMPLES[0]
        ),
    conda:
        "../env/PsiCLASS.yaml"
    log:
        "logs/infer_experiment/{sample}.log".format(sample=SAMPLES[0]),
    shell:
        """
        infer_experiment.py \
            -r {input.bed} \
            -i {input.bam} \
            > {output.txt} 2>{log}
        """


rule annotate_assembly:
    input:
        gtf="PsiCLASS/combined_vote.gtf",
        ref="Resources/gencode46_Hprefixed.annotation.gtf",
    output:
        annotated="gffcmp_noTI.annotated.gtf",
        stats="gffcmp_noTI.stats",
        tracking="gffcmp_noTI.tracking",
        loci="gffcmp_noTI.loci",
    params:
        prefix="gffcmp_noTI",
    shell:
        """
        gffcompare -r {input.ref} -o {params.prefix} {input.gtf} 
        """


rule filter_kallisto_gtf:
    input:
        annotated_gtf="gffcmp_noTI.annotated.gtf",
        original_gtf="PsiCLASS/combined_vote.gtf"  # <--- Crucial: Need the original for the counts!
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

        # PASS 2: Find "j" transcripts in the annotated GTF that pass the threshold
        with open(input.annotated_gtf, "r") as f_in:
            for line in f_in:
                if line.startswith("#") or "\ttranscript\t" not in line:
                    continue
                attrs = line.split("\t")[8]
                code = get_attr(attrs, "class_code")
                
                # gffcompare moves the original PsiCLASS ID to 'oId'
                orig_id = get_attr(attrs, "oId") or get_attr(attrs, "transcript_id")
                curr_tid = get_attr(attrs, "transcript_id")

                # We ONLY want 'j' transcripts for your strategy
                if code == "j" and orig_id in sample_counts:
                    if sample_counts[orig_id] >= params.min_sample_cnt:
                        trusted_transcripts.add(curr_tid)

        # PASS 3: Write out the transcripts and their exons to the final GTF
        with open(input.annotated_gtf, "r") as f_in, open(output.final_gtf, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    f_out.write(line)
                    continue
                
                curr_tid = get_attr(line.split("\t")[8], "transcript_id")
                if curr_tid and curr_tid in trusted_transcripts:
                    f_out.write(line)

        print(
            f"Kept {len(trusted_transcripts)} novel 'j' transcripts "
            f"(sample_cnt >= {params.min_sample_cnt})"
        )



rule build_kallisto_index_combined:
    input:
        novel_fa="kallisto_output/novel_transcripts.fa",
        # Your pre-existing GENCODE transcriptome FASTA (with H prefixes)
        ref_tx_fa="Resources/gencode.v46.Hchr_transcripts.fa" 
    output:
        combined_fa="kallisto_output/combined_transcriptome.fa",
        idx="kallisto_output_combined/kallisto.idx"
    threads: 4
    shell:
        """
        # 1. Concatenate the reference transcriptome with your novel transcripts
        cat {input.ref_tx_fa} {input.novel_fa} > {output.combined_fa}

        # 2. Build the index
        kallisto index -i {output.idx} {output.combined_fa}
        """


rule build_kallisto_index:
    input:
        gtf="kallisto_output/cohort_kallisto_reference.gtf",  
        genome=GENCODE_GENOME_FASTA,
    output:
        fasta="kallisto_output/cohort_transcriptome.fasta",
        index="kallisto_output/cohort_transcriptome.idx",
    threads: 4
    run:
        import os
        import re
        import subprocess

        # Make sure output directory exists
        os.makedirs(os.path.dirname(output.fasta), exist_ok=True)

        # 1. run gffread into a temp fasta
        temp_fasta = output.fasta + ".tmp"
        cmd = f"sed 's/^Hchr/chr/' {input.gtf} | gffread -w {temp_fasta} -g {input.genome} -"
        subprocess.check_call(cmd, shell=True)

        if not os.path.exists(temp_fasta) or os.path.getsize(temp_fasta) == 0:
            raise Exception("ERROR: gffread failed to produce FASTA")

        # 2. Extract annotations from GTF
        tx_info = {}
        with open(input.gtf, "r") as f:
            for line in f:
                if line.startswith("#") or "\\ttranscript\\t" not in line:
                    continue
                parts = line.strip().split("\\t")
                attrs = parts[8]

                def get_attr(key):
                    m = re.search(f'{key} "([^"]+)"', attrs)
                    return m.group(1) if m else "-"

                tid = get_attr("transcript_id")
                gene_id = get_attr("gene_id")
                ref_gene_id = get_attr("ref_gene_id")
                cmp_ref = get_attr("cmp_ref")
                gene_name = get_attr("gene_name")
                class_code = get_attr("class_code")

                # Use ref_gene_id if available, otherwise fallback to the assembled gene_id
                g_id = ref_gene_id if ref_gene_id != "-" else gene_id

                # Format mimic GENCODE:
                # transcript_id|gene_id|HAVANA_gene|HAVANA_transcript|transcript_name|gene_name|length|biotype|
                piped = f"{tid}|{g_id}|-|{cmp_ref}|{tid}|{gene_name}|-|{class_code}|"
                tx_info[tid] = piped

        # 3. Rewrite FASTA headers
        with open(temp_fasta, "r") as fin, open(output.fasta, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    # gffread headers look like >transcript_id gene_id
                    tid = line.strip().split()[0][1:]
                    if tid in tx_info:
                        fout.write(f">{tx_info[tid]}\\n")
                    else:
                        fout.write(line)
                else:
                    fout.write(line)

        # clean up temp fasta
        os.remove(temp_fasta)

        # 4. Build the Kallisto index
        index_cmd = f"kallisto index -i {output.index} {output.fasta}"
        subprocess.check_call(index_cmd, shell=True)


rule all_kaPSI:
    input:
        expand("kallisto_PsiCLASS_idx_combined/{sample}/abundance.tsv", sample=SAMPLES),


rule kallisto_quant_PsiCLASS:
    input:
        fq=lambda wc: f"fastq/{wc.sample}_R1.fastq.gz",
        index="kallisto_output_combined/kallisto.idx",
    output:
        "kallisto_PsiCLASS_idx_combined/{sample}/abundance.tsv",
    threads: 4
    conda:
        "transcript_env.yaml"
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
        "transcripts_tpm_psiclass.tsv",
        "transcripts_counts_psiclass.tsv",
    params:
        names=",".join(SAMPLES),
        script="merge_kallisto.R",
        files=lambda wc, input: ",".join(map(str, input)),
    shell:
        """
        Rscript {params.script} \
            --input "{params.files}" \
            --names "{params.names}" \
            --output . \
            --verbose
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