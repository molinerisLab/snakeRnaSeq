import os

def build_star_flags(d):
    """Turn a dict of {star_flag_name: value} into a CLI string.
    - value None or "" (or key absent)  -> flag omitted -> STAR uses its default
    - value is a list/tuple             -> space-joined (multi-value flags)
    - value is True                     -> bare flag (rare; only for true no-arg flags)
    Flag names are given WITHOUT the leading '--'.
    """
    if not d:
        return ""
    parts = []
    for flag, val in d.items():
        if val is None or val == "":
            continue                                  # unset -> STAR default
        if val is True:
            parts.append(f"--{flag}")                 # genuine no-arg flag
        elif isinstance(val, (list, tuple)):
            parts.append(f"--{flag} " + " ".join(str(x) for x in val))
        else:
            parts.append(f"--{flag} {val}")
    return " ".join(parts)


#############
# Ruleorder #
#############

_star_mode = config.get("STAR_MODE", "two_pass_manual")

if config["LAYOUT"] == "SINGLE":

    ruleorder: generate_unmapped_single > generate_unmapped_R1 > generate_unmapped_pairs
    ruleorder: generate_unmapped_single > generate_unmapped_R2 > generate_unmapped_pairs
    ruleorder: star_twopass_basic_se > star_twopass_basic_pe

elif config["LAYOUT"] == "PAIRED":

    ruleorder: generate_unmapped_pairs > generate_unmapped_R1 > generate_unmapped_single
    ruleorder: generate_unmapped_pairs > generate_unmapped_R2 > generate_unmapped_single
    ruleorder: star_twopass_basic_pe > star_twopass_basic_se

# link_unmapped only serves STAR-emitted FASTQ (single_pass / two_pass_manual).
# Everything else (two_pass_basic, bwa, SAVE_UNMAPPED=BAM) extracts from the BAM.
if (
    config["aligner"] == "star"
    and config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
    and _star_mode in {"single_pass", "two_pass_manual"}
):
    ruleorder: generate_unmapped_pairs > generate_unmapped_R1 > link_unmapped
    ruleorder: generate_unmapped_pairs > generate_unmapped_R2 > link_unmapped
    ruleorder: generate_unmapped_pairs > generate_unmapped_single > link_unmapped
else:
    ruleorder: generate_unmapped_R1 > link_unmapped
    ruleorder: generate_unmapped_R2 > link_unmapped
    ruleorder: generate_unmapped_single > link_unmapped


######
# LINK
######

_bam_sources = {
    "single_pass":      lambda wc: f"Results/star/{wc.sample}/Aligned.sortedByCoord.out.bam",
    "two_pass_manual":  lambda wc: f"Results/pass2/{wc.sample}/Aligned.sortedByCoord.out.bam",
    "two_pass_basic":   lambda wc: f"star_2pass/{wc.sample}/Aligned.sortedByCoord.out.bam",
}

# Where STAR writes its unmapped FASTQ per mode (two_pass_basic emits none —
# it is served by BAM extraction instead, see generate_unmapped_*).
_unmapped_sources = {
    "single_pass":     "Results/star/unmapped/{sample}_unmapped_R{mate}.fastq.gz",
    "two_pass_manual": "Results/pass2/unmapped/{sample}_unmapped_R{mate}.fastq.gz",
}

rule link_bam:
    input:
        _bam_sources[_star_mode],
    output:
        "Results/bam/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "Results/bam/{sample}/link_bam.log",
    conda:
        "transcript_env.yaml"
    shell:
        "ln -srf {input} {output} > {log} 2>&1"

##############
# STAR RULES #
##############


rule star_align_se:
    input:
        fq="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        aln="Results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="Results/star/{sample}/Log.out",
        sj="Results/star/{sample}/SJ.out.tab",
        unmapped=(
            "Results/star/unmapped/{sample}_unmapped_R1.fastq.gz"
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
        log_final="Results/star/{sample}/Log.final.out",
    log:
        "Results/star/{sample}/star.log",
    threads: config["CORES"]["star"]
    conda:
        "transcript_env.yaml"
    params:
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        tmpdir=lambda wc, output: os.path.dirname(output.aln),
        save_unmapped=config["STAR"]["SAVE_UNMAPPED"],
        extra=lambda wc: build_star_flags(config["STAR"].get("EXTRA_PASS2", {})),
    shell:
        """
        mkdir -p {params.tmpdir}

        STAR \
            --runThreadN {threads} \
            --limitBAMsortRAM 10000000000 \
            --genomeLoad LoadAndKeep \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --sjdbOverhang {params.sjdbOver} \
            --outFileNamePrefix {params.tmpdir}/ \
            --quantMode {params.quant_mode} \
            --outSAMtype {params.out_samtype} \
            --outSAMstrandField intronMotif \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --outSAMattributes All \
            {params.extra} \
            $([ "{params.save_unmapped}" = "FASTQ" ] && echo "--outReadsUnmapped Fastx --outSAMunmapped Within" || echo "") \
            > {log} 2>&1

        if [ "{params.save_unmapped}" = "FASTQ" ]; then
            mkdir -p Results/star/unmapped
            gzip -c {params.tmpdir}/Unmapped.out.mate1 > Results/star/unmapped/{wildcards.sample}_unmapped_R1.fastq.gz
        fi
        """


rule star_align_pe:
    input:
        fq1="fastq/fastq_trimmed/{sample}_R1.fastq.gz",
        fq2="fastq/fastq_trimmed/{sample}_R2.fastq.gz",
        idx=STAR_GENOME_DIR,
    output:
        aln="Results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="Results/star/{sample}/Log.out",
        sj="Results/star/{sample}/SJ.out.tab",
        unmapped=(
            [
                "Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
                "Results/star/unmapped/{sample}_unmapped_R2.fastq.gz",
            ]
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
        log_final="Results/star/{sample}/Log.final.out",
    log:
        "Results/star/{sample}/star.log",
    threads: config["CORES"]["star"]
    params:
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        tmpdir=lambda wc, output: os.path.dirname(output.aln),
        save_unmapped=config["STAR"]["SAVE_UNMAPPED"],
        extra=lambda wc: build_star_flags(config["STAR"].get("EXTRA_PASS2", {})),
        strand_field="--outSAMstrandField intronMotif",
        unmapped_flags=lambda wc: (
            "--outReadsUnmapped Fastx --outSAMunmapped Within"
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else ""
        ),
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        mkdir -p {params.tmpdir}
        mkdir -p Results/star/unmapped
        rm -rf {params.tmpdir}/_STARtmp

        STAR \
            --runThreadN {threads} \
            --genomeLoad NoSharedMemory \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --outFileNamePrefix {params.tmpdir}/ \
            --quantMode {params.quant_mode} \
            --outSAMtype {params.out_samtype} \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --outSAMattributes All \
            {params.strand_field} \
            {params.extra} \
            {params.unmapped_flags}

        if [ "{params.save_unmapped}" = "FASTQ" ]; then
            filter_unmapped() {{
                # keep only records whose trailing mate-status tag is 00
                paste - - - - < "$1" \
                  | awk -F'\t' '{{n=split($1,a," "); if (a[n]=="00") print $1"\n"$2"\n"$3"\n"$4}}' \
                  | gzip -c > "$2"
            }}

            filter_unmapped {params.tmpdir}/Unmapped.out.mate1 Results/star/unmapped/{wildcards.sample}_unmapped_R1.fastq.gz
            filter_unmapped {params.tmpdir}/Unmapped.out.mate2 Results/star/unmapped/{wildcards.sample}_unmapped_R2.fastq.gz

            # ---- depletion QC: one row per sample ----
            unmapped_records=$(( $(wc -l < {params.tmpdir}/Unmapped.out.mate1) / 4 ))
            kept=$(( $(gzip -dc Results/star/unmapped/{wildcards.sample}_unmapped_R1.fastq.gz | wc -l) / 4 ))
            kept2=$(( $(gzip -dc Results/star/unmapped/{wildcards.sample}_unmapped_R2.fastq.gz | wc -l) / 4 ))
            input_pairs=$(awk -F'\t' '/Number of input reads/ {{gsub(/ /,"",$2); print $2}}' \
                            {params.tmpdir}/Log.final.out)

            {{
              printf 'sample\tinput_pairs\tunmapped_records\tpairs_to_classifier\thalf_mapped_discarded\tpct_retained\n'
              printf '%s\t%s\t%s\t%s\t%s\t%.3f\n' \
                  "{wildcards.sample}" "$input_pairs" "$unmapped_records" "$kept" \
                  "$(( unmapped_records - kept ))" \
                  "$(awk -v a="$kept" -v b="$input_pairs" 'BEGIN{{print (b>0)? 100*a/b : 0}}')"
            }} > Results/star/unmapped/{wildcards.sample}_depletion_qc.tsv

            # kraken2 --paired requires equal counts in identical order.
            test "$kept" -eq "$kept2" \
                || {{ echo "ERROR: mate count mismatch R1=$kept R2=$kept2" >&2; exit 1; }}

            # A broken filter (awk syntax error, escape mangling) writes valid-but-EMPTY
            # .gz files and exits 0; the mate-count test above passes trivially (0 == 0).
            # This is the guard that actually catches it.
            test "$unmapped_records" -eq 0 -o "$kept" -gt 0 \
                || {{ echo "ERROR: filter kept 0 of $unmapped_records unmapped records" >&2; exit 1; }}
        fi

        rm -rf {params.tmpdir}/_STARtmp
        """


ruleorder: generate_unmapped_pairs > link_unmapped

rule link_unmapped:
    input:
        _unmapped_sources.get(_star_mode, _unmapped_sources["single_pass"]),
    output:
        "fastq/unmapped/{sample}_unmapped_R{mate}.fastq.gz",
    log:
        "fastq/unmapped/{sample}_unmapped_R{mate}.log",
    conda:
        "transcript_env.yaml"
    shell:
        "ln -srf $(realpath {input}) {output}"


rule generate_unmapped_single:
    input:
        lambda wildcards: (
            f"Results/bam/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}.aligned.bam"
        ),
    output:
        "fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    conda:
        "transcript_env.yaml"
    log:
        "fastq/unmapped/{sample}_unmapped_single.log",
    shell:
        """
        samtools view -f 4 {input} | awk '{{print "@"$1; print $10; print "+"; print $11}}' | gzip > {output}
        """


rule generate_unmapped_R1:
    input:
        lambda wildcards: (
            f"Results/bam/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
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

rule all_unmapped_pe: 
    input: 
        expand("fastq/unmapped/{sample}_unmapped_R1.fastq.gz", sample=SAMPLES),
        expand("fastq/unmapped/{sample}_unmapped_R2.fastq.gz", sample=SAMPLES),


rule generate_unmapped_R2:
    input:
        lambda wildcards: (
            f"Results/bam/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
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
rule all_unmap_paired: 
    input: 
        expand("fastq/unmapped/{sample}_unmapped_R1.fastq.gz", sample=SAMPLES),
        expand("fastq/unmapped/{sample}_unmapped_R2.fastq.gz", sample=SAMPLES),

rule generate_unmapped_pairs:
    input:
        bam=lambda wildcards: (
            f"Results/bam/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
            if config["aligner"] == "star"
            else f"aligned_bwa/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
        )
    output:
        r1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
        r2="fastq/unmapped/{sample}_unmapped_R2.fastq.gz"
    threads: 4
    conda:
        "transcript_env.yaml"
    log:
        "fastq/unmapped/{sample}_unmapped.log"
    shell:
        r"""
        set -euo pipefail

        samtools view \
            -u \
            -f 12 \
            -F 2304 \
            {input.bam} 2>> {log} \
        | samtools collate \
            -u -O - 2>> {log} \
        | samtools fastq \
            -@ {threads} \
            -n \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            - 2>> {log}
        """

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
    threads: config["CORES"]["star"]
    conda:
        "transcript_env.yaml"
    log:
        "Results/pass1/{sample}/{sample}_star_pass1.log",
    params:
        tmpdir=lambda wc, output: os.path.dirname(output.sj),
        read_cmd=config["STAR"]["readFilesCommand"],
        limitSjdb=config["STAR"]["limitSjdbInsertNsj"],
        extra=lambda wc: build_star_flags(config["STAR"].get("EXTRA_PASS1", {})),
    shell:
        """
        mkdir -p {params.tmpdir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --limitSjdbInsertNsj {params.limitSjdb} \
            --outFileNamePrefix {params.tmpdir}/{wildcards.sample}_ \
            --outSAMtype None \
            {params.extra} \
            > {log} 2>&1
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
        sj="Results/pass1/merged_filtered_SJ.out.tab",
    output:
        bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.bam",
        gene_counts="Results/pass2/{sample}/ReadsPerGene.out.tab",
        log_final="Results/pass2/{sample}/Log.final.out",
        unmapped=(
            [
                "Results/pass2/unmapped/{sample}_unmapped_R1.fastq.gz",
                "Results/pass2/unmapped/{sample}_unmapped_R2.fastq.gz",
            ]
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ" and config["LAYOUT"] == "PAIRED"
            else (
                ["Results/pass2/unmapped/{sample}_unmapped_R1.fastq.gz"]
                if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
                else []
            )
        ),
    threads: config["CORES"]["star"]
    conda:
        "transcript_env.yaml"
    log:
        "Results/pass2/{sample}/{sample}_star_pass2.log",
    params:
        out_samtype=config["STAR"]["OUT_SAM_TYPE"],
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        quant_mode=config["STAR"]["quantMode"],
        sjdbOver=config["STAR"]["sjdbOverhang"],
        read_cmd=config["STAR"]["readFilesCommand"],
        limitSjdb=config["STAR"]["limitSjdbInsertNsj"],
        tmpdir=lambda wc, output: os.path.dirname(output.bam),
        save_unmapped=config["STAR"]["SAVE_UNMAPPED"],
        extra=lambda wc: build_star_flags(config["STAR"].get("EXTRA_PASS2", {})),
    shell:
        """
        mkdir -p {params.tmpdir}
        mkdir -p Results/pass2/unmapped

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --outSAMstrandField intronMotif \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            --sjdbFileChrStartEnd {input.sj} \
            --sjdbOverhang {params.sjdbOver} \
            --limitSjdbInsertNsj {params.limitSjdb} \
            --genomeLoad NoSharedMemory \
            --outFileNamePrefix {params.tmpdir}/ \
            --outSAMtype {params.out_samtype} \
            --quantMode {params.quant_mode} \
            $([ "{params.save_unmapped}" = "FASTQ" ] && echo "--outReadsUnmapped Fastx --outSAMunmapped Within" || echo "") \
            {params.extra} \
            > {log} 2>&1

        if [ "{params.save_unmapped}" = "FASTQ" ]; then
            if [ -f "{params.tmpdir}/Unmapped.out.mate1" ]; then
                gzip -c {params.tmpdir}/Unmapped.out.mate1 > Results/pass2/unmapped/{wildcards.sample}_unmapped_R1.fastq.gz
            fi
            if [ -f "{params.tmpdir}/Unmapped.out.mate2" ]; then
                gzip -c {params.tmpdir}/Unmapped.out.mate2 > Results/pass2/unmapped/{wildcards.sample}_unmapped_R2.fastq.gz
            fi
        fi
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
    threads: config["CORES"]["star"]
    conda:
        "transcript_env.yaml"
    log:
        "star_2pass/{sample}/star_2pass_basic_se.log",
    params:
        read_cmd=config["STAR"]["readFilesCommand"],
        out_dir=lambda wc, output: os.path.dirname(output.bam),
        gtf=GENCODE_ANNOTATION_GTF,
        sjdbOverhang=config["STAR"]["sjdbOverhang"],
        outsamtype = config["STAR"]["OUT_SAM_TYPE"], 
        quantmode= config["STAR"]["quantMode"]
    shell:
        """
        mkdir -p {params.out_dir}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq} \
            --readFilesCommand {params.read_cmd} \
            --outSAMstrandField intronMotif \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_dir}/ \
            --outSAMtype {params.outsamtype} \
            --outSAMunmapped Within \
            --quantMode {params.quantmode}
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
    threads: config["CORES"]["star"]
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
        outsamtype = config["STAR"]["OUT_SAM_TYPE"], 
        quantmode= config["STAR"]["quantMode"], 
    shell:
        """
        mkdir -p {params.out_dir}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --outSAMstrandField intronMotif \
            --twopassMode Basic \
            --twopass1readsN {params.twopass1readsN} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --outFileNamePrefix {params.out_dir}/ \
            --outSAMtype {params.outsamtype} \
            --outSAMunmapped Within \
            --quantMode {params.quantmode}
        """

rule all_star_chm: 
    input: 
        expand("Results/star_CHM/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)


rule star_align_pe_CHM:
    input:
        fq1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
        fq2="fastq/unmapped/{sample}_unmapped_R2.fastq.gz",
        idx=config["STAR"]["GENOME_DIR_CHM"],
    output:
        aln="Results/star_CHM/{sample}/Aligned.sortedByCoord.out.bam",
        star_log="Results/star_CHM/{sample}/Log.out",
        sj="Results/star_CHM/{sample}/SJ.out.tab",
        log_final="Results/star_CHM/{sample}/Log.final.out",
        unmapped=(
            [
                "Results/star_CHM/unmapped/{sample}_unmapped_R1.fastq.gz",
                "Results/star_CHM/unmapped/{sample}_unmapped_R2.fastq.gz",
            ]
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
        depletion_qc=(
            ["Results/star_CHM/unmapped/{sample}_depletion_qc.tsv"]
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else []
        ),
    log:
        "Results/star_CHM/{sample}/star.log",
    threads:
        config["CORES"]["star"]
    params:
        outdir=lambda wc: f"Results/star_CHM/{wc.sample}",
        outfiltermultimapnmax=config["STAR"]["OUT_FILTER_MULTIMAP_NMAX"],
        read_cmd=config["STAR"]["readFilesCommand"],
        save_unmapped=config["STAR"]["SAVE_UNMAPPED"],
        extra=lambda wc: build_star_flags(
            config["STAR"].get("EXTRA_PASS2", {})
        ),
        unmapped_flags=lambda wc: (
            "--outSAMunmapped Within"
            if config["STAR"]["SAVE_UNMAPPED"] == "FASTQ"
            else ""
        ),
    conda:
        "transcript_env.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="{params.outdir}"
        unmapped_dir="Results/star_CHM/unmapped"
        star_tmp="${{outdir}}/_STARtmp"

        mkdir -p "$outdir" "$unmapped_dir"
        exec > {log} 2>&1

        rm -rf -- "$star_tmp"
        trap 'rm -rf -- "$star_tmp"' EXIT

        STAR \
            --runThreadN {threads} \
            --genomeLoad NoSharedMemory \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand {params.read_cmd} \
            --outFileNamePrefix "$outdir/" \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax {params.outfiltermultimapnmax} \
            {params.extra} \
            {params.unmapped_flags}

        samtools quickcheck -v {output.aln}

        if [ "{params.save_unmapped}" = "FASTQ" ]; then
            r1="$unmapped_dir/{wildcards.sample}_unmapped_R1.fastq.gz"
            r2="$unmapped_dir/{wildcards.sample}_unmapped_R2.fastq.gz"
            qc="$unmapped_dir/{wildcards.sample}_depletion_qc.tsv"

            samtools view -u -f 12 -F 2304 {output.aln} \
            | samtools collate -u -O - \
            | samtools fastq \
                -@ {threads} \
                -n \
                -1 "$r1" \
                -2 "$r2" \
                -0 /dev/null \
                -s /dev/null \
                -

            gzip -t "$r1"
            gzip -t "$r2"

            input_pairs=$(
                awk -F'\t' \
                    '/Number of input reads/ {{
                        gsub(/ /, "", $2)
                        print $2
                    }}' \
                    {output.log_final}
            )

            both_unmapped_records=$(
                samtools view -c -f 12 -F 2304 {output.aln}
            )

            test $(( both_unmapped_records % 2 )) -eq 0

            expected_pairs=$(( both_unmapped_records / 2 ))
            kept_r1=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
            kept_r2=$(( $(gzip -dc "$r2" | wc -l) / 4 ))

            test -n "$input_pairs"
            test "$kept_r1" -eq "$kept_r2"
            test "$kept_r1" -eq "$expected_pairs"
            test "$kept_r1" -le "$input_pairs"

            pct_retained=$(
                awk \
                    -v retained="$kept_r1" \
                    -v input="$input_pairs" \
                    'BEGIN {{
                        printf "%.3f", (input > 0 ? 100 * retained / input : 0)
                    }}'
            )

            printf \
                'sample\tinput_pairs\tpairs_to_classifier\tpct_retained_after_chm13\n%s\t%s\t%s\t%s\n' \
                "{wildcards.sample}" \
                "$input_pairs" \
                "$kept_r1" \
                "$pct_retained" \
                > "$qc"
        fi
        """