##########################################################################
# Helper function to fetch assembly accessions for download based on TaxID
##########################################################################

import pandas as pd

import os

if os.path.exists("../../local/bin/GTDB_NCBI_bact_table.tsv"):
    TAX_MAP = pd.read_csv("../../local/bin/GTDB_NCBI_bact_table.tsv", sep="\t", dtype=str).set_index("kraken_taxid")
else:
    TAX_MAP = pd.DataFrame()

# Create a helper function to dynamically fetch the accession
def get_accession_for_download(wildcards):
    taxid_str = str(wildcards.taxid)
    if taxid_str in TAX_MAP.index:
        return TAX_MAP.loc[taxid_str, "assembly_accession"]
    else:
        raise ValueError(f"CRITICAL FAULT: TaxID {taxid_str} is missing from kraken_master_map.tsv")

if config["LAYOUT"] == "PAIRED":
    ruleorder: kraken_pe_pass1 > kraken_se_pass1
    ruleorder: remove_taxid_reads_pe > remove_taxid_reads_se
    ruleorder: kraken_pe_pass2 > kraken_se_pass2
else:
    ruleorder: kraken_se_pass1 > kraken_pe_pass1
    ruleorder: remove_taxid_reads_se > remove_taxid_reads_pe
    ruleorder: kraken_se_pass2 > kraken_pe_pass2

ruleorder: minimap2_merge > minimap2_index

_unclassified_se = (
    "fastq/unclassified/{sample}_unclassified_R1.fastq.gz"
    if config.get("SAVE_UNCLASSIFIED", False)
    else temp("fastq/unclassified/{sample}_unclassified_R1.fastq.gz")
)
_unclassified_pe_r1 = (
    "fastq/unclassified/{sample}_unclassified_R1.fastq.gz"
    if config.get("SAVE_UNCLASSIFIED", False)
    else temp("fastq/unclassified/{sample}_unclassified_R1.fastq.gz")
)
_unclassified_pe_r2 = (
    "fastq/unclassified/{sample}_unclassified_R2.fastq.gz"
    if config.get("SAVE_UNCLASSIFIED", False)
    else temp("fastq/unclassified/{sample}_unclassified_R2.fastq.gz")
)


##############################################################################
## 0. TARGETS
##############################################################################

rule all_metagenome:
    input:
        expand("alignments_merged/{taxid}/merged_all_samples.bam.bai", taxid=config["kraken_extract_taxid"])


##############################################################################
## 1. K2 DAEMON LIFECYCLE
##############################################################################

rule k2_daemon_start:
    """Load the kraken2 database into RAM"""
    input:
        expand("fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz", sample=SAMPLES),
    output:
        temp(".k2_daemon.sentinel")
    params:
        db=config["kraken_db"],
        use_daemon=config.get("USE_K2_DAEMON", False),
    shell:
        """
        if [ "{params.use_daemon}" = "True" ]; then
            k2 clean --stop-daemon 2>/dev/null || true
            echo "" | k2 classify --db {params.db} --use-daemon \
                --output /dev/null --report /dev/null - 2>/dev/null || true
        fi
        touch {output}
        """

rule k2_daemon_stop:
    """Stop the daemon after all pass2 classifications finish."""
    input:
        kreports=expand("kreports_filtered/{sample}.k2report", sample=SAMPLES),
        sentinel=".k2_daemon.sentinel",
    output:
        temp(".k2_daemon_stopped.sentinel")
    params:
        use_daemon=config.get("USE_K2_DAEMON", False),
    shell:
        """
        if [ "{params.use_daemon}" = "True" ]; then
            k2 clean --stop-daemon || true
        fi
        touch {output}
        """


##############################################################################
## 2. KRAKEN2 CLASSIFICATION  (pass1 -> host/contaminant depletion -> pass2)
##############################################################################

rule kraken_pe_pass1:
    input:
        R1="Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
        R2="Results/star/unmapped/{sample}_unmapped_R2.fastq.gz",
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2",
        unclassified1="unclassified/{sample}_unclassified_1.fq",
        unclassified2="unclassified/{sample}_unclassified_2.fq",
    threads: config["CORES"]["kraken2"]
    shell:
        """
        k2 classify --db {config[kraken_db_pass1]} {config[kraken_options]}\
            --threads {threads} \
            --report {output.report} \
            --output {output.out} \
            --paired {input.R1} {input.R2} \
            --unclassified-out unclassified/{wildcards.sample}_unclassified#.fq
    """

rule kraken_se_pass1:
    input:
        R1="Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2",
    threads: config["CORES"]["kraken2"]
    shell:
        """
        k2 classify --db {config[kraken_db_pass1]} {config[kraken_options]} \
            --threads {threads} \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            {input.R1}
        """

rule remove_taxid_reads_se:
    input:
        kraken2_output="koutputs/{sample}.kraken2",
        kraken2_report="kreports/{sample}.k2report",
        fastq_r1="Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        fastq_r1_nonhost="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    params:
        taxids=lambda wc: " ".join(
            [
                str(t)
                for t in config.get("CONTAMINANT_INFO", {}).get("humanID", [])
                + config.get("CONTAMINANT_INFO", {}).get("contaminant", [])
            ]
        ),
        exclude_opt="--exclude",
        children_opt="--include-children",
    shell:
        r"""
        mkdir -p fastq/fastq_taxid_depleted
        out="{output.fastq_r1_nonhost}"
        tmp=$(printf '%s\n' "$out" | sed 's/\.gz$//')

        extract_kraken_reads.py \
            -k {input.kraken2_output} \
            -s1 {input.fastq_r1} \
            -t {params.taxids} \
            -r {input.kraken2_report} \
            {params.exclude_opt} \
            {params.children_opt} \
            -o "$tmp" \
            --fastq-output

        gzip -f "$tmp"
        """

rule remove_taxid_reads_pe:
    input:
        kraken2_output="koutputs/{sample}.kraken2",
        kraken2_report="kreports/{sample}.k2report",
        fastq_r1="Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
        fastq_r2="Results/star/unmapped/{sample}_unmapped_R2.fastq.gz", 
    output:
        fastq_r1_nonhost="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
        fastq_r2_nonhost="fastq/fastq_taxid_depleted/{sample}_R2.fastq.gz", 
    params:
        taxids=lambda wc: " ".join(
            [
                str(t)
                for t in config.get("CONTAMINANT_INFO", {}).get("humanID", [])
                + config.get("CONTAMINANT_INFO", {}).get("contaminant", [])
            ]
        ),
        exclude_opt="--exclude",
        children_opt="--include-children",
    shell:
        r"""
        mkdir -p fastq/fastq_taxid_depleted
        
        # Strip the .gz extension to create temporary uncompressed targets
        out1="{output.fastq_r1_nonhost}"
        tmp1=$(printf '%s\n' "$out1" | sed 's/\.gz$//')
        
        out2="{output.fastq_r2_nonhost}"
        tmp2=$(printf '%s\n' "$out2" | sed 's/\.gz$//')

        # Execute extraction for both strands
        extract_kraken_reads.py \
            -k {input.kraken2_output} \
            -s1 {input.fastq_r1} \
            -s2 {input.fastq_r2} \
            -t {params.taxids} \
            -r {input.kraken2_report} \
            {params.exclude_opt} \
            {params.children_opt} \
            -o "$tmp1" \
            -o2 "$tmp2" \
            --fastq-output

        # Compress both files in parallel
        gzip -f "$tmp1" "$tmp2"
        """

rule kraken_se_pass2:
    input:
        R1="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
        daemon=".k2_daemon.sentinel",
    output:
        report="kreports_filtered/{sample}.k2report",
        out="koutput_filtered/{sample}.kraken2",
        unclassified=_unclassified_se,
    params:
        use_daemon=config.get("USE_K2_DAEMON", False),
    threads: config["CORES"]["kraken2"]
    shell:
        """
        mkdir -p fastq/unclassified
        k2 classify --db {config[kraken_db]} {config[kraken_options]} \
            --threads {threads} \
            $([ "{params.use_daemon}" = "True" ] && echo "--use-daemon" || echo "") \
            --report-minimizer-data \
            --unclassified-out fastq/unclassified/{wildcards.sample}_unclassified_R1.fq \
            --report {output.report} \
            --output {output.out} \
            {input.R1}
        gzip -f fastq/unclassified/{wildcards.sample}_unclassified_R1.fq
        """

rule kraken_pe_pass2:
    input:
        R1="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
        R2="fastq/fastq_taxid_depleted/{sample}_R2.fastq.gz",
        daemon=".k2_daemon.sentinel",
    output:
        report="kreports_filtered/{sample}.k2report",
        out="koutput_filtered/{sample}.kraken2",
        unclassified_r1=_unclassified_pe_r1,
        unclassified_r2=_unclassified_pe_r2,
    params:
        use_daemon=config.get("USE_K2_DAEMON", False),
    threads: config["CORES"]["kraken2"]
    shell:
        """
        mkdir -p fastq/unclassified
        k2 classify --db {config[kraken_db]} {config[kraken_options]} \
            --threads {threads} \
            $([ "{params.use_daemon}" = "True" ] && echo "--use-daemon" || echo "") \
            --report-minimizer-data \
            --unclassified-out fastq/unclassified/{wildcards.sample}_unclassified_#.fq \
            --report {output.report} \
            --output {output.out} \
            --paired {input.R1} {input.R2}
        gzip -f fastq/unclassified/{wildcards.sample}_unclassified_1.fq
        gzip -f fastq/unclassified/{wildcards.sample}_unclassified_R2.fq
        """

"""
.META: *.kreport
    1   frag_perc   Percentage of fragments covered by the clade rooted at this taxon
    2   frag_num    Number of fragments covered by the clade rooted at this taxon
    3   frag_num_direct Number of fragments assigned directly to this taxon
    4   level   A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
    5   tax_id  NCBI taxonomic ID number
    6   name    Indented scientific name

https://github.com/DerrickWood/kraken2/wiki/Manual#classification
"""


##############################################################################
## 3. BRACKEN + ABUNDANCE FILTERING & STATS
##############################################################################

rule braken:
    input:
        report="kreports/{sample}.k2report",
    output:
        report="breports/{sample}.breport",
        out="boutputs/{sample}.braken",
    params:
        db=config["kraken_db"],
        read_len=config["BRACKEN"]["braken_read_len"],
        level=config["BRACKEN"]["braken_level"],
        min_reads=config["BRACKEN"]["braken_min_reads"],
    shell:
        """
        bracken -d {params.db} -i {input.report} \
            -r {params.read_len} -l {params.level} -t {params.min_reads} \
            -o {output.out} -w {output.report}
        """

rule bracken_merged:
    input:
        reports=expand("breports/{sample}.breport", sample=SAMPLES),
        outputs=expand("boutputs/{sample}.braken", sample=SAMPLES),
    output:
        "bracken_not_filtered_merged_abbundances.txt",
    log:
        "log/Bracken_merged.log",
    shell:
        """
        combine_bracken_outputs.py --files {input.outputs} -o {output} 2> {log} 
        """

rule filbraken:
    input:
        report="kreports_filtered/{sample}.k2report",
        daemon_stopped=".k2_daemon_stopped.sentinel",  
    output:
        report="breports_filtered/{sample}.breport",
        out="boutputs_filtered/{sample}.braken",
    params:
        db=config["kraken_db"],
        read_len=config["BRACKEN"]["braken_read_len"],
        level=config["BRACKEN"]["braken_level"],
        min_reads=config["BRACKEN"]["braken_min_reads"],
    shell:
        """
        mkdir -p breports_filtered boutputs_filtered
        bracken -d {params.db} -i {input.report} \
            -r {params.read_len} -l {params.level} -t {params.min_reads} \
            -o {output.out} -w {output.report}
        """

rule filbracken_merged:
    input:
        reports=expand("breports_filtered/{sample}.breport", sample=SAMPLES),
        outputs=expand("boutputs_filtered/{sample}.braken", sample=SAMPLES),
    output:
        "bracken_merged_abbundances.txt",
    log:
        "log/bracken_merged.log",
    shell:
        """
        combine_bracken_outputs.py --files {input.outputs} -o {output} 2> {log} 
        """

rule split_merged:
    input:
        "bracken_merged_abbundances.txt",
    output:
        num="bracken_merged_abbundances.num.txt",
        frac="bracken_merged_abbundances.frac.txt",
    shell:
        """
        grep_columns -k 1,2,3 braken_num  < {input} | perl -pe '$.==1; s/.braken_num//g'  > {output.num};
        grep_columns -k 1,2,3 braken_frac < {input} | perl -pe '$.==1; s/.braken_frac//g' > {output.frac}
        """

rule feature_filter:
    input:
        abundances="bracken_merged_abbundances.num.txt",
        metadata="metadata.txt",
    output:
        filtered="bracken_merged_abbundances.num.filtered.txt",
        tmp=temp("bracken_merged_abbundances.num.tmp"),
    params:
        condition=config["feature_filter"]["condition"],
        g1=config["feature_filter"]["g1"],
        g2=config["feature_filter"]["g2"],
        use_raw_counts=config["feature_filter"]["use_raw_counts"],
        min_exp=config["feature_filter"]["min_exp"],
        min_samples_ratio=config["feature_filter"]["min_samples_ratio"],
    shell:
        """
        perl -pe 's/\t/;/; s/\t/;/' {input.abundances} > {output.tmp}
        echo -en "name\\ttaxonomy_id\\ttaxonomy_lv\\t" > {output.filtered}
        feature_filter {output.tmp} {input.metadata} \
            --condition={params.condition} --g1 {params.g1} --g2 {params.g2} \
            $([ "{params.use_raw_counts}" = "True" ] && echo "--use_raw_counts" || echo "") \
            --min_exp {params.min_exp} --min_samples_ratio {params.min_samples_ratio} \
        | perl -pe 's/;/\t/; s/;/\t/' >> {output.filtered}
        """

rule filter_frac:
    input:
        num_filter="bracken_merged_abbundances.num.filtered.txt",
        frac="bracken_merged_abbundances.frac.txt",
    output:
        "bracken_merged_abbundances.frac.filtered.txt",
    shell:
        "filter_1col 3 <(cut -f 3 {input.num_filter}) < {input.frac} > {output}"

rule collapse_taxid:
    input:
        "bracken_merged_abbundances.num.filtered.txt",
    output:
        "bracken_merged_abbundances.num.filtered.taxid_collapsed.txt",
    shell:
        "perl -pe 's/\t/;/; s/\t/;/' {input} > {output}"

rule degw:
    input:
        abundances="bracken_merged_abbundances.num.filtered.taxid_collapsed.txt",
        metadata="metadata.txt",
    output:
        "bracken_merged_abbundances.num.filtered.taxid_collapsed.degw.txt",
    params:
        condition=config["feature_filter"]["condition"],
        g1=config["feature_filter"]["g1"],
        g2=config["feature_filter"]["g2"],
        use_raw_counts=config["feature_filter"]["use_raw_counts"],
        min_exp=config["feature_filter"]["min_exp"],
        min_samples_ratio=config["feature_filter"]["min_samples_ratio"],
    shell:
        """
        DEGWilcox.R {input.abundances} {input.metadata} \
            --condition={params.condition} --g1 {params.g1} --g2 {params.g2} \
            $([ "{params.use_raw_counts}" = "True" ] && echo "--use_raw_counts" || echo "") \
            --min_exp {params.min_exp} --min_samples_ratio {params.min_samples_ratio} \
        | bawk 'NR==1 {{$1="name\ttaxonomy_id\tlevel"; print}} NR>1{{gsub(/;/, "\t", $1); print}}' > {output}
    """


##############################################################################
## 4. TAXID EXTRACTION -> REFERENCE FETCH -> MINIMAP2 MAPPING
##############################################################################

rule extract_kraken_reads:
    input:
        kraken2_output="koutput_filtered/{sample}.kraken2",
        kraken2_report="kreports_filtered/{sample}.k2report",
        fastq_r1="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
        fastq_r2="fastq/fastq_taxid_depleted/{sample}_R2.fastq.gz" if config["LAYOUT"] == "PAIRED" else []
    output:
        fastq_r1_extracted="fastq_idmapped/{taxid}/{sample}_R1.fastq.gz",
        fastq_r2_extracted="fastq_idmapped/{taxid}/{sample}_R2.fastq.gz" if config["LAYOUT"] == "PAIRED" else []
    params:
        layout=config["LAYOUT"]
    shell:
        r"""
        mkdir -p fastq_idmapped/{wildcards.taxid}
        
        out1="{output.fastq_r1_extracted}"
        tmp1=$(printf '%s\n' "$out1" | sed 's/\.gz$//')

        
        if [ "{params.layout}" = "PAIRED" ]; then
            out2="{output.fastq_r2_extracted}"
            tmp2=$(printf '%s\n' "$out2" | sed 's/\.gz$//')
            
            extract_kraken_reads.py \
                -k {input.kraken2_output} \
                --taxid {wildcards.taxid} \
                --include-children \
                -s {input.fastq_r1} \
                -s2 {input.fastq_r2} \
                --report {input.kraken2_report} \
                --fastq-output \
                -o "$tmp1" \
                -o2 "$tmp2"
                
            gzip -f "$tmp1" "$tmp2"
        else
            extract_kraken_reads.py \
                -k {input.kraken2_output} \
                --taxid {wildcards.taxid} \
                --include-children \
                -s {input.fastq_r1} \
                --report {input.kraken2_report} \
                --fastq-output \
                -o "$tmp1"
                
            gzip -f "$tmp1"
        fi
        """

rule fetch_gtdb_representative:
    output:
        fasta="Resources/genomes/{taxid}/genome.fna",
        fai="Resources/genomes/{taxid}/genome.fna.fai",
        gff="Resources/genomes/{taxid}/annotation.gff" 
    params:
        accession=get_accession_for_download
    shell:
        """
        datasets download genome accession {params.accession} --include genome,gff3 --filename {params.accession}.zip
        unzip -p {params.accession}.zip "ncbi_dataset/data/{params.accession}/*.fna" > {output.fasta}
        unzip -p {params.accession}.zip "ncbi_dataset/data/{params.accession}/genomic.gff" > {output.gff}
        rm {params.accession}.zip
        samtools faidx {output.fasta}
        """

rule minimap2_align:
    """Align extracted reads with minimap2 short-read preset to the dynamic reference."""
    input:
        fq1="fastq_idmapped/{taxid}/{sample}_R1.fastq.gz",
        fq2="fastq_idmapped/{taxid}/{sample}_R2.fastq.gz" if config["LAYOUT"] == "PAIRED" else [],
        ref="Resources/genomes/{taxid}/genome.fna"
    output:
        bam="alignments/{taxid}/{sample}.bam",
    params:
        layout=config["LAYOUT"]
    threads: 8
    shell:
        """
        if [ "{params.layout}" = "PAIRED" ]; then
            # Feed both fq1 and fq2 to minimap2
            minimap2 -ax sr -t {threads} --secondary=no \
                {input.ref} {input.fq1} {input.fq2} \
                | samtools view -bS -F 4 \
                | samtools sort -o {output.bam}
        else
            # Feed only fq1
            minimap2 -ax sr -t {threads} --secondary=no \
                {input.ref} {input.fq1} \
                | samtools view -bS -F 4 \
                | samtools sort -o {output.bam}
        fi
        """

rule minimap2_index:
    """Universal indexer for single-sample BAMs."""
    input:
        "alignments/{taxid}/{sample}.bam"
    output:
        "alignments/{taxid}/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule minimap2_merge:
    """Merge all per-sample minimap2 BAMs into a single file for IGV/coverage."""
    input:
        bams=expand("alignments/{{taxid}}/{sample}.bam", sample=SAMPLES),
        bais=expand("alignments/{{taxid}}/{sample}.bam.bai", sample=SAMPLES)
    output:
        bam="alignments_merged/{taxid}/merged_all_samples.bam",
        bai="alignments_merged/{taxid}/merged_all_samples.bam.bai"
    shell:
        """
        mkdir -p alignments_merged/{wildcards.taxid}
        samtools merge -f {output.bam} {input.bams}
        samtools index {output.bam}
        """


##############################################################################
## 5. MEGAHIT CO-ASSEMBLY + DOWNSTREAM ANNOTATION
##############################################################################

rule megahit_coassembly:
    input:
        r1=lambda wc: expand(
            "fastq_idmapped/{taxid}/{sample}_R1.fastq.gz",
            taxid=wc.taxid,
            sample=SAMPLES
        ),
        r2=lambda wc: (
            expand(
                "fastq_idmapped/{taxid}/{sample}_R2.fastq.gz",
                taxid=wc.taxid,
                sample=SAMPLES
            )
            if LAYOUT == "PAIRED"
            else []
        )
    output:
        assembly=(
            "Megahit/{taxid}/combined_assembly/final.contigs.fa"
        ),
        split_dir=directory(
            "Megahit/{taxid}/combined_assembly/split_fasta"
        )
    params:
        layout=LAYOUT,
        outdir="Megahit/{taxid}/combined_assembly",
        r1_csv=lambda wc, input: ",".join(map(str, input.r1)),
        r2_csv=lambda wc, input: (
            ",".join(map(str, input.r2))
            if input.r2
            else ""
        )
    threads: 8
    resources:
        mem_mb=64000
    log:
        "logs/megahit/{taxid}/combined_assembly.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})"

        # MEGAHIT requires a non-existing output directory.
        rm -rf {params.outdir:q}

        if [[ "{params.layout}" == "PAIRED" ]]; then
            megahit \
                -1 {params.r1_csv:q} \
                -2 {params.r2_csv:q} \
                -o {params.outdir:q} \
                -t {threads} \
                --presets meta-sensitive \
                > {log:q} 2>&1
        else
            megahit \
                -r {params.r1_csv:q} \
                -o {params.outdir:q} \
                -t {threads} \
                --presets meta-sensitive \
                > {log:q} 2>&1
        fi

        mkdir -p {output.split_dir:q}

        seqkit split \
            {output.assembly:q} \
            --by-part 20 \
            --out-dir {output.split_dir:q} \
            >> {log:q} 2>&1
        """

rule blastn:
    input:
        fasta="Megahit/{sample}_assembly/split_fasta/{contig}.fa",
    output:
        out="blastn/{sample}/{contig}.out",
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames",
    shell:
        """
        mkdir -p blastn/{wildcards.sample}
        blastn -query {input.fasta} -db refseq_rna -out {output.out} -outfmt "{params.outfmt}"
        """

rule diamond_blastx:
    input:
        fasta="results/trinity_coassembly/150191/trinity_out/Trinity.fasta"
    output:
        tsv="results/diamond_megahit_nr/{taxid}_assembly_diamond.tsv"
    threads: 12
    log:
        "logs/diamond/{taxid}_diamond.log"
    params:
        db="/home/reference_data/bioinfotree/task/blast/nr/nr.gz", 
        outfmt="100"
    shell:
        """
        diamond blastx \
            --db {params.db} \
            --query {input.fasta} \
            --out {output.tsv} \
            --outfmt {params.outfmt} \
            --threads {threads} \
            --max-target-seqs 5 \
            --evalue 1e-5 > {log} 2>&1
        """

rule run_transdecoder:
    input:
        fasta="results/trinity_coassembly/{taxid}/trinity_out/Trinity.fasta"
    output:
        pep="results/transdecoder/{taxid}/Trinity.fasta.transdecoder.pep",
        bed="results/transdecoder/{taxid}/Trinity.fasta.transdecoder.bed"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.pep)
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        # 1. Estrazione delle ORF lunghe (almeno 100 amminoacidi di default)
        TransDecoder.LongOrfs -t ../../../{input.fasta}
        
        # 2. Predizione delle ORF più probabili
        TransDecoder.Predict -t ../../../{input.fasta}
        """


##############################################################################
## 6. ALTERNATIVE PROFILERS  (MetaPhlAn4 / Kaiju + Kaiju-Kraken merge)
##############################################################################

rule metaphlan4:
    input:
        "fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    output:
        profile="metaphlan4/{sample}.profile.txt",
        bowtie2out="metaphlan4/{sample}.bowtie2.bz2",
    params:
        db= config["metaphlan_db"],
        idx = config["metaphlan_index"]
    threads: 4
    shell:
        """
        mkdir -p metaphlan4
        metaphlan {input} --input_type fastq --bowtie2db {params.db} --index {params.idx} --bowtie2out {output.bowtie2out} --nproc {threads}  -o {output.profile} --offline
        """

rule kaiju:
    input:
        "Results/star/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        report="kaiju/{sample}.kaiju.out",
    log:
        "log/{sample}_kaiju.log",
    params:
        threads=25,
        db_Nodes=config.get("kaiju", {}).get("db_Nodes", ""),
        db_FMI=config.get("kaiju", {}).get("db_FMI", ""),
    shell:
        """
        kaiju -z {params.threads} -t {params.db_Nodes} -f {params.db_FMI} -i {input} -o {output.report}  2> {log}
        """

rule kaiju_multi:
    input:
        expand("Results/star/unmapped/{sample}_unmapped_R1.fastq.gz", sample=SAMPLES),
    output:
        reports=expand("kaiju/{sample}.kaiju.out", sample=SAMPLES),
    log:
        "log/kaiju_multi.log",
    params:
        threads=config.get("kaiju", {}).get("threads", 12),
        db_Nodes=config.get("kaiju", {}).get("db_Nodes", ""),
        db_FMI=config.get("kaiju", {}).get("db_FMI", ""),
        inputs=lambda wc, input: ",".join(input),
        outputs=lambda wc, output: ",".join(output.reports),
    shell:
        """
        kaiju-multi -z {params.threads} \
            -t {params.db_Nodes} \
            -f {params.db_FMI} \
            -i {params.inputs} \
            -o {params.outputs} \
            2> {log}
        """

rule kaiju_report:
    input:
        "kaiju/{sample}.kaiju.out",
    output:
        "kaiju/report/{sample}.kaiju.tsv",
    params:
        threads=12,
        db_Nodes=config.get("kaiju", {}).get("db_Nodes", ""),
        db_FMI=config.get("kaiju", {}).get("db_FMI", ""),
        db_Names=config.get("kaiju", {}).get("db_Names", ""),
        taxa_level="species",
    shell:
        """
        mkdir -p kaiju/report
        kaiju2table -t {params.db_Nodes} -n {params.db_Names} -r {params.taxa_level} -o {output} {input} -u -p
        """

rule all_krakaiju_filtered:
    input:
        expand("kaiju/filtered/{sample}.kaiju.out", sample=SAMPLES),

rule filter_krakaiju_output:
    input:
        "kaiju/merged/{sample}.merged.out",
    output:
        "kaiju/filtered/{sample}.kaiju.out",
    params:
        contaminants=config.get("CONTAMINANT_INFO", {}).get("contaminant"),
        human=config.get("CONTAMINANT_INFO", {}).get("humanID"),
    run:
        contaminant_conditions = " && ".join(
            [f"$3 != {c}" for c in params.contaminants]
        )
        human_conditions = " && ".join([f"$3 != {c}" for c in params.human])
        shell(
            f"mkdir -p kaiju/filtered && awk '{contaminant_conditions} && {human_conditions}' {{input}} > {{output}}"
        )

rule all_kaiju_kraken_merge:
    input:
        expand("kaiju/merged/{sample}.merged.out", sample=SAMPLES),

rule kaiju_kraken_merge:
    input:
        kaiju="kaiju/{sample}.kaiju.out",
        kraken="koutput_filtered/{sample}.kraken2",
    output:
        "kaiju/merged/{sample}.merged.out",
    params:
        conflict=config.get("kaiju_merge", "lca"),
        nodes=config.get("kaiju", {}).get("db_Nodes", ""),
    shell:
        r"""
        mkdir -p kaiju/merged

        if [ "{params.conflict}" = "lca" ] || [ "{params.conflict}" = "lowest" ]; then
            TAXOPT="-t {params.nodes}"
        else
            TAXOPT=""
        fi

        kaiju-mergeOutputs \
            -i <(sort -k2,2 {input.kaiju}) \
            -j <(sort -k2,2 {input.kraken}) \
            -o {output} \
            -c {params.conflict} \
            $TAXOPT \
            -v
        """

rule kaiju_kraken_merge_report:
    input:
        merged="kaiju/filtered/{sample}.kaiju.out",
        db=config["kraken_k2d"],
    output:
        report="kaiju_kraken_merged/{sample}.k2report",
    shell:
        """
        k2report {input.db} {input.merged} {output.report}
        """


##############################################################################
## 7. QC & DIAGNOSTICS  (FastQC on depleted reads -> overrepresented BLAST)
##############################################################################

rule all_unmapped_fastqc:
    input:
        expand("fastqc/fastq_taxid_depleted/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("fastqc/fastq_taxid_depleted/{sample}_R1_fastqc.html", sample=SAMPLES),

rule unmapped_fastqc:
    input:
        fastq="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    output:
        zip_out="fastqc/fastq_taxid_depleted/{sample}_R1_fastqc.zip",
        html_out="fastqc/fastq_taxid_depleted/{sample}_R1_fastqc.html",
    threads: 8
    shell:
        """
        mkdir -p fastqc/fastq_taxid_depleted/
        fastqc --threads {threads} --quiet --outdir fastqc/fastq_taxid_depleted/ {input.fastq}
        """

rule extract_overrepresented:
    input:
        zips=expand("fastqc/fastq_taxid_depleted/{sample}_R1_fastqc.zip", sample=SAMPLES),
    output:
        fasta="diagnostics/overrepresented_sequences.fasta",
    shell:
        """
        tmp_fastqc=$(mktemp)
        tmp_block=$(mktemp)
        trap 'rm -f "$tmp_fastqc" "$tmp_block"' EXIT

        # Initialize an empty FASTA file
        > {output.fasta}
        
        echo "Extracting overrepresented sequences across the cohort..."
        
        for zip_file in {input.zips}; do
            # Extract sample name safely
            sample=$(basename "$zip_file" _R1_fastqc.zip)
            
            # Stream the internal fastqc_data.txt without extracting the whole folder to disk
            if unzip -p "$zip_file" "*/fastqc_data.txt" > "$tmp_fastqc" 2>/dev/null; then
                
                # Isolate the exact block containing the overrepresented sequences
                sed -n '/>>Overrepresented sequences/,/>>END_MODULE/p' "$tmp_fastqc" | \
                grep -v ">>" | grep -v "^#Sequence" > "$tmp_block" || true
                
                if [ -s "$tmp_block" ]; then
                    count=1
                    while read -r line; do
                        seq=$(echo "$line" | cut -f1)
                        pct=$(echo "$line" | cut -f3)
                        source=$(echo "$line" | cut -f4)
                        
                        # Write the FASTA header and sequence
                        echo ">${{sample}}_ovr_${{count}} | pct:${{pct}}% | source:${{source}}" >> {output.fasta}
                        echo "$seq" >> {output.fasta}
                        
                        count=$((count + 1))
                    done < "$tmp_block"
                fi
            fi
        done
        
        : # temp files auto-cleaned via trap
        echo "Extraction complete."
        """

rule blast_overrepresented_t2t:
    input:
        fasta="diagnostics/overrepresented_sequences.fasta",
    output:
        out="diagnostics/blastn_T2T_overrepresented.tsv",
    threads: 8
    params: 
        db = config["blast_db"]
    shell:
        """

        # Safety check: ensure the FASTA is not empty before launching BLAST
        if [ -s {input.fasta} ]; then
            blastn -task blastn \
                -query {input.fasta} \
                -db {params.db} \
                -outfmt "6 qseqid sseqid staxids sscinames pident length evalue bitscore stitle" \
                -max_target_seqs 1 \
                -evalue 10 \
                -word_size 11 \
                -dust no \
                -perc_identity 70 \
                -num_threads {threads} \
                -out {output.out}
        else
            touch {output.out}
        fi
        """

rule blast_overrepresented_nt:
    input:
        fasta="diagnostics/overrepresented_sequences.fasta",
    output:
        out="diagnostics/blastn_nt_overrepresented.tsv",
    shell:
        """
        # Safety check: ensure the FASTA is not empty
        if [ -s {input.fasta} ]; then
            # Using -remote to hit the complete NCBI database
            blastn -task megablast \
                -query {input.fasta} \
                -db nt \
                -remote \
                -outfmt "6 qseqid sseqid staxids sscinames pident length evalue bitscore stitle" \
                -max_target_seqs 1 \
                -out {output.out}
        else
            touch {output.out}
        fi
        """


##############################################################################
## 8. PROJECT-SPECIFIC  (Bowtie2 / S. aureus region hunt; paracoccus targets)
##    TODO: hardcoded refs & filenames -- slated to move to a dataset include
##############################################################################

rule bowtie2:
    input:
        fastq="fastq_unclassified/{sample}_R1.fa",
        index="Resources/S_aureus/S_aureus_index.1.bt2",
    output:
        sam="bowtie2/{sample}.sam",
    params:
        prefix="Resources/S_aureus/S_aureus_index",
    threads: 8
    shell:
        """
        bowtie2 -x {params.prefix} -U {input.fastq} -S {output.sam} -p {threads} -f --no-unal
        """

rule all_samtools_bam:
    input:
        expand("bowtie2/{sample}.bam", sample=SAMPLES),

rule samtools_bam:
    input:
        "bowtie2/{sample}.sam",
    output:
        "bowtie2/{sample}.bam",
    shell:
        """
        samtools view -bS {input} | samtools sort -o {output}
        """

rule all_samtools_index:
    input:
        expand("bowtie2/{sample}.bam.bai", sample=SAMPLES),

rule samtools_index:
    input:
        "bowtie2/{sample}.bam",
    output:
        "bowtie2/{sample}.bam.bai",
    shell:
        """
        samtools index {input}
        """

rule all_merged_beds:
    input:
        expand("bed_regions/{sample}_merged.bed", sample=SAMPLES),

rule bam_to_merged_bed:
    input:
        bam="bowtie2/{sample}.bam",
    output:
        bed="bed_regions/{sample}_merged.bed",
    shell:
        """
        bedtools bamtobed -i {input.bam} | bedtools merge -d 1000 > {output.bed}
        """

rule combine_all_beds:
    input:
        expand("bed_regions/{sample}_merged.bed", sample=SAMPLES),
    output:
        "bed_regions/tutte_le_regioni_sospette.bed",
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        """

rule bam_2_fasta:
    input:
        "bowtie2/{sample}.bam",
    output:
        "bowtie2/fasta/{sample}.fa",
    shell:
        """
        mkdir -p bowtie2/fasta
        samtools fasta -F 4 {input} > {output}
        """

rule blast:
    input:
        fasta="fastq_idmapped/{sample}_R1.fa",
    output:
        out="blastn_T2T/{sample}.blastn.out",
    threads: 8
    params: 
        db = config["blast_db"]
    shell:
        """
        mkdir -p blastn_T2T

        if [ -s {input.fasta} ]; then
            blastn -task blastn \
                -query {input.fasta} \
                -db {params.db} \
                -outfmt "6 qseqid sseqid staxids sscinames pident length evalue bitscore" \
                -evalue 10 \
                -word_size 11 \
                -dust no \
                -perc_identity 70 \
                -num_threads {threads} \
                -out {output.out}
        else
            touch {output.out}
        fi
        """

rule samtools_contigs_2_fasta: 
    input:
        contigs = "contigs_to_blast_paracoccus.txt", 
        ref = "Resources/genomes/150191/GCA_030161055.1/GCA_030161055.1_ASM3016105v1_genomic.fna"
    output: 
        fasta = "targets_to_blast_paracoccus_ref.fasta"
    shell:
        """
        samtools faidx {input.ref} $(cat {input.contigs} | tr '\n' ' ') > {output.fasta}
        """
