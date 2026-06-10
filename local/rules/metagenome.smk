##########################
### Rules for Kraken2  ###
##########################
if config["LAYOUT"] == "PAIRED":

    ruleorder: kraken_pe_pass1 > kraken_se_pass1

else:

    ruleorder: kraken_se_pass1 > kraken_pe_pass1


rule kraken_pe_pass1:
    input:
        R1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
        R2="fastq/unmapped/{sample}_unmapped_R2.fastq.gz",
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2",
        unclassified1="unclassified/{sample}_unclassified_1.fq",
        unclassified2="unclassified/{sample}_unclassified_2.fq",
    threads: 6
    shell:
        """
        kraken2 --db {config[kraken_db]} \
            --threads {threads} \
            --report {output.report} \
            --output {output.out} \
            --paired {input.R1} {input.R2} \
            --unclassified-out unclassified/{wildcards.sample}_unclassified_#.fq
    """


rule kraken_nr_paired_ends:
    input:
        R1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
        R2="fastq/unmapped/{sample}_unmapped_R2.fastq.gz",
    output:
        report="kreports_nr/{sample}.k2report",
        out="koutputs_nr/{sample}.kraken2",
    threads: 6
    shell:
        """
        k2 --db {config[kraken_db_nr]} {config[kraken_options]} \
            --threads {threads} \
            --report {output.report}\
            --output {output.out} \
            --paired {input.R1} {input.R2}
    """


rule kraken_se_pass1:
    input:
        R1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2",
    threads: 6
    shell:
        """
        kraken2 --db {config[kraken_db]} {config[kraken_options]} \
            --threads {threads} \
            --report-minimizer-data \
            --memory-mapping \
            --report {output.report} \
            --output {output.out} \
            {input.R1}
        """


rule remove_taxid_reads:
    input:
        kraken2_output="koutputs/{sample}.kraken2",
        kraken2_report="kreports/{sample}.k2report",
        fastq_r1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        fastq_r1_nonhost="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    params:
        taxids=lambda wc: " ".join(
            [
                str(t)
                for t in config.get("CONTAMINANT_INFO", {}).get("humanID", ["9606"])
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


rule kraken_se_pass2:
    input:
        R1="fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    output:
        report="kreports_filtered/{sample}.k2report",
        out="koutput_filtered/{sample}.kraken2",
    threads: 6
    shell:
        """
        k2 --db {config[kraken_db]} {config[kraken_options]} \
            --threads {threads} \
            --report-minimizer-data \
            --memory-mapping \
            --report {output.report} \
            --output {output.out} \
            {input.R1}
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


#########################
### Rules for Bracken ###
#########################
rule braken:
    input:
        "kreports/{sample}.k2report",
    output:
        report="breports/{sample}.breport",
        out="boutputs/{sample}.braken",
    shell:
        """
        bracken -d {config[kraken_db]} -i {input} -r {config[BRACKEN][braken_read_len]} -l {config[BRACKEN][braken_level]} -t {config[BRACKEN][braken_min_reads]} -o {output.out} -w {output.report}
        """


rule bracken_merged:
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


rule filbraken:
    input:
        "kreports_filtered/{sample}.k2report",
    output:
        report="breports_filtered/{sample}.breport",
        out="boutputs_filtered/{sample}.braken",
    shell:
        """
        mkdir -p breports_filtered boutputs_filtered
        bracken -d {config[kraken_db]} -i {input} -r {config[BRACKEN][braken_read_len]} -l {config[BRACKEN][braken_level]} -t {config[BRACKEN][braken_min_reads]} -o {output.out} -w {output.report}
        """


########################
### Rules for krona  ###
########################


rule krona_txt:
    input:
        "breports_filtered/{sample}.breport",
    output:
        "b_krona_txt/{sample}.b.krona.txt",
    shell:
        """
        kreport2krona.py -r {input} -o {output} --no-intermediate-ranks
    """


rule krona_html:
    input:
        "b_krona_txt/{sample}.b.krona.txt",
    output:
        "krona_html/{sample}.krona.html",
    shell:
        """
        ktImportText {input} -o {output}
        """


rule spit_merged:
    input:
        "bracken_merged_abbundances.txt",
    output:
        num="bracken_merged_abbundances.num.txt",
        frac="bracken_merged_abbundances.frac.txt",
    shell:
        """
        "grep_columns -k 1,2,3 braken_num  < {input} | perl -pe '$.==1; s/.braken_num//g'  > {output.num};"
        "grep_columns -k 1,2,3 braken_frac < {input} | perl -pe '$.==1; s/.braken_frac//g' > {output.frac}"
        """


rule feature_filter:
    input:
        abundances="bracken_merged_abbundances.num.txt",
        metadata="metadata.txt",
    output:
        filtered="bracken_merged_abbundances.num.filtered.txt",
    params:
        condition=config["feature_filter"]["condition"],
        g1=config["feature_filter"]["g1"],
        g2=config["feature_filter"]["g2"],
        use_raw_counts=config["feature_filter"]["use_raw_counts"],
        min_exp=config["feature_filter"]["min_exp"],
        min_samples_ratio=config["feature_filter"]["min_samples_ratio"],
    shell:
        """
        perl -pe 's/\t/;/; s/\t/;/' {input.abundances} > {input.abundances}.tmp;\
        echo -en "name\\ttaxonomy_id\\ttaxonomy_lv\\t" > {output.filtered}
        feature_filter {input.abundances}.tmp {input.metadata} \
            --condition={params.condition} --g1 {params.g1} --g2 {params.g2} \
            --use_raw_counts --min_exp {params.min_exp} --min_samples_ratio {params.min_samples_ratio} \
        | perl -pe 's/;/\t/; s/;/\t/' >> {output.filtered};\
        rm {input.abundances}.tmp
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
            --use_raw_counts --min_exp {params.min_exp} --min_samples_ratio {params.min_samples_ratio} \
        | bawk 'NR==1 {{$1="name\ttaxonomy_id\tlevel"; print}} NR>1{{gsub(/;/, "\t", $1); print}}' > {output}
    """


rule extract_unclassified_id_paired:
    input:
        "koutput_filtered/{sample}.kraken2",
    output:
        "fastq_unclassified/{sample}.id",
    shell:
        """
        awk '{{print $2,$1}}' {input} | collapsesets 2 | bawk '$2=="U"' > {output}
    """


rule extract_kraken_reads:
    input:
        kraken2_output="koutputs/{sample}.kraken2",
        kraken2_report="kreports/{sample}.k2report",
        fastq_r1="fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        fastq_r1_unclassified="fastq_idmapped/{sample}_R1.fastq.gz",
    params:
        taxid=config.get("kraken_extract_taxid"),
    shell:
        """
        mkdir -p fastq_idmapped
        extract_kraken_reads.py \
            -k {input.kraken2_output} \
            --taxid {params.taxid} \
            --include-children \
            -s {input.fastq_r1} \
            --report {input.kraken2_report} \
            --fastq-output \
            -o >(gzip > {output.fastq_r1_unclassified}) 
        """


#########################
### MEGAHIT assembly ####
#########################


rule megahit_assembly:
    input:
        "fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
    output:
        assembly="Megahit/{sample}_assembly/final.contigs.fa",
        split_dir=directory("Megahit/{sample}_assembly/split_fasta"),
    params:
        outdir="Megahit/{sample}_assembly",
    threads: 8
    shell:
        """
        megahit -r {input}  -o {params.outdir} -t {threads} --preset meta-sensitive --keep-tmp
        mkdir -p {output.split_dir}
        seqkit split {output.assembly} -p 20 -O {output.split_dir}
        """


#################################
### BLASTN against RefSeq RNA ###
#################################


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
    shell:
        """
        mkdir -p blastn_T2T
        export BLASTDB=/mnt/nobackup/home/reference_data/bioinfotree/task/blast/T2T-CHM13

        if [ -s {input.fasta} ]; then
            blastn -task blastn \
                -query {input.fasta} \
                -db T2T-CHM13 \
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


##################
### Metaphlan4 ###
##################


rule metaphlan4:
    input:
        "fastq/fastq_taxid_depleted/{sample}_R1.fastq.gz",
    output:
        profile="metaphlan4/{sample}.profile.txt",
        bowtie2out="metaphlan4/{sample}.bowtie2.bz2",
    params:
        db="/home/reference_data/bioinfotree/task/metaphlan/dataset/nobackup/mpa_vJan25_CHOCOPhlAnSGB_202503",
    threads: 4
    shell:
        """
        mkdir -p metaphlan4
        metaphlan {input} --input_type fastq --bowtie2db {params.db} --index mpa_vJan25_CHOCOPhlAnSGB_202503 --bowtie2out {output.bowtie2out} --nproc {threads}  -o {output.profile} --offline
        """


##########Metagenomic Analysis and Plot#######################


rule common_taxa_in_samples:
    input:
        expand("boutputs_filtered/{sample}.braken", sample=SAMPLES),
    output:
        "results/common_taxa.csv",
        directory("results"),
    params:
        samples=SAMPLES,
    shell:
        """
        mkdir -p results 
        Rscript ../../local/src/common_species.R
        """


rule reads_human_contam_classified:
    output:
        "classified_vs_human_contaminant_barplot_normalized.png",
        "classified_vs_human_contaminant_barplot_percentage.png",
    shell:
        """
        Rscript /home/molinerislab/NeriMetagenome/workflow/src/plot_reads.R
        """


rule kaiju:
    input:
        "fastq/unmapped/{sample}_unmapped_R1.fastq.gz",
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
        expand("fastq/unmapped/{sample}_unmapped_R1.fastq.gz", sample=SAMPLES),
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
        kaiju2table -t {params.db_Nodes} -n {params.db_Names} -r {params.taxa_level} -o {output} {input} -r species -u -p 
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


rule all_kkreport:
    input:
        expand("kaiju_kraken_merged/{sample}.k2report", sample=SAMPLES),


rule kaiju_kraken_merge_report:
    input:
        merged="kaiju/filtered/{sample}.kaiju.out",
        db=config["kraken_k2d"],
    output:
        report="kaiju_kraken_merged/{sample}.k2report",
    shell:
        """
        /home/molinerislab/NeriMetagenome/workflow/kraken2/src/k2report {input.db} {input.merged} {output.report}
        """


rule all_bowtie2:
    input:
        expand("bowtie2/{sample}.sam", sample=SAMPLES),


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


# =============================================================================
# minimap2 alignment of extracted reads vs generic references
# =============================================================================

MINIMAP2 = "/home/molinerislab/IsellaIsoforms/local/env/conda/bin/minimap2"

# Define available references here. You can add more species as needed.
MINIMAP2_REFS = "Resources/GCF_022494545.1_ASM2249454v1_genomic.fna"
# {
#     "aureus": "Resources/GCF_022494545.1_ASM2249454v1_genomic.fna",
#     "cerus": "Resources/b_cerus/ncbi_dataset/data/GCF_030518615.1/GCF_030518615.1_ASM3051861v1_genomic.fna"
# }


ruleorder: minimap2_merge > minimap2_index


rule all_minimap2:
    """Align FASTA reads against all defined references and merge them."""
    input:
        expand("minimap2_{species}/merged_all_samples.bam.bai", species="aureus"),


rule minimap2_align:
    """Align extracted reads (FASTA) with minimap2 short-read preset to a specific species."""
    input:
        fq="fastq_idmapped/{sample}_R1.fastq.gz",
    output:
        bam="minimap2_{species}/{sample}.bam",
    params:
        ref=MINIMAP2_REFS,
    threads: 8
    shell:
        """
        {MINIMAP2} -ax sr -t {threads} --secondary=no \
            {params.ref} {input.fq} \
            | samtools view -bS -F 4 \
            | samtools sort -o {output.bam}
        """


rule minimap2_index:
    input:
        "minimap2_{species}/{sample}.bam",
    output:
        "minimap2_{species}/{sample}.bam.bai",
    shell:
        "samtools index {input}"


rule minimap2_merge:
    """Merge all per-sample minimap2 BAMs into a single file for IGV/coverage."""
    input:
        bams=expand("minimap2_{{species}}/{sample}.bam", sample=SAMPLES),
        bais=expand("minimap2_{{species}}/{sample}.bam.bai", sample=SAMPLES),
    output:
        bam="minimap2_{species}/merged_all_samples.bam",
        bai="minimap2_{species}/merged_all_samples.bam.bai",
    shell:
        """
        samtools merge -f {output.bam} {input.bams}
        samtools index {output.bam}
        """


rule all_unmapped_fastqc:
    input:
        expand("fastqc/fastq_h_depleted/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("fastqc/fastq_h_depleted/{sample}_R1_fastqc.html", sample=SAMPLES),


rule unmapped_fastqc:
    input:
        fastq="fastq/fastq_h_depleted/{sample}_R1.fastq.gz",
    output:
        zip_out="fastqc/fastq_h_depleted/{sample}_R1_fastqc.zip",
        html_out="fastqc/fastq_h_depleted/{sample}_R1_fastqc.html",
    threads: 8
    shell:
        """
        mkdir -p fastqc/fastq_h_depleted/
        fastqc --threads {threads} --quiet --outdir fastqc/fastq_h_depleted/ {input.fastq}
        """


rule extract_overrepresented:
    input:
        zips=expand("fastqc/fastq_h_depleted/{sample}_R1_fastqc.zip", sample=SAMPLES),
    output:
        fasta="diagnostics/overrepresented_sequences.fasta",
    shell:
        """
        # Initialize an empty FASTA file
        > {output.fasta}
        
        echo "Extracting overrepresented sequences across the cohort..."
        
        for zip_file in {input.zips}; do
            # Extract sample name safely
            sample=$(basename "$zip_file" _R1_fastqc.zip)
            
            # Stream the internal fastqc_data.txt without extracting the whole folder to disk
            if unzip -p "$zip_file" "*/fastqc_data.txt" > temp_fastqc.txt 2>/dev/null; then
                
                # Isolate the exact block containing the overrepresented sequences
                sed -n '/>>Overrepresented sequences/,/>>END_MODULE/p' temp_fastqc.txt | \
                grep -v ">>" | grep -v "^#Sequence" > temp_block.txt || true
                
                if [ -s temp_block.txt ]; then
                    count=1
                    while read -r line; do
                        seq=$(echo "$line" | cut -f1)
                        pct=$(echo "$line" | cut -f3)
                        source=$(echo "$line" | cut -f4)
                        
                        # Write the FASTA header and sequence
                        echo ">${{sample}}_ovr_${{count}} | pct:${{pct}}% | source:${{source}}" >> {output.fasta}
                        echo "$seq" >> {output.fasta}
                        
                        count=$((count + 1))
                    done < temp_block.txt
                fi
            fi
        done
        
        rm -f temp_fastqc.txt temp_block.txt
        echo "Extraction complete."
        """


rule blast_overrepresented_t2t:
    input:
        # We point directly to the single aggregated FASTA we just created
        fasta="diagnostics/overrepresented_sequences.fasta",
    output:
        # A single summary TSV file for the entire cohort
        out="diagnostics/blastn_T2T_overrepresented.tsv",
    threads: 8
    shell:
        """
        export BLASTDB=/mnt/nobackup/home/reference_data/bioinfotree/task/blast/T2T-CHM13

        # Safety check: ensure the FASTA is not empty before launching BLAST
        if [ -s {input.fasta} ]; then
            blastn -task blastn \
                -query {input.fasta} \
                -db T2T-CHM13 \
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
        # Taking the FASTA we extracted from FastQC
        fasta="diagnostics/overrepresented_sequences.fasta",
    output:
        # A new summary TSV for the global hits
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
