##########################
### Rules for Kraken2  ###
##########################

if config["LAYOUT"] == "PAIRED":
    ruleorder: kraken_paired_ends > kraken_single_end
else:
    ruleorder: kraken_single_end > kraken_paired_ends


def kraken_input_r1(wc):
    if HUMAN_FILTER_ENABLED:
        return f"fastq/unmapped_human_filtered/{wc.sample}_unmapped_R1.fastq.gz"
    return f"fastq/unmapped/{wc.sample}_unmapped_R1.fastq.gz"


def kraken_input_r2(wc):
    if HUMAN_FILTER_ENABLED:
        return f"fastq/unmapped_human_filtered/{wc.sample}_unmapped_R2.fastq.gz"
    return f"fastq/unmapped/{wc.sample}_unmapped_R2.fastq.gz"


rule kraken_paired_ends:
    input:
        R1=kraken_input_r1,
        R2=kraken_input_r2
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2",
        unclassified1="unclassified/{sample}_unclassified_1.fq",
        unclassified2="unclassified/{sample}_unclassified_2.fq"
    threads: 4
    shell:
        r"""
        k2 classify \
          --db {config[kraken_db]} \
          --threads {threads} \
          --report {output.report} \
          --output {output.out} \
          --paired {input.R1} {input.R2} \
          --unclassified-out unclassified/{wildcards.sample}_unclassified#.fq \
          --minimum-hit-groups {config[hitgrp]}
        """


rule kraken_single_end:
    input:
        R1=lambda wc: kraken_input_r1(wc)
    output:
        report="kreports/{sample}.k2report",
        out="koutputs/{sample}.kraken2"
    threads: 4
    shell:
        r"""
        k2 classify \
          --db {config[kraken_db]} \
          --minimum-hit-groups {config[hitgrp]} \
          --threads {threads} \
          --report {output.report} \
          {input.R1} \
          > {output.out}
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

rule bracken:
    input:
        "kreports/{sample}.k2report"
    output:
        report="breports/{sample}.breport",
        out="boutputs/{sample}.bracken"
    shell:
        r"""
        bracken \
          -d {config[kraken_db]} \
          -i {input} \
          -r {config[BRACKEN][bracken_read_len]} \
          -l {config[BRACKEN][bracken_level]} \
          -t {config[BRACKEN][bracken_min_reads]} \
          -o {output.out} \
          -w {output.report}
        """


rule bracken_merged:
    input:
        outputs=expand("boutputs_filtered/{sample}.bracken", sample=SAMPLES)
    output:
        "bracken_merged_abundances.tsv"
    shell:
        r"""
        combine_bracken_outputs.py --files {input.outputs} -o {output} 2> log.txt
        """

#######################################
### Rules for krona and other rules ###
#######################################
rule krona_txt:
    input:
        "breports_filtered/{sample}.breport"
    output:
        "b_krona_txt/{sample}.b.krona.txt"
    shell:"""
        kreport2krona.py -r {input} -o {output} --no-intermediate-ranks
    """

rule krona_html:
    input: "b_krona_txt/{sample}.b.krona.txt"
    output: "krona_html/{sample}.krona.html"
    shell: "ktImportText {input} -o {output}"
    
rule spit_merged:
    input:
        "bracken_merged_abundances.txt"
    output:
        num="bracken_merged_abundances.num.txt",
        frac="bracken_merged_abundances.frac.txt"
    shell:
        "grep_columns -k 1,2,3 bracken_num  < {input} | perl -pe '$.==1; s/.bracken_num//g'  > {output.num};"
        "grep_columns -k 1,2,3 bracken_frac < {input} | perl -pe '$.==1; s/.bracken_frac//g' > {output.frac}"
    
rule feature_filter:
    input:
        abundances="bracken_merged_abundances.num.txt",
        metadata="metadata.txt"
    output:
        filtered="bracken_merged_abundances.num.filtered.txt"
    params:
        condition=config['feature_filter']['condition'],
        g1=config['feature_filter']["g1"],
        g2=config['feature_filter']["g2"],
        use_raw_counts=config['feature_filter']["use_raw_counts"],
        min_exp=config['feature_filter']['min_exp'],
        min_samples_ratio=config['feature_filter']['min_samples_ratio']
    shell: """
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
        num_filter="bracken_merged_abundances.num.filtered.txt",
        frac="bracken_merged_abundances.frac.txt"
    output:
        "bracken_merged_abundances.frac.filtered.txt"
    shell:
        "filter_1col 3 <(cut -f 3 {input.num_filter}) < {input.frac} > {output}"

rule collapse_taxid:
    input: "bracken_merged_abundances.num.filtered.txt"
    output: "bracken_merged_abundances.num.filtered.taxid_collapsed.txt"
    shell: "perl -pe 's/\t/;/; s/\t/;/' {input} > {output}"

rule degw:
    input: 
        abundances="bracken_merged_abundances.num.filtered.taxid_collapsed.txt",
        metadata="metadata.txt"
    output:
        "bracken_merged_abundances.num.filtered.taxid_collapsed.degw.txt",
    params:
        condition=config['feature_filter']['condition'],
        g1=config['feature_filter']["g1"],
        g2=config['feature_filter']["g2"],
        use_raw_counts=config['feature_filter']["use_raw_counts"],
        min_exp=config['feature_filter']['min_exp'],
        min_samples_ratio=config['feature_filter']['min_samples_ratio']
    shell:"""
        DEGWilcox.R {input.abundances} {input.metadata} \
            --condition={params.condition} --g1 {params.g1} --g2 {params.g2} \
            --use_raw_counts --min_exp {params.min_exp} --min_samples_ratio {params.min_samples_ratio} \
        | bawk 'NR==1 {{$1="name\ttaxonomy_id\tlevel"; print}} NR>1{{gsub(/;/, "\t", $1); print}}' > {output}
    """

rule extract_unclassified_id_paired:
    input: "koutput_filtered/{sample}.kraken2"
    output: "fastq_unclassified/{sample}.id"
    shell: """"
        awk '{{print $2,$1}}' {input} | collapsesets 2 | bawk '$2=="U"' > {output}
    """

rule extract_kraken_unclassified_reads:
    input:
        kraken2_output = "koutput_filtered/{sample}.kraken2",
        kraken2_report = "kreports_filtered/{sample}.k2report",
        fastq_r1 = "fastq/{sample}_R1.fastq.gz",
        fastq_r2 = "fastq/{sample}_R2.fastq.gz"
    output:
        fastq_r1_unclassified = "fastq_unclassified/{sample}_R1.fastq.gz",
        fastq_r2_unclassified = "fastq_unclassified/{sample}_R2.fastq.gz"
    params:
        taxid = 1,
        # Reference the config file here
        exclude_opt = "--exclude" if config["kraken"]["exclude_classified"] else "",
        children_opt = "--include-children" if config["kraken"]["include_children"] else ""
    shell:
        """
        extract_kraken_reads.py \
            -k {input.kraken2_output} \
            --taxid {params.taxid} \
            {params.exclude_opt} \
            {params.children_opt} \
            -s1 {input.fastq_r1} \
            -s2 {input.fastq_r2} \
            --report {input.kraken2_report} \
            -o >(gzip > {output.fastq_r1_unclassified}) \
            -o2 >(gzip > {output.fastq_r2_unclassified})
        """

###########################
### Remove contaminant ####
###########################

rule filter_kraken_output:
    input:
        "koutputs/{sample}.kraken2"
    output:
        "koutput_filtered/{sample}.kraken2"
    params:
        contaminants = config.get("CONTAMINANT_INFO", {}).get("contaminant"),  
        human = config.get("CONTAMINANT_INFO", {}).get("humanID")
    run:
        contaminant_conditions = " && ".join([f"$3 != {c}" for c in params.contaminants])
        human_conditions = " && ".join([f"$3 != {c}" for c in params.human])
        shell(f"mkdir -p koutput_filtered && awk '{contaminant_conditions} && {human_conditions}' {{input}} > {{output}}")


rule filt_k2report:
    input:
        kraken2_filtered = "koutput_filtered/{sample}.kraken2", 
        db = config["kraken_k2d"]
    output:
        report = "kreports_filtered/{sample}.k2report"
    shell:
        """
        mkdir -p kreports_filtered
        /home/molinerislab/NeriMetagenome/workflow/kraken2/src/k2report {input.db} {input.kraken2_filtered} {output.report}
        """

rule filbracken:
    input:
        "kreports_filtered/{sample}.k2report"  
    output:
        report="breports_filtered/{sample}.breport",
        out="boutputs_filtered/{sample}.bracken"
    shell:"""
        mkdir -p breports_filtered boutputs_filtered
        bracken -d {config[kraken_db]} -i {input} -r {config[BRACKEN][bracken_read_len]} -l {config[BRACKEN][bracken_level]} -t {config[BRACKEN][bracken_min_reads]} -o {output.out} -w {output.report}
    """

######################### 
### MEGAHIT assembly ####
#########################

rule megahit_assembly:
    input:
        "fastq/unmapped/{sample}_unmapped.fastq.gz"
    output:
        assembly="Megahit/{sample}_assembly/final.contigs.fa",
        split_dir=directory("Megahit/{sample}_assembly/split_fasta")
    params:
        outdir="Megahit/{sample}_assembly",
    threads: 8
    shell:
        """
        rm -rf {params.outdir}
        megahit -r {input}  -o {params.outdir} -t {threads} --preset meta-sensitive --keep-tmp
        mkdir -p {output.split_dir}
        seqkit split {output.assembly} -p 20 -O {output.split_dir}
        """

#################################
### BLASTN against RefSeq RNA ###
#################################

rule blastn:
    input:
        fasta="Megahit_meta/{sample}_assembly/split_fasta/{contig}.fa"
    output:
        out="blastn/{sample}/{contig}.out"
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames"
    shell:
        """
        blastn -query {input.fasta} -db refseq_rna -out {output.out} -outfmt "{params.outfmt}"
        """


######################
### Structural QC ###
#####################

rule bracken_sample_totals:
    input:
        expand("boutputs_filtered/{sample}.bracken", sample=SAMPLES)
    output:
        "qc/sample_totals.txt"
    shell:
        r"""
        echo "Sample Total_Reads" > {output}
        for f in {input}; do
            sample=$(basename "$f" .bracken)
            total=$(awk -F'\t' '
                NR==1 {{
                    for(i=1;i<=NF;i++) {{
                        if($i=="new_est_reads") col=i
                    }}
                }}
                NR>1 {{
                    sum += $col
                }}
                END {{
                    print sum
                }}
            ' "$f")
            echo "$sample $total" >> {output}
        done
        """

rule qc_abundance_matrix:
    input:
        matrix="bracken_merged_abundances.tsv",
        sample_totals="qc/sample_totals.txt"
    output:
        checked_num_matrix="qc/abundance_num.tsv",
        checked_frac_matrix="qc/abundance_frac.tsv",
        sample_name_map="qc/sample_name_verification.tsv",
        missing_summary="qc/missing_value_summary.tsv",
        sample_summary="qc/sample_qc_summary.tsv",
        taxa_summary="qc/taxa_qc_summary.tsv",
        notes="qc/filtering_notes.txt"
    params:
        min_sample_presence=2,
        min_total_abundance=10
    script:
        "../../local/src/qc_abundance_matrix.py"


rule filter_taxa:
    input:
        matrix="bracken_merged_abundances.tsv"
    output:
        filtered="bracken_merged_abundances.filtered.tsv",
        summary="filter_summary.txt"
    params:
        min_prevalence=config["METAGENOMICS"]["TAXA_FILTERING"]["min_prevalence"],
        min_total=config["METAGENOMICS"]["TAXA_FILTERING"]["min_total"]
    log:
        "logs/filter_taxa.log"
    script:
        "../../local/src/filter_taxa.py"


rule normalize_abundance:
    input:
        matrix="bracken_merged_abundances.filtered.tsv"
    output:
        relative="abundances.filtered.relative.tsv",
        clr="abundances.filtered.clr.tsv"
    log:
        "logs/normalize_abundance.log"
    script:
        "../../local/src/normalize_abundance.py"
        