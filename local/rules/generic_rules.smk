#%.xlsx: %
#    cat $< | tab2xlsx > $@
rule tab2xlsx:
    input: 
        "{file}"
    output: 
        "{file}.xlsx"
    shell: 
        "cat < {input} | tab2xlsx > {output}"

rule gz2xlsx:
    input: 
        "{file}.gz"
    output: 
        "{file}.xlsx"
    shell: 
        "zcat < {input} | tab2xlsx > {output}"

##TODO add BIGWIG_BIN_SIZE = config["BIGWIG_BIN_SIZE"] in snakefile or config
##TODO add GFF3_MIN_WINDOW_READS, BEDGRAPH_FILTER, FASTQ_FILTERING
BEDGRAPH_FILTER = 1

##NEW##
#%.xls: %
#    zcat $< | tab2xls > $@
rule tab2xls:
    input: 
        "{file}"
    output: 
        "{file}.xls"
    shell: 
        "tab2xls < {input} > {output}"

        
#Non so se tab2xlsx e tab2xls accettino anche gz in input, in caso le regole qui sotto sono gestite sopra    
#%.header_added.xls: %.header_added.gz
#    zcat $< | tab2xls > $@

#%.xls: %.gz
#    zcat $< | tab2xls > $@

#%.header_added.xls: %.header_added
#    cat $< | tab2xls > $@

#Se tab2xlsx e tab2xls non accettano gz in input serve un workaround. Chat GPT suggerisce di fare cosi'
# rule unzip:
#     input:
#         "data/{file}.gz"
#     output:
#         "data/{file}"
#     shell:
#         "gunzip -c {input} > {output}"

# Include a wildcard constraint to match both the original and unzipped filenames
# wildcard_constraints:
#     file="[^.]+",  # Matches any file name without an extension
#     ext="(?:\\.gz)?",  # Matches optional .gz extension

#TODO check the rule and change the name, too generic!
# rule process:
#     input:
#         "data/{file}{ext}"
#     output:
#         xlsx="data/{file}.xlsx",
#         xls="data/{file}.xls"
#     run:
#         if input.ext == ".gz":
#             unzipped_file = f"data/{wildcards.file}"
#             shell("gunzip -c {input} > {unzipped_file}")
#             shell("tab2xlsx < {unzipped_file} > {output.xlsx}")
#             shell("tab2xls < {unzipped_file} > {output.xls}")
#         else:
#             shell("tab2xlsx < {input} > {output.xlsx}")
#             shell("tab2xls < {input} > {output.xls}")

rule fastqc:
    input:
        "{path}.fastq.gz"
    output:
        fastqc_html="fastqc/{path}_fastqc.html",
        fastqc_zip="fastqc/{path}_fastqc.zip"
    threads: config["CORES"]
    container: 
        "docker://quay.io/biocontainers/fastqc:0.11.3--0"
    shell:"""
        mkdir -p $(dirname {output.fastqc_html})
        fastqc -t {threads} -o `dirname {output.fastqc_html}` {input}
    """

# #%.fastq.gz: %.fqz
# #    fqz_comp -d < $< | gzip > $@
# rule fqz2fastqgz:
#     input: 
#         "{file}.fqz"
#     output: 
#         "{file}.fastq.gz"
#     shell: 
#         "fqz_comp -d < {input} | tab2xlsx > {output}"

#%.cram: %.bam
#    samtools view -@ $(CORES) -T $(GENCODE_GENOME_FASTA) -C -o $@ $<
rule bam2cram:
    input: 
        "{file}.bam"
    output: 
        "{file}.cram"
    params:
        threads = config["CORES"],
        genome = GENCODE_GENOME_FASTA
    shell: 
        "samtools view -@ {params.threads} -T {params.genome} -C -o {output} {input}"

# ALL_cram: $(addprefix STAR/$(FASTQ_FILTERING)/,$(addsuffix .STAR/Aligned.sortedByCoord.out.cram, $(SAMPLES)))
#     mkdir -p $(FASTQ_FILTERING); cd $(FASTQ_FILTERING);\
#     for s in $(SAMPLES); do\
#          ln ../STAR/$(FASTQ_FILTERING)/$$s.STAR/Aligned.sortedByCoord.out.cram $$s.cram;\
#     done;
#rule ALL_cram:
#    input:
#        expand("STAR/{filter}/{samples}.STAR/Aligned.sortedByCoord.out.cram", filter=config["FASTQ_FILTERING"], samples = SAMPLES)
#    params:
#        filter=config["FASTQ_FILTERING"]
#    shell:
#        """
#        mkdir -p {params.filter}; cd {params.filter};\
#        for s in {params.filter}; do\
#            ln -s ../STAR/{params.filter}/$s.STAR/Aligned.sortedByCoord.out.cram $s.cram;\
#        done;
#        """

#%.header_added.gz: %.gz
#    (bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
rule gz_header_add:
    input: 
        "{path}.gz"
    output: 
        "{path}.header_added.gz"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; zcat {input} ) | gzip > {output}"
        
#%.header_added: %
#    (bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
rule header_add:
    input: 
        "{file}"
    output: 
        "{file}.header_added"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; cat {input} ) > {output}"

#%.bam.bai: %.bam
#    samtools index $<
rule get_bai:
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.bai"
    shell: 
        "samtools index {input}"

#%.bam.id: %.bam
#    samtools view $< | cut -f 1 | bsort -S8% | uniq > $@
rule get_bam_id:
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.id"
    shell: 
        "samtools view {input} | cut -f 1 | bsort -S8% | uniq > {output}"

#%.bw: %.bam %.bam.bai
    #bamCoverage --binSize=$(BIGWIG_BIN_SIZE) -b $< -o $@ --numberOfProcessors=$(CORES)
rule get_bw:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.bw"
    params:
        threads = config["CORES"],
        bin_size = "BIGWIG_BIN_SIZE"
    shell: 
        "bamCoverage --binSize={params.bin_size} -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

#%.norm.bw: %.bam %.bam.bai
#    bamCoverage --binSize=$(BIGWIG_BIN_SIZE) --normalizeUsing=CPM -b $< -o $@ --numberOfProcessors=$(CORES)
rule get_norm_bw:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.norm.bw"
    params:
        threads = config["CORES"],
        bin_size = "BIGWIG_BIN_SIZE"
    shell: 
        "bamCoverage --binSize={params.bin_size}  --normalizeUsing=CPM -b {input.bam} -o {output} --numberOfProcessors={params.threads}"


#Le due regole sotto sono uguali a quelle sopra solo che aggiungono nell'output la bin_size

#%.$(BIGWIG_BIN_SIZE).bw: %.bam %.bam.bai
#    bamCoverage --binSize=$(BIGWIG_BIN_SIZE) -b $< -o $@ --numberOfProcessors=$(CORES)
rule get_bw_binsize:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.{bin_size}.bw"
    params:
        threads = config["CORES"],
        bin_size = BIGWIG_BIN_SIZE
    shell: 
        "bamCoverage --binSize={params.bin_size} -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

#%.norm.$(BIGWIG_BIN_SIZE).bw: %.bam %.bam.bai
#    bamCoverage --binSize=$(BIGWIG_BIN_SIZE) --normalizeUsing=CPM -b $< -o $@ --numberOfProcessors=$(CORES)
rule get_norm_bw_binsize:
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.{bin_size}.norm.bw"
    params:
        threads = config["CORES"],
        bin_size = "BIGWIG_BIN_SIZE"
    shell: 
        "bamCoverage --binSize={params.bin_size}  --normalizeUsing=CPM -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

# %.bed: %.bam
#     bedtools bamtobed -splitD < $< | bsort -k1,1V -k2,2n > $@
rule bam2bed:
    input: 
        "{file}.bam"
    output: 
        "{file}.bed"
    shell: 
        "bedtools bamtobed -splitD < {input} | bsort -k1,1V -k2,2n > {output}"

# %.bedGraph.ranges: %.bedGraph
#     bawk '{print $$1,$$2; print $$1,$$3}' $< | stat_base -o -g -b | bawk '{print $$1,0,$$2}' > $@
rule get_bedgraph_ranges:
    input: 
        "{file}.bedGraph"
    output: 
        "{file}.bedGraph.ranges"
    shell: 
        "bawk '{{print $1,$2; print $1,$3}}' {input} | stat_base -o -g -b | bawk {{print $1,0,$2}}' > {output}"


# %.bedGraph.gff3: %.bedGraph %.bedGraph.ranges
#     bawk '$$4>$(GFF3_MIN_WINDOW_READS)' $< \
#     | repeat_group_pipe '\
#         grep -w $$1 $^2 | append_each_row -B "##sequence-region";\
#                grep -w $$1 $^2 | bawk '\''{print $$$$1,".","chromosome",$$$$2,$$$$3,$$$$4,".","+-",".","."}'\'';\
#                bawk '\''{print $$$$1,".","reads",$$$$2,$$$$3,$$$$4,".","+-",".","."}'\''\
#     ' 1 > $@
rule get_gff3:
    input: 
        bedGraph = "{file}.bedGraph",
        ranges = "{file}.bedGraph.ranges"
    output: 
        "{file}.bedGraph.gff3"
    params:
        min_window_reads = GFF3_MIN_WINDOW_READS
    shell: 
        """
        bawk '$4 > {params.min_window_reads}' {input.bedGraph} \
        | repeat_group_pipe '\
            grep -w $1 {input.ranges} | append_each_row -B "##sequence-region";\
               grep -w $1 {input.ranges} | bawk '\''{{print $$1,".","chromosome",$$2,$$3,$$4,".","+-",".","."}}'\'';\
               bawk '\''{{print $$1,".","reads",$$2,$$3,$$4,".","+-",".","."}}'\''\
        ' 1 > {output}
        """

# %.bedGraph.$(BEDGRAPH_FILTER).merged.bed: %.bedGraph
#     bawk '$$4>$(BEDGRAPH_FILTER)' $< | union -s --allow-duplicates | bawk '{print $$1";"$$2";"$$3, $$4}' | expandsets 2 | stat_base -o -g -b | tr ";" "\t" > $@
rule get_filtered_bed:
    input: 
        "{file}.bedGraph"
    output: 
        "{file}.bedGraph.{BEDGRAPH_FILTER}.merged.bed"
    params:
        bedgraph_filter = BEDGRAPH_FILTER
    shell:
        """
        "bawk '$4 > {params.bedgraph_filter} {input} | union -s --allow-duplicates | bawk '{{print $1";"$2";"$3, $4}}' \
        | expandsets 2 | stat_base -o -g -b | tr ";" "\t" > {output}
        """

# %.fa.gz: %.fastq.gz
#     zcat $< | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip >$@
rule get_fa:
    input: 
        "{file}.fastq.gz"
    output: 
        "{file}.fa.gz"
    shell: 
        "zcat {input} | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip > {output}"




#Queste non so se servono

#BW/%.bw:STAR/$(FASTQ_FILTERING)/%.STAR/Aligned.sortedByCoord.out.bw
#    mkdir -p BW
#ln $< $@
# BW/%.norm.bw: STAR/$(FASTQ_FILTERING)/%.STAR/Aligned.sortedByCoord.out.ribo.ex.norm.bw
#     mkdir -p BW
#     ln $< $@
# BW/pool: metadata.txt $(addprefix BW/, $(addsuffix .norm.bw,$(SAMPLES)))
#     unhead $< | cut -f -2 | collapsesets 1 | perl -pe 'print "cd BW; bigWigMerge "; s/;/.norm.bw /; s/\t/.norm.bw\t/; s/\n/.pool.norm.bw\n/' | parallel -j $(CORES)
# ALL_BW: $(addprefix BW/, $(addsuffix .bw,$(SAMPLES))) $(addprefix BW/, $(addsuffix .norm.bw,$(SAMPLES))) BW/pool
#     @echo done