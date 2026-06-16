
# =============================================================================
# 1. BAM/CRAM/SAM MANIPULATION
# =============================================================================

rule bam2cram:
    """Convert BAM to CRAM using a reference genome for better compression."""
    input: 
        "{file}.bam"
    output: 
        "{file}.cram"
    params:
        threads = config["CORES"],
        genome = GENCODE_GENOME_FASTA
    shell: 
        "samtools view -@ {params.threads} -T {params.genome} -C -o {output} {input}"

rule all_bai:
    input: 
        expand("Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.H.unique.bam.bai", sample=SAMPLES)

rule get_bai:
    """Index a BAM file to create a .bai file."""
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.bai"
    wildcard_constraints:
        file=".*Results/pass2.*|.*star_2pass.*"
    shell: 
        "samtools index {input}"

rule get_bam_id:
    """Extract unique read IDs from a BAM file."""
    input: 
        "{file}.bam"
    output: 
        "{file}.bam.id"
    shell: 
        "samtools view {input} | cut -f 1 | bsort -S8% | uniq > {output}"

rule bam2bed:
    """Convert BAM alignments to BED format using bedtools."""
    input: 
        "{file}.bam"
    output: 
        "{file}.bed"
    shell: 
        "bedtools bamtobed -splitD < {input} | bsort -k1,1V -k2,2n > {output}"

rule bam_unique:
    input:
        bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.bam",
    output:
        unique_bam="Results/pass2/{sample}/Aligned.sortedByCoord.out.ribo.ex.unique.bam",
    shell:
        """
        samtools view -h {input.bam} \
        | awk '$1 ~ /^@/ || $0 ~ /NH:i:1/' \
        | samtools view -b -o {output.unique_bam}
        samtools index {output.unique_bam}
        """

# =============================================================================
# 2. BIGWIG AND BEDGRAPH (VISUALIZATION)
# =============================================================================

if not 'BIGWIG_BIN_SIZE' in globals():
    BIGWIG_BIN_SIZE = 5

rule get_bw:
    """Generate a BigWig track from a BAM file."""
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

rule get_norm_bw:
    """Generate a CPM-normalized BigWig track."""
    input: 
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output: 
        "{file}.norm.bw"
    params:
        threads = config["CORES"],
        bin_size = BIGWIG_BIN_SIZE
    shell: 
        "bamCoverage --binSize={params.bin_size}  --normalizeUsing=CPM -b {input.bam} -o {output} --numberOfProcessors={params.threads}"

rule get_bw_binsize:
    """Generate a BigWig track with variable bin size."""
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

rule get_norm_bw_binsize:
    """Generate a CPM-normalized BigWig track with variable bin size."""
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

rule get_bedgraph_ranges:
    """Extract start-end ranges from a bedGraph file."""
    input: 
        "{file}.bedGraph"
    output: 
        "{file}.bedGraph.ranges"
    shell: 
        "bawk '{{print $1,$2; print $1,$3}}' {input} | stat_base -o -g -b | bawk {{print $1,0,$2}}' > {output}"

if not "GFF3_MIN_WINDOW_READS" in globals():
    GFF3_MIN_WINDOW_READS = 0

rule get_gff3:
    """Convert bedGraph to GFF3 format, filtering by min window reads."""
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
               grep -w $1 {input.ranges} | bawk '\''{{print $1,".","chromosome",$2,$3,$4,".","+-",".","."}}'\'';\
               bawk '\''{{print $1,".","reads",$2,$3,$4,".","+-",".","."}}'\''\
        ' 1 > {output}
        """

if not "BEDGRAPH_FILTER" in globals():
    BEDGRAPH_FILTER = 100

rule get_filtered_bed:
    """Filter bedGraph entries and merge adjacent regions."""
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

# =============================================================================
# 3. NORMALIZATION AND R-BASED SCRIPTS (EDGER)
# =============================================================================

rule tmm:
    """Calculate TMM normalization factors and CPM counts using edgeR."""
    input:
        "{file}.gz"
    output: 
        "{file}.tmm.gz"
    shell: """
        r -e 'library(edgeR);\
            x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\
            y <- DGEList(counts=x);\
            y <- calcNormFactors(y,method="TMM");\
            write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\
            write.table(cpm(y, normalized.lib.sizes=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)';
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

rule ltmm:
    """Calculate TMM normalization factors and Log2-CPM counts using edgeR."""
    input: 
        "{file}.gz"
    output: 
        "{file}.ltmm.gz"
    shell: """
        r -e 'library(edgeR);\
            x <- read.table("{input}", header=T,check.names=FALSE,row.names=1);\
            y <- DGEList(counts=x);\
            y <- calcNormFactors(y,method="TMM");\
            write.table(y$samples,"{output}.factors", sep="\t", quote=F, col.names=NA, row.names=T);\
            write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "{output}.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
        (echo -n "Geneid"; cat {output}.tmp) | gzip > {output}
        rm {output}.tmp
    """

# =============================================================================
# 4. DOWNSTREAM CLEANUP
# =============================================================================



rule clean_bam_downstream:
    shell:
        """
        if [[ $(ls star/*.srt.mapq.*bam) ]]; then
            echo ""
            echo "clean downstream processed bam"
            echo ""
            for i in $(ls bam/*.srt.mapq.*bam); do
                echo rm $i
                rm $i
            done
        else
            echo "no downstream processed bam files found"
        fi

        if [[ $(ls bam/*.srt.mapq.*bam.bai) ]]; then
            echo ""
            echo "clean downstream processed bai"
            echo ""
            for i in $(ls bam/*.srt.mapq.*bam.bai); do
                echo rm $i
                rm $i
            done
        else
            echo "no downstream processed bai files found"
        fi
    """
rule symlink_star_bai_short:
    input:
        "Results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        "star/{sample}.bai"
    shell:
        "mkdir -p star; ln -sfr {input} {output}"
