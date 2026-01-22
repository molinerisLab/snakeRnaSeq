import os

configfile: "config.yaml"

if os.path.exists("prj_Snakefile"):
	include: "prj_Snakefile"

FTP_PREFIX="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_" + config["SPECIES"] + "/release_"+ config["VERSION"] +"/"


#rule all:
#	input:
#		config["GENOME_ASSEMBLY"]+".primary_assembly.genome.fa",
#		"primary_assembly.annotation.gtf",
#		"basic.annotation.gtf",
#		"rseqc.HouseKeepingGenes.bed.gz",
#		config["VERSION"]+"primary_assembly.rRNA_complete.bed",
#		config["VERSION"]+"basic.rRNA_complete.bed",
#		config["VERSION"]+"/repeat_rmsk.bed.gz"

rule all_reference:
  input:
    config["GENOME_ASSEMBLY"]+".primary_assembly.genome.fa",
    "basic.annotation.gtf.gz",
    "primary_assembly.annotation.gtf.gz",
    "rseqc.HouseKeepingGenes.bed.gz",
    "primary_assembly.annotation.rRNA_complete.bed",
    config["GENOME_ASSEMBLY"]+".primary_assembly.star_index/"+config["STAR_VERSION"]+"/SA"


#rule print_config:
#	output: "print_config"
#	shell:


"""
.META: {file}.txt.gz
	1	bin
	2	swScore
	3	milliDiv
	4	milliDel
	5	milliIns
	6	genoName
	7	genoStart
	8	genoEnd
	9	genoLeft
	10	strand
	11	repName
	12	repClass
	13	repFamily
	14	repStart
	15	repEnd
	16	repLeft
	17	id

"""


# Generate this file : GRCh38.primary_assembly.genome.fa and GRCm38.primary_assembly.genome.fa
# Generate the zipped file
rule wget_fa_gz:
    output:
        config["GENOME_ASSEMBLY"]+".primary_assembly.genome.fa.gz"
    shell:
        "wget -O {output} -c "+FTP_PREFIX+config["GENOME_ASSEMBLY"]+".primary_assembly.genome.fa.gz"

# Unzip the file
rule gunzip_fa:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    shell:
        "gunzip {input} > {output}"


# Generate this file : $(HSAPIENS_VERSION)/basic.annotation.gtf.gz and $(MMUSCULUS_VERSION)/basic.annotation.gtf.gz
rule annotation:
    output:
        "{annotation}.gtf.gz"
    shell:
        "wget -O {output} -c "+FTP_PREFIX+"gencode.v"+config["VERSION"]+".{wildcards.annotation}.gtf.gz"
# Generate this file :$(HSAPIENS_VERSION)/rseqc.HouseKeepingGenes.bed.gz and $(MMUSCULUS_VERSION)/rseqc.HouseKeepingGenes.bed.gz

rule reseqc_HouseKeepingGenes:
    input:
        GENCODE_DIR+"/rseqc.HouseKeepingGenes.bed.gz"
    output:
        "rseqc.HouseKeepingGenes.bed.gz"
    shell:
        "cp {input} {output}"


# Generate this file : $(HSAPIENS_VERSION)/primary_assembly.annotation.rRNA_complete.bed and $(MMUSCULUS_VERSION)/primary_assembly.annotation.rRNA_complete.bed
# repeat_rmsk.txt.gz:
# 	wget -O $@ http://hgdownload.soe.ucsc.edu/goldenPath/$(UCSC_VERSION)/database/rmsk.txt.gz
#
# repeat_rmsk.ribosomal.bed: repeat_rmsk.txt.gz
# 	bawk '$$13=="rRNA" {print $$6,$$7,$$8,$$11,$$2,$$10}' $< \
# 	| perl -lane 'if($$F[0]=~m/_/){$$chr=shift(@F); $$chr=~m/chr[^_]+_([^_]+)_?/; $$chr=$$1; $$chr=~s/v(\\d+)$$/.\1/; $$_=$$chr."\t".join("\t",@F);} print ' \
# 	| bawk '{print $$0,$$2,$$3,0,1,$$3-$$2",",0","}' > $@                  * to bed12*
#
# %.annotation.rRNA.bed: %.annotation.gtf
# 	perl -lne 'print if m/gene_type "rRNA"/ or m/gene_type "Mt_rRNA"/ or m/gene_type "Mt_tRNA"/' $< | gtf2bed_Aronesty - > $@
#
# %.annotation.rRNA_complete.bed: %.annotation.rRNA.bed repeat_rmsk.ribosomal.bed
# 	bsort -k1,1 -k2,2n $^ >$@

rule repeat_rms_txt_gz: #This rule work
    output:
        "repeat_rmsk.txt.gz"
    shell:
        "wget -O {output} http://hgdownload.soe.ucsc.edu/goldenPath/"+config["UCSC_VERSION"]+"/database/rmsk.txt.gz"


rule repeat_rmsk_ribosomal_bed:
    input:
        "repeat_rmsk.txt.gz"
    output:
        "repeat_rmsk.ribosomal.bed"
    shell:"""
        bawk '$13=="rRNA" {{print $6,$7,$8,$11,$2,$10}}' {input}  \
	    | perl -lane 'if($F[0]=~m/_/){{$chr=shift(@F); $chr=~m/chr[^_]+_([^_]+)_?/; $chr=$$1; $chr=~s/v(\\d+)$/.\1/; $_=$chr."\t".join("\t",@F);}} print ' \
	    | bawk '{{print $0,$2,$3,0,1,$3-$2",",0","}}' > {output}
    """

rule annotation_gtf_unzip: #This rule work
    input:
        "primary_assembly.annotation.gtf.gz"
    output:
        "primary_assembly.annotation.gtf"
    shell:
        "gunzip < {input} > {output}"

rule annotation_rRNA_bed:
    input:
        "primary_assembly.annotation.gtf"
    output:
        "primary_assembly.annotation.rRNA.bed"
    shell:"""
        perl -lne 'print if m/gene_type "rRNA"/ or m/gene_type "Mt_rRNA"/ or m/gene_type "Mt_tRNA"/' {input} |
        gtf2bed_Aronesty - > {output}
        """

rule annotation_rRNA_complete_bed:
    input:
        p1="primary_assembly.annotation.rRNA.bed",
        p2="repeat_rmsk.ribosomal.bed"
    output:
        "primary_assembly.annotation.rRNA_complete.bed"
    shell:
        "bsort -k1,1 -k2,2n {input.p1} {input.p2} > {output}"



# Generate this file :$(BIOINFO_REFERENCE_ROOT)/ucsc/hsapiens/hg38_20180212/genome_chr and $(BIOINFO_REFERENCE_ROOT)/ucsc/mmusculus/mm10_20180212/genome_chr
#rule genome_chr:

# Generate this file :$(HSAPIENS_VERSION)/repeat_rmsk.bed.gz and $(MMUSCULUS_VERSION)/repeat_rmsk.bed.gz
# repeat_rmsk.bed.gz: repeat_rmsk.txt.gz
# 	bawk '{print $$genoName,$$genoStart,$$genoEnd,$$repName ";" $$repFamily ";" $$repClass,$$swScore,$$strand}' $< | gzip > $@

rule repeat_rmsk_bed_gz: #This rule wor with awk
    input:
        "repeat_rmsk.txt.gz"
    output:
        "repeat_rmsk.bed.gz"
    shell:"""
		bawk '{{print genoName,genoStart,genoEnd,repName ";" repFamily ";" repClass,swScore,strand}}' {input}  \
		| gzip > {output}
	"""

rule star_index:
	input:
		genome=config["GENOME_ASSEMBLY"]+".primary_assembly.genome.fa",
		annot="primary_assembly.annotation.gtf"
	output:
		config["GENOME_ASSEMBLY"]+".primary_assembly.star_index/{STAR_VERSION}/SA"
	threads: 12
	shell:"""
		#STAR_VERSION=$(STAR --version 2>/dev/null);\
		#REQUIRED_VERSION={wildcards.STAR_VERSION};\
		#if [ -z "$STAR_VERSION" ]; then\
		#  echo "STAR is not installed or not in PATH"\
		#  exit 1\
		#fi\
		#if [ "$STAR_VERSION" == "$REQUIRED_VERSION" ]; then\
		#  echo "STAR version is exactly $REQUIRED_VERSION"\
		#else\
		#  echo "STAR version is $STAR_VERSION, required is $REQUIRED_VERSION"\
		#  exit 1\
		#fi\
		mkdir -p $(dirname {output});\
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir $(dirname {output}) --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annot} --sjdbOverhang {config[STAR][sjdbOverhang]}
	"""

#baw/GRCh38.p14.primary_assembly.genome.bwtsw.sa
rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("bwa/{genome}.{alg}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "bwa/{genome}.{alg}.log",
    #conda:
    #    config['PRJ_ROOT']+"/local/env/bwa.yaml"
    params:
        extra=lambda w: f"-a {w.alg}",
    wrapper:
        "v4.3.0/bio/bwa/index"


### Rule kallisto index ### TODO: modify the rule to use the same logic as above

rule build_kallisto_index:
    input:
        transcripts=transcriptome_fasta_path,
        downloaded=f"{REFERENCE_DIR}/.transcriptome_fasta_downloaded"
    output:
        f"{REFERENCE_DIR}/kallisto_index/transcripts.idx"
    conda: "transcript_env.yaml"
    shell:
        """
        mkdir -p {REFERENCE_DIR}/kallisto_index
        kallisto index \
            -i {output} \
            {input.transcripts}
        """

VERSION = "GRCh38"
def get_rseqc_url(version):
    base = "https://sourceforge.net/projects/rseqc/files/BED"
    if version in ["GRCh38", "hg38"]:
        return f"{base}/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download"
    elif version in ["GRCm38", "mm10"]:
        return f"{base}/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed.gz/download"
    else:
        raise ValueError(f"Unknown version: {version}")

rule download_rseqc_housekeeping:
    output:
        f"{GENCODE_DIR}/rseqc.HouseKeepingGenes.bed.gz"
    params:
        url = get_rseqc_url(VERSION)
    shell:
        """
        # 3. Use curl with -L (follow redirects) and -k (insecure/skip SSL check)
        curl -L -k -o {output} "{params.url}"
        """