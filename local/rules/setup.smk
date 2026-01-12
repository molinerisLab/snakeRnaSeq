#############
# Genecode  #
#############

GENCODE_DIR=REFERENCE_ROOT+"/bioinfotree/task/gencode/dataset/"+config["GENCODE"]["SPECIES"]+"/"+config["GENCODE"]["VERSION"]
GENCODE_GENOME_FASTA=GENCODE_DIR+"/"+config["GENCODE"]["GENOME"]+".genome.fa"
GENCODE_ANNOTATION_GTF=GENCODE_DIR+"/"+config["GENCODE"]["ANNOTATION"]+".gtf"
GENCODE_ANNOTATION_BED=GENCODE_DIR+"/"+config["GENCODE"]["ANNOTATION"]+".bed"
RSEQC_REF_BED=GENCODE_DIR+"/basic.annotation.bed"
GENOME_KEY = "GRCh"

##########
# Star  #
##########

STAR_VERSION = shell("STAR --version", read=True).strip()
STAR_GENOME_DIR = GENCODE_DIR+"/"+config["GENCODE"]["GENOME"]+".star_index/"+STAR_VERSION
FEATURECOUNTS_PARAM = config["FEATURECOUNTS_PARAM"] + " -s" + config["STRANDED"]
LCPM_PRIOR_COUNT = 0.25 
LTPM_PRIOR_COUNT = 0.0001


##############
reference_fasta_path     = f"{config['REFERENCE_DIR']}/{config['GENOME_ASSEMBLY']}.genome.fa"
annotation_gtf_path      = f"{config['REFERENCE_DIR']}/46/primary_assembly.annotation.gtf" #TODO: make dynamic 46
transcriptome_fasta_path = f"{config['REFERENCE_DIR']}/gencode.v{config['GENCODE_RELEASE']}.transcripts.fa"