#!/usr/bin/env Rscript
suppressMessages(library("optparse"))
suppressMessages(library("limma"))
suppressMessages(library("GEOquery"))

option_list=list(
  make_option(c("-m", "--medip"), action="store_true", default=FALSE, help="If count_table comes from our Medip pipeline or the counts_table contains the first three columns as chr, start, end [default \"%default\"]"),
   make_option(c("-l", "--labels"), action="store", default="", dest="labels", metavar= "LABELS", help="Names of the samples used to create the peak table. On otpion [default \"%default\"]")
  )


parser<-OptionParser(usage = "%prog counts_table input_info > eset.rda
                     'count_table' is the path of a delimited file (.txt) containing numeric data to be adjusted. Each column must represent a different sample, each row must correspond to a gene.
                     'metadata.txt' is the path of a delimited file (.txt) containing informations about samples classification. The first column must contain sample names and the following ones must represent batch, covariates and other metadata
The colnames of the counts should be in the sample column of the metadata.txt and should be unique.

Obtains an eset.rda file with an eSet object containing counts from count_table and the phenoData of samples coming from metadata.txt file. 
", option_list=option_list)

# I run with also the *.merge.bam.bai hence I hard remove them here. To be fixed
arguments <- parse_args(parser, positional_arguments = 2)
opt<-arguments$options

sink(stderr())#the stdout ahs to be kept clean form messages

args <- arguments$args

stdin = file(args[2])
pData=read.table(stdin,h=T,sep="\t",quote="",stringsAsFactors = F)


if (!"sample"%in%colnames(pData)){
	stop("There is no column named 'sample' in your 'metadata.txt' file")
}


pData$sample = make.names(pData$sample)



counts = read.delim(file=file(args[1]),header=T,row.names=1);

if(!all(pData$sample%in%colnames(counts))){
  stop("Not all the samples in the metadata are in the counts, check the sample names")
}

if(any(duplicated(pData$sample))){
  stop("There are duplicated samples in the metadata file")
}

counts = counts[,pData$sample]

metadata <- data.frame(labelDescription= colnames(pData),row.names=colnames(pData))

phenoData <- new("AnnotatedDataFrame", data=pData[which(pData$sample%in%colnames(counts)),], varMetadata=metadata) # In this way the pData is reduced to what you have in counts.
# The metadata.txt should have 'samples' unique and hopefully as the first column. Only if 'samples' is unique I can apply the row.names. 
# I have to apply the row.names otherwise the cols of counts and the rows of pData(phenoData) are not equivalent
row.names(pData(phenoData))=pData(phenoData)$sample



experimentData <- new("MIAME",
name="Jose M. Garcia Manteiga",
lab="CTGB Lab",
contact="ctgb_lab@hsr.it",
title="Title of the experiment",
abstract="A description of the experiment",
url="www.lab.not.exist.com",
other=list(
notes="Created from text files"
))


fData=as.data.frame(row.names=row.names(counts),stringsAsFactors=F,counts[,1:3])

fData.metadata <- data.frame(labelDescription= c(
"Chromosome",
"Start position of Peak",
"End position of Peak"))

annotation = new("AnnotatedDataFrame",data=fData,varMetadata=fData.metadata)
## The assayData expression forces counts to have the columns in the same order given by the metadata

# Eventually we could add a table to annotate genes to include additonal annotation in fData/annotation

if (opt$medip){
eSet <- ExpressionSet(assayData=as.matrix(counts[,4:ncol(counts)]),
	phenoData=phenoData,
	experimentData=experimentData,
	featureData=annotation)#,
	#annotation="Ensembl_GTF_v2")
} else {
eSet <- ExpressionSet(assayData=as.matrix(counts),
	phenoData=phenoData,
	experimentData=experimentData,
	featureData=annotation)#,
	}

save(eSet,file="/dev/stdout")

#OUTPUT Warnings in STDERR
w=warnings()
sink(stderr())
if(!is.null(w)){
print(w)
        }
sink()

