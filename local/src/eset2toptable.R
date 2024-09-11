#!/usr/bin/env Rscript
#options(error = function() traceback(2))
suppressMessages(library("optparse"))
suppressMessages(library("GenomeInfoDb"))
#suppressMessages(library("GEOquery"))#required by pData
suppressMessages(library("limma"))
suppressMessages(library("edgeR"))
suppressMessages(library("statmod"))
#suppressMessages(library("WriteXLS"))
suppressMessages(library("DESeq2"))
#suppressMessages(library("ChIPpeakAnno")) # To uncomment whenever the library will be installed in the cluster

option_list <- list(
        make_option(c("-c", "--conditions_colname"), action="store", default="", dest="conditions_colname", metavar = "CONDITIONS_COLNAME", help="The name of the column containing the condition factor [default \"%default\"]"),
	make_option(c("-d", "--check_reference_condition"), action="store_true", default=FALSE, help="Print the names of the factors in condition column to choose reference condition for contrasts. A conditions column must be also specified with -c option [default \"%default\"]"),
	make_option(c("-m", "--microarray"), action="store_true", default=FALSE, help="Analyze microarray data rather than RNA-Seq counts. For RNA-Seq counts, the voom function transform reads in the eset dataset obtained with featureCounts and normalize with calcNormFactors prior to limma linear fit analysis using edgeR and limma functions [default \"%default\"]"),
	make_option(c("-l","--min_expr"),action="store",default=1, dest="min_expr", metavar = "MIN_EXPR", help="The minimun cpm value to consider a gene present in a single sample [default \"%default\"]"),
	make_option(c("-n","--number_replicates"),action="store",default="", dest="number_replicates", metavar = "NUMBER_REPLICATES", help="The number of replicates per group of samples to be used in the design. This number will be used to set the minimum amount of samples in which a gene have to pass MIN_EXPR cpm to be considered as expressed [default \"%default\"]"),
	make_option(c("-f","--factorial"),action="store", default="", dest="factorial", metavar= "FACTORIAL", help="Factorial 2x2 design. FACTORIAL contains the name of the column containing the factors to be combined in the design separated by a space [default \"%default\"]"),
	make_option(c("-t","--tool"), action="store", default="limma", dest="tool",metavar="TOOL", help="Statistical tool to be used for differential expression analysis amongst 'limma','deseq2','edger' [default \"%default\"]"),
	make_option(c("-e","--medip"), action="store_true", default=FALSE, help="If eset.rda comes from our counts_table2eset with -m (medip), eset2toptable will work with edgeR and ajust output to annotate peaks with 'ChIPpeakAnno' [default \"%default\"]"),
	make_option(c("-r","--trend"), action="store_true", default=FALSE, help="trend option of limma eBayes [default \"%default\"]"),
	make_option(c("-z","--normalization"), action="store", default="TMM", help="edger normalization tyope, available values: TMM, upperquartile, RLE, none [default \"%default\"]"),
	make_option(c("-a","--anova_like"),action="store_true",default=FALSE, dest="anova_like",  help="Save the result in the object 'anova_like' in output .RData file. Available only in edgeR mode. [default \"%default\"]"),
	make_option(c("-s","--robust"),action="store_true",default=FALSE,dest="robust",help="robust TRUE in QLF and estimate dispersion, suitable to mitigate outlier [default \"%default\"]")
)

usage = "%prog [options] eset.rda FORMULA CONTRASTS
Transform an eSet object (Expression Dataset) to produce an RData object containing a top table data.frame for each contrast in a list (with toptable functions or similar) with either of the three tools, limma, deseq2 or edger. The results are chosen from a single condition. The -m microarray option is classically applying a reference level(e.g. normal, untreated, healthy...) defined by the user. [DEPRECATED]
FORMULA is a classical R formula that limma, edger or deseq will use to produce a design. Typical formulas include:	
~ condition ( The reference level of the factor 'condition' is taken as intercept and coeficients of this fitting are contrasts of the other levels against it. Used with microarray and refernce condition. [DEPRECATED]
~ 0+condition ( There is no reference level and specific contrasts are to be defined ) 
~ 0+condition+treatment ( Two factors design)
~ batch+condition (A condition corrected by batch effect, or a simple paired experiment)
~ 0+condition + covariate1 + covariate2 ... ( Covariates are introduced to correct for them in the design)
CONTRASTS are a list of expressions of the kind 'A-B', '(A+B)-(C+D)' defining the contrasts of the linear model that limma, edger or deseq2 will include in the function makeContrasts(). The expression should contain valid levels of the condition column.
For factorial designs(-f,--factorial) limma is to be chosen as tool (edger and deseq2 compatibility under development), whereby a combination of both factors (TS=paste(condition1,condition2)) and then all the desired contrasts are to be specified.
...

.DOC: stdout
	RData
"

parser<-OptionParser(usage = usage, option_list=option_list)


sink(stderr())
#### Assign arguments to objects
#### Define a function to take the first element of a list and remove it from it
shift_fn <- function(x) {
  if(length(x) == 0) {return(NA)}
  shiftret <- x[1]
  assign(as.character(substitute(x)), x[2:(length(x))], parent.frame())
  return(shiftret)
}

arguments <- parse_args(parser, positional_arguments = c(1,Inf))
args <- arguments$args
opt<-arguments$options

if(opt$microarray && opt$tool!="limma"){
  stop("ERROR: -m option requires -t limma")
}

if(opt$medip){opt$tool="edger"}

if(nchar(opt$factorial)>0 && opt$tool!="limma"){
  stop("ERROR: -f option requires -t limma")
}



### I take the first argument, the eset file and assign to eset object
stdin <- file(shift_fn(args))
open(stdin, blocking=TRUE)
eset = get(load(stdin))

if (nchar(opt$number_replicates) == 0) {
	stop("Please specify the -n value")
}
opt$number_replicates=as.integer(opt$number_replicates)

if(nchar(opt$factorial)==0){
  ### I take now the design FORMULA and the remaning will be the contrasts (a list of all of them)
  formula = shift_fn(args)
  formula.t = formula # I keep the formula in text form to strip out the '0' for Deseq2 afterwards
  ## Deciding how many 0s to add to the contrast vector depending on the number of factors in design
  # First I use  stringr library
  suppressMessages(library("stringr"))
  # to calculate the number of factors
  nfactors = str_count(as.character(formula.t),fixed("+")) 
  ## I transform it in a lm formula expression
  tryCatch( {
          formula <- as.formula(formula)
  }, error = function(err) {
          stop("ERROR: You have not given a suitable formula as design!")
  })
}
	
	
	
### The remaining arguments are already the contrasts
contrasts.args<-args
####

## If run just for check the reference condition, I only open it and write out the levels of the conditions colname for the user to choose afterwards which one is the reference
if (opt$check_reference_condition){
  if (nchar(opt$conditions_colname) <= 0) {
	  stop("Please specify the --conditions_colname")
  }
  write.table(levels(pData(eset)[,opt$conditions_colname]), col.names = FALSE, row.names =FALSE, quote = FALSE, file="", sep="\t")
  quit(save="no",status=0)
}

## Checking that the conditions colname and the reference level are given


## If opt$rna_seq is not chosen it means it is an eset from an array then the voom phase is skipped

#condition = factor(pData(eset)[,opt$conditions_colname],levels=unique(as.character(pData(eset)[,opt$conditions_colname])))

if (nchar(opt$factorial) > 0){
  .f=unlist(strsplit(opt$factorial,split=" +"))
  TS = factor(c(paste(pData(eset)[,.f[1]],pData(eset)[,.f[2]],sep=".")))
  design=model.matrix(~0+TS)
  colnames(design)=levels(TS)
} else {
  TS=c()
  design = model.matrix(data=model.frame(na.action=na.fail,pData(eset)),object=formula)
  colnames(design) = gsub(paste0("^",strsplit(arguments$args[2],"\\+")[[1]][2],"\\s*"),"",colnames(design))
  #colnames(design)[1:length(levels(get_all_vars(formula,data=pData(eset))[1][,1]))]=levels(get_all_vars(formula,data=pData(eset))[1][,1]) # Puts as colnames of design the name of the levels in the first condition factor to extract the contrasts easily using the names of the levels
}
print(design)
###########
print(contrasts.args)
###########

write(contrasts.args)
contrasts = makeContrasts(contrasts=contrasts.args, levels=design)

if (exists("nfactors")){
 if (nfactors == 1) {
	add.after=c()
	add.before = 0
	} else if (nfactors == 2) {
		add.after = 0
		add.before = 0
 } else {
	if(opt$tool=="deseq2"){
		stop("eset2toptable cannot yet manage more than 2 factors or you have not put the 0 on the formula")
	}
 }
}


if (!opt$microarray && (opt$tool=="limma" || opt$tool=="edger")) {
	y.all<-DGEList(counts=exprs(eset),genes=fData(eset))
	isexpr= rowSums(cpm(y.all)>opt$min_expr) >= opt$number_replicates
	y <- y.all[isexpr,,keep.lib.sizes=FALSE]
	y <- calcNormFactors(y,method=opt$normalization)
}else{
	e <- exprs(eset)
	isexpr <- rowSums(e>opt$min_expr) >= opt$number_replicates
	y <- e[isexpr,]
}

anova_like=NA

if (opt$tool=="limma"){
	if(!opt$microarray){
		v = voomWithQualityWeights(y,design,plot=FALSE)
		voom = v$E
                my_rpkm=""
		if("Length"%in%colnames(fData(eset))){
		        my_rpkm = rpkm(y,log=T,gene.length=y$genes$Length)
		        write.table(my_rpkm,"rpkm.log.expression.values.txt",row.names=T,sep='\t')
                }
		write.table(voom,"voom.log.expression.values.txt",row.names=T,sep="\t")
		fit_init = lmFit(object=v,design=design)
	}else{
		fit_init = lmFit(y, design=design)
	}
	fit = contrasts.fit(fit_init,contrasts)
	fit = eBayes(fit,robust=T,trend=opt$trend)
	# results = decideTests(fit,method="separate",adjust.method="BH",p.value=0.05)
	top.list = lapply(colnames(contrasts),function(x){
	  	.top = topTable(fit, coef=x,number = Inf,adjust.method = "BH",sort.by = "p")
	})
	names(top.list)=colnames(contrasts)
} else if(opt$tool=="edger" || opt$medip==TRUE){
	ey = estimateDisp(y,design,robust=opt$robust) # equivalente a v = voom(y), sort of... to estimate disp.
	pdf('eset.plotBVC.pdf')
	plotBCV(ey)
	dev.off()
	efit = glmQLFit(ey,design,robust=opt$robust) # equivalente a fit = lmFit(v,design)
        if(opt$anova_like){
                anova_like = as.data.frame(topTags(glmQLFTest(efit, contrast=contrasts),adjust.method = "BH",sort.by = "p",n = Inf)$table)
        }
	top.list = apply(contrasts,2,function(x) {
	  #.elrt = glmQLFTest(efit,contrast=x)
	  .elrt = glmQLFTest(efit,contrast=x)
	  .top = as.data.frame(topTags(.elrt,adjust.method = "BH",sort.by = "p",n = Inf)$table)
	  })
	#elrt = glmLRT(efit,contrast=contrasts[,2]) # equivalente a fit2 = contrasts.fit(fit,cmat)
	#topedge = lapply(topTags(elrt) # con defaults, per essere veloce...
} else if(opt$tool=="deseq2"){
	# Build the dds object for deseq2 analysis
	# The formula to the dds object has to be changed to include the intercept, otherwise the betaPrior cannot be TRUE.
	formula.t = gsub(formula.t,pattern="~0\\+",replacement="~",fixed=F,perl=T)
	formula <- as.formula(formula.t)
	dds=DESeqDataSetFromMatrix(countData=as.matrix(exprs(eset)),colData=pData(eset),design=formula)
	deseq=DESeq(object = dds,test = "Wald", betaPrior = TRUE, minReplicatesForReplace = Inf,parallel = T)
  top.list = apply(contrasts,2,function(x) {
    .res = DESeq2::results(deseq,contrast=c(add.before,x,add.after),independentFiltering = TRUE,pAdjustMethod = "BH") # the contrasts were made with the original formula containing the 0, hence without the intercept. And so the design object made with model.matrix. We need an intercept value for the contrasts and we set it always to 0. By doing so, the contrast vector of 0,1, and -1 made withouth the intercept can be used by adding the intercept factor (0) with DESeq2 results function. I have tested that 0,0,1,-1 (intercept, A, B, C) is equivalent to B-C contrast also if you call it by contrast=c("condition","B","C") in a intercept included betaPrior=T design in Deseq2. I still need to test and implement the use of factorial designs with deseq2. For it you'd better still use edgerR and/or limma. For batch designs one should put ~0+condition+batch) you can use this formula as design and then call the contrast in the same way. However, then I have to add, not only a 0 before the vector of -1,1 contrasts, but also after, hence I test with an if the number of factors in the design, and if it is only one, I just add 1 before, if it is 2, I add also another 0 at the end, hence adding the intercept and the level of the batch condition (second factor) removed.
    .res = .res[order(.res$pvalue,decreasing=F),]
    .top = as.data.frame(.res)
  })
} else{
	stop(paste("Invalid value for option -t (",opt$tool,")"))
}
			
sink()
# This line is removed. Rather, I will save an object with all the objects needed to produce either the toptable.txt (with the write.fit line below) or the other objects needed to produce quality plots and other
#write.fit(fit, results, "", digits=5, adjust="BH", method="separate", F.adjust="BH", sep="\t")
print(opt$tool)

if(!opt$microarray){
  if(opt$tool=="limma"){
    save(list=c("TS","y","v","voom","my_rpkm","fit_init","fit","top.list"),file="/dev/stdout")
  } else if(opt$tool=="edger"){
    save(list=c("TS","y","efit","top.list","anova_like"),file="/dev/stdout")
  } else if(opt$tool=="deseq2"){
    save(list=c("TS","dds","deseq","top.list"),file="/dev/stdout")
  } else if(opt$medip==TRUE){
    save(list=c("TS","y","efit","top.list","opt$number_replicates"),file="/dev/stdout")
  }
}else{
    save(list=c("TS","y","fit_init","fit","top.list","anova_like"),file="/dev/stdout")
}


w=warnings()
sink(stderr())
if(!is.null(w)){
print(w)
}
sink()



