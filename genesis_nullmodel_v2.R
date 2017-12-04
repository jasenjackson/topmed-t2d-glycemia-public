args<-commandArgs(trailingOnly=T)
#===mandatory parameters
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
# genotype.files <- args[5]
sample.file = args[5]
label <- args[6]

#==optional parameters
kinship.matrix <- args[7]
pheno.id <- args[8]


library(plyr)
# install.packages("gap",repos='http://cran.us.r-project.org')
# library(gap)



# This function provides a basic check that the phenotype file is a suitable
# for seqMeta.  It checks that the required column names are in the phenotype
# data frame. 
#
checkPhenotype <- function(p, outcome, covariates, id.col=NULL, gender.col=NULL) {
  if (!is.null(id.col)) {
    if (anyDuplicated(p[ , id.col])) {
      stop("Duplicated phenotype ids.")
    }
  }
  
  missing.covariates <- !(covariates %in% colnames(p))
  if (any(missing.covariates)) {
    msg <- paste("Covariates:", covariates[missing.covariates], "not found in phenotype file.\n", sep=" ")
    print(colnames(p))
    print(covariates %in% colnames(p))
    print(covariates[covariates %in% colnames(p)])
    stop(msg)
  } 
  return(invisible(NULL)) 
}



# This function reduces a data set to only the variables used in a model
# removing subjects with missing data.  Also, it makes the row names of
# the resulting data fram the subject identifier
#
# JAB addition: subsets to complete cases (i.e. no NAs in outcome or covariates)
#
# p: a data frame containing the variables in the model
#
# formula: a character vector which can be coered to an object of class 
#          "formula" (with as.formula): a symbolic description of the model to
#          be fitted. The details of model specification are given under 
#          'Details' in the "lm" help file.
#
# id: (optional) colunm name identifier of the subjects
#
# gender: (optional) colunm name identifier for the gender classification of
#         the subjects.
#
# returns: data frame with only the columns specified in the formula and with
#          the (optional) row names as the subject identifier.
#

reducePheno <- function(pheno.data, 
                        outcome, 
                        covariates = NULL, 
                        hetvars = NULL, 
                        id=NULL, 
                        gender=NULL) {
  checkPhenotype(pheno.data, outcome, covariates, id.col=id, gender.col=gender)   
  if (!is.null(id)) {
    rownames(pheno.data) <- pheno.data[ ,id]
  }
  
  all.terms <- unique(c(outcome, covariates, hetvars, gender))
  cat('all terms',print(all.terms),'\n')
  pheno.data <- as.data.frame(pheno.data) 
  pheno <- na.omit(pheno.data[, all.terms, drop=F])
  return(pheno)
}

# Calculate MAF
#
# dose: matrix with dosages rows individuals, columns are variants


Maf <- function(dose){
  aaf <- colMeans(dose,na.rm=T)/2
  return(min(1-aaf, aaf))
}


split.by.comma <- function(cur.string){
  cur.string <- gsub('"', '', cur.string)
  out <- unlist(strsplit(cur.string, ","))
  if (length(out) == 0){
    out = NULL
  }
  return(out)
}

GetFamilyDistribution <- function(response.type) {
  if (response.type == "continuous"){
    family = "gaussian"
  } else if (response.type == "dichotomous"){
    family = "binomial"
  } else {
    msg = paste("Don't know how to deal with response type", response.type)
    stop(msg)
  }
  return(family)
}

GetKinshipMatrix <- function(kinship.matrix){
  cat('Loading Kinship Matrix:',kinship.matrix,'\n')
  if(grepl('Rda',kinship.matrix,ignore.case=TRUE)){
    kmatr = get(load(kinship.matrix))
  }
  else{
    kmatr = as.matrix(read.csv(kinship.matrix,as.is=T,check.names=F,row.names=1))
  }
  
  cat('Loaded Kinship NROW:',NROW(kmatr),' NCOL:',NCOL(kmatr),'\n')
  kmatr
}

cat('label',label,'\n')
cat('kinship.matrix',kinship.matrix,'\n')

cat('outcome.type',outcome.type,'\n')



suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
# suppressMessages(library(gap))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
#suppressMessages(library(parallel))
## Setup

covariates <- split.by.comma(covariate.string)  


## phenotype 
# phenotype.data <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)
phenotype.data <- fread(phenotype.file,sep="\t",header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
phenotype.data <- phenotype.data[!duplicated(phenotype.data[,1]),]
if(NCOL(phenotype.data) < 2){
  
  msg = paste("Is the phenotype file a CSV?  Too few columns from read.csv()")
  warning(msg)
}


cat('Input pheno N=',nrow(phenotype.data),'\n')
pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates, id=pheno.id)

cat('Output pheno N=',nrow(pheno),'\n')

## Report dropped individuals
dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
  cat("Dropped because of incomplete cases:", length(dropped.ids) )
}


# load sample ids
sample.ids <- unique(readLines(sample.file))
#For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
all.terms <- unique(c(outcome.name, covariates))
print(covariates)
print(row.names(pheno))
pheno <- pheno[row.names(pheno) %in% sample.ids,na.omit(all.terms,drop=F)]
cat('Output pheno after mergeing with Genos N=',nrow(pheno),'\n')
if(nrow(pheno) == 0){
  msg = paste("Phenotype ID column doesn't match IDs in GDS")
  stop(msg)
}

full.sample.ids <- sample.ids 

#subset to phenotyped samples
# seqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
# sample.ids <- seqGetData(f, "sample.id")
# pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]

## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = GetKinshipMatrix(kinship.matrix)
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
if(nrow(pheno) == 0){
  msg = paste("Phenotype ID column doesn't match IDs in Kinship Matrix")
  stop(msg)
}


# Get sample ids to check order 
# seqSetFilter(f,sample.id = row.names(pheno))
# sample.ids <- seqGetData(f, "sample.id")

# if (!(identical(sample.ids,row.names(pheno)) && identical(row.names(kmatr),row.names(pheno)))){
# stop("Something is off problem with re-ordering")
# }

# seqClose(f)

sample.data <- data.frame(scanID = row.names(pheno),  
                          pheno, 
                          stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
modified.pheno = pheno[full.sample.ids,,drop=FALSE]
row.names(modified.pheno) <- full.sample.ids

sample.data.for.annotated <- data.frame(sample.id = full.sample.ids,
                                        modified.pheno,
                                        stringsAsFactors=F)
rm(modified.pheno)

annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)

###################
## NULL MODEL
##################
# Should depend on response type

cat('start fit....\n')
kmatr = as.matrix(kmatr)
cat('Fitting model ')
nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr)

write(full.sample.ids, file =paste(label,"_sample_ids.txt",sep=""),
      ncolumns = 1,
      append = FALSE, sep = " ")

save(nullmod,annotated.frame,file=paste(label,".RDa",sep=""))






####################
# Summary Stats
####################

## Generate summary statistics & plots (if outcome.type == dichotomous)
## Future versions of this script will generate summary statsitics & plots regardless of outcome type
library("ggplot2")
if (outcome.type == "dichotomous"){ 
  
  ## identify quantitative covariates (Future versions of this script may have more robust selection method)
  cov_to_remove <- c("study_ancestry","STUDY_ANCESTRY","sex")
  quant_covars <- covariates[! covariates %in% cov_to_remove]
  
  ## Create dataframes including: all subjects (pheno),
  #  males (m), females (f), cases (case), controls (control), 
  #  male cases (m_case), male controls (m_control), 
  #  female cases(f_case) and female controls (f_control)
  
  pheno <- pheno[!is.na(pheno[,outcome.name]),]
  case = pheno[pheno[,outcome.name] == 1,]
  control = pheno[pheno[,outcome.name] == 0,]
  male.val <- "M" ## look for initial M/F coding
  female.val <- "F"
  m = pheno[pheno$sex==male.val,]
  f = pheno[pheno$sex==female.val,]
  if (NROW(m)==0 && NROW(f)==0){ ## otherwise, look for 1/2 coding
    male.val <- 1
    m = pheno[pheno$sex==male.val,]
    female.val <- 2
    f = pheno[pheno$sex==female.val,]
  }
  m_case = case[case$sex==male.val,]
  m_control = control[control$sex==male.val,]
  f_case = case[case$sex==female.val,]
  f_control = control[control$sex==female.val,]
  
  ## initiate data_frame for final stats table (include first row with counts, n, for each group)
  stats_df <- data.frame(matrix(NA, nrow= 1, ncol = 9))
  colnames(stats_df) <- c("all", "m", "f", "case", "control", "m_case", "m_control", "f_case", "f_control")
  rownames(stats_df) <- "n"
  stats_df[1,] <- c(nrow(pheno),nrow(m),nrow(f),nrow(case), nrow(control), nrow(m_case),nrow(m_control),nrow(f_case),nrow(f_control))
 
  ## look for "STUDY_ANCESTRY" coding..  
  if ("STUDY_ANCESTRY" %in% covariates){
    study_ancestry_val <- "STUDY_ANCESTRY"
  } else {
    study_ancestry_val <- "study_ancestry"
  }
  
  ## prepare "plot_pheno" dataframe for facetting 
  plot_pheno <- pheno
  colnames(plot_pheno)[which(names(plot_pheno) == outcome.name)]  <- "outcome"## rename outcome column to generic "outcome" 
  plot_pheno[,"outcome"] <- as.factor(plot_pheno[,"outcome"]) 
  levels(plot_pheno[,"outcome"]) <- c("Controls", "Cases")
  
  ##Iterate through each quantitative covariate, 
  # generate boxplot+violinplot for both and add
  # values to statistics dataframe

  pdf(paste(label,"_plots.pdf",sep=""), width=11)
  for (i in quant_covars){
    print(i)
    
    ## generate table for stats of i trait
    i_stats_df <- data.frame(matrix(NA, nrow = 3, ncol = ncol(stats_df)))
    i_stats_df[,1] <- c(mean(pheno[,i]), median(pheno[,i]), sd(pheno[,i])) ##fill out i stats for ALL
    i_stats_df[,2] <- c(mean(m[,i]), median(m[,i]), sd(m[,i])) ##fill out i stats for m
    i_stats_df[,3] <- c(mean(f[,i]), median(f[,i]), sd(f[,i])) ##fill out i stats for f
    i_stats_df[,4] <- c(mean(case[,i]), median(case[,i]), sd(case[,i])) ##fill out i stats for f
    i_stats_df[,5] <- c(mean(control[,i]), median(control[,i]), sd(control[,i])) ##fill out i stats for m_case
    i_stats_df[,6] <- c(mean(m_case[,i]), median(m_case[,i]), sd(m_case[,i])) ##fill out i stats for f
    i_stats_df[,7] <- c(mean(m_control[,i]), median(m_control[,i]), sd(m_control[,i])) ##fill out i stats for m_control
    i_stats_df[,8] <- c(mean(f_case[,i]), median(f_case[,i]), sd(f_case[,i])) ##fill out i stats for f_case
    i_stats_df[,9] <- c(mean(f_control[,i]), median(f_control[,i]), sd(f_control[,i])) ##fill out i stats for m_case
    
    ## merge table with stats_df
    names(i_stats_df) <- names(stats_df)
    row.names(i_stats_df) <- c(paste("mean_",i,sep=""),paste("median_",i,sep=""), paste("sd_of_",i,sep=""))
    stats_df <- rbind(stats_df, i_stats_df)
    
    ## plot boxplots of  i distribution in cases/controls split by ancestry & sex 
    plot_title <- paste(i," By Ancestry (Boxplot)", sep="")
    plot <- ggplot(plot_pheno, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1))
    print(plot + geom_boxplot(aes(fill=factor(sex))) + scale_fill_discrete(name = "", labels=c("Female", "Male")) + labs(title = plot_title, x="") + facet_grid(outcome ~ .))
    
    ## plot violin plots of same data 
    plot_title <- paste(i," By Ancestry (violin plot)", sep="")
    plot <- ggplot(plot_pheno, aes_string(study_ancestry_val, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1))
    print(plot + geom_violin(aes(fill=factor(sex))) + scale_fill_discrete(name = "", labels=c("Female", "Male")) + labs(title = plot_title, x="") + facet_grid(outcome ~ .))
    
  }
  dev.off()
  stats_df <- data.matrix(stats_df)
  write.table(x=stats_df, file=paste(label,"_stats.tsv",sep=""), col.names = TRUE, row.names = TRUE, sep="\t")
}
