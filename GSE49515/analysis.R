#GSE49515 analysis (only HCC vs Healthy)

library(BiocManager)
library(GEOquery)
library(affy)
library(oligo)
library(hgu133plus2.db)
library(org.Hs.eg.db)
#==============================================================================
#Retrieving the whole GSE as a list of ExpressionSets
gse33146 <- getGEO('GSE49515')
gse66417 <- gse66417[[1]]
#Reading CEL data using the oligo package (Raw data)
library(oligoClasses)
gse33146_celdata <- read.celfiles(list.celfiles('GSE49515',full.names=TRUE,listGzipped=TRUE))
#===============================================================================
library(tidyverse)
pd <- pData(gse33146)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse33146_celdata <- read.celfiles(paste0('GSE49515/',pd$cel_file),phenoData=phenoData(gse49515))
pData(gse33146_celdata)[,c("geo_accession","cell line:ch1","culture medium:ch1")]
varLabels(gse66417)
##upstream analysis of microarray data
#accessing the expression data
#the rma function returns the function basicRMA() from the oligo package takes in an ExpressionFeatureSet object. 
##Hence, the following one line of code performs RMA normalisation, 
#and returns an ExpressionSet object containing the background corrected, normalised, and summarised expression data.
gse33146_eset <- rma(gse33146_celdata)
expression_data <- exprs(gse33146_eset)
class(expression_data)
expression_data <- as.data.frame(expression_data)
write.csv(expression_data, "expression.csv")
#===============================================================================
expression_data <- read.csv("expression.csv")
pheno_data <- read.csv("pheno_data.csv")
nm1 <- c("GSM1200294", "GSM1200295", "GSM1200296", "GSM1200297", "GSM1200298", "GSM1200299", "GSM1200300", "GSM1200301", "GSM1200302", "GSM1200303", "GSM1200304", "GSM1200305", "GSM1200306", "GSM1200307", "GSM1200308", "GSM1200309", 
         "GSM1200310", "GSM1200311", "GSM1200312", "GSM1200313")

design <- model.matrix(~0 + status, data = pheno_data)
colnames(design) <- c("healthy", "HCC")
design

contrast_matrix <- makeContrasts(HCC - healthy, levels=design)
contrast_matrix

library(limma)
eset <- read.csv("agg.csv")
fit <- lmFit(eset,design)
fit2 <- contrasts.fit(fit,contrasts=contrast_matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2,lfc=1))
#===============================================================================
RES <- limma::topTable(fit2, number = Inf) 
DEGS <- limma::topTable(fit2, number = Inf , p.value = 0.05, lfc = 1.5) 
write.csv(RES, "RES.csv")
write.csv(DEGS, "DEGS_lfc1.5.csv")
#===============================================================================




