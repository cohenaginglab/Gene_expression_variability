# Author: Frederik Dufour
# Date: Mon Feb 17 10:23:51 2020
# Title: Process ComBat 
# Description: Process data with the ComBat algorithm
#------------------------------------------------------
# init
#------------------------------------------------------
## libraries
libs = list("data.table","ggpubr","tidyverse","beadarray","broom","sva")
lapply(libs,require,character.only = T)

## Set workdir 
## Change to reuse the script
WORKDIR = setwd("~/Documents/data/")

## Files and directories    
dat_Filename = "AdNeuro/data/GSE63063_rawData.txt"
meta_Filename = "AdNeuro/metadata/GSE63063_metadata.RData"
DATASETNAME = "AdNeuro_Combined"

## Load data
dat = data.frame(fread(dat_Filename,data.table = F),row.names = T)
N = ncol(dat)
gc()

## Load metadata
obj_name = load(meta_Filename)
assign("metadat",get(obj_name))

## Filter missing age or sex
metadat = filter(metadat, !is.na(age) & !is.na(gender)) 

## Check sample are in the same order in metadata and data
idx = na.omit(match(metadat$fileID,colnames(dat)))
dat = dat[,idx]
idx = na.omit(match(colnames(dat),metadat$fileID))
metadat = metadat[idx,]

if(!all(colnames(dat) == metadat$fileID)) stop("Metadat and dat ids don't match!!!")

## Remove probes with NAs
dat = na.omit(dat)
print(paste0("DIM: RAW: ",paste0(dim(dat),collapse = ", ")))

## Background subtract:  RMA 
dat = limma::backgroundCorrect.matrix(dat,method = "normexp",normexp.method = "mle")

## Transform data
dat = apply(dat,2,log2)

## Adjust data for platform
dat = sva::ComBat(dat = dat,batch = metadat$platform,mod = NULL)
gc()

## Adjust data for age and sex
mod0 = model.matrix(~splines::bs(metadat$age,df = 3) + as.factor(gender),data = metadat)
fit = limma::lmFit(dat,design = mod0)
dat = residuals(fit,dat)
gc()

## Adjust data for batch effects 
dat = sva::ComBat(dat = dat,batch = metadat$batchid,mod = NULL)
gc()

## Calculate detection p-value using unmatch probes
probes = row.names(dat)
require(illuminaHumanv4.db)
qual=unlist(mget(probes, illuminaHumanv4.db::illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
## Create annotatedDataFrame for pheno data
row.names(metadat) = metadat$fileID
phenoData <- new("AnnotatedDataFrame",data=metadat)
exprSet = ExpressionSet(assayData =as.matrix(dat),phenoData = phenoData,annotation = "GPL10558")
## Convert to ExpressionSetIllumina
summaryData = as(exprSet, "ExpressionSetIllumina")
## Calculate detection p-values with negative probes
fData(summaryData)$Status = ifelse(fData(summaryData)$PROBEQUALITY=="No match","negative","regular")
detectionPVal <- calculateDetection(summaryData,status=fData(summaryData)$Status)
## Clean environment
rm(phenoData,exprSet,summaryData)

## Keep probes with at least one sample called detected (p>0.05) based on estimated detection p values  
to_rm = as_tibble(t(detectionPVal)) %>% reshape2::melt(.) %>% 
  group_by(variable) %>% filter(value < 0.05) %>% 
  tally() %>% mutate(percent = n/N) %>%  filter(percent >0)
dat = dat[row.names(dat) %in% as.character(to_rm$variable),]
print(paste0("DIM: Detection p values: ",paste0(dim(dat),collapse = ", ")))

#### Filter bad quality probes
probes = row.names(dat)
qual=unlist(mget(probes, illuminaHumanv4.db::illuminaHumanv4PROBEQUALITY, ifnotfound=NA))
table(qual)
rem = qual == "No match" | qual == "Bad" | is.na(qual)
dat = dat[!rem,]
print(paste0("DIM: Bad quality probes: ",paste0(dim(dat),collapse = ", ")))

#### Filter X,Y chr probes
probes = row.names(dat)
chrom = mget(probes, illuminaHumanv4.db::illuminaHumanv4CHR, ifnotfound=NA)
table(lengths(chrom))
dat = dat[-which(lengths(chrom)>1),]
probes = row.names(dat)
chrom = unlist(mget(probes, illuminaHumanv4.db::illuminaHumanv4CHR, ifnotfound=NA))
rem = chrom == "X" | chrom == "Y" | chrom == "Un" | is.na(chrom)
table(rem)
dat = dat[!rem,]
detach(package:illuminaHumanv4.db)

print(paste0("DIM: XY filter: ",paste0(dim(dat),collapse = ", ")))

## Save 
saveName = paste0("AdNeuro/data/",format(Sys.Date(),"%d%m%Y"),"_",DATASETNAME,"_ComBat_age_sex_resid_filtered.txt")
fwrite(x = as.data.frame(dat), file = saveName,row.names = T,col.names = T)
