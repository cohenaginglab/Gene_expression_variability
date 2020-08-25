# Author: Frederik Dufour
# Title:  Step 1 - Sample Quality Control: clustering
# Description: Cluster ids to find problematic samples
## This file was written by combining all scripts for: parsing the metadata and doing sample quality control

###-------------------------------------------------------------------
## Libraries
library(data.table)
library(tidyverse)


###-------------------------------------------------------------------
## Get working directory
WORKDIR = getwd()

###-------------------------------------------------------------------
## Leipzig

## Read data 
FILENAME = paste0(WORDIR,"/Leipzig/data/GSE65907_non_normalized_47323probes.txt")
dat = fread(FILENAME,data.table = F) 
row.names(dat) = dat$ID_REF
dat = dat[,-1]

DetectionPval = dat[,seq(from = 2, to = ncol(dat), by = 2)]
dat = dat[,seq(from = 1, to = ncol(dat), by = 2)]


## Parse metadata
METAFILENAME = paste0(WORKDIR, "/Leipzig/metadata/LLHS_SampleDescription.txt")
Leipzig_carac = fread(METAFILENAME,data.table = F)

## Parse metadata
sample_title_order = Leipzig_carac$`Comment [Sample_title]`
sample_title_order = str_replace(string = sample_title_order,pattern = ".*[0-9]_",replacement = "")
sample_title_order = as.numeric(sample_title_order)
anyNA(sample_title_order)
COLNAMES = colnames(Leipzig_carac)
idx = str_detect(string = tolower(COLNAMES),pattern = "characteristics")
idx[1] = TRUE
COLNAMES = COLNAMES[idx]
Leipzig_carac = Leipzig_carac[,idx]
COLNAMES = str_replace(string = COLNAMES, pattern = ".*\\[",replacement = "") 
COLNAMES = str_replace(string = COLNAMES, pattern = "\\]",replacement = "") 
COLNAMES = str_replace(string = COLNAMES, pattern = " ",replacement = "_") 
colnames(Leipzig_carac)= COLNAMES
colnames(Leipzig_carac)[grep("age",COLNAMES)] = "age"
colnames(Leipzig_carac)[grep("sex",COLNAMES)] = "gender"
Leipzig_carac = Leipzig_carac[sample_title_order,]
Leipzig_carac$Source_Name = str_replace(string =Leipzig_carac$Source_Name ,pattern = " 1", "")
colnames(Leipzig_carac)[1] = "fileID"
Leipzig_carac$batchid = Leipzig_carac$sentrix_id


## I compared values on GEO and the sample order are the same between metadata and data
colnames(dat) = Leipzig_carac$fileID

## Check for missing age, sex or batch 
anyNA(Leipzig_carac$age)
anyNA(Leipzig_carac$gender)
anyNA(Leipzig_carac$sentrix_id)

## Check if there are any sample with age < 20
any(Leipzig_carac$age<20)

## Save metadata to RData
save(Leipzig_carac,file = "Leipzig/metadata/18022020_Leipzig_metadata_ordored_reanalysis.RData")

## QC data
## Scale data
dat_scaled = apply(dat[,-1],1,scale,T,T)

## Cluster
cl = hclust(dist(dat_scaled))
plot(cl)  

## Plot randomly selected probes
par(mfrow=c(2,5))    
for(i in sample(2:ncol(dat_scaled),10)){
  plot(as.numeric(dat[i,-1]), col = as.factor(Leipzig_carac$processing_batch))
}

####  Nothing really problematic.

## save expression data
fwrite(dat,file = paste0(WORDIR,"/Leipzig/data/GSE65907_non_normalized_47323probes_exprs.txt"),row.names = T,col.names = T)


###-------------------------------------------------------------------
### AddNeuro

## Clear working environment
rm(list = ls())
gc()
WORKDIR = getwd()
## Parse metadata
METAFILENAME = paste0(WORKDIR, "/AdNeuro/metadata/E-GEOD-63063.sdrf.txt")
meta = fread(METAFILENAME,data.table = F)
## Parse colnames
colnames(meta) = make.unique(colnames(meta))
COLNAMES = colnames(meta)
meta  = meta %>% select(one_of(c("Source Name",grep("Characteristics",COLNAMES,value = T),grep("Comment",COLNAMES,value = T)))) %>%
  mutate(`Source Name` = str_extract(`Source Name`,pattern = "GSM\\d*"),
         batchid = str_extract(`Comment [Sample_title]`,"\\d+")) %>% select(-`Comment [Derived ArrayExpress FTP file]`)
colnames(meta) = gsub("Characteristics|Comment|\\[|\\]","",colnames(meta))
colnames(meta) = gsub("^ ","",colnames(meta))
colnames(meta) = gsub(" ","_",colnames(meta))
colnames(meta) = tolower(colnames(meta))
colnames(meta) = gsub("-","",colnames(meta))
colnames(meta) = gsub("sex","gender",colnames(meta))
colnames(meta) = gsub("source_name","fileID",colnames(meta))

## Check missing and age minimum
any(meta$age<20)
anyNA(meta$age)
anyNA(meta$gender)
anyNA(meta$batchid)
any(meta$gender == "")
any(meta$batchid == "")

## Data downloaded from 2 GSE files measured on 2 platform
## Combine the data
dat1 = fread("data/GSE63060_non-normalized.txt",data.table = F)
dat2 = fread("data/GSE63061_non-normalized.txt",data.table = F)

dat = dat1 %>% left_join(dat2,.,by="IND_ID")
colnames(dat)[2:ncol(dat)]  = str_extract(colnames(dat)[2:ncol(dat)],"\\d+_*[:upper:]*") 
colnames(dat1)[2:ncol(dat1)]  = str_extract(colnames(dat1)[2:ncol(dat1)],"\\d+_*[:upper:]*") 
colnames(dat2)[2:ncol(dat2)]  = str_extract(colnames(dat2)[2:ncol(dat2)],"\\d+_*[:upper:]*") 


## Reorder data to match metadata
idx = match(meta$sample_title,colnames(dat))
idx = na.omit(idx)
dat = dat[,c(1,idx)]

idx1 = match(meta$sample_title,colnames(dat1))
idx1 = na.omit(idx1)
dat1 = dat1[,c(1,idx1)]

idx2 = match(meta$sample_title,colnames(dat2))
idx2 = na.omit(idx2)
dat2 = dat2[,c(1,idx2)]

## Reorder meta to match data
idx = match(colnames(dat),meta$sample_title)
idx = na.omit(idx)
meta = meta[idx,]

idx1 = match(colnames(dat1),meta$sample_title)
idx1 = na.omit(idx1)
meta1 = meta[idx1,]

idx2 = match(colnames(dat2),meta$sample_title)
idx2 = na.omit(idx2)
meta2 = meta[idx2,]


all.equal(meta$sample_title,colnames(dat)[2:ncol(dat)])
all.equal(meta1$sample_title,colnames(dat1)[2:ncol(dat1)])
all.equal(meta2$sample_title,colnames(dat2)[2:ncol(dat2)])

colnames(dat)[2:ncol(dat)] = meta$fileID
colnames(dat1)[2:ncol(dat1)] = meta1$fileID
colnames(dat2)[2:ncol(dat2)] = meta2$fileID

## QC data: Clustering
## Scale data
dat_scaled = apply(dat[,-1],1,scale,T,T)

## Cluster
cl = hclust(dist(dat_scaled))
plot(cl)  

## Two ids are possibly outliers: idx are 85 and 195 
COLORS = rep("black",nrow(dat_scaled))
COLORS[c(85,195)] = "red"

par(mfrow=c(2,6))    
for(i in sample(2:ncol(dat_scaled),10)){
  plot(as.numeric(dat[i,]), col = COLORS)
} 

## Clustering suggest they are outliers but looking at randomly selected probes, they are not obvious outliers in the data
## I will keep them in the analysis.
## There are 2 platforms in this dataset V3 and  V4 of the ilmn-HT12. The processing step needs to be adjusted to remove the effects of platform. 

## Create platform variable
meta = meta %>% mutate(platform = case_when(
  fileID %in% meta1$fileID ~ "V3",
  fileID %in% meta2$fileID ~ "V4",
))


## Save meta
fwrite(meta,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63063_metadata.txt"),col.names = T,row.names = F)
save(meta,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63063_metadata.RData"))
fwrite(dat,file = paste0(WORKDIR,"/AdNeuro/data/GSE63063_rawData.txt"),col.names = T,row.names = F)

fwrite(meta1,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63060_metadata.txt"),col.names = T,row.names = F)
save(meta1,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63060_metadata.RData"))
fwrite(dat1,file = paste0(WORKDIR,"/AdNeuro/data/GSE63060_rawData.txt"),col.names = T,row.names = F)

fwrite(meta2,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63061_metadata.txt"),col.names = T,row.names = F)
save(meta2,file = paste0(WORKDIR,"/AdNeuro/metadata/GSE63061_metadata.RData"))
fwrite(dat2,file = paste0(WORKDIR,"/AdNeuro/data/GSE63061_rawData.txt"),col.names = T,row.names = F)

## Check age distrivution
meta1 %>% ggplot(.,aes(age)) + geom_histogram()
meta2 %>% ggplot(.,aes(age)) + geom_histogram()




###-------------------------------------------------------------------
### GTP
## Clear working environment
rm(list = ls())
gc()
WORKDIR = getwd()

## read metadata
METAFILENAME = paste0(WORKDIR, "/TranscriptAgingHuman/metadata/E-GEOD-58137_metadata.txt") 
meta = fread(METAFILENAME,data.table = F)

## Parse metadata
colnames(meta) = make.unique(colnames(meta))
COLNAMES = colnames(meta)
meta = meta %>% select(one_of(c("Source Name",grep("Characteristics",COLNAMES,value = T),grep("Comment",COLNAMES,value = T),grep("FactorValue",COLNAMES,value = T)))) %>%
  mutate(`Source Name` = str_extract(`Source Name`,pattern = "GSM\\d*"),
         batchID = paste0(meta$`FactorValue [PLATEIDAMPLIFICATION]`,meta$`FactorValue [PLATEIDHYBRIDIZATION]`)) %>% 
  select(one_of(c("Source Name","batchID",grep("Characteristics",COLNAMES,value = T),grep("Comment",COLNAMES,value = T)))) %>% select(-`Comment [Derived ArrayExpress FTP file]`)

colnames(meta) = gsub("Characteristics|Comment|\\[|\\]","",colnames(meta))
colnames(meta) = gsub("^ ","",colnames(meta))
colnames(meta) = gsub(" ","_",colnames(meta))
colnames(meta) = tolower(colnames(meta))
colnames(meta) = gsub("-","",colnames(meta))
colnames(meta) = gsub("sex","gender",colnames(meta))
colnames(meta) = gsub("source_name","fileID",colnames(meta))

## Check for missing meta
any(meta$age<20)
anyNA(meta$age)
anyNA(meta$gender)
anyNA(meta$batchid)
any(meta$gender == "")
any(meta$batchid == "")
meta = meta %>% filter(age > 20)

## data 
dat  =  fread(paste0(WORKDIR, "/TranscriptAgingHuman/data/GSE58137_Raw_279_samplesremoved.csv"),data.table = F)
pval_matrix = dat %>% select(one_of(c("Probe_ID",grep("Detection",colnames(dat),value = T))))
dat = dat %>% select(-one_of(grep("Detection",colnames(dat),value = T)))

colnames(dat)[2:ncol(dat)] = str_extract(colnames(dat)[2:ncol(dat)],"\\d+")

## Keep only the matching samples in meta and data
idx = match(as.character(meta$sample_title),colnames(dat))
idx = na.omit(idx)
dat = dat[,c(1,idx)]


idx = match(colnames(dat),as.character(meta$sample_title))
idx = na.omit(idx)
meta = meta[idx,]

all.equal(as.character(meta$sample_title),colnames(dat)[2:ncol(dat)])
colnames(dat)[2:ncol(dat)] = meta$fileID

## QC clustering
dat_scaled = apply(dat[,-1],1,scale,T,T)
cl = hclust(dist(dat_scaled),method = "single")
plot(cl)  


COLORS = rep("black",nrow(dat_scaled))
COLORS[c(61)] = "red"

for(i in sample(2:ncol(dat_scaled),10)){
  temp_dat=as.numeric(dat[i,-1])
  plot(temp_dat, col = COLORS)
  
} 

par(mfrow=c(2,5))    
for(i in sample(2:ncol(dat_scaled),10)){
  temp_dat=as.numeric(dat[i,-1])
  plot(temp_dat[which(!is.na(temp_dat))], col = as.factor(meta$plateidhybridization)[which(!is.na(temp_dat))])
  
} 

## Nothing really strange

## Write data
fwrite(meta,file = paste0(WORKDIR, "/TranscriptAgingHuman/metadata/GSE58137_metadata.txt"),col.names = T,row.names = F)
save(meta,file = paste0(WORKDIR, "/TranscriptAgingHuman/metadata/GSE58137_metadata.RData"))

fwrite(dat,file = paste0(WORKDIR, "/TranscriptAgingHuman/data/GSE58137_RawData.txt"),col.names = T,row.names = F)



###-------------------------------------------------------------------
####### CHDWB
## Clear  
rm(list = ls())
gc()
WORKDIR = getwd()

METAFILENAME = paste0(WORKDIR,"/CHDWB/metadata/GSE61672_metadata.txt")
meta = fread(METAFILENAME)

## Count number of samples with age
meta %>% select(age = grep("age",colnames(meta))) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% nrow(.)

## N age >= 20
meta %>% select(age = grep("age",colnames(meta))) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% filter(age>=20) %>% nrow(.)

## Age distribution
meta %>% select(age = grep("age",colnames(meta))) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% filter(age>=20) %>% ggplot(.,aes(x = age)) + geom_histogram()


## parse and filter meta
meta = meta %>% filter(age>=20)
colnames(meta) = tolower(colnames(meta))
colnames(meta) = gsub(" ","_",colnames(meta))
colnames(meta) = gsub("sex","gender",colnames(meta))
colnames(meta) = gsub("gsm","fileID",colnames(meta))
colnames(meta) = gsub("hybridization_batch","batchid",colnames(meta))

## Check missing and age minimum
any(meta$age<20)
anyNA(meta$age)
anyNA(meta$gender)
anyNA(meta$batchid)
any(meta$gender == "")
any(meta$batchid == "")


## Some samples with missing sex variable
meta %>% group_by(gender) %>% tally()
meta = meta %>% filter(gender != "")

## get data 
FILENAME = paste0(WORKDIR,"/CHDWB/data/GSE61672_non-normalized.txt")
dat = fread(FILENAME,data.table = F)
p_values = dat[,grep("Detection",colnames(dat))]
dat_expr = dat[,grep("Detection",colnames(dat),invert = T)]
row.names(dat_expr) = dat$ID_REF
row.names(p_values) = dat$ID_REF
dat_expr = dat_expr[,-1]

META_COLNAMES_ID = str_extract(meta$title,"GG\\d_\\d+") %>% as.character(.)

idx = match(META_COLNAMES_ID, colnames(dat_expr))
dat_expr = dat_expr[,idx]
p_values = p_values[,idx]

colnames(dat_expr) = meta$fileID
colnames(p_values) = meta$fileID

all.equal(meta$fileID,colnames(dat_expr))

## QC: Clustering
## Scale data
dat_scaled = apply(dat[,-1],1,scale,T,T)

## Cluster
cl = hclust(dist(dat_scaled),method = "single")
plot(cl)  

## Maybe Idx 193 is problematic? 
COLORS = rep("black",nrow(dat_scaled))
COLORS[c(193)] = "red"

par(mfrow=c(2,5))    
for(i in sample(2:ncol(dat_scaled),10)){
  plot(as.numeric(dat[i,-1]), col = COLORS)
  
} 
## DATA is ok


## Save data
fwrite(meta,file = paste0(WORKDIR,"/CHDWB/metadata/GSE61672_metadata.txt"),col.names = T,row.names = F)
save(meta,file = paste0(WORKDIR,"/CHDWB/metadata/GSE61672_metadata.RData"))

fwrite(dat_expr,file = paste0(WORKDIR,"/CHDWB/data/GSE61672_expr.txt"),col.names = T,row.names = T)
fwrite(p_values,file = paste0(WORKDIR,"/CHDWB/data/GSE61672_pvalues.txt"),col.names = T,row.names = T)


###-------------------------------------------------------------------
### UKPSSR
rm(list = ls())
gc()
WORKDIR = getwd()

## Read meta
METAFILENAME = paste0(WORKDIR,"/UKPSSR/metadata/E-MTAB-8277.sdrf.txt")
meta = fread(METAFILENAME)
colnames(meta) = make.unique(colnames(meta))


## Count N samples
meta %>% select(age = grep("age",colnames(meta))[2]) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% nrow(.)

## N age >= 20
meta %>% select(age = grep("age",colnames(meta))[2]) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% filter(age>=20) %>% nrow(.)

## Age distributions
meta %>% select(age = grep("age",colnames(meta))[2]) %>% mutate(age = str_extract(age,"\\d+")) %>% mutate_all(as.numeric) %>% filter(age>=20) %>% ggplot(.,aes(x = age)) + geom_histogram()

## Filter and parse meta
COLNAMES = colnames(meta)
meta  = meta %>% select(one_of(c("Source Name",grep("Characteristics",COLNAMES,value = T),grep("Comment",COLNAMES,value = T)))) %>% select(-`Comment [ArrayExpress FTP file]`) %>%
  mutate(batchid = str_extract(`Source Name`,"\\d+")) 

colnames(meta) = gsub("Characteristics|Comment|\\[|\\]","",colnames(meta))
colnames(meta) = gsub("^ ","",colnames(meta))
colnames(meta) = gsub(" ","_",colnames(meta))
colnames(meta) = tolower(colnames(meta))
meta = meta %>% mutate(age = as.numeric(age))

colnames(meta) = gsub("sex","gender",colnames(meta))
colnames(meta) = gsub("source_name","fileID",colnames(meta))


## Check missing and age minimum
any(meta$age<20)
anyNA(meta$age)
anyNA(meta$gender)
anyNA(meta$batchid)
any(meta$gender == "")
any(meta$batchid == "")


## get data from idat files
Annotation_bgx_file = paste0(WORKDIR,"Data_Analysis/HT12_probeInformation/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt")

FILEDIR = paste0(WORKDIR, "/data/E-MTAB-8277/")
Files_to_load = paste0(FILEDIR, list.files(FILEDIR,pattern = "idat"))
dat = limma::read.idat(Files_to_load,bgxfile = Annotation_bgx_file)
dat_E = dat$E
dat = dat_E
rm(dat_E)
gc()


idx = na.omit(match(row.names(dat), dat$genes$Array_Address_Id))
ilmn_probes_id = dat$genes$Probe_Id[idx]
row.names(dat) = ilmn_probes_id

colnames(dat) = str_extract(colnames(dat),"\\d+_[[:alpha:]]")

idx = match(meta$fileID,colnames(dat))
idx = na.omit(idx)
dat = dat[,idx]
all.equal(meta$fileID,colnames(dat))

## Modification done after trying to process ComBat
## Reading the file adds an X to sample names because the names are numbers
meta = meta %>% mutate(fileID = paste0("X",.$fileID))

## QC: Clustering
dat_scaled = apply(dat[,-1],1,scale,T,T)
cl = hclust(dist(dat_scaled),method = "single")
plot(cl)  

## These 2 ids seem odd.
COLORS = rep("black",nrow(dat_scaled))
COLORS[c(125,165)] = "red"

par(mfrow=c(2,5))    
for(i in sample(2:ncol(dat_scaled),10)){
  temp_dat=as.numeric(dat[i,-1])
  plot(temp_dat, col = COLORS)
  
} 
## Signal is 0 for many probes. Sometime there is signal but it seems to be background signal
## I will remove these 2 samples.

dat[,c(126,166)] = NULL
dat_scaled = dat_scaled[-c(125,165),]
meta = meta[-c(125,165),]

## Replot to look at data
par(mfrow=c(2,5))    
for(i in sample(2:ncol(dat_scaled),10)){
  temp_dat=as.numeric(dat[i,-1])
  plot(temp_dat, col = as.factor(meta$batchid))
  
} 

## Save data
fwrite(meta,file = paste0(WORKDIR,"/UKPSSR/metadata/EMTAB8277_metadata.txt"),col.names = T,row.names = F)
save(meta,file = paste0(WORKDIR,"/UKPSSR/metadata/EMTAB8277_metadata.RData"))

fwrite(as.data.frame(dat),file = paste0(WORKDIR,"/UKPSSR/data/EMTAB8277_expr.txt"),col.names = T,row.names = T)


