  # Author: Frederik Dufour
# Date: Mon Feb 17 10:27:37 2020
# Title: Calculate DM  
# Description:  After preprocessing, calculate DM and correlate to age
#------------------------------------------------------
# init
#------------------------------------------------------
if(length(commandArgs(TRUE))!=2) stop("########################### 
ERROR: Command should be: Rscript Script.R Parameter_file SystemFile.
############################
The parameter file is a CSV file that contains the following:
col1 = data;
col2 = reference; 
col3 = metadata[data]
The System file is a tab separated file that is read line-by-line and the first entry in the line is the GO term or the system's name. The rest of the line lists the ILMN id in the system.
Ex: GO:0000001  ILMN1010101010  ILMN267262626 ILMN6263537
    GO:000002 ILMN2263636 ILMN378383874 ILMN237843  ILMN27237637  ILMN437387")

## libraries
libs = list("data.table","ggpubr","tidyverse","beadarray","broom")
lapply(libs,require,character.only = T)

## define functions
### Define reference population: returns vector of means, inverse covariance matrix, and vector of SD
define_ref_population = function(expression_data,keepAll = T ,condition = "", N_sd_call_outlier = 3, idx = 1:N) {
  reference_population = expression_data
  # remove outliers N_sd from mean
  if(!keepAll){
    if(condition == "deviation_call"){
      for(i in 1:NCOL(expression_data)){
        probe_data =  expression_data[,i]
        median_probe_data = median(probe_data,na.rm = T)
        sd_probe_data =  sd(probe_data,na.rm = T)
        outliers = which(probe_data > (median_probe_data + (N_sd_call_outlier * sd_probe_data)) | probe_data < (median_probe_data - (N_sd_call_outlier * sd_probe_data)))
        reference_population[outliers,i] = NA
      }
    }
    
    if(condition == "idx"){
      reference_population = reference_population[idx,]
    }
  }
  return(reference_population)
}

#function will center and reduce based on reference population, and return a mahalanobis distance per ID
calculate.mahalanobis.distance = function(system_data,reference_population){
  #Standerdized data by population reference population's mean and sd
  col_medians = apply(reference_population,2, function(xx) median(xx,na.rm = T))
  col_sd = apply(reference_population,2,function(xx) sd(xx,na.rm = T))
  reference_population = t(apply(reference_population,1,function(dat) {(dat - col_medians) / col_sd ** as.logical(col_sd)}))
  system_data = t(apply(system_data,1,function(dat) {(dat - col_medians) / col_sd ** as.logical(col_sd)}))
  #Define covariance matrix with pairwise.complete and invert the matrix from reference population
  COV = cov(reference_population,use = "pairwise.complete")
  inv_ref_covariance_matrix = chol2inv(chol(COV))
  #calculate distance returns NA if any missing value by default
  col_medians = apply(reference_population,2, function(xx) median(xx,na.rm = T))
  D2 = mahalanobis(system_data,col_medians, inv_ref_covariance_matrix, inverted = T)
  D = sqrt(D2)
  return(D)
}

## Read systems
systems_Filename = commandArgs(TRUE)[2]
systems = scan(file = systems_Filename,what = "character",sep = "\n")
systems = lapply(systems, function(xx) return(unlist(strsplit(xx,split = "\t"))))
for(i in 1:length(systems)){
  names(systems)[i] = systems[[i]][1]
}
systems = lapply(systems, function(xx) return(xx[-1]))


main = function(Line){
  ## Files and directories
  dat_Filename = Line[,1]
  meta_Filename = Line[,3]
 
  ## read data
  dat = data.frame(fread(file = dat_Filename,data.table = F),row.names = 1)
  # ----------------------------------------------------------------------------
  
  ## Load metadata
  obj_name = load(meta_Filename)
  assign("metadat",get(obj_name))
  
  idx = na.omit(match(metadat$fileID,colnames(dat)))
  dat = dat[,idx]
  idx = na.omit(match(colnames(dat),metadat$fileID))
  metadat = metadat[idx,]
  
  probes = row.names(dat)
  N = ncol(dat)
  
  # female_idx = grep("[f,F]", metadat$gender)
  # n_females = length(female_idx)
  # n_males = N - n_females
  
  DMs = vector("list",length(systems))
  # DMs_females = vector("list",length(systems)) 
  # DMs_males = vector("list",length(systems))
  
  for(group in names(systems)){
    i = grep(paste0("^",group,"$"),names(systems))
    sys_probes = systems[[i]]
    ## match probes to data
    idx = match(sys_probes,table = probes)
    idx = idx[!is.na(idx)]
    
    ## set temporary data
    dat2 = t(dat[idx,])
    if(length(idx)!=1){
      ###########################
      # Reference population
      ###########################
      ## Define reference population based on all dat2a
      reference_pop = define_ref_population(expression_data = dat2,keepAll = T)
      # reference_pop_females = define_ref_population(expression_data = dat2[female_idx,],keepAll = T)
      # reference_pop_males = define_ref_population(expression_data = dat2[-female_idx,],keepAll = T)
      
      ##########################
      # Mahalanobis distances
      ##########################
      distance = try(calculate.mahalanobis.distance(dat2,reference_pop),silent = T)
      # distance_females = try(calculate.mahalanobis.distance(dat2[female_idx,],reference_pop_females),silent = T)
      # distance_males = try(calculate.mahalanobis.distance(dat2[-female_idx,],reference_pop_males),silent = T)
      
      # If distance cannot be calculated because of a singular covariance matrix, return list of NAs
      if(class(distance)=="try-error"){
        distance = rep(NA,N)
      }
      # 
      # if(class(distance_females)=="try-error"){
      #   distance_females = rep(NA,n_females)
      # }
      # 
      # if(class(distance_males)=="try-error"){
      #   distance_males = rep(NA,n_males)
      # }
      
      names(DMs)[i] = group
      # names(DMs_females)[i] = group
      # names(DMs_males)[i] = group
      
      DMs[[i]] = log(distance)
      # DMs_females[[i]] = log(distance_females)
      # DMs_males[[i]] = log(distance_males)
    }else{
      distance = rep(NA,N)
      # distance_females = rep(NA,n_females)
      # distance_males = rep(NA,n_males)
      
      names(DMs)[i] = group
      # names(DMs_females)[i] = group
      # names(DMs_males)[i] = group
      
      DMs[[i]] = distance
      # DMs_females[[i]] = distance_females
      # DMs_males[[i]] = distance_males
      
    }
    
    if(i%%500 == 0) print(i)
    
  }
  
  DMs = do.call(rbind,DMs)
  # DMs_females = do.call(rbind,DMs_females)
  # DMs_males = do.call(rbind,DMs_males)
  
  ## Remove AllNAs
  N_NAs = apply(DMs,1,function(xx) sum(is.na(xx)))
  DMs = DMs[which(N_NAs != N),]
  
  # N_NAs = apply(DMs_females,1,function(xx) sum(is.na(xx)))
  # DMs_females = DMs_females[which(N_NAs != n_females),]
  # 
  # N_NAs = apply(DMs_males,1,function(xx) sum(is.na(xx)))
  # DMs_males = DMs_males[which(N_NAs != n_males),]
  
  
  ### Correlations
  Corrdat = as_tibble(t(DMs)) %>% melt(value.name = "DM" , variable.name = "GO") %>%
    group_by(GO) %>% mutate(age = metadat$age) %>% do(Fit = cor.test(x = .$age,.$DM)) %>% tidy(Fit,conf.int = T) %>% 
    ungroup() 
  
  
  # Corrdat_females = as_tibble(t(DMs_females)) %>% melt(value.name = "DM" , variable.name = "GO") %>%
  #   group_by(GO) %>% mutate(age = metadat$age[female_idx]) %>% do(Fit = cor.test(x = .$age,.$DM)) %>% tidy(Fit,conf.int = T) %>% 
  #   ungroup()
  # 
  # Corrdat_males = as_tibble(t(DMs_males)) %>% melt(value.name = "DM" , variable.name = "GO") %>%
  #   group_by(GO) %>% mutate(age = metadat$age[-female_idx]) %>% do(Fit = cor.test(x = .$age,.$DM)) %>% tidy(Fit,conf.int = T) %>% 
  #   ungroup()
  
  
  DM_saveName = gsub(pattern = "/data/\\d+" ,replacement = paste0("/mahalanobis_distances/",format(Sys.Date(),"%d%m%Y")),dat_Filename)
  DM_saveName = gsub("filtered","ref_all_median_mahalanobis_distances",DM_saveName)
  
  Corr_saveName = gsub(pattern = "/data/\\d+" ,replacement = paste0("/correlation_age/",format(Sys.Date(),"%d%m%Y")),dat_Filename)
  Corr_saveName = gsub("filtered","ref_all_median_age_dm_correlations",Corr_saveName)
  
  # DM_saveName_females = gsub("sex_resid_ref_all_median_mahalanobis_distances.txt","females_ref_all_median_mahalanobis_distances.txt",DM_saveName)
  # DM_saveName_males = gsub("sex_resid_ref_all_median_mahalanobis_distances.txt","males_ref_all_median_mahalanobis_distances.txt",DM_saveName)
  # 
  # Corr_saveName_females = gsub("sex_resid_ref_all_median_age_dm_correlations.txt","females_ref_all_median_age_dm_correlations.txt",Corr_saveName)
  # Corr_saveName_males = gsub("sex_resid_ref_all_median_age_dm_correlations.txt","males_ref_all_median_age_dm_correlations.txt",Corr_saveName)
  
  
  fwrite(as.data.frame(DMs),file = DM_saveName,col.names = T,row.names = T)
  fwrite(Corrdat,file = Corr_saveName,col.names = T,row.names = F)
  
  # fwrite(as.data.frame(DMs_females),file = DM_saveName_females,col.names = T,row.names = T)
  # fwrite(Corrdat_females,file = Corr_saveName_females,col.names = T,row.names = F)
  # 
  # fwrite(as.data.frame(DMs_males),file = DM_saveName_males,col.names = T,row.names = T)
  # fwrite(Corrdat_males,file = Corr_saveName_males,col.names = T,row.names = F)
}

params = fread(file = commandArgs(T)[1], data.table = F,header = T)

for(i in 1:nrow(params)){
  main(params[i,])
}