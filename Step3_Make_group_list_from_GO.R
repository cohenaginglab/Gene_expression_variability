# Author: Frederik Dufour
# Date: Fri Dec 20 11:02:35 2019
# Title:  Reduce the number of GO terms tested
# Description: Let's try to reduce the number of systems tested in order to have more independant systems tested.
#------------------------------------------------------
# init
#------------------------------------------------------
rm(list = ls())
## set working dir
setwd("~/")
## libraries
libs = list("data.table","ggpubr","tidyverse","mgsa","ontologyIndex", "illuminaHumanv4.db","RBGL","graph","GO.db")
lapply(libs,require,character.only = T)
#------------------------------------------------------
## functions
transform.line = function(linn){
  probes = strsplit(linn,split = "\t")
  probes = probes[[1]][grep(pattern = "[[:digit:]+]",x = probes[[1]])]
  probes = unique(tolower(probes))
}
## read obo file
GO_index = get_ontology("24102019_go_basic.obo",extract_tags = "everything")
## read gaf file with mappings
Mapping = readGAF("20191007_release_goa_human.gaf.gz",evidence = c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP","ISS","ISA","ISM","ISO","IBA","IBD"))
## Set Illumina to symbol mapping
symbol2ILMN  = as.list(illuminaHumanv4SYMBOLREANNOTATED)
#####################################################################################

## read in available probes in datasets
pathToData = "~/path_to_data"
extension = "_ComBat_age_sex_resid_filtered.txt"

  datasets = c(paste0(pathToData,"AdNeuro/data/09062020_AddNeuro_Combined", extension),
               paste0(pathToData,"CHDWB/data/09062020_CHDWB", extension),
               paste0(pathToData,"Leipzig/data/09062020_Leipzig", extension),
               paste0(pathToData,"TranscriptAgingHuman/data/09062020_TranscriptAgingHuman", extension),
               paste0(pathToData,"UKPSSR/data/09062020_UKPSSR", extension))
  
  ## read in available probes in datasets
  for(i in 1:length(datasets)){
    assign(paste0("impossible_name",i),fread(datasets[i],select = 1)$V1)
    assign(paste0("impossible_name_2",i), fread(datasets[i],nrows = 1))
  }
  

## Min-max group size
min_size = 10
max_size = min(unlist(lapply(mget(ls(pattern = "impossible_name_2\\d+")), length)))-1
rm(list = ls(pattern = "impossible_name_2\\d+"))

probe_data = Reduce(intersect,mget(ls(pattern = "impossible_name")))

rm(list = ls(pattern = "impossible_name")); gc()

COLNAMES= c("GO", "keep","N")
GO_group_info = as.data.frame(matrix(NA,length(GO_index$id),length(COLNAMES)))
colnames(GO_group_info) = COLNAMES
GO_group_info$GO = GO_index$id
GO_group_info$keep = F


## Map genes to GO_terms and filter on size of group between 10-30
i = 0
for(go_id in GO_index$id){
  #go_id = GO_index$id[2] # troubleshoot
  i = i+1
  ## map genes to probes 
  genes_symbol = try(as.character(Mapping@itemAnnotations[unlist(Mapping@sets[[grep(pattern = go_id, x = names(Mapping@sets))]]),1]),silent = T)
  ## Filter out GO with less than 3 genes in the group
  if(length(genes_symbol) < 3 | class(genes_symbol) == "try-error") next
  GO_probes = names(symbol2ILMN)[which(symbol2ILMN %in% genes_symbol)]
  ## keep only probes available in data
  GO_probes = GO_probes[GO_probes %in% probe_data]
  ## IF N probes greater than min(obs_Leipzig,obs_RS3)
  if(length(GO_probes) > max_size | length(GO_probes) < min_size) next
  ## Backtransform to gene and filter if less than 3 genes included
  genes_symbol2 = unlist(symbol2ILMN[which(names(symbol2ILMN) %in% GO_probes)])
  if(length(genes_symbol2) < 3) next
  ## write to file 
  GO_group_info$keep[i] = T
  GO_group_info$N[i] = length(GO_probes)
  ##### END loop 1
  if(i%%1000==0) print(i / length(GO_index$id) *100)
}

## Stats 
summary(GO_group_info)
go_keep_ids = GO_group_info$GO[which(GO_group_info$keep == T)]

## Keep biological processes only
getDistanceToRoot = function(GO_term,GOPARENTS_graph,GOANCESTOR){
  SUBGRAPH = subGraph(c(GO_term,get(GO_term, GOANCESTOR)), GOPARENTS_graph)
  DISTANCE = dijkstra.sp(SUBGRAPH, GO_term)$distances["all"]
  return(DISTANCE)
}

#biological process
BIOproc_parentTable = toTable(GOBPPARENTS)
BIOproc_GOtermGraph = ftM2graphNEL(as.matrix(BIOproc_parentTable[, 1:2]), W=rep(1,dim(BIOproc_parentTable)[1]))

## get distance to root
DISTANCE_to_ROOT = lapply(go_keep_ids, function(GO_term) tryCatch(getDistanceToRoot(GO_term ,BIOproc_GOtermGraph,GOBPANCESTOR)[[1]],error = function(e) return(NA)))
DISTANCE_to_ROOT = unlist(DISTANCE_to_ROOT)

## set keep to false for all GOs with distance = NA
GO_group_info$keep[which(GO_group_info$keep)][is.na(DISTANCE_to_ROOT)] = F
DISTANCE_to_ROOT = na.omit(DISTANCE_to_ROOT)
go_keep_ids = GO_group_info$GO[which(GO_group_info$keep)][order(DISTANCE_to_ROOT,decreasing = T)]
DISTANCE_to_ROOT[order(DISTANCE_to_ROOT,decreasing = T)]
    

## Function to get probes in groups
get_probes <- function(go_id) {
  genes_symbol = try(as.character(Mapping@itemAnnotations[unlist(Mapping@sets[[grep(pattern = go_id, x = names(Mapping@sets))]]),1]),silent = T)
  ## Filter out GO with less than 3 genes in the group
  if(length(genes_symbol) < 3 | class(genes_symbol) == "try-error") return(NA)
  GO_probes = names(symbol2ILMN)[which(symbol2ILMN %in% genes_symbol)]
  ## keep only probes available in data
  GO_probes = GO_probes[GO_probes %in% probe_data]
  return(GO_probes)
}  

  
## Open a connection to a file to save lists of probes
writer_file = "GOterms2probeListFromComBat_V3.txt"
writer = file(writer_file,open="w")
for(go_id in GO_group_info$GO[GO_group_info$keep]){
  cat(go_id,paste(get_probes(go_id),sep = "\t"),"\n",file =writer,sep ="\t",append = T)
}
close(writer)
        
