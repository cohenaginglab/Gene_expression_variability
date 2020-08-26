  # Author: Frederik Dufour
# Date: Sat Apr 11 10:50:28 2020
# Title: Data comparaison with 3 datasets processed with raw
# Description:
# Compare Leipzig, AdNeuro, RS and Grady Trauma Project: log2 transformed data with residualization for sex
#------------------------------------------------------
# init
#------------------------------------------------------
rm(list = ls())
## libraries
libs = list("data.table", "ggpubr", "tidyverse", "grid")
lapply(libs, require, character.only = T)
#------------------------------------------------------
#Working directories
## set working dir
setwd("~/")
# directory to save the figures
figure_save_directory = "~/"
# Directory where all data is stored  
data_directory = "~/Documents/"

#------------------------------------------------------
## Load data
#------------------------------------------------------
## Dataset Short names
shortName = c("Leipzig",
              "AddNeuro",
              "CHDWB",
              "TranscriptAgingHuman",
              "UKPSSR")

extension = function(shortName)
  return(
    paste0(
      "_age_sex_resid_ref_",
      shortName,
      "_median_age_dm_correlations.txt"
    )
  )


ComBat_datasets_files = c(
  paste0(
    data_directory,
    "Leipzig/correlation_age/19082020_Leipzig_ComBat",
    extension("all")
  ),
  paste0(
    data_directory,
    "AdNeuro/correlation_age/19082020_AddNeuro_Combined_ComBat",
    extension("all")
  ),
  paste0(
    data_directory,
    "CHDWB/correlation_age/19082020_CHDWB_ComBat",
    extension("all")
  ),
  paste0(
    data_directory,
    "TranscriptAgingHuman/correlation_age/19082020_TranscriptAgingHuman_ComBat",
    extension("all")
  ),
  paste0(
    data_directory,
    "UKPSSR/correlation_age/19082020_UKPSSR_ComBat",
    extension("all")
  )
)


ComBat_datasets = gsub(
  paste0(
    data_directory,
    "[[:alpha:]]+_*[[:alpha:]]*/correlation_age/\\d+_"
  ),
  "",
  ComBat_datasets_files
)
ComBat_datasets = gsub(extension("[[:alnum:]]+_*[[:alnum:]]*"), "", ComBat_datasets)
ComBat_datasets = gsub("TranscriptAgingHuman", "GTP", ComBat_datasets)
ComBat_datasets = gsub("_Combined", "", ComBat_datasets)


#### UKPSSR change name
ComBat_datasets_files[5] = str_replace(ComBat_datasets_files[5], ".txt", "_2id_rm.txt")

## Load correlations data
ComBat_Cor_objects_names = paste0(ComBat_datasets, "_Cor")
for (i in 1:length(ComBat_Cor_objects_names)) {
  assign(ComBat_Cor_objects_names[i],
         fread(ComBat_datasets_files[i], data.table = F))
}

## Joins correlations to JoinCorrelations object
JoinCorrelations = Leipzig_ComBat_Cor  %>% 
  dplyr::select(GO, estimate, p.value) %>%
  setNames(nm = c(
    "GO",
    paste0(ComBat_datasets[1], "_estimate"),
    paste0(ComBat_datasets[1], "_p.value")
  )) %>%
  full_join(
    AddNeuro_ComBat_Cor %>% #filter(term == "age") %>% 
      dplyr::select(GO, estimate, p.value) %>% setNames(nm = c(
      "GO",
      paste0(ComBat_datasets[2], "_estimate"),
      paste0(ComBat_datasets[2], "_p.value")
    )),
    by = "GO"
  )


for (i in 3:length(ComBat_Cor_objects_names)) {
  JoinCorrelations = JoinCorrelations %>%
    full_join(
      get(ComBat_Cor_objects_names[i]) %>% #filter(term == "age") %>% 
        dplyr::select(GO, estimate, p.value) %>% setNames(nm = c(
        "GO",
        paste0(ComBat_datasets[i], "_estimate"),
        paste0(ComBat_datasets[i], "_p.value")
      )),
      by = "GO"
    )
}

summary(JoinCorrelations)

#------------------------------------------------------
## Load metadata 
#------------------------------------------------------
metadata_filenames = c(
  paste0(data_directory,"Leipzig/metadata/18022020_Leipzig_metadata_ordored_reanalysis.RData"),
  paste0(data_directory,"AdNeuro/metadata/GSE63063_metadata.RData"),
  paste0(data_directory,"CHDWB/metadata/GSE61672_metadata.RData"),
  paste0(data_directory,"TranscriptAgingHuman/metadata/GSE58137_metadata.RData"),
  paste0(data_directory,"UKPSSR/metadata/EMTAB8277_metadata.RData")
)

Carac_objects_names = gsub("ComBat", "Carac", ComBat_datasets)
for (i in 1:length(metadata_filenames)) {
  obj_savename = load(metadata_filenames[i])
  assign(Carac_objects_names[i], value = get(obj_savename))
  rm(list = obj_savename)
}

#------------------------------------------------------
### Table 1
#------------------------------------------------------
get_carac = function(xx) {
  get(xx) %>%
    dplyr::select(age, gender) %>%
    filter(!is.na(age) & !is.na(gender)) %>%
    summarise(
      N = n(),
      median_age = round(median(age), 1),
      age_min = floor(min(age)),
      age_max = floor(max(age)),
      n_females = sum(grepl("f", tolower(gender))),
      females_percent = round(n_females / N * 100, 2)
    )
}

# Create the table
carac_table = Carac_objects_names %>%
  map(get_carac) %>%
  reduce(bind_rows) %>%
  mutate(dataset = str_remove(Carac_objects_names, "_Carac"))
# View
carac_table
## Write to file
write.table(
  carac_table,
  file = paste0(figure_save_directory, "../tables/table1_dataset_caracteristics.txt"),
  row.names = T,
  col.names = T,
  quote = F,
  sep = ","
)

#------------------------------------------------------
### Supp Figure 1: age distributions 
#------------------------------------------------------
## Join metadata
JoinAgeSex = get(Carac_objects_names[1]) %>% dplyr::select(fileID, age, gender) %>%
  mutate(
    dataset = str_remove(Carac_objects_names[1], "_Carac"),
    sex = ifelse(
      grepl("[f,F]", gender),
      "females",
      ifelse(grepl("[m,M]", gender), "males", "NA")
    )
  )

for (Carac_objects_name in Carac_objects_names[-1]) {
  JoinAgeSex = get(Carac_objects_name) %>% dplyr::select(fileID, age, gender) %>%
    mutate(
      dataset = str_remove(Carac_objects_name, "_Carac"),
      sex = ifelse(
        grepl("[f,F]", gender),
        "females",
        ifelse(grepl("[m,M]", gender), "males", "NA")
      )
    ) %>%
    bind_rows(JoinAgeSex, .)
}

## PLot
age_sex_distributions = JoinAgeSex %>% 
  ggplot(., aes(age, fill = sex)) + 
  geom_histogram(bins = 50) + 
  facet_wrap( ~ factor(dataset, levels = str_remove(Carac_objects_names, "_Carac")), scales = "free_y") + 
  theme_bw() + 
  scale_fill_brewer("Sex", palette = "Set1") + 
  theme(text = element_text(size = 15)) + 
  labs(x = "Age", y = "Count")

## Export to pdf and png
png_savename = paste0(figure_save_directory, "distributions_age_sex.png")
pdf_savename = paste0(figure_save_directory, "distributions_age_sex.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename, width = fig_width, height = fig_height)
age_sex_distributions
dev.off()

png(file = png_savename, width = fig_width, height = fig_height, units = "in", res = 300)
age_sex_distributions
dev.off()


#------------------------------------------------------
## Plot correlation distributions
#------------------------------------------------------
corr_density_plot = JoinCorrelations %>% select_at(vars(contains("ComBat_estimate"))) %>%
  setNames(stringr::str_remove_all(names(.), '_ComBat_estimate')) %>% reshape2::melt(.) %>%
  ggplot(., aes(value, fill = variable)) + geom_density(alpha = 0.7) + scale_fill_brewer("Dataset", palette = "Set2") +
  labs(y = "Density", x = "Correlation coefficients") +
  theme(axis.title = element_text(size = 15))
corr_density_plot

## Export to pdf and png
png_savename = paste0(figure_save_directory, "correlation_distributions_density.png")
pdf_savename = paste0(figure_save_directory, "correlation_distributions_density.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename, width = fig_width, height = fig_height)
corr_density_plot
dev.off()

png(file = png_savename, width = fig_width, height = fig_height, units = "in", res = 300)
corr_density_plot
dev.off()


#------------------------------------------------------
## Plot p-value distributions
#------------------------------------------------------
pval_hist_plot = JoinCorrelations %>% 
  select_at(vars(contains("ComBat_p.value"))) %>%
  setNames(stringr::str_remove_all(names(.), '_ComBat_p.value')) %>% 
  reshape2::melt(.) %>%
  ggplot(., aes(value)) + geom_histogram(bins = 20) + labs(y = "Counts", x = "p value") + facet_wrap( ~variable)

## Export to pdf and png
png_savename = paste0(figure_save_directory, "pval_distributions_histogram.png")
pdf_savename = paste0(figure_save_directory, "pval_distributions_histogram.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename, width = fig_width, height = fig_height)
pval_hist_plot
dev.off()

png(file = png_savename, width = fig_width, height = fig_height, units = "in", res = 300)
pval_hist_plot
dev.off()


#------------------------------------------------------
## Adjust p values for multiple comparaison
JoinCorrelations = JoinCorrelations %>% select_at(vars(contains("ComBat_p.value"))) %>% 
  transmute_all(p.adjust, method = "BH", n = nrow(Leipzig_ComBat_Cor)) %>%
  setNames(stringr::str_replace_all(names(.), 'p.value', 'q.value')) %>% 
  mutate(GO = JoinCorrelations$GO) %>%
  left_join(JoinCorrelations, ., by = "GO") %>% 
  as_tibble(.)

## Get GO terms and ontology
JoinCorrelations$ONTOLOGY = AnnotationDbi::select(GO.db::GO.db, keys = JoinCorrelations$GO, columns = "ONTOLOGY")$ONTOLOGY
JoinCorrelations$TERM = AnnotationDbi::select(GO.db::GO.db, keys = JoinCorrelations$GO, columns = "TERM")$TERM

## Make sure only the BP were analysed
all(JoinCorrelations$ONTOLOGY == "BP")


## Count significant correlations with p values and ComBat
JoinCorrelations %>% select_at(vars(contains("ComBat_p.value"))) %>% setNames(str_remove_all(names(.), "_p.value")) %>% reshape2::melt(.) %>% filter(value < 0.05) %>% group_by(variable) %>% tally()

JoinCorrelations %>% select_at(vars(contains("ComBat_q.value"))) %>% setNames(str_remove_all(names(.), "_q.value")) %>% reshape2::melt(.) %>% filter(value < 0.05) %>% group_by(variable) %>% tally()

## What are the systems that replicate with q values < 0.05 ??
JoinCorrelations %>% select_at(vars(c("GO", "TERM",  contains("ComBat_q.value")))) %>% setNames(str_remove_all(names(.), "_q.value")) %>% reshape2::melt(.) %>% filter(value < 0.05) %>% group_by(GO, TERM) %>% tally() %>% filter(n>1) %>% view(.)

#------------------------------------------------------
## Figure 1. Barplot of significant correlations
#------------------------------------------------------
## Get null distribution
null_dist = replicate(1000, sum(runif(nrow(JoinCorrelations)) < 0.025))
lwr05 = floor(quantile(null_dist, 0.05))
upp95 = ceiling(quantile(null_dist, 0.95))
xstart = 0
xend = length(ComBat_datasets) + 1
rdm_annotation = data.frame(
  CorrelationSign = "Negative",
  x = 2,
  y = round(mean(null_dist)) + 50,
  label = paste0("Expected at random (p < 0.05): ", round(mean(null_dist)))
)

### Count significant correlation by p values q values and correlation sign
pvals = JoinCorrelations %>% select_at(vars(contains("ComBat_p.value"))) %>% setNames(str_remove_all(names(.), "_ComBat_p.value")) %>% reshape2::melt(., value.name = "pvalue")
qvals = JoinCorrelations %>% select_at(vars(contains("ComBat_q.value"))) %>% setNames(str_remove_all(names(.), "_ComBat_q.value")) %>% reshape2::melt(., value.name = "qvalue")
COUNTS = JoinCorrelations %>% select_at(vars(contains("ComBat_estimate"))) %>%
  setNames(str_remove_all(names(.), "_ComBat_estimate")) %>%
  reshape2::melt(., value.name = "estimate") %>%
  mutate(CorrelationSign = factor(ifelse(estimate > 0, "Positive", "Negative"))) %>%
  bind_cols(., pvals %>% dplyr::select(pvalue)) %>%
  bind_cols(., qvals %>% dplyr::select(qvalue)) %>%
  mutate(p_significant = pvalue < 0.05,
         q_significant = qvalue < 0.05) %>%
  group_by(variable, CorrelationSign) %>% summarise(p_count = sum(p_significant),
                                                    q_count = sum(q_significant)) %>%
  ungroup()

## Define color pal
COLORS = RColorBrewer::brewer.pal(6, "Paired")[c(1, 2, 5, 6)]

## Bar plot split by corr sign
g_bar_positive = ggplot(COUNTS) +
  ## Main
  geom_bar(data = filter(COUNTS, CorrelationSign == "Positive"), 
           aes(x = variable, y = p_count, fill = "p-value < 0.05"), 
           stat = "identity") +
  ## Uncomment these to include q values in figure
  # geom_bar(
  #   data = filter(COUNTS, CorrelationSign == "Positive"),
  #   aes(x = variable, y = q_count, fill = "q-value < 0.05"),
  #   stat = "identity",
  #   alpha = 0.7
  # ) +
  ## Colors
  ## Change "values = COLORS[2]" to "values = COLORS[c(1,2)]" include q values in figure
  scale_fill_manual("Significance threshold", values = COLORS[2]) +
  ## Null distribution
  geom_rect(
    aes(xmin = xstart, xmax = xend , ymin = lwr05, ymax = upp95),
    fill = "gray60",
    alpha = 0.1 ) +
  geom_hline(yintercept = round(mean(null_dist)),
             col = "gray30",
             lty = 2) +
  geom_text(data = rdm_annotation,
    aes(x = x, y = y, label = label),
    col = "gray30",
    size = 4 ) +
  ## Theme
  ylim(c(0, 2200)) +
  theme_classic()  +
  theme(
    title = element_text(size = 15, hjust = 0),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.title = element_text(
      size = 15,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  ## Labels
  labs(x = "Dataset", y = "Number of significant systems", title = "Positive associations")   +
  geom_text(
    data = filter(COUNTS, CorrelationSign == "Positive"),
    aes(x = variable, y = p_count + 50, label = p_count),
    size = 4
  ) #+
  ## Uncomment following and the "+" line above to add qvalues
  # geom_text(
  #   data = filter(COUNTS, CorrelationSign == "Positive"),
  #   aes(x = variable, y = q_count, label = q_count),
  #   size = 4,
  ## Switch "COLORS[2]" to "COLORS[c(1,2)]" to add q values
  #   col = COLORS[2]
  # )


g_bar_negative = ggplot(NULL) +
  ## Main
  geom_bar(
    data = filter(COUNTS, CorrelationSign == "Negative"),
    aes(x = variable, y = p_count, fill = "p-value < 0.05"),
    stat = "identity"
  ) +
  ## Uncomment these to include q values in figure
  # geom_bar(
  #   data = filter(COUNTS, CorrelationSign == "Negative"),
  #   aes(x = variable, y = q_count, fill = "q-value < 0.05"),
  #   stat = "identity",
  #   alpha = 75
  # ) +
  
  ## Colors
  ## Change "values = COLORS[4]" to "values = COLORS[c(3,4)]" include q values in figure
  scale_fill_manual("Significance threshold", values = COLORS[4]) +
  ## Null distribution
  geom_rect(
    data = COUNTS,
    aes(
      xmin = xstart,
      xmax = xend ,
      ymin = lwr05,
      ymax = upp95
    ),
    fill = "gray60",
    alpha = 0.1
  ) +
  geom_hline(yintercept = round(mean(null_dist)),
             col = "gray30",
             lty = 2) +
  geom_text(
    data = rdm_annotation,
    aes(x = x, y = y, label = label),
    col = "gray30",
    size = 4
  ) +
  ## Theme
  ylim(c(0, 2200)) +
  theme_classic()  +
  theme(
    title = element_text(size = 15, hjust = 0),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.title = element_text(
      size = 15,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "none"
  ) +
  ## Labels
  labs(x = "Dataset", y = "Number of significant systems", title = "Negative associations")   +
  geom_text(
    data = filter(COUNTS, CorrelationSign == "Negative"),
    aes(x = variable, y = p_count + 50, label = p_count),
    size = 4
  ) #+
  ## Uncomment following and the "+" line above to add qvalues
  # geom_text(
  #   data = filter(COUNTS, CorrelationSign == "Negative"),
  #   aes(x = variable, y = q_count + 20, label = q_count),
  #   size = 4,
  ## Switch "COLORS[2]" to "COLORS[c(1,2)]" to add q values
  #   col = COLORS[4]
  # )

## plot
fig1 = ggarrange(
  g_bar_positive,
  g_bar_negative,
  labels = LETTERS[1:2],
  font.label = list(size = 20)
)

fig1

## Export to pdf and png
png_savename = paste0(figure_save_directory, "number_significant_correlations_barplot.png")
pdf_savename = paste0(figure_save_directory, "number_significant_correlations_barplot.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename, width = fig_width, height = fig_height)
fig1
dev.off()

png(file = png_savename, width = fig_width, height = fig_height, units = "in", res = 300)
fig1
dev.off()


#--------------------------------
### SuperExactTest : p.value and correlation sign must be equivalent
#--------------------------------
library(SuperExactTest)

## Prep data into list for the test
SuperExactData = lapply(ComBat_Cor_objects_names[length(ComBat_Cor_objects_names):1], function(COR_OBJ) {
  get(COR_OBJ) %>% filter(p.value < 0.05) %>%
    mutate(uniq_GO_id = paste0(GO, ifelse(.$estimate > 0, "pos", "neg"))) %>%
    dplyr::select(uniq_GO_id)
})
names(SuperExactData) = c(gsub(pattern = "_ComBat", "", ComBat_datasets[length(ComBat_Cor_objects_names):1]))
SuperExactData = lapply(SuperExactData, function(xx)
  return(xx$uniq_GO_id))

## Super exact test
SuperTestRes = supertest(SuperExactData, n = nrow(JoinCorrelations))

#------------------------------------------------------
## Super test results table
#------------------------------------------------------
## Get summary table
SuperTestTable = summary(SuperTestRes)
# Filter for intersections only and adjust p values with BH
SuperTestTable = as_tibble(SuperTestTable$Table) %>% filter(Degree > 1) %>%  mutate(Q.value = p.adjust(P.value, method = "fdr"))
## write to file
fwrite(
  SuperTestTable,
  file = paste0(figure_save_directory,"../tables/SuperTest_table.txt"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

#------------------------------------------------------
## Figure 2. Barplot of overlap
#------------------------------------------------------
p_threshold_for_figure = SuperTestTable %>% filter(Q.value < 0.05) %>% summarise(max(P.value)) 

plot(
  SuperTestRes,
  Layout = "landscape",
  sort.by = "degree",
  degree = 2:length(ComBat_datasets),
  show.set.size = T,
  minMinusLog10PValue = (-log10(p_threshold_for_figure$`max(P.value)`)),
  keep.empty.intersections = T,
  margin = c(.5, 9, 1.5, 2),
  show.expected.overlap = T,
  expected.overlap.style = "horizBar",
  color.expected.overlap = "blue",
  ylab = "Number of overlapping systems"
)

## Export to pdf and png
png_savename = paste0(figure_save_directory, "overlap_supertest_barplot.png")
pdf_savename = paste0(figure_save_directory, "overlap_supertest_barplot.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename, width = fig_width, height = fig_height)
plot(
  SuperTestRes,
  Layout = "landscape",
  sort.by = "degree",
  degree = 2:length(ComBat_datasets),
  show.set.size = T,
  minMinusLog10PValue = (-log10(p_threshold_for_figure$`max(P.value)`)),
  keep.empty.intersections = T,
  margin = c(.5, 9, 1.5, 2),
  show.expected.overlap = T,
  expected.overlap.style = "horizBar",
  color.expected.overlap = "blue",
  ylab = "Number of overlapping systems"
)
dev.off()

png(file = png_savename, width = fig_width, height = fig_height, units = "in", res = 300)
plot(
  SuperTestRes,
  Layout = "landscape",
  sort.by = "degree",
  degree = 2:length(ComBat_datasets),
  show.set.size = T,
  minMinusLog10PValue = (-log10(p_threshold_for_figure$`max(P.value)`)),
  keep.empty.intersections = T,
  margin = c(.5, 9, 1.5, 2),
  show.expected.overlap = T,
  expected.overlap.style = "horizBar",
  color.expected.overlap = "blue",
  ylab = "Number of overlapping systems"
)
dev.off()

#------------------------------------------------------
## TABLE 2 : replicated in all 5 datasets
#------------------------------------------------------
## What are the five terms that are replicated between datasets? 
replicated_in_all = SuperTestTable %>% dplyr::filter(Degree == 5) %>% dplyr::select(Elements) %>% str_remove_all(string = .,pattern = "neg| ") %>% str_split(.,",")
JoinCorrelations %>% dplyr::select(GO,TERM) %>% dplyr::filter(GO %in% replicated_in_all[[1]])


## Group size of replicated systems
system_list_idx = match(replicated_in_all[[1]], systems_names)
lengths(systems)[system_list_idx]

## Number of unique genes
lengths(systems_symbols)[system_list_idx]

#------------------------------------------------------
## pvalue of the overlap ?
#------------------------------------------------------
SuperTestTable %>% dplyr::filter(Degree == 5) %>% dplyr::select(P.value)


#------------------------------------------------------
### Get overlapping results
#------------------------------------------------------
Significant_overlap = as_tibble(summary(SuperTestRes)$Table) %>%
  mutate(Q.value = p.adjust(P.value, method = "BH")) %>%
  filter(Q.value < 0.05 & Degree > 2)

## Get list of all the system that overlap
sys_with_overlap = Significant_overlap$Elements %>% str_split(., ",") %>% unlist(.) %>% str_remove_all(., " ") %>% unique(.) %>% str_remove_all(., "neg|pos")
length(sys_with_overlap)

#------------------------------------------------------
## Supp File
#------------------------------------------------------
## Save correlation results replicated to file
JoinCorrelations %>% filter(GO %in% sys_with_overlap) %>% 
  fwrite(.,
  file = paste0(figure_save_directory, "../tables/Systems_overlapping_3_sets.txt"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)


#------------------------------------------------------
## Overlapping systems: analyses
#------------------------------------------------------
## Sign of correlations
sys_sign = Significant_overlap$Elements %>% str_split(., ",") %>% unlist(.) %>% str_remove_all(., " ") %>% unique(.) %>% str_extract_all(., "neg|pos") %>% table(unlist(.))

## Names of the systems
GOTERM_overlap = AnnotationDbi::select(GO.db::GO.db, keys = sys_with_overlap, columns = "TERM")$TERM
names(GOTERM_overlap) = sys_with_overlap

### To use in external tools (such as revigo)
cat(names(GOTERM_overlap), sep = ",")

#------------------------------------------------------
#### Check overlap with Brinkmeyer-Langford et al.
### Check if an ancestor of these terms is in the list 216 replicated GO terms
require(ontologyIndex)
require(GO.db)

Brinkmeyer_etal = c(
  "GO:0001309",
  "GO:0032214",
  "GO:0003691",
  "GO:0044062",
  "GO:0009756",
  "GO:0003840",
  "GO:0036374",
  "GO:2000146"
)

# read obo file for ontology structure
go_obo = get_ontology(
  "25032020_go-basic.obo",
  extract_tags = "everything",
  propagate_relationships = c("is_a", "is_part", "regulates")
)

## get ancestors
Brinkmeyer_ancestors = lapply(Brinkmeyer_etal, function(xx)
  get_ancestors(go_obo, terms = xx))

## Check overlap with our study
brinkmeyer_overlap_idx = which(unlist(lapply(Brinkmeyer_ancestors, function(xx)
  any(xx %in% sys_with_overlap))))

dufour_al_sys_with_overlap = GOTERM_overlap[sys_with_overlap %in% Brinkmeyer_ancestors[[brinkmeyer_overlap_idx]]]
select(GO.db, keys = Brinkmeyer_etal[brinkmeyer_overlap_idx], columns = "TERM")

#------------------------------------------------------
## Supp Figure 2  
#------------------------------------------------------
png_savename = paste0(figure_save_directory, "Spearman_of_corr_corplot.png")
pdf_savename = paste0(figure_save_directory, "Spearman_of_corr_corplot.pdf")
# lab_x_pos = -0.5
# lab_y_pos = 6.5

pdf(pdf_savename)
#par(mfrow = c(1, 2))
corrplot_data = JoinCorrelations %>% select_at(vars(contains("ComBat_estimate"))) %>%
  setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)

COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
corrplot::corrplot(
  COR$r,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  p.mat = COR$P,
  #title = "All systems tested"#,
  #mar = c(0, 0, 0, 0) 
)
dev.off()
# 
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[1],
#   cex = 1.5
# )

### Overlapping
# corrplot_data = JoinCorrelations %>% filter(GO %in% sys_with_overlap) %>% select_at(vars(contains("ComBat_estimate"))) %>%
#   setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)
# 
# COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
# 
# corrplot::corrplot(
#   COR$r,
#   method = "color",
#   type = "upper",
#   addCoef.col = "black",
#   p.mat = COR$P,
#   #main = "Systems replicated in 3 datasets"#,
#   #mar = c(0, 0, 0, 0) 
# )
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[2],
#   cex = 1.5
# )

# ### NON Overlapping
# corrplot_data = JoinCorrelations %>% filter(!GO %in% sys_with_overlap) %>% select_at(vars(contains("ComBat_estimate"))) %>%
#   setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)
# 
# COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
# 
# corrplot::corrplot(
#   COR$r,
#   method = "color",
#   type = "upper",
#   addCoef.col = "black",
#   p.mat = COR$P,
#   main = "Systems not replicated in 3 datasets",
#   mar = c(1, 0, 2, 0) + 0.5
# )
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[3],
#   cex = 1.5
# )
#dev.off()

###
png(png_savename)
# par(mfrow = c(1, 2))
corrplot_data = JoinCorrelations %>% select_at(vars(contains("ComBat_estimate"))) %>%
setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)

COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
corrplot::corrplot(
COR$r,
method = "color",
type = "upper",
addCoef.col = "black",
p.mat = COR$P,
#title = "All systems tested"#,
#mar = c(0, 0, 0, 0) 
)
dev.off()
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[1],
#   cex = 1.5
# )

### Overlapping
# corrplot_data = JoinCorrelations %>% filter(GO %in% sys_with_overlap) %>% select_at(vars(contains("ComBat_estimate"))) %>%
#   setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)
# 
# COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
# 
# corrplot::corrplot(
#   COR$r,
#   method = "color",
#   type = "upper",
#   addCoef.col = "black",
#   p.mat = COR$P,
#   #main = "Systems replicated in 3 datasets"#,
#   #mar = c(0, 0, 0, 0) 
# )
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[2],
#   cex = 1.5
# )

# ### NON Overlapping
# corrplot_data = JoinCorrelations %>% filter(!GO %in% sys_with_overlap) %>% select_at(vars(contains("ComBat_estimate"))) %>%
#   setNames(str_replace_all(names(.), "_ComBat_estimate", "")) %>% as.matrix(.)
# 
# COR = Hmisc::rcorr(na.omit(corrplot_data), type = "spearman")
# 
# corrplot::corrplot(
#   COR$r,
#   method = "color",
#   type = "upper",
#   addCoef.col = "black",
#   p.mat = COR$P,
#   main = "Systems not replicated in 3 datasets",
#   mar = c(1, 0, 2, 0) + 0.5
# )
# text(
#   x = lab_x_pos ,
#   y = lab_y_pos,
#   labels = LETTERS[3],
#   cex = 1.5
# )
#dev.off()


#------------------------------------------------------
### Figure 3: GOslim map and enrichment
#------------------------------------------------------
#####  Get go_slims for grouping values
require(ontologyIndex)
require(GO.db)
# read obo file for ontology structure
#"is_part", "regulates"
  
go_obo = get_ontology(
  "/media/fdufour/Seagate Expansion Drive/Data/Systems/GO/25032020_go-basic.obo",
  extract_tags = "everything", propagate_relationships = c("is_a","is_part", "regulates")
)
goslim_obo = get_ontology(
  "/media/fdufour/Seagate Expansion Drive/Data/Systems/GO/25032020_goslim_generic.obo",
  extract_tags = "minimal"
)
goslim_ids = na.omit(str_extract(goslim_obo$id, "GO:\\d+"))

go_id = "GO:0007045" ## for troubleshooting purposes


get_goslim = function(GO_id, go_obo = go_obo, goslim_ids = goslim_ids) {
  #GO_id = go_id ## for troubleshooting purposes
  list_idx = grep(GO_id, go_obo$id)
  
  if (length(list_idx) == 0){
    return(data.frame(
      GO_id = GO_id,
      GO_slim = NA))
  }
  
  GO_ancestors = go_obo$ancestors[[list_idx]] 
  GO_slims = GO_ancestors[which(GO_ancestors %in% unlist(goslim_ids))]
  
  if (length(GO_slims) == 0){
    return(data.frame(
      GO_id = GO_id,
      GO_slim = NA))
  }
  
  return(data.frame(
    GO_id = GO_id,
    GO_slim = GO_slims))
}

## Test
get_goslim(GO_id = go_id,go_obo = go_obo, goslim_ids = goslim_ids)

## Map replicated systems to GO slims
goslimmap_list = lapply(X = sys_with_overlap,
                        FUN = get_goslim,
                        go_obo = go_obo,
                        goslim_ids = goslim_ids)
goslimmap = do.call(what = "rbind", goslimmap_list)

## Save the mappings to go slim 
as.data.frame(goslimmap) %>%  
  na.omit(.) %>% 
  mutate(GOTerm = AnnotationDbi::select(GO.db::GO.db, keys = as.character(GO_id),columns = "TERM")$TERM,
        GOslimTerm = AnnotationDbi::select(GO.db::GO.db, keys = as.character(GO_slim),columns = "TERM")$TERM) %>%    
  arrange(GOslimTerm) %>% 
  fwrite(.,
         file = paste0(figure_save_directory, "../tables/GOslim_map_non_summerised.txt"),
         row.names = F, 
         col.names = T)

## Map all systems tested to GO slims to calculate expected proportions
goslimmap_list_expected = lapply(X = JoinCorrelations$GO,
  FUN = get_goslim,
  go_obo = go_obo,
  goslim_ids = goslim_ids
)
goslimmap_expected = do.call(what = "rbind", goslimmap_list_expected)


### Prepare object for plotting
Expected = na.omit(goslimmap_expected)  %>%  group_by(GO_slim) %>% tally() %>%
  dplyr::rename(Expected_Count = n) %>%
  mutate(
    TERM = AnnotationDbi::select(GO.db::GO.db, keys = as.character(.$GO_slim) , columns = "TERM")$TERM,
    TOTAL = case_when(TERM == "biological_process" ~ Expected_Count)
  )  %>%
  fill(TOTAL, .direction = "downup") %>%
  na.omit() %>%
  filter(TERM != "biological_process") %>%
  mutate(Expected_Proportion = Expected_Count / TOTAL) %>% dplyr::select(-TOTAL)

Observed = na.omit(goslimmap)  %>%  group_by(GO_slim) %>% tally() %>%
  dplyr::rename(Observed_Count = n) %>%
  mutate(
    TERM = AnnotationDbi::select(GO.db::GO.db, keys = as.character(.$GO_slim) , columns = "TERM")$TERM,
    TOTAL = case_when(TERM == "biological_process" ~ Observed_Count)
  )  %>%
  fill(TOTAL, .direction = "downup") %>%
  na.omit() %>%
  filter(TERM != "biological_process") %>%
  mutate(Observed_Proportion = Observed_Count / TOTAL)

## Test enrichment
z_test =  Observed %>%
  left_join(., Expected, by = c("GO_slim", "TERM")) %>%
  mutate(result = pmap(., ~ prop.test(
    x = ..2,
    n = ..4,
    p = ..7,
    correct = F
  ))) %>%
  mutate(result = map(result, broom::tidy)) %>%
  unnest(result) %>%
  mutate(qvalue = p.adjust(.$p.value, method = "fdr")) %>%
  dplyr::select(c(1, 8:16))


## Combine and save to a file
Observed %>% left_join(., Expected, by = c("GO_slim", "TERM"))  %>% 
  left_join(., z_test, by = "GO_slim") %>%
  mutate(Expected_Count_2 = floor(Expected_Proportion * TOTAL),
         p.value_2 = case_when(p.value < 0.05 ~ p.value)) %>% 
  fwrite(.,
         file = paste0(figure_save_directory, "../tables/GOslim_map.txt"),
         row.names = F,
         col.names = T)



## Color based on p value, If pval > 0.05 set to white
goslim_map_pval = Observed %>% 
  left_join(., Expected, by = c("GO_slim", "TERM"))  %>% 
  left_join(., z_test, by = "GO_slim") %>%
  mutate(Expected_Count_2 = floor(Expected_Proportion * TOTAL),
         p.value = case_when(p.value < 0.05 ~ p.value)) %>%
  ggplot(., aes(Observed_Count, fct_reorder(TERM, Observed_Count))) + 
  geom_bar(aes(fill = -log10(p.value)), stat = "identity", col = "gray50") + 
  scale_x_continuous(n.breaks = 20) + 
  labs(y ="GO slim term", x = "Count") + 
  geom_errorbar(aes(x = Expected_Count_2, xmin = Expected_Count_2, xmax = Expected_Count_2, col = "Expected")) +
  scale_colour_manual(
    name = '',
    values = c("Expected" = "blue"),
    guide = "none"
  ) +
  theme_bw() + 
  theme(text = element_text(size = 15), legend.position = "right") + 
  scale_fill_gradientn(expression("-Log"[10]*"(P)"), colours = heat.colors(100, rev = T), na.value = "white")


## Color based on q value (BH), If qval > 0.05 set to white
goslim_map_qval = Observed %>% 
  left_join(., Expected, by = c("GO_slim", "TERM"))  %>% 
  left_join(., z_test, by = "GO_slim") %>%
  mutate(Expected_Count_2 = floor(Expected_Proportion * TOTAL),
         p.value = case_when(qvalue < 0.05 ~ p.value)) %>%
  ggplot(., aes(Observed_Count, fct_reorder(TERM, Observed_Count))) + 
  geom_bar(aes(fill = -log10(p.value)), stat = "identity", col = "gray50") + 
  scale_x_continuous(n.breaks = 20) + 
  labs(y ="GO slim term", x = "Count") + 
  geom_errorbar(aes(x = Expected_Count_2, xmin = Expected_Count_2, xmax = Expected_Count_2, col = "Expected")) +
  scale_colour_manual(
    name = '',
    values = c("Expected" = "blue"),
    guide = "none"
  ) +
  theme_bw() + 
  theme(text = element_text(size = 15), legend.position = "right") + 
  scale_fill_gradientn(expression("-Log"[10]*"(P)"), colours = heat.colors(100, rev = T), na.value = "white")


## Export plots
png_savename_pval = paste0(figure_save_directory, "GO_slim_enrichment_FILTER_pvalues.png")
pdf_savename_pval = paste0(figure_save_directory, "GO_slim_enrichment_FILTER_pvalues.pdf")
png_savename_qval = paste0(figure_save_directory, "GO_slim_enrichment_FILTER_qvalues.png")
pdf_savename_qval = paste0(figure_save_directory, "GO_slim_enrichment_FILTER_qvalues.pdf")
fig_width = 10
fig_height = 8

pdf(file = pdf_savename_pval,
    width = fig_width,
    height = fig_height)
goslim_map_pval
dev.off()

png(file = png_savename_pval,
    width = fig_width,
    height = fig_height,
    units = "in",
    res = 300)
goslim_map_pval
dev.off()

pdf(file = pdf_savename_qval,
    width = fig_width,
    height = fig_height)
goslim_map_qval
dev.off()

png(file = png_savename_qval,
    width = fig_width,
    height = fig_height,
    units = "in",
    res = 300)
goslim_map_qval
dev.off()

### -------------------------------------------------------
## Fig supp: System size distribution
## Read systems
systems_Filename = "GOterms2probeListFromComBat.txt"
systems = scan(file = systems_Filename, what = "character", sep = "\n")
systems = lapply(systems, function(xx)
  return(unlist(strsplit(xx, split = "\t"))))
systems_names = lapply(systems, function(xx)
  xx[1])
systems = lapply(systems, function(xx)
  return(xx[-1]))


## Get probe symbols
map2symbol = as.list(illuminaHumanv4.db::illuminaHumanv4SYMBOL)
systems_symbols = lapply(systems, function(xx) return(unique(unlist(map2symbol[match(xx,names(map2symbol))]))))

## PLot group size distributions
group_size_scatter = ggplot(data.frame(Sys_size = unlist(lengths(systems)), uniq_symbols = unlist(lengths(systems_symbols))), aes(Sys_size, uniq_symbols)) + 
  geom_point() + 
  labs(x = "Group size", y = "Number of unique genes") + 
  theme_bw()+
  theme(text = element_text(size = 15))

group_size_distributions = ggExtra::ggMarginal(group_size_scatter,type = "histogram")
  

## Export figure to pdf and png
png_savename = paste0(figure_save_directory, "distributions_group_size.png")
pdf_savename = paste0(figure_save_directory, "distributions_group_size.pdf")
fig_width = 10
fig_height = 8


pdf(file = pdf_savename,
    width = fig_width,
    height = fig_height)
group_size_distributions
dev.off()

png(
  file = png_savename,
  width = fig_width,
  height = fig_height,
  units = "in",
  res = 300
)
group_size_distributions
dev.off()

