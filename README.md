# Analysis scripts for Dufour **et al**. Physiological System Dysregulation in Gene Expression Correlates Negatively with Age

The manuscript was submitted to [eLife's special issue on aging, geroscience and longevity.](https://elifesciences.org/inside-elife/4f706531/special-issue-call-for-papers-in-aging-geroscience-and-longevity?utm_source=civimail&utm_medium=email&utm_campaign=aging_special_issue_announcement) 

The "script" directory includes scripts for:
1. [sample quality control](./Step1_QC_Clustering.R) 
2. [data processing](./Step2_ProcessDataComBat.R)
3. [create groups from gene ontology](./Step3_Make_group_list_from_GO.R)
4. [DM estimation and correlation analyses](./Step4_DM_Estimation_Correlation2Age.R)
5. [results comparaison and figures](./Step5_Analyse_correlations_generate_figures.R)


## To reuse the scripts: 

### List of packages to install: 
1. "tidyverse"
2. "data.table"
3. "ggpubr"
4. "beadarray"
5. "broom"
6. "sva"
7. "illuminaHumanv4.db"
8. "limma"
9. "mgsa"
10. "ontologyIndex"
11. "RBGL"
12. "graph"
13. "GO.db"



### Sample quality control
1. Change the file names to match your own.
2. Variable names included with the metadata of datasets were changed by F.Dufour. Change your variable names (or the variable names in the script) to match your own data. 

### Data processing.
You can run the script from bash, using:  
```
Rscript dir\_to\_script/Step2_ProcessDataComBat.R dir\_to\_param_file/parameters.txt working_directory
```
The parameter file is a CSV file \(with a header\) that contains the following: 
1. col1 = dir\_to\_data/data.txt 
2. col2 = dir\_to\_metadata/metadata.txt

The working directory includes the ending "/". Example: "~/user/documents/"

### create groups from gene ontology
You need to download the [ontology](http://geneontology.org/docs/download-ontology/) file \(.obo\) and the [go association file](http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/) from the Gene Ontology Consortium website and change the directories and file names in the script.

You will also need the data files from the datasets in the analyses to set the minimum group size and create groups based on the common probes among all datasets after preprocessing. 

### DM estimation and correlation analyses.
You can run the script from bash, using: 
```
Rscript dir\_to\_script/Step3_DM_Estimation_Correlation2Age.R dir\_to\_param_file/parameters.txt dir\_to\_systems_file/systems.txt
```

The parameter file is a CSV file \(with a header\) that contains the following:
1. col1 = dir\_to\_data/data.txt
2. col2 = dir\_to\_metadata/metadata.txt


The System file is a tab separated file that is read line-by-line and the first entry in the line is the GO term or the system's name. The rest of the line lists the ILMN id in the system.  
Ex:  
GO:0000001  ILMN1010101010  ILMN267262626 ILMN6263537   
GO:000002 ILMN2263636 ILMN378383874 ILMN237843  ILMN27237637  ILMN437387   

### Results comparaison and figures
1. Change the file names to match your own.

## Raw data for datasets included in the paper can be downloaded from: 
- Leipzig LIFE Heart Study (Leipzig) [Accession: GSE65907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65907)
- AddNeuroMed Consortium (AddNeuro) [Accession: GSE63063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63063)
- Grady Trauma Project (GTP) [Accession: GSE58137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58137)
- Emory University Center for Health Discovery and Well-being (CHDWB) [Accession: GSE61672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61672)
- UK Primary Sjögren’s Syndrome Registry (UKPSSR)  [Accession: E-MTAB-8277](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8277/)


## AddNeuro dataset: 
This dataset includes samples measured on both the illumina human-HT12-v3 and illumina human-HT12-v4 platforms. 
An additional step was added to process this data. The data was adjusted for platform version using the ComBat algorithm implemented in the sva package see [vignette](https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf). 

The data preprocessing script specific to AddNeuro can be access [here](./Step2_ProcessDataComBat_AddNeuro.R). 
To use the script, change the working directory in the file.

