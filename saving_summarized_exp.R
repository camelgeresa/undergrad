library(readxl)
library(flowCore)
library(CATALYST)
library(diffcyt)
library(ggplot2)
library(dplyr)
library(HDF5Array)
library(stringr)

# Diffcyt edgeR only works with summarized experiment.
# Need to convert singlecellexperiment to Summarized Experiment. 

# LOAD CLUSTERED DATASET
load("Workspace.RData")
merging_table2 <- "./Excel_files/TE94_cluster_merging_2.xlsx"
merging_table2 <- read_excel(merging_table2)
merging_table2
head(data.frame(merging_table2))
sce <- mergeClusters(sce, k = "meta66", table = merging_table2, id = "Cluster2", overwrite = TRUE)

merge_table3 <- './Excel_files/map3.xlsx'
merge_table3 <- read_excel(merge_table3)
head(data.frame(merge_table3))
sce <- mergeClusters(sce, k = "Cluster2", table = merge_table3, id = "Cluster3", overwrite = T)

set.seed(2024)


clustering_to_use <- 'Cluster2'

# Summarized cell counts only for samples found in the blood:
blood_sce <- filterSCE(sce, grepl('Blood', colData(sce)$condition),!grepl('00h',colData(sce)$condition))
code_id <- colData(blood_sce)$cluster_id
cluster_id <- metadata(blood_sce)$cluster_codes[, clustering_to_use][code_id]
stopifnot(length(cluster_id) == nrow(colData(blood_sce)), 
          length(code_id) == nrow(colData(blood_sce)))
colData(blood_sce)$code_id <- code_id
colData(blood_sce)$cluster_id <- cluster_id
# clustering identifier to store in metadata
clustering_name <- clustering_to_use
cs_by_s <- split(seq_len(ncol(blood_sce)), colData(blood_sce)$sample_id)

# Design matrix
experiment_info <- metadata(blood_sce)$experiment_info
ei <- experiment_info %>%
  mutate(
    drug_type  = as.factor(
      case_when(
    grepl('drug', condition)~'drug',
    grepl('naproxen',condition)~'naproxen',
    grepl('control',condition)~'control')),
    
    time = str_extract(sample_id,"\\b\\d+h\\b"), #as.numeric(str_extract(sample_id,"(?<=\\s)\\d+(?=h)")),
    condition = paste(drug_type, time, sep='.'))
  
design <-  model.matrix(~0+ei$condition) # slight chnage from source to remove intercept, since we do not have a reference.
colnames(design) <- sub("ei$condition", "", colnames(design), fixed=T)
write.csv(design, 'design_blood.csv', row.names = F)

# re-order according to experiment_info
cs <- unlist(cs_by_s[as.character(experiment_info$sample_id)])
es <- t(assays(blood_sce)[["exprs"]])[cs, , drop = FALSE]
# create SummarizedExperiment (in transposed format compared to SingleCellExperiment)
d_se <- SummarizedExperiment(
  assays = list(exprs = es), 
  rowData = colData(blood_sce)[cs, ], 
  colData = rowData(blood_sce), 
  metadata = metadata(blood_sce)
)

d_counts <- calcCounts(d_se) 
#dir <- tempfile("h5_blood_")
blood_dir <- 'blood_sum_exp'
h5_se0 <- saveHDF5SummarizedExperiment(d_counts, blood_dir) # we can now reuse this for diffcyt
list.files(dir)

# For tissue: 

tissue_sce <- filterSCE(sce, grepl('Peritoneum', colData(sce)$condition),!grepl('00h',colData(sce)$condition))
code_id <- colData(tissue_sce)$cluster_id
cluster_id <- metadata(tissue_sce)$cluster_codes[, clustering_to_use][code_id]
stopifnot(length(cluster_id) == nrow(colData(tissue_sce)), 
          length(code_id) == nrow(colData(tissue_sce)))
colData(tissue_sce)$code_id <- code_id
colData(tissue_sce)$cluster_id <- cluster_id
# clustering identifier to store in metadata
clustering_name <- clustering_to_use
cs_by_s <- split(seq_len(ncol(tissue_sce)), colData(tissue_sce)$sample_id)

experiment_info <- metadata(tissue_sce)$experiment_info
ei <- experiment_info %>%
  mutate(
    drug_type  = as.factor(
      case_when(
        grepl('drug', condition)~'drug',
        grepl('naproxen',condition)~'naproxen',
        grepl('control',condition)~'control')),
    
    time = str_extract(sample_id,"\\b\\d+h\\b"), #as.numeric(str_extract(sample_id,"(?<=\\s)\\d+(?=h)")),
    condition = paste(drug_type, time, sep='.'))

design <-  model.matrix(~0+ei$condition) # slight chnage from source to remove intercept, since we do not have a reference.
colnames(design) <- sub("ei$condition", "", colnames(design), fixed=T)
write.csv(design, 'design_tissue.csv', row.names = F)

# re-order according to experiment_info
cs <- unlist(cs_by_s[as.character(experiment_info$sample_id)])
es <- t(assays(tissue_sce)[["exprs"]])[cs, , drop = FALSE]
# create SummarizedExperiment (in transposed format compared to SingleCellExperiment)
d_se <- SummarizedExperiment(
  assays = list(exprs = es), 
  rowData = colData(tissue_sce)[cs, ], 
  colData = rowData(tissue_sce), 
  metadata = metadata(tissue_sce)
)

d_counts <- calcCounts(d_se) 
tissue_dir <- 'tissue_sum_exp'
#dir <- tempfile("h5_tissue_")
h5_se0 <- saveHDF5SummarizedExperiment(d_counts, tissue_dir) # we can now reuse this for diffcyt
list.files(dir)

