library(diffcyt)
library(ggplot2)
library(dplyr)
library(HDF5Array)


tissue_counts <- loadHDF5SummarizedExperiment(tissue_dir)
blood_counts <- loadHDF5SummarizedExperiment(blood_dir)

design <- read.csv('design_tissue.csv') # blood and tissue have the same design matrix

contrasts = makeContrasts(
  DrugvsControl04h = drug.4h- control.4h,
  DrugvsControl24h = drug.24h- control.24h,
  DrugvsControl48h = drug.48h -control.48h,
  DrugvsControl72h = drug.72h - control.72h,
  DrugvsNaproxen04h = drug.4h - naproxen.4h,
  DrugvsNaproxen24h = drug.24h - naproxen.24h,
  DrugvsNaproxen48h = drug.48h - naproxen.48h,
  DrugvsNaproxen72h = drug.72h - naproxen.72h,
  levels = design)

edgeR_function <- function(counts){
  contrasts_list <- c()
  results <- data.frame()
  for (i in 1:8){
    temp_contrast <- contrasts[,i]
    contrast_name <- colnames(contrasts)[i]
    temp_contrast_val <- as.numeric(temp_contrast)
    
    res <- testDA_edgeR(counts, design, temp_contrast_val)
    top_res <- data.frame(diffcyt::topTable(res, all = T, format_vals = TRUE))
    contrast_names_rep <- rep(contrast_name, nrow(top_res))
    contrasts_list <- c(contrasts_list, contrast_names_rep) # 31 clusters
    
    results <- rbind(results, top_res)
  }
  results$contrast <- contrasts_list
  
  return(results)
}


tissue_results <- edgeR_function(tissue_counts)
blood_results <- edgeR_function(blood_counts)

write.csv(tissue_results, 'tissue_results_DA.csv', row.names = F)
write.csv(blood_results, 'blood_results_DA.csv', row.names = F)

