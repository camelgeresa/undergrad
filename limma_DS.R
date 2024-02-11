library(readxl)
library(flowCore)
library(CATALYST)
library(diffcyt)
library(ggplot2)
library(stringr)
library(limma)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(cowplot)
library(ggpubr)

# LOAD CLUSTERED DATASET
# load("Workspace.RData")
# merging_table2 <- "./Excel_files/TE94_cluster_merging_2.xlsx"
# merging_table2 <- read_excel(merging_table2)
# merging_table2
# head(data.frame(merging_table2))
# sce <- mergeClusters(sce, k = "meta66", table = merging_table2, id = "Cluster", overwrite = TRUE)

# merge_table3 <- './Excel_files/map3.xlsx'
# merge_table3 <- read_excel(merge_table3)
# head(data.frame(merge_table3))
# sce <- mergeClusters(sce, k = "Cluster2", table = merge_table3, id = "Cluster3", overwrite = T)

#p_cutoff <- 0.01
#markers <- rownames(sce)
#markers <- gsub("[/-]", ".", markers)

#exprs_df <- data.frame(t(assay(sce,'exprs')))
#exprs_df$sample_id <- colData(sce)$sample_id
#exprs_df$cluster <- gsub("^\\d+\\s*-\\s*", "", cluster_ids(sce, k= 'Cluster2'))
#exprs_df$cluster <- cluster_ids(sce, k= 'Cluster3')

# summary df contains the median per sample_id & cluster combo
# summary_df <- exprs_df %>%
#   group_by(sample_id, cluster) %>%
#   summarise_all(.funs = median) 

# counts <- exprs_df %>%
#   group_by(sample_id, cluster) %>%
#   count() 

# summary_df %>%
#   select(-sample_id) %>%
#   group_by(cluster) %>%
#   summarise_all(.funs = median)
#rm(exprs_df)
# 
count_data <- read_excel('count_data.xlsx')
summary_df <- read_excel('summary_df.xlsx')
markers <- colnames(summary_df)[-c(1,2)]

count_data$cluster <- gsub("NaÃ¯ve", "Naïve", count_data$cluster)
count_data$cluster <- gsub("naÃ¯ve",'naïve',count_data$cluster)

limma_df <-summary_df %>%
  mutate(time = as.numeric(str_extract(sample_id,"(?<=\\s)\\d+(?=h)"))) %>%
  filter(time > 0, cluster != 'Innate lymphoid cells')

clusters <- c("Blood neutrophils (immature)",
              "Blood neutrophils (mature)",                
              "Blood neutrophils (aged)",                   
              "Tissue neutrophils (first)",                
              "Tissue neutrophils (second)",  
              "Tissue neutrophils (third)",
              "Tissue neutrophils (fourth)",               
              "Tissue neutrophils (fifth)",                 
              "Blood naïve B cells",                        
              "Blood memory B cells",                      
              "Tissue B1 cells",                            
              "Tissue activated B cells",                  
              "Naïve CD4 T cells",                         
              "Blood effector CD4 T cells",                
              "Tissue effector CD4 T cells",                
              "Naïve CD8 T cells",                         
              "Blood effector CD8 T cells",                
              "Tissue effector CD8 T cells",    
              "Conventional blood dendritic cells",         
              "Conventional type 1 tissue dendritic cells",
              "Conventional type 2 tissue dendritic cells", 
              "Monocyte-derived tissue dendritic cells" ,  
              "Naïve NK cells" ,                            
              "Activated NK cells",       
              "Tissue-resident macrophages" ,               
              "Infiltrating tissue macrophages (first)", 
              "Infiltrating tissue macrophages (second)" , 
              "Infiltrating tissue macrophages (third)" ,  
              "Classical blood monocytes" ,                 
              "Non-classical blood monocytes",  
              "Eosinophils")

limma_df$cluster <- gsub("NaÃ¯ve", "Naïve", limma_df$cluster)
limma_df$cluster <- gsub("naÃ¯ve",'naïve',limma_df$cluster)

limma_df$cluster <- factor(limma_df$cluster, levels = clusters)
limma_df$loc <- sub("(\\w+).*", "\\1", limma_df$sample_id)
limma_df$cluster_loc <- paste(limma_df$cluster, limma_df$loc)


# count_data %>%
#   mutate(cluster_type = case_when(
#     grepl('NK',cluster) ~ 'NK',
#     grepl('B cell|B1', cluster) ~ 'B cell',
#     grepl('CD4', cluster) ~ 'CD4 T cell',
#     grepl('neutrophil', cluster) ~ 'Aneutrophil',
#     grepl('dend', cluster) ~ 'DC',
#     grepl('CD8', cluster) ~ 'CD8 T cell',
#     grepl('macro', cluster) ~ 'macrophaghes',
#     grepl('monocyt',cluster) ~ 'monocytes',
#     grepl('Eosi',cluster) ~ 'xeosinophils',
#     TRUE ~ NA
#   )) %>%
#   arrange(cluster_type) %>%
#   distinct(cluster) %>%
#   pull()
  
# Plot showing why we can't use raw counts as weights in limma.
# count_data %>%
#   group_by(time) %>%
#   summarise(total_per_hour = sum(n)) %>%
#   ggplot(mapping = aes(x = factor(time), y=total_per_hour))+
#   geom_col() +
#   labs(title = 'Total cell counts over time', x = 'time (h)', y='cell count')
  
clusters_loc <- unique(limma_df$cluster_loc)

# can't use raw counts, so normalizing for total condition 
count_data$cluster <- factor(count_data$cluster, levels = clusters)
count_data <- count_data %>% 
  filter(time > 0) %>%
  arrange(cluster) %>%
  group_by(time, drug_type, location) %>%
  mutate(total_cells_per_condition = sum(n)) %>%
  ungroup() %>%
  mutate(sample_id_condition_prop = sample_id_sum/ total_cells_per_condition)  

#any(count_data$sample_id_condition_prop >1) # FALSE - good.
  
weights <- count_data %>%
  dplyr::select(all_of(c('cluster','sample_id_condition_prop', 'sample_id')))
  
limma_df <- merge(limma_df, weights, by= c('sample_id','cluster'))
limma_df <- limma_df %>% arrange(cluster)

top <- data.frame()
for (cl in clusters_loc){ 
  print(cl)
  temp_df <- limma_df %>%
    filter(cluster_loc == cl) 
  times <- as.factor(as.numeric(str_extract(temp_df$sample_id,"(?<=\\s)\\d+(?=h)")))
  drugs <- sub(".*\\b(\\w+)\\s+\\d+", "\\1", temp_df$sample_id)
  Treat <- paste(drugs,times,sep ='.')
  design <-  model.matrix(~0+Treat)
  colnames(design) <- sub("^Treat", "", colnames(design))
  limma_df2 <- t(temp_df[markers])
  colnames(limma_df2) <- Treat
  
  # Computing weights, should theoretically be the same between loops. 
  limma_weights <- limma_df %>%
    filter(cluster_loc == cl) %>%
    select(sample_id_condition_prop) %>%
    t()
  colnames(limma_weights) <- Treat
  
  cm <- makeContrasts(
    DrugvsControl04h = drug.4- control.4,
    DrugvsControl24h = drug.24- control.24,
    DrugvsControl48h = drug.48 -control.48,
    DrugvsControl72h = drug.72 - control.72,
    DrugvsNaproxen04h = drug.4 - naproxen.4,
    DrugvsNaproxen24h = drug.24 - naproxen.24,
    DrugvsNaproxen48h = drug.48 - naproxen.48,
    DrugvsNaproxen72h = drug.72 - naproxen.72,
    ControlvsNaproxen04h = control.4 - naproxen.4,
    ControlvsNaproxen24h = control.24 - naproxen.24,
    ControlvsNaproxen48h = control.48 - naproxen.48,
    ControlvsNaproxen72h = control.72 - naproxen.72,
    # interaction effects -> what markers are DE at 24h that are not at 0h?
    # DrugvsControl24vs4 = (drug.24-control.24) - (drug.4-control.4),
    # # DrugvsControl24h - DrugvsControl04h
    # DrugvsControl48vs24 = (drug.48-control.48) - (drug.24-control.24),
    # DrugvsControl72vs48 = (drug.72-control.72) - (drug.48-control.48),
    # DrugvsNaproxen24vs4 = (drug.24 - naproxen.24) - (drug.4-naproxen.4),
    # DrugvsNaproxen48vs24 = (drug.48 - naproxen.48) - (drug.24 - naproxen.24),
    # DrugvsNaproxen72vs48 = (drug.72 - naproxen.72) - (drug.48 - naproxen.48),
    levels = design
  )
  
  
  fit <- lmFit(limma_df2, design, weights = limma_weights) 
  fit <- contrasts.fit(fit, cm)
  efit <- eBayes(fit, trend = T)
  cm_cols <- colnames(cm)
  
  # we have 12 contrasts 
  for (coef in 1:12){
    temp_table <- limma::topTable(efit, coef = coef, number = Inf, adjust.method = "BH", sort.by = "none")
    temp_table$contrast = cm_cols[coef]
    temp_table$cluster = cl
    if (coef == 1){
      temp_markers <- rownames(temp_table)
      temp_table$marker = temp_markers
    }
    else{
      temp_table$marker = temp_markers # without the numbers after
    }
    #top <- rbind(top, topTable(efit, coef = coef, number = Inf, adjust.method = "BH", sort.by = "none"))
    top <- rbind(top, temp_table)
    }
}

write.csv(top, 'limma_results_with_weights2.csv', row.names = F)

for (cl in clusters_loc){
  print(cl)
  if (!grepl('./newDE_plots', getwd())){
    setwd('./newDE_plots')
  }
  #p <- medians_df %>%
  p <- limma_df %>%
    filter(cluster_loc == cl) %>%
    select(all_of(c(markers,'sample_id'))) %>%
    pivot_longer(markers, values_to = 'median_exprs', names_to = 'marker') %>%
    mutate(times = as.numeric(str_extract(sample_id,"(?<=\\s)\\d+(?=h)")),
           drugs = sub(".*\\b(\\w+)\\s+\\d+", "\\1", sample_id)) %>%
    filter(times  > 0) %>%
    ggplot(mapping = aes(x= as.factor(times), y= median_exprs, fill = drugs))+
    geom_boxplot() +
    facet_wrap(~marker, scales = 'free_y')+
    labs(title = cl, x = 'time (h)', y= 'expression', fill = 'treatment')
  print(p)
  ggsave(paste(cl, 'boxplot2.png',sep='.'), p, width=1920, height=1009, units='px',dpi=96)

}


lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.01")
count_data$cluster_loc <- paste(count_data$cluster, count_data$location)


plot_heatmap <- function(cl_list, dir_to_save, var_to_plot = 'adj.P.val'){
  # cl_lis : list of clusters
  # dir to save: what dir to save the plots in 
  # var_to_plot: either p_val or t_val (which tells us if up or down)
  
  if (!grepl(dir_to_save, getwd())){
    setwd(dir_to_save)
  }
  
  for (cl in cl_list){
    heatmap_df <- top %>% 
      filter(cluster == cl) %>%
      select(all_of(c('marker',var_to_plot, 'contrast'))) %>%
      pivot_wider(values_from = .data[[var_to_plot]], names_from = marker) %>%
      column_to_rownames(var = 'contrast')
    
    heatmap_data <- as.matrix(heatmap_df)
    col_fun = colorRamp2(c(1, 0.5, 0), c("red", "white", "green"))
    
    if (var_to_plot == 'adj.P.val'){
      ht = Heatmap(heatmap_data, name = "p value", column_title = cl,
                   col = col_fun, 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(heatmap_data[i, j] < 0.01)
                       grid.text('*', x, y, gp = gpar(fontsize = 10))})
    }
    
    else if (var_to_plot == 't'){
      ht = Heatmap(heatmap_data, name = "t value", column_title = cl,
                   col = col_fun)
    }
    
    
    png(file = paste(cl,'heatmap.png', sep='.'), width = 1920, height = 1009, units = "px", res = 96)  # Set the file name and other parameters
    draw(ht, annotation_legend_list = lgd_sig)
    dev.off()
  }
  setwd('~/')
}


blood_clusters <- character(0)
tissue_clusters <- character(0)

for (string in clusters_loc) {
  words <- strsplit(string, "\\s+")[[1]]
  last_word <- tail(words, 1)
  
  if (last_word == "Blood") {
    blood_clusters <- c(blood_clusters, string)
  } else if (last_word == "Peritoneum") {
    tissue_clusters <- c(tissue_clusters, string)
  }
}

print(blood_clusters)
print(tissue_clusters)

plot_heatmap(blood_clusters, './newDE_plots')
plot_heatmap(tissue_clusters,'./peritoneum_plots')

# top %>%
#   filter(adj.P.Val < p_cutoff) %>%
#   group_by(cluster, marker) %>%
#   count() %>%
#   filter(n>2) %>%
#   arrange(desc(n)) %>%
#   print(n= (Inf))
# 
# nrow(top %>% filter(adj.P.Val < 0.01)) # 587
# 
# top2 <- top %>%
#   mutate(loc = word(cluster, -1),
#          cluster = str_remove(cluster, "\\s+\\w+$")) %>%
#   filter(adj.P.Val < p_cutoff) 
# 
# cd8_plot <- limma_df %>%
#   filter(cluster_loc == 'Naïve CD8 T cells Blood') %>%
#   mutate(drug_type = sub(".*\\b(\\w+)\\s+\\d+", "\\1", sample_id)) %>%
#   ggplot(mapping = aes(x = factor(time), y = CD8a, color= drug_type))+
#   geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), 
#            aes(fill = drug_type), width=0.8) +
#   geom_point(position = position_dodge(0.8))+
#   scale_color_manual(values = c('black','black','black'), guide='none')+
#   labs(title = 'CD8a expression in naïve CD8+ T cells from the blood', x = 'time (h)', y = 'CD8a expression', fill = 'Treatment')+
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  


