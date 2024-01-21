library(readxl)
library(flowCore)
library(CATALYST)
#library(cowplot)
library(diffcyt)
#library(corrplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(edgeR)
library(FlowSOM)
library(ConsensusClusterPlus)
library(purrr)

# SET WORKING DIRECTORY

setwd('C:/Users/Camel/OneDrive/Desktop/Gilroy')

load("Workspace.RData")
merging_table2 <- "./Excel_files/TE94_cluster_merging_2.xlsx"
merging_table2 <- read_excel(merging_table2)
sce <- mergeClusters(sce, k = "meta66", table = merging_table2, id = "Cluster2", overwrite = TRUE)

# merge_table3 <- './Excel_files/map3.xlsx'
# merge_table3 <- read_excel(merge_table3)
# sce <- mergeClusters(sce, k = "Cluster2", table = merge_table3, id = "Cluster3", overwrite = T)

cd4_df <- read.csv('cd4_df.csv')
cd4_df$location <- sub("(\\w+).*", "\\1", cd4_df$sample_id)

count_data <- read.csv('count_data.csv')

set.seed(2024)

cd4_sce <- filterSCE(sce, grepl('CD4',cluster_id), k ='Cluster2')
plotDR(cd4_sce, color_by = 'Cluster2')

# plot the expression of these 3 CD4+ clusters:
cd4_markers <- c('CD127', 'CD25')

exprs <- t(assay(cd4_sce, 'exprs'))[,t_cell_markers]
cd4_df$CD25 <- cd25_exprs

cd4_df <- cd4_df %>%
  group_by(sample_id, cluster) %>%
  mutate(median_cd25_per_sample_id = median(CD25)) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(total_sample_id_cell_count = n()) %>%
  ungroup() 

t_cell_activation_markers <- c('CD44','CD25', 'CD69','PD-1')
cd4_df <- cbind(cd4_df, exprs[,-c(2,3)])

# Blood vs tissue effector?
cd4_df %>%
  mutate(across(all_of(t_cell_markers), scale)) %>%
  group_by(cluster) %>%
  summarise_at(t_cell_markers, median) %>%
  pivot_longer(!cluster, names_to = 'marker',values_to = 'median_exprs') %>%
  ggplot(mapping = aes(x = marker,y=cluster, fill = median_exprs))+
  geom_tile()+
  geom_text(aes(label = round(median_exprs,1)))+
  scale_fill_gradient(low="white", high="blue")+
  theme_minimal()

cd4_df %>%
  group_by(cluster,location, condition) %>%
  mutate(cell_count = n()) %>%
  ggplot(mapping = aes(x=time, y= cell_count, color = drug_type))+
  geom_point()+
  geom_line()+
  facet_grid(location ~ cluster, scales='free_y') +
  labs(title = 'T cell count over time', x= 'time (h)', y = 'cell count')+
  scale_x_discrete(
       breaks = c(0,4,24,48,72),
       labels = c('0','4','24','48','72'))+
  theme_minimal()


# ggplot(data = cd4_df, mapping = aes(x = factor(time, levels = unique(cd4_df$time)), y = median_cd25_per_sample_id, 
#                                     fill = drug_type, size = total_sample_id_cell_count)) +
#   geom_point(shape=21, position = position_dodge(0.75)) +
#   facet_grid(location ~ cluster) +
#   labs(title = 'Median CD25 expression over time', x= 'time (h)', y = 'CD25 expression')+
#   scale_x_discrete(
#     breaks = c(0,4,24,48,72),
#     labels = c('0','4','24','48','72'))+
#   theme_minimal()




# Can't find a cluster which has a Treg population. Reclustering....----
t_cell_markers <- c('CD4','CD25','CD127','CD44','CD45','LAG3','PD-1','CD69','CD62L')
rowData(cd4_sce)$used_for_clustering <- rownames(cd4_sce) %in% t_cell_markers

fsom <- ReadInput(flowFrame(t(assay(cd4_sce, "exprs"))))
som <- BuildSOM(fsom, colsToUse=t_cell_markers, 
                silent=TRUE, xdim=10, ydim=10)
mc <- suppressWarnings(suppressMessages(
  ConsensusClusterPlus(t(som$map$codes), 
                       maxK = 50, reps = 100, 
                       distance = "euclidean", 
                       seed = 2024, plot = NULL)))

k <- 100
mcs <- seq_len(50)[-1]
# construct data.frame of clustering codes
codes <- data.frame(seq_len(k), map(mc[-1], "consensusClass"))
codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", mcs))
#cd4_sce <- runDR(cd4_sce, dr=c("UMAP"), cells = 1000, features = "state", seed = 2024)
flowsom_clusters <- factor(som$map$mapping[, 1])

.triangle <- function(m) {
  n <- ncol(m)
  nm <- matrix(0, ncol=n, nrow=n)
  fm <- m
  nm[upper.tri(nm)] <- m[upper.tri(m)]
  fm <- t(nm) + nm
  diag(fm) <-  diag(m)
  nm <- fm
  nm[upper.tri(nm)] <- NA
  diag(nm) <- NA
  m[lower.tri(nm)]
}

.plot_delta_area <- function(mc) {
  # empirical CDF distribution
  maxK <- length(mc)
  v <- lapply(mc[seq_len(maxK)[-1]], function(x) .triangle(x$ml))
  h <- lapply(v, function(x) {
    h <- graphics::hist(x, breaks=seq(0, 1, .01), plot=FALSE)
    h$counts <- cumsum(h$counts) / sum(h$counts)
    return(h)
  })
  # calculate area under CDF curve, by histogram method &
  # calculate proportional increase relative to prior k
  areaK <- vapply(h, function(x) cumsum(x$counts * .01)[100], numeric(1))
  deltaK <- c(areaK[1], diff(areaK) / areaK[seq_len(maxK-2)])
  
  df <- data.frame(k=seq_len(maxK)[-1], y=deltaK)
  y_max <- ceiling(max(df$y)*2)/2
  ggplot(df, aes_string(x="k", y="y")) + 
    theme_classic() + geom_line(color="steelblue", lty=2) + 
    geom_point(size=2.5, color="navy") + coord_fixed(4) +
    scale_x_continuous(breaks=seq(2, 20, 2), expand=c(0,.5)) +
    scale_y_continuous(limits=c(0, y_max), expand=c(0,.125), 
                       breaks=function(x) seq(x[1]+.125, x[2], .5)) +
    ylab("Relative change\nin area under CDF curve") +
    theme(plot.title=element_text(face="bold"),
          axis.text=element_text(color="black"),
          panel.grid.major=element_line(color="grey", size=.2))
}

.plot_delta_area(mc) 

# PLOTS with the new clusters 

#cluster_4 = codes[flowsom_clusters,] #rearranges codes in order of the clusters
#cols_to_normalise <- c('median_cd25','median_cd25')

# cluster_df <- cd4_df %>%
#   mutate(cluster10 = cluster_4$meta20) %>%
#   select(all_of(c('cluster10', 'CD25', 'CD127'))) %>%
#   group_by(cluster10) %>%
#   summarise(median_cd25 = median(CD25),
#          median_cd127 = median(CD127)) %>%
#   mutate(across(all_of(cols_to_normalise), ~(.-min(.))/ (max(.) - min(.))))
# Normalise the medians ACROSS the clusters (or normalise the data first and then find medians)?

# cd4_df %>%
#   mutate(cluster10 = cluster_4$meta20) %>%
#   mutate(scaled_cd25 = scale(CD25),
#          scaled_cd127 = scale(CD127)) %>%
#   select(all_of(c('scaled_cd25','scaled_cd127'))) %>%
#   pivot_longer(everything(),names_to = 'marker', values_to = 'exprs') %>%
#   ggplot(mapping = aes(x=exprs, fill = marker))+
#   geom_histogram(bins=60, alpha = 0.5, position = 'identity')+
#   labs(title = 'Marker expression for CD4 T cells') +
#   theme_minimal()
# 
# cd4_df %>%
#   mutate(cluster20 = cluster_4$meta40) %>%
#   mutate(scaled_cd25 = scale(CD25),
#          scaled_cd127 = scale(CD127)) %>%
#   group_by(cluster20)%>%
#   summarise(median_scaled_cd25 = median(scaled_cd25),
#             median_scaled_cd127 = median(scaled_cd127),
#             num_cells_per_cluster  = n()/916452) %>%
#   #select(all_of(c('cluster20','scaled_cd25','scaled_cd127'))) %>%
#   pivot_longer(!c('cluster20','num_cells_per_cluster'),names_to = 'marker', values_to = 'median_exprs') %>%
#   ggplot(mapping = aes(x=cluster20,y= marker, fill = median_exprs))+
#   geom_tile()+
#   geom_text(aes(label = round(num_cells_per_cluster,2)))+ 
#   scale_fill_distiller(palette = "RdPu")+
#   labs(title = 'Marker expression for CD4 T cells') +
#   theme_minimal()
#   
# cd4_df %>%
#   select(all_of(c('CD25','CD127'))) %>%
#   pivot_longer(everything(),names_to = 'marker', values_to = 'exprs') %>%
#   ggplot(mapping = aes(x=exprs, fill = marker)) + 
#   geom_histogram(bins=60, alpha = 0.5, position = 'identity')+
#   labs(title = 'Marker expression for CD4 T cells') +
#   theme_minimal()



# PLOT to confirm the values match with the CATALYST plot ----
# plotPbExprs(cd4_sce, features= 'CD25',fun='median', facet_by = 'cluster_id', k ='Cluster2')
# # Try to recreate this boxplot:
# cond_levels <- cd4_df %>% 
#   mutate(cond = paste(location,condition,sep='.')) %>% 
#   distinct(cond) %>% pull()
# 
# cd4_df %>%
#   mutate(cond = paste(location,condition,sep='.')) %>%
#   group_by(cluster, sample_id) %>%
#   mutate(median_per_cluster_cd25 = median(CD25)) %>%
#   ungroup() %>%
#   ggplot(mapping = aes(x = factor(cond, levels = cond_levels) , y =median_per_cluster_cd25, color = factor(cond, levels = cond_levels))) +
#   geom_boxplot()+
#   facet_wrap(~factor(cluster, levels = unique(cd4_df$cluster)), scales = 'free_y')+
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# ----
