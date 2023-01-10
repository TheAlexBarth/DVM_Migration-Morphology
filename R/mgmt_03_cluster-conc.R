####
# Binned density by cluster #############
####

# |- Set up ----------------------------

rm(list = ls())
library(EcotaxaTools)
source('./R/tools.R')

cluster_cope <- readRDS('./data/02_cluster-copepods.rds')
uvp_data <- readRDS('./data/01_uvp-trim-final.rds') |> 
  trim_to_cope()

####
# Assign Clusters to individuals ########
####

for(i in 1:length(uvp_data$zoo_files)) {
  idx <- which(cluster_cope$orig_id %in% uvp_data$zoo_files[[i]]$orig_id)
  uvp_data$zoo_files[[i]]$cluster <- cluster_cope$cluster[idx]
  uvp_data$zoo_files[[i]]$PC1 <- cluster_cope$PC1[idx]
  uvp_data$zoo_files[[i]]$PC2 <- cluster_cope$PC2[idx]
}

###
# Calculate Concentrations #######
###

# this process can be slow
# run in 20m bins
start = Sys.time()
cluster_bins <- uvp_data |> uvp_zoo_conc(breaks = seq(0,1200,20),
                                         cat_col = 'cluster',
                                         func_col = 'cluster')
end = Sys.time()
end-start


saveRDS(cluster_bins, './data/03_cluster-conc.rds')
