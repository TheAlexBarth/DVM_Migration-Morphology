####
# Alternative morphogroup assessment #########
####

rm(list = ls())
library(EcotaxaTools)
source('./R/tools.R')


pca_cope <- readRDS('./data/02_PCA-copepods.rds')
uvp_data <- readRDS('./data/01_uvp-trim-final_large.rds') |> 
  trim_to_cope()

# Morphogroup definition ####

# |- Assign percentiles of morphogroups -----------------
pca_cope$PC1_percentile <- sapply(pca_cope$PC1, function(x) ecdf(pca_cope$PC1)(x))
pca_cope$PC2_percentile <- sapply(pca_cope$PC2, function(x) ecdf(pca_cope$PC2)(x))

# |- Assign morphogroups --------

morpho_switch <- function(x) {
  if(x < .25) {
    return('low')
  } else if (x < .75) {
    return('mid')
  } else {
    return('high')
  }
}
pca_cope$pc1_class <- sapply(pca_cope$PC1_percentile, morpho_switch)
pca_cope$pc2_class <- sapply(pca_cope$PC2_percentile, morpho_switch)


for(i in 1:length(uvp_data$zoo_files)) {
  idx <- which(pca_cope$orig_id %in% uvp_data$zoo_files[[i]]$orig_id)
  uvp_data$zoo_files[[i]]$pc1_class <- pca_cope$pc1_class[idx]
  uvp_data$zoo_files[[i]]$pc2_class <- pca_cope$pc2_class[idx]
}

# |- trim out based on PC1 ---------------------

mid_copes <- uvp_data |> 
  mod_zoo(func = function(x) x[which(x$pc1_class == 'mid'),])

# Find any casts with no-copepods
# Need to move to etxtools
drop_casts <- mid_copes$zoo_files |>
  sapply(nrow) |> 
  sapply(function(x) x == 0) |> 
  which() |> 
  names()

if(length(drop_casts) > 0) {
  mid_copes$par_files <- mid_copes$par_files[which(!(names(mid_copes$par_files) %in% drop_casts))]
  mid_copes$zoo_files <- mid_copes$zoo_files[which(!(names(mid_copes$zoo_files) %in% drop_casts))]
  mid_copes$meta <- mid_copes$meta[which(!(mid_copes$meta$profileid%in%drop_casts)),]
}
# restore class structure
mid_copes <- as_ecopart_obj(mid_copes)


# |- Trim out for large copepods based on pc2 ---------------------
big_copes <- uvp_data |> 
  mod_zoo(func = function(x) x[which(x$pc1_class == 'high'),])

# Find any casts with no-copepods
drop_casts <- big_copes$zoo_files |>
  sapply(nrow) |> 
  sapply(function(x) x == 0) |> 
  which() |> 
  names()

if(length(drop_casts) > 0) {
  big_copes$par_files <- big_copes$par_files[which(!(names(big_copes$par_files) %in% drop_casts))]
  big_copes$zoo_files <- big_copes$zoo_files[which(!(names(big_copes$zoo_files) %in% drop_casts))]
  big_copes$meta <- big_copes$meta[which(!(big_copes$meta$profileid%in%drop_casts)),]
}
# restore class structure
big_copes <- as_ecopart_obj(big_copes)

# |- Trim out PC2 based on small copepods -------------------------
sm_copes <- uvp_data |> 
  mod_zoo(func = function(x) x[which(x$pc1_class == 'low'),])

# Find any casts with no-copepods
drop_casts <- sm_copes$zoo_files |>
  sapply(nrow) |> 
  sapply(function(x) x == 0) |> 
  which() |> 
  names()

if(length(drop_casts) > 0) {
  sm_copes$par_files <- sm_copes$par_files[which(!(names(sm_copes$par_files) %in% drop_casts))]
  sm_copes$zoo_files <- sm_copes$zoo_files[which(!(names(sm_copes$zoo_files) %in% drop_casts))]
  sm_copes$meta <- sm_copes$meta[which(!(sm_copes$meta$profileid%in%drop_casts)),]
}
# restore class structure
sm_copes <- as_ecopart_obj(sm_copes)


# Calculate group concentration #####

start = Sys.time()
pc1_bins <- uvp_data |> uvp_zoo_conc(breaks = seq(0,1200,20),
                                         cat_col = 'pc1_class',
                                         func_col = 'pc1_class')

pc2_mid_bins <- mid_copes |> uvp_zoo_conc(breaks = seq(0,1200,20),
                                          cat_col = 'pc2_class',
                                          func_col = 'pc2_class')

pc2_lg_bins <- big_copes |> uvp_zoo_conc(breaks = seq(0,1200,20),
                                         cat_col = 'pc2_class',
                                         func_col = 'pc2_class')

pc2_sm_bins <- sm_copes |> uvp_zoo_conc(breaks = seq(0,1200,20),
                                        cat_col = 'pc2_class',
                                        func_col = 'pc2_class')


end = Sys.time()
end-start

saveRDS(pc1_bins, './data/03b_pc1-bins.rds')
saveRDS(pc2_mid_bins, './data/03b_pc2-bins.rds')
saveRDS(pc2_lg_bins, './data/03b_pc2-lg-bins.rds')
saveRDS(pc2_sm_bins, './data/03b_pc2-sm-bins.rds')
