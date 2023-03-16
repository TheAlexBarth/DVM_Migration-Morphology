####
# Weighted PCA on Copepod Morphology #####################
####

library(EcotaxaTools)
library(FactoMineR)
source('./R/tools.R')


###
# Set up ###############
###

# |- Import and clean data ---------------------------------

uvp_data <- readRDS('./data/01_uvp-trim-final_large.rds')


# # What's the relative abundance of copepods across multicelluar zoops
# rel_cope <- uvp_data |>
#   merge_casts(list(all_casts = names(uvp_data$zoo_files))) |>
#   mod_zoo(names_drop, c('not-living', 'Trichodesmium','Rhizaria', 'darksphere',
#                         'house','temp circle','duplicate'),
#           drop_children = T) |>
#   add_zoo(names_to, col_name = 'name',
#           new_names = c('Copepoda', 'Eumalacostraca',
#                                   'Chaetognatha','Ostracoda','living',
#                                   'Hydrozoa','Annelida','Pteropoda',
#                                   'Actinopterygii','Crustacea')) |>
#   rel_taxa()
# rel_cope$all_casts[order(rel_cope$all_casts$rel_abundance, decreasing = T),]
# # we see copepods to be 60% of all observed taxa, 14.2% living misc, 
# # and 6%ish of Eumalacostraca, Ostracods, Chaetognaths
# 
# crust_data <- trim_to_crust(uvp_data)
uvp_data <- trim_to_cope(uvp_data)

# |- Pair weights to each observation --------------------

# we're going to weight by volume sampled to account for any duplicates
# weights will be applied only in bins where it is oversampled
# frame height is 3.11cm so max volume sampled is 35.36977
max_vol <- (1/.0311) * 1.1

#get the sampled volumes for each cast-bin wise
sampled_vols <- get_ecopart_vol(uvp_data)
# crust_vols <- get_ecopart_vol(crust_data)

# wrap in an if to safety check - should never be an issue unless files are whack
if(all(names(sampled_vols) == names(uvp_data$zoo_files))) {
  
  for(i in 1:length(uvp_data$zoo_files)) {
    #round observed depth to volume sampled depths
    obs_depths <- uvp_data$zoo_files[[i]]$depth_including_offset |> 
      sapply(nearest, sampled_vols[[i]]$depth)
    
    #index to the sampled volumes
    d_idx <- match(obs_depths, sampled_vols[[i]]$depth)
    
    # get volumed sampled in each location
    v_sampled <- sampled_vols[[i]]$vol_sampled[d_idx]
    v_sampled[which(v_sampled < max_vol)] <- max_vol #replace less than max_vol
    
    # assign weights based on volume sampled
    uvp_data$zoo_files[[i]]['weight'] <- 1/v_sampled
  }
}

# # wrap in an if to safety check - should never be an issue unless files are whack
# if(all(names(crust_vols) == names(crust_data$zoo_files))) {
#   
#   for(i in 1:length(crust_data$zoo_files)) {
#     #round observed depth to volume sampled depths
#     obs_depths <- crust_data$zoo_files[[i]]$depth_including_offset |> 
#       sapply(nearest, crust_vols[[i]]$depth)
#     
#     #index to the sampled volumes
#     d_idx <- match(obs_depths, crust_vols[[i]]$depth)
#     
#     # get volumed sampled in each location
#     v_sampled <- crust_vols[[i]]$vol_sampled[d_idx]
#     v_sampled[which(v_sampled < max_vol)] <- max_vol #replace less than max_vol
#     
#     # assign weights based on volume sampled
#     crust_data$zoo_files[[i]]['weight'] <- 1/v_sampled
#   }
# }

# |- Format into a long df ------------------------
all_copes <- uvp_data$zoo_files |> 
  list_to_tib('profileid')

# all_crust <- crust_data$zoo_files |> 
#   list_to_tib('profileid')

###
# PCA ####
###

# |- Features ----

# Description:
# --SIZE--
# major - the major axis of fitted ellipse
# minor - the minor axis of fitted ellipse
# area - the area of pixels
# feret - max distance between two points on perimeter
# perim. - objects perimeter
# --SHAPE--
# circ. - circularity measure (4pi*area / perim)
# elongation - major/minor
# thickr - thickness ratio; max thickness/mean thickness
# symetriev - bilateral symmetry
# symetrievc - bilateral symmetry for >25th percentile of grey pixels
# --TRANS--
# mean - mean grey value
# median - median grey value
# histcum1 - 25th percentile of grey value
# stddev - stddev of grey value
# skew - skewness of histrgram of grey values
# --COMPLEX--
# perimmajor - perim/major : higher for complex
# perimferet - perim/feret
# fractal - fractal dimension of perimeter

features <- c('major','minor', 'area','feret','perim.',
              'circ.','elongation', 'thickr', 'symetriev',
              'symetrievc','mean','median','histcum1',
              'stddev','skew','perimmajor','perimferet',
              'fractal')

# get a df of just interested features
cope_features <- all_copes[,which(names(all_copes) %in% features)]

# |- Run PCA --------------------------------------------

set.seed(20230314)
cope_pca <- PCA(cope_features, scale.unit = TRUE,
                row.w = all_copes$weight, graph = FALSE)


# |- Loading scores visualize ---------------------------------
# cope_pca$var$coord[order(cope_pca$var$coord[,1], decreasing = T),]
# dim.1 is characterized by size and and anti-correlatedd with circularity:: size
# dim.2 is positive with skewed grey values and anti- with darker values:: transparency
# dim.3 is positive with complexity and anti- with elongation:: appendage visibility
# dim.4 is anti- with symmetry and positive with minor:: orientation (weak)

# |- Save PCA -----------------------------------------
saveRDS(cope_pca, './data/02_cope-pca-res.rds')


# |- Bringing it all together ----
results <- data.frame(orig_id = all_copes$orig_id,
                      profileid = all_copes$profileid,
                      PC1 = cope_pca$ind$coord[,1],
                      PC2 = cope_pca$ind$coord[,2])

# |- Quick plot for fun ----
# plot(results$PC1, results$PC2, col = results$cluster)

# |- Save Results
saveRDS(results,'./data/02_PCA-copepods.rds')
