####
# Occupancy model approach ########
####

rm(list = ls())

library(EcotaxaTools)
library(jagsUI)
library(dplyr)
library(readr)
library(lubridate)
library(wiqid)
library(tidyr)
source('./R/tools.R')

uvp_data <- readRDS('./data/01_uvp-trim-final.rds') |> 
  trim_to_cope()

uvp_data$meta$cruise_id <- paste0(uvp_data$meta$stationid, 
                                  uvp_data$meta$programid)

####
# All Copepod Full Column Model ######
####

# |- Copepod Organization ------------------
# Get copepod counts
copepod_counts <- uvp_data$zoo_files |> 
  lapply(bin_taxa,
         depth_breaks = seq(0,1200,20), 
         force_bins = T)


copepod_counts <- copepod_counts |> 
  lapply(was_detected)



# |- Volume Data-----------

# this function returns in L
uvp_vol <- uvp_data |>
  get_ecopart_vol() |>
  lapply(ecopart_vol_bin,depth_breaks = seq(0,1200,20)) |> 
  list_to_tib('profileid')


# |- Full Organization -------------

# copepod count data frame
cc_df <- copepod_counts |> 
  list_to_tib('profileid') |> 
  left_join(uvp_vol) |> 
  left_join(uvp_data$meta[,c('cruise_id','tod','profileid')]) |> 
  bin_format()

# set db to factor
cc_df$db <- cc_df$db |> factor(levels=unique(cc_df$db))

#trim no-sampled regions
cc_df <- cc_df[-which(is.na(cc_df$vol_sampled)),]



# |- Full Column Model ---------------

# set up data into matrix format

cc_df_split <- cc_df |>
  split(f = cc_df$tod) |> 
  lapply(function(x) split(x, x$db)) |> 
  lapply(function(x) lapply(x,surv_assign))

# organize by matrices

detection_list <- list()
detection_list[['day']] <- cc_df_split$day |> 
  lapply(listify_matrix)
detection_list[['night']] <-  cc_df_split$night |> 
  lapply(listify_matrix)

# |- organizing for jags layout -------------

detection_matrix <- list()
detection_matrix$day <- detection_list$day |> 
  lapply(`[[`,'detections') |> 
  lapply(max_fill,19) |> 
  do.call(what = rbind,)
detection_matrix$night <- detection_list$night |> 
  lapply(`[[`,'detections') |> 
  lapply(max_fill, 9) |> 
  do.call(what = rbind,)

# This creates a db as a factor for each row of the detection matrix
# the approach is a bit of hard-coding because you need to verify 
# that that replication numbers are equivalent to the rows
get_db <- function(tod) {
  row_len <- detection_list[[tod]]|> 
    lapply(`[[`, 'detections') |> 
    sapply(nrow)
  
  db_out <- NULL
  for(j in 1:length(row_len)) {
    db_out <- c(db_out, rep(names(detection_list[[tod]])[j], row_len[j]))
  }
  return(db_out)
}

db <- list(
  day = get_db('day'),
  night = get_db('night')
)


vol_matrix <- list()
vol_matrix$day <- detection_list$day |> 
  lapply(`[[`,'vol_s') |> 
  lapply(max_fill, 19) |> 
  do.call(what = rbind)
vol_matrix$night <- detection_list$night |> 
  lapply(`[[`,'vol_s') |> 
  lapply(max_fill, 9) |> 
  do.call(what = rbind)

# |- Run for different section of column ------------------------------------


get_model_data <- function(tod,min_d, max_d, db = db) {
  
  db_code <- code_db(tod, min_d, max_d, db = db)

  mod_data <- list(
    Y = detection_matrix[[tod]][db_code$idx,],
    n = rowSums(!is.na(detection_matrix[[tod]][db_code$idx,])),
    nSites = nrow(detection_matrix[[tod]][db_code$idx,]),
    z = ifelse(rowSums(detection_matrix[[tod]][db_code$idx,], 
                       na.rm = T)>0,
               1,NA),
    vs = standardize(vol_matrix[[tod]][db_code$idx,]),
    db = db_code$code - (min(db_code$code) - 1),
    dbLength = max(db_code$code - (min(db_code$code) - 1))
  )
  
  return(mod_data)
}

epi_mod_day <- get_model_data('day',0,200, db = db)
epi_mod_night <- get_model_data('night', 0, 200, db = db)

param_out <- c('p0','a0','a_vs','psi','N')

day_post_epi <- jags(epi_mod_day,
                          parameters.to.save = param_out,
                          model.file = './R/jags_04_occ-detection-only.jags',
                          DIC = F, n.chains = 3, n.iter = 10000,
                          n.thin = 2, parallel = T)

night_post_epi <- jags(epi_mod_night,
                   parameters.to.save = param_out,
                   model.file = './R/jags_04_occ-detection-only.jags',
                   DIC = F, n.chains = 3, n.iter = 10000,
                    n.thin = 2, parallel = T)


saveRDS(
  list(db = db,day = day_post_epi,
           night = night_post_epi),
  file = './data/04_full-col-model-euphotic.rds'
)

topmeso_day <- get_model_data('day',200,600,db = db)
topmeso_night <- get_model_data('night',200, 600, db = db)

param_out <- c('p0','a0','a_vs','psi','N')

day_post_top <- jags(topmeso_day,
                 parameters.to.save = param_out,
                 model.file = './R/jags_04_occ-detection-only.jags',
                 DIC = F, n.chains = 3, n.iter = 10000,
                 n.thin = 2, parallel = T)

night_post_top <- jags(topmeso_night,
                   parameters.to.save = param_out,
                   model.file = './R/jags_04_occ-detection-only.jags',
                   DIC = F, n.chains = 3, n.iter = 10000,
                   n.thin = 2, parallel = T)


saveRDS(
  list(db = db,day = day_post_top,
       night = night_post_top),
  file = './data/04_full-col-model-top-meso.rds'
)


botmeso_day <- get_model_data('day',600,1200, db = db)
botmeso_night <- get_model_data('night', 600, 1200, db = db)

param_out <- c('p0','a0','a_vs','psi','N')

day_post_bot <- jags(botmeso_day,
                 parameters.to.save = param_out,
                 model.file = './R/jags_04_occ-detection-only.jags',
                 DIC = F, n.chains = 3, n.iter = 10000,
                 n.thin = 2, parallel = T)

night_post_bot <- jags(botmeso_night,
                   parameters.to.save = param_out,
                   model.file = './R/jags_04_occ-detection-only.jags',
                   DIC = F, n.chains = 3, n.iter = 10000,
                   n.thin = 2, parallel = T)


saveRDS(
  list(db = db,day = day_post_bot,
       night = night_post_bot),
  file = './data/04_full-col-model-bot-meso.rds'
)

