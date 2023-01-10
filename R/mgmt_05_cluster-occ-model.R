###
# Occurence probability management for clusters #########
###

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

cluster_counts <- readRDS('./Data/03_cluster-conc.rds') |> 
  list_to_tib('profileid')

# |- Set-up -----------------

cluster_counts$detected <- ifelse(cluster_counts$conc_m3 > 0,
                                  1, 0)


uvp_vol <- uvp_data |>
  get_ecopart_vol() |>
  lapply(ecopart_vol_bin,depth_breaks = seq(0,1200,20)) |> 
  list_to_tib('profileid')

# cluster counts merger with volS
ck_df <- cluster_counts |> 
  left_join(uvp_vol) |> 
  left_join(uvp_data$meta[,c('cruise_id','tod','profileid')]) |> 
  bin_format()

# trim out non-sampled regions
ck_df <- ck_df[!is.na(ck_df$vol_sampled),]


# split by cluster
cluster_list <- vector('list', length(unique(ck_df$group)))
for(i in 1:length(cluster_list)) {
  temp_list <- ck_df |> 
    filter(group == i) |> 
    split(ck_df$tod[ck_df$group == i]) |> 
    lapply(function(x) split(x, x$db)) |> 
    lapply(function(x) lapply(x,surv_assign))
  
  cluster_list[[i]][['day']] <- temp_list[['day']] |> 
    lapply(listify_matrix)
  cluster_list[[i]][['night']] <- temp_list[['night']] |> 
    lapply(listify_matrix)
  
  rm(temp_list,i)
}

# |- Model Data Formatting -------------
# The depth-bin specific approach requires that there are separate models for each
# time of day and each cluster

# Here, they must be formatting in a big list

detection_matrix_list <- vector('list', 4)
vol_matrix_list <- vector('list', 4)
db_list <- vector('list', 4)

# for later rbind, all cols must be same legth
# to know what to ask for in max fill:
get_mcl <- function(tod, i) {
  x <- cluster_list[[i]][[tod]] |> 
    lapply(`[[`, 'detections') |> 
    sapply(ncol) |> 
    max()
  return(x)
}

get_detect <- function(tod, i) {
  x <- cluster_list[[i]][[tod]] |> 
    lapply(`[[`, 'detections') |> 
    lapply(max_fill, max_col_length[[tod]]) |> 
    do.call(what = rbind,)
  return(x)
}

get_vol <- function(tod, i) {
  x <- cluster_list[[i]][[tod]] |> 
    lapply(`[[`, 'vol_s') |> 
    lapply(max_fill, max_col_length[[tod]]) |> 
    do.call(what = rbind,)
  return(x)
}


get_db <- function(tod, i) {
  row_len <- cluster_list[[i]][[tod]] |> 
    lapply(`[[`, 'detections') |> 
    sapply(nrow)
  
  db_out <- NULL
  for(j in 1:length(row_len)) {
    db_out <- c(db_out, rep(names(row_len)[j], row_len[j]))
  }
  return(db_out)
}


# |-|- Loop for holder lists ---------------------  
# loop through to create container matrices for each big list

for(i in 1:length(db_list)) {
  for(tod in names(cluster_list[[i]])) {
    max_col_length <- list(day = NA, night = NA)
    max_col_length[[tod]] <- get_mcl(tod,i)
    detection_matrix_list[[i]][[tod]] <- get_detect(tod,i)
    vol_matrix_list[[i]][[tod]] <- get_vol(tod, i)
    db_list[[i]][[tod]] <- get_db(tod, i)
  }
}

# |-|- Model Data Finalize -------------------------

# To keep code concise, I'm creating 3 lists for each db region
# a list will hold the model data for each cluster which each tod

get_model_data <- function(i, tod, min_d, max_d) {
  
  db_code <- code_db(tod, min_d, max_d, db_list[[i]])
  
  mod_data <- list(
    Y = detection_matrix_list[[i]][[tod]][db_code$idx,],
    n = rowSums(!is.na(detection_matrix_list[[i]][[tod]][db_code$idx,])),
    nSites = nrow(detection_matrix_list[[i]][[tod]][db_code$idx,]),
    z = ifelse(rowSums(detection_matrix_list[[i]][[tod]][db_code$idx,], 
                       na.rm = T)>0,
               1,NA),
    vs = standardize(vol_matrix_list[[i]][[tod]][db_code$idx,]),
    db = db_code$code - (min(db_code$code) - 1),
    dbLength = max(db_code$code - (min(db_code$code) - 1))
  )
  
  return(mod_data)
}

epi_mod_data <- vector('list', 4)
for(i in 1:length(epi_mod_data)) {
  for(tod in c('day','night')) {
    epi_mod_data[[i]][[tod]] <- get_model_data(i, tod, 0, 200)
  }
}

tmeso_mod_data <- vector('list',4)
for(i in 1:length(tmeso_mod_data)) {
  for(tod in c('day','night')) {
    tmeso_mod_data[[i]][[tod]] <- get_model_data(i, tod, 200, 600)
  }
}

bmeso_mod_data <- vector('list', 4)
for(i in 1:length(bmeso_mod_data)) {
  for(tod in c('day', 'night')) {
    bmeso_mod_data[[i]][[tod]] <- get_model_data(i, tod, 600, 1200)
  }
}

###
# Run JAGS MODELS ##################
###

# |- Format iterative way to run it --------------
param_out <- c('p0','a0','a_vs','psi','N')

run_model <- function(mod_data) {
  mod <- jags(mod_data,
              parameters.to.save = param_out,
              model.file = './R/jags_04_occ-detection-only.jags',
              DIC = F, n.chains = 3, n.iter = 10000,
              n.thin = 2, parallel = T)
  return(mod)
}

epi_mod_output <- vector('list',4)
for(i in 1:length(epi_mod_data)) {
  for(tod in c('day','night')) {
    epi_mod_output[[i]][[tod]] <- run_model(epi_mod_data[[i]][[tod]])
  }
}

tmeso_mod_output <- vector('list',4)
for(i in 1:length(tmeso_mod_data)) {
  for(tod in c('day','night')) {
    tmeso_mod_output[[i]][[tod]] <- run_model(tmeso_mod_data[[i]][[tod]])
  }
}

bmeso_mod_output <- vector('list',4)
for(i in 1:length(bmeso_mod_data)) {
  for(tod in c('day','night')) {
    bmeso_mod_output[[i]][[tod]] <- run_model(bmeso_mod_data[[i]][[tod]])
  }
}

# |- Save Model Output ------------------------------

saveRDS(db_list, './data/05_db-list.rds')
saveRDS(epi_mod_output, './data/05_cluster-occ-mod-epi.rds')
saveRDS(tmeso_mod_output, './data/05_cluster-occ-mod-tmeso.rds')
saveRDS(bmeso_mod_output, './data/05_cluster-occ-mod-bmeso.rds')
