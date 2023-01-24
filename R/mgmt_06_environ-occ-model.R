###
# Occurrence probability by environmental factors #########
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

uvp_data <- readRDS('./data/01_uvp-trim-final_large.rds') |> 
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

####
# Environmental Data Set-up ##############
####


ctd_data <- readRDS('./data/01_ctd-trim-final.rds')

# |-  DCM Location -------------

dcm_max <- function(df) {
  df$Depth[which.max(df$RFU)]
}

dcm_depth <- tibble(
  ctd_origfilename = as.integer(names(ctd_data)),
  dcm_d = ctd_data |> sapply(dcm_max)
) |>
  left_join(uvp_data$meta[,c("ctd_origfilename", 'profileid', 'sampledate')])

# |-Integrated particle abundance ----------------------------

par_conc <- uvp_data |> uvp_par_conc(max_esd = .45)

intg_up250 <- function(df) {
  trim_df <- df[df$depth <= 250,]
  return(sum(trim_df$par_conc))
}

intg_particle <- tibble(
  profileid = names(par_conc),
  intg_part = par_conc |> sapply(intg_up250)
)

cast_metrics <- full_join(dcm_depth,
                          intg_particle,
                          by = 'profileid')
cast_metrics$sampledate <- as.Date(cast_metrics$sampledate)

# |- Merging MODIS and CTD data ------------------
modis_raw <- readRDS('./data/00_modis-data.rds')

final_metrics <- cast_metrics |>
  left_join(modis_raw$par, by = c('sampledate' = 'UTC')) |>
  left_join(modis_raw$dac, by = c('sampledate' = 'UTC'))

# |- Correction for modis data ---------------
# For some sets of modis data, there are missing days
# In these cases, we will map the nearest available day
nan_idx <- which(is.nan(final_metrics$dac))

for(i in 1:length(nan_idx)) {
  
  available_dates_par <- modis_raw$par$UTC[which(!(is.nan(modis_raw$par$par)))]
  
  
  #get nearest availabe dates
  nearest_available_par <- nearest(final_metrics$sampledate[nan_idx[i]],
                                   available_dates_par)
  
  if(abs(nearest_available_par - final_metrics$sampledate[nan_idx[i]]) < 3) {
    final_metrics$par[nan_idx[i]] <- modis_raw$par$par[which(modis_raw$par$UTC ==
                                                               nearest_available_par)]
    
  }
  
  available_dates_dac <- modis_raw$dac$UTC[which(!(is.nan(modis_raw$dac$dac)))]
  
  nearest_available_dac <- nearest(final_metrics$sampledate[nan_idx[i]],
                                   available_dates_dac)
  
  if(abs(nearest_available_dac - final_metrics$sampledate[nan_idx[i]]) < 3) {
    final_metrics$dac[nan_idx[i]] <- modis_raw$dac$dac[which(modis_raw$dac$UTC ==
                                                               nearest_available_dac)]
  }
}

# remove one instance of bad ctd cast and merge in cruise date then summarize
environ_summary <- final_metrics[!is.na(final_metrics$profileid),] |> 
  left_join(uvp_data$meta[,c('cruise_id','tod','profileid')]) |> 
  group_by(cruise_id, tod) |> 
  summarise(dcm_d = mean(dcm_d), intg_part = mean(intg_part), 
            par = mean(par), dac = mean(dac))


####
# Formatting Detection and Environmental ####
####

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

# add environmental vars to cluster list
# This is definitely not the cleanest way to do this but it is just stacking
# a house of cards to this stupid list structure I've made.
# In future, this might be cleaner but here - I need to repeat the data for the
# matching dbs and not all dbs have the same cruises visiting that "site"

add_environ <- function(x, tod) {
  x[['env_data']] <- data.frame(cruise_id = x$cruise_id) |> 
    left_join(environ_summary[environ_summary$tod == tod,]) |> 
    suppressMessages()
  
  x$cruise_id <- NULL
  return(x)
}

for(i in 1:length(cluster_list)) {
  for(tod in c('day','night')) {
    cluster_list[[i]][[tod]] <- cluster_list[[i]][[tod]] |> 
      lapply(add_environ, tod)
  }
}

# |- Model Data Formatting -------------

# The depth-bin specific approach requires that there are separate models for each
# time of day and each cluster

# Here, they must be formatting in a big list

detection_matrix_list <- vector('list', 4)
vol_matrix_list <- vector('list', 4)
db_list <- vector('list', 4)
env_data_list <- vector('list', 4)

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

get_site_covars <- function(tod, i) {
  x <- cluster_list[[i]][[tod]] |> 
    lapply(`[[`, 'env_data') |> 
    do.call(what = rbind,)
  return(x)
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
    env_data_list[[i]][[tod]] <- get_site_covars(tod,i)
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
    dcm = standardize(env_data_list[[i]][[tod]]$dcm_d[db_code$idx]),
    prey = standardize(env_data_list[[i]][[tod]]$intg_part[db_code$idx]),
    par = standardize(env_data_list[[i]][[tod]]$par[db_code$idx]),
    dac = standardize(env_data_list[[i]][[tod]]$dac[db_code$idx]),
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
param_out <- c('p0','a0','a_vs',
               'b0','b_dcm','b_prey','b_par','b_dac')

run_model <- function(mod_data) {
  mod <- jags(mod_data,
              parameters.to.save = param_out,
              model.file = './R/jags_04_occ-model.jags',
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

saveRDS(db_list, './data/06_dblist.rds')
saveRDS(epi_mod_output, './data/06_occ-model-output_epi_large.rds')
saveRDS(tmeso_mod_output, './data/06_occ-model-output_tmeso_large.rds')
saveRDS(bmeso_mod_output, './data/06_occ-model-output_bmeso_large.rds')
saveRDS(environ_summary, './data/06_enviornmental-summary.rds')
