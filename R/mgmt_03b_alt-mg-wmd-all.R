####
# Alternative MG wmd analysis #######
####


rm(list = ls())
library(EcotaxaTools)
library(ggplot2)
library(ggpubr)
source('./R/tools.R')

pc1_bins <- readRDS('./data/03b_pc1-bins.rds')
uvp_data  <- readRDS('./data/01_uvp-trim-final_large.rds') |> trim_to_cope()

cast_times <- list(day = uvp_data$meta$profileid[which(uvp_data$meta$tod == 'day')],
                   night = uvp_data$meta$profileid[which(uvp_data$meta$tod == 'night')])



# |- Assign cruise_id ----------------- 

assign_cruise_id <- function(bin_group) {
  for(profileid in names(bin_group)) {
    uvp_meta_subset <- uvp_data$meta[uvp_data$meta$profileid == profileid,]
    bin_group[[profileid]]['cruise_id'] <- paste0(uvp_meta_subset$stationid,
                                                  uvp_meta_subset$programid)
  }
  return(bin_group)
}

pc1_bins <- pc1_bins |> assign_cruise_id()

# |- Bin Contrained Bootstrapping --------------------------
# bcBoot will randomly select across all casts an observation from each depth bin
# then these observations are used to construct a profile which can then be used
# to calculate the wmd. First, what must be fed-in is a dataframe of concentrations
# from casts. This can be done with merge_cast for criteria matching


bin_booter <- function(bin_group) {
  data_pool <- bin_group |> 
    lapply(bin_format) |> 
    merge_casts(cast_times) |> 
    list_to_tib('tod')
  
  wmd_boots <- list(
    day = list(),
    night = list()
  )
  
  data_pool <- data_pool[data_pool$max_d <= 600,]
  
  for(i in unique(data_pool$group)) {
    wmd_boots$day[[i]] <- bc_boot_wmd(data = data_pool[data_pool$tod == 'day' &
                                                         data_pool$group == i,],
                                      iter = 999)
    
    wmd_boots$night[[i]] <- bc_boot_wmd(data = data_pool[data_pool$tod == 'night' &
                                                           data_pool$group == i,],
                                        iter = 999)
  }
  
  for(i in unique(data_pool$group)) {
    print(i)
    print(wmd_boots$day[[i]] |> quantile(c(0.025,0.975)))
    print(wmd_boots$night[[i]] |> quantile(c(0.025,0.975)))
    print('------------------')
  }
  
  return(wmd_boots)
}

set.seed(20230314)
pc1_wmd_boots <- bin_booter(pc1_bins)



# Same run for PC2 #######################

# |- Run for small PC2 ----------------------------------
pc2_bin_sm <- readRDS('./data/03b_pc2-sm-bins.rds')

pc2_sm_wmd_boots <- pc2_bin_sm |> 
  assign_cruise_id() |> 
  bin_booter()


# |- Run for mid PC2 ---------------------------------
pc2_bin_mid <- readRDS('./data/03b_pc2-bins.rds')

pc2_mid_wmd_boots <- pc2_bin_mid |> 
  assign_cruise_id() |> 
  bin_booter()

# |- Run for lg pc2 ----------------------------------
pc2_bin_lg <- readRDS('./data/03b_pc2-lg-bins.rds')

pc2_lg_wmd_boots <- pc2_bin_lg |>
  assign_cruise_id() |> 
  bin_booter()


# Save All Data ######################
saveRDS(list(
  pc1 = pc1_wmd_boots,
  pc2_sm = pc2_sm_wmd_boots,
  pc2_mid = pc2_mid_wmd_boots,
  pc2_lg = pc2_lg_wmd_boots
),'./data/03b_alt-mg-wmd-bins.rds')
