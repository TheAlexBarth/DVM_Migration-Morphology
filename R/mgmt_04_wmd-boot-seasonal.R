###
# Calcuating WMD bootstrapped for each group #####
###

# |- Set up ---------------------

rm(list = ls())
library(EcotaxaTools)
library(ggplot2)
library(ggpubr)
library(lubridate)
source('./R/tools.R')

cluster_cope <- readRDS('./data/03_cluster-conc.rds')
uvp_data  <- readRDS('./data/01_uvp-trim-final.rds') |> trim_to_cope()

cast_times <- list(day = uvp_data$meta$profileid[which(uvp_data$meta$tod == 'day')],
                   night = uvp_data$meta$profileid[which(uvp_data$meta$tod == 'night')])

# |- Assign cruise_id  & season ----------------- 

uvp_data$meta$season <- rep(NA, nrow(uvp_data$meta))
for(i in 1:nrow(uvp_data$meta)) {
  if(month(uvp_data$meta$sampledate[i]) %in% c(12,1,2)) {
    uvp_data$meta$season[i] <- 'winter'
  } else if (month(uvp_data$meta$sampledate[i]) %in% c(3,4,5)) {
    uvp_data$meta$season[i] <- 'spring'
  } else if (month(uvp_data$meta$sampledate[i]) %in% c(6,7,8)) {
    uvp_data$meta$season[i] <- 'summer'
  } else if (month(uvp_data$meta$sampledate[i]) %in% c(9,10,11)) {
    uvp_data$meta$season[i] <- 'fall'
  }
}


for(i in 1:length(cluster_cope)) {
  cluster_cope[[i]]['cruise_id'] <- paste0(uvp_data$meta$stationid[i], 
                                           uvp_data$meta$programid[i])
  
  cluster_cope[[i]]['season'] <- uvp_data$meta$season[i]
}


# |- Bin Contrained Bootstrapping --------------------------
# bcBoot will randomly select across all casts an observation from each depth bin
# then these observations are used to construct a profile which can then be used
# to calculate the wmd. First, what must be fed-in is a dataframe of concentrations
# from casts. This can be done with merge_cast for criteria matching

data_pool <- cluster_cope |> 
  lapply(bin_format) |> 
  merge_casts(cast_times) |> 
  list_to_tib('tod')

wmd_boots <- list(
  day = list(),
  night = list()
)

data_pool <- data_pool[data_pool$max_d <= 600,]

set.seed(20221110)
for(i in 1:length(unique(data_pool$group))) {
  
  wmd_boots$day[[i]] <- bc_boot_wmd(data = data_pool[data_pool$tod == 'day' &
                                                     data_pool$group == i,],
                                    iter = 999)
  
  wmd_boots$night[[i]] <- bc_boot_wmd(data = data_pool[data_pool$tod == 'night' &
                                                         data_pool$group == i,],
                                      iter = 999)
}



for(i in 1:4) {
  print(wmd_boots$day[[i]] |> quantile(c(0.025,0.975)))
  print(wmd_boots$night[[i]] |> quantile(c(0.025,0.975)))
  print('------------------')
}

###
# Set up season-specific ###################
###

# init shell
cruise_boot <- vector('list',length(unique(data_pool$season)))
names(cruise_boot) <- unique(data_pool$season)

fill_boot <- function(x) {
  x <- vector('list',4)
  names(x) <- c(paste0('c',1:4))
  for(i in 1:length(x)) {
    x[[i]][['wmd_d']] <- NA
    x[[i]][['wmd_n']] <- NA
    x[[i]][['dn_diff']] <- NA
  }
  return(x)
}
cruise_boot <- lapply(cruise_boot, fill_boot)

# |- Fill out bootstrapped data for each ----------------

start = Sys.time()
for(i in 1:length(cruise_boot)) {
  season <- names(cruise_boot)[i]
  for(j in 1:length(unique(data_pool$group))) {
    sub_data <- data_pool[data_pool$season == season &
                            data_pool$group == j,]
    
    if(length(unique(sub_data$tod)) <=1) {
      next
    }
    
    cruise_boot[[i]][[j]]$wmd_d <- bc_boot_wmd(data = sub_data[sub_data$tod == 'day',],
                                               iter = 999)
    cruise_boot[[i]][[j]]$wmd_n <- bc_boot_wmd(data = sub_data[sub_data$tod == 'night',],
                                               iter = 999)
    cruise_boot[[i]][[j]]$dn_diff <- cruise_boot[[i]][[j]]$wmd_d - cruise_boot[[i]][[j]]$wmd_n
  }
}
end = Sys.time()
end-start


# |- Save the data -----------------
saveRDS(wmd_boots, './data/04_wmd-all-merged.rds')
saveRDS(cruise_boot, './data/04_wmd-by-cruise-seasonal.rds')
