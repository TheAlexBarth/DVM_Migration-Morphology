###
# Organize casts from ecopart project #######
###

# Requires desktop connection

library(EcotaxaTools)
library(batsFtpReadr)
library(suncalc)
library(lubridate)
# read in data
# trim extra par files
raw_data <- ecopart_import("~/BATS_data/export_all/export_raw_2023", trim_to_zoo = T)

# |- Meta data management -----------------------

# need to edit some stationId by hand
corrected_meta <- read.csv('~/BATS_data/export_all/cast_names_raw.csv') # need to update
raw_data$meta$stationid <- corrected_meta$stationid
raw_data$meta$ctd_origfilename <- corrected_meta$ctdref
raw_data$meta$programid <- corrected_meta$proj_id


# apply tod calcs
# There's sometimes I wrote 6 instead of 3 on the latitude
for(i in 1:length(raw_data$meta$latitude)) {
  if(raw_data$meta$latitude[i] > 60) {
    raw_data$meta$latitude[i] <- raw_data$meta$latitude[i] - 30
  }
}

# |- Assign Time of Day -----
raw_data$meta$tod <- rep(NA, nrow(raw_data$meta))
for(i in 1:length(raw_data$meta$tod)) {
  casttime <- raw_data$meta$sampledate[i]
  suntimes <- getSunlightTimes(
    date = as.Date(casttime),
    lat = raw_data$meta$latitude[i],
    lon = raw_data$meta$longitude[i],
    tz = 'UTC'
  )
  
  if(casttime < suntimes$nauticalDawn | casttime > suntimes$nauticalDusk) {
    raw_data$meta$tod[i] <- 'night'
  } else if (casttime > suntimes$nauticalDawn & casttime < suntimes$nauticalDusk) {
    raw_data$meta$tod[i] <- 'day'
  } else {
    raw_data$meta$tod[i] <- 'twilight'
  }
}

# |- Organizing and removing HS casts -------------------------

# Removing all hs casts then also
hs_index <- which(raw_data$meta$stationid %in% c('HS','hs'))

# need to sort file lists
raw_data$par_files <- raw_data$par_files[order(names(raw_data$par_files), 
                                               raw_data$meta$profileid)]

raw_data$zoo_files <- raw_data$zoo_files[order(names(raw_data$zoo_files),
                                               raw_data$meta$profileid)]

#slim down the ecopart_obj
raw_data$par_files <- raw_data$par_files[-hs_index]
raw_data$zoo_files <- raw_data$zoo_files[-hs_index]
raw_data$meta <- raw_data$meta[-hs_index,]

#remove the 389 project which is not available in ctd data
raw_data_trim <- list()
raw_data_trim$par_files <- raw_data$par_files[which(raw_data$meta$programid != 389)]
raw_data_trim$zoo_files <- raw_data$zoo_files[which(raw_data$meta$programid != 389)]
raw_data_trim$meta <- raw_data$meta[which(raw_data$meta$programid != 389),]

raw_data_trim <- as_ecopart_obj(raw_data_trim)

# |- Getting the CTD data --------------------------

# first for just the bats data
bats_proj <- raw_data_trim$meta$programid[which(raw_data_trim$meta$stationid == 'gf')]

bats_ctd <- read_bats_ctd(unique(bats_proj), 'gf')

# and the bloom data
bloom_proj <- raw_data_trim$meta$programid[which(raw_data_trim$meta$stationid == 'bloom')]
bloom_ctd <- read_bats_ctd(unique(bloom_proj), 'bloom')
bats_ctd$gf379b <- bloom_ctd

#split out based on cast
for(i in 1:length(bats_ctd)){
  bats_ctd[[i]] <- split(bats_ctd[[i]], f = bats_ctd[[i]]$ctd_id)
}
bats_ctd <- unlist(bats_ctd, recursive = F)
names(bats_ctd) <- names(bats_ctd) |> 
  strsplit('.',fixed = T) |> 
  sapply(`[[`,2)

# |- More Meta Management ------------------------------


# confirm the ctd-profile matching
verify <- data.frame(profileid = raw_data_trim$meta$profileid,
                      ctd_origfilename = raw_data_trim$meta$ctd_origfilename,
                      sampledate = raw_data_trim$meta$sampledate,
                      matcheddate = rep(NA, nrow(raw_data_trim$meta)),
                      diff = rep(NA, nrow(raw_data_trim$meta)),
                      flag = rep(F, nrow(raw_data_trim$meta)))

matcheddate <- vector(length = nrow(raw_data_trim$meta))
class(matcheddate) <- 'POSIXct'
matcheddate <- as_datetime(matcheddate)
for(i in 1:nrow(verify)) {
  matcheddate[i] <- bats_ctd[[as.character(verify$ctd_origfilename[i])]]$Date[1]
  verify$diff[i] <- difftime(matcheddate[i], verify$sampledate[i], units = 'mins')
  if(abs(verify$diff[i]) > 10){
    verify$flag[i] <- T
  }
}
verify$matcheddate <- matcheddate


# trim to only have the ctd casts that we need

#this is a bit of ugly code that gets the job done:
bats_ctd <- bats_ctd[which(names(bats_ctd) %in% raw_data_trim$meta$ctd_origfilename)]

# Format for script source
ctd_data <- bats_ctd
uvp_data <- raw_data_trim


# |- Save all the data --------------------------
# 
# saveRDS(bats_ctd, './Data/00_ctd-filtered.rds', compress = TRUE)
# saveRDS(raw_data_trim, file = gzfile('./Data/00_raw-trim-casts_large.rds.gz'), compress = TRUE)
