######
# Importing k490 downwelling data #######
######

rm(list = ls())
library(readr)
library(dplyr)
library(EcotaxaTools)

par_df <- read_csv('~/BATS_data/AquaMODIS_490/modis_par-raw.csv', skip = 1)
dwi_df <- read_csv('~/BATS_data/AquaMODIS_490/modis_490-raw.csv', skip = 1)

dwi_df$dwi <- 1/dwi_df$`m-1`


# |- Cleaning and formatting data ------------------
par_df <- par_df[-which(is.na(par_df$`einstein m-2 day-1`)),]

par_sum <- par_df |> 
  group_by(UTC) |> 
  summarise(mean(`einstein m-2 day-1`))
names(par_sum) <- c('Date','par')

dwi_sum <- dwi_df |> 
  group_by(UTC) |> 
  summarise(mean(dwi, na.rm = T))
names(dwi_sum) <- c('Date','dwi')

modis_mosum <- left_join(par_sum, dwi_sum)

# |- Save Data -----------------------------
saveRDS(modis_mosum, './data/05_modis-month-sum.rds')


#####
# Calculating Chlorophll metrics ##########
#####

# Interested in the DCM location
rm(list = ls())
library(EcotaxaTools)


# |- ctd data format -------------------------------
ctd_data <- readRDS('./data/01_ctd-trim-final.rds')
uvp_data <- readRDS('./data/01_uvp-trim-final.rds')

# |- getting DCM Location -------------

dcm_max <- function(df) {
  df$Depth[which.max(df$RFU)]
}

dcm_depth <- tibble(
  ctd_origfilename = as.integer(names(ctd_data)),
  dcm_d = ctd_data |> sapply(dcm_max)
) |> 
  left_join(uvp_data$meta[,c("ctd_origfilename", 'profileid', 'sampledate')])

# |- getting integrated particle abundance ----------------------------

#need to fix particle concentration
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
saveRDS(cast_metrics,'./data/05_cast-metrics.rds')
