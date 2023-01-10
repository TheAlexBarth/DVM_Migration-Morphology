###
# Importing 490 data #######
###

rm(list = ls())

library(dplyr)
library(readr)

modis_490_day <- read_csv('~/BATS_data/AquaMODIS_490/modis_490_day.csv', skip = 1)
modis_490_day <- modis_490_day |> 
  group_by(UTC) |> 
  summarise(dwi_490 = mean(`m-1`, na.rm = T))

modis_par_day <- read.csv('~/BATS_data/AquaMODIS_490/modis_490_day.csv', skip = 1) |> 
  group_by(UTC) |> 
  summarize(par = mean(m.1, na.rm = T))

modis_490_day$UTC <- as.Date(modis_490_day$UTC)
modis_490_day$dac <- 1/modis_490_day$dwi_490 #calculate the diffuse attenuation coef
modis_par_day$UTC <- as.Date(modis_par_day$UTC)

saveRDS(list(
  par = modis_par_day,
  dac = modis_490_day
),
"./Data/00_modis-data.rds", compress = TRUE)
