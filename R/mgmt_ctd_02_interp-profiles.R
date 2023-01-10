###
# Plot CTD data ######
###

rm(list = ls())
library(EcotaxaTools)
library(MBA)
library(ggplot2)
library(lubridate)
library(tidyr)
library(suncalc)
library(ggpubr)

###
# Load in data #########
###


# |- ctd data format -------------------------------
ctd_data <- readRDS('./data/01_ctd-trim-final.rds')

#mush all casts together
all_casts <- ctd_data|> list_to_tib('ctd_orgfilename')
all_casts <- all_casts[,-1]

#trim mess from the data
all_casts <- all_casts[-which(is.na(all_casts$Depth)),]
all_casts <- all_casts[-which(is.na(all_casts$Sal)),]
all_casts <- all_casts[-which(is.na(all_casts$DO)),]
all_casts <- all_casts[-which(is.na(all_casts$RFU)),]

#keep only shallow casts
all_casts <- all_casts[which(all_casts$Depth < 1200),]

#set cruise name
all_casts$cruise <- paste0(all_casts$proj_id, all_casts$cruise_id)

####
# Surface interpolation #############################################
####

# I'm using b-spline interpolation to try and get pretty figures

# |- worker fx ---------------------

surf_interp <- function(df,
                        value,
                        res = 300,
                        ...) {
  xyz <- cbind(decimal_date(df$Date),
               df$Depth,
               df[[value]])
  interp_output <- mba.surf(xyz, no.X = res, no.Y = res, ...)
  long_interp <- interp_output$xyz.est$z |> 
    as.data.frame() |> 
    pivot_longer(cols = everything())
  
  outdf <- data.frame(Date = date_decimal(interp_output$xyz.est$x),
                      Depth = interp_output$xyz.est$y) |> 
    expand.grid()
  outdf[[value]] <- c(interp_output$xyz.est$z)
  return(outdf)
}

# set up
value_class <- function(cruise){
  ret_list <- vector('list',4)
  names(ret_list) <- c('Sal',"Temperature",'RFU','DO')
  return(ret_list)
}


# |- Loop for individual cruises -------------

# shell storage
cruises <- unique(all_casts$cruise)
plot_values <- lapply(cruises, value_class)
names(plot_values) <- cruises

# loop it out
for(i in 1:length(plot_values)) {
  for(j in 1:length(plot_values[[i]])) {
    df <- all_casts[all_casts$cruise == names(plot_values)[i],]
    value <- names(plot_values[[i]])[j]
    plot_values[[i]][[j]] <- surf_interp(df = df, value = value)
  }
}



####
#  Creating All Cruise profiles ############
####

# create sub-set frames
early_cruises <- all_casts[year(all_casts$Date) < 2020,]
later_cruises <- all_casts[year(all_casts$Date) >= 2020,]

full_proj_values <- list()
full_proj_values$early <- value_class()
full_proj_values$late <- value_class()

for(i in 1:length(full_proj_values)) {
  if (i < 2) {
    df <- early_cruises
  } else {
    df <- later_cruises
  }
  
  for(j in 1:length(full_proj_values[[i]])) {
    value <- names(full_proj_values[[i]])[j]
    full_proj_values[[i]][[j]] <- surf_interp(df = df, value = value)
  }
}

###
# Save Data #####
###
saveRDS(list(by_cruise = plot_values,
          full_proj = full_proj_values),
     './data/ctd_02_interp-grid.rds')

