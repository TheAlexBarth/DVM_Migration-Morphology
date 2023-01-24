####
# Filter uvp and ctd geographically ####
####





# |- get data in data -----------
source('./R/mgmt_00_profile-org.R')
rm(list = ls()[-which(ls() %in% c('ctd_data','uvp_data'))])


# |- Get ctd lat lons ---------------
ctd_loc <- data.frame(ctdid = names(ctd_data),
                      lat = rep(NA, length(ctd_data)),
                      lon = rep(NA, length(ctd_data)),
                      tod = rep(NA, length(ctd_data)),
                      progid = rep(NA,length(ctd_data)))

for(i in 1:nrow(ctd_loc)){
  ctd_loc$lat[i] <- ctd_data[[i]]$Lat[1]
  ctd_loc$lon[i] <- ctd_data[[i]]$Lon[1]
  ctd_loc$tod[i] <- uvp_data$meta$tod[which(uvp_data$meta$ctd_origfilename == ctd_loc$ctdid[i])]
  ctd_loc$progid[i] <- uvp_data$meta$programid[which(uvp_data$meta$ctd_origfilename == ctd_loc$ctdid[i])]
}


# # |- All cast map ---------------------
# all_cast_map <- ggplot()+
#   geom_sf(data = world, fill = 'green')+
#   geom_point(aes(x = -ctd_loc$lon, y = ctd_loc$lat,
#                  shape = ctd_loc$tod, col = ctd_loc$tod),
#              position = position_jitter(width = .05, height = 0.05))+
#   coord_sf(xlim = c(-65,-62), ylim = c(30,33))+
#   geom_segment(aes(x = -65.1, xend = -62.9,
#                    y = 31.8, yend = 32.75),
#                size = 2, col = 'red')+
#   theme_bw()+
#   labs(x = "",y = "", shape = "")+
#   scale_color_manual(values = c('beige','black'))+
#   theme(panel.background = element_rect(fill = 'cornflowerblue'),
#         axis.text.x = element_text(angle = 45, hjust = c(1,1)),
#         legend.position = 'bottom')


# |- Removing casts outside range ------------------------
drop_casts <- which(
  ctd_loc$lat < 30.55 |
    ctd_loc$lat > 32.1 & ctd_loc$lon > 64.25
)

trim_ctd <- ctd_data[which(names(ctd_data) %in% ctd_loc$ctdid[-drop_casts])]

# |- Create for ctd map ---------------------------------------
# only run once, not on sources
trim_loc <- ctd_loc[-drop_casts,]
saveRDS(trim_loc, './data/01_ctd-locations.rds')

# |- trim and save data ------------------------------------

trim_uvp <- list()
trim_uvp$meta <- uvp_data$meta[which(uvp_data$meta$ctd_origfilename %in% names(trim_ctd)),]
trim_uvp$par_files <- uvp_data$par_files[which(names(uvp_data$par_files) %in% trim_uvp$meta$profileid)]
trim_uvp$zoo_files <- uvp_data$zoo_files[which(names(uvp_data$zoo_files) %in% trim_uvp$meta$profileid)]

trim_uvp <- trim_uvp |> as_ecopart_obj()

saveRDS(trim_ctd, './Data/01_ctd-trim-final.rds')
saveRDS(trim_uvp, './Data/01_uvp-trim-final_large.rds')
