###
# Plotting Occupancy of different clusters ##########
###

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(wiqid)
source('./R/tools.R')

mod_output <- readRDS('./data/05_cluster-occ-model-output.rds')

# |- Extact credible intervals from models -------------

psi_extract <- function(model, tod, min_d, max_d, db = db) {
  db_length <- ncol(model[[tod]]$sims.list$psi)
  
  out_df <- data.frame(
    db = db_decode(c(1:db_length), tod, min_d, max_d, db = db),
    mean = rep(NA, db_length),
    low = rep(NA, db_length),
    high = rep(NA, db_length)
  )
  
  for(i in 1:db_length) {
    out_df[i,-1] <- c(mean(model[[tod]]$sims.list$psi[,i]),
                      hdi(model[[tod]]$sims.list$psi[,i]))
  }
  return(out_df)
}


occ_pred <- vector('list', 4)
for(i in 1:4) {
  for(tod in c('day', 'night')) {
    occ_pred[[i]][[tod]] <- do.call(rbind,list(
      psi_extract(mod_output$epi[[i]],
                  tod, 0, 200,
                  db = mod_output$db[[i]]),
      psi_extract(mod_output$tmeso[[i]],
                  tod, 200, 600,
                  db = mod_output$db[[i]]),
      psi_extract(mod_output$bmeso[[i]],
                  tod, 600, 1200,
                  db = mod_output$db[[i]]))) |> 
      bin_format()
  }
}

names(occ_pred) <- c(1:4)
occ_pred <- occ_pred |> 
  lapply(list_to_tib, 'tod') |> 
  list_to_tib('cluster')


###
# Plotting ######
###
ribbon_plot <- function(cluster) {
  plot_data <- occ_pred[occ_pred$cluster == cluster,]
  
  out_plot <- ggplot(plot_data) +
    geom_ribbon(aes(x = mp,
                        y = mean,
                        ymin = low,
                        ymax = high,
                        fill = tod,
                    color = tod),
                alpha = 0.5) +
    geom_line(aes(x = mp,
                  y = mean,
                  color = tod),
              size = 1,
              alpha = .75) +
    facet_grid(~tod)+
    scale_x_reverse() +
    scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
    scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
    coord_flip()+
    theme_pubr() +
    labs(x = 'Depth [m]', y = "Occupancy") +
    theme(legend.position = 'none')
  return(out_plot)
}

point_plot <- function(cluster) {
  plot_data <- occ_pred[occ_pred$cluster == cluster,]
  
  out_plot <- ggplot(plot_data) +
    geom_pointrange(aes(x = mp,
                        y = mean,
                        ymin = low,
                        ymax = high,
                        color = tod),
                    position = position_dodge(width = 5)) +
    scale_x_reverse() +
    scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
    scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
    coord_flip()+
    theme_pubr() +
    labs(x = 'Depth [m]', y = "Occupancy") +
    theme(legend.position = 'none')
  return(out_plot)
}

ribbon_occupancy <- vector('list',4)
point_occupancy <- vector('list', 4)
for(i in 1:4) {
  ribbon_occupancy[[i]] <- ribbon_plot(i)
  point_occupancy[[i]] <- point_plot(i)
}

ggarrange(ribbon_occupancy[[1]],
          ribbon_occupancy[[2]] + theme(axis.title.y = element_blank()),
          ribbon_occupancy[[3]] + theme(axis.title.y = element_blank()),
          ribbon_occupancy[[4]] + theme(axis.title.y = element_blank()),
          ncol = 4)

ggarrange(point_occupancy[[1]],
          point_occupancy[[2]] + theme(axis.title.y = element_blank()),
          point_occupancy[[3]] + theme(axis.title.y = element_blank()),
          point_occupancy[[4]] + theme(axis.title.y = element_blank()),
          ncol = 4)
