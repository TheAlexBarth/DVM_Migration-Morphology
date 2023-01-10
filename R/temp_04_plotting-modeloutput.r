rm(list = ls())
library(wiqid)
library(ggpubr)
library(ggplot2)
library(EcotaxaTools)
source('./R/tools.R')


full_col_mod_epi <- readRDS('./data/04_full-col-model-euphotic.rds')
full_col_mod_tmeso <- readRDS('./data/04_full-col-model-top-meso.rds')
full_col_mod_bmeso <- readRDS('./data/04_full-col-model-bot-meso.rds')
count_data <- readRDS('./data/04_count-data-25mBins.rds')
####
# Detection Probability ########
####

# |- Volume Sampled Slope ----------

# vs range 2.2 to 1288.1
vs_range <- seq(0,1300, length = 100)
vs_rangeS <- standardize2match(vs_range, count_data$vol_sampled)

vol_range_sims <- function(tod, model) {
  

  out_data <- data.frame(
    mean = rep(NA, length(vs_range)),
    low = rep(NA, length(vs_range)),
    high = rep(NA, length(vs_range)),
    vs = vs_range
  )

  
  for(i in 1:length(vs_range)) {
    temp_vs <- with(model[[tod]]$sims.list,
                    a0 + a_vs * vs_rangeS[i])
    p_temp <- plogis(temp_vs)
    out_data[i, 1:3] <- c(mean(p_temp), hdi(p_temp))
  }
  
  return(out_data)
}

day_pred_epi <- vol_range_sims('day', full_col_mod_epi)
night_pred_epi <- vol_range_sims('night', full_col_mod_epi)

# |- Epipelagic Plot --------------------------------------
epi_vs_plot <- ggplot() + 
  geom_line(data = day_pred_epi,
            aes(x = vs, y = mean,
                col = 'day')) +
  geom_ribbon(data = day_pred_epi,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'day'), alpha = .25) + 
  geom_line(data = night_pred_epi,
            aes(x = vs, y = mean, col = 'night')) +
  geom_ribbon(data = night_pred_epi,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'night'), alpha = .25)+
  scale_fill_manual(values = dn_cols[c('day','night')]) + 
  scale_color_manual(values = dn_cols[c('day','night')]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Detection Probability', x = 'Volume Sampled', fill = "") +
  guides(color = 'none')+
  theme_pubr()
  




# |- tmesopelagic Plot --------------------------------------
day_pred_tmeso <- vol_range_sims('day', full_col_mod_tmeso)
night_pred_tmeso <- vol_range_sims('night', full_col_mod_tmeso)


tmeso_vs_plot <- ggplot() + 
  geom_line(data = day_pred_tmeso,
            aes(x = vs, y = mean,
                col = 'day')) +
  geom_ribbon(data = day_pred_tmeso,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'day'), alpha = .25) + 
  geom_line(data = night_pred_tmeso,
            aes(x = vs, y = mean, col = 'night')) +
  geom_ribbon(data = night_pred_tmeso,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'night'), alpha = .25)+
  scale_fill_manual(values = dn_cols[c('day','night')]) + 
  scale_color_manual(values = dn_cols[c('day','night')]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Detection Probability', x = 'Volume Sampled', fill = "") +
  guides(color = 'none')+
  theme_pubr()


# |- Bot Mesopelagic Plot --------------------------------------------

day_pred_bmeso <- vol_range_sims('day', full_col_mod_bmeso)
night_pred_bmeso <- vol_range_sims('night', full_col_mod_bmeso)


bmeso_vs_plot <- ggplot() + 
  geom_line(data = day_pred_bmeso,
            aes(x = vs, y = mean,
                col = 'day')) +
  geom_ribbon(data = day_pred_bmeso,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'day'), alpha = .25) + 
  geom_line(data = night_pred_bmeso,
            aes(x = vs, y = mean, col = 'night')) +
  geom_ribbon(data = night_pred_bmeso,
              aes(x = vs, ymin = low, ymax = high,
                  fill = 'night'), alpha = .25)+
  scale_fill_manual(values = dn_cols[c('day','night')]) + 
  scale_color_manual(values = dn_cols[c('day','night')]) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = 'Detection Probability', x = 'Volume Sampled', fill = "") +
  guides(color = 'none')+
  theme_pubr()

ggarrange(epi_vs_plot, tmeso_vs_plot, bmeso_vs_plot,
          nrow = 1, common.legend = T)

###
# Occurrence Probability #####
###

# |- Model Extraction ------------------

psi_extract <- function(model, tod, min_d, max_d) {
  db_length <- ncol(model[[tod]]$sims.list$psi)
  
  out_df <- data.frame(
    db = db_decode(c(1:db_length), tod, min_d, max_d, db = model$db),
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



day_psi <- do.call(rbind, list(
  psi_extract(full_col_mod_epi, 'day', 0, 200),
  psi_extract(full_col_mod_tmeso, 'day', 200, 600),
  psi_extract(full_col_mod_bmeso, 'day', 600, 1200)
)) |> 
  bin_format()


night_psi <- do.call(rbind, list(
  psi_extract(full_col_mod_epi, 'night', 0, 200),
  psi_extract(full_col_mod_tmeso, 'night', 200, 600),
  psi_extract(full_col_mod_bmeso, 'night', 600, 1200)
)) |> 
  bin_format()

day_psi$tod <- 'day'
night_psi$tod <- 'night'

psi_df <- rbind(day_psi, night_psi)

# |- Occurrence plot ---------------------

ggplot(psi_df) + 
  geom_pointrange(aes(y = mp,
                      x = mean,
                      xmin = low,
                      xmax = high,
                      color = tod)) +
  facet_grid(~tod) +
  scale_y_reverse() +
  scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  theme_pubr() +
  theme(legend.position = 'none')

