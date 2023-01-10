###
# Plotting output for environmental output ####
###

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(wiqid)
source('./R/tools.R')

db_list <- readRDS('./data/06_dblist.rds')
epi_mod <- readRDS('./data/06_occ-model-output_epi.rds')
tmeso_mod <- readRDS('./data/06_occ-model-output_tmeso.rds')
bmeso_mod <- readRDS('./data/06_occ-model-output_bmeso.rds')
environ_summary <- readRDS('./data/06_enviornmental-summary.rds')


# |- Extracting the parameters ------------------------

# to make the slope interpretable as the effect of 1 unit increase of that variable
# I need to match a value of 1 to the scaled value of the variable

conv_to_plogis <- function(x, var) {
  1/(1+exp(-x * standardize2match(1,var)))
}



# function to extract the features
feature_extract <- function(model, tod, min_d, max_d, db) {
  
  b_dcm_hdi <- model[[tod]]$sims.list$b_dcm |> 
    apply(2,hdi)
  b_dcm_mean <- model[[tod]]$sims.list$b_dcm |> 
    apply(2,mean)
  b_prey_hdi <- model[[tod]]$sims.list$b_prey |> 
    hdi()
  b_prey_mean <- model[[tod]]$sims.list$b_prey |> 
    apply(2,mean)
  b_par_hdi <- model[[tod]]$sims.list$b_par |> 
    hdi()
  b_par_mean <- model[[tod]]$sims.list$b_par |> 
    apply(2,mean)
  b_dac_hdi <- model[[tod]]$sims.list$b_dac |> 
    hdi()
  b_dac_mean <- model[[tod]]$sims.list$b_dac |> 
    apply(2,mean)
  
  db_out <- db_decode(1:length(b_dac_mean), tod, min_d, max_d, db = db)
  
  out_df <- data.frame(
    db = db_out,
    dcm_low = b_dcm_hdi[1,],
    dcm_high = b_dcm_hdi[2,],
    dcm = b_dcm_mean,
    prey_low = b_prey_hdi[1,],
    prey_high = b_prey_hdi[2,],
    prey = b_prey_mean,
    dac_low = b_dac_hdi[1,],
    dac_high = b_dac_hdi[2,],
    dac = b_dac_mean,
    par_low = b_par_hdi[1,],
    par_high = b_par_hdi[2,],
    par = b_par_mean
  )
  
  return(out_df |> bin_format())
}

epi_list <- vector('list',4)
tmeso_list <- vector('list',4)
bmeso_list <- vector('list',4)
for(i in 1:4) {
  for(tod in c('day','night')) {
    epi_list[[i]][[tod]] <- feature_extract(epi_mod[[i]], tod, 0, 200, db_list[[i]])
    tmeso_list[[i]][[tod]] <- feature_extract(tmeso_mod[[i]], tod, 200, 600, db_list[[i]])
    bmeso_list[[i]][[tod]] <- feature_extract(bmeso_mod[[i]], tod, 600, 1200, db_list[[i]])
  }
}

# mash into one big data frame 

names(epi_list) <- 1:4
epi_df <- epi_list |> 
  lapply(list_to_tib,'tod') |> 
  list_to_tib('cluster')
names(tmeso_list) <- 1:4
tmeso_df <- tmeso_list |> 
  lapply(list_to_tib,'tod') |> 
  list_to_tib('cluster')
names(bmeso_list) <- 1:4
bmeso_df <- bmeso_list |> 
  lapply(list_to_tib,'tod') |> 
  list_to_tib('cluster')

feature_df <- do.call(rbind, list(epi_df, tmeso_df, bmeso_df))


####
# Plotting Effects #######
####

dcm_plot <- ggplot(feature_df) + 
  geom_ribbon(aes(x = mp,
                  ymin = dcm_low,
                  ymax = dcm_high,
                  fill = tod),
              alpha = 0.5) +
  geom_line(aes(x = mp,
                  y = dcm,
                  color = tod)) +
  geom_hline(yintercept = 0) +
  facet_grid(~cluster+tod) +
  scale_x_reverse() +
  scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  coord_flip()+
  theme_pubr() +
  labs(x = 'Depth [m]', y = "Effect on LogOdds(occupancy)") +
  theme(legend.position = 'none')

prey_plot <- ggplot(feature_df) + 
  geom_ribbon(aes(x = mp,
                  ymin = prey_low,
                  ymax = prey_high,
                  fill = tod),
              alpha = 0.5) +
  geom_line(aes(x = mp,
                y = prey,
                color = tod)) +
  geom_hline(yintercept = 0) +
  facet_grid(~cluster+tod) +
  scale_x_reverse() +
  scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  coord_flip()+
  theme_pubr() +
  labs(x = 'Depth [m]', y = "Effect on LogOdds(occupancy)") +
  theme(legend.position = 'none')

par_plot <- ggplot(feature_df) + 
  geom_ribbon(aes(x = mp,
                  ymin = par_low,
                  ymax = par_high,
                  fill = tod),
              alpha = 0.5) +
  geom_line(aes(x = mp,
                y = par,
                color = tod)) +
  geom_hline(yintercept = 0) +
  facet_grid(~cluster+tod) +
  scale_x_reverse() +
  scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  coord_flip()+
  theme_pubr() +
  labs(x = 'Depth [m]', y = "Effect on LogOdds(occupancy)") +
  theme(legend.position = 'none')

dac_plot <- ggplot(feature_df) + 
  geom_ribbon(aes(x = mp,
                  ymin = dac_low,
                  ymax = dac_high,
                  fill = tod),
              alpha = 0.5) +
  geom_line(aes(x = mp,
                y = dac,
                color = tod)) +
  geom_hline(yintercept = 0) +
  facet_grid(~cluster+tod) +
  scale_x_reverse() +
  scale_fill_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  scale_color_manual(values = c(dn_cols['day'],dn_cols['night'])) +
  coord_flip()+
  theme_pubr() +
  labs(x = 'Depth [m]', y = "Effect on LogOdds(occupancy)") +
  theme(legend.position = 'none')

ggarrange(dcm_plot + theme(axis.title.x = element_blank(),
                           axis.line.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_blank()),
          prey_plot+ theme(axis.title.x = element_blank(),
                           axis.line.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_blank(),
                           strip.background = element_blank(),
                           strip.text = element_blank()),
          par_plot+ theme(axis.title.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_blank(),
                          strip.background = element_blank(),
                          strip.text = element_blank()),
          dac_plot+ theme(strip.background = element_blank(),
                          strip.text = element_blank()),
          nrow = 4)



