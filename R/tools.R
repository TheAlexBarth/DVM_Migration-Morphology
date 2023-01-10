###
# Tools #######
###

# These are code to run functions / actions which are used across multiple
# Scripts but I want to avoid saving too many large data files

#' trim_to_cope
#' 
#' Trim the uvp_data object to only include copepods

trim_to_cope <- function(uvp_data) {
  uvp_data <- uvp_data |> 
    add_zoo(names_to, col_name = 'name', 
            new_names = c('living','not-living','Copepoda'), 
            suppress_print = T) |> 
    mod_zoo(names_keep, 'Copepoda')
  
  # Find any casts with no-copepods
  drop_casts <- uvp_data$zoo_files |>
    sapply(nrow) |> 
    sapply(function(x) x == 0) |> 
    which() |> 
    names()
  
  if(length(drop_casts) > 0) {
    uvp_data$par_files <- uvp_data$par_files[which(names(uvp_data$par_files) != drop_casts)]
    uvp_data$zoo_files <- uvp_data$zoo_files[which(names(uvp_data$zoo_files) != drop_casts)]
    uvp_data$meta <- uvp_data$meta[which(uvp_data$meta$profileid != drop_casts),]
  }
  # restore class structure
  uvp_data <- as_ecopart_obj(uvp_data)
  
  # trim UVP to only have above 1200m
  uvp_data <- uvp_data |> trim_ecopart_depth(1200)
  
  return(uvp_data)
}


trim_to_crust <- function(uvp_data) {
  uvp_data <- uvp_data |> 
    add_zoo(names_to, col_name = 'name', 
            new_names = c('living','not-living','Crustacea'), 
            suppress_print = T) |> 
    mod_zoo(names_keep, 'Crustacea')
  
  # Find any casts with no-copepods
  drop_casts <- uvp_data$zoo_files |>
    sapply(nrow) |> 
    sapply(function(x) x == 0) |> 
    which() |> 
    names()
  
  if(length(drop_casts) > 0) {
    uvp_data$par_files <- uvp_data$par_files[which(names(uvp_data$par_files) != drop_casts)]
    uvp_data$zoo_files <- uvp_data$zoo_files[which(names(uvp_data$zoo_files) != drop_casts)]
    uvp_data$meta <- uvp_data$meta[which(uvp_data$meta$profileid != drop_casts),]
  }
  # restore class structure
  uvp_data <- as_ecopart_obj(uvp_data)
  
  # trim UVP to only have above 1200m
  uvp_data <- uvp_data |> trim_ecopart_depth(1200)
  
  return(uvp_data)
}

#' Calculate weighted mean depth
#' 
#' based on mid-points of depth bins and observed concentrations of a profile
wmd <- function(conc_vect, mp_vect) {
  total_conc <- sum(conc_vect)
  wd <- (conc_vect * mp_vect) / total_conc
  return(sum(wd))
}


#' Bin-Constrained bootstrapping
#' 
#' 
bc_boot_wmd <- function(data, iter) {
  
  out_wmd <- rep(NA, iter)
  for(i in 1:iter){
    
    boot_data <- data |> 
      split(f = data$mp) |> 
      lapply(bc_sample) |> 
      list_to_tib('mp')
    
    out_wmd[i] <- wmd(boot_data$conc_m3, as.numeric(boot_data$mp))
  
  }
  return(out_wmd)
}

#' Helper function for bc_boot
#' 
#' sample the data field at a given mp
bc_sample <- function(sub_data) {
  
  boot_conc <- sample(sub_data$conc_m3,replace = T)
  return(data.frame(conc_m3 = boot_conc,
                    mp = sub_data$mp))
  
}

#' day night colors
dn_cols <- c(
  `day` = '#F9C687',
  `midday` = '#F0983F',
  `darknight` = '#2A394C',
  `night` = '#436871'
)

#' messy conversion function to keep it clean later
quant_to_df <- function(quant) {
  rdf <- data.frame(
    low = rep(NA,1),
    mid = rep(NA,1),
    high = rep(NA,1)
  )
  rdf$low <- quant[1]
  rdf$mid <- quant[2]
  rdf$high <- quant[3]
  return(rdf)
}

#' wmd_df_cruise
# honestly, at this point I've created a ridiculous data structure
# the list is too deep - it is a messy tree
# there has to be a better way but this is going to just be an ugly solution
# I have re-sampled distributions in a tree-like list. I then need to get the quantiles
# and force it all into a top-level dataframe
wmd_df_construct <- function(wmd_data, method) {
  
  wmd_list <- vector('list',length(wmd_data))
  names(wmd_list) <- names(wmd_data)
  for(i in 1:length(wmd_list)) {
    #get quantiles recursively
    wmd_list[[i]] <- wmd_data[[i]] |> 
      lapply(function(x) lapply(x, quantile, probs = c(0.025,.5,0.975), na.rm = T)) 
    
    
    wmd_list[[i]] <- wmd_list[[i]] |> 
      lapply(function(x) lapply(x, quant_to_df)) |> 
      lapply(list_to_tib, 'metric') |> 
      list_to_tib('cluster')
  }
  
  if(method == 'cruise'){
    wmd_df <- wmd_list |> list_to_tib('cruise_id')
    
    wmd_df$date <- rep(NA, nrow(wmd_df))
    for(i in 1:length(unique(wmd_df$cruise_id))) {
      cruise_temp <- unique(wmd_df$cruise_id)[i]
      meta_idx <- which(uvp_data_meta$cruise_id == unique(wmd_df$cruise_id)[i])
      wmd_idx <- which(wmd_df$cruise_id == unique(wmd_df$cruise_id)[i])
      wmd_df$date[wmd_idx] <- unique(uvp_data_meta$dateabv[meta_idx])
    }
  } else if(method == 'season') {
    wmd_df <- wmd_list |> list_to_tib('season')
  } else {
    stop('Method must be cruise or season')
  }
  return(wmd_df)
}



season_assign <- function(date) {
  mo = lubridate::month(date)
  if(mo %in% c(12,1,2)) {
    return('winter')
  } else if(mo %in% c(3,4,5)) {
    return('spring')
  } else if(mo %in% c(6,7,8)) {
    return('summer')
  } else if(mo %in% c(9,10,11)) {
    return('fall')
  }
}

# |- Useful for jags set-up and 04 model scripts ---------------

#function to get db into integer format
code_db <- function(tod,min_d,max_d, db = db){
  
  bin_limits <- get_bin_limtis(db[[tod]])
  d_idx <- which(bin_limits$min_d >= min_d &
                          bin_limits$max_d <= max_d)
  
  db_code <- as.factor(db[[tod]][d_idx]) |> 
    as.integer()
  
  return(list(
    code = db_code,
    idx = d_idx
    ))
}

#function to get db out of integer format
db_decode <- function(db_idx, tod, min_d, 
                      max_d, db = db) {
  
  bin_limits <- get_bin_limtis(db[[tod]])
  d_idx <- which(bin_limits$min_d >= min_d &
                   bin_limits$max_d <= max_d)
  
  #very ulgly if else patch
  if(is.factor(db[[tod]][d_idx])) {
    min_db <- as.factor(db[[tod]][d_idx]) |>
      as.integer() |> 
      min()
    
    code <- db_idx + (min_db - 1)
    
    
    return(levels(db[[tod]])[code])
  } else {
    
    return(unique(db[[tod]][d_idx]))
    
  }
}


# detect presence absence 
# function to add detection to the bins
was_detected <- function(df) {
  df$detected <- ifelse(df$x > 0, 1, 0)
  return(df)
}

# assign a survey/profile number
# For each cruise, a different profile is a visit to a site
# A f
surv_assign <- function(df) { 
  df$survey_id <- rep(NA, nrow(df))
  for(i in 1:length(unique(df$cruise_id))) {
    cruise_idx <- which(df$cruise_id == unique(df$cruise_id)[i])
    for(j in 1:length(unique(df$profileid[cruise_idx]))) {
      cast_idx <- which(df$profileid[cruise_idx] == unique(df$profileid[cruise_idx])[j])
      df$survey_id[cruise_idx][cast_idx] = j
    }
  }
  return(df)
}


listify_matrix <- function(df) {
  out_list <- list()
  
  out_list[['cruise_id']] <- unique(df$cruise_id)
  
  out_list[['detections']] <- df |> 
    select(c(cruise_id, detected, survey_id)) |> 
    pivot_wider(names_from = 'survey_id',values_from = detected) |> 
    select(-cruise_id) |> 
    as.matrix()
  
  out_list[['vol_s']] <- df |> 
    select(c(cruise_id, vol_sampled, survey_id)) |> 
    pivot_wider(names_from = survey_id, values_from = vol_sampled) |> 
    select(-cruise_id) |> 
    as.matrix()
  
  return(out_list)
}



max_fill<- function(x, max_cols) {
  if(ncol(x) > max_cols) {
    error('wrong max_cols')
  } else if (ncol(x) == max_cols) {
    return(x)
  } else {
    return(cbind(x, matrix(nrow = nrow(x), ncol = max_cols - ncol(x))))
  }
}

