require(sp)
require(tidyverse)
require(fs)
require(lubridate)

FINAL_DIR <- '/media/sf_3TBBU/UNSWStats/masters_data/final'

ESSENTIAL_SHAPES_DIR <- '../../essential-data/shape'

DEFAULT.COLS <- c('suburb', 'loc', 'tx', 'ty', 'longitude', 'latitude')

load_predictions <- function(path_to_exp){
  preds_df <- NULL
  pred_assets <- dir_ls(path_to_exp, regexp = '*.tsv', recurse = T)
  for(asset in pred_assets){
    preds_df <- rbind(preds_df, read_tsv(asset))  
  }
  preds_df$residuals <- preds_df$y - preds_df$mu
  return(preds_df)
}

calc_rmses <- function(path_to_exp){
  preds_df <- load_predictions(path_to_exp)
  
  rmse_df <- preds_df %>% select(lga, residuals) %>% 
    group_by(lga) %>% 
    mutate(rmse = sqrt(mean(residuals^2))) %>% 
    select(lga, rmse) %>% distinct() %>% arrange(lga)  
  return(rmse_df)
}

load_all_data <- function(){
  all_data_sf <- st_read(paste0(FINAL_DIR, '/', 'all_data.gpkg')) %>% 
    rename(area.t3.loss = area.delta32, area.t2.loss = area.delta21) %>% 
    select(-ends_with('.500'), -ends_with('.250'), -area.t1, -area.t3) %>% 
    select(-starts_with('cons.act.')) %>%  # only 129 across Syd basin
    mutate(
      lga_id = as.numeric(factor(lga)),
      fold_id = ifelse(lga_id %in% c(2,8), 18, lga_id),
      loc = as.character(loc)
      ) %>%
    select(loc, fold_id, lga_id, lga, suburb, area.t3.loss, everything())
  return(all_data_sf)
}

extract_actual_loss <- function(full_sf){
  # load buffered polygons at 35m radius
  neg.loss_buff35_sf <- st_read(paste0(ESSENTIAL_SHAPES_DIR, '/negative-loss-35-buf.gpkg'))
  # extract spatial points
  actual_loss_sf <- full_sf %>% filter(loc %in% neg.loss_buff35_sf$loc)
  return(actual_loss_sf)
}

#######################################################################

load_canopy_loss <- function(path_to_vectors){
  c_loss_sf <- st_read(path_to_vectors)
  c_loss_sf <- c_loss_sf %>%
    mutate(
      loc = as.character(loc),
      lga = as.character(lga)
    )
  c_loss_sf <- c_loss_sf %>% select(area.t3.loss, loc, lga_id, lga, suburb, everything())
  return(c_loss_sf)
}

# load_canopy_loss_as_spatial <- function(){
#   return(as_Spatial(load_canopy_loss()))
# }

build_formula <- function(data_sp){
  covs <- paste(names(data_sp@data %>% 
                        select(-area.t3.loss, -fold_id, -loc, -lga_id, 
                               -lga, -suburb, -canopy.loss, -loss)), 
                sep = '', collapse = ' + ')
  f <- as.formula(paste('area.t3.loss ~', covs))
  return(f)
}

to_loss_dataframe <- function(c_loss_sf){
  c_loss_df <- cbind(as.data.frame(c_loss_sf) %>% select(-geom), st_coordinates(c_loss_sf))  
  return(c_loss_df)
}

load_suburb_wide <- function(path_to_file){
  sub_wide_df <- read_tsv(path_to_file, col_types = cols(tx = col_integer(), ty = col_integer())) %>% 
    mutate(loc = paste0(tx, '_', ty)) %>% 
    select(suburb, loc, everything())
  return(sub_wide_df)
}

# general purpose pivot full wide frame to long format
suburb_to_long <- function(wide_df){
  long_df <- wide_df %>% select(all_of(DEFAULT.COLS), date1, t1, area.prop1) %>% rename(date = date1, t = t1, area.prop = area.prop1)
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date2, t2, area.prop2) %>% rename(date = date2, t = t2, area.prop = area.prop2))
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date3, t3, area.prop3) %>% rename(date = date3, t = t3, area.prop = area.prop3))
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date4, t4, area.prop4) %>% rename(date = date4, t = t4, area.prop = area.prop4))
  return(long_df)
}

# long data frame at times t1 and t2 only
suburb_to_diffs_long <- function(wide_df){
  long_df <- wide_df %>% select(all_of(DEFAULT.COLS), date1, t1, area.prop1) %>% rename(date = date1, t = t1, area.prop = area.prop1)
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date2, t2, area.prop2) %>% rename(date = date2, t = t2, area.prop = area.prop2))
  return(long_df)
}

# long data frame at times t1 and t3 only
suburb_to_prediction_long <- function(wide_df){
  long_df <- wide_df %>% select(all_of(DEFAULT.COLS), date1, t1, area.prop1) %>% rename(date = date1, t = t1, area.prop = area.prop1)
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date3, t3, area.prop3) %>% rename(date = date3, t = t3, area.prop = area.prop3))
  return(long_df)
}

# long data frame w/o t4
suburb_to_forecast_long <- function(wide_df){
  long_df <- wide_df %>% select(all_of(DEFAULT.COLS), date1, t1, area.prop1) %>% rename(date = date1, t = t1, area.prop = area.prop1)
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date2, t2, area.prop2) %>% rename(date = date2, t = t2, area.prop = area.prop2))
  long_df <- rbind(long_df, wide_df %>% select(all_of(DEFAULT.COLS), date3, t3, area.prop3) %>% rename(date = date3, t = t3, area.prop = area.prop3))
  return(long_df)
}

build_aoi_dataset <- function(a_sf){
  a_sf %>%
    rename(area.t3.loss = area.delta32, area.t2.loss = area.delta21) %>% 
    select(-ends_with('.500'), -ends_with('.250'), -lga, -suburb, 
           -loc, -area.t1, -area.t3) %>% 
    select(-starts_with('cons.act.')) %>% # only 129 across Syd basin
    select(area.t3.loss, everything())  
}

sample_nonloss <- function(a_df, size, thres){
  non_loss_df <- a_df %>% filter(area.delta32 >= -thres)
  sample_df <- sample_n(non_loss_df, size=size, replace=F)
  return(sample_df)
}

trim_fold <- function(f_df){
  f_df %>%
    rename(area.t3.loss = area.delta32, area.t2.loss = area.delta21) %>% 
    select(-ends_with('.500'), -ends_with('.250'), -suburb, 
           -loc, -area.t1, -area.t3) %>% 
    select(-starts_with('cons.act.')) %>% # only 129 across Syd basin
    select(lga, area.t3.loss, canopy.loss, everything())  
}