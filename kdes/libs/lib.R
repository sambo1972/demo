require(tidyverse)
require(lubridate)
require(sf)
require(fs)

LON_LAT_PRECISION = 6
TILE_SIZE = 256
ZOOM = 21
REZ = 0.074645535411904
PIX_AREA = REZ**2

epsg.mga.zone56 <- 28356
SURVEY_DATES <- c(as_date("2017-02-11"), as_date("2018-12-27"), as_date("2019-10-22"))

POOLS_NEW <- 'pools.new'
SOLAR_NEW <- 'solar.new'
CONSTRUCTION_NEW <- 'cons.new'
CONSTRUCTION_COMPLETED <- 'cons.comp'
CONSTRUCTION_ACTIVE <- 'cons.act'
ROOFS_NEW <- 'roofs.new'
ROOFS_DEMOLISHED <- 'demolitions'
ROOFS_ACTIVE <- 'roofs'

COVARIATE_GROUPS <- c(CONSTRUCTION_COMPLETED, CONSTRUCTION_NEW, ROOFS_DEMOLISHED, 
                      POOLS_NEW, ROOFS_NEW, ROOFS_ACTIVE, SOLAR_NEW)


save_vectors <- function(vectors_sf, path_to_file){
  if(file_exists(path_to_file)){
    file_delete(path_to_file)
  }
  st_write(vectors_sf, path_to_file)
}

save_data <- function(data_df, path_to_file){
  if(file_exists(path_to_file)){
    file_delete(path_to_file)
  }
  write_tsv(data_df, path_to_file)
}

safe_identifier <- function(s){
  str_replace_all(s, ' ', '_')  
}

save_model <- function(mdl, target_dir, identifier){
  path_to_model <- paste0(target_dir, '/', safe_identifier(identifier), '.RDS')
  if(file_exists(path_to_model)){
    file_delete(path_to_model)
  }
  saveRDS(mdl, file = path_to_model)
  return(path_to_model)
}

get_months <- function(d1, d2){
  p <- as.period(interval(d1, d2), unit = 'day')
  months <- round(p@day/30, digits=1)
  return(months)
}

date_to_text <- function(sd){
  dtext <- as_date(sd) %>% format('%Y-%m-%d')
  return(dtext)
}

get_state_postcode <- function(path){
  m <- str_split(str_trim(path), '/', simplify = T)
  state_postcode = m[1,8]
  m <- str_split(state_postcode, '_', simplify = T)
  state <- m[1,1]
  postcode <- m[1,2]
  return(c(state, postcode))
}

# function to map months to seasons
season <- function(sd){
  return(ifelse(month(sd) %in% c(9,10,11), 'spring', 
                ifelse(month(sd) %in% c(12,1,2), 'summer', 
                       ifelse(month(sd) %in% c(3,4,5), 'autumn', 'winter'))))
}
