require(tidyverse)
require(lubridate)

source('libs/lib.R')

union_suburb_parcels <- function(suburb_parcel_list, all_parcel_ids){
  suburb_parcels_union_df <- NULL
  for(suburb_parcels_df in suburb_parcel_list){
    if(is.null(suburb_parcels_union_df)){
      suburb_parcels_union_df <- suburb_parcels_df
    }
    else{
      suburb_parcels_extra_df <- suburb_parcels_df %>% filter(!(parcel_id %in% suburb_parcels_union_df$parcel_id))
      suburb_parcels_union_df <- rbind(suburb_parcels_union_df, suburb_parcels_extra_df)
    }
  }
  suburb_parcels_union_df <- suburb_parcels_union_df %>% filter(parcel_id %in% all_parcel_ids)
  suburb_parcels_union_sf <- st_as_sf(suburb_parcels_union_df, crs=4326, wkt=c('wkt'))
  return(suburb_parcels_union_sf)
}

# build set of canonical (unique) parcel ids that are within the same spatial extent over all time points
get_common_parcel_ids <- function(all_parcel_ids, suburb_all_vectors_utm_sf, survey_dates){
  for(sdate in survey_dates){
    df <- suburb_all_vectors_utm_sf %>% filter(survey_date == sdate)
    pids <- unique(df$parcel_id)
    sd <- base::setdiff(all_parcel_ids, pids)
    if(length(sd) != 0){
      all_parcel_ids <- base::setdiff(all_parcel_ids, sd)
    }
  }
  return(all_parcel_ids)
}

# build set of canonical (unique) parcel ids across the set
get_all_parcel_ids <- function(suburb_parcel_list){
  all_parcel_ids <- NULL
  for(suburb_parcels_df in suburb_parcel_list){
    all_parcel_ids <- c(all_parcel_ids, suburb_parcels_df$parcel_id)  
  }  
  all_parcel_ids <- unique(all_parcel_ids)
  return(all_parcel_ids)
}

make_parcel_vectors <- function(parcels.df){
  parcels.df$wkt <- paste0('POINT (', parcels.df$longitude, ' ', parcels.df$latitude, ')')
  
  parcels.df <- parcels.df %>%
    select(parcel_id, property_id, survey_date, address, longitude, latitude, map_browser, wkt)
    
  # create spatial from data frame
  parcels_sf <- st_as_sf(parcels.df, crs=4326, wkt=c('wkt'))
  # MGA94 /zone 56
  parcels_utm_sf <- st_transform(parcels_sf, crs = epsg.mga.zone56)
  return(parcels_utm_sf)
}

make_class_vectors <- function(vectors.df){
  vectors.df$wkt <- paste0('POINT (', vectors.df$longitude, ' ', vectors.df$latitude, ')')
  vectors.df <- vectors.df %>%
    select(parcel_id, property_id, survey_date, cls, type, address, longitude, latitude, confidence, map_browser, wkt) %>%
    arrange(desc(confidence))
  # create spatial from data frame
  vectors_sf <- st_as_sf(vectors.df, crs=4326, wkt=c('wkt'))
  # MGA94 /zone 56
  vectors_utm_sf <- st_transform(vectors_sf, crs = epsg.mga.zone56)
  return(vectors_utm_sf)
}

extract_class_vectors <- function(shape_name, sdate, roofs.df, objects.df){
  vector_list <- list()
  vector_list[['error']] <- NA
  
  # pools
  clz <- 1001
  pool.vectors.df <- objects.df %>% filter(feature_class == clz) %>%
    mutate(type = 'pool', cls = clz) %>%
    select(parcel_id, property_id, survey_date, cls, type, address, longitude, latitude, confidence, map_browser, wkt)
  if(nrow(pool.vectors.df) == 0){
    vector_list[[as.character(clz)]] <- NULL
  }
  else{
    pool.vectors.df <- make_class_vectors(pool.vectors.df)
    vector_list[[as.character(clz)]] <- pool.vectors.df
  }
  
  # roofs
  clz <- 1002
  roof.vectors.df <- roofs.df %>% 
    mutate(type = 'roof', cls = clz) %>% 
    select(parcel_id, property_id, survey_date, cls, type, address, longitude, latitude, confidence, map_browser, wkt)
  if(nrow(roof.vectors.df) == 0){
    vector_list[[as.character(clz)]] <- NULL
  }
  else{
    roof.vectors.df <- make_class_vectors(roof.vectors.df)
    vector_list[[as.character(clz)]] <- roof.vectors.df
  }
  
  # solar
  clz <- 1003
  solar.vectors.df <- objects.df %>% filter(feature_class == clz) %>%
    mutate(type = 'solar', cls = clz) %>%
    select(parcel_id, property_id, survey_date, cls, type, address, longitude, latitude, confidence, map_browser, wkt)
  if(nrow(solar.vectors.df) == 0){
    vector_list[[as.character(clz)]] <- NULL
  }
  else{
    solar.vectors.df <- make_class_vectors(solar.vectors.df)
    vector_list[[as.character(clz)]] <- solar.vectors.df
  }
  
  # construction
  clz <- 1004
  cons.vectors.df <- objects.df %>% filter(feature_class == clz) %>%
    mutate(type = 'construction', cls = clz) %>%
    select(parcel_id, property_id, survey_date, cls, type, address, longitude, latitude, confidence, map_browser, wkt)
  if(nrow(cons.vectors.df) == 0){
    vector_list[[as.character(clz)]] <- NULL
  }
  else{
    cons.vectors.df <- make_class_vectors(cons.vectors.df)
    vector_list[[as.character(clz)]] <- cons.vectors.df
  }
  
  if(nrow(roof.vectors.df) == 0 && nrow(cons.vectors.df) == 0){
    vector_list[['error']] <- paste('No vectors for', shape_name, 'for', sdate)
  }
  return(vector_list)
}

create_timeline_datasets <- function(shape_name, parcels.df, objects.df, roofs.df, survey_dates){
  timeline <- list()
  timeline[['error']] <- NA
  
  for(sdate in survey_dates){
    time_datasets <- list()
    dtext <- date_to_text(sdate)
    
    # parcels
    parcels.sd.df <- parcels.df %>% filter(survey_date == sdate) %>% split(.$parcel_id) %>% map_dfr(get_latest_survey)
    if(nrow(parcels.sd.df) == 0){
      timeline[['error']] <- paste(shape_name, ': no parcels found for date', dtext)
      return(timeline)
    }
    parcel.addrs.df <- parcels.sd.df %>% split(.$parcel_id) %>% map_dfr(ex.property.au)
    # join address and prop info to main frame
    parcels.sd.df <- parcels.sd.df %>% left_join(parcel.addrs.df, by="parcel_id")
    rm(parcel.addrs.df)
    
    # objects
    objects.sd.df <- objects.df %>% filter(survey_date == sdate) %>% split(.$parcel_id) %>% map_dfr(get_latest_survey)
    if(nrow(objects.sd.df) == 0){
      timeline[['error']] <- paste(shape_name, ': no objects found for date', dtext)
      return(timeline)
    }
    
    # roofs
    roofs.sd.df <- roofs.df %>% filter(survey_date == sdate) %>% split(.$parcel_id) %>% map_dfr(get_latest_survey)
    if(nrow(roofs.sd.df) == 0){
      timeline[['error']] <- paste(shape_name, ': no roofs found for date', dtext)
      return(timeline)
    }
    
    # all data present for this date so join parcel data
    objects.sd.df <- parcels.sd.df %>% 
      select(parcel_id, property_id, survey_date, survey_id, address) %>%
      inner_join(objects.sd.df, by=c("parcel_id", "survey_date", "survey_id"))
    
    roofs.sd.df <- parcels.sd.df %>% 
      select(parcel_id, property_id, survey_date, survey_id, address) %>%
      inner_join(roofs.sd.df, by=c("parcel_id", "survey_date", "survey_id"))
    # list for this date
    time_datasets[['parcels']] <- parcels.sd.df; time_datasets[['objects']] <- objects.sd.df; time_datasets[['roofs']] <- roofs.sd.df
    
    timeline[[dtext]] <- time_datasets
  }
  return(timeline)
}

load_datasets <- function(aoi_sf_albers, shape_name, data_dir){
  data_sets <- list()
  data_sets[['error']] <- NA
  
  # parcels
  parcels.df <- read.parcels(paste0(data_dir, '/parcels.csv.gz'))
  if(nrow(parcels.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no parcels found')
    return(data_sets)
  }
  parcels.df <- clip.vectors.to.shape(parcels.df, aoi_sf_albers, albers.crs)
  if(nrow(parcels.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no parcels found after clipping')
    return(data_sets)
  }
  
  # objects
  objects.df <- read.objects(paste0(data_dir, '/objects.csv.gz'))
  if(nrow(objects.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no object data found')
    return(data_sets)
  }
  objects.df <- clip.vectors.to.shape(objects.df, aoi_sf_albers, albers.crs)
  if(nrow(objects.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no object data found after clipping')
    return(data_sets)
  }
  objects.df <- objects.df %>% mutate(confidence = round(100.0*raw_confidence, digits=1))
  
  # roofs
  roofs.df <- read.roofs(paste0(data_dir, '/roofs.csv.gz'))
  if(nrow(roofs.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no roof data found')
    return(data_sets)
  }
  roofs.df <- clip.vectors.to.shape(roofs.df, aoi_sf_albers, albers.crs, src.geom = 'enhanced_wkt')
  if(nrow(roofs.df) == 0){
    data_sets[['error']] <- paste(shape_name, ': no roof data found after clipping')
    return(data_sets)
  }
  roofs.df <- roofs.df %>% mutate(confidence = round(100.0*confidence, digits=1)) %>% rename(wkt = enhanced_wkt)
  
  data_sets[['parcels']] <- parcels.df; data_sets[['objects']] <- objects.df; data_sets[['roofs']] <- roofs.df
  return(data_sets)
}