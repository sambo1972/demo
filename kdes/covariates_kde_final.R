rm(list=ls())
library(sp)
library(sf)
library(gstat)
library(fields)
library(tidyverse)
library(fs)
library(lubridate)
library(parallel)

source('libs/lib.R')
source('libs/dataframes.R')
source('libs/covariates.R')

####################################################################################
# PATHS
####################################################################################

canopy_data_dir <- '/media/sf_3TBBU/UNSWStats/masters_data/canopy'
change_data_dir <- '/media/sf_3TBBU/UNSWStats/masters_data/change_vectors'
target_dir <- '/media/sf_3TBBU/UNSWStats/masters_data/kde_covariates'
if(!dir_exists(target_dir)){
  dir_create(target_dir, recurse = T)
}

# load change datasets
solar_new_sf <- st_read(paste0(change_data_dir, '/new_solar.gpkg'))
pools_new_sf <- st_read(paste0(change_data_dir, '/new_pools.gpkg'))
construction_new_sf <- st_read(paste0(change_data_dir, '/new_construction.gpkg'))
construction_completed_sf <- st_read(paste0(change_data_dir, '/completed_construction.gpkg'))
construction_active_sf <- st_read(paste0(change_data_dir, '/active_construction.gpkg'))
roofs_new_sf <- st_read(paste0(change_data_dir, '/new_roofing.gpkg'))
demolitions_sf <- st_read(paste0(change_data_dir, '/demolitions.gpkg'))
roofs_sf <- st_read(paste0(change_data_dir, '/existing_roofing.gpkg'))

canopy_assets <- dir_ls(canopy_data_dir, regexp = '*.tsv.gz')

canopy_assets <- canopy_assets[4:5]

cluster <- makeCluster(detectCores())
clusterEvalQ(cluster, {
  library(tidyverse)
  library(sf)
  library(fields)
})
clusterExport(cluster, c('epsg.mga.zone56', 'BANDWIDTHS', 'kde.gauss', 'BUFFER.DIST.METRES'))
# extract kdes over all suburbs
for(canopy_asset in canopy_assets){
  # load suburb data
  suburb <- path_ext_remove(path_ext_remove(path_file(canopy_asset)))
  suburb_target_dir <- paste0(target_dir, '/', suburb)
  if(!dir_exists(suburb_target_dir)){
    dir_create(suburb_target_dir, recurse = T)
  }
  suburb_df_wide <- load_suburb_wide(canopy_asset) %>% 
    mutate(
      area.delta = round(area.prop3 - area.prop1, digits=4)
    ) %>% 
    select(suburb, loc, longitude, latitude, area.delta)
  suburb_sf_4326 <- st_as_sf(suburb_df_wide, crs=4326, coords = c('longitude', 'latitude'))
  suburb_sf_utm <- st_transform(suburb_sf_4326, crs=epsg.mga.zone56)
  rm(suburb_sf_4326)
  rm(suburb_df_wide)
  # create spatially localised subsets of above
  suburb_buf_boundary <- st_union(st_buffer(suburb_sf_utm, BUFFER.DIST.METRES))
  sub_solar_new_sf <- st_intersection(solar_new_sf, suburb_buf_boundary)
  sub_pools_new_sf <- st_intersection(pools_new_sf, suburb_buf_boundary)
  sub_cons_new_sf <- st_intersection(construction_new_sf, suburb_buf_boundary)
  sub_cons_comp_sf <- st_intersection(construction_completed_sf, suburb_buf_boundary)
  sub_cons_act_sf <- st_intersection(construction_active_sf, suburb_buf_boundary)
  sub_roofs_new_sf <- st_intersection(roofs_new_sf, suburb_buf_boundary)
  sub_roofs_dem_sf <- st_intersection(demolitions_sf, suburb_buf_boundary)
  sub_roofs_sf <- st_intersection(roofs_sf, suburb_buf_boundary) 
  # extract canopy coordinates
  suburb_sf_coords <- st_coordinates(suburb_sf_utm)
  # generate and save the kdes
  print(paste(suburb, ': new pools'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=POOLS_NEW, cov_sf_utm=sub_pools_new_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', POOLS_NEW, '.tsv.gz'))
  print(paste(suburb, ': new solar'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=SOLAR_NEW, cov_sf_utm=sub_solar_new_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', SOLAR_NEW, '.tsv.gz'))
  print(paste(suburb, ': new construction'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=CONSTRUCTION_NEW, cov_sf_utm=sub_cons_new_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', CONSTRUCTION_NEW, '.tsv.gz'))
  print(paste(suburb, ': completed construction'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=CONSTRUCTION_COMPLETED, cov_sf_utm=sub_cons_comp_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', CONSTRUCTION_COMPLETED, '.tsv.gz'))
  print(paste(suburb, ': active construction'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=CONSTRUCTION_ACTIVE, cov_sf_utm=sub_cons_act_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', CONSTRUCTION_ACTIVE, '.tsv.gz'))
  print(paste(suburb, ': new roofs'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=ROOFS_NEW, cov_sf_utm=sub_roofs_new_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', ROOFS_NEW, '.tsv.gz'))
  print(paste(suburb, ': demolitions'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=ROOFS_DEMOLISHED, cov_sf_utm=sub_roofs_dem_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', ROOFS_DEMOLISHED, '.tsv.gz'))
  print(paste(suburb, ': active roofs'))
  kde.mtrx <- parApply(cluster, suburb_sf_coords, MARGIN=1, gen_point_kde_covariates, cov.name=ROOFS_ACTIVE, cov_sf_utm=sub_roofs_sf)
  save_data(as.data.frame(t(kde.mtrx)), paste0(suburb_target_dir, '/', ROOFS_ACTIVE, '.tsv.gz'))
  # bind together and save as vector asset
  suburb_covariates_df <- map_dfc(dir_ls(suburb_target_dir, regexp = '*.tsv.gz'), read_tsv)
  suburb_sf_utm <- st_as_sf(cbind(as.data.frame(suburb_sf_utm), suburb_covariates_df))
  path_to_vectors <- paste0(suburb_target_dir, '/', suburb, '.gpkg')
  save_vectors(suburb_sf_utm, path_to_vectors)
  print(paste(suburb, ': saved to', path_to_vectors))
}
stopCluster(cluster)
print('Done')