require(fields)
require(tidyverse)

BANDWIDTHS <- c(10, 50, 100, 250, 500)
BUFFER.DIST.METRES <- 500

# for each row (response location) sum across the columns (covariates)
kde.gauss <- function(dist.vec, bndw){
  return(sum(exp(-(1/(2*bndw^2)) * dist.vec^2)))
}

gen_point_kde_covariates <- function(coords_utm, cov.name, cov_sf_utm, bandwidths=BANDWIDTHS, 
                                     kde.function = kde.gauss, buffer.dist=BUFFER.DIST.METRES){
  X <- coords_utm['X']; Y <- coords_utm['Y']
  
  # create point and a buffer
  point_canopy_sf_utm <- st_as_sf(st_sfc(st_point(c(X, Y))), crs=epsg.mga.zone56)
  point_buf_sf <- st_buffer(point_canopy_sf_utm, buffer.dist)
  point_covariates_sf <- st_intersection(cov_sf_utm, point_buf_sf)
  
  # calculate distance matrix
  canopy.mtrx <- matrix(c(X, Y), nrow=1, byrow = T)
  cov.mtrx <- st_coordinates(point_covariates_sf)
  d.mtrx <- rdist(canopy.mtrx, cov.mtrx)
  
  kde.vals <- NULL
  for(b in bandwidths){
    kde.val <- apply(d.mtrx, MARGIN = 1, FUN = kde.function, bndw=b)
    kde.vals <- c(kde.vals, kde.val)
  }
  names(kde.vals) <- paste0(cov.name, ".", bandwidths)
  
  return(kde.vals)
}

gen_kde_covariates <- function(cov.name, canopy_df_utm, cov_df_utm, bandwidths, kde.function = kde.gauss){
  # calculate distance matrix
  canopy.mtrx <- st_coordinates(canopy_df_utm)
  cov.mtrx <- st_coordinates(cov_df_utm)
  d.mtrx <- rdist(canopy.mtrx, cov.mtrx)
  
  # local function to map across columns
  kde_over_bandwidths <- function(b, dist.mx, kde.fun){
    f.hat <- apply(dist.mx, MARGIN = 1, FUN = kde.fun, bndw=b)
    return(f.hat)
  }
  
  # map over bandwidths and column bind into resulting data frame
  kde_cov_df <- map_dfc(bandwidths, kde_over_bandwidths, dist.mx=d.mtrx, kde.fun=kde.function)
  names(kde_cov_df) <- paste0(cov.name, ".", bandwidths)  
  
  return(kde_cov_df)
}

gen_construction_covariates <- function(suburb_sf_utm, construction_new_sf, construction_completed_sf, construction_active_sf){
  cons_list <- list()
  cons_list[[CONSTRUCTION_NEW]] <- gen_kde_covariates(cov.name = CONSTRUCTION_NEW, suburb_sf_utm, construction_new_sf, bandwidths = BANDWIDTHS)
  cons_list[[CONSTRUCTION_COMPLETED]] <- gen_kde_covariates(cov.name = CONSTRUCTION_COMPLETED, suburb_sf_utm, construction_completed_sf, bandwidths = BANDWIDTHS)
  cons_list[[CONSTRUCTION_ACTIVE]] <- gen_kde_covariates(cov.name = CONSTRUCTION_ACTIVE, suburb_sf_utm, construction_active_sf, bandwidths = BANDWIDTHS)
  return(cons_list)
}

gen_roof_covariates <- function(suburb_sf_utm, roofs_new_sf, demolitions_sf, roofs_sf){
  roofs_list <- list()
  roofs_list[[ROOFS_NEW]] <- gen_kde_covariates(cov.name = ROOFS_NEW, suburb_sf_utm, roofs_new_sf, bandwidths = BANDWIDTHS)
  roofs_list[[ROOFS_DEMOLISHED]] <- gen_kde_covariates(cov.name = ROOFS_DEMOLISHED, suburb_sf_utm, demolitions_sf, bandwidths = BANDWIDTHS)
  roofs_list[[ROOFS_ACTIVE]] <- gen_kde_covariates(cov.name = ROOFS_ACTIVE, suburb_sf_utm, roofs_sf, bandwidths = BANDWIDTHS)
  return(roofs_list)
}

#######################################################################