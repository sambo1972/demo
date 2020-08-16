require(sp)
require(tidyverse)
require(fs)
require(lubridate)
require(FRK)
require(glmnet)

select_variables <- function(dataset_df, lambda.tst = exp(seq(-10,0,0.2)), min.magnitude = NA, apply.limits=F){
  y <- dataset_df$area.t3.loss; x <- as.matrix(dataset_df[,6:51])
  # fit Lasso model using lga id as fold identifier
  model.cv.lasso.fit <- cv.glmnet(x, y, alpha=1, 
                                  foldid = dataset_df$lga_id,
                                  parallel = T,
                                  family='gaussian',
                                  lambda = lambda.tst)
  
  # re-fit using entire set and use predict function to get coefficients
  model.lasso.full <- glmnet(x, y, family='gaussian', alpha=1)
  model.lasso.coeffs <- predict(model.lasso.full, type='coefficients', 
                                s=model.cv.lasso.fit$lambda.min)
  # make a dataframe out of coefficients
  coeffs <- model.lasso.coeffs[2:nrow(model.lasso.coeffs),1]
  coeffs_df <- data.frame(variable=as.character(names(coeffs)), coefficient=coeffs)
  
  # select variables subject to 1) no more than one present in each covariate group
  # and 2) if present, at least min.magnitude in absolute magnitude
  if(apply.limits){
    variables <- c('area.t2', 'area.t2.loss')
    for(cov_grp in COVARIATE_GROUPS){
      group_df <- coeffs_df %>% filter(str_starts(variable, cov_grp) & coefficient != 0) %>% 
        top_n(n=1, wt=abs(coefficient))
      if(!is.na(min.magnitude)){
        group_df <- group_df %>% filter(abs(coefficient) >= min.magnitude)  
      }
      if(nrow(group_df) == 1){
        variables <- c(variables, as.character(group_df$variable[1]))
      }
    }
  }else{
    non_zero_df <- coeffs_df %>% filter(coefficient != 0)
    print(non_zero_df)
    variables <- as.character(non_zero_df$variable)
  }
  return(unique(variables))
}

add_poly_squares <- function(closs_sf){
  area.t2.mtrx <- poly(closs_sf$area.t2, 2, simple=T)
  area.t2.loss.mtrx <- poly(closs_sf$area.t2.loss, 2, simple=T)
  
  cons.comp.10.mtrx <- poly(closs_sf$cons.comp.10, 2, simple=T)
  cons.comp.50.mtrx <- poly(closs_sf$cons.comp.50, 2, simple=T)
  cons.comp.100.mtrx <- poly(closs_sf$cons.comp.100, 2, simple=T)
  
  cons.new.10.mtrx <- poly(closs_sf$cons.new.10, 2, simple=T)
  cons.new.50.mtrx <- poly(closs_sf$cons.new.50, 2, simple=T)
  cons.new.100.mtrx <- poly(closs_sf$cons.new.100, 2, simple=T)
  
  demolitions.10.mtrx <- poly(closs_sf$demolitions.10, 2, simple=T)
  demolitions.50.mtrx <- poly(closs_sf$demolitions.50, 2, simple=T)
  demolitions.100.mtrx <- poly(closs_sf$demolitions.100, 2, simple=T)
  
  pools.new.10.mtrx <- poly(closs_sf$pools.new.10, 2, simple=T)
  pools.new.50.mtrx <- poly(closs_sf$pools.new.50, 2, simple=T)
  pools.new.100.mtrx <- poly(closs_sf$pools.new.100, 2, simple=T)
  
  roofs.new.10.mtrx <- poly(closs_sf$roofs.new.10, 2, simple=T)
  roofs.new.50.mtrx <- poly(closs_sf$roofs.new.50, 2, simple=T)
  roofs.new.100.mtrx <- poly(closs_sf$roofs.new.100, 2, simple=T)
  
  roofs.10.mtrx <- poly(closs_sf$roofs.10, 2, simple=T)
  roofs.50.mtrx <- poly(closs_sf$roofs.50, 2, simple=T)
  roofs.100.mtrx <- poly(closs_sf$roofs.100, 2, simple=T)
  
  solar.new.10.mtrx <- poly(closs_sf$solar.new.10, 2, simple=T)
  solar.new.50.mtrx <- poly(closs_sf$solar.new.50, 2, simple=T)
  solar.new.100.mtrx <- poly(closs_sf$solar.new.100, 2, simple=T)
  
  closs_sf <- closs_sf %>% 
    mutate(
      area.t2 = area.t2.mtrx[,1], area.t2.sq = area.t2.mtrx[,2],
      area.t2.loss = area.t2.loss.mtrx[,1], area.t2.loss.sq = area.t2.loss.mtrx[,2],
      
      cons.comp.10 = cons.comp.10.mtrx[,1], cons.comp.10.sq = cons.comp.10.mtrx[,2],
      cons.comp.50 = cons.comp.50.mtrx[,1], cons.comp.50.sq = cons.comp.50.mtrx[,2],
      cons.comp.100 = cons.comp.100.mtrx[,1], cons.comp.100.sq = cons.comp.100.mtrx[,2],
      
      cons.new.10 = cons.new.10.mtrx[,1], cons.new.10.sq = cons.new.10.mtrx[,2],
      cons.new.50 = cons.new.50.mtrx[,1], cons.new.50.sq = cons.new.50.mtrx[,2],
      cons.new.100 = cons.new.100.mtrx[,1], cons.new.100.sq = cons.new.100.mtrx[,2],
      
      demolitions.10 = demolitions.10.mtrx[,1], demolitions.10.sq = demolitions.10.mtrx[,2],
      demolitions.50 = demolitions.50.mtrx[,1], demolitions.50.sq = demolitions.50.mtrx[,2],
      demolitions.100 = demolitions.100.mtrx[,1], demolitions.100.sq = demolitions.100.mtrx[,2],
      
      pools.new.10 = pools.new.10.mtrx[,1], pools.new.10.sq = pools.new.10.mtrx[,2],
      pools.new.50 = pools.new.50.mtrx[,1], pools.new.50.sq = pools.new.50.mtrx[,2],
      pools.new.100 = pools.new.100.mtrx[,1], pools.new.100.sq = pools.new.100.mtrx[,2],
      
      roofs.new.10 = roofs.new.10.mtrx[,1], roofs.new.10.sq = roofs.new.10.mtrx[,2],
      roofs.new.50 = roofs.new.50.mtrx[,1], roofs.new.50.sq = roofs.new.50.mtrx[,2],
      roofs.new.100 = roofs.new.100.mtrx[,1], roofs.new.100.sq = roofs.new.100.mtrx[,2],
      
      roofs.10 = roofs.10.mtrx[,1], roofs.10.sq = roofs.10.mtrx[,2],
      roofs.50 = roofs.50.mtrx[,1], roofs.50.sq = roofs.50.mtrx[,2],
      roofs.100 = roofs.100.mtrx[,1], roofs.100.sq = roofs.100.mtrx[,2],
      
      solar.new.10 = solar.new.10.mtrx[,1], solar.new.10.sq = solar.new.10.mtrx[,2],
      solar.new.50 = solar.new.50.mtrx[,1], solar.new.50.sq = solar.new.50.mtrx[,2],
      solar.new.100 = solar.new.100.mtrx[,1], solar.new.100.sq = solar.new.100.mtrx[,2]
    )
  return(closs_sf)
}

# sample observations in each LGA within +/- measurement error 
# which are effectively zero loss events. 
sample.zeros <- function(full_sf, measurement.error, sample_factor=10){
  zero_sample_locs <- NULL
  for(l in unique(full_sf$lga)){
    print(paste('Sampling lga', l))
    
    lga_actual_loss_sf <- full_sf %>% filter(lga == l & area.t3.loss < -measurement.error)
    
    print(paste('Actual loss', nrow(lga_actual_loss_sf)))
    
    lga_zero_loss_sf <- full_sf %>% filter(lga == l & area.t3.loss <= measurement.error & area.t3.loss >= -measurement.error)
    
    print(paste('Zero loss rows', nrow(lga_zero_loss_sf)))
    
    # number of actual loss obs times factor
    sample.size <- min(nrow(lga_actual_loss_sf) * sample_factor, nrow(lga_zero_loss_sf))
    print(paste('Sampling', sample.size, 'rows'))
    
    lga_zero_sample_sf <- lga_zero_loss_sf %>% sample_n(size=sample.size, replace=F)
    
    zero_sample_locs <- c(zero_sample_locs, lga_zero_sample_sf$loc)
  }
  zero_sample_sf <- full_sf %>% filter(loc %in% zero_sample_locs)
  return(zero_sample_sf)
}

fit.sre.scv.intra.lga <- function(exp_dir, mdl_form, test_fold_id, data_sp, res, fn='bisquare',
                                  fs_obs=F, estimate.err=T, std = 1){

  # test fold
  test_fold_sp <- subset(data_sp, scv_fold == test_fold_id)
  test_identifier <- safe_identifier(as.character(test_fold_sp[1,]$lga))

  # training data is the response at all locations (except test fold)
  training_folds <- base::setdiff(unique(data_sp@data$scv_fold), test_fold_id)
  training_sp <- subset(data_sp, scv_fold %in% training_folds)[, c('area.t3.loss')]

  fit.sre.fold.delegate(exp_dir, mdl_form, test_fold_id, test_identifier, test_fold_sp,
                        training_sp, data_sp, res, fn, fs_obs, by.lga=F, estimate.err, std)
}

fit.sre.scv.by.lga <- function(exp_dir, mdl_form, test_fold_id, data_sp, res, fn='bisquare', 
                               fs_obs=F, estimate.err=T, std = 1){
  
  # test fold
  test_fold_sp <- subset(data_sp, lga_id == test_fold_id)
  test_identifier <- safe_identifier(as.character(test_fold_sp[1,]$lga))
  
  # training data is the response at all locations (except test fold)
  training_folds <- base::setdiff(unique(data_sp@data$lga_id), test_fold_id)
  training_sp <- subset(data_sp, lga_id %in% training_folds)[, c('area.t3.loss')]
  
  fit.sre.fold.delegate(exp_dir, mdl_form, test_fold_id, test_identifier, test_fold_sp,
                        training_sp, data_sp, res, fn, fs_obs, by.lga=T, estimate.err, std)
}

fit.sre.fold.delegate <- function(exp_dir, mdl_form, test_fold_id, test_identifier, test_fold_sp, 
                                  training_sp, data_sp, res, fn, fs_obs, by.lga, estimate.err, std){
  # initialise logging
  basicConfig(level='INFO')
  
  path_to_log <- paste0(exp_dir, '/', test_identifier, '.log')
  new.log <- !file_exists(path_to_log)
  addHandler(writeToFile, file=path_to_log, level='INFO')
  
  if(new.log){
    loginfo('# *************** Model setup ***************')
    loginfo(paste0('Formula\t', mdl_form))
    loginfo(paste0('Basis function\t', fn))
    loginfo(paste0('Resolutions\t', res))
    loginfo(paste0('Est error\t', estimate.err)) 
    if(estimate.err == F){
      loginfo(paste0('Error\t', std))
    }
  }
  loginfo(paste0('\nTest identifier\t', test_identifier))
  loginfo(paste0('Test fold\t', test_fold_id))
  loginfo(paste0('Test rows\t', nrow(test_fold_sp)))
  loginfo(paste0('Training rows\t', nrow(training_sp)))
  
  if(estimate.err == F){
    training_sp$std <- std
  }
  
  # build BAUs - include all covariates but not the response
  aoi_BAUs <- data_sp[,c(-1)]
  aoi_BAUs$fs <- 1 # fine scale variation at BAU level
  
  # construct basis
  aoi_basis <- auto_basis(manifold = plane(), data=data_sp, nres = res, type=fn, regular=0)
  
  dt_start <- lubridate::now()
  S <- SRE(as.formula(mdl_form), data = list(training_sp), 
           BAUs = aoi_BAUs, basis = aoi_basis, 
           est_error = estimate.err, average_in_BAU = F)
  elapsed <- lubridate::now() - dt_start
  loginfo(paste0('Binning & error estimation elapsed\t', round(elapsed, 2), ' ', units(elapsed)))
  
  dt_start <- lubridate::now()
  S <- SRE.fit(SRE_model = S, n_EM = 5000, tol = 0.01, print_lik=F)
  elapsed <- lubridate::now() - dt_start
  loginfo(paste0('Model fit elapsed\t', round(elapsed, 2), ' ', units(elapsed)))
  
  # predict at all locations
  aoi_preds_df <- as(predict(S, obs_fs = fs_obs), "data.frame")
  
  # extract predictions for test fold ordered by yhat (mu)
  if(by.lga){
    test_preds_df <- aoi_preds_df %>% filter(lga_id == test_fold_id)
  }
  else{
    test_preds_df <- aoi_preds_df %>% filter(scv_fold == test_fold_id)
  }
  test_preds_df <- test_preds_df %>% 
    mutate(y = test_fold_sp$area.t3.loss) %>% 
    select(y, mu, sd, everything()) %>% 
    arrange(y)
  
  loginfo(paste0('Fold prediction rows\t', nrow(test_preds_df)))
  
  # calculate RMS prediction error
  rmspe <- sqrt(mean((test_fold_sp$area.t3.loss - test_preds_df$mu)^2))
  
  loginfo(paste0('RMSPE\t', round(rmspe, 2)))
  
  # return list
  ret.lst = list("identifier" = test_identifier, "fold_id" = test_fold_id, "model" = S, 
                 "predictions" = test_preds_df, "rmspe" = rmspe)
  return(ret.lst)
}

fit.sre.full <- function(exp_dir, mdl_form, data_sp, res, fn='bisquare', fs_obs=F, estimate.err=T, std = 1){
  # training data is the response at all locations
  training_sp <- data_sp[, c('area.t3.loss')]
  if(estimate.err == F){
    training_sp$std <- std
  }
  # build BAUs - include all covariates but not the response
  aoi_BAUs <- data_sp[,c(-1)]
  aoi_BAUs$fs <- 1 # fine scale variation at BAU level
  
  # construct basis
  aoi_basis <- auto_basis(manifold = plane(), data=data_sp, nres = res, type=fn, regular=0)
  
  S <- SRE(as.formula(mdl_form), data = list(training_sp), 
           BAUs = aoi_BAUs, basis = aoi_basis, 
           est_error = estimate.err, average_in_BAU = F)
  S <- SRE.fit(SRE_model = S, n_EM = 5000, tol = 0.01, print_lik=F)
  return(S)
}

############################################################
# Models
############################################################

# fit a simple temporal linear model at each location
fit_lm_at_location <- function(data){
  mdl.lm <- lm(area.prop ~ 1 + t, data = data)
}

predict_lm_at_location <- function(mdl, t_pred){
  predict(mdl, newdata = t_pred, interval = 'prediction') %>% 
    data.frame() %>% 
    mutate(se = (upr - lwr)/(2 * 1.96)) %>% 
    select(fit, se)
}

predict_with_simple_linear <- function(sub_df){
  # fit simple linear model over time
  locations_lm <- sub_df %>% 
    group_by(longitude, latitude) %>% 
    nest() %>% 
    mutate(model = purrr::map(data, fit_lm_at_location)) %>% 
    mutate(model_df = purrr::map(model, tidy)) 
  # predict missing time at each location
  predictions_t2 <- locations_lm %>% 
    mutate(preds = purrr::map(model, predict_lm_at_location, t_pred=data.frame(t=13.8))) %>% 
    unnest(preds)
  return(predictions_t2$fit)
}

##############################################################################################

fit.sre <- function(a_sf, basis_type='bisquare', num_resolutions=3, max_em_iter=500, verbose=F){
  aoi_sp <- as_Spatial(aoi_sf)
  # build formula
  covs <- paste(names(aoi_sp@data %>% select(-area.t3.loss)), sep = '', collapse = ' + ')
  f <- as.formula(paste('area.t3.loss ~', covs))
  # basis
  aoi_basis <- auto_basis(manifold = plane(), data=aoi_sp, nres = num_resolutions,
                          type=basis_type, regular=0)
  # build BAUs - include all covariates but not the response
  aoi_BAUs <- aoi_sp[,c(-1)]
  aoi_BAUs$fs <- 1 # fine scale variation at BAU level

  # construct spatial random effects model
  S <- SRE(f = f,
           data = list(aoi_sp[,c(1)]), # data frame should include response and locations only
           BAUs = aoi_BAUs,
           basis = aoi_basis,
           est_error = T,
           average_in_BAU = F)
  # fit
  S <- SRE.fit(SRE_model = S, n_EM = max_em_iter, tol = 0.01, print_lik=verbose)
  return(S)
}