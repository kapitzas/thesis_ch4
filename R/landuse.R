rm(list = ls())

#------------------------#
#### 0. Load packages ####
#------------------------#

source("./R/99_helperfunctions.R")
#detach("package:flutes", unload=TRUE)
#devtools::install_github("kapitzas/flutes", force = TRUE)
require(flutes)
require(devtools)
require(WorldClimTiles)
require(rgdal)
require(gdaltools)
require(raster)
require(stringr)
require(sf)
require(flutes)
require(tidyverse)
require(nnet)

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")
temp_path <- file.path("/Volumes", "external", "c3 processing", "temp")
out_path <- file.path(".", "output")
data_path <- file.path(".", "data")

#-----------------------#
#### 1. Prepare data ####
#-----------------------#

# 1. a load map containing the 30 GTAP regions of our aggregation
gtap_borders <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
table_gtap <- table(gtap_borders[])

# 1. b load land use and present and future environmental variables
lu_global <- stack(list.files(processed_path, pattern = "lu_5min", full.names = TRUE, recursive = TRUE))
lu_global <- lu_global[[which(grepl(pattern = paste0(c("^urb", "^crop", "^pas", "^pri", "^sec"), collapse = "|"), names(lu_global)))]]
covs_files <- list.files(processed_path, pattern = "5min", full.names = TRUE, recursive = TRUE)
covs_global <- stack(covs_files[which(!grepl(pattern = paste0(c("urb", "crop", "pas", "pri", "sec"), collapse = "|"), covs_files))])
# 1. c Calculate time step increments in covariates that change in the future (bioclim and pop density)

# subset the present time steps
covs_pres <- covs_global[[-which(grepl("rcp|ssp|biorealms", names(covs_global)))]]
covs_pres <- covs_pres[[which(grepl("bio|pop", names(covs_pres)))]]

# subset future time steps (bio under rcp45 and rcp85 and pop under ssp2 and ssp5)
covs_fut45 <- covs_global[[which(grepl("rcp45|ssp2", names(covs_global)))]]
covs_fut85 <- covs_global[[which(grepl("rcp85|ssp5", names(covs_global)))]]
covs_fut85 <- covs_fut85[[c(which(grepl("pop", names(covs_fut85))), which(grepl("bio", names(covs_fut85))))]] # reod
covs_fut45 <- covs_fut45[[c(which(grepl("pop", names(covs_fut45))), which(grepl("bio", names(covs_fut45))))]]

# calculate increments in the changing time steps (2020-2100 means we have 8 10-year time steps)
cov_fut_increment_45 <- (covs_fut45 - covs_pres)/8
cov_fut_increment_85 <- (covs_fut85 - covs_pres)/8
names(cov_fut_increment_85) <- names(cov_fut_increment_45) <- names(covs_pres)

covs_incr <- c(cov_fut_increment_45, cov_fut_increment_85, cov_fut_increment_85) #2x 85, bc validation scenario (LUH1 demand) is also 85

# 1.d load mask
mask <- covs_global[[grep(paste(c("mask"), collapse = "|"), names(covs_global))]]

# 1.e Calculate initial neighbourhood rasters
# weights (see ?flutes::neighbourhood)
weights <- list(matrix(1/9, 3, 3, byrow= TRUE)) #size of window
weights <- rep(weights, length.out = 5)

# get non-na cells
inds_global <- which(!is.na(mask[]))
neigh_global <- neighbourhood(lu_global, weights = weights, mask = mask, suffix = "neigh", cols = 1:5, format = "stack")

#----------------------------#
#### 2. Suitability model ####
#----------------------------#
# 2. a subset covs to include the ones we need in suitability model
covs_suitmodel <- covs_global[[-grep(paste(c("biorealms", "rcp45", "rcp85", "lu", "ssp", "gtap", "unsubregions", "gls", "pa", "mask"), collapse = "|"), names(covs_global))]]

# 2. b get protected areas
pa_global <- covs_global[[grep(paste(c("pa_"), collapse = "|"), names(covs_global))]]


# 2 c. Estimate suit model for each region

# get gtap regions (names and our id)
country_code <- read.csv(file.path(data_path, "GTAP_regions.csv"))

country_code <- distinct(country_code[,c(3,4)])

# iterate through regions
suitmodels <- list()
i <- 1
for(i in 1:30){
  print(i)
  # current region name and code
  regions <- country_code$GTAP_aggregation[i]
  map_code <- country_code$GTAP_code[i]

  # get regional data from global data set:
  # indices of current region
  inds <- which(gtap_borders[]==map_code)

  # land use data
  landuse <- lu_global[inds]
  colnames(landuse) <- unlist(lapply(strsplit(colnames(landuse), "_"), FUN = function(x) x[[1]]))

  # env covariates
  covs <- covs_suitmodel[inds]
  neigh <- neigh_global[inds]
  data <- cbind(covs, neigh)

  ssize <- 10000
  if(nrow(covs) < ssize) ssize <- nrow(covs)

  # Determine correlations in data and reduce predictor set
  preds <- colnames(correlations(data, sub = sample(1:nrow(covs), ssize)))
  data <- data[,preds]

  # Run Suitability model
  suitmodels[[i]] <- suitmodel(form = paste(colnames(data), collapse = "+"),
                               lu = landuse,
                               data = data,
                               sub = sample(1:nrow(covs), ssize),
                               resolution = 10000,
                               maxit = 1000,
                               model = TRUE)
}
suitmodels[[1]]$terms

str(suitmodels[[1]])
#-------------------------------#
#### 3. Land use allocations ####
#-------------------------------#

# 3. a prepare data and initialise variables
suitmodels <- readRDS(file.path(out_path, "suitmodels.rds"))

# we have three scenarios we need to run land use simulations for: rcp45 (SSP2 baseline), rcp85 (SSP5 baseline) and dluh
# rcp45: this is mapped to SSP2
# rcp85: this is mapped to SSP5
# dluh: scenario that uses the LUH1 demands

# rcp 45: projecting 
scens <- c("rcp45", "rcp85",  "dluh")

# names of land use classes
lu_names <- c("cropland", "pasture", "primary", "secondary", "urban")

# get demands calculated from GTAP INT predictions
d1 <- readRDS(file.path(out_path, "demand_RCP45.rds"))
d2 <- readRDS(file.path(out_path, "demand_RCP85.rds"))
d1 <- d1[,match(country_code$GTAP_aggregation, colnames(d1))]
d2 <- d2[,match(country_code$GTAP_aggregation, colnames(d2))]
demand <- list(d1, d2)

# get demands estimated from LUH predictions for validation scenario rluh
dmd_ch <- readRDS(file.path(out_path, "demand_rluh.rds"))

# make some lists to store the final demands. Need to do this because demand for pri and sec is only caclulated during the allocation process
demands <- demands_byts <- final_demands <- list()
inds_global <- which(!is.na(mask[]))
# load parameters for our model constraint
lu_newestablishment <- readRDS(file.path("output", "lu_newestablishment.rds"))
lu_newestablishment <- lu_newestablishment[match(lu_newestablishment[,1], country_code$GTAP_code),]
ts <- 1:8

# save.image("ch3_workspace.RData")

# load("ch3_workspace.RData")

# 3. b allocate

# initialise output data frames
covs_global <- as.matrix(covs_suitmodel)

lu_out <- lu_scen <- list()

for(k in 1:length(scens)){
  
  # initialise neighbourhood and land use variables in the beginning of each scenario
  neigh_ts <- as.matrix(neigh_global)
  lu_ts <- as.matrix(lu_global)
  incr_global <- as.matrix(stack(covs_incr[k]))
  lu_global_ts <- matrix(NA, nrow = ncell(lu_global), ncol = 5)
  colnames(lu_global_ts) <- names(lu_global)
  
  for(j in ts){
    
    # create empty matrix to store time step 
    
    # iterate through regions
    
    for(i in 1:30){
      
      # current region name and code
      regions <- country_code$GTAP_aggregation[i]
      map_code <- country_code$GTAP_code[i]
      inds <- which(gtap_borders[]==map_code)
      
      mod <- suitmodels[[i]]
      preds <- mod$coefnames[-1]
      
      # subset land use to current region
      landuse <- lu_ts[inds,]
      colnames(landuse) <- unlist(lapply(strsplit(colnames(landuse), "_"), FUN = function(x) x[[1]]))
      
      # subset covs and neighbourhood data to current region
      covs <- covs_global[inds,]
      neigh <- neigh_ts[inds,]
      data <- cbind(covs, neigh)
      
      colMeans(landuse)
      # add timestep-specific increment to model covariates that change into future 
      incr <- incr_global[inds, which(colnames(incr_global)%in%preds)]
      changing_covs <- which(colnames(data)%in%colnames(incr))
      data[,changing_covs] <- data[,changing_covs] + incr * j
      
      # predict suitability models
      
      suitpred <- predict(mod, newdata = data,  type = "probs")
      colnames(suitpred) <- colnames(landuse)
      
      # specify parameters for model
      params <- list(
        max_dev = 0.1, # how much can the prediction deviate from the demand target (%)
        resolution = 10000, # the reoslution at which the multinom model works
        growth = lu_newestablishment[i,-1] * ((2100-2020)/length(ts)), # the model constraint
        no_change = NULL, # land use classes not allowed to decrease (urban, becuase of high inittial investment)
        max_iter = 200
      )
      
      # just making sure that the constraint parameters are correctly matched with the region
      if(map_code != lu_newestablishment[i,1]) {stop("constraint parameters matched to the wrong region")}
      
      # Estimate demand for current scenario, region and time step
      # demand change for 10 years (1/8 of demand change 2020-2100)
      if (scens[k]%in%c("rcp45", "rcp85")){
      dmd <- demand[[k]][,regions]/length(ts)
      
      # make sure the land use classes and demand values match correctly
      dmd <- dmd[match(colnames(landuse), tolower(names(dmd)))]
      
      # calculate demand for primary and secondary types (weighting by relative suitability)
      rel_suits <- colMeans(suitpred[,c("primary", "secondary")])
      rel_suits <- rel_suits/sum(rel_suits)
      
      # absolute demand for the current time step
      dmd_abs <- colMeans(landuse) * (1+dmd/100)
      dmd_abs[c("primary", "secondary")] <- (1- sum(dmd_abs, na.rm = TRUE)) * rel_suits
      dmd_abs <- rbind(colMeans(landuse), dmd_abs)
      }
      
      if(scens[k] == "dluh"){
        dmd_abs <- dmd_ch[map_code,]
        dmd_abs <- rbind(colMeans(landuse), colMeans(landuse) + dmd_abs)
      }
      
      if(any(colSums(landuse) == 0)){
        params$no_change <- c(params$no_change, which(colSums(landuse) == 0))
      }
      
      
      # allocate
      lu_local <- NA
      while(!is.matrix(lu_local)){
        lu_local <- allocation(
          lu = landuse, 
          sm = suitpred, 
          ln = neigh,
          dmd = dmd_abs, 
          params = params, 
          constraint = TRUE,
          pa = NULL)
      }
      
      # store matrix
      lu_global_ts[inds,] <- as.matrix(lu_local)
      final_demands[[i]] <- dmd_abs
      print(i)
    }
    
    # Calculate new global neighbourhood rasters for the next time step
    neigh_ts <- as.matrix(neighbourhood(lu_global_ts[inds_global,], weights = weights, mask = mask, suffix = "neigh", cols = 1:5, format = "stack"))
    
    # update time step variable (lu_ts is input land use for next time step)
    lu_ts <- lu_global_ts
    
    # store predicted maps without NA to save space
    lu_scen[[j]] <- lu_global_ts[inds_global,]
    
    # store final demands, because primary and secondary land demand get calculated in each time step
    demands_byts[[j]] <- final_demands
  }
  lu_out[[k]] <- lu_scen
  demands[[k]] <- demands_byts
}

# write land use predictions into rasters
for(i in 1:3){
  lustack <- stack(apply(lu_out[[i]][[8]], 2, FUN = function(x) {r <- mask; r[inds_global] <- x; r}))
  writeRaster(lustack, filename = file.path(processed_path, paste0(scens[i], "_", names(lustack), ".tif")), format = "GTiff", bylayer = TRUE, overwrite = TRUE)
}

# save land use predictions and demands to disk
saveRDS(lu_out, file.path(out_path, "lu_simulations.rds"))
saveRDS(demands, file.path(out_path, "final_demands.rds"))


#--------------------------------------------------#
#### 4. Prediction of fine res lu intensity map ####
#--------------------------------------------------#
# needs to be run on HPC
# 4. a Make prediction data frame with high res data

# Load 0.5 arcmin mask and extract indices
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))
inds_5min <- which(!is.na(getValues(mask_5min)))
mask_30min <- raster(file.path(processed_path, "mask_30min.tif"))
inds_30min <- which(!is.na(getValues(mask_30min)))

# Since there are 2E10^8 data points, we have to predict the models in digestable chunks to avoid memory issues. 
# Create chunks
length(inds_5min)
inds_5min_breaks <- chunk(inds_5min, n.chunks = 10)
inds_30min_breaks <- chunk(inds_30min, n.chunks = 10)

# Load fine res data
files_pred <- list()

files_5min <- list.files(processed_path, pattern = "5min", full.names = TRUE)
files_5min <- files_5min[which(grepl(paste(c("lu", "pop", "unsubregions") , collapse = "|"), files_5min))]
files_30min <- list.files(processed_path, pattern = "30min", full.names = TRUE)
files_30min <- files_30min[which(grepl(paste(c("lu", "pop", "unsubregions") , collapse = "|"), files_30min))]

files_pred[[1]] <- files_5min[-which(grepl(paste(c("ssp2", "ssp5", "rcp45", "rcp85", "dluh") , collapse = "|"), files_5min))]
files_pred[[2]] <- files_5min[which(grepl(paste(c("rcp45", "ssp2", "unsubregions") , collapse = "|"), files_5min))]
files_pred[[3]] <- files_5min[which(grepl(paste(c("rcp85", "ssp5", "unsubregions") , collapse = "|"), files_5min))]
files_pred[[4]] <- files_5min[which(grepl(paste(c("dluh", "ssp5", "unsubregions") , collapse = "|"), files_5min))]
files_pred[[5]] <- files_30min[which(grepl(paste(c("repr", "ssp5", "unsubregions") , collapse = "|"), files_30min))]

# Write very large files into the smaller chunks so we can read them in from disk step by step

scens <- c("pres", "rcp45", "rcp85", "dluh", "repr") # need to have same order as above
landuses <- c("cropland", "pasture", "primary", "secondary", "urban")
intensities <- c("minimal", "light", "intense")
result <- readRDS(file.path(data_path, "lumodels.rds"))
logfile <- file.path(data_path, paste0("gam_log.txt"))

writeLines(c(""), logfile)
cat(paste("scen", scens[k], "\n"), file = logfile, append = T)
k <- 5
for(k in 1:length(scens)){
  
  files_5min <- files_pred[[k]]
  
  for(i in 1:length(files_5min)){
    r <- raster(files_5min[[i]])
    nm <- names(r)
    for(j in 1:length(inds_5min_breaks)){
      if(k%in%1:4){
        saveRDS(r[inds_5min_breaks[[j]]], file = file.path(temp_path, paste0(nm, "_chunk", j, "_nona.rds")))
      }
      if(k == 5){
        saveRDS(r[inds_30min_breaks[[j]]], file = file.path(temp_path, paste0(nm, "_chunk", j, "_nona.rds")))
      }
      print(paste(i, j))
    }
    rm(r)
    removeTmpFiles(h = 0)
  }
  
  
  # Load models
  # 4. b Predict, processing chunks on parallel cores
  
  
  # create cluster
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  # system("ps")
  # system("pkill -f R")
  
  lu_int_names <- as.vector(t(outer(landuses, intensities, paste, sep="_")))[-c(4,14)]
  
  foreach(j = 1:length(inds_5min_breaks), .packages = c("mgcv", "foreach")) %dopar% {
    
    # Load covariates of current chunk into data frame
    gam_data <- lapply(list.files(temp_path, pattern = paste0("chunk", j, "_nona.rds"), full.names = TRUE), FUN = function(x) {readRDS(x)})
    gam_data <- as.data.frame(do.call("cbind", gam_data))
    
    if(k%in%c(1)) { 
      colnames(gam_data) <- c("base_pop_30min", "cropland", "pasture", "primary", "secondary", "unsubregions_30min", "urban")
    }
    
    if(k%in%c(2,3,4,5)) { 
      colnames(gam_data) <- c("cropland", "pasture", "primary", "secondary", "urban", "base_pop_30min", "unsubregions_30min")
    }
    chunk_inds <- inds_5min_breaks[[j]]
    
    gam_data$unsubregions_30min[gam_data$unsubregions_30min%in%c(57, 54, 61,29,155)] <- 99
    gam_data$unsubregions_30min <- as.factor(gam_data$unsubregions_30min)
    
    # Predict the models (as in validation)
    
    results <- foreach(i = 1:5, .packages = c("mgcv")) %do% {
      
      cat(paste("chunk", j," ", landuses[i], "\n"), file = logfile, append = T)
      
      # Load estimated model objects
      m_minimal_or_not <- result[[i]][[1]]
      m_light_or_intense <- result[[i]][[2]]
      colnames(gam_data)
      
      # get predictions from first model (probability of being low)
      minimal <- predict(m_minimal_or_not, type = "response", newdata = gam_data)
      not_minimal <- 1 - minimal
      
      # and second model (probability of being medium *if not low*)
      if(length(m_light_or_intense) != 1){
        light_if_not_minimal <- predict(m_light_or_intense, type = "response", newdata = gam_data)
        intense_if_not_minimal <- 1 - light_if_not_minimal
        light <- light_if_not_minimal * not_minimal
        intense <- intense_if_not_minimal * not_minimal
        predicted <- cbind(minimal, light, intense)
      } else {
        predicted <- cbind(minimal, not_minimal)
      }
      
      # Apply modelled intensity probabilities to land use type fractions
      predicted_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
      colnames(predicted_out) <- paste0(landuses[i], "_", colnames(predicted_out))
      predicted_out
    }
    
    predicted_out <- as.data.frame(do.call('cbind', results))
    colnames(predicted_out) <- lu_int_names
    
    # Save each predicted chunk indiividually to avoid memory issues, we can load them back in and merge into map later
    saveRDS(predicted_out, file = file.path(out_path, paste0("predicted_chunk", j, "_nona.rds")))
  }
  unlink(list.files(temp_path, full.names = TRUE))
  stopCluster(cl)
  
  # Combine predicted chunks to single maps, one for each scen and lnad use intensity type
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  logfile <- file.path(data_path, paste0("log.txt"))
  writeLines(c(""), logfile)
  j <- 1
  foreach(j = 1:length(lu_int_names)) %do% {
    out <- numeric()
    lu_int <- lu_int_names[j]
    
    for(i in 1:10){
      cat(paste(paste(j, i), "\n"), file = logfile, append = T)
      predicted_chunk <- readRDS(list.files(out_path, pattern = paste0("predicted_chunk",i, "_nona.rds"), full.names = TRUE))
      out <- c(out, predicted_chunk[,which(colnames(predicted_chunk) == lu_int)])
    }
    
    if(scens[k] != "repr"){
      out_ras <- mask_5min
      out_ras[inds_5min] <- out
    }
    
    if(scens[k] == "repr"){
      out_ras <- mask_30min
      out_ras[inds_30min] <- out
    }
    writeRaster(out_ras, file.path(out_path, paste0(scens[k], "_", lu_int, "_5min", ".tif")), format = "GTiff", overwrite = TRUE)
  }
  unlink(list.files(out_path, pattern = paste0("predicted_chunk"), full.names = TRUE))
  stopCluster(cl)
}
