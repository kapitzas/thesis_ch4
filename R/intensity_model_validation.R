# Preproccesing

# This script estimates parameters of the relationship between land use type and intensity at 0.5 degree resolution with population intensity and fractional land use cover (following Newbold et al., 2015). The modelling is done using gam. The parameters are then applied to 0.5 arcmin resolution population density and fractional cover, to obtain a global map of land use type and intensity at that resoltuion. 



rm(list = ls())

#----------------------------#
#### 0. packages, folders ####
#----------------------------#


# Packages
# devtools::install_github("kapitzas/gdaltools") # package to directly link R to rgdal
# devtools::install_github("kapitzas/flutes") # integerify function

require(devtools)
require(WorldClimTiles)
require(rgdal)
require(gdaltools)
require(raster)
require(stringr)
require(rgeos)
require(sf)
require(mgcv)
library(BBmisc)
require(data.table)
library("doParallel")
require(flutes)

# Folders
# setwd("/home/student.unimelb.edu.au/kapitzas/ch3") only on boab
processed_path <- file.path(getwd(), "preprocessed layers")
data_path <- file.path(getwd(), "data")
temp_path <- file.path(data_path, "temp")
out_path <- file.path(data_path, "output")
int_path <- file.path(data_path, "lu_intensity")

#-----------------------#
#### 1. Build models ####
#----------------------## 

# 1. a Create modelling dataframe
# this creates a data frame with the land use + intensity classes estimated from GLS data set as response and un subregion and pop density as covariates at 30 min res.


# Load mask, extract non-NA indices

mask_30min <- raster(file.path(processed_path, "mask_30min.tif"))
inds_30min <- which(!is.na(mask_30min[]))


# Load 30 min data for modelling into data frame
files_30min <- list.files(processed_path, pattern = "30min", full.names = TRUE)
files_30min <- grep(paste(c("predicts", "unsubregions", "base_pop"),collapse="|"), files_30min, value=TRUE)

for(i in 1:length(files_30min)){
  r <- raster(files_30min[i])
  plot(r, main = names(r))
}

gam_data <- as.data.frame(stack(files_30min))[inds_30min,]

landuses <- c("cropland", "pasture", "primary", "secondary", "urban")
intensities <- c("minimal", "light", "intense")

# turn response columns (fractions of land use type and intensity classes) into integers for GAM
# to run GAM, the fractional data have to be in integer representaiton

for(i in 1:length(landuses)){
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  gam_data[landuses[i]] <- rowSums(gam_data[,lu_inds])
  gam_data[,lu_inds] <- integerify(gam_data[,lu_inds], resolution = 1000)
}

# Combine pacific micro states into one region (otherwise there is not enough data to build models on some un subregions)
gam_data[is.na(gam_data)] <- 0
gam_data$unsubregions_30min[gam_data$unsubregions_30min%in%c(57, 54, 61,29,155)] <- 99
gam_data$unsubregions_30min <- as.factor(gam_data$unsubregions_30min) # this is random effect, need to trun into factor

# 1. b Build models

# Prepare CPU cluster
cl <- makeCluster(5)
registerDoParallel(cl)

# Log file to keep track
logfile <- file.path(data_path, paste0("gam_log.txt"))
writeLines(c(""), logfile)

results <- foreach(i = 1:5, .packages = c("mgcv")) %dopar% {
  
  # write to logfile
  cat(paste(landuses[i],"\n"), file = logfile, append = T)
  
  lus <- landuses[i] # land use classes types
  
  # Two different formulas, one for minimal or not, one for light or intense, if not minimal
  f_minimal_or_not <- as.formula(paste("minimal_or_not ~", paste(c("unsubregions_30min", paste0("s(", c("base_pop_30min", lus), ")"), paste0("s(", c(lus, "base_pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  f_light_or_intense <- as.formula(paste("light_or_intense ~", paste(c("unsubregions_30min", paste0("s(", c("base_pop_30min", lus), ")"), paste0("s(", c(lus, "base_pop_30min"), ",by = unsubregions_30min)")), collapse = "+")))
  
  # Get indices of gam_data columns matching land use type and intensity classes
  lu_inds <- grep(paste(paste(landuses[i], intensities, sep = "_"), collapse = "|"), colnames(gam_data))
  
  
  # Modelling probability for minimial or not minimal (first model)
  
  cat(paste("minimal or not: ",landuses[i],"\n"), file = logfile, append = T)
  
  # Response data matrix: When we only two intensities for a lnad use type, get minimal and not minimal (whichever it is).
  if(length(lu_inds) == 2){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], gam_data[,lu_inds][, 2])
  }
  
  # When we have all three intensities for a lnad use type, get minimal and the sum of those that are not not minimal to model against
  if(length(lu_inds) == 3){
    minimal_or_not <- cbind(gam_data[,lu_inds][,1], rowSums(gam_data[,lu_inds][,-1]))
  }
  # run "minimal or not" GAM
  m_minimal_or_not <- tryCatch(mgcv::gam(f_minimal_or_not, family = binomial, data = gam_data), error = function(e) NA)
  if(is.na(m_minimal_or_not)) {cat(paste("returned NA - minimal or not: ",landuses[i],"\n"), file = logfile, append = T)}
  
  # Modelling probability for light or intense, when not minimal (second model)
  
  # only when we have all three classes in a type we need to do this
  if(length(lu_inds)==3){
    cat(paste("light or intense: ",landuses[i],"\n"), file = logfile, append = T)
    
    # Response data matrix: all intensities except minimal
    light_or_intense <- as.matrix(gam_data[,lu_inds][,-1])
    
    # run "light or intense" GAM
    m_light_or_intense <- tryCatch(mgcv::gam(f_light_or_intense, family = binomial, data = gam_data), error = function(e) NA)
    if(is.na(m_light_or_intense)) {cat(paste("returned NA - light or intense: ",landuses[i],"\n"), file = logfile, append = T)}
  } else {
    m_light_or_intense <- NA 
  }
  
  # Return the intensity models for each land use type
  list(m_minimal_or_not, m_light_or_intense, landuses[i])
}

# Save models, so we can load them to predict and do validation below
saveRDS(results, file.path(data_path, "lumodels.rds"))

stopCluster(cl)

#---------------------#
#### 2. Validation ####
#---------------------#

# 2. a Load intensity models
result <- readRDS(file.path(data_path, "lumodels.rds"))
valid_df <- list()

# 2. b Predict intensity for each land use type
for(i in 1:length(landuses)){
  print(i)
  # Load models for current land use type
  m_minimal_or_not <- result[[i]][[1]] # light or intense for pasture
  m_light_or_intense <- result[[i]][[2]]
  
  # predict first model (probability of being minimal)
  minimal <- predict(m_minimal_or_not, type = "response", newdata = gam_data)
  not_minimal <- 1 - minimal
  
  # predict second model (probability of being light, *if not minimal*)
  if(length(m_light_or_intense) != 1){
    light_if_not_minimal <- predict(m_light_or_intense, type = "response", newdata = gam_data)
    intense_if_not_minimal <- 1 - light_if_not_minimal
    light <- light_if_not_minimal * not_minimal
    intense <- intense_if_not_minimal * not_minimal
    predicted <- cbind(minimal, light, intense)
  } else {
    predicted <- cbind(minimal, not_minimal)
  }
  
  # Now we have the probability of being minimal, light or intense, given our model covariates pop dens and fractional land use type (wihtout intnesity)
  
  # multiply each of the land use type columns with the according two or three predicted probability columns to get land use type and intensity map
  valid_out <- sweep(predicted, 1, as.matrix(gam_data[landuses[i]]), FUN = "*")
  colnames(valid_out) <- paste0(landuses[i], "_", colnames(valid_out))
  valid_df[[i]] <- valid_out
}


# combine to data frame
valid_df <- do.call("cbind", valid_df)

# 2. c Map errors

# Get observed data (the data set we created from gls and harmonized lu data)
obs_df <- stack(list.files(processed_path, pattern = "predicts_30min", full.names = TRUE))
obs_df <- as.data.frame(obs_df)[inds_30min,]


files_bla <- list.files(processed_path, full.names= TRUE, pattern = "min.tif")
for (i in 1:length(files_bla)){
  print(length(which(!is.na(raster(files_bla[[i]])[]))))
}

dim(obs_df) == dim(valid_df)


# Caluclate and map RMSE
error <- as.matrix(obs_df) - as.matrix(valid_df)
error <- sqrt(rowMeans(error^2))
error_map <- mask_30min
error_map[] <- NA
error_map[inds_30min] <- error
spplot(error_map)
