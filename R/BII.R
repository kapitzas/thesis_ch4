rm(list = ls())

# from https://adrianadepalma.github.io/BII_tutorial/bii_example.html

#----------------------------#
#### 0. Packages, folders ####
#----------------------------#

library(dplyr) # for easy data manipulation
library(tidyr) # ditto
library(magrittr) # for piping
library(lme4) # for mixed effects models
library(car) # for logit transformation with adjustment
library(raster) # for working with raster data
library(geosphere) # calculating geographic distance between sites
library(foreach) # running loops
library(doParallel)
require(tidyverse)
source("./R/99_helperfunctions.R")

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
data_path <- file.path(".", "data")
temp_path <- file.path(data_path, "temp")
out_path <- file.path("output")
int_path <- file.path(data_path, "lu_intensity")
processed_path <- file.path(raw_path, "Global", "processed rasters")

#----------------------------------#
#### 1. Clean Biodiversity Data ####
#----------------------------------#

# 1. a Load and filter

diversity <- readRDS(file.path(data_path, "PREDICTS", "database.rds"))

# Filter out samples with use land intensity and predominant land use
diversity <- diversity %>%
  filter(Use_intensity != "Cannot decide") %>%
  filter(Predominant_land_use != "Cannot decide")


diversity$LandUse <- NA

# 1. b combine land use intensity classes to new variable

# Primary vegetation
diversity$LandUse[which(diversity$Predominant_land_use == "Primary vegetation" & diversity$Use_intensity == "Minimal use")] <- "Primary minimal"
diversity$LandUse[which(diversity$Predominant_land_use == "Primary vegetation" & diversity$Use_intensity == "Light use")] <- "Primary light"
diversity$LandUse[which(diversity$Predominant_land_use == "Primary vegetation" & diversity$Use_intensity == "Intense use")] <- "Primary intense"

# Secondary vegetation
diversity$LandUse[which(grepl("secondary", tolower(diversity$Predominant_land_use)))] <- "Secondary vegetation"
diversity$LandUse[which(diversity$LandUse == "Secondary vegetation" & diversity$Use_intensity == "Minimal use")] <- "Secondary minimal"
diversity$LandUse[which(diversity$LandUse == "Secondary vegetation" & diversity$Use_intensity == "Light use")] <- "Secondary light"
diversity$LandUse[which(diversity$LandUse == "Secondary vegetation" & diversity$Use_intensity == "Intense use")] <- "Secondary intense"

# Cropland
diversity$LandUse[which(diversity$Predominant_land_use == "Cropland" & diversity$Use_intensity == "Minimal use")] <- "Cropland minimal"
diversity$LandUse[which(diversity$Predominant_land_use == "Cropland" & diversity$Use_intensity == "Light use")] <- "Cropland light"
diversity$LandUse[which(diversity$Predominant_land_use == "Cropland" & diversity$Use_intensity == "Intense use")] <- "Cropland intense"

# Pasture
diversity$LandUse[which(diversity$Predominant_land_use == "Pasture" & (diversity$Use_intensity == "Light use" | diversity$Use_intensity == "Minimal use" ))] <- "Pasture light"
diversity$LandUse[which(diversity$Predominant_land_use == "Pasture" & diversity$Use_intensity == "Intense use")] <- "Pasture intense"

# Urban
diversity$LandUse[which(diversity$Predominant_land_use == "Urban" & (diversity$Use_intensity == "Intense use" | diversity$Use_intensity == "Light use" ))] <- "Urban intense"
diversity$LandUse[which(diversity$Predominant_land_use == "Urban" & diversity$Use_intensity == "Minimal use")] <- "Urban minimal"

# Plantatation - we have no map to inform this, so it's going to be secondary intense
diversity$LandUse[which(diversity$Predominant_land_use == "Plantation forest")] <- "Secondary intense"

diversity <- diversity %>%  mutate(
  LandUse = factor(LandUse),
  LandUse = relevel(LandUse, ref = "Primary minimal")
)

table(diversity$LandUse, diversity$Predominant_land_use)

#---------------------------------------------#
#### 2. Calculate abundance and similarity ####
#---------------------------------------------#

# code adopted from https://adrianadepalma.github.io/BII_tutorial/bii_example.html

# 2. a Total abundance
abundance_data <- diversity %>%
  
  # pull out just the abundance measures
  filter(Diversity_metric_type == "Abundance") %>%
  
  # group by SSBS (each unique value corresponds to a unique site)
  group_by(SSBS) %>%
  
  # now add up all the abundance measurements within each site
  mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  
  # ungroup
  ungroup() %>%
  
  # pull out unique sites
  distinct(SSBS, .keep_all = TRUE) %>%
  
  # now group by Study ID
  group_by(SS) %>%
  
  # pull out the maximum abundance for each study
  mutate(MaxAbundance = max(TotalAbundance)) %>%
  
  # ungroup
  ungroup() %>%
  
  # now rescale total abundance, so that within each study, abundance varies from 0 to 1.
  mutate(RescaledAbundance = TotalAbundance/MaxAbundance)


# 2. b Compositional Similarity (asymmetric Jaccard Index)

# Prepare data for compositional similarity estimation
cd_data_input <- diversity %>%
  
  # pull out only the abundance data
  filter(Diversity_metric_type == "Abundance") %>%
  
  # group by Study
  group_by(SS) %>%
  
  # calculate the number of unique sampling efforts within that study
  mutate(n_sample_effort = n_distinct(Sampling_effort)) %>%
  
  # calculate the number of unique species sampled in that study
  mutate(n_species = n_distinct(Taxon_name_entered)) %>%
  
  # check if there are any Primary minimal sites in the dataset
  mutate(n_primin_records = sum(LandUse == "Primary minimal")) %>%
  
  # ungroup
  ungroup() %>%
  
  # now keep only the studies with one unique sampling effort
  filter(n_sample_effort == 1) %>%
  
  # and keep only studies with more than one species 
  # as these studies clearly aren't looking at assemblage-level diversity
  filter(n_species > 1) %>%
  
  # and keep only studies with at least some Primary minimal data
  filter(n_primin_records > 0) %>%
  
  # drop empty factor levels
  droplevels()


# Calcuate dissimilarity in parallel

studies <- distinct(cd_data_input, SS) %>%
  pull()

registerDoParallel(cores = 2)

logfile <- file.path("dissimilarity_log.txt")
writeLines(c(""), logfile)

cd_data <- foreach(s = studies,
                   .combine = rbind,
                   .packages = c("dplyr", "magrittr", "geosphere")) %dopar% {
                     cat(paste(which(studies%in%s), "/", length(studies), "\n"), file = logfile, append = T)

                     # filter out the given study
                     data_ss <- filter(cd_data_input, SS == s)

                     # pull out the SSBS and LandUse information (we'll use this later to assign a land use contrast to each pair of site
                     site_data <- data_ss %>%
                       dplyr::select(SSBS, LandUse) %>%
                       distinct(SSBS, .keep_all = TRUE)

                     # pull out the sites that are Primary minimal (we only want to use comparisons with the baseline)
                     baseline_sites <- site_data %>%
                       filter(LandUse == "Primary minimal") %>%
                       pull(SSBS)

                     # pull out all the sites
                     site_list <- site_data %>%
                       pull(SSBS)


                     # get all site x site comparisons for this study
                     site_comparisons <- expand.grid(baseline_sites, site_list) %>%

                       # rename the columns so they will be what the compositional similarity function expects for ease
                       rename(s1 = Var1, s2 = Var2) %>%

                       # remove the comparisons where the same site is being compared to itself
                       filter(s1 != s2)


                     # apply the compositional similarity function over each site combination in the dataset
                     sor <- apply(site_comparisons, 1, function(y) getJacAbSym(data = data_ss, s1 = y['s1'], s2 = y['s2']))

                     # calculate the geographic distance between sites
                     # first pull out the lat and longs for each site combination
                     s1LatLong <- as.matrix(data_ss[match(site_comparisons$s1, data_ss$SSBS), c('Longitude','Latitude')])
                     s2LatLong <- as.matrix(data_ss[match(site_comparisons$s2, data_ss$SSBS), c('Longitude','Latitude')])

                     # then calculate the distance between sites
                     dist <- distHaversine(s1LatLong, s2LatLong)

                     # pull out the land-use contrast for those site combinations
                     Contrast <- paste(site_data$LandUse[match(site_comparisons$s1, site_data$SSBS)],
                                       site_data$LandUse[match(site_comparisons$s2, site_data$SSBS)],
                                       sep = "-")

                     # put all the information into a single dataframe

                     study_results <- data.frame(site_comparisons,
                                                 sor,
                                                 dist,
                                                 Contrast,
                                                 SS = s,
                                                 stringsAsFactors = TRUE)



                   }

registerDoSEQ()

cd_data <- cd_data %>%

  # Firstly, we only care about comparisons where Primary minimal is the first site
  # so pull the contrast apart
  separate(Contrast, c("s1_LandUse", "s2_LandUse"), sep = "-", remove = FALSE) %>%

  # filter sites where s1_LandUse is the baseline site
  filter(s1_LandUse == "Primary minimal") %>%

  # logit transform the compositional similarity
  mutate(logitCS = logit(sor, adjust = 0.001, percents = FALSE)) %>%

  # log10 transform the geographic distance between sites
  mutate(log10geo = log10(dist + 1)) %>%

  # make primary minimal-primary minimal the baseline again
  mutate(Contrast = factor(Contrast),
         Contrast = relevel(Contrast, ref = "Primary minimal-Primary minimal"))

saveRDS(cd_data, file.path("output", "cd_data.rds"))

cd_data <- readRDS(file.path(out_path, "cd_data.rds"))

#-----------------------#
#### 3. Build models ####
#-----------------------#

# 3.a Model compositional similarity as a function of the land-use contrast and the geographic distance between sites

cd_m <- lmer(logitCS ~ Contrast + log10geo + (1|SS), data = cd_data)

saveRDS(cd_m, file.path(out_path, "cd_model.rds"))

# 3. b Model abundance as a function of the land-use level
ab_m <- lmer(sqrt(RescaledAbundance) ~ LandUse + (1|SS), data = abundance_data)

saveRDS(ab_m, file.path(out_path, "ab_model.rds"))


#-------------------------#
#### 4. Predict models ####
#-------------------------#

# 4.a predict abundance
# set up data frame with lan-use levels
newdata_ab <- data.frame(LandUse = levels(abundance_data$LandUse)) %>%
  
  mutate(ab_m_preds = predict(ab_m, ., re.form = NA) ^ 2)

# 4. b predict compositional similarity
# set up data frame with land-use contrasts
newdata_cd <- data.frame(Contrast = levels(cd_data$Contrast),
                         log10geo = 0) %>%
  mutate(cd_m_preds = predict(cd_m, ., re.form = NA) %>%
           inv_logit(a = 0.001))

# 4.c Calculate mean comp similarity and abundance for each class/contrast (for figure)
dat <- 
  newdata_cd %>% 
  mutate(Contrast = stringr::str_remove(pattern = "Primary minimal-", string = Contrast)) %>% 
  left_join(newdata_ab, by = c("Contrast" = "LandUse")) %>% 
  select(-log10geo) %>% 
  column_to_rownames('Contrast') %>% 
  apply(2, function(x) {x[-1]/x[1]}) %>% 
  data.frame %>% 
  mutate("BII" = cd_m_preds * ab_m_preds)

write.csv(dat, file = file.path(out_path, "cd_ab_BII.csv"))

# 4.d Predict globally at 5min res
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))
inds_sub_5min <- which(!is.na(mask_5min[]))
mask_30min <- raster(file.path(processed_path, "mask_30min.tif"))
inds_sub_30min <- which(!is.na(mask_30min[]))

unsub_5min <- raster(list.files(processed_path, pattern = "unsubregions_5min", full.names = TRUE))
unsub_30min <- raster(list.files(processed_path, pattern = "unsubregions_30min", full.names = TRUE))

scens <- c("pres", "rcp45", "rcp85", "dluh", "repr")
for (k in 1:length(scens)){
  
  if(k%in%c(1:4)){
  mask_sub <- mask_5min
  inds_sub <- which(!is.na(mask_5min[]))
  }
  if(k==5){
    mask_sub <- mask_30min
    inds_sub <- which(!is.na(mask_30min[]))
  }
  
  files_5min <- list.files(processed_path, full.names = TRUE)
  files_5min <- files_5min[grepl(paste(c("light", "intense", "minimal"), collapse = "|"), files_5min)]
  files_predict <- files_5min[grepl(scens[k], files_5min)]
  lus <- as.data.frame(stack(files_predict))[inds_sub,]
  colnames(lus) <- colnames(lus) %>% str_replace("_5min", "") %>% str_replace("pres_", "")
  
  ras_list <- list()
  i <- 1
  for(i in 1:ncol(lus)){
    ras <- mask_sub
    ras[inds_sub] <- lus[,i]
    ras_list[[i]] <- ras
  }
  
  ras <- stack(ras_list)
  names(ras) <- colnames(lus)
  newdata_ab <- newdata_ab[order(newdata_ab$LandUse),]
  
  ab_raster <- ras[[1]]
  ab_raster[inds_sub] <- 0
  
  for(i in 1:nrow(newdata_ab)){
    ab_raster <- ab_raster + (newdata_ab[i,2] * ras[[i]])
  }
  
  # devide by reference value
  ab_raster <- ab_raster /  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary minimal']
  # plot(ab_raster)
  
  newdata_cd <- newdata_cd[order(newdata_cd$Contrast),]
  
  cd_raster <- ras[[1]]
  cd_raster[inds_sub] <- 0
  
  for(i in 1:nrow(newdata_ab)){
    cd_raster <- cd_raster + (newdata_cd[i,3] * ras[[i]])
  }
  
  cd_raster <- cd_raster / newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary minimal']
  
  bii <- ab_raster * cd_raster
  writeRaster(bii, filename = file.path(out_path, paste0("bii_", scens[k], ".tif")), format = "GTiff", overwrite = TRUE)
  writeRaster(cd_raster, filename = file.path(out_path, paste0("cd_", scens[k], ".tif")), format = "GTiff", overwrite = TRUE)
  writeRaster(ab_raster, filename = file.path(out_path, paste0("ab_", scens[k], ".tif")), format = "GTiff", overwrite = TRUE)
  print(scens[k])
}
