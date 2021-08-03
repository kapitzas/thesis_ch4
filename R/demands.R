#----------------------------#
#### 0. Packages, folders ####
#----------------------------#

rm(list = ls())
require(rgdal)
require(raster)
require(raster)

data_path <- file.path(".", "data")
temp_path <- file.path("/Volumes", "external", "c3 processing", "temp")
out_path <- file.path(".", "output")
demand_path <- file.path(data_path, "demand calculation")
raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")

#-----------------------------------------------------#
#### 1. Match UN countries to our GTAP aggregation ####
#-----------------------------------------------------#
# Load data
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))

gtap_regions <- read.csv(file.path(data_path, "GTAP_regions.csv"))
world_borders <- readOGR(file.path(raw_path, "Global", "TM_WORLD_BORDERS_SIMPL-0", "TM_WORLD_BORDERS_SIMPL-0.3.shp"))
regions_table <- world_borders@data
regions_table$ISO3 <- tolower(regions_table$ISO3)


# match countries to our GTAP aggregation, depending on available info

matching_list<- list(
  "by_ISO2" = list(
    "e25" = c("AT", "EE", "IT", "PT", "BE", "FI", "LV", "FR", "LT", "SK", "DE", "LU", 'SI', 'CY', 'EL', 'MT', 'ES', 'CZ', 'HU', 'NL', 'SE', 'DK', 'IE', 'PL'),
    "gbr" = "UK"
  ),
  
  "by_UN" = list(
    "xsu" = c(031, 051, 112, 233, 268, 398, 417, 428, 440, 498, 762, 795, 804, 860)
  ),
  
  "by_subregion" = list(
    "xoc" = c(54,57,61),
    "xea" = 30,
    "xse" = 35,
    "xsa" = 34,
    "xsm" = c(5, 13, 29),
    'xer' = c(39, 151, 154, 830, 155),
    'xws' = 145,
    'xnf' = 15,
    'ssa' = c(11, 14, 17, 18)
  )
)

regions_table['GTAP'] <- NA
gtap_mapping <- read.csv(file.path(demand_path, "GTAP_countrymappingraw.csv"))[1:30,1:2]
gtap_mapping$Code.30 <- trimws(gtap_mapping$Code.30)
# Allocate countries we retained in GTAP aggregation

existing <- which(regions_table$ISO3%in%gtap_mapping$Code.30)
regions_table$GTAP[existing] <- regions_table$ISO3[existing]

# Match EU25 countries (ISO2 info)
matches <- matching_list[[1]]
i <- 1
for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$ISO2%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match former soviet union countries (UN code info)
matches <- matching_list[[2]]

for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$UN%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match other countries (UN SUBREGION code info)
matches <- matching_list[[3]]
for(i in 1:length(matches)){
  regions_table$GTAP[which(regions_table$SUBREGION%in%matches[[i]] & is.na(regions_table$GTAP))] <- names(matches[i])
}

# Match the unallocated rest
regions_table$GTAP[is.na(regions_table$GTAP)] <- "xtw" # This includes Taiwan.

# add our own GTAP code column (so we can make a numeric world raster with our regions)
regions_table$GTAP_code <- gtap_regions$GTAP_code[match(regions_table$GTAP, gtap_regions$GTAP_aggregation)]
world_borders@data <- regions_table

# Rasterize
writeOGR(world_borders, dsn = file.path(temp_path), layer = "GTAP_aggregation_borders", driver = "ESRI Shapefile", overwrite_layer = TRUE)
output <- file.path(temp_path, "gtap_aggregation_5min.tif")
input <- file.path(temp_path, "GTAP_aggregation_borders.shp")
gdaltools::rasterize_shp(input, output, res = 0.083, c(-180, 180, -90, 90), attribute = "GTAP_code")

#--------------------------#
#### 2. land use demand ####
#--------------------------#

# we need to match up the various data sources first
# a) FAO: harvested area in each FAO commodity in 2019 by country. Need to match the ~ 200 FAO crop types to the 8 GTAP groups, and the FAO countries to our GTAP aggregation. Result is area harvested in each GTAP commodity in 2019 wihtin each of our GTAP regions.
# b) WORLDBANK: urban pop growth rate. Future data (until 2050) is only available for worldbank aggregations of the world (i.e. low income coutnries, south asia, etc.). These groups don't overlap with our GTAP regions, so we need to make an educated guess. Result is the annual growth rate of the urban population for each of our GTAP aggregations.

# 2. a 2019 harvested area in each GTAP commodity and region.
# Load downloaded FAO raw data for 2019, GTAP classes and GTAP regions look-up tables
fao_crops <- read.csv(file.path(demand_path, "FAOSTAT_harvestedraw.csv"))
gtap_classes <- read.csv(file.path(demand_path, "GTAP_classes.csv"))
areas <- read.csv(file.path(data_path, "GTAP_regions.csv"))

# Match GTAP classes to FAO crop types
fao_crops$GTAP <- NA
fao_crops$GTAP <- gtap_classes$GTAP.sector[match(fao_crops$Item, unique(gtap_classes$Item))]
fao_crops$GTAP_code <- areas$GTAP_code[match(fao_crops$Area, areas$fao_areas)]

# Aggregate area harvested by commodity and region
crops_by_region <- aggregate(Value ~ GTAP + GTAP_code, data = fao_crops, FUN = "sum")
crops_by_region$GTAP_country <- areas$GTAP_aggregation[match(crops_by_region$GTAP_code, areas$GTAP_code)]
write.csv(crops_by_region, file.path(demand_path, "GTAP_harvestedbyregion.csv"), row.names = TRUE)


scens <- c("RCP45", "RCP85")
# Estimate demand trajectory 

for(j in 1:length(scens)){
  crops_by_region <- read.csv(file.path(demand_path, "GTAP_harvestedbyregion.csv"))
  
  gtap_traj <- read.csv(file.path(demand_path, "GTAP trajectories", paste0("qfe_", scens[j], "_2100.csv")))
  gtap_traj <- gtap_traj[c(which(gtap_traj$qfe%in%unique(crops_by_region$GTAP)),9), ]
  gtap_regions <- unique(crops_by_region$GTAP_country)
  
  # cropland and pasture
  colnames(crops_by_region)
  change_by_region <- matrix(NA, ncol = 30, nrow = 5)
  colnames(change_by_region) <- unique(gtap_regions)
  rownames(change_by_region) <- c("Cropland", "Pasture", "Urban", "Primary", "Secondary")
  
  for(i in 1:length(unique(gtap_regions))){
    last_ts <- gtap_traj[1:8, c("qfe", gtap_regions[i])]
    weights <- crops_by_region[which(crops_by_region$GTAP_country%in%gtap_regions[i]), c("Value", "GTAP")]
    missing <- last_ts$qfe[which(!last_ts$qfe%in%weights$GTAP)]
    if(length(missing) > 0){
      weights <- rbind(weights, data.frame("Value" = rep(mean(weights$Value), length(missing)), "GTAP" = missing))
    }
    change_by_region[1,i] <- weighted.mean(last_ts[,2], w = as.numeric(weights[match(last_ts[,1], weights$GTAP), 1]))
    change_by_region[2,i] <- gtap_traj[9, c(gtap_regions[i])]
  }
  
  # Urban land
  # global_popgrowth <- read.csv(file.path(demand_path, "WORLDBANK_urbanpopraw.csv"))
  # global_popgrowth <- global_popgrowth[which(global_popgrowth$Indicator.Code == "SP.URB.GROW"),]
  # global_popgrowth3 <- na.omit(cbind(global_popgrowth[,1:2], rowMeans(global_popgrowth[,paste0("X", 2020:2050)])))
  
  # now we copy the values over manually to WORLDBANK_urbanpopbyregion, because they don't really match up and we have to make educated guesses.
  
  urban_popgrowth <- read.csv(file.path(demand_path,  "WORLDBANK_urbanpopbyregion.csv"))
  
  
  popgrowth <- (1 + urban_popgrowth$pop.growth/100)^80
  change_by_region[3,] <- popgrowth[match(colnames(change_by_region), trimws(urban_popgrowth$Code.30))]
  saveRDS(change_by_region, file.path(out_path, paste0("demand_", scens[j], ".rds")))
}

#---------------------------------------------------------------------#
#### 3. Get demand from HURTT projections for M8.5 demand scenario ####
#---------------------------------------------------------------------#

# get demand from HURTT predictions
lu_global <- stack(list.files(processed_path, pattern = "lu_5min", full.names = TRUE, recursive = TRUE))
lu_global <- lu_global[[-which(grepl("rcp", names(lu_global)))]]
gtap_borders <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
table_gtap <- table(gtap_borders[])

lu_files <- list.files("/Volumes/external/c3 processing/LUHa_u2.v1_message.v1/updated_states", full.names = TRUE)
lu_files <- lu_files[grep(pattern = paste(c("gcrop", "gothr", "gpast", "gsecd", "gurbn"), collapse = "|"), lu_files)]

lu_hurtt <-stack(lu_files[grep(pattern = "2100", lu_files)])
cell_sums <- sum(lu_hurtt)
lu_hurtt[which(cell_sums[] < 0.999)] <- NA
lu_hurtt <- as.data.frame(lu_hurtt)

# Estimate demand by GTAP region RCP85
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))
gtap_aggregation <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
lu <- raster(lu_files[[1]])
crs(lu) <- crs(mask_5min)
mask_05degree <- projectRaster(mask_5min, lu, method = "ngb")
gtap_aggregation_05degree <- projectRaster(gtap_aggregation, mask_05degree, method = "ngb")[]

lu_global_df <- as.data.frame(lu_global)
colnames(lu_global_df) <- unlist(lapply(strsplit(colnames(lu_global_df), "_"), FUN = function(x) x[[1]]))

dmds <- list()
lus <- cbind( "gtap_borders" = gtap_borders[], lu_global_df)
dmds[[1]] <- aggregate(.~ gtap_borders, FUN = mean, data = na.omit(lus))
rowSums(dmds[[1]][,-1])

lus2 <- cbind("gtap_borders" = gtap_aggregation_05degree, lu_hurtt)
dmds[[2]] <- aggregate(.~ gtap_borders, FUN = mean, data = na.omit(lus2))
colnames(dmds[[2]])[2:6] <- c("cropland", "primary", "pasture", "secondary", "urban")
dmds[[2]] <- dmds[[2]][,match(colnames(lu_global_df), tolower(names(dmds[[2]])))]

dmd_ch <- (dmds[[2]] - dmds[[1]][,-1])/8
saveRDS(dmd_ch, file.path(out_path, "demand_rluh.rds"))
