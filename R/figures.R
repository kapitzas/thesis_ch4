rm(list = ls())
source("./R/99_helperfunctions.R")
#detach("package:flutes", unload=TRUE)
#devtools::install_github("kapitzas/flutes", force = TRUE)

#-------------------------#
#### Packages, folders ####
#-------------------------#
require(lme4)
require(raster)
require(rasterVis)
require(viridis)
require(scales)
require(effectsize)
require(tidyverse)
require(ggstance)
require(cowplot)
require(flutes)

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")
temp_path <- file.path("/Volumes", "external", "c3 processing", "temp")
out_path <- file.path(".", "output")
data_path <- file.path(".", "data")
thesis_path <- file.path("/Users", "simon", "ownCloud", "PhD", "writing", "thesis")
figure_thesis <- file.path(thesis_path, "chapters", "figures", "chapter4")
figure_chapter <- file.path(".", "figures")

# Load some data we need for the figures
mask_5min <- raster(file.path(processed_path, "mask_5min.tif"))
mask_30min <- raster(file.path(processed_path, "mask_30min.tif"))

ecoregs_30min <- raster(file.path(processed_path, "biorealms_30min.tif"))
bioreg_names <- read.csv(file.path(out_path, "ecoreglevels.csv"),  na.strings = "-9999")

single_width <- 3.25 #(82.5mm, single column)
medium_width <- 4.33 #110mm, 1.5 x page)
full_width <- 6 #165 mm

#---------------------------------#
#### Figure: Map of 30 regions ####
#---------------------------------#

gtap <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
classes <- read.csv("/Users/simon/ownCloud/PhD/chapter3/data/GTAP_regions.csv")
classes <- unique(classes[,c(4,3)])
classes <- classes[order(classes$GTAP_code),]
colnames(classes)[1] <- "ID"
tar <- levels(gtap)[[1]]
tar <- classes
levels(gtap) <- as.data.frame(tar)


l1 <- levelplot(gtap, margin = F, colorkey = list(space='right'),
                xlab="", ylab="",
                par.settings=list(
                  strip.border=list(col='transparent'),
                  strip.background=list(col='transparent'),
                  axis.line=list(col='transparent')),
                scales=list(draw=FALSE),
                col.regions = magma(100),
                maxpixels = 1e6
)


#-------------------------------------------#
#### Figure: Map of biogeographic realms ####
#-------------------------------------------#

bioreg_names$Name <- c("Australasia", "Antarctic", "Afrotropic", "Indo-Malay", "Nearctic", "Neotropical", "Oceania", "Palearctic")
colnames(bioreg_names)[3] <- "ID"

levels(ecoregs_30min) <- bioreg_names[,c(3,4)]

l2 <- levelplot(ecoregs_30min, margin = F,
                colourkey = "bottom",
                xlab="", ylab="",
                par.settings=list(
                  strip.border=list(col='transparent'),
                  strip.background=list(col='transparent'),
                  axis.line=list(col='transparent')),
                scales=list(draw=FALSE),
                col.regions = magma(100),
                maxpixels = 1e6
)

l <- plot_grid(l1, l2, ncol = 1, labels = c("(a)", "(b)"), label_size = 10, label_fontface = "plain")

fig_paths <- c(file.path(figure_thesis, "fig_regionmap.pdf"),
               file.path(figure_chapter, "fig_regionmap.pdf"))
multi_pdf(x = l, paths = fig_paths, height = 7, width = full_width)

#-------------------------------------#
#### Table: Harvested area in 2019 ####
#-------------------------------------#

harvested <- read.csv("/Users/simon/ownCloud/PhD/chapter3/data/demand calculation/GTAP_harvestedbyregion.csv")
regions <- harvested %>% select(GTAP_country) %>% unique()
out <- matrix(nrow = 30, ncol = 9)
colnames(out) <- c("region", "c_b", "gro", "ocr", "osd", "pdr", "pfb", "v_f", "wht")
out <- as.data.frame(out)

for(i in 1:nrow(regions)){
  harv_coun <- 
    harvested %>% 
    filter(GTAP_country == regions[i,]) %>%  
    select(c(GTAP, Value))
  out[i,1] <- regions[i,]
  out[i, match(harv_coun$GTAP, colnames(out))] <- harv_coun$Value
}
out[,-1] <- round(out[,-1]/rowSums(out[,-1], na.rm = TRUE)  * 100, 2)
write.csv(out, "/Users/simon/ownCloud/PhD/writing/thesis/tables/ch4_fao_estimates.csv")

#-------------------------------#
#### Figure: CGE projections ####
#-------------------------------#

qfe85 <- read.csv(file.path(data_path, "demand calculation", "GTAP trajectories", "qfe_RCP85_2100.csv"))[c(1:9), ]
qfe60 <- read.csv(file.path(data_path, "demand calculation", "GTAP trajectories", "qfe_RCP45_2100.csv"))[c(1:9), ]
gtap_classes <- qfe85$qfe
ctl_85 <- qfe85[9,]
ctl_60 <- qfe60[9,]
qfe60 <- qfe60[-9,]
qfe85 <- qfe85[-9,]
qfe85[,-1] <- round(qfe85[,-1], 2)
qfe60[,-1] <- round(qfe60[,-1], 2)

colnames(qfe60)[-1] == t(out)[1,]
qfe60[,1] == rownames(t(out))[-1]

out_t <- t(out)
colnames(out_t) <- out_t[1,]
out_t2 <- apply(out_t[-1,], 2, FUN = function(x) {as.numeric(x)})
rownames(out_t2) <- rownames(out_t)[-1]

nas <- which(is.na(out_t2), arr.ind = TRUE)
nas[,1] <- rownames(out_t2)[nas[,1]]
nas[,2] <- colnames(out_t2)[as.numeric(nas[,2])]

na_rows <- match(nas[,1], qfe60[,1])
na_cols <- match(nas[,2], colnames(qfe60))
nas2 <- cbind(na_rows, na_cols)

for(i in 1:nrow(nas2)){
  qfe60[nas2[i,1], nas2[i,2]] <- NA
  qfe85[nas2[i,1], nas2[i,2]] <- NA
}

qfe85 <- rbind(qfe85, ctl_85)
qfe60 <- rbind(qfe60, ctl_60)
qfe85 <- qfe85[,-1] %>% 
  gather() %>% 
  mutate(classes = rep(gtap_classes, 240/8))
qfe85$rcp <- 'SSP5'

qfe60 <- qfe60[,-1] %>% 
  gather() %>% 
  mutate(classes = rep(gtap_classes, 240/8))
qfe60$rcp <- 'SSP2'

# get total ranges by SSP for text
qfe60 %>% filter(value == min(value, na.rm = TRUE))
qfe60 %>% filter(value == max(value, na.rm = TRUE))

qfe85 %>% filter(value == min(value, na.rm = TRUE))
qfe85 %>% filter(value == max(value, na.rm = TRUE))

# get ranges for sectors:

# wheat
qfe85 %>% filter(classes == "wht") %>% arrange(value)
qfe85 %>% filter(classes == "ctl") %>% arrange(value)
qfe60 %>% filter(classes == "ctl") %>% arrange(value)

qfe <- rbind(qfe85, qfe60)

fig_gtap <- qfe %>% 
  ggplot(aes(y = key, x = classes, fill= value)) + 
  facet_wrap(~rcp, strip.position="bottom") +
  geom_tile() + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position="bottom", 
        legend.key.width= unit(1.2, 'cm')) +
  xlab("GTAP commodity") +
  ylab("GTAP INT region") + 
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue", name = "Change [%]")


fig_paths <- c(file.path(figure_thesis, "fig_gtap.pdf"),
               file.path(figure_chapter, "fig_gtap.pdf"))
multi_pdf(x = fig_gtap, paths = fig_paths, width = medium_width + 0.3, height = 5)

#---------------------------#
#### Figure: BII results ####
#---------------------------#

bio_measures <- c("bii_", "ab_", "cd_")
ylabs <- c("Intactness", "Abundance", "Similarity")
p_list <- list()
ecoregs_5min <- raster(file.path(processed_path, "biorealms_5min.tif"))
ecoregs_30min <- raster(file.path(processed_path, "biorealms_30min.tif"))
bioreg_names <- read.csv(file.path(out_path, "ecoreglevels.csv"),  na.strings = "-9999")
bioreg_names$Name <- c("Australasia", "Antarctic", "Afrotropic", "Indo-Malay", "Nearctic", "Neotropical", "Oceania", "Palearctic")
colnames(bioreg_names)[3] <- "ID"

levels(ecoregs_30min) <- bioreg_names[,c(3,4)]
levels(ecoregs_5min) <- bioreg_names[,c(3,4)]
cols <- c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB")

for(i in 1:length(bio_measures)){
  bii <- list.files(out_path, pattern = bio_measures[i], full.names = TRUE)
  bii <- bii[grepl("tif", bii)]
  
  bii_ras <- stack(bii[-which(grepl("repr", bii))])
  
  inds_5min <- which(!is.na(mask_5min[]))
  inds_30min <- which(!is.na(mask_30min[]))
  
  bii_5min <- as.data.frame(stack(bii_ras, ecoregs_5min)[inds_5min])
  bii_30min <- as.data.frame(stack(bii[[which(grepl("repr", bii))]], ecoregs_30min)[inds_30min])
  bii_diff <- cbind(bii_5min[c(1,3,4)] -  bii_5min[,2])
  
  bii_diff <- bii_diff[,1:3] %>% gather()
  bii_diff$bioregions <- bii_5min$biorealms_5min
  colnames(bii_diff) <- c("scenario", "diff",  "biorealms")
  
  bii_30min <- as.data.frame(stack(bii[[5]], ecoregs_30min)[inds_30min])
  bii_30min <- cbind(bii_30min, "pres" = projectRaster(raster(bii[2]), mask_30min, method = "bilinear")[inds_30min])
  bii_30min <- data.frame( "scenario" = "repr", "biorealms" = bii_30min$biorealms_30min, "diff" = bii_30min[paste0(bio_measures[i], "repr")] - bii_30min$pres)
  colnames(bii_30min) <- c("scenario", "biorealms", "diff")
  
  bii_plot <- rbind(bii_diff, bii_30min[,match(colnames(bii_diff), colnames(bii_30min))])
  
  bii_plot <- na.omit(bii_plot)
  
  bii_plot$biorealms <- as.factor(bii_plot$biorealms)
  bii_plot$scenario <- as.factor(bii_plot$scenario)
  levels(bii_plot$scenario) <- c("M8.5 demand", "SSP2", "SSP5", "M8.5")
  levels(bii_plot$biorealms) <- bioreg_names$Name[-2] #not antarctic
  bii_plot$scenario <- factor(bii_plot$scenario,levels = levels(bii_plot$scenario)[c(2,3,1,4)],ordered = TRUE)
  
  p <- bii_plot %>% 
    ggplot(aes(x = scenario, y = diff, colour = biorealms)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    # geom_boxplot(outlier.shape = NA, fatten = NULL) + 
    geom_linerange(stat = "summary",
                   fun.min = function(z) {quantile(z,0.25)},
                   fun.max = function(z) {quantile(z,0.75)},
                   fun = mean, 
                   linetype = 1,
                   size = 0.8,
                   alpha = 0.5,
                   position = position_dodge(0.75)) +
    geom_point(stat = "summary", fun = mean, position = position_dodge(0.75), mapping = aes(fill = biorealms), size = 3, shape = 19) +
    #stat_summary(fun=mean, geom="point", shape=19, color="black", aes(group = biorealms), position = position_dodge(0.75)) +
    scale_y_continuous(limits = c(-0.1, 0.005)) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols)
  
  p <- p + 
    theme_minimal() +
    theme(panel.border = element_blank(), 
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position= "none" ,
          axis.title.x = element_blank(),
          legend.title=element_blank()) +
    ylab(ylabs[i])
  
  if(i%in%c(1,2)){
    p <- p + theme(axis.text.x = element_blank())
  }
  
  if(i == 2){
    p <- p + theme(legend.position= c(0.5, 0.1)) +
      guides(colour = guide_legend(ncol = 3))
  }
  
  p_list[[i]] <- p
  print(i)
}

fig <- plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], ncol = 1, rel_heights = c(1,1,1), labels = c("(a)", "(b)", "(c)") , label_size = 10, label_fontface = "plain")

fig_paths <- c(file.path(figure_thesis, "fig_biodiv.pdf"),
               file.path(figure_chapter, "fig_biodiv.pdf"))


multi_pdf(x = fig, paths = fig_paths, width = medium_width, height = 8)

#-----------------------------------------------------#
#### Figure: BII change vs commodity demand change ####
#-----------------------------------------------------#

# Load bii prediction
bio_measures <- c("bii", "ab", "cd")
gtap <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))
inds_5min <- which(!is.na(mask_5min[]))

# RCP8.5
ecoregs_5min <- raster(file.path(processed_path, "biorealms_5min.tif"))
qfe85 <- read.csv(file.path(data_path, "demand calculation", "GTAP trajectories", "qfe_RCP85_2100.csv"))[c(1:9), ]
qfe85[,-1] <- round(qfe85[,-1], 2)
gtap_classes <- qfe85$qfe

qfe_85 <- as.data.frame(t(qfe85[,-1]))
n <- qfe85$qfe
colnames(qfe_85) <- n
i <- 1
for(i in 1:3){
  bii <- list.files(out_path, pattern = bio_measures[i], full.names = TRUE)
  bii <- bii[grepl("tif", bii)]
  bii_ras <- stack(bii[-which(grepl(paste0(c("repr", "dluh"), collapse = "|"), bii))])
  bii_5min <- as.data.frame(stack(bii_ras[[c(1, 3)]], ecoregs_5min)[inds_5min])
  bii_diff <- cbind(bii_5min[,2] -  bii_5min[,1])
  gtap_bii <- data.frame("bii" = bii_diff, "gtap" = gtap[inds_5min])
  colnames(gtap_bii) <- c(bio_measures[i], "gtap")
  aggr <- aggregate(as.formula(paste0(colnames(gtap_bii)[1],"~","gtap")), data = gtap_bii, FUN = mean)
  qfe_85[,9 + i] <- aggr[,2]
  colnames(qfe_85)[9 + i] <- colnames(gtap_bii)[1]
}

commodities <- colnames(qfe_85)[1:9]
qfe_85$countries <- row.names(qfe_85)
qfe_85$bii <- qfe_85$bii * 100

# RCP4.5
qfe45 <- read.csv(file.path(data_path, "demand calculation", "GTAP trajectories", "qfe_RCP45_2100.csv"))[c(1:9), ]
qfe45[,-1] <- round(qfe45[,-1], 2)
gtap_classes <- qfe45$qfe

qfe_45 <- as.data.frame(t(qfe45[,-1]))
n <- qfe45$qfe
colnames(qfe_45) <- n

for(i in 1:3){
  bii <- list.files(out_path, pattern = bio_measures[i], full.names = TRUE)
  bii <- bii[grepl("tif", bii)]
  bii_ras <- stack(bii[-which(grepl(paste0(c("repr", "dluh"), collapse = "|"), bii))])
  bii_5min <- as.data.frame(stack(bii_ras[[c(1, 2)]], ecoregs_5min)[inds_5min])
  bii_diff <- cbind(bii_5min[,2] -  bii_5min[,1])
  gtap_bii <- data.frame("bii" = bii_diff, "gtap" = gtap[inds_5min])
  colnames(gtap_bii) <- c(bio_measures[i], "gtap")
  aggr <- aggregate(as.formula(paste0(colnames(gtap_bii)[1],"~","gtap")), data = gtap_bii, FUN = mean)
  qfe_45[,9 + i] <- aggr[,2]
  colnames(qfe_45)[9 + i] <- colnames(gtap_bii)[1]
}

commodities <- colnames(qfe_45)[1:9]
qfe_45$countries <- row.names(qfe_45)
qfe_45$bii <- qfe_45$bii * 100

qfe1 <- qfe_85 %>%  
  gather(commodity, change , all_of(commodities)) %>% 
  gather(type, measured, "bii", "ab", "cd")

qfe2 <- qfe_45 %>%  
  gather(commodity, change , all_of(commodities)) %>% 
  gather(type, measured, "bii", "ab", "cd")

qfe <- rbind(qfe1, qfe2)

qfe_df <- rbind(qfe_45, qfe_85)

coef_out <- list()

for(j in 1:3){
  coeffs <- data.frame("effect" = numeric(), "ci_low" = numeric(), "ci_high"= numeric(), "type" = character())
  for(i in 1:9){
    tmodel <- lm(paste0(bio_measures[j], "~", commodities[i]), data = qfe_df)
    coefs <- effectsize(tmodel)
    coeffs[i,1:3] <- c(coefs$Std_Coefficient[2], coefs$CI_low[2], coefs$CI_high[2])
    coeffs[i, 4] <- gsub('_', '', bio_measures[j])
  }
  coeffs$commodity <- commodities
  coeffs <- coeffs[order(coeffs$effect, decreasing = TRUE),]
  coeffs$commodity <- factor(coeffs$commodity, levels = coeffs$commodity)
  coef_out[[j]] <- coeffs
}

coefs <- do.call("rbind", coef_out)

coefs$type <- factor(coefs$type, levels = rev(c("bii", "ab", "cd")))
levels(coefs$type) <- c("similarity", "abundance", "BII")
cols <- c("black", "grey", "grey")
group.colors <- c(BII = "black", abundance = "grey", similarity = "grey")
group.lty <- c(BII = "solid", abundance = "dashed", similarity = "dotted")

qfe$commodity <- factor(qfe$commodity, levels = rev(coef_out[[1]]$commodity))

p1 <- qfe %>% 
  filter(type == "bii") %>% 
  ggplot(aes(x = change, y = measured)) +
  facet_wrap(~commodity, scales = "free_x") +
  geom_point(size = 0.7) +
  geom_smooth(method="lm", color = "black", size=0.8) +
  scale_x_continuous(breaks= pretty_breaks()) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("BII difference [%]") +
  xlab("Change in land endowment [%]")

# Russia
qfe %>% 
  filter(type == "bii") %>% 
  filter(measured < -5)

p2 <- coefs %>%
  ggplot(aes(y = commodity, x = effect, group = type, col = type)) + 
  geom_errorbar(aes(xmax=ci_low, xmin=ci_high, linetype = type), width=0, position = position_dodgev(height=0.6)) +
  geom_point(size=2, position = position_dodgev(height=0.6)) +
  geom_vline(xintercept = 0,  linetype="dashed") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank(),
        legend.position= c(0.8, 0.8) ,
        panel.grid.minor = element_blank()) +
  
  xlab(expression(paste("Effect" %+-% "95% CI"))) +
  scale_colour_manual(values=group.colors, guide = guide_legend(reverse = TRUE)) +
  scale_linetype_manual(values=group.lty, guide = guide_legend(reverse = TRUE))


fig <- plot_grid(p2, p1, labels = c("(a)", "(b)"), rel_widths = c(0.8, 1), label_size = 10, label_fontface = "plain")

fig_paths <- c(file.path(figure_thesis, "fig_biigtap.pdf"),
               file.path(figure_chapter, "fig_biigtap.pdf"))

multi_pdf(x = fig, paths = fig_paths, width = full_width, height = full_width/2+0.3)


#---------------------------------------------#
#### Figure: effect sizes of biodiv models ####
#---------------------------------------------#

ab_model <- readRDS(file.path(out_path, "ab_model.rds"))
cd_model <- readRDS(file.path(out_path, "cd_model.rds"))
ab_coefs <- effectsize(ab_model, two_sd = TRUE)
cd_coefs <- effectsize(cd_model, two_sd = TRUE)

ab_coefs <- data.frame("parameter" = ab_coefs$Parameter, "effect" = ab_coefs$Std_Coefficient, "ci_low" = ab_coefs$CI_low, "ci_high" = ab_coefs$CI_high)

ab_coefs$parameter <- c("Primary minimal", "Cropland intense", "Cropland light", "Cropland minimal", "Pasture intense", "Pasture light", "Primary intense", "Primary light", "Secondary intense", "Secondary light", "Secondary minimal", "Urban intense", "Urban minimal")
ab_coefs <- ab_coefs[order(ab_coefs$effect, decreasing = TRUE),]
ab_coefs$parameter <- factor(ab_coefs$parameter, levels = ab_coefs$parameter)

cd_coefs <- data.frame("parameter" = cd_coefs$Parameter, "effect" = cd_coefs$Std_Coefficient, "ci_low" = cd_coefs$CI_low, "ci_high" = cd_coefs$CI_high)

cd_coefs$parameter <- c("Primary minimal", "Primary light", "Urban minimal", "Primary intense", "Urban intense", "Cropland minimal", "Secondary intense", "Pasture light", "Pasture intense", "Secondary light", "Secondary minimal", "Cropland intense", "Cropland light", "Log distance")
cd_coefs <- cd_coefs[order(cd_coefs$effect, decreasing = TRUE),]
cd_coefs$parameter <- factor(cd_coefs$parameter, levels = cd_coefs$parameter)

p1 <- ab_coefs %>% 
  ggplot(aes(y = parameter, x = effect)) + 
  geom_errorbar(aes(xmax=ci_low, xmin=ci_high), width = 0.3) +
  geom_point(size=1.5) +
  geom_vline(xintercept = 0,  linetype="dashed") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position= c(0.9, 0.6) ,
        panel.grid.minor = element_blank()) +
  xlab(expression(paste("Effect on abundance" %+-% "95% CI"))) +
  scale_colour_manual(values=group.colors, guide = guide_legend(reverse = TRUE)) +
  scale_linetype_manual(values=group.lty, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(-1, 1.5))


p2 <- cd_coefs %>% 
  ggplot(aes(y = parameter, x = effect)) + 
  geom_errorbar(aes(xmax=ci_low, xmin=ci_high), width = 0.3) +
  geom_point(size=1.5) +
  geom_vline(xintercept = 0,  linetype="dashed") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position= c(0.9, 0.6) ,
        panel.grid.minor = element_blank()) +
  xlab(expression(paste("Effect on similarity" %+-% "95% CI"))) +
  scale_colour_manual(values=group.colors, guide = guide_legend(reverse = TRUE)) +
  scale_linetype_manual(values=group.lty, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(-1, 1.5))

fig <- plot_grid(p1, p2, labels = c("(a)", "(b)"), ncol = 1, rel_heights = c(1 - 1/13, 1), label_size = 10, label_fontface = "plain")


fig_paths <- c(file.path(figure_thesis, "fig_bdeffects.pdf"),
               file.path(figure_chapter, "fig_bdeffects.pdf"))

multi_pdf(x = fig, paths = fig_paths, width = full_width, height = full_width - 0.5)

#---------------------------------------------------#
#### Figure: Biodiversity metrics for each class ####
#---------------------------------------------------#

bii <- read.csv(file.path(out_path, "cd_ab_BII.csv"))
bii <- bii[order(bii$BII, decreasing = TRUE),]
p1 <- 
  bii %>% 
  arrange(desc(BII)) %>% 
  mutate("class" = factor(X, levels = X),
         "similarity" = cd_m_preds,
         "abundance" = ab_m_preds,
         .keep = "unused") %>% 
  gather(key = key, value = value, similarity, abundance, BII) %>% 
  
  #make plot
  ggplot(aes(x = class, y = value, col = key, shape = key)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position= c(0.25, 0.2)) +
  coord_flip() +
  scale_colour_manual(values=c("grey", "black", "grey")) +
  scale_shape_manual(values= c(15, 19, 17))

fig_paths <- c(file.path(figure_thesis, "fig_biodivbyclass.pdf"),
               file.path(figure_chapter, "fig_biodivbyclass"))
multi_pdf(x = p1, paths = fig_paths, width = single_width, height = single_width)

#-----------------------------#
#### Figure: land use maps ####
#-----------------------------#

l <- list()
rasters <- list()
landuses <- c("Cropland", "Pasture", "Primary", "Secondary")
for(i in 1:length(landuses)){
  pm_dluh <- raster(file.path(processed_path, paste0("dluh_",tolower(landuses[i]), "_lu_5min.tif")))
  pm_rcp85 <- raster(file.path(processed_path, paste0("rcp85_", tolower(landuses[i]), "_lu_5min.tif")))
  rasters[[i]] <- stack(pm_dluh, pm_rcp85)
}
rasters <- stack(rasters)
panel_names <- c("M8.5 Cropland", "SSP5 Cropland", "M8.5 Pasture", "SSP5 Pasture", "M8.5 Primary", "SSP5 Primary", "M8.5 Secondary", "SSP5 Secondary")
l <- levelplot(rasters, margin = F, colorkey = list(space='bottom',
                                                    labels=list(cex=0.8, font=1, width = 0.7),
                                                    height=0.5,width=1),
               xlab = list(label = "Fractional cover", vjust = -.2, cex = 1),
               ylab="",
               par.settings=list(
                 strip.border=list(col='transparent'),
                 strip.background=list(col='transparent'),
                 axis.line=list(col='transparent')),
               scales=list(draw=FALSE),
               col.regions = viridis(100),
               maxpixels = 5e4, layout=c(2, 4), names.attr = panel_names)

l
fig_paths <- c(file.path(figure_thesis, "fig_landusemaps.pdf"),
               file.path(figure_chapter, "fig_landusemaps.pdf"))

multi_pdf(x = l, paths = fig_paths, width = full_width, height = full_width * 1.2)

#--------------------------------------------------------------#
#### Figure: land use type/intensity changes (SSP5 vs M8.5) ####
#--------------------------------------------------------------#

raw_path <-  file.path("/Volumes", "external", "OneDrive - The University of Melbourne", "PhD - Large Files", "raw data")
processed_path <- file.path(raw_path, "Global", "processed rasters")

mask_5min <- raster(list.files(processed_path, pattern = "mask_5min", full.names = TRUE))
inds_5min <- which(!is.na(mask_5min[]))
gtap_aggregation_5min <- raster(file.path(processed_path, "gtap_aggregation_5min.tif"))

luh_demand <- 
  # load data into stack
  processed_path %>%  
  list.files(pattern = "dluh|gtap", full.names = TRUE) %>% 
  stack %>% 
  
  # turn to df and subset
  as.data.frame %>% 
  select(matches("minimal|light|intense|gtap")) %>% 
  slice(inds_5min) %>% 
  # calculate mean by region
  group_by(gtap_aggregation_5min) %>% 
  summarise(across(everything(), list(mean))) %>% 
  mutate(scenario = "M8.5 demand")

class_names <- c("Cropland intense", "Cropland light", "Cropland minimal", "Pasture intense", "Pasture light", "Primary intense", "Primary light", "Primary minimal", "Secondary intense", "Secondary light", "Secondary minimal", "Urban intense", "Urban minimal")

colnames(luh_demand) <- c("GTAP", class_names, "scenario")

ssp5_list <- 
  processed_path %>%  
  list.files(pattern = "rcp85|gtap", full.names = TRUE)

ssp5_list <- ssp5_list[which(grepl(pattern = "minimal|light|intense|gtap", ssp5_list))]

ssp5_demand <- 
  # load data into stack
  ssp5_list %>% 
  stack %>% 
  
  # turn to df and subset
  as.data.frame %>% 
  slice(inds_5min) %>% 
  
  # calculate mean by region
  group_by(gtap_aggregation_5min) %>% 
  summarise(across(everything(), list(mean))) %>% 
  mutate(scenario = "SSP5")


pres_list <- 
  processed_path %>%  
  list.files(pattern = "pres|gtap", full.names = TRUE)

pres_list <- pres_list[which(grepl(pattern = "minimal|light|intense|gtap", pres_list))]

pres_demand <- 
  # load data into stack
  pres_list %>% 
  stack %>% 
  
  # turn to df and subset
  as.data.frame %>% 
  slice(inds_5min) %>% 
  
  # calculate mean by region
  group_by(gtap_aggregation_5min) %>% 
  summarise(across(everything(), list(mean)))  %>% 
  mutate(scenario = "present")

colnames(ssp5_demand) <- colnames(pres_demand) <- colnames(luh_demand)

demand <- rbind(ssp5_demand, luh_demand, pres_demand) %>% 
  gather("land_use", "demand", -one_of("GTAP", "scenario")) %>% 
  separate(land_use, into = c("type", "intensity"), sep = " ", remove = FALSE) %>% 
  mutate(land_use = factor(land_use, levels = rev(class_names)),
         type = factor(type, levels = rev(unique(type))))
demand$scenario

p1 <- ggplot(demand, aes(x = land_use, y = demand,  colour = scenario)) +
  geom_point(stat = "summary", fun = mean, position = position_dodge(0.75), size = 3, shape = 19) +
  geom_point(position = position_dodge(0.75), size = 1, shape = 20) +
  coord_flip(ylim = c(0, 0.4)) + 
  scale_colour_manual(values=viridis::cividis(3, begin = 0.1, end = 0.9)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position= "none") +
  ylab("Fraction of landscape")

p1
p2 <- ggplot(demand, aes(x = type, y = demand, colour = scenario)) +
  geom_point(stat = "summary", fun = mean, position = position_dodge(0.75), size = 3, shape = 19) +
  geom_point(position = position_dodge(0.75), size = 1, shape = 20) +
  coord_flip(ylim = c(0, 0.4)) + 
  scale_colour_manual(values=viridis::cividis(3, begin = 0.1, end = 0.9)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.7, 0.1)) +
  ylab("Fraction of landscape")

fig <- plot_grid(p2, p1, ncol =2, rel_widths = c(0.85, 1), labels = c("(a)", "(b)"), label_size = 10, label_fontface = "plain")

fig_paths <- c(file.path(figure_thesis, "fig_landusecomp.pdf"),
               file.path(figure_chapter, "fig_landusecomp.pdf"))

multi_pdf(x = fig, paths = fig_paths, width = full_width, height = full_width)


#----------------------------------------#
#### Table: land-use model parameters ####
#----------------------------------------#

# land use data
# get regional data from global data set:

# indices of current region
lu_global <- stack(list.files(processed_path, pattern = "lu_5min", full.names = TRUE, recursive = TRUE))
lu_global <- lu_global[[which(grepl(pattern = paste0(c("^urb", "^crop", "^pas", "^pri", "^sec"), collapse = "|"), names(lu_global)))]]
covs_files <- list.files(processed_path, pattern = "5min", full.names = TRUE, recursive = TRUE)
covs_global <- stack(covs_files[which(!grepl(pattern = paste0(c("urb", "crop", "pas", "pri", "sec"), collapse = "|"), covs_files))])

weights <- list(matrix(1/9, 3, 3, byrow= TRUE)) #size of window
weights <- rep(weights, length.out = 5)

# get non-na cells
inds <- which(!is.na(mask_5min[]))
neigh_global <- neighbourhood(lu_global, weights = weights, mask = mask_5min, suffix = "neigh", cols = 1:5, format = "stack")

covs_suitmodel <- covs_global[[-grep(paste(c("biorealms", "rcp45", "rcp85", "lu", "ssp", "gtap", "unsubregions", "gls", "pa", "mask"), collapse = "|"), names(covs_global))]]

# land use data
landuse <- lu_global[inds]
colnames(landuse) <- unlist(lapply(strsplit(colnames(landuse), "_"), FUN = function(x) x[[1]]))
landuse <-landuse[,c(3,4,1,2,5)]
# env covariates
covs <- covs_suitmodel[inds]
neigh <- neigh_global[inds]
data <- cbind(covs, neigh)

# Determine correlations in data and reduce predictor set
preds <- colnames(correlations(data, sub = 100000))
data <- data[,preds]
subs <- sample(1:nrow(landuse), 100000)
counts <- integerify(x = landuse[subs, ], resolution = 10000)
colnames(counts) <- colnames(landuse)
data_sub <- as.data.frame(data[subs, ])
form <- paste(colnames(data), collapse = "+")
f <- as.formula(paste("counts", "~", form))
suit_model <- nnet::multinom(f, data = data_sub, model = TRUE, maxit = 1000)
ef <- effectsize(suit_model)
ef <- as.data.frame(ef)
ef$Response <- str_to_title(ef$Response)
ef$Parameter <- rep(c("intercept", "pop. density", "bio13", "bio15", "bio19", "bio2", "bio7", "bio8", "dist built-up", "dist lakes", "dist rivers", "soil phosphorus", "slope", "elevation", "wilting point", "Primary neigh", "Secondary neigh", "Urban neigh") ,4)
ef$Parameter <- as.factor(ef$Parameter)
ef$Parameter <- relevel(ef$Parameter, "intercept")
class(suit_model)
p1 <- ef %>% 
  ggplot(aes(y = Parameter, x = Std_Coefficient, shape = Response)) + 
  geom_point(size=1.5, position = position_dodge(0.5)) +
  geom_vline(xintercept = 0,  linetype="dashed") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.minor = element_line(colour = "grey",size=0.1),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position= c(0.2, 0.5)) +
  xlab(expression(paste("Standardized effect")))


fig_paths <- c(file.path(figure_thesis, "fig_landuseeffects.pdf"),
               file.path(figure_chapter, "fig_landuseeffects.pdf"))

multi_pdf(x = p1, paths = fig_paths, width = full_width, height = full_width)

# DIscussion figure

# files <- list.files("/Volumes/external/OneDrive - The University of Melbourne/PhD - Large Files/raw data/Global/LUHa_u2t1.v1/updated_states/", full.names = TRUE)
# pri_list <- files[which(grepl("gothr", files))]
# sec_list <- files[which(grepl("gsec", files))]
# 
# ma <- matrix(ncol = 2, nrow = length(pri_list))
# for(i in 1:length(pri_list)){
#   ma[i,1] <- mean(raster(pri_list[[i]])[])
#   ma[i,2] <- mean(raster(sec_list[[i]])[])
#   print(i)
# }
# 
# ma <- as.data.frame(ma)
# ma$year <- 1700:2005
# colnames(ma) <- c("primary", "secondary", "year")
# ma %>% 
#   gather(key = "type", value = "coverage", primary, secondary) %>% 
#   ggplot(aes(x = year, y = coverage, col = type)) +
#            geom_line()
# tail(ma, 20)
#          