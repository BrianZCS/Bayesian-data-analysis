###############Initialize##################

rm(list=ls())

library(tidyverse)
library(RColorBrewer)
library(rstan)

set.seed(12345)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("posterior_data.Rdata")
load("med_df.Rds")

med.df <- med.df %>% droplevels()
type.lvl <- levels(med.df$name)
y_post_samples <- rstan::extract(npp_fit, par="post_pred")[["post_pred"]]

# loads in type.lvl - vector of species in order of their identifier/index 
# loads in y_post_samples - 3D list/data-frame structure of the posterior distribution.

# vector of zone order from south to north based on median latitude coordinate.
# zone order is determined by median of latitude coordinates for pixels in 
# a region.
zone_order <- c(15, 13, 14, 11, 12, 8, 5, 6, 2, 3, 4, 1, 0, 7, 10, 9)

# data frame to store zone/region, zone order, forest species, 
# median NPP (50% quantile) from posterior for species-zone,
# lower NPP (2.5% quantile), and upper NPP (97.5% quantile) 
df <- data.frame() 

# vector to store zone order of Mixed Deciduous/Coniferous Forest
zone_order_vec <- c()
# vector to store posterior predictions for Mixed Deciduous/Coniferous Forest
# across all years and iterations.
mdcf_vec <- c()

# loop through species/forest types
for (si in 1:10) {
  species_name <- type.lvl[si]
  # loop through zones
  for (zi in 1:16) {
    zo <- zone_order[zi] # get zone order from south to north
    s_z_post <- c(y_post_samples[,si,zi,]) # get posterior predictions for species-zone across all years and iterations
    if (species_name == "Mixed Deciduous/Coniferous Forested Wetland") {
        zo_vec <- rep(zo, length(s_z_post))
        mdcf_vec <- c(mdcf_vec, s_z_post)
        zone_order_vec <- c(zone_order_vec, zo_vec)
    }
    s_z_l95 <- as.numeric(quantile(s_z_post, probs=0.025, na.rm=T)[1])
    s_z_u95 <- as.numeric(quantile(s_z_post, probs=0.975, na.rm=T)[1])
    s_z_median <- as.numeric(quantile(s_z_post, probs=0.5, na.rm=T)[1])
    
    tmp.df <- data.frame(Median_NPP=s_z_median, Low_NPP=s_z_l95, Upper_NPP=s_z_u95, Zone_ID=zi, Zone_Order=zo, Species_Name=species_name) 
    df <- rbind(df, tmp.df)
    }
}

# color for region zones (gradient palette because coloring assigned based on zone order)
region_colors_num <- 16
region_colors <- colorRampPalette(brewer.pal(8, "PuBu"))(region_colors_num)

# create first plot of zone median bar plots for species (grids) with bars representing zones ordered south to north
png("Zone_Facet_Grid.png", units='in', height=5, width=15, res=600)
dodge <- position_dodge(width=0.9)
ggplot(df, aes(x=reorder(Zone_ID, Zone_Order), y=Median_NPP, fill=as.factor(Zone_Order))) + geom_bar(color='black', stat='identity') + 
  geom_errorbar(aes(ymin = Low_NPP, ymax = Upper_NPP), position = dodge, width = 0.25) + scale_fill_manual(values=region_colors) + facet_wrap(~Species_Name) + theme_classic()
dev.off()

# create second plot for posterior-prediction of Mixed Deciduous/Coniferous Forest" across regions (gridds)  
png("Zone_Giant_Histogram.png", units='in', height=3, width=10, res=600)
mdcf.df <- data.frame(Zone_Order=zone_order_vec, MDCF_Posterior=mdcf_vec) 
ggplot(mdcf.df, aes(x=MDCF_Posterior, fill=as.factor(Zone_Order))) + geom_histogram(color='black', bins=50) + 
  theme_classic() + facet_wrap(~Zone_Order) + scale_fill_manual(values=region_colors) + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())
dev.off()