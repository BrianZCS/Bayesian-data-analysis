###########SETUP##############
library(tidyverse)
rm(list=ls())

# Your working directory goes here
setwd("C:/Users/caleb/Desktop/479/Project/full_dat")

load("tidy_all_forest_data.Rds")

#############FILTER DATA#################

# 102 million observations to start
lvl <- levels(tidy.df$name); lvl

# filter to common private forest types
filt.df <- tidy.df %>% 
  filter(name %in% lvl[c(1,7,9,10,11,12,14,15,17,18)]) %>% 
  mutate(year = as.factor(year)) %>% 
  droplevels()

levels(filt.df$name)

rm(tidy.df)
str(filt.df)

count.df <- filt.df %>% 
  group_by(zone, name) %>% 
  summarise(count = n()/21)

count.df %>% 
  ggplot(aes(x=name, y=log(count), fill=name)) +
    geom_col() +
    facet_wrap(~zone) +
    labs(title="Forest Type Log(count) by Zone") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

# remove zone-type categories with log(count) < 5 -> 148 observations - 22 categories in total
to.remove <- which(log(count.df$count) < 5)
l <- length(to.remove)
to.rm.df <- count.df[to.remove,]

for(i in 1:l) {
  cat(paste0("\nremoving zone ", to.rm.df$zone[i], "\t type: ", to.rm.df$name[i]))
  filt.df <- filt.df %>% 
    filter(!(zone==to.rm.df$zone[i] & name==to.rm.df$name[i]))
}

save(filt.df, file="filt_df.Rds")
load("count_df.Rds")
################CALC GRAND MEAN#####################
rm(list=ls())
load("filt_df.Rds")
mean(filt.df$npp, na.rm=T)

##############CHECK COUNTS#################
# recalculate counts after filtering  - down to 78.5 million observations
rm(count.df)

count.df <- filt.df %>% 
  group_by(zone, name) %>% 
  summarise(count = n()/21)

save(count.df, file="count_df.Rds")

load("count_df.Rds")

count.df %>% 
  ggplot(aes(x=name, y=log(count), fill=name)) +
  geom_col() +
  facet_wrap(~zone) +
  labs(title="Forest Type Log(n) by Zone (Filtered)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

head(filt.df)

############IQR vs COUNT - RAW DATA#################


###########IQR vs COUNT - AGGREGATED DATA################



###############AGGREGATE MEDIAN AND SD##################
# calculate median for each zone-type-year combination
load("med_df.Rds")

med.df <- filt.df %>% 
  group_by(zone, name, year) %>% 
  summarise(med = median(npp, na.rm=T))

# plot medians - points over time
# something weird going on with aspen forested wetland, deciduous/coniferous wetland 
# in years 2020 and 2021
# in zones 5, 7, 8, 9, 10, 12, 14, 15
# see below - zero.df

med.df %>% 
  ggplot(aes(x=year, y=med, color=name)) +
  geom_point() +
  facet_wrap(~zone) +
  labs(title="Median NPP vs. Year: by Type, Zone", y="Median NPP")

med.df %>% 
  ggplot(aes(x=year, y=med, color=zone)) +
  geom_point() +
  facet_wrap(~name) +
  labs(title="Median NPP vs. Year: by Zone, Type", y="Median NPP")

# medians - lines over time - by zone
med.df %>% 
  mutate(year = as.numeric(year)) %>% 
  ggplot(aes(x=year, y=med, color=name)) +
  geom_line() +
  facet_wrap(~zone) +
  labs(title="Median NPP vs. Year: by Type, Zone", y="Median NPP")

med.df %>% 
  mutate(year = as.numeric(year)) %>% 
  ggplot(aes(x=year, y=med, color=zone)) +
  geom_line() +
  facet_wrap(~name) +
  labs(title="Median NPP vs. Year: by Zone, Type", y="Median NPP")

# boxplot - wrapped by zone
med.df %>% 
  ggplot(aes(x=name,y=med, color=name)) +
  geom_boxplot() +
  facet_wrap(~zone) +
  labs(title="Median NPP (all years): by Zone, Type", y="Median NPP", x="") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# boxplot - wrapped by type
med.df %>% 
  ggplot(aes(x=zone,y=med, color=zone)) +
  geom_boxplot() +
  facet_wrap(~name) +
  labs(title="Median NPP (all years): by Type, Zone", y="Median NPP", x="zone")

save(med.df, file="med_df.Rds")

# some years have a median of 0 - we may want to remove those years, but for now they're in the data
zero.ind <- which(med.df$med == 0); zero.ind
zero.df <- med.df[zero.ind,]; zero.df
unique(zero.df$name)
unique(zero.df$zone)
unique(zero.df$year)

med.var.df <- med.df %>% 
  group_by(name, zone) %>% 
  summarise(med_var = var(med, na.rm=T)) %>% 
  mutate(med_sd = sqrt(med_var)) %>% 
  mutate(log_sd = log(med_sd))

save(med.var.df, file="med_var_df.Rds")

rm(count.df); rm(filt.df); rm(to.rm.df); rm(zero.df)


###############FIT MLR###################
setwd("C:/Users/caleb/Desktop/479/Project/full_dat/share_final_data")
load("med_df.Rds")
load("med_var_df.Rds")

mean(med.df$med)
min(med.var.df$med_sd)
levels(med.df$name)

mean(med.var.df$log_sd)

med.lm <- lm(med ~ name + zone + year, data=med.df); summary(med.lm)

mean(med.df$med)

# median coefficients (scalars fed to Bayesian model as params)
names(med.lm$coefficients)
mu.intercept <- med.lm$coefficients[1]; mu.intercept
mu.type.coefs <- med.lm$coefficients[2:10]; mu.type.coefs
mu.zone.coefs <- med.lm$coefficients[11:25]; mu.zone.coefs
mu.year.coefs <- med.lm$coefficients[26:45]; mu.year.coefs

mean(mu.type.coefs); sd(mu.type.coefs)
mean(mu.zone.coefs); sd(mu.zone.coefs)
mean(mu.year.coefs); sd(mu.year.coefs)

# histograms
hist(mu.type.coefs)
hist(mu.zone.coefs)
hist(mu.year.coefs)

# save
setwd("C:/Users/caleb/Desktop/479/Project/full_dat/share_final_data/stan_dat")
getwd()
save(mu.intercept, file="mu_intercept.Rds")
save(mu.type.coefs, file="mu_type_coefs.Rds")
save(mu.zone.coefs, file="mu_zone_coefs.Rds")
save(mu.year.coefs, file="mu_year_coefs.Rds")

# sd coefficients (scalars fed to Bayesian model as params)
sd.lm <- lm(log_sd ~ 0 + zone + name, data=med.var.df); summary(sd.lm)

names(sd.lm$coefficients)
#sd.intercept <- sd.lm$coefficients[1]; sd.intercept
sd.zone.coefs <- sd.lm$coefficients[1:16]; sd.zone.coefs
sd.type.coefs <- sd.lm$coefficients[17:25]; sd.type.coefs

# mean and sd of zone and name sd coefs
mean(sd.zone.coefs); sd(sd.zone.coefs)

mean(sd.type.coefs); sd(sd.type.coefs)

# hists
hist(sd.zone.coefs, breaks=10)
hist(sd.type.coefs, breaks=10)

# save
save(sd.intercept, file="sd_intercept.Rds")
save(sd.zone.coefs, file="sd_zone_coefs.Rds")
save(sd.type.coefs, file="sd_type_coefs.Rds")
