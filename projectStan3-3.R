###############Initialize##################

rm(list=ls())

library(rstan)

# YOUR WORKING DIRECTORY HERE
working.dir <- "C:/Users/caleb/Desktop/479/Project/full_dat/share_final_data/stan_dat/" 

setwd(working.dir)

#################Load Data#####################

# median data matrix - aggregated data points being predicted
load("median/dat_mat.Rds")

#  mu data
#  intercept:
mu.intercept = 7212 # grand mean of median values from med.df - in med_df.Rds

#  type scale - hyperparameters elicited using OLS multiple regression
type.coefs.scale <- 420 
#  zone scale:
zone.coefs.scale <- 371
# year scale:
year.coefs.scale <- 807

# sd data - just intercept
#  intercept:
sd.intercept <- 6.854 # grand mean of log(sd) values from med.var.df - in med_var_df.Rds

J = 10 # number species
K = 16 # number zones
L = 21 # number years

################Construct Model Object####################
# construct model object
npp_model <- stan_model(file = paste0(working.dir,"npp_model3-3.stan"))

################Fit Stan Model#####################
npp_data <- list(J=J, K=K, L=L,
                 dat_mat = dat.mat,
                 mu_intercept = mu.intercept,
                 type_coefs_scale = type.coefs.scale,
                 zone_coefs_scale = zone.coefs.scale,
                 year_coefs_scale = year.coefs.scale,
                 sd_intercept = sd.intercept
                 )

npp_fit <- sampling(object = npp_model,
                    data = npp_data,
                    iter=5000, 
                    chains=4)


###########Summary################
npp_summary <- summary(npp_fit)

str(npp_summary)

npp_summary$summary

summary.df <- as.data.frame(npp_summary$summary)

min(summary.df$Rhat[!is.na(summary.df$Rhat)])
max(summary.df$Rhat[!is.na(summary.df$Rhat)])

min(summary.df$n_eff[!is.na(summary.df$n_eff)])
max(summary.df$n_eff[!is.na(summary.df$n_eff)])
low.n.ind <- which((summary.df$n_eff-1440.089)<200)
summary.df[low.n.ind,]
#npp_summary$c_summary

############Extract Samples#################
library(tidyverse)

# load type.df - associates type numbers with names
load("type_df.Rds")
str(type.df)

# choose your own adventure - samples of any parameter can be meaningful
# e.g. in identifying whether one zone > another zone, one species > another species

# see possible variables to extract
# std_*_coefs are standardized - just in terms of sd
str(rstan::extract(npp_fit))

#################Posterior Predictive - Standardized Coefficients###########################
# extract samples for a particular param:

# standardized type coefficients
std_type_coefs <- rstan::extract(npp_fit, par="std_type_coefs")[["std_type_coefs"]]
str(std_type_coefs)
# convert to df
std.type.df <- as.data.frame(std_type_coefs)
# rename columns
colnames(std.type.df) <- 1:10
# make tidy
library(tidyverse)
std.type.df <- std.type.df %>% 
  pivot_longer(names_to="type", cols=1:10) %>%
  mutate(type = as.numeric(type)) %>% 
  arrange(type) %>% 
  mutate(type = as.factor(type)) %>% 
  right_join(type.df, by="type")
# check
str(std.type.df)
levels(std.type.df$name)
# ggplot
std.type.df %>% 
  ggplot(aes(x=value, color=name)) +
  geom_density() +
  facet_wrap(~name) +
  labs(title="Standardized Type Coefficient Estimates")
  
# standardized zone coefficients
std_zone_coefs <- rstan::extract(npp_fit, par="std_zone_coefs")[["std_zone_coefs"]]
str(std_zone_coefs)
# to df
std.zone.df <- as.data.frame(std_zone_coefs)
str(std.zone.df)
# rename cols
colnames(std.zone.df) <- 1:16
# make tidy
std.zone.df <- std.zone.df %>% 
  pivot_longer(names_to="zone", cols=1:16) %>% 
  mutate(zone = as.numeric(zone)) %>% 
  arrange(zone) %>% 
  mutate(zone = as.factor(zone))
# ggplot
std.zone.df %>% 
  ggplot(aes(x=value, color=zone)) +
  geom_density() +
  facet_wrap(~zone) +
  labs(title="Standardized Zone Coefficient Estimates")

# standardized year coefficients
std_year_coefs <- rstan::extract(npp_fit, par="std_year_coefs")[["std_year_coefs"]]
str(std_year_coefs)
std.year.df <- as.data.frame(std_year_coefs)
# rename cols
colnames(std.year.df) <- 1:21
# make tidy
std.year.df <- std.year.df %>% 
  pivot_longer(names_to="year", cols=1:21) %>% 
  mutate(year = as.numeric(year)) %>% 
  arrange(year) %>% 
  mutate(year = as.factor(year))
# ggplot
std.year.df %>% 
  ggplot(aes(x=value, color=year)) + 
  geom_density() +
  facet_wrap(~year) +
  labs(title="Standardized Year Coefficient Estimates")


##################Posterior Predictive######################
y_post_samples <- rstan::extract(npp_fit, par="post_pred")[["post_pred"]]
str(y_post_samples)
save(y_post_samples, file="observed_y_post_samples.Rds")

# posterior predictive
# Aspen Zone 1, year 2000
hist(y_post_samples[,1,1,1], main="Posterior Predictive NPP: Aspen in Zone 1, Year 2000")

# Aspen Zone 1
hist(y_post_samples[,1,1,], main="Posterior Predictive NPP: Aspen in Zone 1")
# Aspen Zone 6
hist(y_post_samples[,1,6,], main="Posterior Predictive NPP: Aspen in Zone 6, All Years")
# Aspen Zone 7
hist(y_post_samples[,1,7,], main="Posterior Predictive NPP: Aspen in Zone 7, All Years")

quantile(y_post_samples[,1,6,], probs=c(0.1,0.9))
quantile(y_post_samples[,1,7,], probs=c(0.1,0.9))

# Mixed Deciduous
hist(y_post_samples[,3,,], main="Posterior Predictive NPP: Mixed Deciduous, All Zones + Years")
# Northern Hardwoods
hist(y_post_samples[,4,,], main="Posterior Predictive NPP: Northern Hardwoods, All Zones + Years")

# just in zones 8 and 9
hist(y_post_samples[,3,8:9,], main="Posterior Predictive NPP: Mixed Deciduous, Zones 8,9")
# Northern Hardwoods
hist(y_post_samples[,4,8:9,], main="Posterior Predictive NPP: Northern Hardwoods, Zones 8, 9")

quantile(y_post_samples[,3,8:9,] , probs=c(0.1,0.9))
quantile(y_post_samples[,4,8:9,] , probs=c(0.1,0.9))

# Aspen all zones
hist(y_post_samples[,1,,], main="Posterior Predictive NPP: Aspen in All Zones")

# zone 1 all types
hist(y_post_samples[,,1,], main="Posterior Predictive NPP: All Types, Zone 1")

# year 13 all types
hist(y_post_samples[,,,13], main="Posterior Predictive NPP: All Types, Year 13")


