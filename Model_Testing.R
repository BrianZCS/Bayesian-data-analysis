###############Initialize##################

rm(list=ls())

library(rstan)
library(tidyverse)
library(hash)
# YOUR WORKING DIRECTORY HERE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#################Load Data#####################
J = 10 # number species
K = 16 # number zones
L = 21 # number years

load("med_df.Rds")
med.df <- med.df %>% droplevels()

# randomly partition dataset for model parameter imputation + modeling
# and testing. 80% used for modeling and 20% used for testing.
random_selection <- sample(seq(1, dim(med.df)[1], by=1), 0.8*dim(med.df)[1])
model.df <- med.df[random_selection,]
test.df <- med.df[-random_selection,]

# Create data matrix of the subsampled data
dat.mat <- array(-1, dim=c(J, K, L))
type.lvl <- levels(model.df$name)
for(l in 1:L) { # year loop
  mat.tmp <- array(-1, dim=c(J,K))
  year.tmp <- model.df %>% filter(year == (l-1))
  for(j in 1:J) { # species/type loop
    type.tmp <- year.tmp %>% filter(name == type.lvl[j])
    ind.tmp <- as.numeric(type.tmp$zone)
    mat.tmp[j,ind.tmp] <- as.numeric(type.tmp$med)
  }
  # for each year - assign dat.mat layer to year.tmp
  dat.mat[,,l] <- mat.tmp
}

# Mu Related Inputs to Model
mu.intercept <- mean(model.df$med)

med.lm <- lm(med ~ name + zone + year, data=model.df)

mu.type.coefs <- med.lm$coefficients[2:10]
mu.zone.coefs <- med.lm$coefficients[11:25]
mu.year.coefs <- med.lm$coefficients[26:45]

zone.coefs.scale <- sd(mu.zone.coefs)
type.coefs.scale <- sd(mu.type.coefs)
year.coefs.scale <- sd(mu.year.coefs)

# SD Related Inputs to Model
med.var.df <- model.df %>% 
  group_by(name, zone) %>% 
  summarise(med_var = var(med, na.rm=T)) %>% 
  mutate(med_sd = sqrt(med_var)) %>% 
  mutate(log_sd = log(med_sd))

sd.intercept <- mean(med.var.df$log_sd)

################Construct Posterior Model Object####################
# construct model object
npp_model <- stan_model(file = "npp_model3-3.stan")

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


###########Summary and Rhat Testing################
npp_summary <- summary(npp_fit)
print(quantile(npp_summary[[1]][,"Rhat"], na.rm=T))

############Test Model Predictions vs. Observed Medians#################

y_post_samples <- rstan::extract(object=npp_fit, pars="post_pred")[["post_pred"]]

# create hash/dictionaries mapping years to year indicies, species to species indices, etc.
years_vector <- seq(0, 20, by=1)
species_vector <- c(type.lvl)

# [n, type, zone, year]
posterior_vs_observed <- data.frame()
for (test_index in 1:dim(test.df)[1]){
  test_year <- test.df$year[test_index]
  test_zone <- test.df$zone[test_index]
  test_name <- test.df$name[test_index]
  observed <- test.df$med[test_index]
  test_year_index <- which(years_vector == test_year)[1]
  test_zone_index <- test_zone
  test_species_index <- which(species_vector == test_name)[1]
  predictions <- y_post_samples[,test_species_index, test_zone_index, test_year_index]
  pred_l95 <- as.numeric(quantile(predictions, probs = 0.025, na.rm=T))
  pred_u95 <- as.numeric(quantile(predictions, probs = 0.975, na.rm=T))
  pred_range <- pred_u95 - pred_l95
  pred_mean <- mean(predictions)
  obs_in_interval <- (observed >= pred_l95) & (observed <= pred_u95)
  test_case <- data.frame(test_index, test_year, test_zone, test_name, observed, pred_l95, pred_u95, pred_mean, obs_in_interval, pred_range)
  posterior_vs_observed <- rbind(posterior_vs_observed, test_case)
}

png("Observation_vs_Posterior_Predictions.png", units='in', height=5, width=12, res=600)
ggplot(posterior_vs_observed, aes(x=reorder(as.factor(test_index), observed), xend=reorder(as.factor(test_index), observed))) + geom_segment(aes(y=pred_l95, yend=pred_u95), color='grey') +
  geom_point(aes(y=pred_mean), color='black') + geom_point(aes(y=observed, color=obs_in_interval), shape=2) + theme_classic() + ylab("NPP") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

################Construct Prior Model Object####################
# construct model object
npp_model <- stan_model(file = "prior_pred_cal.stan")


################Fit Stan Model#####################
npp_data <- list(
  J=J, K=K, L=L,
  mu_intercept = mu.intercept,
  type_coefs_scale = type.coefs.scale,
  zone_coefs_scale = zone.coefs.scale,
  year_coefs_scale = year.coefs.scale,
  sd_intercept = sd.intercept
)

npp_prior_fit <- sampling(object = npp_model, algorithm = "Fixed_param",
                    data = npp_data,
                    iter=5000, 
                    chains=4)

y_prior_samples <- rstan::extract(npp_prior_fit, par="y_prior")[["y_prior"]]

# [n, type, zone, year]
prior_vs_observed <- data.frame()
for (test_index in 1:dim(test.df)[1]){
  test_year <- test.df$year[test_index]
  test_zone <- test.df$zone[test_index]
  test_name <- test.df$name[test_index]
  observed <- test.df$med[test_index]
  test_year_index <- which(years_vector == test_year)[1]
  test_zone_index <- test_zone
  test_species_index <- which(species_vector == test_name)[1]
  predictions <- y_prior_samples[,test_species_index, test_zone_index, test_year_index]
  pred_l95 <- as.numeric(quantile(predictions, probs = 0.025))
  pred_u95 <- as.numeric(quantile(predictions, probs = 0.975))
  pred_mean <- mean(predictions)
  obs_in_interval <- (observed >= pred_l95) & (observed <= pred_u95)
  pred_range <- pred_u95 - pred_l95
  test_case <- data.frame(test_index, test_year, test_zone, test_name, observed, pred_l95, pred_u95, pred_mean, obs_in_interval, pred_range)
  prior_vs_observed <- rbind(prior_vs_observed, test_case)
}

png("Observation_vs_Prior_Predictions.png", units='in', height=5, width=12, res=600)
ggplot(prior_vs_observed, aes(x=reorder(as.factor(test_index), observed), xend=reorder(as.factor(test_index), observed))) + geom_segment(aes(y=pred_l95, yend=pred_u95), color='grey') +
  geom_point(aes(y=pred_mean), color='black') + geom_point(aes(y=observed, color=obs_in_interval), shape=2) + theme_classic() + ylab("NPP") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

interval.range.df <- data.frame(V1=prior_vs_observed$pred_range, V2=rep('0_Prior-Predictive', 530))
interval.range.df <- rbind(interval.range.df, data.frame(V1=posterior_vs_observed$pred_range, V2=rep('1_Posterior-Predictive', 530)))
colnames(interval.range.df) <- c("Prediction_Range", "Prediction")

png("Prior_vs_Posterior_Range.png", units='in', height=3, width=4, res=600)
ggplot(interval.range.df, aes(x=Prediction, y =Prediction_Range)) + geom_boxplot() + theme_classic() + ylab("95% Credible Interval Range") + xlab("Model") # + theme(text = element_text(size = 20)) 
dev.off() 