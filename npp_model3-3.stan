// Stan code for model:

// data block - includes prior hyperparameters
data {
  // data:
  // n params
  int<lower=0> J; // number of types - expecting 10
  int<lower=0> K; // number of zones - expecting 16
  int<lower=0> L; // number of years - expecting 21
  
  // median data - response to fit
  real dat_mat[J,K,L]; // median data
  
  // MU DATA
  // mu intercept
  real mu_intercept;
  // mu TYPE scale coefficient
  real type_coefs_scale;
  // mu ZONE scale coefficient
  real zone_coefs_scale;
  // mu YEAR scale coefficient
  real year_coefs_scale;
  
  // SD DATA
  // sd intercept
  real sd_intercept;
}

parameters {
  // MU
  real std_type_coefs[J];
  real std_zone_coefs[K];
  real std_year_coefs[L];
  // SD
  real std_sd_type_coefs[J];
  real std_sd_zone_coefs[K];
}

// transformed parameters block - happens AFTER model block in TIME
// but comes first in stan code

// we draw standardized parameters in model block - need to transform them
// using scale params and add intercept values before we can use them
// to estimate aggregated median NPP vals
transformed parameters {
  // mu coefs
  real year_coefs[L];
  real zone_coefs[K];
  real type_coefs[J];
  // mu matrix
  real mu_mat[J, K, L];
  // sd matrix
  real sd_mat[J, K];
  
  // zone loop
  for(k in 1:K) {
    // transform norm_zone_coefs
    zone_coefs[k] = (std_zone_coefs[k]*zone_coefs_scale);
    // type loop
    for(j in 1:J) {
      // transform norm_type_coefs
      type_coefs[j] = std_type_coefs[j]*type_coefs_scale;
      // transform norm_sd_type_coefs
      sd_mat[j, k] = exp(sd_intercept + std_sd_zone_coefs[k] + std_sd_type_coefs[j]);
      // year loop
      for(l in 1:L) {
        // transform norm_year_coefs
        year_coefs[l] = std_year_coefs[l]*year_coefs_scale;
        mu_mat[j,k,l] = mu_intercept + zone_coefs[k] + type_coefs[j] + year_coefs[l];
      }
    }
  }
}

// model block - steps to model eventual response of interest
model {
  // zone loop
  for(k in 1:K) {
    std_zone_coefs[k] ~ normal(0, 1);
    std_sd_zone_coefs[k] ~ normal(0, 1);
    // type loop
    for(j in 1:J) {
      std_type_coefs[j] ~ normal(0, 1);
      std_sd_type_coefs[j] ~ normal(0, 1);
      // year loop
      for(l in 1:L) {
        std_year_coefs[l] ~ normal(0, 1);
        if(dat_mat[j,k,l] < 0) {continue;}
        dat_mat[j, k, l] ~ normal(mu_mat[j,k,l], sd_mat[j, k]);
      }
    }
  }
}

generated quantities {
  // initialize generated quantities :
  real post_pred[J, K, L]; // array to store posterior predictions of observed data
  // loop for years
  for(l in 1:L) {
    // loop for zones
    for(k in 1:K) {
      // loop for types
      for(j in 1:J) {
        // if we don't actually have data, continue
        if(dat_mat[j,k,l] < 0) {continue;}
        post_pred[j,k,l] = normal_rng(mu_mat[j,k,l], sd_mat[j,k]);
      }
    }
  }
}

