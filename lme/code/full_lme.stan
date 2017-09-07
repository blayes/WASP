data {
  // data for model fitting
  int<lower=0> nobs;                          // total no. of individuals
  int<lower=0> nfixef;                        // total no. of fixed effects
  int<lower=0> nranef;                        // total no. of random effects
  int<lower=0> ngroup;                        // total no. of clusters
  matrix[nobs, nfixef] xmat;                  // fixed effects design matrix
  matrix[nobs, nranef] zmat;                  // random effects design matrix
  int group[nobs];                            // cluster ids
  real yvec[nobs];                            // observations
}

transformed data {
  // both fix eff. and rand. eff. are apriori centered at 0
  vector[nranef] meanRanef;
  vector[nfixef] meanFixef;
  
  meanRanef <- rep_vector(0.0, nranef);
  meanFixef <- rep_vector(0.0, nfixef);
}

parameters {
  corr_matrix[nranef] corrRanef;       // correlation matrix of rand. eff.
  vector<lower=0>[nranef] sclRanef;    // scale matrix of rand. eff.
  vector[nranef] ranef[ngroup];        // population level rand. eff.

  vector[nfixef] fixef;                // population level fix. eff.
  real<lower=0> stdErrFixef;           // std err. in pop. level. fix. eff.

  real<lower=0> stdErr;                // population level std. err.
}

transformed parameters {
  real yHat[nobs];
  matrix[nranef, nranef] covRanef;

  for (ii in 1:nobs) {
    yHat[ii] <- xmat[ii] * fixef + zmat[ii] * ranef[group[ii]]; // individual level mean 
  }

  covRanef <- quad_form_diag(corrRanef, sclRanef);
}

model {
  stdErr ~ cauchy(0, 2.5);
  
  // sample rand. eff.
  sclRanef ~ cauchy(0, 2.5);
  corrRanef ~ lkj_corr(2);
  for (ii in 1:ngroup) {
    ranef[ii] ~ multi_normal(meanRanef, covRanef); 
  }
  
  // sample fix. eff.
  stdErrFixef ~ cauchy(0, 2.5);  
  fixef ~ normal(meanFixef, stdErrFixef);

  // sample data
  yvec ~ normal(yHat, stdErr); 
}
