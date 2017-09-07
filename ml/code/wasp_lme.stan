functions {
  // takes care of the modified likelihood
  real stoc_approx_log (vector y, vector mu, real stdErr, matrix covMat, matrix zmat, real nrep) {
    matrix[num_elements(y), num_elements(y)] covY; 
    matrix[num_elements(y), num_elements(y)] L;
    
    covY <- quad_form(covMat, zmat');
    for (kk in 1:(rows(zmat)))
      covY[kk, kk] <- covY[kk, kk] + stdErr;

    L <- cholesky_decompose(covY);
    
    return (nrep * multi_normal_cholesky_log(y, mu, L));
  } 
}

data {
  int<lower=0> nobs;                          // total no. of obs
  int<lower=0> nfixef;                        // total no. of fixed effects
  int<lower=0> nranef;                        // total no. of random effects
  int<lower=0> ngroup;                        // total no. of clusters
  real<lower=0> nrep;                         // total no. of resamples
  matrix[nobs, nfixef] xmat;                  // fixed effects design matrix
  matrix[nobs, nranef] zmat;                  // random effects design matrix
  int group[nobs];                            // cluster ids
  vector[nobs] yvec;                          // observations
  int pos1[ngroup];                           // database indices ...
  int pos2[ngroup];                           // to handle ragged arrays
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

  vector[nfixef] fixef;                // population level fix. eff.
  real<lower=0> stdErrFixef;           // std err. in pop. level. fix. eff.
  
  real<lower=0> stdErr;                // population level std. err.
}

transformed parameters {
  matrix[nranef, nranef] covRanef;
  vector[nobs] mu;

  mu <- xmat * fixef;

  covRanef <- quad_form_diag(corrRanef, sclRanef);
}

model {  

  stdErr ~ cauchy(0, 2.5);

  // prior for fix. eff.
  stdErrFixef ~ cauchy(0, 2.5);  
  fixef ~ normal(meanFixef, stdErrFixef);
  
  // prior for rand. eff.
  sclRanef ~ cauchy(0, 2.5);
  corrRanef ~ lkj_corr(2);
  
  for (ii in 1:ngroup) {
    segment(yvec, pos1[ii], pos2[ii] - pos1[ii] + 1) ~ stoc_approx(segment(mu, pos1[ii], pos2[ii] - pos1[ii] + 1), 
                                                                   stdErr, covRanef,
                                                                   block(zmat, pos1[ii], 1, pos2[ii] - pos1[ii] + 1, nranef),
                                                                   nrep);
  }
}
