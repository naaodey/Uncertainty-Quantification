data {
  int<lower=1> N;                     // number of observations
  vector[N] y;                        // observed test measures
  real<lower=0> alpha_0;
  real<lower=0> beta_0;
  real<lower=0> alpha_1;
  real<lower=0> beta_1;
}

parameters {
  real<lower=0, upper=1> p;        // prevalence
}

model {
  
  
  p ~ uniform (0,1);
  
  for (i in 1:N) {
    target += log_mix(
      p,
      gamma_lpdf(y[i] | alpha_1, beta_1),
      gamma_lpdf(y[i] | alpha_0, beta_0)
    );
  }
}


generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  for (i in 1:N) {
    real lp_pos = log(p)   + gamma_lpdf(y[i] | alpha_1, beta_1);
    real lp_neg = log1m(p) + gamma_lpdf(y[i] | alpha_0, beta_0);

    // log-likelihood for LOO / WAIC
    log_lik[i] = log_sum_exp(lp_pos, lp_neg);

    // posterior predictive draw
    if (bernoulli_rng(p) == 1)
      y_rep[i] = gamma_rng(alpha_1, beta_1);
    else
      y_rep[i] = gamma_rng(alpha_0, beta_0);
  }
}

