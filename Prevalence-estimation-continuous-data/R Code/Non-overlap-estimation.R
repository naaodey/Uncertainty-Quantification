library(rstan)
library(ggplot2)
library(loo)
library(tidyr)
library(dplyr)
library(bayesplot)
options(mc.cores = parallel::detectCores())

# spin = shortest posterior interval
# equivalent to  highest posterior density interval for unimodal posteriors

spin = function(x, lower = NULL, upper = NULL, conf = 0.95){
  x = sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n = length(x)
  gap = round(conf*n)
  width = x[(gap+1):n] - x[1:(n-gap)]
  index = min(which(width == min(width)))
  x[c(index, index + gap)]
}


#' @param N number individuals enrolled for the study
#' @param mean shape parameter for unhealthy individuals
#' @param beta_1 scale parameter for unhealthy individuals
#' @param alpha_0 shape parameter for healthy individuals
#' @param beta_0 scale parameter for healthy individuals
#' @param p  prevalence 



# Data Simulation, Non-Overlapping or Minimal lap

data_sim_n = function(N, p, mean = 2, beta_1 = 1, alpha_0 = 10, beta_0 = 1) {
  z = rbinom(N, 1, p)
  
  x = rep(0, N)
  
  for(i in seq(1, N, 1)){
    z = rbinom(1, 1, p)
    
    if(z == 1){
      x[i] = rgamma(1, shape = mean, scale = beta_1 )
    }else{
      x[i] = rgamma(1, shape = alpha_0, scale = beta_0)
    }
  }
  return(x)
}


##################### Simple Expectation-Maximation Algorithm (Classical Approach) To Estimate Prevalence ##########################

em_gam_al_n = function(data_x_n, mean = 2, beta_1 = 1,  alpha_0 = 10, beta_0 = 1, n_iter = 100, tol = 1e-6){
  n = length(data_x_n)
  
  #Setting initial guess for Prevalence
  p = sum(data_x_n > (max(data_x_n) / 3)) / length(data_x_n)
  
  loglik = -Inf
  loglik_trace = numeric(n_iter)
  
  
  for (iter in 1:n_iter){
    # E-step: The likelihood components (Responsiblities)
    f1 = dgamma(data_x_n, shape = mean, scale = beta_1 )
    f2 = dgamma(data_x_n, shape = alpha_0, scale = beta_0)
    
    gamma = p * f1 / (p * f1 + (1 - p) * f2)
    
    # M-step: Updating prevalence
    
    p = mean(gamma)
    
    
    ## Log-likelihood and Convergence Check
    loglik_p = sum(log(p * f1 + (1 - p) * f2))
    
    loglik_trace[iter] = loglik_p
    
    if (abs(loglik_p - loglik) < tol) {
      loglik_trace = loglik_trace[1:iter]
      break
    } 
    loglik = loglik_p
    
  }
  
  return(list(pi = p, likelihood = loglik_p, trace = loglik_trace))
}


# Bootsrap Method to Estimate Matrices

bootstrap_em_gam_n = function(data_x_n, n_sim, p = 0.08) {
  
  n = length(data_x_n)
  
  # Initializing vector for bootstrapped prevalence estimates
  p_boot = numeric(n_sim)
  
  for (j in 1:n_sim) {
    # Bootstrapping: sampling with replacement
    boot_x = sample(data_x_n, n, replace = TRUE)  
    
    # Fitting the model 
    sim_fit = em_gam_al_n(boot_x)
    
    # Storing the prevalence estimate
    p_boot[j] = sim_fit$pi
  }
  
  # Calculating standard error, bias, and confidence intervals
  se = sd(p_boot) 
  bias = mean(p_boot) - p
  
  ci = cbind(spin(p_boot, lower = 0, upper = 1, conf = 0.95))
  
  return(list(
    pi = p_boot,
    se = se,
    bias = bias,
    con_int = ci
  ))
}

# Density plot of mixture model

den_mix_n = function(data_x_n, p, mean = 2, beta_1 = 1, alpha_0 = 10, beta_0 = 1){
  
  
  # First gamma distribution for unhealthy persons
  f1 = dgamma(data_x_n, shape = mean, scale = beta_1)
  
  # Second gamma distribution for healthy persons
  f2 = dgamma(data_x_n, shape = alpha_0, scale = beta_0)
  
  # Mixture density
  
  mix_den = p * f1 + (1 - p) * f2
  
  den_data = data.frame(x = data_x_n, f1 = f1, f2 = f2, mix = mix_den)
  
  # Density Plot
  ggplot(den_data, aes(x = x)) +
    geom_line(aes(y = f1, color = "Component 1"), linewidth = 1) +
    geom_line(aes(y = f2, color = "Component 2"), linewidth = 1.5) +
    geom_line(aes(y = mix, color = "Mixture Density"), linewidth = 1.2) +
    labs(title = "Density of Mixture Model",
         x = "Value",
         y = "Density") +
    scale_color_manual(values = c("Component 1" = "#DAA520",
                                  "Component 2" = "#FF7F50",
                                  "Mixture Density" = "#6A5ACD" )) +
    theme_minimal()
}



################################## Analysis ###############################

set.seed(123)

# Simulated data
data_x_n = data_sim_n(5000, mean = 2, p = 0.08) 

# Fit the model
fit_n = em_gam_al_n(data_x) 

# Fit Bootstrap method
x_boot_n = bootstrap_em_gam_n(data_x_n, 10000)

cat("Prevalence:", mean(x_boot_n$pi), "\n")
cat("Standard Error:", x_boot_n$se, "\n")
cat("Bias:", x_boot_n$bias, "\n")
cat("Confidence Interval:", x_boot_n$con_int[1], "-", x_boot_n$con_int[2], "\n")

# Plot the mixture density
pdf("../img/mixture_density_n.pdf", height = 3.5, width = 4.5)
par(mar=c(3,3,0,1), mgp=c(2, .7, 0), tck=-.02)
den_mix(data_x_n, p = 0.08)  
dev.off()



############################### Simple Bayesian Estimation for Prevalence ###############

# Initial value for Prevalence
ini_val_n = function(){
  list(p = sum(data_x_n > (max(data_x_n) / 3)) / length(data_x_n))
}

# Stan data
bay_data_n = list(N  = length(data_x_n), y = data_x_n,
                alpha_1 = 2, beta_1 = 1,
                alpha_0 = 10, beta_0 = 1)



# Uniform Prior: As original prior
uni_model = stan_model("Uniform.stan")

uni_fit_n = sampling(uni_model, data = bay_data_n, init = ini_val_n, chains = 4,
                   iter = 10000, refresh = 0)

print(uni_fit_n, pars = 'p', digits_summary = 3)

# Posterior draws 
uni_draw_n = rstan::extract(uni_fit_n, pars = 'p')

# 95% highest posterior density bounds
print(spin(uni_draw_n[["p"]], lower = 0, upper = 1, conf = 0.95))



# Jefferys' Prior: First alternative Prior
jef_model = stan_model("Jeff.stan")

jef_fit_n = sampling(jef_model, data = bay_data_n, init = ini_val_n, chains = 4,
                   iter = 10000, refresh = 0)

print(jef_fit_n, pars = 'p', digits_summary = 3)

# Posterior draws 
jef_draw_n = rstan::extract(jef_fit_n, pars = 'p')


# 95% highest posterior density bounds
print(spin(jef_draw_n[["p"]], lower = 0, upper = 1, conf = 0.95))



# MDIP: Second alternative Prior
mdip_model = stan_model("MDIP.stan")

mdip_fit_n = sampling(mdip_model, data = bay_data_n, init = ini_val_n, chains = 4,
                    iter = 10000, refresh = 0)

print(mdip_fit_n, pars = "p", digits_summary = 3)

# Posterior draws 
mdip_draw_n = rstan::extract(mdip_fit_n, pars = "p")

# 95% highest posterior density bounds
print(spin(mdip_draw_n[["p"]], lower = 0, upper = 1, conf = 0.95))


