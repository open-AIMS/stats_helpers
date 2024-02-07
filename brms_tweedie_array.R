library(brms)
library(tidyverse)
library(mgcv)
tweedie <- brms::custom_family(
  "tweedie", dpars = c("mu", "phi", "mtheta", "M"),
  links = rep("identity", 4), lb = c(0, 0, 1, 0),
  ub = c(NA, NA, 2, NA), type = "real"
)

stan_funs <- "
  int num_non_zero_fun(array[] real y) {
    int A = 0;
    int N = num_elements(y);
    
    for (n in 1 : N) {
      if (y[n] != 0) {
        A += 1;
      }
    }
    return A;
  }
  
  array[] int non_zero_index_fun(array[] real y, int A) {
    int N = num_elements(y);
    array[A] int non_zero_index;
    int counter = 0;
    for (n in 1 : N) {
      if (y[n] != 0) {
        counter += 1;
        non_zero_index[counter] = n;
      }
    }
    return non_zero_index;
  }
  
  array[] int zero_index_fun(array[] real y, int Z) {
    int N = num_elements(y);
    array[Z] int zero_index;
    int counter = 0;
    for (n in 1 : N) {
      if (y[n] == 0) {
        counter += 1;
        zero_index[counter] = n;
      }
    }
    return zero_index;
  }
  
  void check_tweedie(real mu, real phi, real mtheta) {
    if (mu < 0) {
      reject(\"mu must be >= 0; found mu =\", mu);
    }
    if (phi < 0) {
      reject(\"phi must be >= 0; found phi =\", phi);
    }
    if (mtheta < 1 || mtheta > 2) {
      reject(\"mtheta must be in [1, 2]; found mtheta =\", mtheta);
    }
  }
  
  //void check_tweedie(vector mu, real phi, real mtheta) {
  //  int N = num_elements(mu);
  //  if (phi < 0) {
  //    reject(\"phi must be >= 0; found phi =\", phi);
  //  }
  //  if (mtheta < 1 || mtheta > 2) {
  //    reject(\"mtheta must be in [1, 2]; found mtheta =\", mtheta);
  //  }
  //  
  //  for (n in 1 : N) {
  //    if (mu[n] < 0) {
  //      reject(\"mu must be >= 0; found mu =\", mu[n], \"on element\", n);
  //    }
  //  }
  //}
  
  real tweedie_lpdf(array[] real y, real mu, real phi, real mtheta, int M) {
    check_tweedie(mu, phi, mtheta);
    int N = num_elements(y);
    int N_non_zero = num_non_zero_fun(y);
    int N_zero = N - N_non_zero;
    array[N_non_zero] int non_zero_index = non_zero_index_fun(y, N_non_zero);
    int A = num_elements(non_zero_index);
    int NmA = N - A;
    real lambda = 1 / phi * mu ^ (2 - mtheta) / (2 - mtheta);
    real alpha = (2 - mtheta) / (mtheta - 1);
    real beta = 1 / phi * mu ^ (1 - mtheta) / (mtheta - 1);
    real lp = -NmA * lambda;
    
    for (n in 1 : A) {
      vector[M] ps;
      for (m in 1 : M) {
        ps[m] = poisson_lpmf(m | lambda) + gamma_lpdf(y[non_zero_index[n]] | m * alpha, beta);
      }
      lp += log_sum_exp(ps);
    }
    return lp;
  }
  
  // vector version for mu is untested
  //real tweedie_lpdf(real y, vector mu, real phi, real mtheta, int M) {
  //  check_tweedie(mu, phi, mtheta);
  //  int N = num_elements(y);
  //  int N_non_zero = num_non_zero_fun(y);
  //  int N_zero = N - N_non_zero;
  //  array[N_zero] int zero_index = zero_index_fun(y, N_zero);
  //  array[N_non_zero] int non_zero_index = non_zero_index_fun(y, N_non_zero);
  //  int A = num_elements(non_zero_index);
  //  int NmA = N - A;
  //  vector[N] lambda = 1 / phi * mu ^ (2 - mtheta) / (2 - mtheta);
  //  real alpha = (2 - mtheta) / (mtheta - 1);
  //  vector[N] beta = 1 / phi * mu ^ (1 - mtheta) / (mtheta - 1);
  //  real lp = -sum(lambda[zero_index]);
  //  
  //  for (n in 1 : A) {
  //    vector[M] ps;
  //    for (m in 1 : M) {
  //      ps[m] = poisson_lpmf(m | lambda[n]) + gamma_lpdf(y[non_zero_index[n]] //| m * alpha, beta[n]);
  //    }
  //    lp += log_sum_exp(ps);
  //  }
  //  return lp;
  //}
  
  real tweedie_rng(real mu, real phi, real mtheta) {
    check_tweedie(mu, phi, mtheta);
    
    real lambda = 1 / phi * mu ^ (2 - mtheta) / (2 - mtheta);
    real alpha = (2 - mtheta) / (mtheta - 1);
    real beta = 1 / phi * mu ^ (1 - mtheta) / (mtheta - 1);
    
    int N = poisson_rng(lambda);
    real tweedie_val;
    
    if (mtheta == 1) {
      return phi * poisson_rng(mu / phi);
    }
    if (mtheta == 2) {
      return gamma_rng(1 / phi, beta);
    }
    if (N * alpha == 0) {
      return 0.;
    }
    
    return gamma_rng(N * alpha, beta);
  }
"

stanvars <- brms::stanvar(scode = stan_funs, block = "functions")
set.seed(10)
x <- rnorm(100)
mu <- 20 + x * 5
set.seed(10)
y <- mgcv::rTweedie(mu, p = 1.5, phi = 1)
testdf <- data.frame(x, y)
brms::get_prior(
  y ~ x, data = testdf, family = tweedie, stanvars = stanvars
)
brms::make_stancode(
  y ~ x, data = testdf, family = tweedie, stanvars = stanvars,
  backend = "cmdstanr"
)
fit <- brms::brm(
  y ~ x, data = testdf, family = tweedie, stanvars = stanvars,
  backend = "cmdstanr"
)

# brms::expose_functions(fit, vectorize = TRUE)
