tweedie <- brms::custom_family(
  "tweedie", dpars = c("mu", "phi", "mtheta"),
  links = rep("identity", 3), lb = c(0, 0, 1),
  ub = c(NA, NA, 2), type = "real", loop = FALSE, var = "M"
)

stan_funs <- "
  int num_non_zero_fun(vector y) {
    int A = 0;
    int N = num_elements(y);
    for (n in 1 : N) {
      if (y[n] != 0) {
        A += 1;
      }
    }
    return A;
  }
  
  array[] int non_zero_index_fun(vector y, int A) {
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
  
  array[] int zero_index_fun(vector y, int Z) {
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
  
  void check_tweedie(vector mu, real phi, real mtheta) {
    int N = num_elements(mu);
    if (phi < 0) {
      reject(\"phi must be >= 0; found phi =\", phi);
    }
    if (mtheta < 1 || mtheta > 2) {
      reject(\"mtheta must be in [1, 2]; found mtheta =\", mtheta);
    }
    for (n in 1 : N) {
      if (mu[n] < 0) {
        reject(\"mu must be >= 0; found mu =\", mu[n], \"on element\", n);
      }
    }
  }
  
  real tweedie_lpdf(vector y, vector mu, real phi, real mtheta, int M) {
    check_tweedie(mu, phi, mtheta);
    int N = num_elements(y);
    int N_non_zero = num_non_zero_fun(y);
    int N_zero = N - N_non_zero;
    array[N_zero] int zero_index = zero_index_fun(y, N_zero);
    array[N_non_zero] int non_zero_index = non_zero_index_fun(y, N_non_zero);
    int A = num_elements(non_zero_index);
    int NmA = N - A;
    vector[N] lambda = 1 / phi * mu ^ (2 - mtheta) / (2 - mtheta);
    real alpha = (2 - mtheta) / (mtheta - 1);
    vector[N] beta = 1 / phi * mu ^ (1 - mtheta) / (mtheta - 1);
    real lp = -sum(lambda[zero_index]);
    for (n in 1 : A) {
      vector[M] ps;
      for (m in 1 : M) {
        ps[m] = poisson_lpmf(m | lambda[n]) + gamma_lpdf(y[non_zero_index[n]] | m * alpha, beta[n]);
      }
      lp += log_sum_exp(ps);
    }
    return lp;
  }
  
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

gen_qt <- "
  vector[N] mu = rep_vector(0.0, N);
  vector[N] r_tweedie = rep_vector(0.0, N);
  // generate values for each mu
  mu += Intercept + Xc * b;
  for (n in 1 : N) {
    r_tweedie[n] = tweedie_rng(mu[n], phi, mtheta);
  }
"

tweedie_lpdf <- function(y, mu, phi, mtheta) {
  out <- numeric(length = length(mu))
  for (i in seq_along(out)) {
    out[i] <- tweedie::dtweedie(
      y = y, mu = mu[i], phi = phi[i], power = mtheta[i]
    )
  }
  out
}

log_lik_tweedie <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  mtheta <- prep$dpars$mtheta
  y <- prep$data$Y[i]
  tweedie_lpdf(y, mu, phi, mtheta)
}

tweedie_rng <- function(mu, phi, mtheta) {
  out <- numeric(length = length(mu))
  for (i in seq_along(out)) {
    out[i] <- tweedie::rtweedie(
      n = 1, mu = mu[i], power = mtheta[i], phi = phi[i]
    )
  }
  out
}

posterior_predict_tweedie <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  mtheta <- prep$dpars$mtheta
  tweedie_rng(mu, phi, mtheta)
}

posterior_epred_tweedie <- function(prep) {
  prep$dpars$mu
}

stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
  brms::stanvar(30L, "M", scode = "int<lower=1> M;", block = "data")
