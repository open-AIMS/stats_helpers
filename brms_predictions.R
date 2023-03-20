# This script shows a simple example to illustrate how posterior_predict works
#   under the hood.
source("R/packages.R")

set.seed(1)
dat <- DHARMa::createData(
  sampleSize = 500, intercept = 6, fixedEffects = 0, numGroups = 250,
 randomEffectVariance = 0.2, family = gaussian(), scale = 0.1
)

## BRMS
mod <- brms::brm(
  observedResponse ~ 1 + (1 | group), family = gaussian(), data = dat,
  cores = 3, chains = 3, warmup = 500, thin = 5, iter = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 17)
)

# from brms
options(mc.cores = 1)
response <- dat$observedResponse
npreds <- nrow(as_draws_df(mod))
manual_preds_brms <- matrix(0, npreds, nrow(mod$data))
mean_fit_brms <- as_draws_df(mod) %>%
  dplyr::pull(b_Intercept)
sigma_brms <- as_draws_df(mod) %>%
  dplyr::pull(sigma)
set.seed(10)
for (i in seq_len(ncol(manual_preds_brms))) {
  manual_preds_brms[, i] <- rnorm(
    nrow(manual_preds_brms), mean_fit_brms, sigma_brms
  )
}

options(mc.cores = 1)
set.seed(10)
sim_fit_brms <- brms::posterior_predict(
  mod, re_formula = NA
)
dimnames(sim_fit_brms) <- NULL
all.equal(manual_preds_brms, sim_fit_brms)
