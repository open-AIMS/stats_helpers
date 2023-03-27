rm(list=ls())

library(dplyr)
library(INLA)
source("R/dharma.R")

set.seed(1)
dat <- DHARMa::createData(
  sampleSize = 500, intercept = 6, fixedEffects = 0, numGroups = 250,
  randomEffectVariance = 0.2, family = gaussian(), scale = 0.1
)


prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))

form <- observedResponse ~ 1 + f(group, model='iid', hyper = prior.prec) 

mod <- inla(form,
                 data=dat,
                 control.predictor=list(link=1, compute=TRUE),
                 control.compute = list(return.marginals.predictor=TRUE,
                                        config = TRUE, dic= TRUE), 
                 verbose = FALSE)

####################################################### Posterior draws
n_draws <- 1000
sim <- inla.posterior.sample(n_draws, result=mod, seed=123)

response <- dat$observedResponse
manual_preds_inla <- matrix(0, n_draws, nrow(dat))

means = sapply(sim, function(x) x[['latent']])

i.mod <- sapply(c('Intercept', '^Predictor'),
                function(x) grep(x, sim[[1]]$latent %>% rownames))

mean_fit_inla <- means[i.mod[[1]],] %>%  as.data.frame

# INLA uses precision not standard deviation

sigma_inla <- unlist(lapply(sim,
                            function(s) 1/sqrt(sim[[1]]$hyperpar["Precision for the Gaussian observations"])))

set.seed(10)
for (i in seq_len(ncol(manual_preds_inla))) {
  manual_preds_inla[, i] <- rnorm(
    nrow(manual_preds_inla), mean_fit_inla$., sigma_inla
  )
}

set.seed(10)
sim_fit_inla <- means[i.mod[[2]],] 

dimnames(sim_fit_inla) <- NULL
all.equal(t(manual_preds_inla), sim_fit_inla)

#### pp_check
pp_check.foo <- function(object, type = c("multiple", "overlaid"), ...) {
  type <- match.arg(type)
  y <- object[["y"]]
  yrep <- object[["yrep"]]

  if (type == "overlaid") {
    ppc_dens_overlay(y, yrep, ...) 
  } else {
    ppc_hist(y, yrep, ...)
  }
}

#convert to class foo for pp_check
x <- list(y = response, yrep = t(sim_fit_inla))
class(x) <- "foo"
pp_check.foo(x, type = "overlaid")

# DHARMa using posterior predictive distributions 

fitted_median_inla <- apply(sim_fit_inla, 1, median)

DHARMa_inla_posterior <- DHARMa::createDHARMa(
  simulatedResponse = sim_fit_inla,
  observedResponse = response,
  fittedPredictedResponse = fitted_median_inla,
)

a <- gg_dharma(DHARMa_inla_posterior, form = factor(rep(1, nrow(dat))))
b <- gg_disp_hist(DHARMa_inla_posterior)
a / b

# DHARMa using manual posterior predictive distributions 

fitted_median_manual_inla <- apply(t(manual_preds_inla), 1, median)

DHARMa_inla_posterior_manual <- DHARMa::createDHARMa(
  simulatedResponse = t(manual_preds_inla),
  observedResponse = response,
  fittedPredictedResponse = fitted_median_manual_inla,
)

a <- gg_dharma(DHARMa_inla_posterior_manual, form = factor(rep(1, nrow(dat))))
b <- gg_disp_hist(DHARMa_inla_posterior_manual)
a / b
