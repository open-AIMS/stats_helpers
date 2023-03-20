source("R/packages.R")
source("R/dharma.R")

set.seed(1)
dat <- DHARMa::createData(
  sampleSize = 500, intercept = 6, fixedEffects = 0, numGroups = 250,
 randomEffectVariance = 0.2, family = gaussian(), scale = 0.1
)

## BRMS
gaus_mod_brms <- brms::brm(
  observedResponse ~ 1 + (1 | group), family = gaussian(), data = dat,
  cores = 3, chains = 3, warmup = 500, thin = 5, iter = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 17)
)

brms_dharma_res <- make_brms_dharma_res(
  gaus_mod_brms, integerResponse = FALSE
)

plot(brms_dharma_res, form = factor(rep(1, nrow(dat))))

a <- gg_dharma(brms_dharma_res, form = factor(rep(1, nrow(dat))))
b <- gg_disp_hist(brms_dharma_res)
a / b

## glmmTMB
gaus_mod_glmmTMB <- glmmTMB::glmmTMB(
  observedResponse ~ 1 + (1 | group), family = gaussian(), data = dat
)
glmmTMB_dharma_res <- DHARMa::simulateResiduals(gaus_mod_glmmTMB)
plot(glmmTMB_dharma_res)
a <- gg_dharma(glmmTMB_dharma_res, form = factor(rep(1, nrow(dat))))
b <- gg_disp_hist(glmmTMB_dharma_res)
a / b

## glmmTMB
gaus_mod_lme4 <- lme4::lmer(
  observedResponse ~ 1 + (1 | group), data = dat
)
lmer_dharma_res <- DHARMa::simulateResiduals(gaus_mod_lme4, use.u = FALSE)
plot(lmer_dharma_res)
a <- gg_dharma(lmer_dharma_res, form = factor(rep(1, nrow(dat))))
b <- gg_disp_hist(lmer_dharma_res)
a / b
