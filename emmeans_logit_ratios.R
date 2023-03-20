# This script shows how to obtain fold change between probabilities in emmeans
#   in a binomial GLM with a logit link from brms.
source("R/packages.R")

prob_means <- c(
  0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.46, 0.46, 0.46, 0.46, 0.46,
  0.46, 0.46, 0.46, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.45,
  0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51,
  0.51, 0.51, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.40, 0.40,
  0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46,
  0.46, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52, 0.52, 0.54, 0.54, 0.54, 0.54,
  0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60,
  0.60, 0.60, 0.54, 0.60, 0.64
)

set.seed(1)
dat <- expand.grid(
  trials = 50, cat_a = letters[1:12], cat_b = LETTERS[1:8]
) |>
  dplyr::mutate(
    twowaypred = paste0(cat_a, "_", cat_b), probability = prob_means,
    success = rbinom(n(), trials, probability)
  )
binom_mod <- brms::brm(
  success | trials(trials) ~ twowaypred, data = dat,
  family = binomial(link = "logit"), cores = 4, chains = 4, iter = 1e4
)

# The column `ratio` illustrates the desired target calculation.
#   It gives a fold-change between parameters on the data (probability) scale.
# Compare "a" + "A" (intercept) with "a" + "B".
as_draws_df(binom_mod) %>%
  dplyr::transmute(
    eff_1 = plogis(b_Intercept),
    eff_2 = plogis(b_Intercept + b_twowaypreda_B),
    logit_diff = -1 * b_twowaypreda_B,
    pairs_regrid = exp(logit_diff),
    ratio = eff_1 / eff_2
  ) %>%
  purrr::map_dfr(ggdist::median_hdci, .id = "par") %>%
  dplyr::mutate(across(where(is.numeric), \(x)round(x, 3)))

# Now achieve the same outcomes with emmeans
# equivalent to `eff_1` and `eff_2` above
emmeans::emmeans(binom_mod, ~ twowaypred) %>%
  emmeans::regrid() %>%
  data.frame() %>%
  slice(1:2) %>%
  dplyr::mutate(across(where(is.numeric), \(x)round(x, 3)))

# equivalent to `logit_diff` above
emmeans::emmeans(binom_mod, ~ twowaypred) %>%
  pairs %>%
  data.frame() %>%
  slice(1) %>%
  dplyr::mutate(across(where(is.numeric), \(x)round(x, 3)))

# equivalent to `pairs_regrid` above
emmeans::emmeans(binom_mod, ~ twowaypred) %>%
  pairs %>%
  emmeans::regrid() %>%
  data.frame() %>%
  slice(1) %>%
  dplyr::mutate(across(where(is.numeric), \(x)round(x, 3)))

# equivalent to `ratio` above
emmeans::emmeans(binom_mod, ~ twowaypred) %>%
  emmeans::regrid(transform = "log") %>%
  pairs %>%
  emmeans::regrid() %>%
  data.frame() %>%
  slice(1) %>%
  dplyr::mutate(across(where(is.numeric), \(x)round(x, 3)))
