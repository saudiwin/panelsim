library(tidyverse)
library(panelsim)
require(ggplot2)

# see if DiD works

# loop 1,000 times
# canonical DiD: two time points, interact with time point

true_did_eff <- 5.5
main_coef_x <- -4

# assume treatment effects are homogeneous

first_run <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,did.eff.sd = 0,
                     cross.eff.mean = main_coef_x,case.int.sd = .5,
                     cross.int.sd = .5,cross.eff.sd = 0,
                     case.eff.sd = 0,case.eff.mean = .5,
                     cross.int.mean = .5,case.int.mean = 2)

calc_did <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,
                           cross.eff.mean = main_coef_x,
                           cross.int.mean = .5,case.int.mean = 2,
                           prior_true_vals = first_run$fixed_params)

  did_est <- lm(y ~ x*factor(time), data=did_canonical$data)

  tibble(ATE_DID=did_est$coefficients['x:factor(time)2'],
         ATE_int=did_est$coefficients['factor(time)2'],
         overall_int=did_est$coefficients['(Intercept)'],
         coef_x=did_est$coefficients['x'])

},mc.cores=8) %>% bind_rows

summary(calc_did$ATE_DID)

summary(calc_did$ATE_int)

summary(calc_did$overall_int)

summary(calc_did$coef_x)

did_canonical <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,
                         cross.eff.mean = main_coef_x,
                         cross.int.mean = .5,case.int.mean = 2,
                         prior_true_vals = first_run$fixed_params)

did_est <- lm(y ~ x*factor(time), data=did_canonical$data)

y_pred <- predict(did_est)

summary(y_pred)
summary(did_canonical$data$y)

# matches simulated values

# we can do this with a constant "treatment" effect with more T
# recover coefficient if we still find the coef with a given time point

true_did_eff <- 1.5
main_coef_x <- 5

# run once to generate true values

first_run <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,did.eff.sd = 0,
                     cross.eff.mean = main_coef_x,case.eff.mean = 1,
                     cross.eff.sd=0,case.eff.sd=0,
                     cross.int.mean = .5,case.int.mean = 2)

calc_did <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,
                           cross.eff.mean = main_coef_x,case.eff.sd=0,
                           cross.eff.sd=0,
                           cross.int.mean = .5,case.int.mean = 2,
                           prior_true_vals = first_run$fixed_params)

  did_est <- lm(y ~ x*factor(time), data=did_canonical$data)

  did_est$coefficients['x:factor(time)2']

},mc.cores=8) %>% unlist

summary(calc_did)

# but if we switch to TWFE, no dice:

true_did_eff <- -1.5
main_coef_x <- 5

calc_did_twfe <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=10,did.eff.mean = true_did_eff,
                           cross.eff.mean = main_coef_x,
                           cross.int.mean = .5,case.int.mean = 2)

  did_est <- lm(y ~ x + factor(time) + factor(case), data=did_canonical$data)

  did_est$coefficients['x']

},mc.cores=8) %>% unlist

summary(calc_did_twfe)

## ATE now equal to +2

## Try with heterogenous DiD

true_did_eff <- 2
true_did_sd <- .5
main_coef_x <- -2

# assume treatment effects are homogeneous

calc_did <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=10,did.eff.mean = true_did_eff,
                           did.eff.sd = true_did_sd,
                           cross.eff.mean = main_coef_x,
                           cross.int.mean = .5,case.int.mean = 2)

  did_est <- lm(y ~ x*factor(time), data=did_canonical$data)

  tibble(ATE_DID=did_est$coefficients[grepl(x=names(did_est$coefficients),
                                            pattern='x:factor')],
         coef_x=did_est$coefficients['x'])

},mc.cores=8) %>% bind_rows

# still same average ATE but more variance

summary(calc_did$ATE_DID)

# main effect is he same

summary(calc_did$coef_x)

# try with chatgpt code

compute_xy <- function(a, b, c, d, e, f) {
  denom <- a*c + b - d*f - e
  if (abs(denom) < .Machine$double.eps) {
    stop("Denominator is (nearly) zero — can't compute x,y")
  }
  x <- (d - a) / denom
  y <- (a*c*d - a*d*f - a*e + b*d) / denom
  return(tibble(x = x, y = y))
}

