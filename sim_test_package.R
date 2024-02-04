library(tidyverse)
library(panelsim)
require(ggplot2)

set.seed(22902)

# see if DiD works

# loop 1,000 times
# canonical DiD: two time points, interact with time point

true_did_eff <- 2
main_coef_x <- -2

# assume treatment effects are homogeneous

calc_did <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=2,did.eff.mean = true_did_eff,
                           cross.eff.mean = main_coef_x,
                           cross.int.mean = .5,case.int.mean = 2)

  did_est <- lm(y ~ x*factor(time), data=did_canonical$data)

  tibble(ATE_DID=did_est$coefficients['x:factor(time)2'],
         coef_x=did_est$coefficients['x'])

},mc.cores=8) %>% bind_rows

summary(calc_did$ATE_DID)

summary(calc_did$coef_x)

# matches simulated values

# we can do this with a constant "treatment" effect with more T
# recover coefficient if we still find the coef with a given time point

true_did_eff <- -1.5
main_coef_x <- 5

calc_did <- parallel::mclapply(1:1000, function(i) {

  did_canonical <- tw_data(N=1000,T=10,did.eff.mean = true_did_eff,
                           cross.eff.mean = main_coef_x,
                           cross.int.mean = .5,case.int.mean = 2)

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

