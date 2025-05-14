library(tidyverse)
library(panelsim)
require(ggplot2)

# see if DiD works

# loop 1,000 times
# canonical DiD: two time points, interact with time point

true_did_eff <- 0.2
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
                           prior_true_vals = first_run,
                           binary_outcome=F,binary_x=F)

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

# try with treatment effects, look at what happens to the 1-way case FE


