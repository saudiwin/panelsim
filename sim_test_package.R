library(tidyverse)
library(panelsim)
require(ggplot2)
set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- tw_sim(iter=100, parallel=T, arg=c('omm.x.case',
                                          'case.int.sd'),at1=-2:2,
               at2=1:5,cores = 2)

sim1 <- tw_sim(N=250,T=50,iter=30, parallel=T, arg="case.eff.mean", range1=c(-2,5),cores=6,cross.eff.mean=0,cross.eff.sd=0)

## plot simulation

tw_plot(sim1,facet_scales = "free") + geom_point(alpha=0.5)

sim_bin <- tw_sim(N=50,T=50,iter=30, parallel=T, arg="case.eff.mean", range1=c(0,0.5),cores=6,cross.eff.mean=0,cross.eff.sd=0,
                  models=c("randomfx","wfe"),
                  treat_effect=T,
                  binary_outcome=T)

test_bin <- tw_data(N=100,T=100,case.eff.mean = -0.1,treat_effect = T,binary_outcome = T)

## plot simulation

tw_plot(sim1,facet_scales = "free_y") + geom_point(alpha=0.5)
