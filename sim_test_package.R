library(tidyverse)
library(panelsim)
require(ggplot2)
set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- tw_sim(iter=100, parallel=T, arg=c('omm.x.case',
                                          'case.int.sd'),at1=-2:2,
               at2=1:5,cores = 2)

sim1 <- tw_sim(N=150,T=150,iter=30, parallel=T, arg="case.eff.mean", at1=-2:5,cores=4,cross.eff.mean=0,cross.eff.sd=0)

## plot simulation

tw_plot(sim1) + geom_point(alpha=0.5)
