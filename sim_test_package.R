library(tidyverse)
library(twowaysim)
set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- tw_sim(iter=10, parallel=T, arg=c('omm.x.case'),at=-2:2,cores = 3)

## plot simulation

tw_plot(sim1)
