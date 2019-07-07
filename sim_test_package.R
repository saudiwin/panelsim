library(tidyverse)
library(twowaysim)
set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- tw_sim(iter=10, parallel=T, arg=c('omm.x.case',
                                          'case.int.sd'),at1=-2:2,
               at2=1:5,cores = 6)

## plot simulation

tw_plot(sim1)
