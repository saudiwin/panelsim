library(tidyverse)
library(twowaysim)
set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- tw_sim(iter=30, parallel=T, arg="case.eff.mean", at=-2:5,cores = 3)
