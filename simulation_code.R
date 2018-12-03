library(tidyr)
library(twowaysim)
setwd("~/Box Sync/Between Effects/Simulation")
source("twsim.R")

set.seed(22902)

## Changing the mean of the within-time slopes
sim1 <- twsim(iter=30, parallel=FALSE, arg="beta.mean", at=-2:5)
sim1 <- gather(sim1, `Two-way FE`:`RE (v_t)`, key="Model", value="Coefficient",
              factor_key=TRUE)

write_csv(sim1, path="sim1.csv")

## Changing the variance of the within-time slopes
sim2 <- twsim(iter=30, parallel=FALSE, arg="beta.sd", at=seq(0,1,by=.1))
sim2 <- gather(sim2, `Two-way FE`:`RE (v_t)`, key="Model", value="Coefficient",
              factor_key=TRUE)

write_csv(sim2, path="sim2.csv")

## Changing the variance of the within-case slopes
sim3 <- twsim(iter=30, parallel=FALSE, arg="gamma.sd", at=seq(0,1,by=.1))
sim3 <- gather(sim3, `Two-way FE`:`RE (v_t)`, key="Model", value="Coefficient",
              factor_key=TRUE)

write_csv(sim3, path="sim3.csv")

## Changing N and T, temporal autocorrelation
iterations <- 500
parallel <- TRUE
sim4 <- data.frame()

simtemp <- twsim(iter=iterations, parallel=parallel, arg="time.ac", at=c(0,.25,.75,.95), N=30, T=30, re_vt=FALSE)
simtemp <- gather(simtemp, `Two-way FE`:`RE (u_i)`, key="Model", value="Coefficient",
              factor_key=TRUE)
simtemp <- mutate(simtemp, type="(30,30)")
sim4 <- bind_rows(sim4, simtemp)

simtemp <- twsim(iter=iterations, parallel=parallel, arg="time.ac", at=c(0,.25,.75,.95), N=100, T=10, re_vt=FALSE)
simtemp <- gather(simtemp, `Two-way FE`:`RE (u_i)`, key="Model", value="Coefficient",
                  factor_key=TRUE)
simtemp <- mutate(simtemp, type="(100,10)")
sim4 <- bind_rows(sim4, simtemp)

simtemp <- twsim(iter=iterations, parallel=parallel, arg="time.ac", at=c(0,.25,.75,.95), N=1000, T=3, re_vt=FALSE)
simtemp <- gather(simtemp, `Two-way FE`:`RE (u_i)`, key="Model", value="Coefficient",
                  factor_key=TRUE)
simtemp <- mutate(simtemp, type="(1000,3)")
sim4 <- bind_rows(sim4, simtemp)

sim4 <- sim4 %>%
  mutate(ac="AC =",
         type = factor(type, levels=c("(30,30)","(100,10)","(1000,3)"))) %>%
  unite(time.ac, ac, time.ac, sep=" ")

write_csv(sim4, path="sim4.csv")


source("sim_august2017_graphs.R")