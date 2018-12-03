library(tidyverse)

#box_sync_path <- normalizePath('C:/Users/bobku/Box Sync/')

#setwd(paste0(box_sync_path,'Between Effects/Simulation'))

sim1 <- read_csv("sim1.csv") %>% mutate(beta.mean=as.numeric(beta.mean))
sim2 <- read_csv("sim2.csv") %>% mutate(beta.mean=as.numeric(beta.mean))
sim3 <- read_csv("sim3.csv") %>% mutate(beta.mean=as.numeric(beta.mean))

sim1 <- mutate(sim1, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))
sim2 <- mutate(sim2, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))
sim3 <- mutate(sim3, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))


pdf("simresults_august2017.pdf", width=8, height=5)

## Changing the mean of the within-time slopes
g1 <- ggplot(sim1, aes(x=beta.mean, y=Coefficient)) +
  geom_jitter() +
  geom_smooth(method="loess") +
  geom_line(aes(x = beta.mean, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = beta.mean, y = beta.mean, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Mean of the within-time slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (varying)     ", "Mean slope within cases (fixed at -3)")) 
g1


## Changing the variance of the within-time slopes
g2 <- ggplot(sim2, aes(x=beta.sd, y=Coefficient)) +
  geom_jitter() +
  geom_smooth(method="loess") +
  geom_line(aes(x = beta.sd, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = beta.sd, y = 3, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Standard deviation of the within-time slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (fixed at 3)     ", "Mean slope within cases (fixed at -3)"))
g2

## Changing the variance of the within-case slopes
g3 <- ggplot(sim3, aes(x=gamma.sd, y=Coefficient)) +
  geom_jitter() +
  geom_smooth(method="loess") +
  geom_line(aes(x = gamma.sd, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = gamma.sd, y = 3, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Standard deviation of the within-case slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (fixed at 3)     ", "Mean slope within cases (fixed at -3)"))
g3

dev.off()

sim1 <- read_csv("sim1big.csv")
sim2 <- read_csv("sim2big.csv")
sim3 <- read_csv("sim3big.csv")

sim1 <- mutate(sim1, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))
sim2 <- mutate(sim2, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))
sim3 <- mutate(sim3, Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)", "RE (v_t)")))


pdf(paste0(box_sync_path,"Between Effects/Paper/Manuscript/simresults_august2017big.pdf"), width=8, height=5)

## Changing the mean of the within-time slopes
g1 <- ggplot(sim1, aes(x=beta.mean, y=Coefficient)) +
  geom_jitter(size=.5,alpha=.5) +
  geom_smooth(method="loess") +
  geom_line(aes(x = beta.mean, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = beta.mean, y = beta.mean, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Mean of the within-time slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (varying)     ", "Mean slope within cases (fixed at -3)")) 
g1


## Changing the variance of the within-time slopes
g2 <- ggplot(sim2, aes(x=beta.sd, y=Coefficient)) +
  geom_jitter(size=.5,alpha=.5) +
  geom_smooth(method="loess") +
  geom_line(aes(x = beta.sd, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = beta.sd, y = 3, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Standard deviation of the within-time slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (fixed at 3)     ", "Mean slope within cases (fixed at -3)"))
g2

## Changing the variance of the within-case slopes
g3 <- ggplot(sim3, aes(x=gamma.sd, y=Coefficient)) +
  geom_jitter(size=.5,alpha=.5) +
  geom_smooth(method="loess") +
  geom_line(aes(x = gamma.sd, y = -3, linetype="solid"), color="red") +
  geom_line(aes(x = gamma.sd, y = 3, linetype="dashed"), color="red") +
  facet_wrap(~ Model) +
  xlab("Standard deviation of the within-case slopes in the DGP") +
  theme(legend.position='bottom') +
  scale_linetype_manual(name = "",
                        values = c("dashed" = 2, "solid" = 1), 
                        labels = c("Mean slope within time points (fixed at 3)     ", "Mean slope within cases (fixed at -3)"))
g3

dev.off()

sim4 <- read_csv("sim4.csv")
sim4 <- mutate(sim4, type = factor(type, 
                                   levels=c("(30,30)","(100,10)","(1000,3)")),
               Model = factor(Model, 
                                    levels=c("Two-way FE", "Case FE",
                                             "Time FE", "Pooled OLS",
                                             "RE (u_i)")))

pdf(paste0(box_sync_path,"Between Effects/Paper/Manuscript/simresults_august2017_2.pdf"), width=9, height=7)
## Changing N and T, temporal autocorrelation
g4 <- ggplot(sim4, aes(y=Coefficient, x=factor(type))) +
  geom_boxplot() + 
  facet_grid(time.ac~Model) + 
  theme_bw() +
  xlab("Sample Sizes (N, T)") + 
  ylab("Coefficients") +
  scale_y_continuous(breaks=c(-3,0,3))
g4

dev.off()

