###################################################
# Simulation files for "Why You Can’t Control Both Space and Time: 
# Reassessing Fixed Effects"
#
# Jonathan Kropko
# jkropko@virginia.edu
#
# Robert Kubinec
# rmk7xy@virginia.edu
# 
# gen_errormat() produces a correlation matrix for the errors. Arguments:
# data: a data frame with variables case and time containing case and time IDs
# time.ac: correlation between errors in adjacent time points for the same case (default = 0) 
# spatial.ac: correlation between errors in adjacent cases at the same time (default = 0) 
#
# twsimdata() generates the simulated data. Arguments:
# N: Number of cases (default = 30)
# T: Number of time points (default = 30)
# alpha.i.mean: Mean of case-specific intercepts (default = 0)
# alpha.i.sd: Standard deviation of case-specific intercepts (default = 1)
# alpha.t.mean: Mean of time-specific intercepts (default = 0)
# alpha.t.sd: Standard deviation of time-specific intercepts (default = 1)
# beta.mean: Mean of the (cross-sectional) slopes of the best fit lines within time points (default = 3)
# beta.sd: SD of the (cross-sectional) slopes of the best fit lines within time points (default = .5)
# gamma.mean: Mean of the (over-time) slopes of the best fit lines within cases (default = -3)
# gamma.sd: SD of the (over-time) slopes of the best fit lines within cases (default = .5)
# noise.sd: SD of the error applied after generation of y (default = 1)
# eta1.i: coefficient for how an omitted variable z_i depends on x (default = 0)
# eta1.t: coefficient for how an omitted variable z_t depends on x (default = 0)
# eta2.i: coefficient for how an omitted variable z_i depends on y (default = 0)
# eta2.t: coefficient for how an omitted variable z_t depends on y (default = 0)
# binary: whether to make the predictor X binary by cutting at the median (default = FALSE)
# unbalance: whether to generate unbalanced panels with time frames 1-K, 
#     where k is drawn from ceiling(runif(min=(T/3), max=T)) (default = FALSE)
# time.ac: correlation between errors in adjacent time points for the same case (default = 0) 
# spatial.ac: correlation between errors in adjacent cases at the same time (default = 0) 
#
# twreg() runs regressions within each case and time to compare
# DGP coefficients to realized coefficients. 
# Arguments: $data element of output of twsimdata()
#
# twmodel() runs pooled OLS, one-way FEs, RE, and two-way FE models
# Arguments: output of twsimdata()
#
# twtransform() creates a second version of the simulated data that
# contains case-level, time-level, and overall means of x and y and 
# the twoway tranformations. 
# Arguments: output of twsimdata()
#
# twtrans2() multiplies the various transformed variables from twtransform()
# for the purpose of investigating the formal one-way and two-way FE estimators.
# Arguments: output of twtransform()
#
# twsim() runs the simulation. Arguments:
# iter: number of iterations
# cores: the number of processing cores to devote to the simulation (defaults to detectCores())
# parallel: whether the simulation should run with parallel processing (default = FALSE)
# ... : arguments passed to twsimdata()
###################################################


#' @import dplyr tidyr
#' @importFrom MASS mvrnorm
#' @export
gen_errormat <- function(data, time.ac, spatial.ac){

      eps <- dplyr::select(data, case, time)
      eps <- mutate(eps, id = 1:n(), one=1)
      eps <- full_join(eps, eps, by="one")
      eps <- eps %>%
            mutate(corr = time.ac^(abs(time.x - time.y))*(case.x==case.y) +
                         spatial.ac^(abs(case.x - case.y))*(time.x==time.y)) %>%
            dplyr::select(id.x, id.y, corr) %>%
            spread(key=id.y, value=corr) %>%
            dplyr::select(-id.x)
      eps <- as.matrix(eps)
      diag(eps) <- 1
      return(eps)
}

#' @export
twsimdata <- function(N = 30, T = 30, alpha.i.mean = 0, alpha.i.sd = 1, 
                      alpha.t.mean = 0, alpha.t.sd = 1, beta.mean = 3, 
                      beta.sd = .5, gamma.mean = -3, gamma.sd = .5, noise.sd = 1, 
                      eta1.i = 0, eta1.t = 0, eta2.i = 0, eta2.t = 0, 
                      binary = FALSE, unbalance = FALSE, time.ac = 0, spatial.ac = 0){

      bet.data <- data.frame(case=1:N, 
                             alpha.i=rnorm(N, mean=alpha.i.mean, sd=alpha.i.sd), 
                             gamma=rnorm(N, mean=gamma.mean, sd=gamma.sd), 
                             zi=rnorm(N), 
                             unbal=ifelse(unbalance, ceiling(runif(N, min=(T/3), max=T)), T),
                             one=1)
      with.data <- data.frame(time=1:T, 
                              alpha.t=rnorm(T, mean=alpha.t.mean, sd=alpha.t.sd), 
                              beta=rnorm(T, mean=beta.mean, sd=beta.sd), 
                              zt=rnorm(T), 
                              one=1)
      data <- full_join(bet.data, with.data, by="one")
      if(time.ac==0 && spatial.ac==0) {
        # if no autocorrelation, much faster to just use rnorm
        data <- data %>%
              mutate(noise=rnorm(n=N*T,sd=noise.sd),
                     x=(alpha.i-alpha.t)/(beta - gamma) + eta1.i*zi + eta1.t*zt, 
                     y=(beta*alpha.i-gamma*alpha.t)/(beta - gamma) + eta2.i*zi + eta2.t*zt + noise,
                     tokeep = time <= unbal) %>%
              dplyr::select(case, time, y, x, gamma, beta, alpha.i, alpha.t, noise, 
                            zi, zt, tokeep) %>%
              arrange(case, time)
      } else {
        eps <- gen_errormat(data, time.ac, spatial.ac)
        data <- data %>%
          mutate(noise=mvrnorm(n=1, mu=rep(0, N*T), Sigma=eps)*noise.sd,
                 x=(alpha.i-alpha.t)/(beta - gamma) + eta1.i*zi + eta1.t*zt, 
                 y=(beta*alpha.i-gamma*alpha.t)/(beta - gamma) + eta2.i*zi + eta2.t*zt + noise,
                 tokeep = time <= unbal) %>%
          dplyr::select(case, time, y, x, gamma, beta, alpha.i, alpha.t, noise, 
                        zi, zt, tokeep) %>%
          arrange(case, time)
      }
      if(binary) data <- mutate(data, x=(x>=median(x)))
      
      pars <- data.frame(N = N, T = T, alpha.i.mean = alpha.i.mean, alpha.i.sd = alpha.i.sd, 
                      alpha.t.mean = alpha.t.mean, alpha.t.sd = alpha.t.sd, beta.mean = beta.mean, 
                      beta.sd = beta.sd, gamma.mean = gamma.mean, gamma.sd = gamma.sd, noise.sd = noise.sd, 
                      eta1.i = eta1.i, eta1.t = eta1.t, eta2.i = eta2.i, eta2.t = eta2.t, 
                      binary = binary, unbalance = unbalance, time.ac = time.ac, spatial.ac = spatial.ac)
      toreturn <- list(data, pars)
      names(toreturn) <- c("data", "pars")
      return(toreturn)
}

#' @export
twreg <- function(d){
      data <- d$data
      pars <- d$pars
      
      #Cross-sectional effects
      data.xs <- split(data, data$time)
      realized <- sapply(data.xs, FUN=function(d){
            reg <- lm(y ~ x, data=d)
            coef(reg)[2]
      }, simplify=TRUE)
      DGP <- dplyr::summarize(group_by(data,time), dgp=mean(beta))
      xs.coef <- cbind(DGP, realized)
      xs.cor <- cor(xs.coef)[2,3]
      xs.plot <- function(x=NULL){
            plot(xs.coef[,2], xs.coef[,3], xlab="DGP", ylab="Realized", pch=16) 
            abline(lm(xs.coef[,3] ~ xs.coef[,2]))
      }
      xs <- list(xs.coef, xs.cor, xs.plot)
      names(xs) <- c("coefs", "cor", "plot")
      
      #Over-time effects
      data.time <- split(data, data$case)
      realized <- sapply(data.time, FUN=function(d){
            reg <- lm(y ~ x, data=d)
            coef(reg)[2]
      }, simplify=TRUE)
      DGP <- dplyr::summarize(group_by(data,case), dgp=mean(gamma))
      time.coef <- cbind(DGP, realized)
      time.cor <- cor(time.coef)[2,3]
      time.plot <- function(x=NULL){
            plot(time.coef[,2], time.coef[,3], xlab="DGP", ylab="Realized", pch=16)
            abline(lm(time.coef[,3] ~ time.coef[,2]))
      }
      time <- list(time.coef, time.cor, time.plot)
      names(time) <- c("coefs", "cor", "plot")
      
      toreturn <- list(xs, time)
      names(toreturn) <- c("cross_section", "over_time")
      return(toreturn)
}

#' @export
twmodel <- function(d, re_vt=TRUE){
      data <- d$data
      pars <- d$pars
      
      require(plm)
      pdata <- data[,apply(data, 2, sd)!=0]
      pdata.time <- pdata.frame(pdata, index=c("case","time"), row.names=FALSE)
      pdata.xs <- pdata.frame(pdata, index=c("time","case"), row.names=FALSE)
      
      #1) Pooled OLS
      pooled.ols <- lm(y ~ x, data=data)
      
      #2) One way FEs
      #case.fe.dummy <- lm(y ~ x + factor(case), data=data)
      #time.fe.dummy <- lm(y ~ x + factor(time), data=data)
      case.fe.plm <- plm(y ~ x, model="within", data=pdata.time)
      time.fe.plm <- plm(y ~ x, model="within", data=pdata.xs)
      
      #3) Random effects
      if(re_vt) random.xs <- plm(y ~ x, model="random", data=pdata.xs)
      random.time <- plm(y ~ x, model="random", data=pdata.time)
      
      #4) Twoway FE
      twoway.fe <- lm(y ~ x + factor(case) + factor(time), data=data)
      identified <- all(!is.na(coef(twoway.fe)))
      twcoef <- ifelse(identified, summary(twoway.fe)$coefficients[2,1], NA)
      twse <- ifelse(identified, summary(twoway.fe)$coefficients[2,2], NA)
      
      ### Compile saved results in a data frame  
      results <- data.frame(
            model = c("Two-way FE", "Case FE", "Time FE", "Pooled OLS",
                      "RE (u_i)"),
            coef = c(twcoef,
                     #summary(case.fe.dummy)$coefficients[2,1],
                     summary(case.fe.plm)$coefficients[1,1],
                     #summary(time.fe.dummy)$coefficients[2,1],
                     summary(time.fe.plm)$coefficients[1,1],
                     summary(pooled.ols)$coefficients[2,1],
                     summary(random.time)$coefficients[2,1]),
            se = c(twse,
                  # summary(case.fe.dummy)$coefficients[2,2],
                   summary(case.fe.plm)$coefficients[1,2],
                   #summary(time.fe.dummy)$coefficients[2,2],
                   summary(time.fe.plm)$coefficients[1,2],
                   summary(pooled.ols)$coefficients[2,2],
                   summary(random.time)$coefficients[2,2])
      )
      if(re_vt){
        results_vt <- data.frame(model = "RE (v_t)", 
                                 coef = summary(random.xs)$coefficients[2,1],
                                 se = summary(random.xs)$coefficients[2,2])
        results <- bind_rows(results, results_vt)
      }
      return(results)
}

#' @export
twtransform <- function(d){
      require(dplyr)
      pars <- d$pars
      data <- d$data
      if(pars$unbalance) data <- filter(data, tokeep)
      data <- dplyr::select(data, case, time, x, y)
      data.case <- data %>%
            group_by(case) %>%
            summarize(x.bar.i=mean(x), y.bar.i=mean(y)) 
      data <- full_join(data, data.case, by="case")
      data.time <- data %>%
            group_by(time) %>%
            summarize(x.bar.t=mean(x), y.bar.t=mean(y)) 
      data <- full_join(data, data.time, by="time") 
      data2 <- data %>%
            mutate(x.bar=mean(x), y.bar=mean(y),
                  x.pool = x - x.bar,
                  y.pool = y - y.bar,
                  x.caseFE = x - x.bar.i,
                  y.caseFE = y - y.bar.i,
                  x.timeFE = x - x.bar.t,
                  y.timeFE = y - y.bar.t,
                  x.tw = x - x.bar.i - x.bar.t + x.bar,
                  y.tw = y - y.bar.i - y.bar.t + y.bar) 
      return(data2)
}

#' @export
twtrans2 <- function(data){
      yc.xc <- with(data, sum(y.caseFE*x.caseFE))
      yc.xt <- with(data, sum(y.caseFE*x.timeFE))
      yc.xp <- with(data, sum(y.caseFE*x.pool))
      yt.xc <- with(data, sum(y.timeFE*x.caseFE))
      yt.xt <- with(data, sum(y.timeFE*x.timeFE))
      yt.xp <- with(data, sum(y.timeFE*x.pool))
      yp.xc <- with(data, sum(y.pool*x.caseFE))
      yp.xt <- with(data, sum(y.pool*x.timeFE))
      yp.xp <- with(data, sum(y.pool*x.pool))
      tw.num <- yc.xc + yc.xt - yc.xp + yt.xc + yt.xt -
            yt.xp - yp.xc - yp.xt + yp.xp
      
      xc.xc <- with(data, sum(x.caseFE*x.caseFE))
      xt.xt <- with(data, sum(x.timeFE*x.timeFE))
      xp.xp <- with(data, sum(x.pool*x.pool))      
      xc.xt <- with(data, sum(x.caseFE*x.timeFE))
      xc.xp <- with(data, sum(x.caseFE*x.pool))
      xt.xp <- with(data, sum(x.timeFE*x.pool))
      tw.denom <- xc.xc + xt.xt + xp.xp + 2*xc.xt - 2*xc.xp - 2*xt.xp
      
      tw <- tw.num/tw.denom

      p <- c("yc.xc"=yc.xc, "yc.xt"=yc.xt, "yc.xp"=yc.xp, 
             "yt.xc"=yt.xc, "yt.xt"=yt.xt, "yt.xp"=yt.xp, "yp.xc"=yp.xc, "yp.xt"=yp.xt, 
             "yp.xp"=yp.xp, "tw.num"=tw.num, "xc.xc"=xc.xc, "xt.xt"=xt.xt, "xp.xp"=xp.xp, 
             "xc.xt"=xc.xt, "xc.xp"=xc.xp, "xt.xp"=xt.xp, "tw.denom"=tw.denom, "tw"=tw )
      return(p)
}

#' @export
twsim <- function(iter=1000, cores=NULL, parallel=FALSE, arg, at, verbose=TRUE, re_vt=TRUE, ...){
      el <- expand.grid(iteration = 1:iter, arg = at)
      if(parallel){
            if(file.exists('output.txt')) file.remove('output.txt')
            loopfun <- ifelse(verbose, pbmclapply::pbmclapply, parallel::mclapply)
            if(is.null(cores)) cores <- parallel::detectCores()
            p <- loopfun(1:nrow(el), FUN=function(x){
                  sink('output.txt',append = T)
                  print(paste0('Now on row ',x, ' out of ',nrow(el),' rows in the simulation.'))
                  sink()
                  z <- list(foo = el[x,2])
                  names(z) <- arg
                  extra.args <- as.list(substitute(list(...)))[-1L]
                  z<-c(z,extra.args)
                  d <- do.call(twsimdata, args=z)
                  m <- twmodel(d, re_vt=re_vt)
                  m <- sapply(m, FUN=function(x){x}, simplify=TRUE)
                  return(c(el[x,2], m[,2]))
            }, mc.cores=cores) 
            p <- sapply(p, FUN=function(x){x}, simplify=TRUE)
      } else {
            p <- sapply(1:nrow(el), FUN=function(x){
                  if(verbose) print(el[x,])
                  z <- list(foo = el[x,2])
                  extra.args <- as.list(substitute(list(...)))[-1L]
                  z<-c(z,extra.args)
                  names(z) <- arg
                  d <- do.call(twsimdata, args=z)
                  m <- twmodel(d, re_vt=re_vt)
                  m <- sapply(m, FUN=function(x){x}, simplify=TRUE)
                  return(c(el[x,2], m[,2]))      
            })
      }
      p <- data.frame(t(p))
      if(re_vt){
        colnames(p) <- c(arg, "Two-way FE", "Case FE", "Time FE", "Pooled OLS","RE (u_i)", "RE (v_t)")
      } else {
        colnames(p) <- c(arg, "Two-way FE", "Case FE", "Time FE", "Pooled OLS","RE (u_i)")
      }
      return(p)
}


