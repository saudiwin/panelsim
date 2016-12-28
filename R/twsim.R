###################################################
# Simulation files for "On the Simultaneous Use
# of Fixed Effects on Cases and Time Points"
#
# Jonathan Kropko
# jkropko@virginia.edu
#
# Robert Kubinec
# rmk7xy@virginia.edu
#
# Goals of the simulation:
# (1) Demonstrate that while the twoway FE estimator
# does in fact remove omitted variables that are fixed
# across cases or time points, it reintroduces a model
# misspecification bias regarding the dimension (cases
# or time) that the effect should be present on, a bias
# that oneway FEs correct
#
# (2) Demonstrate the the twoway FE estimator is
# sensitive to balance in the panels while oneway
# FEs are not.
#
# twsimdata() generates the simulated data. Arguments:
# N: Number of cases
# T: Number of time points
# betsd: Standard deviation of case unit effects
# withsd: Standard deviation of time unit effects
# totalsd: Standard deviation of the overall residual
# alpha: Intercept
# betaBmean: Mean of the cross-sectional effects
# betaBsd: Standard deviation over time of cross-sectional effects
# betaWmean: Mean of the over-time effects
# betaWsd: Standard deviation across cases of over-time effects
# binary: whether to make the predictor X binary by cutting at the median
# unbalance: whether to generate unbalanced panels with time frames 1-K,
#     where k is drawn from ceiling(runif(min=(T/3), max=T))
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
# twsim() runs the simulation. Arguments:
# iter: number of iterations
# cores: the number of processing cores to devote to the simulation (defaults to detectCores())
# ... : arguments passed to twsimdata()
###################################################

#' @import dplyr
#' @export
twsimdata <- function(N=30, T=30, ui.mean=0, ui.sd=1, vt.mean=0, vt.sd=1,
                      ep.ui.sd=.1, ep.vt.sd=.1, betaBmean=3, betaBsd=0, betaWmean=-3, betaWsd=0,
                      binary=FALSE, unbalance=FALSE, gamma.ux=0, gamma.vx=0, gamma.uy=0, gamma.vy=0){



      bet.data <- data.frame(case=1:N, ui=rnorm(N, mean=ui.mean, sd=ui.sd),
                             betaW=rnorm(N, mean=betaWmean, sd=betaWsd), zi=rnorm(N))
      with.data <- data.frame(time=1:T, vt=rnorm(T, mean=vt.mean, sd=vt.sd),
                              betaB=rnorm(T, mean=betaBmean, sd=betaBsd), zt=rnorm(T))
      if(unbalance){
            bet.data$unbal <- ceiling(runif(N, min=(T/3), max=T))
            data <- merge(bet.data, with.data)
      } else data <- merge(bet.data, with.data)
      data$ep.ui <- rnorm(nrow(data), mean=0, sd=ep.ui.sd)
      data$ep.vt <- rnorm(nrow(data), mean=0, sd=ep.vt.sd)

      data <- mutate(data, x=(vt-ui)/(betaW - betaB) +
                           (ep.vt - ep.ui)/(betaW-betaB) + gamma.ux*zi + gamma.vx*zt,
                     y=(betaW*vt-betaB*ui)/(betaW - betaB) +
                           (betaW*ep.vt - betaB*ep.ui)/(betaW-betaB) + gamma.uy*zi + gamma.vy*zt)
      if(binary) data$x <- data$x >= median(data$x)
      if(unbalance) data$I <- (data$time <= data$unbal)
      data <- arrange(data, case, time)
      pars <- data.frame(N=N, T=T, ui.mean=ui.mean, ui.sd=ui.sd,
                         vt.mean=vt.mean, vt.sd=vt.sd,
                         ep.ui.sd=ep.ui.sd, ep.vt.sd=ep.vt.sd,
                         betaBmean=betaBmean, betaBsd=betaBsd,
                         betaWmean=betaWmean, betaWsd=betaWsd,
                         binary=binary, unbalance=unbalance)
      toreturn <- list(data,pars)
      names(toreturn) <- c("data", "pars")
      return(toreturn)
}
#' @export
twreg <- function(d){
      data <- d$data
      pars <- d$pars

      #Cross-sectional effects
      data.bet <- split(data, data$time)
      realized <- sapply(data.bet, FUN=function(d){
            reg <- lm(y ~ x, data=d)
            coef(reg)[2]
      }, simplify=TRUE)
      DGP <- summarize(group_by(data,time), dgp=mean(betaB))
      bet.coef <- cbind(DGP, realized)
      bet.cor <- cor(bet.coef)[2,3]
      bet.plot <- function(x=NULL){
            plot(bet.coef[,2], bet.coef[,3], xlab="DGP", ylab="Realized", pch=16)
            abline(lm(bet.coef[,3] ~ bet.coef[,2]))
      }
      bet <- list(bet.coef, bet.cor, bet.plot)
      names(bet) <- c("coefs", "cor", "plot")

      #Over-time effects
      data.with <- split(data, data$case)
      realized <- sapply(data.with, FUN=function(d){
            reg <- lm(y ~ x, data=d)
            coef(reg)[2]
      }, simplify=TRUE)
      DGP <- summarize(group_by(data,case), dgp=mean(betaW))
      with.coef <- cbind(DGP, realized)
      with.cor <- cor(with.coef)[2,3]
      with.plot <- function(x=NULL){
            plot(with.coef[,2], with.coef[,3], xlab="DGP", ylab="Realized", pch=16)
            abline(lm(with.coef[,3] ~ with.coef[,2]))
      }
      with <- list(with.coef, with.cor, with.plot)
      names(with) <- c("coefs", "cor", "plot")

      toreturn <- list(bet, with)
      names(toreturn) <- c("cross_section", "over_time")
      return(toreturn)
}

#' @import plm
#' @export
twmodel <- function(d){

      data <- d$data
      pars <- d$pars

      #1) Pooled OLS
      pooled.ols <- lm(y ~ x, data=data)

      #2) One way FEs
      case.fe.dummy <- lm(y ~ x + factor(case), data=data)
      time.fe.dummy <- lm(y ~ x + factor(time), data=data)
      #case.fe.plm <- plm(y ~ x, model="within", data=pdata.with)
      #time.fe.plm <- plm(y ~ x, model="within", data=pdata.bet)

      #3) Random effects
      pdata.with <- pdata.frame(data, index=c("case","time"), row.names=FALSE)
      pdata.bet <- pdata.frame(data, index=c("time","case"), row.names=FALSE)
      random.ui <- plm(y ~ x, model="random", data=pdata.with)
      #random.vt <- plm(y ~ x, model="random", data=pdata.bet)

      #4) Twoway FE
      ##Manually leaving out a case (combine first and last cases)
      case2 <- factor(data$case)
      levels(case2)[pars$N] <- levels(case2)[1]
      tw.dropcase <- lm(y ~ x + case2 + factor(time), data=data)

      ##Manually leaving out a time (combine first and last time points)
      time2 <- factor(data$time)
      levels(time2)[pars$T] <- levels(time2)[1]
      tw.droptime <- lm(y ~ x + factor(case) + time2, data=data)

      ##Twoway FE, without dropping any factor categories
      twoway.fe <- lm(y ~ x + factor(case) + factor(time), data=data)

      #Check if Twoway FE is identified
      tw.singular <- any(summary(twoway.fe)[[5]])

      #5) SD and Variance weighted average of between effects and within effects
      ## Calculating from theoretical DGP parameters
      dgp.avg <- (pars$betaBmean*pars$betaWsd+pars$betaWmean*pars$betaBsd)/
            (pars$betaWsd+pars$betaBsd)
      dgp.avg.var <- (pars$betaBmean*pars$betaWsd^2+pars$betaWmean*pars$betaBsd^2)/
            (pars$betaWsd^2+pars$betaBsd^2)

      ## Calculating from realized cross-sectional coefficients
      real.case <- as.matrix(summarize(group_by(data,case), betaWmean=mean(betaW)))
      real.time <- as.matrix(summarize(group_by(data,time), betaWmean=mean(betaB)))
      coefBmean <- mean(real.time[,2])
      coefBsd <- sd(real.time[,2])
      coefWmean <- mean(real.case[,2])
      coefWsd <- sd(real.case[,2])
      real.avg <- (coefBmean*coefWsd+coefWmean*coefBsd)/(coefWsd+coefBsd)
      real.avg.var <- (coefBmean*coefWsd^2+coefWmean*coefBsd^2)/(coefWsd^2+coefBsd^2)

      #6) Constructing the Imai and Kim result (proposition 4)
      ##Extracting the relevant coefficients
      beta.fe <- coef(case.fe.dummy)[2]
      beta.fe.time <- coef(time.fe.dummy)[2]
      beta.pool <- coef(pooled.ols)[2]

      ##Calculating the weights
      #wgt.pool
      wgt.pool <- sum((data$x - mean(data$x))^2)/(nrow(data)-1)

      #wgt.fe
      data.case <- split(data, data$case)
      S.i.squared <- sapply(data.case, FUN=function(d){
            sum((d$x - mean(d$x))^2)/(nrow(data)-1)
      }, simplify=TRUE)
      wgt.fe <- mean(S.i.squared)

      #wgt.fe.time
      data.time <- split(data, data$time)
      S.t.squared <- sapply(data.time, FUN=function(d){
            sum((d$x - mean(d$x))^2)/(nrow(data)-1)
      }, simplify=TRUE)
      wgt.fe.time <- mean(S.t.squared)

      ## Calculate the weighted average
      imai.kim <-  (wgt.fe*beta.fe + wgt.fe.time*beta.fe.time - wgt.pool*beta.pool)/
            (wgt.fe + wgt.fe.time - wgt.pool)

      ### Compile saved results in a list
      return(list(pooled.ols=summary(pooled.ols)$coefficients[2,1],
            pooled.ols.se=summary(pooled.ols)$coefficients[2,2],

            case.fe.dummy=summary(case.fe.dummy)$coefficients[2,1],
            case.fe.dummy.se=summary(case.fe.dummy)$coefficients[2,2],
            time.fe.dummy=summary(time.fe.dummy)$coefficients[2,1],
            time.fe.dummy.se=summary(time.fe.dummy)$coefficients[2,2],

            random.ui=summary(random.ui)$coefficients[2,1],
            random.ui.se=summary(random.ui)$coefficients[2,2],
            #random.vt=summary(random.vt)$coefficients[2,1],
            #random.vt.se=summary(random.vt)$coefficients[2,2],

            twoway.fe=summary(twoway.fe)$coefficients[2,1],
            twoway.fe.se=summary(twoway.fe)$coefficients[2,2],
            tw.coef.dc=summary(tw.dropcase)$coefficients[2,1],
            tw.se.dc=summary(tw.dropcase)$coefficients[2,2],
            tw.coef.dt=summary(tw.droptime)$coefficients[2,1],
            tw.se.dt=summary(tw.droptime)$coefficients[2,2],
            tw.singular=tw.singular,

            dgp.avg=dgp.avg, real.avg=real.avg,
            dgp.avg.var=dgp.avg.var, real.avg.var=real.avg.var,

            imai.kim=imai.kim))
}

#' @import dplyr
#' @export
twtransform <- function(d){

      pars <- d$pars
      if(pars$unbalance){
            data <- d$data[,which(names(d$data) %in% c("case", "time", "x", "y", "I"))]
            data <- mutate(data, xI=x*I, yI=y*I)
            data.case <- group_by(data, case)
            data.case <- summarize(data.case, sum.xI.case=sum(xI),
                 sum.yI.case=sum(yI), sum.I.case=sum(I))
            data.time <- group_by(data, time)
            data.time <- summarize(data.time, sum.xI.time=sum(xI),
                 sum.yI.time=sum(yI), sum.I.time=sum(I))
            data2 <- merge(data, data.case, by="case")
            data2 <- merge(data2, data.time, by="time")
            data2 <- mutate(data2,
                            x.bar.i = sum.xI.case/sum.I.case,
                            x.bar.t = sum.xI.time/sum.I.time,
                            x.bar = sum(xI)/(sum.I.case*sum.I.time),
                            y.bar.i = sum.yI.case/sum.I.case,
                            y.bar.t = sum.yI.time/sum.I.time,
                            y.bar = sum(yI)/(sum.I.case*sum.I.time))
      } else {
            data <- d$data[,which(names(d$data) %in% c("case", "time", "x", "y"))]
            data.case <- group_by(data, case)
            data.case <- summarize(data.case, x.bar.i=mean(x), y.bar.i=mean(y))
            data.time <- group_by(data, time)
            data.time <- summarize(data.time, x.bar.t=mean(x), y.bar.t=mean(y))
            data2 <- merge(data, data.case, by="case")
            data2 <- merge(data2, data.time, by="time")
            data2 <- mutate(data2, x.bar=mean(x), y.bar=mean(y))
      }

      data2 <- mutate(data2,
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

#' @import parallel
#' @export
twsim <- function(iter=1000, cores=NULL, ...){

      z <- list(...)
      if(is.null(cores)) cores <- detectCores()
      p <- mclapply(1:iter, FUN=function(x){
            d <- do.call(twsimdata, args=z)
            m <- twmodel(d)
            m$iter <- x
            return(m)
      }, mc.cores=cores)
}

# From the dimension confusion bias simulation
#' @export
dc.sim <- function(iter, zbeta, betaWmean, argnum){

  dta <- twsimdata(N=30, T=30, ui.mean=0, ui.sd=1, vt.mean=0, vt.sd=1,
                   ep.ui.sd=.1, ep.vt.sd=.1, betaBmean=3, betaBsd=0, betaWmean=betaWmean, betaWsd=0,
                   binary=FALSE, unbalance=FALSE, gamma.ux=zbeta, gamma.vx=zbeta,
                   gamma.uy=.5, gamma.vy=.5)

  data <- dta$data

  xi <- data.frame(summarize(group_by(data,case), dgp=mean(x)))[,2]
  xt <- data.frame(summarize(group_by(data,time), dgp=mean(x)))[,2]
  yi <- data.frame(summarize(group_by(data,case), dgp=mean(y)))[,2]
  yt <- data.frame(summarize(group_by(data,time), dgp=mean(y)))[,2]
  zi <- data.frame(summarize(group_by(data,case), dgp=mean(zi)))[,2]
  zt <- data.frame(summarize(group_by(data,time), dgp=mean(zt)))[,2]

  corr.xz.i <- cor(cbind(xi,yi,zi))[1,3]
  corr.xz.t <- cor(cbind(xt,yt,zt))[1,3]

  reg1.1 <- lm(y ~ x + factor(case), data=data)
  reg2.1 <- lm(y ~ x + factor(time), data=data)
  reg3.1 <- lm(y ~ x + factor(case) + factor(time), data=data)

  beta.caseFE <- coef(reg1.1)[2]
  beta.timeFE <- coef(reg2.1)[2]
  beta.twowayFE <- coef(reg3.1)[2]

  reg1.2 <- lm(y ~ x + factor(case) + zi + zt, data=data)
  reg2.2 <- lm(y ~ x + factor(time) + zi + zt, data=data)
  reg3.2 <- lm(y ~ x + factor(case) + factor(time) + zi + zt, data=data)
  beta.caseFE2 <- coef(reg1.2)[2]
  beta.timeFE2 <- coef(reg2.2)[2]
  beta.twowayFE2 <- coef(reg3.2)[2]

  argnum <- argnum

  res <- c(iter, argnum, zbeta, betaWmean, corr.xz.i, corr.xz.t,
           beta.caseFE, beta.timeFE, beta.twowayFE,
           beta.caseFE2, beta.timeFE2, beta.twowayFE2)
  names(res) <- c("iter", "argnum", "zbeta", "betaWmean", "corr.xz.i", "corr.xz.t",
                  "beta.caseFE", "beta.timeFE", "beta.twowayFE",
                  "beta.caseFE2", "beta.timeFE2", "beta.twowayFE2")
  return(res)
}

#' @import parallel
#' @export
twsim <- function(iterations=1000, cores=NULL, ...){

  if(is.null(cores)) cores <- detectCores()
  iter <- data.frame(one=1,iter=1:iterations)
  zbeta <- data.frame(one=1,zbeta=seq(0,.5,by=.05))
  betaWmean <- data.frame(one=1, betaWmean=seq(-3,3,by=.25)[-length(seq(-3,3,by=.25))])
  args <- merge(iter, zbeta, by="one")
  args <- merge(args, betaWmean, by="one")
  args <- data.frame(args[,-1], argnum=1:nrow(args))
  p <- mcmapply(1:nrow(args), FUN=function(x){
    d <- do.call(dc.sim, args=as.list(args[x,]))
    return(d)
  }, mc.cores=cores, SIMPLIFY=TRUE)
}


# From Empirical Examples -------------------------------------------------

# used to iterate over V-DEM posterior samples
over_posterior <- function(x=NULL,y=NULL,modelfunc=NULL,merged_data=NULL,select_vars=NULL,...) {
  to_analyze <- merged_data[,c(x,select_vars)]
  names(to_analyze) <- c('v2x_polyarchy',select_vars)
  # Need to avoid an evaluation error by resetting modelfunc to current namespace
  modelfunc <- force(modelfunc)
  model1 <- modelfunc(formula = y,data = to_analyze,...)
  # Use the sandwich estimator to adjust variances
  # Need to drop coefs that come out as NA
  # This happens with the country of the Republic of Vietnam in the within-between model (model 5)
  # The sandwich estimator will automatically drop NA coefficients from the VCOV matrix, causing an error
  # With MASS
  coefs <- coef(model1)[!is.na(coef(model1))]
  sds <- sandwich::vcovHC(model1,type='HC0')

  samples <- MASS::mvrnorm(mu=coefs,Sigma=sds)
  return(samples)
}

#' VDEM download function
#' @export
download_vdem <- function() {
  if(!dir.exists('data/')) {
    dir.create('data/')
  }
  print("Note, this file is 350MB and will take some time to download.")
  download.file('https://virginia.box.com/shared/static/20lo0j7pmjpi0i6iyezbo4ewbmdb2zbj.sqlite',
                destfile='data/vdem_data.sqlite',mode='wb')
}


#' @export
run_vdem <- function(varnames=NULL,full_formula=NULL,modelfunc=lm,select_vars=NULL,num_cores=1,
                     num_iters=900,dbcon=NULL,...) {

  #If SQLite database does not exist, then make user run download function
  if(!file.exists('data/vdem_data.sqlite')) {
      stop('You must first run the download_vdem() function before running VDEM models to download the vdem data.')
  }


  # Load analysis variables from SQLITE database
  dbcon <- RSQLite::dbConnect(RSQLite::SQLite(), dbcon)
  select_vars <- c(select_vars,'country_text_id','year')
  vdem_data <- RSQLite::dbGetQuery(dbcon,paste0('SELECT ',paste0(select_vars,collapse=','),' FROM vdem_data')) %>%
    as_tibble %>% mutate(year_factor=factor(year))

# if rest than the full subset of the posterior samples are used, use a specific number of samples
    if(is.null(num_iters) || num_iters==900) {
    num_iters <- 1:length(varnames)
  } else {
    num_iters <- sample(1:length(varnames),num_iters)
  }
  var_collapse <- c(varnames[num_iters],'country_text_id','year')
  var_collapse <- paste0(var_collapse,collapse=",")
  vdem_pos <- RSQLite::dbGetQuery(dbcon,paste0("SELECT ",var_collapse," FROM vdem_pos"))
  # merge together for analysis
  merged_data <- dplyr::left_join(vdem_data,vdem_pos,by=c('country_text_id','year'))
  model1 <- parallel::mclapply(varnames[num_iters],over_posterior,y=full_formula,select_vars=select_vars,modelfunc=modelfunc,merged_data=merged_data,...,mc.cores=num_cores)
  beta_names <- names(model1[[1]])
  model1 <- matrix(unlist(model1),nrow=length(num_iters),dimnames=list(num_iters,beta_names),byrow = TRUE)
  results <- tibble::data_frame(betas=colnames(model1),coef=apply(model1,2,mean),
             sd=apply(model1,2,sd),upper=coef + 1.96*sd,lower=coef - 1.96*sd)


  RSQLite::dbDisconnect(dbcon)

  return(results)
}

#' @import ggplot2
#' @import plotly
#' @export
panel_balance <- function(dbcon=NULL,select_vars=NULL,use_plotly=TRUE) {

   dbcon <- RSQLite::dbConnect(RSQLite::SQLite(), dbcon)

  if(!(("year_factor" %in% select_vars) && ("country_name"%in% select_vars))) {
    select_vars <- c(select_vars,'year_factor','country_name')
  }
  select_vars <- select_vars[!duplicated(select_vars)]

  merged_data <- RSQLite::dbGetQuery(dbcon,paste0('SELECT ',paste0(select_vars,collapse=','),' FROM vdem_data')) %>% as.tbl
  if(use_plotly==TRUE) {
    countries <- plot_ly(merged_data,x=~country_name) %>% add_histogram()
    years <-  plot_ly(merged_data,x=~year_factor) %>% add_histogram()
  } else {
    countries <- ggplot(data=merged_data,aes(reorder(country_name,country_name,function(x)-length(x)))) + geom_bar(alpha=0.5) + ylab("") + xlab("")
    years <- ggplot(data=merged_data,aes(reorder(year_factor,year_factor,function(x)-length(x)))) + geom_bar(alpha=0.5) + ylab("") + xlab("")
  }


 RSQLite::dbDisconnect(dbcon)

 return(list(countries=countries,years=years))
}


# Helper Functions --------------------------------------------------------


#' @export
to_round <- function(colvar,sigdif) {
  outputvar <- format(round(colvar,digits=sigdif),nsmall=sigdif)
  return(outputvar)
}

#' @export
prepare_for_display <- function(x) {
  x <- filter(x, betas %in% control_vars) %>% mutate_if(is.numeric,to_round,sigdif=3) %>% mutate(TCI=paste0("(",lower,",",upper,")"))
  x <-  x %>%  select(betas,coef,TCI) %>% gather(value_type,to_display,-betas) %>% unite(row_value,betas,value_type)
  return(x)
}

