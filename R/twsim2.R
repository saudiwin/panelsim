###################################################
# Simulation files for "Interpretation and Identification of Within-unit and Cross-sectional Variation in Panel Data Models"
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
# tw_data() generates the simulated data. Arguments:
# N: Number of cases (default = 30)
# T: Number of time points (default = 30)
# case.int.mean: Mean of case-specific intercepts (default = 0)
# case.int.sd: Standard deviation of case-specific intercepts (default = 1)
# cross.int.mean: Mean of time-specific intercepts (default = 0)
# cross.int.sd: Standard deviation of time-specific intercepts (default = 1)
# cross.eff.mean: Mean of the (cross-sectional) slopes of the best fit lines within time points (default = 3)
# cross.eff.sd: SD of the (cross-sectional) slopes of the best fit lines within time points (default = .5)
# case.eff.mean: Mean of the (over-time) slopes of the best fit lines within cases (default = -3)
# case.eff.sd: SD of the (over-time) slopes of the best fit lines within cases (default = .5)
# noise.sd: SD of the error applied after generation of y (default = 1)
# omm.x.case: coefficient for how an omitted variable z_i depends on x (default = 0)
# omm.x.cross: coefficient for how an omitted variable z_t depends on x (default = 0)
# omm.y.case: coefficient for how an omitted variable z_i depends on y (default = 0)
# omm.y.cross: coefficient for how an omitted variable z_t depends on y (default = 0)
# binary: whether to make the predictor X binary by cutting at the median (default = FALSE)
# unbalance: whether to generate unbalanced panels with time frames 1-K,
#     where k is drawn from ceiling(runif(min=(T/3), max=T)) (default = FALSE)
# time.ac: correlation between errors in adjacent time points for the same case (default = 0)
# spatial.ac: correlation between errors in adjacent cases at the same time (default = 0)
#
# tw_reg() runs regressions within each case and time to compare
# DGP coefficients to realized coefficients.
# Arguments: $data element of output of tw_data()
#
# tw_model() runs pooled OLS, one-way FEs, RE, and two-way FE models
# Arguments: output of tw_data()
#
# twtransform() creates a second version of the simulated data that
# contains case-level, time-level, and overall means of x and y and
# the twoway tranformations.
# Arguments: output of tw_data()
#
# twtrans2() multiplies the various transformed variables from twtransform()
# for the purpose of investigating the formal one-way and two-way FE estimators.
# Arguments: output of twtransform()
#
# tw_sim() runs the simulation. Arguments:
# iter: number of iterations
# cores: the number of processing cores to devote to the simulation (defaults to detectCores())
# parallel: whether the simulation should run with parallel processing (default = FALSE)
# ... : arguments passed to tw_data()
###################################################


#' Function to generate error matrix
#'
#' This function generates a distribution of errors given a set of
#' over-time or between-cases coefficients. Errors can be affected by
#' either spatial or time autocorrelation.
#'
#' @import dplyr tidyr
#' @importFrom MASS mvrnorm
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

#' Generate One-way and Two-way Fixed Effects Panel Data
#'
#' This function will produce panel data where variation can exist in the cross-section, over time
#' or in both dimensions simultaneously. Furthermore, effect heterogeneity by case or cross section
#' is also allowed, along with interactive effects such as differences-in-differences.
#'
#' The \code{tw_data} function is the workhorse of the \code{twowaysim} package. It accepts as
#' input the dimensions of the panel/TSCS data to be generated, and also parameters that
#' determine the extent of variance and heterogeneity in either the cross-sectional or
#' over-time effects in the data or the interaction thereof. The parameter \code{N} determines how many observations
#' exist for each case or unit in the panel, while \code{T} determines how many time points exist
#' per case or unit. To create a model with a homogenous (static) within-unit over-time (case) effect,
#' simply set \code{case.eff.mean} to a non-zero number and set \code{case.eff.sd} to zero. Similarly,
#' setting \code{cross.eff.mean} to a non-zero number and \code{cross.eff.sd} to zero will produce a
#' panel dataset with a cross-sectional effect of X on Y where the effect of X does not vary across
#' countries (no effect heterogeneity). Increasing \code{cross.eff.sd} and \code{case.eff.sd} will result in more
#' effect heterogeneity across countries and time points. If both \code{case.eff.mean} and \code{cross.eff.mean} are
#' non-zero, then Y will have both dimensions of variance. A 1-way fixed effects model with intercepts on
#' cases will return the \code{case.eff.mean} coefficient and a model with intercepts on time points will return
#' the \code{cross.eff.mean} estimate, whereas a 2-way model (intercepts on cases and time points) will return
#' a difficult-to-characterize weighted average.
#'
#' To estimate a difference-in-differences (time interaction) effect, set \code{did.eff.mean} to a non-zero number and set \code{did.eff.sd} to zero if the DiD effects are supposed to be homogenous. Note that to estimate a standard "canonical" DiD setup, the time points \code{T} should be no more than 2. For more information about generating DiD specifications with `tw_data`, see the blog post \url{https://www.robertkubinec.com/post/did_dnd/index.html}.
#'
#' We refer you to Kropko and Kubinec (2020) for more information
#' on the difference between these models:
#' \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0231349}.
#'
#' The parameters \code{case.int} and \code{cross.int} represent the values of the intercepts for the
#' cases or time points. Changing these parameters will increase or decrease the amount of unexplained
#' variance (random noise) in the dataset.
#'
#' The additional parameters in the function allow the user to create unbalanced panels (varying numbers
#' of observations per case or time point if \code{unbalance=TRUE}), auto-correlation in the effects and
#' omitted variables. Autocorrelation can exist either in the over time dimension or the cross-sectional
#' dimension. To increase time autocorrelation, set \code{time.ac} to a value between 0 and 1 where values
#' closer to one signal higher autocorrelation. To increase spatial (cross-sectional) autocorrelation,
#' set \code{spatial.ac} to a value between 0 and 1.
#'
#' Finally, to include omitted variables, set one of the \code{omm} parameters to a non-zero value.
#' Omitted variables that vary within cases (over time) can be included by setting an \code{omm} parameter
#' subscripted with \code{case} to a non-zero value, and the same is possible for variables that vary
#' in the cross-section \code{cross}. The analyst can also decide whether the omitted variable is correlated
#' with the independent variable of interest \code{x} or the dependent variable \code{y} by choosing the
#' subscript of \code{omm}.
#'
#' @return The function returns a named list where \code{object$data} is a \code{data.frame} and
#' \code{object$pars} are the original parametes used to generate the data. The value of generated coefficients is returned in a list as \code{object$fixed_params}. For repeated sampling, the object can be given to the `prior_true_vals` argument to allow for fixed population parameters for intercepts and vectors of effects.
#'
#' @examples
#'
#' # case (over-time) effect with no effect heterogeneity
#'
#' case1 <- tw_data(case.eff.mean=-1,case.eff.sd=0)
#'
#' # case (over-time) effect with substantial effect heterogeneity across countries
#'
#' case2 <- tw_data(case.eff.mean=-1,case.eff.sd=1)
#'
#' # cross-section effect with no effect heterogeneity
#'
#' cross1 <- tw_data(cross.eff.mean=-1,cross.eff.sd=0)
#'
#' # cross-section effect with substantial effect heterogeneity across countries
#'
#' cross2 <- tw_data(cross.eff.mean=-1,cross.eff.sd=1)
#'
#' # panel data with a cross-sectional effect of 3 and a case (over-time) effect of -1
#'
#' both_case_cross <- tw_data(cross.eff.mean=3,
#'                              case.eff.mean=-1)
#'
#' @seealso \code{\link{tw_model}} for running linear models on the data and and
#' \code{\link{tw_sim}} function for running Monte Carlo simulations on panel data.
#'
#' @param N The number of observations for each case/unit.
#' @param T The number of time points per observation.
#' @param case.int.mean The mean of the case/unit intercepts/fixed effects
#' @param case.int.sd The SD of the case/unit intercepts/fixed effects
#' @param cross.int.mean The mean of the cross-sectional intercepts/fixed effects
#' @param cross.int.sd The SD of the cross-sectional intercepts/fixed effects
#' @param cross.eff.mean The mean of the cross-sectional effect of X on Y
#' @param cross.eff.sd The SD of the cross-sectional effect of X on Y
#' @param case.eff.mean The mean of the case (over-time) effect of X on Y
#' @param case.eff.sd The SD of the case (over-time) effect of X on Y
#' @param did.eff.mean The mean of the difference-in-difference effect of X on Y
#' @param did.eff.sd The SD of the difference-in-difference effect of X on Y
#' @param wid.eff.mean The mean of the difference-in-cases effect of X on Y
#' @param wid.eff.sd The SD of the difference-in-cases effect of X on Y
#' @param noise.sd The residual variance of the data
#' @param omm.x.case The value of an omitted variable correlated with X that varies across cases/units
#' @param omm.x.cross The value of an omitted variable correlated with X that varies cross-sectionally
#' @param omm.y.case The value of an omitted variable correlated with Y that varies across cases/units
#' @param omm.y.cross The value of an omitted variable correlated with Y that varies cross-sectionally
#' @param treat_effect A vector of length 1 or N*T that is equal to 1 for assignment to treatment and 0 for assignment to control. If this argument is not NULL, `tw_data` will generate data for y that takes x as fixed to these values (i.e., treatment is being assigned/manipulated).
#' @param binary_outcome Whether the Y (outcome) variable should be converted to 0/1. Note that this can lead to some measurement bias as Y is simulated as continuous.
#' @param binary_x Whether X should be converted to 0/1. Note that doing so may lead to some bias in estimation as X is simulated as continuous.
#' @param unbalance Whether to simulate varying numbers of observations by cases or time points.
#' @param time.ac A value between 0 and 1 giving the over-time autocorrelation in effect of X on Y
#' @param spatial.ac A value between 0 and 1 giving the cross-sectional (spatial)
#' autocorrelation in the effect of X on Y
#' @param prior_true_vals A fitted `tw_data` object with generated coefficients that can be used to keep vectors of intercepts/effects fixed over repeated sampling
#'
#'
#' @export
tw_data <- function(N = 30, T = 30, case.int.mean = .5, case.int.sd = .5,
                      cross.int.mean = -1, cross.int.sd = .5, cross.eff.mean = 0,
                      did.eff.mean = 1, did.eff.sd = 0.25,
                      wid.eff.mean=0,wid.eff.sd=0,
                      cross.eff.sd = .5, case.eff.mean = 0.5, case.eff.sd = .5, noise.sd = 1,
                      omm.x.case = 0, omm.x.cross = 0, omm.y.case = 0, omm.y.cross = 0,
                      treat_effect = NULL,this_case=NULL,this_time=NULL,
                    binary_outcome=FALSE,binary_x=FALSE,unbalance = FALSE,
                    time.ac = 0, spatial.ac = 0,
                    prior_true_vals=NULL) {

      # if treatment effects are used, need to exclude negative numbers
  if(case.eff.sd>0 && is.null(prior_true_vals)) {

    gamma <- rnorm(N, mean=case.eff.mean, sd=case.eff.sd)

  } else if(!is.null(prior_true_vals)) {

    gamma <- prior_true_vals$fixed_params$gamma

  }  else {
    gamma <- rep(case.eff.mean,N)
  }

  if(cross.eff.sd>0 && is.null(prior_true_vals)) {

    beta <- rnorm(T, mean=cross.eff.mean, sd=cross.eff.sd)

  } else if(!is.null(prior_true_vals)) {

    beta <- prior_true_vals$fixed_params$beta

    } else {
    beta <- rep(cross.eff.mean,T)
  }

  if(did.eff.sd>0 && is.null(prior_true_vals)) {

    did <- rnorm(T, mean=did.eff.mean, sd=did.eff.sd)

  } else if(!is.null(prior_true_vals)) {

    did <- prior_true_vals$fixed_params$did

  } else {

    did <- rep(did.eff.mean,T)

  }

  if(wid.eff.sd>0 && is.null(prior_true_vals)) {

    wid <- rnorm(N, mean=wid.eff.mean, sd=wid.eff.sd)

  } else if(!is.null(prior_true_vals)) {

    wid <- prior_true_vals$fixed_params$wid

  } else {

    wid <- rep(wid.eff.mean,N)

  }

  # need global intercepts

  if(is.null(prior_true_vals)) {

    i.intercept <- rnorm(1)
    t.intercept <- rnorm(1)

  } else {

    i.intercept=prior_true_vals$fixed_params$i.intercept
    t.intercept=prior_true_vals$fixed_params$t.intercept

  }


      if(binary_outcome) {
        if(case.eff.mean>.9999 || cross.eff.mean>.9999 || omm.x.cross>.9999 ||
           omm.y.cross>.9999 || omm.x.case>.9999 || omm.y.case>.9999) {
          stop("You specified a value for a treatment/omitted variable that is greater than .9999 (i.e. not a probability).\n
               Change treat_effect to false to calculate on all real numbers instead of probabilities.")
        }

        case.eff.mean <- suppressWarnings(qlogis(case.eff.mean+0.5))
        cross.eff.mean <- suppressWarnings(qlogis(cross.eff.mean+0.5))
        omm.x.cross <- suppressWarnings(qlogis(omm.x.cross+0.5))
        omm.x.case <- suppressWarnings(qlogis(omm.x.case+0.5))
        omm.y.case <- suppressWarnings(qlogis(omm.y.case+0.5))
        omm.y.cross <- suppressWarnings(qlogis(omm.y.cross+0.5))

        # case.eff.mean <- qlogis(case.eff.mean)
        # cross.eff.mean <- qlogis(cross.eff.mean)
        # omm.x.cross <- qlogis(omm.x.cross)
        # omm.x.case <- qlogis(omm.x.case)
        # omm.y.case <- qlogis(omm.y.case)
        # omm.y.cross <- qlogis(omm.y.cross)

      }

      # we need different data if we are going to simulate gsynt

        if(case.int.sd>0 & is.null(prior_true_vals)) {

          alpha.i <- rnorm(N, mean=case.int.mean, sd=case.int.sd)

        } else if(!is.null(prior_true_vals)) {

          alpha.i <- prior_true_vals$fixed_params$alpha.i

        } else {

          alpha.i <- rep(case.int.mean, N)

        }

          bet.data <- data.frame(case=1:N,
                                 alpha.i=alpha.i,
                                 gamma=gamma,
                                 zi=rnorm(N),
                                 wid=wid,
                                 unbal=ifelse(unbalance, ceiling(runif(N, min=(T/3), max=T)), T),
                                 one=1) %>%
            mutate(alpha.i=ifelse(case==1,0,alpha.i),
                   wid=ifelse(case==1,0,wid))

        if(cross.int.sd>0 & is.null(prior_true_vals)) {

          alpha.t <-  rnorm(T, mean=cross.int.mean, sd=cross.int.sd)

        } else if(!is.null(prior_true_vals)) {

          alpha.t <- prior_true_vals$fixed_params$alpha.t

        } else {

          alpha.t <- rep(cross.int.mean, N)

        }

          with.data <- data.frame(time=1:T,
                                  alpha.t=alpha.t,
                                  beta=beta,
                                  did=did,
                                  zt=rnorm(T),
                                  one=1) %>%
            mutate(alpha.t=ifelse(time==1,0,alpha.t),
                   did=ifelse(time==1,0,did))




      data <- full_join(bet.data, with.data, by="one")


      if(time.ac==0 && spatial.ac==0) {

        # if no autocorrelation, much faster to just use rnorm

        if(is.null(treat_effect)) {

          # generate both X and Y

          data <- data %>%
            mutate(noise=rnorm(n=N*T,sd=noise.sd),
                   denom=(beta + did - gamma - wid),
                   x=((alpha.i + i.intercept - alpha.t - t.intercept)/denom) + omm.x.case*zi + omm.x.cross*zt,
                   y= ((beta * alpha.i + did*alpha.i + beta * i.intercept + did * i.intercept - alpha.t * gamma - wid * alpha.t - gamma * t.intercept - wid * t.intercept)/denom) + omm.y.case*zi + omm.y.cross*zt + noise,
                   tokeep = time <= unbal) %>%
            dplyr::select(case, time, y, x, gamma, beta, alpha.i, alpha.t, noise,
                          zi, zt, wid, did, tokeep) %>%
            arrange(case, time)

        } else {

          # X is assigned, generate Y

          # make the treat effect longer if it needs to be

          if(length(treat_effect)==1) treat_effect <- rep(treat_effect, N*`T`)

          data <- data %>%
            mutate(noise=rnorm(n=N*T,sd=noise.sd),
                   x=treat_effect,
                   y= ifelse(treat_effect,
                             t.intercept + alpha.t + beta + did + omm.y.case*zi + omm.y.cross*zt + noise,
                             t.intercept + alpha.t + noise),
                   tokeep = time <= unbal) %>%
            dplyr::select(case, time, y, x, gamma, beta, alpha.i, alpha.t, noise,
                          zi, zt, wid, did, tokeep) %>%
            arrange(case, time)

        }

        if(binary_outcome) data$y <- as.numeric(data$y > median(data$y))
        if(binary_x) data$x <- as.numeric(data$x > median(data$x))

      } else {
        eps <- gen_errormat(data, time.ac, spatial.ac)
        data <- data %>%
          mutate(noise=mvrnorm(n=1, mu=rep(0, N*T), Sigma=eps)*noise.sd,
                 x=(alpha.i-alpha.t)/(beta - gamma) + omm.x.case*zi + omm.x.cross*zt,
                 y=(beta*alpha.i-gamma*alpha.t)/(beta - gamma) + omm.y.case*zi + omm.y.cross*zt + noise,
                 tokeep = time <= unbal) %>%
          dplyr::select(case, time, y, x, gamma, beta, alpha.i, alpha.t, noise,
                        zi, zt, tokeep) %>%
          arrange(case, time)
      }
      #if(binary) data <- mutate(data, x=(x>=median(x)))

      pars <- data.frame(N = N, T = T, case.int.mean = case.int.mean, case.int.sd = case.int.sd,
                      cross.int.mean = cross.int.mean, cross.int.sd = cross.int.sd, cross.eff.mean = cross.eff.mean,
                      cross.eff.sd = cross.eff.sd, case.eff.mean = case.eff.mean, case.eff.sd = case.eff.sd, noise.sd = noise.sd,
                      omm.x.case = omm.x.case, omm.x.cross = omm.x.cross, omm.y.case = omm.y.case, omm.y.cross = omm.y.cross,
                      binary_x=binary_x,
                      binary_outcome=binary_outcome,
                      unbalance = unbalance, time.ac = time.ac, spatial.ac = spatial.ac)

      fixed_params <- list(gamma=gamma,
                           wid=wid,
                           beta=beta,
                           did=did,
                           alpha.i=bet.data$alpha.i,
                           alpha.t=with.data$alpha.t,
                           i.intercept=i.intercept,
                           t.intercept=t.intercept)

      toreturn <- list(data, pars, fixed_params)
      names(toreturn) <- c("data", "pars", "fixed_params")
      return(toreturn)
}

#' @export
tw_reg <- function(d){
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

#' Fit Panel Data Given Simulated Data
#'
#' This function will take the result of data generated by the
#' \code{\link{tw_data}} function and fit panel data models
#' using the \code{\link[plm]{plm}} model-fitting function.
#' Available models include standard fixed effects models (1-way
#' and 2-way), random effects models and pooled OLS.
#'
#' @param gendata A data frame of simulated panel data procued by
#' \code{\link{tw_data}}
#' @param models Whether to include random effects estimates in the
#' returned data frame (defaults to true).
#'
#' @return The function returns a \code{data.frame} with the coefficients
#' of all of the requested models (includes 1-way and 2-way FE models by default).
#'
#' @examples
#'
#' \dontrun{
#' # generate panel data
#' gendata <- tw_data()
#' # estimate panel models
#' estimates <- tw_model(gendata)
#' }
#'
#'
#' @import gsynth
#' @import wfe
#' @export
tw_model <- function(gendata, models=c("randomfx","wfe"),binary_outcome=FALSE){

      data <- gendata$data
      pars <- gendata$pars

      pdata <- data[,apply(data, 2, sd)!=0]
      pdata <- select(pdata,time_id="time",everything())
      pdata.time <- pdata.frame(as.data.frame(pdata), index=c("case","time_id"), row.names=FALSE,drop.index = T)
      pdata.xs <- pdata.frame(as.data.frame(pdata), index=c("time_id","case"), row.names=FALSE,drop.index = T)

      # need to add binary x for wfe/gsynth
      if(!(all(pdata$x %in% c(0,1)))) {
        pdata$bin_x <- as.numeric(pdata$x>median(pdata$x))
      } else {
        pdata$bin_x <- pdata$x
      }

      #1) Pooled OLS
      pooled.ols <- lm(y ~ x, data=data)

      #2) One way FEs

      case.fe.plm <- plm(y ~ x, model="within", data=pdata.time)
      time.fe.plm <- plm(y ~ x, model="within", data=pdata.xs)



      #3) Random effects
      if("randomfx" %in% models) {
        random.xs <- plm(y ~ x, model="random", data=pdata.xs)
        random.time <- plm(y ~ x, model="random", data=pdata.time)

        random_time_coef <- summary(random.time)$coefficients[2,1]
        random_time_se <- summary(random.time)$coefficients[2,2]
        random_time_pval <- summary(random.time)$coefficients[2,4]

        random_xs_coef <- summary(random.xs)$coefficients[2,1]
        random_xs_se <- summary(random.xs)$coefficients[2,2]
        random_xs_pval <- summary(random.xs)$coefficients[2,4]
      }

      #4) Twoway FE
      twoway.fe <- lm(y ~ x + factor(case) + factor(time), data=data)
      identified <- all(!is.na(coef(twoway.fe)))
      twcoef <- summary(twoway.fe)$coefficients[2,1]
      twse <- summary(twoway.fe)$coefficients[2,2]
      tw_pval <- summary(twoway.fe)$coefficients[2,4]

      # 5) Generalized Synthetic Control

      if("gsynth" %in% models) {

        gsynth_out <- try(suppressMessages(gsynth(formula = y~bin_x,index=c("case","time_id"),
                                                  EM=T,
                                                  inference="parametric",
                                                  parallel=F,
                                                  data=pdata,se=T)))

        if(!('try-error' %in% class(gsynth_out))) {
          coef_gsynth <- gsynth_out$est.avg[,"ATT.avg"]
          p_gsynth <- gsynth_out$est.avg[,"p.value"]
          se_gsynth <- gsynth_out$est.avg[,"S.E."]
          r_gsynth <- gsynth_out$r.cv
        } else {
          print(paste("R package gsynth failed to estimate because of ",gsynth_out))
          coef_gsynth <- NA
          p_gsynth <- NA
          se_gsynth <- NA
          r_gsynth <- NA
        }


      } else {
        coef_gsynth <- NA
        p_gsynth <- NA
        se_gsynth <- NA
        r_gsynth <- NA
      }

      # 6) Weighted Fixed Effects

      if("wfe" %in% models) {
        wfe_out <- try(suppressMessages(wfe(y~bin_x,data=pdata,treat="bin_x",unit.index="case",time.index="time_id")))

        if(!('try-error' %in% class(wfe_out))) {
          sum_wfe <- summary(wfe_out)
          coef_wfe <- sum_wfe$coefficients[,"Estimate"]
          p_wfe <- sum_wfe$coefficients[,"p.value"]
          se_wfe <- sum_wfe$coefficients[,"Std.Error"]
        } else {
          print(paste("R package wfe failed to estimate because of ",wfe_out))
          coef_wfe <- NA
          p_wfe <- NA
          se_wfe <- NA
        }
        coef_wfe <- NA
        p_wfe <- NA
        se_wfe <- NA

      } else {
        coef_wfe <- NA
        p_wfe <- NA
        se_wfe <- NA
      }

      # 7) Logit GLM (if binary_outcome)

      if(binary_outcome) {
        bin_out_case <- try(suppressMessages(glm(y~x + factor(case),family="binomial",data=pdata)))
        bin_out_time <- try(suppressMessages(glm(y~x + factor(time_id),family="binomial",data=pdata)))
        bin_out_twoway <- try(suppressMessages(glm(y~x + factor(time_id) + factor(case),family="binomial",data=pdata)))

        if(!('try-error') %in% bin_out_case) {
          bin_out_case_coef <- summary(bin_out_case)$coefficients[2,1]
          bin_out_case_se <- summary(bin_out_case)$coefficients[2,2]
          bin_out_case_pval <- summary(bin_out_case)$coefficients[2,4]
        } else {
          bin_out_case_coef <- NA
          bin_out_case_se <- NA
          bin_out_case_pval <- NA
        }

        if(!('try-error') %in% bin_out_time) {
          bin_out_time_coef <- summary(bin_out_time)$coefficients[2,1]
          bin_out_time_se <- summary(bin_out_time)$coefficients[2,2]
          bin_out_time_pval <- summary(bin_out_time)$coefficients[2,4]
        } else {
          bin_out_time_coef <- NA
          bin_out_time_se <- NA
          bin_out_time_pval <- NA
        }

        if(!('try-error') %in% bin_out_twoway) {
          bin_out_twoway_coef <- summary(bin_out_twoway)$coefficients[2,1]
          bin_out_twoway_se <- summary(bin_out_twoway)$coefficients[2,2]
          bin_out_twoway_pval <- summary(bin_out_twoway)$coefficients[2,4]
        } else {
          bin_out_twoway_coef <- NA
          bin_out_twoway_se <- NA
          bin_out_twoway_pval <- NA
        }

      } else {
        bin_out_case_coef <- NA
        bin_out_case_se <- NA
        bin_out_case_pval <- NA
        bin_out_time_coef <- NA
        bin_out_time_se <- NA
        bin_out_time_pval <- NA
        bin_out_twoway_coef <- NA
        bin_out_twoway_se <- NA
        bin_out_twoway_pval <- NA
      }

      ### Compile saved results in a data frame
      results <- tibble(
            model = c("Two-way FE", "Case FE", "Time FE", "Pooled OLS",
                      "RE (u_i)","RE (v_t)","Generalized\nSynthetic\nControl","Weighted\nFixed Effects",
                      "Case Logit","Time Logit","Two-way Logit"),
            coef = c(twcoef,
                     #summary(case.fe.dummy)$coefficients[2,1],
                     summary(case.fe.plm)$coefficients[1,1],
                     #summary(time.fe.dummy)$coefficients[2,1],
                     summary(time.fe.plm)$coefficients[1,1],
                     summary(pooled.ols)$coefficients[2,1],
                     random_time_coef,
                     random_xs_coef,
                     coef_gsynth,
                     coef_wfe,
                     bin_out_case_coef,
                     bin_out_time_coef,
                     bin_out_twoway_coef),
            se = c(twse,
                  # summary(case.fe.dummy)$coefficients[2,2],
                   summary(case.fe.plm)$coefficients[1,2],
                   #summary(time.fe.dummy)$coefficients[2,2],
                   summary(time.fe.plm)$coefficients[1,2],
                   summary(pooled.ols)$coefficients[2,2],
                  random_time_se,
                  random_xs_se,
                  se_gsynth,
                  se_wfe,
                  bin_out_case_se,
                  bin_out_time_se,
                  bin_out_twoway_se),
            p_val = c(tw_pval,
                      summary(case.fe.plm)$coefficients[1,4],
                      #summary(time.fe.dummy)$coefficients[2,2],
                      summary(time.fe.plm)$coefficients[1,4],
                      summary(pooled.ols)$coefficients[2,4],
                      random_time_pval,
                      random_xs_pval,
                      p_gsynth,
                      p_wfe,
                      bin_out_case_pval,
                      bin_out_time_pval,
                      bin_out_twoway_pval),
            tw_id=identified,
            gsynth_factor=r_gsynth
      )

      return(results)
}

#' Sumarize results of \code{tw_sim} simulations
twtransform <- function(d){

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

#' Function to calculate coefficients for simulations
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

#' Run Monte Carlo Simulations on Panel Data
#'
#' Use this function to run your own Monte Carlo experiments on panel data.
#'
#' The \code{tw_sim} function acts as a frontend to the \code{link{tw_data}} function that
#' generates panel data given a set of parameters. The arguments are fairly self-explanatory, and
#' any underlying panel data parameters are passed on to \code{link{tw_data}}, so review that help
#' site for more information. The simulation works by passing a character value of one of the
#' parameters from the \code{link{tw_data}} function to the \code{arg} option, such as \code{cross.eff.mean}
#' for varying cross-sectional effects. The \code{at} option must take a numeric vector giving the range of
#' values to iterate over, such as \code{seq(-1,1,by=.1)} to test values of \code{cross.eff.mean} from
#' -1 to 1 in increments of .1.
#'
#' @return This function returns a \code{data.frame} with the two-way fixed effects coefficient,
#' both case (over-time) and cross-sectional 1-way fixed effects estimates, pooled OLS model
#' estimates, and a random-effects estimate. Each row of the \code{data.frame} represents
#' one sample of panel data drawn from \code{link{tw_data}}. This data can then be plotted with the
#' \code{\link{tw_plot}} function to see results of the simulation.
#'
#' @examples
#' \dontrun{
#'  ex_sim <- tw_sim(iter=150,arg='case.eff.mean',at=seq(-5,5,by=.1))
#'
#' }
#'
#' @param iter The number of Monte Carlo simulations to run. This effects the precision of the estimates.
#' @param cores The number of cores to run the simulations in parallel.
#' @param parallel whether to use parallel processing if \code{cores>1}.
#' @param arg A character value of which model parameter from \code{link{tw_data}} to vary across simulations.
#' @param at A numeric sequence giving the range of values of the \code{arg} parameter to vary across.
#' @param models A list of additional models besides fixed effects models to fit to the data.
#' @param ... All additional parameters to pass on to \code{link{tw_data}}, such as
#' average effects in case or cross-sectional dimensions of the data.
#'
#' @importFrom plm plm pdata.frame
#' @export
tw_sim <- function(iter=1000, cores=1, parallel=FALSE, arg='cross.eff.mean',
                   range1=c(-1,1),
                   range2=NULL, models=c("randomfx","wfe"), ...){


      # sample from uniform bounds

      range1 <- runif(iter,min=range1[1],max=range1[2])

      if(!is.null(range2)) range2 <- runif(iter,min=range2[1],max=range2[2])


      if(is.null(range2)) {
        el <- expand.grid(iteration = 1:iter, input1 = range1, arg=arg)
        # split to help with parallel processing
        split_el <- split(el,list(el$input,el$iteration))
      } else {
        el <- expand.grid(iteration = 1:iter, input1 = range1, input2=range2,arg=arg)
        # split to help with parallel processing
        split_el <- split(el,list(el$input1,el$input2,el$iteration))
      }

      # check for binary outcomes

      check_args <- list(...)

      if(!is.null(check_args$binary_outcome) &&
                  check_args$binary_outcome==101) {
        binary_outcome <- T
      } else {
        binary_outcome <- F
      }


      if(parallel){
            #if(file.exists('output.txt')) file.remove('output.txt')
            #loopfun <- ifelse(verbose, pbmclapply::pbmclapply, parallel::mclapply)
            if(is.null(cores)) cores <- parallel::detectCores()
            p <- mclapply(split_el, FUN=function(el){
                  sink('output.txt',append = T)
                    print(paste0('Now on iteration ',el$iteration[1], ' of input ',
                           el$input1[1], ' for parameter ',el$arg[1],collapse=" "))
                  sink()

                  z <- list(foo = el$input1[1])
                  names(z) <- el$arg[1]
                  extra.args <- as.list(substitute(list(...)))[-1L]
                  if(length(el$arg)>1) {
                    other_z <- list(foo = el$input2[1])
                    names(other_z) <- el$arg[2]
                    z<-c(z,other_z,extra.args)
                  } else {
                    z<-c(z,extra.args)
                  }
                  #z$gsynth <- "gsynth" %in% models
                  d <- do.call(tw_data, args=z)
                  m <- tw_model(d, models=models,binary_outcome=binary_outcome)
                  m$iteration <- el$iteration

                  out <- data_frame(input1=unique(el$input1),
                                    iteration=unique(el$iteration))
                  names(out) <- c(as.character(el$arg[1]),'iteration')

                  if(length(el$arg)>1) {
                    out <- mutate(out,
                                  input2=unique(el$input2))
                    names(out)[names(out)=='input2'] <- as.character(el$arg[2])
                  }
                  out <- left_join(out,m,by='iteration')

            }, mc.cores=cores)  %>% bind_rows

      } else {

        p <- lapply(split_el, FUN=function(el){

          print(paste0('Now on iteration ',el$iteration[1], ' of input ',
                       el$input1[1], ' for parameter ',el$arg[1],collapse=" "))


          z <- list(foo = el$input1[1])
          names(z) <- el$arg[1]
          extra.args <- as.list(substitute(list(...)))[-1L]
          if(length(el$arg)>1) {
            other_z <- list(foo = el$input2[1])
            names(other_z) <- el$arg[2]
            z<-c(z,other_z,extra.args)
          } else {
            z<-c(z,extra.args)
          }

          #z$gsynth <- "gsynth" %in% models

          d <- do.call(tw_data, args=z)
          m <- tw_model(d, models=models,binary_outcome=binary_outcome)
          m$iteration <- el$iteration

          out <- data_frame(input1=unique(el$input1),
                            iteration=unique(el$iteration))
          names(out) <- c(as.character(el$arg[1]),'iteration')

          if(length(el$arg)>1) {
            out <- mutate(out,
                          input2=unique(el$input2))
            names(out)[names(out)=='input2'] <- as.character(el$arg[2])
          }
          out <- left_join(out,m,by='iteration')

        })  %>% bind_rows
      }

      return(p)
}

#' Plot Results of \code{tw_sim} simulations
#'
#' This function will generate plots comparing the performance of different estimators
#' at returning \code{\link{tw_sim}} Monte Carlo simulations.
#'
#' @param gensim A set of simulations produced by \code{\link{tw_sim}}
#' @param display_est Whether to display estimated coefficients of
#' independent variables from the linear model (\code{"coef"}) or the
#' standard error of these estimates (\code{"se"})
#' @param use_ci Whether to add confidence intervals around LOESS curves.
#' Turn off if the CIs are too big relative to plot.
#' @param xvar In case of multiple attributes varied in the simulation, which
#' attribute should be on the x axis? See \code{\link{tw_data}} for possible
#' attribute names. Should be an unquoted variable name like
#' \code{case.int.sd}.
#' @param yvar In case of multiple attributes varied in the simulation, which
#' attribute should be on the y axis? See \code{\link{tw_data}} for possible
#' attribute names. Should be an unquoted variable name like
#' \code{case.int.sd}.
#' @param facet_scales The value of \code{scales} option passed to the
#' \code{ggplot2} \code{facet_wrap} function to control the scales.
#' Changes to \code{'free_x'},\code{'free_y'} or \code{'free'} if the
#' scales obscure comparison
#' @param ... Currently ignored
#'
#' @import ggplot2
#' @export
tw_plot <- function(gensim,
                    display_est='coef',
                    use_ci=T,
                    xvar=NULL,
                    yvar=NULL,
                    facet_scales="fixed",
                    ...) {


  att_data <- select(gensim,-iteration,-model,-coef,-se,-p_val,-tw_id,-gsynth_factor)
  att_names <- names(att_data)


  # different plots for whether we have one or two contrasts

  if(length(att_names)>1) {

    if(is.null(xvar) || is.null(yvar)) {
      gensim %>%
        ggplot(aes(y=!! sym(att_names[1]),
                    x=!! sym(att_names[2]))) +
          geom_raster(aes(fill=!! sym(display_est))) +
          facet_wrap(~model,scales=facet_scales)
    } else {
      gensim %>%
        ggplot(aes_(y=!! enquo(yvar),
                    x=!! enquo(xvar))) +
        geom_raster(aes(fill=coef)) +
        facet_wrap(~model,scales=facet_scales)
    }
  } else {

    gensim %>%  ggplot(aes(y= !! sym(display_est),
                                          x=!! sym(att_names))) +
      geom_abline(slope=1,intercept=0,linetype=2,alpha=0.5) +
      geom_smooth(se=use_ci,aes(colour=model,linetype=model)) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face="bold")) +
      ylab("Estimated Coefficient") +
      xlab(paste0("Fixed Value for ",att_names)) +
      facet_wrap(~model,scales=facet_scales)
  }

}

#' Would like to have a function that can generate power curves
tw_power <- function(gensim,
                     ...) {

}

#' Function to determine whether to sample from binary outcome or the real line
det_outcome <- function(var,treat_effect=NULL,x=NULL) {
  if(treat_effect && x) {
    as.numeric(var>median(var))
  } else if(treat_effect && !x) {
    as.numeric(runif(length(var))<plogis(var))
  } else {
    var
  }
}

