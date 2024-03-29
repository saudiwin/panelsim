% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twsim2.R
\name{tw_model}
\alias{tw_model}
\title{Fit Panel Data Given Simulated Data}
\usage{
tw_model(gendata, models = c("randomfx", "wfe"), binary_outcome = FALSE)
}
\arguments{
\item{gendata}{A data frame of simulated panel data procued by
\code{\link{tw_data}}}

\item{models}{Whether to include random effects estimates in the
returned data frame (defaults to true).}
}
\value{
The function returns a \code{data.frame} with the coefficients
of all of the requested models (includes 1-way and 2-way FE models by default).
}
\description{
This function will take the result of data generated by the
\code{\link{tw_data}} function and fit panel data models
using the \code{\link[plm]{plm}} model-fitting function.
Available models include standard fixed effects models (1-way
and 2-way), random effects models and pooled OLS.
}
\examples{

\dontrun{
# generate panel data
gendata <- tw_data()
# estimate panel models
estimates <- tw_model(gendata)
}


}
