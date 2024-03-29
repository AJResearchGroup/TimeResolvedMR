% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_all_in_one.R
\name{analyze_all_in_one}
\alias{analyze_all_in_one}
\title{Convenience function to estimate time-dependent effects on exposure and
outcome}
\usage{
analyze_all_in_one(
  pgs,
  exposure,
  outcome,
  exposure_age,
  outcome_age,
  covariates,
  exposure_modeltype = c("loess", "glm", "aalen"),
  interaction_exponents = 1,
  age_range = NULL,
  age_step = 1
)
}
\arguments{
\item{pgs}{Numeric vector of olygenic scores of exposure}

\item{exposure}{Numeric vector of exposure measurements or status indicator}

\item{outcome}{Numeric vector of outcome measurements or status indicator}

\item{exposure_age}{Numeric vector of age at exposure assessment or
follow-up for exposure}

\item{outcome_age}{Numeric vector of age at outcome assessment or follow-up
of outcome}

\item{covariates}{Covariates to include in the model}

\item{exposure_modeltype}{Type of the model to use for regression of exposure.
Either \code{loess} (default), \code{glm} or \code{aalen}}

\item{interaction_exponents}{Exponents to use for PGS-time interaction terms
in \code{glm} model. Otherwise unused.}

\item{age_range}{Range between which ages genetic effects on the exposure are
to be calculated for LOESS model. Otherwise unused. \code{NULL} means the range
should be inferred from the data. (Default: \code{NULL})}

\item{age_step}{Step size between ages for LOESS model. Otherwise unused.
(Default: 1)}
}
\value{
A named list containing the models for exposure and outcome as
elements named \code{exposure} and \code{outcome}. Depending on the requested model,
these will be of type \code{\link[=AalenTimeDependentModel]{AalenTimeDependentModel()}}, \code{\link[=GlmTimeDependentModel]{GlmTimeDependentModel()}}
or \code{\link[=LoessTimeDependentModel]{LoessTimeDependentModel()}}.
}
\description{
Convenience function to estimate time-dependent effects on exposure and
outcome
}
\details{
The function performs all analyses internally which can also be
carried out by the other functions in this package.
It is intended for quick one-off analyses as the other functions allow you
to reuse models.
The current implementation only allows for binary outcomes and has been
tested on continuous exposures, only.
Individuals with a recorded outcome event before the age at exposure
assessment will be removed from the model.
}
\examples{
# Generate random example data
library(stats)
pgs <- rnorm(1000)
exposure <- pgs * 3
exposure_age <- rep(40:70, length.out = 1000)
covariates <- data.frame(covar = rnorm(1000))
outcome <- sample(c(TRUE,FALSE), 1000, replace=TRUE)
outcome_age <- rnorm(1000, mean=60, sd = 3) |> pmax(40)

# Default: LOESS model for exposure and Aalen for outcome. Age range is
# inferred from data.
analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
  covariates)

# Only estimate effects every other year between ages 50 and 65
analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
  covariates, age_range = c(50,65), age_step = 2)

# Use GLM with linear pgs-age interaction instead for exposure
analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
  covariates, exposure_modeltype = "glm")

# GLM for exposure with linear and quartic interactions
analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
  covariates, exposure_modeltype = "glm", interaction_exponents = c(1,4))

}
\seealso{
\code{\link[=time_dependent_loess]{time_dependent_loess()}}, \code{\link[=time_dependent_glm]{time_dependent_glm()}},
\code{\link[=time_dependent_aalen]{time_dependent_aalen()}}
}
