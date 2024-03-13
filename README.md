
# TimeResolvedMR

<!-- badges: start -->
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/AJResearchGroup/TimeResolvedMR/r.yml)
![GitHub License](https://img.shields.io/github/license/AJResearchGroup/TimeResolvedMR)
<!-- badges: end -->

Estimate time-varying effects of an exposure on an outcome using genetic
instrumental variables through Mendelian1 Randomization

## Installation

You can install the development version of TimeResolvedMR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AJResearchGroup/TimeResolvedMR")
```

## Example

The following is the basic workflow we used in our upcoming paper.
This estimates the time-dependent effect of a continuous exposure (BMI) on a
binary outcome (Type 2 diabetes). The time-dependent genetic effects are
estimated using a generalized linear model including interactions PGS:age and
PGS:(age^4) and an Aalen additive hazards model.

``` r
library(TimeResolvedMR)

# Assuming you got polygenic scores and all exposure,outcome and covariate
# measurements you need

exposure_model <- time_dependent_glm(
  pgs = pgs,
  pheno = bmi,
  age =  age_at_assessment,
  covariates = covariates,
  exponents = c(1,4)
)

outcome_model <- time_dependent_aalen(
  pgs = pgs,
  pheno = t2_diabetes,
  event_age = age_at_diagnose,
  covariates = covariates
)

time_dependent_MR(exposure_model, outcome_model)
```

