#' Estimate time-dependent effects on exposure and effects on outcome
#'
#' @param data Data frame to be used for regression, includes PGS, outcome,
#'   exposure and covariate data.
#' @param pgs Column name of the polygenic score
#' @param exposure Column name of the exposure measurements
#' @param outcome Column name of the outcome measurements
#' @param exposure_age Column name of age at assessment of exposure
#' @param outcome_age Column name of time at event (outcome), if any
#' @param covariates Column names of covariates
#' @param interaction_exponents Exponents to use for PGS-time interaction
#' @param age_range Range between which ages genetic effects on the exposure are
#'   to be calculated
#' @param age_step Step size between ages
#'
#' @return A list containing estimates
#' @export
#'
#' @examples
calculate_genetic_effects <- function(
    data, pgs = "pgs", exposure = "exposure", outcome = "outcome",
    exposure_age = "exposure_age", outcome_age = "outcome_age",
    covariates = setdiff(colnames(data), c(exposure, outcome, outcome_age, exposure_age)),
    interaction_exponents = 1, age_range = NULL, age_step = 1) {

  # Make sure we are dealing with column names
  stopifnot(
    "pgs must be a column name in data" = pgs %in% colnames(data),
    "exposure must be a column name in data" = exposure %in% colnames(data),
    "outcome must be a column name in data" = outcome %in% colnames(data),
    "outcome_age must be a column name in data" = outcome_age %in% colnames(data),
    "exposure_age must be a column name in data" = exposure_age %in% colnames(data)
  )
  data_clean <- remove_outcome_cases(
    data = data, outcome = outcome, outcome_age = outcome_age,
    exposure_age = exposure_age
  )
  time_dependent_coefficients <- get_time_dependent_effects(
    data = data_clean, pgs = pgs, exposure = exposure, covariates = covariates,
    exposure_age = exposure_age, exponents = interaction_exponents
  )
  age_models <- get_age_stratified_effects(
    data = data_clean, pgs = pgs, exposure = exposure, covariates = covariates,
    exposure_age = exposure_age, age_range = age_range,
    age_step = age_step
  )
  outcome_model <- regress_outcome(data = data_clean, pgs = pgs, covariates = covariates)

  list(
    timedep = time_dependent_coefficients,
    byage = age_models,
    outcome = outcome_model
  )
}

#' Calculate time-dependent effects of a PGS on an exposure
#'
#' @param data Individual-level data in a data frame
#' @param pgs Column name of PGS
#' @param covariates Column names of covariates to include
#' @param exposure Column name of exposure
#' @param exposure_age Column name of age when exposure was measured
#' @param exponents A numeric vector of exponents for the PGS-Age interaction term
#'   (Default 1)
#'
#' @details The relationship is always modeled as
#'   \code{exposure ~ pgs + age + covariates + interaction_terms} where
#'   \code{interaction_terms} is an arbitrary number of PGS-Age interactions.
#'   These interactions take the form \code{pgs:(age^exponent)} for each
#'   exponent in \code{exponents}. Therefore, it is possible to include multiple
#'   interaction terms to model, e.g., non-linear interactions. A value of
#'   \code{NULL} causes no interaction term to be included.
#'
#' @return A linear model of type \code{\link[stats]{glm}} with the exposure
#'   as response. The predictors are the PGS, covariates, age and an arbitrary
#'   number of PGS-Age interaction terms.
#' @export
#'
#' @examples
get_time_dependent_effects <- function(data, pgs = "pgs", exposure = "exposure",
                                       exposure_age = "exposure_age",
                                       covariates = setdiff(colnames(data), c(exposure, outcome, exposure_age)),
                                       exponents = 1) {
  # Include PGS, PGS-Age interaction with arbitrary exponents and all covariates
  # Exponent = NULL is equivalent to no age interaction. This is on purpose
  interaction_terms <- stringr::str_glue("{pgs}:I({exposure_age}^{exponents})")
  model_formula <- reformulate(
    termlabels = c(pgs, interaction_terms, covariates),
    response = exposure
  )
  glm(model_formula, family = "gaussian", data = data)
}

#' Calculate age-stratified effects of a polygenic score on an exposure
#'
#' @param data Individual-level data in a data frame
#' @param pgs Column name of PGS
#' @param covariates Column names of covariates to include
#' @param exposure Column name of exposure
#' @param exposure_age Column name of age when exposure was measured
#' @param age_range Age range (min, max) for effect estimation or \code{NULL} for
#'   automatic detection
#' @param age_step Number of years between each age stratum (default 1)
#'
#' @details
#'   By default, \code{age_range} is inferred to be \code{range(data[, exposure_age])},
#'   i.e. span from the minimum to the maximum age. The range is start- and
#'   end-\emph{inclusive}. If the value does not represent a valid range, i.e.
#'   element 2 is smaller than element 1, an error will be emitted. If the age
#'   range only includes one age (i.e. element 1 = element 2), you will be warned
#'   but the regression is still performed.
#' @return A named list of two elements. \describe {
#'   \item{pgs_by_age}{
#'     A named list of \link[stats:glm]{linear models}
#'     stratified by \code{exposure_age}, One model for each age as given by
#'     \code{age_range} and \code{age_step}. The first model includes all
#'     individuals with \code{exposure_age} at or below the minimum age and the
#'     last model includes all individuals at or above the maximum age. Element
#'     names correspond to the age group.
#'   }
#'   \item{age_effect}{
#'     A \link[stats:loess]{LOESS} model of \code{beta ~ age}, where
#'     \code{beta} is the effect estimates of the PGS on the exposure at a
#'     certain age and \code{age} is the age at exposure assessment.
#'   }}
#' @export
#'
#' @examples
get_age_stratified_effects <- function(data, pgs = "pgs",
                                       exposure = "exposure",
                                       exposure_age = "exposure_age",
                                       covariates = setdiff(colnames(data), c(pgs, exposure, exposure_age)),
                                       age_range = NULL, age_step = 1) {
  # Guess Age range if not given
  if (is.null(age_range)) {
    age_range <- range(data[, exposure_age], na.rm = TRUE)
  }
  if (any(is.na(age_range)) || (age_range[1] > age_range[2])) {
    stop(stringr::str_glue("Invalid age range: {age_range}"))
  }
  if (age_range[1] == age_range[2]) warning("Age range is only one age")

  age_seq <- seq(from = age_range[1], to = age_range[2], by = age_step)
  # Regress min age and before
  below_min_age <- list(data[, exposure_age] <= age_seq[1])
  # Each age inbetween
  at_specific_age <- lapply(age_seq[-c(1, length(age_seq))], \(age) {
    data[, exposure_age] == age
  })
  # Max age and after
  above_max_age <- list(data[, exposure_age] >= age_seq[length(age_seq)])

  age_filters <- c(below_min_age, at_specific_age, above_max_age)

  age_formula <- reformulate(
    termlabels = c(pgs, covariates),
    response = exposure
  )

  model_list <- lapply(age_filters, \(age_filter) {
    glm(
      age_formula,
      family = "gaussian",
      data = data[age_filter, ]
    )
  })
  model_names <- age_seq
  model_names[length(model_names)] <- paste0(">=", model_names[length(model_names)])
  model_names[1] <- paste0("<=", model_names[1])
  names(model_list) <- model_names

  age_coefficients <- model_list |> lapply(\(model) coef(model)[pgs]) |> unlist()
  coef_and_age <- data.frame(beta = age_coefficients, age = age_seq)
  # time-dependent (disease-specific) effects are estimated from a LOESS fit to age-stratified effects
  loess_model <- loess(
    beta ~ age,
    control = loess.control(surface = "direct"),
    family = "gaussian",
    span = 0.50,
    degree = 1,
    data = coef_and_age
  )
  list(pgs_by_age = model_list, age_effect = loess_model)
}

#' Removes individuals with prior outcome events
#'
#' @param data A data frame containing the data set to be analyzed
#' @param outcome Column name of outcome status (i.e. had event or not)
#' @param outcome_age Column name of the age when the event occured
#' @param exposure_age Column name of the age when the exposure was measured
#'
#' @return A subset data frame where all observations where the outcome appeared
#'   before the exposure were removed
#'
#' @export
#'
#' @examples
remove_outcome_cases <- function(data, outcome = "outcome",
                                 outcome_age = "outcome_age",
                                 exposure_age = "exposure_age") {
  prior_event <-  # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
    as.logical(data[, outcome]) &
    (data[, outcome_age] < data[, exposure_age])
  data[!prior_event, ]
}

#' Subsets and renames columns to default values of analysis functions
#'
#' @param data A data frame containing the data to be analyzed
#' @param pgs Column name of polygenic scores
#' @param exposure Column name of exposure measurements.
#'   Can be case/control or continuous.
#' @param outcome Column name of outcome measurements.
#'   At this point only case/control
#' @param exposure_age Column name of age at which the exposure was measured/happened
#' @param outcome_age Column name of age at which the outcome happened
#' @param covariates Column names of covariates
#'
#' @return A renamed data frame with only columns specified by the function's
#'   arguments
#'
#' @details This function does not modify the original data frame.
#'   By default \code{covariates} is set to all column names that are not set as
#'   values of the other arguments.
#'
#' @export
#'
#' @examples
rename_columns <- function(data, pgs, exposure, outcome,
                           exposure_age, outcome_age,
                           covariates = setdiff(colnames(data), c(pgs, exposure, outcome, exposure_age, outcome_age))) {
  transformed <- data[, c(pgs, exposure, exposure_age, outcome_age, covariates)]
  colnames(transformed)[1:4] <- c("pgs", "exposure", "exposure_age", "outcome_age")
  transformed
}

#' Calculate cumulative PGS effect on outcome
#'
#' @param data The data to be analyzed
#' @param pgs Column name of PGS
#' @param outcome Column name of outcome status
#' @param outcome_age Column name of age when outcome happened
#' @param covariates Column names of covariates to include in the model
#'
#' @return An Aalen model of type \code{\link[timereg]{aalen}}
#' @export
#'
#' @examples
regress_outcome <- function(data, pgs = "pgs", outcome = "outcome",
                            outcome_age = "outcome_age",
                            covariates = setdiff(colnames(data), pgs)) {
  # PGS effect on outcome
  outcome_formula <- reformulate(
    termlabels = c(pgs, covariates),
    response = stringr::str_glue("survival::Surv({outcome_age}, {outcome})")
  )
  outcome_model <- timereg::aalen(outcome_formula, data = data)
}
