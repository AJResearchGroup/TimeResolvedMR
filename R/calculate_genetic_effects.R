#' Internal function for model selection
#'
#' @param model_type The model to run
#' @param ... Arguments which are passed to the model
#'
#' @return The desired model
#' @keywords internal
select_model <- function(model_type, pgs, pheno, age, covariates, ...) {
  dots <- list(...)
  switch(model_type,
    loess = time_dependent_loess(
      pgs = pgs, pheno = pheno, age = age, covariates = covariates,
      age_range <- dots[["age_range"]], age_step = dots[["age_step"]]),
    glm = time_dependent_glm(
      pgs = pgs, pheno = pheno, age = age, covariates = covariates,
      exponents = dots[["exponents"]]
    ),
    aalen = time_dependent_aalen(
      pgs = pgs, event = pheno, event_age = age, covariates = covariates,
    ),
    rlang::abort(paste0("Unknown model type ", model_type))
  )
}

#' Convenience function to estimate time-dependent effects on exposure and
#' outcome
#'
#' @param pgs Numeric vector of olygenic scores of exposure
#' @param exposure Numeric vector of exposure measurements or status indicator
#' @param outcome Numeric vector of outcome measurements or status indicator
#' @param exposure_age Numeric vector of age at exposure assessment or
#'   follow-up for exposure
#' @param outcome_age Numeric vector of age at outcome assessment or follow-up
#'   of outcome
#' @param covariates Covariates to include in the model
#' @param exposure_modeltype Type of the model to use for regression of exposure.
#'   Either `loess` (default), `glm` or `aalen`
#' @param exposure_modeltype Type of the model to use for regression of outcome.
#'   Either `aalen` (default), `loess` or `glm`
#' @param interaction_exponents Exponents to use for PGS-time interaction terms
#'   in `glm` model. Otherwise unused
#' @param age_range Range between which ages genetic effects on the exposure are
#'   to be calculated for LOESS model. Otherwise unused. `NULL` means the range
#'   should be inferred from the data. (Default: `NULL`)
#' @param age_step Step size between ages for LOESS model. Otherwise unused.
#'   (Default: 1)
#'
#' @return A named list containing the models for exposure and outcome as
#'   elements named `exposure` and `outcome`. Depending on the requested model,
#'   these will be of type [AalenTimeDependentModel()], [GlmTimeDependentModel()]
#'   or [LoessTimeDependentModel()].
#' @seealso
#'   [time_dependent_loess()], [time_dependent_glm()],
#'   [time_dependent_aalen()]
#' @export
#'
#' @examples
calculate_genetic_effects <- function(
    pgs, exposure, outcome, exposure_age, outcome_age, covariates,
    exposure_modeltype = c("loess", "glm", "aalen"),
    outcome_modeltype = c("aalen", "loess", "glm"),
    interaction_exponents = 1, age_range = NULL, age_step = 1) {

  exposure_modeltype <- match.arg(exposure_modeltype)
  outcome_modeltype <- match.arg(outcome_modeltype)

  exposure_model <- select_model(
    model_type = exposure_model, pgs = pgs, pheno = exposure, age = exposure_age,
    covariates = covariates, interaction_exponents = interaction_exponents,
    age_range = age_range, age_step = age_step
  )

  outcome_model <- select_model(
    model_type = outcome_model, pgs = pgs, pheno = outcome, age = outcome_age,
    covariates = covariates, interaction_exponents = interaction_exponents,
    age_range = age_range, age_step = age_step
  )
  list(
    exposure = exposure_time_dependent_model,
    outcome = outcome_time_dependent_model
  )
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

