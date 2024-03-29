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
      pgs = pgs, event = pheno, event_age = age, covariates = covariates
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
#' @param interaction_exponents Exponents to use for PGS-time interaction terms
#'   in `glm` model. Otherwise unused.
#' @param age_range Range between which ages genetic effects on the exposure are
#'   to be calculated for LOESS model. Otherwise unused. `NULL` means the range
#'   should be inferred from the data. (Default: `NULL`)
#' @param age_step Step size between ages for LOESS model. Otherwise unused.
#'   (Default: 1)
#'
#' @details The function performs all analyses internally which can also be
#'   carried out by the other functions in this package.
#'   It is intended for quick one-off analyses as the other functions allow you
#'   to reuse models.
#'   The current implementation only allows for binary outcomes and has been
#'   tested on continuous exposures, only.
#'   Individuals with a recorded outcome event before the age at exposure
#'   assessment will be removed from the model.
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
#' # Generate random example data
#' library(stats)
#' pgs <- rnorm(1000)
#' exposure <- pgs * 3
#' exposure_age <- rep(40:70, length.out = 1000)
#' covariates <- data.frame(covar = rnorm(1000))
#' outcome <- sample(c(TRUE,FALSE), 1000, replace=TRUE)
#' outcome_age <- rnorm(1000, mean=60, sd = 3) |> pmax(40)
#'
#' # Default: LOESS model for exposure and Aalen for outcome. Age range is
#' # inferred from data.
#' analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
#'   covariates)
#'
#' # Only estimate effects every other year between ages 50 and 65
#' analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
#'   covariates, age_range = c(50,65), age_step = 2)
#'
#' # Use GLM with linear pgs-age interaction instead for exposure
#' analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
#'   covariates, exposure_modeltype = "glm")
#'
#' # GLM for exposure with linear and quartic interactions
#' analyze_all_in_one(pgs, exposure, outcome, exposure_age, outcome_age,
#'   covariates, exposure_modeltype = "glm", interaction_exponents = c(1,4))
#'
analyze_all_in_one <- function(
    pgs, exposure, outcome, exposure_age, outcome_age, covariates,
    exposure_modeltype = c("loess", "glm", "aalen"),
    interaction_exponents = 1, age_range = NULL, age_step = 1) {

  exposure_modeltype <- match.arg(exposure_modeltype)
  outcome_modeltype <- "aalen"

  exposure_model <- select_model(
    model_type = exposure_modeltype, pgs = pgs, pheno = exposure, age = exposure_age,
    covariates = covariates, interaction_exponents = interaction_exponents,
    age_range = age_range, age_step = age_step
  )

  outcome_model <- select_model(
    model_type = outcome_modeltype, pgs = pgs, pheno = outcome, age = outcome_age,
    covariates = covariates, interaction_exponents = interaction_exponents,
    age_range = age_range, age_step = age_step
  )
  age_seq <- if(!is.null(age_range))
    seq(age_range[1], age_range[2], by = age_step)
    else 0:max(outcome_age)
  time_dependent_MR(
    age_seq = age_seq, exposure_model = exposure_model,
    outcome_model = outcome_model, method = "trapezoidal")
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
#' # Random example data
#' library(stats)
#' data <- data.frame(
#'   pgs = rnorm(50),
#'   outcome = sample(c(TRUE,FALSE), 50, replace=TRUE),
#'   outcome_age = rnorm(50, mean = 50, sd = 5) |> pmax(20),
#'   expoure_age = rnorm(50, mean = 30, sd = 2) |> pmax(15)
#' )
#'
#' # Columns are already named the same as default values
#' remove_outcome_cases <- function(data)
#'
#' # Case where one column has a different name
#' data$event <- data$outcome
#' data$outcome <- NULL
#' remove_outcome_cases <- function(data, outcome = "event")
remove_outcome_cases <- function(data, outcome = "outcome",
                                 outcome_age = "outcome_age",
                                 exposure_age = "exposure_age") {
  prior_event <-  # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
    as.logical(data[, outcome]) &
    (data[, outcome_age] < data[, exposure_age])
  data[!prior_event, ]
}
