#Parameters
#TODO replace with Roxygen
# pgs : Polygenic risk scores (hopefully a data frame)
# exp.data: some sort of exposure data. Measurements maybe?
# data: outcome data (see above)
# Merged PGS: if TRUE, compute merged PGS, otherwise subcohort (??)
# w.outtype: Either binary or continuous
# age.step: Step size for age sequence
calculate_genetic_effects <- function(
    data, pgs, exposure, outcome, time_at_event, age_at_assessment,
    covariates = setdiff(colnames(data), c(exposure, outcome, time_at_event, age_at_assessment)),
    interaction_exponents = 1, age_range = NULL, age.step = 1,
    merged_pgs = TRUE) {

  # Make sure we are dealing with column names
  stopifnot(
    "pgs must be a column name in data" = pgs %in% colnames(data),
    "exposure must be a column name in data" = exposure %in% colnames(data),
    "outcome must be a column name in data" = outcome %in% colnames(data),
    "time_at_event must be a column name in data" = time_at_event %in% colnames(data),
    "age_at_assessment must be a column name in data" = age_at_assessment %in% colnames(data)
  )

  # Guess Age range if not given
  if (is.null(age_range)) {
    age_range <- range(data[, age_at_assessment], na.rm = TRUE)
  }
  if (is.na(age_range) || (age_range[0] > age_range[1])) {
    stop(stringr::str_glue("Invalid age range: {age_range}"))
  }
  if (age_range[0] == age_range[1]) warning("Age range is only one age")

  # Set up age bounds
  prior_event <-  # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
    (data[, exposure] == 1) &
    (data[, time_at_event] < data[, age_at_assessment])

  data_clean <- na.omit(data)
  data_clean$outcome_controls <- !dis

  healthy <- data_clean$outcome_controls

  # Include PGS, PGS-Age interaction with arbitrary exponents and all covariates
  # Exponent = NULL is equivalent to no age interaction. This is on purpose
  interaction_terms <- stringr::str_glue("{pgs}:I({age_at_assessment}^{exponents})")
  model_formula <- reformulate(
    termlabels = c(pgs, interaction_terms, covariates),
    response = exposure
  )
  model_exposure <- glm(
      model_formula,
      family = "gaussian",
      subset = healthy,
      data = data_clean
  )

  time_dependent_coefficients <- coef(model_exposure)[c(pgs, interaction_terms)]

  # local regression (LOESS) fitting to age-stratified data:
  age_seq <- seq(from = age_range[1], to = age_range[2], by = age.step)
  coef.exp.age <-
    matrix(rep(NA, length(age_seq) * 3), length(age_seq), 3)
  # Regress min age and before
  below_min_age <- list(data_clean[, age_at_assessment] <= age_seq[1])
  # Each age inbetween
  at_specific_age <- lapply(age_seq[c(-1, length(age_seq))], \(age) {
    data_clean[, age_at_assessment] == age
  })
  # Max age and after
  above_max_age <- list(data_clean[, age_at_assessment] >= age_seq[length(age_seq)])

  age_filters <- c(below_min_age, at_specific_age, above_max_age)

  age_formula <- reformulate(
    termlabels = c(pgs, covariats),
    response = exposure
  )

  coefficients_list <- lapply(age_filters, \(filter) {
    model <- glm(
      age_formula,
      family = "gaussian",
      subset = healthy & filter,
      data = exp.data.age
    )
    summary(model)$coef[pgs, c(1,2,4)] |> as.data.frame()
  })
  age_coefficients <- do.call(rbind, coefficients_list)
  names(age_coefficients) <- c("beta", "se", "p")
  age_coefficients$age <- age_seq

  # time-dependent (disease-specific) effects are estimated from a LOESS fit to age-stratified effects (NOTE: outcome cases are removed prior to effect estimation)
  loess_model <- loess(
    beta ~ age,
    control = loess.control(surface = "direct"),
    family = "gaussian",
    span = 0.50,
    degree = 1,
    data = age_coefficients
  )   # linear extrapolation

  # PGS effect on outcome
  outcome_formula <- reformulate(
    termlabels = c(pgs, covariates),
    response = stringr::str_glue("Surv({time_at_event}, {outcome})")
  )
  outcome_model <- timereg::aalen(outcome_formula, data = covars)

  list(
    timedep = time_dependent_coefficients,
    byage = age_coefficients,
    outcome = outcome_model
  )
}
