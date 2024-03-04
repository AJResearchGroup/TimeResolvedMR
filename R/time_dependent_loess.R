#' Calculate time-dependent effects of a polygenic score on a phenotype using
#' a LOESS model
#'
#' @param pgs Numeric vector of polygenic scores
#' @param pheno Numeric vector of phenotype measurements
#' @param age Numeric vector of age at assessment
#' @param covariates Data frame of covariates
#' @param age_range Age range (min, max) for effect estimation or \code{NULL} for
#'   automatic detection
#' @param age_step Number of years between each age stratum (default 1)
#'
#' @details
#'   By default, \code{age_range} is inferred to be \code{range(data[, age])},
#'   i.e. span from the minimum to the maximum age. The range is start- and
#'   end-\emph{inclusive}. If the value does not represent a valid range, i.e.
#'   element 2 is smaller than element 1, an error will be emitted. If the age
#'   range only includes one age (i.e. element 1 = element 2), you will be warned
#'   but the regression is still performed.
#' @return An instance of \code{\link{LoessTimeDependentModel}} containing a
#'   model of the relationship between age and genetic effect and all
#'   age-stratified models generated in the calculation
#' @seealso [LoessTimeDependentModel()], [stats::loess()], [time_dependent_glm()]
#' @export
#'
#' @examples
#' # Generate random example data
#' library(stats)
#' pgs <- rnorm(1000)
#' pheno <- pgs * 3
#' age <- rep(40:70, length.out = 1000)
#' covariates <- data.frame(covar = rnorm(1000))
#'
#' # Infer age range from data
#' time_dependent_loess(pgs, pheno, age, covariates)
#'
#' # Use explicit age range from 50 to 65
#' time_dependent_loess(pgs, pheno, age, covariates, age_range = c(50,65))
#'
#' # Only estimate effects for every other year
#' time_dependent_loess(pgs, pheno, age, covariates, age_step = 2)
time_dependent_loess <- function(pgs, pheno, age, covariates,
                                 age_range = NULL, age_step = 1) {
  # Guess Age range if not given
  if (is.null(age_range)) {
    age_range <- range(age, na.rm = TRUE)
  }
  if (any(is.na(age_range)) || (age_range[1] > age_range[2])) {
    rlang::abort(stringr::str_glue("Invalid age range: {age_range}"))
  }
  if (age_range[1] == age_range[2]) rlang::warn("Age range is only one age")

  age_seq <- seq(from = age_range[1], to = age_range[2], by = age_step)
  # Regress min age and before
  below_min_age <- list(age <= age_seq[1])
  # Each age inbetween
  at_specific_age <- lapply(age_seq[-c(1, length(age_seq))], \(a) age == a)
  # Max age and after
  above_max_age <- list(age >= age_seq[length(age_seq)])

  age_filters <- c(below_min_age, at_specific_age, above_max_age)

  age_formula <- reformulate(
    termlabels = c("pgs", colnames(covariates)),
    response = "pheno"
  )
  data <- covariates
  data$pgs <- pgs
  data$pheno <- pheno
  data$age <- age

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

  age_coefficients <- model_list |> vapply(\(model) coef(model)["pgs"], numeric(1))
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
  LoessTimeDependentModel(
    model = loess_model,
    stratifiedModels = model_list
  )
}
