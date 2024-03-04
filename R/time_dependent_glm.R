#' Calculate time-dependent effects of a PGS on an exposure using a glm with
#' interaction terms
#'
#' @param pgs Numeric vector of polygenic scores
#' @param pheno Numeric vector of phenotype measurements
#' @param age Numeric vector of age at assessment
#' @param covariates Data frame of covariates
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
#' @return An instance of \code{\link{GlmTimeDependentModel}} containing a
#' \code{\link[stats]{glm}} modeling the relationship and extracted coefficients
#'
#' @export
#'
#' @examples
time_dependent_glm <- function(pgs, pheno, age, covariates, exponents = 1) {
  # Include PGS, PGS-Age interaction with arbitrary exponents and all covariates
  # Exponent = NULL is equivalent to no age interaction. This is on purpose
  interaction_terms <- stringr::str_glue("pgs:I(age^{exponents})")
  model_formula <- reformulate(
    termlabels = c("pgs", interaction_terms, colnames(covariates)),
    response = "pheno"
  )
  model <- glm(model_formula, family = "gaussian", data = covariates)
  GlmTimeDependentModel(
    model = model,
    exponents = exponents,
    fixedEffect = coef(model)["pgs"],
    interactionEffects = coef(model)[interaction_terms]
  )
}
