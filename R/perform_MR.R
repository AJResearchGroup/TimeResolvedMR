###########################################################################################
# Script 5: Causal estimation of the time-dependent, life-course effect of BMI on outcome #
###########################################################################################



integrate_estimates <- function(age, beta, method = c("trapezoidal", "midpoint")){
  stopifnot(
    is.numeric(age),
    is.numeric(beta),
  )
  method <- match.arg(method)
  effects <- switch(method,
    trapezoidal = pracma::cumtrapz(age, beta),
    midpoint =  c(0, cumsum(diff(age) * zoo::rollmean(beta, 2))),
    "You should not be here"
  )
  c(age, effects) |> matrix(ncol = 2, dimnames = list(NULL, c("age", "effect")))
}

calculate_mr_variance_trapz <- function(effect, variance) {
  stopifnot(
    is.numeric(effect),
    is.numeric(variance),
    length(effect) == length(variance),
    length(effect) > 0
  )
  n <- length(effect)
  var1 <- variance[1] / (16 * effect[1] ^ 2)
  if (n < 2) return(var1)

  var2 <- var1 + variance[2] / (4 * effect[2] ^ 2)
  if (n < 3) return (c(var1,var2))

  # At this point, we know that n > 2, so the following is safe:
  point_variances <- c(
    var1,
    vapply(2:n, \(k){ variance[k] / (4 * effect[k] ^2) }, numeric(1))
  )

  c(
    var1,
    var2,
    vapply(3:n, \(k) { var1 + point_variances[k-1] + point_variances[k] }, numeric(1))
  )
}

time_dependent_MR <- function(age.seq, exposure_model, outcome_model,
                              method = c("trapezoidal", "midpoint")) {
  method <- match.arg(method)
  midpoints <-  zoo::rollmean(age.seq, 2)
  total_exposure_effects <- totalEffect(exposure_model, age.seq)
  midpoint_exposure_effects <- totalEffect(exposure_mode, midpoints)
  outcome_effects <- totalEffect(outcome_model, age.seq)
  dBeta <- diff(outcome_effects) / diff(age.seq)
  wald_ratio <- dBeta / midpoint_effects
  Gamma <- integrate_estimates(midpoints, wald_ratio)
  Gamma_variance <- switch(method,
    trapezoidal = calculate_mr_variance_trapz(effect = total_effects, variance = beta[, "var"]),
    midpoint =  beta[, "var"] / total_effects ^ 2,
    "Method is neither trapezoidal nor midpoint. There should've been an error earlier."
  )
  data.frame(
    Gamma = Gamma,
    var = Gamma_variance,
    L95 = Gamma - 1.96 * sqrt(Gamma_variance),
    L95 = Gamma + 1.96 * sqrt(Gamma_variance)
  )

}

#Glin <- time_dependent_MR(age.seq = age.seq, model_type = "glm", betas = coef.exp.lin["pgs:age_at_assessment"], exponents = 1, fixed_effect = coef.exp.lin["pgs"])
#Gqdr <- time_dependent_MR(age.seq = age.seq, model_type = "glm",  betas = coef.exp.lin[c("pgs:age_at_assessment","pgs:I(age_at_assessment^4)")], exponents = c(1,4), fixed_effect = coef.exp.lin["pgs"])
#Gloe <- time_dependent_MR(age.seq = age.deq, model_type = "loess", loess_model = loess.m)
