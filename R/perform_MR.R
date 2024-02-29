###########################################################################################
# Script 5: Causal estimation of the time-dependent, life-course effect of BMI on outcome #
###########################################################################################

load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/BMI_DATA_TRMR/output_data_bmi_t2dm_aalen_steiger_15152.RData",
  verbose = TRUE
)

# [,1]: PGS effect, [,2]: Interaction effect
coef.exp.lin          <- as.numeric(coef.exp.lin)
coef.exp.qdr          <- as.numeric(coef.exp.qdr)

pgs <- "pgs"

bcum.raw <-
  cbind(p.out$cum[, c("time", pgs)], p.out$var.cum[, pgs])   # time, B(t) and Var(B(t))
colnames(bcum.raw) <- c("time", "pgs", "pgs.var")

# interpolate B(t) and Var(B(t)) for non-existing t (ages), i.e., times with no events:
pgs_interpolated <- approx(x = floor(bcum.raw[, "time"]), y = bcum.raw[, "pgs"],
                           xout = age.seq, method = "constant", rule = 2)$y
pgs.var_interpolated <- approx(x = floor(bcum.raw[, "time"]), y = bcum.raw[, "pgs.var"],
                               xout = age.seq, method = "constant", rule = 2)$y
bcum <- c(age.seq, pgs_interpolated, pgs.var_interpolated) |>
  matrix(ncol = 3, dimnames = list(NULL, colnames(bcum.raw)))

# calculate b(t) = dB/dt and corresponding standard errors:
# Numeric differentiation, dx = diff(x) in this case
dPGS <- c(0, diff(bcum[, "pgs"]) / diff(bcum[, "time"]))
dPGS.var <- c(0, diff(bcum[, "pgs.var"]) / diff(bcum[, "time"]))
dPGS.low <- dPGS - 1.96 * sqrt(dPGS.var)
dPGS.high <- dPGS + 1.96 * sqrt(dPGS.var)
beta <- c(age.seq, dPGS, dPGS.var, dPGS.low, dPGS.high) |>
  matrix(ncol = 5, dimnames = list(NULL, c("time", "beta", "var", "l95", "h95")))

dBeta <- diff(beta[, "beta"]) / diff(beta[,"time"])    # numerical derivative (central difference), calculated in midpoint t = k - 1/2

Glin <- time_dependent_MR(age.seq = age.seq, model_type = "glm", betas = coef.exp.lin["pgs:age_at_assessment"], exponents = 1, fixed_effect = coef.exp.lin["pgs"])
Gqdr <- time_dependent_MR(age.seq = age.seq, model_type = "glm",  betas = coef.exp.lin[c("pgs:age_at_assessment","pgs:I(age_at_assessment^4)")], exponents = c(1,4), fixed_effect = coef.exp.lin["pgs"])
Gloe <- time_dependent_MR(age.seq = age.deq, model_type = "loess", loess_model = loess.m)

integrate_estimates <- function(age, beta, q.trap = TRUE){
  stopifnot(
    is.numeric(age),
    is.numeric(beta),
    is.logical(q.trap)
  )
  effects <- if (q.trap) pracma::cumtrapz(age, beta)
    else c(0, cumsum(diff(age) * zoo::rollmean(beta, 2)))
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

get_total_genetic_effect <- function(ages, betas, exponents, fixed_effect) {
  stopifnot(
    length(betas) == length(exponents),
    length(betas) > 0, length(exponents) > 0, length(ages) > 0
  )
  result <- vapply(ages, \(age) {
    sum(betas * age ^ exponents) + fixed_effect
  }, numeric(1))
  # Add 0 in beginning bc integral at time 0 is 0
  c(0,ages,0, result) |> matrix(ncol=2)
}

time_dependent_MR <- function(age.seq, model_type = c("loess", "glm"), loess_model, betas, exponents, fixed_effect,) {
  model_type <- match.arg(model_type)
  midpoints <-  zoo::rollmean(age.seq, 2)
  total_effects <- NULL
  midpoint_effects <- NULL
  if(model_type == "loess"){
    total_effects <- predict(loess_model, age.se, se = TRUE)
    midpoint_effects <- predict(loess_model, midpoints, se = TRUE)
  } else {
    total_effects <- get_total_genetic_effect(ages = age.seq, betas = betas,
                                              exponents = exponents,
                                              fixed_effect = fixed_effect)
    midpoint_effects <- get_total_genetic_effect(ages = midpoints, betas = betas,
                                              exponents = exponents,
                                              fixed_effect = fixed_effect)
  }
  wald_ratio <- dBeta / midpoint_effects
  Gamma <- integrate_estimates(midpoints, wald_ratio)
  Gamma_variance <-
    if (q.trap) calculate_mr_variance_trapz(effect = total_effects, variance = beta[, "var"])
    else beta[, "var"] / total_effects ^ 2
  data.frame(
    Gamma = Gamma,
    var = Gamma_variance,
    L95 = Gamma - 1.96 * sqrt(Gamma_variance),
    L95 = Gamma + 1.96 * sqrt(Gamma_variance)
  )

}
