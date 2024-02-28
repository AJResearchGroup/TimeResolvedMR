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

# calcualte central difference (midpoint) to estimate dbdt and gamma(t) = dbdt/beta_G(t), then use trapezoidal rule for integration to obtain the MR estimate Gamma(t):
dbdt <- matrix(rep(0, (length(age.seq) - 1) * 2), (length(age.seq) - 1), 2)
gmat <- matrix(rep(0, (length(age.seq) - 1) * 4), (length(age.seq) - 1), 4)

gmat[, 1] <- dbdt[, 1] <- zoo::rollmean(age.seq, 2)

dbdt[, 2] <- diff(beta[, "beta"])    # numerical derivative (central difference), calculated in midpoint t = k - 1/2

gmat[, 2] <- dbdt[, 2] /
  (coef.exp.lin["pgs"] + coef.exp.lin["pgs:age_at_assessment"] * gmat[, 1])   # time-resolved wald ratio at midpoint t = k - 1/2; linear time-dependence
gmat[, 3] <- dbdt[, 2] / (
  coef.exp.qdr["pgs"] + coef.exp.qdr["pgs:age_at_assessment"] * gmat[, 1] +
  coef.exp.qdr["pgs:I(age_at_assessment^4)"] * gmat[, 1] ^ 4   # quartic time-dependence
)

loess.f  <- predict(loess.m, gmat[, 1], se = TRUE)

gmat[, 4] <-
  dbdt[, 2] / loess.f$fit                                        # LOESS time-dependence (with linear extrapolation)




# integration of gamma(t) to Gamma(t):
q.trap <-
  FALSE   # q.trap = TRUE -> trapezoidal rule, q.trap = FALSE -> midpoint rule (Riemann sum)

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

Glin <- integrate_estimates(gmat[,1], gmat[,2], q.trap)
Gqdr <- integrate_estimates(gmat[,1], gmat[,3], q.trap)
Gloe <- integrate_estimates(gmat[,1], gmat[,4], q.trap)

loess.f   <- predict(loess.m, beta[, 1], se = TRUE)

total_effect_linear <- get_total_genetic_effect(betas = coef.exp.lin["pgs:age_at_assessment"], exponents = 1, fixed_effect = coef.exp.lin["pgs"])
total_effect_qdr <- get_total_genetic_effect(betas = coef.exp.lin[c("pgs:age_at_assessment","pgs:I(age_at_assessment^4)")], exponents = c(1,4), fixed_effect = coef.exp.lin["pgs"])
total_effect_loess <- loess.f$fit

if (q.trap) {
  Glin[, 3] <- calculate_mr_variance_trapz(effect = total_effect_linear, variance = beta[, "var"])
  Gqdr[, 3] <- calculate_mr_variance_trapz(effect = total_effect_qdr, variance = beta[, "var"])
  Gloe[, 3] <- calculate_mr_variance_trapz(effect = total_effect_loess, variance = beta[, "var"])
} else {
  # 95% CIs using the variance for the midpoint rule (-> uncorrelated errors; wider CIs than for the trapezoidal rule)
  Glin[, 3] <- beta[, "var"] / total_effect_linear ^ 2
  Gqdr[, 3] <- beta[, "var"] / total_effect_qdr ^ 2
  Gloe[, 3] <- beta[, "var"] / total_effect_loess ^ 2
}
Glin[, 4] <- Glin[, 2] - 1.96 * sqrt(Glin[, 3])
Glin[, 5] <- Glin[, 2] + 1.96 * sqrt(Glin[, 3])

Gqdr[, 4] <- Gqdr[, 2] - 1.96 * sqrt(Gqdr[, 3])
Gqdr[, 5] <- Gqdr[, 2] + 1.96 * sqrt(Gqdr[, 3])

Gloe[, 4] <- Gloe[, 2] - 1.96 * sqrt(Gloe[, 3])
Gloe[, 5] <- Gloe[, 2] + 1.96 * sqrt(Gloe[, 3])

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
  vapply(ages, \(age) {
    sum(betas * age ^ exponents) + fixed_effect
  }, numeric(1))
}

# beta <- beta
# Glinear <- Glin
# Gquartic <- Gqdr
# Gloess <- Gloe
