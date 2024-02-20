###########################################################################################
# Script 5: Causal estimation of the time-dependent, life-course effect of BMI on outcome #
###########################################################################################

load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/BMI_DATA_TRMR/output_data_bmi_t2dm_aalen_steiger_15152.RData",
  verbose = TRUE
)

coef.exp.lin          <- as.numeric(coef.exp.lin)
coef.exp.qdr          <- as.numeric(coef.exp.qdr)

coef.exp.lin
coef.exp.qdr


bcum.raw <-
  cbind(p.out$cum[, c(1, 3)], p.out$var.cum[, 3])   # time, B(t) and Var(B(t))
colnames(bcum.raw) <- c("time", "pgs", "pgs.var")

# interpolate B(t) and Var(B(t)) for non-existing t (ages), i.e., times with no events:
w.last <- 1
bcum <- {
}
for (j in 1:length(age.seq))
{
  w.tmp <- which(trunc(bcum.raw[, 1]) == age.seq[j])
  n.tmp <- length(w.tmp)

  if (n.tmp > 0)
  {
    bcum <- rbind(bcum, bcum.raw[w.tmp, ])
    w.last <- w.tmp[length(w.tmp)]
  } else
  {
    insert.tmp <- c(age.seq[j], as.numeric(bcum.raw[w.last, c(2, 3)]))
    bcum <- rbind(bcum, insert.tmp)
  }
}
# for (j in 1:nrow(bcum$cum))   # if bcum$cum[,3] = NA, then assume that beta (=dB/dt) is zero and interpolate to the next non-NA time-point
# {
#     if (is.na(bcum$cum[j,3])) bcum$cum[j,3] <- bcum$cum[j-1,3]
# }


# calculate b(t) = dB/dt and corresponding standard errors:
dt       <-
  1.0   # for now, suppose that dt = 1 (this should be implemented in the future to allow for non-unitary dt)
beta     <-
  matrix(rep(0, length(age.seq) * 5), length(age.seq), 5)   # dB/dt
beta[, 1] <- age.seq
for (j in 2:length(age.seq))
{
  w.Bt1 <-
    max(which(trunc(bcum[, 1]) == age.seq[j]))     # identify maximum time-point for age = j
  w.Bt0 <-
    max(which(trunc(bcum[, 1]) == age.seq[j - 1]))   # identify maximum time-point for age = j-1

  beta[j, 2] <-
    bcum[w.Bt1, 2] - bcum[w.Bt0, 2]         # beta(t)dt = B(t) - B(t-), suppose that dt = 1
  beta[j, 3] <-
    bcum[w.Bt1, 3] - bcum[w.Bt0, 3]         # variance of beta(t); it appears that var.cum is more stable than robvar.cum...?
  #    beta[j,3] <- bcum$robvar.cum[w.Bt1,3] - bcum$robvar.cum[w.Bt0,3]   # variance of beta(t); robvar.cum (less stable)
  beta[j, 4] <- beta[j, 2] - 1.96 * sqrt(beta[j, 3])   # 95% lower CI
  beta[j, 5] <- beta[j, 2] + 1.96 * sqrt(beta[j, 3])   # 95% upper CI
}

# calcualte central difference (midpoint) to estimate dbdt and gamma(t) = dbdt/beta_G(t), then use trapezoidal rule for integration to obtain the MR estimate Gamma(t):
dbdt <- matrix(rep(0, (length(age.seq) - 1) * 2), (length(age.seq) - 1), 2)
gmat <- matrix(rep(0, (length(age.seq) - 1) * 4), (length(age.seq) - 1), 4)

dbdt[, 1] <-
  age.seq[2:length(age.seq)] - 0.5   # midpoint time t = k - 1/2, db/dt(t) is estimated in the interior of [0,age_max]
gmat[, 1] <-
  dbdt[, 1]                           # midpoint time t = k - 1/2, gamma(t) is estimated in the interior of [0,age_max]

dbdt[, 2] <-
  beta[2:length(age.seq), 2] - beta[1:(length(age.seq) - 1), 2]   # numerical derivative (central difference), calculated in midpoint t = k - 1/2

gmat[, 2] <-
  dbdt[, 2] / (coef.exp.lin[1] + coef.exp.lin[2] * gmat[, 1])   # time-resolved wald ratio at midpoint t = k - 1/2; linear time-dependence
gmat[, 3] <-
  dbdt[, 2] / (coef.exp.qdr[1] + coef.exp.qdr[2] * gmat[, 1] + coef.exp.qdr[3] *
                 gmat[, 1] ^ 4)   # quartic time-dependence

loess.f  <- predict(loess.m, gmat[, 1], se = TRUE)

gmat[, 4] <-
  dbdt[, 2] / loess.f$fit                                        # LOESS time-dependence (with linear extrapolation)




# integration of gamma(t) to Gamma(t):
q.trap <-
  FALSE   # q.trap = TRUE -> trapezoidal rule, q.trap = FALSE -> midpoint rule (Riemann sum)

Glin <- matrix(rep(0, length(age.seq) * 5), length(age.seq), 5)
Gqdr <- matrix(rep(0, length(age.seq) * 5), length(age.seq), 5)
Gloe <- matrix(rep(0, length(age.seq) * 5), length(age.seq), 5)

if (q.trap)
  # trapezoidal rule
{
  Glin[1, 1] <- 0
  Glin[2:length(age.seq), 1] <- age.seq[2:length(age.seq)] - 0.5
  Gqdr[1, 1] <- 0
  Gqdr[2:length(age.seq), 1] <- age.seq[2:length(age.seq)] - 0.5
  Gloe[1, 1] <- 0
  Gloe[2:length(age.seq), 1] <- age.seq[2:length(age.seq)] - 0.5

  Glin[1, 2] <-
    0                                                               # integral at t=0 is zero
  Glin[2, 2] <-
    0.5 * gmat[1, 2] / 2                                                 # approximation of integral at the first time-step t=1/2 (= half the default time-step)
  Gqdr[1, 2] <-
    0                                                               # integral at t=0 is zero
  Gqdr[2, 2] <-
    0.5 * gmat[1, 3] / 2                                                 # approximation of integral at the first time-step t=1/2 (= half the default time-step)
  Gloe[1, 2] <-
    0                                                               # integral at t=0 is zero
  Gloe[2, 2] <-
    0.5 * gmat[1, 4] / 2                                                 # approximation of integral at the first time-step t=1/2 (= half the default time-step)
  for (j in 2:nrow(gmat))
  {
    if (j == 2)
    {
      Glin[3, 2] <-
        Glin[2, 2] + 0.5 * (gmat[1, 2] + gmat[2, 2])   # approximation of integral at the second time-step (no separate calc. needed, but anyway)
      Gqdr[3, 2] <-
        Gqdr[2, 2] + 0.5 * (gmat[1, 3] + gmat[2, 3])   # approximation of integral at the second time-step (no separate calc. needed, but anyway)
      Gloe[3, 2] <-
        Gloe[2, 2] + 0.5 * (gmat[1, 4] + gmat[2, 4])   # approximation of integral at the second time-step (no separate calc. needed, but anyway)
    }
    if (j  > 2)
      # approximation of integral at remaining time-steps, using the trapezoidal rule
    {
      # Gamma(t) for linear time-dependence in beta_G(t):
      int.tmp <- Glin[2, 2] + 0.5 * gmat[1, 2]
      for (k in 2:(j - 1))
      {
        int.tmp <- int.tmp + gmat[k, 2]
      }
      int.tmp <- int.tmp + 0.5 * gmat[j, 2]
      Glin[j + 1, 2] <- int.tmp

      # Gamma(t) for quartic time-dependence in beta_G(t):
      int.tmp <- Gqdr[2, 2] + 0.5 * gmat[1, 3]
      for (k in 2:(j - 1))
      {
        int.tmp <- int.tmp + gmat[k, 3]
      }
      int.tmp <- int.tmp + 0.5 * gmat[j, 3]
      Gqdr[j + 1, 2] <- int.tmp

      # Gamma(t) for nonlinear (LOESS) time-dependence in beta_G(t):
      int.tmp <- Gloe[2, 2] + 0.5 * gmat[1, 4]
      for (k in 2:(j - 1))
      {
        int.tmp <- int.tmp + gmat[k, 4]
      }
      int.tmp <- int.tmp + 0.5 * gmat[j, 4]
      Gloe[j + 1, 2] <- int.tmp
    }
  }
} else
  # midpoint rule
{
  Glin[, 1] <- age.seq
  Gqdr[, 1] <- age.seq
  Gloe[, 1] <- age.seq

  Glin[1, 2] <-
    0                                                               # integral at t=0 is zero
  Glin[2:nrow(Glin), 2] <-
    cumsum(gmat[, 2])                                     # approximation of integral using the midpoint rule
  Gqdr[1, 2] <-
    0                                                               # integral at t=0 is zero
  Gqdr[2:nrow(Gqdr), 2] <-
    cumsum(gmat[, 3])                                     # approximation of integral using the midpoint rule
  Gloe[1, 2] <-
    0                                                               # integral at t=0 is zero
  Gloe[2:nrow(Gloe), 2] <-
    cumsum(gmat[, 4])                                     # approximation of integral using the midpoint rule
}



# finally, calculate 95% CIs, corresponding to each effect estimate Gamma(t) using the expression for Var(Gamma(t)):
Glin[1, 3] <- 0
Glin[1, 4] <- 0
Glin[1, 5] <- 0
Gqdr[1, 3] <- 0
Gqdr[1, 4] <- 0
Gqdr[1, 5] <- 0
Gloe[1, 3] <- 0
Gloe[1, 4] <- 0
Gloe[1, 5] <- 0
Gprx <- Gqdr
loess.f   <- predict(loess.m, beta[, 1], se = TRUE)

if (q.trap)
  # 95% CIs using the variance for the trapezoiudal-rule integration (-> correlated errors)
{
  for (j in 2:nrow(Glin))
  {
    if (j == 2)
      # derived MR variance:
    {
      Glin[j, 3] <-
        beta[2, 3] / (coef.exp.lin[1] + coef.exp.lin[2] * beta[2, 1]) ^ 2 / 16        # NOTE: variance in Gamma(k-1/2) is estimated from Var[Eta(k)]/beta_G(k)^2/16
      Gqdr[j, 3] <-
        beta[2, 3] / (coef.exp.qdr[1] + coef.exp.qdr[2] * beta[2, 1] + coef.exp.qdr[3] *
                        beta[2, 1] ^ 4) ^ 2 / 16
      Gloe[j, 3] <- beta[2, 3] / (loess.f$fit[2]) ^ 2 / 16

      Glin[j, 4] <- Glin[j, 2] - 1.96 * sqrt(Glin[j, 3])
      Glin[j, 5] <- Glin[j, 2] + 1.96 * sqrt(Glin[j, 3])

      Gqdr[j, 4] <- Gqdr[j, 2] - 1.96 * sqrt(Gqdr[j, 3])
      Gqdr[j, 5] <- Gqdr[j, 2] + 1.96 * sqrt(Gqdr[j, 3])

      Gloe[j, 4] <- Gloe[j, 2] - 1.96 * sqrt(Gloe[j, 3])
      Gloe[j, 5] <- Gloe[j, 2] + 1.96 * sqrt(Gloe[j, 3])

      Gvar <-
        beta[j, 3] / (2 * (
          coef.exp.qdr[1] + coef.exp.qdr[2] * beta[j, 1] + coef.exp.qdr[3] * beta[j, 1] ^
            4
        ) ^ 2)
      Gprx[j, 3] <- Gvar
      Gprx[j, 4] <- Gprx[j, 2] - 1.96 * sqrt(Gvar)
      Gprx[j, 5] <- Gprx[j, 2] + 1.96 * sqrt(Gvar)
    } else if (j == 3)
    {
      Glin[j, 3] <-
        Glin[2, 3] + beta[3, 3] / (coef.exp.lin[1] + coef.exp.lin[2] * beta[3, 1]) ^
        2 / 4
      Gqdr[j, 3] <-
        Gqdr[2, 3] + beta[3, 3] / (coef.exp.qdr[1] + coef.exp.qdr[2] * beta[3, 1] + coef.exp.qdr[3] *
                                     beta[3, 1] ^ 4) ^ 2 / 4
      Gloe[j, 3] <- Gloe[2, 3] + beta[3, 3] / (loess.f$fit[3]) ^ 2 / 4

      Glin[j, 4] <- Glin[j, 2] - 1.96 * sqrt(Glin[j, 3])
      Glin[j, 5] <- Glin[j, 2] + 1.96 * sqrt(Glin[j, 3])

      Gqdr[j, 4] <- Gqdr[j, 2] - 1.96 * sqrt(Gqdr[j, 3])
      Gqdr[j, 5] <- Gqdr[j, 2] + 1.96 * sqrt(Gqdr[j, 3])

      Gloe[j, 4] <- Gloe[j, 2] - 1.96 * sqrt(Gloe[j, 3])
      Gloe[j, 5] <- Gloe[j, 2] + 1.96 * sqrt(Gloe[j, 3])

      Gvar <-
        beta[j, 3] / (2 * (
          coef.exp.qdr[1] + coef.exp.qdr[2] * beta[j, 1] + coef.exp.qdr[3] * beta[j, 1] ^
            4
        ) ^ 2)
      Gprx[j, 3] <- Gvar
      Gprx[j, 4] <- Gprx[j, 2] - 1.96 * sqrt(Gvar)
      Gprx[j, 5] <- Gprx[j, 2] + 1.96 * sqrt(Gvar)
    } else
    {
      Glin[j, 3] <-
        Glin[2, 3] + beta[j - 1, 3] / (coef.exp.lin[1] + coef.exp.lin[2] * beta[j -
                                                                                  1, 1]) ^ 2 / 4 + beta[j, 3] / (coef.exp.lin[1] + coef.exp.lin[2] * beta[j, 1]) ^
        2 / 4
      Gqdr[j, 3] <-
        Gqdr[2, 3] + beta[j - 1, 3] / (coef.exp.qdr[1] + coef.exp.qdr[2] * beta[j -
                                                                                  1, 1] + coef.exp.qdr[3] * beta[j - 1, 1] ^ 4) ^ 2 / 4 + beta[j, 3] / (coef.exp.qdr[1] + coef.exp.qdr[2] *
                                                                                                                                                          beta[j, 1] + coef.exp.qdr[3] * beta[j, 1] ^ 4) ^ 2 / 4
      Gloe[j, 3] <-
        Gloe[2, 3] + beta[j - 1, 3] / (loess.f$fit[j - 1]) ^ 2 / 4 + beta[j, 3] /
        (loess.f$fit[j]) ^ 2 / 4

      Glin[j, 4] <- Glin[j, 2] - 1.96 * sqrt(Glin[j, 3])
      Glin[j, 5] <- Glin[j, 2] + 1.96 * sqrt(Glin[j, 3])

      Gqdr[j, 4] <- Gqdr[j, 2] - 1.96 * sqrt(Gqdr[j, 3])
      Gqdr[j, 5] <- Gqdr[j, 2] + 1.96 * sqrt(Gqdr[j, 3])

      Gloe[j, 4] <- Gloe[j, 2] - 1.96 * sqrt(Gloe[j, 3])
      Gloe[j, 5] <- Gloe[j, 2] + 1.96 * sqrt(Gloe[j, 3])

      Gvar <-
        beta[j, 3] / (2 * (
          coef.exp.qdr[1] + coef.exp.qdr[2] * beta[j, 1] + coef.exp.qdr[3] * beta[j, 1] ^
            4
        ) ^ 2)
      Gprx[j, 3] <- Gvar
      Gprx[j, 4] <- Gprx[j, 2] - 1.96 * sqrt(Gvar)
      Gprx[j, 5] <- Gprx[j, 2] + 1.96 * sqrt(Gvar)
    }
  }
} else
  # 95% CIs using the variance for the midpoint rule (-> uncorrelated errors; wider CIs than for the trapezoidal rule)
{
  for (j in 2:nrow(Glin))
  {
    # these variances must be verified!!!
    Glin[j, 3] <-
      beta[j, 3] / (coef.exp.lin[1] + coef.exp.lin[2] * beta[j, 1]) ^ 2
    Gqdr[j, 3] <-
      beta[j, 3] / (coef.exp.qdr[1] + coef.exp.qdr[2] * beta[j, 1] + coef.exp.qdr[3] *
                      beta[j, 1] ^ 4) ^ 2
    Gloe[j, 3] <- beta[j, 3] / (loess.f$fit[j]) ^ 2

    Glin[j, 4] <- Glin[j, 2] - 1.96 * sqrt(Glin[j, 3])
    Glin[j, 5] <- Glin[j, 2] + 1.96 * sqrt(Glin[j, 3])

    Gqdr[j, 4] <- Gqdr[j, 2] - 1.96 * sqrt(Gqdr[j, 3])
    Gqdr[j, 5] <- Gqdr[j, 2] + 1.96 * sqrt(Gqdr[j, 3])

    Gloe[j, 4] <- Gloe[j, 2] - 1.96 * sqrt(Gloe[j, 3])
    Gloe[j, 5] <- Gloe[j, 2] + 1.96 * sqrt(Gloe[j, 3])
  }
}

# beta <- beta
# Glinear <- Glin
# Gquartic <- Gqdr
# Gloess <- Gloe
