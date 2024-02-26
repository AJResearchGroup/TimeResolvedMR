#Parameters
#TODO replace with Roxygen
# pgs : Polygenic risk scores (hopefully a data frame)
# exp.data: some sort of exposure data. Measurements maybe?
# data: outcome data (see above)
# Merged PGS: if TRUE, compute merged PGS, otherwise subcohort (??)
# w.outtype: Either binary or continuous
# age.step: Step size for age sequence
perform_mr <- function(
    data, pgs, exposure, outcome, time_at_event, age_at_assessment,
    covariates = setdiff(colnames(data), c(exposure, outcome, time_at_event, age_at_assessment)),
    merged_pgs = TRUE, age.step = 1) {

  # Set up age bounds
  age.start <- 0
  age.stop <- 81 # identify maximum age here
  age.seq <- seq(rom = age.start, to = age.stop,  by = age.step)
  prior_event <-  # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
    data[, exposure] == 1 &
    (data[, time_at_event] < data[, age_at_assessment])

  data_clean <- na.omit(data)
  data_clean$outcome_controls <- !dis

  healthy <- which(data_clean$outcome_controls)

  # Include PGS, PGS-Age interaction and all covariates
  linear_formula <- reformulate(
    termlabels = c(
      pgs,
      stringr::str_glue("{pgs}:{age_at_assessment}"),
      covariates
    ),
    response = exposure
  )
  linear_model_exposure <- glm(
      linear_formula,
      family = "gaussian",
      subset = healthy,
      data = data_clean
  )

  # Same as above but include the age^4. I is necessary to prevent special
  # treatment of ^ in formula
  quadratic_formula <- reformulate(
    termlabels = c(
      pgs,
      stringr::str_glue("{pgs}:{age_at_assessment}"),
      stringr::str_glue("{pgs}:I({age_at_assessment}^4)"),
      covariates
    ),
    response = exposure
  )
  quadratic_model_exposure <- glm(
    quadratic_formula,
    family = "gaussian",
    subset = healthy,
    data = data_clean
  )

  coef.exp.lin <-
    summary(linear_model_exposure)$coef[c(2, 51), 1]      # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!
  coef.exp.qdr <-
    summary(quadratic_model_exposure)$coef[c(2, 51, 52), 1]   # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!

  # local regression (LOESS) fitting to age-stratified data:
  age.ukb <- seq(40, 70, by = 1)
  coef.exp.age <-
    matrix(rep(NA, length(age.ukb) * 3), length(age.ukb), 3)
  for (j in 1:length(age.ukb))
  {
    if (j == 1)
    {
      exp.data.age <-
        exp.data.omit[which(exp.data.omit$age_at_assessment <= age.ukb[j]), ]
    } else if (j == length(age.ukb))
    {
      exp.data.age <-
        exp.data.omit[which(exp.data.omit$age_at_assessment >= age.ukb[j]), ]
    } else
    {
      exp.data.age <-
        exp.data.omit[which(exp.data.omit$age_at_assessment == age.ukb[j]), ]
    }

    healthy <- which(exp.data.age$outcome_controls)
    p.exp.age <-
      glm(
        bmi.ztf ~ pgs + sex + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 +
          pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        family = "gaussian",
        subset = healthy,
        data = exp.data.age
      )

    if (dim(summary(p.exp.age)$coef)[1] > 1)
      coef.exp.age[j, ] <- summary(p.exp.age)$coef[2, c(1, 2, 4)]
  }

  # time-dependent (disease-specific) effects are estimated from a LOESS fit to age-stratified effects (NOTE: outcome cases are removed prior to effect estimation)
  loess.m <-
    loess(
      coef.exp.age[, 1] ~ age.ukb,
      control = loess.control(surface = "direct"),
      family = "gaussian",
      span = 0.50,
      degree = 1
    )   # linear extrapolation
  # loess.fit <- predict(loeess.m, t (insert specific variable name), se = TRUE)
  # dEta_dt/loess.fit$fit

  date()

  # PGS -> outcome: estimation of effect of PGS on outcome, using Aalen's additive hazard model
    data$event <- data[, 2]
    data$tstop <- data[, 3]

  dim(exp.data)   # log size of exp.data to check that it is compatible with data and original exp.data (= exp.data.trmr)
  covars <-
    cbind(
      exp.data$sex,
      data$array,
      data$centre,
      data$pc1,
      data$pc2,
      data$pc3,
      data$pc4,
      data$pc5,
      data$pc6,
      data$pc7,
      data$pc8,
      data$pc9,
      data$pc10,
      data$pc11,
      data$pc12,
      data$pc13,
      data$pc14,
      data$pc15,
      data$pc16,
      data$pc17,
      data$pc18,
      data$pc19,
      data$pc20,
      data$pc21,
      data$pc22,
      data$pc23,
      data$pc24,
      data$pc25
    )
  covars <- as.data.frame(covars)

  colnames(covars) <-
    c(
      "sex",
      "array",
      "centre",
      "pc1",
      "pc2",
      "pc3",
      "pc4",
      "pc5",
      "pc6",
      "pc7",
      "pc8",
      "pc9",
      "pc10",
      "pc11",
      "pc12",
      "pc13",
      "pc14",
      "pc15",
      "pc16",
      "pc17",
      "pc18",
      "pc19",
      "pc20",
      "pc21",
      "pc22",
      "pc23",
      "pc24",
      "pc25"
    )

    covars$event <- data$event
    covars$tstop <- outcome_measurement$tstop
  covars$pgs    <- pg
  covars$centre <- as.factor(covars$centre)
  covars        <- na.omit(covars)

  p.out <-
    aalen(
      Surv(tstop, event) ~ pgs + sex + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
      data = covars
    )

}
