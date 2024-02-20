##################################################################################################
# Script 4: Estimation of time-dependent betaG(t) and Eta(t) (= effect of instrument on outcome) #
##################################################################################################

# R_code_for_trmr_aalen_xxxx_15152.R   // generic name, xxxx is replaced by disease-outcome name
# Created by Torgny Karlsson on 2022-11-21

# T2DM, all individuals:

date()

load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_t2dm.RData",
  verbose = TRUE
)   # T2DM: stanardized PGS (IV)
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/data_for_trmr.RData",
  verbose = TRUE
)                          # exposure: BMI
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_t2dm.RData",
  verbose = TRUE
)             # outcome:  disease, please change!!!

library(survival)
library(sandwich)
library(timereg)


q.aalen      <- TRUE    # Aalen's additive hazard model
q.merged.pgs <-
  TRUE    # q.merged.pgs = TRUE -> estimate effect on outcome of merged PGS; q.merged.pgs = FALSE -> estimate effect of PGS in each subsample, separately
w.outtype    <-
  "bin"   # w.outtype = "bin" (binary) -> disease; w.outtype = "cnt" (continuous) -> continuous response

# first of all, rename data-frames (NOTE: disease specific, change for different outcomes):
pgs      <-
  pgs.dm           # PGS <--- NOTE: disease-specific, please change!!!
exp.data <- exp.data.trmr    # exposure data-frame
out.data <-
  t2dm.data.trmr   # disease outcome data-frame, please change!!!


# some initializations:
age.start <-  0
age.stop  <- 81
age.step  <-  1

set.oa  <- c(1:32)
set.ps  <- c(1:32)
set.ad  <- c(1:32)
set.gc  <- c(1:4, 8:35)
set.rc  <- c(1, 5:35)
set.uc  <- c(1:32)
set.af  <-
  c(1, 48:50, 5:32)   # order is important for downstream analyses (event is set to out.data[,2], etc.)
set.cad <- c(1:32)
set.t2d <- c(1:32)
set.dth <- c(1:32)

set.dis <-
  set.t2d             # <--- define outcome, please change!!!

if (q.aalen)
{
  if (q.merged.pgs)
  {
    age.seq          <- seq(age.start, age.stop, age.step)

    # PGS -> exposure: estimation of effect of PGS on exposure BMI
    out.data <- out.data[, set.dis]
    if (w.outtype == "bin")
    {
      # dis      <- which(out.data[,2] == 1)           # <--- NOTE: remove outcome cases for exposure to minimize reversed causation/weak-instrument bias for time-fixed betaG
      dis      <-
        which((out.data[, 2] == 1) &
                (out.data[, 3] < exp.data$age_at_assessment)) # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
      # dis      <- {}   # removed none, assuming that the binary response is just a relatization of an underlying continuous response
    } else
    {
      dis      <-
        {
        }                                 # <--- NOTE: for continuous response, no removal of cases to minimize weak-instrument bias can be performed
    }
    non.dis  <-
      !(seq(1, nrow(exp.data), by = 1) %in% dis)   # healthy controls

    exp.data$pgs <- pgs
    exp.data$outcome_controls <- non.dis
    exp.data.omit <- na.omit(exp.data)

    healthy <- which(exp.data.omit$outcome_controls)
    p.exp.lin <-
      glm(
        bmi.ztf ~ pgs + pgs:age_at_assessment + sex + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit
      )
    p.exp.qdr <-
      glm(
        bmi.ztf ~ pgs + pgs:age_at_assessment + pgs:I(age_at_assessment ^ 4) + sex + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit
      )

    coef.exp.lin <-
      summary(p.exp.lin)$coef[c(2, 51), 1]      # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!
    coef.exp.qdr <-
      summary(p.exp.qdr)$coef[c(2, 51, 52), 1]   # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!

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
    if (w.outtype == "bin")
    {
      out.data$event <- out.data[, 2]
      out.data$tstop <- out.data[, 3]
    }

    dim(exp.data)   # log size of exp.data to check that it is compatible with out.data and original exp.data (= exp.data.trmr)
    covars <-
      cbind(
        exp.data$sex,
        out.data$array,
        out.data$centre,
        out.data$pc1,
        out.data$pc2,
        out.data$pc3,
        out.data$pc4,
        out.data$pc5,
        out.data$pc6,
        out.data$pc7,
        out.data$pc8,
        out.data$pc9,
        out.data$pc10,
        out.data$pc11,
        out.data$pc12,
        out.data$pc13,
        out.data$pc14,
        out.data$pc15,
        out.data$pc16,
        out.data$pc17,
        out.data$pc18,
        out.data$pc19,
        out.data$pc20,
        out.data$pc21,
        out.data$pc22,
        out.data$pc23,
        out.data$pc24,
        out.data$pc25
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

    if (w.outtype == "bin")
    {
      covars$event <- out.data$event
      covars$tstop <- out.data$tstop
    } else
    {
      #            covars$out <- out.data$response
    }
    covars$pgs    <- pgs
    covars$centre <- as.factor(covars$centre)
    covars        <- na.omit(covars)

    p.out <-
      aalen(
        Surv(tstop, event) ~ pgs + sex + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        data = covars
      )

  } else
  {
    # first subsample:
    # second subsample:
    p.out <- 0

  }
}


date()

save(
  age.seq,
  age.ukb,
  coef.exp.lin,
  coef.exp.qdr,
  coef.exp.age,
  p.exp.lin,
  p.exp.qdr,
  p.exp.age,
  loess.m,
  p.out,
  q.aalen,
  q.merged.pgs,
  w.outtype,
  set.dis,
  exp.data,
  out.data,
  covars,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/BMI_DATA_TRMR/output_data_bmi_t2dm_aalen_steiger_15152.RData"
)

date()



# T2DM, females:
# R_code_for_trmr_aalen_xxxx_15152.R   // generic name, xxxx is replaced by disease-outcome name
# Created by Torgny Karlsson on 2022-11-21

date()

load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_t2dm.RData",
  verbose = TRUE
)   # T2DM: stanardized PGS (IV)
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/data_for_trmr.RData",
  verbose = TRUE
)                          # exposure: BMI
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_t2dm.RData",
  verbose = TRUE
)             # outcome:  disease, please change!!!

library(survival)
library(sandwich)
library(timereg)


q.aalen      <- TRUE    # Aalen's additive hazard model
q.merged.pgs <-
  TRUE    # q.merged.pgs = TRUE -> estimate effect on outcome of merged PGS; q.merged.pgs = FALSE -> estimate effect of PGS in each subsample, separately
w.outtype    <-
  "bin"   # w.outtype = "bin" (binary) -> disease; w.outtype = "cnt" (continuous) -> continuous response
w.sex        <- 0       # (females,males) = (0,1)

# first of all, rename data-frames (NOTE: disease specific, change for different outcomes):
pgs      <-
  pgs.dm           # PGS <--- NOTE: disease-specific, please change!!!
exp.data <- exp.data.trmr    # exposure data-frame
out.data <-
  t2dm.data.trmr   # disease outcome data-frame, please change!!!


# some initializations:
age.start <-  0
age.stop  <- 81
age.step  <-  1

set.oa  <- c(1:32)
set.ps  <- c(1:32)
set.ad  <- c(1:32)
set.gc  <- c(1:4, 8:35)
set.rc  <- c(1, 5:35)
set.uc  <- c(1:32)
set.cad <- c(1:32)
set.t2d <- c(1:32)
set.dth <- c(1:32)

set.dis <-
  set.t2d             # <--- define outcome, please change!!!

if (q.aalen)
{
  if (q.merged.pgs)
  {
    age.seq          <- seq(age.start, age.stop, age.step)

    # PGS -> exposure: estimation of effect of PGS on exposure BMI
    out.data <- out.data[, set.dis]
    if (w.outtype == "bin")
    {
      # dis      <- which(out.data[,2] == 1)   # <--- NOTE: remove outcome cases for exposure to minimize reversed causation/weak-instrument bias for time-fixed betaG
      dis      <-
        which((out.data[, 2] == 1) &
                (out.data[, 3] < exp.data$age_at_assessment)) # remove cases with T_event < T_assessment for estimation of exposure effect betaG = betaG(t)
      # dis      <- {}   # removed none, assuming that the binary response is just a relatization of an underlying continuous response
    } else
    {
      dis      <-
        {
        }                                 # <--- NOTE: for continuous response, no removal of cases to minimize weak-instrument bias can be performed
    }
    non.dis  <-
      !(seq(1, nrow(exp.data), by = 1) %in% dis)   # healthy controls

    exp.data$pgs <- pgs
    exp.data$outcome_controls <- non.dis
    exp.data.omit <- na.omit(exp.data)

    healthy <-
      which((exp.data.omit$outcome_controls) &
              (exp.data.omit$sex == w.sex))
    p.exp.lin <-
      glm(
        bmi.ztf ~ pgs + pgs:age_at_assessment + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit
      )
    p.exp.qdr <-
      glm(
        bmi.ztf ~ pgs + pgs:age_at_assessment + pgs:I(age_at_assessment ^ 4) + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit
      )

    coef.exp.lin <-
      summary(p.exp.lin)$coef[c(2, 50), 1]      # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!
    coef.exp.qdr <-
      summary(p.exp.qdr)$coef[c(2, 50, 51), 1]   # this is the time-dependent PGS effect on BMI to be used in the MR estimation!!!

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

      healthy <-
        which((exp.data.age$outcome_controls) &
                (exp.data.age$sex == w.sex))
      p.exp.age <-
        glm(
          bmi.ztf ~ pgs + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
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
    if (w.outtype == "bin")
    {
      out.data$event <- out.data[, 2]
      out.data$tstop <- out.data[, 3]
    }

    dim(exp.data)   # log size of exp.data to check that it is compatible with out.data and original exp.data (= exp.data.trmr)
    covars <-
      cbind(
        exp.data$sex,
        out.data$array,
        out.data$centre,
        out.data$pc1,
        out.data$pc2,
        out.data$pc3,
        out.data$pc4,
        out.data$pc5,
        out.data$pc6,
        out.data$pc7,
        out.data$pc8,
        out.data$pc9,
        out.data$pc10,
        out.data$pc11,
        out.data$pc12,
        out.data$pc13,
        out.data$pc14,
        out.data$pc15,
        out.data$pc16,
        out.data$pc17,
        out.data$pc18,
        out.data$pc19,
        out.data$pc20,
        out.data$pc21,
        out.data$pc22,
        out.data$pc23,
        out.data$pc24,
        out.data$pc25
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

    if (w.outtype == "bin")
    {
      covars$event <- out.data$event
      covars$tstop <- out.data$tstop
    } else
    {
      #            covars$out <- out.data$response
    }
    covars$pgs    <- pgs
    covars$centre <- as.factor(covars$centre)
    covars        <- na.omit(covars)

    p.out <-
      aalen(
        Surv(tstop, event) ~ pgs + array + centre + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25,
        data = subset(covars, sex == w.sex)
      )

  } else
  {
    # first subsample:
    # second subsample:

    p.out <- 0

  }
}

date()

save(
  age.seq,
  age.ukb,
  coef.exp.lin,
  coef.exp.qdr,
  coef.exp.age,
  p.exp.lin,
  p.exp.qdr,
  p.exp.age,
  loess.m,
  p.out,
  q.aalen,
  q.merged.pgs,
  w.outtype,
  set.dis,
  exp.data,
  out.data,
  covars,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/BMI_DATA_TRMR/output_data_bmi_t2dm_aalen_steiger_females_15152.RData"
)
