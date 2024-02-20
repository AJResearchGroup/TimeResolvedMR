############################################
# Script 1: Steiger filtering, subsample 1 #
############################################
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/data_for_trmr.RData",
  verbose = TRUE
)   # needed for w.half.1 and w.half.2
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/GENETIC_INSTR/BMI/genetic_instr_gts_bmi_12_no_filtering.Rata",
  verbose = TRUE
)

# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_all_cause_mortality.RData", verbose=TRUE)   # all-cause mortality
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_osteoarthritis.RData", verbose=TRUE)        # OA
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_t2dm.RData", verbose=TRUE)                  # T2DM
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cancers.RData", verbose=TRUE)               # Gastric, Kidney, Esophagus, and Colorectal cancers
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_esophagus.RData", verbose=TRUE)             # Esphagus, Esophagus adeno, Esophagus squamous cancer
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_ulcerative_colitis.RData", verbose=TRUE)    # UC
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_psoriasis.RData", verbose=TRUE)             # Psoriasis
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cad.RData", verbose=TRUE)                   # CAD + Atrial fibrillation
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cad_2019.RData", verbose=TRUE)              # CAD + Atrial fibrillation
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_stroke.RData", verbose=TRUE)                # Ischaemic stroke (from algorithmically-defined outcomes)
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_alzd.RData", verbose=TRUE)                  # Alzheimer's disease
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_multiple_sclerosis.RData", verbose=TRUE)    # Multiple sclerosis
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_myocardial_infarction.RData", verbose=TRUE) # Myocardial infarction
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_angina_pectoris.RData", verbose=TRUE)       # Angina pectoris
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_female_cancers.RData",
  verbose = TRUE
)          # BC, EC, or OC (female cancers)


library(sandwich)

gts   <- gts.1
snpid <- cbind(colnames(gts)[7:ncol(gts)])


# output (remember to change):

# out.data <-  dth.data.trmr[,c(1:32)]                   # all-cause mortality
# out.data <- arth.data.trmr[,c(1:32)]                   # OA
# out.data <- t2dm.data.trmr[,c(1:32)]                   # T2DM
# out.data <-  cnr.data.trmr[,c(1,2:4,14:41)]            # Gastric cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,5:7,14:41)]            # Kidney cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,8:10,14:41)]           # Esophagus cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,11:13,14:41)]          # Colorectal cancer <- double-check to extract correct cancer type!
# out.data <-  esc.data.trmr[,c(1,2:4,11:38)]            # Esophagus cancer (combined; same as above)
# out.data <-  esc.data.trmr[,c(1,5:7,11:38)]            # Esophagus adenocarcinoma
# out.data <-  esc.data.trmr[,c(1,8:10,11:38)]           # Esophagus squamous cell carcinoma
# out.data <- ulcc.data.trmr[,c(1:32)]                   # Ulcerative colitis
# out.data <- psor.data.trmr[,c(1:32)]                   # Psoriasis
# out.data <-  cad.data.trmr[,c(1:32)]                   # CAD
# out.data <-  strk.data.trmr[,c(1:32)]                  # Ischaemic stroke
# out.data <-  cad.data.trmr[,c(1,48:50,5:32)]           # Atrial fibrillation
# out.data <- myoc.data.trmr[,c(1:32)]                   # Myocardial infarction
# out.data <- apec.data.trmr[,c(1:32)]                   # Angina pectoris
# out.data <- alzd.data.trmr[,c(1:32)]                   # Alzheimer's diseasse
# out.data <- mscl.data.trmr[,c(1:32)]                   # Multiple sclerosis
out.data <-
  cnr.data.trmr[, c(1, 2:4, 11:38)]             # Breast cancer
# out.data <- cnr.data.trmr[,c(1,5:7,11:38)]             # Endometrial cancer
# out.data <- cnr.data.trmr[,c(1,8:10,11:38)]            # Ovarian cancer

colnames(out.data)

snp.out <- matrix(rep(0, nrow(snpid) * 7), nrow(snpid), 7)

out.data$event <- out.data[, 2]
out.data$tstop <- out.data[, 3]

out.data$subset <-
  (seq(1, nrow(gts), by = 1) %in% w.half.1)   # same subsample x in which the variants of gts.x were identified
length(which(out.data$subset))
for (j in 7:ncol(gts))
{
  k <- j - 6
  out.data$gts  <- gts[, j]             # SNP k = j - 6
  out.data.omit <- na.omit(out.data)   # remove all "NA"s

  healthy <- which(out.data.omit$subset)

  # remember to change event.xx and tstop.xx in both p.full, p.reduced, and p.corr:

  p.full    <-
    glm(
      event ~ gts + sex * poly(age_at_eof, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )
  p.reduced <-
    glm(
      event ~ sex * poly(age_at_eof, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )
  p.corr    <-
    glm(
      event ~ gts + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )

  p.model <- p.full
  snp.out[k, 1] <- summary(p.model)$coef[2, 1]
  cov.p <-
    vcovHC(p.model, type = "HC0")   # accounting for heteroscedasticity using White's estimator
  std.err <- sqrt(diag(cov.p))
  snp.out[k, 2] <- std.err[2]
  snp.out[k, 3] <-
    2 * pnorm(abs(snp.out[k, 1] / snp.out[k, 2]), 0, 1, lower.tail = FALSE)
  snp.out[k, 4] <-
    (summary(p.reduced)$deviance - summary(p.full)$deviance) / summary(p.full)$null.deviance   # squared semipartial correlation
  snp.out[k, 5] <-
    1 - summary(p.corr)$deviance / summary(p.corr)$null.deviance                               # squared correlation
  snp.out[k, 6] <-
    0.5 * log((1 + sqrt(snp.out[k, 4])) / (1 - sqrt(snp.out[k, 4])))                                   # Fisher's z-transformation
  snp.out[k, 7] <-
    length(healthy)                                                                          # sample size
  n.individuals <- length(p.corr$residuals)
}
snp.out[1:10, ]
snp.out.bc.rsq.1 <- snp.out   # change name!



# exposure:
exp.data <- exp.data.trmr
exp.data$tstop <-
  out.data$tstop   # this extra column is there in order to enable estimation of the exposure effect in the correct sample (e.g., female-specific cancers)

snp.exp <- matrix(rep(0, nrow(snpid) * 7), nrow(snpid), 7)

exp.data$subset <-
  (seq(1, nrow(gts), by = 1) %in% w.half.1)   # NOTE: estimated in the same (selection) sample, to avoid bias by selecting SNPs in the target sample
for (j in 7:ncol(gts))
{
  k <- j - 6
  exp.data$gts  <- gts[, j]             # SNP k = j - 6
  exp.data.omit <- na.omit(exp.data)   # remove all "NA"s

  healthy <- which(exp.data.omit$subset)

  p.full <-
    glm(
      bmi.ztf ~ gts + sex * poly(age_at_assessment, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre,
      family = "gaussian",
      subset = healthy,
      data = exp.data.omit
    )
  p.reduced <-
    glm(
      bmi.ztf ~ sex * poly(age_at_assessment, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre,
      family = "gaussian",
      subset = healthy,
      data = exp.data.omit
    )
  p.corr <-
    glm(bmi.ztf ~ gts,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit)

  p.model <- p.full
  snp.exp[k, 1:3] <-
    summary(p.model)$coef[2, c(1, 2, 4)]    # main effect
  snp.exp[k, 4]   <-
    (summary(p.reduced)$deviance - summary(p.full)$deviance) / summary(p.full)$null.deviance   # squared semipartial correlation
  snp.exp[k, 5]   <-
    1 - summary(p.corr)$deviance / summary(p.corr)$null.deviance                               # squared correlation
  snp.exp[k, 6]   <-
    0.5 * log((1 + sqrt(snp.exp[k, 4])) / (1 - sqrt(snp.exp[k, 4])))                                   # Fisher's z-transformation
  snp.exp[k, 7]   <-
    length(healthy)                                                                          # sample size

}
snp.exp[1:10, ]
snp.exp.bc.rsq.1 <- snp.exp   # change name!

# test of equality beetween correlated correlations:
rho.xy <- median(sqrt(snp.out[, 4] / snp.exp[, 4]))
rm2 <- 0.5 * (snp.exp[, 4] + snp.out[, 4])
f   <- 0.5 * (1 - rho.xy) / (1 - rm2)
h   <- (1 - f * rm2) / (1 - rm2)
zstat.spc <-
  (snp.exp[, 6] - snp.out[, 6]) * sqrt(snp.out[, 7] - 3) / sqrt(2 * h * (1 - rho.xy))

rho.xy <- median(sqrt(snp.out[, 5] / snp.exp[, 5]))
rm2 <- 0.5 * (snp.exp[, 5] + snp.out[, 5])
f   <- 0.5 * (1 - rho.xy) / (1 - rm2)
h   <- (1 - f * rm2) / (1 - rm2)
zstat.cor <-
  (atanh(sqrt(snp.exp[, 5])) - atanh(sqrt(snp.out[, 5]))) * sqrt(snp.out[, 7] -
                                                                   3) / sqrt(2 * h * (1 - rho.xy))

pval.spc <- pnorm(zstat.spc)
pval.cor <- pnorm(zstat.cor)

fct.a <- sqrt(2 * h * (1 - rho.xy))

psi   <-
  rho.xy * (1 - snp.exp[, 5] - snp.out[, 5]) - 0.5 * sqrt(snp.out[, 5] *
                                                            snp.exp[, 5]) * (1 - snp.exp[, 5] - snp.out[, 5] - rho.xy ^ 2)
chat  <- psi / ((1 - snp.exp[, 5]) * (1 - snp.out[, 5]))
fct.b <- sqrt(2 - 2 * chat)

save(
  snp.exp.bc.rsq.1,
  snp.out.bc.rsq.1,
  zstat.spc,
  zstat.cor,
  pval.spc,
  pval.cor,
  fct.a,
  fct.b,
  n.individuals,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_bc_1.RData"
)   # change name!


############################################
# Script 1: Steiger filtering, subsample 2 #
############################################
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/data_for_trmr.RData",
  verbose = TRUE
)   # needed for w.half.1 and w.half.2
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/GENETIC_INSTR/BMI/genetic_instr_gts_bmi_12_no_filtering.Rata",
  verbose = TRUE
)

# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_all_cause_mortality.RData", verbose=TRUE)   # all-cause mortality
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_osteoarthritis.RData", verbose=TRUE)        # OA
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_t2dm.RData", verbose=TRUE)                  # T2DM
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cancers.RData", verbose=TRUE)               # Gastric, Kidney, Esophagus, and Colorectal cancers
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_esophagus.RData", verbose=TRUE)             # Esphagus, Esophagus adeno, Esophagus squamous cancer
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_ulcerative_colitis.RData", verbose=TRUE)    # UC
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_psoriasis.RData", verbose=TRUE)             # Psoriasis
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cad.RData", verbose=TRUE)                   # CAD + Atrial fibrillation
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_cad_2019.RData", verbose=TRUE)              # CAD + Atrial fibrillation
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_stroke.RData", verbose=TRUE)                # Ischaemic stroke (from algorithmically-defined outcomes)
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_alzd.RData", verbose=TRUE)                  # Alzheimer's disease
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_multiple_sclerosis.RData", verbose=TRUE)    # Multiple sclerosis
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_myocardial_infarction.RData", verbose=TRUE) # Myocardial infarction
# load("/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_angina_pectoris.RData", verbose=TRUE)       # Angina pectoris
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/outcome_data_for_trmr_female_cancers.RData",
  verbose = TRUE
)          # BC, EC, or OC (female cancers)


library(sandwich)

gts   <- gts.2
snpid <- cbind(colnames(gts)[7:ncol(gts)])


# output (remember to change):

# out.data <-  dth.data.trmr[,c(1:32)]                   # all-cause mortality
# out.data <- arth.data.trmr[,c(1:32)]                   # OA
# out.data <- t2dm.data.trmr[,c(1:32)]                   # T2DM
# out.data <-  cnr.data.trmr[,c(1,2:4,14:41)]            # Gastric cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,5:7,14:41)]            # Kidney cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,8:10,14:41)]           # Esophagus cancer <- double-check to extract correct cancer type!
# out.data <-  cnr.data.trmr[,c(1,11:13,14:41)]          # Colorectal cancer <- double-check to extract correct cancer type!
# out.data <-  esc.data.trmr[,c(1,2:4,11:38)]            # Esophagus cancer (combined; same as above)
# out.data <-  esc.data.trmr[,c(1,5:7,11:38)]            # Esophagus adenocarcinoma
# out.data <-  esc.data.trmr[,c(1,8:10,11:38)]           # Esophagus squamous cell carcinoma
# out.data <- ulcc.data.trmr[,c(1:32)]                   # Ulcerative colitis
# out.data <- psor.data.trmr[,c(1:32)]                   # Psoriasis
# out.data <-  cad.data.trmr[,c(1:32)]                   # CAD
# out.data <-  strk.data.trmr[,c(1:32)]                  # Ischaemic stroke
# out.data <-  cad.data.trmr[,c(1,48:50,5:32)]           # Atrial fibrillation
# out.data <- myoc.data.trmr[,c(1:32)]                   # Myocardial infarction
# out.data <- apec.data.trmr[,c(1:32)]                   # Angina pectoris
# out.data <- alzd.data.trmr[,c(1:32)]                   # Alzheimer's diseasse
# out.data <- mscl.data.trmr[,c(1:32)]                   # Multiple sclerosis
out.data <-
  cnr.data.trmr[, c(1, 2:4, 11:38)]             # Breast cancer
# out.data <- cnr.data.trmr[,c(1,5:7,11:38)]             # Endometrial cancer
# out.data <- cnr.data.trmr[,c(1,8:10,11:38)]            # Ovarian cancer

colnames(out.data)

snp.out <- matrix(rep(0, nrow(snpid) * 7), nrow(snpid), 7)

out.data$event <- out.data[, 2]
out.data$tstop <- out.data[, 3]

out.data$subset <-
  (seq(1, nrow(gts), by = 1) %in% w.half.2)   # same subsample x in which the variants of gts.x were identified
length(which(out.data$subset))
for (j in 7:ncol(gts))
{
  k <- j - 6
  out.data$gts  <- gts[, j]             # SNP k = j - 6
  out.data.omit <- na.omit(out.data)   # remove all "NA"s

  healthy <- which(out.data.omit$subset)

  # remember to change event.xx and tstop.xx in both p.full, p.reduced, and p.corr:

  p.full    <-
    glm(
      event ~ gts + sex * poly(age_at_eof, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )
  p.reduced <-
    glm(
      event ~ sex * poly(age_at_eof, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )
  p.corr    <-
    glm(
      event ~ gts + offset(log(tstop)),
      family = "poisson",
      subset = healthy,
      data = out.data.omit
    )

  p.model <- p.full
  snp.out[k, 1] <- summary(p.model)$coef[2, 1]
  cov.p <-
    vcovHC(p.model, type = "HC0")   # accounting for heteroscedasticity using White's estimator
  std.err <- sqrt(diag(cov.p))
  snp.out[k, 2] <- std.err[2]
  snp.out[k, 3] <-
    2 * pnorm(abs(snp.out[k, 1] / snp.out[k, 2]), 0, 1, lower.tail = FALSE)
  snp.out[k, 4] <-
    (summary(p.reduced)$deviance - summary(p.full)$deviance) / summary(p.full)$null.deviance   # squared semipartial correlation
  snp.out[k, 5] <-
    1 - summary(p.corr)$deviance / summary(p.corr)$null.deviance                               # squared correlation
  snp.out[k, 6] <-
    0.5 * log((1 + sqrt(snp.out[k, 4])) / (1 - sqrt(snp.out[k, 4])))                                   # Fisher's z-transformation
  snp.out[k, 7] <-
    length(healthy)                                                                          # sample size
  n.individuals <- length(p.corr$residuals)
}
snp.out[1:10, ]
snp.out.bc.rsq.2 <- snp.out   # change name!



# exposure:
exp.data <- exp.data.trmr
exp.data$tstop <-
  out.data$tstop   # this extra column is there in order to enable estimation of the exposure effect in the correct sample (e.g., female-specific cancers)

snp.exp <- matrix(rep(0, nrow(snpid) * 7), nrow(snpid), 7)

exp.data$subset <-
  (seq(1, nrow(gts), by = 1) %in% w.half.2)   # NOTE: estimated in the same (selection) sample, to avoid bias by selecting SNPs in the target sample
for (j in 7:ncol(gts))
{
  k <- j - 6
  exp.data$gts  <- gts[, j]             # SNP k = j - 6
  exp.data.omit <- na.omit(exp.data)   # remove all "NA"s

  healthy <- which(exp.data.omit$subset)

  p.full <-
    glm(
      bmi.ztf ~ gts + sex * poly(age_at_assessment, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre,
      family = "gaussian",
      subset = healthy,
      data = exp.data.omit
    )
  p.reduced <-
    glm(
      bmi.ztf ~ sex * poly(age_at_assessment, 3) + array + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + centre,
      family = "gaussian",
      subset = healthy,
      data = exp.data.omit
    )
  p.corr <-
    glm(bmi.ztf ~ gts,
        family = "gaussian",
        subset = healthy,
        data = exp.data.omit)

  p.model <- p.full
  snp.exp[k, 1:3] <-
    summary(p.model)$coef[2, c(1, 2, 4)]    # main effect
  snp.exp[k, 4]   <-
    (summary(p.reduced)$deviance - summary(p.full)$deviance) / summary(p.full)$null.deviance   # squared semipartial correlation
  snp.exp[k, 5]   <-
    1 - summary(p.corr)$deviance / summary(p.corr)$null.deviance                               # squared correlation
  snp.exp[k, 6]   <-
    0.5 * log((1 + sqrt(snp.exp[k, 4])) / (1 - sqrt(snp.exp[k, 4])))                                   # Fisher's z-transformation
  snp.exp[k, 7]   <-
    length(healthy)                                                                          # sample size

}
snp.exp[1:10, ]
snp.exp.bc.rsq.2 <- snp.exp   # change name!

# test of equality beetween correlated correlations:
rho.xy <- median(sqrt(snp.out[, 4] / snp.exp[, 4]))
rm2 <- 0.5 * (snp.exp[, 4] + snp.out[, 4])
f   <- 0.5 * (1 - rho.xy) / (1 - rm2)
h   <- (1 - f * rm2) / (1 - rm2)
zstat.spc <-
  (snp.exp[, 6] - snp.out[, 6]) * sqrt(snp.out[, 7] - 3) / sqrt(2 * h * (1 - rho.xy))

rho.xy <- median(sqrt(snp.out[, 5] / snp.exp[, 5]))
rm2 <- 0.5 * (snp.exp[, 5] + snp.out[, 5])
f   <- 0.5 * (1 - rho.xy) / (1 - rm2)
h   <- (1 - f * rm2) / (1 - rm2)
zstat.cor <-
  (atanh(sqrt(snp.exp[, 5])) - atanh(sqrt(snp.out[, 5]))) * sqrt(snp.out[, 7] -
                                                                   3) / sqrt(2 * h * (1 - rho.xy))

pval.spc <- pnorm(zstat.spc)
pval.cor <- pnorm(zstat.cor)

fct.a <- sqrt(2 * h * (1 - rho.xy))

psi   <-
  rho.xy * (1 - snp.exp[, 5] - snp.out[, 5]) - 0.5 * sqrt(snp.out[, 5] *
                                                            snp.exp[, 5]) * (1 - snp.exp[, 5] - snp.out[, 5] - rho.xy ^ 2)
chat  <- psi / ((1 - snp.exp[, 5]) * (1 - snp.out[, 5]))
fct.b <- sqrt(2 - 2 * chat)

save(
  snp.exp.bc.rsq.2,
  snp.out.bc.rsq.2,
  zstat.spc,
  zstat.cor,
  pval.spc,
  pval.cor,
  fct.a,
  fct.b,
  n.individuals,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_bc_2.RData"
)   # change name!



