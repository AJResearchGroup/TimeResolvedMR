###########################################################################
# Script 2: Remove invalid SNPs, identified in the Steiger filtering step #
###########################################################################
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/GENETIC_INSTR/BMI/genetic_instr_gts_bmi_12_no_filtering.Rata",
  verbose = TRUE
)   # to extract SNP rs-numbers

# T2DM, with more stringent Steiger filtering (default):
load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_t2dm_1.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.1 <- which(pval.fdr > 0.05)

load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_t2dm_2.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.2 <- which(pval.fdr > 0.05)

w.rmv.dm.1 <- w.dis.1
w.rmv.dm.2 <- w.dis.2

snp.rmv.dm.1 <- snpid.1[w.dis.1]
snp.rmv.dm.2 <- snpid.2[w.dis.2]

w.rmv.all.dm.1 <-
  w.rmv.dm.1   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms
w.rmv.all.dm.2 <-
  w.rmv.dm.2   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms

snp.rmv.all.dm.1 <- snpid.1[w.rmv.all.dm.1]
snp.rmv.all.dm.2 <- snpid.2[w.rmv.all.dm.2]

w.slct.dm.1 <-
  which(!(seq(1, nrow(snpid.1), by = 1) %in% w.rmv.all.dm.1))   # variants that are selected from the 1st sample and applied to the 2nd sample (PGS.2)
w.slct.dm.2 <-
  which(!(seq(1, nrow(snpid.2), by = 1) %in% w.rmv.all.dm.2))   # variants that are selected from the 2st sample and applied to the 1nd sample (PGS.1)

save(
  w.slct.dm.1,
  w.slct.dm.2,
  w.rmv.all.dm.1,
  w.rmv.all.dm.2,
  snp.rmv.all.dm.1,
  snp.rmv.all.dm.2,
  w.rmv.dm.1,
  w.rmv.dm.2,
  snp.rmv.dm.1,
  snp.rmv.dm.2,
  w.rmv.intr.1,
  w.rmv.intr.2,
  snp.rmv.intr.1,
  snp.rmv.intr.2,
  file = "DATA/snps_selected_for_PGS_steiger_t2dm.RData"
)


# Osteoarthritis, wcith more stringent Steiger filtering:
load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_osteoarthritis_1.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.1 <- which(pval.fdr > 0.05)

load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_osteoarthritis_2.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.2 <- which(pval.fdr > 0.05)

w.rmv.oa.1 <- w.dis.1
w.rmv.oa.2 <- w.dis.2

snp.rmv.oa.1 <- snpid.1[w.dis.1]
snp.rmv.oa.2 <- snpid.2[w.dis.2]

w.rmv.all.oa.1 <-
  w.rmv.oa.1   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms
w.rmv.all.oa.2 <-
  w.rmv.oa.2   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms

snp.rmv.all.oa.1 <- snpid.1[w.rmv.all.oa.1]
snp.rmv.all.oa.2 <- snpid.2[w.rmv.all.oa.2]

w.slct.oa.1 <-
  which(!(seq(1, nrow(snpid.1), by = 1) %in% w.rmv.all.oa.1))   # variants that are selected from the 1st sample and applied to the 2nd sample (PGS.2)
w.slct.oa.2 <-
  which(!(seq(1, nrow(snpid.2), by = 1) %in% w.rmv.all.oa.2))   # variants that are selected from the 2st sample and applied to the 1nd sample (PGS.1)

save(
  w.slct.oa.1,
  w.slct.oa.2,
  w.rmv.all.oa.1,
  w.rmv.all.oa.2,
  snp.rmv.all.oa.1,
  snp.rmv.all.oa.2,
  w.rmv.oa.1,
  w.rmv.oa.2,
  snp.rmv.oa.1,
  snp.rmv.oa.2,
  w.rmv.intr.1,
  w.rmv.intr.2,
  snp.rmv.intr.1,
  snp.rmv.intr.2,
  file = "DATA/snps_selected_for_PGS_steiger_osteoarthritis.RData"
)


# CAD, with more stringent Steiger filtering:
load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_cad_1.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.1 <- which(pval.fdr > 0.05)

load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_cad_2.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.2 <- which(pval.fdr > 0.05)

w.rmv.cad.1 <- w.dis.1
w.rmv.cad.2 <- w.dis.2

snp.rmv.cad.1 <- snpid.1[w.dis.1]
snp.rmv.cad.2 <- snpid.2[w.dis.2]

w.rmv.all.cad.1 <-
  w.rmv.cad.1   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms
w.rmv.all.cad.2 <-
  w.rmv.cad.2   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms

snp.rmv.all.cad.1 <- snpid.1[w.rmv.all.cad.1]
snp.rmv.all.cad.2 <- snpid.2[w.rmv.all.cad.2]

w.slct.cad.1 <-
  which(!(seq(1, nrow(snpid.1), by = 1) %in% w.rmv.all.cad.1))   # variants that are selected from the 1st sample and applied to the 2nd sample (PGS.2)
w.slct.cad.2 <-
  which(!(seq(1, nrow(snpid.2), by = 1) %in% w.rmv.all.cad.2))   # variants that are selected from the 2st sample and applied to the 1nd sample (PGS.1)

save(
  w.slct.cad.1,
  w.slct.cad.2,
  w.rmv.all.cad.1,
  w.rmv.all.cad.2,
  snp.rmv.all.cad.1,
  snp.rmv.all.cad.2,
  w.rmv.cad.1,
  w.rmv.cad.2,
  snp.rmv.cad.1,
  snp.rmv.cad.2,
  w.rmv.intr.1,
  w.rmv.intr.2,
  snp.rmv.intr.1,
  snp.rmv.intr.2,
  file = "DATA/snps_selected_for_PGS_steiger_cad.RData"
)



# Atrial fibrillation, with more stringent Steiger filtering:
load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_atrial_fibrillation_1.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.1 <- which(pval.fdr > 0.05)

load(
  "DATA/STEIGER_SELECTION_SAMPLE/snps_rsq_for_steiger_selection_atrial_fibrillation_2.RData",
  verbose = TRUE
)
pval.fdr <- p.adjust(pnorm(-zstat.cor), method = "fdr")
w.dis.2 <- which(pval.fdr > 0.05)

w.rmv.af.1 <- w.dis.1
w.rmv.af.2 <- w.dis.2

snp.rmv.af.1 <- snpid.1[w.dis.1]
snp.rmv.af.2 <- snpid.2[w.dis.2]

w.rmv.all.af.1 <-
  w.rmv.af.1   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms
w.rmv.all.af.2 <-
  w.rmv.af.2   # SNP(s) with the strongest time-tependences are NOT removed due to potential correlations with specific molecular sub-mechanisms

snp.rmv.all.af.1 <- snpid.1[w.rmv.all.af.1]
snp.rmv.all.af.2 <- snpid.2[w.rmv.all.af.2]

w.slct.af.1 <-
  which(!(seq(1, nrow(snpid.1), by = 1) %in% w.rmv.all.af.1))   # variants that are selected from the 1st sample and applied to the 2nd sample (PGS.2)
w.slct.af.2 <-
  which(!(seq(1, nrow(snpid.2), by = 1) %in% w.rmv.all.af.2))   # variants that are selected from the 2st sample and applied to the 1nd sample (PGS.1)

save(
  w.slct.af.1,
  w.slct.af.2,
  w.rmv.all.af.1,
  w.rmv.all.af.2,
  snp.rmv.all.af.1,
  snp.rmv.all.af.2,
  w.rmv.af.1,
  w.rmv.af.2,
  snp.rmv.af.1,
  snp.rmv.af.2,
  w.rmv.intr.1,
  w.rmv.intr.2,
  snp.rmv.intr.1,
  snp.rmv.intr.2,
  file = "DATA/snps_selected_for_PGS_steiger_atrial_fibrillation.RData"
)
