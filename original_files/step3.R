############################################
# Script 3: Generate disease-specific PGSs #
############################################

# beta effects for BMI, from PLINK:
beta.plink.1 <-
  read.table(
    "DATA/SNP_LISTS/lead_snps_all_data_bmi_1_ordered.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )   # beta effect is found in column #11
beta.plink.2 <-
  read.table(
    "DATA/SNP_LISTS/lead_snps_all_data_bmi_2_ordered.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )   # beta effect is found in column #11

load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/GENETIC_INSTR/BMI/genetic_instr_gts_bmi_12_no_filtering.Rata",
  verbose = TRUE
)   # BMI SNPs with dosage values
load(
  "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/data_for_trmr.RData",
  verbose = TRUE
)   # needed for w.half.1 and w.half.2



# T2DM, with stringent Steiger filtering (now default):
load("DATA/snps_selected_for_PGS_steiger_t2dm.RData", verbose = TRUE)

gts.tmp.1 <- gts.1[, 7:dim(gts.1)[2]]
gts.tmp.2 <- gts.2[, 7:dim(gts.2)[2]]

pgs.dos.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.slct.dm.2]) %*% cbind(beta.plink.2[w.slct.dm.2, 11])   # selection/estimation performed in sample 2 while PGS calculated in sample 1
pgs.dos.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.slct.dm.1]) %*% cbind(beta.plink.1[w.slct.dm.1, 11])   # selection/estimation performed in sample 1 while PGS calculated in sample 2

pgs.dm.1 <-
  (pgs.dos.1 - mean(pgs.dos.1, na.rm = TRUE)) / sd(pgs.dos.1, na.rm = TRUE)
pgs.dm.2 <-
  (pgs.dos.2 - mean(pgs.dos.2, na.rm = TRUE)) / sd(pgs.dos.2, na.rm = TRUE)

pgs.dm   <- exp.data.trmr[, 1]
pgs.dm[w.half.1] <- pgs.dm.1
pgs.dm[w.half.2] <- pgs.dm.2

save(pgs.dm,
     pgs.dm.1,
     pgs.dm.2,
     pgs.dos.1,
     pgs.dos.2,
     w.half.1,
     w.half.2,
     file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_t2dm.RData")


# T2DM with stringent Steiger filtering, cluster splitting (MR-AHC) using independent GWAS meta-analysed results: Weighted allele score (now default)
load("DATA/snps_selected_for_PGS_AHC_clustering_steiger_t2dm.RData",
     verbose = TRUE)

gts.tmp.1 <-
  gts.1[, 7:dim(gts.1)[2]]   # 122 SNPs, identified in sample 1
gts.tmp.2 <-
  gts.2[, 7:dim(gts.2)[2]]   # 106 SNps, identified in sample 2

# the last figure .x in the name, refers to in which subsample the PGS is calculated (w.half.x),
# not in which subsample (y) the SNPs are identified and effects are estimated (w.snp.t2dm.cluster.y)
# the part: "pos.a", "neg.a", etc. denotes cluster classification:

pgs.dos.pos.2.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.t2dm.pos.2.2]) %*% cbind(beta.plink.2[w.snp.t2dm.pos.2.2, 11])   # cluster 1: slct/est in sample 2, PGS calc in sample 1
pgs.dos.pos.2.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.t2dm.pos.2.1]) %*% cbind(beta.plink.1[w.snp.t2dm.pos.2.1, 11])   # cluster 1: slct/est in sample 1, PGS calc in sample 2

pgs.dos.pos.1.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.t2dm.pos.1.2]) %*% cbind(beta.plink.2[w.snp.t2dm.pos.1.2, 11])   # cluster 2: slct/est in sample 2, PGS calc in sample 1
pgs.dos.pos.1.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.t2dm.pos.1.1]) %*% cbind(beta.plink.1[w.snp.t2dm.pos.1.1, 11])   # cluster 2: slct/est in sample 1, PGS calc in sample 2

pgs.dos.neg.1.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.t2dm.neg.1.2]) %*% cbind(beta.plink.2[w.snp.t2dm.neg.1.2, 11])   # cluster 3: slct/est in sample 2, PGS calc in sample 1
pgs.dos.neg.1.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.t2dm.neg.1.1]) %*% cbind(beta.plink.1[w.snp.t2dm.neg.1.1, 11])   # cluster 3: slct/est in sample 1, PGS calc in sample 2


pgs.dm.pos.2.1 <-
  (pgs.dos.pos.2.1 - mean(pgs.dos.pos.2.1, na.rm = TRUE)) / sd(pgs.dos.pos.2.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 1
pgs.dm.pos.2.2 <-
  (pgs.dos.pos.2.2 - mean(pgs.dos.pos.2.2, na.rm = TRUE)) / sd(pgs.dos.pos.2.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 1

pgs.dm.pos.1.1 <-
  (pgs.dos.pos.1.1 - mean(pgs.dos.pos.1.1, na.rm = TRUE)) / sd(pgs.dos.pos.1.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 2
pgs.dm.pos.1.2 <-
  (pgs.dos.pos.1.2 - mean(pgs.dos.pos.1.2, na.rm = TRUE)) / sd(pgs.dos.pos.1.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 2

pgs.dm.neg.1.1 <-
  (pgs.dos.neg.1.1 - mean(pgs.dos.neg.1.1, na.rm = TRUE)) / sd(pgs.dos.neg.1.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 3
pgs.dm.neg.1.2 <-
  (pgs.dos.neg.1.2 - mean(pgs.dos.neg.1.2, na.rm = TRUE)) / sd(pgs.dos.neg.1.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 3


pgs.dm.pos.2   <- exp.data.trmr[, 1]
pgs.dm.pos.2[w.half.1] <- pgs.dm.pos.2.1
pgs.dm.pos.2[w.half.2] <- pgs.dm.pos.2.2

pgs.dm.pos.1   <- exp.data.trmr[, 1]
pgs.dm.pos.1[w.half.1] <- pgs.dm.pos.1.1
pgs.dm.pos.1[w.half.2] <- pgs.dm.pos.1.2

pgs.dm.neg.1   <- exp.data.trmr[, 1]
pgs.dm.neg.1[w.half.1] <- pgs.dm.neg.1.1
pgs.dm.neg.1[w.half.2] <- pgs.dm.neg.1.2

save(
  pgs.dm.pos.2,
  pgs.dm.pos.2.1,
  pgs.dm.pos.2.2,
  pgs.dos.pos.2.1,
  pgs.dos.pos.2.2,
  pgs.dm.pos.1,
  pgs.dm.pos.1.1,
  pgs.dm.pos.1.2,
  pgs.dos.pos.1.1,
  pgs.dos.pos.1.2,
  pgs.dm.neg.1,
  pgs.dm.neg.1.1,
  pgs.dm.neg.1.2,
  pgs.dos.neg.1.1,
  pgs.dos.neg.1.2,
  w.half.1,
  w.half.2,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_AHC_clustering_steiger_t2dm.RData"
)



# Osteoarthritis, with stringent Steiger filtering (now default, from: June 2023):
load("DATA/snps_selected_for_PGS_steiger_osteoarthritis.RData",
     verbose = TRUE)

gts.tmp.1 <- gts.1[, 7:dim(gts.1)[2]]
gts.tmp.2 <- gts.2[, 7:dim(gts.2)[2]]

pgs.dos.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.slct.oa.2]) %*% cbind(beta.plink.2[w.slct.oa.2, 11])   # selection/estimation performed in sample 2 while PGS calculated in sample 1
pgs.dos.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.slct.oa.1]) %*% cbind(beta.plink.1[w.slct.oa.1, 11])   # selection/estimation performed in sample 1 while PGS calculated in sample 2

pgs.oa.1 <-
  (pgs.dos.1 - mean(pgs.dos.1, na.rm = TRUE)) / sd(pgs.dos.1, na.rm = TRUE)
pgs.oa.2 <-
  (pgs.dos.2 - mean(pgs.dos.2, na.rm = TRUE)) / sd(pgs.dos.2, na.rm = TRUE)

pgs.oa   <- exp.data.trmr[, 1]
pgs.oa[w.half.1] <- pgs.oa.1
pgs.oa[w.half.2] <- pgs.oa.2

save(pgs.oa,
     pgs.oa.1,
     pgs.oa.2,
     pgs.dos.1,
     pgs.dos.2,
     w.half.1,
     w.half.2,
     file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_osteoarthritis.RData")




# CAD, with stringent Steiger filtering (now default, from June 2023):
load("DATA/snps_selected_for_PGS_steiger_cad.RData", verbose = TRUE)

gts.tmp.1 <- gts.1[, 7:dim(gts.1)[2]]
gts.tmp.2 <- gts.2[, 7:dim(gts.2)[2]]

pgs.dos.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.slct.cad.2]) %*% cbind(beta.plink.2[w.slct.cad.2, 11])   # selection/estimation performed in sample 2 while PGS calculated in sample 1
pgs.dos.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.slct.cad.1]) %*% cbind(beta.plink.1[w.slct.cad.1, 11])   # selection/estimation performed in sample 1 while PGS calculated in sample 2

pgs.cad.1 <-
  (pgs.dos.1 - mean(pgs.dos.1, na.rm = TRUE)) / sd(pgs.dos.1, na.rm = TRUE)
pgs.cad.2 <-
  (pgs.dos.2 - mean(pgs.dos.2, na.rm = TRUE)) / sd(pgs.dos.2, na.rm = TRUE)

pgs.cad   <- exp.data.trmr[, 1]
pgs.cad[w.half.1] <- pgs.cad.1
pgs.cad[w.half.2] <- pgs.cad.2

save(pgs.cad,
     pgs.cad.1,
     pgs.cad.2,
     pgs.dos.1,
     pgs.dos.2,
     w.half.1,
     w.half.2,
     file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_cad.RData")



# CAD with stringent Steiger filtering, cluster splitting (MR-AHC) using independent GWAS meta-analysed results: Weighted allele score (now default)
# "..._cad_24Oct_2023.RData" denotes clustering using ebi-a-GCST005194, suspected to be the correct GWAS form MR Base (this one has the lowest SE)

### NOTE!!: The three SNP sets/clusters saved in "snps_selected_for_PGS_AHC_clustering_steiger_cad_24Oct_2023.RData" are, unfortunately, switch so that the
### most harmful cluster is called "positive 1", while the moderately harmful cluster is called "positive 2"!!!
load(
  "DATA/snps_selected_for_PGS_AHC_clustering_steiger_cad_24Oct_2023.RData",
  verbose = TRUE
)


gts.tmp.1 <-
  gts.1[, 7:dim(gts.1)[2]]   # 122 SNPs, identified in sample 1
gts.tmp.2 <-
  gts.2[, 7:dim(gts.2)[2]]   # 106 SNps, identified in sample 2

# the last figure .x in the name, refers to in which subsample the PGS is calculated (w.half.x),
# not in which subsample (y) the SNPs are identified and effects are estimated (w.snp.cad.cluster.y)
# the part: "pos.a", "neg.a", etc. denotes cluster classification:

pgs.dos.pos.2.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.cad.pos.2.2]) %*% cbind(beta.plink.2[w.snp.cad.pos.2.2, 11])   # cluster 2: slct/est in sample 2, PGS calc in sample 1
pgs.dos.pos.2.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.cad.pos.2.1]) %*% cbind(beta.plink.1[w.snp.cad.pos.2.1, 11])   # cluster 2: slct/est in sample 1, PGS calc in sample 2

pgs.dos.pos.1.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.cad.pos.1.2]) %*% cbind(beta.plink.2[w.snp.cad.pos.1.2, 11])   # cluster 2: slct/est in sample 2, PGS calc in sample 1
pgs.dos.pos.1.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.cad.pos.1.1]) %*% cbind(beta.plink.1[w.snp.cad.pos.1.1, 11])   # cluster 2: slct/est in sample 1, PGS calc in sample 2

pgs.dos.neg.1.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.snp.cad.neg.1.2]) %*% cbind(beta.plink.2[w.snp.cad.neg.1.2, 11])   # cluster 3: slct/est in sample 2, PGS calc in sample 1
pgs.dos.neg.1.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.snp.cad.neg.1.1]) %*% cbind(beta.plink.1[w.snp.cad.neg.1.1, 11])   # cluster 3: slct/est in sample 1, PGS calc in sample 2


pgs.cad.pos.2.1 <-
  (pgs.dos.pos.2.1 - mean(pgs.dos.pos.2.1, na.rm = TRUE)) / sd(pgs.dos.pos.2.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 2
pgs.cad.pos.2.2 <-
  (pgs.dos.pos.2.2 - mean(pgs.dos.pos.2.2, na.rm = TRUE)) / sd(pgs.dos.pos.2.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 2

pgs.cad.pos.1.1 <-
  (pgs.dos.pos.1.1 - mean(pgs.dos.pos.1.1, na.rm = TRUE)) / sd(pgs.dos.pos.1.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 2
pgs.cad.pos.1.2 <-
  (pgs.dos.pos.1.2 - mean(pgs.dos.pos.1.2, na.rm = TRUE)) / sd(pgs.dos.pos.1.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 2

pgs.cad.neg.1.1 <-
  (pgs.dos.neg.1.1 - mean(pgs.dos.neg.1.1, na.rm = TRUE)) / sd(pgs.dos.neg.1.1, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 1, cluster 3
pgs.cad.neg.1.2 <-
  (pgs.dos.neg.1.2 - mean(pgs.dos.neg.1.2, na.rm = TRUE)) / sd(pgs.dos.neg.1.2, na.rm =
                                                                 TRUE)     # standardized PGS, in subsample 2, cluster 3


pgs.cad.pos.2   <- exp.data.trmr[, 1]
pgs.cad.pos.2[w.half.1] <- pgs.cad.pos.2.1
pgs.cad.pos.2[w.half.2] <- pgs.cad.pos.2.2

pgs.cad.pos.1   <- exp.data.trmr[, 1]
pgs.cad.pos.1[w.half.1] <- pgs.cad.pos.1.1
pgs.cad.pos.1[w.half.2] <- pgs.cad.pos.1.2

pgs.cad.neg.1   <- exp.data.trmr[, 1]
pgs.cad.neg.1[w.half.1] <- pgs.cad.neg.1.1
pgs.cad.neg.1[w.half.2] <- pgs.cad.neg.1.2

save(
  pgs.cad.pos.2,
  pgs.cad.pos.2.1,
  pgs.cad.pos.2.2,
  pgs.dos.pos.2.1,
  pgs.dos.pos.2.2,
  pgs.cad.pos.1,
  pgs.cad.pos.1.1,
  pgs.cad.pos.1.2,
  pgs.dos.pos.1.1,
  pgs.dos.pos.1.2,
  pgs.cad.neg.1,
  pgs.cad.neg.1.1,
  pgs.cad.neg.1.2,
  pgs.dos.neg.1.1,
  pgs.dos.neg.1.2,
  w.half.1,
  w.half.2,
  file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_AHC_clustering_steiger_cad_24Oct_2023.RData"
)




# Atrial fibrillation, with stringent Steiger filtering (now default; from June 2023):
load("DATA/snps_selected_for_PGS_steiger_atrial_fibrillation.RData",
     verbose = TRUE)

gts.tmp.1 <- gts.1[, 7:dim(gts.1)[2]]
gts.tmp.2 <- gts.2[, 7:dim(gts.2)[2]]

pgs.dos.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.slct.af.2]) %*% cbind(beta.plink.2[w.slct.af.2, 11])   # selection/estimation performed in sample 2 while PGS calculated in sample 1
pgs.dos.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.slct.af.1]) %*% cbind(beta.plink.1[w.slct.af.1, 11])   # selection/estimation performed in sample 1 while PGS calculated in sample 2

pgs.af.1 <-
  (pgs.dos.1 - mean(pgs.dos.1, na.rm = TRUE)) / sd(pgs.dos.1, na.rm = TRUE)
pgs.af.2 <-
  (pgs.dos.2 - mean(pgs.dos.2, na.rm = TRUE)) / sd(pgs.dos.2, na.rm = TRUE)

pgs.af   <- exp.data.trmr[, 1]
pgs.af[w.half.1] <- pgs.af.1
pgs.af[w.half.2] <- pgs.af.2

save(pgs.af,
     pgs.af.1,
     pgs.af.2,
     pgs.dos.1,
     pgs.dos.2,
     w.half.1,
     w.half.2,
     file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_PGS_steiger_atrial_fibrillation.RData")



# Generic weighted alleles score, including all SNPs:
gts.tmp.1 <- gts.1[, 7:dim(gts.1)[2]]
gts.tmp.2 <- gts.2[, 7:dim(gts.2)[2]]

w.slct.all.1 <- seq(1, 122, by = 1)
w.slct.all.2 <- seq(1, 106, by = 1)

pgs.dos.1 <-
  as.matrix(gts.tmp.2[w.half.1, w.slct.all.2]) %*% cbind(beta.plink.2[w.slct.all.2, 11])   # selection/estimation performed in sample 2 while PGS calculated in sample 1
pgs.dos.2 <-
  as.matrix(gts.tmp.1[w.half.2, w.slct.all.1]) %*% cbind(beta.plink.1[w.slct.all.1, 11])   # selection/estimation performed in sample 1 while PGS calculated in sample 2

pgs.all.1 <-
  (pgs.dos.1 - mean(pgs.dos.1, na.rm = TRUE)) / sd(pgs.dos.1, na.rm = TRUE)
pgs.all.2 <-
  (pgs.dos.2 - mean(pgs.dos.2, na.rm = TRUE)) / sd(pgs.dos.2, na.rm = TRUE)

pgs.all   <- exp.data.trmr[, 1]
pgs.all[w.half.1] <- pgs.all.1
pgs.all[w.half.2] <- pgs.all.2

save(pgs.all,
     pgs.all.1,
     pgs.all.2,
     pgs.dos.1,
     pgs.dos.2,
     w.half.1,
     w.half.2,
     file = "/proj/sens2017538/nobackup/torgny/TIME_RESOLVED_MR/DATA/input_data_for_TRMR_WAS_with_all_snps.RData")
