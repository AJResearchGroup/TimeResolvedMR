# Main_R_scripts_to_be_cleaned_for_publication.R
# Created by Torgny Karlsson on 2024-02-19.


# Brief work flow:
# 0) GWAS and clumping (to find lead SNPs associated with BMI) in two independent subsamples of UKB is performed prior to this point
# 1) Steiger filtering is performed within the separate subsamples to identify BMI-SNPs that are suspected to be primarily associated with the outcome disease, not BMI
# 2) SNPs suspected to primarily associated with outcome are removed from the SNP sets, making use of the Steiger filtering in step 1
# 3) Generate disease-specific polygenic scores, making use of the filtered SNP sets identified in step 2
#### HERE IS THE INTERESTING STUFF:
# 4) Estimate effects of PGS on BMI (betaG[t]) and on outcome (Eta[t]), making use of the PGSs generated in step 3
# 5) Estimate causal life-course effects of BMI on outcome (Gamma[t]), making use of betaG and Eta from step 4


source("step1.R")
source("step2.R")
source("step3.R")
source("step5.R")
