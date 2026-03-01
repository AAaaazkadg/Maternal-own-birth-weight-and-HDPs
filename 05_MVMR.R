# 03_MVMR
library(TwoSampleMR)
library(MVMR)
library(MendelianRandomization)
childhood_BMI_significant <- subset(childhood_BMI_data,p_value < 5E-08)
write.csv(BW_exp,"BW_exp_mediation.csv")
write.csv(childhood_BMI_significant_exposure_clumped,"childhood_BMI_exp_mediation.csv")
BW_GH_exp_MVMR <- read.csv("BW_exp_mediation.csv")
childhood_BMI_GH_exp_MVMR <- read.csv("childhood_BMI_exp_mediation.csv")

# Two sets of exposure IVs are ready----BW_GH_exp_MVMR, childhood_BMI_GH_exp_MVMR
# Merge BW_GH_exp_MVMR and childhood_BMI_GH_exp_MVMR by columns, keeping the same columns
# Extract the SNP column and p-value column from the exposure factors
SNP_BW_GH = BW_GH_exp_MVMR[,c("SNP","pval.exposure"), drop = FALSE]
SNP_childhood_BMI_GH = childhood_BMI_GH_exp_MVMR[,c("SNP","pval.exposure"), drop = FALSE]
# Merge exposure data and remove duplicates
SNP_all_GH = rbind(SNP_BW_GH, SNP_childhood_BMI_GH)
SNP_uni_GH = SNP_all_GH[!duplicated(SNP_all_GH$SNP),]
write.csv(SNP_uni_GH,"SNP_uni_GH.csv")
SNP_uni_GH <- read.csv("SNP_uni_GH.csv")
# After merging the exposure data, perform linkage disequilibrium clumping again (only p.val and SNP required).
exp_MVMR_GH_clumped <- clump_data(SNP_uni_GH,
                                  clump_kb = 10000,
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,pop="EUR")
# Remove 2 SNPs (20 + 17 - 2 = 35)

# Add a Markername(37)/(38) column to SNP_uni_GH
# Merge SNP with exposure data
# After merging, retain SNP and exposure's ea/nea/other data
merBW_GH <- merge(exp_MVMR_GH_clumped,BW_fetal_data, by.x = "SNP",by.y = "RSID")
merchildhood_BMI_GH <- merge(exp_MVMR_GH_clumped,bmi_data, by.x = "SNP",by.y = "variant_id")
# Change the name of exposure
merBW_GH$id.exposure = "BW"
merchildhood_BMI_GH$id.exposure = "Childhood_BMI"
write.csv(merBW_GH,file = "merBW_GH.csv")
write.csv(merchildhood_BMI_GH,file = "merchildhood_BMI_GH.csv")
merBW_GH <- read.csv("merBW_GH.csv")
merchildhood_BMI_GH <- read.csv("merchildhood_BMI_GH.csv")
# Merge SNPs and outcome
merGH <- merge(SNP_uni_GH,out_hdp_meta,by.x = 'MarkerName', by.y = "MarkerName")
write.csv(merGH,"merGH.csv")

out_dat_GH = format_data(dat = merGH,
                         type = "outcome",
                         snps = SNP_uni_GH$SNP,
                         snp_col = "SNP",
                         beta_col = "BETA",
                         pval_col = "P",
                         se_col = "SE",
                         eaf_col = "FRQ",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         ncase_col = NA,
                         ncontrol_col = NA)
out_dat_GH$id.outcome = "GH"
write.csv(out_dat_GH,"out_dat_GH.csv")
out_dat_GH <- read.csv("out_dat_GH.csv")

# Bind exposure data
merBW_GH <- read.csv("merBW_GH.csv")
merchildhood_BMI_GH <- read.csv("merchildhood_BMI_GH.csv")
expo_dat_GH = rbind(merBW_GH, merchildhood_BMI_GH)
expo_dat_GH <- read.csv("expo_dat_GH.csv")
write.csv(expo_dat_GH,"expo_dat_GH.csv")

expo_dat_GH <- read.csv("expo_dat_GH.csv")
expo_dat_GH$effect_allele.exposure <- toupper(expo_dat_GH$effect_allele.exposure)
expo_dat_GH$other_allele.exposure <- toupper(expo_dat_GH$other_allele.exposure)
mvmr_dat_GH = mv_harmonise_data(expo_dat_GH, out_dat_GH)
save(mvmr_dat_GH, file="mvmr_dat_GH.Rdata")

load("mvmr_dat_GH.Rdata")
mv_multiple(mvmr_dat_GH)
# Not sugnificant

SummaryStats_childhood_BMI = cbind(mvmr_dat_GH[["outcome_beta"]],
                                   mvmr_dat_GH[["exposure_beta"]][,1],
                                   mvmr_dat_GH[["exposure_beta"]][,2],
                                   mvmr_dat_GH[["exposure_se"]][,1],
                                   mvmr_dat_GH[["exposure_se"]][,2],
                                   mvmr_dat_GH[["outcome_se"]])
SummaryStats_childhood_BMI = data.frame(SummaryStats_childhood_BMI)

MVMR_Input_childhood_BMI = mr_mvinput(bx=cbind(SummaryStats_childhood_BMI$X2, SummaryStats_childhood_BMI$X3),
                                      bxse=cbind( SummaryStats_childhood_BMI$X4, SummaryStats_childhood_BMI$X5),
                                      by = SummaryStats_childhood_BMI$X1,
                                      byse = SummaryStats_childhood_BMI$X6)
ivw = mr_mvivw(MVMR_Input_childhood_BMI,
               model = "default",
               correl = FALSE,
               distribution = "normal",
               alpha = 0.05)
ivw
egger = mr_mvegger(MVMR_Input_childhood_BMI,
                   orientate = 1,
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05)
egger

r_input = format_mvmr(BXGs = mvmr_dat_GH[["exposure_beta"]],
                      BYG = mvmr_dat_GH[["outcome_beta"]],
                      seBXGs = mvmr_dat_GH[["exposure_se"]],
                      seBYG = mvmr_dat_GH[["outcome_se"]],
                      RSID = rownames(mvmr_dat_GH[["exposure_beta"]]))

mvmr_results <- ivw_mvmr(r_input)

pleiotropy_test <- pleiotropy_mvmr(r_input)
Fz = strength_mvmr(r_input = r_input, gencov = 0)

presso_results_BMI <-mr_presso(BetaOutcome = "X1",
                               BetaExposure = c("X2", "X3"), 
                               SdOutcome = "X6", 
                               SdExposure = c("X4", "X5"),
                               OUTLIERtest = TRUE, 
                               DISTORTIONtest = TRUE, 
                               data = SummaryStats_childhood_BMI,
                               NbDistribution = 2000, 
                               SignifThreshold = 0.05)
# rs17817449 rs4144829 rs61765651 rs75844534 
expo_dat_GH_presso <- read.csv("expo_dat_GH_presso.csv")
mvmr_dat_GH_presso = mv_harmonise_data(expo_dat_GH_presso, out_dat_GH)
mv_multiple(mvmr_dat_GH_presso)

SummaryStats_childhood_BMI_presso = cbind(mvmr_dat_GH_presso[["outcome_beta"]],
                                          mvmr_dat_GH_presso[["exposure_beta"]][,1],
                                          mvmr_dat_GH_presso[["exposure_beta"]][,2],
                                          mvmr_dat_GH_presso[["exposure_se"]][,1],
                                          mvmr_dat_GH_presso[["exposure_se"]][,2],
                                          mvmr_dat_GH_presso[["outcome_se"]])
SummaryStats_childhood_BMI_presso = data.frame(SummaryStats_childhood_BMI_presso)

MVMR_Input_childhood_BMI_presso = mr_mvinput(bx=cbind(SummaryStats_childhood_BMI_presso$X2, SummaryStats_childhood_BMI_presso$X3),
                                             bxse=cbind( SummaryStats_childhood_BMI_presso$X4, SummaryStats_childhood_BMI_presso$X5),
                                             by = SummaryStats_childhood_BMI_presso$X1,
                                             byse = SummaryStats_childhood_BMI_presso$X6)
ivw = mr_mvivw(MVMR_Input_childhood_BMI_presso,
               model = "default",
               correl = FALSE,
               distribution = "normal",
               alpha = 0.05)
ivw
egger = mr_mvegger(MVMR_Input_childhood_BMI_presso,
                   orientate = 1,
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05)
egger
r_input = format_mvmr(BXGs = mvmr_dat_GH_presso[["exposure_beta"]],
                      BYG = mvmr_dat_GH_presso[["outcome_beta"]],
                      seBXGs = mvmr_dat_GH_presso[["exposure_se"]],
                      seBYG = mvmr_dat_GH_presso[["outcome_se"]],
                      RSID = rownames(mvmr_dat_GH_presso[["exposure_beta"]]))

mvmr_results <- ivw_mvmr(r_input)
qhet_mvmr(r_input)
pleiotropy_test <- pleiotropy_mvmr(r_input)
Fz = strength_mvmr(r_input = r_input, gencov = 0)

# delta method calculate se and p value
res <- read.csv("Total_BW_childhood_BMI_GH.csv")

# Extract beta and se
# UVMR: BW -> CBMI 
b1 <- res[3, "beta"]    # BW to childhood_BMI beta
se1 <- res[3, "se"]     # BW to childhood_BMI se

# MVMR: childhood_BMI -> GH 
b2 <- res[7, "beta"]    # childhood_BM to GH beta(direct)
se2 <- res[7, "se"]     # childhood_BMI to GH se(direct)

# Total: BW -> GH 
total_beta <- res[2, "beta"]  # BW to GH total beta
total_se <- res[2, "se"]      # BW to GH total se

# Calculate indeirect effect
indirect_effect <- b1 * b2

# Calculate indirect se (delta method)
indirect_se <- sqrt(b1^2 * se2^2 + b2^2 * se1^2)

# Calculate z and p value
z_value <- indirect_effect / indirect_se
p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)

# Calculate 95% CI of direct effect
ci_lower <- indirect_effect - 1.96 * indirect_se
ci_upper <- indirect_effect + 1.96 * indirect_se

# Calculate mediation proportion
mediation_proportion <- indirect_effect / total_beta

# Calculate se of mediation proportion (delta method)
prop_se <- abs(mediation_proportion) * sqrt((indirect_se/indirect_effect)^2 + (total_se/total_beta)^2)

# Calculate 95% CI of mediation proportion
prop_ci_lower <- mediation_proportion - 1.96 * prop_se
prop_ci_upper <- mediation_proportion + 1.96 * prop_se

results_table <- data.frame(
  Analysis = c("indirect effect", "mediation proportion"),
  beta = c(round(indirect_effect, 6), round(mediation_proportion, 4)),
  se = c(round(indirect_se, 6), round(prop_se, 4)),
  lo_ci = c(round(ci_lower, 6), round(prop_ci_lower, 4)),
  up_ci = c(round(ci_upper, 6), round(prop_ci_upper, 4)),
  P = c(format.pval(p_value, digits = 3), NA)
)

print(results_table)

