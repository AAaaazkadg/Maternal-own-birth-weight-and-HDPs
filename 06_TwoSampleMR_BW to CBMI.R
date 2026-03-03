# 01 BW → Childhood_BMI
BW_exp_20 <- read.csv("BW_exp_20.csv", header = T)
bmi_data <-  fread("CBMI_GWAS_summary_stat_2020.tsv", sep = "\t",header = TRUE)
Merge_BW_childhood_BMI <- merge(BW_exp_20, bmi_data, by.x = "SNP", by.y = "variant_id")
write.csv(Merge_BW_childhood_BMI, "outcome_aftermerge_BW_childhood_BMI.csv")
Merge_BW_childhood_BMI <- read.csv("outcome_aftermerge_BW_childhood_BMI.csv")
# The outcome is missing the EAF column and needs to be manually completed.
outcome_dat_BW_childhood_BMI <- read_outcome_data(
  snps = BW_exp_20$SNP,
  filename = "outcome_aftermerge_BW_childhood_BMI.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  eaf_col = "EAF"
)

outcome_dat_BW_childhood_BMI$outcome <- "childhood_BMI"
BW_childhood_BMI_dat <- harmonise_data(
  exposure_dat = BW_exp,
  outcome_dat = outcome_dat_BW_childhood_BMI,
  action = 2
)

write.csv(BW_childhood_BMI_dat, "Harmonization_BW_childhood_BMI.csv")
BW_childhood_BMI_dat <- read.csv("Harmonization_BW_childhood_BMI.csv")
BW_childhood_BMI_dat_presso <- read.csv("Harmonization_BW_childhood_BMI_presso.csv",header = T)

set.seed(20260201)
res_1 <- mr(BW_childhood_BMI_dat)
res_2 <- mr(BW_childhood_BMI_dat_presso)
mr_scatter_plot(mr_results = res_1, BW_childhood_BMI_dat)
mr_heterogeneity(BW_childhood_BMI_dat)
mr_pleiotropy_test(BW_childhood_BMI_dat)
mr_leaveoneout_plot(mr_leaveoneout(BW_childhood_BMI_dat))
BW_childhood_BMI_presso <- mr_presso(BetaOutcome = "beta.outcome",
                                     BetaExposure = "beta.exposure",
                                     SdOutcome = "se.outcome",
                                     SdExposure = "se.exposure",
                                     OUTLIERtest = TRUE,
                                     DISTORTIONtest = TRUE,
                                     data = BW_childhood_BMI_dat,
                                     NbDistribution = 1500,
                                     SignifThreshold = 0.05)
# No outlier were identified
