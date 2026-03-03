library(TwoSampleMR)
library(data.table)
# BW → PE
# IVs with SNP
BW_significiant_SEM_exposure_clumped <- read.csv("BW_significiant_SEM_exposure_clumped.csv", header = T)
# IVs with Chr:Pos(GRCh38)
BW_exp_20 <- read.csv("BW_exp_20.csv")
# Description of the columns in the meta-analysis summary file _ Outcome
# Chromosome  - chromosome name
# Position    - chromosomal position
# Marker      - this is the marker name
# Allele1     - the first allele for this marker in the first file where it occurs
# Allele2     - the second allele for this marker in the first file where it occurs
# Freq1       - weighted average of frequency for allele 1 across all studies
# FreqSE      - corresponding standard error for allele frequency estimate
# MinFreq     - minimum frequency for allele 1 across all studies
# MaxFreq     - maximum frequency for allele 1 across all studies
# Effect      - overall estimated effect size for allele1
# StdErr      - overall standard error for effect size estimate
# P-value     - meta-analysis p-value
# Direction   - summary of effect direction for each study, with one '+' or '-' per study
# HetISq      - I^2 statistic which measures heterogeneity on scale of 0-100%
# HetChiSq    - chi-squared statistic in simple test of heterogeneity
# df          - degrees of freedom for heterogeneity statistic
# HetPVal     - P-value for heterogeneity statistic
out_hdp_meta <-  fread("metal_geshtn_European_allBiobanks_omitNone_1.txt", sep = "\t",header = TRUE)
Merge_hdp <- merge(BW_exp_20, out_hdp_meta, by.x = "Markername.GRCh38.", by.y = "MarkerName")
write.csv(Merge_hdp, "outcome_aftermerge_hdp.csv")
outcome_dat_gh <- read_outcome_data(
  snps = BW_exp$SNP,
  filename = "outcome_aftermerge_hdp.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P.value",
  eaf_col = "EAF"
)
outcome_dat_gh$outcome <- "GH"
gh_dat <- harmonise_data(
  exposure_dat = BW_significiant_SEM_exposure_clumped,
  outcome_dat = outcome_dat_gh,
  action = 2
)
write.csv(gh_dat, "Harmonization_gh.csv")
gh_dat<- read.csv("Harmonization_gh.csv",header = T)

# MR estimation 
set.seed(20250929)
res_01 <- mr(gh_dat)
mr(gh_dat, 
   parameters =default_parameters(),
   method_list=c('mr_ivw_mre',
                 'mr_egger_regression',
                 'mr_weighted_median'))
generate_odds_ratios(mr_res = res_01)
res_single <- mr_singlesnp(gh_dat,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR1 <- mr_forest_plot(res_single)
FOR1[[1]]
mr_scatter_plot(mr_results = res_01, gh_dat)
mr_heterogeneity(gh_dat)
mr_pleiotropy_test(gh_dat)
mr_leaveoneout_plot(mr_leaveoneout(gh_dat))
library(MRPRESSO)
gh_presso <- mr_presso(BetaOutcome = "beta.outcome",
                       BetaExposure = "beta.exposure",
                       SdOutcome = "se.outcome",
                       SdExposure = "se.exposure",
                       OUTLIERtest = TRUE,
                       DISTORTIONtest = TRUE,
                       data = gh_dat,
                       NbDistribution = 2000,
                       SignifThreshold = 0.05)
# Identify rs4144829 rs80278614 as outlier
ex1 <- which(gh_dat$SNP=="rs4144829")
gh_dat_presso <- gh_dat[-ex1, ]
ex2 <- which(gh_dat$SNP=="rs75844534")
gh_dat_presso <- gh_dat_presso[-ex2, ]
gh_dat_presso <- write.csv("Harmonization_gh_presso.csv",header = T)
#############################################
out_PE_meta <-  fread("metal_preec_European_allBiobanks_omitNone_1.txt", sep = "\t",header = TRUE)
Merge_PE <- merge(BW_exp_20, out_PE_meta, by.x = "Markername.GRCh38.", by.y = "MarkerName")
write.csv(Merge_PE, "outcome_aftermerge_PE.csv")
outcome_dat_PE <- read_outcome_data(
  snps = BW_exp$SNP,
  filename = "outcome_aftermerge_PE.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P.value",
  eaf_col = "EAF"
)
outcome_dat_PE$outcome <- "PE"
PE_dat <- harmonise_data(
  exposure_dat = BW_significiant_SEM_exposure_clumped,
  outcome_dat = outcome_dat_PE,
  action = 2
)
write.csv(PE_dat, "Harmonization_PE.csv")
PE_dat<- read.csv("Harmonization_PE.csv",header = T)
PE_dat<- read.csv("Harmonization_PE_presso.csv",header = T)
# MR estimation 
res_02 <- mr(PE_dat)
generate_odds_ratios(mr_res = res_02)
mr_scatter_plot(mr_results = res_02, PE_dat)
mr_heterogeneity(PE_dat)
mr_funnel_plot(res_02)
mr_pleiotropy_test(PE_dat)
mr_leaveoneout_plot(mr_leaveoneout(PE_dat))
res_single <- mr_singlesnp(PE_dat,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR2 <- mr_forest_plot(res_single)
FOR2[[1]]

PE_presso <- mr_presso(BetaOutcome = "beta.outcome",
                       BetaExposure = "beta.exposure",
                       SdOutcome = "se.outcome",
                       SdExposure = "se.exposure",
                       OUTLIERtest = TRUE,
                       DISTORTIONtest = TRUE,
                       data = PE_dat,
                       NbDistribution = 2000,
                       SignifThreshold = 0.05)
# Identify rs11698914 rs35261542 rs7076938 as outlier
ex1 <- which(PE_dat$SNP=="rs11698914")
PE_dat_presso <- PE_dat[-ex1, ]
ex2 <- which(PE_dat$SNP=="rs35261542")
PE_dat_presso <- PE_dat_presso[-ex2, ]
ex3 <- which(PE_dat$SNP=="rs7076938")
PE_dat_presso <- PE_dat_presso[-ex3, ]
PE_dat_presso <- write.csv("Harmonization_PE_presso.csv",header = T)
# F-statistic calculation
k <- 1
beta_values <- BW_exp_20$Beta     
se_values <- BW_exp_20$SE         
eaf_values <- BW_exp_20$EAF        
n_values <- BW_exp_20$N  
# Calculate R² of every SNP
BW_exp_20$R2 <- (2 * BW_exp_20$EAF * (1 - BW_exp_20$EAF) * BW_exp_20$Beta^2) / 
  (2 * BW_exp_20$EAF * (1 - BW_exp_20$EAF) * BW_exp_20$Beta^2 + 
     2 * BW_exp_20$EAF * (1 - BW_exp_20$EAF) * BW_exp_20$N * BW_exp_20$SE^2)
# Calculate F of every SNP
BW_exp_20$F_statistic <- BW_exp_20$R2 * (BW_exp_20$N - 2) / (1 - BW_exp_20$R2)
# R² = [2 × EAF × (1 - EAF) × β²] / [2 × EAF × (1 - EAF) × β² + 2 × EAF × (1 - EAF) × N × SE²]
# F = [(N - k - 1) / k] × [R² / (1 - R²)]
