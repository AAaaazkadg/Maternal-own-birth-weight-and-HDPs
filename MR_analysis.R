# MR analysis
# Place this script and all input files in the same working directory before running it.

# 0. Load all required R packages ----------------------------------------------------
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(plinkbinr)
  library(ieugwasr)
  library(MRPRESSO)
  library(mrclust)
  library(ggplot2)
  library(dplyr)
  library(MVMR)
  library(MendelianRandomization)
})


# 1. Prepare and clump the birth weight instruments -------------------------------
assign_many <- function(...) {
  list2env(list(...), envir = parent.frame())
  invisible(NULL)
}
# Import birth weight data from the article supplementary material as the exposure.
BW_SEM <- read.csv("BW_2019_SEM.csv", header = TRUE)
BW_significiant_SEM_exposure <- subset(BW_SEM, Fetal.P.value < 5e-08)

# Identify the variants excluded in the original analysis in one command.
assign_many(
  ex1 = which(BW_significiant_SEM_exposure$SNP == "rs138715366"), # Rare variant: MAF < 0.01
  ex2 = which(BW_significiant_SEM_exposure$SNP == "rs11042596"),  # Imprinted-gene region
  ex3 = which(BW_significiant_SEM_exposure$SNP == "rs1801253"),  # Unclassified effect
  ex4 = which(BW_significiant_SEM_exposure$SNP == "rs10872678"), # MTA effect
  ex5 = which(BW_significiant_SEM_exposure$SNP == "rs560887") , # MNTA effect
  ex6 = which(BW_significiant_SEM_exposure$SNP == "rs4444073") # Non-significant with WLM P = 2.129 × 10⁻⁷
)
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[
  -unique(c(ex1, ex2, ex3, ex4, ex5, ex6)),
]
write.csv(
  BW_significiant_SEM_exposure,
  "BW_significiant_SEM_exposure.csv",
  row.names = FALSE
)

# Format the exposure data for TwoSampleMR.
BW_significiant_SEM_exposure <- read_exposure_data(
  "BW_significiant_SEM_exposure.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Fetal.Beta..SDs.",
  se_col = "Fetal.SE",
  effect_allele_col = "Effect.allele",
  other_allele_col = "Other.allele",
  pval_col = "Fetal.P.value"
)
BW_significiant_SEM_exposure$exposure <- "Birth weight"

# Perform online LD clumping.
BW_significiant_SEM_exposure_clumped <- clump_data(
  BW_significiant_SEM_exposure,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  pop = "EUR"
)

# Perform local LD clumping as a reproducible alternative.
# Keep the PLINK executable and EUR reference files in the working directory,
# or change only these two arguments for the local computer being used.
clumped <- ld_clump(
  dplyr::tibble(
    rsid = BW_significiant_SEM_exposure$SNP,
    pval = BW_significiant_SEM_exposure$pval.exposure,
    id = BW_significiant_SEM_exposure$id.exposure
  ),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  plink_bin = get_plink_exe(),
  bfile = "EUR"
)
# Four of 23 variants were removed because of LD or absence from the reference panel.

BW_exposure_selected <- BW_significiant_SEM_exposure_clumped %>%
  select(SNP = rsid) %>%
  left_join(
    BW_significiant_SEM_exposure,
    by = "SNP"
  )


write.csv(
  BW_exposure_selected,
  "BW_exposure_selected.csv",
  row.names = FALSE
)

BW_exposure_selected <- read.csv("BW_exposure_selected.csv")

# 2. Two-sample MR: birth weight to GH and PE --------------------------------------

# Read all primary files needed in this section.
assign_many(
  BW_exp = read.csv("BW_exp_19.csv", header = TRUE),
  out_GH_meta = fread(
    "metal_geshtn_European_allBiobanks_omitNone_1.txt",
    sep = "\t",
    header = TRUE
  ),
  out_PE_meta = fread(
    "metal_preec_European_allBiobanks_omitNone_1.txt",
    sep = "\t",
    header = TRUE
  )
)

# 2.1 Birth weight to gestational hypertension ------------------------------------

Merge_GH <- merge(
  BW_exp,
  out_GH_meta,
  by.x = "Markername.GRCh38.",
  by.y = "MarkerName"
)
write.csv(Merge_GH, "outcome_aftermerge_GH.csv", row.names = FALSE)

outcome_dat_GH <- read_outcome_data(
  snps = BW_exp$SNP,
  filename = "outcome_aftermerge_GH.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  eaf_col = "Freq1"
)
outcome_dat_GH$outcome <- "GH"

GH_dat <- harmonise_data(
  exposure_dat = BW_exposure_selected,
  outcome_dat = outcome_dat_GH,
  action = 2
)
write.csv(GH_dat, "Harmonization_GH.csv", row.names = FALSE)

# Estimate the main MR effects and perform sensitivity analyses.
set.seed(20250929)
res_01 <- mr(GH_dat)
mr(
  GH_dat,
  parameters = default_parameters(),
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
generate_odds_ratios(mr_res = res_01)
res_single <- mr_singlesnp(
  GH_dat,
  all_method = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
FOR1 <- mr_forest_plot(res_single)
FOR1[[1]]
mr_scatter_plot(mr_results = res_01, GH_dat)
mr_heterogeneity(GH_dat)
mr_pleiotropy_test(GH_dat)
mr_leaveoneout_plot(mr_leaveoneout(GH_dat))
set.seed(20250929)
GH_presso <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = GH_dat,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)
# Identify rs4144829 rs75844534 as outliers.
# Remove the 2 MR-PRESSO outliers identified in the original analysis.
assign_many(
  ex1 = which(GH_dat$SNP == "rs4144829"),
  ex2 = which(GH_dat$SNP == "rs75844534")
)
GH_dat_presso <- GH_dat[-unique(c(ex1, ex2)), ]
write.csv(GH_dat_presso, "Harmonization_GH_presso.csv", row.names = FALSE)

set.seed(20250929)
res_03 <- mr(GH_dat_presso)
mr(
  GH_dat_presso,
  parameters = default_parameters(),
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
generate_odds_ratios(mr_res = res_03)
res_single <- mr_singlesnp(
  GH_dat_presso,
  all_method = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
FOR1 <- mr_forest_plot(res_single)
FOR1[[1]]
mr_scatter_plot(mr_results = res_03, GH_dat_presso)
mr_heterogeneity(GH_dat_presso)
mr_pleiotropy_test(GH_dat_presso)
mr_leaveoneout_plot(mr_leaveoneout(GH_dat_presso))


# 2.2 Birth weight to preeclampsia ------------------------------------------------

Merge_PE <- merge(
  BW_exp,
  out_PE_meta,
  by.x = "Markername.GRCh38.",
  by.y = "MarkerName"
)
write.csv(Merge_PE, "outcome_aftermerge_PE.csv", row.names = FALSE)

outcome_dat_PE <- read_outcome_data(
  snps = BW_exp$SNP,
  filename = "outcome_aftermerge_PE.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  eaf_col = "Freq1"
)
outcome_dat_PE$outcome <- "PE"

PE_dat <- harmonise_data(
  exposure_dat = BW_exposure_selected,
  outcome_dat = outcome_dat_PE,
  action = 2
)
write.csv(PE_dat, "Harmonization_PE.csv", row.names = FALSE)

# Estimate the main MR effects and perform sensitivity analyses.
set.seed(20250929)
res_02 <- mr(PE_dat)
generate_odds_ratios(mr_res = res_02)
mr_scatter_plot(mr_results = res_02, PE_dat)
mr_heterogeneity(PE_dat)
mr_pleiotropy_test(PE_dat)
mr_leaveoneout_plot(mr_leaveoneout(PE_dat))
res_single <- mr_singlesnp(
  PE_dat,
  all_method = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
FOR2 <- mr_forest_plot(res_single)
FOR2[[1]]

PE_presso <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = PE_dat,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)

# Remove the three MR-PRESSO outliers identified in the original analysis.
assign_many(
  ex1 = which(PE_dat$SNP == "rs11698914"),
  ex2 = which(PE_dat$SNP == "rs35261542"),
  ex3 = which(PE_dat$SNP == "rs7076938")
)
PE_dat_presso <- PE_dat[-unique(c(ex1, ex2, ex3)), ]
write.csv(PE_dat_presso, "Harmonization_PE_presso.csv", row.names = FALSE)

set.seed(20250929)
res_04 <- mr(PE_dat_presso)
mr(
  PE_dat_presso,
  parameters = default_parameters(),
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
generate_odds_ratios(mr_res = res_04)
res_single <- mr_singlesnp(
  PE_dat_presso,
  all_method = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)
FOR1 <- mr_forest_plot(res_single)
FOR1[[1]]
mr_scatter_plot(mr_results = res_04, PE_dat_presso)
mr_heterogeneity(PE_dat_presso)
mr_pleiotropy_test(PE_dat_presso)
mr_leaveoneout_plot(mr_leaveoneout(PE_dat_presso))

# Calculate the F statistic for each birth-weight instrument.
assign_many(
  k = 1,
  beta_values = BW_exp$Beta,
  se_values = BW_exp$SE,
  eaf_values = BW_exp$EAF,
  n_values = BW_exp$N
)
BW_exp$R2 <- (
  2 * BW_exp$EAF * (1 - BW_exp$EAF) * BW_exp$Beta^2
) / (
  2 * BW_exp$EAF * (1 - BW_exp$EAF) * BW_exp$Beta^2 +
    2 * BW_exp$EAF * (1 - BW_exp$EAF) * BW_exp$N * BW_exp$SE^2
)
BW_exp$F_statistic <- BW_exp$R2 * (BW_exp$N - 2) / (1 - BW_exp$R2)


# 3. MR-Clust analysis: birth weight to gestational hypertension -------------------

# Reuse the harmonised object created in Section 2.
harmonised_data <- GH_dat
assign_many(
  exposure_abbr = "BW",
  outcome_abbr = "GH"
)

# Extract the variables required for MR-Clust.
mr_clust_data_BW_GH <- data.frame(
  SNP = harmonised_data$SNP,
  betaX_BW = harmonised_data$beta.exposure,
  betaY_GH = harmonised_data$beta.outcome,
  seX_BW = harmonised_data$se.exposure,
  seY_GH = harmonised_data$se.outcome,
  pvalX_BW = harmonised_data$pval.exposure,
  pvalY_GH = harmonised_data$pval.outcome
)
assign_many(
  theta_BW_GH = mr_clust_data_BW_GH$betaY_GH / mr_clust_data_BW_GH$betaX_BW,
  theta_se_BW_GH = mr_clust_data_BW_GH$seY_GH /
    abs(mr_clust_data_BW_GH$betaX_BW)
)
mr_clust_data_BW_GH$theta_BW_GH <- theta_BW_GH
mr_clust_data_BW_GH$theta_se_BW_GH <- theta_se_BW_GH

# Run the MR-Clust expectation-maximisation analysis.
set.seed(12345)
res_em <- mr_clust_em(
  theta = mr_clust_data_BW_GH$theta_BW_GH,
  theta_se = mr_clust_data_BW_GH$theta_se_BW_GH,
  bx = mr_clust_data_BW_GH$betaX_BW,
  by = mr_clust_data_BW_GH$betaY_GH,
  bxse = mr_clust_data_BW_GH$seX_BW,
  byse = mr_clust_data_BW_GH$seY_GH,
  obs_names = mr_clust_data_BW_GH$SNP
)
head(res_em$results$all, n = 24)
head(res_em$results$best, n = 24)

# Create and display the cluster plot.
cluster_plot <- res_em$plots$two_stage +
  ggplot2::ggtitle(
    "MR-Clust Analysis: Birth Weight → Gestational Hypertension"
  ) +
  ggplot2::xlim(
    0,
    max(abs(mr_clust_data_BW_GH$betaX_BW) +
          2 * mr_clust_data_BW_GH$seX_BW)
  ) +
  ggplot2::xlab("Genetic association with Birth Weight (betaX)") +
  ggplot2::ylab(
    "Genetic association with Gestational Hypertension (betaY)"
  )
cluster_plot

write.csv(mr_clust_data_BW_GH, "mr_clust_data_BW_GH.csv", row.names = FALSE)
write.csv(res_em$results$best, "cluster_result_BW_GH.csv", row.names = FALSE)

# Apply the conservative assignment criterion: probability >= 0.7.
filtered_data <- pr_clust(
  dta = res_em$results$best,
  prob = 0.7,
  min_obs = 1
)
write.csv(
  filtered_data,
  "cluster_result_BW_GH_conservation.csv",
  row.names = FALSE
)

keep_indices <- which(
  res_em$results$best$observation %in% filtered_data$observation
)
assign_many(
  filtered_bx = mr_clust_data_BW_GH$betaX_BW[keep_indices],
  filtered_by = mr_clust_data_BW_GH$betaY_GH[keep_indices],
  filtered_bxse = mr_clust_data_BW_GH$seX_BW[keep_indices],
  filtered_byse = mr_clust_data_BW_GH$seY_GH[keep_indices],
  filtered_rsid = mr_clust_data_BW_GH$SNP[keep_indices]
)

conservative_plot <- two_stage_plot(
  res = filtered_data,
  bx = filtered_bx,
  by = filtered_by,
  bxse = filtered_bxse,
  byse = filtered_byse,
  obs_names = filtered_rsid
)

# Define the cluster colour scheme used for both plots.
my_colors <- c(
  "junk" = "black",
  "null" = "gray50",
  "1" = "#E18727FF",
  "2" = "#0072B5FF"
)

cluster_plot_custom <- res_em$plots$two_stage +
  ggplot2::ggtitle(
    "MR-Clust Analysis: Birth Weight → Gestational Hypertension"
  ) +
  ggplot2::xlim(
    0,
    max(abs(mr_clust_data_BW_GH$betaX_BW) +
          2 * mr_clust_data_BW_GH$seX_BW)
  ) +
  ggplot2::xlab("Genetic association with Birth Weight (betaX)") +
  ggplot2::ylab(
    "Genetic association with Gestational Hypertension (betaY)"
  ) +
  ggplot2::scale_color_manual(values = my_colors) +
  ggplot2::scale_fill_manual(values = my_colors)
print(cluster_plot_custom)

conservative_plot <- conservative_plot +
  ggplot2::ggtitle(
    "MR-Clust Analysis: Birth Weight → Gestational Hypertension"
  ) +
  ggplot2::xlim(
    0,
    max(abs(mr_clust_data_BW_GH$betaX_BW) +
          2 * mr_clust_data_BW_GH$seX_BW)
  ) +
  ggplot2::xlab("Genetic association with Birth Weight") +
  ggplot2::ylab(
    "Genetic association with Gestational Hypertension"
  ) +
  ggplot2::scale_color_manual(values = my_colors) +
  ggplot2::scale_fill_manual(values = my_colors)
print(conservative_plot)


# 4. Two-sample MR: birth weight to childhood BMI ---------------------------------
childhood_BMI_data <- fread(
  "BMI.summary_stat_EBI_valid.tsv",
  sep = "\t",
  header = TRUE
)
BW_fetal_data <- fread(
  "Fetal_Effect_European_meta_NG2019.txt",
  sep = "\t",
  header = TRUE
)
Merge_BW_childhood_BMI <- merge(
  BW_exp,
  childhood_BMI_data,
  by.x = "SNP",
  by.y = "variant_id"
)
write.csv(
  Merge_BW_childhood_BMI,
  "outcome_aftermerge_BW_CBMI.csv",
  row.names = FALSE
)

# The original outcome file did not contain EAF; complete this column before use.

outcome_aftermerge_BW_CBMI <- read_outcome_data(
  snps = BW_exp$SNP,
  filename = "outcome_aftermerge_BW_CBMI.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  eaf_col = "EAF"
)
outcome_aftermerge_BW_CBMI$outcome <- "childhood_BMI"

BW_childhood_BMI_dat <- harmonise_data(
  exposure_dat = BW_exposure_selected,
  outcome_dat = outcome_aftermerge_BW_CBMI,
  action = 2
)
write.csv(
  BW_childhood_BMI_dat,
  "Harmonization_BW_childhood_BMI.csv",
  row.names = FALSE
)
BW_childhood_BMI_dat <- read.csv("Harmonization_BW_childhood_BMI.csv")
set.seed(20260201)
res_1 <- mr(BW_childhood_BMI_dat)
generate_odds_ratios(mr_res = res_1)
mr_scatter_plot(mr_results = res_1, BW_childhood_BMI_dat)
mr_heterogeneity(BW_childhood_BMI_dat)
mr_pleiotropy_test(BW_childhood_BMI_dat)
mr_leaveoneout_plot(mr_leaveoneout(BW_childhood_BMI_dat))

BW_childhood_BMI_presso <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = BW_childhood_BMI_dat,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)
#  No outlier were identified, therefore the results for the outlier-corrected MR are set to NA


# 5. Multivariable MR and mediation analysis ---------------------------------------

# Prepare and clump the childhood-BMI instruments.
childhood_BMI_significant <- subset(childhood_BMI_data, p_value < 5e-08)
write.csv(
  childhood_BMI_significant,
  "childhood_BMI_significant.csv",
  row.names = FALSE
)
childhood_BMI_exposure <- read_exposure_data(
  "childhood_BMI_significant.csv",
  sep = ",",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value"
)
childhood_BMI_exposure$exposure <- "Childhood BMI"
childhood_BMI_significant_clumped <- clump_data(
  childhood_BMI_exposure,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  pop = "EUR"
)
# Removing 1336 of 1353 variants due to LD with other variants or absence from LD reference panel

write.csv(BW_exp, "BW_exp_mediation.csv", row.names = FALSE)
write.csv(
  childhood_BMI_significant_clumped,
  "childhood_BMI_exp_mediation.csv",
  row.names = FALSE
)
assign_many(
  BW_GH_exp_MVMR = BW_exp,
  childhood_BMI_GH_exp_MVMR = childhood_BMI_significant_clumped
)

# Combine the two sets of exposure instruments.
assign_many(
  SNP_BW_GH = BW_GH_exp_MVMR[, c("SNP", "pval.exposure"), drop = FALSE],
  SNP_childhood_BMI_GH = childhood_BMI_GH_exp_MVMR[
    , c("SNP", "pval.exposure"), drop = FALSE
  ]
)
SNP_all_GH <- rbind(SNP_BW_GH, SNP_childhood_BMI_GH)
SNP_uni_GH <- SNP_all_GH[!duplicated(SNP_all_GH$SNP), ]
write.csv(SNP_uni_GH, "SNP_uni_GH.csv", row.names = FALSE)
# Clump the combined instrument set again.
exp_MVMR_GH_clumped <- clump_data(
  SNP_uni_GH,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  pop = "EUR"
)
# Removing 2 of 36 variants due to LD with other variants or absence from LD reference panel

# MANUAL DATA CHECKPOINT:
# 1. Add the GRCh37/GRCh38 MarkerName column to SNP_uni_GH.csv.
# 2. Ensure that BW_fetal_data is available and contains the RSID column used below.
# 3. Restart this section after reading the completed SNP_uni_GH.csv if required.
SNP_uni_GH <- read.csv("SNP_uni_GH.csv", header = TRUE)
if (!"MarkerName" %in% names(SNP_uni_GH)) {
  stop(
    "Add the MarkerName column to SNP_uni_GH.csv before continuing Section 5."
  )
}
if (!exists("BW_fetal_data")) {
  stop(
    paste(
      "Create BW_fetal_data before continuing Section 5;",
      "it must contain the RSID column used in the original analysis."
    )
  )
}

# Merge the selected SNPs with the two exposure datasets in one command.
assign_many(
  merBW_GH = merge(
    exp_MVMR_GH_clumped,
    BW_fetal_data,
    by.x = "SNP",
    by.y = "RSID"
  ),
  merchildhood_BMI_GH = merge(
    exp_MVMR_GH_clumped,
    childhood_BMI_data,
    by.x = "SNP",
    by.y = "variant_id"
  )
)
merBW_GH$id.exposure <- "BW"
merchildhood_BMI_GH$id.exposure <- "Childhood_BMI"
write.csv(merBW_GH, "merBW_GH.csv", row.names = FALSE)
write.csv(merchildhood_BMI_GH, "merchildhood_BMI_GH.csv", row.names = FALSE)

# Merge the combined SNP list with the GH outcome data.
merGH <- merge(
  SNP_uni_GH,
  out_GH_meta,
  by.x = "MarkerName",
  by.y = "MarkerName"
)
write.csv(merGH, "merGH.csv", row.names = FALSE)

out_dat_GH <- format_data(
  dat = merGH,
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
  ncontrol_col = NA
)
out_dat_GH$id.outcome <- "GH"
write.csv(out_dat_GH, "out_dat_GH.csv", row.names = FALSE)

# Bind and harmonise the exposure datasets.
expo_dat_GH <- rbind(merBW_GH, merchildhood_BMI_GH)
write.csv(expo_dat_GH, "expo_dat_GH.csv", row.names = FALSE)
expo_dat_GH$effect_allele.exposure <- toupper(
  expo_dat_GH$effect_allele.exposure
)
expo_dat_GH$other_allele.exposure <- toupper(
  expo_dat_GH$other_allele.exposure
)
mvmr_dat_GH <- mv_harmonise_data(expo_dat_GH, out_dat_GH)
save(mvmr_dat_GH, file = "mvmr_dat_GH.Rdata")
load("mvmr_dat_GH.Rdata")

# Estimate the direct multivariable MR effects.
mv_multiple(mvmr_dat_GH)
SummaryStats_childhood_BMI <- cbind(
  mvmr_dat_GH[["outcome_beta"]],
  mvmr_dat_GH[["exposure_beta"]][, 1],
  mvmr_dat_GH[["exposure_beta"]][, 2],
  mvmr_dat_GH[["exposure_se"]][, 1],
  mvmr_dat_GH[["exposure_se"]][, 2],
  mvmr_dat_GH[["outcome_se"]]
)
SummaryStats_childhood_BMI <- data.frame(SummaryStats_childhood_BMI)

MVMR_Input_childhood_BMI <- mr_mvinput(
  bx = cbind(
    SummaryStats_childhood_BMI$X2,
    SummaryStats_childhood_BMI$X3
  ),
  bxse = cbind(
    SummaryStats_childhood_BMI$X4,
    SummaryStats_childhood_BMI$X5
  ),
  by = SummaryStats_childhood_BMI$X1,
  byse = SummaryStats_childhood_BMI$X6
)
ivw <- mr_mvivw(
  MVMR_Input_childhood_BMI,
  model = "default",
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
ivw
egger <- mr_mvegger(
  MVMR_Input_childhood_BMI,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
egger

r_input <- format_mvmr(
  BXGs = mvmr_dat_GH[["exposure_beta"]],
  BYG = mvmr_dat_GH[["outcome_beta"]],
  seBXGs = mvmr_dat_GH[["exposure_se"]],
  seBYG = mvmr_dat_GH[["outcome_se"]],
  RSID = rownames(mvmr_dat_GH[["exposure_beta"]])
)
mvmr_results <- ivw_mvmr(r_input)
pleiotropy_test <- pleiotropy_mvmr(r_input)
Fz <- strength_mvmr(r_input = r_input, gencov = 0)

presso_results_BMI <- mr_presso(
  BetaOutcome = "X1",
  BetaExposure = c("X2", "X3"),
  SdOutcome = "X6",
  SdExposure = c("X4", "X5"),
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = SummaryStats_childhood_BMI,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)
# Outliers identified in the original analysis:
# rs17817449, rs4144829, rs61765651 and rs75844534.
# The following file is the exposure dataset after removing those outliers.
expo_dat_GH_presso <- read.csv("expo_dat_GH_presso.csv", header = TRUE)
mvmr_dat_GH_presso <- mv_harmonise_data(expo_dat_GH_presso, out_dat_GH)
mv_multiple(mvmr_dat_GH_presso)

SummaryStats_childhood_BMI_presso <- cbind(
  mvmr_dat_GH_presso[["outcome_beta"]],
  mvmr_dat_GH_presso[["exposure_beta"]][, 1],
  mvmr_dat_GH_presso[["exposure_beta"]][, 2],
  mvmr_dat_GH_presso[["exposure_se"]][, 1],
  mvmr_dat_GH_presso[["exposure_se"]][, 2],
  mvmr_dat_GH_presso[["outcome_se"]]
)
SummaryStats_childhood_BMI_presso <- data.frame(
  SummaryStats_childhood_BMI_presso
)

MVMR_Input_childhood_BMI_presso <- mr_mvinput(
  bx = cbind(
    SummaryStats_childhood_BMI_presso$X2,
    SummaryStats_childhood_BMI_presso$X3
  ),
  bxse = cbind(
    SummaryStats_childhood_BMI_presso$X4,
    SummaryStats_childhood_BMI_presso$X5
  ),
  by = SummaryStats_childhood_BMI_presso$X1,
  byse = SummaryStats_childhood_BMI_presso$X6
)
ivw <- mr_mvivw(
  MVMR_Input_childhood_BMI_presso,
  model = "default",
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
ivw
egger <- mr_mvegger(
  MVMR_Input_childhood_BMI_presso,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05
)
egger

r_input <- format_mvmr(
  BXGs = mvmr_dat_GH_presso[["exposure_beta"]],
  BYG = mvmr_dat_GH_presso[["outcome_beta"]],
  seBXGs = mvmr_dat_GH_presso[["exposure_se"]],
  seBYG = mvmr_dat_GH_presso[["outcome_se"]],
  RSID = rownames(mvmr_dat_GH_presso[["exposure_beta"]])
)
mvmr_results <- ivw_mvmr(r_input)
qhet_mvmr(r_input)
pleiotropy_test <- pleiotropy_mvmr(r_input)
Fz <- strength_mvmr(r_input = r_input, gencov = 0)


# 6. Delta-method mediation analysis ------------------------------------------------

res <- read.csv("Total_BW_childhood_BMI_GH.csv", header = TRUE)

# Extract all coefficients and standard errors in one command.
assign_many(
  b1 = res[3, "beta"],       # Birth weight to childhood BMI
  se1 = res[3, "se"],
  b2 = res[7, "beta"],       # Childhood BMI to GH, direct effect
  se2 = res[7, "se"],
  total_beta = res[2, "beta"], # Birth weight to GH, total effect
  total_se = res[2, "se"]
)

# Calculate the indirect effect, its uncertainty and mediation proportion.
indirect_effect <- b1 * b2
indirect_se <- sqrt(b1^2 * se2^2 + b2^2 * se1^2)
z_value <- indirect_effect / indirect_se
p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
ci_lower <- indirect_effect - 1.96 * indirect_se
ci_upper <- indirect_effect + 1.96 * indirect_se
mediation_proportion <- indirect_effect / total_beta
prop_se <- abs(mediation_proportion) * sqrt(
  (indirect_se / indirect_effect)^2 + (total_se / total_beta)^2
)
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
