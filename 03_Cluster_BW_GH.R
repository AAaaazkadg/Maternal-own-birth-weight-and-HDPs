# 1. Install and load necessary packages --------------------------------------------------------
library(devtools)
library(mrclust)    
library(ggplot2)   
library(dplyr)      
library(readr)     
library(ISLR)

# 2. Read and prepare data ------------------------------------------------------------
# Read harmonised data
harmonised_data <- read.csv("Harmonization_GH.csv", header = TRUE)

# Define abbreviations for exposure and outcome
exposure_abbr <- "BW"  
outcome_abbr <- "GH"   

# 3. Extract variables required for MR-Clust analysis ------------------------------------------------
# Extract necessary variables from harmonised data
mr_clust_data_BW_GH <- data.frame(
  SNP = harmonised_data$SNP,
  betaX_BW = harmonised_data$beta.exposure,    
  betaY_GH = harmonised_data$beta.outcome,     
  seX_BW = harmonised_data$se.exposure,        
  seY_GH = harmonised_data$se.outcome,         
  pvalX_BW = harmonised_data$pval.exposure,    
  pvalY_GH = harmonised_data$pval.outcome      
)

# Calculate Wald ratio estimates and their standard errors
mr_clust_data_BW_GH$theta_BW_GH <- mr_clust_data_BW_GH$betaY_GH / mr_clust_data_BW_GH$betaX_BW  
mr_clust_data_BW_GH$theta_se_BW_GH <- mr_clust_data_BW_GH$seY_GH / abs(mr_clust_data_BW_GH$betaX_BW)  

# 4. Run MR-Clust clustering analysis -----------------------------------------------------
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
head(res_em$results$all,n=24)   
head(res_em$results$best,n=24)  
# 5. Results visualization  ----------------------------------------------------------------
# Create cluster annotation scatter plot
cluster_plot <- res_em$plots$two_stage +
  ggplot2::ggtitle("MR-Clust Analysis: Birth Weight → Gestational Hypertension") +
  ggplot2::xlim(0, max(abs(mr_clust_data_BW_GH$betaX_BW) +2*mr_clust_data_BW_GH$seX_BW)) + 
  ggplot2::xlab("Genetic association with Birth Weight (betaX)") +
  ggplot2::ylab("Genetic association with Gestational Hypertension (betaY)")
# Display plot
cluster_plot

write.csv(res_em$results$best, "cluster_result_BW_GH_2019.csv")

# Implement method version (B): Probability ≥ 0.7
# 1. First use pr_clust to filter data meeting conservative criteria
filtered_data <- pr_clust(
  dta = res_em$results$best,
  prob = 0.7,    # Assignment probability ≥ 0.7
  min_obs = 1    # Cluster member count ≥ 1
)
write.csv(filtered_data,"cluster_result_BW_GH_2019_conservation.csv")

# 2. Obtain indices of filtered variants
keep_indices <- which(res_em$results$best$observation %in% filtered_data$observation)

# 3. Extract filtered data
filtered_bx <- mr_clust_data_BW_GH$betaX_BW[keep_indices]
filtered_by <- mr_clust_data_BW_GH$betaY_GH[keep_indices]
filtered_bxse <- mr_clust_data_BW_GH$seX_BW[keep_indices]
filtered_byse <- mr_clust_data_BW_GH$seY_GH[keep_indices]
filtered_rsid <- mr_clust_data_BW_GH$SNP[keep_indices]

# 4. Plot two_stage_plot for conservative method
conservative_plot <- two_stage_plot(
  res = filtered_data,          
  bx = filtered_bx,             
  by = filtered_by,             
  bxse = filtered_bxse,         
  byse = filtered_byse,         
  obs_names = filtered_rsid 
)

# 5. Add custom styling
conservative_plot <- conservative_plot +
  ggplot2::ggtitle("MR-Clust Analysis: Birth Weight → Gestational Hypertension") +
  ggplot2::xlim(0, max(abs(mr_clust_data_BW_GH$betaX_BW) +2*mr_clust_data_BW_GH$seX_BW)) + 
  ggplot2::xlab("Genetic association with Birth Weight") +
  ggplot2::ylab("Genetic association with Gestational Hypertension") 

# Display plot
print(conservative_plot)

# Define color scheme
# Adjust colors based on number of clusters
my_colors <- c(
  "junk" = "black",        
  "null" = "gray50",         
  "1" = "#E18727FF",          
  "2" = "#0072B5FF"        
)

# Modify original cluster plot
cluster_plot_custom <- res_em$plots$two_stage +
  ggplot2::ggtitle("MR-Clust Analysis: Birth Weight → Gestational Hypertension") +
  ggplot2::xlim(0, max(abs(mr_clust_data_BW_GH$betaX_BW) + 2*mr_clust_data_BW_GH$seX_BW)) + 
  ggplot2::xlab("Genetic association with Birth Weight (betaX)") +
  ggplot2::ylab("Genetic association with Gestational Hypertension (betaY)") +
  # Apply custom colors
  ggplot2::scale_color_manual(values = my_colors) +
  ggplot2::scale_fill_manual(values = my_colors)

print(cluster_plot_custom)

conservative_plot <- conservative_plot +
  ggplot2::ggtitle("MR-Clust Analysis: Birth Weight → Gestational Hypertension") +
  ggplot2::xlim(0, max(abs(mr_clust_data_BW_GH$betaX_BW) +2*mr_clust_data_BW_GH$seX_BW)) + 
  ggplot2::xlab("Genetic association with Birth Weight") +
  ggplot2::ylab("Genetic association with Gestational Hypertension") + 
# Apply custom colors
  ggplot2::scale_color_manual(values = my_colors) +
  ggplot2::scale_fill_manual(values = my_colors)

print(conservative_plot)