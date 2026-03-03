# Maternal-own-birth-weight-and-HDPs
This repository contains the analysis code for the study investigating the relationship between maternal own birth weight and hypertensive disorders of pregnancy (HDPs), including gestational hypertension (GH) and preeclampsia (PE).
## Overview
We employed an integrated analytical framework combining:
- **Nested case-control studies** 
- **Two-sample Mendelian randomization (MR)** 
- **MR-Clust analysis** 
- **MR mediation analysis** 
## Repository Structure:
```
Maternal-own-birth-weight-and-HDPs/
├── scripts/
│ ├── 01_analysis_UKB.R
│ ├── 02_analysis_SBC.R
│ ├── 03_clump and remove SNPs_BW.R
│ ├── 04_TwoSampleMR_BW to GH and PE.R
│ ├── 05_Cluster_BW_GH.R
│ ├── 06_TwoSampleMR_BW to CBMI.R
│ └── 07_MVMR.R
├── data/
│ ├── DataAvailability.md # Instructions for obtaining GWAS data
└── README.md
```
