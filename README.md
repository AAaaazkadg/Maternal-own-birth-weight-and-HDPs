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
│ ├── 01_nested_case_control.R # Observational analysis (conditional logistic regression, RCS)
│ ├── 02_mr_analysis.R # Two-sample MR (IVW, weighted median, MR-Egger)
│ ├── 03_mr_clust.R # MR-Clust analysis for heterogeneous effects
│ ├── 04_mr_mediation.R # Mediation analysis (childhood BMI)
│ ├── 05_sensitivity.R # Sensitivity analyses (MR-PRESSO, leave-one-out)
│ └── 06_visualization.R # Generate figures (RCS plots, forest plots)
├── data/
│ ├── DataAvailability.md # Instructions for obtaining GWAS data
└── README.md
```
