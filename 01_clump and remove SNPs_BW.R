library(TwoSampleMR)
BW_significiant_SEM <- read.csv("Significiant_BW_2019_SEM.csv",header = T)
BW_significiant_SEM_exposure<- subset(BW_significiant_SEM,Fetal.P.value<5e-08)
write.csv(BW_significiant_SEM_exposure,"BW_significiant_SEM_exposure.csv")
# Remove rs138715366 rare variants (minor allele frequency <0.01)
ex1 <- which(BW_significiant_SEM_exposure$SNP=="rs138715366")
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[-ex1,]
# Remove rs11042596 in the range of imprinted genes
ex2 <- which(BW_significiant_SEM_exposure$SNP=="rs11042596")
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[-ex2,]
# Remove rs1801253 Unclassified effect
ex3 <- which(BW_significiant_SEM_exposure$SNP=="rs1801253")
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[-ex3,]
# Remove rs10872678 MTA effect
ex4 <- which(BW_significiant_SEM_exposure$SNP=="rs10872678")
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[-ex4,]
# Remove rs560887 MNTA effect
ex5 <- which(BW_significiant_SEM_exposure$SNP=="rs560887")
BW_significiant_SEM_exposure <- BW_significiant_SEM_exposure[-ex5,]
write.csv(BW_significiant_SEM_exposure,"BW_significiant_SEM_exposure.csv")
BW_significiant_SEM_exposure <- read_exposure_data("BW_significiant_SEM_exposure.csv" ,
                                                   sep=",",
                                                   snp_col = "SNP",
                                                   beta_col="Fetal.Beta..SDs.",
                                                   se_col="Fetal.SE",
                                                   effect_allele_col = "Effect.allele",
                                                   other_allele_col = "Other.allele",
                                                   pval_col = "Fetal.P.value")
BW_significiant_SEM_exposure$exposure <- "Birth weight"

BW_significiant_SEM_exposure_clumped <- clump_data(BW_significiant_SEM_exposure,
                                                   clump_kb = 10000,
                                                   clump_r2 = 0.001, 
                                                   clump_p1 = 1,
                                                   clump_p2 = 1,pop="EUR")
write.csv(BW_significiant_SEM_exposure_clumped,"BW_significiant_SEM_exposure_clumped.csv")