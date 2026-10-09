# UK Biobank nested case-control analysis
# Maternal birth weight in relation to preeclampsia and gestational hypertension.

# 0. Load all required R packages ---------------------------------------------

required_packages <- c(
  "readxl", "dplyr", "stringr", "MatchIt", "tableone", "survival",
  "Hmisc", "rms", "ggplot2", "grid", "checkmate", "abind", "forestplot"
)
invisible(lapply(required_packages, library, character.only = TRUE))

# 1. Import and merge the UK Biobank source data -------------------------------
# Import participant-level data from two UK Biobank extracts.
data1<-read_excel("data1_participant.xlsx")
data2<-read_excel("data2_participant.xlsx") 
# Merge the two extracts by the unique UK Biobank participant identifier.
data<-merge(data1,data2,by="eid",all=TRUE)

# 2. Clean and derive the analysis variables ----------------------------------
# Columns 150-307 contain diagnosis-date fields whose order corresponds to the ICD-10 codes in field 41270.
names(data)[150:307]<-paste0("p41280_a",101:258)
# Treat the UK Biobank response "Not known" as missing across character fields.
data<-data%>%
  mutate(across(where(is.character),~ifelse(.=="Not known",NA,.)))
# Use the first available self-reported birth-weight value across assessment instances (UK Biobank field 20022).
data<-data%>%
  mutate(BirthWeight=coalesce(p20022_i0,p20022_i1,p20022_i2))
summary(data$BirthWeight)
# Replace missing birth weight with 3.32 kg, the mean value used for imputation in the analysis.
data$BirthWeight[is.na(data$BirthWeight)]<-3.32
# Convert birth weight from kilograms to grams.
data$BirthWeight<-data$BirthWeight*1000
summary(data$BirthWeight)
data<-data%>%
  rename(ICDcode=p41270)
# Z37.0 and Z37.1 indicate singleton live birth and singleton stillbirth, respectively.
# Exclude participants without either code.
pregnancy<-"Z37\\.(0|1)"
data<-data%>%
  filter(str_detect(ICDcode,pregnancy))
# Collapse the detailed UK Biobank ethnicity responses into five categories.
get_race_category<-function(x){
  case_when(
    x%in%c("Any other white background","British","Irish","White")~"white",
    x%in%c("white and black caribbean","White and Black African","White and Asian","Any other mixed background")~"mix",
    x%in%c("Chinese","Bangladeshi","Asian or Asian British","Pakistani","Any other Asian background")~"Asian",
    x%in%c("African","Caribbean","Black or Black British","Any other Black background")~"black",
    x=="Other ethnic group"~"other",
    TRUE~NA_character_
  )
}
# Use the first available ethnicity response across assessment instances.
data<-data%>%
  mutate(
    race=case_when(
      !(p21000_i0%in%c("Do not know","Prefer not to answer"))~get_race_category(p21000_i0),
      !(p21000_i1%in%c("Do not know","Prefer not to answer"))~get_race_category(p21000_i1),
      !(p21000_i2%in%c("Do not know","Prefer not to answer"))~get_race_category(p21000_i2),
      !(p21000_i3%in%c("Do not know","Prefer not to answer"))~get_race_category(p21000_i3),
      TRUE~"none"
    )
  )
table(data$race,useNA="always")
# Assign missing or unclassified ethnicity to the White category.
data$race<-ifelse(is.na(data$race),"white",data$race)
table(data$race,useNA="always")

# 3. Define the PE and GH case cohorts -----------------------------------------
# Define preeclampsia using ICD-10 O14 and gestational hypertension using ICD-10 O13.
data_PE<-data%>%
  filter(grepl("O14",ICDcode))
data_GH<-data%>%
  filter(grepl("O13",ICDcode))
# Identify participants with both PE and GH records.
data_GH$match<-data_GH$eid%in%data_PE$eid
data_both<-data_GH[data_GH$match==TRUE,]

# 3.1 Exclude prior GH from the PE cohort --------------------------------------
# For participants with both O13 and O14 records, compare diagnosis dates to identify a history of hypertensive disorders of pregnancy before the index diagnosis.
# In the PE cohort, calculate the interval as the index PE diagnosis date (O14) minus the earlier GH diagnosis date (O13), using it as a proxy for the interpregnancy interval.
# Classify participants with an interval greater than 182 days as having a history of GH and flag them for exclusion from the PE cohort.
# The position of each code in ICDcode corresponds to the suffix of its diagnosis-date variable p41280_a*.
target_icd_PE<-"O14"
unwanted_icd_PE<-"O13"
days_threshold_PE<-182 
keep_rows_PE<-logical(nrow(data_both))
for(i in 1:nrow(data_both)){  
  icds<-unlist(strsplit(data_both$ICDcode[i],"\\|"))  
  target_pos_PE<-which(startsWith(icds,target_icd_PE))  
  unwanted_pos_PE<-which(startsWith(icds,unwanted_icd_PE))  
  if(length(target_pos_PE)>0&&length(unwanted_pos_PE)>0){  
    target_date_vars<-paste0("p41280_a",target_pos_PE-1)  
    unwanted_date_vars<-paste0("p41280_a",unwanted_pos_PE-1)  
    target_dates_PE<-as.Date(as.character(data_both[i,target_date_vars]),format="%Y-%m-%d")  
    unwanted_dates_PE<-as.Date(as.character(data_both[i,unwanted_date_vars]),format="%Y-%m-%d")  
    target_dates_PE[is.na(target_dates_PE)]<-as.Date(NA)  
    unwanted_dates_PE[is.na(unwanted_dates_PE)]<-as.Date(NA)  
    date_diffs<-outer(target_dates_PE,unwanted_dates_PE,"-")  
    if (any(date_diffs>days_threshold_PE,na.rm=TRUE)){  
      keep_rows_PE[i]<-TRUE 
    }  
  }  
}  
filtered_data_PE<-data_both[keep_rows_PE,]  

# 3.2 Exclude prior PE from the GH cohort --------------------------------------
# In the GH cohort, calculate the interval as the index GH diagnosis date (O13) minus the earlier PE diagnosis date (O14).
# Classify participants with an interval greater than 224 days as having a history of PE and flag them for exclusion from the GH cohort.
target_icd_GH<-"O13"
unwanted_icd_GH<-"O14"
days_threshold_GH<-224 
keep_rows_GH<-logical(nrow(data_both))
for(i in 1:nrow(data_both)){  
  icds<-unlist(strsplit(data_both$ICDcode[i],"\\|"))  
  target_pos_GH<-which(startsWith(icds,target_icd_GH))  
  unwanted_pos_GH<-which(startsWith(icds,unwanted_icd_GH))  
  if(length(target_pos_GH)>0&&length(unwanted_pos_GH)>0){  
    target_date_vars<-paste0("p41280_a",target_pos_GH-1)  
    unwanted_date_vars<-paste0("p41280_a",unwanted_pos_GH-1)  
    target_dates_GH<-as.Date(as.character(data_both[i,target_date_vars]),format="%Y-%m-%d")  
    unwanted_dates_GH<-as.Date(as.character(data_both[i,unwanted_date_vars]),format="%Y-%m-%d")  
    target_dates_GH[is.na(target_dates_GH)]<-as.Date(NA)  
    unwanted_dates_GH[is.na(unwanted_dates_GH)]<-as.Date(NA)  
    date_diffs<-outer(target_dates_GH,unwanted_dates_GH,"-")  
    if (any(date_diffs>days_threshold_GH,na.rm=TRUE)){  
      keep_rows_GH[i]<-TRUE 
    }  
  }  
}  
filtered_data_GH<-data_both[keep_rows_GH,]  

# Exclude participants classified as having a previous PE or GH diagnosis.
data_PE$match<-data_PE$eid%in%filtered_data_PE$eid
data_PE<-data_PE[data_PE$match==FALSE,]
data_GH$match<-data_GH$eid%in%filtered_data_GH$eid
data_GH<-data_GH[data_GH$match==FALSE,]

# 3.3 Apply the final GH eligibility list --------------------------------------
# Apply the derived eligibility list containing participants included in the GH cohort after the preceding cohort-selection procedures.
data_02<-read.csv("~/intermediate_cohort_after_6_2.csv")
data_02_gh<-data_02[data_02$Group=="GH",]
data_GH$match<-data_GH$eid%in%data_02_gh$Participant.ID
table(data_GH$match)
data_GH<-data_GH[data_GH$match==TRUE,]
# Code cases as 1 for the matched case-control analyses.
data_PE$case<-1
data_GH$case<-1

# 4. Define controls and derive the matching variables -------------------------
# Define eligible controls as participants without ICD-10 O10-O16, covering hypertensive disorders in pregnancy.
data_control<-data%>% 
  filter(!grepl("O10|O11|O12|O13|O14|O15|O16",ICDcode))
data_control$case<-0

# 4.1 Derive matching age for PE cases -----------------------------------------
# Match four controls to each case on age at the recorded singleton delivery.
find_z37_position<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions)
  }
}
data_PE$Z37_position<-sapply(data_PE$ICDcode,find_z37_position)
find_z37_position1<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions[1])
  }
}
data_PE$Z37_position1<-sapply(data_PE$ICDcode,find_z37_position1)
# Extract the diagnosis date associated with the first Z37 code for each PE case.
# The corresponding date-field suffix is the zero-based code position.
data_PE$age<-NA
for(i in 1:nrow(data_PE)){
  pos<-data_PE$Z37_position1[i]
  if(!is.na(pos)){
    var_name<-paste0("p41280_a",pos-1) 
    if(var_name%in%names(data_PE)){
      data_PE$age[i]<-data_PE[[var_name]][i]
    }else{
      data_PE$age[i]<-NA
      warning(paste("Variable",var_name,"not found in data_control for row",i))
    }
  }
}
# Convert diagnosis dates stored as seconds since 1970-01-01 to calendar years.
# Subtract the year of birth (UK Biobank field 34) to estimate age at delivery.
data_PE$age<-data_PE$age/31536000
data_PE$age<-data_PE$age+1970
data_PE$age<-round(data_PE$age)
table(data_PE$age,useNA="always")
data_PE$age<-data_PE$age-data_PE$p34

# 4.2 Derive matching age for GH cases -----------------------------------------
# Repeat the same first-Z37 date and age derivation for GH cases.
find_z37_position<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions)
  }
}
data_GH$Z37_position<-sapply(data_GH$ICDcode,find_z37_position)
find_z37_position1<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions[1])
  }
}
data_GH$Z37_position1<-sapply(data_GH$ICDcode,find_z37_position1)
data_GH$age<-NA
for(i in 1:nrow(data_GH)){
  pos<-data_GH$Z37_position1[i]
  if(!is.na(pos)){
    var_name<-paste0("p41280_a",pos-1) 
    if(var_name%in%names(data_GH)){
      data_GH$age[i]<-data_GH[[var_name]][i]
    }else{
      data_GH$age[i]<-NA
      warning(paste("Variable",var_name,"not found in data_control for row",i))
    }
  }
}
data_GH$age<-data_GH$age/31536000
data_GH$age<-data_GH$age+1970
data_GH$age<-round(data_GH$age)
table(data_GH$age,useNA="always")
data_GH$age<-data_GH$age-data_GH$p34

# 4.3 Derive matching age for eligible controls --------------------------------
# Repeat the same first-Z37 date and age derivation for eligible controls.
find_z37_position<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions)
  }
}
data_control$Z37_position<-sapply(data_control$ICDcode,find_z37_position)
find_z37_position1<-function(icd_string){
  codes<-unlist(strsplit(icd_string,"\\|"))
  positions<-grep("Z37",codes)
  if(length(positions)== 0){
    return(NA)
  }else{
    return(positions[1])
  }
}
data_control$Z37_position1<-sapply(data_control$ICDcode,find_z37_position1)
data_control$age<-NA
for(i in 1:nrow(data_control)){
  pos<-data_control$Z37_position1[i]
  if(!is.na(pos)){
    var_name<-paste0("p41280_a",pos-1) 
    if(var_name%in%names(data_control)){
      data_control$age[i]<-data_control[[var_name]][i]
    }else{
      data_control$age[i]<-NA
      warning(paste("Variable",var_name,"not found in data_control for row",i))
    }
  }
}
data_control$age<-data_control$age/31536000
data_control$age<-data_control$age+1970
data_control$age<-round(data_control$age)
table(data_control$age,useNA="always")
na_rows<-which(is.na(data_control$age))
print(na_rows)
data_control<-data_control[!is.na(data_control$age),]
data_control$age<-data_control$age-data_control$p34
# Remove temporary overlap indicators before combining cases and controls.
data_PE$match<-NULL
data_GH$match<-NULL
data_pe<-rbind(data_PE,data_control)
data_gh<-rbind(data_GH,data_control)

# 4.4 Perform 1:4 nearest-neighbour matching -----------------------------------
# Perform separate 1:4 nearest-neighbour matching for the PE and GH analyses.
# Use age as the matching variable and a caliper of 2 on the matching scale.
data_pe_match<-matchit(
  case~age,
  data=data_pe,
  method="nearest",
  ratio=4,
  caliper=2,
  replace=FALSE 
)
data_pe_matched<-match.data(data_pe_match)
summary(data_pe_matched)
# Extract matched PE case records for cohort accounting.
data_pe_check<-data_pe_matched[data_pe_matched$case==1,]
# Compare age distributions after matching.
t.test(age~case,data=data_pe_matched)
data_gh_match<-matchit(
  case~age,
  data=data_gh,
  method="nearest",
  ratio=4,
  caliper=2,
  replace=FALSE 
)
data_gh_matched<-match.data(data_gh_match)
summary(data_gh_matched)
# Extract matched GH case records for cohort accounting.
data_gh_check<-data_gh_matched[data_gh_matched$case==1,]
# Compare age distributions after matching.
t.test(age~case,data=data_gh_matched)

# 5. Descriptive analysis ------------------------------------------------------
# Create datasets for unadjusted models and models adjusted for ethnicity.
data_pe0<-data_pe_matched%>%
  select(eid,BirthWeight,case,age,subclass)
data_gh0<-data_gh_matched%>%
  select(eid,BirthWeight,case,age,subclass)
data_pe1<-data_pe_matched%>%
  select(eid,BirthWeight,case,age,race,subclass)
data_gh1<-data_gh_matched%>%
  select(eid,BirthWeight,case,age,race,subclass)
vars<-c("age","BirthWeight")
catVars<-NULL
matched_table<-CreateTableOne(
  vars=vars,
  strata="case",
  data=data_pe0,
  factorVars=catVars
)
print(matched_table,smd=TRUE,showAllLevels=TRUE)
matched_table<-CreateTableOne(
  vars=vars,
  strata="case",
  data=data_gh0,
  factorVars=catVars
)
print(matched_table,smd=TRUE,showAllLevels=TRUE)

# 6. Conditional logistic regression analysis --------------------------------
# Evaluate the linear association between maternal birth weight and hypertensive disorders of pregnancy while accounting for matched subclasses.
# PE model with birth weight entered in grams.
logistics_pe0<-clogit(
  case~BirthWeight+strata(subclass),
  data=data_pe0
)
summary(logistics_pe0)
# Rescale birth weight so that the effect estimate represents a 500-g increase.
data_pe0$BirthWeight_500g<-data_pe0$BirthWeight/500
logistics_pe0000<-clogit(
  case~BirthWeight_500g+strata(subclass),
  data=data_pe0
)
summary(logistics_pe0000)
# Ethnicity-adjusted PE model per 500-g increase in birth weight.
data_pe1$BirthWeight_500g<-data_pe1$BirthWeight/500
logistics_pe1000<-clogit(
  case~BirthWeight_500g+race+strata(subclass),
  data=data_pe1
)
summary(logistics_pe1000)
# GH model with birth weight entered in grams.
logistics_gh0<-clogit(
  case~BirthWeight+strata(subclass),
  data=data_gh0
)
summary(logistics_gh0)
data_gh0$BirthWeight_500g<-data_gh0$BirthWeight/500
logistics_gh0000<-clogit(
  case~BirthWeight_500g+strata(subclass),
  data=data_gh0
)
summary(logistics_gh0000) 
# Ethnicity-adjusted GH model per 500-g increase in birth weight.
data_gh1$BirthWeight_500g<-data_gh1$BirthWeight/500
logistics_gh1000<-clogit(
  case~BirthWeight_500g+race+strata(subclass),
  data=data_gh1
)
summary(logistics_gh1000)

# 7. Restricted cubic spline analysis -----------------------------------------
# Explore potential non-linear associations between maternal birth weight and each hypertensive disorder. Models include ethnicity and matched subclass.
# Inspect the observed range and interquartile range before fitting PE splines.
range(data_pe1$BirthWeight)
quantile(data_pe1$BirthWeight,probs=c(0.25,0.75))
# datadist stores predictor distributions used by rms for prediction and plotting.
dd<-datadist(data_pe1)
options(datadist="dd")
# Fit and plot PE models using three, four, and five spline knots.
rcs_pe3<-lrm(case~rcs(BirthWeight,3)+race+strat(subclass),data=data_pe1)
summary(rcs_pe3)
anova(rcs_pe3)
pred_pe<-Predict(rcs_pe3,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_pe,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Preeclampsia Risk(3)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
rcs_pe4<-lrm(case~rcs(BirthWeight,4)+race+strat(subclass),data=data_pe1)
summary(rcs_pe4)
anova(rcs_pe4)
pred_pe<-Predict(rcs_pe4,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_pe,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Preeclampsia Risk(4)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
rcs_pe5<-lrm(case~rcs(BirthWeight,5)+race+strat(subclass),data=data_pe1)
summary(rcs_pe5)
anova(rcs_pe5)
pred_pe<-Predict(rcs_pe5,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_pe,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Preeclampsia Risk(5)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
# Compare model fit across knot specifications.
# A lower AIC indicates better relative fit, while lrtest compares nested specifications.
AIC(rcs_pe3, rcs_pe4, rcs_pe5)
lrtest(rcs_pe3,rcs_pe4)
lrtest(rcs_pe4,rcs_pe5)
# Repeat the spline analysis for GH using the same sequence of knot numbers.
options(datadist="dd")
dd<-datadist(data_gh1)
range(data_gh1$BirthWeight)
quantile(data_gh1$BirthWeight,probs=c(0.25,0.75))
rcs_gh3<-lrm(case~rcs(BirthWeight,3)+race+strat(subclass),data=data_gh1)
summary(rcs_gh3)
anova(rcs_gh3)
pred_gh<-Predict(rcs_gh3,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_gh,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Gestational hypertension Risk(3)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
rcs_gh4<-lrm(case~rcs(BirthWeight,4)+race+strat(subclass),data=data_gh1)
summary(rcs_gh4)
anova(rcs_gh4)
pred_gh<-Predict(rcs_gh4,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_gh,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Gestational hypertension Risk(4)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
rcs_gh5<-lrm(case~rcs(BirthWeight,5)+race+strat(subclass),data=data_gh1)
summary(rcs_gh5)
anova(rcs_gh5)
pred_gh<-Predict(rcs_gh5,BirthWeight,ref.zero=TRUE,fun=exp)
ggplot(pred_gh,aes(x=BirthWeight,y=yhat))+ 
  geom_line(size=1.5,color="#AD002AFF")+ 
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.2,fill="#925E9FFF")+
  geom_hline(yintercept=1,linetype="dashed")+ 
  labs(x="Mother's Birth Weight (g)", 
       y="Odds Ratio (95% CI)", 
       title="Association between Mother's Birth Weight and Gestational hypertension Risk(5)")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
AIC(rcs_gh3, rcs_gh4, rcs_gh5)
lrtest(rcs_gh3,rcs_gh4)
lrtest(rcs_gh4,rcs_gh5)

# 8. Categorical birth-weight analysis ----------------------------------------
# Divide maternal birth weight into five clinically interpretable categories:
# <2500, 2500-2999, 3000-3499, 3500-3999, and >=4000 g.
# Use the 3000-3499 g group as the reference category.
data_preeclampsia<-data_pe0%>%
  mutate(BirthWeight=case_when(
    BirthWeight<2500~"low",
    BirthWeight>=2500&BirthWeight<=2999~"lownormal",
    BirthWeight>=3000&BirthWeight<=3499~"normal",
    BirthWeight>=3500&BirthWeight<=3999~"highnormal",
    BirthWeight>=4000~"high"
  ))
table(data_preeclampsia$BirthWeight)
data_preeclampsia$BirthWeight<-factor(data_preeclampsia$BirthWeight,ordered=FALSE)
data_preeclampsia$BirthWeight<-relevel(data_preeclampsia$BirthWeight,ref="normal")
logistics_preeclampsia<-clogit(
  case~BirthWeight+strata(subclass),
  data=data_preeclampsia
)
summary(logistics_preeclampsia)
# Repeat the PE categorical model with adjustment for ethnicity.
data_preeclampsia1<-data_pe1%>%
  mutate(BirthWeight=case_when(
    BirthWeight<2500~"low",
    BirthWeight>=2500&BirthWeight<=2999~"lownormal",
    BirthWeight>=3000&BirthWeight<=3499~"normal",
    BirthWeight>=3500&BirthWeight<=3999~"highnormal",
    BirthWeight>=4000~"high"
  ))
table(data_preeclampsia1$BirthWeight)
data_preeclampsia1$BirthWeight<-factor(data_preeclampsia1$BirthWeight,ordered=FALSE)
data_preeclampsia1$BirthWeight<-relevel(data_preeclampsia1$BirthWeight,ref="normal")
logistics_preeclampsia1<-clogit(
  case~BirthWeight+race+strata(subclass),
  data=data_preeclampsia1
)
summary(logistics_preeclampsia1)
# Fit the corresponding unadjusted categorical model for GH.
data_gestationalhypertension<-data_gh0%>%
  mutate(BirthWeight=case_when(
    BirthWeight<2500~"low",
    BirthWeight>=2500&BirthWeight<=2999~"lownormal",
    BirthWeight>=3000&BirthWeight<=3499~"normal",
    BirthWeight>=3500&BirthWeight<=3999~"highnormal",
    BirthWeight>=4000~"high"
  ))
table(data_gestationalhypertension$BirthWeight)
data_gestationalhypertension$BirthWeight<-factor(data_gestationalhypertension$BirthWeight,ordered=FALSE)
data_gestationalhypertension$BirthWeight<-relevel(data_gestationalhypertension$BirthWeight,ref="normal")
logistics_gestationalhypertension<-clogit(
  case~BirthWeight+strata(subclass),
  data=data_gestationalhypertension
)
summary(logistics_gestationalhypertension)
# Repeat the GH categorical model with adjustment for ethnicity.
data_gestationalhypertension1<-data_gh1%>%
  mutate(BirthWeight=case_when(
    BirthWeight<2500~"low",
    BirthWeight>=2500&BirthWeight<=2999~"lownormal",
    BirthWeight>=3000&BirthWeight<=3499~"normal",
    BirthWeight>=3500&BirthWeight<=3999~"highnormal",
    BirthWeight>=4000~"high"
  ))
table(data_gestationalhypertension1$BirthWeight)
data_gestationalhypertension1$BirthWeight<-factor(data_gestationalhypertension1$BirthWeight,ordered=FALSE)
data_gestationalhypertension1$BirthWeight<-relevel(data_gestationalhypertension1$BirthWeight,ref="normal")
logistics_gestationalhypertension1<-clogit(
  case~BirthWeight+race+strata(subclass),
  data=data_gestationalhypertension1
)
summary(logistics_gestationalhypertension1)

# 9. Forest plots --------------------------------------------------------------
# Display odds ratios and 95% confidence intervals from the categorical analyses.
# The values below correspond to the final categorical model results.
box_colors <- c("black", "#0099B4FF", "#AD002AFF", "black", "#925E9FFF", "#42B540FF")
styles <- fpShapesGp(
  box = lapply(box_colors, function(col) gpar(fill = col, col = col, shape = "circle")),
  lines = lapply(box_colors, function(col) gpar(col = col, lwd = 2))
)

# Draw the four coloured legend lines with one function call.
legend_line_x <- list(
  c(0.500, 0.530),
  c(0.582, 0.612),
  c(0.730, 0.760),
  c(0.883, 0.913)
)
legend_line_colors <- c("#0099B4FF", "#AD002AFF", "#925E9FFF", "#42B540FF")
draw_legend_lines <- function(y) {
  invisible(Map(
    function(x, col) {
      grid.lines(
        x = unit(x, "npc"),
        y = unit(y, "npc"),
        gp = gpar(col = col)
      )
    },
    legend_line_x,
    legend_line_colors
  ))
}

data_forest_pe <- data.frame(
  Group = c("low","low-normal","normal","high-normal","high"),
  n = c(115,272,999,380,124),
  OddsRatio = c(1.61,1.57,1,0.90,0.90),
  LowerCI = c(1.04,1.14,NA,0.66,0.54),
  UpperCI = c(2.50,2.14,NA,1.23,1.49),
  P_value = c(0.03,0.009,NA,0.50,0.68),
  stringsAsFactors = FALSE
)
data_forest_pe$P_formatted <- ifelse(
  data_forest_pe$P_value < 0.01,
  "<0.01",
  sprintf("%.2f", data_forest_pe$P_value)
)
data_forest_pe$LowerCI[data_forest_pe$Group == "normal"] <- NA
data_forest_pe$UpperCI[data_forest_pe$Group == "normal"] <- NA
data_forest_pe$P_formatted[data_forest_pe$Group == "normal"] <- "-"
data_forest_pe$OddsRatio_formatted <- ifelse(
  data_forest_pe$Group == "normal",
  "1",
  sprintf("%.2f", data_forest_pe$OddsRatio)
)
tabletext_pe <- cbind(
  c("Group", data_forest_pe$Group),
  c("n", sprintf("%d", data_forest_pe$n)),
  c("OR (95% CI)", ifelse(data_forest_pe$Group == "normal", 
                          "1 (Reference)",
                          paste0(sprintf("%.2f", data_forest_pe$OddsRatio), 
                                 " (", sprintf("%.2f", data_forest_pe$LowerCI), 
                                 "-", sprintf("%.2f", data_forest_pe$UpperCI), ")"))),
  c("P.value", ifelse(data_forest_pe$Group == "normal", " - ", data_forest_pe$P_formatted))
)
forestplot(tabletext_pe,
           mean = c(NA, data_forest_pe$OddsRatio),
           lower = c(NA, data_forest_pe$LowerCI),
           upper = c(NA, data_forest_pe$UpperCI),
           is.summary = c(TRUE, rep(FALSE, nrow(data_forest_pe))),
           graph.pos = 5,
           title = "",
           xlab = "Odds Ratio (with 95% CI)",
           txt_gp = fpTxtGp(label = gpar(cex=1.3),
                            ticks = gpar(cex=0.8),
                            xlab = gpar(cex = 1.0),
                            title = gpar(cex = 1.9),
                            summary = gpar(cex=1.3)),
           lineheight = unit(12,"mm"),
           colgap = unit(6,"mm"),
           lwd.ci = 2,
           ci.vertices = TRUE,
           box.size = 0.3,
           zero = 1,
           xticks = seq(0, 3, by = 0.5),
           graphwidth = unit(8,"cm"),
           clip = c(0.1, 4),
           mar = unit(c(5, 1, 7, 1), "lines"),
           shapes_gp = styles
)
grid.legend(
  labels = c("low","low-normal","high-normal", "high"),
  pch = 15,
  ncol = 4,
  gp = gpar(
    col  = c("#0099B4FF", "#AD002AFF", "#925E9FFF", "#42B540FF"),
    fill = c("#0099B4FF", "#AD002AFF", "#925E9FFF", "#42B540FF"),
    cex  = 1
  ),
  vp = viewport(
    x = 0.73,
    y = 0.2,
    just = c("center", "center")
  )
)
draw_legend_lines(0.2)
data_forest_gh <- data.frame(
  Group = c("low","low-normal","normal","high-normal","high"),
  n = c(167,400,1429,537,192),
  OddsRatio = c(1.08,1.31,1,0.88,1.35),
  LowerCI = c(0.73,1.00,NA,0.68,0.95),
  UpperCI = c(1.61,1.72,NA,1.15,1.93),
  P_value = c(0.70,0.05,NA,0.36,0.10),
  stringsAsFactors = FALSE
)
data_forest_gh$P_formatted <- ifelse(
  data_forest_gh$P_value < 0.01,
  "<0.01",
  sprintf("%.2f", data_forest_gh$P_value)
)
data_forest_gh$LowerCI[data_forest_gh$Group == "normal"] <- NA
data_forest_gh$UpperCI[data_forest_gh$Group == "normal"] <- NA
data_forest_gh$P_formatted[data_forest_gh$Group == "normal"] <- "-"
data_forest_gh$OddsRatio_formatted <- ifelse(
  data_forest_gh$Group == "normal",
  "1",
  sprintf("%.2f", data_forest_gh$OddsRatio)
)
tabletext_gh <- cbind(
  c("Group", data_forest_gh$Group),
  c("n", sprintf("%d", data_forest_gh$n)),
  c("OR (95% CI)", ifelse(data_forest_gh$Group == "normal", 
                          "1 (Reference)",
                          paste0(sprintf("%.2f", data_forest_gh$OddsRatio), 
                                 " (", sprintf("%.2f", data_forest_gh$LowerCI), 
                                 "-", sprintf("%.2f", data_forest_gh$UpperCI), ")"))),
  c("P.value", ifelse(data_forest_gh$Group == "normal", " - ", data_forest_gh$P_formatted))
)
forestplot(tabletext_gh,
           mean = c(NA, data_forest_gh$OddsRatio),
           lower = c(NA, data_forest_gh$LowerCI),
           upper = c(NA, data_forest_gh$UpperCI),
           is.summary = c(TRUE, rep(FALSE, nrow(data_forest_gh))),
           graph.pos = 5,
           title = "",
           xlab = "Odds Ratio (with 95% CI)",
           txt_gp = fpTxtGp(label = gpar(cex=1.3),
                            ticks = gpar(cex=0.8),
                            xlab = gpar(cex = 1.0),
                            title = gpar(cex = 1.9),
                            summary = gpar(cex=1.3)),
           lineheight = unit(12,"mm"),
           colgap = unit(6,"mm"),
           lwd.ci = 2,
           ci.vertices = TRUE,
           box.size = 0.3,
           zero = 1,
           xticks = seq(0, 3, by = 0.5),
           graphwidth = unit(8,"cm"),
           clip = c(0.1, 4),
           mar = unit(c(5, 1, 7, 1), "lines"),
           shapes_gp = styles
)
grid.legend(
  labels = c("low","low-normal","high-normal", "high"),
  pch = 15,
  ncol = 4,
  gp = gpar(
    col  = c("#0099B4FF", "#AD002AFF", "#925E9FFF", "#42B540FF"),
    fill = c("#0099B4FF", "#AD002AFF", "#925E9FFF", "#42B540FF"),
    cex  = 1
  ),
  vp = viewport(
    x = 0.73,
    y = 0.2,
    just = c("center", "center")
  )
)
draw_legend_lines(0.3)
