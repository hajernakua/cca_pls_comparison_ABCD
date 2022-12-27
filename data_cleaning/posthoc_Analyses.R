###loading any necessary packages 
library(ggplot2)
library(plotly)
library(sparklyr)
library(magrittr)
library(lubridate)
library(dplyr)
library(tidyverse)
library(plyr)
library(standardize)


##reading in the required csv files
CBCL_CT_ALLmerge <- read.csv("CBCL_CT_ALLmerge.csv")
CBCL_Data <- read.csv("abcd_cbcls01.csv")

Site_Family_Data <- read.csv("Site_Family_Data.csv")
Family_Twin_Data <- read.csv("Family_Twin_ABCD_Data.csv")

Baseline_SiteFam_Data <- Site_Family_Data %>% 
  filter(eventname == "baseline_year_1_arm_1")

Baseline_FamTwin_Data <- Family_Twin_Data %>% 
  filter(event_name == "baseline_year_1_arm_1")

Baseline_MRI_Info <- MRI_Info %>% 
  filter(eventname == "baseline_year_1_arm_1")

Baseline_Imaging_Info <- Imaging_Info %>% 
  filter(eventname == "baseline_year_1_arm_1")

Baseline_exclusion_criteria <- exclusion_criteria %>% 
  filter(event_name == "baseline_year_1_arm_1")
names(Baseline_exclusion_criteria)[names(Baseline_exclusion_criteria) == "src_subject_id"] <- "subjectkey"

####################################################################################################
#######       Cleaning CBCL Data
####################################################################################################

Baseline_CBCLData <- CBCLData %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  drop_na(cbcl_scr_syn_totprob_r, cbcl_scr_syn_anxdep_r) %>%
  select(ends_with("_r"), "subjectkey", "sex", "interview_age") %>%
  select(-cbcl_scr_syn_internal_r, -cbcl_scr_syn_external_r, -cbcl_scr_syn_totprob_r,
         -starts_with("cbcl_scr_dsm5"))


####################################################################################################
#######  Removing a sibling from each family 
####################################################################################################

###merging the data frames
names(Baseline_SiteFam_Data)[names(Baseline_SiteFam_Data) == "src_subject_id"] <- "subjectkey"
names(Baseline_FamTwin_Data)[names(Baseline_FamTwin_Data) == "ID"] <- "subjectkey"

Baseline_CBCL_CT_Site_Data <- merge(CBCL_CT_ALLmerge, Baseline_SiteFam_Data, by = "subjectkey")
Baseline_CBCL_CT_SiteFam_Data <- merge(Baseline_CBCL_CT_Site_Data, Baseline_FamTwin_Data, by = "subjectkey")


Baseline_CBCL_CT_SiteFam_Data_NoSib <- Baseline_CBCL_CT_SiteFam_Data %>%
  group_by(rel_family_id.x) %>%
  filter(row_number(rel_relationship == "sibling") == 1)


####################################################################################################
#######  Symptom Endorsement Subset  
####################################################################################################

## removing participants with at least one 0 across their CBCL subscale scores 

Baseline_CBCL_CT_SiteFam_Data_NoSib$cbcl_rowsum <- rowSums(Baseline_CBCL_CT_SiteFam_Data_NoSib[, c(2:12)])
Baseline_CBCL_CT_symptomEndorsement <- filter(Baseline_CBCL_CT_SiteFam_Data_NoSib,cbcl_rowsum > 11)


####################################################################################################
#######  Subclinical Subset 
####################################################################################################

##merging the total CBCL T-score from the uncleaned baseline CBCL with the cleaned participant data 

Baseline_CBCL_CT_SiteFam_Data_NoSib_wTotalCBCL <- merge(Baseline_CBCL_CT_SiteFam_Data_NoSib, Baseline_CBCLData[, c("subjectkey", "cbcl_scr_syn_totprob_t")], by = "subjectkey")
Baseline_CBCL_CT_subclinical <- filter(Baseline_CBCL_CT_SiteFam_Data_NoSib_wTotalCBCL, cbcl_scr_syn_totprob_t > 60)



####################################################################################################
#######  Residualizing Subclinical Subset 
####################################################################################################

######### regressing the CBCL data 
Baseline_CBCL_CT_subclinical_cbclScores <- Baseline_CBCL_CT_subclinical[, c('cbcl_scr_syn_anxdep_r', 'cbcl_scr_syn_withdep_r', 'cbcl_scr_syn_somatic_r', 'cbcl_scr_syn_social_r', 
                                                                                                           'cbcl_scr_syn_thought_r', 'cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r', 'cbcl_scr_syn_aggressive_r',  
                                                                                                            'cbcl_scr_07_sct_r', 'cbcl_scr_07_ocd_r', 'cbcl_scr_07_stress_r', 'interview_age', 'sex', 'abcd_site', "smri_vol_cdk_total", "mri_info_manufacturersmn")]

#### found a problem with the CBCL dataset in which some people have an empty cell for the CBCL scores so I need to remove those
Baseline_CBCL_CT_subclinical_cbclScores <- Baseline_CBCL_CT_subclinical_cbclScores %>% drop_na(cbcl_scr_syn_anxdep_r)

###getting the variables to be the right structure
Baseline_CBCL_CT_subclinical_cbclScores$abcd_site <- as.factor(Baseline_CBCL_CT_subclinical_cbclScores$abcd_site)
Baseline_CBCL_CT_subclinical_cbclScores$sex <- as.factor(Baseline_CBCL_CT_subclinical_cbclScores$sex)
Baseline_CBCL_CT_subclinical_cbclScores$interview_age <- as.numeric(Baseline_CBCL_CT_subclinical_cbclScores$interview_age)
Baseline_CBCL_CT_subclinical_cbclScores$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_subclinical_cbclScores$mri_info_manufacturersmn)
Baseline_CBCL_CT_subclinical_cbclScores$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_subclinical_cbclScores$smri_vol_cdk_total)
Baseline_CBCL_CT_subclinical_cbclScores$smri_vol_cdk_total <- Baseline_CBCL_CT_subclinical_cbclScores$smri_vol_cdk_total /1000

Baseline_CBCL_CT_subclinical_cbclScores[,c(1:11)] <- as.data.frame(lapply(Baseline_CBCL_CT_subclinical_cbclScores[, c(1:11)], function(x) as.numeric(x)))

CBCL_final_subclinical_res <- lapply(Baseline_CBCL_CT_subclinical_cbclScores, function(x) lm(x ~Baseline_CBCL_CT_subclinical_cbclScores$abcd_site + 
                                                                                                Baseline_CBCL_CT_subclinical_cbclScores$interview_age + 
                                                                                                Baseline_CBCL_CT_subclinical_cbclScores$sex +
                                                                                                Baseline_CBCL_CT_subclinical_cbclScores$mri_info_manufacturersmn +
                                                                                                Baseline_CBCL_CT_subclinical_cbclScores$smri_vol_cdk_total)$residuals)


##removing the factors from the list so I can create a data frame 
CBCL_final_subclinical_res <- CBCL_final_subclinical_res[-c(12:19)]
CBCL_final_subclinical_res <- as.data.frame(CBCL_final_subclinical_res)


######### regressing the cortical thickness data
Baseline_CBCL_CT_subclinical_CT <- Baseline_CBCL_CT_subclinical %>%
  select(starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total", "abcd_site", "sex", "interview_age")

###getting the variables to be the right structure
Baseline_CBCL_CT_subclinical_CT$abcd_site <- as.factor(Baseline_CBCL_CT_subclinical_CT$abcd_site)
Baseline_CBCL_CT_subclinical_CT$sex <- as.factor(Baseline_CBCL_CT_subclinical_CT$sex)
Baseline_CBCL_CT_subclinical_CT$interview_age <- as.numeric(Baseline_CBCL_CT_subclinical_CT$interview_age)
Baseline_CBCL_CT_subclinical_CT$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_subclinical_CT$mri_info_manufacturersmn)
Baseline_CBCL_CT_subclinical_CT$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_subclinical_CT$smri_vol_cdk_total)
Baseline_CBCL_CT_subclinical_CT$smri_vol_cdk_total <- Baseline_CBCL_CT_subclinical_CT$smri_vol_cdk_total /1000



Baseline_CBCL_CT_subclinical_CT[, c(1:68)] <- as.data.frame(lapply(Baseline_CBCL_CT_subclinical_CT[, c(1:68)], function(x) as.numeric(x)))


### 4. Regressing out age, sex, abcd_site, total brain volume, and MRI manufacturer

CT_final_subclinical_res <- lapply(Baseline_CBCL_CT_subclinical_CT, function(x) lm(x ~Baseline_CBCL_CT_subclinical_CT$abcd_site + 
                                                                                 Baseline_CBCL_CT_subclinical_CT$interview_age + 
                                                                                 Baseline_CBCL_CT_subclinical_CT$sex +
                                                                                 Baseline_CBCL_CT_subclinical_CT$mri_info_manufacturersmn +
                                                                                 Baseline_CBCL_CT_subclinical_CT$mri_vol_cdk_total)$residuals)
###removing the variables that are not cortical thickness
CT_final_subclinical_res <- CT_final_subclinical_res[-c(69:79)]
CT_final_subclinical_res <- as.data.frame(CT_final_subclinical_res)

####################################################################################################
#######  Residualizing the Symptom Endorsement Subset
####################################################################################################

### 1. getting CBCL data 
Baseline_CBCL_CT_symptomEndorsement_cbclScores <- Baseline_CBCL_CT_symptomEndorsement[, c('cbcl_scr_syn_anxdep_r', 'cbcl_scr_syn_withdep_r', 'cbcl_scr_syn_somatic_r', 'cbcl_scr_syn_social_r', 
                                                                                                        'cbcl_scr_syn_thought_r', 'cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r', 'cbcl_scr_syn_aggressive_r',  
                                                                                                        'cbcl_scr_07_sct_r', 'cbcl_scr_07_ocd_r', 'cbcl_scr_07_stress_r', 'interview_age', 'sex', 'abcd_site', "smri_vol_cdk_total", "mri_info_manufacturersmn")]

#### found a problem with the CBCL dataset in which some people have an empty cell for the CBCL scores so I need to remove those
Baseline_CBCL_CT_symptomEndorsement_cbclScores <- Baseline_CBCL_CT_symptomEndorsement_cbclScores %>% drop_na(cbcl_scr_syn_anxdep_r)


###getting the variables to be the right structure
Baseline_CBCL_CT_symptomEndorsement_cbclScores$abcd_site <- as.factor(Baseline_CBCL_CT_symptomEndorsement_cbclScores$abcd_site)
Baseline_CBCL_CT_symptomEndorsement_cbclScores$sex <- as.factor(Baseline_CBCL_CT_symptomEndorsement_cbclScores$sex)
Baseline_CBCL_CT_symptomEndorsement_cbclScores$interview_age <- as.numeric(Baseline_CBCL_CT_symptomEndorsement_cbclScores$interview_age)
Baseline_CBCL_CT_symptomEndorsement_cbclScores$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_symptomEndorsement_cbclScores$mri_info_manufacturersmn)
Baseline_CBCL_CT_symptomEndorsement_cbclScores$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_symptomEndorsement_cbclScores$smri_vol_cdk_total)
Baseline_CBCL_CT_symptomEndorsement_cbclScores$smri_vol_cdk_total <- Baseline_CBCL_CT_symptomEndorsement_cbclScores$smri_vol_cdk_total /1000

Baseline_CBCL_CT_symptomEndorsement_cbclScores[,c(1:11)] <- as.data.frame(lapply(Baseline_CBCL_CT_symptomEndorsement_cbclScores[, c(1:11)], function(x) as.numeric(x)))


CBCL_final_symptomEndorsement_res <- lapply(Baseline_CBCL_CT_symptomEndorsement_cbclScores, function(x) lm(x ~Baseline_CBCL_CT_symptomEndorsement_cbclScores$abcd_site + 
                                                                                                        Baseline_CBCL_CT_symptomEndorsement_cbclScores$interview_age + 
                                                                                                        Baseline_CBCL_CT_symptomEndorsement_cbclScores$sex +
                                                                                                        Baseline_CBCL_CT_symptomEndorsement_cbclScores$mri_info_manufacturersmn +
                                                                                                        Baseline_CBCL_CT_symptomEndorsement_cbclScores$smri_vol_cdk_total)$residuals)


##removing the factors from the list so I can create a data frame 
CBCL_final_symptomEndorsement_res <- CBCL_final_symptomEndorsement_res[-c(12:19)]
CBCL_final_symptomEndorsement_res <- as.data.frame(CBCL_final_symptomEndorsement_res)


#### residualizing cortical thickness variables 
Baseline_CBCL_CT_symptomEndorsement_CT <- Baseline_CBCL_CT_symptomEndorsement %>%
  select(starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total", "abcd_site", "sex", "interview_age")

###getting the variables to be the right structure
Baseline_CBCL_CT_symptomEndorsement_CT$abcd_site <- as.factor(Baseline_CBCL_CT_symptomEndorsement_CT$abcd_site)
Baseline_CBCL_CT_symptomEndorsement_CT$sex <- as.factor(Baseline_CBCL_CT_symptomEndorsement_CT$sex)
Baseline_CBCL_CT_symptomEndorsement_CT$interview_age <- as.numeric(Baseline_CBCL_CT_symptomEndorsement_CT$interview_age)
Baseline_CBCL_CT_symptomEndorsement_CT$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_symptomEndorsement_CT$mri_info_manufacturersmn)
Baseline_CBCL_CT_symptomEndorsement_CT$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_symptomEndorsement_CT$smri_vol_cdk_total)
Baseline_CBCL_CT_symptomEndorsement_CT$smri_vol_cdk_total <- Baseline_CBCL_CT_symptomEndorsement_CT$smri_vol_cdk_total /1000


Baseline_CBCL_CT_symptomEndorsement_CT[, c(1:68)] <- as.data.frame(lapply(Baseline_CBCL_CT_symptomEndorsement_CT[, c(1:68)], function(x) as.numeric(x)))


### 4. Regressing out age, sex, abcd_site, total brain volume, and MRI manufacturer

CT_final_symptomEndorsement_res <- lapply(Baseline_CBCL_CT_symptomEndorsement_CT, function(x) lm(x ~Baseline_CBCL_CT_symptomEndorsement_CT$abcd_site + 
                                                                                              Baseline_CBCL_CT_symptomEndorsement_CT$interview_age + 
                                                                                              Baseline_CBCL_CT_symptomEndorsement_CT$sex +
                                                                                              Baseline_CBCL_CT_symptomEndorsement_CT$mri_info_manufacturersmn +
                                                                                              Baseline_CBCL_CT_symptomEndorsement_CT$smri_vol_cdk_total)$residuals)
###removing the variables that are not cortical thickness
CT_final_symptomEndorsement_res <- CT_final_symptomEndorsement_res[-c(69:79)]
CT_final_symptomEndorsement_res <- as.data.frame(CT_final_symptomEndorsement_res)
