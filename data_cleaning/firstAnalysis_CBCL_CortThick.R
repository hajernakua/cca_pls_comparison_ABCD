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


######## getting all the data frames we will need
CBCL_Data <- read.csv("abcd_cbcls01.csv")
CorticalThickness_Data <- read.csv("abcd_smrip10201.csv")
Imaging_Info <- read.csv("abcd_imgincl01.csv")
MRI_Info <- read.csv("abcd_mri01.csv")
HeadInjury_Info <- read.csv("abcd_mx01.csv")
exclusion_criteria <- read.csv("exclusion_nback_sst.csv")
Site_Family_Data <- read.csv("Site_Family_Data.csv")
Family_Twin_Data <- read.csv("Family_Twin_ABCD_Data.csv")

####################################################################################################
#######   Extracting Baseline Data and variables of interest 
####################################################################################################


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

Baseline_CBCLData <- CBCLData %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  drop_na(cbcl_scr_syn_totprob_r, cbcl_scr_syn_anxdep_r) %>%
  select(ends_with("_r"), "subjectkey", "sex", "interview_age") %>%
  select(-cbcl_scr_syn_internal_r, -cbcl_scr_syn_external_r, -cbcl_scr_syn_totprob_r,
         -starts_with("cbcl_scr_dsm5"))

Baseline_CorticalThickness_Data <- CorticalThickness_Data %>% 
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total") %>% 
  drop_na(smri_area_cdk_banksstslh)


####################################################################################################
#######       CLEANING CORTICAL THICKNESS DATA
####################################################################################################

## 1. merging cortical thickness with descriptive data to be used to exclude participants 
CorticalThicknessData_wMRIData <- merge(Baseline_CorticalThickness_Data, Baseline_MRI_Info[, c('subjectkey','mri_info_manufacturersmn')], by = "subjectkey")
CorticalThicknessData_wMRIData_Var_QC <- merge(CorticalThicknessData_wMRIData, Baseline_Imaging_Info[, c('subjectkey','imgincl_t1w_include')], by = "subjectkey")
CorticalThicknessData_wMRIData_Var_QC <- merge(CorticalThicknessData_wMRIData_Var_QC, Baseline_exclusion_criteria[, c('subjectkey', 'mrif_score', "fsqc_qc", "iqc_t1_ok_ser")], by = "subjectkey")


###removing participants with incidental findings, poor raw T1w scans, poor freesurfer QC, and failed DAIRC:
CortThick_complete_QC_baseline <- CorticalThicknessData_wMRIData_Var_QC %>% 
  filter(!mrif_score %in% c("Consider clinical referral", "Consider immediate clinical referral")) %>%
  filter(iqc_t1_ok_ser != 0) %>%
  filter(fsqc_qc != "reject") %>%
  filter(imgincl_t1w_include == 1) 


####################################################################################################
#######  GETTING ALL THE PARTICIPANTS FOR BOTH MATRICES 
####################################################################################################

## 1. Merge the CBCL and the Cortical Thickness Matrices
CBCL_CortThick_merge <- merge(Baseline_CBCLData, CortThick_complete_QC_baseline, by = "subjectkey")

####################################################################################################
#######  REMOVING A SIBLING FROM EACH FAMILY 
####################################################################################################

###merging the data frames
names(Baseline_SiteFam_Data)[names(Baseline_SiteFam_Data) == "src_subject_id"] <- "subjectkey"
names(Baseline_FamTwin_Data)[names(Baseline_FamTwin_Data) == "ID"] <- "subjectkey"

Baseline_CBCL_CT_Site_Data <- merge(CBCL_CortThick_merge, Baseline_SiteFam_Data, by = "subjectkey")
Baseline_CBCL_CT_SiteFam_Data <- merge(Baseline_CBCL_CT_Site_Data, Baseline_FamTwin_Data, by = "subjectkey")

Baseline_CBCL_CT_SiteFam_Data_NoSib <- Baseline_CBCL_CT_SiteFam_Data %>%
  group_by(rel_family_id.x) %>%
  filter(row_number(rel_relationship == "sibling") == 1)


####################################################################################################
#######  Residualizing CBCL Scores  
####################################################################################################

##covariates: age, sex, head size, site, and MRI manufacturer 

### 1. getting CBCL data 
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores <- Baseline_CBCL_CT_SiteFam_Data_NoSib[, c('cbcl_scr_syn_anxdep_r', 'cbcl_scr_syn_withdep_r', 'cbcl_scr_syn_somatic_r', 'cbcl_scr_syn_social_r', 
                                                                                                           'cbcl_scr_syn_thought_r', 'cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r', 'cbcl_scr_syn_aggressive_r',  
                                                                                                            'cbcl_scr_07_sct_r', 'cbcl_scr_07_ocd_r', 'cbcl_scr_07_stress_r', 'interview_age', 'sex', 'abcd_site', "smri_vol_cdk_total", "mri_info_manufacturersmn")]


###getting the variables to be the right structure
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$abcd_site <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$abcd_site)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$sex <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$sex)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$interview_age <- as.numeric(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$interview_age)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$mri_info_manufacturersmn)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$smri_vol_cdk_total)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$smri_vol_cdk_total <- Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$smri_vol_cdk_total /1000


Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores[,c(1:11)] <- as.data.frame(lapply(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores[, c(1:11)], function(x) as.numeric(x)))


CBCL_final_res <- lapply(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores, function(x) lm(x ~Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$abcd_site + 
                                                                                          Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$interview_age + 
                                                                                          Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$sex +
                                                                                          Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$mri_info_manufacturersmn +
                                                                                          Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores$smri_vol_cdk_total)$residuals)


##removing the factors from the list so I can create a data frame 
CBCL_final_res <- CBCL_final_res[-c(12:16)]
CBCL_final_res <- as.data.frame(CBCL_final_res)

####################################################################################################
#######  REGRESSING OUT VARIABLES FOR CT SCORES 
####################################################################################################


Baseline_CBCL_CT_SiteFam_Data_NoSib_CT <- Baseline_CBCL_CT_SiteFam_Data_NoSib %>%
  select(starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total", "abcd_site", "sex", "interview_age")

Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$rel_family_id.x <- NULL

###getting the variables to be the right structure
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$abcd_site <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$abcd_site)
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$sex <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$sex)
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn <- as.factor(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn)
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total <- as.numeric(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total)
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total <- Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total/1000


Baseline_CBCL_CT_SiteFam_Data_NoSib_CT[, c(1:68, 76)] <- as.data.frame(lapply(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT[, c(1:68, 76)], function(x) as.numeric(x)))


### 4. Regressing out age, sex, abcd_site, total brain volume, and MRI manufacturer

CT_final_res <- lapply(Baseline_CBCL_CT_SiteFam_Data_NoSib_CT, function(x) lm(x ~Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$abcd_site + 
                                                                                        Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$interview_age + 
                                                                                        Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$sex +
                                                                                        Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn +
                                                                                        Baseline_CBCL_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total)$residuals)
###removing the variables that are not cortical thickness
CT_final_res <- CT_final_res[-c(69:73)]
CT_final_res <- as.data.frame(CT_final_res)

### the CBCL_final_res and CT_final_res will be used in the first group of CCA and PLS analyses 
