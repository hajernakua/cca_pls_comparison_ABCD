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
NIH_Data <- read.csv("abcd_tbss01.csv")
CorticalThickness_Data <- read.csv("abcd_smrip10201.csv")
Imaging_Info <- read.csv("abcd_imgincl01.csv")
MRI_Info <- read.csv("abcd_mri01.csv")
HeadInjury_Info <- read.csv("abcd_mx01.csv")
exclusion_criteria <- read.csv("exclusion_nback_sst.csv")

##these two csv files were extracted from DEAP 
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

### only reading in prepared DFs (cleaned but not residualized)
Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores <- read.csv("/Users/hajernakua/Documents/CCA_PLS_Project_V2_2021/Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores.csv")
Baseline_CBCL_CT_SiteFam_Data_NoSib_CT <- read.csv("/Users/hajernakua/Documents/CCA_PLS_Project_V2_2021/Baseline_CBCL_CT_SiteFam_Data_NoSib_CT.csv")

####################################################################################################
#######       CLEANING NIH DATA
####################################################################################################

NIHData[NIHData == ""] <- NA

Baseline_NIHData <- NIHData %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  drop_na() %>%
  select(starts_with("nihtbx"), "subjectkey", "sex", "interview_age") %>%
  select(ends_with("uncorrected"), -nihtbx_totalcomp_uncorrected, -nihtbx_fluidcomp_uncorrected, -nihtbx_cryst_uncorrected)

####################################################################################################
#######       CLEANING CORTICAL THICKNESS DATA
####################################################################################################

##extracting baseline data
Baseline_CorticalThickness_Data <- CorticalThickness_Data %>% 
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total") %>% 
  drop_na(smri_area_cdk_banksstslh)

## merging cortical thickness with descriptive data to be used to exclude participants 
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
NIH_CortThick_merge <- merge(Baseline_NIHData, CortThick_complete_QC_baseline, by = "subjectkey")

####################################################################################################
#######  REMOVING A SIBLING FROM EACH FAMILY 
####################################################################################################

###merging the data frames
names(Baseline_SiteFam_Data)[names(Baseline_SiteFam_Data) == "src_subject_id"] <- "subjectkey"
names(Baseline_FamTwin_Data)[names(Baseline_FamTwin_Data) == "ID"] <- "subjectkey"

Baseline_NIH_CT_Site_Data <- merge(NIH_CortThick_merge, Baseline_SiteFam_Data, by = "subjectkey")
Baseline_NIH_CT_SiteFam_Data <- merge(Baseline_NIH_CT_Site_Data, Baseline_FamTwin_Data, by = "subjectkey")

Baseline_NIH_CT_SiteFam_Data_NoSib <- Baseline_NIH_CT_SiteFam_Data %>%
  group_by(rel_family_id.x) %>%
  filter(row_number(rel_relationship == "sibling") == 1)


####################################################################################################
#######  Residualizing NIH scores  
####################################################################################################

### 1. getting NIH data 
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores <- Baseline_NIH_CT_SiteFam_Data_NoSib[, c('nihtbx_picvocab_uncorrected', 'nihtbx_flanker_uncorrected', 'nihtbx_list_uncorrected', 'nihtbx_cardsort_uncorrected', 
                                                                                                           'nihtbx_pattern_uncorrected', 'nihtbx_picture_uncorrected','nihtbx_reading_uncorrected', 
                                                                                                           'interview_age', 'sex', 'abcd_site', "smri_vol_cdk_total", "mri_info_manufacturersmn")]


### 2. Regressing out age, sex, & abcd_site from the NIH scores

###getting the variables to be the right structure
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$abcd_site <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$abcd_site)
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$sex <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$sex)
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$interview_age <- as.numeric(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$interview_age)
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$mri_info_manufacturersmn <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$mri_info_manufacturersmn)
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$smri_vol_cdk_total <- as.numeric(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$smri_vol_cdk_total)
Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$smri_vol_cdk_total <- Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$smri_vol_cdk_total /1000


Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores[,c(1:7)] <- as.data.frame(lapply(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores[, c(1:7)], function(x) as.numeric(x)))


NIH_final_res <- lapply(Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores, function(x) lm(x ~Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$abcd_site + 
                                                                                          Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$interview_age + 
                                                                                          Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$sex +
                                                                                          Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$mri_info_manufacturersmn +
                                                                                          Baseline_NIH_CT_SiteFam_Data_NoSib_NIHScores$smri_vol_cdk_total)$residuals)


##removing the factors from the list so I can create a data frame 
NIH_final_res <- NIH_final_res[-c(8:12)]
NIH_final_res <- as.data.frame(NIH_final_res)

####################################################################################################
#######  REGRESSING OUT VARIABLES FOR CT SCORES 
####################################################################################################
Baseline_NIH_CT_SiteFam_Data_NoSib_CT <- Baseline_NIH_CT_SiteFam_Data_NoSib %>%
  select(starts_with("smri_thick_cdk_"), "mri_info_manufacturersmn", "smri_vol_cdk_total", "abcd_site", "sex", "interview_age")

###getting the variables to be the right structure
Baseline_NIH_CT_SiteFam_Data_NoSib_CT$abcd_site <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_CT$abcd_site)
Baseline_NIH_CT_SiteFam_Data_NoSib_CT$sex <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_CT$sex)
Baseline_NIH_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn <- as.factor(Baseline_NIH_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn)
Baseline_NIH_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total <- as.numeric(Baseline_NIH_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total)
Baseline_NIH_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total <- Baseline_NIH_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total/1000


Baseline_NIH_CT_SiteFam_Data_NoSib_CT[, c(1:68, 76)] <- as.data.frame(lapply(Baseline_NIH_CT_SiteFam_Data_NoSib_CT[, c(1:68, 76)], function(x) as.numeric(x)))


### 4. Regressing out age, sex, abcd_site, total brain volume, and MRI manufacturer
CTnih_final_res <- lapply(Baseline_NIH_CT_SiteFam_Data_NoSib_CT, function(x) lm(x ~Baseline_NIH_CT_SiteFam_Data_NoSib_CT$abcd_site + 
                                                                                        Baseline_NIH_CT_SiteFam_Data_NoSib_CT$interview_age + 
                                                                                        Baseline_NIH_CT_SiteFam_Data_NoSib_CT$sex +
                                                                                        Baseline_NIH_CT_SiteFam_Data_NoSib_CT$mri_info_manufacturersmn +
                                                                                        Baseline_NIH_CT_SiteFam_Data_NoSib_CT$smri_vol_cdk_total)$residuals)
###removing the variables that are not cortical thickness
CTnih_final_res <- CTnih_final_res[-c(69:76)]
CTnih_final_res <- as.data.frame(CTnih_final_res)


### the NIH_final_res and CTnih_final_res will be used in the second group of CCA and PLS analyses 
