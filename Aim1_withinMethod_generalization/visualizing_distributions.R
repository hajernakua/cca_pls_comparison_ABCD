##this script plots the distributions of the various datasets

####################################################################################################
#######  Distributions of residualized phenotypic measures
####################################################################################################

######## main analysis - CBCL scores 

#renaming the columns
colnames(CBCL_final_res) <- c("anxiety/depression", "withdrawal/depression", "somatic symptoms",
                              "social problems", "thought problems", "attention problems", "rule breaking",
                              "aggressive behaviour", "sluggish-cognitive-tempo", "OCD symptoms", "stress symptoms")

library(reshape2)
CBCL_final_res_melted <- melt(CBCL_final_res)
CBCL_final_res_melted_plot <- CBCL_final_res_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  theme(axis.text = element_text(size = 20))


######## post-hoc analysis - subclinical subset - CBCL scores 
colnames(CBCL_final_subclinical_res) <- c("anxiety/depression", "withdrawal/depression", "somatic symptoms",
                              "social problems", "thought problems", "attention problems", "rule breaking",
                              "aggressive behaviour", "sluggish-cognitive-tempo", "OCD symptoms", "stress symptoms")


CBCL_final_subclinical_res_melted <- melt(CBCL_final_subclinical_res)
CBCL_final_subclinical_res_melted_plot <- CBCL_final_subclinical_res_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  #geom_histogram(aes(y=..density.., fill= variable), bins = 60) +
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  #theme(text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 20))


######## post-hoc analysis - symptom endorsement subset - CBCL scores
colnames(CBCL_final_symptomEndorsment_res) <- c("anxiety/depression", "withdrawal/depression", "somatic symptoms",
                                           "social problems", "thought problems", "attention problems", "rule breaking",
                                           "aggressive behaviour", "sluggish-cognitive-tempo", "OCD symptoms", "stress symptoms")


CBCL_final_symptomEndorsment_res_melted <- melt(CBCL_final_symptomEndorsment_res)
CBCL_final_symptomEndorsment_res_melted_plot <- CBCL_final_symptomEndorsment_res_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  #geom_histogram(aes(y=..density.., fill= variable), bins = 60) +
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  #theme(text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 20))

######## post-hoc analysis - main analysis - NIH scores
colnames(NIH_final_res) <- c("picture vocabulary task", "flanker task", "list sorting task", "card sorting task",
                             "pattern attention task", "picture sequence task", "oral reading task")

NIH_final_res_melted <- melt(NIH_final_res)
NIH_final_res_melted_plot <- NIH_final_res_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  #geom_histogram(aes(y=..density.., fill= variable), bins = 60) +
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  #theme(text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 20))


####################################################################################################
#######  Distributions of split-half resampling analysis
####################################################################################################
library(ggridges)
### this is an example of one of the functions used to plot the distribution of correlations when conducting the split-half resampling

SingVecV_splitHalf_melted <- melt(SingVecV_splitHalf)
SingVecV_splitHalf_distribution_figure <- SingVecV_splitHalf_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  theme(axis.text = element_text(size = 20))

tiff("SingVecV_splitHalf_distribution_figure.tiff", width = 24, height = 20, units = "cm", res = 200)
print(SingVecV_splitHalf_distribution_figure)
dev.off()


####################################################################################################
#######  Distributions of train-test resampling analysis
####################################################################################################
### this is an example of one of the functions used to plot the distribution of correlations when conducting the train-test resampling

singVal_trainTest_melted <- melt(singVal_trainTest)
singVal_trainTest_distribution_figure <- singVal_trainTest_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  theme(axis.text = element_text(size = 20))

tiff("singVal_trainTest_distribution_figure.tiff", width = 24, height = 20, units = "cm", res = 200)
print(singVal_trainTest_distribution_figure)
dev.off()
