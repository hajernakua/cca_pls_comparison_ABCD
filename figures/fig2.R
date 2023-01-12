# Ju-Chi Yu (2nd author), Nov. 17th, 2022

# Libraries
library(tidyverse)
library(ggseg)
library(ggsegGlasser)
library(broom)
library(PTCA4CATA)
library(plotly)
library(colorspace)
library(tableone)
library("RColorBrewer")
library(reshape2)
library(ggridges)

#getting the loading data
source("functions/cca_pls_analysis.R")

## Reading the U and V
CCA.u <- X_beta_weights_spearman
CCA.v <- Y_beta_weights_spearman
rownames(CCA.v) <- NULL
rownames(CCA.u) <- c("anxdep", "withdep", "somatic", "social", 
                     "thought", "attention", "rulebreak", "aggressive",
                     "sct", "ocs", "stress")

## Barplots for behaviors

# Load function for the barplot
source("Aim1_withinMethod_generalization/Barplot_LVs.R")

CCA_loadings <- Barplot_LV(CCA.u, column2plot = 1, 
              threshold = 0,
              title = "LV1 Behavioural Loadings in CCA", 
              sort = TRUE,
              color4pos = "#612975", color4neg = "#556b2f",
              font.size = 4, horizontal = FALSE,
              ylim.min = NULL, ylim.max = NULL, # these will be set automatically if NULL
              xaxis.textsize = 10, xaxis.lwd = 1)

##getting the PLS figures 
PLS.u <- PLS_svd_scaled_spearman$u
PLS.v <- PLS_svd_scaled_spearman$v

rownames(PLS.u) <- c("anxdep", "withdep", "somatic", "social", 
                     "thought", "attention", "rulebreak", "aggressive",
                     "sct", "ocs", "stress")

PLS_loadings <- Barplot_LV(PLS.u, column2plot = 1, 
                              threshold = 0,
                              title = "LV1 Behavioural Loadings in PLS", 
                              sort = TRUE,
                              color4pos = "#612975", color4neg = "#556b2f",
                              font.size = 4, horizontal = FALSE,
                              ylim.min = -0.5, ylim.max = 0, # these will be set automatically if NULL
                              xaxis.textsize = 10, xaxis.lwd = 1)+ 
                              theme(axis.text.y = element_text(size=11))

## Reading latent scores for CCA
CCA.scores <- UX_VY_mat_CCA_spearman
colnames(CCA.scores) <- c("Lx.1", "Ly.1")

## scatter pot for Lx and Ly

# Load function for the scatter plot
source("figures/LatentScore_plot.R")

# Example code to plot it
CCA_latentScores <- LxLyplot(CCA.scores[,1], CCA.scores[,2], 
               column2plot.Lx = 1,
               column2plot.Ly = 1,
               title = "Latent scores for CCA \n(19.3% of variance)",
               Name4X = "X",
               Name4Y = "Y",
               textsize = 10,
               col.line = "navy",
               col.points = "royalblue3",
               line.width = 1.5, # line width of the regression line
               alpha.points = 0.2 # transparency of the dots
               )

## Reading latent scores for PLS
PLS.scores <- UX_VY_mat_PLS_spearman
colnames(PLS.scores) <- c("Lx.1", "Ly.1")

PLS_latentScores <- LxLyplot4Hajer(PLS.scores[,1], PLS.scores[,2], 
                                   column2plot.Lx = 1,
                                   column2plot.Ly = 1,
                                   title = "Latent scores for PLS \n(81.6% of covariance)",
                                   Name4X = "X",
                                   Name4Y = "Y",
                                   textsize = 10,
                                   col.line = "navy",
                                   col.points = "royalblue3",
                                   line.width = 1.5, # line width of the regression line
                                   alpha.points = 0.2 # transparency of the dots
)

## brain plot for loadings

# Load function for brainplot
source("figures/BrainPlot.R")

# Example code to plot it
dev.new()
CCA_brainLoadings <- Brainplot(CCA.v, column2plot = 1,
                palette = "RdBu", 
                limits = c(-0.6, 0.6), 
                values = c(0, 0.47, 0.5, 0.53, 1),
                title = "LV1 Brain Loadings for CCA",
                brain.position = position_brain(side ~ hemi))



PLS_brainLoadings <- Brainplot(PLS.v, column2plot = 1,
                                     palette = "RdBu", 
                                     limits = c(-0.6, 0.6), 
                                     values = c(0, 0.47, 0.5, 0.53, 1),
                                     title = "LV1 Brain Loadings for PLS",
                                     brain.position = position_brain(side ~ hemi))

#### getting the train-test data
pls_trainTest <- singVal_trainTest
cca_trainTest <- singVal_trainTest

##melting the split-half DFs
pls_trainTest_melted <- melt(pls_trainTest)
cca_trainTest_melted <- melt(cca_trainTest)

##getting the individual train-test distribution figures
pls_trainTest_fig <- pls_trainTest_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  ggtitle(paste0("PLS Test-set Stability")) +
  labs(y="", x="Predicted Singular Values") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "none")

cca_trainTest_fig <- cca_trainTest_melted %>% 
  ggplot(aes(y=variable,x=value,fill=variable)) +
  geom_density_ridges()+
  theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  ggtitle(paste0("CCA Test-set Stability")) +
  labs(y="", x="Predicted Singular Values") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "none")


###merging it all together to get one figure 
Figure2_plot_final <- gridExtra::grid.arrange(
  grobs = list(CCA_loadings, CCA_brainLoadings, 
               PLS_loadings, PLS_brainLoadings,
               CCA_latentScores, PLS_latentScores,
               pls_trainTest_fig, cca_trainTest_fig),
  widths = c(1, 1, 1, 1, 1, 1,1,1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1,1,2,2,5,5,7,7),
                        c(3,3,4,4,6,6,8,8)))
