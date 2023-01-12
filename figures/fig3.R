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
CCA.nih.u <- X_beta_weights_nih_spearman
CCA.nih.v <- Y_beta_weights_nih_spearman
rownames(CCA.nih.v) <- NULL
rownames(CCA.nih.u) <- c("picvocab", "flanker", "list", "cardsort", 
                         "pattern", "picture", "reading")

## Barplots for behaviors
# Load function for the barplot
source("Aim1_withinMethod_generalization/Barplot_LVs.R")

CCA.nih_loadings <- Barplot_LV(CCA.nih.u, column2plot = 1, 
              threshold = 0,
              title = "LV1 Behavioural Loadings in CCA", 
              sort = TRUE,
              color4pos = "#612975", color4neg = "#556b2f",
              font.size = 4, horizontal = FALSE,
              ylim.min = NULL, ylim.max = NULL, # these will be set automatically if NULL
              xaxis.textsize = 10, xaxis.lwd = 1)

##getting the PLS figures 
PLS.nih.u <- PLS_svd_scaled_nih_spearman$u
PLS.nih.v <- PLS_svd_scaled_nih_spearman$v
rownames(PLS.nih.u) <- c("picvocab", "flanker", "list", "cardsort", 
                         "pattern", "picture", "reading")

PLS.nih_loadings <- Barplot_LV(PLS.nih.u, column2plot = 1, 
                              threshold = 0,
                              title = "LV1 Behavioural Loadings in PLS", 
                              sort = TRUE,
                              color4pos = "#612975", color4neg = "#556b2f",
                              font.size = 4, horizontal = FALSE,
                              ylim.min = -0.5, ylim.max = 0, # these will be set automatically if NULL
                              xaxis.textsize = 10, xaxis.lwd = 1)+ 
                              theme(axis.text.y = element_text(size=11))

## Reading latent scores for CCA
CCA.nih.scores <- UX_VY_mat_CCA_nih_spearman
colnames(CCA.nih.scores) <- c("Lx.1", "Ly.1")

## scatter pot for Lx and Ly

# Load function for the scatter plot
source("figures/LatentScore_plot.R")

# Example code to plot it
CCA.nih_latentScores <- LxLyplot(CCA.nih.scores[,1], CCA.nih.scores[,2], 
               column2plot.Lx = 1,
               column2plot.Ly = 1,
               title = "Latent scores for CCA \n(41.6% of variance)",
               Name4X = "X",
               Name4Y = "Y",
               textsize = 10,
               col.line = "navy",
               col.points = "royalblue3",
               line.width = 1.5, # line width of the regression line
               alpha.points = 0.2 # transparency of the dots
               )

## Reading latent scores for PLS
PLS.nih.scores <- UX_VY_mat_PLS_nih_spearman
colnames(PLS.nih.scores) <- c("Lx.1", "Ly.1")

PLS.nih_latentScores <- LxLyplot(PLS.nih.scores[,1], PLS.nih.scores[,2], 
                                   column2plot.Lx = 1,
                                   column2plot.Ly = 1,
                                   title = "Latent scores for PLS \n(75.5% of covariance)",
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
CCA.nih_brainLoadings <- Brainplot(CCA.nih.v, column2plot = 1,
                palette = "RdBu", 
                limits = c(-0.6, 0.6), 
                values = c(0, 0.47, 0.5, 0.53, 1),
                title = "LV1 Brain Loadings for CCA",
                brain.position = position_brain(side ~ hemi))

PLS.nih_brainLoadings <- Brainplot(PLS.nih.v, column2plot = 1,
                                     palette = "RdBu", 
                                     limits = c(-0.6, 0.6), 
                                     values = c(0, 0.47, 0.5, 0.53, 1),
                                     title = "LV1 Brain Loadings for PLS",
                                     brain.position = position_brain(side ~ hemi))

#### getting the train-test data
pls.nih_trainTest <- singVal_nih_trainTest
cca.nih_trainTest <- singVal_nih_trainTest

##melting the split-half DFs
pls.nih_trainTest_melted <- melt(pls.nih_trainTest)
cca.nih_trainTest_melted <- melt(cca.nih_trainTest)

##getting the individual train-test distribution figures
pls.nih_trainTest_fig <- pls.nih_trainTest_melted %>% 
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

cca.nih_trainTest_fig <- cca.nih_trainTest_melted %>% 
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
Figure3_plot_final <- gridExtra::grid.arrange(
  grobs = list(CCA.nih_loadings, CCA.nih_brainLoadings, 
               PLS.nih_loadings, PLS.nih_brainLoadings,
               CCA.nih_latentScores, PLS.nih_latentScores,
               pls.nih_trainTest_fig, cca.nih_trainTest_fig),
  widths = c(1, 1, 1, 1, 1, 1,1,1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1,1,2,2,5,5,7,7),
                        c(3,3,4,4,6,6,8,8)))
