#creating more descriptive names 
beh_cbcl <- CBCL_final_res
brain_cortThick <- CT_final_res

##getting the z-scores
brain_cortThick_res_z <- scale(brain_cortThick, center = TRUE, scale = TRUE)
beh_cbcl_res_z <- scale(beh_cbcl, center = TRUE, scale = TRUE)

####################################################################################################
#######  Preparing the correlation matrices 
####################################################################################################

###preparing both the Pearson and Spearman correlation matrices

beh_cbcl_res_z_cor <- cor(beh_cbcl_res_z)
beh_cbcl_res_z_cor_spearman <- cor(beh_cbcl_res_z, method = "spearman")

brain_cortThick_res_z_cor <- cor(brain_cortThick_res_z)
brain_cortThick_res_z_cor_spearman <- cor(brain_cortThick_res_z, method = "spearman")

#Rxy matrices
cbcl_CT_Rxy <- cor(beh_cbcl_res_z, brain_cortThick_res_z)
cbcl_CT_Rxy_spearman <- cor(beh_cbcl_res_z, brain_cortThick_res_z, method = "spearman")

Rxy_plot <- plot(cbcl_CT_Rxy, key = list(side=3))
Rxy_plot_spearman <- plot(cbcl_CT_Rxy_spearman, key = list(side=3))

#Omega matrices
cbcl_CT_omega <- t(solve(chol(beh_cbcl_res_z_cor))) %*% cbcl_CT_Rxy %*% solve(chol(brain_cortThick_res_z_cor))
cbcl_CT_omega_spearman <- t(solve(chol(beh_cbcl_res_z_cor_spearman))) %*% cbcl_CT_Rxy_spearman %*% solve(chol(brain_cortThick_res_z_cor_spearman))


###visualizing the correlation matrices
library(plot.matrix)
Rxy_plot <- plot(cbcl_CT_Rxy, key = list(side=3))
Rxy_plot_spearman <- plot(cbcl_CT_Rxy_spearman, key = list(side=3))

Omega_plot <- plot(cbcl_CT_omega, key = list(side=3))
Omega_plot_spearman <- plot(cbcl_CT_omega_spearman, key = list(side=3))


##################################################################
#### RUNNING THE CCA
##################################################################

##Performing the decomposition on the cross-product matrices
cca_svd_scaled <- svd(cbcl_CT_omega)
cca_svd_scaled_spearman <- svd(cbcl_CT_omega_spearman)

##getting the variance explained of each singular value
cca_svd_scaled$d^2/sum(cca_svd_scaled$d^2)
cca_svd_scaled_spearman$d^2/sum(cca_svd_scaled_spearman$d^2)

###re-weighting the singular vectors  
X_beta_weights <- solve(chol(beh_cbcl_res_z_cor)) %*% cca_svd_scaled$u 
X_beta_weights_spearman <- solve(chol(beh_cbcl_res_z_cor_spearman)) %*% cca_svd_scaled_spearman$u 

Y_beta_weights <- solve(chol(brain_cortThick_res_z_cor)) %*% cca_svd_scaled$v 
Y_beta_weights_spearman <- solve(chol(brain_cortThick_res_z_cor_spearman)) %*% cca_svd_scaled_spearman$v 

###getting the structural coefficients 
Sxx_StrucCoeff <- beh_cbcl_res_z_cor %*% X_beta_weights
Sxx_StrucCoeff_spearman <- beh_cbcl_res_z_cor_spearman %*% X_beta_weights_spearman

Syy_StrucCoeff <- brain_cortThick_res_z_cor %*% Y_beta_weights
Syy_StrucCoeff_spearman <- brain_cortThick_res_z_cor_spearman %*% Y_beta_weights_spearman


###running the parametric null hypothesis test for CCA (this test does not exist for PLS)
library(CCA)
library(CCP)
CCA_analysis <- cc(beh_cbcl_res_z, brain_cortThick_res_z)
analysis.sig <- p.asym(CCA_analysis$cor,9191,11, 68,tstat = "Wilks")


##################################################################
#### RUNNING THE PLS
##################################################################

##Performing the decomposition on the cross-product matrices
PLS_svd_scaled <- svd(cbcl_CT_Rxy)
PLS_svd_scaled_spearman <- svd(cbcl_CT_Rxy_spearman)

##getting the variance explained of each singular value
PLS_svd_scaled$d^2/sum(PLS_svd_scaled$d^2)
PLS_svd_scaled_spearman$d^2/sum(PLS_svd_scaled_spearman$d^2)


##################################################################
#### Getting the Latent Scores - PLS
##################################################################

#Standardizing the V singular values by the singular values
S_diag_PLS_spearman <- as.data.frame(diag(PLS_svd_scaled_spearman$d))
S_diag_PLS_spearman <- as.matrix(S_diag_PLS_spearman)
V_S_PLS_spearman <- PLS_svd_scaled_spearman$v %*% S_diag_PLS_spearman

#projecting the standardized V singular values by the residualized and normalized brain matrix 
V_Ymat_PLS_spearman <- brain_cortThick_res_z %*% V_S_PLS_spearman
V_Ymat_PLS_spearman <- as.data.frame(V_Ymat_PLS_spearman)

#Standardizing the U singular values by the singular values
U_S_PLS_spearman <- PLS_svd_scaled_spearman$u %*% S_diag_PLS_spearman

#projecting the standardized U singular values by the residualized and normalized behaviour matrix
U_Xmat_PLS_spearman <- beh_cbcl_res_z %*% U_S_PLS_spearman
U_Xmat_PLS_spearman <- as.data.frame(U_Xmat_PLS_spearman)

##creating a DF with the latent scores of the first LV
latentScores_LV1_mat_PLS_spearman <- cbind.data.frame(U_Xmat_PLS_spearman[, c(1)], V_Ymat_PLS_spearman[, c(1)])

##visualizing the relationship
Scatterplot_latentScores_LV1_PLS_spearman <- ggplot(data = latentScores_LV1_mat_PLS_spearman, aes(x = U_Xmat_PLS_spearman[, c(1)], y = V_Ymat_PLS_spearman[, c(1)])) +
  geom_smooth(method = lm, aes(), size = 1) +
  geom_point(colour = "lightpink", size = 0.2) +
  theme(legend.position = "none")  + 
  theme_classic() +
  #theme(text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 20))


##################################################################
#### Getting the Latent Scores - CCA
##################################################################

#Standardizing the V singular values by the singular values
S_diag_cca_spearman <- as.data.frame(diag(cca_svd_scaled_spearman$d))
S_diag_cca_spearman <- as.matrix(S_diag_cca_spearman)
V_S_cca_spearman <- cca_svd_scaled_spearman$v %*% S_diag_cca_spearman

#projecting the standardized V singular values by the residualized and normalized brain matrix 
V_Ymat_cca_spearman <- brain_cortThick_res_z %*% V_S_cca_spearman
V_Ymat_cca_spearman <- as.data.frame(V_Ymat_cca_spearman)

#Standardizing the U singular values by the singular values
U_S_cca_spearman <- cca_svd_scaled_spearman$u %*% S_diag_cca_spearman

#projecting the standardized U singular values by the residualized and normalized behaviour matrix
U_Xmat_cca_spearman <- beh_cbcl_res_z %*% U_S_cca_spearman
U_Xmat_cca_spearman <- as.data.frame(U_Xmat_cca_spearman)

##creating a DF with the latent scores of the first LV
latentScores_LV1_mat_cca_spearman <- cbind.data.frame(U_Xmat_cca_spearman[, c(1)], V_Ymat_cca_spearman[, c(1)])

##visualizing the relationship
Scatterplot_latentScores_LV1_cca_spearman <- ggplot(data = latentScores_LV1_mat_cca_spearman, aes(x = U_Xmat_cca_spearman[, c(1)], y = V_Ymat_cca_spearman[, c(1)])) +
  geom_smooth(method = lm, aes(), size = 1) +
  geom_point(colour = "lightpink", size = 0.2) +
  theme(legend.position = "none")  + 
  theme_classic() +
  #theme(text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 20))

