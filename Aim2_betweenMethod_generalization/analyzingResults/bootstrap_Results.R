#### this script will analyze the bootstrap singular vectors to see which variables are stable in each LV

####################################################################################################
#######      1. SIGN FLIPPING THE SINGULAR VECTORS
####################################################################################################

#1a. performing the sign flipping 
sign_flip_function <- function(bootstrap_singVec, empirical_singVec, numboot) {
  
  sign_flipped_Umat <- vector("list", numboot)
  sign_flipped_Vmat <- vector("list", numboot)
  
  for(idx in 1:numboot) {
    ##multiply each bootstrapped U matrix with the original/empirical U matrix
    boot_singVec_discov <- lapply(1:numboot, function(i) t(bootstrap_singVec$U_mat[[i]]) %*% empirical_singVec$u)
    
    #create a matrix which identifies the negative diagonals 
    sign_flips <- lapply(1:numboot, function(i) diag(boot_singVec_discov[[i]]) / abs(diag(boot_singVec_discov[[i]])))
    sign_mat <- matrix(diag(sign_flips[[idx]]), ncol=11)
    
    #multiplies the corresponding vectors from the U and V matrices by 1 or -1 (if the vector needs to be flipped, they're multiplied by -1 which is stored in the sign_mat matrix)
    flipped_Umat <- bootstrap_singVec$U_mat[[idx]] %*% sign_mat
    flipped_Vmat <- bootstrap_singVec$V_mat[[idx]] %*% sign_mat
    
    sign_flipped_Umat[[idx]] <- flipped_Umat
    sign_flipped_Vmat[[idx]] <- flipped_Vmat
  }
  return(list("sign_flipped_Umat" = sign_flipped_Umat, "sign_flipped_Vmat" = sign_flipped_Vmat))
}

bootstrap_singVec_cca <- sign_flip_function(bootstrap_CCA_singvectors_spearman, cca_svd_scaled_spearman, 1000)



####################################################################################################
  #######      2. GETTING CONFIDENCE INTERVALS FOR BEHAVIOUR
####################################################################################################


##getting the values for LV1
SingVect1_UMat_full <- function(X, numboot) {
  
  X <- X
  numboot <- numboot
  
  SingVect1_BehCog_Umat_CCA <- vector("list", numboot)
  
  for(idx in 1:11) {
    SingVect1_BehCog_Umat_CCA[[idx]] <- sapply(1:numboot, function(i) unlist(X[[i]][idx,1]))
  }
  return(SingVect1_BehCog_Umat_CCA)
}

SingVect1_Umat_CCA_spearman <- SingVect1_UMat_full(bootstrap_singVec_CCA_spearman$sign_flipped_Umat, 1000)


####### getting the confidence intervals 
SingVect1_Umat_CCA_spearman_CI_full <- sapply(1:11, function(i) unlist(SingVect1_Umat_CCA_spearman[[i]]))
SingVect1_Umat_CCA_spearman_CI_full <- as.data.frame(SingVect1_Umat_CCA_spearman_CI_full)

SingVect1_Umat_CCA_spearman_CIval_full <- as.data.frame(apply(SingVect1_Umat_CCA_spearman_CI_full, 2, function(x) quantile(x, c(.05, .95)))) 
SingVect1_Umat_CCA_spearman_CIval_full


##getting the values for LV2
SingVect2_UMat_full <- function(X, numboot) {
  
  X <- X
  numboot <- numboot
  
  SingVect2_BehCog_Umat_PLS <- vector("list", numboot)
  
  for(idx in 1:11) {
    SingVect2_BehCog_Umat_PLS[[idx]] <- sapply(1:numboot, function(i) unlist(X[[i]][idx,2]))
  }
  return(SingVect2_BehCog_Umat_PLS)
}

SingVect2_BehCog_Umat_full <- SingVect2_UMat_full(bootstrap_singVec_PLS_pearson$sign_flipped_Umat, 1000)

####### getting the confidence intervals 
SingVect2_BehCog_Umat_CI_full <- sapply(1:11, function(i) unlist(SingVect2_BehCog_Umat_full[[i]]))
SingVect2_BehCog_Umat_CI_full <- as.data.frame(SingVect2_BehCog_Umat_CI_full)


SingVect2_BehCog_Umat_CIval_full <- as.data.frame(apply(SingVect2_BehCog_Umat_CI_full, 2, function(x) quantile(x, c(.05, .95)))) 
SingVect2_BehCog_Umat_CIval_full



####################################################################################################
#######      3. THE BRAIN LATENT VARIABLE - SPEARMAN
####################################################################################################


##getting the values for LV1
SingVect1_VMat_full <- function(X, numboot) {
  
  X <- X
  numboot <- numboot
  
  SingVect1_BehCog_Vmat_CCA <- vector("list", numboot)
  
  for(idx in 1:68) {
    SingVect1_BehCog_Vmat_CCA[[idx]] <- sapply(1:numboot, function(i) unlist(X[[i]][idx,1]))
  }
  return(SingVect1_BehCog_Vmat_CCA)
}

SingVect1_Vmat_CCA_spearman <- SingVect1_VMat_full(bootstrap_singVec_CCA_spearman$sign_flipped_Vmat, 1000)


####### getting the confidence intervals 
SingVect1_Vmat_CCA_spearman_CI_full <- sapply(1:68, function(i) unlist(SingVect1_Vmat_CCA_spearman[[i]]))
SingVect1_Vmat_CCA_spearman_CI_full <- as.data.frame(SingVect1_Vmat_CCA_spearman_CI_full)


SingVect1_Vmat_CCA_spearman_CIval_full <- as.data.frame(apply(SingVect1_Vmat_CCA_spearman_CI_full, 2, function(x) quantile(x, c(.05, .95)))) 
SingVect1_Vmat_CCA_spearman_CIval_full
