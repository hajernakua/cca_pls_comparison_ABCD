set.seed(22598)
bootstrap_resampling_PLS_SingularVectors <- function(X, Y, numboot, nsub) {

  numboot <- numboot 
  ##creating a data container to hold all the matrices
  U_mat <- vector("list", numboot)
  V_mat <- vector("list", numboot)
  
  Xrow <- nrow(X)
  Yrow <- nrow(Y)
  y <- nsub
  
  ##getting the proper structure for the variables 
  X$abcd_site <- as.factor(X$abcd_site)
  X$sex <- as.factor(X$sex)
  X$interview_age <- as.numeric(X$interview_age)
  X$mri_info_manufacturersmn <- as.factor(X$mri_info_manufacturersmn)
  X$smri_vol_cdk_total <- as.numeric(X$smri_vol_cdk_total)
  X$smri_vol_cdk_total <- X$smri_vol_cdk_total /1000
  X[,c(1:11)] <- as.data.frame(lapply(X[, c(1:11)], function(x) as.numeric(x)))
  
  Y$abcd_site <- as.factor(Y$abcd_site)
  Y$sex <- as.factor(Y$sex)
  Y$mri_info_manufacturersmn <- as.factor(Y$mri_info_manufacturersmn)
  Y$smri_vol_cdk_total <- as.numeric(Y$smri_vol_cdk_total)
  Y$smri_vol_cdk_total <- Y$smri_vol_cdk_total/1000
  Y[, c(1:68, 76)] <- as.data.frame(lapply(Y[, c(1:68, 76)], function(x) as.numeric(x)))
  
  for(i in 1:numboot) {
    print(i) ##permutation number
    resamples_X <- lapply(1:numboot, function(idx) X[sample(Xrow,y, replace = T),])
    resamples_Y <- lapply(1:numboot, function(idx) Y[sample(Xrow,y, replace = T),])
    
    ##residualizing the behavioural data 
    resamples_X[[i]] <- lapply(resamples_X[[i]], function(x) lm(x ~resamples_X[[i]]$abcd_site + 
                                                                  resamples_X[[i]]$interview_age + 
                                                                  resamples_X[[i]]$sex +
                                                                  resamples_X[[i]]$mri_info_manufacturersmn +
                                                                  resamples_X[[i]]$smri_vol_cdk_total)$residuals)
    resamples_X[[i]] <- resamples_X[[i]][-c(12:16)]
    resamples_X[[i]] <- as.data.frame(resamples_X[[i]])
    ##residualizing the brain data
    resamples_Y[[i]] <- lapply(resamples_Y[[i]], function(x) lm(x ~resamples_Y[[i]]$abcd_site + 
                                                                  resamples_Y[[i]]$interview_age + 
                                                                  resamples_Y[[i]]$sex +
                                                                  resamples_Y[[i]]$mri_info_manufacturersmn +
                                                                  resamples_Y[[i]]$smri_vol_cdk_total)$residuals)
    resamples_Y[[i]] <- resamples_Y[[i]][-c(69:76)]
    resamples_Y[[i]] <- as.data.frame(resamples_Y[[i]])
    
    ##getting the z-scores (split half 1)
    resamples_X[[i]] <- scale(resamples_X[[i]], center = TRUE, scale = TRUE)
    resamples_Y[[i]] <- scale(resamples_Y[[i]], center = TRUE, scale = TRUE)
    
    ##getting the Rxy matrix
    Rxy <- cor(resamples_X[[i]], resamples_Y[[i]], method="spearman")
    
    ##doing the svd
    svd_Rxy <- svd(Rxy)
    
    ##getting the singular vector matrices
    U_matrix <- svd_Rxy$u
    V_matrix <- svd_Rxy$v
    
    U_mat[[i]] <- U_matrix
    V_mat[[i]] <- V_matrix
  }
  
  ##this returns a data container with a series of variables which each have a list of matrices that can be extracted
  return(list("U_mat" = U_mat, "V_mat" = V_mat))
}


bootstrap_PLS_singvectors_spearman <-  bootstrap_resampling_PLS_SingularVectors(phenotypic_data, brain_data, 1000, 9191)
save(bootstrap_PLS_singvectors_spearman, file = "bootstrap_PLS_singvectors_spearman.RData")
