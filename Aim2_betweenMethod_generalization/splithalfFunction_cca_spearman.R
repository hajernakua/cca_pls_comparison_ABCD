set.seed(1312206748567)
CCA_splithalf_spearman <- function(X, Y, numboot){
  
  X <- X ##X matrix
  Y <- Y ##Y matrix
  numboot <- numboot # number of permutations 

  
  # Initialize empty list of length numboot for the different outputs we want
  SingVectors_V1 <- vector("list", numboot)
  SingVectors_U1 <- vector("list", numboot)
  SingVectors_V2 <- vector("list", numboot)
  SingVectors_U2 <- vector("list", numboot)
  SingVectors_V3 <- vector("list", numboot)
  SingVectors_U3 <- vector("list", numboot)
  SingVectors_V4 <- vector("list", numboot)
  SingVectors_U4 <- vector("list", numboot)
  SingVectors_V5 <- vector("list", numboot)
  SingVectors_U5 <- vector("list", numboot)
  SingVectors_V6 <- vector("list", numboot)
  SingVectors_U6 <- vector("list", numboot)
  SingVectors_V7 <- vector("list", numboot)
  SingVectors_U7 <- vector("list", numboot)
  SingVectors_V8 <- vector("list", numboot)
  SingVectors_U8 <- vector("list", numboot)
  SingVectors_V9 <- vector("list", numboot)
  SingVectors_U9 <- vector("list", numboot)
  SingVectors_V10 <- vector("list", numboot)
  SingVectors_U10 <- vector("list", numboot)
  SingVectors_V11 <- vector("list", numboot)
  SingVectors_U11 <- vector("list", numboot)
  
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
    print(i) # permutation number
    
    split_half1 <- sort(sample(nrow(X), floor(nrow(X)/2), replace=F))
    split_half2 <- setdiff(1:nrow(X),split_half1)
    
    ##splitting the X matrices based on the indices of the previous steps
    new_sampleX1 <- lapply(1:numboot, function(i) X[split_half1,])
    new_sampleX2 <- lapply(1:numboot, function(i) X[split_half2,]) 
    
    ##splitting the Y matrices based on the indices of the previous steps (this ensures that the same rows are indexed in the X and Y matrices)
    new_sampleY1 <- lapply(1:numboot, function(i) Y[split_half1,])
    new_sampleY2 <- lapply(1:numboot, function(i) Y[split_half2,]) 
    
    ##residualizing the behavioural data (split half 1)
    new_sampleX1[[i]] <- lapply(new_sampleX1[[i]], function(x) lm(x ~new_sampleX1[[i]]$abcd_site + 
                                                                    new_sampleX1[[i]]$interview_age + 
                                                                    new_sampleX1[[i]]$sex +
                                                                    new_sampleX1[[i]]$mri_info_manufacturersmn +
                                                                    new_sampleX1[[i]]$smri_vol_cdk_total)$residuals)
    new_sampleX1[[i]] <- new_sampleX1[[i]][-c(12:16)]
    new_sampleX1[[i]] <- as.data.frame(new_sampleX1[[i]])
    
    ##residualizing the brain data (split half 1)
    new_sampleY1[[i]] <- lapply(new_sampleY1[[i]], function(x) lm(x ~new_sampleY1[[i]]$abcd_site + 
                                                                    new_sampleY1[[i]]$interview_age + 
                                                                    new_sampleY1[[i]]$sex +
                                                                    new_sampleY1[[i]]$mri_info_manufacturersmn +
                                                                    new_sampleY1[[i]]$smri_vol_cdk_total)$residuals)
    new_sampleY1[[i]] <- new_sampleY1[[i]][-c(69:76)]
    new_sampleY1[[i]] <- as.data.frame(new_sampleY1[[i]])
    
    ##getting the z-scores (split half 1)
    new_sampleX1[[i]] <- scale(new_sampleX1[[i]], center = TRUE, scale = TRUE)
    new_sampleY1[[i]] <- scale(new_sampleY1[[i]], center = TRUE, scale = TRUE)
    
    
    ##residualizing the behavioural data (split half 2)
    new_sampleX2[[i]] <- lapply(new_sampleX2[[i]], function(x) lm(x ~new_sampleX2[[i]]$abcd_site + 
                                                                    new_sampleX2[[i]]$interview_age + 
                                                                    new_sampleX2[[i]]$sex +
                                                                    new_sampleX2[[i]]$mri_info_manufacturersmn +
                                                                    new_sampleX2[[i]]$smri_vol_cdk_total)$residuals)
    new_sampleX2[[i]] <- new_sampleX2[[i]][-c(12:16)]
    new_sampleX2[[i]] <- as.data.frame(new_sampleX2[[i]])
    
    ##residualizing the brain data (split half 1)
    new_sampleY2[[i]] <- lapply(new_sampleY2[[i]], function(x) lm(x ~new_sampleY2[[i]]$abcd_site + 
                                                                    new_sampleY2[[i]]$interview_age + 
                                                                    new_sampleY2[[i]]$sex +
                                                                    new_sampleY2[[i]]$mri_info_manufacturersmn +
                                                                    new_sampleY2[[i]]$smri_vol_cdk_total)$residuals)
    new_sampleY2[[i]] <- new_sampleY2[[i]][-c(69:76)]
    new_sampleY2[[i]] <- as.data.frame(new_sampleY2[[i]])
    
    ##getting the z-scores (split half 1)
    new_sampleX2[[i]] <- scale(new_sampleX2[[i]], center = TRUE, scale = TRUE)
    new_sampleY2[[i]] <- scale(new_sampleY2[[i]], center = TRUE, scale = TRUE)
    
    ##calculating the omega matrix
    Sxx_sample1 <- cor(new_sampleX1[[i]], method = "spearman")
    Sxx_sample2 <- cor(new_sampleX2[[i]], method = "spearman")
    Syy_sample1 <- cor(new_sampleY1[[i]], method = "spearman")
    Syy_sample2 <- cor(new_sampleY2[[i]], method = "spearman")
    Rxy_1 <- cor(new_sampleX1[[i]], new_sampleY1[[i]], method = "spearman")
    Rxy_2 <- cor(new_sampleX2[[i]], new_sampleY2[[i]], method = "spearman")
    
    Omega1 <- t(solve(chol(Sxx_sample1))) %*% Rxy_1 %*% solve(chol(Syy_sample1))
    Omega2 <- t(solve(chol(Sxx_sample2))) %*% Rxy_2 %*% solve(chol(Syy_sample2))
    
    ##performing the decomposition 
    svd1 <- svd(Omega1)
    svd2 <- svd(Omega2)
    
    ##getting the beta weights for sample 1
    X_beta_weights1 <- solve(chol(Sxx_sample1)) %*% svd1$u 
    Y_beta_weights1 <- solve(chol(Syy_sample1)) %*% svd1$v 
    
    ##getting the beta weights for sample 2
    X_beta_weights2 <- solve(chol(Sxx_sample2)) %*% svd2$u 
    Y_beta_weights2 <- solve(chol(Syy_sample2)) %*% svd2$v 
    
    ##getting the correlation matrices for beta weights between sample 1 and 2
    BetaWeight.X_Corr <- cor(X_beta_weights1, X_beta_weights2) 
    BetaWeight.Y_Corr <- cor(Y_beta_weights1, Y_beta_weights2)
    
    ####BETA WEIGHTS
    diagonal_LV_Y_Corr_betaWeight <- as.data.frame(diag(abs(BetaWeight.Y_Corr)))
    diagonal_LV_X_Corr_betaWeight <- as.data.frame(diag(abs(BetaWeight.X_Corr)))

    LV1_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[1,])
    LV1_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[1,])

    LV2_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[2,])
    LV2_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[2,])

    LV3_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[3,])
    LV3_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[3,])

    LV4_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[4,])
    LV4_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[4,])

    LV5_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[5,])
    LV5_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[5,])

    LV6_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[6,])
    LV6_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[6,])

    LV7_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[7,])
    LV7_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[7,])

    LV8_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[8,])
    LV8_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[8,])

    LV9_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[9,])
    LV9_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[9,])

    LV10_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[10,])
    LV10_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[10,])

    LV11_corr_elements_Y_betaWeight <- rbind(diagonal_LV_Y_Corr_betaWeight[11,])
    LV11_corr_elements_X_betaWeight <- rbind(diagonal_LV_X_Corr_betaWeight[11,])
  
  
    
    SingVectors_V1[[i]] <- LV1_corr_elements_Y_betaWeight
    SingVectors_U1[[i]] <- LV1_corr_elements_X_betaWeight
    SingVectors_V2[[i]] <- LV2_corr_elements_Y_betaWeight
    SingVectors_U2[[i]] <- LV2_corr_elements_X_betaWeight
    SingVectors_V3[[i]] <- LV3_corr_elements_Y_betaWeight
    SingVectors_U3[[i]] <- LV3_corr_elements_X_betaWeight
    SingVectors_V4[[i]] <- LV4_corr_elements_Y_betaWeight
    SingVectors_U4[[i]] <- LV4_corr_elements_X_betaWeight
    SingVectors_V5[[i]] <- LV5_corr_elements_Y_betaWeight
    SingVectors_U5[[i]] <- LV5_corr_elements_X_betaWeight
    SingVectors_V6[[i]] <- LV6_corr_elements_Y_betaWeight
    SingVectors_U6[[i]] <- LV6_corr_elements_X_betaWeight
    SingVectors_V7[[i]] <- LV7_corr_elements_Y_betaWeight
    SingVectors_U7[[i]] <- LV7_corr_elements_X_betaWeight
    SingVectors_V8[[i]] <- LV8_corr_elements_Y_betaWeight
    SingVectors_U8[[i]] <- LV8_corr_elements_X_betaWeight
    SingVectors_V9[[i]] <- LV9_corr_elements_Y_betaWeight
    SingVectors_U9[[i]] <- LV9_corr_elements_X_betaWeight
    SingVectors_V10[[i]] <- LV10_corr_elements_Y_betaWeight
    SingVectors_U10[[i]] <- LV10_corr_elements_X_betaWeight
    SingVectors_V11[[i]] <- LV11_corr_elements_Y_betaWeight
    SingVectors_U11[[i]] <- LV11_corr_elements_X_betaWeight

  } 
  return(list("SingVectors_V1"=SingVectors_V1, "SingVectors_V2"=SingVectors_V2, "SingVectors_V3"=SingVectors_V3, "SingVectors_V4"=SingVectors_V4,
              "SingVectors_V5"=SingVectors_V5, "SingVectors_V6"=SingVectors_V6, "SingVectors_V7"=SingVectors_V7, "SingVectors_V8"=SingVectors_V8,
              "SingVectors_V9"=SingVectors_V9, "SingVectors_V10"=SingVectors_V10, "SingVectors_V11"=SingVectors_V11,
              "SingVectors_U1"=SingVectors_U1, "SingVectors_U2"=SingVectors_U2, "SingVectors_U3"=SingVectors_U3, "SingVectors_U4"=SingVectors_U4,
              "SingVectors_U5"=SingVectors_U5, "SingVectors_U6"=SingVectors_U6, "SingVectors_U7"=SingVectors_U7, "SingVectors_U8"=SingVectors_U8,
              "SingVectors_U9"=SingVectors_U9, "SingVectors_U10"=SingVectors_U10, "SingVectors_U11"=SingVectors_U11))
}

##running this function in parallel 
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


CCA_splithalf_spearman_output <- foreach(i=1:100) %dopar% {
  CCA_splithalf_spearman_output <- CCA_splithalf_spearman(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores, Baseline_CBCL_CT_SiteFam_Data_NoSib_CT, 100) #calling a function
}
save(CCA_splithalf_spearman_output, file = "CCA_splithalf_spearman_output.RData")
