set.seed(1312206748567)
PLS_splithalf_spearman <- function(X, Y, numboot){
  
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
    
    ##getting the correlation and cross-correlation matrices for each split half 
    Rxy_1 <- cor(new_sampleX1[[i]], new_sampleY1[[i]], method = "spearman")
    Rxy_2 <- cor(new_sampleX2[[i]], new_sampleY2[[i]], method = "spearman")
    
    ##performing the decomposition 
    svd1 <- svd(Rxy_1)
    svd2 <- svd(Rxy_2)
    
    ##getting the correlation matrices of the u and v decomposition 
    LV_X_Corr <- cor(svd1$u, svd2$u) 
    LV_Y_Corr <- cor(svd1$v, svd2$v)
    
    ####getting the diagonal and absolute value of the correlation matrices (the sign is not relevant here)
    diagonal_LV_Y_Corr <- as.data.frame(diag(abs(LV_Y_Corr)))
    diagonal_LV_X_Corr <- as.data.frame(diag(abs(LV_X_Corr)))

    ### extracting the first singular vector from each iteration in a single variable 
    LV1_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[1,])
    LV1_corr_elements_X <- rbind(diagonal_LV_X_Corr[1,])

    LV2_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[2,])
    LV2_corr_elements_X <- rbind(diagonal_LV_X_Corr[2,])

    LV3_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[3,])
    LV3_corr_elements_X <- rbind(diagonal_LV_X_Corr[3,])

    LV4_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[4,])
    LV4_corr_elements_X <- rbind(diagonal_LV_X_Corr[4,])

    LV5_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[5,])
    LV5_corr_elements_X <- rbind(diagonal_LV_X_Corr[5,])

    LV6_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[6,])
    LV6_corr_elements_X <- rbind(diagonal_LV_X_Corr[6,])
    
    LV7_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[7,])
    LV7_corr_elements_X <- rbind(diagonal_LV_X_Corr[7,])
    
    LV8_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[8,])
    LV8_corr_elements_X <- rbind(diagonal_LV_X_Corr[8,])
    
    LV9_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[9,])
    LV9_corr_elements_X <- rbind(diagonal_LV_X_Corr[9,])
    
    LV10_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[10,])
    LV10_corr_elements_X <- rbind(diagonal_LV_X_Corr[10,])
    
    LV11_corr_elements_Y <- rbind(diagonal_LV_Y_Corr[11,])
    LV11_corr_elements_X <- rbind(diagonal_LV_X_Corr[11,])
  
  
    
    SingVectors_V1[[i]] <- LV1_corr_elements_Y
    SingVectors_U1[[i]] <- LV1_corr_elements_X
    SingVectors_V2[[i]] <- LV2_corr_elements_Y
    SingVectors_U2[[i]] <- LV2_corr_elements_X
    SingVectors_V3[[i]] <- LV3_corr_elements_Y
    SingVectors_U3[[i]] <- LV3_corr_elements_X
    SingVectors_V4[[i]] <- LV4_corr_elements_Y
    SingVectors_U4[[i]] <- LV4_corr_elements_X
    SingVectors_V5[[i]] <- LV5_corr_elements_Y
    SingVectors_U5[[i]] <- LV5_corr_elements_X
    SingVectors_V6[[i]] <- LV6_corr_elements_Y
    SingVectors_U6[[i]] <- LV6_corr_elements_X
    SingVectors_V7[[i]] <- LV7_corr_elements_Y
    SingVectors_U7[[i]] <- LV7_corr_elements_X
    SingVectors_V8[[i]] <- LV8_corr_elements_Y
    SingVectors_U8[[i]] <- LV8_corr_elements_X
    SingVectors_V9[[i]] <- LV9_corr_elements_Y
    SingVectors_U9[[i]] <- LV9_corr_elements_X
    SingVectors_V10[[i]] <- LV10_corr_elements_Y
    SingVectors_U10[[i]] <- LV10_corr_elements_X
    SingVectors_V11[[i]] <- LV11_corr_elements_Y
    SingVectors_U11[[i]] <- LV11_corr_elements_X

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


PLS_splithalf_spearman_output <- foreach(i=1:100) %dopar% {
  PLS_splithalf_spearman_output <- PLS_splithalf_spearman(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores, Baseline_CBCL_CT_SiteFam_Data_NoSib_CT, 100) #calling a function
}

save(PLS_splithalf_spearman_output, file = "PLS_splithalf_spearman_output.RData")
