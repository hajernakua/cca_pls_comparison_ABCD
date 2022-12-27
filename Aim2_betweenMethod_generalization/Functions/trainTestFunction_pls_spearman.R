set.seed(126757489)
PLS_traintest_spearman <- function(X, Y, numboot){
  
  # Initialize empty list of length numboot for the different outputs we want
  SV1_test_pred <- vector("list", numboot)
  SV2_test_pred <- vector("list", numboot)
  SV3_test_pred <- vector("list", numboot)
  SV4_test_pred <- vector("list", numboot)
  SV5_test_pred <- vector("list", numboot)
  SV6_test_pred <- vector("list", numboot)
  SV7_test_pred <- vector("list", numboot)
  SV8_test_pred <- vector("list", numboot)
  SV9_test_pred <- vector("list", numboot)
  SV10_test_pred <- vector("list", numboot)
  SV11_test_pred <- vector("list", numboot)

  
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
    
    ##setting up the indices each of the split halfs
    train <- sort(sample(nrow(X), round(0.8*nrow(X)), replace=F))
    test <- setdiff(1:nrow(X),train)
    
    ##getting 80% (training data) of the X and Y matrix
    training_X1 <- lapply(1:numboot, function(i) X[train,])
    training_Y1 <- lapply(1:numboot, function(i) Y[train,])
    
    ##getting the 20% (training data) of the X and Y matrix
    testing_X2 <- lapply(1:numboot, function(i) X[test,])
    testing_Y2 <- lapply(1:numboot, function(i) Y[test,])
    
    
    ###residualizing the training data (X - behaviour)
    training_X1[[i]] <- lapply(training_X1[[i]], function(x) lm(x ~training_X1[[i]]$abcd_site + 
                                                                  training_X1[[i]]$interview_age + 
                                                                  training_X1[[i]]$sex +
                                                                  training_X1[[i]]$mri_info_manufacturersmn +
                                                                  training_X1[[i]]$smri_vol_cdk_total)$residuals)
    training_X1[[i]] <- training_X1[[i]][-c(12:16)]
    training_X1[[i]] <- as.data.frame(training_X1[[i]])
    
    ##residualizing the training (Y - brain)
    training_Y1[[i]] <- lapply(training_Y1[[i]], function(x) lm(x ~training_Y1[[i]]$abcd_site + 
                                                                  training_Y1[[i]]$interview_age + 
                                                                  training_Y1[[i]]$sex +
                                                                  training_Y1[[i]]$mri_info_manufacturersmn +
                                                                  training_Y1[[i]]$smri_vol_cdk_total)$residuals)
    training_Y1[[i]] <- training_Y1[[i]][-c(69:76)]
    training_Y1[[i]] <- as.data.frame(training_Y1[[i]])
    
    ##z-transforming training data
    training_X1[[i]] <- scale(training_X1[[i]], center = TRUE, scale = TRUE)
    training_Y1[[i]] <- scale(training_Y1[[i]], center = TRUE, scale = TRUE)
    
    
    ###residualizing the testing data (X - behaviour)
    testing_X2[[i]] <- lapply(testing_X2[[i]], function(x) lm(x ~testing_X2[[i]]$abcd_site + 
                                                                testing_X2[[i]]$interview_age + 
                                                                testing_X2[[i]]$sex +
                                                                testing_X2[[i]]$mri_info_manufacturersmn +
                                                                testing_X2[[i]]$smri_vol_cdk_total)$residuals)
    testing_X2[[i]] <- testing_X2[[i]][-c(12:16)]
    testing_X2[[i]] <- as.data.frame(testing_X2[[i]])
    
    ##residualizing the testing data (Y - brain)
    testing_Y2[[i]] <- lapply(testing_Y2[[i]], function(x) lm(x ~testing_Y2[[i]]$abcd_site + 
                                                                testing_Y2[[i]]$interview_age + 
                                                                testing_Y2[[i]]$sex +
                                                                testing_Y2[[i]]$mri_info_manufacturersmn +
                                                                testing_Y2[[i]]$smri_vol_cdk_total)$residuals)
    testing_Y2[[i]] <- testing_Y2[[i]][-c(69:76)]
    testing_Y2[[i]] <- as.data.frame(testing_Y2[[i]])

    ##z-transforming testing data
    testing_X1[[i]] <- scale(testing_X1[[i]], center = TRUE, scale = TRUE)
    testing_Y1[[i]] <- scale(testing_Y1[[i]], center = TRUE, scale = TRUE)
    
    ##getting Rxy matrix of the training subset
    Rxy_train <- cor(training_Y1[[i]], training_X1[[i]], method = "spearman")
    
    #decomposition 
    svd_Rxy <- svd(Rxy_train)
    
    ##getting Rxy matrix of the testing subset
    Rxy_test <- cor(testing_Y2[[i]], testing_X2[[i]], method = "spearman")
    
    ##arranging the equation to solve for the singular values of the test set
    sing_val_test_prediction <- t(svd_Rxy$u) %*% Rxy_test %*% svd_Rxy$v
    sing_val_test_prediction <- as.data.frame(diag(abs(sing_val_test_prediction)))
    
    ##extracting the individual singular value elements of each LV
    SV1_test <- rbind(sing_val_test_prediction[1,])
    SV2_test <- rbind(sing_val_test_prediction[2,])
    SV3_test <- rbind(sing_val_test_prediction[3,])
    SV4_test <- rbind(sing_val_test_prediction[4,])
    SV5_test <- rbind(sing_val_test_prediction[5,])
    SV6_test <- rbind(sing_val_test_prediction[6,])
    SV7_test <- rbind(sing_val_test_prediction[7,])
    SV8_test <- rbind(sing_val_test_prediction[8,])
    SV9_test <- rbind(sing_val_test_prediction[9,])
    SV10_test <- rbind(sing_val_test_prediction[10,])
    SV11_test <- rbind(sing_val_test_prediction[11,])
    
    SV1_test_pred[[i]] <- SV1_test
    SV2_test_pred[[i]] <- SV2_test
    SV3_test_pred[[i]] <- SV3_test
    SV4_test_pred[[i]] <- SV4_test
    SV5_test_pred[[i]] <- SV5_test
    SV6_test_pred[[i]] <- SV6_test
    SV7_test_pred[[i]] <- SV7_test
    SV8_test_pred[[i]] <- SV8_test
    SV9_test_pred[[i]] <- SV9_test
    SV10_test_pred[[i]] <- SV10_test
    SV11_test_pred[[i]] <- SV11_test
    
    
  } 
  return(list("SV1_test_pred"=SV1_test_pred, "SV2_test_pred"=SV2_test_pred, "SV3_test_pred"=SV3_test_pred, "SV4_test_pred"=SV4_test_pred,
              "SV5_test_pred"=SV5_test_pred, "SV6_test_pred"=SV6_test_pred, "SV7_test_pred"=SV7_test_pred, "SV8_test_pred"=SV8_test_pred, "SV9_test_pred"=SV9_test_pred, "SV10_test_pred"=SV10_test_pred,
              "SV11_test_pred"=SV11_test_pred))
}

##running this function in parallel
library(foreach)
library(doParallel)

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, type="FORK") 
registerDoParallel(cl)

PLS_testTrain_spearman_output <- foreach(i=1:100) %dopar% {
  PLS_testTrain_spearman_output <- PLS_traintest_spearman(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores, Baseline_CBCL_CT_SiteFam_Data_NoSib_CT, 100, 9191, 0.5) #calling a function
}

save(PLS_testTrain_spearman_output, file = "PLS_testTrain_spearman_output.RData")
