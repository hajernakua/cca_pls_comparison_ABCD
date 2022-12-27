set.seed(12348)
PLS_permutation_spearman <- function(X, Y, numboot) {
  
  ##creating data container that we will use to input the singular values later on
  sum_of_squares1 <- vector("list", numboot)
  sum_of_squares2 <- vector("list", numboot)
  sum_of_squares3 <- vector("list", numboot)
  sum_of_squares4 <- vector("list", numboot)
  sum_of_squares5 <- vector("list", numboot)
  sum_of_squares6 <- vector("list", numboot)
  sum_of_squares7 <- vector("list", numboot)
  sum_of_squares8 <- vector("list", numboot)
  sum_of_squares9 <- vector("list", numboot)
  sum_of_squares10 <- vector("list", numboot)
  sum_of_squares11 <- vector("list", numboot)
  
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
  Y[, c(1:68, 73)] <- as.data.frame(lapply(Y[, c(1:68, 73)], function(x) as.numeric(x)))
  
  
  for(i in 1:numboot) {
    
    print(i) 
    
    resamplesX <- sample(nrow(X), replace=F)
    resamples_X <- lapply(1:numboot, function(i) X[resamplesX,])
    resamples_Y <- lapply(1:numboot, function(i) Y)
    
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
    resamples_Y[[i]] <- resamples_Y[[i]][-c(69:73)]
    resamples_Y[[i]] <- as.data.frame(resamples_Y[[i]])
    
    ##getting the z-scores
    resamples_X[[i]] <- scale(resamples_X[[i]], center = TRUE, scale = TRUE)
    resamples_Y[[i]] <- scale(resamples_Y[[i]], center = TRUE, scale = TRUE)
    
    ##steps for the CCA analysis prior to decomposition 
    Rxy <- cor(resamples_X[[i]], resamples_Y[[i]], method = "spearman")
    
    ##decomposition
    decomposition <- svd(Rxy)
    
    ##getting the singular values
    singular_values <- as.data.frame(decomposition$d)
    
    ##getting the sum of square singular values 
    sing_values_squared <- as.data.frame(apply(singular_values, 1, function(x) x^2))
    sing_values_squared_summed1 <- sum(sing_values_squared)

    singular_values2_11 <- as.data.frame(singular_values[-c(1),])
    sing_values_squared2_11 <- as.data.frame(apply(singular_values2_11, 1, function(x) x^2))
    sing_values_squared_summed2 <- sum(sing_values_squared2_11)

    singular_values3_11 <- as.data.frame(singular_values[-c(1, 2),])
    sing_values_squared3_11 <- as.data.frame(apply(singular_values3_11, 1, function(x) x^2))
    sing_values_squared_summed3 <- sum(sing_values_squared3_11)

    singular_values4_11 <- as.data.frame(singular_values[-c(1, 2, 3),])
    sing_values_squared4_11 <- as.data.frame(apply(singular_values4_11, 1, function(x) x^2))
    sing_values_squared_summed4 <- sum(sing_values_squared4_11)

    singular_values5_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4),])
    sing_values_squared5_11 <- as.data.frame(apply(singular_values5_11, 1, function(x) x^2))
    sing_values_squared_summed5 <- sum(sing_values_squared5_11)

    singular_values6_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5),])
    sing_values_squared6_11 <- as.data.frame(apply(singular_values6_11, 1, function(x) x^2))
    sing_values_squared_summed6 <- sum(sing_values_squared6_11)

    singular_values7_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5, 6),])
    sing_values_squared7_11 <- as.data.frame(apply(singular_values7_11, 1, function(x) x^2))
    sing_values_squared_summed7 <- sum(sing_values_squared7_11)

    singular_values8_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5, 6, 7),])
    sing_values_squared8_11 <- as.data.frame(apply(singular_values8_11, 1, function(x) x^2))
    sing_values_squared_summed8 <- sum(sing_values_squared8_11)

    singular_values9_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5, 6, 7, 8),])
    sing_values_squared9_11 <- as.data.frame(apply(singular_values9_11, 1, function(x) x^2))
    sing_values_squared_summed9 <- sum(sing_values_squared9_11)

    singular_values10_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5, 6, 7, 8, 9),])
    sing_values_squared10_11 <- as.data.frame(apply(singular_values10_11, 1, function(x) x^2))
    sing_values_squared_summed10 <- sum(sing_values_squared10_11)

    singular_values11_11 <- as.data.frame(singular_values[-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),])
    sing_values_squared11_11 <- as.data.frame(apply(singular_values11_11, 1, function(x) x^2))
    sing_values_squared_summed11 <- sum(sing_values_squared11_11)
    
    
    sum_of_squares1[[i]] <- sing_values_squared_summed1
    sum_of_squares2[[i]] <- sing_values_squared_summed2
    sum_of_squares3[[i]] <- sing_values_squared_summed3
    sum_of_squares4[[i]] <- sing_values_squared_summed4
    sum_of_squares5[[i]] <- sing_values_squared_summed5
    sum_of_squares6[[i]] <- sing_values_squared_summed6
    sum_of_squares7[[i]] <- sing_values_squared_summed7
    sum_of_squares8[[i]] <- sing_values_squared_summed8
    sum_of_squares9[[i]] <- sing_values_squared_summed9
    sum_of_squares10[[i]] <- sing_values_squared_summed10
    sum_of_squares11[[i]] <- sing_values_squared_summed11


  }
  return(list("sum_of_squares1"=sum_of_squares1, "sum_of_squares2"=sum_of_squares2, "sum_of_squares3"=sum_of_squares3, "sum_of_squares4"=sum_of_squares4,
              "sum_of_squares5"=sum_of_squares5, "sum_of_squares6"=sum_of_squares6, "sum_of_squares7"=sum_of_squares7, "sum_of_squares8"=sum_of_squares8, "sum_of_squares9"=sum_of_squares9, "sum_of_squares10"=sum_of_squares10,
              "sum_of_squares11"=sum_of_squares11))
}

##running this function in parallel 
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


PLS_permutation_spearman_output <- foreach(i=1:100) %dopar% {
  PLS_permutation_spearman_output <- PLS_permutation_spearman(Baseline_CBCL_CT_SiteFam_Data_NoSib_cbclScores, Baseline_CBCL_CT_SiteFam_Data_NoSib_CT, 100) #calling a function
}
