#This script analyzes an Rdata file that is an output of the trainTestFunction_cca/pls_spearman.R script 


#####################################################
#### Extracting predicted singular values
#####################################################

####X LATENT VARIABLES
singVal1_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV1_test_pred))
singVal1_trainTest <- as.vector(singVal1_trainTest)

singVal2_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV2_test_pred))
singVal2_trainTest <- as.vector(singVal2_trainTest)

singVal3_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV3_test_pred))
singVal3_trainTest <- as.vector(singVal3_trainTest)

singVal4_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV4_test_pred))
singVal4_trainTest <- as.vector(singVal4_trainTest)

singVal5_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV5_test_pred))
singVal5_trainTest <- as.vector(singVal5_trainTest)

singVal6_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV6_test_pred))
singVal6_trainTest <- as.vector(singVal6_trainTest)

singVal7_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV7_test_pred))
singVal7_trainTest <- as.vector(singVal7_trainTest)

singVal8_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV8_test_pred))
singVal8_trainTest <- as.vector(singVal8_trainTest)

singVal9_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV9_test_pred))
singVal9_trainTest <- as.vector(singVal9_trainTest)

singVal10_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV10_test_pred))
singVal10_trainTest <- as.vector(singVal10_trainTest)

singVal11_trainTest <- sapply(1:100, function(i) unlist(trainTestFunction_Rdata_file[[i]]$SV11_test_pred))
singVal11_trainTest <- as.vector(singVal11_trainTest)

singVal_trainTest <- cbind.data.frame(singVal1_trainTest, singVal2_trainTest, singVal3_trainTest, singVal4_trainTest, singVal5_trainTest, singVal6_trainTest, 
                                      singVal7_trainTest, singVal8_trainTest, singVal9_trainTest, singVal10_trainTest, singVal11_trainTest)





### calculating the z-score of the distribution of predicted singular values
apply(singVal_trainTest, 2, function(x) (mean(x)) / sd(x))

