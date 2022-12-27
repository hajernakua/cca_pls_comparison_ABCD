#This script analyzes an Rdata file that is an output of the splithalfFunction_cca/pls_spearman.R script 


#####################################
#### Extracting singular vectors
#####################################

##for PLS analyses, we extract the loadings from the U and V matrices
##for CCA analyses, we extract the loadings from the beta weights (re-weighted U and V matrices)

####X LATENT VARIABLES
SingVecU1_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec1_U))
SingVecU1_splitHalf <- as.data.frame(as.vector(SingVecU1_splitHalf))

SingVecU2_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec2_U))
SingVecU2_splitHalf <- as.data.frame(as.vector(SingVecU2_splitHalf))

SingVecU3_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec3_U))
SingVecU3_splitHalf <- as.data.frame(as.vector(SingVecU3_splitHalf))

SingVecU4_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec4_U))
SingVecU4_splitHalf <- as.data.frame(as.vector(SingVecU4_splitHalf))

SingVecU5_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec5_U))
SingVecU5_splitHalf <- as.data.frame(as.vector(SingVecU5_splitHalf))

SingVecU6_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec6_U))
SingVecU6_splitHalf <- as.data.frame(as.vector(SingVecU6_splitHalf))

SingVecU7_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec7_U))
SingVecU7_splitHalf <- as.data.frame(as.vector(SingVecU7_splitHalf))

SingVecU8_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec8_U))
SingVecU8_splitHalf <- as.data.frame(as.vector(SingVecU8_splitHalf))

SingVecU9_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec9_U))
SingVecU9_splitHalf <- as.data.frame(as.vector(SingVecU9_splitHalf))

SingVecU10_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec10_U))
SingVecU10_splitHalf <- as.data.frame(as.vector(SingVecU10_splitHalf))

SingVecU11_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec11_U))
SingVecU11_splitHalf <- as.data.frame(as.vector(SingVecU11_splitHalf))

SingVecU_splitHalf <- cbind.data.frame(SingVecU1_splitHalf, SingVecU2_splitHalf, SingVecU3_splitHalf, SingVecU4_splitHalf, SingVecU5_splitHalf, SingVecU6_splitHalf, 
                                       SingVecU7_splitHalf, SingVecU8_splitHalf, SingVecU9_splitHalf, SingVecU10_splitHalf, SingVecU11_splitHalf)


##Y VARIABLES
SingVecV1_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec1_V))
SingVecV1_splitHalf <- as.data.frame(as.vector(SingVecV1_splitHalf))

SingVecV2_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec2_V))
SingVecV2_splitHalf <- as.data.frame(as.vector(SingVecV2_splitHalf))

SingVecV3_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec3_V))
SingVecV3_splitHalf <- as.data.frame(as.vector(SingVecV3_splitHalf))

SingVecV4_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec4_V))
SingVecV4_splitHalf <- as.data.frame(as.vector(SingVecV4_splitHalf))

SingVecV5_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec5_V))
SingVecV5_splitHalf <- as.data.frame(as.vector(SingVecV5_splitHalf))

SingVecV6_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec6_V))
SingVecV6_splitHalf <- as.data.frame(as.vector(SingVecV6_splitHalf))

SingVecV7_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec7_V))
SingVecV7_splitHalf <- as.data.frame(as.vector(SingVecV7_splitHalf))

SingVecV8_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec8_V))
SingVecV8_splitHalf <- as.data.frame(as.vector(SingVecV8_splitHalf))

SingVecV9_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec9_V))
SingVecV9_splitHalf <- as.data.frame(as.vector(SingVecV9_splitHalf))

SingVecV10_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec10_V))
SingVecV10_splitHalf <- as.data.frame(as.vector(SingVecV10_splitHalf))

SingVecV11_splitHalf <- sapply(1:100, function(i) unlist(splithalfFunction_Rdata_file[[i]]$singVec11_V))
SingVecV11_splitHalf <- as.data.frame(as.vector(SingVecV11_splitHalf))

SingVecV_splitHalf <- cbind.data.frame(SingVecV1_splitHalf, SingVecV2_splitHalf, SingVecV3_splitHalf, SingVecV4_splitHalf, SingVecV5_splitHalf, SingVecV6_splitHalf, 
                                       SingVecV7_splitHalf, SingVecV8_splitHalf, SingVecV9_splitHalf, SingVecV10_splitHalf, SingVecV11_splitHalf)



###################################
#### Calculating reproducibility
###################################

###calculating the z-score for the correlations of singular vectors across the split-half iterations 

apply(SingVecU_splitHalf, 2, function(x) (mean(x)) / sd(x))
apply(SingVecV_splitHalf, 2, function(x) (mean(x)) / sd(x))

## calculating the mean of the distributions 
apply(SingVecU_splitHalf, 2, function(x) (mean(x)))
apply(SingVecV_splitHalf, 2, function(x) (mean(x)))
