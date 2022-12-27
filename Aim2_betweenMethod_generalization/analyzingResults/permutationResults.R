#This script analyzes an Rdata file that is an output of the permutationFunction_cca/pls_spearman.R script 



##################################################################
#### getting the sum of squares for each LV
##################################################################

singVal_spearman <- as.data.frame(svd_scaled_spearman$d)

## this next stretch of code is computing the sum of squares for each latent variable
## for the first LV sum of squares, we're computing the some of squares of all singular values
## then, for each subsequent one, we're removing the preceding singular values 

####1 THROUGH 11
singVal_squared <- as.data.frame(apply(singVal_spearman, 1, function(x) x^2))
singVal_squared_summed1 <- sum(singVal_squared)
#####2 THROUGH 11
singVal2_11 <- as.data.frame(singVal_spearman[-c(1),])
singVal_squared2_11 <- as.data.frame(apply(singVal2_11, 1, function(x) x^2))
singVal_squared_summed2 <- sum(singVal_squared2_11)
######3 THROUGH 11
singVal3_11 <- as.data.frame(singVal_spearman[-c(1, 2),])
singVal_squared3_11 <- as.data.frame(apply(singVal3_11, 1, function(x) x^2))
singVal_squared_summed3 <- sum(singVal_squared3_11)
######4 THROUGH 11
singVal4_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3),])
singVal_squared4_11 <- as.data.frame(apply(singVal4_11, 1, function(x) x^2))
singVal_squared_summed4 <- sum(singVal_squared4_11)
######5 THROUGH 11
singVal5_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4),])
singVal_squared5_11 <- as.data.frame(apply(singVal5_11, 1, function(x) x^2))
singVal_squared_summed5 <- sum(singVal_squared5_11)
######6 THROUGH 11
singVal6_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5),])
singVal_squared6_11 <- as.data.frame(apply(singVal6_11, 1, function(x) x^2))
singVal_squared_summed6 <- sum(singVal_squared6_11)
######7 THROUGH 11
singVal7_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5, 6),])
singVal_squared7_11 <- as.data.frame(apply(singVal7_11, 1, function(x) x^2))
singVal_squared_summed7 <- sum(singVal_squared7_11)
######8 THROUGH 11
singVal8_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5, 6, 7),])
singVal_squared8_11 <- as.data.frame(apply(singVal8_11, 1, function(x) x^2))
singVal_squared_summed8 <- sum(singVal_squared8_11)
######9 THROUGH 11
singVal9_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5, 6, 7, 8),])
singVal_squared9_11 <- as.data.frame(apply(singVal9_11, 1, function(x) x^2))
singVal_squared_summed9 <- sum(singVal_squared9_11)
######10 THROUGH 11
singVal10_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5, 6, 7, 8, 9),])
singVal_squared10_11 <- as.data.frame(apply(singVal10_11, 1, function(x) x^2))
singVal_squared_summed10 <- sum(singVal_squared10_11)
######11 THROUGH 11
singVal11_11 <- as.data.frame(singVal_spearman[-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),])
singVal_squared11_11 <- as.data.frame(apply(singVal11_11, 1, function(x) x^2))
singVal_squared_summed11 <- sum(singVal_squared11_11)



##################################################################
#### Creating data frames of the sum of square singular values 
##################################################################


singVal1_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares1))
singVal1_permutation_results <- as.data.frame(as.vector(singVal1_permutation_results))

singVal2_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares2))
singVal2_permutation_results <- as.data.frame(as.vector(singVal2_permutation_results))

singVal3_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares3))
singVal3_permutation_results <- as.data.frame(as.vector(singVal3_permutation_results))

singVal4_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares4))
singVal4_permutation_results <- as.data.frame(as.vector(singVal4_permutation_results))

singVal5_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares5))
singVal5_permutation_results <- as.data.frame(as.vector(singVal5_permutation_results))

singVal6_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares6))
singVal6_permutation_results <- as.data.frame(as.vector(singVal6_permutation_results))

singVal7_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares7))
singVal7_permutation_results <- as.data.frame(as.vector(singVal7_permutation_results))

singVal8_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares8))
singVal8_permutation_results <- as.data.frame(as.vector(singVal8_permutation_results))

singVal9_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares9))
singVal9_permutation_results <- as.data.frame(as.vector(singVal9_permutation_results))

singVal10_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares10))
singVal10_permutation_results <- as.data.frame(as.vector(singVal10_permutation_results))

singVal11_permutation_results <- sapply(1:100, function(i) unlist(permutation_Rdata_file[[i]]$sum_of_squares11))
singVal11_permutation_results <- as.data.frame(as.vector(singVal11_permutation_results))

##################################################################
#### getting the p-value 
##################################################################

##filtering the number of permuted singular values that are greater or equal to the empirical value
##a singular value (and thereby, LV) will be considered significant if the number of permuted singular values
##that are greater/equal to the empirical singular value is less than 0.05%

pval_singVal1 <- filter(singVal1_permutation_results, singVal1_permutation_results >= singVal_squared_summed1)
pval_singVal2 <- filter(singVal2_permutation_results, singVal2_permutation_results >= singVal_squared_summed2)
pval_singVal3 <- filter(singVal3_permutation_results, singVal3_permutation_results >= singVal_squared_summed3)
pval_singVal4 <- filter(singVal4_permutation_results, singVal4_permutation_results >= singVal_squared_summed4)
pval_singVal5 <- filter(singVal5_permutation_results, singVal5_permutation_results >= singVal_squared_summed5)
pval_singVal6 <- filter(singVal6_permutation_results, singVal6_permutation_results >= singVal_squared_summed6)
pval_singVal7 <- filter(singVal7_permutation_results, singVal7_permutation_results >= singVal_squared_summed7)
pval_singVal8 <- filter(singVal8_permutation_results, singVal8_permutation_results >= singVal_squared_summed8)
pval_singVal9 <- filter(singVal9_permutation_results, singVal9_permutation_results >= singVal_squared_summed9)
pval_singVal10 <- filter(singVal10_permutation_results, singVal10_permutation_results >= singVal_squared_summed10)
pval_singVal11 <- filter(singVal11_permutation_results, singVal11_permutation_results >= singVal_squared_summed11)
