##scripts/05_coexistence_mechanisms_with_facilitation
source("functions/toolbox.R")

library(mvtnorm)
library(parallel)
library(doParallel)
library(dplyr)
library(foreach)
library(lme4)
library(lmerTest)
library(ggplot2)
library(Rmisc)

#matrices for computing, here I did not appliy rm_facilitation to remove the positive interactions
alpha_full_range <- list()
lambda <- list()
for(i in 1:8){
  alpha_full_range[[i]] <- -all_matrices[[i]] #pay attention to the sign of competition coefficients
  lambda[[i]] <- as.data.frame(all_intrinsic[[i]]) %>% dplyr::select( -lwr, -upr)
}

##calculation of coexistence mechanisms based on the sown composition including facilitation
## Use foreach to parallelize the calculation across datasets
# Initialize a cluster for parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

set.seed(123)
# computationally intensive, I provide the code but also load the intermediate results after this calculation
# load the intermediate results: observed coexistence mechanisms based on the sown composition
dataset_results_ini_fac <- readRDS("intermediate_results/dataset_results_ini_fac.rds")

dataset_results_ini_fac <- foreach(i = 1:8, .combine = 'c', .packages = c('mvtnorm','dplyr')) %dopar% {
  
  # Temporary list to store results for this dataset
  result_list <- list()
  
  # Loop over each plot for the current dataset
  for (j in unique(richness_ini[[i]]$plot)) {
    
    # Extract richness for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == j, 2]
    
    # Extract species combination for the current plot (as a character vector)
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == j, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- combination1[combination1 %in% rownames(alpha_full_range[[i]])]
    combination1 <- rownames(alpha_full_range[[i]])[rownames(alpha_full_range[[i]]) %in% combination1]  # reorder to match alpha
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_full_range[[i]],          # Alpha matrix for dataset i
        intrinsic = lambda[[i]],     # Intrinsic matrix for dataset i
        n = n1,                      # Richness for the current plot
        combination = combination1,  # Species combination for the current plot
        list_names = names(all_matrices)[i],  # Dataset name
        plot = j                     # The specific plot being analyzed
      )
      
      # Store the result for this plot
      result_list[[paste("Dataset", i, "Plot", j)]] <- result1
    } else {
      warning(paste("No valid species found for plot", j, "in dataset", i))
    }
  }
  
  return(result_list)
}
# Stop the cluster after completion
#stopCluster(cl)
dataset_results_ini1_fac <- bind_rows(dataset_results_ini_fac)
richness_ini_1_fac<- bind_rows(richness_ini)
comb_str_ini_fac <-dataset_results_ini1_fac[,c(1:8, 15:19)]
comb_str_ini_1_fac<-data.frame(comb_str_ini_fac, richness_ini_1_fac$richness)
comb_str_ini_1_fac$season <- c(rep("June", 136), rep("August", 136))
str(comb_str_ini_1_fac)

#load the saved intermediate results of 99 times randomized_matrices_new_fac and of 99 times coexistence mechanisms calculated based on full range of competition including facilitation
# The code was same as I used in script 3, just replaced the alpha to alpha_full_range, and I don't repeat here
randomized_matrices_new_fac <- readRDS("intermediate_results/randomized_matrices_new_fac.rds")

dataset_results_ini_rand1_fac <- readRDS("intermediate_results/dataset_results_ini_rand1_fac.rds")
binded_results_ini_rand_C_june_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand1_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_C_june_fac[[j]] <- selected_columns
}
comb_str_ini_rand1_fac <-bind_rows(binded_results_ini_rand_C_june_fac)
summary(as.numeric(subset(comb_str_ini_rand1_fac, plot %in% c(21,117,247,290))$theta))#for 20 species plot Omega and theta are exactly the same

dataset_results_ini_rand2_fac <- readRDS("intermediate_results/dataset_results_ini_rand2_fac.rds")
binded_results_ini_rand_N_june_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand2_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_N_june_fac[[j]] <- selected_columns
}
comb_str_ini_rand2_fac <-bind_rows(binded_results_ini_rand_N_june_fac)
summary(as.numeric(subset(comb_str_ini_rand2_fac, plot %in% c(43,59,99,308))$Omega))

dataset_results_ini_rand3_fac <- readRDS("intermediate_results/dataset_results_ini_rand3_fac.rds")
# Combine all result_list entries into a single data frame
binded_results_ini_rand_F_june_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand3_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_F_june_fac[[j]] <- selected_columns
}
comb_str_ini_rand3_fac <-bind_rows(binded_results_ini_rand_F_june_fac)
summary(as.numeric(subset(comb_str_ini_rand3_fac, plot==158)$Omega))##all 20 plots are same

dataset_results_ini_rand4_fac <- readRDS("intermediate_results/dataset_results_ini_rand4_fac.rds")
# Combine all result_list entries into a single data frame
binded_results_ini_rand_NF_june_fac<- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand4_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_NF_june_fac[[j]] <- selected_columns
}
comb_str_ini_rand4_fac <-bind_rows(binded_results_ini_rand_NF_june_fac)
summary(as.numeric(subset(comb_str_ini_rand4_fac, plot == 83)$Omega))#all the same

dataset_results_ini_rand5_fac <- readRDS("intermediate_results/dataset_results_ini_rand5_fac.rds")
binded_results_ini_rand_C_August_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand5_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_C_August_fac[[j]] <- selected_columns
}
comb_str_ini_rand5_fac <-bind_rows(binded_results_ini_rand_C_August_fac)
summary(as.numeric(subset(comb_str_ini_rand5_fac, plot %in% c(21, 117, 247, 290))$Omega))

dataset_results_ini_rand6_fac <- readRDS("intermediate_results/dataset_results_ini_rand6_fac.rds")
# Combine all result_list entries into a single data frame
binded_results_ini_rand_N_August_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand6_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_N_August_fac[[j]] <- selected_columns
}
comb_str_ini_rand6_fac <-bind_rows(binded_results_ini_rand_N_August_fac)
summary(as.numeric(subset(comb_str_ini_rand6_fac, plot %in% c(59, 99, 308) )$Omega))#all the same

dataset_results_ini_rand7_fac <- readRDS("intermediate_results/dataset_results_ini_rand7_fac.rds")
# Combine all result_list entries into a single data frame
binded_results_ini_rand_F_August_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand7_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_F_August_fac[[j]] <- selected_columns
}
comb_str_ini_rand7_fac <-bind_rows(binded_results_ini_rand_F_August_fac)
summary(as.numeric(subset(comb_str_ini_rand7_fac, plot %in% c(158, 171, 295, 328) )$Omega))

dataset_results_ini_rand8_fac <- readRDS("intermediate_results/dataset_results_ini_rand8_fac.rds")
# Combine all result_list entries into a single data frame
binded_results_ini_rand_NF_August_fac <- list()
# Loop over each iteration (j = 1 to 99)
for (j in 1:num_iterations) {
  # Extract the result_list for the current iteration and combine rows
  result_list_df <- bind_rows(dataset_results_ini_rand8_fac[[j]]$result_list)
  selected_columns <- result_list_df[, c("Omega", "Omega1", "Omega2", "theta", "theta1", "differential", "overlap","feasibility", "coex_rate","treatment", "theta2", "distance_to_exclusion", "plot")]
  # Store the selected columns for the current iteration
  binded_results_ini_rand_NF_August_fac[[j]] <- selected_columns
}
comb_str_ini_rand8_fac <-bind_rows(binded_results_ini_rand_NF_August_fac)
summary(as.numeric(subset(comb_str_ini_rand8_fac, plot==233)$Omega))#all the same

##
comb_str_ini_rand_fac<-rbind(comb_str_ini_rand1_fac, comb_str_ini_rand2_fac, comb_str_ini_rand3_fac, comb_str_ini_rand4_fac,
                            comb_str_ini_rand5_fac, comb_str_ini_rand6_fac, comb_str_ini_rand7_fac, comb_str_ini_rand8_fac)

comb_str_ini_rand_fac$season <- c(rep("June", 13464), rep("August", 13464))
comb_str_ini_rand_fac$Treatment <- c(rep("C", 3366), rep("N", 3366), rep("F", 3366), rep("NF", 3366),
                                     rep("C", 3366), rep("N", 3366), rep("F", 3366), rep("NF", 3366))
comb_str_ini_rand_fac$year<-2016
comb_str_ini_rand_fac$Omega<-as.numeric(comb_str_ini_rand_fac$Omega)
comb_str_ini_rand_fac$theta<-as.numeric(comb_str_ini_rand_fac$theta)
comb_str_ini_rand_fac$distance_to_exclusion<-as.numeric(comb_str_ini_rand_fac$distance_to_exclusion)
comb_str_ini_rand_fac$differential<-as.numeric(comb_str_ini_rand_fac$differential)

comb_str_ini_rand_fac_plot.info<-merge(comb_str_ini_rand_fac, plot_descriptions, all.x = F, by = c("plot") )
summary(subset(comb_str_ini_rand_fac_plot.info, plot %in% c(43, 59, 99, 308) & season=="June"))#20 speceis plots have the same omega and theta, differential in 20 species plot was 0.

## check the variance of 99 times randomization
## Calculate mean of 'value' based on 'group'
## Aggregate to calculate the mean for multiple columns
mean_values_rand_ini_fac <- aggregate(cbind(Omega, theta, distance_to_exclusion, differential) ~ plot + season, 
                                      data = comb_str_ini_rand_fac, 
                                      FUN = function(x)  mean(x, na.rm = TRUE),  # Custom NA handling
                                      na.action = na.pass) # Prevent dropping of NA rows
names(mean_values_rand_ini_fac)[3:6] <-c("mean_omega_rand","mean_theta_rand","mean_distance_to_exclusion", "mean_differential") 

sd_values_rand_ini_fac <- aggregate(cbind(Omega, theta, distance_to_exclusion, differential) ~ plot + season, 
                                    data = comb_str_ini_rand_fac, 
                                    FUN = function(x) sd(x, na.rm = TRUE),  # Calculate sd
                                    na.action = na.pass)
names(sd_values_rand_ini_fac)[3:6] <-c("sd_omega_rand","sd_theta_rand","sd_distance_to_exclusion", "sd_differential") 
str(mean_values_rand_ini_fac)#272 mean values for 99 random

merged_data_ini_rand_fac <- comb_str_ini_1_fac %>%
  left_join(mean_values_rand_ini_fac, by = c("plot","season"))%>%left_join(sd_values_rand_ini_fac, by = c("plot","season"))

merged_data_ini_rand_fac$Omega<-as.numeric(merged_data_ini_rand_fac$Omega)
merged_data_ini_rand_fac$theta<-as.numeric(merged_data_ini_rand_fac$theta)
merged_data_ini_rand_fac$differential<-as.numeric(merged_data_ini_rand_fac$differential)
merged_data_ini_rand_fac$year<-2016

merged_data_ini_rand_fac$Treatment<- c(rep("C", 34), rep("N", 34), rep("F", 34), rep("NF", 34),
                                       rep("C", 34), rep("N", 34), rep("F", 34), rep("NF", 34))
combined_sown_coex_fac<-merge(merged_data_ini_rand_fac, plot_descriptions, all.x = F, by = c("plot") )

combined_sown_coex_fac$Omega<-as.numeric(combined_sown_coex_fac$Omega)
combined_sown_coex_fac$theta<-as.numeric(combined_sown_coex_fac$theta)
combined_sown_coex_fac$differential<-as.numeric(combined_sown_coex_fac$differential)

##sta.Omega: predicted niche differences with facilitation; 
##sta.theta: predicted competitive ability differences with facilitation; 
##sta.differential: predicted indirect interactions with facilitation; 
combined_sown_coex_fac$sta.Omega<-(combined_sown_coex_fac$Omega-combined_sown_coex_fac$mean_omega_rand)/combined_sown_coex_fac$sd_omega_rand
combined_sown_coex_fac$sta.theta<-(combined_sown_coex_fac$theta-combined_sown_coex_fac$mean_theta_rand)/combined_sown_coex_fac$sd_theta_rand
combined_sown_coex_fac$sta.differential<-(combined_sown_coex_fac$differential-combined_sown_coex_fac$mean_differential)

