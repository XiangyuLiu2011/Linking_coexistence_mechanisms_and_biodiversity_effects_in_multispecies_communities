##scripts/03_coexistence_mechanisms_without_facilitation
#This is to calculate coexistence mechanisms without facilitation with competition coefficients and intrinsic growth rates

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

#1. load competition coefficients (alpha), intrinsic growth rates (lambda)
#open matrices and set dataframe
alpha_june_control<-read.table("data/biommatrix_june_2020_control.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_june_control<-read.table("data/SE_biommatrix_june_2020_control.txt",header = T, sep = "\t", row.names = 1)
alpha_june_control<-alpha_june_control/SE_alpha_june_control
alpha_june_control<-as.data.frame(alpha_june_control)

alpha_june_nitrogen<-read.table("data/biommatrix_june_2020_nitrogen.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_june_nitrogen<-read.table("data/SE_biommatrix_june_2020_nitrogen.txt",header = T, sep = "\t", row.names = 1)
alpha_june_nitrogen<-alpha_june_nitrogen/SE_alpha_june_nitrogen
alpha_june_nitrogen<-as.data.frame(alpha_june_nitrogen)

alpha_june_fungicide<-read.table("data/biommatrix_june_2020_fungicide.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_june_fungicide<-read.table("data/SE_biommatrix_june_2020_fungicide.txt",header = T, sep = "\t", row.names = 1)
alpha_june_fungicide<-alpha_june_fungicide/SE_alpha_june_fungicide
alpha_june_fungicide<-as.data.frame(alpha_june_fungicide)

alpha_june_combined<-read.table("data/biommatrix_june_2020_combined.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_june_combined<-read.table("data/SE_biommatrix_june_2020_combined.txt",header = T, sep = "\t", row.names = 1)
alpha_june_combined<-alpha_june_combined/SE_alpha_june_combined
alpha_june_combined<-as.data.frame(alpha_june_combined)

alpha_june_control<-as.matrix(alpha_june_control)
alpha_june_nitrogen<-as.matrix(alpha_june_nitrogen)
alpha_june_fungicide<-as.matrix(alpha_june_fungicide)
alpha_june_combined<-as.matrix(alpha_june_combined)

alpha_august_control<-read.table("data/biommatrix_august_2020_control.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_august_control<-read.table("data/SE_biommatrix_august_2020_control.txt",header = T, sep = "\t", row.names = 1)
alpha_august_control<-alpha_august_control/SE_alpha_august_control
alpha_august_control<-as.data.frame(alpha_august_control)

alpha_august_nitrogen<-read.table("data/biommatrix_august_2020_nitrogen.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_august_nitrogen<-read.table("data/SE_biommatrix_august_2020_nitrogen.txt",header = T, sep = "\t", row.names = 1)
alpha_august_nitrogen<-alpha_august_nitrogen/SE_alpha_august_nitrogen
alpha_august_nitrogen<-as.data.frame(alpha_august_nitrogen)

alpha_august_fungicide<-read.table("data/biommatrix_august_2020_fungicide.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_august_fungicide<-read.table("data/SE_biommatrix_august_2020_fungicide.txt",header = T, sep = "\t", row.names = 1)
alpha_august_fungicide<-alpha_august_fungicide/SE_alpha_august_fungicide
alpha_august_fungicide<-as.data.frame(alpha_august_fungicide)

alpha_august_combined<-read.table("data/biommatrix_august_2020_combined.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_august_combined<-read.table("data/SE_biommatrix_august_2020_combined.txt",header = T, sep = "\t", row.names = 1)
alpha_august_combined<-alpha_august_combined/SE_alpha_august_combined
alpha_august_combined<-as.data.frame(alpha_august_combined)

alpha_august_control<-as.matrix(alpha_august_control)
alpha_august_nitrogen<-as.matrix(alpha_august_nitrogen)
alpha_august_fungicide<-as.matrix(alpha_august_fungicide)
alpha_august_combined<-as.matrix(alpha_august_combined)

lambda_june_control <- read.table("data/biomintrinsic_june_2020_control.txt",  sep = "\t", header = T)
lambda_june_nitrogen <- read.table( "data/biomintrinsic_june_2020_nitrogen.txt",  sep = "\t", header = T)
lambda_august_control <- read.table( "data/biomintrinsic_august_2020_control.txt",  sep = "\t", header = T)
lambda_august_nitrogen <- read.table( "data/biomintrinsic_august_2020_nitrogen.txt",  sep = "\t", header = T)
lambda_june_fungicide <- read.table( "data/biomintrinsic_june_2020_fungicide.txt",  sep = "\t", header = T)
lambda_august_fungicide <- read.table( "data/biomintrinsic_august_2020_fungicide.txt",  sep = "\t", header = T)
lambda_june_combined <- read.table( "data/biomintrinsic_june_2020_combined.txt",  sep = "\t", header = T)
lambda_august_combined <- read.table( "data/biomintrinsic_august_2020_combined.txt",  sep = "\t", header = T)
species<-rownames(alpha_june_control)
##
#make a list of interaction matrices
all_matrices<-list(alpha_june_control,alpha_june_nitrogen,alpha_june_fungicide,alpha_june_combined,alpha_august_control,alpha_august_nitrogen,alpha_august_fungicide,alpha_august_combined)
names(all_matrices)<-c("June C", "June N", "June F", "June NF", "August C", "August N", "August F", "August NF")

all_intrinsic<-list(lambda_june_control,lambda_june_nitrogen,lambda_june_fungicide,lambda_june_combined,lambda_august_control,lambda_august_nitrogen,lambda_august_fungicide,lambda_august_combined)
names(all_intrinsic)<-c("June C", "June N", "June F", "June NF", "August C", "August N", "August F", "August NF")
all_intrinsic<-lapply(all_intrinsic, function(x){exp(x)})
all_intrinsic<-lapply(all_intrinsic, function(x){x/1000})

#run a function that replaces facilitation with a random number between the minimum competition and half of the minimum of competition
all_matrices_competition<-lapply(all_matrices, rm_facilitation)

#matrices for computing
alpha <- list()
lambda <- list()
for(i in 1:8){
  alpha[[i]] <- -all_matrices_competition[[i]] #pay attention to the sign of competition coefficients, positive sign means competition
  lambda[[i]] <- as.data.frame(all_intrinsic[[i]]) %>% dplyr::select( -lwr, -upr)
}

##use the sown species composition (sown composition in PaNDiv)
##Initialize empty lists for combination, comp, and richness
sown_composition<-read.table("data/Plot_Descriptions.txt", header=T, fill = TRUE)
plot_descriptions<-read.csv("data/Plot.info.csv")
sown_composition$Species <- as.factor(sown_composition$Species)
levels(sown_composition$Species)<-c("Am", "Ao", "As", "Be","Cj","Cb","Dg","Dc","Fr","Ga","Hp","Hs",
                                    "Hl", "Lp", "Pm", "Pt", "Pg", "Ra", "Sp", "To")#keep species same naming style
sown_composition_VIP<-sown_composition[sown_composition$VIP==1,]
##exclude As and Hs 
sown_composition_VIP1<-sown_composition_VIP[!sown_composition_VIP$Species %in% c("As","Hs"),]
table(sown_composition_VIP1$Species)
sown_composition_VIP1$treatment[sown_composition_VIP1$N==0&sown_composition_VIP1$F==0] <- "control"
sown_composition_VIP1$treatment[sown_composition_VIP1$N==1&sown_composition_VIP1$F==0] <- "nitrogen"
sown_composition_VIP1$treatment[sown_composition_VIP1$N==0&sown_composition_VIP1$F==1] <- "fungicide"
sown_composition_VIP1$treatment[sown_composition_VIP1$N==1&sown_composition_VIP1$F==1] <- "NF"

df_june <- sown_composition_VIP1
df_june$season <- "June"
df_august <- sown_composition_VIP1
df_august$season <- "August"
df_full <- rbind(df_june, df_august)
df_full$I <- with(df_full, ifelse(treatment == "control" & season == "June", 1,
                           ifelse(treatment == "nitrogen" & season == "June", 2,
                           ifelse(treatment == "fungicide" & season == "June", 3,
                           ifelse(treatment == "NF" & season == "June", 4,
                           ifelse(treatment == "control" & season == "August", 5,
                           ifelse(treatment == "nitrogen" & season == "August", 6,
                           ifelse(treatment == "fungicide" & season == "August", 7,
                           ifelse(treatment == "NF" & season == "August", 8, NA)))))))))


head(df_full)
table(df_full$season)  # Should be 1004 June, 1004 August
table(df_full$I)    
combination_ini <- list()
comp_ini <- list()
richness_ini <- list()

for (i in sort(unique(df_full$I))) {
  
  # Initialize empty dataframes for each index
  combination_ini[[i]] <- data.frame()
  comp_ini[[i]] <- data.frame()
  richness_ini[[i]] <- data.frame()
  
  for (j in unique(df_full[df_full$I==i,]$Plot_Nr)) {
    
    # Subset the data where plot matches and relative percentage >= 0
    subset_data <- subset(df_full[df_full$I==i,], Plot_Nr == j )
    richness_val <- nrow(subset_data)
    
    # Only consider cases where richness is greater than 1
    if (richness_val > 1) {
      
      # Append to combination list
      combination_ini[[i]] <- rbind(combination_ini[[i]], subset_data)
      
      # Store species from combination for comp
      comp_species <- subset(combination_ini[[i]], Plot_Nr == j)$Species
      comp_ini[[i]] <- rbind(comp_ini[[i]], data.frame(plot = j, species = comp_species))
      
      # Store the richness value
      richness_ini[[i]] <- rbind(richness_ini[[i]], data.frame(plot = j, richness = richness_val))
    }
  }
}

##the run is slow, therefore I use foreach
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

set.seed(123)
# Use foreach to parallelize the calculation across datasets
# computationally intensive, I provide the code but also load the intermediate results after this calculation
# load the intermediate results: observed coexistence mechanisms based on the sown composition
dataset_results_ini <- readRDS("intermediate_results/dataset_results_ini.rds")

# computationally intensive and suggest to skip
dataset_results_ini <- foreach(i = 1:8, .combine = 'c', .packages = c('mvtnorm','dplyr')) %dopar% {
  
  # Temporary list to store results for this dataset
  result_list <- list()
  
  # Loop over each plot for the current dataset
  for (j in unique(richness_ini[[i]]$plot)) {
    
    # Extract richness for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == j, 2]
    
    # Extract species combination for the current plot (as a character vector)
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == j, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- combination1[combination1 %in% rownames(alpha[[i]])]
    combination1 <- rownames(alpha[[i]])[rownames(alpha[[i]]) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha[[i]],          # Alpha matrix for dataset i
        intrinsic = lambda[[i]],     # Intrinsic matrix for dataset i
        n = n1,                      # Richness for the current plot
        combination = combination1,  # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
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
#saveRDS(dataset_results_ini, file = "intermediate_results/dataset_results_ini.rds")
#Stop the cluster after completion
#stopCluster(cl)

#calculate Omega_all_pairs separately (because in previous code I did not calculate it, redo everything requires too much time)
n_samples <- 10000L
dataset_names <- names(all_matrices_competition)

overlap_results_ini <- foreach(i = seq_along(alpha), .combine = "c") %dopar% {
  result_list <- list()
  A <- as.matrix(alpha[[i]])
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    j <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == j])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(A)) ||!all(requested_species %in% colnames(A))
    ) {
      stop(
        "Invalid species composition: dataset ", i, ", plot ", j
      )
    }
    combination1 <- rownames(A)[rownames(A) %in% requested_species]
    n1 <- length(combination1)
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == j]))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i,", plot ", j)
    }
    
    alpha_community <- A[combination1, combination1, drop = FALSE]
    
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))) {
      stop(
        "Invalid interaction matrix: dataset ", i, ", plot ", j)
    }
    
    seed_plot <- as.integer(123L + i * 100000L + plot_index)
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_,Omega_all = NA_real_,overlap = NA_real_)
      
    } else {set.seed(seed_plot)
            y <- compute_overlap(alpha_community,Nrand = n_samples)
    }
    result_list[[paste("Dataset", i, "Plot", j)]] <-
      data.frame(
        plot = as.character(j),
        treatment = dataset_names[i],
        season = sub("\\s.*$", "", dataset_names[i]),
        richness = n1,
        Omega_overlap = as.numeric(y$Omega),
        Omega_all_pairs = as.numeric(y$Omega_all),
        differential_new = as.numeric(y$Omega - y$Omega_all),
        overlap_new = as.numeric(y$overlap),
        overlap_n_samples = n_samples,
        overlap_seed = seed_plot,
        stringsAsFactors = FALSE
      )
  }
  result_list
}
overlap_results_ini1 <- dplyr::bind_rows(overlap_results_ini)

dataset_results_ini1 <- bind_rows(dataset_results_ini)
comb_str_ini <-dataset_results_ini1[,c(1:8, 15:19)]

comb_str_ini_1<-merge(comb_str_ini, overlap_results_ini1, , by=c("plot", "treatment"))
comb_str_ini_1$Omega<-as.numeric(comb_str_ini_1$Omega)
comb_str_ini_1$theta<-as.numeric(comb_str_ini_1$theta)
comb_str_ini_1$Omega1<-as.numeric(comb_str_ini_1$Omega1)
comb_str_ini_1$theta1<-as.numeric(comb_str_ini_1$theta1)
comb_str_ini_1$Omega2<-as.numeric(comb_str_ini_1$Omega2)
comb_str_ini_1$theta2<-as.numeric(comb_str_ini_1$theta2)
comb_str_ini_1$differential<-as.numeric(comb_str_ini_1$differential)
comb_str_ini_1$feasibility<-as.numeric(comb_str_ini_1$feasibility)
comb_str_ini_1$distance_to_exclusion<-as.numeric(comb_str_ini_1$distance_to_exclusion)
head(comb_str_ini_1)

######
##99 times of randomization for alpha and lambda (details see methods) 
##Next, calculation of coexistence mechanisms based on the randomized alpha and lambda  
##calculation was too slow, that why we decided to run for each alpha[[i]] and lambda[[i]] separately.
##Intermediate results after randomization were saved. So to save time and show what we have done, you could skip these randomization and load these intermediate results

# Step 1: Randomize alpha and lambda matrices (i.e., June C matrices) 99 times
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
num_iterations <- 99
name_vec <- colnames(alpha[[1]])
set.seed(123)

#load the saved intermediate results
randomized_matrices_new <- readRDS("intermediate_results/randomized_matrices_new_99.rds")

##only if you want to calculate it by yourself, otherwise load the intermediate results
randomized_matrices_new <- foreach(i = 1:8, .combine = 'list',.multicombine = TRUE, .packages = 'dplyr') %dopar% {
  rand_alpha_list_ad <- vector("list", num_iterations)
  rand_lambda_list_ad <- vector("list", num_iterations)
  
  for (j in 1:num_iterations) {
    # Create randomized alpha and lambda matrices
    alpha_rand_ad <- alpha[[i]]
    lambda_rand_ad <- lambda[[i]]
    
    # Randomize row and column names
    names_rand_ad <- sample(name_vec)
    rownames(alpha_rand_ad) <- names_rand_ad
    colnames(alpha_rand_ad) <- names_rand_ad#this makes sure that intras and inters are randomized within its pool
    rownames(lambda_rand_ad) <- names_rand_ad
    
    # Store randomized matrices
    rand_alpha_list_ad[[j]] <- alpha_rand_ad
    rand_lambda_list_ad[[j]] <- lambda_rand_ad
  }
  
  # Return lists of randomized matrices for the current dataset
  list(alpha = rand_alpha_list_ad, lambda = rand_lambda_list_ad)
}

#the calculation goes too slow, that why I decide to run for each set of alpha[[i]] and lambda[[i]] separately
#load the saved intermediate results or run the calculation
dataset_results_ini_rand1_ad <- readRDS("intermediate_results/dataset_results_ini_rand1_ad.rds")

##only if you want to calculate it by yourself, otherwise load the intermediate results
i=1
dataset_results_ini_rand1_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}
#saveRDS(dataset_results_ini_rand1_ad,"intermediate_results/dataset_results_ini_rand1_ad.rds")

#rerun compute_overlap() for the 99 iterations
overlap_results_ini_rand1 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop("Invalid species composition: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    # Match the species order in this randomized matrix.
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) || n_recorded != n1) {
      stop("Richness mismatch: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    alpha_community <- alpha_rand[combination1, combination1,drop = FALSE]
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i,", iteration ", j, ", plot ", k)
    }
    # Reproducible seed for dataset × iteration × plot.
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      y <- list(Omega = NA_real_,Omega_all = NA_real_,overlap = NA_real_)
      
    } else {set.seed(seed_plot)
      y <- compute_overlap(alpha_community,Nrand = n_samples)
    }
    plot_key <- paste("Plot", k)
    iteration_result_list[[plot_key]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand1,"intermediate_results/overlap_results_ini_rand1.rds")
columns_to_update <- c("Omega_overlap", "Omega_all_pairs", "differential", "overlap","randomization", "overlap_n_samples", "overlap_seed")

comb_str_ini_rand1_new <- combine_one_dataset(dataset_results_ini_rand1_ad, 
                                              overlap_results_ini_rand1, "June", "C")
comb_str_ini_rand1_new <-comb_str_ini_rand1_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand1_new)
summary(as.numeric(subset(comb_str_ini_rand1_new, plot %in% c(21, 117, 247, 290))$Omega))#for 20 species plot Omega and theta are exactly same

##
dataset_results_ini_rand2_ad <- readRDS("intermediate_results/dataset_results_ini_rand2_ad.rds")

i=2#the calculation goes too slow, that why I decide to run for each alpha[[i]] and lambda[[i]]
dataset_results_ini_rand2_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand2 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop(
        "Invalid species composition: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]
    ))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i,", iteration ", j, ", plot ", k
      )
    }
    
    alpha_community <- alpha_rand[combination1, combination1,drop = FALSE
    ]
    
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    # Reproducible seed for dataset × iteration × plot.
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index
    )
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_,Omega_all = NA_real_, overlap = NA_real_)
      
    } else {set.seed(seed_plot)
      y <- compute_overlap(alpha_community,Nrand = n_samples)}
    
    plot_key <- paste("Plot", k)
    
    iteration_result_list[[plot_key]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand2,"intermediate_results/overlap_results_ini_rand2.rds")

comb_str_ini_rand2_new <- combine_one_dataset(dataset_results_ini_rand2_ad,
                                              overlap_results_ini_rand2, "June", "N")
comb_str_ini_rand2_new <-comb_str_ini_rand2_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand2_new)
# Overwrite the intermediate results with the updated version.
#saveRDS(dataset_results_ini_rand2_ad,"intermediate_results/dataset_results_ini_rand2_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand2_new, plot %in% c(43,59,99,308))$Omega))

##
dataset_results_ini_rand3_ad <- readRDS("intermediate_results/dataset_results_ini_rand3_ad.rds")

i=3
dataset_results_ini_rand3_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
}
#saveRDS(overlap_results_ini_rand3,"intermediate_results/overlap_results_ini_rand3.rds")

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand3 <- foreach(j = seq_len(num_iterations),.packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop("Invalid species composition: dataset ", i,", iteration ", j, ", plot ", k)
    }
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]
    ))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1) {
      stop(
        "Richness mismatch: dataset ", i,", iteration ", j, ", plot ", k
      )
    }
    
    alpha_community <- alpha_rand[combination1, combination1, drop = FALSE]
    
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i,", iteration ", j, ", plot ", k
      )
    }
    
    # Reproducible seed for dataset × iteration × plot.
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_,Omega_all = NA_real_,overlap = NA_real_)
      
    } else {
      
      set.seed(seed_plot)
      y <- compute_overlap(alpha_community,Nrand = n_samples)
    }
    
    plot_key <- paste("Plot", k)
    
    iteration_result_list[[plot_key]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
comb_str_ini_rand3_new <- combine_one_dataset(dataset_results_ini_rand3_ad, overlap_results_ini_rand3, "June", "F")
comb_str_ini_rand3_new <-comb_str_ini_rand3_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand3_new)
# Overwrite the intermediate file with the updated results.
saveRDS(dataset_results_ini_rand3_ad, "intermediate_results/dataset_results_ini_rand3_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand3_new, plot==158)$theta))##all 20 plots are same

##
dataset_results_ini_rand4_ad <- readRDS("intermediate_results/dataset_results_ini_rand4_ad.rds")

i=4
dataset_results_ini_rand4_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand4 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop("Invalid species composition: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]))
    
    if (n1 < 2L || length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i, ", iteration ", j, ", plot ", k
      )
    }
    alpha_community <- alpha_rand[combination1, combination1, drop = FALSE]
    
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_, Omega_all = NA_real_, overlap = NA_real_)
      
    } else {
      set.seed(seed_plot)
      y <- compute_overlap(alpha_community, Nrand = n_samples
      )
    }
    iteration_result_list[[paste("Plot", k)]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand4,"intermediate_results/overlap_results_ini_rand4.rds")
comb_str_ini_rand4_new <- combine_one_dataset(dataset_results_ini_rand4_ad, overlap_results_ini_rand4, "June", "NF")
comb_str_ini_rand4_new <-comb_str_ini_rand4_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand4_new)

# Overwrite the intermediate file with the updated results.
saveRDS(dataset_results_ini_rand4_ad, "intermediate_results/dataset_results_ini_rand4_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand4_new, plot == 83)$Omega))#all the same

##
dataset_results_ini_rand5_ad <- readRDS("intermediate_results/dataset_results_ini_rand5_ad.rds")

i=5
dataset_results_ini_rand5_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand5 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop(
        "Invalid species composition: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    
    alpha_community <- alpha_rand[combination1, combination1, drop = FALSE]
    
    if (!is.numeric(alpha_community) || any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_, Omega_all = NA_real_, overlap = NA_real_)
      
    } else {
      
      set.seed(seed_plot)
      y <- compute_overlap(alpha_community, Nrand = n_samples)
    }
    iteration_result_list[[paste("Plot", k)]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand5,"intermediate_results/overlap_results_ini_rand5.rds")

comb_str_ini_rand5_new <- combine_one_dataset(dataset_results_ini_rand5_ad, overlap_results_ini_rand5, "August", "C")
comb_str_ini_rand5_new <-comb_str_ini_rand5_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand5_new)

# Overwrite the intermediate file with the updated results.
saveRDS(dataset_results_ini_rand5_ad, "intermediate_results/dataset_results_ini_rand5_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand5_new, plot %in% c(21, 117, 247, 290))$Omega))

##
dataset_results_ini_rand6_ad <- readRDS("intermediate_results/dataset_results_ini_rand6_ad.rds")

i=6
dataset_results_ini_rand6_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand6 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) || anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) || !all(requested_species %in% colnames(alpha_rand))
    ) {
      stop("Invalid species composition: dataset ", i, ", iteration ", j, ", plot ", k)
    }
    
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]))
    
    if (n1 < 2L ||length(n_recorded) != 1L ||is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i,", iteration ", j, ", plot ", k
      )
    }
    
    alpha_community <- alpha_rand[combination1, combination1,drop = FALSE]
    
    if (!is.numeric(alpha_community) ||any(!is.finite(alpha_community))) {
      stop(
        "Invalid matrix: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      
      y <- list(Omega = NA_real_,Omega_all = NA_real_, overlap = NA_real_)
      
    } else {
      
      set.seed(seed_plot)
      y <- compute_overlap(alpha_community,Nrand = n_samples)
    }
    
    iteration_result_list[[paste("Plot", k)]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand6,"intermediate_results/overlap_results_ini_rand6.rds")
comb_str_ini_rand6_new <- combine_one_dataset(dataset_results_ini_rand6_ad, overlap_results_ini_rand6, "August", "N")
comb_str_ini_rand6_new <-comb_str_ini_rand6_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand6_new)
# Overwrite the intermediate file with the updated results.
saveRDS(dataset_results_ini_rand6_ad,"intermediate_results/dataset_results_ini_rand6_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand6_new, plot %in% c(59, 99, 308) )$Omega))#all the same

##
dataset_results_ini_rand7_ad <- readRDS("intermediate_results/dataset_results_ini_rand7_ad.rds")

i=7
dataset_results_ini_rand7_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}

# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand7 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (anyNA(requested_species) ||anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||!all(requested_species %in% colnames(alpha_rand))
    ) {
      stop("Invalid species composition: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% requested_species]
    n1 <- length(combination1)
    
    n_recorded <- as.numeric(as.character(richness_ini[[i]]$richness[richness_ini[[i]]$plot == k]))
    
    if (n1 < 2L ||length(n_recorded) != 1L || is.na(n_recorded) ||n_recorded != n1
    ) {
      stop("Richness mismatch: dataset ", i,", iteration ", j, ", plot ", k)
    }
    
    alpha_community <- alpha_rand[combination1, combination1, drop = FALSE]
    
    if (
      !is.numeric(alpha_community) ||
      any(!is.finite(alpha_community))
    ) {
      stop("Invalid matrix: dataset ", i, ", iteration ", j, ", plot ", k
      )
    }
    # Reproducible seed for dataset × iteration × plot.
    seed_plot <- as.integer(i * 10000000L + j * 10000L + plot_index)
    
    if (all(alpha_community == 0)) {
      y <- list(Omega = NA_real_, Omega_all = NA_real_, overlap = NA_real_)
      
    } else {
      
      set.seed(seed_plot)
      y <- compute_overlap(alpha_community, Nrand = n_samples)
    }
    
    iteration_result_list[[paste("Plot", k)]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand7,"intermediate_results/overlap_results_ini_rand7.rds")
comb_str_ini_rand7_new <- combine_one_dataset(dataset_results_ini_rand7_ad, overlap_results_ini_rand7, "August", "F")
comb_str_ini_rand7_new <-comb_str_ini_rand7_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand7_new)

# Overwrite the intermediate file with the updated results.
saveRDS(dataset_results_ini_rand7_ad, "intermediate_results/dataset_results_ini_rand7_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand7_new, plot %in% c(158, 171, 295, 328) )$theta))#

##
dataset_results_ini_rand8_ad <- readRDS("intermediate_results/dataset_results_ini_rand8_ad.rds")

i=8
dataset_results_ini_rand8_ad <- foreach(j = 1:num_iterations, .combine = 'list', .multicombine = TRUE, .packages = c('mvtnorm', 'dplyr'))  %dopar% {
  
  if (j %% 10 == 0) {
    cat(sprintf("Dataset %d, Iteration %d completed\n", i, j), file = log_file, append = TRUE)
  }
  
  # Create a copy of the original alpha and lambda matrices
  alpha_rand <- randomized_matrices_new[[i]]$alpha[[j]]
  lambda_rand <- randomized_matrices_new[[i]]$lambda[[j]]
  
  # List to store results for this iteration
  iteration_result_list <- list()
  iteration_alpha_results <- list()
  iteration_lambda_results <- list()
  
  # Loop over each plot for the current dataset
  for (k in unique(richness_ini[[i]]$plot)) {
    # Extract richness and composition for the current plot
    n1 <- richness_ini[[i]][richness_ini[[i]]$plot == k, 2]
    combination1 <- as.character(comp_ini[[i]][comp_ini[[i]]$plot == k, 2])
    # ensure that the species in the subset (4, 8 and 20 species) are matched and always ordered according to alpha_rand
    combination1 <- rownames(alpha_rand)[rownames(alpha_rand) %in% combination1]
    
    # Proceed only if combination1 is not empty
    if (length(combination1) == 0) {
      warning(paste("No valid species found for plot", k, "of dataset", i, ". Skipping this plot."))
      next
    }
    # Run the structural_coex function for the current plot if combination is valid
    if (length(combination1) > 0) {
        result1 <- structural_coex_new3(
        alpha = alpha_rand,         # Alpha matrix for dataset i
        intrinsic = lambda_rand,    # Intrinsic matrix for dataset i
        n = n1,                     # Richness for the current plot
        combination = combination1, # Species combination for the current plot
        list_names = names(all_matrices_competition)[i],  # Dataset name
        plot = k                    # The specific plot being analyzed
      )
      # Store alpha and lambda matrices and results
      iteration_alpha_results[[paste("Plot", k)]] <- alpha_rand[combination1, combination1, drop = FALSE]
      iteration_lambda_results[[paste("Plot", k)]] <- lambda_rand[combination1, , drop = FALSE]
      iteration_result_list[[paste("Plot", k)]] <- result1
      
    } 
  }
  list(
    result_list = iteration_result_list,
    alpha_results = iteration_alpha_results,
    lambda_results = iteration_lambda_results
  )
  
}
# Recalculate ONLY the overlap-related metrics.
overlap_results_ini_rand8 <- foreach(j = seq_len(num_iterations), .packages = "dplyr") %dopar% {
  
  alpha_rand <- as.matrix(randomized_matrices_new[[i]]$alpha[[j]])
  iteration_result_list <- list()
  plot_ids <- unique(richness_ini[[i]]$plot)
  
  for (plot_index in seq_along(plot_ids)) {
    
    k <- plot_ids[plot_index]
    requested_species <- as.character(comp_ini[[i]]$species[comp_ini[[i]]$plot == k])
    
    if (
      anyNA(requested_species) ||
      anyDuplicated(requested_species) ||
      !all(requested_species %in% rownames(alpha_rand)) ||
      !all(requested_species %in% colnames(alpha_rand))
    ) {
      stop(
        "Invalid species composition: dataset ", i,
        ", iteration ", j, ", plot ", k
      )
    }
    
    combination1 <- rownames(alpha_rand)[
      rownames(alpha_rand) %in% requested_species
    ]
    
    n1 <- length(combination1)
    
    n_recorded <- as.numeric(as.character(
      richness_ini[[i]]$richness[
        richness_ini[[i]]$plot == k
      ]
    ))
    
    if (
      n1 < 2L ||
      length(n_recorded) != 1L ||
      is.na(n_recorded) ||
      n_recorded != n1
    ) {
      stop(
        "Richness mismatch: dataset ", i,
        ", iteration ", j, ", plot ", k
      )
    }
    
    alpha_community <- alpha_rand[
      combination1, combination1,
      drop = FALSE
    ]
    
    if (
      !is.numeric(alpha_community) ||
      any(!is.finite(alpha_community))
    ) {
      stop(
        "Invalid matrix: dataset ", i,
        ", iteration ", j, ", plot ", k
      )
    }
    
    # Reproducible seed for dataset × iteration × plot.
    seed_plot <- as.integer(
      i * 10000000L + j * 10000L + plot_index
    )
    
    if (all(alpha_community == 0)) {
      
      y <- list(
        Omega = NA_real_,
        Omega_all = NA_real_,
        overlap = NA_real_
      )
      
    } else {
      
      set.seed(seed_plot)
      
      y <- compute_overlap(
        alpha_community,
        Nrand = n_samples
      )
    }
    
    iteration_result_list[[paste("Plot", k)]] <- data.frame(
      plot = as.character(k),
      treatment = names(all_matrices_competition)[i],
      randomization = j,
      Omega_overlap = as.numeric(y$Omega),
      Omega_all_pairs = as.numeric(y$Omega_all),
      differential = as.numeric(y$Omega - y$Omega_all),
      overlap = as.numeric(y$overlap),
      overlap_n_samples = n_samples,
      overlap_seed = seed_plot,
      stringsAsFactors = FALSE
    )
  }
  
  iteration_result_list
}
#saveRDS(overlap_results_ini_rand8,"intermediate_results/overlap_results_ini_rand8.rds")
comb_str_ini_rand8_new <- combine_one_dataset(dataset_results_ini_rand8_ad, overlap_results_ini_rand8, "August", "NF")
comb_str_ini_rand8_new <-comb_str_ini_rand8_new[,c("plot", "Omega", "theta", "differential", "feasibility", "distance_to_exclusion", "differential_old", "Omega_overlap", "Omega_all_pairs", "season", "Treatment", "year")]
head(comb_str_ini_rand8_new)

# Overwrite the intermediate file with the updated results.
#saveRDS(dataset_results_ini_rand8_ad, "intermediate_results/dataset_results_ini_rand8_ad.rds")
summary(as.numeric(subset(comb_str_ini_rand8_new, plot==233)$Omega))#all the same
##
comb_str_ini_rand_new <- dplyr::bind_rows(
  comb_str_ini_rand1_new,
  comb_str_ini_rand2_new,
  comb_str_ini_rand3_new,
  comb_str_ini_rand4_new,
  comb_str_ini_rand5_new,
  comb_str_ini_rand6_new,
  comb_str_ini_rand7_new,
  comb_str_ini_rand8_new
)

metric_names <- c("Omega", "theta", "distance_to_exclusion", "differential", "differential_old", "Omega_all_pairs")

comb_str_ini_rand_new <- comb_str_ini_rand_new %>%
  dplyr::mutate(plot = as.character(plot),
    dplyr::across(dplyr::all_of(metric_names),
      ~ as.numeric(as.character(.x))
    )
  )

mean_values_rand_ini_new <- comb_str_ini_rand_new %>%
  dplyr::group_by(plot, season) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(metric_names),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

names(mean_values_rand_ini_new)[3:8] <- c(
  "mean_omega_rand", "mean_theta_rand",
  "mean_distance_to_exclusion", "mean_differential",
  "mean_differential_old", "mean_Omega_all_pairs"
)

sd_values_rand_ini_new <- comb_str_ini_rand_new %>%
  dplyr::group_by(plot, season) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(metric_names),
      ~ sd(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

names(sd_values_rand_ini_new)[3:8] <- c(
  "sd_omega_rand", "sd_theta_rand",
  "sd_distance_to_exclusion", "sd_differential",
  "sd_differential_old", "sd_Omega_all_pairs"
)

names(plot_descriptions)[2] <- "plot"

comb_str_ini_1$plot<-as.character(comb_str_ini_1$plot)
merged_data_ini_rand_new <- comb_str_ini_1 %>%left_join(mean_values_rand_ini_new, by = c("plot","season"))%>%left_join(sd_values_rand_ini_new, by = c("plot","season"))
combined_sown_coex_new<-merge(merged_data_ini_rand_new, plot_descriptions, all.x = F, by = c("plot") )
View(subset(combined_sown_coex_new, SD%in%c(4,8)))
summary(subset(combined_sown_coex_new, SD%in%c(4,8))$Omega_all_pairs)

##sta.Omega: predicted niche differences without facilitation; 
##sta.theta: predicted competitive ability differences without facilitation; 
##sta.differential: predicted indirect interactions without facilitation; 
combined_sown_coex_new$sta.Omega<-(combined_sown_coex_new$Omega-combined_sown_coex_new$mean_omega_rand)/combined_sown_coex_new$sd_omega_rand
combined_sown_coex_new$sta.theta<-(combined_sown_coex_new$theta-combined_sown_coex_new$mean_theta_rand)/combined_sown_coex_new$sd_theta_rand
combined_sown_coex_new$sta.differential_new<-(combined_sown_coex_new$differential_new-combined_sown_coex_new$mean_differential)/combined_sown_coex_new$sd_differential
combined_sown_coex_new$sta.differential<-(combined_sown_coex_new$differential-combined_sown_coex_new$mean_differential_old)/combined_sown_coex_new$sd_differential_old

combined_sown_coex_new <- combined_sown_coex_new %>%
  dplyr::mutate(
    zero_sd_Omega_all_pairs =!is.na(sd_Omega_all_pairs) & sd_Omega_all_pairs == 0,
    
    sta.Omega_all_pairs = (Omega_all_pairs - mean_Omega_all_pairs
    ) / dplyr::if_else(
      is.finite(sd_Omega_all_pairs) & sd_Omega_all_pairs > 0,
      sd_Omega_all_pairs,
      NA_real_
    )
  )

combined_sown_coex_new %>%
  dplyr::mutate(
    observed_zero = Omega_all_pairs == 0,
    null_mean_zero = mean_Omega_all_pairs == 0,
    null_sd_zero = sd_Omega_all_pairs == 0,
    standardized_missing = is.na(sta.Omega_all_pairs)
  ) %>% dplyr::count(
    observed_zero,
    null_mean_zero,
    null_sd_zero,
    standardized_missing,
    name = "n_rows"
  ) %>% print()

combined_sown_coex_new <- combined_sown_coex_new %>%
  dplyr::mutate(
    zero_null_match = dplyr::coalesce(
      Omega_all_pairs == 0 &
        mean_Omega_all_pairs == 0 &
        sd_Omega_all_pairs == 0,
      FALSE
    ),
    sta.Omega_all_pairs_zero = dplyr::if_else(
      zero_null_match,
      0,
      sta.Omega_all_pairs
    )
  )
combined_sown_coex_new <- combined_sown_coex_new %>% dplyr::mutate(delta.Omega_all_pairs = Omega_all_pairs - mean_Omega_all_pairs)