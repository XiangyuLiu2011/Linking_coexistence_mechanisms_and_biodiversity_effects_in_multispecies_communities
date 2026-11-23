
##functions for normalized structural niche difference and fitness difference
# load necessary packages
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
library(dplyr)
library(anisoFun)

# Make the internal helper available to the custom theta1() function.
check_A_int <- getFromNamespace("check_A_int", "anisoFun")
check_center <- getFromNamespace("check_center", "anisoFun")
incenter_inradius_isoprob_calculation<- getFromNamespace("incenter_inradius_isoprob_calculation", "anisoFun")
cone_vertices_director_vertices<- getFromNamespace("cone_vertices_director_vertices", "anisoFun")
# loading function from anisoFun package manually
# Set the path to the folder containing the .R files
folder_path <- "D:/Postdoc/Linking_coexistence_mechanisms_and_biodiversity_effects/coexistence_BEF_code/scripts/code distance to exclusion"
# Get a list of all .R files in the folder
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)

# Load each file using source
lapply(files, source)

#function: change facilitation to very low competition in a competition matrix
rm_facilitation<- function (all) {
  seq<-runif(n= length(all[all>0]), min = min(abs(all))/10, max = min(abs(all))/9)
  seq<- -(seq)
  all[all>0]<-seq
  return(all)
}

# structural niche difference (output on a log scale)
Omega <- function(alpha) {
  n <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha, tol = 1e-40) # tolerance modified due to singularities
  d <- pmvnorm(lower = rep(0, n), upper = rep(Inf, n), mean = rep(0, n), sigma = Sigma,  seed = 123)##seed part added by Xiangyu
  out <- log10(d[1]) + n * log10(2)
  return(out)
}
# function that computes the normalized feasibility from an interaction matrix
# inputs: alpha = interaction matrix
# output: out = the normalized feasibility
# The range is [0,0.5]
OmegaNorm <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  
  # Function to check if the matrix is invertible
  f <- function(m) {
    result <- try(solve(t(m) %*% m), silent = TRUE)
    return(!inherits(result, "try-error")) # # Ensure a single logical value is returned
  }
  
  if (!f(alpha)) {
    return(0)
  } else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

# structural fitness difference (in degree)
theta <- function(alpha, r) {
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c * r) / (sqrt(sum(r^2)) * sqrt(sum(r_c^2)))) * 180 / pi
  # print(r_c*r)
  return(out)
}

##structural fitness difference (in degree) with normalization, function received from Oscar.
thetaNorm <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r, na.rm = TRUE)/(sqrt(sum(r^2, na.rm = TRUE))*sqrt(sum(r_c^2, na.rm = TRUE))))*90/pi
  return(out)
}

#vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}

#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha, r) > 0)
  return(out)
}

#test which pairs in a system (alpha and r) are feasible (output 1 = feasible, 0 = not feasible)
test_feasibility_pairs <- function(alpha,r){
  n <- length(r)
  c <- combn(n,2)
  nc <- dim(c)[2]
  f <- rep(NA,nc)
  for (i in 1:nc){
    f[i] <- prod(solve(alpha[c[,i],c[,i]],r[c[,i]])>0)
  }
  out <- list(pairs = c, feasibility = f)
  return(out)
}

#compute the feasiblity domain, the feasibility domain of all pairs, and their overlap (Nrand = number of randomization)
compute_overlap <- function(alpha,Nrand){
  
  n <- dim(alpha)[1]
  
  counter_f <- 0
  counter_overlap <- 0
  counter_all <- 0
  
  for (i in 1:Nrand){
    
    r_rand <- abs(rnorm(n))
    r_rand <- r_rand/sqrt(sum(r_rand^2))
    
    f1 <- test_feasibility(alpha,r_rand)  
    f2 <- test_feasibility_pairs(alpha,r_rand)$feasibility  
    
    counter_f <- counter_f + f1
    counter_all <- counter_all + prod(f2)
    counter_overlap <- counter_overlap + f1*prod(f2)
    
  }
  
  Omega <- counter_f/Nrand
  Omega_all <- counter_all/Nrand
  overlap <- counter_overlap/Nrand
  
  out <- list(Omega = Omega, Omega_all = Omega_all, overlap = overlap)
  return(out)
  
}

##I adapt normalization of structural ND and FD from the methods of Luo et al. 2022 PNAS
##The range is [0,1]
OmegaNorm2 <- function(alpha) {
    S <- nrow(alpha)
    omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- (d[1]* 2^S)^(1 / (S-1))
    return(out)
  }
  
  # Function to check if the matrix is invertible
  f <- function(m) {
    result <- try(solve(t(m) %*% m), silent = TRUE)
    return(!inherits(result, "try-error")) # # Ensure a single logical value is returned
  }
  
  if (!f(alpha)) {
    return(0)
  } else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

r_centroid2 <- function(alpha) {
  n <- nrow(alpha)  # Number of species
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))  # Scaling matrix
  alpha_n <- alpha %*% D  # Normalized alpha
  r_c <- rowSums(alpha_n) / n  # Centroid calculation
  return(r_c)  # Return as a vector
}

##Both OmegaNorm2 and thetaNorm2 are translated based the matlab code from Luo et al. 2022
thetaNorm2 <- function(alpha, lambda) {
  S <- nrow(alpha)  # Number of species
  
  # Calculate RR = (c * lambda_M - e)
  c <- diag(alpha)  # Extract c from the diagonal of alpha
  e <- rowSums(alpha)  # Sum of rows to get e
  RR <- c * lambda - e
  
  # Calculate r_c using the custom r_centroid function
  r_c <- r_centroid2(alpha)
  
  # Avoid division by zero by checking the norms of RR and r_c
  norm_RR <- sqrt(sum(RR^2))
  norm_r_c <- sqrt(sum(r_c^2))
  
  if (norm_RR == 0 || norm_r_c == 0) {
    warning("One of the norms is zero, returning NA.")
    return(NA)
  }
  
  # Calculate theta(S)
  theta_S <- acos(sum(RR * r_c) / (norm_RR * norm_r_c))
  
  # Calculate theta1(S) = theta(S) * sqrt(S)
  theta1_S <- theta_S * sqrt(S)
  
  # Return the theta1 value
  return(theta1_S)
}

#combination <- Filter(function(x) nrow(x) > 0, combination)
#comp <- Filter(function(x) nrow(x) > 0, comp)
#richness <- Filter(function(x) nrow(x) > 0, richness)

##Normalized theta1 function received from Oscar, based on anisoFun, however, the output are not the same compare to thetaNrom2
##I understand that this is the distance from the incenter to all the edges of the FD (along the surface of the unit ball) multiply sqrt(number of species)
##this function also received from Oscar.
theta1 <- function(A_int){
  
  center <- incenter_inradius_isoprob_calculation(A_int)[[1]] # the intrinsic growth rates of the feasibility domain's incenter  
  
  test_result_A_int <- check_A_int(A_int)
  test_result_center <- check_center(center, A_int)
  
  if(test_result_A_int[[1]]==T & test_result_center[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    if(test_result_center[[2]] == "Turning the initial growth rates into a unit vector."){
      
      center <- center/sqrt(sum(center*center))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    dimensions <- ncol(A_int)
    curves_data <- cone_vertices_director_vertices(A_int)
    
    col_number <- ncol(curves_data)
    
    dir_ver_mat <- as.matrix(curves_data[,(col_number-(dimensions-1)):col_number])
    
    # we check the direction of the director vectors
    correction_direction <- NULL
    
    for(i in 1: dimensions){
      
      if(all((dir_ver_mat %*% as.numeric(curves_data[i,(1:dimensions)]) %>% round(12))>=0)){
        correction_direction <- c(correction_direction,1)
      }else{
        correction_direction <- c(correction_direction,-1)
      }
      
    }
    
    acr_distance <- NULL
    
    for(i in 1:nrow(curves_data)){
      
      normal_vector_i <- correction_direction[i] * curves_data[i,c((col_number-(dimensions-1)):col_number)]
      sin_theta <- sum(normal_vector_i*center)
      acr_distance <- c(acr_distance,asin(sin_theta))
      
    }
    
    if(length(unique(acr_distance)) == 1){
      
      result <- abs(acr_distance[1])*sqrt(dimensions)
      
    }else{
      
      result <- tibble(species = colnames(A_int),
                       arc_distance = abs(acr_distance)*sqrt(dimensions))
    }
    
    return(result)
    
    
  }else{
    
    if(test_result_A_int[[2]]!=""){
      error_messages <- paste0(test_result_A_int[[2]],"\n",
                               test_result_center[[2]],"\n")
    }else{
      error_messages <- paste0(test_result_center[[2]],"\n")
    }
    
    cat(error_messages)
    
  }
  
}

##
structural_coex_new <- function(alpha, intrinsic, n, combination, list_names, plot) {
  
  pair_names <- apply(combn(n, 2), 2, paste, collapse = "_") #name all the possible pairs
  
  # Directly use the species combination for this plot
  combos <- as.matrix(combination)
  colnames(combos) <- NULL
  
  # Initialize results matrix
  results_combos <- matrix(nrow = 1, ncol = (5 + ncol(combn(n, 2)) + 2))
  row.names(results_combos) <- paste(combination, collapse = "_") # Combo name for this plot
  colnames(results_combos) <- c("Omega", "theta", "differential", "overlap", "feasibility", pair_names, "coex_rate", "treatment")
  
  # Set up the alpha and intrinsic matrices for the selected species in this plot
  zeroes <- matrix(data = 0, nrow = n, ncol = n)
  intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combination))
  alpha2 <- as.matrix(alpha[combination, combination])
  
  # Check for zero matrix condition
  if (isFALSE(FALSE %in% (alpha2 == zeroes))) {
    results_combos[1, ] <- NA
  } else {
    # omega
    results_combos[1, 1] <- OmegaNorm(alpha2)
    # theta
    results_combos[1, 2] <- thetaNorm(alpha2, intrinsic2)
    # differential and overlap
    y <- compute_overlap(alpha2, 10000)
    results_combos[1, 3] <- y$Omega - y$Omega_all
    results_combos[1, 4] <- y$overlap
    # feasibility (all)
    results_combos[1, 5] <- test_feasibility(alpha2, intrinsic2)
    x <- test_feasibility_pairs(alpha2, intrinsic2)
    # feasibility (pair by pair)
    results_combos[1, 6:(5 + ncol(combn(n, 2)))] <- x$feasibility
    # coexistence rate
    results_combos[1, ncol(results_combos)-1] <- sum(x$feasibility) / ncol(combn(n, 2))
    # treatment (list name or other identifier)
    results_combos[1, ncol(results_combos)] <- list_names
  }
  
  results_combos <- as.data.frame(results_combos)
  
  # Include plot name as part of the output
  results_combos$plot <- plot
  
  return(results_combos)
}


##function theta1 and sp_distance_to_exclusion cannot handle the 2 species plots.
structural_coex_new2 <- function(alpha, intrinsic, n, combination, list_names, plot) {
  
  pair_names <- apply(combn(n, 2), 2, paste, collapse = "_") #name all the possible pairs
  
  # Directly use the species combination for this plot
  combos <- as.matrix(combination)
  colnames(combos) <- NULL
  
  # Initialize results matrix
  results_combos <- matrix(nrow = 1, ncol = (5 + ncol(combn(n, 2)) + 3))
  row.names(results_combos) <- paste(combination, collapse = "_") # Combo name for this plot
  colnames(results_combos) <- c("Omega", "theta", "differential", "overlap", "feasibility", pair_names, "coex_rate", "treatment", "distance_to_exclusion")
  
  # Set up the alpha and intrinsic matrices for the selected species in this plot
  zeroes <- matrix(data = 0, nrow = n, ncol = n)
  intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combination))
  alpha2 <- as.matrix(alpha[combination, combination])
  
  # Check for zero matrix condition
  if (isFALSE(FALSE %in% (alpha2 == zeroes))) {
    results_combos[1, ] <- NA
  } else {
    # omega
    results_combos[1, 1] <- OmegaNorm2(alpha2)
    # theta
    results_combos[1, 2] <- theta1(alpha2)$arc_distance[1]
    # differential and overlap
    y <- compute_overlap(alpha2, 10000)
    results_combos[1, 3] <- y$Omega - y$Omega_all
    results_combos[1, 4] <- y$overlap
    # feasibility (all)
    results_combos[1, 5] <- test_feasibility(alpha2, intrinsic2)
    x <- test_feasibility_pairs(alpha2, intrinsic2)
    # feasibility (pair by pair)
    results_combos[1, 6:(5 + ncol(combn(n, 2)))] <- x$feasibility
    # coexistence rate
    results_combos[1, ncol(results_combos)-1] <- sum(x$feasibility) / ncol(combn(n, 2))
    # treatment (list name or other identifier)
    results_combos[1, ncol(results_combos)] <- list_names
    min_distance <- sp_distance_to_exclusion(intrinsic2, -alpha2)
    results_combos[1, ncol(results_combos)] <- min(min_distance$arc_distance)
  }
  
  results_combos <- as.data.frame(results_combos)
  
  # Include plot name as part of the output
  results_combos$plot <- plot
  
  return(results_combos)
}

## in thisstructural_coex_new3, I put omega1 (OmegaNorm) (Song et al. 2018, range [0,0.5]) and omega2 (OmegaNorm2) (Luo et al. 2022, range [0,1])
## and I put theta1 (thetaNorm) (theta1 is from Song et al.) and theta2 (theta1) ((Luo et al. 2022) as well.
## theta and distance to exclusion of 2 species plot will be assign to NA
structural_coex_new3 <- function(alpha, intrinsic, n, combination, list_names, plot) {
  
  pair_names <- apply(combn(n, 2), 2, paste, collapse = "_") #name all the possible pairs
  
  # Directly use the species combination for this plot
  combos <- as.matrix(combination)
  colnames(combos) <- NULL
  
  # Initialize results matrix
  results_combos <- matrix(nrow = 1, ncol = (8 + ncol(combn(n, 2)) + 4))
  row.names(results_combos) <- paste(combination, collapse = "_") # Combo name for this plot
  colnames(results_combos) <- c("Omega", "Omega1", "Omega2","theta", "theta1", "differential", "overlap", "feasibility", pair_names, "coex_rate", "treatment", "theta2", "distance_to_exclusion")
  
  # Set up the alpha and intrinsic matrices for the selected species in this plot
  zeroes <- matrix(data = 0, nrow = n, ncol = n)
  intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combination))
  alpha2 <- as.matrix(alpha[combination, combination])
  
  # Check for zero matrix condition
  if (isFALSE(FALSE %in% (alpha2 == zeroes))) {
    results_combos[1, ] <- NA
  } else {
    # omega without normalization
    results_combos[1, 1] <- 10^Omega(alpha2) 
    # omega1 normalized based on Song et al. 2018 (0~0.5)
    results_combos[1, 2] <- OmegaNorm(alpha2)
    # omega2 normalized based on Luo et al. 2022 (0~1)
    results_combos[1, 3] <- OmegaNorm2(alpha2)
    # theta without normalization
    results_combos[1, 4] <- theta(alpha2, intrinsic2)
    # theta normalized based on Song et al. 2018 
    results_combos[1, 5] <- thetaNorm(alpha2, intrinsic2)
    # differential and overlap
    y <- compute_overlap(alpha2, 10000)
    results_combos[1, 6] <- y$Omega - y$Omega_all
    results_combos[1, 7] <- y$overlap
    # feasibility (all)
    results_combos[1, 8] <- test_feasibility(alpha2, intrinsic2)
    x <- test_feasibility_pairs(alpha2, intrinsic2)
    # feasibility (pair by pair)
    results_combos[1, 9:(8 + ncol(combn(n, 2)))] <- x$feasibility
    # coexistence rate
    results_combos[1, ncol(results_combos)-3] <- sum(x$feasibility) / ncol(combn(n, 2))
    # treatment (list name or other identifier)
    results_combos[1, ncol(results_combos)-2] <- list_names
    
    # Handle special case for two species
    if (n == 2) {
      # Provide a simple fallback for two-species communities, as arc distance might not apply well
      results_combos[1, ncol(results_combos)-1] <- NA
      results_combos[1, ncol(results_combos)] <- NA  # Or provide a different simple metric for 2D cases
    } else {
      # For n > 2, use the usual sp_distance_to_exclusion calculation
      # theta 2
      results_combos[1, ncol(results_combos)-1] <- theta1(alpha2)$arc_distance[1]
      min_distance <- sp_distance_to_exclusion(intrinsic2, -alpha2)
      results_combos[1, ncol(results_combos)] <- min(min_distance$arc_distance)
    }
  }
  
  results_combos <- as.data.frame(results_combos)
  
  # Include plot name as part of the output
  results_combos$plot <- plot
  
  return(results_combos)
}


#I save it as toolbox.R

## function to predict from multigroup
## provide model fit, pred.data (names same as in model), level of the x variables in pred.data (as a vector with length = ncol (pred.data)), intercept.level = FALSE to use 0 as intercept, otherwise level 1 or 2
predict.mgroup <- function(model_fit, pred.data, yvar, level, intercept.level =FALSE){
  
  pp <- parameterEstimates(model_fit)
  ngroups <- lavInspect(model_fit,"ngroups")
  preds <- names(pred.data[[1]])
  
  ## select the required coefficients based on variables in pred.data and specified yvar
  coef.mat <- pp[pp$op=="~"&pp$rhs%in%preds&pp$lhs%in%yvar,]
  ## subset correct level
  coef.mat <- coef.mat[coef.mat$level%in%level,]
  
  # split into list based on groups
  coef.list <- split(coef.mat, coef.mat[,"group"])
  
  ## if discount intercepts use 0
  if(intercept.level==FALSE){
    intercepts <- rep(0, ngroups)
  }
  
  ## otherwise take from model
  else{
    ii <- pp[pp$op=="~1"&pp$lhs%in%yvar&pp$level%in%intercept.level,]
    intercepts <- ii[,"est"]
    ii.groups <- ii[,"group"]
  }
  ## do the prediction
  pred.list <- list()
  
  for(i in 1:ngroups){
    xx <- coef.list[[i]][,"est"]
    cc <- coef.list[[i]][,"rhs"]
    ## sort to same order as pred data columns
    xx2 <- xx[match(names(pred.data[[i]]),cc)]
    
    xy <- t(t(pred.data[[i]])*xx2)
    xyi <- cbind(intercepts[i], xy)
    yy <- apply(xyi, 1, sum)
    predd <- cbind(pred.data[[i]], yy)
    names(predd)[which(names(predd)=="yy")] <- yvar
    pred.list[[i]] <- predd 
  }
  
  return(pred.list)
  
}

library(feasibilityR)
norm1 <- function(a) {
  if (is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sum(x))
  } else {
    a <- a / sum(a)
  }
  return(a)
}

norm2 <- function(a) {
  if (is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    a <- a / sqrt(sum(a^2))
  }
  return(a)
}

feasibility_resistance_full <- function(matA, r, norm = "l2", nsample = 100, all = FALSE) {
  simplex_sampling <- function(m, n) {
    r <- list()
    for (j in 1:m) {
      dist <- c(sort(runif(n - 1, 0, 1)), 1)
      r[[j]] <- c(dist[1], diff(dist))
    }
    return(r)
  }
  euclidean_distance <- function(a, b) {
    sqrt(sum((a - b)^2))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }
  
  vertices <- combn(seq_len(nrow(matA)), nrow(matA) - 1, simplify = FALSE)
  distances_list <- list()
  for (i in seq_len(length(vertices))) {
    vertex <- vertices[[i]]
    if (norm == "l1") {
      r <- norm1(r); matA <- norm1(matA)
      distances_list[[i]] <- 1:nsample %>%
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA) - 1))
          border_point <- c(matA[, vertex] %*% t)
          euclidean_distance(r, border_point)
        })
    }
    if (norm == "l2") {
      r <- norm2(r); matA <- norm2(matA)
      distances_list[[i]] <- 1:nsample %>%
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA) - 1))
          border_point <- c(matA[, vertex] %*% t)
          border_point_norm <- border_point / sqrt(sum(border_point^2))
          arc_length(r, border_point_norm)
        })
    }
  }
  names(distances_list) <- unlist(lapply(vertices, paste, collapse = "_"))
  if (all) {
    return(distances_list)
  } else {
    return(sapply(distances_list, min))
  }
}

plot.res.sem <- function(fit_data, pred_term, residual_list, cols, alpha, xlab, ylab, xlim, ylim){
  
  plot(1, xlim=xlim, ylim =ylim, xlab =xlab, ylab =ylab, col ="white", las =1)
  
  cols2 <- alpha(cols, alpha=alpha) ## for the points
  
  for(i in 1:4){
    
    x <- fit_data[[i]][,pred_term]
    
    points(x, residual_list[[1]][[i]], col = cols2[i], pch =16)
    
    cx <- residual_list[[2]][names(residual_list[[2]])%in%pred_term][i]
    a <- residual_list[[3]][i]
    y <- a + x*cx 
    
    xf <- x[order(x)]
    yf <- y[order(x)] ## so it draws a smooth lines not joining all dots
    
    lines(xf, yf, lwd = 2, col = cols[i])
    
  }
}

library(dplyr)

metric_names <- c("Omega", "theta", "distance_to_exclusion","differential", "differential_old", "Omega_all_pairs")

# Combine one saved dataset with its corresponding overlap results.
combine_one_dataset <- function(saved_results, overlap_results,
                                season_name, treatment_name) {
  
  overlap_columns <- c(
    "Omega_overlap", "Omega_all_pairs", "differential", "overlap",
    "overlap_n_samples", "overlap_seed"
  )
  
  stopifnot(
    length(saved_results) == 99L,
    length(overlap_results) == length(saved_results)
  )
  
  combined_iterations <- vector("list", length(saved_results))
  
  for (j in seq_along(saved_results)) {
    
    original <- dplyr::bind_rows(saved_results[[j]]$result_list) %>%
      dplyr::mutate(plot = as.character(plot))
    
    updated <- dplyr::bind_rows(overlap_results[[j]]) %>%
      dplyr::mutate(plot = as.character(plot))
    
    stopifnot(
      !anyNA(original$plot),
      !anyNA(updated$plot),
      !anyDuplicated(original$plot),
      !anyDuplicated(updated$plot),
      setequal(original$plot, updated$plot),
      all(overlap_columns %in% names(updated))
    )
    
    # Match plots within the same randomization.
    matched_rows <- match(original$plot, updated$plot)
    
    # Confirm treatment and iteration where recorded.
    if ("treatment" %in% names(original) &&
        "treatment" %in% names(updated)) {
      stopifnot(
        as.character(original$treatment) ==
          as.character(updated$treatment[matched_rows])
      )
    }
    
    if ("randomization" %in% names(updated)) {
      stopifnot(all(updated$randomization == j))
    }
    
    # Preserve previous estimates before replacing them.
    if (!"differential_old" %in% names(original)) {
      original$differential_old <- original$differential
    }
    
    if (!"overlap_old" %in% names(original)) {
      original$overlap_old <- original$overlap
    }
    
    # Update only the overlap-related columns.
    for (column in overlap_columns) {
      original[[column]] <- updated[[column]][matched_rows]
    }
    
    original$randomization <- j
    original$season <- season_name
    original$Treatment <- treatment_name
    original$year <- 2016
    
    numeric_columns <- intersect(
      c(
        "Omega", "Omega1", "Omega2",
        "theta", "theta1", "theta2",
        "distance_to_exclusion", "feasibility", "coex_rate",
        "differential_old", "overlap_old",
        overlap_columns
      ),
      names(original)
    )
    
    original <- original %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(numeric_columns),
          ~ as.numeric(as.character(.x))
        )
      )
    
    combined_iterations[[j]] <- original
  }
  
  combined <- dplyr::bind_rows(combined_iterations)
  
  stopifnot(
    !anyDuplicated(combined[c("plot", "randomization")])
  )
  combined
}
