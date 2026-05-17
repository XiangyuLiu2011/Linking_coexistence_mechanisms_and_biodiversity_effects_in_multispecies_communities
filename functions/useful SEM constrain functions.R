##function that extract paths for the multigroup model
library(lavaan)
library(pbapply) # For progress bar

test_stepwise_constraints_cfa <- function(fit, data) {
  fit_configural <- fit
  group_var <- lavInspect(fit, "group")
  
  path.list <- extract_paths_cfa(fit_configural)
  
  # Remove covariance paths (~~)
  path.list <- lapply(path.list, function(x) x[!grepl("~~", x)])
  
  fit_list <- list()
  model_compare <- data.frame(Path = unlist(path.list), 
                              ChiSquare = NA, pValue = NA, stringsAsFactors = FALSE)
  
  generate_coef_name <- function(index) paste0("coef_", index)
  
  ngroups <- length(levels(data[, group_var]))
  groups <- seq(1, ngroups, 1)
  
  total_steps <- sum(sapply(path.list, length))
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0
  
  for (i in seq_along(path.list)) {
    for (j in seq_along(path.list[[i]])) {
      step <- step + 1
      setTxtProgressBar(pb, step)
      
      path <- path.list[[i]][j]
      coef_name <- generate_coef_name(step)
      
      constrained_model_string <- constrain_path_cfa(fit_configural, path, coef_name, groups_constrain = groups)
      
      fit_constrained <- sem(constrained_model_string, data = data, group = group_var, fixed.x = FALSE)
      fit_list[[step]] <- fit_constrained
      
      compare <- anova(fit_configural, fit_constrained)
      chi_square <- compare[2, "Chisq diff"]
      p_value <- compare[2, "Pr(>Chisq)"]
      
      model_compare[step, "ChiSquare"] <- chi_square
      model_compare[step, "pValue"] <- p_value
    }
  }
  close(pb)
  
  model_compare <- model_compare[order(-model_compare$pValue),]
  model_compare2 <- model_compare
  new_fit_list <- list()
  latest_fit_configural <- fit_configural
  
  pb2 <- txtProgressBar(min = 0, max = nrow(model_compare2), style = 3)
  
  while (nrow(model_compare2) > 0) {
    setTxtProgressBar(pb2, nrow(model_compare2))
    
    path <- model_compare2$Path[1]
    coef_name <- generate_coef_name(nrow(model_compare2))
    
    next_constrain_ms <- constrain_path_cfa(latest_fit_configural, path, coef_name, groups_constrain = groups)
    
    fit_constrained <- sem(next_constrain_ms, data = data, group = group_var, fixed.x = FALSE)
    compare <- anova(latest_fit_configural, fit_constrained)
    
    chi_square <- compare[2, "Chisq diff"]
    p_value <- compare[2, "Pr(>Chisq)"]
    
    # Check if p_value is NA or if degrees of freedom are the same
    if (is.na(p_value)) {
      warning("p_value is NA for path: ", path)
      model_compare2 <- model_compare2[-1, , drop = FALSE]
      next
    }
    
    if (compare[2, "Df diff"] == 0) {
      warning("Degrees of freedom are the same for path: ", path)
      model_compare2 <- model_compare2[-1, , drop = FALSE]
      next
    }
    
    if (p_value > 0.05) { 
      new_fit_list[[paste0("Step_", nrow(model_compare2))]] <- fit_constrained
      latest_fit_configural <- fit_constrained
      final_model_string <- next_constrain_ms
      model_compare2 <- model_compare2[-1, , drop = FALSE]
    } else {
      model_compare2 <- model_compare2[-1, , drop = FALSE]
    }
    
    print(coef_name)
  }
  close(pb2)
  
  return(list(
    fit_list = new_fit_list,
    final_fit = latest_fit_configural,
    single_path_constrain = model_compare,
    final_model_string = final_model_string
  ))
}

constrain_path_cfa <- function(fit, path, coef_name, groups_constrain, intercepts = FALSE) {
  
  model_paths <- extract_paths_cfa(fit, intercepts = intercepts)
  model_pathsc <- model_paths
  
  con_ms <- list()
  
  for (i in groups_constrain) {
    
    mp <- model_pathsc[[i]]
    
    con_paths <- gsub(
      pattern = paste0("^", path, "$"),
      replacement = paste0(
        sub("~\\s*", "~ ", path),
        coef_name,
        " * ",
        sub(".*~\\s*", "", path)
      ),
      x = mp
    )
    
    con_ms[[i]] <- paste(
      "Group: ", i, "\n",
      paste(con_paths, collapse = "\n")
    )
  }
  
  total_groups <- lavInspect(fit, "ngroups")
  unconstrained_groups <- setdiff(seq_len(total_groups), groups_constrain)
  
  uncon_ms <- lapply(unconstrained_groups, function(g) {
    paste(
      "Group: ", g, "\n",
      paste(model_paths[[g]], collapse = "\n")
    )
  })
  
  fin_ms <- vector("list", total_groups)
  fin_ms[groups_constrain] <- con_ms[groups_constrain]
  fin_ms[unconstrained_groups] <- uncon_ms
  
  final_model_string <- paste(fin_ms, collapse = "\n")
  
  return(final_model_string)
}

sort.pars <- function(model, response="All", op.type = "regression", sig=FALSE, group.names, level = "All"){
  
  pp <- parameterEstimates(model)
  pp$group <- factor(pp$group)
  levels(pp$group) <- group.names
  
  op.type[op.type=="regression"] <- "~"
  op.type[op.type=="intercepts"] <- "~1"
  op.type[op.type=="covariances"] <- "~~"
  
  if(any(op.type=="All")){
    pp.reg <- pp
  }
  else{
    pp.reg <- pp[pp$op%in%op.type,]
  }
  
  oo <- order(pp.reg$lhs, pp.reg$rhs)
  pp.reg2 <- pp.reg[oo,]
  
  if(sig==TRUE)
  {
    pp.reg2 <- pp.reg2[which(pp.reg2$pvalue<0.1),]  
  }
  else{}
  if(response =="All"){
    pp.return <- pp.reg2
  }
  else{
    pp.return <- pp.reg2[pp.reg2$lhs%in%response,]
  }
  
  if(level == "All"){
    pp.return <- pp.return
  }
  else{
    pp.return <- pp.return[pp.return$level == level, ]
  }
  
  return(pp.return)
}

## ============================================================
## Extract lavaan model paths by group
## ============================================================

extract_paths_cfa <- function(fit, intercepts = FALSE) {
  
  pp <- parTable(fit)
  
  ## Select directed paths and covariances
  pp2 <- pp[pp$op %in% c("~", "~~"), ]
  
  ## Remove non-free parameters, such as automatically estimated
  ## exogenous covariances and variances
  pp3 <- pp2[pp2$user == 1, ]
  
  lhs <- pp3$lhs
  rhs <- pp3$rhs
  ops <- pp3$op
  groups <- pp3$group
  
  ## Keep existing parameter labels if some paths are already constrained
  if ("label" %in% names(pp3)) {
    label <- pp3$label
    cons <- rhs[label != ""]
    conlab <- label[label != ""]
    cons2 <- paste(conlab, cons, sep = "*")
    rhs[label != ""] <- cons2
  }
  
  ## Reconstruct model syntax
  all.paths <- paste(lhs, ops, rhs, sep = " ")
  
  ## Split model syntax by group
  path.list <- split(all.paths, groups)
  
  ## Optionally include intercepts
  if (intercepts) {
    
    ii <- pp[pp$op == "~1", ]
    
    lhsi <- ii$lhs
    grii <- ii$group
    rhsi <- rep("~1", length(lhsi))
    
    ## Keep existing parameter labels for intercepts
    if ("label" %in% names(ii)) {
      label <- ii$label
      conlab <- label[label != ""]
      rhsi[label != ""] <- paste("~", conlab, "*1", sep = "")
    }
    
    ints <- paste(lhsi, rhsi, sep = " ")
    
    ## Add intercepts to each group
    for (j in unique(groups)) {
      path.list[[as.character(j)]] <- c(
        path.list[[as.character(j)]],
        ints[grii == j]
      )
    }
  }
  
  return(path.list)
}