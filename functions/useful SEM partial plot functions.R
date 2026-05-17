## ============================================================
## Calculate partial residuals from lavaan SEM
## ============================================================

part.res.sem <- function(
    fit,
    pred_term,
    response,
    data,
    groups,
    level,
    correct_group = FALSE
) {
  
  group_var <- lavInspect(fit, "group")
  
  ppfull <- parameterEstimates(fit)
  
  ## Select parameters from the requested level and groups
  pp <- ppfull[ppfull$level == level, ]
  pp <- pp[pp$group %in% groups, ]
  
  ## Select all effects on the response variable
  pp2 <- pp[pp$lhs == response, ]
  
  ## Remove variance and covariance terms
  pp3 <- pp2[pp2$op != "~~", ]
  
  estimates <- pp3[, "est"]
  term_names <- pp3[, "rhs"]
  ggl <- pp3[, "group"]
  
  ## Separate slopes and intercepts
  slope_names <- term_names[-which(pp3$op == "~1")]
  intercept_names <- term_names[which(pp3$op == "~1")]
  
  ggl <- ggl[-which(pp3$op == "~1")]
  
  slopes <- estimates[-which(pp3$op == "~1")]
  intercepts <- estimates[which(pp3$op == "~1")]
  
  ## Remove the predictor of interest
  slope_names2 <- slope_names[-which(slope_names %in% pred_term)]
  slopes2 <- slopes[-which(slope_names %in% pred_term)]
  ggl2 <- ggl[-which(slope_names %in% pred_term)]
  
  ss_mat <- data.frame(
    slope_names2 = slope_names2,
    slopes2 = slopes2,
    ggl2 = ggl2
  )
  
  part_res <- list()
  
  ## Data should be organised as a list, with one element per group
  dd <- data
  
  for (i in seq_along(groups)) {
    
    ## Select all predictors except the predictor of interest
    dd2 <- dd[[i]][, match(unique(slope_names2), colnames(dd[[i]]))]
    
    ## For level-2 variables, remove intercepts from other predictors
    if (level == 2) {
      
      ii <- pp[pp$op == "~1", "est"]
      ii_names <- pp[pp$op == "~1", "lhs"]
      ggs <- pp[pp$op == "~1" & pp$level == 2, "group"]
      
      ii2_names <- ii_names[-which(ii_names %in% pred_term)]
      ii2 <- ii[-which(ii_names %in% pred_term)]
      
      ii2 <- ii2[-which(ii2_names %in% response)]
      ggs2 <- ggs[-which(ii_names %in% pred_term)]
      ggs2 <- ggs2[-which(ii2_names %in% response)]
      ii2_names <- ii2_names[-which(ii2_names %in% response)]
      
      ii_mat <- data.frame(
        ggs2 = ggs2,
        ii2 = ii2,
        ii2_names = ii2_names
      )
      
      ## Remove intercepts from selected variables
      dd2a <- dd2[, match(unique(ii2_names), colnames(dd2))]
      dd2b <- t(t(dd2a) - ii_mat[ggs2 == i, "ii2"])
      dd2 <- dd2b
    }
    
    ## Effect of all other predictors: slope x value
    ddeff <- t(t(dd2) * ss_mat[ggl2 == i, "slopes2"])
    
    ## Total effect of all other predictors
    ddsum <- apply(ddeff, 1, sum)
    
    ## Extract response
    resp <- dd[[i]][, match(response, colnames(dd[[i]]))]
    
    ## Partial residuals
    if (correct_group == TRUE) {
      part_res[[i]] <- as.vector(resp) - ddsum - intercepts[i]
    } else {
      part_res[[i]] <- as.vector(resp) - ddsum
    }
  }
  
  ## Return the slope of the focal predictor
  return_slopes <- slopes[which(slope_names %in% pred_term)]
  names(return_slopes) <- slope_names[which(slope_names %in% pred_term)]
  
  result <- list(
    partial_residuals = part_res,
    predictor_slope = return_slopes,
    intercept = intercepts
  )
  
  return(result)
}

plot.res.sem <- function(fit_data, pred_term, residual_list, cols,
                         alpha = 0.1, xlab = "", ylab = "",
                         xlim = NULL, ylim = NULL) {
  
  # number of groups (works for 2,3,4,...)
  ngrp <- length(fit_data)
  
  ## create sensible defaults for xlim if missing
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
    xr <- diff(xlim)
    if (xr == 0) xr <- ifelse(abs(xlim[1]) > 0, abs(xlim[1]) * 0.1, 0.1)
    xlim <- xlim + c(-0.05, 0.05) * xr
  }
  
  ## create sensible defaults for ylim if missing
  if (is.null(ylim)) {
    # residual values
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    # predicted lines (use same formula you already used)
    ys_pred <- unlist(lapply(seq_len(ngrp), function(i) {
      x <- fit_data[[i]][, pred_term]
      cx <- residual_list[[2]][ names(residual_list[[2]]) %in% pred_term ][i]
      a  <- residual_list[[3]][i]
      a + x * cx
    }))
    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
    yr <- diff(ylim)
    if (yr == 0) yr <- ifelse(abs(ylim[1]) > 0, abs(ylim[1]) * 0.1, 0.1)
    ylim <- ylim + c(-0.05, 0.05) * yr
  }
  
  ## start plot
  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
  
  cols2 <- alpha(cols, alpha = alpha) ## for the points (requires scales::alpha)
  
  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    
    ## points: partial residuals
    points(x, residual_list[[1]][[i]], col = cols2[i], pch = 16)
    
    ## slope and intercept from residual_list (keeps your original logic)
    cx <- residual_list[[2]][ names(residual_list[[2]]) %in% pred_term ][i]
    a  <- residual_list[[3]][i]
    y  <- a + x * cx
    
    xf <- x[order(x)]
    yf <- y[order(x)]
    lines(xf, yf, lwd = 2, col = cols[i])
  }
}


plot.res.sem <- function(fit_data, pred_term, residual_list, cols,
                         alpha = 0.1, xlab = "", ylab = "",
                         xlim = NULL, ylim = NULL,
                         group_labels = NULL) {
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # x limits
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
  }
  
  # y limits - include both residuals and predicted regression lines
  if (is.null(ylim)) {
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    
    # Pre-calculate regression lines to determine y-range
    ys_pred <- unlist(lapply(seq_len(ngrp), function(i) {
      x <- fit_data[[i]][, pred_term]
      y_res <- residual_list[[1]][[i]]
      
      # Calculate linear regression
      lm_fit <- lm(y_res ~ x)
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      
      return(a + x * cx)
    }))
    
    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
  }
  
  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
  
  cols2 <- alpha(cols, alpha = alpha)
  
  # Store regression coefficients for potential later use
  regression_coefs <- list()
  
  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    y_res <- residual_list[[1]][[i]]
    
    # Plot partial residuals
    points(x, y_res, col = cols2[i], pch = 16)
    
    # Calculate linear regression
    lm_fit <- lm(y_res ~ x)
    a <- coef(lm_fit)[1]  # intercept
    cx <- coef(lm_fit)[2] # slope
    
    # Store coefficients
    regression_coefs[[i]] <- list(intercept = a, slope = cx, 
                                  r_squared = summary(lm_fit)$r.squared)
    
    # Calculate predicted values
    y_pred <- a + x * cx
    
    # Plot regression line
    xf <- x[order(x)]
    yf <- y_pred[order(x)]
    lines(xf, yf, lwd = 2, col = cols[i])
    
    # Add group label at the end of the regression line
    text(max(xlim), tail(yf, 1), group_labels[i], col = cols[i], pos = 4, cex = 0.8)
    
    # Optional: Add R² value to the plot
    r2_text <- sprintf("R² = %.2f", summary(lm_fit)$r.squared)
    text(mean(xlim), min(ylim) + (i * diff(ylim)/15), r2_text, 
         col = cols[i], cex = 0.7, pos = 4)
  }
  
  # Return regression coefficients invisibly for potential further analysis
  names(regression_coefs) <- group_labels
  invisible(regression_coefs)
}

plot.res.sem2 <- function(fit_data, pred_term, residual_list, cols,
                         alpha = 0.1, xlab = "", ylab = "",
                         xlim = NULL, ylim = NULL,
                         group_labels = NULL,
                         group_combinations = NULL) {
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # If group_combinations is provided, use it to combine groups
  if (!is.null(group_combinations)) {
    # Combine data for the specified groups
    combined_fit_data <- list()
    combined_residuals <- list()
    combined_labels <- c()
    combined_cols <- c()
    
    for (comb_idx in seq_along(group_combinations)) {
      groups_to_combine <- group_combinations[[comb_idx]]
      
      # Combine predictor data
      combined_pred <- unlist(lapply(groups_to_combine, function(i) fit_data[[i]][, pred_term]))
      
      # Combine residuals
      combined_res <- unlist(lapply(groups_to_combine, function(i) residual_list[[1]][[i]]))
      
      combined_fit_data[[comb_idx]] <- data.frame(x = combined_pred)
      colnames(combined_fit_data[[comb_idx]]) <- pred_term
      
      combined_residuals[[comb_idx]] <- combined_res
      combined_labels[comb_idx] <- paste(group_labels[groups_to_combine], collapse = " + ")
      combined_cols[comb_idx] <- cols[groups_to_combine[1]] # Use first group's color
    }
    
    # Replace original data with combined data
    fit_data <- combined_fit_data
    residual_list[[1]] <- combined_residuals
    group_labels <- combined_labels
    cols <- combined_cols
    ngrp <- length(fit_data)
  }
  
  # x limits
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
  }
  
  # y limits - include both residuals and predicted regression lines
  if (is.null(ylim)) {
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    
    # Pre-calculate regression lines to determine y-range
    ys_pred <- unlist(lapply(seq_len(ngrp), function(i) {
      x <- fit_data[[i]][, pred_term]
      y_res <- residual_list[[1]][[i]]
      
      # Calculate linear regression
      lm_fit <- lm(y_res ~ x)
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      
      return(a + x * cx)
    }))
    
    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
  }
  
  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
  
  cols2 <- alpha(cols, alpha = alpha)
  
  # Store regression coefficients for potential later use
  regression_coefs <- list()
  
  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    y_res <- residual_list[[1]][[i]]
    
    # Plot partial residuals
    points(x, y_res, col = cols2[i], pch = 16)
    
    # Calculate linear regression
    lm_fit <- lm(y_res ~ x)
    a <- coef(lm_fit)[1]  # intercept
    cx <- coef(lm_fit)[2] # slope
    
    # Store coefficients
    regression_coefs[[i]] <- list(intercept = a, slope = cx, 
                                  r_squared = summary(lm_fit)$r.squared)
    
    # Calculate predicted values
    y_pred <- a + x * cx
    
    # Plot regression line
    xf <- x[order(x)]
    yf <- y_pred[order(x)]
    lines(xf, yf, lwd = 2, col = cols[i])
    
    # Add group label at the end of the regression line
    text(max(xlim), tail(yf, 1), group_labels[i], col = cols[i], pos = 4, cex = 0.8)
    
    # Optional: Add R² value to the plot
    r2_text <- sprintf("R² = %.2f", summary(lm_fit)$r.squared)
    #text(mean(xlim), min(ylim) + (i * diff(ylim)/15), r2_text, 
        # col = cols[i], cex = 0.7, pos = 4)
  }
  
  # Return regression coefficients invisibly for potential further analysis
  names(regression_coefs) <- group_labels
  invisible(regression_coefs)
}


plot.res.sem2 <- function(fit_data, pred_term, residual_list, cols,
                          alpha = 0.1, xlab = "", ylab = "",
                          xlim = NULL, ylim = NULL,
                          group_labels = NULL,
                          group_combinations = NULL) {
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # If group_combinations is provided, use it to combine groups
  if (!is.null(group_combinations)) {
    # Combine data for the specified groups
    combined_fit_data <- list()
    combined_residuals <- list()
    combined_labels <- c()
    combined_cols <- c()
    
    for (comb_idx in seq_along(group_combinations)) {
      groups_to_combine <- group_combinations[[comb_idx]]
      
      # Combine predictor data
      combined_pred <- unlist(lapply(groups_to_combine, function(i) fit_data[[i]][, pred_term]))
      
      # Combine residuals
      combined_res <- unlist(lapply(groups_to_combine, function(i) residual_list[[1]][[i]]))
      
      combined_fit_data[[comb_idx]] <- data.frame(x = combined_pred)
      colnames(combined_fit_data[[comb_idx]]) <- pred_term
      
      combined_residuals[[comb_idx]] <- combined_res
      combined_labels[comb_idx] <- paste(group_labels[groups_to_combine], collapse = " + ")
      combined_cols[comb_idx] <- cols[groups_to_combine[1]] # Use first group's color
    }
    
    # Replace original data with combined data
    fit_data <- combined_fit_data
    residual_list[[1]] <- combined_residuals
    group_labels <- combined_labels
    cols <- combined_cols
    ngrp <- length(fit_data)
  }
  
  # x limits
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
  }
  
  # y limits - include both residuals and predicted regression lines
  if (is.null(ylim)) {
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    
    # Pre-calculate regression lines to determine y-range (only for groups with R² > 0)
    ys_pred <- unlist(lapply(seq_len(ngrp), function(i) {
      x <- fit_data[[i]][, pred_term]
      y_res <- residual_list[[1]][[i]]
      
      # Calculate linear regression
      lm_fit <- lm(y_res ~ x)
      r_squared <- summary(lm_fit)$r.squared
      
      # Only include in y-range calculation if R² > 0
      if (r_squared > 0) {
        a <- coef(lm_fit)[1]  # intercept
        cx <- coef(lm_fit)[2] # slope
        return(a + x * cx)
      } else {
        return(numeric(0)) # Return empty vector
      }
    }))
    
    # Filter out empty vectors
    ys_pred <- ys_pred[lengths(ys_pred) > 0]
    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
  }
  
  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
  
  cols2 <- alpha(cols, alpha = alpha)
  
  # Store regression coefficients for potential later use
  regression_coefs <- list()
  
  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    y_res <- residual_list[[1]][[i]]
    
    # Plot partial residuals
    points(x, y_res, col = cols2[i], pch = 16)
    
    # Calculate linear regression
    lm_fit <- lm(y_res ~ x)
    r_squared <- summary(lm_fit)$r.squared
    
    # Store coefficients
    regression_coefs[[i]] <- list(intercept = coef(lm_fit)[1], 
                                  slope = coef(lm_fit)[2], 
                                  r_squared = r_squared)
    
    # Only plot regression line if R² > 0
    if (r_squared > 0) {
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      
      # Calculate predicted values
      y_pred <- a + x * cx
      
      # Plot regression line
      xf <- x[order(x)]
      yf <- y_pred[order(x)]
      lines(xf, yf, lwd = 2, col = cols[i])
      
      # Add group label at the end of the regression line
      text(max(xlim), tail(yf, 1), group_labels[i], col = cols[i], pos = 4, cex = 0.8)
      
      # Add R² value to the plot
      r2_text <- sprintf("R² = %.2f", r_squared)
      text(mean(xlim), min(ylim) + (i * diff(ylim)/15), r2_text, 
           col = cols[i], cex = 0.7, pos = 4)
    } else {
      # Add group label without regression line (positioned differently)
      text(mean(xlim), min(ylim) + (i * diff(ylim)/15), 
           paste0(group_labels[i], " (R² = 0)"), 
           col = cols[i], cex = 0.7, pos = 4)
    }
  }
  
  # Return regression coefficients invisibly for potential further analysis
  names(regression_coefs) <- group_labels
  invisible(regression_coefs)
}

plot.res.sem3 <- function(fit_data, pred_term, residual_list, cols,
                         alpha = 0.1, xlab = "", ylab = "",
                         xlim = NULL, ylim = NULL,
                         group_labels = NULL,
                         group_combinations = NULL,
                         show_regression = NULL) {
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # If group_combinations is provided, use it to combine groups
  if (!is.null(group_combinations)) {
    # Combine data for the specified groups
    combined_fit_data <- list()
    combined_residuals <- list()
    combined_labels <- c()
    combined_cols <- c()
    
    for (comb_idx in seq_along(group_combinations)) {
      groups_to_combine <- group_combinations[[comb_idx]]
      
      # Combine predictor data
      combined_pred <- unlist(lapply(groups_to_combine, function(i) fit_data[[i]][, pred_term]))
      
      # Combine residuals
      combined_res <- unlist(lapply(groups_to_combine, function(i) residual_list[[1]][[i]]))
      
      combined_fit_data[[comb_idx]] <- data.frame(x = combined_pred)
      colnames(combined_fit_data[[comb_idx]]) <- pred_term
      
      combined_residuals[[comb_idx]] <- combined_res
      combined_labels[comb_idx] <- paste(group_labels[groups_to_combine], collapse = " + ")
      combined_cols[comb_idx] <- cols[groups_to_combine[1]] # Use first group's color
    }
    
    # Replace original data with combined data
    fit_data <- combined_fit_data
    residual_list[[1]] <- combined_residuals
    group_labels <- combined_labels
    cols <- combined_cols
    ngrp <- length(fit_data)
  }
  
  # Set default: show all regression lines if not specified
  if (is.null(show_regression)) {
    show_regression <- rep(TRUE, ngrp)
  } else if (length(show_regression) != ngrp) {
    stop("show_regression must have the same length as the number of groups")
  }
  
  # x limits
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
  }
  
  # Pre-calculate all regressions first
  regression_results <- list()
  valid_regression_indices <- c() # Track which groups have valid regressions to show
  
  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    y_res <- residual_list[[1]][[i]]
    
    # Calculate linear regression
    lm_fit <- lm(y_res ~ x)
    r_squared <- summary(lm_fit)$r.squared
    
    regression_results[[i]] <- list(
      x = x,
      y_res = y_res,
      lm_fit = lm_fit,
      r_squared = r_squared,
      has_valid_regression = (r_squared > 0 & !is.na(r_squared)),
      show_line = show_regression[i] # Whether to show this regression line
    )
    
    if (regression_results[[i]]$has_valid_regression && regression_results[[i]]$show_line) {
      valid_regression_indices <- c(valid_regression_indices, i)
    }
  }
  
  # y limits - include both residuals and predicted regression lines (only for valid regressions to show)
  if (is.null(ylim)) {
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) regression_results[[i]]$y_res))
    
    # Only include valid regression lines that should be shown
    ys_pred <- unlist(lapply(valid_regression_indices, function(i) {
      x <- regression_results[[i]]$x
      lm_fit <- regression_results[[i]]$lm_fit
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      return(a + x * cx)
    }))
    
    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
  }
  
  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
  
  cols2 <- alpha(cols, alpha = alpha)
  
  # Store regression coefficients for potential later use
  regression_coefs <- list()
  
  for (i in seq_len(ngrp)) {
    x <- regression_results[[i]]$x
    y_res <- regression_results[[i]]$y_res
    lm_fit <- regression_results[[i]]$lm_fit
    r_squared <- regression_results[[i]]$r_squared
    has_valid_regression <- regression_results[[i]]$has_valid_regression
    show_line <- regression_results[[i]]$show_line
    
    # Plot partial residuals
    points(x, y_res, col = cols2[i], pch = 16)
    
    # Store coefficients
    regression_coefs[[i]] <- list(
      intercept = coef(lm_fit)[1], 
      slope = coef(lm_fit)[2], 
      r_squared = r_squared,
      has_valid_regression = has_valid_regression,
      show_line = show_line
    )
    
    # Only plot regression line if R² > 0, not NA, AND we want to show it
    if (has_valid_regression && show_line) {
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      
      # Calculate predicted values
      y_pred <- a + x * cx
      
      # Plot regression line
      xf <- x[order(x)]
      yf <- y_pred[order(x)]
      lines(xf, yf, lwd = 2, col = cols[i])
      
      # Add group label at the end of the regression line
      text(max(xlim), tail(yf, 1), group_labels[i], col = cols[i], pos = 4, cex = 0.8)
      
      # Add R² value to the plot
      r2_text <- sprintf("R² = %.2f", r_squared)
      #text(mean(xlim), min(ylim) + (i * diff(ylim)/15), r2_text, 
           #col = cols[i], cex = 0.7, pos = 4)
    } else {
      # Add group label without regression line
      #text(mean(xlim), min(ylim) + (i * diff(ylim)/15), 
           #paste0(group_labels[i], ifelse(has_valid_regression, sprintf(" (R² = %.2f)", r_squared), " (R² = 0)")), 
           #col = cols[i], cex = 0.7, pos = 4)
    }
  }
  
  # Return regression coefficients invisibly for potential further analysis
  names(regression_coefs) <- group_labels
  invisible(regression_coefs)
}


plot.res.sem3_auto <- function(fit_data, pred_term, residual_list, cols,
                               alpha = 0.1, xlab = "", ylab = "",
                               xlim = NULL, ylim = NULL,
                               group_labels = NULL,
                               group_combinations = NULL,
                               show_regression = NULL) {
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # --- Detect if group_combinations means "all groups together" ---
  combine_all <- !is.null(group_combinations) && 
    length(group_combinations) == 1 &&
    length(group_combinations[[1]]) == ngrp
  
  if (combine_all) {
    message("Detected 'combine all groups' mode: fitting one regression across all points.")
    
    # Combine all data
    all_x <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    all_y <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    
    # Fit one regression across all
    lm_fit <- lm(all_y ~ all_x)
    r2 <- summary(lm_fit)$r.squared
    slope <- coef(lm_fit)[2]
    intercept <- coef(lm_fit)[1]
    
    # Define xlim and ylim
    if (is.null(xlim)) xlim <- range(all_x, na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(all_y, na.rm = TRUE)
    
    # Plot base
    plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
    
    # --- Determine nitrogen color group ---
    # Assuming groups 1–2 = no N, 3–4 = N
    noN_idx <- c(1, 2)
    N_idx <- c(3, 4)
    
    cols2 <- alpha(cols, alpha = alpha)
    
    # Plot points
    for (i in seq_len(ngrp)) {
      x <- fit_data[[i]][, pred_term]
      y_res <- residual_list[[1]][[i]]
      points(x, y_res, col = cols2[i], pch = 16)
    }
    
    # Plot one regression line
    x_seq <- seq(min(all_x, na.rm = TRUE), max(all_x, na.rm = TRUE), length.out = 100)
    y_pred <- intercept + slope * x_seq
    lines(x_seq, y_pred, col = "black", lwd = 2)
    #text(max(xlim), max(y_pred), sprintf("R² = %.2f", r2), pos = 2, cex = 0.8)
    
    # Add legend for nitrogen treatments
    legend("right",
           legend = c("Withou nitrogen", "With nitrogen"),
           col = c(unique(cols[noN_idx])[1], unique(cols[N_idx])[1]),
           pch = 16,
           bty = "n")
    
    # Return model info
    return(invisible(list(
      slope = slope,
      intercept = intercept,
      R2 = r2
    )))
  }
  
  # Otherwise, default to your original behavior
  message("Using standard plot.res.sem3() behavior.")
  plot.res.sem3(
    fit_data, pred_term, residual_list, cols,
    alpha = alpha, xlab = xlab, ylab = ylab,
    xlim = xlim, ylim = ylim,
    group_labels = group_labels,
    group_combinations = group_combinations,
    show_regression = show_regression
  )
}


overlay_sem_lines <- function(residual_list,           # e.g. pr_indirect_SEPS (same structure as plot.res.sem3 expects)
                              Bio_list_unscaled,      # list of unscaled data.frames (names match groups)
                              sem_model,              # lavaan fit object, e.g. fc1_sown_N_indirect_14.0
                              response = "SEPS",
                              predictor = "sta.differential",
                              group_combinations = NULL, # same list of combos used in plot.res.sem3 (or NULL)
                              group_names = NULL,     # optional; default = names(Bio_list_unscaled) or combined names
                              line_col = "red",
                              lwd = 2) {
  
  # Load lavaan if needed
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan required")
  pe <- lavaan::parameterEstimates(sem_model)
  
  # Prepare groups / combined groups exactly like plot.res.sem3 does
  if (is.null(group_names)) group_names <- names(Bio_list_unscaled)
  
  if (!is.null(group_combinations)) {
    combined_fit_data <- list()
    combined_residuals <- list()
    combined_labels <- c()
    for (comb_idx in seq_along(group_combinations)) {
      groups_to_combine <- group_combinations[[comb_idx]]
      # pooled x (unscaled) for combined group:
      pooled_x <- unlist(lapply(groups_to_combine, function(i) Bio_list_unscaled[[i]][[predictor]]))
      # pooled y (raw response, unscaled)
      pooled_y_raw <- unlist(lapply(groups_to_combine, function(i) Bio_list_unscaled[[i]][[response]]))
      combined_fit_data[[comb_idx]] <- data.frame(x = pooled_x, y_raw = pooled_y_raw)
      combined_residuals[[comb_idx]] <- unlist(lapply(groups_to_combine, function(i) residual_list[[1]][[i]]))
      combined_labels[comb_idx] <- paste(group_names[groups_to_combine], collapse = " + ")
    }
    # Replace with combined
    Bio_comb <- combined_fit_data
    res_comb <- combined_residuals
    groups_to_process <- seq_along(Bio_comb)
    final_group_names <- combined_labels
  } else {
    # no combination: one panel per list element
    Bio_comb <- lapply(seq_along(Bio_list_unscaled), function(i) {
      data.frame(x = Bio_list_unscaled[[i]][[predictor]],
                 y_raw = Bio_list_unscaled[[i]][[response]])
    })
    res_comb <- lapply(seq_along(Bio_list_unscaled), function(i) residual_list[[1]][[i]])
    groups_to_process <- seq_along(Bio_comb)
    final_group_names <- if (!is.null(group_names)) group_names else names(Bio_list_unscaled)
  }
  
  # For each displayed group, compute SEM raw slope and overlay line
  for (i in groups_to_process) {
    # Determine which lavaan group index to read coefficient from.
    # If group_combinations used we pick the first group's index inside that combination,
    # which matches plot.res.sem3 behavior (it used first group's color).
    if (!is.null(group_combinations)) {
      source_groups <- group_combinations[[i]]
      lavaan_group_index <- source_groups[1]
    } else {
      lavaan_group_index <- i
    }
    
    # Extract standardized slope from lavaan (op == "~")
    beta_std_idx <- which(pe$lhs == response & pe$rhs == predictor & pe$op == "~" & pe$group == lavaan_group_index)
    if (length(beta_std_idx) == 0) {
      warning(sprintf("No SEM coefficient found for group %s (group index %s). Skipping line.", 
                      final_group_names[i], lavaan_group_index))
      next
    }
    beta_std <- pe$est[beta_std_idx]
    
    # Raw x and raw response (unscaled) for this displayed group
    df_pair <- Bio_comb[[i]]
    x_un <- df_pair$x
    y_raw <- df_pair$y_raw
    
    # Partial residuals for this displayed group
    y_res <- res_comb[[i]]
    # Safety: ensure lengths align - if not, try to match by length or warn
    if (length(y_res) != length(x_un)) {
      # If lengths differ, we still compute sd on pooled raw data for converting slope,
      # but intercept calculation is computed from y_res mean and x_un mean (works if order differs).
      warning(sprintf("Length mismatch for group '%s': length(x)=%d, length(residuals)=%d. Continuing but check ordering.",
                      final_group_names[i], length(x_un), length(y_res)))
    }
    
    # Compute sd_x and sd_y (use raw/unscaled y from Bio_list_unscaled)
    sd_x <- stats::sd(x_un, na.rm = TRUE)
    sd_y <- stats::sd(y_raw, na.rm = TRUE)
    if (is.na(sd_x) || sd_x == 0) {
      warning(sprintf("sd(x) is NA or 0 for group %s; cannot compute raw slope. Skipping.", final_group_names[i]))
      next
    }
    beta_raw <- beta_std * (sd_y / sd_x)
    
    # Compute intercept so line goes through mean of partial residuals
    intercept <- mean(y_res, na.rm = TRUE) - beta_raw * mean(x_un, na.rm = TRUE)
    
    # Build x range and add line to current plot
    xseq <- seq(min(x_un, na.rm = TRUE), max(x_un, na.rm = TRUE), length.out = 200)
    yline <- intercept + beta_raw * xseq
    
    # Draw the line on CURRENT device / plot
    lines(xseq, yline, col = line_col, lwd = lwd)
    
    # Optionally add text annotation of slope in the top-right corner of the group's plotting area
    # compute reasonable coords:
    x_pos <- median(xseq, na.rm = TRUE)
    y_pos <- median(yline, na.rm = TRUE)
    label <- sprintf("β_std=%.3f\nβ_raw=%.3f", beta_std, beta_raw)
    text(x = x_pos, y = y_pos, labels = label, adj = c(1.1, -0.2), cex = 0.7, col = line_col)
  }
  
  invisible(TRUE)
}

##updated on 06/02/2026
plot.res.sem3_auto1 <- function(fit_data, pred_term, residual_list, cols,
                               alpha = 0.1, xlab = "", ylab = "",
                               xlim = NULL, ylim = NULL,
                               group_labels = NULL,
                               group_combinations = NULL,
                               show_regression = NULL,
                               pch = 21,            ## CHANGED: new arg
                               cex = 1.2,           ## CHANGED: new arg
                               lwd = 0.6) {         ## CHANGED: new arg
  library(scales) # for alpha()
  
  ngrp <- length(fit_data)
  
  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }
  
  # --- Detect if group_combinations means "all groups together" ---
  combine_all <- !is.null(group_combinations) &&
    length(group_combinations) == 1 &&
    length(group_combinations[[1]]) == ngrp
  
  if (combine_all) {
    message("Detected 'combine all groups' mode: fitting one regression across all points.")
    
    # Combine all data
    all_x <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    all_y <- unlist(lapply(seq_len(ngrp), function(i) residual_list[[1]][[i]]))
    
    # Fit one regression across all
    lm_fit <- lm(all_y ~ all_x)
    r2 <- summary(lm_fit)$r.squared
    slope <- coef(lm_fit)[2]
    intercept <- coef(lm_fit)[1]
    
    # Define xlim and ylim
    if (is.null(xlim)) xlim <- range(all_x, na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(all_y, na.rm = TRUE)
    
    # Plot base
    plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)
    
    ## CHANGED: remove cols2 <- alpha(cols, alpha=alpha)
    ## We now use bg=alpha(color, alpha) and col=color for ggplot-like points
    
    # Plot points (ggplot-like)
    for (i in seq_len(ngrp)) {
      x <- fit_data[[i]][, pred_term]
      y_res <- residual_list[[1]][[i]]
      points(x, y_res,
             pch = pch, cex = cex, lwd = lwd,                 ## CHANGED
             bg  = scales::alpha(cols[i], alpha),             ## CHANGED
             col = cols[i])                                   ## CHANGED
    }
    
    # Plot one regression line
    x_seq <- seq(min(all_x, na.rm = TRUE), max(all_x, na.rm = TRUE), length.out = 100)
    y_pred <- intercept + slope * x_seq
    lines(x_seq, y_pred, col = "black", lwd = 2)
    
    return(invisible(list(
      slope = slope,
      intercept = intercept,
      R2 = r2
    )))
  }
  
  # Otherwise, default to your original behavior
  message("Using standard plot.res.sem3() behavior.")
  plot.res.sem3_1(
    fit_data, pred_term, residual_list, cols,
    alpha = alpha, xlab = xlab, ylab = ylab,
    xlim = xlim, ylim = ylim,
    group_labels = group_labels,
    group_combinations = group_combinations,
    show_regression = show_regression,
    pch = pch, cex = cex, lwd = lwd      ## CHANGED: pass through
  )
}

plot.res.sem3_1 <- function(fit_data, pred_term, residual_list, cols,
                          alpha = 0.1, xlab = "", ylab = "",
                          xlim = NULL, ylim = NULL,
                          group_labels = NULL,
                          group_combinations = NULL,
                          show_regression = NULL,
                          pch = 21,          ## CHANGED: new arg
                          cex = 1.2,         ## CHANGED: new arg
                          lwd = 0.6) {       ## CHANGED: new arg
  library(scales) # for alpha()

  ngrp <- length(fit_data)

  if (is.null(group_labels)) {
    group_labels <- names(fit_data)
    if (is.null(group_labels)) {
      group_labels <- paste0("Group", seq_len(ngrp))
    }
  }

  # If group_combinations is provided, use it to combine groups
  if (!is.null(group_combinations)) {
    # Combine data for the specified groups
    combined_fit_data <- list()
    combined_residuals <- list()
    combined_labels <- c()
    combined_cols <- c()

    for (comb_idx in seq_along(group_combinations)) {
      groups_to_combine <- group_combinations[[comb_idx]]

      # Combine predictor data
      combined_pred <- unlist(lapply(groups_to_combine, function(i) fit_data[[i]][, pred_term]))

      # Combine residuals
      combined_res <- unlist(lapply(groups_to_combine, function(i) residual_list[[1]][[i]]))

      combined_fit_data[[comb_idx]] <- data.frame(x = combined_pred)
      colnames(combined_fit_data[[comb_idx]]) <- pred_term

      combined_residuals[[comb_idx]] <- combined_res
      combined_labels[comb_idx] <- paste(group_labels[groups_to_combine], collapse = " + ")
      combined_cols[comb_idx] <- cols[groups_to_combine[1]] # Use first group's color
    }

    # Replace original data with combined data
    fit_data <- combined_fit_data
    residual_list[[1]] <- combined_residuals
    group_labels <- combined_labels
    cols <- combined_cols
    ngrp <- length(fit_data)
  }

  # Set default: show all regression lines if not specified
  if (is.null(show_regression)) {
    show_regression <- rep(TRUE, ngrp)
  } else if (length(show_regression) != ngrp) {
    stop("show_regression must have the same length as the number of groups")
  }

  # x limits
  if (is.null(xlim)) {
    xs <- unlist(lapply(seq_len(ngrp), function(i) fit_data[[i]][, pred_term]))
    xlim <- range(xs, na.rm = TRUE)
  }

  # Pre-calculate all regressions first
  regression_results <- list()
  valid_regression_indices <- c() # Track which groups have valid regressions to show

  for (i in seq_len(ngrp)) {
    x <- fit_data[[i]][, pred_term]
    y_res <- residual_list[[1]][[i]]

    # Calculate linear regression
    lm_fit <- lm(y_res ~ x)
    r_squared <- summary(lm_fit)$r.squared

    regression_results[[i]] <- list(
      x = x,
      y_res = y_res,
      lm_fit = lm_fit,
      r_squared = r_squared,
      has_valid_regression = (r_squared > 0 & !is.na(r_squared)),
      show_line = show_regression[i] # Whether to show this regression line
    )

    if (regression_results[[i]]$has_valid_regression && regression_results[[i]]$show_line) {
      valid_regression_indices <- c(valid_regression_indices, i)
    }
  }

  # y limits - include both residuals and predicted regression lines (only for valid regressions to show)
  if (is.null(ylim)) {
    ys_res <- unlist(lapply(seq_len(ngrp), function(i) regression_results[[i]]$y_res))

    # Only include valid regression lines that should be shown
    ys_pred <- unlist(lapply(valid_regression_indices, function(i) {
      x <- regression_results[[i]]$x
      lm_fit <- regression_results[[i]]$lm_fit
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope
      return(a + x * cx)
    }))

    ylim <- range(c(ys_res, ys_pred), na.rm = TRUE)
  }

  plot(1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = "white", las = 1)

  ## CHANGED: remove cols2 <- alpha(cols, alpha=alpha)
  ## We'll use bg = alpha(cols[i], alpha) for fill instead.

  # Store regression coefficients for potential later use
  regression_coefs <- list()

  for (i in seq_len(ngrp)) {
    x <- regression_results[[i]]$x
    y_res <- regression_results[[i]]$y_res
    lm_fit <- regression_results[[i]]$lm_fit
    r_squared <- regression_results[[i]]$r_squared
    has_valid_regression <- regression_results[[i]]$has_valid_regression
    show_line <- regression_results[[i]]$show_line

    # Plot partial residuals (ggplot-like)
    points(x, y_res,
           pch = pch, cex = cex, lwd = lwd,                  ## CHANGED
           bg  = scales::alpha(cols[i], alpha),              ## CHANGED
           col = cols[i])                                     ## CHANGED

    # Store coefficients
    regression_coefs[[i]] <- list(
      intercept = coef(lm_fit)[1],
      slope = coef(lm_fit)[2],
      r_squared = r_squared,
      has_valid_regression = has_valid_regression,
      show_line = show_line
    )

    # Only plot regression line if R² > 0, not NA, AND we want to show it
    if (has_valid_regression && show_line) {
      a <- coef(lm_fit)[1]  # intercept
      cx <- coef(lm_fit)[2] # slope

      # Calculate predicted values
      y_pred <- a + x * cx

      # Plot regression line
      xf <- x[order(x)]
      yf <- y_pred[order(x)]
      lines(xf, yf, lwd = 2, col = cols[i])

      # Add group label at the end of the regression line
      text(max(xlim), tail(yf, 1), group_labels[i], col = cols[i], pos = 4, cex = 0.8)
    }
  }

  # Return regression coefficients invisibly for potential further analysis
  names(regression_coefs) <- group_labels
  invisible(regression_coefs)
}


