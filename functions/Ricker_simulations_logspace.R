# Log-space Ricker simulations with irreversible extinction.
# Uses prepared alpha_sim and lambda_sim without further transformations.
# Species falling below the extinction threshold are permanently set to zero.
# Zero abundance is represented by -Inf in log space.
# Surviving species follow the same Ricker update as cxr RK.
# Extremely large growth may still exceed numerical limits.

library(cxr)
library(dplyr)
library(tidyr)

ricker_log_step <- function(log_abundance, C, lambda) {

  species <- names(lambda)

  stopifnot(
    !is.null(species),
    all(is.finite(log_abundance) | log_abundance == -Inf),
    all(is.finite(C)),
    all(is.finite(lambda)),
    all(lambda > 0),
    identical(names(log_abundance), species),
    identical(rownames(C), species),
    identical(colnames(C), species)
  )

  alive <- is.finite(log_abundance)

  # Extinct species retain log abundance -Inf, equivalent to abundance 0.
  z_next <- stats::setNames(rep(-Inf, length(species)), species)
  pressure <- stats::setNames(rep(0, length(species)), species)

  for (i in which(alive)) {

    active <- alive & C[i, ] != 0

    if (any(active)) {
      terms <- log(abs(C[i, active])) + log_abundance[active]
      pressure[i] <- sum(sign(C[i, active]) * exp(terms))
    }

    z_next[i] <- log_abundance[i] + log(lambda[i]) - pressure[i]
  }

  status <- if (any(!is.finite(pressure[alive]))) {
    "Nonfinite interaction pressure"
  } else if (any(!is.finite(z_next[alive]))) {
    "Nonfinite log abundance"
  } else if (any(z_next[alive] > log(.Machine$double.xmax))) {
    "Abundance exceeds numerical range"
  } else {
    "Valid"
  }

  list(
    log_abundance = z_next,
    interaction_pressure = pressure,
    status = status
  )
}

# Check agreement with cxr where ordinary-scale calculations are safe.
# Check log-space calculations against:
# 1. The Ricker equation.
# 2. The actual cxr Ricker projection function.
local({

  sp <- c("a", "b")

  C <- matrix(
    c(0.2, -0.03,
      0.04, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(sp, sp)
  )

  lam <- setNames(c(2, 1.5), sp)
  n <- setNames(c(1, 0.5), sp)

  # Direct mathematical calculation.
  reference_equation <- as.numeric(
    lam * n * exp(-drop(C %*% n))
  )

  # Actual package function, preserving coefficient names.
  reference_cxr <- vapply(sp, function(focal) {

    other_species <- setdiff(sp, focal)

    result <- cxr::RK_project_alpha_pairwise_lambdacov_none_alphacov_none(
      lambda = lam[focal],

      alpha_intra = setNames(
        as.numeric(C[focal, focal]),
        focal
      ),

      alpha_inter = setNames(
        as.numeric(C[focal, other_species]),
        other_species
      ),

      lambda_cov = NULL,
      alpha_cov = NULL,
      abundance = n,
      covariates = NULL
    )

    as.numeric(result)

  }, numeric(1))

  # Log-space calculation.
  actual <- ricker_log_step(log(n), C, lam)
  actual_abundance <- as.numeric(exp(actual$log_abundance))

  print(data.frame(
    species = sp,
    equation = reference_equation,
    cxr = as.numeric(reference_cxr),
    log_space = actual_abundance
  ))

  stopifnot(
    actual$status == "Valid",

    isTRUE(all.equal(
      actual_abundance,
      reference_equation,
      tolerance = 1e-12
    )),

    isTRUE(all.equal(
      actual_abundance,
      as.numeric(reference_cxr),
      tolerance = 1e-12
    ))
  )

  # Check that an extinct species stays absent and has no interaction effect.
z <- c(a = -Inf, b = log(0.5))

update <- ricker_log_step(
  log_abundance = z,
  C = C,
  lambda = lam
)

expected_b <- lam["b"] * 0.5 * exp(-C["b", "b"] * 0.5)

stopifnot(
  update$status == "Valid",
  update$log_abundance["a"] == -Inf,
  isTRUE(all.equal(
    unname(exp(update$log_abundance["b"])),
    unname(expected_b),
    tolerance = 1e-12
  ))
)

# Check that a completely extinct community stays at zero.
update_all_extinct <- ricker_log_step(
  log_abundance = c(a = -Inf, b = -Inf),
  C = C,
  lambda = lam
)

stopifnot(
  update_all_extinct$status == "Valid",
  all(update_all_extinct$log_abundance == -Inf)
)

message(
  "Checks passed: agreement with the Ricker equation and cxr; ",
  "extinct species remain at zero."
)
})

run_extinction_simulation_ricker_log <- function(
    alpha_sim,
    lambda_sim,
    comp_ini,
    richness_ini,
    dataset_names,
    remove_facilitation = TRUE,
    n_steps = 200L,
    initial_abundance = 1,
    extinction_threshold = 0.001,
    replacement_competition = 0.00001,
    richness_to_use = c(2, 3, 4, 6, 7, 8),
    start_rule = c("per_species", "equal_total"),
    threshold_rule = c("fixed", "fraction_of_start")
) {
  
  start_rule <- match.arg(start_rule)
  threshold_rule <- match.arg(threshold_rule)
  # initial_abundance is per species or per plot according to start_rule.
  # extinction_threshold is an absolute value or a fraction according to threshold_rule.
  stopifnot(
    length(alpha_sim) == length(dataset_names),
    length(lambda_sim) == length(dataset_names),
    length(comp_ini) == length(dataset_names),
    length(richness_ini) == length(dataset_names),
    n_steps >= 1,
    length(initial_abundance) == 1L, is.finite(initial_abundance), initial_abundance > 0,
    length(extinction_threshold) == 1L, is.finite(extinction_threshold),
    threshold_rule != "fraction_of_start" || extinction_threshold < 1,
    extinction_threshold > 0,
    replacement_competition >= 0
  )
  
  scenario <- if (remove_facilitation) {
    "Without facilitation"
  } else {
    "With facilitation"
  }
  
  summary_list <- list()
  trajectory_list <- list()
  species_list <- list()
  failure_list <- list()
  failure_examples <- list()
  result_id <- 0L
  
  for (i in seq_along(dataset_names)) {
    
    dataset <- dataset_names[i]
    season <- sub(" .*", "", dataset)
    treatment <- sub("^[^ ]+ ", "", dataset)
    
    # Align full matrix rows and columns.
    A <- as.matrix(alpha_sim[[i]])
    storage.mode(A) <- "double"
    
    stopifnot(
      !is.null(rownames(A)),
      !is.null(colnames(A)),
      !anyDuplicated(rownames(A)),
      !anyDuplicated(colnames(A)),
      setequal(rownames(A), colnames(A))
    )
    
    A <- A[rownames(A), rownames(A), drop = FALSE]
    
    L <- as.data.frame(lambda_sim[[i]])
    
    stopifnot(
      "intrinsic_growth" %in% names(L),
      !anyDuplicated(rownames(L))
    )
    
    composition_data <- as.data.frame(comp_ini[[i]])
    richness_data <- as.data.frame(richness_ini[[i]])
    
    stopifnot(
      "plot" %in% names(composition_data),
      "plot" %in% names(richness_data),
      ncol(composition_data) >= 2,
      ncol(richness_data) >= 2
    )
    
    richness_values <- as.numeric(
      as.character(richness_data[[2]])
    )
    
    eligible_plots <- unique(as.character(
      richness_data$plot[
        !is.na(richness_values) &
          richness_values %in% richness_to_use
      ]
    ))
    
    for (plot_id in eligible_plots) {
      
      richness_this <- unique(richness_values[
        as.character(richness_data$plot) == plot_id
      ])
      stopifnot(length(richness_this) == 1L)
      
      species <- unique(as.character(
        composition_data[
          as.character(composition_data$plot) == plot_id,
          2
        ]
      ))
      
      stopifnot(
        length(species) > 0,
        !anyNA(species),
        all(species %in% rownames(A)),
        all(species %in% rownames(L))
      )
      
      C <- A[species, species, drop = FALSE]
      
      lambda_plot <- stats::setNames(
        as.numeric(as.character(
          L[species, "intrinsic_growth"]
        )),
        species
      )
      
      stopifnot(
        all(is.finite(C)),
        all(is.finite(lambda_plot)),
        all(lambda_plot > 0)
      )
      
      n_facilitative <- sum(C < 0)
      
      if (remove_facilitation) {
        C[C < 0] <- replacement_competition
      }
      
      n_species <- length(species)
      start_per_species <- if (start_rule == "equal_total") initial_abundance / n_species else initial_abundance
      initial_plot_total <- if (start_rule == "equal_total") initial_abundance else n_species * initial_abundance
      plot_threshold <- if (threshold_rule == "fraction_of_start") extinction_threshold * start_per_species else extinction_threshold
      stopifnot(is.finite(initial_plot_total), plot_threshold > 0, start_per_species >= plot_threshold)
      first_above_total_step <- NA_integer_
      first_above_total_species <- NA_character_
      n_loss_half <- ceiling(n_species / 2)
      
      # Positive equilibrium check; not a stability test.
      equilibrium <- tryCatch(
        solve(C, log(lambda_plot)),
        error = function(e) rep(NA_real_, n_species)
      )
      
      feasible_Ricker <- if (all(is.finite(equilibrium))) {
        all(equilibrium > 0)
      } else {
        NA
      }
      
      # Row 1 = step 0.
      N <- matrix(
        NA_real_,
        nrow = n_steps + 1L,
        ncol = n_species,
        dimnames = list(NULL, species)
      )
      N[1, ] <- log(start_per_species)
      
      first_crossing <- stats::setNames(
        rep(NA_integer_, n_species),
        species
      )
      
      time_first_observed <- NA_integer_
      time_half_observed <- NA_integer_
      first_species <- NA_character_
      
      simulation_status <- "Completed"
      failure_step <- NA_integer_
      last_valid_step <- 0L
      
      for (step in seq_len(n_steps)) {
        
        log_previous <- N[step, ]
        names(log_previous) <- species
        update <- ricker_log_step(log_previous, C, lambda_plot)
        log_next <- update$log_abundance

        # A finite log prediction can be compared even if exp(log_next) would overflow.
        # This descriptive upper crossing has no effect on the population update.
        if (is.na(first_above_total_step) &&
            update$status %in% c("Valid", "Abundance exceeds numerical range")) {
          above_total <- log_next > log(initial_plot_total)
          if (any(above_total)) {
            first_above_total_step <- step
            first_above_total_species <- paste(species[above_total], collapse = "; ")
          }
        }
        if (update$status != "Valid") {
          simulation_status <- update$status
          failure_step <- step
          failure_list[[length(failure_list) + 1L]] <- tibble::tibble(
            start_rule = start_rule, threshold_rule = threshold_rule,
            scenario = scenario, dataset = dataset, plot = plot_id,
            season = season, treatment = treatment, failure_step = step,
            species = species,
            log_abundance_previous = as.numeric(log_previous),
            interaction_pressure = as.numeric(update$interaction_pressure),
            log_abundance_next = as.numeric(log_next),
            simulation_status = simulation_status
          )
          failure_examples[[paste(dataset, plot_id, sep = " | ")]] <- list(
            start_rule = start_rule, threshold_rule = threshold_rule,
            scenario = scenario, dataset = dataset, plot = plot_id,
            failure_step = step, alpha_matrix = C, lambda = lambda_plot,
            log_initial_abundances = log_previous, update = update
          )
          break
        }

        # Record first downward crossing for each species.
        crossed_now <- log_previous >= log(plot_threshold) &
          log_next < log(plot_threshold)
        
        newly_crossed <- crossed_now & is.na(first_crossing)
        first_crossing[newly_crossed] <- step
        
        if (
          is.na(time_first_observed) &&
          any(newly_crossed)
        ) {
          time_first_observed <- step
          first_species <- paste(
            species[newly_crossed],
            collapse = "; "
          )
        }
        
        # Cumulative distinct species that have crossed at least once.
        if (
          is.na(time_half_observed) &&
          sum(!is.na(first_crossing)) >= n_loss_half
        ) {
          time_half_observed <- step
        }
        
        # Apply irreversible extinction after recording threshold crossings.
        # All species that have crossed the threshold remain at abundance zero.
        log_next[!is.na(first_crossing)] <- -Inf

        N[step + 1L, ] <- log_next
        last_valid_step <- step
      }
      
      completed <- simulation_status == "Completed"
      
      event_first <- if (!is.na(time_first_observed)) {
        TRUE
      } else if (completed) {
        FALSE
      } else {
        NA
      }
      
      event_half <- if (!is.na(time_half_observed)) {
        TRUE
      } else if (completed) {
        FALSE
      } else {
        NA
      }
      
      time_first <- if (!is.na(time_first_observed)) {
        time_first_observed
      } else if (completed) {
        n_steps
      } else {
        NA_integer_
      }
      
      time_half <- if (!is.na(time_half_observed)) {
        time_half_observed
      } else if (completed) {
        n_steps
      } else {
        NA_integer_
      }
      
      last_log_abundance <- N[last_valid_step + 1L, ]
      last_abundance <- exp(last_log_abundance)
      below_last <- last_log_abundance < log(plot_threshold)
      
      result_id <- result_id + 1L
      
      summary_list[[result_id]] <- tibble::tibble(
        scenario = scenario,
        dataset = dataset,
        plot = plot_id,
        season = season,
        treatment = treatment,
        richness = richness_this,
        n_species_simulated = n_species,
        model_family = "RK",
        numerical_method = "log abundance",
        feasible_Ricker = feasible_Ricker,
        n_facilitative_coefficients = n_facilitative,
        initial_abundance = start_per_species,
        extinction_threshold = plot_threshold,
        start_rule = start_rule,
        threshold_rule = threshold_rule,
        initial_plot_total = initial_plot_total,
        first_above_total_step = first_above_total_step,
        first_above_total_species = first_above_total_species,
        max_steps = n_steps,
        n_loss_half = n_loss_half,
        
        # Times of irreversible extinction under the specified threshold rule.
        time_to_first_extinction = time_first_observed,
        time_to_half_extinction = time_half_observed,
        first_species = first_species,
        time_first = time_first,
        event_first = event_first,
        time_half = time_half,
        event_half = event_half,
        right_censored_first = !event_first,
        right_censored_half = !event_half,
        n_ever_crossed = sum(!is.na(first_crossing)),
        n_below_threshold_last = sum(below_last),
        proportion_remaining_last = mean(!below_last),
        proportion_remaining_end = if (completed) {
          mean(!below_last)
        } else {
          NA_real_
        },
        simulation_status = simulation_status,
        failure_step = failure_step,
        last_valid_step = last_valid_step
      )
      
      valid_N <- N[
        seq_len(last_valid_step + 1L),
        ,
        drop = FALSE
      ]
      
      trajectory_list[[result_id]] <-
        as.data.frame(valid_N) %>%
        dplyr::mutate(step = 0:last_valid_step) %>%
        tidyr::pivot_longer(
          cols = -step,
          names_to = "species",
          values_to = "log_abundance"
        ) %>%
        dplyr::mutate(
          scenario = scenario,
          dataset = dataset,
          plot = plot_id,
          season = season,
          treatment = treatment,
          richness = richness_this,
          abundance = exp(log_abundance),
          extinct = log_abundance == -Inf,
          abundance_rounded_to_zero = is.finite(log_abundance) & abundance == 0,
          below_threshold = log_abundance < log(plot_threshold),
          start_rule = start_rule, threshold_rule = threshold_rule,
          simulation_completed = completed
        )
      
      species_list[[result_id]] <- tibble::tibble(
        scenario = scenario,
        dataset = dataset,
        plot = plot_id,
        season = season,
        treatment = treatment,
        richness = richness_this,
        species = species,
        log_abundance_last = as.numeric(last_log_abundance),
        log_abundance_end = if (completed) as.numeric(last_log_abundance) else rep(NA_real_, n_species),
        abundance_end = if (completed) as.numeric(last_abundance) else rep(NA_real_, n_species),
        start_rule = start_rule, threshold_rule = threshold_rule,
        first_crossing_step = as.integer(first_crossing),
        ever_below_threshold = !is.na(first_crossing),
        abundance_last = as.numeric(last_abundance),
        below_threshold_last = as.logical(below_last),
        recovered_by_last_step =
          !is.na(first_crossing) & !below_last,
        simulation_completed = completed,
        last_valid_step = last_valid_step
      )
    }
  }
  
  list(
    summary = dplyr::bind_rows(summary_list),
    trajectories = dplyr::bind_rows(trajectory_list),
    species_extinctions = dplyr::bind_rows(species_list),
    failure_diagnostics = dplyr::bind_rows(failure_list),
    failure_examples = failure_examples,
    package_version = as.character(utils::packageVersion("cxr"))
  )
}


# Run after loading existing five simulation inputs.
# Produces named runs with summary, trajectories, species records and diagnostics.
run_ricker_start_comparison <- function(alpha_sim, lambda_sim, comp_ini, richness_ini,
                                        dataset_names, threshold_rule = "fraction_of_start") {
  threshold_rule <- match.arg(threshold_rule, c("fixed", "fraction_of_start"))
  runs <- list()
  for (start in c("per_species", "equal_total")) {
    for (remove in c(TRUE, FALSE)) {
      key <- paste(start, if (remove) "without" else "with", sep = "_")
      runs[[key]] <- run_extinction_simulation_ricker_log(
        alpha_sim, lambda_sim, comp_ini, richness_ini, dataset_names,
        remove_facilitation = remove, n_steps = 200L,
        initial_abundance = 1, extinction_threshold = 0.001,
        start_rule = start, threshold_rule = threshold_rule
      )
    }
  }
  runs
}
