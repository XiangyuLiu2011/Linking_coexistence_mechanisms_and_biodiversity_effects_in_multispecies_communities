##scripts/09_Time_to_extinction_analysis

# install.packages("cxr")
library(cxr)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(tidyr)
library(ggplot2)
library(lme4)
library(patchwork)
library(sjPlot)
library(performance)
library(tibble)
source("functions/Ricker_simulations_logspace.R")

# Use the inputs you have already prepared:
#
# alpha_sim: Positive coefficients = competition. Negative coefficients = facilitation.
# lambda_sim: No additional transformation is applied here.
# comp_ini, richness_ini:Community composition and richness, in the same dataset order.
dataset_names <- c("June C", "June N", "June F", "June NF", "August C", "August N", "August F", "August NF")

ricker_start_runs <- run_ricker_start_comparison(
  alpha_sim = alpha_sim,
  lambda_sim = lambda_sim,
  comp_ini = comp_ini,
  richness_ini = richness_ini,
  dataset_names = dataset_names,
  threshold_rule = "fraction_of_start"
)

ricker_start_runs_no_recovery <- run_ricker_start_comparison(
  alpha_sim = alpha_sim,
  lambda_sim = lambda_sim,
  comp_ini = comp_ini,
  richness_ini = richness_ini,
  dataset_names = dataset_names,
  threshold_rule = "fraction_of_start"
)

summary_no_recovery <- dplyr::bind_rows(lapply(ricker_start_runs_no_recovery, function(x) x$summary))
species_no_recovery <- dplyr::bind_rows(lapply(ricker_start_runs_no_recovery, function(x) x$species_extinctions))
summary_no_recovery %>% dplyr::count(start_rule, scenario, simulation_status, season) %>% print(n = Inf)

# 1. Combine plot-level results from the four runs.
ricker_start_summary <- bind_rows(
  ricker_start_runs$equal_total_without$summary, ricker_start_runs$equal_total_with$summary,
  ricker_start_runs$per_species_without$summary, ricker_start_runs$per_species_with$summary
)
stopifnot(all(ricker_start_summary$threshold_rule == "fraction_of_start"),
          !anyDuplicated(ricker_start_summary[c("start_rule", "scenario", "dataset", "plot")]))
print(count(ricker_start_summary, start_rule, scenario, simulation_status), n = Inf)

# 2. Choose the starting condition ----
# For the sensitivity analysis change to "per_species" and repeat sections 2–8.
rule <- "equal_total"
plot_summary <- filter(summary_no_recovery, start_rule == rule)
stopifnot(nrow(plot_summary) > 0L)
scenario_levels <- c("Without facilitation", "With facilitation")
endpoint_levels <- c("First species below threshold", "Half the species ever below threshold")

# 3. Summarize observed, censored and unresolved endpoints ----
# Check completion and feasibility separately for each scenario.
plot_summary %>% dplyr::count(scenario, feasible_Ricker, simulation_status)

event_summary_both <- plot_summary %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    first_observed = sum(event_first %in% TRUE),
    first_censored = sum(event_first %in% FALSE),
    first_unusable = sum(is.na(event_first)),
    half_observed = sum(event_half %in% TRUE),
    half_censored = sum(event_half %in% FALSE),
    half_unusable = sum(is.na(event_half)),
    .groups = "drop"
  )

crossing_data_both <- plot_summary %>%
  dplyr::select(
    scenario, dataset, plot, season, treatment, richness,
    max_steps, simulation_status,
    time_first, event_first, time_half, event_half
  ) %>%
  tidyr::pivot_longer(
    cols = c(time_first, event_first, time_half, event_half),
    names_to = c(".value", "endpoint"),
    names_pattern = "(time|event)_(first|half)"
  ) %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = scenario_levels),
    Endpoint = factor(
      endpoint,
      levels = c("first", "half"),
      labels = endpoint_levels
    ),
    Season = factor(season, levels = c("June", "August")),
    plot = as.character(plot),
    season = as.character(season)
  )
simulation_limit_both <- unique(crossing_data_both$max_steps)
stopifnot(length(simulation_limit_both) == 1L)

crossing_summary_both <- crossing_data_both %>%
  dplyr::group_by(scenario, Endpoint, Season) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    n_observed = sum(event %in% TRUE),
    n_censored = sum(event %in% FALSE),
    n_unusable = sum(is.na(event)),
    percent_censored = if (sum(!is.na(event)) > 0) {
      100 * sum(event %in% FALSE) / sum(!is.na(event))
    } else {
      NA_real_
    },
    mean_observed = if (any(event %in% TRUE)) {mean(time[event %in% TRUE])
    } else {
      NA_real_
    },
    median_observed = if (any(event %in% TRUE)) {median(time[event %in% TRUE])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )
print(crossing_summary_both, width = Inf)

# 4. Plot crossing times for each scenario ----
plot_annotations_both <- crossing_summary_both %>%
  dplyr::mutate(
    median_text = dplyr::if_else(
      is.na(median_observed),
      "NA",
      format(median_observed, trim = TRUE)
    ),
    label = paste0(
      "Observed: ", n_observed,
      "\nMedian (observed): ", median_text, " steps",
      "\nNo event by ", simulation_limit_both, ": ",
      n_censored,
      "\nUnresolved: ", n_unusable
    )
  )

# Retain your intervals, adapting the last boundary if needed.
time_breaks_both <- sort(unique(c(
  c(0, 5, 10, 20, 50, 100, 150, 200)[c(0, 5, 10, 20, 50, 100, 150, 200) < simulation_limit_both],
  simulation_limit_both
)))

time_labels_both <- paste0(head(time_breaks_both, -1), "–", tail(time_breaks_both, -1))

crossing_binned_both <- crossing_data_both %>%
  dplyr::filter(event %in% TRUE) %>%
  dplyr::mutate(
    Time_interval = cut(
      time,
      breaks = time_breaks_both,
      labels = time_labels_both,
      include.lowest = TRUE,
      right = TRUE
    )
  )

# Include zero-count bins, even in panels with no observed events.
crossing_counts_both <- crossing_binned_both %>%
  dplyr::count(scenario, Endpoint, Season, Time_interval) %>%
  tidyr::complete(
    scenario, Endpoint, Season, Time_interval,
    fill = list(n = 0L)
  )

crossing_distribution_without <- ggplot(
  dplyr::filter(
    crossing_counts_both,
    scenario == "Without facilitation"
  ),
  aes(x = Time_interval, y = n, fill = Endpoint)
) +
  geom_col(
    width = 0.85,
    colour = "white",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  geom_text(
    data = dplyr::filter(
      plot_annotations_both,
      scenario == "Without facilitation"
    ),
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    size = 3,
    lineheight = 1.15
  ) +
  facet_grid(
    Endpoint ~ Season,
    scales = "free_y",
    drop = FALSE,
    labeller = labeller(
      Endpoint = c(
        "First species below threshold" = "First species\nbelow threshold",
        "Half the species ever below threshold" = "Half the species ever\nbelow threshold"
      )
    )
  )  +
  scale_fill_manual(
    values = c(
      "First species below threshold" = "#4477AA",
      "Half the species ever below threshold" = "#AA3377"
    ),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.45))
  ) +
  labs(
    title = "Without facilitation",
    x = "Time (simulation steps)",
    y = "Number of plots"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    panel.spacing = grid::unit(1, "lines")
  )
crossing_distribution_without

crossing_distribution_with <- ggplot(
  dplyr::filter(
    crossing_counts_both,
    scenario == "With facilitation"
  ),
  aes(x = Time_interval, y = n, fill = Endpoint)
) +
  geom_col(
    width = 0.85,
    colour = "white",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  geom_text(
    data = dplyr::filter(
      plot_annotations_both,
      scenario == "With facilitation"
    ),
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    size = 3,
    lineheight = 1.15
  ) +
  facet_grid(
    Endpoint ~ Season,
    scales = "free_y",
    drop = FALSE,
    labeller = labeller(
      Endpoint = c(
        "First species below threshold" = "First species\nbelow threshold",
        "Half the species ever below threshold" = "Half the species ever\nbelow threshold"
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "First species below threshold" = "#4477AA",
      "Half the species ever below threshold" = "#AA3377"
    ),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.45))
  ) +
  labs(
    title = "With facilitation",
    x = "Time (simulation steps)",
    y = "Number of plots"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    panel.spacing = grid::unit(1, "lines")
  )
crossing_distribution_with

# 5. Match the predictors and inspect the model sample ----
predictor_columns <- c("plot", "season", "N", "F", "sta.Omega", "SD", "Month", "Block", "composition", "CE")
predictors_without <- data.scaled.both.seasons %>%
  select(all_of(predictor_columns)) %>%
  mutate(plot = as.character(plot), season = as.character(season), scenario = "Without facilitation")
predictors_with <- data.scaled.both.seasons_fac %>%
  select(all_of(predictor_columns)) %>%
  mutate(plot = as.character(plot), season = as.character(season), scenario = "With facilitation")
extinction_predictors_both <- bind_rows(predictors_without, predictors_with)
join_keys <- c("scenario", "plot", "season")

unmatched_predictors <- ricker_start_summary %>%
  mutate(plot = as.character(plot), season = as.character(season)) %>%
  anti_join(extinction_predictors_both, by = join_keys) %>% distinct(scenario, plot, season)
if (nrow(unmatched_predictors)) print(unmatched_predictors, n = Inf)

df_extinction_both <- crossing_data_both %>%
  dplyr::mutate(scenario = as.character(scenario)) %>%
  dplyr::left_join(
    extinction_predictors_both,
    by = c("scenario", "plot", "season")
  ) %>%
  dplyr::mutate(
    Block = factor(Block),
    composition = factor(composition),
    plot = factor(plot)
  )

# Identify records with usable outcomes and model predictors.
numeric_predictors <- c("N", "F", "sta.Omega", "SD", "Month", "CE")

df_extinction_both <- df_extinction_both %>%
  dplyr::mutate(
    usable_model =
      !is.na(event) &
      is.finite(time) & time > 0 &
      dplyr::if_all(
        dplyr::all_of(numeric_predictors),
        ~ is.finite(.x)
      ) &
      !is.na(Block) & !is.na(composition) & !is.na(plot)
  )

df_extinction_both %>%
  dplyr::group_by(scenario, Endpoint) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_usable = sum(usable_model),
    n_observed_usable = sum(usable_model & event %in% TRUE),
    n_censored_usable = sum(usable_model & event %in% FALSE),
    .groups = "drop"
  ) %>% print()

lmm_extinction_data <- df_extinction_both %>% dplyr::filter(usable_model, !is.na(event)) %>% dplyr::mutate(log_time = log(time))
table(lmm_extinction_data$scenario, lmm_extinction_data$endpoint)
# 6. Fit and inspect each observed-crossing LMM ----
dat_first_without <- lmm_extinction_data %>% filter(scenario == "Without facilitation", endpoint == "first") %>% droplevels()
m_first_without <- lmerTest::lmer(log(time) ~ sta.Omega + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                                  data = dat_first_without, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_first_without)
Anova(m_first_without, test="F")
qqnorm(resid(m_first_without))
qqline(resid(m_first_without))
check_overdispersion(m_first_without)#No overdispersion
check_normality(m_first_without)#not normally distributed
m_first_without1 <- lmerTest::lmer(log(time) ~ CE + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                                  data = dat_first_without, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_first_without1)
ND_persistence_plot1 <- plot_model(
  m_first_without, type = "pred", terms = "sta.Omega [all]",
  pred.type = "fe", show.data = FALSE,
  ci.lvl = 0.95, transform = "exp", colors = "black"
) + labs(x = "Standardized niche differences",
         y = "Predicted time to first species loss\n(simulation steps)",
         title = NULL
) + theme_classic() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    plot.caption = element_text(size = 9, hjust = 0)
  )
ND_persistence_plot1
##
dat_half_without <- lmm_extinction_data %>% filter(scenario == "Without facilitation", endpoint == "half") %>% droplevels()
hist(log(dat_half_without$time))
m_half_without <- lmerTest::lmer(log(time) ~ sta.Omega  + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                                 data = dat_half_without, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_half_without)
Anova(m_half_without, test="F")
qqnorm(resid(m_half_without))
qqline(resid(m_half_without))
check_overdispersion(m_half_without)#No overdispersion
check_normality(m_half_without)#not normally distributed
ND_persistence_plot2 <- plot_model(
  m_half_without, type = "pred", terms = "sta.Omega [all]",
  pred.type = "fe", show.data = FALSE,
  ci.lvl = 0.95, transform = "exp", colors = "black"
) + labs(x = "Standardized niche differences",
         y = "Predicted time to half species loss\n(simulation steps)",
         title = NULL
) + theme_classic() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    plot.caption = element_text(size = 9, hjust = 0)
  )
ND_persistence_plot2

m_half_without1 <- lmerTest::lmer(log(time) ~ CE + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                                 data = dat_half_without, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_half_without1)
##
dat_first_with <- lmm_extinction_data %>% filter(scenario == "With facilitation", endpoint == "first") %>% droplevels()
hist(dat_first_with$time)
m_first_with <- lmerTest::lmer(log(time) ~ sta.Omega + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                               data = dat_first_with, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_first_with)
Anova(m_first_with, test="F")
qqnorm(resid(m_first_with))
qqline(resid(m_first_with))
check_overdispersion(m_first_with)#No overdispersion
check_normality(m_first_with)#not normally distributed
m_first_with1 <- lmerTest::lmer(log(time) ~ CE + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                               data = dat_first_with, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_first_with1)

dat_half_with <- lmm_extinction_data %>% filter(scenario == "With facilitation", endpoint == "half") %>% droplevels()
hist(log(dat_half_with$time))
m_half_with <- lmerTest::lmer(log(time) ~ sta.Omega + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                              data = dat_half_with, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_half_with)
Anova(m_half_with, test="F")
qqnorm(resid(m_half_with))
qqline(resid(m_half_with))
check_overdispersion(m_half_with)#No overdispersion
check_normality(m_half_with)

m_half_with1 <- lmerTest::lmer(log(time) ~ CE + Month + (1 | Block) + (1 | composition) + (1 | plot), 
                              data = dat_half_with, control = lme4::lmerControl(optimizer = "bobyqa"))
summary(m_half_with1)

combined_crossing_distribution <- (crossing_distribution_without / (crossing_distribution_with)) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain")  )
combined_crossing_distribution

##This is fig. S5
ggsave(filename = "figures/combined_crossing_distribution.tiff",
       plot = combined_crossing_distribution, width = 8, height = 12,
       units = "in", dpi = 300, compression = "lzw")

