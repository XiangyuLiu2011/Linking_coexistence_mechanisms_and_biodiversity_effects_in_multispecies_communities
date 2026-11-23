##scripts/08_Visualization_of_the_contribution_of_indirect_interactions
##As the reviewer requested, sta.differential and observed community-pair differential (the contribution of indirect interactions) are visualized.
##This is code for figure S11
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# I followed the Figure 2 that used for diversity effects
error_type <- "SD"  
season_colours <- c("June" = "#009E73","August" = "#E69F00")

# Define the eight panels and their order.
panel_layout <- tibble::tribble(
  ~Panel, ~Metric,        ~SD, ~Facilitation,
  "A",    "Observed",       4, "Without facilitation",
  "B",    "Observed",       8, "Without facilitation",
  "C",    "Standardized",   4, "Without facilitation",
  "D",    "Standardized",   8, "Without facilitation",
  "E",    "Observed",       4, "With facilitation",
  "F",    "Observed",       8, "With facilitation",
  "G",    "Standardized",   4, "With facilitation",
  "H",    "Standardized",   8, "With facilitation"
)

# ------------------------------------------------------------
# Prepare individual observations of community-pair differential
D_observed_plot_data <- two_coexistence_dfs %>%
  dplyr::transmute(
    plot, season = factor(as.character(season), levels = c("June", "August")),
    SD = as.numeric(as.character(SD)),
    differential = as.numeric(as.character(differential)),
    differential_fac = as.numeric(as.character(differential_fac))
  ) %>%
  dplyr::filter(SD %in% c(4, 8)) %>%
  tidyr::pivot_longer(
    cols = c(differential, differential_fac),
    names_to = "Facilitation",
    values_to = "D"
  ) %>%
  dplyr::filter(is.finite(D)) %>%
  dplyr::mutate(
    Facilitation = factor(
      Facilitation,
      levels = c("differential", "differential_fac"),
      labels = c("Without facilitation", "With facilitation")
    ),
    Richness = factor(
      SD, levels = c(4, 8),
      labels = c("4 species", "8 species")
    ),
    Panel = factor(
      paste(Facilitation, Richness, sep = " | "),
      levels = c(
        "Without facilitation | 4 species",
        "Without facilitation | 8 species",
        "With facilitation | 4 species",
        "With facilitation | 8 species"
      )
    )
  )

D_points <- dplyr::bind_rows(D_observed_plot_data %>%
    dplyr::transmute(plot, season, SD, Facilitation, D,
      Metric = "Observed"),
  D_plot_data %>%
    dplyr::transmute(plot, season, SD, Facilitation, D,
      Metric = "Standardized")
) %>% dplyr::mutate(
    SD = as.numeric(as.character(SD)),
    Facilitation = as.character(Facilitation),
    season = factor(as.character(season), levels = c("June", "August"))
  ) %>% dplyr::filter(SD %in% c(4, 8),
    is.finite(D), !is.na(season)
  ) %>% dplyr::left_join(
    panel_layout, by = c("Metric", "SD", "Facilitation")
  ) %>% dplyr::mutate(
    # Offset raw observations and summaries slightly.
    x_raw = as.numeric(season) - 0.08,
    x_summary = as.numeric(season) + 0.12
  )

# ------------------------------------------------------------
# 3. Calculate mean, SD and SE.
# ------------------------------------------------------------
D_summary <- D_points %>%
  dplyr::group_by(Panel, season) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_D = mean(D),
    sd_D = sd(D),
    se_D = sd_D / sqrt(n),
    percent_positive = 100 * mean(D > 0),
    x_summary = dplyr::first(x_summary),
    .groups = "drop"
  ) %>% dplyr::mutate(
    error_D = if (error_type == "SE") se_D else sd_D,
    positive_label = sprintf("%.1f%% > 0", percent_positive)
  )

# ------------------------------------------------------------
# 4. Shared plotting style.
# ------------------------------------------------------------
D_point_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, face = "plain", hjust = 0),
    plot.tag = element_text(size = 10, face = "plain"),
    legend.position = "none",
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 12, 8, 8)
  )

# ------------------------------------------------------------
# Shared plot layers.
# The summary data are supplied separately to each plot.
# ------------------------------------------------------------
D_common_layers <- list(
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4),
  geom_errorbar(
    aes(x = x_summary, ymin = mean_D - error_D, ymax = mean_D + error_D
    ), width = 0.08, linewidth = 0.6, colour = "black"),
  geom_point(
    aes(x = x_summary, y = mean_D, fill = season),
    shape = 23,size = 3,stroke = 0.6,colour = "black"),
  geom_text(aes(x = as.numeric(season), y = Inf, label = positive_label),
    vjust = 1.2, size = 3.2, colour = "black", show.legend = FALSE),
  scale_x_continuous(breaks = c(1, 2),
    labels = c("June", "August"), limits = c(0.6, 2.4), expand = expansion(mult = 0)),
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))  ),
  scale_fill_manual(values = season_colours,
    breaks = c("June", "August"), drop = FALSE, name = NULL),
  D_point_theme
)

# ------------------------------------------------------------
# A: Observed, 4 species, without facilitation.
# ------------------------------------------------------------
D_A <- ggplot(D_summary %>% dplyr::filter(Panel == "A")) +
  geom_point(data = D_points %>% dplyr::filter(Panel == "A"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(
    title = "4 species plots without facilitation",
    x = NULL,
    y = "Observed community-pair differential"
  )
D_A

# ------------------------------------------------------------
# B: Observed, 8 species, without facilitation.
# ------------------------------------------------------------
D_B <- ggplot(D_summary %>% dplyr::filter(Panel == "B")
) + geom_point(
    data = D_points %>% dplyr::filter(Panel == "B"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) + D_common_layers +
  labs(
    title = "8 species plots without facilitation",
    x = NULL, y = ""
  )
D_B

# ------------------------------------------------------------
# C: Standardized, 4 species, without facilitation.
# ------------------------------------------------------------
D_C <- ggplot(D_summary %>% dplyr::filter(Panel == "C")) +
  geom_point(
    data = D_points %>% dplyr::filter(Panel == "C"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(title = "4 species plots without facilitation",
    x = NULL, y = "Standardized community-pair differential"
  )
D_C

# ------------------------------------------------------------
# D: Standardized, 8 species, without facilitation.
# ------------------------------------------------------------
D_D <- ggplot(
  D_summary %>% dplyr::filter(Panel == "D")
) +
  geom_point(
    data = D_points %>% dplyr::filter(Panel == "D"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(title = "8 species plots without facilitation",
    x = NULL, y = "")
D_D

# ------------------------------------------------------------
# E: Observed, 4 species, with facilitation.
# ------------------------------------------------------------
D_E <- ggplot(D_summary %>% dplyr::filter(Panel == "E")
) + geom_point(
    data = D_points %>% dplyr::filter(Panel == "E"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(
    title = "4 species plots with facilitation",
    x = NULL,
    y = "Observed community-pair differential"
  )

D_E

# ------------------------------------------------------------
# F: Observed, 8 species, with facilitation.
# ------------------------------------------------------------
D_F <- ggplot(D_summary %>% dplyr::filter(Panel == "F")
) + geom_point(
    data = D_points %>% dplyr::filter(Panel == "F"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(
    title = "8 species plots with facilitation",
    x = NULL, y = ""  )

D_F

# ------------------------------------------------------------
# G: Standardized, 4 species, with facilitation.
# ------------------------------------------------------------
D_G <- ggplot(
  D_summary %>% dplyr::filter(Panel == "G")
) +
  geom_point(
    data = D_points %>% dplyr::filter(Panel == "G"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  D_common_layers +
  labs(
    title = "4 species plots with facilitation",
    x = NULL, y = "Standardized community-pair differential"
  )
D_G

# ------------------------------------------------------------
# H: Standardized, 8 species, with facilitation.
# ------------------------------------------------------------
D_H <- ggplot(D_summary %>% dplyr::filter(Panel == "H")) +
  geom_point(
    data = D_points %>% dplyr::filter(Panel == "H"),
    aes(x = x_raw, y = D, fill = season),
    position = position_jitter(width = 0.065, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    colour = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) + D_common_layers +
  labs(
    title = "8 species plots with facilitation",
    x = NULL,  y = "")
D_H

# ------------------------------------------------------------
# Combine the eight independently editable plots.
# ------------------------------------------------------------
D_distribution_combined <- (
  patchwork::wrap_plots(
    D_A, D_B,
    D_C, D_D,
    D_E, D_F,
    D_G, D_H,
    ncol = 2, widths = c(1, 1), heights = c(1, 1, 1, 1), guides = "collect") + 
    patchwork::plot_annotation(tag_levels = "A"))
D_distribution_combined
##This figure S11
ggsave(filename = "figures/Community_pair_differential_8panels.pdf",plot = D_distribution_combined, width = 12, height = 14,  units = "in")
