##scripts/10.Extra_Feasibility_domain_all_pairs

##Results for Table S4 and Fig.S12, S13
source("functions/useful SEM constrain functions.R")
source("functions/useful SEM partial plot functions.R")
library(lavaan)
library(corrplot)
library(ggpubr)
##To answer the reviewer's question on the contribution of indirect interactions
##The contribution of indirect interactions D = Omega - Omega_all, Omega_all is the feasibility domian for all pairs
##we decided to use Omega_all in the SEM, if ND~CE relationship disappears, indicating indirect interactions may be important for CE
head(data.scaled.both.seasons)
names(data.scaled.both.seasons)[17]<-"sta.Omega_all_pairs"
##SEM model without indirect interactions (the SEM for the major revision)
model_string_both_seasons_N_Omega_all <- '
CE ~ F + SD + sownSLA + sownMPD + sta.Omega_all_pairs + sta.theta 
SE_inter ~  F + SD + sownSLA + sownMPD + sta.Omega_all_pairs + sta.theta 
SE_intra ~ F + SD + sownSLA + sownMPD + sta.Omega_all_pairs + sta.theta 
sta.Omega_all_pairs ~ F + SD + sownSLA + sownMPD 
sta.theta ~ F + SD + sownSLA + sownMPD 

sownSLA~~sownMPD
sta.Omega_all_pairs~~sta.theta 
'
fituw_both_seasons_N_Omega_all <- sem(model_string_both_seasons_N_Omega_all, data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_Omega_all, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.027
summary(fituw_both_seasons_N_Omega_all, standardized =TRUE, rsquare=TRUE)

##now I need to constrain the SEM model
constrain_model_both_seasons_N_Omega_all <- test_stepwise_constraints_cfa(fituw_both_seasons_N_Omega_all, data.scaled.both.seasons)
constrain_model_both_seasons_N_Omega_all$final_fit
best_string_both_seasons_N_Omega_all <- constrain_model_both_seasons_N_Omega_all$final_model_string

confituw_both_seasons_N_Omega_all <- sem(unlist(best_string_both_seasons_N_Omega_all), data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_Omega_all, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_Omega_all, c("pvalue","cfi","rmsea","srmr"))##p:0.13, cfi:0.971, rmsea:0.059 ,srmr0.075
sort.pars(confituw_both_seasons_N_Omega_all, response = "CE", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_Omega_all, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_Omega_all, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_Omega_all, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_Omega_all, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#Constrain model mannually
#sownSLA on sta.theta N vs. NoN
mcs1_sown_N_Omega_all_1.0 <- constrain_path_cfa(confituw_both_seasons_N_Omega_all, "sta.theta ~ sownSLA", coef_name = "coef244N", groups_constrain = c(1,3))
fc1a_sown_N_Omega_all_1.0 <- sem(mcs1_sown_N_Omega_all_1.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_1.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_1.0, "sta.theta ~ sownSLA", coef_name = "coef244NoN", groups_constrain = c(2,4))
fc1_sown_N_Omega_all_1.0 <- sem(mcs2_sown_N_Omega_all_1.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_1.0)
lavTestLRT(fc1_sown_N_Omega_all_1.0, confituw_both_seasons_N_Omega_all)#accepted
sort.pars(fc1_sown_N_Omega_all_1.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_1.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_1.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_1.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_1.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))
fitMeasures(fc1_sown_N_Omega_all_1.0, c("pvalue","cfi","rmsea","srmr"))## p = 0.093, cfi:0.965, rmsea:0.063 ,srmr0.079

#sownMPD on CE N vs. NoN
mcs1_sown_N_Omega_all_2.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_1.0, "CE ~ sownMPD", coef_name = "coef245N", groups_constrain = c(1,3))
fc1a_sown_N_Omega_all_2.0 <- sem(mcs1_sown_N_Omega_all_2.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_2.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_2.0, "CE ~ sownMPD", coef_name = "coef245NoN", groups_constrain = c(2,4))
fc1_sown_N_Omega_all_2.0 <- sem(mcs2_sown_N_Omega_all_2.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_2.0)
lavTestLRT(fc1_sown_N_Omega_all_1.0, fc1_sown_N_Omega_all_2.0)#rejected

#sownMPD on SE_inter N vs. NoN
mcs1_sown_N_Omega_all_3.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_1.0, "SE_inter ~ sownMPD", coef_name = "coef245N", groups_constrain = c(1,3))
fc1a_sown_N_Omega_all_3.0 <- sem(mcs1_sown_N_Omega_all_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_3.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_3.0, "SE_inter ~ sownMPD", coef_name = "coef245NoN", groups_constrain = c(2,4))
fc1_sown_N_Omega_all_3.0 <- sem(mcs2_sown_N_Omega_all_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_3.0)
lavTestLRT(fc1_sown_N_Omega_all_1.0, fc1_sown_N_Omega_all_3.0)#accepted
sort.pars(fc1_sown_N_Omega_all_3.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_3.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_3.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_3.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_3.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#F on CE June vs. August
mcs1_sown_N_Omega_all_4.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_3.0, "CE ~ F", coef_name = "coef246Aug", groups_constrain = c(1,2))
fc1a_sown_N_Omega_all_4.0 <- sem(mcs1_sown_N_Omega_all_4.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_4.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_4.0, "CE ~ F", coef_name = "coef246Jun", groups_constrain = c(3,4))
fc1_sown_N_Omega_all_4.0 <- sem(mcs2_sown_N_Omega_all_4.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_4.0)
lavTestLRT(fc1_sown_N_Omega_all_4.0, fc1_sown_N_Omega_all_3.0)#accepted
sort.pars(fc1_sown_N_Omega_all_4.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#SD on CE June vs. August
mcs1_sown_N_Omega_all_5.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_4.0, "CE ~ SD", coef_name = "coef247N", groups_constrain = c(1,3))
fc1a_sown_N_Omega_all_5.0 <- sem(mcs1_sown_N_Omega_all_5.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_5.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_5.0, "CE ~ SD", coef_name = "coef247NoN", groups_constrain = c(2,4))
fc1_sown_N_Omega_all_5.0 <- sem(mcs2_sown_N_Omega_all_5.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_5.0)
lavTestLRT(fc1_sown_N_Omega_all_4.0, fc1_sown_N_Omega_all_5.0)#accepted
sort.pars(fc1_sown_N_Omega_all_5.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_5.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_5.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_5.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_5.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#F on SE_inter June vs. August
mcs1_sown_N_Omega_all_6.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_5.0, "SE_inter ~ F", coef_name = "coef248Aug", groups_constrain = c(1,2))
fc1a_sown_N_Omega_all_6.0 <- sem(mcs1_sown_N_Omega_all_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_6.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_6.0, "SE_inter ~ F", coef_name = "coef248Jun", groups_constrain = c(3,4))
fc1_sown_N_Omega_all_6.0 <- sem(mcs2_sown_N_Omega_all_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_6.0)
lavTestLRT(fc1_sown_N_Omega_all_6.0, fc1_sown_N_Omega_all_5.0)#accepted
sort.pars(fc1_sown_N_Omega_all_6.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_6.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_6.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_6.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_6.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#sta.Omega_all_pairs on N vs. NoN
mcs1_sown_N_Omega_all_7.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_6.0, "SE_inter ~ sta.Omega_all_pairs", coef_name = "coef249N", groups_constrain = c(1,3))
fc1a_sown_N_Omega_all_7.0 <- sem(mcs1_sown_N_Omega_all_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_Omega_all_7.0 <- constrain_path_cfa(fc1a_sown_N_Omega_all_7.0, "SE_inter ~ sta.Omega_all_pairs", coef_name = "coef249NoN", groups_constrain = c(2,4))
fc1_sown_N_Omega_all_7.0 <- sem(mcs2_sown_N_Omega_all_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_Omega_all_7.0)
lavTestLRT(fc1_sown_N_Omega_all_7.0, fc1_sown_N_Omega_all_6.0)#accepted
sort.pars(fc1_sown_N_Omega_all_7.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_7.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_7.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_7.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_Omega_all_7.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#SD on SE_intra all four groups
mcs1_sown_N_Omega_all_8.0 <- constrain_path_cfa(fc1_sown_N_Omega_all_7.0, "SE_intra ~ SD", coef_name = "coef250", groups_constrain = c(1,2,3,4))
fc1a_sown_N_Omega_all_8.0 <- sem(mcs1_sown_N_Omega_all_8.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1a_sown_N_Omega_all_8.0)
lavTestLRT(fc1a_sown_N_Omega_all_8.0, fc1_sown_N_Omega_all_7.0)#accepted
sort.pars(fc1a_sown_N_Omega_all_8.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1a_sown_N_Omega_all_8.0, response = "sta.Omega_all_pairs", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1a_sown_N_Omega_all_8.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1a_sown_N_Omega_all_8.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1a_sown_N_Omega_all_8.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#the best model is fc1a_sown_N_Omega_all_8.0
fitMeasures(fc1a_sown_N_Omega_all_8.0, c("pvalue","cfi","rmsea","srmr"))## p = 0.052, cfi:0.951, rmsea:0.068 ,srmr0.086
summary(fc1a_sown_N_Omega_all_8.0, standardized =TRUE, rsquare=TRUE)
##make a table in word and a figure in ppt to report this results

library(dplyr)
library(ggplot2)
library(smatr)
# 1. Prepare the paired observations from the model without facilitation.
# Use the null-standardized metrics before the additional SEM scaling.
# Keep the 4- and 8-species plots without facilitation.
omega_pair_full_data <- df_ini_coexistence_bio.effects1 %>%
  dplyr::filter(
    SD %in% c(4, 8),
    is.finite(sta.Omega_all_pairs_zero),
    is.finite(sta.Omega)
  ) %>%
  dplyr::mutate(
    richness_group = factor(
      SD, levels = c(4, 8), labels = c("4 species plots", "8 species plots")
    )
  )

# Fit the orthogonal regression to the same observations shown in the plot.
omega_pair_full_MA <- smatr::sma(sta.Omega ~ sta.Omega_all_pairs_zero, method = "MA", data = omega_pair_full_data)
omega_MA_4 <- smatr::sma(sta.Omega ~ sta.Omega_all_pairs_zero, data = subset(omega_pair_full_data, SD == 4), method = "MA")
omega_MA_8 <- smatr::sma(sta.Omega ~ sta.Omega_all_pairs_zero, data = subset(omega_pair_full_data, SD == 8), method = "MA")

summary(omega_pair_full_MA)
summary(omega_MA_4)
summary(omega_MA_8)
omega_pair_full_coef <- coef(omega_pair_full_MA)
coef_4 <- coef(omega_MA_4)
coef_8 <- coef(omega_MA_8)

# A single-group MA fit returns the elevation (intercept), then the slope.
coef_4 <- as.numeric(unlist(coef(omega_MA_4)))
coef_8 <- as.numeric(unlist(coef(omega_MA_8)))

omega_lines_SD <- tibble::tibble(
  SD = c(4, 8), intercept = c(coef_4[1], coef_8[1]), slope = c(coef_4[2], coef_8[2])
)

# Calculate each fitted line over the observed x-range for its richness group.
omega_segments_SD <- omega_pair_full_data %>%
  dplyr::group_by(SD, richness_group) %>%
  dplyr::summarise(
    x_start = min(sta.Omega_all_pairs_zero),
    x_end   = max(sta.Omega_all_pairs_zero),
    .groups = "drop"
  ) %>%
  dplyr::left_join(omega_lines_SD, by = "SD") %>%
  dplyr::mutate(
    y_start = intercept + slope * x_start,
    y_end   = intercept + slope * x_end
  )

omega_pair_full_plot_facet <- ggplot(
  omega_pair_full_data, aes(x = sta.Omega_all_pairs_zero, y = sta.Omega)) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed", linewidth = 0.7
  ) + geom_point(colour = "black", size = 2.5, alpha = 0.5) +
  geom_segment(
    data = omega_segments_SD,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    colour = "black", linewidth = 0.9
  ) +
  facet_wrap(~ richness_group, nrow = 1, scales = "fixed") +
  coord_equal() +
  labs(
    x = "Standardized pairwise feasibility\nwithout facilitation",
    y = "Standardized full-community feasibility\nwithout facilitation"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    strip.text = element_text(size = 11, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing.x = grid::unit(1, "lines")
  )
omega_pair_full_plot_facet

omega_equations_SD <- omega_segments_SD %>%
  dplyr::transmute(
    richness_group,
    equation = sprintf("y = %.3f x %+.3f", slope, intercept)
  )
print(omega_equations_SD)

# 1. Choose a common x-axis window for both richness groups.
omega_x_upper <- unname(quantile(
  omega_pair_full_data$sta.Omega_all_pairs_zero,
  probs = 0.99,
  na.rm = TRUE
))

omega_x_limits <- c(
  min(omega_pair_full_data$sta.Omega_all_pairs_zero, na.rm = TRUE),
  omega_x_upper
)

# Check how many observations the zoom hides in each panel.
omega_pair_full_data %>%
  dplyr::group_by(richness_group) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_above_window = sum(sta.Omega_all_pairs_zero > omega_x_upper),
    .groups = "drop"
  )

# 2. Zoom the display while retaining the lines fitted to all observations.
omega_pair_full_plot_zoom <- omega_pair_full_plot_facet +
  coord_cartesian(xlim = omega_x_limits) +
  theme(
    aspect.ratio = 0.9,
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 13, colour = "black"),
    panel.spacing.x = grid::unit(1, "lines")
  )

omega_pair_full_plot_zoom

omega_pair_full_plot_zoom_eq <- omega_pair_full_plot_zoom +
  geom_text(
    data = omega_equations_SD,
    aes(x = -Inf, y = Inf, label = equation),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.5,
    size = 3.5,
    colour = "black"
  )
# 3. Export at a size suitable for a full-width supplementary figure S13.
ggsave(filename = "figures/omega_pair_full_by_richness_zoom.pdf",
  plot = omega_pair_full_plot_zoom_eq, width = 190, height = 110, units = "mm"
)
ggsave(filename = "figures/omega_pair_full_by_richness_zoom.tif",
       plot = omega_pair_full_plot_zoom_eq, width = 190, height = 110, units = "mm"
)
