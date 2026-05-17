##scripts/06_Linking_coexistence_mechanisms_with_biodiversity_effects_including_facilitation

source("functions/useful SEM constrain functions.R")
source("functions/useful SEM partial plot functions.R")
library(lavaan)
library(corrplot)
library(ggpubr)
library(broom)
#merge two datasets

df_ini_coexistence_bio.effects_fac <- merge(combined_sown_coex_fac, df_biodiversity.effects_c[,c(1:2,8:15)], by=c("season", "plot", "Block"))
df_ini_coexistence_bio.effects_fac$Month<-1*(df_ini_coexistence_bio.effects_fac$season=="August")

##exclude SD==20, because all 20 species plots have the same predicted coexistence_mechanisms
df_ini_coexistence_bio.effects1_fac<-subset(df_ini_coexistence_bio.effects_fac, SD!=20)
head(df_ini_coexistence_bio.effects1_fac)

Scalevars_both_seasons_fac <-df_ini_coexistence_bio.effects1_fac[,c("Month", "N", "F", "SD", "sownSLA", "sownMPD", "sta.Omega", "sta.theta", "sta.differential", "NE", "CE", "SE", "SE_inter", "SE_intra")]
# do the scaling
Scalevars_both_seasons2_fac <- scale(Scalevars_both_seasons_fac)
corrs.2_fac<-cor(Scalevars_both_seasons2_fac, use = "pairwise.complete.obs")
corrplot(corrs.2_fac, method = "square", type = "upper", order = "hclust",
  col = colorRampPalette(c("#336699", "white", "#990000"))(200),
  addCoef.col = "black", number.cex = 0.65, tl.col = "black",
  tl.cex = 0.8, tl.srt = 45, diag = FALSE, mar = c(0, 0, 1, 0) )

cat.vars.both.seasons_fac <- df_ini_coexistence_bio.effects1_fac[,c("plot", "Block", "composition", "season", "N")]
names(cat.vars.both.seasons_fac)[5] <- c("Nitrogen")
cat.vars.both.seasons_fac$group <- interaction(cat.vars.both.seasons_fac$season, cat.vars.both.seasons_fac$Nitrogen)
data.scaled.both.seasons_fac <- cbind(cat.vars.both.seasons_fac, Scalevars_both_seasons2_fac)
head(data.scaled.both.seasons_fac)
data.scaled.both.seasons_fac$group<-as.factor(data.scaled.both.seasons_fac$group)

fituw_both_seasons_N_indirect_fac <- sem(model_string_both_seasons_N_indirect, data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect_fac, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect_fac, standardized =TRUE, rsquare=TRUE)
##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect_fac <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect_fac, data.scaled.both.seasons_fac)
constrain_model_both_seasons_N_indirect_fac$final_fit
best_string_both_seasons_N_indirect_fac <- constrain_model_both_seasons_N_indirect_fac$final_model_string

confituw_both_seasons_N_indirect_fac <- sem(unlist(best_string_both_seasons_N_indirect_fac), data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect_fac, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect_fac, c("cfi","rmsea","srmr"))## cfi:0.98, rmsea:0.05 ,srmr0.08
sort.pars(confituw_both_seasons_N_indirect_fac, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "SE_intra", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.theta", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.differential", group.names = unique(data.scaled.both.seasons_fac$group))
#Constrain model mannually
#F on SEAS (dominance effects) August vs. June
mcs1_sown_N_indirect_1.0_fac <- constrain_path_cfa(confituw_both_seasons_N_indirect_fac, "SE_inter ~ F", coef_name = "coef244Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_1.0_fac <- sem(mcs1_sown_N_indirect_1.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_1.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_1.0_fac, "SE_inter ~ F", coef_name = "coef244Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_1.0_fac <- sem(mcs2_sown_N_indirect_1.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_1.0_fac)
lavTestLRT(fc1_sown_N_indirect_1.0_fac, confituw_both_seasons_N_indirect_fac)#accepted
sort.pars(fc1_sown_N_indirect_1.0_fac, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))

#sta.Omega on SEAS (dominance effects) August vs. June
mcs1_sown_N_indirect_2.0.1_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "SE_inter ~ sta.Omega", coef_name = "coef245Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_2.0.1_fac <- sem(mcs1_sown_N_indirect_2.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0.1_fac <- constrain_path_cfa(fc1a_sown_N_indirect_2.0.1_fac, "SE_inter ~ sta.Omega", coef_name = "coef245Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_2.0.1_fac <- sem(mcs2_sown_N_indirect_2.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0.1_fac)
lavTestLRT(fc1_sown_N_indirect_1.0_fac, fc1_sown_N_indirect_2.0.1_fac)#rejected!

#sownMPD on SEAS (dominance effects) N vs. noN
mcs1_sown_N_indirect_2.0.2_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "SE_inter ~ sownMPD", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_2.0.2_fac <- sem(mcs1_sown_N_indirect_2.0.2_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0.2_fac <- constrain_path_cfa(fc1a_sown_N_indirect_2.0.2_fac, "SE_inter ~ sownMPD", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_2.0.2_fac <- sem(mcs2_sown_N_indirect_2.0.2_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0.2_fac)
lavTestLRT(fc1_sown_N_indirect_1.0_fac, fc1_sown_N_indirect_2.0.2_fac)#rejected!

#sownMPD on CE (complementarity effects) Aug vs. June
mcs1_sown_N_indirect_3.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "CE ~ sownMPD", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0_fac <- sem(mcs1_sown_N_indirect_3.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_3.0_fac, "CE ~ sownMPD", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0_fac <- sem(mcs2_sown_N_indirect_3.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0_fac)
lavTestLRT(fc1_sown_N_indirect_3.0_fac, fc1_sown_N_indirect_1.0_fac)#rejected

#SD on CE (complementarity effects) Aug vs. June
mcs1_sown_N_indirect_4.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "CE ~ SD", coef_name = "coef247N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_4.0_fac <- sem(mcs1_sown_N_indirect_4.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_4.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_4.0_fac, "CE ~ SD", coef_name = "coef247noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_4.0_fac <- sem(mcs2_sown_N_indirect_4.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_4.0_fac)
lavTestLRT(fc1_sown_N_indirect_4.0_fac, fc1_sown_N_indirect_1.0_fac)#accepted

#sownMPD on SEPS June vs. August
mcs1_sown_N_indirect_5.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_4.0_fac, "SE_intra ~ sownMPD", coef_name = "coef248Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_5.0_fac <- sem(mcs1_sown_N_indirect_5.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_5.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_5.0_fac, "SE_intra ~ sownMPD", coef_name = "coef248Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_5.0_fac <- sem(mcs2_sown_N_indirect_5.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_5.0_fac)
lavTestLRT(fc1_sown_N_indirect_5.0_fac, fc1_sown_N_indirect_1.0_fac)#rejected

#sownSLA on theta N vs. noN
mcs1_sown_N_indirect_6.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_4.0_fac, "sta.theta ~ sownSLA", coef_name = "coef249N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_6.0_fac <- sem(mcs1_sown_N_indirect_6.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_6.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_6.0_fac, "sta.theta ~ sownSLA", coef_name = "coef249noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_6.0_fac <- sem(mcs2_sown_N_indirect_6.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_6.0_fac)
lavTestLRT(fc1_sown_N_indirect_6.0_fac, fc1_sown_N_indirect_4.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_6.0_fac, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownSLA on omega N vs. noN
mcs1_sown_N_indirect_7.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_6.0_fac, "sta.Omega ~ sownSLA", coef_name = "coef251N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_7.0_fac <- sem(mcs1_sown_N_indirect_7.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_7.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_7.0_fac, "sta.Omega ~ sownSLA", coef_name = "coef251noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_7.0_fac <- sem(mcs2_sown_N_indirect_7.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_7.0_fac)
lavTestLRT(fc1_sown_N_indirect_6.0_fac, fc1_sown_N_indirect_7.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on sta.differential June vs. Aug
mcs1_sown_N_indirect_8.0.1_fac <- constrain_path_cfa(fc1_sown_N_indirect_7.0_fac, "sta.differential ~ sownMPD", coef_name = "coef252Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_8.0.1_fac <- sem(mcs1_sown_N_indirect_8.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_8.0.1_fac <- constrain_path_cfa(fc1a_sown_N_indirect_8.0.1_fac, "sta.differential ~ sownMPD", coef_name = "coef252June", groups_constrain = c(3,4))
fc1_sown_N_indirect_8.0.1_fac <- sem(mcs2_sown_N_indirect_8.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_8.0.1_fac)
lavTestLRT(fc1_sown_N_indirect_8.0.1_fac, fc1_sown_N_indirect_7.0_fac)#accepted

##fc1_sown_N_8.0.1_fac is the best model, fig.S4 made in PowerPoint based on the best SEM model
summary(fc1_sown_N_indirect_8.0.1_fac, standardized =TRUE, rsquare=TRUE)
fitMeasures(fc1_sown_N_indirect_8.0.1_fac, c("cfi","rmsea","srmr"))## cfi:0.98, rmsea:0.04 ,srmr0.08

sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "sta.Omega", sig=T, group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "sta.theta", sig=T, group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "sta.differential",sig=T, group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "CE", sig=T, group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "SE_inter", sig=T, group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_8.0.1_fac, response = "SE_intra", sig=T, group.names = unique(data.scaled.both.seasons_fac$group))

######Fig.S5
df_coex_fac <- combined_sown_coex_fac[,c(1:2,5,15,33:34)]
names(df_coex_fac)[c(2:3,5:6)] <- c("Omega_fac", "theta_fac","sta.Omega_fac", "sta.theta_fac")
two_coexistence_dfs <- merge(combined_sown_coex_new, df_coex_fac, by=c("season", "plot"))

two_coexistence_dfs %>% #group_by(season) %>%
  do({model <- lm(theta_fac ~ theta, data = .)
  data.frame(tidy(model) %>% filter(term == "theta"),
             glance(model) %>% dplyr::select(r.squared, adj.r.squared)
  ) }) %>% arrange(p.value)

two_coexistence_dfs %>% #group_by(season) %>%
  do({model <- lm(Omega_fac ~ Omega, data = .)
  data.frame(tidy(model) %>% filter(term == "Omega"),
             glance(model) %>% dplyr::select(r.squared, adj.r.squared)
  ) }) %>% arrange(p.value)

two_coexistence_dfs %>% #group_by(season) %>%
  do({model <- lm(sta.theta_fac ~ sta.theta, data = .)
  data.frame(tidy(model) %>% filter(term == "sta.theta"),
             glance(model) %>% dplyr::select(r.squared, adj.r.squared)
  ) }) %>% arrange(p.value)

two_coexistence_dfs %>% #group_by(season) %>%
  do({model <- lm(sta.Omega_fac ~ sta.Omega, data = .)
    data.frame(tidy(model) %>% filter(term == "sta.Omega"),
      glance(model) %>% dplyr::select(r.squared, adj.r.squared)
    ) }) %>% arrange(p.value)

plot.info.sla.vip_c<-plot.info.sla.vip#make a copy of traits
names(plot.info.sla.vip_c)[3] <- "plot"
two_coexistence_dfs<-merge(two_coexistence_dfs, plot.info.sla.vip_c[,c(2:3,10:11)] , by=c("plot", "Block") )

omega_corr <- ggplot(two_coexistence_dfs, aes(x = sta.Omega_fac, y = sta.Omega, color = sownSLA)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01), size = 3) +
  scale_color_gradientn(colors = c("darkgreen", "green", "yellow", "pink", "white"),
                        limits = c(15, 25),
                        name = "Sown CWM SLA") + geom_abline(intercept = 0, slope = 1) + theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold"),
                     strip.background = element_blank(), axis.title = element_text(size = 12, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"), panel.grid = element_blank(),
                     legend.position = "",panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                     axis.line = element_line(color = "black", linewidth = 0.5),
                     axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  labs(x = "Predicted ND with facilitation", y = "Predicted ND without facilitation", color = "Treatment")

omega_corr1 <- ggplot(two_coexistence_dfs, aes(x = Omega_fac, y = Omega, color = sownSLA)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01), size = 3) +
  scale_color_gradientn(colors = c("darkgreen", "green", "yellow", "pink", "white"),
                        limits = c(15, 25),
                        name = "Sown CWM SLA") + geom_abline(intercept = 0, slope = 1) + theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(), axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"), panel.grid = element_blank(),
        legend.position = "",panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  labs(x = "Observed ND with facilitation", y = "Observed ND without facilitation", color = "Treatment")

theta_corr <- ggplot(two_coexistence_dfs, aes(x = sta.theta_fac, y = sta.theta, color = sownSLA)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01), size = 3) +
  scale_color_gradientn(colors = c("darkgreen", "green", "yellow", "pink", "white"),
                        limits = c(15, 25),
                        name = "Sown CWM SLA") + geom_abline(intercept = 0, slope = 1) + theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(), axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"), panel.grid = element_blank(),
        legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  labs(x = "Predicted FD with facilitation", y = "Predicted FD without facilitation", color = "Treatment")

theta_corr1 <- ggplot(two_coexistence_dfs, aes(x = theta_fac, y = theta, color = sownSLA)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.01), size = 3) +
  scale_color_gradientn(colors = c("darkgreen", "green", "yellow", "pink", "white"),
                        limits = c(15, 25),
                        name = "Sown CWM SLA") + geom_abline(intercept = 0, slope = 1) + theme_minimal() + 
  theme(strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(), axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"), panel.grid = element_blank(),
        legend.position = "",panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  labs(x = "Observed FD with facilitation", y = "Observed FD without facilitation", color = "Treatment")

correlation_plots_sown_coexistence <- (omega_corr1 + theta_corr1)/(omega_corr + theta_corr) + plot_annotation(tag_levels = "A") & theme(plot.title = element_text(hjust = 0.5)) 
#ggsave("figures/FigS5.pdf", correlation_plots_sown_coexistence , width = 9, height =8, dpi = 300)
