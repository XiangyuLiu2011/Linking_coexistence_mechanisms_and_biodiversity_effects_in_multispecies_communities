##scripts/04_Linking coexistence_mechanisms_and_biodiversity_effects_without_facilitation
source("functions/useful SEM constrain functions.R")
source("functions/useful SEM partial plot functions.R")
library(lavaan)
library(corrplot)
library(ggpubr)
#make sure two datasets have the same names of Plot and harvest season
#make a copy of df_biodiversity.effects
df_biodiversity.effects_c <- df_biodiversity.effects
names(df_biodiversity.effects_c)[c(1,10)] <- c("plot", "season")

df_ini_coexistence_bio.effects <- merge(combined_sown_coex_new, df_biodiversity.effects_c[,c(1:2,8:15)], by=c("season", "plot", "Block"))
df_ini_coexistence_bio.effects$Month<-1*(df_ini_coexistence_bio.effects$season=="August")

##exclude SD==20, because all 20 species plots have the same predicted coexistence_mechanisms
df_ini_coexistence_bio.effects1<-subset(df_ini_coexistence_bio.effects, SD!=20)
head(df_ini_coexistence_bio.effects1)

Scalevars_both_seasons <-df_ini_coexistence_bio.effects1[,c("Month", "N", "F", "SD", "sownSLA", "sownMPD", "sta.Omega", "sta.theta", "sta.differential", "NE", "CE", "SE", "SE_inter", "SE_intra")]
# do the scaling
Scalevars_both_seasons2 <- scale(Scalevars_both_seasons)
corrs.2<-cor(Scalevars_both_seasons2, use = "pairwise.complete.obs")
corrplot(corrs.2, method = "square", type = "upper", order = "hclust",
  col = colorRampPalette(c("#336699", "white", "#990000"))(200),
  addCoef.col = "black", number.cex = 0.65, tl.col = "black",
  tl.cex = 0.8, tl.srt = 45, diag = FALSE, mar = c(0, 0, 1, 0) )

cat.vars.both.seasons <- df_ini_coexistence_bio.effects1[,c("plot", "Block", "composition", "season", "N")]
names(cat.vars.both.seasons)[5] <- c("Nitrogen")
cat.vars.both.seasons$group <- interaction(cat.vars.both.seasons$season, cat.vars.both.seasons$Nitrogen)
data.scaled.both.seasons <- cbind(cat.vars.both.seasons, Scalevars_both_seasons2)
head(data.scaled.both.seasons)
data.scaled.both.seasons$season<-as.factor(data.scaled.both.seasons$season)
data.scaled.both.seasons$group<-as.factor(data.scaled.both.seasons$group)

##SEM model including indirect interactions
model_string_both_seasons_N_indirect <- '
CE ~ F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
SE_inter ~  F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
SE_intra ~ F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
sta.Omega ~ F + SD + sownSLA + sownMPD + sta.differential
sta.theta ~ F + SD + sownSLA + sownMPD 
sta.differential ~ F + SD + sownSLA + sownMPD

sta.Omega~~sta.theta 
sta.theta~~sta.differential
sownSLA~~sownMPD
'

fituw_both_seasons_N_indirect <- sem(model_string_both_seasons_N_indirect, data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect, standardized =TRUE, rsquare=TRUE)
sort.pars(fituw_both_seasons_N_indirect, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect, data.scaled.both.seasons)
constrain_model_both_seasons_N_indirect$final_fit
best_string_both_seasons_N_indirect <- constrain_model_both_seasons_N_indirect$final_model_string

confituw_both_seasons_N_indirect <- sem(unlist(best_string_both_seasons_N_indirect), data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect, c("cfi","rmsea","srmr"))## cfi:0.99, rmsea:0.02 ,srmr0.07
sort.pars(confituw_both_seasons_N_indirect, response = "CE", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect, response = "sta.differential", group.names = unique(data.scaled.both.seasons$group))

#Constrain model mannually
#F on SEAS (dominance effects) August vs. June
mcs1_sown_N_indirect_1.0 <- constrain_path_cfa(confituw_both_seasons_N_indirect, "SE_inter ~ F", coef_name = "coef244Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_1.0 <- sem(mcs1_sown_N_indirect_1.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_1.0 <- constrain_path_cfa(fc1a_sown_N_indirect_1.0, "SE_inter ~ F", coef_name = "coef244Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_1.0 <- sem(mcs2_sown_N_indirect_1.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_1.0)
lavTestLRT(fc1_sown_N_indirect_1.0, confituw_both_seasons_N_indirect)#accepted
sort.pars(fc1_sown_N_indirect_1.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_1.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#Sown MPD on CE (complementarity effects) N vs. No N
mcs1_sown_N_indirect_2.0 <- constrain_path_cfa(fc1_sown_N_indirect_1.0, "CE ~ sownMPD", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_2.0 <- sem(mcs1_sown_N_indirect_2.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0 <- constrain_path_cfa(fc1a_sown_N_indirect_2.0, "CE ~ sownMPD", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_2.0 <- sem(mcs2_sown_N_indirect_2.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0)
lavTestLRT(fc1_sown_N_indirect_1.0, fc1_sown_N_indirect_2.0)#rejected

#F on CE (complementarity effects) August vs. June
mcs1_sown_N_indirect_2.0.1 <- constrain_path_cfa(fc1_sown_N_indirect_1.0, "CE ~ F", coef_name = "coef247Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_2.0.1 <- sem(mcs1_sown_N_indirect_2.0.1, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0.1 <- constrain_path_cfa(fc1a_sown_N_indirect_2.0.1, "CE ~ F", coef_name = "coef247Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_2.0.1 <- sem(mcs2_sown_N_indirect_2.0.1, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0.1)
lavTestLRT(fc1_sown_N_indirect_1.0, fc1_sown_N_indirect_2.0.1)#accepted
sort.pars(fc1_sown_N_indirect_2.0.1, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SE_intra N vs. noN
mcs1_sown_N_indirect_3.0 <- constrain_path_cfa(fc1_sown_N_indirect_2.0.1, "SE_intra ~ sownMPD", coef_name = "coef248N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0 <- sem(mcs1_sown_N_indirect_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0 <- constrain_path_cfa(fc1a_sown_N_indirect_3.0, "SE_intra ~ sownMPD", coef_name = "coef248noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0 <- sem(mcs2_sown_N_indirect_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0)
lavTestLRT(fc1_sown_N_indirect_2.0.1, fc1_sown_N_indirect_3.0)#rejected

#SD on omega N vs. noN
mcs1_sown_N_indirect_6.0 <- constrain_path_cfa(fc1_sown_N_indirect_2.0.1, "sta.Omega ~ SD", coef_name = "coef250Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_6.0 <- sem(mcs1_sown_N_indirect_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_6.0 <- constrain_path_cfa(fc1a_sown_N_indirect_6.0, "sta.Omega ~ SD", coef_name = "coef250Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_6.0 <- sem(mcs2_sown_N_indirect_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_6.0)
lavTestLRT(fc1_sown_N_indirect_2.0.1, fc1_sown_N_indirect_6.0)#accepted 
sort.pars(fc1_sown_N_indirect_6.0, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownSLA on theta N vs. noN
mcs1_sown_N_indirect_7.0 <- constrain_path_cfa(fc1_sown_N_indirect_6.0, "sta.theta ~ sownSLA", coef_name = "coef251N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_7.0 <- sem(mcs1_sown_N_indirect_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_7.0 <- constrain_path_cfa(fc1a_sown_N_indirect_7.0, "sta.theta ~ sownSLA", coef_name = "coef251noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_7.0 <- sem(mcs2_sown_N_indirect_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_7.0)
lavTestLRT(fc1_sown_N_indirect_7.0, fc1_sown_N_indirect_6.0)#accepted (p=0.08)
sort.pars(fc1_sown_N_indirect_7.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on theta June vs. August
mcs1_sown_N_indirect_8.0 <- constrain_path_cfa(fc1_sown_N_indirect_7.0, "sta.theta ~ sownMPD", coef_name = "coef252Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_8.0 <- sem(mcs1_sown_N_indirect_8.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_8.0 <- constrain_path_cfa(fc1a_sown_N_indirect_8.0, "sta.theta ~ sownMPD", coef_name = "coef252June", groups_constrain = c(3,4))
fc1_sown_N_indirect_8.0 <- sem(mcs2_sown_N_indirect_8.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_8.0)
lavTestLRT(fc1_sown_N_indirect_7.0, fc1_sown_N_indirect_8.0)#accepted
sort.pars(fc1_sown_N_indirect_8.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))
fitmeasures(fc1_sown_N_indirect_8.0, c("cfi","rmsea","srmr"))## cfi:0.99, rmsea:0.03 ,srmr0.08
summary(fc1_sown_N_indirect_8.0, standardized =TRUE, rsquare=TRUE)
sort.pars(fc1_sown_N_indirect_8.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SEAS (dominance effects) N vs. No N
mcs1_sown_N_indirect_11.0 <- constrain_path_cfa(fc1_sown_N_indirect_8.0, "SE_inter ~ sownMPD", coef_name = "coef255N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_11.0 <- sem(mcs1_sown_N_indirect_11.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_11.0 <- constrain_path_cfa(fc1a_sown_N_indirect_11.0, "SE_inter ~ sownMPD", coef_name = "coef255noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_11.0 <- sem(mcs2_sown_N_indirect_11.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_11.0)
lavTestLRT(fc1_sown_N_indirect_8.0, fc1_sown_N_indirect_11.0)#accepted
sort.pars(fc1_sown_N_indirect_11.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))

#sta.Omega on SEAS (dominance effects) N vs. No N
mcs1_sown_N_indirect_12.0 <- constrain_path_cfa(fc1_sown_N_indirect_11.0, "SE_inter ~ sta.Omega", coef_name = "coef256N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_12.0 <- sem(mcs1_sown_N_indirect_12.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_12.0 <- constrain_path_cfa(fc1a_sown_N_indirect_12.0, "SE_inter ~ sta.Omega", coef_name = "coef256noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_12.0 <- sem(mcs2_sown_N_indirect_12.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_12.0)
lavTestLRT(fc1_sown_N_indirect_12.0, fc1_sown_N_indirect_11.0)#accepted
sort.pars(fc1_sown_N_indirect_12.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))

##the best model is fc1_sown_N_indirect_12.0, fig3 was made based on this SEM in PowerPoint
fitMeasures(fc1_sown_N_indirect_12.0, c("cfi","rmsea","srmr"))## cfi:0.99, rmsea:0.04 ,srmr0.08
summary(fc1_sown_N_indirect_12.0, standardized =TRUE, rsquare=TRUE)
sort.pars(fc1_sown_N_indirect_12.0, response = "sta.Omega", sig=T, group.names = unique(data.scaled.both.seasons$group))

##make partial plots
final_model_fit_final <- fc1_sown_N_indirect_12.0
lavInspect(final_model_fit_final, "group.label")#the group order is: August N, August no N, June N and June no N
# split data by group (treatment)
Bio_list <- split(data.scaled.both.seasons, data.scaled.both.seasons$group)
new_names <- names(Bio_list)
new_names <- gsub("August\\.0", "August no N", new_names)
new_names <- gsub("June\\.0", "June no N", new_names)
new_names <- gsub("August\\.1", "August N", new_names)
new_names <- gsub("June\\.1", "June N", new_names)
names(Bio_list) <- new_names
Bio_list <- Bio_list[c("August N", "August no N", "June N", "June no N")]#reordered as final SEM

Bio_list_unscaled <- split(Scalevars_both_seasons, interaction(Scalevars_both_seasons$Month, Scalevars_both_seasons$N, sep = "_"))
names(Bio_list_unscaled) <- c("June no N", "August no N", "June N", "August N")
Bio_list_unscaled <- Bio_list_unscaled[c("August N", "August no N", "June N", "June no N")]
# --- partial residuals for your predictors ---
## Effects of SD on Complementarity (CE)
pr_SLA_theta <- part.res.sem(final_model_fit_final, c("sownSLA"), "sta.theta", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_Omega_CE <- part.res.sem(final_model_fit_final, c("sta.Omega"), "CE", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_theta_CE <- part.res.sem(final_model_fit_final, c("sta.theta"), "CE", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_indirect_CE <- part.res.sem(final_model_fit_final, c("sta.differential"), "CE", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_Omega_SEAS <- part.res.sem(final_model_fit_final, c("sta.Omega"), "SE_inter", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_theta_SEAS <- part.res.sem(final_model_fit_final, c("sta.theta"), "SE_inter", Bio_list, groups = 1:length(Bio_list), level = 1)
	
pr_indirect_SEAS <- part.res.sem(final_model_fit_final, c("sta.differential"), "SE_inter", Bio_list, groups = 1:length(Bio_list), level = 1)
	
pr_Omega_SEPS <- part.res.sem(final_model_fit_final, c("sta.Omega"), "SE_intra", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_theta_SEPS <- part.res.sem(final_model_fit_final, c("sta.theta"), "SE_intra", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_indirect_SEPS <- part.res.sem(final_model_fit_final, c("sta.differential"), "SE_intra", Bio_list, groups = 1:length(Bio_list), level = 1)

cols <- c("#336699", "#99CCFF", "#336699", "#99CCFF")
group_labels <- c("August N", "August no N", "June N", "June no N")

##fig4 in the paper
# Now create the PDF with layout and legend on the right
pdf("figures/fig4.pdf", width = 11, height = 9)
layout_matrix <- matrix(c(
  1,  2,  3,  10,
  4,  5,  6,  10,
  7,  8,  9,  10
), nrow = 3, byrow = TRUE)

layout(layout_matrix, widths = c(1, 1, 1, 0.55))
par(mar = c(4.5, 4.5, 2, 1))

# Define labels
panel_labels <- letters[1:9]

# Plot 1
plot.res.sem3_auto1(Bio_list_unscaled, "sta.Omega", residual_list = pr_Omega_CE, cols = cols,
                   alpha = 0.8, xlab = "", ylab = "Complementarity effects",
                   group_labels = group_labels, group_combinations = list(c(1, 2, 3, 4)), # triggers auto "combine all" mode
                   pch = 21, cex = 1.1, lwd = 0.5
                   )
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("A", side = 3, adj = -0.12, line = 0.2, font = 2)
# Plot 2
plot.res.sem3_1(Bio_list_unscaled, "sta.theta", residual_list = pr_theta_CE,
              cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
              xlab = "",  ylab = "", 
              group_combinations = list(c(1, 3),  c(2, 4)), show_regression = c(FALSE,FALSE))

abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("B", side = 3, adj = -0.12, line = 0.2, font = 2)
# Plot 3
plot.res.sem3_1(Bio_list_unscaled, "sta.differential", residual_list = pr_indirect_CE,
              cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
              xlab = "", 
              ylab = "",
              group_combinations = list(c(1, 3),  c(2, 4)), 
              show_regression = c(FALSE, FALSE))
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("C", side = 3, adj = -0.12, line = 0.2, font = 2)

# Plot 4
plot.res.sem3_1(Bio_list_unscaled, "sta.Omega", residual_list = pr_Omega_SEPS, cols = cols,
                   alpha = 0.8, xlab = "", ylab = "Intraspecific selection effects",
                   group_labels = group_labels, group_combinations = list(c(1, 3),  c(2, 4)), # triggers auto "combine all" mode
                   pch = 21, cex = 1.1, lwd = 0.5, show_regression = c(FALSE, FALSE)
                   )
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("D", side = 3, adj = -0.12, line = 0.2, font = 2)
# Plot 5
plot.res.sem3_1(Bio_list_unscaled, "sta.theta", residual_list = pr_theta_SEPS,
              cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
              xlab = "", 
              ylab = "",
              group_combinations = list(c(1, 3),  c(2, 4)), 
              show_regression = c(FALSE, FALSE))
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("E", side = 3, adj = -0.12, line = 0.2, font = 2)
# Plot 6
plot.res.sem3_1(Bio_list_unscaled, "sta.differential", residual_list = pr_indirect_SEPS, cols = cols,
                   alpha = 0.8, xlab = "", ylab = "",
                   group_labels = group_labels, group_combinations = list(c(1, 3),  c(2, 4)), # triggers auto "combine all" mode
                   pch = 21, cex = 1.1, lwd = 0.5, show_regression = c(FALSE, FALSE)
                   )
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("F", side = 3, adj = -0.12, line = 0.2, font = 2)
# Plot 7
plot.res.sem3_1(Bio_list_unscaled, "sta.Omega", residual_list = pr_Omega_SEAS,
              cols = cols, alpha = 0.8,
              xlab = "Predicted niche differences", 
              ylab = "Interspecific selection effects", 
              group_labels = group_labels,
              group_combinations = list(c(1, 3),  c(2, 4)), show_regression = c(TRUE, TRUE),
              pch = 21, cex = 1.1, lwd = 0.5)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("G", side = 3, adj = -0.12, line = 0.2, font = 2)


# Plot 8
plot.res.sem3_auto1(Bio_list_unscaled, "sta.theta", residual_list = pr_theta_SEAS, cols = cols,
                   alpha = 0.8, xlab = "Predicted competitive ability differences", ylab = "",
                   group_labels = group_labels, group_combinations = list(c(1, 2, 3, 4)),
                   pch = 21, cex = 1.1, lwd = 0.5) 
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("H", side = 3, adj = -0.12, line = 0.2, font = 2)


# Plot 9
plot.res.sem3_1(Bio_list_unscaled, "sta.differential", residual_list = pr_indirect_SEAS,
              cols = cols, alpha = 0.8,
              xlab = "Predicted indirect interactions", 
              ylab = "",
              group_combinations = list(c(1, 3),  c(2, 4)), 
              show_regression = c(FALSE, FALSE), pch = 21, cex = 1.1, lwd = 0.5)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("I", side = 3, adj = -0.12, line = 0.2, font = 2)
# Legend in the right column
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")
legend("center", pch = 19, lty = 1, col = c("#99CCFF", "#336699"),
       legend = c("Without nitrogen", "With nitrogen"),
       bty = "n", pt.cex = 1.2, cex = 0.95, , xpd = NA)
dev.off()

#######
pr_SLA_theta <- part.res.sem(final_model_fit_final, c("sownSLA"), "sta.theta", Bio_list,
                            groups = 1:length(Bio_list), level = 1)

pr_MPD_omega <- part.res.sem(final_model_fit_final, c("sownSLA"), "sta.Omega", Bio_list,
                             groups = 1:length(Bio_list), level = 1)

## ---- FigS7 ---------------
pdf("figures/FigS7.pdf", width = 9, height = 4)
# 1 row: A | B | legend
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.38))
# Global settings for the two panels
par(oma = c(0.5, 0.5, 0.2, 0.2))  # outer margins (small)
par(mar = c(4.2, 4.4, 2.2, 0.8))  # inner margins for panels

# ---- Panel A ----
plot.res.sem3_auto1(
  Bio_list_unscaled, "sownMPD", residual_list = pr_MPD_omega,
  cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
  xlab = "sownMPD",
  ylab = "Predicted niche differences",
  group_labels = group_labels,
  group_combinations = list(c(1, 2, 3, 4))
)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("A", side = 3, adj = -0.08, line = 0.2, font = 2)

# ---- Panel B ----
plot.res.sem3_1(
  Bio_list_unscaled, "sownSLA", residual_list = pr_SLA_theta,
  cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
  xlab = "sownSLA",
  ylab = "Predicted fitness differences",
  group_labels = group_labels,
  group_combinations = list(c(1, 3),  c(2, 4)),
  show_regression = c(TRUE, TRUE)
)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("B", side = 3, adj = -0.08, line = 0.2, font = 2)

# ---- Legend column ----
par(mar = c(0, 0, 0, 0))
plot.new()
par(xpd = NA)
legend(
  "center",
  pch = 19, lty = 1,
  col = cols[c(1, 3)],
  legend = c("Without nitrogen", "With nitrogen"),
  bty = "n",
  pt.cex = 1.4, cex = 1
)
dev.off()

######
##Fig.S3
Fig.indirect.fitness.diff <- ggplot(df_ini_coexistence_bio.effects1,
       aes(x = sta.differential, y = sta.theta, color = factor(N))) +
  geom_point(size = 2, shape = 16, alpha = 0.8,
             position = position_dodge(width = 0.5)) +
  #geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2, color = "black") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", show.legend = FALSE, color = "black") + 
  #facet_wrap(~ season, scales = "free_y", nrow = 1) +
  scale_color_manual(
    values = c("0" = "#99CCFF", "1" = "#336699"),
    breaks = c("0", "1"),
    labels = c("Without nitrogen", "With nitrogen"),
    name = ""  # 
  ) + theme_minimal() +
  theme(
    strip.background = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    panel.grid = element_blank(),
    legend.position = "right",   
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(x = "Indirect interactions", y = "Predicted competitive ability differences")
ggsave("figures/FigS3.pdf", Fig.indirect.fitness.diff , width = 5, height = 4, dpi = 300)
