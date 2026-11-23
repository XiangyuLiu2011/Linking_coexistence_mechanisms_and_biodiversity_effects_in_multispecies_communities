##scripts/06_Linking_coexistence_mechanisms_with_biodiversity_effects_including_facilitation

source("functions/useful SEM constrain functions.R")
source("functions/useful SEM partial plot functions.R")
library(lavaan)
library(corrplot)
library(ggpubr)
library(broom)
library(patchwork)
library(dplyr)
#merge two datasets
df_ini_coexistence_bio.effects_fac <- merge(combined_sown_coex_fac, df_biodiversity.effects_c[,c(1:2,8:17)], by=c("season", "plot", "Block"))
df_ini_coexistence_bio.effects_fac$Month<-1*(df_ini_coexistence_bio.effects_fac$season=="August")

##exclude SD==20, because all 20 species plots have the same predicted coexistence_mechanisms
df_ini_coexistence_bio.effects1_fac<-subset(df_ini_coexistence_bio.effects_fac, SD!=20)
head(df_ini_coexistence_bio.effects1_fac)

Scalevars_both_seasons_fac <-df_ini_coexistence_bio.effects1_fac[,c("Month", "N", "F", "SD", "sownSLA", "sownMPD", "sta.Omega", "sta.theta", "sta.differential", "NE", "CE_sqrt", "SE", "SE_inter", "SE_intra_sqrt")]
names(Scalevars_both_seasons_fac)[c(11,14)]<-c("CE", "SE_intra")
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
##
fituw_both_seasons_N_indirect_fac_1st <- sem(model_string_both_seasons_N_indirect_1st, data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect_fac_1st, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect_fac_1st, standardized =TRUE, rsquare=TRUE)
##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect_fac_1st <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect_fac_1st, data.scaled.both.seasons_fac)
constrain_model_both_seasons_N_indirect_fac_1st$final_fit
best_string_both_seasons_N_indirect_fac_1st <- constrain_model_both_seasons_N_indirect_fac_1st$final_model_string

confituw_both_seasons_N_indirect_fac_1st <- sem(unlist(best_string_both_seasons_N_indirect_fac_1st), data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect_fac_1st, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect_fac_1st, c("pvalue","cfi","rmsea","srmr"))##p:0.074, cfi:0.959, rmsea:0.061 ,srmr0.074
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "SE_intra", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac_1st, response = "sta.differential", group.names = unique(data.scaled.both.seasons_fac$group))

#Constrain model mannually
#sownMPD on sta.differential August vs. June
mcs1_sown_N_indirect_1.0_fac_1st <- constrain_path_cfa(confituw_both_seasons_N_indirect_fac_1st, "sta.differential ~ sownMPD", coef_name = "coef244Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_1.0_fac_1st <- sem(mcs1_sown_N_indirect_1.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_1.0_fac_1st <- constrain_path_cfa(fc1a_sown_N_indirect_1.0_fac_1st, "sta.differential ~ sownMPD", coef_name = "coef244Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_1.0_fac_1st <- sem(mcs2_sown_N_indirect_1.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_1.0_fac_1st)
lavTestLRT(fc1_sown_N_indirect_1.0_fac_1st, confituw_both_seasons_N_indirect_fac_1st)#accepted
sort.pars(fc1_sown_N_indirect_1.0_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))

#F on SEAS (dominance effects) June vs. August
mcs1_sown_N_indirect_2.0_fac_1st <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac_1st, "SE_inter ~ F", coef_name = "coef245Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_2.0_fac_1st <- sem(mcs1_sown_N_indirect_2.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0_fac_1st <- constrain_path_cfa(fc1a_sown_N_indirect_2.0_fac_1st, "SE_inter ~ F", coef_name = "coef245Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_2.0_fac_1st <- sem(mcs2_sown_N_indirect_2.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0_fac_1st)
lavTestLRT(fc1_sown_N_indirect_1.0_fac_1st, fc1_sown_N_indirect_2.0_fac_1st)#accepted
sort.pars(fc1_sown_N_indirect_2.0_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_2.0_fac_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))

#sownSLA on sta.Omega N vs. noN
mcs1_sown_N_indirect_3.0_fac_1st <- constrain_path_cfa(fc1_sown_N_indirect_2.0_fac_1st, "sta.Omega ~ sownSLA", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0_fac_1st <- sem(mcs1_sown_N_indirect_3.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0_fac_1st <- constrain_path_cfa(fc1a_sown_N_indirect_3.0_fac_1st, "sta.Omega ~ sownSLA", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0_fac_1st <- sem(mcs2_sown_N_indirect_3.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0_fac_1st)
lavTestLRT(fc1_sown_N_indirect_3.0_fac_1st, fc1_sown_N_indirect_2.0_fac_1st)#accepted
sort.pars(fc1_sown_N_indirect_3.0_fac_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))

#sownSLA on theta N vs. noN
mcs1_sown_N_indirect_4.0_fac_1st <- constrain_path_cfa(fc1_sown_N_indirect_3.0_fac_1st, "sta.theta ~ sownSLA", coef_name = "coef247N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_4.0_fac_1st <- sem(mcs1_sown_N_indirect_4.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_4.0_fac_1st <- constrain_path_cfa(fc1a_sown_N_indirect_4.0_fac_1st, "sta.theta ~ sownSLA", coef_name = "coef247noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_4.0_fac_1st <- sem(mcs2_sown_N_indirect_4.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_4.0_fac_1st)
lavTestLRT(fc1_sown_N_indirect_4.0_fac_1st, fc1_sown_N_indirect_3.0_fac_1st)#accepted
sort.pars(fc1_sown_N_indirect_4.0_fac_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_4.0_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SE_inter N vs. noN
mcs1_sown_N_indirect_5.0_fac_1st <- constrain_path_cfa(fc1_sown_N_indirect_4.0_fac_1st, "SE_inter ~ sownMPD", coef_name = "coef248N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_5.0_fac_1st <- sem(mcs1_sown_N_indirect_5.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_5.0_fac_1st <- constrain_path_cfa(fc1a_sown_N_indirect_5.0_fac_1st, "SE_inter ~ sownMPD", coef_name = "coef248noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_5.0_fac_1st <- sem(mcs2_sown_N_indirect_5.0_fac_1st, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_5.0_fac_1st)
lavTestLRT(fc1_sown_N_indirect_5.0_fac_1st, fc1_sown_N_indirect_4.0_fac_1st)#accepted
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))

##fc1_sown_N_indirect_5.0_fac_1st, is the best model, fig.S4 made in PowerPoint based on the best SEM model
summary(fc1_sown_N_indirect_5.0_fac_1st, standardized =TRUE, rsquare=TRUE)
fitMeasures(fc1_sown_N_indirect_5.0_fac_1st, c("pvalue","cfi","rmsea","srmr"))##p:0.054, cfi:0.952, rmsea:0.063, srmr 0.080
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "sta.differential", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_5.0_fac_1st, response = "SE_intra", group.names = unique(data.scaled.both.seasons_fac$group))

##-These code was not used anymore (could skip this part)------------------------------------
##To revise reviewer's comment on the circularity between sta.Omega and sta.differential: 
##sta.Omega~sta.differential was removed, and sta.Omega~~sta.differential was added.
fituw_both_seasons_N_indirect_fac <- sem(model_string_both_seasons_N_indirect, data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect_fac, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect_fac, standardized =TRUE, rsquare=TRUE)
##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect_fac <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect_fac, data.scaled.both.seasons_fac)
constrain_model_both_seasons_N_indirect_fac$final_fit
best_string_both_seasons_N_indirect_fac <- constrain_model_both_seasons_N_indirect_fac$final_model_string

confituw_both_seasons_N_indirect_fac <- sem(unlist(best_string_both_seasons_N_indirect_fac), data = data.scaled.both.seasons_fac,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect_fac, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect_fac, c("pvalue","cfi","rmsea","srmr"))##P:0.052, cfi:0.954, rmsea:0.065 ,srmr0.074
sort.pars(confituw_both_seasons_N_indirect_fac, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "SE_intra", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.theta", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(confituw_both_seasons_N_indirect_fac, response = "sta.differential", group.names = unique(data.scaled.both.seasons_fac$group))

#Constrain model mannually
#sownMPD on sta.differential August vs. June
mcs1_sown_N_indirect_1.0_fac <- constrain_path_cfa(confituw_both_seasons_N_indirect_fac, "sta.differential ~ sownMPD", coef_name = "coef244Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_1.0_fac <- sem(mcs1_sown_N_indirect_1.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_1.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_1.0_fac, "sta.differential ~ sownMPD", coef_name = "coef244Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_1.0_fac <- sem(mcs2_sown_N_indirect_1.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_1.0_fac)
lavTestLRT(fc1_sown_N_indirect_1.0_fac, confituw_both_seasons_N_indirect_fac)#accepted
sort.pars(fc1_sown_N_indirect_1.0_fac, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))

#sta.Omega on SEAS (dominance effects) August vs. June
mcs1_sown_N_indirect_2.0.1_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "SE_inter ~ sta.Omega", coef_name = "coef245Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_2.0.1_fac <- sem(mcs1_sown_N_indirect_2.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0.1_fac <- constrain_path_cfa(fc1a_sown_N_indirect_2.0.1_fac, "SE_inter ~ sta.Omega", coef_name = "coef245Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_2.0.1_fac <- sem(mcs2_sown_N_indirect_2.0.1_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0.1_fac)
lavTestLRT(fc1_sown_N_indirect_1.0_fac, fc1_sown_N_indirect_2.0.1_fac)#rejected!

#sownMPD on CE (complementarity effects) N vs. NoN
mcs1_sown_N_indirect_3.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "CE ~ sownMPD", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0_fac <- sem(mcs1_sown_N_indirect_3.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_3.0_fac, "CE ~ sownMPD", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0_fac <- sem(mcs2_sown_N_indirect_3.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0_fac)
lavTestLRT(fc1_sown_N_indirect_3.0_fac, fc1_sown_N_indirect_1.0_fac)#rejected

#sownSLA on theta N vs. noN
mcs1_sown_N_indirect_4.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_1.0_fac, "sta.theta ~ sownSLA", coef_name = "coef247N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_4.0_fac <- sem(mcs1_sown_N_indirect_4.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_4.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_4.0_fac, "sta.theta ~ sownSLA", coef_name = "coef247noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_4.0_fac <- sem(mcs2_sown_N_indirect_4.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_4.0_fac)
lavTestLRT(fc1_sown_N_indirect_4.0_fac, fc1_sown_N_indirect_1.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_4.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownSLA on omega N vs. noN
mcs1_sown_N_indirect_5.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_4.0_fac, "sta.Omega ~ sownSLA", coef_name = "coef248N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_5.0_fac <- sem(mcs1_sown_N_indirect_5.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_5.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_5.0_fac, "sta.Omega ~ sownSLA", coef_name = "coef248noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_5.0_fac <- sem(mcs2_sown_N_indirect_5.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_5.0_fac)
lavTestLRT(fc1_sown_N_indirect_4.0_fac, fc1_sown_N_indirect_5.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_5.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#F on SE_inter June vs. August
mcs1_sown_N_indirect_6.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_5.0_fac, "SE_inter ~ F", coef_name = "coef249Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_6.0_fac <- sem(mcs1_sown_N_indirect_6.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_6.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_6.0_fac, "SE_inter ~ F", coef_name = "coef249Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_6.0_fac <- sem(mcs2_sown_N_indirect_6.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_6.0_fac)
lavTestLRT(fc1_sown_N_indirect_5.0_fac, fc1_sown_N_indirect_6.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_6.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SE_inter N vs. NoN
mcs1_sown_N_indirect_7.0_fac <- constrain_path_cfa(fc1_sown_N_indirect_6.0_fac, "SE_inter ~ sownMPD", coef_name = "coef250N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_7.0_fac <- sem(mcs1_sown_N_indirect_7.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_7.0_fac <- constrain_path_cfa(fc1a_sown_N_indirect_7.0_fac, "SE_inter ~ sownMPD", coef_name = "coef250NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_7.0_fac <- sem(mcs2_sown_N_indirect_7.0_fac, data = data.scaled.both.seasons_fac, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_7.0_fac)
lavTestLRT(fc1_sown_N_indirect_7.0_fac, fc1_sown_N_indirect_6.0_fac)#accepted
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

##fc1_sown_N_indirect_7.0_fac, is the best model, fig.S4 made in PowerPoint based on the best SEM model
summary(fc1_sown_N_indirect_7.0_fac, standardized =TRUE, rsquare=TRUE)
fitMeasures(fc1_sown_N_indirect_7.0_fac, c("pvalue","cfi","rmsea","srmr"))##p:0.028, cfi:0.941, rmsea:0.070, srmr 0.082
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "sta.Omega", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "sta.theta", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "sta.differential", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "CE", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "SE_inter", group.names = unique(data.scaled.both.seasons_fac$group))
sort.pars(fc1_sown_N_indirect_7.0_fac, response = "SE_intra", group.names = unique(data.scaled.both.seasons_fac$group))

#--------------------------------------------------
##Plot the direct vs. indirect effects through niche and competitive ability differences
##when facilitation was included
##extract the standardized effect sizes from final model fit
library(lavaan)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(manymome)
library(parallel)
library(tibble)

final_model_fit_final_fac <- fc1_sown_N_indirect_5.0_fac_1st
group_names_fac <- lavInspect(final_model_fit_final_fac, "group.label")
data.frame(group_number = seq_along(group_names_fac), group_label = group_names_fac)

# Inspect how uncertainty was calculated in the original model.
lavInspect(final_model_fit_final_fac, "options")[c("estimator", "se", "test")]
stopifnot(lavInspect(final_model_fit_final_fac, "converged"))

# Monte Carlo Estimates for 'indirect_effects' and 'cond_indirect_effects'
mc_results_fac <- do_mc(fit = final_model_fit_final_fac, R = 10000, seed = 123)

predictors <- c("F", "SD", "sownSLA", "sownMPD")
responses  <- c("CE", "SE_inter", "SE_intra")
mediators  <- c("sta.Omega", "sta.theta", "sta.differential")

get_effect_fac <- function(x, y, g, m = NULL) {
  manymome::indirect_effect(
    x = x,
    y = y,
    m = m,
    fit = final_model_fit_final_fac,
    group = g,
    standardized_x = FALSE,#use F keeps the constrained paths across groups
    standardized_y = FALSE,
    mc_ci = TRUE,
    mc_out = mc_results_fac,
    level = 0.95,
    parallel = FALSE,
    progress = FALSE
  )
}

effect_row_fac <- function(object, g, x, y, type) {
  
  ci <- as.numeric(confint(object))
  stopifnot(length(ci) == 2L)
  tibble::tibble(
    group = group_names_fac[g],
    predictor = x,
    response = y,
    effect_type = type,
    estimate = as.numeric(coef(object)),
    conf.low = ci[1],
    conf.high = ci[2]
  )
}

# One task per group × predictor × response.
tasks_fac <- expand.grid(
  g = seq_along(group_names_fac),
  x = c(predictors, mediators),
  y = responses,
  stringsAsFactors = FALSE
)

calculate_task_fac <- function(i) {
  
  g <- tasks_fac$g[i]
  x <- tasks_fac$x[i]
  y <- tasks_fac$y[i]
  
  direct <- get_effect_fac(x, y, g)
  
  # Coexistence metrics have only direct effects in this SEM.
  if (x %in% mediators) {
    return(effect_row_fac(direct, g, x, y, "Direct"))
  }
  
  via_ND <- get_effect_fac(x, y, g, m = "sta.Omega")
  via_FD <- get_effect_fac(x, y, g, m = "sta.theta")
  via_D  <- get_effect_fac(x, y, g, m = "sta.differential")
  
  # Preserve the joint-draw calculation for uncertainty in the sum.
  total_indirect <- via_ND + via_FD + via_D
  
  dplyr::bind_rows(
    effect_row_fac(direct, g, x, y, "Direct"),
    effect_row_fac(via_ND, g, x, y, "Via niche differences"),
    effect_row_fac(via_FD, g, x, y, "Via competitive ability differences"),
    effect_row_fac(via_D, g, x, y, "Via community-pair differential"),
    effect_row_fac(total_indirect, g, x, y, "Total indirect")
  )
}

# The wrapper ensures workers are stopped on completion or error.
run_parallel_effects_fac <- function(n_workers = 20L) {
  
  cl <- parallel::makePSOCKcluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(manymome)
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "final_model_fit_final_fac",
      "group_names_fac",
      "mc_results_fac",
      "mediators",
      "tasks_fac",
      "get_effect_fac",
      "effect_row_fac",
      "calculate_task_fac"
    ),
    envir = environment()
  )
  
  results <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(tasks_fac)),
    fun = calculate_task_fac
  )
  
  dplyr::bind_rows(results)
}

# Calculate effects for the model including facilitation.
effects_all_fac <- run_parallel_effects_fac(n_workers = 20L) %>%
  dplyr::mutate(
    group = factor(group, levels = group_names_fac),
    supported = conf.low > 0 | conf.high < 0
  )
View(effects_all_fac)#Unstandardized effect sizes

effects_all_fac %>% dplyr::filter(predictor == "sta.Omega", response == "CE") %>% dplyr::select(group, estimate, conf.low, conf.high)#ns
effects_all_fac %>% dplyr::filter(predictor == "sta.Omega", response == "SE_inter") %>% dplyr::select(group, estimate, conf.low, conf.high)#ns
effects_all_fac %>% dplyr::filter(predictor == "sta.Omega", response == "SE_intra") %>% dplyr::select(group, estimate, conf.low, conf.high)#ns

table(effects_all_fac$response)
table(effects_all_fac$effect_type)
table(effects_all_fac$predictor)

##Make a combined figure with and without facilitation
##“selected_paths” is the results from SEM without facilitation
path_key <- selected_paths %>% dplyr::select(response, predictor, effect_type, Pathway, path_id)
dim(path_key)#due to the space in main text, selected paths are shown
# Combine results from both models.
comparison_data <- bind_rows(
  effects_all %>%
    mutate(
      group = as.character(group),
      Facilitation = "Without facilitation"
    ),
  
  effects_all_fac %>%
    mutate(
      group = as.character(group),
      Facilitation = "With facilitation"
    )
) %>%
  inner_join(
    path_key,
    by = c("response", "predictor", "effect_type")
  ) %>%
  mutate(
    Season = case_when(
      grepl("^June\\.", group)   ~ "June",
      grepl("^August\\.", group) ~ "August"
    ),
    
    Nitrogen = case_when(
      grepl("\\.0$", group) ~ "Without nitrogen",
      grepl("\\.1$", group) ~ "With nitrogen"
    ),
    
    Season = factor(Season, levels = c("June", "August")),
    
    Facilitation = factor(
      Facilitation,
      levels = c("Without facilitation", "With facilitation")
    ),
    
    Response = factor(
      response,
      levels = c("CE", "SE_intra", "SE_inter"),
      labels = c(
        "Complementarity effects",
        "Intraspecific selection effects",
        "Interspecific selection effects"
      )
    ),
    
    supported = conf.low > 0 | conf.high < 0,
    
    Significance = factor(
      if_else(supported, "Significant", "Nonsignificant"),
      levels = c("Significant", "Nonsignificant")
    ),
    
    # Requested pathway order, from top to bottom.
    y_base = nrow(path_key) - match(path_id, path_key$path_id) + 1,
    
    # Separate the two model estimates within each pathway.
    y_pos = y_base + if_else(
      Facilitation == "Without facilitation",
      0.25, -0.25
    )
  )

stopifnot(
  !anyNA(comparison_data$Season),
  !anyNA(comparison_data$Nitrogen),
  !anyNA(comparison_data$y_pos)
)

axis_key <- path_key %>% mutate(y_base = nrow(path_key) - row_number() + 1)

# Limits based on the actual confidence-interval range
effect_range <- range(c(0, comparison_data$estimate,
    comparison_data$conf.low, comparison_data$conf.high), na.rm = TRUE)
padding <- max(0.10 * diff(effect_range), 0.01)
common_x_limits <- effect_range + c(-padding, padding)
facilitation_colours <- c("Without facilitation" = "#DDCC77","With facilitation" = "#8E6BB0")

make_nitrogen_panel <- function(nitrogen_level) {
  
  panel_data <- comparison_data %>%
    filter(Nitrogen == nitrogen_level) %>%
    mutate(
      point_fill = if_else(
        supported,
        unname(facilitation_colours[as.character(Facilitation)]),
        "white"
      ))
  
  ggplot(panel_data, aes(x = estimate, y = y_pos)  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.4
    ) +
     geom_segment(
      aes(
        x = conf.low,
        xend = conf.high,
        color = Facilitation
      ),
      linewidth = 0.75,
      alpha = 0.9
    ) +
    
    geom_point(
      aes(colour = Facilitation,fill = point_fill
      ),
      shape = 21,
      size = 2.2,
      stroke =1
    ) +
    
    facet_grid(
      Response ~ Season,
      scales = "free_y",
      space = "free_y",labeller = labeller(
        Response = c(
          "Complementarity effects" = "Complementarity",
          "Intraspecific selection effects" = "Intraspecific selection",
          "Interspecific selection effects" = "Interspecific selection"
        )
      )
    ) +
    
    scale_colour_manual(
      values = facilitation_colours,
      breaks = c("Without facilitation", "With facilitation"),
      name = NULL
    ) +
    
    scale_fill_identity() +
    scale_y_continuous(
      breaks = axis_key$y_base,
      labels = axis_key$Pathway,
      expand = expansion(add = 0.45)
    ) +
    
    coord_cartesian(clip = "on") +
    
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 21,
          fill = unname(facilitation_colours)
        )
      )
    ) +
    labs(
      title = nitrogen_level,
      x = "Direct and indirect effect sizes",
      y = NULL
    ) + theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, colour = "black"),
      strip.text.y.right = element_text(angle = 270, size = 7.5,lineheight = 1.1, colour = "black"),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      legend.position = "bottom",
      panel.spacing.x = grid::unit(0.5, "lines"),
      panel.spacing.y = grid::unit(0.3, "lines")
    )
}

# C: without nitrogen, June left and August right.
panel_C <- make_nitrogen_panel("Without nitrogen") +labs(x = NULL)
# D: with nitrogen, June left and August right.
panel_D <- make_nitrogen_panel("With nitrogen")

Fig3_CD <- ((panel_C / panel_D) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = list(c("C", "D")))
) & theme(legend.position = "bottom",legend.box = "horizontal",
  plot.tag = element_text(size = 11, face = "plain")
)

Fig3_CD
#ggsave(filename = "figures/Fig3_CD.pdf",plot = Fig3_CD, width = 6, height = 12, units = "in")

#-----------------------------------------------------------------
##Based on Reviewer 2, we tested these relationships by orthogonal regressions
##Fig.S8
df_coex_fac <- combined_sown_coex_fac[,c(1:2,5,7:8,15:17,19,33:35)]
names(df_coex_fac)[c(2:5,7:12)] <- c("Omega_fac", "theta_fac", "differential_fac", "overlap_fac","mean_omega_rand_fac", "mean_theta_rand_fac", "mean_differential_rand_fac",
                                    "sta.Omega_fac", "sta.theta_fac", "sta.differential_fac")
two_coexistence_dfs <- merge(combined_sown_coex_new, df_coex_fac, by=c("season", "plot"))
two_coexistence_dfs <- subset(two_coexistence_dfs, SD != 20)

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

library(smatr)
##Panel Fig.S5a: observed Omega
Obs_omega_MA <- sma(Omega ~ Omega_fac, method = "MA", data = two_coexistence_dfs)
Obs_omega_coef <- summary(Obs_omega_MA)#Omega =0.291*Omega_fac + 0.0039

Mean_omega_rand_MA <- sma(mean_omega_rand ~ mean_omega_rand_fac, method = "MA", data = two_coexistence_dfs)
Mean_omega_rand_coef <- summary(Mean_omega_rand_MA)#Mean expected Omega =0.384*Omega_fac - 0.0008

Sta_omega_MA <- sma(sta.Omega ~ sta.Omega_fac, method = "MA", data = two_coexistence_dfs)
Sta_omega_coef <- summary(Sta_omega_MA)##sta.Omega =1.839*sta.Omega_fac + 0.0605

##same analysis for structural competitive ability differences
Sta_theta_MA <- sma(sta.theta ~ sta.theta_fac, method = "MA", data = two_coexistence_dfs)
Sta_theta_coef <- summary(Sta_theta_MA)#sta.theta =1.063*sta.theta_fac + 0.0538

Obs_theta_MA <- sma(theta ~ theta_fac, method = "MA", data = two_coexistence_dfs)
Obs_theta_coef <- summary(Obs_theta_MA)#Obs.theta =0.925*Obs.theta_fac + 2.673

# Matching x and y limits within each panel
obs_omega_limits <- range(c(two_coexistence_dfs$Omega_fac, two_coexistence_dfs$Omega), finite = TRUE)
obs_theta_limits <- range(c(two_coexistence_dfs$theta_fac, two_coexistence_dfs$theta), finite = TRUE)
sta_omega_limits <- range(c(two_coexistence_dfs$sta.Omega_fac, two_coexistence_dfs$sta.Omega), finite = TRUE)
sta_theta_limits <- range(c(two_coexistence_dfs$sta.theta_fac, two_coexistence_dfs$sta.theta), finite = TRUE)

# Equations generated directly from the fitted coefficients
# %+ formats positive and negative intercepts appropriately
obs_omega_equation <- sprintf("y = %.3f x %+.4f", Obs_omega_coef$Slope[1], Obs_omega_coef$Int[1])
obs_theta_equation <- sprintf("y = %.3f x %+.3f", Obs_theta_coef$Slope[1], Obs_theta_coef$Int[1])

sta_omega_equation <- sprintf("y = %.3f x %+.3f", Sta_omega_coef$Slope[1], Sta_omega_coef$Int[1])
sta_theta_equation <- sprintf("y = %.3f x %+.3f", Sta_theta_coef$Slope[1], Sta_theta_coef$Int[1])

# Panel A: observed Omega (structural niche differences)
obs_omega_corr <- ggplot(two_coexistence_dfs, aes(x = Omega_fac, y = Omega, colour = sownSLA) ) + geom_point(
  position = position_jitter(width = 0.01, height = 0.01, seed = 123),
  size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40", linetype = "dashed", linewidth = 0.7
  ) +
  geom_abline(intercept = Obs_omega_coef$Int[1], slope = Obs_omega_coef$Slope[1], colour = "black", linewidth = 0.9) +
  scale_colour_gradientn(
    colours = c("darkgreen", "green", "yellow", "pink", "white"),
    limits = c(15, 25),
    name = "Sown CWM SLA"
  ) +
  labs(
    x = "Observed niche differences\nwith facilitation",
    y = "Observed niche differences\nwithout facilitation"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )
obs_omega_corr

obs_omega_corr <- obs_omega_corr + coord_equal(xlim = obs_omega_limits, ylim = obs_omega_limits) +
  annotate("text", x = -Inf, y = Inf, label = obs_omega_equation,
    hjust = -0.08, vjust = 1.5, size = 3.5)

# Panel B: observed theta (structural competitive ability differences)
obs_theta_corr <- ggplot(two_coexistence_dfs, aes(x = theta_fac, y = theta, colour = sownSLA) ) + geom_point(
  position = position_jitter(width = 0.01, height = 0.01, seed = 123),
  size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40", linetype = "dashed", linewidth = 0.7
  ) + geom_abline(intercept = Obs_theta_coef$Int[1], slope = Obs_theta_coef$Slope[1], colour = "black", linewidth = 0.9) +
  scale_colour_gradientn(
    colours = c("darkgreen", "green", "yellow", "pink", "white"),
    limits = c(15, 25),
    name = "Sown CWM SLA"
  ) +
  labs(
    x = "Observed competitive ability differences\nwith facilitation",
    y = "Observed competitive ability differences\nwithout facilitation"
  ) +
  coord_equal() + theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )
obs_theta_corr

obs_theta_corr <- obs_theta_corr + coord_equal(xlim = obs_theta_limits, ylim = obs_theta_limits) +
  annotate("text", x = -Inf, y = Inf, label = obs_theta_equation,
    hjust = -0.08, vjust = 1.5, size = 3.5)

# Panel C: standardized Omega
sta_omega_corr <- ggplot(two_coexistence_dfs, aes(x = sta.Omega_fac, y = sta.Omega, colour = sownSLA) ) + geom_point(
  position = position_jitter(width = 0.01, height = 0.01, seed = 123),
  size = 3, alpha = 0.8
) +
  geom_abline(
    intercept = 0, slope = 1,
    colour = "grey40", linetype = "dashed", linewidth = 0.7
  ) +
  geom_abline(intercept = Sta_omega_coef$Int[1], slope = Sta_omega_coef$Slope[1], colour = "black", linewidth = 0.9) +
  scale_colour_gradientn(
    colours = c("darkgreen", "green", "yellow", "pink", "white"),
    limits = c(15, 25),
    name = "Sown CWM SLA"
  ) +
  labs(
    x = "Standardized niche differences\nwith facilitation",
    y = "Standardized niche differences\nwithout facilitation"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )
sta_omega_corr

sta_omega_corr <- sta_omega_corr + coord_equal(xlim = sta_omega_limits, ylim = sta_omega_limits) +
  annotate("text", x = -Inf, y = Inf, label = sta_omega_equation,
    hjust = -0.08, vjust = 1.5, size = 3.5)

# Panel D: standardized theta
sta_theta_corr <- ggplot(two_coexistence_dfs, aes(x = sta.theta_fac, y = sta.theta, colour = sownSLA) ) + geom_point(
  position = position_jitter(width = 0.01, height = 0.01, seed = 123),
  size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40", linetype = "dashed", linewidth = 0.7
  ) + geom_abline(intercept = Sta_theta_coef$Int[1], slope = Sta_theta_coef$Slope[1], colour = "black", linewidth = 0.9) +
  scale_colour_gradientn(
    colours = c("darkgreen", "green", "yellow", "pink", "white"),
    limits = c(15, 25),
    name = "Sown CWM SLA"
  ) +
  labs(
    x = "Standardized competitive ability differences\nwith facilitation",
    y = "Standardized competitive ability differences\nwithout facilitation"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )

Sta_differential_MA <- sma(sta.differential ~ sta.differential_fac, method = "MA", data = two_coexistence_dfs)
Sta_differential_coef <- summary(Sta_differential_MA)#sta.theta =1.063*sta.theta_fac + 0.0538

ggplot(two_coexistence_dfs, aes(x = sta.differential_fac, y = sta.differential, colour = sownSLA) ) + geom_point(
  position = position_jitter(width = 0.01, height = 0.01, seed = 123),
  size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40", linetype = "dashed", linewidth = 0.7
  ) + geom_abline(intercept = Sta_differential_coef$Int[1], slope = Sta_differential_coef$Slope[1], colour = "black", linewidth = 0.9) +
  scale_colour_gradientn(
    colours = c("darkgreen", "green", "yellow", "pink", "white"),
    limits = c(15, 25),
    name = "Sown CWM SLA"
  ) +
  labs(
    x = "Standardized community-pair differential\nwith facilitation",
    y = "Standardized community-pair differential\nwithout facilitation"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )

sta_theta_corr <- sta_theta_corr + coord_equal(
    xlim = sta_theta_limits, ylim = sta_theta_limits) +
  annotate("text", x = -Inf, y = Inf, label = sta_theta_equation,
    hjust = -0.08, vjust = 1.5, size = 3.5)

correlation_plots_sown_coexistence <- (wrap_plots(obs_omega_corr, obs_theta_corr,
    sta_omega_corr, sta_theta_corr, ncol = 2, widths = c(1, 1), heights = c(1, 1), guides = "collect") +
    plot_annotation(tag_levels = "A")
) & theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.tag = element_text(size = 14, face = "bold")
  ) &  guides(
    colour = guide_colourbar(
      direction = "vertical",
      title.position = "top",
      barwidth = grid::unit(0.35, "cm"),
      barheight = grid::unit(3.5, "cm")
    )
  )

##This figure S8
ggsave("figures/FigS5_orthogonal_regressions.pdf", correlation_plots_sown_coexistence, width = 9, height = 9, units = "in")
#ggsave("figures/FigS5_orthogonal_regressions.tiff", correlation_plots_sown_coexistence, width = 9, height = 9, units = "in")

##--For the major revision-------------------
##Statistical tests:Does including facilitation increase niche differences more strongly in high-SLA communities
##Fig.S9
library(dplyr)
library(lmerTest)
library(sjPlot)

delta_data_coexistence <- two_coexistence_dfs %>%
  ungroup() %>%
  mutate(
    delta_obs_Omega = Omega_fac - Omega,
    delta_rand_Omega = mean_omega_rand_fac - mean_omega_rand,
    delta_excess_Omega = delta_obs_Omega - delta_rand_Omega,
    delta_sta_Omega = sta.Omega_fac - sta.Omega,
    
    delta_obs_theta = theta_fac - theta,
    delta_rand_theta = mean_theta_rand_fac - mean_theta_rand,
    delta_excess_theta = delta_obs_theta - delta_rand_theta,
    delta_sta_theta = sta.theta_fac - sta.theta) %>%
  filter(
    if_all(
      c(delta_obs_Omega, delta_rand_Omega, delta_excess_Omega, delta_sta_Omega, 
        delta_obs_theta, delta_rand_theta, delta_excess_theta, delta_sta_theta),
      ~ is.finite(.x)
    ) )
summary(delta_data_coexistence)
# Check the number of observations retained per season
table(delta_data_coexistence$season)
delta_data_coexistence$Month <- 1*(delta_data_coexistence$season=="August")

Scalevars.delta_data <- delta_data_coexistence[,c("Month", "N", "F", "SD",  "sownSLA", "sownMPD")]
# do the scaling
Scalevars.delta_data.2 <- scale(Scalevars.delta_data)
cat.vars.delta_data <- delta_data_coexistence[,c("plot", "Block", "composition","Treatment", 
                                                 "delta_obs_Omega", "delta_rand_Omega", "delta_excess_Omega", "delta_sta_Omega", 
                                                 "delta_obs_theta", "delta_rand_theta", "delta_excess_theta", "delta_sta_theta")]
scaled.delta_data <- cbind(cat.vars.delta_data, Scalevars.delta_data.2)
head(scaled.delta_data)

# Observed increase in Omega
hist(scaled.delta_data$delta_obs_Omega)
m_delta_obs_Omega <- lmer(scaled.delta_data$delta_obs_Omega ~(N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | plot) , control = lmerControl(optimizer = "bobyqa"), 
  data = scaled.delta_data, REML = TRUE)
summary(m_delta_obs_Omega)
plot_model(m_delta_obs_Omega, type="pred", terms=c("sownSLA"),show.data=F)#same as the results from Null model before

library(robustlmm)
m_delta_obs_Omega2 <- rlmer(scaled.delta_data$delta_obs_Omega ~(N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | plot), 
                          data = scaled.delta_data)
summary(m_delta_obs_Omega2)
# Compare coefficient estimates on the same response scale
coef_comparison_delta_obs_Omega <- data.frame(
  Term = names(fixef(m_delta_obs_Omega)),
  Gaussian = unname(fixef(m_delta_obs_Omega)),
  Robust = unname(
    fixef(m_delta_obs_Omega2)[names(fixef(m_delta_obs_Omega))]
  )
)
coef_comparison_delta_obs_Omega

#par(mfrow = c(1, 3))
# 1. Normality: approximately follows the straight line?
qqnorm(residuals(m_delta_obs_Omega) / sigma(m_delta_obs_Omega), main = "Residual Q–Q plot")
qqline(residuals(m_delta_obs_Omega) / sigma(m_delta_obs_Omega), col = "red", lwd = 2)

# 2. Model fit: no curve or funnel-shaped spread?
plot(fitted(m_delta_obs_Omega), residuals(m_delta_obs_Omega) / sigma(m_delta_obs_Omega),
  xlab = "Fitted values",
  ylab = "Scaled residuals",
  main = "Residuals versus fitted",
  pch = 16, col = adjustcolor("black", alpha.f = 0.5)
)
abline(h = 0, col = "red", lty = 2)

# 3. Does residual spread differ between seasons?
boxplot(
  (residuals(m_delta_obs_Omega) / sigma(m_delta_obs_Omega)) ~ factor(model.frame(m_delta_obs_Omega)$Month),
  xlab = "Month",
  ylab = "Scaled residuals",
  main = "Residuals by month"
)
##overall, very bad model fit in lmer

# Increase in the randomized expectation Omega
hist(sqrt(scaled.delta_data$delta_rand_Omega))
m_delta_rand_Omega <- lmer(log(delta_rand_Omega) ~(N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | plot) , control = lmerControl(optimizer = "bobyqa"), 
                         data = scaled.delta_data, REML = TRUE)
summary(m_delta_rand_Omega)
# 1. Normality: approximately follows the straight line?
qqnorm(residuals(m_delta_rand_Omega) / sigma(m_delta_rand_Omega), main = "Residual Q–Q plot")
qqline(residuals(m_delta_rand_Omega) / sigma(m_delta_rand_Omega), col = "red", lwd = 2)

# 2. Model fit: no curve or funnel-shaped spread?
plot(fitted(m_delta_rand_Omega), residuals(m_delta_rand_Omega) / sigma(m_delta_rand_Omega),
     xlab = "Fitted values",
     ylab = "Scaled residuals",
     main = "Residuals versus fitted",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5)
)
abline(h = 0, col = "red", lty = 2)

# Increase in standardized Omega 
hist(scaled.delta_data$delta_sta_Omega)
m_delta_sta_Omega <- lmer(delta_sta_Omega ~(N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | plot) , control = lmerControl(optimizer = "bobyqa"), 
                           data = scaled.delta_data, REML = TRUE)
summary(m_delta_sta_Omega)
plot_model(m_delta_sta_Omega, type="pred", terms=c("sownSLA"),show.data=F)#same as the results from Null model before
plot_model(m_delta_sta_Omega, type="pred", terms=c("SD", "N"),show.data=F)#same as the results from Null model before

library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Extract coefficients and matching 95% confidence intervals -------------
delta_coef <- as.data.frame(coef(summary(m_delta_sta_Omega)))

# Fixed effects only; profile intervals, as in your original confint() call
delta_CI <- confint(m_delta_sta_Omega, parm = "beta_", method = "profile", level = 0.95)

# Match intervals by term name rather than row position
model_delta_sta_Result <- data.frame(
  Variable = rownames(delta_coef),
  Estimate = delta_coef$Estimate,
  lower = delta_CI[rownames(delta_coef), "2.5 %"],
  upper = delta_CI[rownames(delta_coef), "97.5 %"],
  P = delta_coef[["Pr(>|t|)"]]
) %>%
  filter(Variable != "(Intercept)") %>%
  mutate(
    Colour = factor(
      case_when(
        lower > 0 ~ "1",
        upper < 0 ~ "-1",
        TRUE ~ "0"
      ),
      levels = c("-1", "0", "1")
    )
  )

# Check before plotting
model_delta_sta_Result

# 2. Set coefficient order and readable labels -----------------------------
interaction_order <- c(
  "N:F:SD",
  "SD:sownSLA",
  "SD:sownMPD",
  "N:F",
  "F:sownSLA",
  "F:sownMPD",
  "F:SD",
  "N:sownSLA",
  "N:sownMPD",
  "N:SD",
  "Month",
  "sownSLA",
  "sownMPD",
  "F",
  "N",
  "SD"
)

# Check that every coefficient has a plotting position
stopifnot(all(model_delta_sta_Result$Variable %in% interaction_order))

model_delta_sta_Result$Variable <- factor(model_delta_sta_Result$Variable,  levels = interaction_order)

# labels: SLA appears as a subscript
term_labels <- c("N:F:SD" = "N:F:SR",
  "SD:sownSLA" = "SR:CWM[SLA]",
  "SD:sownMPD" = "SR:MPD[SLA]",
  "N:F" = "Nitrogen:Fungicide",
  "F:sownSLA" = "Fungicide:CWM[SLA]",
  "F:sownMPD" = "Fungicide:MPD[SLA]",
  "F:SD" = "Fungicide:SR",
  "N:sownSLA" = "Nitrogen:CWM[SLA]",
  "N:sownMPD" = "Nitrogen:MPD[SLA]",
  "N:SD" = "Nitrogen:SR",
  "Month" = "Season",
  "sownSLA" = "'Sown'~CWM[SLA]",
  "sownMPD" = "'Sown'~MPD[SLA]",
  "F" = "Fungicide",
  "N" = "Nitrogen",
  "SD" = "'Sown richness'"
)

# 3. Coefficient plot ------------------------------------------------------
Plot.model_delta_sta_Result <- ggplot(
  model_delta_sta_Result,
  aes(x = Variable, y = Estimate, fill = Colour)
) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.3,
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "-1" = "#d95f02",
      "0" = "gray",
      "1" = "#1b9e77"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = function(x) {
      parse(text = unname(term_labels[x]))
    }
  ) +
  labs(
    x = NULL,
    y = "Effect size",
    title = "Change in standardized niche differences"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(
      colour = "gray90",
      linetype = "dotted"
    ),
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.5
    ),
    plot.margin = margin(5, 8, 5, 5),
    legend.position = "none",
    plot.title = element_text(
      size = 11, hjust = 0.5, face = "plain"
    )
  ) + coord_flip()
Plot.model_delta_sta_Result

library(emmeans)
library(ggplot2)
# Use exactly the data included in the model
model_data_sta_Omega <- model.frame(m_delta_sta_Omega)

sla_center <- mean(delta_data_coexistence$sownSLA, na.rm = TRUE)
sla_scale  <- sd(delta_data_coexistence$sownSLA, na.rm = TRUE)

sla_values <- seq(min(as.numeric(model_data_sta_Omega$sownSLA)), max(as.numeric(model_data_sta_Omega$sownSLA)),
  length.out = 100
)

# Predicted relationship at the means of the other predictors
pred_SLA <- emmeans(
  m_delta_sta_Omega,
  ~ sownSLA,
  at = list(
    sownSLA = sla_values,
    Month = mean(as.numeric(model_data_sta_Omega$Month)),
    N = mean(as.numeric(model_data_sta_Omega$N)),
    F = mean(as.numeric(model_data_sta_Omega$F)),
    SD = mean(as.numeric(model_data_sta_Omega$SD)),
    sownMPD = mean(as.numeric(model_data_sta_Omega$sownMPD))
  ),
  lmer.df = "satterthwaite"
)

pred_SLA_sta_Omega <- as.data.frame(
  confint(pred_SLA, level = 0.95, adjust = "none")
)

# Restore SLA to its original units
stopifnot(length(sla_center) == 1, length(sla_scale) == 1)
pred_SLA_sta_Omega$sownSLA_unscaled <-pred_SLA_sta_Omega$sownSLA * sla_scale + sla_center

names(pred_SLA_sta_Omega)[c(2,5:6)]<-c("predicted_delta", "lower_delta", "upper_delta")

# Plot the predicted relationship and pointwise 95% confidence interval
SLA_delta_plot <- ggplot(
  pred_SLA_sta_Omega,
  aes(x = sownSLA_unscaled, y = predicted_delta)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_ribbon(
    aes(ymin = lower_delta, ymax = upper_delta),
    fill = "gray",
    alpha = 0.20
  ) +
  geom_line(
    colour = "black",
    linewidth = 0.9
  ) +
  labs(
    x = expression("Sown CWM"[SLA]),
    y = paste0(
      "Predicted change in standardized\n",
      "niche differences"
    )
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.5
    ),
    axis.ticks = element_line(
      colour = "black", linewidth = 0.5
    )
  )
SLA_delta_plot

# Values on the scale used for model fitting
SD_values <- sort(unique(as.numeric(model_data_sta_Omega$SD)))
N_values  <- sort(unique(as.numeric(model_data_sta_Omega$N)))

# Check the mapping before plotting
SD_values
N_values

# Original richness labels: change if your model includes other levels
richness_labels <- c(4, 8)

stopifnot(
  length(SD_values) == length(richness_labels),
  length(N_values) == 2
)

# Predictions at each observed richness × nitrogen combination
# Other predictors are held at their means
pred_N_SD_sta_Omega <- emmeans(
  m_delta_sta_Omega,
  ~ N | SD,
  at = list(
    N = N_values,
    SD = SD_values,
    F = mean(as.numeric(model_data_sta_Omega$F)),
    Month = mean(as.numeric(model_data_sta_Omega$Month)),
    sownSLA = mean(as.numeric(model_data_sta_Omega$sownSLA)),
    sownMPD = mean(as.numeric(model_data_sta_Omega$sownMPD))
  ),
  lmer.df = "satterthwaite"
)

pred_N_SD_sta_Omega <- as.data.frame(
  confint(pred_N_SD_sta_Omega, level = 0.95, adjust = "none")
)

# Restore readable labels
pred_N_SD_sta_Omega$Richness <- richness_labels[
  match(pred_N_SD_sta_Omega$SD, SD_values)
]

pred_N_SD_sta_Omega$Nitrogen <- factor(
  pred_N_SD_sta_Omega$N,
  levels = N_values,
  labels = c("Without nitrogen", "With nitrogen")
)

names(pred_N_SD_sta_Omega)[c(3,6:7)]<-c("predicted_delta", "lower_delta", "upper_delta")
N_SD_delta_plot <- ggplot(pred_N_SD_sta_Omega, aes(x = Richness, y = predicted_delta,
    colour = Nitrogen, group = Nitrogen) ) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5
  ) + geom_errorbar(
    aes(ymin = lower_delta, ymax = upper_delta),
    position = position_dodge(width = 0.7),
    width = 0.15, linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  scale_colour_manual(
    values = c("Without nitrogen" = "#99CCFF", "With nitrogen" = "#336699"
    ), name = NULL) +
  scale_x_continuous(breaks = richness_labels, expand = expansion(mult = 0.15)) +
  labs(x = "Sown species richness",
    y = "Predicted change in standardized\nniche differences"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "bottom"
  )

N_SD_delta_plot

# 4. Combine: coefficient plot left; predictions stacked right -------------
figure_layout <- "
AB
AC
"

combined_delta_sta_Omega <- (
  Plot.model_delta_sta_Result +
    SLA_delta_plot +
    N_SD_delta_plot +
    plot_layout(
      design = figure_layout,
      widths = c(1.25, 1),
      heights = c(1, 1)
    ) +
    plot_annotation(tag_levels = "A")
) &
  theme(
    plot.tag = element_text(size = 14, face = "plain")
  )

combined_delta_sta_Omega
##This figure S9
ggsave(filename = "figures/Delta_standardized_Omega.pdf", plot = combined_delta_sta_Omega, width = 11, height = 8,  units = "in")
ggsave(filename = "figures/Delta_standardized_Omega.tiff", plot = combined_delta_sta_Omega, width = 9, height = 8,  units = "in")

# 1. Normality: approximately follows the straight line?
qqnorm(residuals(m_delta_sta_Omega) / sigma(m_delta_sta_Omega), main = "Residual Q–Q plot")
qqline(residuals(m_delta_sta_Omega) / sigma(m_delta_sta_Omega), col = "red", lwd = 2)

# 2. Model fit: no curve or funnel-shaped spread?
plot(fitted(m_delta_sta_Omega), residuals(m_delta_sta_Omega) / sigma(m_delta_sta_Omega),
     xlab = "Fitted values",
     ylab = "Scaled residuals",
     main = "Residuals versus fitted",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5)
)
abline(h = 0, col = "red", lty = 2)

# 3. Does residual spread differ between seasons?
boxplot(
  (residuals(m_delta_sta_Omega) / sigma(m_delta_sta_Omega)) ~ factor(model.frame(m_delta_sta_Omega)$Month),
  xlab = "Month",
  ylab = "Scaled residuals",
  main = "Residuals by month"
)

m_delta_sta_Omega_season <- update(m_delta_sta_Omega, . ~ . + sownSLA:Month)
anova(m_delta_sta_Omega_season, type = 2)
summary(m_delta_sta_Omega_season)
