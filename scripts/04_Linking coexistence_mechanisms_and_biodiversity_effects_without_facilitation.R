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
#signed square-root transformation to improve residual normality and homogeneity of variance
df_biodiversity.effects_c <- df_biodiversity.effects_c %>%
  dplyr::mutate(
    CE_sqrt = sign(CE) * sqrt(abs(CE)),
    SE_intra_sqrt = sign(SE_intra) * sqrt(abs(SE_intra))
)

df_ini_coexistence_bio.effects <- merge(combined_sown_coex_new, df_biodiversity.effects_c[,c(1:2,8:17)], by=c("season", "plot", "Block"))
df_ini_coexistence_bio.effects$Month<-1*(df_ini_coexistence_bio.effects$season=="August")

##exclude SD==20, because all 20 species plots have the same predicted coexistence_mechanisms
df_ini_coexistence_bio.effects1<-subset(df_ini_coexistence_bio.effects, SD!=20)
head(df_ini_coexistence_bio.effects1)

Scalevars_both_seasons <-df_ini_coexistence_bio.effects1[,c("Month", "N", "F", "SD", "sownSLA", "sownMPD", "sta.Omega", "sta.theta", "sta.differential", "sta.differential_new","sta.Omega_all_pairs_zero", "NE", "CE_sqrt", "SE", "SE_inter", "SE_intra_sqrt")]
names(Scalevars_both_seasons)[c(13,16)]<-c("CE", "SE_intra")
# do the scaling
Scalevars_both_seasons2 <- scale(Scalevars_both_seasons)
corrs.2<-cor(Scalevars_both_seasons2, use = "pairwise.complete.obs")
pdf("outputs/correlation_matrix.pdf", width = 12, height = 12)

corrplot::corrplot(
  corrs.2,
  method = "square",
  type = "upper",
  order = "hclust",
  col = colorRampPalette(c("#336699", "white", "#990000"))(200),
  addCoef.col = "black",
  number.cex = 0.65,
  tl.col = "black",
  tl.cex = 0.7,
  tl.srt = 45,
  diag = FALSE,
  mar = c(0, 0, 1, 0)
)
dev.off()

cat.vars.both.seasons <- df_ini_coexistence_bio.effects1[,c("plot", "Block", "composition", "season", "N")]
names(cat.vars.both.seasons)[5] <- c("Nitrogen")
cat.vars.both.seasons$group <- interaction(cat.vars.both.seasons$season, cat.vars.both.seasons$Nitrogen)
data.scaled.both.seasons <- cbind(cat.vars.both.seasons, Scalevars_both_seasons2)
head(data.scaled.both.seasons)
data.scaled.both.seasons$season<-as.factor(data.scaled.both.seasons$season)
data.scaled.both.seasons$group<-as.factor(data.scaled.both.seasons$group)

##SEM model including indirect interactions (the SEM for the first submission)
model_string_both_seasons_N_indirect_1st <- '
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

fituw_both_seasons_N_indirect_1st <- sem(model_string_both_seasons_N_indirect_1st, data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect_1st, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect_1st, standardized =TRUE, rsquare=TRUE)
sort.pars(fituw_both_seasons_N_indirect_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect_1st <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect_1st, data.scaled.both.seasons)
constrain_model_both_seasons_N_indirect_1st$final_fit
best_string_both_seasons_N_indirect_1st <- constrain_model_both_seasons_N_indirect_1st$final_model_string

confituw_both_seasons_N_indirect_1st <- sem(unlist(best_string_both_seasons_N_indirect_1st), data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect_1st, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect_1st, c("pvalue", "cfi","rmsea","srmr"))##p:0.464, cfi:0.998, rmsea:0.01 ,srmr0.072
sort.pars(confituw_both_seasons_N_indirect_1st, response = "CE", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect_1st, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))
sort.pars(confituw_both_seasons_N_indirect_1st, response = "sta.differential", group.names = unique(data.scaled.both.seasons$group))

#Constrain model mannually
#F on SEAS (dominance effects) August vs. June
mcs1_sown_N_indirect_1.0_1st <- constrain_path_cfa(confituw_both_seasons_N_indirect_1st, "SE_inter ~ F", coef_name = "coef244Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_1.0_1st <- sem(mcs1_sown_N_indirect_1.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_1.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_1.0_1st, "SE_inter ~ F", coef_name = "coef244Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_1.0_1st <- sem(mcs2_sown_N_indirect_1.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_1.0_1st)
lavTestLRT(fc1_sown_N_indirect_1.0_1st, confituw_both_seasons_N_indirect_1st)#accepted
sort.pars(fc1_sown_N_indirect_1.0_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_1.0_1st, response = "CE", group.names = unique(data.scaled.both.seasons$group))

#Sown MPD on CE (complementarity effects) N vs. No N
mcs1_sown_N_indirect_2.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_1.0_1st, "CE ~ sownMPD", coef_name = "coef246N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_2.0_1st <- sem(mcs1_sown_N_indirect_2.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_2.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_2.0_1st, "CE ~ sownMPD", coef_name = "coef246noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_2.0_1st <- sem(mcs2_sown_N_indirect_2.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_2.0_1st)
lavTestLRT(fc1_sown_N_indirect_1.0_1st, fc1_sown_N_indirect_2.0_1st)#rejected

#sta.Omega on SE_inter N vs. noN
mcs1_sown_N_indirect_3.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_1.0_1st, "SE_inter ~ sta.Omega", coef_name = "coef248N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0_1st <- sem(mcs1_sown_N_indirect_3.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_3.0_1st, "SE_inter ~ sta.Omega", coef_name = "coef248noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0_1st <- sem(mcs2_sown_N_indirect_3.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0_1st)
lavTestLRT(fc1_sown_N_indirect_1.0_1st, fc1_sown_N_indirect_3.0_1st)#accepted
sort.pars(fc1_sown_N_indirect_3.0_1st, response = "CE", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0_1st, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0_1st, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownSLA on theta N vs. noN
mcs1_sown_N_indirect_4.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_3.0_1st, "sta.theta ~ sownSLA", coef_name = "coef250N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_4.0_1st <- sem(mcs1_sown_N_indirect_4.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_4.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_4.0_1st, "sta.theta ~ sownSLA", coef_name = "coef250NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_4.0_1st <- sem(mcs2_sown_N_indirect_4.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_4.0_1st)
lavTestLRT(fc1_sown_N_indirect_3.0_1st, fc1_sown_N_indirect_4.0_1st)#accepted 
sort.pars(fc1_sown_N_indirect_4.0_1st, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#SD on sta.Omega June vs. August
mcs1_sown_N_indirect_5.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_4.0_1st, "sta.Omega ~ SD", coef_name = "coef251Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_5.0_1st <- sem(mcs1_sown_N_indirect_5.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_5.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_5.0_1st, "sta.Omega ~ SD", coef_name = "coef251Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_5.0_1st <- sem(mcs2_sown_N_indirect_5.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_5.0_1st)
lavTestLRT(fc1_sown_N_indirect_4.0_1st, fc1_sown_N_indirect_5.0_1st)#accepted
sort.pars(fc1_sown_N_indirect_5.0_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on sta.theta June vs. August
mcs1_sown_N_indirect_6.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_5.0_1st, "sta.theta ~ sownMPD", coef_name = "coef252Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_6.0_1st <- sem(mcs1_sown_N_indirect_6.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_6.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_6.0_1st, "sta.theta ~ sownMPD", coef_name = "coef252Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_6.0_1st <- sem(mcs2_sown_N_indirect_6.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_6.0_1st)
lavTestLRT(fc1_sown_N_indirect_6.0_1st, fc1_sown_N_indirect_5.0_1st)#accepted
sort.pars(fc1_sown_N_indirect_6.0_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SE_inter N vs. NoN
mcs1_sown_N_indirect_7.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_6.0_1st, "SE_inter ~ sownMPD", coef_name = "coef253N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_7.0_1st <- sem(mcs1_sown_N_indirect_7.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_7.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_7.0_1st, "SE_inter ~ sownMPD", coef_name = "coef253NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_7.0_1st <- sem(mcs2_sown_N_indirect_7.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_7.0_1st)
lavTestLRT(fc1_sown_N_indirect_6.0_1st, fc1_sown_N_indirect_7.0_1st)#accepted
sort.pars(fc1_sown_N_indirect_7.0_1st, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sta.differential on June vs. August
mcs1_sown_N_indirect_8.0_1st <- constrain_path_cfa(fc1_sown_N_indirect_7.0_1st, "sta.Omega ~ sta.differential", coef_name = "coef254Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_8.0_1st <- sem(mcs1_sown_N_indirect_8.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_8.0_1st <- constrain_path_cfa(fc1a_sown_N_indirect_8.0_1st, "sta.Omega ~ sta.differential", coef_name = "coef254Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_8.0_1st <- sem(mcs2_sown_N_indirect_8.0_1st, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_8.0_1st)
lavTestLRT(fc1_sown_N_indirect_8.0_1st, fc1_sown_N_indirect_7.0_1st)#rejected

##the best model is fc1_sown_N_indirect_7.0_1st, fig3 was made based on this SEM in PowerPoint
fitMeasures(fc1_sown_N_indirect_7.0_1st, c("pvalue", "cfi","rmsea","srmr"))##p:0.283, cfi:0.986, rmsea:0.037 ,srmr0.075
summary(fc1_sown_N_indirect_7.0_1st, standardized =TRUE, rsquare=TRUE)
sort.pars(fc1_sown_N_indirect_7.0_1st, response = "CE", group.names = unique(data.scaled.both.seasons$group))

##-These code was not used anymore (could skip this part)------------------------------------
###To revise reviewer's comment on the circularity between sta.Omega and sta.differential: 
##sta.Omega~sta.differential was removed, and sta.Omega~~sta.differential was added.
model_string_both_seasons_N_indirect <- '
CE ~ F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
SE_inter ~  F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
SE_intra ~ F + SD + sownSLA + sownMPD + sta.Omega + sta.theta + sta.differential
sta.Omega ~ F + SD + sownSLA + sownMPD 
sta.theta ~ F + SD + sownSLA + sownMPD 
sta.differential ~ F + SD + sownSLA + sownMPD

sta.Omega~~sta.theta 
sta.Omega~~sta.differential
sta.theta~~sta.differential
sownSLA~~sownMPD
'

fituw_both_seasons_N_indirect <- sem(model_string_both_seasons_N_indirect, data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
fitMeasures(fituw_both_seasons_N_indirect, c("cfi","rmsea","srmr"))## cfi:1, rmsea:0.00 ,srmr0.025
summary(fituw_both_seasons_N_indirect, standardized =TRUE, rsquare=TRUE)
sort.pars(fituw_both_seasons_N_indirect, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))
show(fituw_both_seasons_N_indirect)

##now I need to constrain the SEM model
constrain_model_both_seasons_N_indirect <- test_stepwise_constraints_cfa(fituw_both_seasons_N_indirect, data.scaled.both.seasons)
constrain_model_both_seasons_N_indirect$final_fit
best_string_both_seasons_N_indirect <- constrain_model_both_seasons_N_indirect$final_model_string

confituw_both_seasons_N_indirect <- sem(unlist(best_string_both_seasons_N_indirect), data = data.scaled.both.seasons,  group = "group", fixed.x = FALSE)
summary(confituw_both_seasons_N_indirect, standardized =TRUE, rsquare=TRUE)
fitmeasures(confituw_both_seasons_N_indirect, c("pvalue","cfi","rmsea","srmr"))##p:0.460, cfi:0.999, rmsea:0.011 ,srmr0.067
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

#sta.Omega on SE_inter N vs. noN
mcs1_sown_N_indirect_3.0 <- constrain_path_cfa(fc1_sown_N_indirect_1.0, "SE_inter ~ sta.Omega", coef_name = "coef248N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_3.0 <- sem(mcs1_sown_N_indirect_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_3.0 <- constrain_path_cfa(fc1a_sown_N_indirect_3.0, "SE_inter ~ sta.Omega", coef_name = "coef248noN", groups_constrain = c(2,4))
fc1_sown_N_indirect_3.0 <- sem(mcs2_sown_N_indirect_3.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_3.0)
lavTestLRT(fc1_sown_N_indirect_1.0, fc1_sown_N_indirect_3.0)#accepted
sort.pars(fc1_sown_N_indirect_3.0, response = "CE", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0, response = "SE_inter", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0, response = "SE_intra", group.names = unique(data.scaled.both.seasons$group))
sort.pars(fc1_sown_N_indirect_3.0, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on sta.Omega N vs. noN
mcs1_sown_N_indirect_4.0 <- constrain_path_cfa(fc1_sown_N_indirect_3.0, "sta.Omega ~ sownMPD", coef_name = "coef250N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_4.0 <- sem(mcs1_sown_N_indirect_4.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_4.0 <- constrain_path_cfa(fc1a_sown_N_indirect_4.0, "sta.Omega ~ sownMPD", coef_name = "coef250NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_4.0 <- sem(mcs2_sown_N_indirect_4.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_4.0)
lavTestLRT(fc1_sown_N_indirect_3.0, fc1_sown_N_indirect_4.0)#accepted 
sort.pars(fc1_sown_N_indirect_4.0, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

#SD on sta.Omega June vs. August
mcs1_sown_N_indirect_5.0 <- constrain_path_cfa(fc1_sown_N_indirect_4.0, "sta.Omega ~ SD", coef_name = "coef251Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_5.0 <- sem(mcs1_sown_N_indirect_5.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_5.0 <- constrain_path_cfa(fc1a_sown_N_indirect_5.0, "sta.Omega ~ SD", coef_name = "coef251Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_5.0 <- sem(mcs2_sown_N_indirect_5.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_5.0)
lavTestLRT(fc1_sown_N_indirect_4.0, fc1_sown_N_indirect_5.0)#accepted
sort.pars(fc1_sown_N_indirect_5.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownSLA on sta.theta N vs. NoN
mcs1_sown_N_indirect_6.0 <- constrain_path_cfa(fc1_sown_N_indirect_5.0, "sta.theta ~ sownSLA", coef_name = "coef252N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_6.0 <- sem(mcs1_sown_N_indirect_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_6.0 <- constrain_path_cfa(fc1a_sown_N_indirect_6.0, "sta.theta ~ sownSLA", coef_name = "coef252NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_6.0 <- sem(mcs2_sown_N_indirect_6.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_6.0)
lavTestLRT(fc1_sown_N_indirect_6.0, fc1_sown_N_indirect_5.0)#accepted
sort.pars(fc1_sown_N_indirect_6.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on SE_inter N vs. NoN
mcs1_sown_N_indirect_7.0 <- constrain_path_cfa(fc1_sown_N_indirect_6.0, "SE_inter ~ sownMPD", coef_name = "coef253N", groups_constrain = c(1,3))
fc1a_sown_N_indirect_7.0 <- sem(mcs1_sown_N_indirect_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_7.0 <- constrain_path_cfa(fc1a_sown_N_indirect_7.0, "SE_inter ~ sownMPD", coef_name = "coef253NoN", groups_constrain = c(2,4))
fc1_sown_N_indirect_7.0 <- sem(mcs2_sown_N_indirect_7.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_7.0)
lavTestLRT(fc1_sown_N_indirect_7.0, fc1_sown_N_indirect_6.0)#accepted
sort.pars(fc1_sown_N_indirect_7.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

#sownMPD on sta.theta N vs. NoN
mcs1_sown_N_indirect_8.0 <- constrain_path_cfa(fc1_sown_N_indirect_7.0, "sta.theta ~ sownMPD", coef_name = "coef254Aug", groups_constrain = c(1,2))
fc1a_sown_N_indirect_8.0 <- sem(mcs1_sown_N_indirect_8.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
mcs2_sown_N_indirect_8.0 <- constrain_path_cfa(fc1a_sown_N_indirect_8.0, "sta.theta ~ sownMPD", coef_name = "coef254Jun", groups_constrain = c(3,4))
fc1_sown_N_indirect_8.0 <- sem(mcs2_sown_N_indirect_8.0, data = data.scaled.both.seasons, group = "group", fixed.x = FALSE)
extract_paths_cfa(fc1_sown_N_indirect_8.0)
lavTestLRT(fc1_sown_N_indirect_8.0, fc1_sown_N_indirect_7.0)#accepted
sort.pars(fc1_sown_N_indirect_8.0, response = "sta.theta", group.names = unique(data.scaled.both.seasons$group))

##the best model is fc1_sown_N_indirect_8.0
fitMeasures(fc1_sown_N_indirect_8.0, c("pvalue", "cfi","rmsea","srmr"))##p:0.327, cfi:0.99, rmsea:0.032 ,srmr0.074
summary(fc1_sown_N_indirect_8.0, standardized =TRUE, rsquare=TRUE)
sort.pars(fc1_sown_N_indirect_8.0, response = "sta.Omega", group.names = unique(data.scaled.both.seasons$group))

##---------
##make partial plots
final_model_fit_final <- fc1_sown_N_indirect_7.0_1st
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
pr_SLA_Omega <- part.res.sem(final_model_fit_final, c("sownSLA"), "sta.Omega", Bio_list, groups = 1:length(Bio_list), level = 1)

pr_MPD_Omega <- part.res.sem(final_model_fit_final, c("sownMPD"), "sta.Omega", Bio_list, groups = 1:length(Bio_list), level = 1)

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
pdf("figures/fig4_new.pdf", width = 9, height = 12)  

# Layout: 2 columns for panels, plus space for legend
layout_matrix <- matrix(c(
  1, 2, 7,  # Row 1: Panel A, Panel B, empty for legend alignment
  3, 4, 7,  # Row 2: Panel C, Panel D
  5, 6, 7   # Row 3: Panel E, Panel F
), nrow = 3, byrow = TRUE)

# Column widths: first column, second column, legend column
layout(layout_matrix, widths = c(1, 1, 0.55))

# Set margins for panels (can adjust as needed)
par(mar = c(4.5, 4.5, 2, 1))

# ---- Panel A ----
plot.res.sem3_auto1(
  Bio_list_unscaled, "sownSLA", residual_list = pr_SLA_Omega,
  cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
  xlab = "Functional composition (sown CWMSLA)",
  ylab = "Niche differences",
  group_labels = group_labels,
  group_combinations = list(c(1, 2, 3, 4))
)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("A", side = 3, adj = -0.08, line = 0.2, font = 2)

# ---- Panel B ----
plot.res.sem3_1(
  Bio_list_unscaled, "sownSLA", residual_list = pr_SLA_theta,
  cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
  xlab = "Functional composition (sown CWMSLA)",
  ylab = "Competitive ability differences",
  group_labels = group_labels,
  group_combinations = list(c(1, 3),  c(2, 4)),
  show_regression = c(TRUE, TRUE)
)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("B", side = 3, adj = -0.08, line = 0.2, font = 2)

# ---- Panel C ----
plot.res.sem3_auto1(Bio_list_unscaled, "sta.Omega", residual_list = pr_Omega_CE, cols = cols,
                    alpha = 0.8, 
                    xlab = "Niche differences",  
                    ylab = "Complementarity effects",
                    group_labels = group_labels, group_combinations = list(c(1, 2, 3, 4)),
                    pch = 21, cex = 1.1, lwd = 0.5
)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("C", side = 3, adj = -0.12, line = 0.2, font = 2)

# ---- Panel D ----
plot.res.sem3_1(Bio_list_unscaled, "sta.theta", residual_list = pr_theta_CE,
                cols = cols, alpha = 0.8, pch = 21, cex = 1.1, lwd = 0.5,
                xlab = "Competitive ability differences",  
                ylab = "Complementarity effects", 
                group_combinations = list(c(1, 3),  c(2, 4)), show_regression = c(FALSE,FALSE))
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("D", side = 3, adj = -0.12, line = 0.2, font = 2)

# ---- Panel E ----
plot.res.sem3_1(Bio_list_unscaled, "sta.Omega", residual_list = pr_Omega_SEAS,
                cols = cols, alpha = 0.8,
                xlab = "Niche differences", 
                ylab = "Interspecific selection effects", 
                group_labels = group_labels,
                group_combinations = list(c(1, 3),  c(2, 4)), show_regression = c(TRUE, TRUE),
                pch = 21, cex = 1.1, lwd = 0.5)
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("E", side = 3, adj = -0.12, line = 0.2, font = 2)

# ---- Panel F ----
plot.res.sem3_auto1(Bio_list_unscaled, "sta.theta", residual_list = pr_theta_SEAS, cols = cols,
                    alpha = 0.8, xlab = "Competitive ability differences", ylab = "Interspecific selection effects",
                    group_labels = group_labels, group_combinations = list(c(1, 2, 3, 4)),
                    pch = 21, cex = 1.1, lwd = 0.5) 
abline(h = 0, lty = 2, col = "grey9", lwd = 1.2)
mtext("F", side = 3, adj = -0.12, line = 0.2, font = 2)

# Legend in the right column (shared across all rows)
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")
legend("center", pch = 19, lty = 1, col = c("#99CCFF", "#336699"),
       legend = c("Without nitrogen", "With nitrogen"),
       bty = "n", pt.cex = 1.2, cex = 0.95, xpd = NA)
dev.off()

#--------------------------------------------------
##As reviewer requested, I start to plot the direct vs. indirect effects through niche and competitive ability differences on diversity effects
##extract the standardized effect sizes from: fc1_sown_N_indirect_7.0_1st (the best SEM without facilitation)
##Note that this visualization will be continued in R script.06
library(lavaan)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(manymome)
library(parallel)
library(tibble)
final_model_fit_final <- fc1_sown_N_indirect_7.0_1st
group_names <- lavInspect(final_model_fit_final, "group.label")
data.frame(group_number = seq_along(group_names), group_label = group_names)

# Inspect how uncertainty was calculated in the original model.
lavInspect(final_model_fit_final, "options")[c("estimator", "se", "test")]
stopifnot(lavInspect(final_model_fit_final, "converged"))

# Monte Carlo Estimates for 'indirect_effects' and 'cond_indirect_effects'
mc_results <- do_mc(fit = final_model_fit_final, R = 10000, seed = 123)
predictors <- c("F", "SD", "sownSLA", "sownMPD")
responses  <- c("CE", "SE_inter", "SE_intra")
mediators  <- c("sta.Omega", "sta.theta", "sta.differential")

get_effect <- function(x, y, g, m = NULL, level = 0.95) {
  manymome::indirect_effect(
    x = x,
    y = y,
    m = m,
    fit = final_model_fit_final,
    group = g,
    standardized_x = FALSE,
    standardized_y = FALSE,
    mc_ci = TRUE,
    mc_out = mc_results,
    level = level,
    parallel = FALSE,
    progress = FALSE
  )
}

effect_row <- function(object, g, x, y, type, level) {
  ci <- as.numeric(stats::confint(object, level = level))
  stopifnot(length(ci) == 2L)
  
  tibble::tibble(
    group = group_names[g],
    predictor = x,
    response = y,
    effect_type = type,
    estimate = as.numeric(coef(object)),
    conf.low = ci[1],
    conf.high = ci[2]
  )
}

# One task per group × predictor × response.
tasks <- expand.grid(
  g = seq_along(group_names),
  x = c(predictors, mediators),
  y = responses,
  stringsAsFactors = FALSE
)

calculate_task <- function(i, level = 0.95) {
  
  g <- tasks$g[i]
  x <- tasks$x[i]
  y <- tasks$y[i]
  
  direct <- get_effect(x, y, g, level = level)
  
  if (x %in% mediators) {
    return(
      effect_row(direct, g, x, y, "Direct", level = level)
    )
  }
  
  via_ND <- get_effect(
    x, y, g, m = "sta.Omega", level = level
  )
  
  via_FD <- get_effect(
    x, y, g, m = "sta.theta", level = level
  )
  
  via_D <- get_effect(
    x, y, g, m = "sta.differential", level = level
  )
  
  total_indirect <- via_ND + via_FD + via_D
  
  dplyr::bind_rows(
    effect_row(
      direct, g, x, y, "Direct", level = level
    ),
    effect_row(
      via_ND, g, x, y,
      "Via niche differences", level = level
    ),
    effect_row(
      via_FD, g, x, y,
      "Via competitive ability differences", level = level
    ),
    effect_row(
      via_D, g, x, y,
      "Via community-pair differential", level = level
    ),
    effect_row(
      total_indirect, g, x, y,
      "Total indirect", level = level
    )
  )
}
# The wrapper ensures workers are stopped on completion or error.
run_parallel_effects <- function(n_workers = 20L, level = 0.95) {
  
  cl <- parallel::makePSOCKcluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(manymome)
    NULL
  })
  
  # Copy shared objects once to each worker.
  parallel::clusterExport(
    cl,
    varlist = c(
      "final_model_fit_final", "group_names", "mc_results",
      "mediators", "tasks",
      "get_effect", "effect_row", "calculate_task"
    ),
    envir = environment()
  )
  
  # Load balancing assigns a new task whenever a worker finishes.
  results <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(tasks)),
    fun = calculate_task,
    level = level
  )
  
  dplyr::bind_rows(results)
}

# Start with multiple workers (time consuming).
effects_all <- run_parallel_effects(n_workers = 20L) %>%
  mutate(
    group = factor(group, levels = group_names),
    supported = conf.low > 0 | conf.high < 0
  )
View(effects_all)#Unstandardized effect sizes

effects_all %>% dplyr::filter(effect_type=="Via niche differences", response == "CE")
effects_all %>% dplyr::filter(effect_type=="Direct",predictor=="sta.Omega", response == "CE")
effects_all %>% dplyr::filter(effect_type=="Direct",predictor=="SD", response == "CE")
effects_all %>% dplyr::filter(effect_type=="Direct",predictor=="sownMPD", response == "CE")
effects_all %>% dplyr::filter(effect_type=="Direct",predictor=="sownSLA", response == "CE")

effects_all %>% dplyr::filter(effect_type=="Via niche differences", response == "SE_intra")
effects_all %>% dplyr::filter(effect_type=="Via competitive ability differences", response == "SE_intra")

table(effects_all$response)
table(effects_all$effect_type)
table(effects_all$predictor)

# Selected pathways, in the requested top-to-bottom order.
# seasonal = FALSE: effect is shared between June and August -> one black point.
# seasonal = TRUE: retain season-specific effects and colors.
# august_only = TRUE: show only August.
# june_only = TRUE: show only June.
selected_paths <- tribble(
  ~response, ~predictor, ~effect_type, ~Pathway, ~seasonal, ~august_only, ~june_only,
  
  "CE", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "CE", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE, FALSE,
  "CE", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE, FALSE,

  "SE_intra", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "SE_intra", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE, FALSE,
  "SE_intra", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE, FALSE,

  "SE_intra", "sta.theta", "Direct", "CAD", FALSE, FALSE, FALSE,

  "SE_inter", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "SE_inter", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE,FALSE,
  "SE_inter", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE,FALSE,

  "SE_inter", "sta.theta", "Direct", "CAD", FALSE, FALSE, FALSE,
  "SE_inter", "sownSLA", "Via competitive ability differences", "Sown CWM SLA via CAD", FALSE, FALSE,FALSE,
  "SE_inter", "sownMPD", "Via competitive ability differences", "Sown MPD SLA via CAD", TRUE, FALSE,TRUE,
  "SE_inter", "F", "Via competitive ability differences", "Fungicide via CAD", TRUE, TRUE, FALSE
) %>%   mutate(path_id = row_number())
##These code will continute with the same code in 06_Lining_coexistence_mechanims.R

##Same thing for the SEM with ND (sta.Omega)~~D (sta.differential)
final_model_fit_final_8 <- fc1_sown_N_indirect_8.0
group_names_8 <- lavaan::lavInspect(final_model_fit_final_8, "group.label")
data.frame(group_number = seq_along(group_names_8), group_label = group_names_8)

lavaan::lavInspect(final_model_fit_final_8, "options")[c("estimator", "se", "test")]
stopifnot(lavaan::lavInspect(final_model_fit_final_8, "converged"))

# Generate Monte Carlo draws from the new model.
mc_results_8 <- manymome::do_mc(fit = final_model_fit_final_8, R = 10000, seed = 123)
predictors_8 <- c("F", "SD", "sownSLA", "sownMPD")
responses_8 <- c("CE", "SE_inter", "SE_intra")
mediators_8 <- c("sta.Omega", "sta.theta", "sta.differential")

# ------------------------------------------------------------
# 2. Functions for calculating and extracting effects.
# ------------------------------------------------------------
get_effect_8 <- function(x, y, g, m = NULL, level = 0.95) {
  manymome::indirect_effect(
    x = x,
    y = y,
    m = m,
    fit = final_model_fit_final_8,
    group = g,
    standardized_x = FALSE,
    standardized_y = FALSE,
    mc_ci = TRUE,
    mc_out = mc_results_8,
    level = level,
    parallel = FALSE,
    progress = FALSE
  )
}

effect_row_8 <- function(object, g, x, y, type, level) {
  
  ci <- as.numeric(stats::confint(object, level = level))
  stopifnot(length(ci) == 2L)
  
  tibble::tibble(
    group = group_names_8[g],
    predictor = x,
    response = y,
    effect_type = type,
    estimate = as.numeric(stats::coef(object)),
    conf.low = ci[1],
    conf.high = ci[2]
  )
}

# One task per group × predictor × response.
tasks_8 <- expand.grid(
  g = seq_along(group_names_8),
  x = c(predictors_8, mediators_8),
  y = responses_8,
  stringsAsFactors = FALSE
)

calculate_task_8 <- function(i, level = 0.95) {
  
  g <- tasks_8$g[i]
  x <- tasks_8$x[i]
  y <- tasks_8$y[i]
  direct <- get_effect_8(x, y, g, level = level)
  
  # For coexistence metrics, extract direct effects only.
  if (x %in% mediators_8) {
    return(
      effect_row_8(
        direct, g, x, y, "Direct", level = level
      )
    )
  }
  
  via_ND <- get_effect_8(x, y, g, m = "sta.Omega", level = level)
  
  via_FD <- get_effect_8(x, y, g, m = "sta.theta", level = level)
  
  via_D <- get_effect_8(x, y, g, m = "sta.differential", level = level)
  
  # Sum the three single-mediator pathways using joint draws.
  total_indirect <- via_ND + via_FD + via_D
  
  dplyr::bind_rows(
    effect_row_8(
      direct, g, x, y,
      "Direct", level = level
    ),
    effect_row_8(
      via_ND, g, x, y,
      "Via niche differences", level = level
    ),
    effect_row_8(
      via_FD, g, x, y,
      "Via competitive ability differences", level = level
    ),
    effect_row_8(
      via_D, g, x, y,
      "Via community-pair differential", level = level
    ),
    effect_row_8(
      total_indirect, g, x, y,
      "Total indirect", level = level
    )
  )
}

# ------------------------------------------------------------
# 3. Parallel calculation.
# ------------------------------------------------------------
run_parallel_effects_8 <- function(
    n_workers = 20L,
    level = 0.95
) {
  
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
      "final_model_fit_final_8",
      "group_names_8",
      "mc_results_8",
      "mediators_8",
      "tasks_8",
      "get_effect_8",
      "effect_row_8",
      "calculate_task_8"
    ),
    envir = environment()
  )
  
  results <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(tasks_8)),
    fun = calculate_task_8,
    level = level
  )
  
  dplyr::bind_rows(results)
}

# ------------------------------------------------------------
# 4. Calculate effects with 95% confidence intervals.
# ------------------------------------------------------------
effects_all_8 <- run_parallel_effects_8(n_workers = 20L, level = 0.95) %>%
  dplyr::mutate(group = factor(group, levels = group_names_8),
               supported = conf.low > 0 | conf.high < 0
)
View(effects_all_8)

effects_all_8 %>% dplyr::filter(effect_type == "Via niche differences", response == "CE")
effects_all_8 %>% dplyr::filter(effect_type == "Direct", predictor == "sta.Omega", response == "CE")
effects_all_8 %>% dplyr::filter(effect_type == "Direct", predictor == "SD", response == "CE"  )
effects_all_8 %>% dplyr::filter(effect_type == "Direct", predictor == "sownMPD", response == "CE")
effects_all_8 %>% dplyr::filter(effect_type == "Direct", predictor == "sownSLA", response == "CE")
effects_all_8 %>% dplyr::filter(effect_type == "Via niche differences", response == "SE_intra")
effects_all_8 %>% dplyr::filter(effect_type == "Via competitive ability differences",response == "SE_intra")

table(effects_all_8$response)
table(effects_all_8$effect_type)
table(effects_all_8$predictor)

# Selected pathways, in the requested top-to-bottom order.
# seasonal = FALSE: effect is shared between June and August -> one black point.
# seasonal = TRUE: retain season-specific effects and colors.
# august_only = TRUE: show only August.
# june_only = TRUE: show only June.
selected_paths_without_D_on_ND <- tribble(
  ~response, ~predictor, ~effect_type, ~Pathway, ~seasonal, ~august_only, ~june_only,
  
  "CE", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "CE", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE, FALSE,
  "CE", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE, FALSE,
  "CE", "sownMPD", "Via niche differences", "Sown MPD SLA via ND", TRUE, FALSE, TRUE,
  
  "SE_intra", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "SE_intra", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE, FALSE,
  "SE_intra", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE, FALSE,
  "SE_intra", "sownMPD", "Via niche differences", "Sown MPD SLA via ND", TRUE, FALSE, TRUE,
  
  "SE_intra", "sta.theta", "Direct", "CAD", FALSE, FALSE, FALSE,
  
  "SE_inter", "sta.Omega", "Direct", "ND", FALSE, FALSE, FALSE,
  "SE_inter", "SD", "Via niche differences", "Sown richness via ND", TRUE, TRUE,FALSE,
  "SE_inter", "sownSLA", "Via niche differences", "Sown CWM SLA via ND", FALSE, FALSE,FALSE,
  "SE_inter", "sownMPD", "Via niche differences", "Sown MPD SLA via ND", TRUE, FALSE, TRUE,
  "SE_inter", "sta.theta", "Direct", "CAD", FALSE, FALSE, FALSE,
  "SE_inter", "sownSLA", "Via competitive ability differences", "Sown CWM SLA via CAD", FALSE, FALSE,FALSE,
  "SE_inter", "sownMPD", "Via competitive ability differences", "Sown MPD SLA via CAD", TRUE, FALSE,TRUE,
  "SE_inter", "F", "Via competitive ability differences", "Fungicide via CAD", TRUE, TRUE, FALSE
) %>%   mutate(path_id = row_number())
##These code will continue with the same code in 06_Lining_coexistence_mechanims.R