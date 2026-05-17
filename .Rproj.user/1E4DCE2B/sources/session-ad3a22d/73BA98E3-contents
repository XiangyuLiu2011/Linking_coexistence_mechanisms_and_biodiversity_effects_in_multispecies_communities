##scripts/02_linear_mixed_effects_models_biodiversity_effects
library(lme4)
library(dplyr)
library(lmerTest)
library(sjPlot)
library(multcomp)
library(ggrepel)
library(corrplot)
library(car)
library(performance)
library(emmeans)
library(ggplot2)
library(ggforce)
library(patchwork)
library(visreg)
##Average diversity effects per plot across years but keep the harvests
Mean.Biomass.Tripartite <- AP_Biomass %>%
  group_by(Plot_Nr, Month) %>%
  dplyr::summarise(across(c("NE", "CE", "SE", "SE_inter", "SE_intra"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
Mean.Biomass.Tripartite <- as.data.frame(Mean.Biomass.Tripartite)
#write.csv(df_biodiversity.effects.1, "outputs/Mean.Biomass.Tripartite17_23_new.csv", row.names = FALSE)
##load CWM_SLA and MPD_SLA
plot.info.sla <- read.csv("data/plot.info.sla.updated.csv")
plot.info.sla.vip <- subset(plot.info.sla, VIP==1)

df_biodiversity.effects <- merge(plot.info.sla.vip[,c(2:8,10:11)], Mean.Biomass.Tripartite, by=c("Plot_Nr") )
summary(df_biodiversity.effects)
df_biodiversity.effects.1 <- subset(df_biodiversity.effects, SD!=20)
df_biodiversity.effects.1$treatment[df_biodiversity.effects.1$N==0&df_biodiversity.effects.1$F==0] <- "control"
df_biodiversity.effects.1$treatment[df_biodiversity.effects.1$N==1&df_biodiversity.effects.1$F==0] <- "nitrogen"
df_biodiversity.effects.1$treatment[df_biodiversity.effects.1$N==0&df_biodiversity.effects.1$F==1] <- "fungicide"
df_biodiversity.effects.1$treatment[df_biodiversity.effects.1$N==1&df_biodiversity.effects.1$F==1] <- "NF"

df_biodiversity.effects.1$Month <- 1*(df_biodiversity.effects.1$Month=="August")
Scalevars.BEF <- df_biodiversity.effects.1[,c("Month", "N", "F", "SD",  "sownSLA", "sownMPD", "NE", "CE", "SE", "SE_inter", "SE_intra")]
# do the scaling
Scalevars.BEF.2 <- scale(Scalevars.BEF)

corrs.BEF<-cor(Scalevars.BEF.2, use = "pairwise.complete.obs")
corrplot(corrs.BEF, order="AOE", method = 'number')
##correlation between sownMPD and sownSLA are relatively high (0.65)
cat.vars.BEF <- df_biodiversity.effects.1[,c("Plot_Nr", "Block", "composition","treatment")]
data.scaled.BEF <- cbind(cat.vars.BEF, Scalevars.BEF.2)
head(data.scaled.BEF)

hist(data.scaled.BEF$NE)
m1_NE<-lmer(NE ~ (N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_NE)#Month (-0.27***)， N:SD:0.13 (p=0.06)
Anova(m1_NE, test="F")
qqnorm(resid(m1_NE))
qqline(resid(m1_NE))
check_overdispersion(m1_NE)#No overdispersion
check_normality(m1_NE)#

##
hist(data.scaled.BEF$CE)
m1_CE<-lmer(CE ~ (N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF )
summary(m1_CE)#SD: 0.14*, N:sownMPD (-0.21*), N:SD (0.16**), Month (-0.40***)
Anova(m1_CE, test="F")
emm_CE_N_SD <- emmeans(m1_CE, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_CE_N_within_SD <- contrast(emm_CE_N_SD, method = "pairwise", adjust = "tukey")
contr_CE_N_within_SD 

plot_model(m1_CE, type="pred", terms=c("sownMPD", "N"),show.data=F)#same as the results from Null model before
plot_model(m1_CE, type="pred", terms=c("SD", "N"),show.data=F)#same as the results from Null model before
qqnorm(resid(m1_CE))
qqline(resid(m1_CE))
check_overdispersion(m1_CE)#No overdispersion
check_normality(m1_CE)#non-normality

##########
m1_SEAS<-lmer(SE_inter ~ (N + F + SD + sownSLA + sownMPD)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_SEAS)#N (0.34***), month (-0.11*)
Anova(m1_SEAS, test="F")
qqnorm(resid(m1_SEAS))
qqline(resid(m1_SEAS))
check_overdispersion(m1_SEAS)#No overdispersion
check_normality(m1_SEAS)#normality
emm_SEAS_N_SD <- emmeans(m1_SEAS, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_SEAS_N_within_SD <- contrast(emm_SEAS_N_SD, method = "pairwise", adjust = "tukey")
contr_SEAS_N_within_SD

hist(data.scaled.BEF$SE_intra)
m1_SEPS<-lmer(SE_intra ~ (N + F + SD  + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_SEPS)#N:sownMPD, N:SD
Anova(m1_SEPS, test="F")
emm_SEPS_N_SD <- emmeans(m1_SEPS, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_SEPS_N_within_SD <- contrast(emm_SEPS_N_SD, method = "pairwise", adjust = "tukey")
contr_SEPS_N_within_SD

plot_model(m1_SEPS, type="pred", terms=c("SD", "N"),show.data=F)#same as the results from Null model before
plot_model(m1_SEPS, type="pred", terms=c("sownMPD", "N"),show.data=F)#same as the results from Null model before

qqnorm(resid(m1_SEPS))
qqline(resid(m1_SEPS))
check_overdispersion(m1_SEPS)#No overdispersion
check_normality(m1_SEPS)#non-normality
##Fig.2---------------------------------------------------
y_all <- c(df_biodiversity.effects.1$NE,
           df_biodiversity.effects.1$CE,
           df_biodiversity.effects.1$SE_inter,
           df_biodiversity.effects.1$SE_intra)
y_all <- y_all[!is.na(y_all)]

y_max <- max(abs(y_all))
y_min <- min(y_all)
y_limits <- c(y_min, y_max)
y_breaks <- c(-100, -75, -50, -25, 0, 30, 60, 90, 120)

NE.sown.SD.fig <- ggplot(
  data = df_biodiversity.effects.1,
  aes(x = factor(SD), y = NE, fill = factor(N), color = factor(N))
) +
  geom_violin(aes(group = interaction(SD, N)), trim = FALSE, alpha = 0.7,
              linewidth = 0.4, position = position_dodge(width = 0.6)) +
  #geom_sina(aes(group = interaction(SD, N)), size = 0.5, alpha = 0.8,
  #maxwidth = 0.4, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("0" = "#99CCFF", "1" = "#336699")) +
  scale_color_manual(values = c("0" = "#99CCFF", "1" = "#336699"),
                     labels = c("0" = "Without nitrogen", "1" = "With nitrogen")) +
  scale_y_continuous(breaks = y_breaks) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "plain"),
    strip.background = element_blank(),
    axis.title = element_text(size = 11, face = "plain"),
    axis.text = element_text(size = 10, face = "plain"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "", y = "Net effect", fill = "Nitrogen")

CE.sown.SD.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(SD), y = CE, fill = factor(N), color = factor(N))) +
  geom_violin(aes(group = interaction(SD, N)), trim = FALSE, alpha = 0.7, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  #geom_sina(aes(group = interaction(SD, N)), size = 0.5, alpha = 0.8, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "#99CCFF", "1" = "#336699")) +
  scale_color_manual(values = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen")) +
  #facet_wrap(~ Month) +
  theme_minimal() + scale_y_continuous(breaks = y_breaks) +
  coord_cartesian(ylim = y_limits) +
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "", y = "Complementarity effects", fill = "Nitrogen")
y_pos_SD4_CE <- max(df_biodiversity.effects.1$CE[df_biodiversity.effects.1$SD == 4], na.rm = TRUE) + 
  0.05 * diff(range(df_biodiversity.effects.1$CE))
CE.sown.SD.fig <- CE.sown.SD.fig + annotate("text", x = "4", y = y_pos_SD4_CE-60, label = "*", size = 5)

SEAS.sown.SD.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(SD), y = SE_inter, fill = factor(N), color = factor(N))) +
  geom_violin(aes(group = interaction(SD, N)), trim = FALSE, alpha = 0.7, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  #geom_sina(aes(group = interaction(SD, N)), size = 0.5, alpha = 0.8, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "#99CCFF", "1" = "#336699")) +
  scale_color_manual(values = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen")) +
  #facet_wrap(~ Month) +
  theme_minimal() + scale_y_continuous(breaks = y_breaks) +
  coord_cartesian(ylim = y_limits) +
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 10, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "Species richness", y = "Interspecific selection effects", fill = "Nitrogen")
y_pos_SD8_SEAS <- max(df_biodiversity.effects.1$SE_inter[df_biodiversity.effects.1$SD == 8], na.rm = TRUE) + 
  0.05 * diff(range(df_biodiversity.effects.1$SE_inter))
y_pos_SD4_SEAS <- max(df_biodiversity.effects.1$SE_inter[df_biodiversity.effects.1$SD == 4], na.rm = TRUE) + 
  0.05 * diff(range(df_biodiversity.effects.1$SE_inter))
SEAS.sown.SD.fig <- SEAS.sown.SD.fig + annotate("text", x = "8", y = y_pos_SD8_SEAS, label = "***", size = 5) + annotate("text", x = "4", y = y_pos_SD4_SEAS, label = "***", size = 5)

SEPS.sown.SD.fig <- ggplot(
  data = df_biodiversity.effects.1,
  aes(x = factor(SD), y = SE_intra, fill = factor(N), color = factor(N))
) + geom_violin(aes(group = interaction(SD, N)), trim = FALSE, alpha = 0.7, linewidth = 0.4,
    position = position_dodge(width = 0.6)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen"),
    name = "Treatment") +
  scale_color_manual(
    values = c("0" = "#99CCFF", "1" = "#336699"),
    labels = c("0" = "Without nitrogen", "1" = "With nitrogen"),
    name = "Treatment"
  ) + theme_minimal() +
  scale_y_continuous(breaks = y_breaks) +
  coord_cartesian(ylim = y_limits) +
  theme(
    strip.text = element_text(size = 11, face = "plain"),
    strip.background = element_blank(),
    axis.title = element_text(size = 11, face = "plain"),
    axis.text = element_text(size = 10, face = "plain"),
    panel.grid = element_blank(),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    x = "Species richness",
    y = "Intraspecific selection effects"
  )

y_pos_SD8_SEPS <- max(df_biodiversity.effects.1$SE_intra[df_biodiversity.effects.1$SD == 8], na.rm = TRUE) + 
  0.05 * diff(range(df_biodiversity.effects.1$SE_intra))
SEPS.sown.SD.fig <- SEPS.sown.SD.fig + annotate("text", x = "8", y = y_pos_SD8_SEPS, label = "*", size = 5) 

data.scaled.BEF$N_raw <- ifelse(data.scaled.BEF$N < 0, 0, 1)
data.scaled.BEF$N_f <- factor(data.scaled.BEF$N_raw)#N as a factor
m1_CE2 <- lmer(CE ~ (N_f + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N_f:F:SD +
                 (1 | Block) + (1 | composition) + (1 | Plot_Nr), data = data.scaled.BEF, control = lmerControl(optimizer = "bobyqa") )
vis1 <- visreg(m1_CE2, xvar = "sownMPD", by = "N_f", partial = TRUE, plot = FALSE)
names(vis1)
partial_resid <- vis1$res
colnames(partial_resid)[11] <- c("CE_partial")
table(partial_resid$N_f)
fit1 <- vis1$fit
fit_long <- data.frame(sownMPD = fit1$sownMPD, CE_fitted = fit1$visregFit,
                       lower = fit1$visregLwr, upper = fit1$visregUpr, N_f = fit1$N_f)
##add the unscaled sownMPD to partial_resid
mpd_mean <- mean(df_biodiversity.effects.1$sownMPD, na.rm = TRUE)
mpd_sd   <- sd(df_biodiversity.effects.1$sownMPD, na.rm = TRUE)
partial_resid$sownMPD_raw <- partial_resid$sownMPD * mpd_sd + mpd_mean
fit_long$sownMPD_raw <- fit_long$sownMPD * mpd_sd + mpd_mean
CE.sown.MPD.visreg.fig <- ggplot() + 
  # partial residuals (raw MPD)
  geom_point(data = partial_resid, aes(x = sownMPD_raw, y = CE_partial, color = N_f),
             alpha = 1, size = 1) + 
  # CI ribbon (raw MPD)
  geom_ribbon(data = fit_long, aes(x = sownMPD_raw, ymin = lower, ymax = upper, fill = N_f),
              alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_line(data = fit_long, aes(x = sownMPD_raw, y = CE_fitted, color = N_f),
            linewidth = 0.9) + 
  scale_color_manual(values = c("0" = "#99CCFF", "1" = "#336699")) +
  scale_fill_manual(values  = c("0" = "#99CCFF", "1" = "#336699")) +
  labs(x = "Functional diversity (sown MPD)",
       y = "Complementarity effects\n(partial residuals)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11)
  )

m1_SEPS2 <- lmer(SE_intra ~ (N_f + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N_f:F:SD +
                   (1 | Block) + (1 | composition) + (1 | Plot_Nr), data = data.scaled.BEF, control = lmerControl(optimizer = "bobyqa") )
#partial plots
vis2 <- visreg(m1_SEPS2, xvar = "sownMPD", by = "N_f", partial = TRUE, plot = FALSE)
names(vis2)
partial_resid2 <- vis2$res
colnames(partial_resid2)[11] <- c("SEPS_partial")
table(partial_resid2$N_f)
fit2 <- vis2$fit
fit_long2 <- data.frame(sownMPD = fit2$sownMPD, CE_fitted = fit2$visregFit,
                        lower = fit2$visregLwr, upper = fit2$visregUpr, N_f = fit2$N_f)
##add the unscaled sownMPD to partial_resid
partial_resid2$sownMPD_raw <- partial_resid2$sownMPD * mpd_sd + mpd_mean
fit_long2$sownMPD_raw <- fit_long2$sownMPD * mpd_sd + mpd_mean
SEPS.sown.MPD.visreg.fig <- ggplot() + 
  # partial residuals (raw MPD)
  geom_point(data = partial_resid2, aes(x = sownMPD_raw, y = SEPS_partial, color = N_f),
             alpha = 1, size = 1) + 
  # CI ribbon (raw MPD)
  geom_ribbon(data = fit_long2, aes(x = sownMPD_raw, ymin = lower, ymax = upper, fill = N_f),
              alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_line(data = fit_long2, aes(x = sownMPD_raw, y = CE_fitted, color = N_f),
            linewidth = 0.9) + 
    scale_color_manual(values = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen")) +
  scale_fill_manual(values  = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen")) +
  labs(x = "Functional diversity (sown MPD)",
       y = "Intraspecific selection effects\n(partial residuals)", fill= "Treatment", color="Treatment") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11)
  )

combined_mean_BEFs_N_interactions <- (NE.sown.SD.fig + CE.sown.SD.fig)/(SEAS.sown.SD.fig + SEPS.sown.SD.fig) / (CE.sown.MPD.visreg.fig + SEPS.sown.MPD.visreg.fig) + plot_layout(heights = c(1, 1, 1)) + plot_annotation(tag_levels = "A")
ggsave("figures/Fig2.pdf", combined_mean_BEFs_N_interactions, width = 7, height = 8.5, dpi = 300)

##Fig.S1---------------------------------------------------
NE.season.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(Month), y = NE, fill = factor(Month), color = factor(Month))) +
  geom_violin(aes(group = factor(Month)), trim = FALSE, alpha = 0.6, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_sina(aes(group = factor(Month)), size = 0.5, alpha = 0.6, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "white", "1" = "white")) +
  scale_color_manual(values = c("0" = "#009E73", "1" = "#E69F00"), labels = c("0" = "June", "1" = "August")) +
  scale_x_discrete(labels = c("0" = "June", "1" = "August")) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 11, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "", y = "Net effect", fill = "") 
y_pos_season_NE <- max(df_biodiversity.effects.1$NE, na.rm = TRUE) + 0.05 * diff(range(df_biodiversity.effects.1$NE, na.rm = TRUE))
NE.season.fig <- NE.season.fig + annotate("text", x = 1.5, y = y_pos_season_NE, label = "***", size = 5)

CE.season.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(Month), y = CE, fill = factor(Month), color = factor(Month))) +
  geom_violin(aes(group = factor(Month)), trim = FALSE, alpha = 0.6, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_sina(aes(group = factor(Month)), size = 0.5, alpha = 0.6, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "white", "1" = "white")) +
  scale_color_manual(values = c("0" = "#009E73", "1" = "#E69F00"), labels = c("0" = "June", "1" = "August")) +
  scale_x_discrete(labels = c("0" = "June", "1" = "August")) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 11, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "", y = "Complementarity effects", fill = "") 
y_pos_season_CE <- max(df_biodiversity.effects.1$CE, na.rm = TRUE) + 0.05 * diff(range(df_biodiversity.effects.1$CE, na.rm = TRUE))
CE.season.fig <- CE.season.fig + annotate("text", x = 1.5, y = y_pos_season_CE, label = "***", size = 5)

SEAS.season.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(Month), y = SE_inter, fill = factor(Month), color = factor(Month))) +
  geom_violin(aes(group = factor(Month)), trim = FALSE, alpha = 0.6, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_sina(aes(group = factor(Month)), size = 0.5, alpha = 0.6, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "white", "1" = "white")) +
  scale_color_manual(values = c("0" = "#009E73", "1" = "#E69F00"), labels = c("0" = "June", "1" = "August")) +
  scale_x_discrete(labels = c("0" = "June", "1" = "August")) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 11, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "Season", y = "Interspecific selection effects", fill = "") 
y_pos_season_SEAS <- max(df_biodiversity.effects.1$SE_inter, na.rm = TRUE) + 0.05 * diff(range(df_biodiversity.effects.1$SE_inter, na.rm = TRUE))
SEAS.season.fig <- SEAS.season.fig + annotate("text", x = 1.5, y = y_pos_season_SEAS, label = "*", size = 5)

SEPS.season.fig <- ggplot(data = df_biodiversity.effects.1, aes(x = factor(Month), y = SE_intra, fill = factor(Month), color = factor(Month))) +
  geom_violin(aes(group = factor(Month)), trim = FALSE, alpha = 0.6, linewidth = 0.4, position = position_dodge(width = 0.6) ) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_sina(aes(group = factor(Month)), size = 0.5, alpha = 0.6, maxwidth = 0.4, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), shape = 16, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "white", "1" = "white")) +
  scale_color_manual(values = c("0" = "#009E73", "1" = "#E69F00"), labels = c("0" = "June", "1" = "August")) +
  scale_x_discrete(labels = c("0" = "June", "1" = "August")) +
  theme_minimal() + 
  theme(strip.text = element_text(size = 11, face = "plain"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11, face = "plain"),
        axis.text = element_text(size = 11, face = "plain"),
        panel.grid = element_blank(),
        legend.position = "",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + labs(x = "Season", y = "Intraspecific selection effects", fill = "") 
y_pos_season_SEPS <- max(df_biodiversity.effects.1$SE_intra, na.rm = TRUE) + 0.05 * diff(range(df_biodiversity.effects.1$SE_intra, na.rm = TRUE))
SEPS.season.fig <- SEPS.season.fig + annotate("text", x = 1.5, y = y_pos_season_SEPS, label = "***", size = 5)
#shows variation in biodiversity effects between two harvest seasons
Fig.S1_combined <- (NE.season.fig + CE.season.fig)/(SEAS.season.fig + SEPS.season.fig) + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
ggsave("figures/Fig.S1.pdf", Fig.S1_combined, width = 7, height = 7, dpi = 300)

##Fig.S2---------------------------------------------------
##get coefficient of NE model
model_m1_NE_Summary <- as.data.frame(coef(summary(m1_NE)))
model_m1_NE_CI <- confint(m1_NE)[-c(1,2),]
model_m1_NE_Result <- as.data.frame(cbind(model_m1_NE_Summary, model_m1_NE_CI[3:19,]))
model_m1_NE_Result$Variable <- row.names(model_m1_NE_Result)
row.names(model_m1_NE_Result)
model_m1_NE_Result

model_m1_NE_Result <- model_m1_NE_Result %>%
  dplyr::select(Variable, Estimate, lower = '2.5 %', upper = '97.5 %', P = 'Pr(>|t|)') %>%
  filter(Variable != '(Intercept)') %>%
  mutate(Y = 'NE') %>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0))) %>%
  mutate(Variable = plyr::mapvalues(Variable, 
                                    from = c("N", "F", "SD", "sownSLA", "sownMPD", "Month", "N:F", "N:SD", "N:sownSLA",
                                             "N:sownMPD", "F:SD",  "F:sownSLA", "F:sownMPD", "SD:sownSLA","SD:sownMPD", "N:F:SD"
                                    ),
                                    to = c("Nitrogen", "Fungicide", "Sown richness", "sownSLA", "sownMPD", "Season", 
                                           "Nitrogen:Fungicide", "Nitrogen:SR", "Nitrogen:sownSLA", "Nitrogen:sownMPD", "Fungicide:SR", "Fungicide:sownSLA",
                                           "Fungicide:sownMPD", "SR:sownSLA", "SR:sownMPD","N:F:SR"))) %>%  mutate(Colour = as.factor(Colour))

interaction_order <- c("N:F:SR","SR:sownSLA", "SR:sownMPD", "Nitrogen:Fungicide" ,"Fungicide:sownSLA", "Fungicide:sownMPD","Fungicide:SR", "Nitrogen:sownSLA","Nitrogen:sownMPD", "Nitrogen:SR", "Season", "sownSLA", "sownMPD", "Fungicide", "Nitrogen", "Sown richness")

model_m1_NE_Result$Variable <- factor(model_m1_NE_Result$Variable, levels = interaction_order)
#plot effect size
# Create the effect size plot
Plot.model_m1_NE_Result <- ggplot(data = model_m1_NE_Result, aes(x = Variable, y = Estimate, fill = Colour)) +
  # Bars: Only significant effects get bars, non-significant remain empty
  geom_bar(stat = "identity", position = "dodge", width = 0.6, aes(alpha = as.factor(Colour))) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_NE_Result %>% filter(Colour != 0), 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "gray", linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'transparent',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
  labs(x = "", y = "Effect size", title = "Net diversity effect") +
  theme(
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2, 0.4, 0.2, 0), "cm"),
    legend.position = "none",
    strip.text = element_text(size = 11),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold")
  ) +  theme(panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"))+
  coord_flip()  # Flip coordinates for better readability
Plot.model_m1_NE_Result

#FigS2.B
#get coefficient of CE model
model_m1_CE_Summary <- as.data.frame(coef(summary(m1_CE)))
model_m1_CE_CI <- confint(m1_CE)[-c(1,2),]
model_m1_CE_Result <- as.data.frame(cbind(model_m1_CE_Summary, model_m1_CE_CI[3:19,]))
model_m1_CE_Result$Variable <- row.names(model_m1_CE_Result)
row.names(model_m1_CE_Result)
model_m1_CE_Result

model_m1_CE_Result <- model_m1_CE_Result %>%
  dplyr::select(Variable, Estimate, lower = '2.5 %', upper = '97.5 %', P = 'Pr(>|t|)') %>%
  filter(Variable != '(Intercept)') %>%
  mutate(Y = 'CE') %>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0))) %>%
  mutate(Variable = plyr::mapvalues(Variable, 
                                    from = c("N", "F", "SD", "sownSLA", "sownMPD", "Month", "N:F", "N:SD", "N:sownSLA",
                                             "N:sownMPD", "F:SD",  "F:sownSLA", "F:sownMPD", "SD:sownSLA","SD:sownMPD", "N:F:SD"
                                    ),
                                    to = c("Nitrogen", "Fungicide", "Sown richness", "sownSLA", "sownMPD", "Season", 
                                           "Nitrogen:Fungicide", "Nitrogen:SR", "Nitrogen:sownSLA", "Nitrogen:sownMPD", "Fungicide:SR", "Fungicide:sownSLA",
                                           "Fungicide:sownMPD", "SR:sownSLA", "SR:sownMPD","N:F:SR"))) %>%  mutate(Colour = as.factor(Colour))

model_m1_CE_Result$Variable <- factor(model_m1_CE_Result$Variable, levels = interaction_order)
#plot effect size
# Create the effect size plot
Plot.model_m1_CE_Result <- ggplot(data = model_m1_CE_Result, aes(x = Variable, y = Estimate, fill = Colour)) +
  # Bars: Only significant effects get bars, non-significant remain empty
  geom_bar(stat = "identity", position = "dodge", width = 0.6, aes(alpha = as.factor(Colour))) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_CE_Result %>% filter(Colour != 0), 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "gray", linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'transparent',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
  labs(x = "", y = "Effect size", title = "Complementarity effects") +
  theme(
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2, 0.4, 0.2, 0), "cm"),
    legend.position = "none",
    strip.text = element_text(size = 11),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold")
  ) +  theme(panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"))+
  coord_flip()  # Flip coordinates for better readability
Plot.model_m1_CE_Result

##get coefficient of SEAS model (inter-specific selection)
model_m1_SEAS_Summary <- as.data.frame(coef(summary(m1_SEAS)))
model_m1_SEAS_CI <- confint(m1_SEAS)[-c(1,2),]
model_m1_SEAS_Result <- as.data.frame(cbind(model_m1_SEAS_Summary, model_m1_SEAS_CI[3:19,]))
model_m1_SEAS_Result$Variable <- row.names(model_m1_SEAS_Result)
row.names(model_m1_SEAS_Result)
model_m1_SEAS_Result

model_m1_SEAS_Result <- model_m1_SEAS_Result %>%
  dplyr::select(Variable, Estimate, lower = '2.5 %', upper = '97.5 %', P = 'Pr(>|t|)') %>%
  filter(Variable != '(Intercept)') %>%
  mutate(Y = 'SEAS') %>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0))) %>%
  mutate(Variable = plyr::mapvalues(Variable, 
                                    from = c("N", "F", "SD", "sownSLA", "sownMPD", "Month", "N:F", "N:SD", "N:sownSLA",
                                             "N:sownMPD", "F:SD",  "F:sownSLA", "F:sownMPD", "SD:sownSLA","SD:sownMPD", "N:F:SD"
                                    ),
                                    to = c("Nitrogen", "Fungicide", "Sown richness", "sownSLA", "sownMPD", "Season", 
                                           "Nitrogen:Fungicide", "Nitrogen:SR", "Nitrogen:sownSLA", "Nitrogen:sownMPD", "Fungicide:SR", "Fungicide:sownSLA",
                                           "Fungicide:sownMPD", "SR:sownSLA", "SR:sownMPD","N:F:SR"))) %>%  mutate(Colour = as.factor(Colour))

model_m1_SEAS_Result$Variable <- factor(model_m1_SEAS_Result$Variable, levels = interaction_order)
#plot effect size
# Create the effect size plot
Plot.model_m1_SEAS_Result <- ggplot(data = model_m1_SEAS_Result, aes(x = Variable, y = Estimate, fill = Colour)) +
  # Bars: Only significant effects get bars, non-significant remain empty
  geom_bar(stat = "identity", position = "dodge", width = 0.6, aes(alpha = as.factor(Colour))) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_SEAS_Result %>% filter(Colour != 0), 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "gray", linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'transparent',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
  labs(x = "", y = "Effect size", title = "Interspecific selection effects") +
  theme(
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2, 0.4, 0.2, 0), "cm"),
    legend.position = "none",
    strip.text = element_text(size = 11),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold")
  ) +  theme(panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"))+
  coord_flip()  # Flip coordinates for better readability
Plot.model_m1_SEAS_Result

##get coefficient of SEPS model (intra-specific selection)
model_m1_SEPS_Summary <- as.data.frame(coef(summary(m1_SEPS)))
model_m1_SEPS_CI <- confint(m1_SEPS)[-c(1,2),]
model_m1_SEPS_Result <- as.data.frame(cbind(model_m1_SEPS_Summary, model_m1_SEPS_CI[3:19,]))
model_m1_SEPS_Result$Variable <- row.names(model_m1_SEPS_Result)
row.names(model_m1_SEPS_Result)
model_m1_SEPS_Result

model_m1_SEPS_Result <- model_m1_SEPS_Result %>%
  dplyr::select(Variable, Estimate, lower = '2.5 %', upper = '97.5 %', P = 'Pr(>|t|)') %>%
  filter(Variable != '(Intercept)') %>%
  mutate(Y = 'SEPS') %>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0))) %>%
  mutate(Variable = plyr::mapvalues(Variable, 
                                    from = c("N", "F", "SD", "sownSLA", "sownMPD", "Month", "N:F", "N:SD", "N:sownSLA",
                                             "N:sownMPD", "F:SD",  "F:sownSLA", "F:sownMPD", "SD:sownSLA","SD:sownMPD", "N:F:SD"
                                    ),
                                    to = c("Nitrogen", "Fungicide", "Sown richness", "sownSLA", "sownMPD", "Season", 
                                           "Nitrogen:Fungicide", "Nitrogen:SR", "Nitrogen:sownSLA", "Nitrogen:sownMPD", "Fungicide:SR", "Fungicide:sownSLA",
                                           "Fungicide:sownMPD", "SR:sownSLA", "SR:sownMPD","N:F:SR"))) %>%  mutate(Colour = as.factor(Colour))

model_m1_SEPS_Result$Variable <- factor(model_m1_SEPS_Result$Variable, levels = interaction_order)
#plot effect size
# Create the effect size plot
Plot.model_m1_SEPS_Result <- ggplot(data = model_m1_SEPS_Result, aes(x = Variable, y = Estimate, fill = Colour)) +
  # Bars: Only significant effects get bars, non-significant remain empty
  geom_bar(stat = "identity", position = "dodge", width = 0.6, aes(alpha = as.factor(Colour))) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_SEPS_Result %>% filter(Colour != 0), 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "gray", linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'transparent',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
  labs(x = "", y = "Effect size", title = "Intraspecific selection effects") +
  theme(
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2, 0.4, 0.2, 0), "cm"),
    legend.position = "none",
    strip.text = element_text(size = 11),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold")
  ) +  theme(panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"))+
  coord_flip()  # Flip coordinates for better readability
Plot.model_m1_SEPS_Result

combined_mean_BEFs_effect_sizes <- (Plot.model_m1_NE_Result + Plot.model_m1_CE_Result)/(Plot.model_m1_SEAS_Result + Plot.model_m1_SEPS_Result) + plot_annotation(tag_levels = "A")

##Fig. S2 in the manuscript
ggsave("figures/Fig.S2.pdf", combined_mean_BEFs_effect_sizes, width = 8, height = 8, dpi = 300)

