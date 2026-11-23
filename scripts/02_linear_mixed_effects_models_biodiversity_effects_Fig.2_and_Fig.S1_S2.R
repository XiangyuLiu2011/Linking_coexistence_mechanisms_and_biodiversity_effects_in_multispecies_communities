##scripts/02_linear_mixed_effects_models_biodiversity_effects
#should first run script.1, as AP_Biomass is needed here.
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

##signed square-root transformation to improve residual normality and homogeneity of variance
df_biodiversity.effects.1 <- df_biodiversity.effects.1 %>%
  dplyr::mutate(
    CE_sqrt = sign(CE) * sqrt(abs(CE)),
    SE_intra_sqrt = sign(SE_intra) * sqrt(abs(SE_intra))
    )

df_biodiversity.effects.1$Month <- 1*(df_biodiversity.effects.1$Month=="August")
Scalevars.BEF <- df_biodiversity.effects.1[,c("Month", "N", "F", "SD",  "sownSLA", "sownMPD", "NE", "CE", "SE", "SE_inter", "SE_intra")]
# do the scaling
Scalevars.BEF.2 <- scale(Scalevars.BEF)

corrs.BEF<-cor(Scalevars.BEF.2, use = "pairwise.complete.obs")
corrplot(corrs.BEF, order="AOE", method = 'number')
##correlation between sownMPD and sownSLA are relatively high (0.65)
cat.vars.BEF <- df_biodiversity.effects.1[,c("Plot_Nr", "Block", "composition", "treatment", "CE_sqrt", "SE_intra_sqrt")]
data.scaled.BEF <- cbind(cat.vars.BEF, Scalevars.BEF.2)
head(data.scaled.BEF)

hist(scale(data.scaled.BEF$NE))
m1_NE<-lmer(NE ~ (N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_NE)#Month (-0.26***; named as season in the manuscript), N:SD:0.13 (p=0.06)
Anova(m1_NE, test="F")
qqnorm(resid(m1_NE))
qqline(resid(m1_NE))
check_overdispersion(m1_NE)#No overdispersion
check_normality(m1_NE)#residuals appear as normally distributed

##
hist(data.scaled.BEF$CE_sqrt)
hist(data.scaled.BEF$CE)
m1_CE<-lmer(CE_sqrt ~ (N + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF )
summary(m1_CE)#SD: 0.40 (p=0.06), N:sownMPD (-0.52*), N:SD (0.40*), Month (-0.64***)
Anova(m1_CE, test="F")
emm_CE_N_SD <- emmeans(m1_CE, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_CE_N_within_SD <- contrast(emm_CE_N_SD, method = "pairwise", adjust = "tukey")
contr_CE_N_within_SD#N weakened complementary in 4-species plots

plot_model(m1_CE, type="pred", terms=c("sownMPD", "N"),show.data=F)#same as the results from Null model before
plot_model(m1_CE, type="pred", terms=c("SD", "N"),show.data=F)#same as the results from Null model before
qqnorm(resid(m1_CE))
qqline(resid(m1_CE))
check_overdispersion(m1_CE)#No overdispersion
check_normality(m1_CE)#residuals were not normally distributed, but this was the best transformation

##########
hist(data.scaled.BEF$SE_inter)
m1_SEAS<-lmer(SE_inter ~ (N + F + SD + sownSLA + sownMPD)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_SEAS)#N (0.33***), month (-0.11, p=0.05)
Anova(m1_SEAS, test="F")
qqnorm(resid(m1_SEAS))
qqline(resid(m1_SEAS))
check_overdispersion(m1_SEAS)#No overdispersion
check_normality(m1_SEAS)#residuals were normally distributed
emm_SEAS_N_SD <- emmeans(m1_SEAS, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_SEAS_N_within_SD <- contrast(emm_SEAS_N_SD, method = "pairwise", adjust = "tukey")
contr_SEAS_N_within_SD

hist(data.scaled.BEF$SE_intra_sqrt)
m1_SEPS<-lmer(SE_intra_sqrt ~ (N + F + SD  + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N:F:SD + (1 | Block) + (1 | composition) + (1 | Plot_Nr) , control = lmerControl(optimizer = "bobyqa"), data = data.scaled.BEF)
summary(m1_SEPS)#SD: -0.37**, Month:0.56***, N:SD: -0.20* 
Anova(m1_SEPS, test="F")
plot_model(m1_SEPS, type="pred", terms=c("SD", "N"),show.data=F)#same as the results from Null model before

emm_SEPS_N_SD <- emmeans(m1_SEPS, ~ N | SD)
# pairwise contrasts: effect of N at each SD
contr_SEPS_N_within_SD <- contrast(emm_SEPS_N_SD, method = "pairwise", adjust = "tukey")
contr_SEPS_N_within_SD

qqnorm(resid(m1_SEPS))
qqline(resid(m1_SEPS))
check_overdispersion(m1_SEPS)#No overdispersion
check_normality(m1_SEPS)#residuals were not normally distributed, but this was the best transformation

##Fig.2--for the major revision-------------------------------------------------
richness_levels <- sort(unique(df_biodiversity.effects.1$SD[!is.na(df_biodiversity.effects.1$SD)]))

BEF_plot_data <- df_biodiversity.effects.1 %>% filter(!is.na(SD), !is.na(N)) %>%
  mutate(SD_group = factor(SD, levels = richness_levels),
    N_group = factor(N, levels = c(0, 1)),
    x_center = as.numeric(SD_group) +
      if_else(N_group == "0", -0.22, 0.22),
    x_raw = x_center - 0.075,
    x_summary = x_center + 0.075
  )

treatment_colors <- c("0" = "#99CCFF", "1" = "#336699")
treatment_labels <- c("0" = "Without nitrogen", "1" = "With nitrogen")

BEF_theme <- theme_minimal() + theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "none")

# 1. Net effect ------------------------------------------------------------
NE_plot_data <- BEF_plot_data %>% filter(is.finite(NE))

NE_summary <- NE_plot_data %>% group_by(SD_group, N_group) %>%
  summarise(
    n = n(),
    mean_NE = mean(NE),
    sd_NE = sd(NE),
    se_NE = sd_NE / sqrt(n),
    x_summary = first(x_summary),
    .groups = "drop"
  )

NE.sown.SD.fig <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed",color = "black", linewidth = 0.4) +
  geom_point(data = NE_plot_data,
    aes(x = x_raw, y = NE, fill = N_group),
    position = position_jitter(width = 0.045, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    color = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) + geom_errorbar(
    data = NE_summary,
    aes(x = x_summary, ymin = mean_NE - sd_NE, ymax = mean_NE + sd_NE),
    width = 0.07, linewidth = 0.6, color = "black"
  ) + geom_point(data = NE_summary, aes(x = x_summary, y = mean_NE, fill = N_group),
    shape = 23, size = 3.8, stroke = 0.6, color = "black"
  ) + scale_x_continuous(breaks = seq_along(richness_levels),
    labels = richness_levels, expand = expansion(add = 0.45)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment",
    drop = FALSE
  ) + BEF_theme + labs(x = NULL, y = "Net effect")
NE.sown.SD.fig

# 2. Complementarity effects ----------------------------------------------
CE_plot_data <- BEF_plot_data %>% filter(is.finite(CE))

CE_summary <- CE_plot_data %>% group_by(SD_group, N_group) %>%
  summarise(
    n = n(),
    mean_CE = mean(CE),
    sd_CE = sd(CE),
    se_CE = sd_CE / sqrt(n),
    x_summary = first(x_summary),
    .groups = "drop"
  )

CE.sown.SD.fig <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_point(data = CE_plot_data, aes(x = x_raw, y = CE, fill = N_group),
    position = position_jitter(width = 0.045, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    color = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_errorbar(data = CE_summary,
    aes(x = x_summary, ymin = mean_CE - sd_CE, ymax = mean_CE + sd_CE),
    width = 0.07, linewidth = 0.6, color = "black"
  ) + geom_point(data = CE_summary, aes(x = x_summary, y = mean_CE, fill = N_group),
    shape = 23, size = 3.8, stroke = 0.6, color = "black"
  ) + scale_x_continuous(
    breaks = seq_along(richness_levels), labels = richness_levels, expand = expansion(add = 0.45)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment",
    drop = FALSE
  ) + BEF_theme + labs(x = NULL, y = "Complementarity effects")
CE.sown.SD.fig

# 3. Interspecific selection effects ---------------------------------------
SE_inter_plot_data <- BEF_plot_data %>% filter(is.finite(SE_inter))

SE_inter_summary <- SE_inter_plot_data %>% group_by(SD_group, N_group) %>%
  summarise(n = n(),
    mean_inter = mean(SE_inter),
    sd_inter = sd(SE_inter),
    se_inter = sd_inter / sqrt(n),
    x_summary = first(x_summary),
    .groups = "drop"
  )

SE_inter.sown.SD.fig <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_point(
    data = SE_inter_plot_data,
    aes(x = x_raw, y = SE_inter, fill = N_group),
    position = position_jitter(width = 0.045, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    color = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = SE_inter_summary,
    aes(x = x_summary,
      ymin = mean_inter - sd_inter,
      ymax = mean_inter + sd_inter
    ), width = 0.07, linewidth = 0.6, color = "black") +
  geom_point(
    data = SE_inter_summary,
    aes(x = x_summary, y = mean_inter, fill = N_group),
    shape = 23, size = 3.8, stroke = 0.6, color = "black"
  ) +
  scale_x_continuous(
    breaks = seq_along(richness_levels),
    labels = richness_levels,
    expand = expansion(add = 0.45)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  scale_fill_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment",
    drop = FALSE
  ) + BEF_theme + labs(x = "Sown species richness", y = "Interspecific selection effects")
SE_inter.sown.SD.fig

# 4. Intraspecific selection effects ---------------------------------------
SE_intra_plot_data <- BEF_plot_data %>% filter(is.finite(SE_intra))

SE_intra_summary <- SE_intra_plot_data %>%
  group_by(SD_group, N_group) %>%
  summarise(
    n = n(),
    mean_intra = mean(SE_intra),
    sd_intra = sd(SE_intra),
    se_intra = sd_intra / sqrt(n),
    x_summary = first(x_summary),
    .groups = "drop"
  )

SE_intra.sown.SD.fig <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4
  ) + geom_point(
    data = SE_intra_plot_data,
    aes(x = x_raw, y = SE_intra, fill = N_group),
    position = position_jitter(width = 0.045, height = 0, seed = 123),
    shape = 21, size = 2, stroke = 0.3,
    color = "grey40", alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = SE_intra_summary,
    aes(x = x_summary, 
      ymin = mean_intra - sd_intra,
      ymax = mean_intra + sd_intra
    ), width = 0.07, linewidth = 0.6, color = "black") +
  geom_point(
    data = SE_intra_summary,
    aes(x = x_summary, y = mean_intra, fill = N_group),
    shape = 23, size = 3.8, stroke = 0.6, color = "black"
  ) +
  scale_x_continuous(
    breaks = seq_along(richness_levels),
    labels = richness_levels,
    expand = expansion(add = 0.45)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  scale_fill_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment",
    drop = FALSE
  ) + theme_minimal() + BEF_theme +
  labs(x = "Sown species richness", y = "Intraspecific selection effects")
SE_intra.sown.SD.fig

# CE: nitrogen marginally reduced CE when richness at 4, but marginally increased it when richness at 8
CE_span <- diff(range(c(0, CE_plot_data$CE, CE_summary$mean_CE - CE_summary$sd_CE,
    CE_summary$mean_CE + CE_summary$sd_CE), na.rm = TRUE))

CE_sig <- CE_plot_data %>%
  filter(SD == 4) %>%
  group_by(SD_group) %>%
  summarise(top_raw = max(CE, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    CE_summary %>%
      group_by(SD_group) %>%
      summarise(
        top_bar = max(mean_CE + sd_CE, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "SD_group"
  ) %>%
  mutate(
    x = match(as.character(SD_group), as.character(richness_levels)),
    y = pmax(top_raw, top_bar) + 0.06 * CE_span,
    label = "*"
  )
CE.sown.SD.fig <- CE.sown.SD.fig + geom_text(data = CE_sig, aes(x = x, y = 0.5*y, label = label), size = 3.5)

# Interspecific selection: nitrogen comparisons at richness 4 and 8
inter_span <- diff(range(c(0, SE_inter_plot_data$SE_inter,
    SE_inter_summary$mean_inter - SE_inter_summary$sd_inter,
    SE_inter_summary$mean_inter + SE_inter_summary$sd_inter),
  na.rm = TRUE
))

inter_sig <- SE_inter_plot_data %>%
  filter(SD %in% c(4, 8)) %>%
  group_by(SD_group) %>%
  summarise(top_raw = max(SE_inter), .groups = "drop") %>%
  left_join(
    SE_inter_summary %>%
      group_by(SD_group) %>%
      summarise(
        top_bar = max(mean_inter + sd_inter, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "SD_group"
  ) %>%
  mutate(
    x = as.numeric(SD_group),
    y = pmax(top_raw, top_bar) + 0.06 * inter_span,
    label = "***"
  )

SE_inter.sown.SD.fig <- SE_inter.sown.SD.fig + geom_text(data = inter_sig, aes(x = x, y = y, label = label), size = 3.5)

# Combine and save ---------------------------------------------------------
combined_mean_BEFs_N_interactions <- ((NE.sown.SD.fig + CE.sown.SD.fig) / (SE_inter.sown.SD.fig + SE_intra.sown.SD.fig) +
    plot_layout(heights = c(1, 1), guides = "collect") + plot_annotation(tag_levels = "A")
) & theme(legend.position = "right", legend.justification = "center")

combined_mean_BEFs_N_interactions
ggsave("figures/Fig2_major_revision.pdf",combined_mean_BEFs_N_interactions, width = 8.5, height = 7, units = "in", dpi = 300)
ggsave("figures/Fig2_major_revision.tiff", combined_mean_BEFs_N_interactions, width = 8.5, height = 7, dpi = 300)

##Fig.S1---------------------------------------------
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
SEAS.season.fig <- SEAS.season.fig + annotate("text", x = 1.5, y = y_pos_season_SEAS, label = "†", size = 4)

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
ggsave("figures/Fig.S1_major_revison.pdf", Fig.S1_combined, width = 7, height = 7, dpi = 300)
ggsave("figures/Fig.S1_major_revison.tiff", Fig.S1_combined, width = 8, height = 8, dpi = 300)

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
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_NE_Result, aes(ymin = lower, ymax = upper), width = 0.3, linewidth = 0.5) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "black", linewidth = 0.5,  linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'gray',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  #scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
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
    plot.title = element_text(size = 11, hjust = 0.5, face = "plain")
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
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_CE_Result, 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "black", linewidth = 0.5, linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'gray',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  #scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
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
    plot.title = element_text(size = 11, hjust = 0.5, face = "plain")
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
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_SEAS_Result, 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "black", linewidth = 0.5, linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'gray',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  #scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
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
    plot.title = element_text(size = 11, hjust = 0.5, face = "plain")
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
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  # Error bars only for significant effects (filter out Colour == 0)
  geom_errorbar(data = model_m1_SEPS_Result, 
                aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  # Reference line at 0
  geom_hline(aes(yintercept = 0), colour = "black", linewidth = 0.5, linetype = "dashed") +
  theme_bw() +
  # Define fill colors for significant effects
  scale_fill_manual(name = NULL,
                    values = c('-1' = '#d95f02',  # Dark orange (negative)
                               '0' = 'gray',
                               '1' = '#1b9e77')   # Teal/blue-green (positive)
  ) +
  # Adjust transparency: Significant effects solid, non-significant fully transparent
  #scale_alpha_manual(values = c('1' = 1, '0' = 0, '-1' = 1), guide = "none") +
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
    plot.title = element_text(size = 11, hjust = 0.5, face = "plain")
  ) +  theme(panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"))+
  coord_flip()  # Flip coordinates for better readability
Plot.model_m1_SEPS_Result

combined_mean_BEFs_effect_sizes <- (Plot.model_m1_NE_Result + Plot.model_m1_CE_Result)/(Plot.model_m1_SEAS_Result + Plot.model_m1_SEPS_Result) + plot_annotation(tag_levels = "A")

##Fig. S2 in the manuscript
ggsave("figures/Fig.S2_major_revision.pdf", combined_mean_BEFs_effect_sizes, width = 8, height = 8, dpi = 300)
ggsave("figures/Fig.S2_major_revision.tiff", combined_mean_BEFs_effect_sizes, width = 8, height = 8, dpi = 300)

##--Fig.S3-------------
data.scaled.BEF$N_raw <- ifelse(data.scaled.BEF$N < 0, 0, 1)
data.scaled.BEF$N_f <- factor(data.scaled.BEF$N_raw)#N as a factor
m1_CE2 <- lmer(CE_sqrt ~ (N_f + F + SD + sownMPD + sownSLA)^2 + Month - sownSLA:sownMPD + N_f:F:SD +
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
  scale_fill_manual(values = c("0" = "#99CCFF", "1" = "#336699"), labels = c("0" = "Without nitrogen", "1" = "With nitrogen"),
                    name = "Treatment") +
  scale_color_manual(
    values = c("0" = "#99CCFF", "1" = "#336699"),
    labels = c("0" = "Without nitrogen", "1" = "With nitrogen"),
    name = "Treatment"
  ) + theme_minimal() +
  labs(x = "Functional diversity (sown MPD)",
       y = "Complementarity effects (partial residuals)") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11)
  )
#save the new fig.S3
CE.sown.MPD.visreg.fig 
ggsave("figures/FigS3_new.pdf", CE.sown.MPD.visreg.fig, width = 4.5, height = 4, dpi = 300)
ggsave("figures/FigS3_new.tiff", CE.sown.MPD.visreg.fig, width = 4.5, height = 4, dpi = 300)

##Fig.S14(linear model assumption diagnosis)----------------------------------------------------------------
lme_model_list <- list(
  NE = m1_NE,
  CE = m1_CE,
  SE_inter = m1_SEAS,
  SE_intra = m1_SEPS
)

response_labels <- c(
  NE = "Net biodiversity effect",
  CE = "Complementarity effect",
  SE_inter = "Interspecific selection effect",
  SE_intra = "Intraspecific selection effect"
)

diagnostic_panels <- list()

diagnostic_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "plain"),
    axis.text = element_text(colour = "black"),
    plot.tag = element_text(size = 12, face = "plain")
  )

for (response_name in names(lme_model_list)) {
  
  fit <- lme_model_list[[response_name]]
  
  # Conditional residuals divided by the estimated residual SD.
  # These are scaled residuals, not leverage-adjusted studentized residuals.
  diagnostic_data <- data.frame(
    fitted = as.numeric(fitted(fit)),
    residual = as.numeric(residuals(fit)) / sigma(fit)
  )
  
  diagnostic_data <- diagnostic_data[
    is.finite(diagnostic_data$fitted) &
      is.finite(diagnostic_data$residual),
  ]
  
  p_qq <- ggplot(diagnostic_data, aes(sample = residual)) +
    stat_qq(size = 1.3, alpha = 0.65) +
    stat_qq_line(colour = "#D55E00", linewidth = 0.6) +
    labs(
      title = response_labels[[response_name]],
      x = "Theoretical normal quantiles",
      y = "Scaled residuals"
    ) +
    diagnostic_theme
  
  p_residual <- ggplot(
    diagnostic_data,
    aes(x = fitted, y = residual)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    geom_point(size = 1.3, alpha = 0.65) +
    labs(
      title = response_labels[[response_name]],
      x = "Fitted values",
      y = "Scaled residuals"
    ) +
    diagnostic_theme
  
  diagnostic_panels[[paste0(response_name, "_qq")]] <- p_qq
  diagnostic_panels[[paste0(response_name, "_residual")]] <- p_residual
}

LMM_diagnostic_figure <- patchwork::wrap_plots(
  diagnostic_panels,
  ncol = 2,
  byrow = TRUE,
  widths = c(1, 1),
  heights = rep(1, 4)
) + patchwork::plot_annotation(tag_levels = "A")
LMM_diagnostic_figure

ggsave(filename = "figures/Fig_S14_LMM_residual_diagnostics.pdf",plot = LMM_diagnostic_figure, width = 10, height = 12, units = "in")
ggsave(filename = "figures/Fig_S14_LMM_residual_diagnostics.tiff",plot = LMM_diagnostic_figure, width = 10, height = 12, units = "in")
