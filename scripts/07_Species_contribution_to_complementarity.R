##scripts/07_Species_contribution_to_complementarity
library(dplyr)
library(broom)
library(ggplot2)
library(ggrepel)
library(patchwork)
ap_Biomass.1$Sp.contri.on.CE <- ap_Biomass.1$Mi*ap_Biomass.1$dRYi 

mean.Sp.contri.to.CE <- subset(ap_Biomass.1, SD!=20) %>%
  group_by(Plot_Nr, Month, N, F, Species) %>%
  summarise(
    Plot_Nr = first(Plot_Nr),
    Block = first(block),
    composition = first(composition),
    Month = first(Month),
    Species = first(Species),
    N = first(N),
    F = first(F),
    mdRYi = mean(dRYi, na.rm = TRUE),
    mMi = mean(Mi, na.rm = TRUE),
    mSp.contri.on.CE = mean(Sp.contri.on.CE, na.rm = TRUE),
    mRel.cover = mean(ncover3, na.rm = TRUE),
    .groups = "drop"
  )

mean.Sp.contri.to.CE <- as.data.frame(mean.Sp.contri.to.CE)

mean.Sp.contri.to.CE$Species <- as.factor(mean.Sp.contri.to.CE$Species)
mean.Sp.contri.to.CE <- mean.Sp.contri.to.CE %>% filter(!Species %in% c("Anthriscus sylvestris", "Heracleum sphondylium"))%>%
  droplevels()

levels(mean.Sp.contri.to.CE$Species)<-c("Am", "Ao", "Be","Cj","Cb","Dg","Dc","Fr","Ga","Hp",
                                        "Hl", "Lp", "Pm", "Pt", "Pg", "Ra", "Sp", "To")#keep species same naming style

mean.Sp.contri.to.CE$Treatment[mean.Sp.contri.to.CE$N==0&mean.Sp.contri.to.CE$F==0] <- "control"
mean.Sp.contri.to.CE$Treatment[mean.Sp.contri.to.CE$N==1&mean.Sp.contri.to.CE$F==0] <- "nitrogen"
mean.Sp.contri.to.CE$Treatment[mean.Sp.contri.to.CE$N==0&mean.Sp.contri.to.CE$F==1] <- "fungicide"
mean.Sp.contri.to.CE$Treatment[mean.Sp.contri.to.CE$N==1&mean.Sp.contri.to.CE$F==1] <- "NF"

mean.Sp.contri.to.CE$Treatment<-factor(mean.Sp.contri.to.CE$Treatment, levels=c("control", "nitrogen", "fungicide", "NF"))
mean.Sp.contri.to.CE$Month<-factor(mean.Sp.contri.to.CE$Month, levels = c("June", "August") )  
summary(mean.Sp.contri.to.CE)


Species.level.means <- mean.Sp.contri.to.CE %>%
  group_by(Month, N, Species) %>%
  summarise(
    mmSp.contri.on.CE = mean(mSp.contri.on.CE, na.rm = TRUE),
    se_mmSp.contri.on.CE = sd(mSp.contri.on.CE, na.rm = TRUE) / sqrt(sum(!is.na(mSp.contri.on.CE))),
    mmMi = mean(mMi, na.rm = TRUE),
    mmdRYi = mean(mdRYi, na.rm = TRUE),
    mmRel.cover = mean(mRel.cover, na.rm = TRUE),
    se_mmRel.cover = sd(mRel.cover, na.rm = TRUE) / sqrt(sum(!is.na(mRel.cover))),
    N_treatment = first(N),
    .groups = "drop")
Species.level.means<-as.data.frame(Species.level.means) 
Species.level.means <- Species.level.means %>% filter(!(Species %in% c("As", "Hs")))
Species.level.means$Season <- 1*(Species.level.means$Month=="August")


#facilitation received by species
facil_received_prop_C_June <- apply(all_matrices[[1]], 2, function(x) mean(x > 0))#proportion of facilitative interactions for each column/species
facil_received_prop_N_June <- apply(all_matrices[[2]], 2, function(x) mean(x > 0))
facil_received_prop_F_June <- apply(all_matrices[[3]], 2, function(x) mean(x > 0))
facil_received_prop_NF_June <- apply(all_matrices[[4]], 2, function(x) mean(x > 0))
facil_received_prop_C_August <- apply(all_matrices[[5]], 2, function(x) mean(x > 0))
facil_received_prop_N_August <- apply(all_matrices[[6]], 2, function(x) mean(x > 0))
facil_received_prop_F_August <- apply(all_matrices[[7]], 2, function(x) mean(x > 0))
facil_received_prop_NF_August <- apply(all_matrices[[8]], 2, function(x) mean(x > 0))

# Create a helper function to reshape the data
create_facil_df <- function(facil_received, treatment, month) {
  data.frame(
    Species = names(facil_received),
    treatment = treatment,
    month = month,
    facil_received_prop = as.numeric(facil_received),
    stringsAsFactors = FALSE
  )
}
# Combine all facilitation data into a single data frame
facilitation_df <- rbind(
  create_facil_df(facil_received_prop_C_June, "control", "June"),
  create_facil_df(facil_received_prop_N_June, "nitrogen", "June"),
  create_facil_df(facil_received_prop_F_June, "fungicide", "June"),
  create_facil_df(facil_received_prop_NF_June, "NF", "June"),
  create_facil_df(facil_received_prop_C_August, "control", "August"),
  create_facil_df(facil_received_prop_N_August, "nitrogen", "August"),
  create_facil_df(facil_received_prop_F_August, "fungicide", "August"),
  create_facil_df(facil_received_prop_NF_August, "NF", "August")
)

facilitation_df <- facilitation_df %>%
  dplyr::mutate(N_treatment = dplyr::if_else(treatment %in% c("control", "fungicide"), 0, 1 ) )
facilitation_df_N <- facilitation_df %>%
  group_by(month, N_treatment, Species) %>%
  summarise(
    mfacil_received_prop = mean(facil_received_prop, na.rm = TRUE),
    .groups = "drop")
names(facilitation_df_N)[1] <- "Month"

combined_species.level.means_facilitation <- merge(Species.level.means, facilitation_df_N, by = c("Month", "N_treatment", "Species"), all = TRUE)
head(combined_species.level.means_facilitation)

combined_species.level.means_facilitation %>% group_by(Month, N_treatment) %>% do({
  df_clean <- na.omit(.)  # Remove rows with NAs
  if(nrow(df_clean) > 1) {  # Check that there's enough data
    model <- lm(mmSp.contri.on.CE ~ mfacil_received_prop, data = df_clean)
    tidy_res <- tidy(model)
    glance_res <- glance(model)
    tidy_res$r.squared <- glance_res$r.squared
    tidy_res$adj.r.squared <- glance_res$adj.r.squared
    tidy_res
  } else { tibble(
    term = NA, estimate = NA, std.error = NA,
    statistic = NA, p.value = NA,
    r.squared = NA, adj.r.squared = NA
  )
  } }) %>%  filter(term == "mfacil_received_prop") %>% arrange(p.value)
#June no N: p < 0.01, N: p < 0.10; August no N: p < 0.10, N: p<0.01

combined_species.level.means_facilitation %>% group_by(Month, N_treatment) %>% do({
  df_clean <- na.omit(.)  # Remove rows with NAs
  if(nrow(df_clean) > 1) {  # Check that there's enough data
    model <- lm( mmRel.cover ~ mfacil_received_prop, data = df_clean)
    tidy_res <- tidy(model)
    glance_res <- glance(model)
    tidy_res$r.squared <- glance_res$r.squared
    tidy_res$adj.r.squared <- glance_res$adj.r.squared
    tidy_res
  } else { tibble(
    term = NA, estimate = NA, std.error = NA,
    statistic = NA, p.value = NA,
    r.squared = NA, adj.r.squared = NA
  )
  } }) %>%  filter(term == "mfacil_received_prop") %>% arrange(p.value)
#all significant

######
##code for Fig.5 and Fig.S6
plot_facil_df <- combined_species.level.means_facilitation %>%
  dplyr::mutate(
    Month = as.character(Month),
    N_label = factor(
      N_treatment,
      levels = c(0, 1),
      labels = c("Without nitrogen", "With nitrogen")
    )
  )

focal_species_selected <- c("Ra", "Ao", "Dc", "To", "Pt")

## common plotting function
make_facil_panel <- function(data, month_select, y_var, se_var, y_lab,
                             show_legend = FALSE,
                             use_linetype = FALSE,
                             reverse_linetype = FALSE) {
  
  plot_df <- data %>%
    dplyr::filter(Month == month_select,
      !is.na(mfacil_received_prop),
      !is.na(.data[[y_var]])
    )
  
  label_df <- plot_df %>%
    dplyr::filter(Species %in% focal_species_selected)
  
  p <- ggplot(
    plot_df,
    aes(
      x = mfacil_received_prop,
      y = .data[[y_var]],
      color = N_label
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.02),
      size = 2,
      alpha = 0.9
    ) +
    geom_errorbar(
      aes(
        ymin = .data[[y_var]] - .data[[se_var]],
        ymax = .data[[y_var]] + .data[[se_var]]
      ),
      position = position_dodge(width = 0.02),
      width = 0.02,
      alpha = 0.6,
      linewidth = 0.5,
      na.rm = TRUE
    )  +
    labs(
      x = "Proportion of facilitation",
      y = y_lab,
      color = "Treatment",
      fill = "Treatment",
      linetype = "Treatment"
    ) +
    scale_color_manual(
      values = c(
        "Without nitrogen" = "#99CCFF",
        "With nitrogen" = "#336699"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Without nitrogen" = "#99CCFF",
        "With nitrogen" = "#336699"
      )
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 11, face = "plain"),
      strip.background = element_blank(),
      axis.title = element_text(size = 11, face = "plain"),
      axis.text = element_text(size = 11, face = "plain"),
      panel.grid = element_blank(),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  if (use_linetype) {
    
    line_values <- if (reverse_linetype) {
      c(
        "Without nitrogen" = "solid",
        "With nitrogen" = "dashed"
      )
    } else {
      c(
        "Without nitrogen" = "dashed",
        "With nitrogen" = "solid"
      )
    }
    
    p <- p +
      geom_smooth(
        aes(fill = N_label, linetype = N_label),
        method = "lm",
        se = TRUE,
        linewidth = 1,
        alpha = 0.1
      ) + 
      scale_linetype_manual(values = line_values)
    
  } else {
    
    p <- p +
      geom_smooth(
        aes(fill = N_label),
        method = "lm",
        linetype = "solid",
        se = TRUE,
        linewidth = 1,
        alpha = 0.1
      ) 
  }
  p <- p +
    geom_text_repel(
      data = label_df,
      aes(label = Species),
      size = 2.8,
      box.padding = 0.2,
      point.padding = 0.3,
      segment.alpha = 0.3,
      force = 2,
      max.iter = 5000,
      min.segment.length = 0,
      color = "black")
  return(p)
}

Fig.4_new_a1 <- make_facil_panel(
  data = plot_facil_df,
  month_select = "June",
  y_var = "mmSp.contri.on.CE",
  se_var = "se_mmSp.contri.on.CE",
  y_lab = "Species' contribution to complementarity",
  show_legend = FALSE,
  use_linetype = TRUE,
  reverse_linetype = TRUE
)

Fig.4_new_b2 <- make_facil_panel(
  data = plot_facil_df,
  month_select = "August",
  y_var = "mmSp.contri.on.CE",
  se_var = "se_mmSp.contri.on.CE",
  y_lab = "Species' contribution to complementarity",
  show_legend = FALSE,
  use_linetype = TRUE,
  reverse_linetype = FALSE
)

Fig.4_new_c3 <- make_facil_panel(
  data = plot_facil_df,
  month_select = "June",
  y_var = "mmRel.cover",
  se_var = "se_mmRel.cover",
  y_lab = "Species' mean relative abundance",
  show_legend = TRUE,
  use_linetype =  FALSE,
  reverse_linetype = FALSE
)
Fig.4_new_d4 <- make_facil_panel(
  data = plot_facil_df,
  month_select = "August",
  y_var = "mmRel.cover",
  se_var = "se_mmRel.cover",
  y_lab = "Species' mean relative abundance",
  show_legend = TRUE,
  use_linetype = FALSE,
  reverse_linetype = FALSE
)

Fig4.new_June <- (Fig.4_new_a1 + Fig.4_new_c3)  +  plot_annotation(tag_levels = 'A')
Fig4.new_August <- (Fig.4_new_b2 +  Fig.4_new_d4) +  plot_annotation(tag_levels = 'A')
ggsave("figures/Fig5.pdf", Fig4.new_June, width = 8, height = 3.5)
ggsave("figures/FigS6.pdf", Fig4.new_August, width = 8, height = 3.5)
