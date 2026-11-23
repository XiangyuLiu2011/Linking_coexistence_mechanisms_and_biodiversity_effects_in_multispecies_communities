##scripts/01_calculate_biodiversity_effects.R

# Biodiversity effects were calculated using the tripartite partitioning approach.
# The net biodiversity effect (NE) was partitioned into complementarity effects (CE),
# interspecific selection effects (SE_inter), and intraspecific selection effects (SE_intra) (methods see https://doi.org/10.1002/ecm.1521).
library(dplyr)
library(tidyr)

#1.load data
cover <- read.csv("data/PaNDiv_cover_2017_2023_cleaned.csv")
biomass <- read.csv("data/biomass.vip.plots_17_23_cleaned.csv")
# Note:Heracleum sphondylium and Anthriscus sylvestris were poorly established
# in the PaNDiv experiment and were already excluded from the cleaned input data.
ap_Biomass <- merge(cover, biomass, by = c("Year", "Month", "Plot_Nr"),all.x = TRUE)
summary(ap_Biomass)
ap_Biomass <- ap_Biomass %>% dplyr::filter(!is.na(meanAB))

# target.biomass.2 = total plot-level aboveground biomass
# ncover3 = relative cover of species
# Yoi = observed biomass of species i in a plot
ap_Biomass <- ap_Biomass %>% dplyr::mutate(Yoi = target.biomass.2 * ncover3)

#In 2017, Daucus carota had very small biomass in control monoculture plot
#This can make the CE extremely high, even replace june DC control monoculture yield with mean from the other treatments do not change this
#Daucus carota in 2017 June was drop when calculating biodiversity effects, otherwise it drives extremely high CE
excluded_Dc_June2017 <- ap_Biomass %>% dplyr::filter(Harvest == "2017.June", Species == "Daucus carota")
# Exclude Daucus carota from BOTH monocultures and mixtures in June 2017.
ap_Biomass <- ap_Biomass %>% dplyr::filter(!(Harvest == "2017.June" & Species == "Daucus carota"))

#2.Extract monoculture biomass Mi
Mi <- ap_Biomass %>%
  dplyr::filter(SD == 1) %>%
  dplyr::group_by(N, F, Harvest, Species) %>%
  dplyr::summarise(
    Mi = if (all(is.na(Yoi))) {
      NA_real_
    } else {
      mean(Yoi, na.rm = TRUE)
    },
    .groups = "drop"
  )

invalid_Mi_keys <- Mi %>% dplyr::filter(is.na(Mi) | !is.finite(Mi) | Mi <= 0)
print(invalid_Mi_keys, n = Inf)

# Remove monocultures and excluded species from mixture data ("Anthriscus sylvestris" and "Heracleum sphondylium" not well established)
ap_Biomass.target <- ap_Biomass %>% dplyr::filter(SD != 1, !Species %in% c("Anthriscus sylvestris", "Heracleum sphondylium"))

#3.Merge Mi into mixture data
ap_Biomass.1 <- ap_Biomass.target %>%
  dplyr::left_join(
    Mi_for_BEF,
    by = c("N", "F", "Harvest", "Species")
  )

# Record species rows excluded because their reference is unavailable.
excluded_due_to_Mi <- ap_Biomass.1 %>% dplyr::filter(is.na(Mi) | !is.finite(Mi) | Mi <= 0)

print(nrow(excluded_due_to_Mi))
print(summary(excluded_due_to_Mi$Yoi))
print(table(excluded_due_to_Mi$Species))
print(table(excluded_due_to_Mi$Harvest))

# Record the previously identified unreliable plot–harvest.
excluded_plot54_June2018 <- ap_Biomass.1 %>% dplyr::filter(Year == 2018, Month == "June", Plot_Nr == 54)


#4. Apply exclusions before recalculating richness.
ap_Biomass.1 <- ap_Biomass.1 %>%
  dplyr::filter(
    is.finite(Mi),
    Mi > 0,
    !(Year == 2018 & Month == "June" & Plot_Nr == 54)
  )

# Each species must have one row per plot–harvest.
stopifnot(!anyDuplicated(ap_Biomass.1[c("Harvest", "Plot_Nr", "Species")]))
stopifnot(all(is.finite(ap_Biomass.1$Yoi)), all(is.finite(ap_Biomass.1$ncover3)))

#5. Recalculate species richness used in partitioning.
ap_Biomass.1 <- ap_Biomass.1 %>%
  dplyr::group_by(Harvest, Plot_Nr) %>%
  dplyr::mutate(
    Num.sp.used = dplyr::n_distinct(Species)
  ) %>% dplyr::ungroup()

species_count_check <- ap_Biomass.1 %>%
  dplyr::distinct(Harvest, Plot_Nr, Num.sp, Num.sp.used) %>%
  dplyr::mutate(
    lost_species = Num.sp - Num.sp.used)

print(sum(species_count_check$lost_species > 0, na.rm = TRUE))
print(table(ap_Biomass.1$SD, ap_Biomass.1$Num.sp.used))

#6. Calculate relative yields and species contributions.
# equal expected proportions among the retained species.
ap_Biomass.1 <- ap_Biomass.1 %>%
  dplyr::group_by(Harvest, Plot_Nr) %>%
  dplyr::mutate(
    RYei = 1 / Num.sp.used,
    RYoi = Yoi / Mi,
    Yei = RYei * Mi,
    dRYi = RYoi - RYei,
    
    mean_Mi = mean(Mi),
    Sp.contri.on.CE = dRYi * Mi,
    Sp.contri.on.NE = Yoi - Yei,
    # Contributions that sum to the community-level CE.
    Sp.contri.on.CE_partition = dRYi * mean_Mi,
    check_species = Sp.contri.on.CE - Sp.contri.on.NE
  ) %>% dplyr::ungroup()
print(summary(ap_Biomass.1$check_species))

# Confirm Dc is absent from June 2017.
stopifnot(!any(ap_Biomass.1$Harvest == "2017.June" & ap_Biomass.1$Species == "Daucus carota"))

# Check retained richness in previously affected plots.
ap_Biomass.1 %>%
  dplyr::filter(
    Harvest == "2017.June",
    Plot_Nr %in% c(21, 107, 125, 150, 155, 258, 284, 290, 305)
  ) %>%
  dplyr::distinct(Plot_Nr, SD, Num.sp.used) %>%
  dplyr::arrange(Plot_Nr) %>% print(n = Inf)

#7. Inspect the later extreme cases without excluding them.
ap_Biomass.1 %>% dplyr::filter(
    (Harvest == "2018.August" & Plot_Nr == 204) |
      (Harvest == "2023.August" & Plot_Nr == 323) |
       (Harvest == "2017.June" & Plot_Nr == 305)
  ) %>% dplyr::select(
    Harvest, Plot_Nr, Species,
    ncover3, target.biomass.2,
    Mi, Yei, Yoi, RYei, RYoi, dRYi,
    Sp.contri.on.CE,
    Sp.contri.on.CE_partition
  ) %>%   print(n = Inf, width = Inf)

#8. Apply tripartite partitioning.
bef_results <- data.frame()

for (harvest_i in unique(ap_Biomass.1$Harvest)) {
  
  harvest_data <- ap_Biomass.1[ap_Biomass.1$Harvest == harvest_i,]
  
  for (plot_i in unique(harvest_data$Plot_Nr)) {
    
    plot_data <- harvest_data[harvest_data$Plot_Nr == plot_i, ]
    
    # Observed and expected biomass of retained species.
    Yo <- sum(plot_data$Yoi)
    Ye <- sum(plot_data$Yei)
    NE <- Yo - Ye
    mdRY <- mean(plot_data$dRYi)
    
    # Abundance and performance shifts.
    mdRYAS <- mean(plot_data$ncover3 - 1 / plot_data$Num.sp.used)
    mdRYPS <- mean(plot_data$Yoi / plot_data$Mi - plot_data$ncover3)
    mMi <- mean(plot_data$Mi)
    N_used <- plot_data$Num.sp.used[1]
    
    # Complementarity.
    CE <- N_used * mdRY * mMi
    
    # Total selection.
    SE <- sum((plot_data$dRYi - mdRY) * (plot_data$Mi - mMi))
    
    # Interspecific selection.
    SE_inter <- sum(((plot_data$ncover3 - 1 / plot_data$Num.sp.used) - mdRYAS) * (plot_data$Mi - mMi))
    
    # Intraspecific selection.
    SE_intra <- sum(((plot_data$Yoi / plot_data$Mi - plot_data$ncover3) - mdRYPS) * (plot_data$Mi - mMi))
    
    plot_partition <- cbind.data.frame(
      plot_data[1, c("Year", "Month", "Block", "composition", "Plot_Nr", "N", "F", "SD", "Num.sp", "Num.sp.used")],
      NE = NE,
      CE = CE,
      SE = SE,
      SE_inter = SE_inter,
      SE_intra = SE_intra
    )
    
    bef_results <- rbind.data.frame(
      bef_results,
      plot_partition
    )
  }
}

#9. Check partition identities.
# ------------------------------------------------------------
partition_checks <- bef_results %>%
  dplyr::summarise(
    max_error_NE = max(
      abs(NE - CE - SE_inter - SE_intra),
      na.rm = TRUE
    ),
    max_error_SE = max(
      abs(SE - SE_inter - SE_intra),
      na.rm = TRUE
    )
  )

print(partition_checks)
summary(bef_results)

# Three plot–harvest records to investigate.
flagged_partition_records <- tibble::tribble(
  ~Year, ~Month,   ~Plot_Nr,
  2017,  "June",       305,
  2018,  "August",     204,
  2018,  "August",     174,
  2022,  "June",       33,
  2023,  "August",     323
)

# Inspect all species contributing to these records.
# Large RYoi values help identify the small-denominator cases.
flagged_species_details <- ap_Biomass.1 %>%
  dplyr::semi_join(
    flagged_partition_records,
    by = c("Year", "Month", "Plot_Nr")
  ) %>%
  dplyr::select(
    Year, Month, Harvest, Plot_Nr, Species, N, F,
    Mi, Yoi, RYoi, dRYi,
    Sp.contri.on.CE_partition
  ) %>%
  dplyr::arrange(Year, Month, Plot_Nr, dplyr::desc(RYoi))

print(flagged_species_details, n = Inf, width = Inf)

AP_Biomass <- bef_results %>%
  dplyr::left_join(
    flagged_partition_records %>%
      dplyr::mutate(flag_partition = TRUE),
    by = c("Year", "Month", "Plot_Nr")
  ) %>%
  dplyr::mutate(
    flag_partition = dplyr::coalesce(flag_partition, FALSE),
    
    # Preserve the calculated values for comparison.
    CE_before_exclusion = CE,
    SE_intra_before_exclusion = SE_intra,
    SE_before_exclusion = SE,
    
    CE = dplyr::if_else(flag_partition, NA_real_, CE),
    SE_intra = dplyr::if_else(
      flag_partition, NA_real_, SE_intra
    ),
    SE = dplyr::if_else(flag_partition, NA_real_, SE)
  )
head(AP_Biomass)

