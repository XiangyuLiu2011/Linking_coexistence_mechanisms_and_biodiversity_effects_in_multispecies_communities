##scripts/01_calculate_biodiversity_effects.R
# Biodiversity effects were calculated using the tripartite partitioning approach.
# The net biodiversity effect (NE) was partitioned into complementarity effects (CE),
# interspecific selection effects (SE_inter), and intraspecific selection effects (SE_intra).
library(dplyr)
library(tidyr)

#1.load data
cover <- read.csv("data/PaNDiv_cover_2017_2023_cleaned.csv")
biomass <- read.csv("data/biomass.vip.plots_17_23_cleaned.csv")
# Note:
# Heracleum sphondylium and Anthriscus sylvestris were poorly established
# in the PaNDiv experiment and were already excluded from the cleaned input data.

ap_Biomass <- merge(cover, biomass, by=c("Year", "Month", "Plot_Nr"), all.x=T)
summary(ap_Biomass)
ap_Biomass <- ap_Biomass[!is.na(ap_Biomass$meanAB),]# Keep rows with valid biomass information

# target.biomass.2 = total plot-level aboveground biomass
# ncover3 = relative cover / relative contribution of species i
# Yoi = observed biomass of species i in a plot
ap_Biomass %>% filter(is.na(target.biomass.2))
ap_Biomass <- ap_Biomass %>% mutate(Yoi = target.biomass.2 * ncover3)

#In 2017, Dc had very small biomass in control monoculture plot, and will make the CE extremely high, even replace june DC control monoculture yield with mean from the other treatments do not change this
#Dc in 2017 June was drop when calculating biodiversity effects, otherwise it drives extremely high CE
ap_Biomass[ap_Biomass$Harvest=="2017.June" & ap_Biomass$SD==1 & ap_Biomass$composition=="Dc", "target.biomass.2"] <- rep(NA,4)

#2.Extract monoculture biomass Mi
Mi <- ap_Biomass %>%
  dplyr::filter(SD == 1) %>%
  dplyr::group_by(N, F, Harvest, Species) %>%
  dplyr::summarise(
    Mi = mean(Yoi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Mi = dplyr::if_else(is.nan(Mi), NA_real_, Mi)
  )

# Check monoculture reference biomass
cat("\nProblematic monoculture reference biomass:\n")
invalid_Mi_keys <- Mi %>%
  dplyr::filter(is.na(Mi) | Mi <= 0)
#there are a few cases where monoculture biomass of species were zero, these species were excluded

Mi_for_BEF <- Mi %>% dplyr::filter(!is.na(Mi), Mi > 0)

# Remove monocultures and excluded species from mixture data ("Anthriscus sylvestris" and "Heracleum sphondylium" not well established)
ap_Biomass.target <- ap_Biomass %>%
  filter(
    SD != 1,
    !Species %in% c("Anthriscus sylvestris", "Heracleum sphondylium")
  )
#3.Merge Mi into mixture data
ap_Biomass.1 <- ap_Biomass.target %>%
  dplyr::left_join(
    Mi_for_BEF,
    by = c("N", "F", "Harvest", "Species") )

excluded_due_to_Mi <- ap_Biomass.1 %>%
  dplyr::filter(is.na(Mi) | Mi <= 0)
cat("\nRows excluded because Mi is missing or invalid:\n")
print(nrow(excluded_due_to_Mi))
summary(excluded_due_to_Mi$Yoi)#their biomasses were quite small

#4. Remove rows with invalid Mi
invalid_mi <- ap_Biomass.1 %>% filter(is.na(Mi) | Mi <= 0)
cat("\nRows with missing or invalid Mi after merging:\n")
print(nrow(invalid_mi))

cat("\nInvalid Mi by species:\n")
print(table(invalid_mi$Species))

cat("\nInvalid Mi by harvest:\n")
print(table(invalid_mi$Harvest))

ap_Biomass.1 <- ap_Biomass.1 %>% dplyr::filter(!is.na(Mi), Mi > 0)
#5. Recalculate species richness used in BEF partitioning
# Num.sp.used is the number of species actually used in the BEF calculation
# after removing rows with unavailable monoculture reference biomass.
ap_Biomass.1 <- ap_Biomass.1 %>%
  group_by(Harvest, Plot_Nr) %>%
  mutate(Num.sp.used = n_distinct(Species)) %>% ungroup()

# Check whether any plot lost species because of invalid Mi
species_count_check <- ap_Biomass.1 %>%
  distinct(Harvest, Plot_Nr, Num.sp, Num.sp.used) %>%
  mutate(lost_species = Num.sp - Num.sp.used)
print(sum(species_count_check$lost_species > 0))

table(ap_Biomass.1$SD, ap_Biomass.1$Num.sp.used)#in 2018 june, plot54 is a 8 species plot but only record with one species To, which is impossible, excluded from calculation
ap_Biomass.1 <- ap_Biomass.1 %>% filter(!(Year == 2018 & Month == "June" & Plot_Nr == 54 & Species == "Taraxacum officinale") )

#6. Calculate relative yields and expected yields
ap_Biomass.1 <- ap_Biomass.1 %>%
  mutate(
    RYei = 1 / Num.sp.used,   # expected relative yield
    RYoi = Yoi / Mi,          # observed relative yield
    Yei  = RYei * Mi,         # expected yield
    dRYi = RYoi - RYei        # deviation from expected relative yield
  )

# Species-level identity check:
# dRYi * Mi should equal Yoi - Yei.
ap_Biomass.1 <- ap_Biomass.1 %>%
  mutate(check_species = dRYi * Mi - (Yoi - Yei))
print(summary(ap_Biomass.1$check_species))#close to zero means they are very close

##7. Apply tripartite partitioning
bef_results <- data.frame()

for (harvest_i in unique(ap_Biomass.1$Harvest)) {
  
  harvest_data <- ap_Biomass.1[ap_Biomass.1$Harvest == harvest_i, ]
  
  for (plot_i in unique(harvest_data$Plot_Nr)) {
    
    plot_data <- harvest_data[harvest_data$Plot_Nr == plot_i, ]
    
    # Observed and expected yield
    Yo <- sum(plot_data$Yoi)
    Ye <- sum(plot_data$Yei)
    NE <- Yo - Ye
    
    # Mean deviation in relative yield
    mdRY <- mean(plot_data$dRYi)
    
    # Components based on abundance shift and performance shift
    mdRYAS <- mean(plot_data$ncover3 - 1 / plot_data$Num.sp.used)
    mdRYPS <- mean(plot_data$Yoi / plot_data$Mi - plot_data$ncover3)
    
    # Mean monoculture yield and species richness used
    mMi <- mean(plot_data$Mi)
    N_used <- plot_data$Num.sp.used[1]
    
    # Complementarity effect
    CE <- N_used * mdRY * mMi
    
    # Total selection effect
    SE <- sum((plot_data$dRYi - mdRY) * (plot_data$Mi - mMi))
    
    # Interspecific selection effect
    SE_inter <- sum(
      ((plot_data$ncover3 - 1 / plot_data$Num.sp.used) - mdRYAS) *
        (plot_data$Mi - mMi)
    )
    
    # Intraspecific selection effect
    SE_intra <- sum(
      ((plot_data$Yoi / plot_data$Mi - plot_data$ncover3) - mdRYPS) *
        (plot_data$Mi - mMi)
    )
    
    plot_partition <- cbind.data.frame(
      plot_data[
        c(
          "Year", "Month", "Block", "composition", "Plot_Nr",
          "N", "F", "SD", "Num.sp", "Num.sp.used"
        )
      ][1, ],
      NE = NE,
      CE = CE,
      SE = SE,
      SE_inter = SE_inter,
      SE_intra = SE_intra
    )
    
    bef_results <- rbind.data.frame(bef_results, plot_partition)
  }
}

AP_Biomass <- bef_results
summary(AP_Biomass)

AP_Biomass$check_NE <- AP_Biomass$CE + AP_Biomass$SE - AP_Biomass$NE
summary(AP_Biomass$check_NE)
AP_Biomass$check_SE <- AP_Biomass$SE_inter + AP_Biomass$SE_intra - AP_Biomass$SE
summary(AP_Biomass$check_SE)

#I remove outliers based on between 6 times of 95% upper and 5% lower quantile (Wagg 2022 Nature Communications)
# Function to set values outside 6 times the 5% and 95% quantiles to NA
remove_outliers_quantile <- function(x, multiplier = 6) {
  
  q_low <- quantile(x, probs = 0.05, na.rm = TRUE)
  q_high <- quantile(x, probs = 0.95, na.rm = TRUE)
  
  lower_limit <- multiplier * q_low
  upper_limit <- multiplier * q_high
  
  x[x <= lower_limit | x >= upper_limit] <- NA
  
  return(x)
}

AP_Biomass <- AP_Biomass %>%
  mutate(
    NE = remove_outliers_quantile(NE, multiplier = 6),
    CE = remove_outliers_quantile(CE, multiplier = 6),
    SE_inter = remove_outliers_quantile(SE_inter, multiplier = 6),
    SE_intra = remove_outliers_quantile(SE_intra, multiplier = 6)
  )
summary(AP_Biomass)
#write.csv(AP_Biomass, "outputs/biodiversity_effects_2017_2023_outliers_removed.csv", row.names = FALSE)
