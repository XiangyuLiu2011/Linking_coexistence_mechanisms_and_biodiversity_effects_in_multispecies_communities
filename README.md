In this folder you may find the data and code for the analysis of "**Linking mechanisms of species coexistence and biodiversity effects on productivity in multispecies communities under environmental variation**".



What can you find in this project?

**a) data**

Competition coefficients and intrinsic growth rates for coexistence mechanisms

In the data folder, "biommatrix" and "biomintrinsic" files containing all competition coefficients and intrinsic growth rates used during the analysis. Text files which start with "biommatrix" contain competition coefficients, and their associated standard errors, which are named "SE\_biommatrix". The species intrinsic growth rates are also present in this folder under the name starting with "biomintrinsic". These matrices are derived over the whole year 2020 in two harvest seasons, June and August, respectively (see Daniel et al. 2024, Ecology Letters, for more information on the method).



Biomass and cover from 2017 to 2023 for biodiversity effects

“PaNDiv\_cover\_2017\_2023\_cleaned.csv” and “biomass.vip.plots\_17\_23\_cleaned.csv” contain species level cover and biomass in each plot in June and August harvest from 2017 to 2023. These two files were used for partitioning biodiversity effects. 



In “PaNDiv\_cover\_2017\_2023\_cleaned.csv”: block,	Plot\_Nr (plot number), unique.composition.id (names of species composition), functional.composition, species.diversity (sown species richness),	 nitrogen (1 == with nitrogen; 0 == without nitrogen),	fungicide (1 == with fungicide; 0 == without fungicide),	Month (June and August),	Year (2017 to 2023),	Species	cover (cover of each species), total.cover (total cover including weeds),	Harvest (a combined information of Year and Month),	 ncover3 (species relative cover with weeds, Hs and As excluded), num\_excluded\_species (number of species excluded As and Hs),Num.sp (number of species without As and Hs). In “biomass.vip.plots\_17\_23\_cleaned.csv”: meanAB (average of two sampled biomass),	total.coverw (total cover with weeds),rel.weeds (relative cover of weeds),target.biomass (total biomass without weeds),	rel.As (relative cover of As),	rel.Hs (relative cover of Hs),target.biomass.2 (total biomass without weeds, As and Hs) \[I don’t repeat the same plot level information].



“Plot.info.csv” and “Plot\_Descriptions.txt” contain information like block, plot number (Plot\_Nr), species richness (SD), N (nitrogen addition: 1 == with nitrogen; 0 == without nitrogen), F (fungicide application: 1 == with fungicide; 0 == without fungicide), functional composition (in category, F: fast community, M: mixed with fast- and slow-species, S: slow community), VIP (1 indicates plots that sown composition was maintained by weeding). “Plot\_Descriptions.txt” also contains species that we had sown in each plot. “plot.info.sla.updated.csv” contains the plot information, as well as the sown CWMSLA (community weighted mean of specific leaf area based on the sown composition) and sown MPD (mean pairwise distance in specific leaf area based on the sown composition, see Pichon et al. 2024 Ecology letters).



**b) Intermediate results**

Coexistence mechanisms

Because we developed a randomization approach to calculate whether a given composition had higher or lower niche and competitive ability differences than expected under a random reshuffling of interaction coefficients and intrinsic growth rates (predicted coexistence mechanims; see methods). Initially we calculated coexistence mechanisms without facilitation (by setting positive interactions to close to zero; see methods), and later tested whether including facilitation altered links between coexistence mechanism and biodiversity effects. To save time for the readers who have interest in our analysis. Intermediate results for coexistence mechanisms were saved in “intermediate\_results” folder. “dataset\_results\_ini.rds” and “dataset\_results\_ini\_fac.rds” are observed coexistence mechanisms without and with facilitation. Files with names like “dataset\_results\_ini\_rand1\_ad.rds” and “dataset\_results\_ini\_rand1\_fac.rds” are coexistence mechanisms of each randomization without and with facilitation (99 times in total). These files are labeled from 1 to 8, aligning the "biommatrix" and “biomintrinsic”, representing control, nitrogen, fungicide and combined treatments in June and August, respectively. 



Biodiversity effects

In the output folder, there are two files containing the biodiversity effects at each harvest season in every year (“biodiversity\_effects\_2017\_2023\_outliers\_removed.csv”) and the averaged biodiversity effects (“Mean.Biomass.Tripartite17\_23.csv”).



**c) functions and code**

This repository contains modified functions from previous studies:



"toolbox" contains functions from the paper: "A structural approach for understanding multispecies coexistence", by Serguei Saavedra et al., Ecology Monographs (2017). This code was modified to calculate structural niche differences (Omega), structural fitness differences (theta) and indirect interactions (differential).



"useful SEM constrain functions" contain a set of functions that were used to make the stepwise constraining of the multigroup SEM model and to check the model output. “useful SEM partial plot functions” contain a set of functions that were used to make the partial plots from best SEM model.



This R project contains a set of original code stored in “scripts” folder. These code files are numbered from (1) to (7) and they contain the main analyses of the study. Each file can be run independently, as they all source all files that are needed to run the analysis. However, the code was build in this specific order and the logic of the analysis might be easier to follow if each file is read in this order.



