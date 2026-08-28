#### Load packages ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(sf)
library(ggplot2)
library(patchwork)
library(colorspace)
library(readxl)
library(parallel)

#### Supporting scripts ####
source("Code/12_Analyses/12_00_LoadData.R")
source("Code/99_Supporting/genetic_algorithm.R")

#### Read in additional data ####
# Plant checklist with additional attributes for scenario filtering
plant_checklist <- read.csv("Data_Clean/Species_Checklists/plant_checklist.csv")

n_plants <- 6

# Filter predicted interactions to only GEOOS + simple for flower visitation
predicted_interactions <- predicted_interactions %>%
  filter(method_type %in% c("host", "Level1_Simple", "Level2_GEOOS"))

#### find_randomMix: random plant selection baseline ####
find_randomMix <- function(ecoregion,
                           interaction_data_subset,
                           candidate_plants,
                           focal_pollinators,
                           checklist,
                           n_plants,
                           sp_by_region,
                           n_iterations = 100) {
  
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  plants <- candidate_plants[candidate_plants %in% checklist_region[checklist_region$taxon == "plants", ]$genus_species] %>%
    unique()
  
  pollinators <- focal_pollinators[focal_pollinators %in% checklist_region[checklist_region$taxon %in% c("bees", "moths", "butterflies", "hoverflies"), ]$genus_species] %>%
    unique()
  
  interaction_matrix_filtered <- build_interaction_matrix(
    ecoregion_data,
    plants = plants,
    pollinators = pollinators
  )
  
  if (nrow(interaction_matrix_filtered) < n_plants) {
    warning(paste("Only", nrow(interaction_matrix_filtered),
                  "plants available in ecoregion", ecoregion))
    return(list(
      random_mix = data.frame(
        ecoregion = ecoregion,
        iteration = NA,
        n_plants = nrow(interaction_matrix_filtered),
        pollinators_supported = NA
      ),
      taxon_summary = NULL
    ))
  }
  
  all_pollinators_in_region <- checklist_region %>%
    filter(taxon %in% c("bees", "moths", "butterflies", "hoverflies"))
  
  results <- lapply(1:n_iterations, function(iter) {
    random_plants <- sample(rownames(interaction_matrix_filtered), n_plants)
    covered_pollinators <- apply(interaction_matrix_filtered[random_plants, , drop = FALSE], 2, any)
    data.frame(
      ecoregion = ecoregion,
      iteration = iter,
      n_plants = n_plants,
      pollinators_supported = sum(covered_pollinators)
    )
  })
  
  random_mix_results <- do.call(rbind, results)
  
  taxon_summary_list <- lapply(1:n_iterations, function(iter) {
    random_plants <- sample(rownames(interaction_matrix_filtered), n_plants)
    covered_pollinators <- apply(interaction_matrix_filtered[random_plants, , drop = FALSE], 2, any)
    supported_pollinator_names <- names(covered_pollinators)[covered_pollinators]
    
    taxon_totals <- all_pollinators_in_region %>%
      group_by(taxon) %>%
      summarise(total_in_region = n(), .groups = "drop")
    
    taxon_supported <- all_pollinators_in_region %>%
      filter(genus_species %in% supported_pollinator_names) %>%
      group_by(taxon) %>%
      summarise(supported_by_mix = n(), .groups = "drop")
    
    taxon_totals %>%
      left_join(taxon_supported, by = "taxon") %>%
      mutate(
        supported_by_mix = replace_na(supported_by_mix, 0),
        percent_supported = round(100 * supported_by_mix / total_in_region, 1),
        ecoregion = ecoregion,
        iteration = iter
      ) %>%
      dplyr::select(ecoregion, iteration, taxon, total_in_region, supported_by_mix, percent_supported)
  })
  
  taxon_summary <- do.call(rbind, taxon_summary_list)
  
  taxon_avg_summary <- taxon_summary %>%
    group_by(ecoregion, taxon, total_in_region) %>%
    summarise(
      mean_supported = round(mean(supported_by_mix), 1),
      sd_supported   = round(sd(supported_by_mix), 1),
      mean_percent   = round(mean(percent_supported), 1),
      sd_percent     = round(sd(percent_supported), 1),
      .groups = "drop"
    )
  
  message(paste("\n=== Average Taxon Group Summary for", ecoregion,
                "(across", n_iterations, "iterations) ==="))
  print(taxon_avg_summary)
  
  return(list(
    random_mix    = random_mix_results,
    taxon_summary = taxon_summary
  ))
}

### Scenario 0: Unconstrained ###
topMix_noConstraints <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_top_mix(
    ecoregion,
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants        = unique(plant_checklist$scientificName),
    focal_pollinators       = unique(checklist_cleaned$scientificName[checklist_cleaned$taxon %in% c("bees", "butterflies", "moths", "hoverflies")]),
    checklist               = checklist_cleaned,
    sp_by_region            = sp_by_region,
    n_plants                = n_plants,
    n_generations           = 200,
    population_size         = 200,
    mutation_rate           = 0.1
  )
}, mc.cores = 12)

scenario0_results      <- do.call(rbind, lapply(topMix_noConstraints, function(x) x$top_mix))
scenario0_taxon_summary <- do.call(rbind, lapply(topMix_noConstraints, function(x) x$taxon_summary))

write.csv(scenario0_results,
          "Data_Clean/Scenario_Results/Scenario0_topMixUnconstrained.csv")
write.csv(scenario0_results %>% group_by(ecoregion) %>% slice_head(n = 1) %>% dplyr::select(-c(rank, scientificName)),
          "Data_Clean/Scenario_Results/Scenario0_pols.csv")

random_scenario0_results_list <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_randomMix(
    ecoregion,
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants        = unique(plant_checklist$scientificName),
    focal_pollinators       = unique(checklist_cleaned$scientificName[checklist_cleaned$taxon %in% c("bees", "butterflies", "moths", "hoverflies")]),
    checklist               = checklist_cleaned,
    n_plants                = n_plants,
    sp_by_region            = sp_by_region,
    n_iterations            = 100
  )
}, mc.cores = 12)

random_scenario0_results      <- do.call(rbind, lapply(random_scenario0_results_list, function(x) x$random_mix))
random_scenario0_taxon_summary <- do.call(rbind, lapply(random_scenario0_results_list, function(x) x$taxon_summary))

### quick summary for redding vs LA
# Redding
scen0_supportByRegion <- scenario0_taxon_summary %>% 
  group_by(ecoregion) %>%
  summarise(total_in_region = sum(total_in_region),
            supported_by_mix = sum(supported_by_mix)) %>%
  mutate(perc_supported = supported_by_mix/total_in_region*100) 


## Plot
prepare_pollinator_plot_data <- function(ecoregion_code, 
                                         optimized_results, 
                                         random_taxon_summary) {
  
  # Get random scenario stats by taxon for this ecoregion
  random_stats <- random_taxon_summary %>%
    filter(ecoregion == ecoregion_code) %>%
    group_by(taxon, total_in_region) %>%
    summarise(
      mean_percent_random = mean(percent_supported),
      sd_percent_random = sd(percent_supported),
      se_percent_random = sd(percent_supported) / sqrt(n()),
      ci_lower = mean(percent_supported) - 1.96 * se_percent_random,
      ci_upper = mean(percent_supported) + 1.96 * se_percent_random,
      .groups = 'drop'
    )
  
  # Combine with your optimized results
  pollinator_data <- optimized_results %>%
    filter(ecoregion == ecoregion_code) %>%
    left_join(random_stats, by = c("taxon", "total_in_region")) %>%
    mutate(
      taxon_group = tools::toTitleCase(taxon),
      proportion = total_in_region / sum(total_in_region),
      end_angle = cumsum(proportion) * 2 * pi,
      start_angle = lag(end_angle, default = 0),
      mid_angle = (start_angle + end_angle) / 2
    )
  
  return(pollinator_data)
}

plot_pollinator_support_simple <- function(pollinator_data, ecoregion_name) {
  
  tick_positions <- seq(0, 100, by = 25)
  tick_labels <- c(0, 25, 50, 75, 100)
  
  colors_optimized <- c(
    "Bees" = "#FFCC00",
    "Butterflies" = "#E69190",
    "Moths" = "#BD9ECD",
    "Hoverflies" = "#F76700",
    "Hummingbirds" = "#52A273"
  )
  
  ggplot(pollinator_data) +
    
    # Radial grid lines
    geom_segment(
      data = data.frame(x = tick_positions / 100),
      aes(x = x, xend = x, y = 0, yend = 2 * pi),
      color = "grey70",
      linewidth = 0.2,
      alpha = 0.5
    ) +
    
    # Optimized scenario (solid colors)
    geom_rect(
      aes(xmin = 0, xmax = percent_supported/100, 
          ymin = start_angle, ymax = end_angle,
          fill = taxon_group),
      alpha = 1,
      color = "white",
      linewidth = 0.3
    ) +
    
    # Percentage labels
    geom_text(
      data = data.frame(
        x = tick_labels / 100,
        label = paste0(tick_labels, "%")
      ),
      aes(x = x, y = 0, label = label),
      size = 4 / .pt,
      color = "black",
      vjust = -0.5,
      family = "Arial"
    ) +
    
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors_optimized) +
    
    theme_void() +
    theme(
      legend.position = "none",
      text = element_text(family = "Arial", size = 4),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    xlim(0, 1.2)
}# Make plots
SCo_data <- prepare_pollinator_plot_data("SCo",
                                         scenario0_taxon_summary,
                                         random_scenario0_taxon_summary)
plot_SCo <- plot_pollinator_support_simple(SCo_data, "South Coast (SCo)")

CaRF_data <- prepare_pollinator_plot_data("CaRF",
                                          scenario0_taxon_summary,
                                          random_scenario0_taxon_summary)
plot_CaRF <- plot_pollinator_support_simple(CaRF_data, "Central Coast Ranges (CaRF)")

# Save
ggsave("Figures/scenario0_pie/scenario0_SCo.png", 
       plot = plot_SCo,
       width = 35, 
       height = 35, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/scenario0_pie/scenario0_SCo.pdf", 
       plot = plot_SCo,
       width = 35, 
       height = 35, 
       units = "mm",
       device = cairo_pdf)

ggsave("Figures/scenario0_pie/scenario0_CaRF.png", 
       plot = plot_CaRF,
       width = 35, 
       height = 35, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/scenario0_pie/scenario0_CaRF.pdf", 
       plot = plot_CaRF,
       width = 35, 
       height = 35, 
       units = "mm",
       device = cairo_pdf)


#### Scenario 1: habitat quality improvement ####
## Objective: Maximize pollinator richness supported
## Constraints: Chaparral habitat, common (remove plants with RPR 1-3)

## Filter candidate plant list
# Read in supporting constraint data
chaparral_plants <- read.csv("Data_Clean/Species_Traits_Attributes/Calflora_Communities.csv") %>%
  dplyr::select(genus_species_harm, community) %>%
  rename(lower = genus_species_harm) %>%
  unique() %>%
  # only chaparral
  filter(community == "chaparral")

# Filter plant_checklist to candidates based on constraints
plantCandidates_scenario1 <- plant_checklist %>%
  # constrain 1: filter out rare plants (Anything with CRPR 1-3)
  filter(!CRPR %in% c("1A", "1B.1", "1B.2", "1B.3", "2A", "2B.1", "2B.2", "2B.3", "3", "3.1", "3.2", "3.3")) %>%
  # constraint 2: found in chaparral habitat
  filter(scientificName %in% chaparral_plants$lower)

focalPollinators_scenario1 <- checklist_cleaned %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies"))

scenario1_results_list <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_top_mix(
    ecoregion, 
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants = unique(plantCandidates_scenario1$scientificName),
    focal_pollinators = unique(focalPollinators_scenario1$scientificName),
    checklist = checklist_cleaned,
    sp_by_region = sp_by_region,
    n_plants = n_plants,
    n_generations = 200, 
    population_size = 200,
    mutation_rate = 0.1
  )
},
mc.cores = 12)


scenario1_results <- do.call(rbind, lapply(scenario1_results_list, function(x) x$top_mix))
scenario1_taxon_summary <- do.call(rbind, lapply(scenario1_results_list, function(x) x$taxon_summary))

# Save scenario1_results
write.csv(scenario1_results, 
          "Data_Clean/Scenario_Results/Scenario1_topmix.csv")
write.csv(scenario1_results %>% group_by(ecoregion) %>% slice_head(n=1) %>% dplyr::select(-c(rank, scientificName)), 
          "Data_Clean/Scenario_Results/Scenario1_pols.csv")



## Compare to random draw
random_scenario1_results_list <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_randomMix(
    ecoregion, 
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants = unique(plantCandidates_scenario1$scientificName),
    focal_pollinators = unique(focalPollinators_scenario1$scientificName),
    checklist = checklist_cleaned,
    n_plants = n_plants,
    sp_by_region = sp_by_region,
    n_iterations = 100
  )
},
mc.cores = 12)

# Extract and combine the tables
random_scenario1_results <- do.call(rbind, lapply(random_scenario1_results_list, function(x) x$random_mix))
random_scenario1_taxon_summary <- do.call(rbind, lapply(random_scenario1_results_list, function(x) x$taxon_summary))

## Summary how many supported in random total vs optimized
scen1_random_propSupported <- random_scenario1_taxon_summary %>% 
  group_by(taxon) %>% 
  summarise(support_prop = mean(percent_supported)) 

scen1_propSupported <- scenario1_taxon_summary %>% 
  group_by(taxon) %>% 
  summarise(support_prop = mean(percent_supported)) 

scen1_propSupported %>% 
  left_join(scen1_random_propSupported, by = "taxon") %>%
  mutate(absolute_inc = support_prop.x-support_prop.y,
         relative_inc = (support_prop.x-support_prop.y)/support_prop.y*100)


## Plot
prepare_pollinator_plot_data <- function(ecoregion_code, 
                                         optimized_results, 
                                         random_taxon_summary) {
  
  # Get random scenario stats by taxon for this ecoregion
  random_stats <- random_taxon_summary %>%
    filter(ecoregion == ecoregion_code) %>%
    group_by(taxon, total_in_region) %>%
    summarise(
      mean_percent_random = mean(percent_supported),
      sd_percent_random = sd(percent_supported),
      se_percent_random = sd(percent_supported) / sqrt(n()),
      ci_lower = mean(percent_supported) - 1.96 * se_percent_random,
      ci_upper = mean(percent_supported) + 1.96 * se_percent_random,
      .groups = 'drop'
    )
  
  # Combine with your optimized results
  pollinator_data <- optimized_results %>%
    filter(ecoregion == ecoregion_code) %>%
    left_join(random_stats, by = c("taxon", "total_in_region")) %>%
    mutate(
      taxon_group = tools::toTitleCase(taxon),
      proportion = total_in_region / sum(total_in_region),
      end_angle = cumsum(proportion) * 2 * pi,
      start_angle = lag(end_angle, default = 0),
      mid_angle = (start_angle + end_angle) / 2
    )
  
  return(pollinator_data)
}

plot_pollinator_support <- function(pollinator_data, ecoregion_name) {
  
  tick_positions <- seq(0, 100, by = 25)
  tick_labels <- c(0, 25, 50, 75, 100)
  
  colors_optimized <- c(
    "Bees" = "#FFCC00",
    "Butterflies" = "#E69190",
    "Moths" = "#BD9ECD",
    "Hoverflies" = "#F76700",
    "Hummingbirds" = "#52A273"
  )
  
  colors_random <- lighten(colors_optimized, amount = 0.4)
  
  # Calculate error bar data with radius-adjusted cap width
  pollinator_data <- pollinator_data %>%
    mutate(
      radius_ci_lower = ci_lower/100,
      radius_ci_upper = ci_upper/100,
      radius_mean = mean_percent_random/100,
      # Make angular offset inversely proportional to radius
      # so visual width stays constant
      cap_offset_lower = 0.03 / pmax(radius_ci_lower, 0.01),  # Adjust 0.05 to change cap size
      cap_offset_upper = 0.03 / pmax(radius_ci_upper, 0.01),
      cap_offset_mean = 0.03 / pmax(radius_mean, 0.01)
    )
  
  ggplot(pollinator_data) +
    
    # Radial grid lines
    geom_segment(
      data = data.frame(x = tick_positions / 100),
      aes(x = x, xend = x, y = 0, yend = 2 * pi),
      color = "grey70",
      linewidth = 0.2,
      alpha = 0.5
    ) +
    
    # Optimized scenario (solid colors)
    geom_rect(
      aes(xmin = 0, xmax = percent_supported/100, 
          ymin = start_angle, ymax = end_angle,
          fill = taxon_group),
      alpha = 1,
      color = "white",
      linewidth = 0.3
    ) +
    
    # Random scenario mean (on top of CI)
    geom_rect(
      aes(xmin = 0, xmax = mean_percent_random/100, 
          ymin = start_angle, ymax = end_angle,
          fill = paste0(taxon_group, "_random")),
      alpha = 1,
      color = "white",
      linewidth = 0.3
    ) +
    
    # Error bars - main horizontal line
    geom_segment(
      aes(x = radius_ci_lower, 
          xend = radius_ci_upper,
          y = mid_angle, 
          yend = mid_angle),
      color = "black",
      linewidth = 0.3,
      alpha = 1
    ) +
    
    # Error bar - left cap (radius-adjusted)
    geom_segment(
      aes(x = radius_ci_lower, 
          xend = radius_ci_lower,
          y = mid_angle - cap_offset_lower, 
          yend = mid_angle + cap_offset_lower),
      color = "black",
      linewidth = 0.3,
      alpha = 1
    ) +
    
    # Error bar - right cap (radius-adjusted)
    geom_segment(
      aes(x = radius_ci_upper, 
          xend = radius_ci_upper,
          y = mid_angle - cap_offset_upper, 
          yend = mid_angle + cap_offset_upper),
      color = "black",
      linewidth = 0.3,
      alpha = 1
    ) +
    
    # Percentage labels
    geom_text(
      data = data.frame(
        x = tick_labels / 100,
        label = paste0(tick_labels, "%")
      ),
      aes(x = x, y = 0, label = label),
      size = 4 / .pt,
      color = "black",
      vjust = -0.5,
      family = "Arial"
    ) +
    
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(
      colors_optimized,
      setNames(colors_random, paste0(names(colors_random), "_random")),
      setNames(colors_random, paste0(names(colors_random), "_random_ci"))
    )) +
    
    theme_void() +
    theme(
      legend.position = "none",
      text = element_text(family = "Arial", size = 4),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    xlim(0, 1.2)
}
# Make plots
SCo_data <- prepare_pollinator_plot_data("SCo",
                                         scenario1_taxon_summary,
                                         random_scenario1_taxon_summary)
plot_SCo <- plot_pollinator_support(SCo_data, "South Coast (SCo)")

CaRF_data <- prepare_pollinator_plot_data("CaRF",
                                          scenario1_taxon_summary,
                                          random_scenario1_taxon_summary)
plot_CaRF <- plot_pollinator_support(CaRF_data, "Central Coast Ranges (CaRF)")

# Save
ggsave("Figures/scenario1_pie/scenario1_SCo.png", 
       plot = plot_SCo,
       width = 35, 
       height = 35, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/scenario1_pie/scenario1_SCo.pdf", 
       plot = plot_SCo,
       width = 35, 
       height = 35, 
       units = "mm",
       device = cairo_pdf)

ggsave("Figures/scenario1_pie/scenario1_CaRF.png", 
       plot = plot_CaRF,
       width = 35, 
       height = 35, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/scenario1_pie/scenario1_CaRF.pdf", 
       plot = plot_CaRF,
       width = 35, 
       height = 35, 
       units = "mm",
       device = cairo_pdf)

## Summary #s
SCo_data %>%
  mutate(percent_increase = percent_supported - mean_percent_random) %>%
  mutate(percent_community_contr = percent_increase*proportion) %>%
  pull(percent_community_contr) %>%
  sum()

CaRF_data %>%
  mutate(percent_increase = percent_supported - mean_percent_random) %>%
  mutate(percent_community_contr = percent_increase*proportion) %>%
  pull(percent_community_contr) %>%
  sum()



#### Scenario 2: drought tolerant specialist bee garden ####
## Objective: Maximize oligolectic bee richness supported
## Constraints: Drought tolerant plants, nursery availability, known/presumed pollen plant for oligoleges

## Constrain plants
calscape_plants_list <- read.csv("Data_Clean/Species_Traits_Attributes/Calscape_Plants_Harmonized.csv")
calscape_plant_data <- read.csv("Data_Raw/Species_Traits_Attributes/plants_table_5-7.csv")

# add harmonized names to work with
calscape_plant_data <- left_join(calscape_plant_data, calscape_plants_list, by = "species")

# Constrain plant list
plantCandidates_scenario2 <- calscape_plant_data %>%
  # Constrain to drought tolerant plants ("low" in drought column = low water need)
  filter(grepl("low", drought, ignore.case = TRUE)) %>%
  # Must be available in at least one nursery
  filter(nurseries != "" & !is.na(nurseries))

# focal pollinators
# we will use Jarrod Fowler's specialist bee list to identify species that are oligolectic, and what their known/presumed host plants are
# but we'll use out spatially explicit flower visitation interaction network ("force" flower visitation to pollen collection for known host plants of oligoleges)
jarrodFowler_oligoleges <- read.csv("Data_Raw/Interactions/fowler/jarrodFowler_scraped_interactions.csv") %>%
  rename(sourceTaxonName = higher,
         targetTaxonName = lower)

## Harmonize scientific Names
# Load harmonization functions
source("Code/99_Supporting/name_harmonization.R")
# gnvpath
gnverifier_path <- "/home/tb625/bin"

# 1. Clean raw names
oligolege_names <- unique(jarrodFowler_oligoleges$sourceTaxonName)
plantInt_names  <- unique(jarrodFowler_oligoleges$targetTaxonName)

names_oligolege_parsed <- bdc_clean_names(sci_names = oligolege_names, save_outputs = FALSE)
names_plantInt_parsed  <- bdc_clean_names(sci_names = plantInt_names, save_outputs = FALSE)

# 2. Harmonize names (bees = Hymenoptera, plants = Plantae)
names_harm_oligolege <- harmonize_names(
  names = unique(names_oligolege_parsed$names_clean),
  higher_tax = "Hymenoptera",
  names_file = "Temp/names_oligoleges.tsv",
  names_harm_file = "Temp/names_harmonized_oligoleges.tsv",
  gnv_path = gnverifier_path
)
names_harm_plantInts <- harmonize_names(
  names = unique(names_plantInt_parsed$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_plantInts.tsv",
  names_harm_file = "Temp/names_harmonized_plantInts.tsv",
  gnv_path = gnverifier_path
)

# 3. Apply Jepson taxonomy (plants only)
# Skip here because jepson corrections are for species level, but data here are only at genus level
# names_harm_plantInts <- apply_jepson_taxonomy(names_harm_plantInts)

# Join back to original Jarrod Fowler list
oligoleges <- jarrodFowler_oligoleges %>%
  # Bees - First, with parsed names
  left_join(names_oligolege_parsed, by = join_by("sourceTaxonName" == "scientificName")) %>%
  rename(sourceTaxonName_parsed = names_clean) %>%
  dplyr::select(sourceTaxonName_parsed, targetTaxonName, interactionTypeName) %>%
  # Bees - Second, with harmonized names
  left_join(names_harm_oligolege, by = join_by("sourceTaxonName_parsed" == "ScientificName")) %>%
  rename(sourceTaxonName_harm = names_harm) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName, interactionTypeName) %>%
  # Plants - First, with parsed names
  left_join(names_plantInt_parsed, by = join_by("targetTaxonName" == "scientificName")) %>%
  rename(targetTaxonName_parsed = names_clean) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_parsed, interactionTypeName) %>%
  # Plants - Second, with harmonized names
  left_join(names_harm_plantInts, by = join_by("targetTaxonName_parsed" == "ScientificName")) %>%
  rename(targetTaxonName_harm = names_harm) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, interactionTypeName)

# remove rows with NA
oligoleges <- na.omit(oligoleges)
# now, clean plant genus column and filter to only CA species 
oligoleges_CA <- oligoleges %>%
  # force all plant names to genus - fowler data lists plant families and genera
  mutate(targetTaxonName_harm = word(targetTaxonName_harm, 1, 1)) %>%
  # filter to only CA pollinator species
  filter(sourceTaxonName_harm %in% checklist_cleaned$genus_species)

find_topMix_oligoleges <- function(ecoregion, 
                                   interaction_data_subset, 
                                   candidate_plants, 
                                   all_pollinators, 
                                   checklist,
                                   n_plants,
                                   oligoleges_pol,
                                   sp_by_region,
                                   optimization_mode = "oligolectic",
                                   n_generations = 100, 
                                   population_size = 200,
                                   mutation_rate = 0.1) {
  
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  # checklist of sp in this specific ecoregion
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  # Get unique plants and pollinators in this ecoregion
  plants <- candidate_plants[candidate_plants %in% checklist_region[checklist_region$taxon == "plants"]$genus_species] %>%
    unique()
  
  pollinators <- all_pollinators[all_pollinators %in% checklist_region[checklist_region$taxon %in% c("bees", "moths", "butterflies", "hoverflies")]$genus_species] %>%
    unique()
  
  # Get oligolectic bees from oligoleges_pol
  oligolectic_bees <- unique(oligoleges_pol$sourceTaxonName_harm)
  
  # Create a plant-pollinator interaction matrix (FULL VISITATION DATA)
  interaction_matrix_full <- ecoregion_data %>%
    dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  # Convert to logical matrix for more efficient operations
  interaction_matrix_full <- as.matrix(interaction_matrix_full) > 0
  
  # Filter FULL interaction matrix to include only plants and pollinators of interest
  interaction_matrix_full_filtered <- interaction_matrix_full[
    rownames(interaction_matrix_full) %in% plants, 
    colnames(interaction_matrix_full) %in% pollinators,
    drop = FALSE
  ]
  
  # Create OLIGOLECTIC POLLEN interaction matrix
  oligolectic_mask <- interaction_matrix_full_filtered
  oligolectic_mask[,] <- FALSE  # Initialize all as FALSE
  
  # For each cell in the matrix, check if it's a valid oligolectic interaction
  for (plant_species in rownames(interaction_matrix_full_filtered)) {
    # Extract genus from plant species name (first word)
    plant_genus <- sub("\\s+.*", "", plant_species)
    
    for (bee_species in colnames(interaction_matrix_full_filtered)) {
      # Check if this bee-plant genus combo exists in oligoleges_pol
      is_oligolectic <- any(
        oligoleges_pol$sourceTaxonName_harm == bee_species & 
          oligoleges_pol$targetTaxonName_harm == plant_genus &
          oligoleges_pol$interactionTypeName == "collectsPollenOf"
      )
      
      if (is_oligolectic) {
        oligolectic_mask[plant_species, bee_species] <- TRUE
      }
    }
  }
  
  # Apply the mask to create oligolectic-only matrix
  interaction_matrix_oligolectic <- interaction_matrix_full_filtered & oligolectic_mask
  
  # Filter to only include oligolectic bees that have at least one pollen interaction
  oligolectic_bees_in_region <- oligolectic_bees[oligolectic_bees %in% colnames(interaction_matrix_oligolectic)]
  interaction_matrix_oligolectic <- interaction_matrix_oligolectic[
    , 
    colnames(interaction_matrix_oligolectic) %in% oligolectic_bees_in_region,
    drop = FALSE
  ]
  
  # Choose which matrix to use for optimization
  if (optimization_mode == "oligolectic") {
    optimization_matrix <- interaction_matrix_oligolectic
    message(paste("Optimizing for oligolectic bee pollen interactions in", ecoregion))
  } else if (optimization_mode == "all_pollinators") {
    optimization_matrix <- interaction_matrix_full_filtered
    message(paste("Optimizing for all pollinator interactions in", ecoregion))
  } else {
    stop("optimization_mode must be either 'oligolectic' or 'all_pollinators'")
  }
  
  # Check if we have enough plants
  if (nrow(optimization_matrix) < n_plants) {
    warning(paste("Only", nrow(optimization_matrix), 
                  "plants available in ecoregion", ecoregion))
    n_plants <- nrow(optimization_matrix)
  }
  
  all_plants <- rownames(optimization_matrix)
  
  message(paste("Using genetic algorithm with", length(all_plants), "candidate plants..."))
  
  # Fitness function
  calculate_fitness <- function(plant_combo) {
    covered_pollinators <- apply(optimization_matrix[plant_combo, , drop = FALSE], 
                                 2, any)
    sum(covered_pollinators)
  }
  
  # Helper function to create a unique key for a combination
  combo_key <- function(plant_combo) {
    paste(sort(plant_combo), collapse = "|||")
  }
  
  # Initialize population with random plant combinations
  population <- replicate(population_size, 
                          sample(all_plants, n_plants), 
                          simplify = FALSE)
  
  # Track best solutions
  best_fitness <- -Inf
  best_solutions <- list()
  fitness_history <- numeric(n_generations)
  
  # Genetic algorithm main loop
  for (gen in 1:n_generations) {
    # Calculate fitness for each individual
    fitness_scores <- sapply(population, calculate_fitness)
    
    # Update best solutions
    current_max <- max(fitness_scores)
    
    if (current_max > best_fitness) {
      best_fitness <- current_max
      best_solutions <- list()
      
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        best_solutions[[key]] <- population[[idx]]
      }
    } else if (current_max == best_fitness) {
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        if (!key %in% names(best_solutions)) {
          best_solutions[[key]] <- population[[idx]]
        }
      }
    }
    
    fitness_history[gen] <- best_fitness
    
    # Selection: tournament selection
    select_parent <- function() {
      tournament_size <- 5
      tournament_indices <- sample(length(population), tournament_size)
      tournament_fitness <- fitness_scores[tournament_indices]
      population[[tournament_indices[which.max(tournament_fitness)]]]
    }
    
    # Create next generation
    new_population <- list()
    
    # Elitism: keep best 10% of population
    n_elite <- max(1, floor(population_size * 0.1))
    elite_indices <- order(fitness_scores, decreasing = TRUE)[1:n_elite]
    new_population[1:n_elite] <- population[elite_indices]
    
    # Generate rest through crossover and mutation
    for (i in (n_elite + 1):population_size) {
      parent1 <- select_parent()
      parent2 <- select_parent()
      
      # Crossover
      n_from_parent1 <- sample(1:(n_plants-1), 1)
      child <- c(sample(parent1, n_from_parent1), 
                 sample(parent2, n_plants - n_from_parent1))
      
      # Ensure no duplicates
      child <- unique(child)
      while (length(child) < n_plants) {
        child <- c(child, sample(setdiff(all_plants, child), 1))
      }
      child <- child[1:n_plants]
      
      # Mutation
      if (runif(1) < mutation_rate) {
        n_mutations <- sample(1:2, 1)
        positions <- sample(n_plants, n_mutations)
        child[positions] <- sample(setdiff(all_plants, child), n_mutations)
      }
      
      new_population[[i]] <- child
    }
    
    population <- new_population
    
    # Progress update every 10 generations
    if (gen %% 10 == 0) {
      message(paste("Generation", gen, "- Best fitness:", best_fitness, 
                    "- Unique solutions found:", length(best_solutions)))
    }
  }
  
  message(paste("Best combination(s) support", best_fitness, "pollinators"))
  message(paste("Found", length(best_solutions), "unique combination(s) with this score"))
  
  # Select best solution (first one if multiple)
  best_solution <- best_solutions[[1]]
  
  # Calculate final coverage using FULL visitation matrix
  covered_pollinators_full_final <- apply(
    interaction_matrix_full_filtered[best_solution, , drop = FALSE], 
    2, 
    any
  )
  
  # Calculate final coverage using OLIGOLECTIC matrix
  covered_oligolectic_final <- apply(
    interaction_matrix_oligolectic[best_solution, , drop = FALSE],
    2,
    any
  )
  
  # Calculate taxon group statistics
  all_pollinators_in_region <- checklist_region %>%
    filter(taxon %in% c("bees", "moths", "butterflies", "hoverflies"))
  
  supported_pollinator_names <- names(covered_pollinators_full_final)[covered_pollinators_full_final]
  
  # Calculate total counts by taxon group in the region
  taxon_totals <- all_pollinators_in_region %>%
    group_by(taxon) %>%
    summarise(
      total_in_region = n(),
      .groups = 'drop'
    )
  
  # Calculate supported counts by taxon group
  taxon_supported <- all_pollinators_in_region %>%
    filter(genus_species %in% supported_pollinator_names) %>%
    group_by(taxon) %>%
    summarise(
      supported_by_mix = n(),
      .groups = 'drop'
    )
  
  # Combine and calculate percentages
  taxon_summary <- taxon_totals %>%
    left_join(taxon_supported, by = "taxon") %>%
    mutate(
      supported_by_mix = replace_na(supported_by_mix, 0),
      percent_supported = round(100 * supported_by_mix / total_in_region, 1),
      ecoregion = ecoregion
    ) %>%
    dplyr::select(ecoregion, taxon, total_in_region, supported_by_mix, percent_supported)
  
  # Calculate oligolectic bee summary
  total_oligolectic_bees <- length(oligolectic_bees_in_region)
  supported_oligolectic_bees <- sum(covered_oligolectic_final)
  percent_oligolectic_supported <- round(100 * supported_oligolectic_bees / total_oligolectic_bees, 1)
  
  oligolectic_summary <- data.frame(
    ecoregion = ecoregion,
    optimization_mode = optimization_mode,
    total_oligolectic_bees = total_oligolectic_bees,
    supported_oligolectic_bees = supported_oligolectic_bees,
    percent_oligolectic_supported = percent_oligolectic_supported
  )
  
  # Print summaries
  message(paste("\n=== Taxon Group Summary for", ecoregion, "(Full Visitation Data) ==="))
  print(taxon_summary)
  
  message(paste("\n=== Oligolectic Bee Summary for", ecoregion, "==="))
  print(oligolectic_summary)
  
  # Create results dataframe
  top_mix <- data.frame(
    ecoregion = ecoregion,
    rank = 1:n_plants,
    scientificName = best_solution,
    cumulative_target_pollinators = best_fitness,
    cumulative_oligolectic_bees = sum(covered_oligolectic_final),
    cumulative_all_pollinators = sum(covered_pollinators_full_final)
  )
  
  # Add attributes
  attr(top_mix, "all_best_combinations") <- best_solutions
  attr(top_mix, "n_best_combinations") <- length(best_solutions)
  attr(top_mix, "fitness_history") <- fitness_history
  attr(top_mix, "final_fitness") <- best_fitness
  
  # Return list with all tables
  return(list(
    top_mix = top_mix,
    taxon_summary = taxon_summary,
    oligolectic_summary = oligolectic_summary
  ))
}

sample_random_mix <- function(ecoregion, 
                              interaction_data_subset, 
                              candidate_plants, 
                              checklist,
                              n_plants,
                              oligoleges_pol,
                              sp_by_region,
                              n_iterations = 100) {
  
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  # checklist of sp in this specific ecoregion
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  # Get oligolectic bees
  oligolectic_bees <- unique(oligoleges_pol$sourceTaxonName_harm)
  
  # Create interaction matrix
  interaction_matrix_full <- ecoregion_data %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  interaction_matrix_full <- as.matrix(interaction_matrix_full) > 0
  
  # Create oligolectic matrix
  oligolectic_mask <- interaction_matrix_full
  oligolectic_mask[,] <- FALSE
  
  for (plant_species in rownames(interaction_matrix_full)) {
    plant_genus <- sub("\\s+.*", "", plant_species)
    for (bee_species in colnames(interaction_matrix_full)) {
      is_oligolectic <- any(
        oligoleges_pol$sourceTaxonName_harm == bee_species & 
          oligoleges_pol$targetTaxonName_harm == plant_genus &
          oligoleges_pol$interactionTypeName == "collectsPollenOf"
      )
      if (is_oligolectic) {
        oligolectic_mask[plant_species, bee_species] <- TRUE
      }
    }
  }
  
  interaction_matrix_oligolectic <- interaction_matrix_full & oligolectic_mask
  
  # Filter to plants that exist in both matrices
  available_plants <- intersect(rownames(interaction_matrix_oligolectic), candidate_plants)
  oligolectic_bees_in_region <- oligolectic_bees[oligolectic_bees %in% colnames(interaction_matrix_oligolectic)]
  
  total_oligolectic_bees <- length(oligolectic_bees_in_region)
  
  # Run iterations
  results <- replicate(n_iterations, {
    selected_plants <- sample(available_plants, min(n_plants, length(available_plants)), replace = FALSE)
    
    # Calculate coverage
    covered_oligo <- apply(
      interaction_matrix_oligolectic[selected_plants, oligolectic_bees_in_region, drop = FALSE], 
      2, 
      any
    )
    
    supported_oligolectic_bees <- sum(covered_oligo)
    percent_oligolectic_supported <- 100 * supported_oligolectic_bees / total_oligolectic_bees
    
    data.frame(
      ecoregion = ecoregion,
      optimization_mode = "random",
      total_oligolectic_bees = total_oligolectic_bees,
      supported_oligolectic_bees = supported_oligolectic_bees,
      percent_oligolectic_supported = percent_oligolectic_supported
    )
  }, simplify = FALSE)
  
  do.call(rbind, results)
}


# V1: optimize oligolege diversity supported (based on presumed pollen collection interactions)
scenario2_results_list <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_topMix_oligoleges(
    ecoregion, 
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants = unique(plantCandidates_scenario2$genus_species_harm),
    all_pollinators = unique(checklist_cleaned["taxon" != "plants",]$scientificName),
    checklist = checklist_cleaned,
    n_plants = n_plants,
    oligoleges_pol = oligoleges_CA,
    sp_by_region = sp_by_region,
    optimization_mode = "oligolectic",
    n_generations = 200, 
    population_size = 200,
    mutation_rate = 0.1
  )
},
mc.cores = 12)

# V2: optimized by all pols
scenario2_results_list_allPols <- mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_topMix_oligoleges(
    ecoregion, 
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants = unique(plantCandidates_scenario2$genus_species_harm),
    all_pollinators = unique(checklist_cleaned["taxon" != "plants",]$scientificName),
    checklist = checklist_cleaned,
    n_plants = n_plants,
    oligoleges_pol = oligoleges_CA,
    sp_by_region = sp_by_region,
    optimization_mode = "all_pollinators",
    n_generations = 200, 
    population_size = 200,
    mutation_rate = 0.1
  )
},
mc.cores = 12)

# V3: 
# Run 100 random iterations for each ecoregion
random_scenario2_oligolege_summary <- do.call(rbind, mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  sample_random_mix(
    ecoregion, 
    interaction_data_subset = filter(predicted_interactions, interactionTypeName == "visitsFlowersOf"),
    candidate_plants = unique(plantCandidates_scenario2$genus_species_harm),
    checklist = checklist_cleaned,
    n_plants = n_plants,
    oligoleges_pol = oligoleges_CA,
    sp_by_region = sp_by_region,
    n_iterations = 100
  )
},
mc.cores = 12))


# Optimizing for oligolege support
scenario2_results <- do.call(rbind, lapply(scenario2_results_list, function(x) x$top_mix))
scenario2_taxon_summary <- do.call(rbind, lapply(scenario2_results_list, function(x) x$taxon_summary))
scenario2_oligolege_summary <- do.call(rbind, lapply(scenario2_results_list, function(x) x$oligolectic_summary))

# Optimizing for all pols
scenario2_results_allPols <- do.call(rbind, lapply(scenario2_results_list_allPols, function(x) x$top_mix))
scenario2_taxon_summary_allPols <- do.call(rbind, lapply(scenario2_results_list_allPols, function(x) x$taxon_summary))
scenario2_oligolege_summary_allPols <- do.call(rbind, lapply(scenario2_results_list_allPols, function(x) x$oligolectic_summary))

# random % supported oligoleges
random_scenario2_oligolege_summary %>% group_by(ecoregion) %>% summarise(perc_random = mean(percent_oligolectic_supported))

# Save plant lists
write.csv(scenario2_results, "Data_Clean/Scenario_Results/Scenario2_oligolegeOpt.csv")
write.csv(scenario2_results_allPols, "Data_Clean/Scenario_Results/Scenario2_allPolOpt.csv")

write.csv(scenario2_results %>% group_by(ecoregion) %>% slice_head(n=1) %>% dplyr::select(-c(rank, scientificName)), 
          "Data_Clean/Scenario_Results/Scenario2_oligolegeOpt_pols.csv")
write.csv(scenario2_results_allPols %>% group_by(ecoregion) %>% slice_head(n=1) %>% dplyr::select(-c(rank, scientificName)), 
          "Data_Clean/Scenario_Results/Scenario2_allPolOpt_pols.csv")


## Plot
plot_oligolege_support <- function(optimized_summary, random_summary, ecoregion_code) {
  
  # Calculate random stats for this ecoregion
  random_stats <- random_summary %>%
    filter(ecoregion == ecoregion_code) %>%
    summarise(
      mean_percent = mean(percent_oligolectic_supported),
      sd_percent = sd(percent_oligolectic_supported),
      se_percent = sd(percent_oligolectic_supported) / sqrt(n()),
      ci_lower = mean(percent_oligolectic_supported) - 1.96 * se_percent,
      ci_upper = mean(percent_oligolectic_supported) + 1.96 * se_percent,
      .groups = 'drop'
    ) %>%
    mutate(
      ecoregion = ecoregion_code,
      optimization_mode = "random"
    )
  
  # Prepare optimized data for this ecoregion
  optimized_data <- optimized_summary %>%
    filter(ecoregion == ecoregion_code) %>%
    mutate(
      mean_percent = percent_oligolectic_supported,
      sd_percent = 0,
      se_percent = 0,
      ci_lower = NA,  # No error bars for optimized
      ci_upper = NA
    ) %>%
    dplyr::select(ecoregion, optimization_mode, mean_percent, sd_percent, 
                  se_percent, ci_lower, ci_upper)
  
  # Combine data
  plot_data <- bind_rows(optimized_data, random_stats) %>%
    mutate(
      optimization_mode = factor(optimization_mode, 
                                 levels = c("random", "all_pollinators", "oligolectic"))
    )
  
  # Create plot
  ggplot(plot_data, aes(x = optimization_mode, y = mean_percent, fill = optimization_mode)) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "random" = "#FFE0A2",
        "all_pollinators" = "#FFCC00",
        "oligolectic" = "#A98A10"
      )
    ) +
    # scale_x_discrete(
    #   labels = c(
    #     "random" = "Random",
    #     "all_pollinators" = "Optimize\nPollinator\nDiversity",
    #     "oligolectic" = "Optimize\nOligolege\nDiversity"
    #   )
    # ) +
    scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
    labs(
      x = NULL,
      # y = "% Oligoleges Supported\n(Pollen Interaction)"
      y = NULL
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", size = 4),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      # axis.text.x = element_text(angle = 0, hjust = 0.5, size = 4),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 4),
      plot.margin = ggplot2::margin(2, 2, 2, 2, "mm")
    )
}

# Create and save plots
plot_SCo <- plot_oligolege_support(
  bind_rows(scenario2_oligolege_summary, scenario2_oligolege_summary_allPols),
  random_scenario2_oligolege_summary,
  "SCo"
)

plot_CaRF <- plot_oligolege_support(
  bind_rows(scenario2_oligolege_summary, scenario2_oligolege_summary_allPols),
  random_scenario2_oligolege_summary,
  "CaRF"
)

# Save
ggsave("Figures/scenario2_bar/scen2_SCo.png", 
       plot_SCo, 
       width = 30, 
       height = 30, 
       units = "mm", 
       dpi = 2000)
ggsave("Figures/scenario2_bar/scen2_SCo.pdf", 
       plot = plot_SCo,
       width = 30, 
       height = 30, 
       units = "mm", 
       device = cairo_pdf)

ggsave("Figures/scenario2_bar/scen2_CaRF.png", 
       plot_CaRF, 
       width = 30, 
       height = 30, 
       units = "mm", 
       dpi = 2000)
ggsave("Figures/scenario2_bar/scen2_CaRF.pdf", 
       plot = plot_CaRF,
       width = 30, 
       height = 30, 
       units = "mm", 
       device = cairo_pdf)


#### Scenario 3: supporting interactions for impereled species ####
## Objective: Maximize butterfly richness supported, among known declining butterflies
## Constraints: butterfly needs host+nectar plant

# We will use Edwards et al (2025) Science to get all declining butterfly species in SW
# and identify planting mix that best supports them
butterfly_trends_SW <- read_excel("Data_Raw/Edwards2025/science.adp4671_tables_s1_to_s10.xlsx", sheet = "Table S9") %>%
  filter(Region == "Southwest")
# butterfly_trends <- read_excel("Data_Raw/Edwards2025/science.adp4671_tables_s1_to_s10.xlsx", sheet = "Table S5")

butterflies_declining_SW <- butterfly_trends_SW %>%
  filter(Trend == "declining") %>%
  # standardize column names
  rename(scientificName = `Scientific name`)

## Harmonize scientific Names
# Load harmonization functions
source("Code/99_Supporting/name_harmonization.R")

# 1. Clean raw names
butterfly_names <- unique(butterflies_declining_SW$scientificName)
names_butterfly_parsed <- bdc_clean_names(sci_names = butterfly_names, save_outputs = FALSE)

# 2. Harmonize names (Lepidoptera)
names_harm_butterfly <- harmonize_names(
  names = unique(names_butterfly_parsed$names_clean),
  higher_tax = "Lepidoptera",
  names_file = "Temp/names_butterflies_declining.tsv",
  names_harm_file = "Temp/names_harmonized_butterflies_declining.tsv",
  gnv_path = gnverifier_path
)

# list of all focal pols
butterflies_declining_SW_list <- unique(pull(names_harm_butterfly, names_harm))

# plant candidates
plant_candidates_scenario3 <- checklist_cleaned %>%
  filter(taxon == "plants") %>%
  pull(genus_species) %>%
  unique()

find_topMix_butterfly <- function(ecoregion, 
                                  interaction_data_subset, 
                                  candidate_plants, 
                                  focal_pollinators, 
                                  n_plants,
                                  checklist,
                                  n_generations = 100,
                                  population_size = 200,
                                  mutation_rate = 0.1) {
  
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  # checklist of sp in this specific ecoregion
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  # Get unique plants and pollinators in this ecoregion
  plants <- candidate_plants[candidate_plants %in% checklist_region[checklist_region$taxon == "plants"]$genus_species] %>%
    unique()
  
  pollinators <- focal_pollinators[focal_pollinators %in% checklist_region[checklist_region$taxon %in% c("butterflies")]$genus_species] %>%
    unique()
  
  # Create separate matrices for host plant and nectar interactions
  host_matrix <- ecoregion_data %>%
    filter(interactionTypeName == "hasHost") %>%
    dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  nectar_matrix <- ecoregion_data %>%
    filter(interactionTypeName == "visitsFlowersOf") %>%
    dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fn = max,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  # Convert to logical matrices
  host_matrix <- as.matrix(host_matrix) > 0
  nectar_matrix <- as.matrix(nectar_matrix) > 0
  
  # Filter matrices to include only plants and pollinators of interest
  host_matrix_filtered <- host_matrix[
    rownames(host_matrix) %in% plants, 
    colnames(host_matrix) %in% pollinators,
    drop = FALSE
  ]
  
  nectar_matrix_filtered <- nectar_matrix[
    rownames(nectar_matrix) %in% plants, 
    colnames(nectar_matrix) %in% pollinators,
    drop = FALSE
  ]
  
  # Ensure both matrices have the same pollinator and plant columns
  all_plants <- union(rownames(host_matrix_filtered), rownames(nectar_matrix_filtered))
  all_pols <- union(colnames(host_matrix_filtered), colnames(nectar_matrix_filtered))
  
  # Expand matrices to include all plants and pollinators
  expand_matrix <- function(mat, all_rows, all_cols) {
    # Add missing columns
    missing_cols <- setdiff(all_cols, colnames(mat))
    if (length(missing_cols) > 0) {
      missing_mat <- matrix(FALSE, nrow = nrow(mat), ncol = length(missing_cols))
      colnames(missing_mat) <- missing_cols
      rownames(missing_mat) <- rownames(mat)
      mat <- cbind(mat, missing_mat)
    }
    
    # Add missing rows
    missing_rows <- setdiff(all_rows, rownames(mat))
    if (length(missing_rows) > 0) {
      missing_mat <- matrix(FALSE, nrow = length(missing_rows), ncol = ncol(mat))
      rownames(missing_mat) <- missing_rows
      colnames(missing_mat) <- colnames(mat)
      mat <- rbind(mat, missing_mat)
    }
    
    # Reorder to match specified order
    mat[all_rows, all_cols, drop = FALSE]
  }
  
  host_matrix_filtered <- expand_matrix(host_matrix_filtered, all_plants, all_pols)
  nectar_matrix_filtered <- expand_matrix(nectar_matrix_filtered, all_plants, all_pols)
  
  # Check plants available for use
  if (length(all_plants) < n_plants) {
    warning(paste("Only", length(all_plants), "plants available in ecoregion", ecoregion))
    n_plants <- length(all_plants)
  }
  
  message(paste("Using genetic algorithm with", length(all_plants), "candidate plants..."))
  
  # Fitness function - requires distinct host and nectar plants for but to be supported
  calculate_fitness <- function(plant_combo) {
    supported <- sapply(all_pols, function(lep) {
      # Get plants that serve as host for this butterfly
      host_plants <- plant_combo[host_matrix_filtered[plant_combo, lep, drop = TRUE]]
      # Get plants that serve as nectar for this butterfly
      nectar_plants <- plant_combo[nectar_matrix_filtered[plant_combo, lep, drop = TRUE]]
      
      # Check if there's at least one host plant AND at least one DIFFERENT nectar plant
      has_host <- length(host_plants) > 0
      has_nectar <- length(nectar_plants) > 0
      has_distinct <- length(setdiff(nectar_plants, host_plants)) > 0
      
      has_host && has_nectar && has_distinct
    })
    sum(supported)
  }
  
  # Helper function to create a unique key for a combination
  combo_key <- function(plant_combo) {
    paste(sort(plant_combo), collapse = "|||")
  }
  
  # Initialize population with random plant combinations
  population <- replicate(population_size, 
                          sample(all_plants, n_plants), 
                          simplify = FALSE)
  
  # Track best solutions - use a named list where names are combo keys
  best_fitness <- -Inf
  best_solutions <- list()  # Will store all combinations with max fitness
  fitness_history <- numeric(n_generations)
  
  # Genetic algorithm main loop
  for (gen in 1:n_generations) {
    # Calculate fitness for each individual
    fitness_scores <- sapply(population, calculate_fitness)
    
    # Update best solutions
    current_max <- max(fitness_scores)
    
    if (current_max > best_fitness) {
      # Found a new best - reset the list
      best_fitness <- current_max
      best_solutions <- list()
      
      # Add all solutions with this fitness
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        best_solutions[[key]] <- population[[idx]]
      }
    } else if (current_max == best_fitness) {
      # Found more solutions with the same max fitness
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        if (!key %in% names(best_solutions)) {
          best_solutions[[key]] <- population[[idx]]
        }
      }
    }
    
    fitness_history[gen] <- best_fitness
    
    # Selection: tournament selection
    select_parent <- function() {
      tournament_size <- 5
      tournament_indices <- sample(length(population), tournament_size)
      tournament_fitness <- fitness_scores[tournament_indices]
      population[[tournament_indices[which.max(tournament_fitness)]]]
    }
    
    # Create next generation
    new_population <- list()
    
    # Elitism: keep best 10% of population
    n_elite <- max(1, floor(population_size * 0.1))
    elite_indices <- order(fitness_scores, decreasing = TRUE)[1:n_elite]
    new_population[1:n_elite] <- population[elite_indices]
    
    # Generate rest through crossover and mutation
    for (i in (n_elite + 1):population_size) {
      parent1 <- select_parent()
      parent2 <- select_parent()
      
      # Crossover: combine parents
      n_from_parent1 <- sample(1:(n_plants-1), 1)
      child <- c(sample(parent1, n_from_parent1), 
                 sample(parent2, n_plants - n_from_parent1))
      
      # Ensure no duplicates
      child <- unique(child)
      while (length(child) < n_plants) {
        child <- c(child, sample(setdiff(all_plants, child), 1))
      }
      child <- child[1:n_plants]
      
      # Mutation: replace random plants
      if (runif(1) < mutation_rate) {
        n_mutations <- sample(1:2, 1)
        positions <- sample(n_plants, n_mutations)
        child[positions] <- sample(setdiff(all_plants, child), n_mutations)
      }
      
      new_population[[i]] <- child
    }
    
    population <- new_population
    
    # Progress update every 10 generations
    if (gen %% 10 == 0) {
      message(paste("Generation", gen, "- Best fitness:", best_fitness, 
                    "- Unique solutions found:", length(best_solutions)))
    }
  }
  
  message(paste("Best combination(s) support", best_fitness, "lepidoptera species"))
  message(paste("Found", length(best_solutions), "unique combination(s) with this score"))
  message(paste("Selecting best mix from", length(best_solutions), "tied combinations..."))
  
  # Function to count total pollinators supported by a plant combination
  count_total_pollinators <- function(plant_combo) {
    # Get all interactions for these plants in this ecoregion
    plant_interactions <- ecoregion_data %>%
      filter(targetTaxonName_harm %in% plant_combo)
    
    # Count unique pollinators (not just focal ones)
    n_pollinators <- length(unique(plant_interactions$sourceTaxonName_harm))
    
    return(n_pollinators)
  }
  
  # Evaluate all best solutions
  pollinator_counts <- sapply(best_solutions, count_total_pollinators)
  
  # Select the combination(s) with most total pollinators
  max_pollinators <- max(pollinator_counts)
  best_idx <- which.max(pollinator_counts)
  best_solution <- best_solutions[[best_idx]]
  
  message(paste("Selected mix supports", max_pollinators, "total pollinator species"))
  message(paste("(vs focal butterflies:", best_fitness, ")"))
  
  
  # Calculate interaction support for butterflies
  butterfly_support <- data.frame(
    butterfly = all_pols,
    has_host = sapply(all_pols, function(lep) {
      any(host_matrix_filtered[best_solution, lep, drop = TRUE], na.rm = TRUE)
    }),
    has_nectar = sapply(all_pols, function(lep) {
      any(nectar_matrix_filtered[best_solution, lep, drop = TRUE], na.rm = TRUE)
    })
  )
  
  # Add check for distinct plants
  butterfly_support <- butterfly_support %>%
    mutate(
      has_distinct_nectar = sapply(all_pols, function(lep) {
        host_plants <- best_solution[host_matrix_filtered[best_solution, lep, drop = TRUE]]
        nectar_plants <- best_solution[nectar_matrix_filtered[best_solution, lep, drop = TRUE]]
        length(setdiff(nectar_plants, host_plants)) > 0
      })
    ) %>%
    mutate(
      support_category = case_when(
        has_host & has_nectar & has_distinct_nectar ~ "both",
        has_host & !has_nectar ~ "host_only",
        has_host & has_nectar & !has_distinct_nectar ~ "host_only",  # Same plant = host only
        !has_host & has_nectar ~ "nectar_only",
        TRUE ~ "none"
      )
    )
  
  # Count butterflies in each category
  n_both <- sum(butterfly_support$support_category == "both")
  n_host_only <- sum(butterfly_support$support_category == "host_only")
  n_nectar_only <- sum(butterfly_support$support_category == "nectar_only")
  
  # Get lists of butterflies in each category
  leps_both <- butterfly_support$butterfly[butterfly_support$support_category == "both"]
  leps_host_only <- butterfly_support$butterfly[butterfly_support$support_category == "host_only"]
  leps_nectar_only <- butterfly_support$butterfly[butterfly_support$support_category == "nectar_only"]
  
  message(paste("\nSupport breakdown (requires distinct host and nectar plants):"))
  message(paste("  Both host and nectar:", n_both))
  message(paste("  Host only:", n_host_only))
  message(paste("  Nectar only:", n_nectar_only))
  
  # Create detailed results
  result <- data.frame(
    ecoregion = ecoregion,
    rank = 1:n_plants,
    scientificName = best_solution,
    total_leps_supported = best_fitness,
    n_both_interactions = n_both,
    n_host_only = n_host_only,
    n_nectar_only = n_nectar_only
  )
  
  # which lepidoptera are supported
  result$supported_leps_list <- list(leps_both)
  
  # Add category-specific lists
  result$leps_both <- list(leps_both)
  result$leps_host_only <- list(leps_host_only)
  result$leps_nectar_only <- list(leps_nectar_only)
  
  # Add all best combinations as attribute
  attr(result, "all_best_combinations") <- best_solutions
  attr(result, "n_best_combinations") <- length(best_solutions)
  attr(result, "fitness_history") <- fitness_history
  attr(result, "final_fitness") <- best_fitness
  attr(result, "butterfly_support_details") <- butterfly_support
  
  return(result)
}

scenario3_results <- do.call(rbind, mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_topMix_butterfly(
    ecoregion, 
    interaction_data_subset = predicted_interactions,
    candidate_plants = plant_candidates_scenario3,
    focal_pollinators = butterflies_declining_SW_list,
    n_plants = n_plants,
    checklist = checklist_cleaned,
    n_generations = 100,
    population_size = 500,
    mutation_rate = 0.1
  )
},
mc.cores = 12))

## Compare to random draw
find_randomMix_butterflies <- function(ecoregion,
                                       interaction_data_subset, 
                                       candidate_plants, 
                                       focal_pollinators, 
                                       n_plants,
                                       checklist,
                                       n_iterations = 100) {
  
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  # checklist of sp in this specific ecoregion
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  # Get unique plants and pollinators in this ecoregion
  plants <- candidate_plants[candidate_plants %in% checklist_region[checklist_region$taxon == "plants"]$genus_species] %>%
    unique()
  
  pollinators <- focal_pollinators[focal_pollinators %in% checklist_region[checklist_region$taxon %in% c("butterflies")]$genus_species] %>%
    unique()
  
  # Create separate matrices for host plant and nectar interactions
  host_matrix <- ecoregion_data %>%
    filter(interactionTypeName == "hasHost") %>%
    dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  nectar_matrix <- ecoregion_data %>%
    filter(interactionTypeName == "visitsFlowersOf") %>%
    dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, potentialOrConfirmed) %>%
    distinct() %>%
    pivot_wider(names_from = sourceTaxonName_harm, 
                values_from = potentialOrConfirmed,
                values_fn = max,
                values_fill = 0) %>%
    tibble::column_to_rownames("targetTaxonName_harm")
  
  # Convert to logical matrices
  host_matrix <- as.matrix(host_matrix) > 0
  nectar_matrix <- as.matrix(nectar_matrix) > 0
  
  # Filter matrices to include only plants and pollinators of interest
  host_matrix_filtered <- host_matrix[
    rownames(host_matrix) %in% plants, 
    colnames(host_matrix) %in% pollinators,
    drop = FALSE
  ]
  
  nectar_matrix_filtered <- nectar_matrix[
    rownames(nectar_matrix) %in% plants, 
    colnames(nectar_matrix) %in% pollinators,
    drop = FALSE
  ]
  
  # Ensure both matrices have the same pollinator and plant columns
  all_plants <- union(rownames(host_matrix_filtered), rownames(nectar_matrix_filtered))
  all_pols <- union(colnames(host_matrix_filtered), colnames(nectar_matrix_filtered))
  
  # Expand matrices to include all plants and pollinators
  expand_matrix <- function(mat, all_rows, all_cols) {
    # Add missing columns
    missing_cols <- setdiff(all_cols, colnames(mat))
    if (length(missing_cols) > 0) {
      missing_mat <- matrix(FALSE, nrow = nrow(mat), ncol = length(missing_cols))
      colnames(missing_mat) <- missing_cols
      rownames(missing_mat) <- rownames(mat)
      mat <- cbind(mat, missing_mat)
    }
    
    # Add missing rows
    missing_rows <- setdiff(all_rows, rownames(mat))
    if (length(missing_rows) > 0) {
      missing_mat <- matrix(FALSE, nrow = length(missing_rows), ncol = ncol(mat))
      rownames(missing_mat) <- missing_rows
      colnames(missing_mat) <- colnames(mat)
      mat <- rbind(mat, missing_mat)
    }
    
    # Reorder to match specified order
    mat[all_rows, all_cols, drop = FALSE]
  }
  
  host_matrix_filtered <- expand_matrix(host_matrix_filtered, all_plants, all_pols)
  nectar_matrix_filtered <- expand_matrix(nectar_matrix_filtered, all_plants, all_pols)
  
  # Check if we have enough plants
  if (length(all_plants) < n_plants) {
    warning(paste("Only", length(all_plants), 
                  "plants available in ecoregion", ecoregion))
    return(data.frame(
      ecoregion = ecoregion,
      iteration = NA,
      n_plants = length(all_plants),
      butterflies_supported = NA,
      n_both_interactions = NA,
      n_host_only = NA,
      n_nectar_only = NA
    ))
  }
  
  # Fitness function - requires distinct host and nectar plants for but to be supported
  calculate_fitness <- function(plant_combo) {
    supported <- sapply(all_pols, function(lep) {
      # Get plants that serve as host for this butterfly
      host_plants <- plant_combo[host_matrix_filtered[plant_combo, lep, drop = TRUE]]
      # Get plants that serve as nectar for this butterfly
      nectar_plants <- plant_combo[nectar_matrix_filtered[plant_combo, lep, drop = TRUE]]
      
      # Check if there's at least one host plant AND at least one DIFFERENT nectar plant
      has_host <- length(host_plants) > 0
      has_nectar <- length(nectar_plants) > 0
      has_distinct <- length(setdiff(nectar_plants, host_plants)) > 0
      
      has_host && has_nectar && has_distinct
    })
    sum(supported)
  }
  
  #  calculate detailed interaction breakdown with distinct plant requirement
  calculate_interaction_breakdown <- function(plant_combo) {
    butterfly_support <- data.frame(
      butterfly = all_pols,
      has_host = sapply(all_pols, function(lep) {
        any(host_matrix_filtered[plant_combo, lep, drop = TRUE], na.rm = TRUE)
      }),
      has_nectar = sapply(all_pols, function(lep) {
        any(nectar_matrix_filtered[plant_combo, lep, drop = TRUE], na.rm = TRUE)
      })
    )
    
    # Add check for distinct plants
    butterfly_support <- butterfly_support %>%
      mutate(
        has_distinct_nectar = sapply(all_pols, function(lep) {
          host_plants <- plant_combo[host_matrix_filtered[plant_combo, lep, drop = TRUE]]
          nectar_plants <- plant_combo[nectar_matrix_filtered[plant_combo, lep, drop = TRUE]]
          length(setdiff(nectar_plants, host_plants)) > 0
        })
      ) %>%
      mutate(
        support_category = case_when(
          has_host & has_nectar & has_distinct_nectar ~ "both",
          has_host & !has_nectar ~ "host_only",
          has_host & has_nectar & !has_distinct_nectar ~ "host_only",  # Same plant = host only
          !has_host & has_nectar ~ "nectar_only",
          TRUE ~ "none"
        )
      )
    
    list(
      n_both = sum(butterfly_support$support_category == "both"),
      n_host_only = sum(butterfly_support$support_category == "host_only"),
      n_nectar_only = sum(butterfly_support$support_category == "nectar_only")
    )
  }
  
  # Run multiple random iterations
  results <- lapply(1:n_iterations, function(iter) {
    # Randomly sample n plants
    random_plants <- sample(all_plants, n_plants)
    
    # Calculate how many butterflies are supported by this random mix
    n_supported <- calculate_fitness(random_plants)
    
    # Calculate interaction breakdown
    breakdown <- calculate_interaction_breakdown(random_plants)
    
    # Return results for this iteration
    data.frame(
      ecoregion = ecoregion,
      iteration = iter,
      n_plants = n_plants,
      butterflies_supported = n_supported,
      n_both_interactions = breakdown$n_both,
      n_host_only = breakdown$n_host_only,
      n_nectar_only = breakdown$n_nectar_only
    )
  })
  
  # Combine all iterations
  combined_results <- do.call(rbind, results)
  
  # Add summary statistics as attributes
  attr(combined_results, "mean_both") <- mean(combined_results$n_both_interactions)
  attr(combined_results, "mean_host_only") <- mean(combined_results$n_host_only)
  attr(combined_results, "mean_nectar_only") <- mean(combined_results$n_nectar_only)
  attr(combined_results, "sd_both") <- sd(combined_results$n_both_interactions)
  attr(combined_results, "sd_host_only") <- sd(combined_results$n_host_only)
  attr(combined_results, "sd_nectar_only") <- sd(combined_results$n_nectar_only)
  
  return(combined_results)
}


scenario3_random_results <- do.call(rbind, mclapply(unique(ecoregions$JEPCODE), function(ecoregion) {
  find_randomMix_butterflies(
    ecoregion, 
    interaction_data_subset = predicted_interactions,
    candidate_plants = plant_candidates_scenario3,
    focal_pollinators = butterflies_declining_SW_list,
    checklist = checklist_cleaned,
    n_plants = n_plants,
    n_iterations = 100
  )
},
mc.cores = 12))

scenario3_random_summary <- scenario3_random_results %>%
  group_by(ecoregion) %>%
  summarise(
    mean_butterflies_supported = mean(butterflies_supported, na.rm = TRUE),
    mean_both_interactions = mean(n_both_interactions, na.rm = TRUE),
    mean_host_only = mean(n_host_only, na.rm = TRUE),
    mean_nectar_only = mean(n_nectar_only, na.rm = TRUE),
    sd_butterflies_supported = sd(butterflies_supported, na.rm = TRUE),
    sd_both_interactions = sd(n_both_interactions, na.rm = TRUE),
    sd_host_only = sd(n_host_only, na.rm = TRUE),
    sd_nectar_only = sd(n_nectar_only, na.rm = TRUE),
    n_iterations = n()
  )

# View the summary
print(scenario3_random_summary)

# Save
write.csv(scenario3_results[,1:7], "Data_Clean/Scenario_Results/Scenario3_opt.csv")
write.csv(scenario3_results[,1:7] %>% group_by(ecoregion) %>% slice_head(n=1) %>% dplyr::select(-c(rank, scientificName)), 
          "Data_Clean/Scenario_Results/Scenario3_pols.csv")

##
print(mean(scenario3_results$n_both_interactions))
print(mean(scenario3_random_results$n_both_interactions))
