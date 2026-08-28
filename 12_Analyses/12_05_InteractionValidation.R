#### Load packages ####
library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggdist)
library(units)
library(data.table)
library(purrr)

# Note: this script loads its own data rather than sourcing 14_00_LoadData.R

#### Read in data ####
# California ecoregions
ecoregions <- read_sf("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid()

## Checklists
checklist_cleaned <- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# Plant checklist
plant_checklist <- checklist_cleaned[checklist_cleaned$taxon == "plants", ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Pollinator checklists
pollinator_checklists <- checklist_cleaned[checklist_cleaned$taxon %in% c("bees", "butterflies", "moths", "hoverflies"), ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Occurrences
occurrences <- fread("Data_Clean/Occurrences/occurrences_clean.csv") %>%
  filter(genus_species %in% checklist_cleaned$genus_species)

# Phenology
phenology <- fread("Data_Clean/Phenology/phenology_best_available.csv") %>% distinct()

## Interactions
interactions <- fread("Data_Clean/Interactions/ints_final_clean_slim.csv")

## Spatially explicit interactions in CA
interactions$decimalLatitude  <- as.numeric(interactions$decimalLatitude)
interactions$decimalLongitude <- as.numeric(interactions$decimalLongitude)

interactions_sf <- st_as_sf(
  interactions[!is.na(interactions$decimalLatitude) &
                 !is.na(interactions$decimalLongitude), ],
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)

interactions_sf_witheco <- st_join(interactions_sf, ecoregions) %>%
  st_drop_geometry()

interactions_regional <- interactions_sf_witheco %>%
  filter(!is.na(JEPCODE))

interactions_witheco <- left_join(interactions, interactions_regional) %>%
  mutate(sourceTaxonSpeciesName = word(sourceTaxonName_harm, 2),
         targetTaxonSpeciesName = word(targetTaxonName_harm, 2))

interactions_regional_visitsFlowersOf <- filter(interactions_regional, interactionTypeName == "visitsFlowersOf")
interactions_regional_hasHost         <- filter(interactions_regional, interactionTypeName == "hasHost")

spatialInteractions_unique_visitsFlowersOf <- interactions_regional_visitsFlowersOf %>%
  filter(str_count(targetTaxonName_harm, "\\S+") >= 2,
         str_count(sourceTaxonName_harm, "\\S+") >= 2,
         !is.na(JEPCODE)) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE, interactionTypeName) %>%
  mutate(
    sourceTaxonName_harm    = word(sourceTaxonName_harm, start = 1, end = 2),
    targetTaxonName_harm    = word(targetTaxonName_harm, start = 1, end = 2),
    sourceTaxonGenusName    = word(sourceTaxonName_harm, 1),
    targetTaxonGenusName    = word(targetTaxonName_harm, 1),
    sourceTaxonSpeciesName  = word(sourceTaxonName_harm, 2),
    targetTaxonSpeciesName  = word(targetTaxonName_harm, 2)
  ) %>%
  filter(sourceTaxonName_harm %in% checklist_cleaned$genus_species,
         targetTaxonName_harm %in% checklist_cleaned$genus_species) %>%
  unique()

#### Nectar Interaction Predictions ####
predict_visitsFlowersOf <- read.csv("Data_Clean/Interactions/predictions_supporting_visitsFlowersOf.csv")
predict_hasHost         <- read.csv("Data_Clean/Interactions/predictions_supporting_hasHost.csv")

#### Normalize Spatial Overlap ####
sdm_example <- rast(list.files("Data_Clean/SDMs/sdm_by_taxon/plants/continuous", full.names = TRUE)[1])
res_value <- res(sdm_example)[1]
crs_value <- crs(sdm_example)

ecoregions_union_full <- ecoregions %>%
  st_union() %>%
  st_make_valid()

ecoregions_vect_full <- vect(as(ecoregions_union_full, "Spatial"))

ones_rast <- rast(ext(ecoregions_vect_full), resolution = res_value, crs = crs_value)
values(ones_rast) <- 1
ones_rast <- mask(ones_rast, ecoregions_vect_full)

ecoregions_union <- ecoregions %>%
  group_by(JEPCODE) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

unique_jepcodes <- ecoregions_union$JEPCODE

pixel_counts <- data.frame(JEPCODE = unique_jepcodes, n_pixels = NA_integer_)

for (i in seq_along(unique_jepcodes)) {
  eco_single <- ecoregions_union %>% filter(JEPCODE == unique_jepcodes[i])
  eco_vect   <- vect(as(eco_single, "Spatial"))
  masked     <- terra::mask(terra::crop(ones_rast, eco_vect), eco_vect)
  pixel_counts$n_pixels[i] <- sum(!is.na(values(masked)))
}

predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

predict_hasHost <- predict_hasHost %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

#### Normalize temporal overlap ####
predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  mutate(temporal_overlap_normalized = temporal_overlap / plant_duration)

#### Predict flower visitation interactions ####
predictions_visitsFlowersOf <- predict_visitsFlowersOf %>%
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    temporal_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized * temporal_overlap_normalized
  ))

predictions_visitsFlowersOf_NA <- predictions_visitsFlowersOf[is.na(predictions_visitsFlowersOf$interactionScore), ]
predictions_visitsFlowersOf    <- predictions_visitsFlowersOf[!is.na(predictions_visitsFlowersOf$interactionScore), ]

predictions_visitsFlowersOf <- predictions_visitsFlowersOf %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1))

predictions_hasHost <- predict_hasHost %>%
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized
  ))

predictions_hasHost_NA <- predictions_hasHost[is.na(predictions_hasHost$interactionScore), ]
predictions_hasHost    <- predictions_hasHost[!is.na(predictions_hasHost$interactionScore), ]

predictions_hasHost <- predictions_hasHost %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1))

###########################################################################
#### Validation 1 - Leave out validation  ####
###########################################################################all_pollinators <- pollinator_checklists$genus_species
predictions_visitsFlowersOf_Validation1 <- predictions_visitsFlowersOf %>%
  filter(method_type == "Level1_Simple") %>%
  unique()

all_plants      <- plant_checklist$genus_species
all_JEPCODE     <- unique(ecoregions$JEPCODE)

predictions_visitsFlowersOf_Validation1$targetTaxonGenusName <- sub(" .*", "", predictions_visitsFlowersOf_Validation1$targetTaxonName_harm)

# Validation function with multiple null runs
run_validation_flowerVisitation <- function(interactions_witheco,
                                            predictions_visitsFlowersOf,
                                            spatialInteractions_unique_visitsFlowersOf,
                                            plant_checklist,
                                            pollinator_checklists,
                                            ecoregions,
                                            validation_prop = 0.2,
                                            n_null_runs = 100,
                                            seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Define potential interaction pairs (pol species - plant genus)
  potential_pairs <- predictions_visitsFlowersOf %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonGenusName) %>%
    distinct()
  
  # 2. Randomly select interactions to remove (validation set)
  spatialInteractions_toRemove <- spatialInteractions_unique_visitsFlowersOf %>%
    slice_sample(prop = validation_prop)
  
  # 3. Remove validation set from training data
  interactions_filteredValidation <- interactions_witheco %>%
    filter(interactionTypeName == "visitsFlowersOf") %>%
    anti_join(spatialInteractions_toRemove, 
              by = c("sourceTaxonGenusName", "targetTaxonGenusName", 
                     "sourceTaxonSpeciesName", "targetTaxonSpeciesName", 
                     "JEPCODE"))
  
  # 4. Identify which potential pairs are still allowed after removal
  ints_pol_plantgenera_filteredValidation <- interactions_filteredValidation %>%
    dplyr::select(targetTaxonGenusName, sourceTaxonName_harm) %>%
    filter(str_count(sourceTaxonName_harm, "\\S+") >= 2) %>%
    mutate(sourceTaxonName_harm = word(sourceTaxonName_harm, 1, 2))
  
  ints_pol_plantgenera_checklist_filteredValidation <- ints_pol_plantgenera_filteredValidation %>%
    filter(sourceTaxonName_harm %in% pollinator_checklists$genus_species &
             targetTaxonGenusName %in% plant_checklist$genus) %>%
    unique()
  
  potential_pairs_validation <- semi_join(potential_pairs, 
                                          ints_pol_plantgenera_checklist_filteredValidation,
                                          by = c("sourceTaxonName_harm", "targetTaxonGenusName"))
  
  # 5. Mark which predictions are still "allowed" (validationPotentialLink)
  predictions_marked <- predictions_visitsFlowersOf %>%
    mutate(validationPotentialLink = ifelse(
      paste(sourceTaxonName_harm, targetTaxonGenusName) %in% 
        paste(potential_pairs_validation$sourceTaxonName_harm, 
              potential_pairs_validation$targetTaxonGenusName),
      1, 0
    ))
  
  # 6. Add observed interactions from filtered training data
  int_observed <- interactions_filteredValidation %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    mutate(observedInteraction = 1) %>%
    distinct()
  
  predictions_with_obs <- predictions_marked %>%
    left_join(int_observed, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    mutate(observedInteraction = replace_na(observedInteraction, 0))
  
  # 7. Calculate threshold from confirmed interactions
  confirmed_probs <- predictions_with_obs %>%
    filter(observedInteraction == 1, validationPotentialLink == 1) %>%
    pull(interactionScore)
  
  lower_threshold <- quantile(confirmed_probs, 0.25, na.rm = TRUE)
  
  # 8. Filter for predicted interactions (above threshold OR observed, AND allowed)
  filtered_interactions <- predictions_with_obs %>%
    filter((interactionScore >= lower_threshold | observedInteraction == 1),
           validationPotentialLink == 1)
  
  # 9. Calculate actual recovery rate
  matches_fromRemoved <- spatialInteractions_toRemove %>%
    semi_join(filtered_interactions, 
              by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
  
  actual_recovery_prop <- nrow(matches_fromRemoved) / nrow(spatialInteractions_toRemove)
  
  # 10. Count reasons for non-recovery
  n_removed_forbidden <- spatialInteractions_toRemove %>%
    anti_join(potential_pairs_validation,
              by = c("sourceTaxonName_harm")) %>%
    nrow()
  
  n_removed_below_threshold <- spatialInteractions_toRemove %>%
    semi_join(potential_pairs_validation, 
              by = c("sourceTaxonName_harm")) %>%
    anti_join(filtered_interactions,
              by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    nrow()
  
  # 10b. Per-interaction reason (for completeness-binned pie charts downstream)
  removed_forbidden_ids <- spatialInteractions_toRemove %>%
    anti_join(potential_pairs_validation, by = c("sourceTaxonName_harm"))
  
  removed_interactions_with_reason <- spatialInteractions_toRemove %>%
    mutate(
      reason = case_when(
        paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
          paste(matches_fromRemoved$sourceTaxonName_harm, matches_fromRemoved$targetTaxonName_harm, matches_fromRemoved$JEPCODE) ~ "Recovered",
        paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
          paste(removed_forbidden_ids$sourceTaxonName_harm, removed_forbidden_ids$targetTaxonName_harm, removed_forbidden_ids$JEPCODE) ~ "Forbidden interaction",
        TRUE ~ "Below threshold"
      )
    )
  
  # 11. Null models - run multiple times
  n_draw <- filtered_interactions %>%
    anti_join(int_observed, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    nrow()
  
  all_pollinators <- pollinator_checklists$genus_species
  all_plants <- plant_checklist$genus_species
  all_JEPCODE <- unique(ecoregions$JEPCODE)
  
  # Initialize vectors
  null_1_props <- numeric(n_null_runs)
  null_2_props <- numeric(n_null_runs)
  
  # For constrained validation (Null 2), what possible pairs?
  potential_for_null2 <- predictions_with_obs %>%
    filter(validationPotentialLink == 1,
           observedInteraction == 0) %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct()
  
  for (null_run in 1:n_null_runs) {
    ## Null model 1 - completely random interactions
    null_1 <- expand.grid(
      sourceTaxonName_harm = sample(all_pollinators, min(length(all_pollinators), 500)),
      targetTaxonName_harm = sample(all_plants, min(length(all_plants), 500)),
      JEPCODE = all_JEPCODE,
      stringsAsFactors = FALSE
    ) %>%
      slice_sample(n = n_draw, replace = FALSE)
    
    null_1_overlap <- spatialInteractions_toRemove %>%
      semi_join(null_1, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
    
    null_1_props[null_run] <- nrow(null_1_overlap) / nrow(spatialInteractions_toRemove)
    
    ## Null model 2 - random from allowed potential interactions
    if (nrow(potential_for_null2) >= n_draw) {
      null_2 <- potential_for_null2 %>%
        slice_sample(n = n_draw, replace = FALSE)
      
      null_2_overlap <- spatialInteractions_toRemove %>%
        semi_join(null_2, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
      
      null_2_props[null_run] <- nrow(null_2_overlap) / nrow(spatialInteractions_toRemove)
    } else {
      null_2_props[null_run] <- NA
    }
  }
  
  # 12. Calculate p-values
  null_1_beats_actual <- mean(null_1_props >= actual_recovery_prop)
  null_2_beats_actual <- mean(null_2_props >= actual_recovery_prop, na.rm = TRUE)
  
  # 13. Return results
  list(
    summary = data.frame(
      validation_prop = validation_prop,
      n_removed = nrow(spatialInteractions_toRemove),
      n_predicted = nrow(filtered_interactions),
      n_observed_in_training = nrow(int_observed),
      threshold = lower_threshold,
      actual_overlap_n = nrow(matches_fromRemoved),
      actual_recovery_prop = actual_recovery_prop,
      n_removed_forbidden = n_removed_forbidden,
      n_removed_below_threshold = n_removed_below_threshold,
      removed_interactions_with_reason = removed_interactions_with_reason,
      null1_mean_prop = mean(null_1_props),
      null1_sd_prop = sd(null_1_props),
      null1_beats_actual_prop = null_1_beats_actual,
      null2_mean_prop = mean(null_2_props, na.rm = TRUE),
      null2_sd_prop = sd(null_2_props, na.rm = TRUE),
      null2_beats_actual_prop = null_2_beats_actual
    ),
    null_1_props = null_1_props,
    null_2_props = null_2_props,
    removed_interactions = spatialInteractions_toRemove,
    recovered_interactions = matches_fromRemoved
  )
}

# wrapper function
run_multiple_validation_flowerVisitation <- function(interactions_witheco,
                                                     predictions_visitsFlowersOf,
                                                     spatialInteractions_unique_visitsFlowersOf,
                                                     plant_checklist,
                                                     pollinator_checklists,
                                                     ecoregions,
                                                     n_runs = 100,
                                                     n_null_runs = 100,
                                                     validation_prop = 0.2) {
  
  # Initialize results list
  results_list <- list()
  all_null_1_props <- list()
  all_null_2_props <- list()
  
  # Run validation in loop with progress
  for (i in 1:n_runs) {
    cat("Completed", i, "of", n_runs, "runs\n")
    
    result <- run_validation_flowerVisitation(
      interactions_witheco = interactions_witheco,
      predictions_visitsFlowersOf = predictions_visitsFlowersOf,
      spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
      plant_checklist = plant_checklist,
      pollinator_checklists = pollinator_checklists,
      ecoregions = ecoregions,
      validation_prop = validation_prop,
      n_null_runs = n_null_runs,
      seed = i
    )
    
    results_list[[i]] <- result$summary
    all_null_1_props[[i]] <- result$null_1_props
    all_null_2_props[[i]] <- result$null_2_props
  }
  
  # Combine results
  list(
    summary_df = do.call(rbind, results_list),
    null_1_distributions = all_null_1_props,
    null_2_distributions = all_null_2_props
  )
}

# Run the validation
validation_results_flowerVisitation <- run_multiple_validation_flowerVisitation(
  interactions_witheco = interactions_witheco,
  predictions_visitsFlowersOf = predictions_visitsFlowersOf_Validation1,
  spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
  plant_checklist = plant_checklist,
  pollinator_checklists = pollinator_checklists,
  ecoregions = ecoregions,
  n_runs = 25,
  n_null_runs = 100,
  validation_prop = 0.2
)

# View summary statistics across all runs
# summary(validation_results_flowerVisitation$summary_df)

saveRDS(validation_results_flowerVisitation, "Data_Clean/Validation/interaction_validation_results.rds")

validation_results_flowerVisitation <- readRDS("Data_Clean/Validation/interaction_validation_results.rds")

#### Plot: Interaction plausibility score histograms ####
predictions_visitsFlowersOf_Validation1 <- predictions_visitsFlowersOf_Validation1 %>%
  mutate(observed = paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
           paste(spatialInteractions_unique_visitsFlowersOf$sourceTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$targetTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$JEPCODE))



# Split data
obs_true_A  <- predictions_visitsFlowersOf_Validation1 %>% filter(observed == TRUE)
obs_false_A <- predictions_visitsFlowersOf_Validation1 %>% filter(observed == FALSE)

scale_factor_A <- max(table(cut(obs_false_A$interactionScore, breaks = 30))) /
  max(table(cut(obs_true_A$interactionScore,  breaks = 30)))

hist_scores <- ggplot() +
  geom_histogram(data = obs_false_A, aes(x = interactionScore), bins = 30,
                 fill = "#FFCC00", color = "#331E0F") +
  geom_histogram(data = obs_true_A,
                 aes(x = interactionScore, y = after_stat(count) * scale_factor_A),
                 bins = 30, fill = "#BAFF76", color = "#3B7004", alpha = 0.5) +
  geom_vline(xintercept = quantile(obs_true_A$interactionScore, 0.25),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(
    name = "Density of unobserved\ninteractions (yellow)",
    sec.axis = sec_axis(~ . / scale_factor_A,
                        name = "Density of observed spatially-\nexplicit interactions (green)")
  ) +
  scale_x_continuous(name = "Interaction Score") +
  theme_cowplot()

## Has Host Interaction Scores Plot
#### Plot: Interaction plausibility score histograms ####
predictions_hasHost <- predictions_hasHost %>%
  mutate(observed = paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
           paste(spatialInteractions_unique_visitsFlowersOf$sourceTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$targetTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$JEPCODE))

# Split data
obs_true_B  <- predictions_hasHost %>% filter(observed == TRUE)
obs_false_B <- predictions_hasHost %>% filter(observed == FALSE)

scale_factor_B <- max(table(cut(obs_false_B$interactionScore, breaks = 30))) /
  max(table(cut(obs_true_B$interactionScore,  breaks = 30)))

hist_scores_hasHost <- ggplot() +
  geom_histogram(data = obs_false_B,
                 aes(x = interactionScore),
                 bins = 30,
                 fill = "#FFCC00", color = "#331E0F") +
  geom_histogram(data = obs_true_B,
                 aes(x = interactionScore, y = after_stat(count) * scale_factor_B),
                 bins = 30,
                 fill = "#BAFF76", color = "#3B7004", alpha = 0.5) +
  geom_vline(xintercept = quantile(obs_true_B$interactionScore, 0.25),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(
    name = "Density of unobserved\ninteractions (yellow)",
    sec.axis = sec_axis(~ . / scale_factor_B,
                        name = "Density of observed spatially-\nexplicit interactions (green)")
  ) +
  scale_x_continuous(name = "Interaction Score") +
  theme_cowplot()

hist_scores <- hist_scores + theme(axis.title = element_text(size = 12),
                                   axis.text = element_text(size = 10))
hist_scores_hasHost <- hist_scores_hasHost + theme(axis.title = element_text(size = 12),
                                                   axis.text = element_text(size = 10))

combined_hist <- (hist_scores | hist_scores_hasHost) +
  plot_annotation(tag_levels = 'A')


ggsave("Figures/Supplementary/prediction_interactionScore_histogram.png", 
       combined_hist, 
       width = 200, 
       height = 100, 
       units = "mm", 
       dpi = 2000)

# corresponding percentile for unobserved
mean(obs_false_A$interactionScore <= quantile(obs_false_A$interactionScore, 0.25))


#### Plot: Recall rate comparison null models vs NECTAR ####
# Prepare data - expand null distributions
plot_data <- bind_rows(
  # Null 1 distributions
  tibble(
    fold = rep(1:25, each = 100),
    model = "Null 1: Completely Random",
    recovery_prop = unlist(validation_results_flowerVisitation$null_1_distributions)
  ),
  # Null 2 distributions
  tibble(
    fold = rep(1:25, each = 100),
    model = "Null 2: Phylogenetically Constrained",
    recovery_prop = unlist(validation_results_flowerVisitation$null_2_distributions)
  )
) %>%
  mutate(fold = factor(fold, levels = 25:1))  # Reverse for top-to-bottom

# Actual recovery as points
actual_data <- validation_results_flowerVisitation$summary_df %>%
  mutate(
    fold = factor(row_number(), levels = 25:1),
    model = "Actual"
  ) %>%
  dplyr::select(fold, model, recovery_prop = actual_recovery_prop)

validation_recall <- ggplot() +
  # show intervals of null distributions
  stat_interval(data = plot_data,
                aes(y = fold, x = recovery_prop, color = model),
                .width = c(.50, .95),
                position = position_dodge(width = 0.4),
                linewidth = 2) +
  # add points for actual recovery (now with shape aesthetic mapped)
  geom_point(data = actual_data,
             aes(y = fold, x = recovery_prop, shape = model),
             size = 3, color = "#FFCC00") +
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.3) +
  scale_color_manual(
    values = c("Null 1: Completely Random" = "gray70", 
               "Null 2: Phylogenetically Constrained" = "gray40"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Actual" = 18),
    name = NULL,
    labels = c("NECTAR")
  ) +
  scale_y_discrete(labels = 1:25) +
  labs(
    x = "Recall Rate (Proportion of withheld interactions recovered)",
    y = "Validation Fold"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

ggsave("Figures/Supplementary/validation_recall.png", 
       validation_recall, 
       width = 200, 
       height = 200, 
       units = "mm", 
       dpi = 2000)

#### Plot: Summary of why non-recovery of interaction ####
reason_summary <- validation_results_flowerVisitation$summary_df %>%
  summarise(
    `Recovered` = sum(actual_overlap_n),
    `Forbidden interaction` = sum(n_removed_forbidden),
    `Below threshold` = sum(n_removed_below_threshold)
  ) %>%
  pivot_longer(everything(), names_to = "reason", values_to = "n_interactions") %>%
  mutate(
    pct = round(n_interactions / sum(n_interactions) * 100, 1),
    reason = factor(reason, levels = c("Recovered", "Below threshold", "Forbidden interaction"))
  )

# Create pie chart
validation_outcome <- ggplot(reason_summary, aes(x = "", y = n_interactions, fill = reason)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y") +
  scale_fill_manual(
    values = c("Recovered" = "#FFCC00", 
               "Below threshold" = "gray60", 
               "Forbidden interaction" = "gray80")
  ) +
  theme_void() +
  theme(legend.position = "right") +
  labs(
    fill = "Outcome of\nWithheld Interaction",
  ) +
  # add percentage labels
  geom_text(aes(label = paste0(pct, "%")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold")

ggsave("Figures/Supplementary/validation_outcomeBreakdown.png", 
       validation_outcome, 
       width = 200, 
       height = 200, 
       units = "mm", 
       dpi = 2000)



###########################################################################
#### Validation 2 - Recall rate vs. data completeness      ####
###########################################################################
## Rerun leave-out validation, keeping per-interaction recovery outcomes ####
n_runs <- 25

per_run_flagged <- map(1:n_runs, function(i) {
  result <- run_validation_flowerVisitation(
    interactions_witheco = interactions_witheco,
    predictions_visitsFlowersOf = predictions_visitsFlowersOf_Validation1,
    spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
    plant_checklist = plant_checklist,
    pollinator_checklists = pollinator_checklists,
    ecoregions = ecoregions,
    validation_prop = 0.2,
    n_null_runs = 1, # only run null models once, not important here
    seed = i
  )
  
  recovered_ints <- result$recovered_interactions %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct() %>%
    mutate(recovered = TRUE)
  
  result$removed_interactions %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    left_join(recovered_ints, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    mutate(recovered = replace_na(recovered, FALSE),
           run = i)
})

## Data completeness from whole interaction dataset
compute_completeness <- function(taxon_col) {
  interactions_witheco %>%
    filter(interactionTypeName == "visitsFlowersOf",
           str_count(.data[[taxon_col]], "\\S+") >= 2) %>%
    mutate(!!taxon_col := word(.data[[taxon_col]], 1, 2)) %>%
    count(.data[[taxon_col]], name = "n_records_total")
}

completeness_by_pollinator <- compute_completeness("sourceTaxonName_harm")
completeness_by_plant      <- compute_completeness("targetTaxonName_harm")

## Bin the data completeness values
# Then compare within a single validation iteration, summarise across iterations
# This is so we can avoid pseudo-replication issues of witholding the same interactions over and over.
add_completeness_bins <- function(completeness_df) {
  completeness_df %>%
    mutate(
      completeness_bin = cut(
        log10(n_records_total),
        breaks = quantile(log10(n_records_total), probs = seq(0, 1, 0.25), na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Low", "Mid-low", "Mid-high", "High")
      )
    )
}

completeness_by_pollinator <- add_completeness_bins(completeness_by_pollinator)
completeness_by_plant      <- add_completeness_bins(completeness_by_plant)

## Recall rate
compute_recall_by_bin <- function(per_run_flagged, completeness_df, join_col) {
  map_dfr(per_run_flagged, function(run_df) {
    run_df %>%
      left_join(completeness_df, by = join_col) %>%
      filter(!is.na(completeness_bin)) %>%
      group_by(run, completeness_bin) %>%
      summarize(
        n_removed   = n(),
        n_recovered = sum(recovered),
        recall_rate = n_recovered / n_removed,
        .groups = "drop"
      )
  })
}

recall_by_run_bin_pollinator <- compute_recall_by_bin(per_run_flagged, completeness_by_pollinator, "sourceTaxonName_harm")
recall_by_run_bin_plant      <- compute_recall_by_bin(per_run_flagged, completeness_by_plant, "targetTaxonName_harm")

## plot data completeness vs. recall rate
p_pollinator <- ggplot(recall_by_run_bin_pollinator, aes(x = completeness_bin, y = recall_rate)) +
  geom_boxplot(fill = "#FFCC00", color = "#331E0F", outlier.size = 1) +
  geom_jitter(width = 0.1, height = 0, color = "#331E0F", alpha = 0.4, size = 1) +
  labs(
    x = "Pollinator data completeness",
    y = "Recall rate"  ) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 14)) +
  theme_cowplot()

p_plant <- ggplot(recall_by_run_bin_plant, aes(x = completeness_bin, y = recall_rate)) +
  geom_boxplot(fill = "#FFCC00", color = "#331E0F", outlier.size = 1) +
  geom_jitter(width = 0.1, height = 0, color = "#331E0F", alpha = 0.4, size = 1) +
  labs(
    x = "Plant data completeness",
    y = NULL) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 14)) +
  theme_cowplot()

## Combine
combined_completeness <- (p_pollinator | p_plant) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_completeness)

ggsave("Figures/Supplementary/recall_vs_completeness.png",
       combined_completeness,
       width = 200,
       height = 130,
       units = "mm",
       dpi = 2000)

###########################################################################
#### Validation 3 - External Sources  ####
###########################################################################
source("Code/99_Supporting/name_harmonization.R")
# Define gnverifier path
gnverifier_path <- "/home/tb625/bin"

### Gather raw data ###
## Dataset 1: Hung 2019 San Diego plant-pollinator surveys in coastal sage scrub remnant habitat
SCo_Hung_2019_plants <- read_excel("Data_Raw/External_Validation/SCo_Hung_2019.xlsx", sheet = "Plant codes")

SCo_Hung_2019_1 <- read_excel("Data_Raw/External_Validation/SCo_Hung_2019.xlsx", sheet = "Interactions 2015") %>%
  left_join(SCo_Hung_2019_plants)
SCo_Hung_2019_2 <- read_excel("Data_Raw/External_Validation/SCo_Hung_2019.xlsx", sheet = "Interactions 2016") %>%
  left_join(SCo_Hung_2019_plants)

SCo_Hung_2019 <- rbind(SCo_Hung_2019_1, SCo_Hung_2019_2) %>%
  dplyr::select('Pollinator', 'Plant species name') %>%
  rename(targetTaxon = 'Plant species name',
         sourceTaxon = 'Pollinator') %>%
  mutate(sourceTaxon = gsub(".", " ", sourceTaxon, fixed = TRUE)) %>%
  mutate(JEPCODE = "SCo",
         dataset = "Hung_2019")

## Dataset 2: Hazlehurst 2023 San Clemente Island Rare Plant surveys of plant-pollinator interactions
SChI_Hazlehurst_2023 <- read_excel("Data_Raw/External_Validation/SChI_Hazlehurst_2023.xlsx") %>%
  mutate(JEPCODE = "sChI",
         dataset = "Hazlehurst_2023")

## Dataset 3: BBNA
BBNA <- read.csv("Data_Raw/External_Validation/BBNA.share.01.08.2025.csv") 

BBNA_sf <- BBNA %>%
  filter(!is.na(longitude),
         !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(ecoregions))

BBNA_with_eco <- st_join(BBNA_sf, ecoregions) %>%
  st_drop_geometry() %>%
  filter(!is.na(JEPCODE))

BBNA_cleaned <- BBNA_with_eco %>%
  filter(plant.host.genus != "",
         plant.host.species != "") %>%
  mutate(targetTaxon = paste(plant.host.genus, plant.host.species),
         sourceTaxon = paste("Bombus", species),
         dataset = "BBNA") %>%
  dplyr::select(targetTaxon, sourceTaxon, JEPCODE, dataset) %>%
  unique()

## Dataset 4: BUGS (Biodiversity in Urban Gardens - CNPS inat project)
BUGS <- read.csv("Data_Raw/External_Validation/BUGS.csv") 

# Filter to only research grade observations, genus is in checklist, and where description is not empty (i.e. potentially has some plant interaction information)
BUGS_filtered <- BUGS %>%
  mutate(genus = word(scientific_name, 1, 1)) %>%
  filter(genus %in% unique(checklist_cleaned$genus),
         quality_grade == "research",
         description != "")

write.csv(BUGS_filtered, "Data_Raw/External_Validation/BUGS_filtered_raw.csv")

## Clean up plant names manually then read in file
BUGS_clean <- read.csv("Data_Raw/External_Validation/BUGS_filtered_clean.csv")

BUGS_clean_sf <- BUGS_clean %>%
  filter(!is.na(longitude),
         !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(ecoregions))

BUGS_with_eco <- st_join(BUGS_clean_sf, ecoregions) %>%
  st_drop_geometry() %>%
  filter(!is.na(JEPCODE))

BUGS_clean_witheco <- BUGS_with_eco %>%
  mutate(targetTaxon = targetTaxonName,
         sourceTaxon = sourceTaxonName,
         dataset = "BUG") %>%
  filter(targetTaxon != "") %>%
  dplyr::select(targetTaxon, sourceTaxon, JEPCODE, dataset) %>%
  unique()



## Combine
ext_validation_data_raw = rbind(SCo_Hung_2019, SChI_Hazlehurst_2023, BBNA_cleaned, BUGS_clean_witheco)

# Be sure to first install gnparser from https://github.com/gnames/gnparser
parse_names_plant <-
  bdc_clean_names(sci_names = ext_validation_data_raw$targetTaxon, save_outputs = FALSE)
parse_names_pol <-
  bdc_clean_names(sci_names = ext_validation_data_raw$sourceTaxon, save_outputs = FALSE)

# Examine errors
wrong_names_plant <- filter(parse_names_plant, is.na(names_clean) | quality == 0)
wrong_names_pol <- filter(parse_names_pol, is.na(names_clean) | quality == 0)

### Name harmonization
names_plants <- parse_names_plant %>% pull(names_clean) %>% unique()
names_animals <- parse_names_pol %>% pull(names_clean) %>% unique()

# 2. Harmonize names
names_harm_plants <- harmonize_names(
  names = names_plants,
  higher_tax = "Plantae",
  names_file = "Temp/names_validation_plants.tsv",
  names_harm_file = "Temp/names_harmonized_validation_plants.tsv",
  gnv_path = gnverifier_path
)

names_harm_animals <- harmonize_names(
  names = names_animals,
  higher_tax = "Animalia",
  names_file = "Temp/names_validation_animals.tsv",
  names_harm_file = "Temp/names_harmonized_validation_animals.tsv",
  gnv_path = gnverifier_path
)

# 3. Apply Jepson taxonomy to plants
names_harm_plants <- apply_jepson_taxonomy(names_harm_plants)

# 4. Join back with occurrence data
ext_validation_data_clean <- left_join(ext_validation_data_raw,
                                       names_harm_plants,
                                       by = c("targetTaxon" = "ScientificName")) %>%
  rename(targetTaxonName_harm = names_harm)
ext_validation_data_clean <- left_join(ext_validation_data_clean,
                                       names_harm_animals,
                                       by = c("sourceTaxon" = "ScientificName")) %>%
  rename(sourceTaxonName_harm = names_harm)

ext_validation_data_clean <- ext_validation_data_clean %>% 
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, JEPCODE, dataset) %>%
  # downscale taxonomy to genus species level to match our data
  mutate(targetTaxonName_harm = word(targetTaxonName_harm, 1, 2),
         sourceTaxonName_harm = word(sourceTaxonName_harm, 1, 2),)

# Filter for only plants and pollinators in our dataset
checklist_cleaned <- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

ext_validation_data_clean <- ext_validation_data_clean %>%
  filter(sourceTaxonName_harm %in% checklist_cleaned$genus_species,
         targetTaxonName_harm %in% checklist_cleaned$genus_species)

interactions_visitsFlowersOf_only <- interactions %>%
  filter(interactionTypeName == "visitsFlowersOf")

# Flag the novel interactions (not in prior observed interaction data)
ext_validation_data_clean <- ext_validation_data_clean %>%
  mutate(novel_eco = !(paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
                         paste(spatialInteractions_unique_visitsFlowersOf$sourceTaxonName_harm,
                               spatialInteractions_unique_visitsFlowersOf$targetTaxonName_harm,
                               spatialInteractions_unique_visitsFlowersOf$JEPCODE)) &
           (paste(sourceTaxonName_harm, targetTaxonName_harm) %in%
              paste(interactions_visitsFlowersOf_only$sourceTaxonName_harm,
                    interactions_visitsFlowersOf_only$targetTaxonName_harm)),
         novel = !(paste(sourceTaxonName_harm, targetTaxonName_harm) %in%
                     paste(interactions_visitsFlowersOf_only$sourceTaxonName_harm,
                           interactions_visitsFlowersOf_only$targetTaxonName_harm))) %>%
  distinct()

# Quick descriptive breakdown by dataset, before running anything
ext_validation_data_clean %>%
  count(dataset, novel) %>%
  print()
ext_validation_data_clean %>%
  count(dataset, novel_eco) %>%
  print()

## Get % of unique interactions that are novel per dataset
# species pairs
ext_validation_data_clean %>%
  dplyr::select(-c(JEPCODE, novel_eco)) %>%
  distinct() %>%
  group_by(dataset) %>%
  summarise(perc_novel_spPair = sum(novel)/n(),
            n = n())

# spatially explicit observations
ext_validation_data_clean %>%
  group_by(dataset) %>%
  filter(novel==0) %>%
  summarise(perc_novel_ecoregion = sum(novel_eco)/n())

ext_validation_novel_withDataset <- ext_validation_data_clean %>%
  filter(novel) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE, dataset) %>%
  distinct()
ext_validation_novel_eco_withDataset <- ext_validation_data_clean %>%
  filter(novel_eco) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE, dataset) %>%
  distinct()

ext_validation_novel <- ext_validation_data_clean %>%
  filter(novel) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  distinct()
ext_validation_novel_eco <- ext_validation_data_clean %>%
  filter(novel_eco) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  distinct()

# Run validation
run_external_validation_flowerVisitation <- function(predictions_visitsFlowersOf,
                                                     ext_validation_novel,
                                                     spatialInteractions_unique_visitsFlowersOf,
                                                     plant_checklist,
                                                     pollinator_checklists,
                                                     ecoregions,
                                                     n_null_runs = 100,
                                                     seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  predictions_visitsFlowersOf <- predictions_visitsFlowersOf %>%
    mutate(targetTaxonGenusName = sub(" .*", "", targetTaxonName_harm))
  
  potential_pairs <- predictions_visitsFlowersOf %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonGenusName) %>%
    distinct()
  
  ext_validation_novel <- ext_validation_novel %>%
    mutate(targetTaxonGenusName = sub(" .*", "", targetTaxonName_harm),
           validationPotentialLink = paste(sourceTaxonName_harm, targetTaxonGenusName) %in%
             paste(potential_pairs$sourceTaxonName_harm, potential_pairs$targetTaxonGenusName))
  
  n_novel_total <- nrow(ext_validation_novel)
  
  observed_pairs <- spatialInteractions_unique_visitsFlowersOf %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct() %>%
    mutate(observedInteraction = 1)
  
  predictions_with_obs <- predictions_visitsFlowersOf %>%
    left_join(observed_pairs, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    mutate(observedInteraction = replace_na(observedInteraction, 0))
  
  confirmed_probs <- predictions_with_obs %>%
    filter(observedInteraction == 1) %>%
    pull(interactionScore)
  
  lower_threshold <- quantile(confirmed_probs, 0.25, na.rm = TRUE)
  
  filtered_interactions <- predictions_with_obs %>%
    filter(interactionScore >= lower_threshold | observedInteraction == 1)
  
  matches_novel <- ext_validation_novel %>%
    semi_join(filtered_interactions, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
  
  actual_recovery_prop <- nrow(matches_novel) / n_novel_total
  
  n_novel_forbidden <- ext_validation_novel %>%
    filter(!validationPotentialLink) %>%
    nrow()
  
  n_novel_below_threshold <- ext_validation_novel %>%
    filter(validationPotentialLink) %>%
    anti_join(filtered_interactions, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    nrow()
  
  n_draw <- filtered_interactions %>%
    filter(observedInteraction == 0) %>%
    nrow()
  
  all_pollinators <- pollinator_checklists$genus_species
  all_plants <- plant_checklist$genus_species
  all_JEPCODE <- unique(ecoregions$JEPCODE)
  
  null_1_props <- numeric(n_null_runs)
  null_2_props <- numeric(n_null_runs)
  
  potential_for_null2 <- predictions_with_obs %>%
    filter(observedInteraction == 0) %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct()
  
  for (null_run in 1:n_null_runs) {
    null_1 <- expand.grid(
      sourceTaxonName_harm = sample(all_pollinators, min(length(all_pollinators), 500)),
      targetTaxonName_harm = sample(all_plants, min(length(all_plants), 500)),
      JEPCODE = all_JEPCODE,
      stringsAsFactors = FALSE
    ) %>%
      slice_sample(n = n_draw, replace = FALSE)
    
    null_1_overlap <- ext_validation_novel %>%
      semi_join(null_1, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
    
    null_1_props[null_run] <- nrow(null_1_overlap) / n_novel_total
    
    if (nrow(potential_for_null2) >= n_draw) {
      null_2 <- potential_for_null2 %>%
        slice_sample(n = n_draw, replace = FALSE)
      
      null_2_overlap <- ext_validation_novel %>%
        semi_join(null_2, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
      
      null_2_props[null_run] <- nrow(null_2_overlap) / n_novel_total
    } else {
      null_2_props[null_run] <- NA
    }
  }
  
  null_1_pvalue <- mean(null_1_props >= actual_recovery_prop)
  null_2_pvalue <- mean(null_2_props >= actual_recovery_prop, na.rm = TRUE)
  
  list(
    summary = data.frame(
      n_novel_total              = n_novel_total,
      n_predicted                = nrow(filtered_interactions),
      n_draw                     = n_draw,
      threshold                  = lower_threshold,
      actual_overlap_n           = nrow(matches_novel),
      actual_recovery_prop       = actual_recovery_prop,
      n_novel_forbidden          = n_novel_forbidden,
      n_novel_below_threshold    = n_novel_below_threshold,
      null1_mean_prop            = mean(null_1_props),
      null1_sd_prop              = sd(null_1_props),
      null1_pvalue               = null_1_pvalue,
      null2_mean_prop            = mean(null_2_props, na.rm = TRUE),
      null2_sd_prop              = sd(null_2_props, na.rm = TRUE),
      null2_pvalue                = null_2_pvalue
    ),
    null_1_props = null_1_props,
    null_2_props = null_2_props,
    novel_interactions = ext_validation_novel,
    recovered_interactions = matches_novel
  )
}

## for truly novel, and for new to ecoregion 

external_validation_results <- run_external_validation_flowerVisitation(
  predictions_visitsFlowersOf = predictions_visitsFlowersOf_Validation1,
  ext_validation_novel = ext_validation_novel,
  spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
  plant_checklist = plant_checklist,
  pollinator_checklists = pollinator_checklists,
  ecoregions = ecoregions,
  n_null_runs = 100,
  seed = 1
)

external_validation_results_eco <- run_external_validation_flowerVisitation(
  predictions_visitsFlowersOf = predictions_visitsFlowersOf_Validation1,
  ext_validation_novel = ext_validation_novel_eco,
  spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
  plant_checklist = plant_checklist,
  pollinator_checklists = pollinator_checklists,
  ecoregions = ecoregions,
  n_null_runs = 100,
  seed = 1
)

print(external_validation_results$summary)
print(external_validation_results_eco$summary)

saveRDS(external_validation_results, "Data_Clean/Validation/external_validation_results.rds")
saveRDS(external_validation_results_eco, "Data_Clean/Validation/external_validation_results_eco.rds")

#### Recall rate summary per external dataset
summarize_recall_by_dataset <- function(results, ext_validation_withDataset) {
  recovered_keys <- results$recovered_interactions %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct() %>%
    mutate(recovered = TRUE)
  
  ext_validation_withDataset %>%
    left_join(recovered_keys, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    mutate(recovered = replace_na(recovered, FALSE)) %>%
    group_by(dataset) %>%
    summarize(
      n_novel     = n(),
      n_recovered = sum(recovered),
      recall_rate = n_recovered / n_novel,
      .groups = "drop"
    )
}

recall_by_dataset <- summarize_recall_by_dataset(external_validation_results, ext_validation_novel_withDataset)
recall_by_dataset_eco <- summarize_recall_by_dataset(external_validation_results_eco, ext_validation_novel_eco_withDataset)

print(recall_by_dataset)
print(recall_by_dataset_eco)

write.csv(recall_by_dataset, "Data_Clean/Validation/external_validation_recall_by_dataset.csv", row.names = FALSE)
write.csv(recall_by_dataset_eco, "Data_Clean/Validation/external_validation_recall_by_dataset_eco.csv", row.names = FALSE)

#### Plot: recall by dataset, species-level vs. ecoregion-level novelty ####

recall_by_dataset_combined <- bind_rows(
  recall_by_dataset %>% mutate(novelty = "Novel (species pair)"),
  recall_by_dataset_eco %>% mutate(novelty = "Novel (ecoregion-specific)")
)

p_recall_by_dataset <- ggplot(recall_by_dataset_combined, aes(x = dataset, y = recall_rate, fill = novelty, color = novelty)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(n_recovered, "/", n_novel)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3, color = "black") +
  scale_fill_manual(values = c("Novel (species pair)" = "#331E0F", "Novel (ecoregion-specific)" = "#FFCC00"),
                    name = NULL) +
  scale_color_manual(values = c("Novel (species pair)" = "#331E0F", "Novel (ecoregion-specific)" = NA),
                     guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  labs(x = "External dataset", y = "Recall rate") +
  theme_cowplot() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))

print(p_recall_by_dataset)

ggsave("Figures/Supplementary/external_validation_recall_by_dataset.png",
       p_recall_by_dataset,
       width = 200,
       height = 140,
       units = "mm",
       dpi = 2000)

#### Plot: recall by dataset, species-level vs. ecoregion-level novelty ####

make_external_validation_plot <- function(results, x_ax) {
  plot_data <- bind_rows(
    tibble(model = "Null 1: Completely Random",
           recovery_prop = results$null_1_props),
    tibble(model = "Null 2: Phylogenetically Constrained",
           recovery_prop = results$null_2_props)
  )
  
  actual_prop <- results$summary$actual_recovery_prop
  
  ggplot(plot_data, aes(x = recovery_prop, fill = model)) +
    geom_histogram(bins = 30, color = "#331E0F", alpha = 0.7, position = "identity") +
    geom_vline(xintercept = actual_prop, linetype = "dashed", color = "#FFCC00", linewidth = 1) +
    scale_fill_manual(
      values = c("Null 1: Completely Random" = "gray70",
                 "Null 2: Phylogenetically Constrained" = "gray40"),
      name = NULL
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      x = x_ax,
      y = "Count (null draws)") +
    theme_cowplot() +
    theme(legend.position = "top")
}

p_ext_species <- make_external_validation_plot(external_validation_results, "Recall Rate of novel external interactions \n(novel species-pair level)")
p_ext_eco     <- make_external_validation_plot(external_validation_results_eco, "Recall Rate of novel external interactions \n(novel ecoregion level)")

panel_theme_tweak <- theme(axis.title = element_text(size = 11),
                           axis.text.x = element_text(size = 10),
                           axis.text.y = element_text(size = 10),
                           plot.title = element_text(size = 12),
                           plot.subtitle = element_text(size = 9))

p_ext_species <- p_ext_species + panel_theme_tweak
p_ext_eco     <- p_ext_eco + panel_theme_tweak

combined_external_validation <- (p_ext_species / p_ext_eco) +
  plot_annotation(tag_levels = 'A')

print(combined_external_validation)

ggsave("Figures/Supplementary/external_validation_recall.png",
       combined_external_validation,
       width = 200,
       height = 150,
       units = "mm",
       dpi = 2000)

###########################################################################
#### Validation 4 - Threshold Sensitivity  ####
###########################################################################
predictions_visitsFlowersOf_Validation2 <- predictions_visitsFlowersOf %>%
  filter(method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>%
  unique()

### Get confirmed_probs (interactionScores of known/observed links)
observed_pairs <- spatialInteractions_unique_visitsFlowersOf %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  distinct() %>%
  mutate(observedInteraction = 1)

predictions_with_obs_full <- predictions_visitsFlowersOf_Validation1 %>%
  left_join(observed_pairs, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
  mutate(observedInteraction = replace_na(observedInteraction, 0))

confirmed_probs <- predictions_with_obs_full %>%
  filter(observedInteraction == 1) %>%
  pull(interactionScore)

### Threshold quantiles to compare
quantile_grid <- seq(0.1, 0.40, by = 0.05)

thresholds_df <- data.frame(
  quantile_level = quantile_grid,
  threshold_value = quantile(confirmed_probs, quantile_grid, na.rm = TRUE)
)

print(thresholds_df)

grid_labels <- paste0("q", quantile_grid * 100)

interactions_by_threshold <- map_dfr(seq_len(nrow(thresholds_df)), function(i) {
  predictions_visitsFlowersOf_Validation2 %>%
    filter(interactionScore >= thresholds_df$threshold_value[i]) %>%
    count(method_type, name = "n_interactions") %>%
    mutate(
      threshold_label = grid_labels[i],
      quantile_level  = thresholds_df$quantile_level[i],
      threshold_value = thresholds_df$threshold_value[i]
    )
})

interactions_by_threshold_summary <- interactions_by_threshold %>%
  group_by(threshold_label, quantile_level, threshold_value) %>%
  mutate(n_total      = sum(n_interactions),
         pct_of_total = n_interactions / n_total * 100) %>%
  ungroup() %>%
  arrange(quantile_level, method_type)

print(interactions_by_threshold_summary)
write.csv(interactions_by_threshold_summary,
          "Tables/InteractionCounts_by_Threshold_MethodType.csv",
          row.names = FALSE)


### Ecoregion pollinator diversity supported using data from each threshold
compute_pols_supported <- function(threshold_value, predictions_visitsFlowersOf) {
  predictions_visitsFlowersOf %>%
    filter(interactionScore >= threshold_value) %>%
    group_by(JEPCODE, targetTaxonName_harm) %>%
    summarize(n_pols_supported = n_distinct(sourceTaxonName_harm), .groups = "drop")
}

support_by_threshold <- map(
  thresholds_df$threshold_value,
  ~ compute_pols_supported(.x, predictions_visitsFlowersOf_Validation2)
)

names(support_by_threshold) <- grid_labels  # grid_labels already derived from quantile_grid above

### Rank plants (per region) by total # pols supported
rank_plants <- function(support_df) {
  support_df %>%
    group_by(JEPCODE) %>%
    arrange(desc(n_pols_supported), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

ranked_by_threshold <- map(support_by_threshold, rank_plants)

### function to get top-N plant sets per ecoregion
get_top_n_sets <- function(ranked_df, n) {
  ranked_df %>%
    filter(rank <= n) %>%
    group_by(JEPCODE) %>%
    summarize(top_plants = list(targetTaxonName_harm), .groups = "drop")
}

### Compare top-50 sets between each threshold and baseline (25% threshold)
baseline_label <- "q25"

comparison_labels <- setdiff(names(ranked_by_threshold), baseline_label)

compare_topN_overlap <- function(n, baseline_label, comparison_labels, ranked_by_threshold) {
  baseline_sets <- get_top_n_sets(ranked_by_threshold[[baseline_label]], n)
  
  map_dfr(comparison_labels, function(lbl) {
    comp_sets <- get_top_n_sets(ranked_by_threshold[[lbl]], n)
    
    joined <- full_join(baseline_sets, comp_sets, by = "JEPCODE",
                        suffix = c("_baseline", "_comparison"))
    
    joined %>%
      rowwise() %>%
      mutate(
        n_baseline   = length(top_plants_baseline),
        n_comparison = length(top_plants_comparison),
        n_shared     = length(intersect(top_plants_baseline, top_plants_comparison)),
        pct_retained = n_shared / n_baseline * 100
      ) %>%
      ungroup() %>%
      mutate(comparison = lbl, top_n = n) %>%
      dplyr::select(JEPCODE, top_n, comparison, n_baseline, n_comparison,
                    n_shared, pct_retained)
  })
}

top50_comparison <- compare_topN_overlap(50, baseline_label, comparison_labels, ranked_by_threshold)
top25_comparison <- compare_topN_overlap(25, baseline_label, comparison_labels, ranked_by_threshold)

topN_comparison_all <- bind_rows(top50_comparison, top25_comparison)

topN_comparison_all %>% 
  group_by(comparison, top_n) %>% 
  summarise(mean_pct_retained = mean(pct_retained)) %>% 
  View()

### Summarise across ecoregions
topN_summary <- topN_comparison_all %>%
  group_by(top_n, comparison) %>%
  summarize(
    mean_pct_retained   = mean(pct_retained, na.rm = TRUE),
    median_pct_retained = median(pct_retained, na.rm = TRUE),
    sd_pct_retained      = sd(pct_retained, na.rm = TRUE),
    min_pct_retained     = min(pct_retained, na.rm = TRUE),
    max_pct_retained     = max(pct_retained, na.rm = TRUE),
    n_ecoregions         = n(),
    .groups = "drop"
  )

print(topN_summary)

#### Rank correlation and rank-biased overlap (RBO) of top 50 plants ####
compute_rank_correlation_n <- function(n, baseline_label, comparison_labels, ranked_by_threshold) {
  baseline_df <- ranked_by_threshold[[baseline_label]] %>%
    filter(rank <= n)
  
  map_dfr(comparison_labels, function(lbl) {
    comp_df <- ranked_by_threshold[[lbl]] %>%
      filter(rank <= n)
    
    merged <- inner_join(
      baseline_df %>%
        dplyr::select(JEPCODE, targetTaxonName_harm, rank_baseline = rank),
      comp_df %>%
        dplyr::select(JEPCODE, targetTaxonName_harm, rank_comparison = rank),
      by = c("JEPCODE", "targetTaxonName_harm")
    )
    
    merged %>%
      group_by(JEPCODE) %>%
      summarize(
        spearman_rho = cor(rank_baseline, rank_comparison, method = "spearman"),
        n_shared_species = n(),
        .groups = "drop"
      ) %>%
      mutate(comparison = lbl, top_n = n)
  })
}

rank_correlation_top25 <- compute_rank_correlation_n(25, baseline_label, comparison_labels, ranked_by_threshold)
rank_correlation_top50 <- compute_rank_correlation_n(50, baseline_label, comparison_labels, ranked_by_threshold)

### Rank correlation among ALL plants 
compute_rank_correlation_all_plants <- function(baseline_label, comparison_labels, ranked_by_threshold) {
  baseline_df <- ranked_by_threshold[[baseline_label]]
  
  map_dfr(comparison_labels, function(lbl) {
    comp_df <- ranked_by_threshold[[lbl]]
    
    merged <- inner_join(
      baseline_df %>%
        dplyr::select(JEPCODE, targetTaxonName_harm, rank_baseline = rank),
      comp_df %>%
        dplyr::select(JEPCODE, targetTaxonName_harm, rank_comparison = rank),
      by = c("JEPCODE", "targetTaxonName_harm")
    )
    
    merged %>%
      group_by(JEPCODE) %>%
      summarize(
        spearman_rho = cor(rank_baseline, rank_comparison, method = "spearman"),
        n_shared_species = n(),
        .groups = "drop"
      ) %>%
      mutate(comparison = lbl, top_n = "All")
  })
}

rank_correlation_all_plants <- compute_rank_correlation_all_plants(baseline_label, comparison_labels, ranked_by_threshold)

rank_correlation_all <- bind_rows(
  rank_correlation_top25 %>%
    mutate(top_n = as.character(top_n)),
  rank_correlation_top50 %>%
    mutate(top_n = as.character(top_n)),
  rank_correlation_all_plants
)

rank_correlation_all %>% 
  group_by(comparison, top_n) %>% 
  summarise(mean_rho = mean(spearman_rho)) %>% 
  View()

#### Plot: Validation 2 - threshold sensitivity ####
threshold_order <- thresholds_df %>%
  mutate(label = grid_labels) %>%
  arrange(quantile_level) %>%
  pull(label)

topN_comparison_all <- topN_comparison_all %>%
  mutate(comparison = factor(comparison, levels = threshold_order),
         top_n = factor(top_n, levels = c(25, 50)))

rank_correlation_all <- rank_correlation_all %>%
  mutate(comparison = factor(comparison, levels = threshold_order),
         top_n = factor(top_n, levels = c("25", "50", "All")))

quantile_lookup <- thresholds_df %>%
  mutate(label = grid_labels) %>%
  dplyr::select(label, quantile_level)

axis_labels <- setNames(
  paste0(quantile_lookup$quantile_level * 100, "%"),
  quantile_lookup$label
)

### Panel A - % of baseline (25%) top-50 plants consistent
p_retention <- ggplot(topN_comparison_all, aes(x = comparison, y = pct_retained, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00"), name = "Top-N") +
  labs(
    x = "Comparison Threshold",
    y = "% of top N plants shared with baseline") +
  theme_cowplot()

# Panel B - Rank correlation of top N plants (plus all plants)
p_rank_corr <- ggplot(rank_correlation_all, aes(x = comparison, y = spearman_rho, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00", "All" = "#B58A00"), name = "Top-N") +
  labs(
    x = "Comparison Threshold",
    y = "Spearman rank correlation with baseline") +
  theme_cowplot()

#### Combine panels ####
panel_theme_tweak <- theme(axis.title = element_text(size = 12),
                           axis.text = element_text(size = 10),
                           plot.title = element_text(size = 11),
                           legend.position = "top")

p_retention <- p_retention + panel_theme_tweak
p_rank_corr <- p_rank_corr + panel_theme_tweak

combined_threshold_sensitivity <- (p_retention | p_rank_corr) +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/Supplementary/threshold_sensitivity.png",
       combined_threshold_sensitivity,
       width = 200,
       height = 130,
       units = "mm",
       dpi = 2000)


####

## Average Spearman rho per ecoregion, across all thresholds tested above
ecoregion_sensitivity <- rank_correlation_all %>%
  filter(top_n == "All") %>%
  group_by(JEPCODE) %>%
  summarize(
    mean_spearman_rho = mean(spearman_rho, na.rm = TRUE),
    sd_spearman_rho    = sd(spearman_rho, na.rm = TRUE),
    n_comparisons      = n(),
    .groups = "drop"
  )

## Join to ecoregion polygons
ecoregions_map_data <- ecoregions %>%
  group_by(JEPCODE) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid() %>%
  left_join(ecoregion_sensitivity, by = "JEPCODE")

## Map
p_map_sensitivity <- ggplot(ecoregions_map_data) +
  geom_sf(aes(fill = mean_spearman_rho), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    high = "white",
    mid = "#F9B753",
    low = "#331E0F",
    midpoint = 0.94,
    name = "Mean Spearman ρ\nvs. baseline",
    na.value = "grey85"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, hjust = 0.5)
  )

ggsave("Figures/Supplementary/ecoregion_threshold_sensitivity_map.png",
       p_map_sensitivity,
       width = 180,
       height = 150,
       units = "mm",
       dpi = 2000)


###########################################################################
#### Validation 5 - Predicted interaction data scope used to construct plant rankings ####
###########################################################################
predictions_default <- predictions_visitsFlowersOf %>%
  filter(method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>%
  unique()

predictions_default_with_obs <- predictions_default %>%
  left_join(observed_pairs, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
  mutate(observedInteraction = replace_na(observedInteraction, 0))

fixed_threshold <- predictions_default_with_obs %>%
  filter(observedInteraction == 1) %>%
  pull(interactionScore) %>%
  quantile(0.25, na.rm = TRUE)

## Compare three method_type filters
# 1. Only interactions reconstructed with phylogenetic constraint
# 2. (default) phylo constraint plus graph embedding out of sample (pollinators with no pre-existing interaction data)
# 3. Phylo constraint plus graph embedding out of sample (pollinators with no pre-existing interaction data) plus graph embedding in sample (addt. potential plant genera inferred for pollinators with existing data)
methodtype_sets <- list(
  Level1_only   = c("Level1_Simple"),
  Level1_Level2 = c("Level1_Simple", "Level2_GEOOS"),  # default / baseline
  All_methods   = NULL                                  # no filter (adds Level3_GEOOSIS)
)

build_ranked_by_methodtype <- function(method_types, predictions_visitsFlowersOf, threshold_value) {
  preds_subset <- if (is.null(method_types)) {
    predictions_visitsFlowersOf %>% unique()
  } else {
    predictions_visitsFlowersOf %>% filter(method_type %in% method_types) %>% unique()
  }
  compute_pols_supported(threshold_value, preds_subset) %>%
    rank_plants()
}

ranked_by_methodtype <- map(
  methodtype_sets,
  ~ build_ranked_by_methodtype(.x, predictions_visitsFlowersOf, fixed_threshold)
)

baseline_label_mt    <- "Level1_Level2"
comparison_labels_mt <- setdiff(names(ranked_by_methodtype), baseline_label_mt)

top50_comparison_mt <- compare_topN_overlap(50, baseline_label_mt, comparison_labels_mt, ranked_by_methodtype)
top25_comparison_mt <- compare_topN_overlap(25, baseline_label_mt, comparison_labels_mt, ranked_by_methodtype)
topN_comparison_mt  <- bind_rows(top50_comparison_mt, top25_comparison_mt)

topN_summary_mt <- topN_comparison_mt %>%
  group_by(top_n, comparison) %>%
  summarize(
    mean_pct_retained   = mean(pct_retained, na.rm = TRUE),
    median_pct_retained = median(pct_retained, na.rm = TRUE),
    sd_pct_retained     = sd(pct_retained, na.rm = TRUE),
    min_pct_retained    = min(pct_retained, na.rm = TRUE),
    max_pct_retained    = max(pct_retained, na.rm = TRUE),
    n_ecoregions        = n(),
    .groups = "drop"
  )
print(topN_summary_mt)

rank_correlation_top25_mt <- compute_rank_correlation_n(25, baseline_label_mt, comparison_labels_mt, ranked_by_methodtype)
rank_correlation_top50_mt <- compute_rank_correlation_n(50, baseline_label_mt, comparison_labels_mt, ranked_by_methodtype)
rank_correlation_all_mt   <- compute_rank_correlation_all_plants(baseline_label_mt, comparison_labels_mt, ranked_by_methodtype)

rank_correlation_methodtype <- bind_rows(
  rank_correlation_top25_mt %>% mutate(top_n = as.character(top_n)),
  rank_correlation_top50_mt %>% mutate(top_n = as.character(top_n)),
  rank_correlation_all_mt
)

rank_correlation_methodtype %>%
  group_by(comparison, top_n) %>%
  summarise(mean_rho = mean(spearman_rho), .groups = "drop") %>%
  print()

#### Plot: method-type sensitivity ####
method_order  <- c("Level1_only", "All_methods")
method_labels <- c(
  Level1_only = "Phylogenetic\nconstraint only",
  All_methods = "All methods\n(+ GE in-sample\nlearning)"
)

topN_comparison_mt <- topN_comparison_mt %>%
  mutate(comparison = factor(comparison, levels = method_order),
         top_n = factor(top_n, levels = c(25, 50)))
rank_correlation_methodtype <- rank_correlation_methodtype %>%
  mutate(comparison = factor(comparison, levels = method_order),
         top_n = factor(top_n, levels = c("25", "50", "All")))

p_retention_mt <- ggplot(topN_comparison_mt, aes(x = comparison, y = pct_retained, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00"), name = "Top-N") +
  labs(x = "Prediction methods included", y = "% of top N plants shared") +
  theme_cowplot()

p_rank_corr_mt <- ggplot(rank_correlation_methodtype, aes(x = comparison, y = spearman_rho, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00", "All" = "#B58A00"), name = "Top-N") +
  labs(x = "Prediction methods included", y = "Spearman rank correlation") +
  theme_cowplot()

panel_theme_tweak_mt <- theme(axis.title = element_text(size = 12),
                              axis.text = element_text(size = 10),
                              legend.position = "top")
p_retention_mt <- p_retention_mt + panel_theme_tweak_mt
p_rank_corr_mt <- p_rank_corr_mt + panel_theme_tweak_mt

combined_methodtype_sensitivity <- (p_retention_mt | p_rank_corr_mt) +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/Supplementary/methodtype_sensitivity.png",
       combined_methodtype_sensitivity,
       width = 200, height = 130, units = "mm", dpi = 2000)

###########################################################################
#### Validation 6 - Method-type sensitivity within each pollinator taxon group ####
###########################################################################

## Taxon lookup for pollinator species (bees, butterflies, moths, hoverflies)
taxon_lookup <- checklist_cleaned %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies")) %>%
  dplyr::select(genus_species, taxon) %>%
  unique()

taxon_groups <- c("bees", "butterflies", "moths", "hoverflies")

## Extend build_ranked_by_methodtype to optionally restrict to one pollinator taxon
build_ranked_by_methodtype_taxon <- function(method_types, predictions_visitsFlowersOf,
                                             threshold_value, taxon_group = NULL,
                                             taxon_lookup = NULL) {
  preds_subset <- if (is.null(method_types)) {
    predictions_visitsFlowersOf %>% unique()
  } else {
    predictions_visitsFlowersOf %>% filter(method_type %in% method_types) %>% unique()
  }
  
  if (!is.null(taxon_group)) {
    taxon_species <- taxon_lookup %>%
      filter(taxon == taxon_group) %>%
      pull(genus_species)
    
    preds_subset <- preds_subset %>%
      filter(sourceTaxonName_harm %in% taxon_species)
  }
  
  compute_pols_supported(threshold_value, preds_subset) %>%
    rank_plants()
}

## Build ranked plant lists for every combination of taxon group x method-type set
ranked_by_methodtype_taxon <- map(taxon_groups, function(tg) {
  map(
    methodtype_sets,
    ~ build_ranked_by_methodtype_taxon(.x, predictions_visitsFlowersOf, fixed_threshold,
                                       taxon_group = tg, taxon_lookup = taxon_lookup)
  )
})
names(ranked_by_methodtype_taxon) <- taxon_groups

## Run the same top-N overlap and rank-correlation comparisons within each taxon group,
## reusing baseline_label_mt / comparison_labels_mt already defined above
topN_comparison_taxon <- map_dfr(taxon_groups, function(tg) {
  ranked_list <- ranked_by_methodtype_taxon[[tg]]
  top50 <- compare_topN_overlap(50, baseline_label_mt, comparison_labels_mt, ranked_list)
  top25 <- compare_topN_overlap(25, baseline_label_mt, comparison_labels_mt, ranked_list)
  bind_rows(top50, top25) %>% mutate(taxon = tg)
})

topN_summary_taxon <- topN_comparison_taxon %>%
  group_by(taxon, top_n, comparison) %>%
  summarize(
    mean_pct_retained   = mean(pct_retained, na.rm = TRUE),
    median_pct_retained = median(pct_retained, na.rm = TRUE),
    sd_pct_retained     = sd(pct_retained, na.rm = TRUE),
    min_pct_retained    = min(pct_retained, na.rm = TRUE),
    max_pct_retained    = max(pct_retained, na.rm = TRUE),
    n_ecoregions        = n(),
    .groups = "drop"
  )
print(topN_summary_taxon)

rank_correlation_taxon <- map_dfr(taxon_groups, function(tg) {
  ranked_list <- ranked_by_methodtype_taxon[[tg]]
  top25 <- compute_rank_correlation_n(25, baseline_label_mt, comparison_labels_mt, ranked_list)
  top50 <- compute_rank_correlation_n(50, baseline_label_mt, comparison_labels_mt, ranked_list)
  allp  <- compute_rank_correlation_all_plants(baseline_label_mt, comparison_labels_mt, ranked_list)
  bind_rows(
    top25 %>% mutate(top_n = as.character(top_n)),
    top50 %>% mutate(top_n = as.character(top_n)),
    allp
  ) %>% mutate(taxon = tg)
})

rank_correlation_taxon %>%
  group_by(taxon, comparison, top_n) %>%
  summarise(mean_rho = mean(spearman_rho), .groups = "drop") %>%
  print()

#### Plot: method-type sensitivity, faceted by taxon group ####
topN_comparison_taxon <- topN_comparison_taxon %>%
  mutate(comparison = factor(comparison, levels = method_order),
         top_n = factor(top_n, levels = c(25, 50)),
         taxon = str_to_title(taxon))

rank_correlation_taxon <- rank_correlation_taxon %>%
  mutate(comparison = factor(comparison, levels = method_order),
         top_n = factor(top_n, levels = c("25", "50", "All")),
         taxon = str_to_title(taxon))

p_retention_taxon <- ggplot(topN_comparison_taxon, aes(x = comparison, y = pct_retained, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  facet_wrap(~ taxon) +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00"), name = "Top-N") +
  labs(x = "Prediction methods included", y = "% of top N plants shared") +
  theme_cowplot()

p_rank_corr_taxon <- ggplot(rank_correlation_taxon, aes(x = comparison, y = spearman_rho, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ taxon) +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00", "All" = "#B58A00"), name = "Top-N") +
  labs(x = "Prediction methods included", y = "Spearman rank correlation") +
  theme_cowplot()

panel_theme_tweak_taxon <- theme(axis.title = element_text(size = 12),
                                 axis.text = element_text(size = 9),
                                 strip.text = element_text(size = 11, face = "bold"),
                                 legend.position = "top")
p_retention_taxon <- p_retention_taxon + panel_theme_tweak_taxon
p_rank_corr_taxon <- p_rank_corr_taxon + panel_theme_tweak_taxon

combined_methodtype_sensitivity_taxon <- (p_retention_taxon / p_rank_corr_taxon) +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/Supplementary/methodtype_sensitivity_by_taxon.png",
       combined_methodtype_sensitivity_taxon,
       width = 220, height = 220, units = "mm", dpi = 2000)

###########################################################################
#### Validation 7 - Top plant genera for both flower visitation and host ####
###########################################################################
predicted_interactions <- fread("Data_Clean/Interactions/Predicted_Interactions.csv") %>%
  filter(sourceTaxonType != "hummingbirds")

## Flower visitation
ranked_visitsFlowersOf <- predicted_interactions %>%
  filter(interactionTypeName == "visitsFlowersOf",
         method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>%
  group_by(JEPCODE, targetTaxonName_harm) %>%
  summarize(n_pols_supported = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
  rank_plants()

ranked_hasHost <- predicted_interactions %>%
  filter(interactionTypeName == "hasHost") %>%
  group_by(JEPCODE, targetTaxonName_harm) %>%
  summarize(n_pols_supported = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
  rank_plants()

## Top 25 per ecoregion
top25_visitsFlowersOf <- ranked_visitsFlowersOf %>%
  filter(rank <= 25) %>%
  mutate(genus = word(targetTaxonName_harm, 1))

top25_hasHost <- ranked_hasHost %>%
  filter(rank <= 25) %>%
  mutate(genus = word(targetTaxonName_harm, 1))

## Genus-level summary
summarize_top_genera <- function(top25_df, n_top_genera = 5) {
  top25_df %>%
    group_by(genus) %>%
    summarize(
      n_total_appearances       = n(),
      n_distinct_species        = n_distinct(targetTaxonName_harm),
      n_ecoregions_represented  = n_distinct(JEPCODE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_total_appearances)) %>%
    slice_head(n = n_top_genera)
}

top_genera_visitsFlowersOf <- summarize_top_genera(top25_visitsFlowersOf, 5) %>%
  mutate(interaction_type = "Flower Visitation")

top_genera_hasHost <- summarize_top_genera(top25_hasHost, 5) %>%
  mutate(interaction_type = "Host Use")

top_genera_summary <- bind_rows(top_genera_visitsFlowersOf, top_genera_hasHost)

print(top_genera_summary)
write.csv(top_genera_summary, "Tables/Top5PlantGenera_FromTop25PlantsPerRegion_Summary.csv", row.names = FALSE)


###########################################################################
#### Validation 8 - Plant ranking sensitivity to data poor taxa removal ####
###########################################################################
species_poorBoyce <- read.csv("Data_Clean/Validation/species_boyce_0.25_or_lower.csv")
species_ultrarare <- read.csv("Data_Clean/Validation/ultra_rare_species.csv")

exclude_poorBoyce <- species_poorBoyce$species
exclude_ultrarare <- species_ultrarare$species

## Default flower visitation network: Level1_Simple + Level2_GEOOS
predictions_visitsFlowersOf_default <- predicted_interactions %>%
  filter(interactionTypeName == "visitsFlowersOf",
         method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>%
  unique()

## Three candidate species pools: default (no exclusion), excluding poor-Boyce
## species, excluding ultra-rare (data-poor) species. Exclusion applies to
## EITHER side of the interaction (plant or pollinator).
species_sets <- list(
  Default           = character(0),
  Exclude_PoorBoyce = exclude_poorBoyce,
  Exclude_UltraRare = exclude_ultrarare
)

build_ranked_by_speciesFilter <- function(exclude_species, predictions_default) {
  preds_subset <- if (length(exclude_species) == 0) {
    predictions_default
  } else {
    predictions_default %>%
      filter(!sourceTaxonName_harm %in% exclude_species,
             !targetTaxonName_harm %in% exclude_species)
  }
  
  preds_subset %>%
    group_by(JEPCODE, targetTaxonName_harm) %>%
    summarize(n_pols_supported = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
    rank_plants()
}

ranked_by_speciesFilter <- map(
  species_sets,
  ~ build_ranked_by_speciesFilter(.x, predictions_visitsFlowersOf_default)
)

## Reuse the same comparison functions from earlier validations
baseline_label_sf    <- "Default"
comparison_labels_sf <- setdiff(names(ranked_by_speciesFilter), baseline_label_sf)

top50_comparison_sf <- compare_topN_overlap(50, baseline_label_sf, comparison_labels_sf, ranked_by_speciesFilter)
top25_comparison_sf <- compare_topN_overlap(25, baseline_label_sf, comparison_labels_sf, ranked_by_speciesFilter)
topN_comparison_sf  <- bind_rows(top50_comparison_sf, top25_comparison_sf)

topN_summary_sf <- topN_comparison_sf %>%
  group_by(top_n, comparison) %>%
  summarize(
    mean_pct_retained   = mean(pct_retained, na.rm = TRUE),
    median_pct_retained = median(pct_retained, na.rm = TRUE),
    sd_pct_retained     = sd(pct_retained, na.rm = TRUE),
    min_pct_retained    = min(pct_retained, na.rm = TRUE),
    max_pct_retained    = max(pct_retained, na.rm = TRUE),
    n_ecoregions        = n(),
    .groups = "drop"
  )
print(topN_summary_sf)

rank_correlation_top25_sf <- compute_rank_correlation_n(25, baseline_label_sf, comparison_labels_sf, ranked_by_speciesFilter)
rank_correlation_top50_sf <- compute_rank_correlation_n(50, baseline_label_sf, comparison_labels_sf, ranked_by_speciesFilter)
rank_correlation_all_sf   <- compute_rank_correlation_all_plants(baseline_label_sf, comparison_labels_sf, ranked_by_speciesFilter)

rank_correlation_speciesFilter <- bind_rows(
  rank_correlation_top25_sf %>% mutate(top_n = as.character(top_n)),
  rank_correlation_top50_sf %>% mutate(top_n = as.character(top_n)),
  rank_correlation_all_sf
)

rank_correlation_speciesFilter %>%
  group_by(comparison, top_n) %>%
  summarise(mean_rho = mean(spearman_rho), .groups = "drop") %>%
  print()

#### Plot ####
sf_order  <- c("Exclude_PoorBoyce", "Exclude_UltraRare")
sf_labels <- c(
  Exclude_PoorBoyce = "Excluding\npoor-Boyce\nSDM species",
  Exclude_UltraRare = "Excluding\nrare species"
)

topN_comparison_sf <- topN_comparison_sf %>%
  mutate(comparison = factor(comparison, levels = sf_order),
         top_n = factor(top_n, levels = c(25, 50)))
rank_correlation_speciesFilter <- rank_correlation_speciesFilter %>%
  mutate(comparison = factor(comparison, levels = sf_order),
         top_n = factor(top_n, levels = c("25", "50", "All")))

p_retention_sf <- ggplot(topN_comparison_sf, aes(x = comparison, y = pct_retained, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = sf_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00"), name = "Top-N") +
  labs(x = NULL, y = "% of top N plants shared with default") +
  theme_cowplot()

p_rank_corr_sf <- ggplot(rank_correlation_speciesFilter, aes(x = comparison, y = spearman_rho, fill = top_n)) +
  geom_boxplot(color = "#331E0F", outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_discrete(labels = sf_labels) +
  scale_fill_manual(values = c("25" = "#FFE999", "50" = "#FFCC00", "All" = "#B58A00"), name = "Top-N") +
  labs(x = NULL, y = "Spearman rank correlation with default") +
  theme_cowplot()

panel_theme_tweak_sf <- theme(axis.title = element_text(size = 12),
                              axis.text = element_text(size = 10),
                              legend.position = "top")
p_retention_sf <- p_retention_sf + panel_theme_tweak_sf
p_rank_corr_sf <- p_rank_corr_sf + panel_theme_tweak_sf

combined_speciesFilter_sensitivity <- (p_retention_sf | p_rank_corr_sf) +
  plot_annotation(tag_levels = 'A')

ggsave("Figures/Supplementary/speciesFilter_sensitivity.png",
       combined_speciesFilter_sensitivity,
       width = 200, height = 130, units = "mm", dpi = 2000)


rank_correlation_speciesFilter %>%
  filter(comparison == "Exclude_PoorBoyce", top_n == "All") %>%
  summarise(median_rho = median(spearman_rho, na.rm = TRUE))

topN_summary_sf %>%
  filter(comparison == "Exclude_PoorBoyce") %>%
  dplyr::select(top_n, mean_pct_retained)

