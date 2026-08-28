#### Load packages ####
library(sf)
library(terra)
library(dplyr)

#### Read in Data ####
checklist_cleaned <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

## Load jepson ecoregions
ecoregions <- read_sf("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid()

## interaction observations
interactions <- read.csv("Data_Clean/Interactions/ints_final_clean_slim.csv")

## Pre-processed spatial/pheno overlap for predictions
predict_visitsFlowersOf <- read.csv("Data_Clean/Interactions/predictions_supporting_visitsFlowersOf.csv")
predict_hasHost <- read.csv("Data_Clean/Interactions/predictions_supporting_hasHost.csv")

#### Data Prep ####
## Clean up checklist for here
checklist_pollinator_simple <- checklist_cleaned %>%
  filter(taxon != "plants") %>%
  dplyr::select(genus_species, taxon) %>%
  rename(sourceTaxonType = taxon) %>%
  unique()

## get spatially explicit interactions observed in CA
# Convert to numeric (non-numeric values will become NA)
interactions$decimalLatitude <- as.numeric(interactions$decimalLatitude)
interactions$decimalLongitude <- as.numeric(interactions$decimalLongitude)

# create sf object
interactions_sf <- st_as_sf(
  interactions[!is.na(interactions$decimalLatitude) & 
                 !is.na(interactions$decimalLongitude),], 
  coords = c("decimalLongitude", "decimalLatitude"), 
  crs = 4326,
  remove = FALSE
)

# Spatial filter - keep only points that intersect with regions
interactions_regional <- st_join(interactions_sf, ecoregions) %>%
  filter(!is.na(JEPCODE)) %>%
  st_drop_geometry()

#### Normalize Spatial Overlap ####
## Create rast template covering all ecoregions
# Get resolution from an example SDM
sdm_example <- rast(list.files("Data_Clean/SDMs/sdm_by_taxon/plants/continuous", full.names = T)[1])
res_value <- res(sdm_example)[1]  # Assuming square pixels
crs_value <- crs(sdm_example)

# merge all ecoregions
ecoregions_union_full <- ecoregions %>%
  st_union() %>%
  st_make_valid()

# Convert to SpatVector
ecoregions_vect_full <- vect(as(ecoregions_union_full, "Spatial"))

# Create a new raster from scratch covering the full extent of all ecoregions
ones_rast <- rast(ext(ecoregions_vect_full), resolution = res_value, crs = crs_value)
values(ones_rast) <- 1

# Mask to ecoregions boundary
ones_rast <- mask(ones_rast, ecoregions_vect_full)

## Count pixels per JEPCODE
# Group by JEPCODE and union geometries
ecoregions_union <- ecoregions %>%
  group_by(JEPCODE) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# Get unique JEPCODES
unique_jepcodes <- ecoregions_union$JEPCODE

# empty pixel counts df
pixel_counts <- data.frame(
  JEPCODE = unique_jepcodes,
  n_pixels = NA_integer_
)

# Count pixels for each JEPCODE
for (i in seq_along(unique_jepcodes)) {
  eco_single <- ecoregions_union %>% filter(JEPCODE == unique_jepcodes[i])
  eco_sp <- as(eco_single, "Spatial")
  eco_vect <- vect(eco_sp)
  
  masked <- terra::mask(terra::crop(ones_rast, eco_vect), eco_vect)
  pixel_counts$n_pixels[i] <- sum(!is.na(values(masked)))
}

## Join and normalize spatial overlap
# For flower visitation predictions
predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

# For host predictions
predict_hasHost <- predict_hasHost %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

#### Normalize temporal overlap ####
## Divide temporal overlap by plant flowering duration
# For flower visitation predictions
predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  mutate(temporal_overlap_normalized = temporal_overlap / plant_duration)

#### Predict flower visitation interactions ####
# Prediction using switch functions through case_when, temporal and spatial overlap
predictions_visitsFlowersOf <- predict_visitsFlowersOf %>% 
  # interaction prediction
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    temporal_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized * temporal_overlap_normalized
  ))

## Remove the handful of cases where interactionScore is NA due to missing spatial overlap value
predictions_visitsFlowersOf_NA <- predictions_visitsFlowersOf[is.na(predictions_visitsFlowersOf$interactionScore),]
predictions_visitsFlowersOf <- predictions_visitsFlowersOf[!is.na(predictions_visitsFlowersOf$interactionScore),]

# observed nectaring interactions
int_observed_visitsFlowersOf <- interactions_regional %>% 
  filter(interactionTypeName == "visitsFlowersOf") %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>% 
  mutate(observedInteraction = 1) %>%
  unique()

# Add in observed interactions to new columns
predictions_visitsFlowersOf <- predictions_visitsFlowersOf %>%
  left_join(int_observed_visitsFlowersOf) %>%
  mutate(observedInteraction = if_else(!is.na(observedInteraction), 1, 0)) %>%
  # Add type of interaction
  mutate(interactionTypeName = "visitsFlowersOf")

# Filter for confirmed observations to help select threshold for "likely" interactions
confirmedScores_visitsFlowersOf <- predictions_visitsFlowersOf$interactionScore[
  predictions_visitsFlowersOf$observedInteraction == 1
]

# Find a reasonable threshold for distinguishing "Potential" interactions - perhaps the 25th percentile of confirmed
lower_threshold_visitsFlowersOf <- quantile(confirmedScores_visitsFlowersOf, 0.25)

# Filter to keep high-probability interactions OR confirmed interactions
predictions_visitsFlowersOf_filtered <- predictions_visitsFlowersOf[
  predictions_visitsFlowersOf$interactionScore >= lower_threshold_visitsFlowersOf |
    predictions_visitsFlowersOf$observedInteraction == 1,
]

predictions_visitsFlowersOf_filtered <- 
  predictions_visitsFlowersOf_filtered[,c("targetTaxonName_harm",
                                          "sourceTaxonName_harm",
                                          "JEPCODE",
                                          "method_type",
                                          "interactionTypeName",
                                          "observedInteraction")]


#### Host Interaction Predictions ####
### For host plants, only consider host plant interactions that have been observed anywhere before
# If observed in that ecoregion, receives value of 1
# If observed elsewhere but not in ecoregion, probability based amount of area of co-occurrence
# If not observed anywhere, recieves value of 0

predictions_hasHost <- predict_hasHost %>% 
  # interaction prediction
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized
  ))

## Remove the handful of cases where interactionScore is NA due to missing spatial overlap value
predictions_hasHost_NA <- predictions_hasHost[is.na(predictions_hasHost$interactionScore),]
predictions_hasHost <- predictions_hasHost[!is.na(predictions_hasHost$interactionScore),]

# observed nectaring interactions
int_observed_hasHost <- interactions_regional %>% 
  filter(interactionTypeName == "hasHost") %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>% 
  mutate(observedInteraction = 1) %>%
  unique()

# Add in observed interactions to new columns
predictions_hasHost <- predictions_hasHost %>%
  left_join(int_observed_hasHost) %>%
  mutate(observedInteraction = if_else(!is.na(observedInteraction), 1, 0)) %>%
  # Add type of interaction
  mutate(interactionTypeName = "hasHost")

# Filter for confirmed observations to help select threshold for "likely" interactions
confirmedScores_hasHost <- predictions_hasHost$interactionScore[
  predictions_hasHost$observedInteraction == 1
]

# Find a reasonable threshold for distinguishing "Potential" interactions - perhaps the 25th percentile of confirmed
lower_threshold_hasHost <- quantile(confirmedScores_hasHost, 0.25)

# Filter to keep high-probability interactions OR confirmed interactions
predictions_hasHost_filtered <- predictions_hasHost[
  predictions_hasHost$interactionScore >= lower_threshold_hasHost |
    predictions_hasHost$observedInteraction == 1,
]

predictions_hasHost_filtered <- 
  predictions_hasHost_filtered[,c("targetTaxonName_harm",
                                  "sourceTaxonName_harm",
                                  "JEPCODE",
                                  "method_type",
                                  "interactionTypeName",
                                  "observedInteraction")]

#### Final aggregation, cleaning ####
# Combine interaction predictions datasets
predicted_interactions <- bind_rows(predictions_visitsFlowersOf_filtered, predictions_hasHost_filtered)

# Join with checklists to get higher taxon
predicted_interactions <- left_join(predicted_interactions, checklist_pollinator_simple, by = join_by("sourceTaxonName_harm" == "genus_species"))

# Add calscape categorical column to distinguish "Confirmed" and "Potential"
predicted_interactions$categoricalInteraction <- ifelse(
  predicted_interactions$observedInteraction == 1, 
  "Confirmed", 
  "Potential"
)

#### Save Predicted interactions ####
write.csv(unique(predicted_interactions), "Data_Clean/Interactions/Predicted_Interactions.csv", row.names = F)
