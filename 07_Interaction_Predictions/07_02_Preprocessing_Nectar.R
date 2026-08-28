#### Load packages ####
library(sf)
library(terra)
library(dplyr)
library(stringr)
library(data.table)
library(parallel)

#### List which species have phenometrics and SDMs ####
### Phenology
phen_best <- read.csv("Data_Clean/Phenology/phenology_best_available.csv") 

species_phen <- phen_best %>%
  filter(!is.na(phen_start)) %>%
  pull(genus_species) %>%
  unique()

### SDMs
# Set base path
base_path <- "Data_Clean/SDMs/sdm_by_taxon"

# Get all taxon folders
taxon_folders <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)

# Initialize empty list to store results
results_list <- list()

# Loop through each taxon folder
for (taxon in taxon_folders) {
  
  # Build path to continuous sdm folder
  continuous_path <- file.path(base_path, taxon, "continuous")
  
  # Get all files in continuous folder
  all_files <- list.files(continuous_path, pattern = "\\.tif$", full.names = FALSE)
  
  # Remove the manifest file
  manifest_pattern <- paste0(taxon, "_continuous_manifest\\.csv")
  files_filtered <- all_files[!grepl(manifest_pattern, all_files)]
  
  # Extract species name
  species_names <- gsub("\\_continuous.tif$", "", files_filtered)
  
  # Replace underscores with spaces for species names
  species_names <- gsub("_", " ", species_names)
  
  # Create data frame for this taxon
  if (length(species_names) > 0) {
    taxon_df <- data.frame(
      species = species_names,
      taxon = taxon,
      stringsAsFactors = FALSE
    )
    results_list[[taxon]] <- taxon_df
  }
}

# Combine
sdms_manifest <- bind_rows(results_list)

### Checklist
checklist_cleaned <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# Filter phenometrics and sdms to only genus_species in checklist (should only remove subspecies from SDMs)
species_phen <- species_phen[species_phen %in% checklist_cleaned$genus_species]
species_sdms <- sdms_manifest$species[sdms_manifest$species %in% checklist_cleaned$genus_species]

checklist_prediction <- checklist_cleaned %>% 
  dplyr::select(genus_species, taxon, genus, specificEpithet) %>%
  unique() %>%
  filter(genus_species %in% species_phen & genus_species %in% species_sdms) 

# Get unique pollinator and plant species
plant_species <- checklist_prediction[checklist_prediction$taxon == "plants",]
pollinator_species <- checklist_prediction[checklist_prediction$taxon != "plants",]

#### Create long df that detail which interactions to try to predict based on GE and phylogenetic constraint ####
### Level 1: Simple - with just a plant genus phylogenetic constraint
# read in interaction observations
ints <- read.csv("Data_Clean/Interactions/ints_final_clean_slim.csv")

### Temporary until we model pollen interactions seperately
ints$interactionTypeName[ints$interactionTypeName == "collectsPollenOf"] <- "visitsFlowersOf"

# get unique pollinator species-plant genus interactions for flower visitation ints
ints_pol_plantgenera <- ints %>%
  # filter for only flower visitation|
  filter(interactionTypeName == "visitsFlowersOf") %>%
  dplyr::select(targetTaxonGenusName, sourceTaxonName_harm) %>%
  # remove rows where sourceTaxonName_harm is only one word (genus only)
  filter(str_count(sourceTaxonName_harm, "\\S+") >= 2) %>%
  # keep only first two words (Genus species)
  mutate(sourceTaxonName_harm = word(sourceTaxonName_harm, 1, 2))

# now filter for only plant genera and pol species in checklist
ints_pol_plantgenera_checklist <- ints_pol_plantgenera %>%
  filter(sourceTaxonName_harm %in% pollinator_species$genus_species &
           targetTaxonGenusName %in% plant_species$genus) %>%
  unique()

# now expand for all plant species in each genus to get pairwise species we want to predict over
# using phylo constraint
predictions_simple <- ints_pol_plantgenera_checklist %>%
  left_join(plant_species, by = c("targetTaxonGenusName" = "genus")) %>%
  # Remove rows where no matching plant species found
  filter(!is.na(genus_species)) %>%
  rename(targetTaxonName_harm = genus_species) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm) %>%
  mutate(method_type = "Level1_Simple")

### Graph embedding
bees_GE <- read.csv("Data_Clean/Graph_Embedding/Results_Long/bees_predicted_net.csv")
butterflies_GE <- read.csv("Data_Clean/Graph_Embedding/Results_Long/butterflies_predicted_net.csv")
hoverflies_GE <- read.csv("Data_Clean/Graph_Embedding/Results_Long/hoverflies_predicted_net.csv")
moths_GE <- read.csv("Data_Clean/Graph_Embedding/Results_Long/moths_predicted_net.csv")

# GE Predictions
GE_predictions <- bind_rows(bees_GE, butterflies_GE, hoverflies_GE, moths_GE)

# split into OOS and IS
GE_predictions_IS <- GE_predictions %>% 
  filter(type_of_learning == "In-sample learning",
         predicted_interaction == 1 & observed_interaction == 0) # only newly predicted ints
GE_predictions_OOS <- GE_predictions %>% 
  filter(type_of_learning == "Out-of-sample learning",
         predicted_interaction == 1)

###  Level 2: Graph Embedding - adding Out of Sample (OOS) predictions only
## OOS 
predictions_GEOOS_new <- GE_predictions_OOS %>%
  rename(sourceTaxonName_harm = genus_species) %>%
  left_join(plant_species, by = c("plant_genus" = "genus")) %>%
  # Remove rows where no matching plant species found
  filter(!is.na(genus_species)) %>%
  rename(targetTaxonName_harm = genus_species) %>%
  # remove pollinators without SDM or pheno data 
  filter(sourceTaxonName_harm %in% pollinator_species$genus_species) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm)  %>%
  mutate(method_type = "Level2_GEOOS")

predictions_GEOOS <- bind_rows(predictions_simple, predictions_GEOOS_new)

###  Level 3: Graph Embedding - adding In Sample (IS) plus OOS
## OOS 
predictions_GEOOSIS_new <- GE_predictions_IS %>%
  rename(sourceTaxonName_harm = genus_species) %>%
  left_join(plant_species, by = c("plant_genus" = "genus")) %>%
  # Remove rows where no matching plant species found
  filter(!is.na(genus_species)) %>%
  rename(targetTaxonName_harm = genus_species) %>%
  # remove pollinators without SDM or pheno data 
  filter(sourceTaxonName_harm %in% pollinator_species$genus_species) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm)  %>%
  mutate(method_type = "Level3_GEOOSIS")

predictions_all <- bind_rows(predictions_GEOOS, predictions_GEOOSIS_new)

#### Expand predictions_all to all regions where the species co-occur based on raw occurrences ####
sp_by_region <- read.csv("Data_Clean/Species_Checklists/ecoregion_checklist.csv")

# Get pollinator regions
pollinator_regions <- sp_by_region %>%
  filter(genus_species %in% pollinator_species$genus_species) %>%
  dplyr::select(genus_species, JEPCODE) %>%
  rename(sourceTaxonName_harm = genus_species)

# Get plant regions
plant_regions <- sp_by_region %>%
  filter(genus_species %in% pollinator_species$genus_species) %>%
  dplyr::select(genus_species, JEPCODE) %>%
  rename(targetTaxonName_harm = genus_species)

# Expand predictions to include JEPCODE where both species co-occur
# Get pollinator regions
pollinator_regions <- sp_by_region %>%
  dplyr::select(genus_species, JEPCODE) %>%
  rename(sourceTaxonName_harm = genus_species)

# Get plant regions  
plant_regions <- sp_by_region %>%
  dplyr::select(genus_species, JEPCODE) %>%
  rename(targetTaxonName_harm = genus_species)

# Expand predictions to include JEPCODE where both species co-occur
predictions_all_regions <- predictions_all %>%
  inner_join(pollinator_regions, by = "sourceTaxonName_harm") %>%
  inner_join(plant_regions, by = c("targetTaxonName_harm", "JEPCODE"))

### Remove spent datasets
rm(list = setdiff(ls(), c("predictions_all_regions", 
                          "checklist_prediction",
                          "phen_best",
                          "plant_species",
                          "pollinator_species")))

#### Load all SDMs at once ####
# Load all plant SDMs
plant_sdm_files <- list.files("Data_Clean/SDMs/sdm_by_taxon/plants/continuous", 
                              pattern = "\\.tif$", 
                              full.names = TRUE)
# Remove manifest if it exists
plant_sdm_files <- plant_sdm_files[!grepl("manifest", plant_sdm_files)]

# Load as a single raster stack
plant_sdms <- rast(plant_sdm_files)

# Get clean names (remove .tif and replace _ with space)
plant_names <- basename(plant_sdm_files) %>%
  gsub("\\_continuous.tif$", "", .) %>%
  gsub("_", " ", .)
names(plant_sdms) <- plant_names

# Load all pollinator SDMs at once
pollinator_folders <- c("bees", "butterflies", "moths", "hoverflies", "hummingbirds")

# Get all tif files from all pollinator folders
all_pollinator_files <- c()
for (taxon in pollinator_folders) {
  taxon_files <- list.files(
    file.path("Data_Clean/SDMs/sdm_by_taxon", taxon, "continuous"),
    pattern = "\\.tif$",
    full.names = TRUE
  )
  taxon_files <- taxon_files[!grepl("manifest", taxon_files)]
  all_pollinator_files <- c(all_pollinator_files, taxon_files)
}

# Load all as a single raster stack
pollinator_sdms <- rast(all_pollinator_files)

# Set clean names
pollinator_names <- basename(all_pollinator_files) %>%
  gsub("\\_continuous.tif$", "", .) %>%
  gsub("_", " ", .)
names(pollinator_sdms) <- pollinator_names

# Load jepson ecoregions
ecoregions <- read_sf("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid() %>%
  st_transform(crs(plant_sdms))  # Match SDM CRS

# Pre-calculate region vectors
jep_codes <- unique(ecoregions$JEPCODE)
region_vects <- list()
for (code in jep_codes) {
  region_geom <- st_geometry(ecoregions[ecoregions$JEPCODE == code,])
  if (length(region_geom) > 0) {
    region_vects[[code]] <- terra::vect(region_geom)
  }
}

#### Calculate spatial overlap ####
# Set number of cores
num_cores <- 25

# Convert to data.table
predictions_dt <- as.data.table(predictions_all_regions)
predictions_dt[, row_id := .I]

cat("Processing", nrow(predictions_dt), "rows using", num_cores, "cores\n")

# Split data into chunks for parallel processing
chunk_size <- ceiling(nrow(predictions_dt) / num_cores)
chunks <- split(predictions_dt$row_id, ceiling(seq_along(predictions_dt$row_id) / chunk_size))

# Function to process a chunk
process_chunk_spatialOverlap <- function(row_ids) {
  library(terra)
  library(data.table)
  
  # Create results vector
  results <- rep(NA_real_, length(row_ids))
  
  # Progress tracking variables
  total <- length(row_ids)
  last_message <- 0
  
  for (idx in seq_along(row_ids)) {
    i <- row_ids[idx]
    
    # Print progress every 10,000 iterations
    if (idx %% 10000 == 0 || idx == total) {
      message(sprintf("Processing: %d / %d (%.1f%%)", 
                      idx, total, (idx/total)*100))
      last_message <- idx
    }
    
    pollinator_sp <- predictions_dt$sourceTaxonName_harm[i]
    plant_sp <- predictions_dt$targetTaxonName_harm[i]
    jep_code <- predictions_dt$JEPCODE[i]
    
    # Check if both species have SDMs
    if (pollinator_sp %in% names(pollinator_sdms) && 
        plant_sp %in% names(plant_sdms)) {
      
      # Get individual SDMs
      pol_sdm <- pollinator_sdms[[pollinator_sp]]
      plant_sdm <- plant_sdms[[plant_sp]]
      
      # Calculate overlap
      overlap <- pol_sdm * plant_sdm
      
      # Get region mask
      if (jep_code %in% names(region_vects)) {
        region_vect <- region_vects[[jep_code]]
        
        # Create mask
        mask <- terra::rasterize(
          x = region_vect,
          y = overlap,
          field = 1,
          background = 0
        )
        
        # Calculate total overlap in region
        masked_overlap <- overlap * mask
        total_overlap <- terra::global(masked_overlap, fun = "sum", na.rm = TRUE)$sum
        
        results[idx] <- total_overlap
        
        # Clean up
        rm(masked_overlap, mask)
      }
      
      rm(overlap)
    }
  }
  
  return(data.table(row_id = row_ids, spatial_overlap = results))
}

# Run parallel processing
results_list <- mclapply(chunks, process_chunk_spatialOverlap, mc.cores = num_cores)

# checkpoint
saveRDS(results_list, "Temp/flowerVisit_int_temp.rds")
results_list <- readRDS("Temp/flowerVisit_int_temp.rds")

# Combine results
results_dt <- rbindlist(results_list)

# Merge back into original data
predictions_dt <- merge(predictions_dt, results_dt, by = "row_id", all.x = TRUE)
predictions_dt[, row_id := NULL]


#### Calculate phenological overlap ####
# Load phenology data
phen_best <- as.data.table(phen_best)
phen_best <- phen_best[, .(genus_species, JEPCODE, phen_start, phen_end, duration)]

# Function to calculate circular temporal overlap
calculate_circular_overlap <- function(x_start, x_end, y_start, y_end) {
  # Function to handle a single pair of sp
  calc_single <- function(xs, xe, ys, ye) {
    # Check if ranges wrap around year
    x_wraps <- xe < xs
    y_wraps <- ye < ys
    
    if (!x_wraps && !y_wraps) {
      # Neither wraps - simple case
      overlap_start <- pmax(xs, ys)
      overlap_end <- pmin(xe, ye)
      return(pmax(0, overlap_end - overlap_start))
      
    } else if (x_wraps && y_wraps) {
      # Both wrap around year
      # Convert to non-wrapping by splitting into two segments
      # x: [0, xe] and [xs, 365]
      # y: [0, ye] and [ys, 365]
      
      overlap1 <- pmin(xe, ye)  # Overlap in [0, min(xe, ye)]
      overlap2 <- pmax(0, 365 - pmax(xs, ys))  # Overlap in [max(xs, ys), 365]
      return(overlap1 + overlap2)
      
    } else if (x_wraps) {
      # Only x wraps: x is [0, xe] and [xs, 365]
      # Check overlap with y [ys, ye]
      
      overlap1 <- pmax(0, pmin(xe, ye) - pmax(0, ys))  # [0, xe] overlap with [ys, ye]
      overlap2 <- pmax(0, pmin(365, ye) - pmax(xs, ys))  # [xs, 365] overlap with [ys, ye]
      return(overlap1 + overlap2)
      
    } else {
      # Only y wraps: y is [0, ye] and [ys, 365]
      # Check overlap with x [xs, xe]
      
      overlap1 <- pmax(0, pmin(xe, ye) - pmax(xs, 0))  # [xs, xe] overlap with [0, ye]
      overlap2 <- pmax(0, pmin(xe, 365) - pmax(xs, ys))  # [xs, xe] overlap with [ys, 365]
      return(overlap1 + overlap2)
    }
  }
  
  # Vectorize the calculation
  mapply(calc_single, x_start, x_end, y_start, y_end)
}

# Prepare phenology data with proper column names
phen_pollinator <- phen_best[, .(genus_species, JEPCODE, phen_start, phen_end)]
setnames(phen_pollinator, 
         c("genus_species", "phen_start", "phen_end"),
         c("sourceTaxonName_harm", "pol_phen_start", "pol_phen_end"))

phen_plant <- phen_best[, .(genus_species, JEPCODE, phen_start, phen_end, duration)]
setnames(phen_plant,
         c("genus_species", "phen_start", "phen_end", "duration"),
         c("targetTaxonName_harm", "plant_phen_start", "plant_phen_end", "plant_duration"))

# Merge phenology data with predictions
predictions_dt <- merge(
  predictions_dt,
  phen_pollinator,
  by = c("sourceTaxonName_harm", "JEPCODE"),
  all.x = TRUE
)

predictions_dt <- merge(
  predictions_dt,
  phen_plant,
  by = c("targetTaxonName_harm", "JEPCODE"),
  all.x = TRUE
)

# Calculate temporal overlap
predictions_dt[, temporal_overlap := calculate_circular_overlap(
  pol_phen_start, pol_phen_end,
  plant_phen_start, plant_phen_end
)]

# Remove temporary phenology columns if desired
predictions_dt[, c("pol_phen_start", "pol_phen_end", 
                   "plant_phen_start", "plant_phen_end") := NULL]

# save!!
write.csv(predictions_dt, "Data_Clean/Interactions/predictions_supporting_visitsFlowersOf.csv", row.names=F)

