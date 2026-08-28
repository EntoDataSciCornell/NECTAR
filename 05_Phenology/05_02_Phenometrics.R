#### Load packages ####
library(lubridate)
library(parallel)
library(sf)
library(stringr)
library(dplyr)

#### Load Data ####
## read in checklists 
checklists <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

## Read in occurrence Data
occurrences_clean <- read.csv("Data_Clean/Occurrences/occurrences_clean.csv")

# Also read in gbif_plant_annotation data
gbif_plant_annotated <- read.csv("Data_Clean/Occurrences/gbif_plant_annotations_clean.csv") %>%
  dplyr::select(-X) %>%
  # phenologyCode is always F, but gets read in as False (logical), force back
  mutate(phenologyCode = "F")

## filter out non-flowering observations from plants
phenol_data <- occurrences_clean %>%
  filter(!(category == "all_plants")) %>%  #remove non-pheno plants
  filter(taxon != "plants" | (taxon == "plants" & phenologyCode == "F"))

## filter out all larval/juvenile life stages for pollinator taxa
# unique(phenol_data$lifeStage)
# Define larval/juvenile life stages to exclude
larval_stages <- c(
  "pupa", "Pupa",
  "Larvae", "Larval", "Larva",
  "Egg", "Egg cluster",
  "cocoon",
  "reared ova ex: 27 Feb 91", "ex ova on endive",
  "Adult w/ pupal case",
  "Larval Case", "Adult, Larval case"
)

# Filter out larval/juvenile stages for non-plant taxa
phenol_data <- phenol_data %>%
  filter(
    taxon == "plants" | 
      !(lifeStage %in% larval_stages) |
      is.na(lifeStage)
  )

# create DOY column using date information
phenol_data$DOY <- yday(make_date(phenol_data$yearClean, 
                                  phenol_data$monthClean, 
                                  phenol_data$dayClean))

# Merge occurrences with gbif_plant_annotated
phenol_data <- bind_rows(phenol_data, gbif_plant_annotated)

# Filter out NAs from DOY
phenol_data <- phenol_data %>%
  filter(!is.na(DOY))

## restrict records to only CA
# Read jepcode delineations shapefile
regions <- st_read("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid() %>%
  st_transform(4326)

# Filter points from phenol_data without jepson region
# Convert to sf object
phenol_data_sf <- st_as_sf(phenol_data, 
                           coords = c("decimalLongitude", "decimalLatitude"), 
                           crs = 4326,
                           remove = FALSE)  # Keep original lat/lon columns

# Spatial filter - keep only points that intersect with regions
phenol_data_sf <- st_filter(phenol_data_sf, regions, .predicate = st_intersects)

# Convert back to regular dataframe (keeping geometry dropped)
phenol_data <- st_drop_geometry(phenol_data_sf)

# remove spent dataset
rm(phenol_data_sf)

## split species into different type groups based on data availability
# summarise data
sp_phenol_sum <- phenol_data %>% 
  group_by(genus_species) %>%
  summarise(n = n(),
            n_unique = length(unique(DOY)))

# Type 1 - at least 10 unique observations on 3 or more unique days of the year
type1_sp <- sp_phenol_sum %>%
  filter(n >= 10,
         n_unique >= 3) %>%
  pull(genus_species)

# Type 2 - less than 10 unique observations on at least 2 or unique days of the year OR seen only on 2 days of the year
type2_sp <- sp_phenol_sum %>%
  filter((n < 10 & n_unique >= 2) | # V1 - not at least 10 records, but at least seen on 2 unique days
           (n >= 10 & n_unique == 2) ) %>% # V2 - more than 10 records, but only seen 2 unique days
  pull(genus_species)

# Type 3 - only seen on one unique day of the year
type3_sp <- sp_phenol_sum %>%
  filter(n_unique == 1) %>%
  pull(genus_species)

## Split the data frame into a list of subsets by species name
subsets_by_species <- split(phenol_data, phenol_data$genus_species)

# then break up by type
subsets_by_species_type1 <- subsets_by_species[names(subsets_by_species) %in% type1_sp]
subsets_by_species_type2 <- subsets_by_species[names(subsets_by_species) %in% type2_sp]
subsets_by_species_type3 <- subsets_by_species[names(subsets_by_species) %in% type3_sp]

rm(subsets_by_species)

#### Type 1: statewide-level Phenology ####
# source("Code/99_Supporting/phenol_functions.R")
source("Code/99_Supporting/phenol_skewedvonmises.R")

# wrapper function
phenol_wrapper_circular_skewed <- function(sp_data, 
                                           sp_name, 
                                           doy_col = "DOY", 
                                           quantiles = c(0.05, 0.95), 
                                           iterations = 100, 
                                           include_ci = FALSE,
                                           year_length = 365,
                                           n_terms = 1000,
                                           min_obs = 10,
                                           min_obs_unique = 3) {
  
  # Convert DOY to theta
  sp_theta <- doy_to_radians(sp_data[[doy_col]])
  
  # Return named NA structure if not enough phenol observations for the species
  if (length(sp_theta) < min_obs | length(unique(sp_theta)) < min_obs_unique) {
    result <- list(genus_species = sp_name)
    for (q in quantiles) {
      q_name <- paste0("q", sprintf("%02d", q * 100))
      result[[q_name]] <- NA
    }
    return(result)
  }
  
  # Nonparametric bootstrap quantile estimates using skewed von Mises
  if (include_ci == TRUE){
    result_list <- bootstrap_circular_quantiles_skewed(sp_theta, 
                                                       species_name = sp_name, 
                                                       quantiles = quantiles, 
                                                       n_boot = iterations, 
                                                       year_length = year_length,
                                                       include_ci = include_ci) 
  } else{
    ## If no CIs requested, no need to do bootstrapping. Just get maximum likelihood parameter values from skewed VM
    # Fit sine skewed von Mises distribution
    svm_params <- fit_skewed_von_mises(sp_theta)
    
    # if returned NULL, end here
    if (is.null(svm_params)) {
      result <- list(genus_species = sp_name)
      for (q in quantiles) {
        q_name <- paste0("q", sprintf("%02d", q * 100))
        result[[q_name]] <- NA
      }
      return(result)
    }
    
    # Calculate quantiles from the fitted distribution
    quants <- qskewed_vonmises(quantiles, 
                               mu = svm_params$mu, 
                               kappa = svm_params$kappa, 
                               lambda = svm_params$lambda,
                               n_terms = n_terms)
    
    # Create named list of results
    result_list <- list()
    
    # Add species name if provided
    if (!is.null(sp_name)) {
      result_list$genus_species <- sp_name
    }
    
    # Add quantile estimates (no CIs)
    for (i in 1:length(quants)) {
      q_name <- paste0("q", sprintf("%02d", quantiles[i] * 100))
      
      # Point estimate - convert from radians to DOY
      result_list[[q_name]] <- radians_to_doy(quants[i], year_length)
    }
  }
  
  return(result_list)
}

# Apply with names preserved
phen_circular <- mcMap(phenol_wrapper_circular_skewed,
                       subsets_by_species_type1,
                       names(subsets_by_species_type1),
                       MoreArgs = list(quantiles = c(0.05, 0.95), 
                                       include_ci = FALSE,
                                       min_obs = 10,
                                       min_obs_unique = 3),
                       mc.cores = 10)

# Combine results
phen_df <- bind_rows(phen_circular) %>%
  rename(phen_start = q05, phen_end = q95) %>%
  mutate(duration = ifelse(phen_start>=phen_end,
                           (365-phen_start + phen_end),
                           phen_end - phen_start)) # length of adult active period

# Save as csv
write.csv(phen_df, "Data_Clean/Phenology/phenology_simple.csv")

#### Type 1: Ecoregion-level Phenology ####
# Regional phenology wrapper function
phenol_wrapper_regional <- function(sp_data, 
                                    sp_name, 
                                    region_cols,
                                    regions_sf,
                                    doy_col = "DOY", 
                                    quantiles = c(0.05, 0.95), 
                                    iterations = 100, 
                                    include_ci = FALSE,
                                    year_length = 365,
                                    n_terms = 1000,
                                    min_obs = 10,
                                    min_obs_unique = 3) {
  
  # Convert to sf object (only this species)
  sp_sf <- st_as_sf(sp_data, 
                    coords = c("decimalLongitude", "decimalLatitude"), 
                    crs = 4326)
  
  # Spatial join to assign regions
  sp_sf <- st_join(sp_sf, regions_sf)
  
  # Drop geometry for faster processing
  sp_df <- st_drop_geometry(sp_sf)
  
  # Loop over each region column
  all_results <- map_dfr(region_cols, function(region_col) {
    
    # Check if column exists
    if (!region_col %in% names(sp_df)) {
      warning(paste("Column", region_col, "not found in data. Skipping."))
      return(NULL)
    }
    
    # Group by region and calculate phenology
    results <- sp_df %>%
      group_by(across(all_of(region_col))) %>%
      group_modify(~{
        n_obs <- nrow(.x)
        
        # Check if we have enough observations
        if (nrow(.x) < min_obs | length(unique(.x[[doy_col]])) < min_obs_unique) {
          # Create result with NA values
          result <- tibble(
            genus_species = sp_name,
            n_observations_used = n_obs
          )          
          for (q in quantiles) {
            q_name <- paste0("q", sprintf("%02d", q * 100))
            result[[q_name]] <- NA
          }
          return(result)
        }
        
        # Convert DOY to theta
        sp_theta <- doy_to_radians(.x[[doy_col]])
        
        # Calculate phenology based on include_ci flag
        if (include_ci == TRUE) {
          # Full bootstrap with CIs
          boot_quantiles <- bootstrap_circular_quantiles_skewed(
            sp_theta, 
            species_name = sp_name, 
            quantiles = quantiles, 
            n_boot = iterations, 
            year_length = year_length,
            include_ci = include_ci
          )
          result <- as_tibble(boot_quantiles)
        } else {
          # Just get maximum likelihood parameter values from skewed VM
          svm_params <- fit_skewed_von_mises(sp_theta)
          
          # stop here if NA values for distribution fit
          if (is.null(svm_params)) {
            result <- tibble(
              genus_species = sp_name,
              n_observations_used = n_obs
            )
            for (q in quantiles) {
              q_name <- paste0("q", sprintf("%02d", q * 100))
              result[[q_name]] <- NA
            }
            return(result)
          }
          
          # Calculate quantiles from the fitted distribution
          quants <- qskewed_vonmises(
            quantiles, 
            mu = svm_params$mu, 
            kappa = svm_params$kappa, 
            lambda = svm_params$lambda,
            n_terms = n_terms
          )
          
          # Create result tibble
          result <- tibble(
            genus_species = sp_name,
            n_observations_used = n_obs
          )
          
          # Add quantile estimates
          for (i in seq_along(quants)) {
            q_name <- paste0("q", sprintf("%02d", quantiles[i] * 100))
            result[[q_name]] <- radians_to_doy(quants[i], year_length)
          }
        }
        
        return(result)
      }) %>%
      # Add a column indicating which regional grouping this is
      mutate(region_type = region_col, .before = 1) %>%
      # Rename the grouping column
      rename(region_name = all_of(region_col)) %>% 
      ungroup()
    
    return(results)
  })
  
  # Clean up
  rm(sp_sf, sp_df)
  
  return(all_results)
}

phen_region_results <- mcMap(
  phenol_wrapper_regional,
  subsets_by_species_type1,
  names(subsets_by_species_type1),
  MoreArgs = list(
    region_cols = c("PROV", "REGION", "JEPCODE"),
    regions_sf = regions,
    quantiles = c(0.05, 0.95),
    iterations = 100,
    include_ci = FALSE,
    min_obs = 10,
    min_obs_unique = 3
  ),
  mc.cores = 10
)

# Combine results
phen_region_results_df <- bind_rows(phen_region_results) %>%
  rename(phen_start = q05, phen_end = q95) %>%
  mutate(duration = ifelse(phen_start >= phen_end,
                           (365 - phen_start + phen_end),
                           phen_end - phen_start))

# Break into each region type
phen_prov <- phen_region_results_df %>%
  filter(region_type == "PROV") %>%
  dplyr::select(-region_type) %>%
  rename(PROV = region_name)
phen_region <- phen_region_results_df %>%
  filter(region_type == "REGION") %>%
  dplyr::select(-region_type) %>%
  rename(REGION = region_name)
phen_jepcode <- phen_region_results_df %>%
  filter(region_type == "JEPCODE") %>%
  dplyr::select(-region_type) %>%
  rename(JEPCODE = region_name)

# Save output files
write.csv(phen_prov, "Data_Clean/Phenology/phenology_eco_prov.csv", row.names = FALSE)
write.csv(phen_region, "Data_Clean/Phenology/phenology_eco_region.csv", row.names = FALSE)
write.csv(phen_jepcode, "Data_Clean/Phenology/phenology_eco_jepcode.csv", row.names = FALSE)

#### Type 2: between min and max ####
# Requirements: less than 10 data points and at least 2 unique days OR 2 unique days
pheno_interval <- function(sp_data,
                           sp_name, 
                           doy_col = "DOY") {
  
  sp_theta <- unique(doy_to_radians(sp_data[[doy_col]]))
  
  # Sort thetas
  theta_sorted <- sort(sp_theta)
  n <- length(theta_sorted)
  
  # Calculate gaps between consecutive points (including wrap-around for change of year)
  gaps <- numeric(n)
  for (i in 1:(n-1)) {
    gaps[i] <- theta_sorted[i+1] - theta_sorted[i]
  }
  # turn of year gap
  gaps[n] <- (2 * pi - theta_sorted[n]) + theta_sorted[1]
  
  # Duration spans from the point after the largest gap to the point before it
  max_gap_idx <- which.max(gaps)
  
  # Start of interval is the point after the largest gap
  start_idx <- (max_gap_idx %% n) + 1
  # End of interval is the point at the largest gap
  end_idx <- max_gap_idx
  
  start_theta <- theta_sorted[start_idx]
  end_theta <- theta_sorted[end_idx]
  
  # Calculate interval length in radians
  if (start_idx <= end_idx) {
    duration_radians <- end_theta - start_theta
  } else {
    # Wraps around the year
    duration_radians <- (2 * pi - start_theta) + end_theta
  }
  
  # Convert back to days
  duration_days <- radians_to_doy(duration_radians)
  
  # Convert start and end back to DOY
  start_doy <- radians_to_doy(start_theta)
  end_doy <- radians_to_doy(end_theta)
  
  list(
    genus_species = sp_name,
    phen_start = start_doy,
    phen_end = end_doy,
    duration = round(duration_days)
  )
}

phen_type2_results <- mcMap(
  pheno_interval,
  subsets_by_species_type2,
  names(subsets_by_species_type2),
  MoreArgs = list(
    doy_col = "DOY"
  ),
  mc.cores = 10
)

# Combine results
phen_type2_results_df <- bind_rows(phen_type2_results)

write.csv(phen_type2_results_df, "Data_Clean/Phenology/phenology_type2.csv", row.names = FALSE)

#### Type 3: buffer around only day ####
# Requirements: only seen on one unique day of the year
pheno_buffer <- function(sp_data,
                         sp_name, 
                         doy_col = "DOY",
                         buffer_days = 7) {
  
  # Check if only one DOY
  if (length(unique(sp_data[[doy_col]])) != 1) {
    warning(paste("More than one unique point, not valid for this method"))
    return(NULL)
  }
  
  # Get unique day
  day_obs <- unique(sp_data[[doy_col]])
  
  # Calculate start and end with wrapping
  phen_start_raw <- day_obs - buffer_days
  phen_end_raw <- day_obs + buffer_days
  
  # Wrap to [1, 365] range
  phen_start <- ((phen_start_raw - 1) %% 365) + 1
  phen_end <- ((phen_end_raw - 1) %% 365) + 1
  
  list(
    genus_species = sp_name,
    phen_start = phen_start,
    phen_end = phen_end,
    duration = 2 * buffer_days + 1 # always 15 due to 7 day buffer
  )
}

phen_type3_results <- mcMap(
  pheno_buffer,
  subsets_by_species_type3,
  names(subsets_by_species_type3),
  MoreArgs = list(
    doy_col = "DOY"
  ),
  mc.cores = 10
)

# Combine results
phen_type3_results_df <- bind_rows(phen_type3_results)

write.csv(phen_type3_results_df, "Data_Clean/Phenology/phenology_type3.csv", row.names = FALSE)

#### Combine and iteratively choose highest resolution estimates ####
# MERGING OLD AND NEW
phen_l3 <- read.csv("Data_Clean/Phenology/phenology_eco_jepcode.csv")
phen_l2 <- read.csv("Data_Clean/Phenology/phenology_eco_region.csv")
phen_l1 <- read.csv("Data_Clean/Phenology/phenology_eco_prov.csv")
phen_simple <- read.csv("Data_Clean/Phenology/phenology_simple.csv")
phen_type2 <- read.csv("Data_Clean/Phenology/phenology_type2.csv")
phen_type3 <- read.csv("Data_Clean/Phenology/phenology_type3.csv")

## Create "best available" pheno dataset
# read in ecoregion checklist
sp_by_region <- read.csv("Data_Clean/Species_Checklists/ecoregion_checklist.csv") %>%
  # only species level, remove subspecies
  dplyr::select(-scientificName) %>%
  unique()

# Species and region
phen_best <- sp_by_region %>%
  # get "best" phenology information
  mutate(
    phen_start = coalesce(
      phen_l3$phen_start[match(paste(genus_species, JEPCODE), 
                               paste(phen_l3$genus_species, phen_l3$JEPCODE))],  # Check ecoregion 3 level
      phen_l2$phen_start[match(paste(genus_species, REGION), 
                               paste(phen_l2$genus_species, phen_l2$REGION))],  # Check ecoregion 2 level
      phen_l1$phen_start[match(paste(genus_species, PROV), 
                               paste(phen_l1$genus_species, phen_l1$PROV))],  # Check ecoregion 1 level
      phen_simple$phen_start[match(genus_species, phen_simple$genus_species)],  # Check species level
      phen_type2$phen_start[match(genus_species, phen_type2$genus_species)],  # Type 2
      phen_type3$phen_start[match(genus_species, phen_type3$genus_species)]  # Type 3
    ),
    
    phen_end = coalesce(
      phen_l3$phen_end[match(paste(genus_species, JEPCODE), 
                             paste(phen_l3$genus_species, phen_l3$JEPCODE))],  # Check ecoregion 3 level
      phen_l2$phen_end[match(paste(genus_species, REGION), 
                             paste(phen_l2$genus_species, phen_l2$REGION))],  # Check ecoregion 2 level
      phen_l1$phen_end[match(paste(genus_species, PROV), 
                             paste(phen_l1$genus_species, phen_l1$PROV))],  # Check ecoregion 1 level
      phen_simple$phen_end[match(genus_species, phen_simple$genus_species)], # Check species level
      phen_type2$phen_end[match(genus_species, phen_type2$genus_species)],  # Type 2
      phen_type3$phen_end[match(genus_species, phen_type3$genus_species)]  # Type 3
      
    ),
    
    duration = coalesce(
      phen_l3$duration[match(paste(genus_species, JEPCODE), 
                             paste(phen_l3$genus_species, phen_l3$JEPCODE))],  # Check ecoregion 3 level
      phen_l2$duration[match(paste(genus_species, REGION), 
                             paste(phen_l2$genus_species, phen_l2$REGION))],  # Check ecoregion 2 level
      phen_l1$duration[match(paste(genus_species, PROV), 
                             paste(phen_l1$genus_species, phen_l1$PROV))],  # Check ecoregion 1 level
      phen_simple$duration[match(genus_species, phen_simple$genus_species)],  # Check species level
      phen_type2$duration[match(genus_species, phen_type2$genus_species)],  # Type 2
      phen_type3$duration[match(genus_species, phen_type3$genus_species)]  # Type 3
    ),
    
    # What scale are the phenology estimates from??
    phen_scale = case_when( # keep track of which scale the estimate came from
      !is.na(phen_l3$phen_start[match(paste(genus_species, JEPCODE), 
                                      paste(phen_l3$genus_species, phen_l3$JEPCODE))]) ~ "JEPCODE",
      !is.na(phen_l2$phen_start[match(paste(genus_species, REGION), 
                                      paste(phen_l2$genus_species, phen_l2$REGION))]) ~ "REGION",
      !is.na(phen_l1$phen_start[match(paste(genus_species, PROV), 
                                      paste(phen_l1$genus_species, phen_l1$PROV))]) ~ "PROV",
      !is.na(phen_simple$phen_start[match(genus_species, phen_simple$genus_species)]) ~ "CALIFORNIA",
      !is.na(phen_type2$phen_start[match(genus_species, phen_type2$genus_species)]) ~ "TYPE2",
      !is.na(phen_type3$phen_start[match(genus_species, phen_type3$genus_species)]) ~ "TYPE3",
      TRUE ~ NA_character_
    ),
    
    # Number of observations at each scale
    n_observations_JEPCODE = phen_l3$n_observations_used[match(paste(genus_species, JEPCODE),
                                                               paste(phen_l3$genus_species, phen_l3$JEPCODE))],
    n_observations_REGION = phen_l2$n_observations_used[match(paste(genus_species, REGION),
                                                              paste(phen_l2$genus_species, phen_l2$REGION))],
    n_observations_PROV = phen_l1$n_observations_used[match(paste(genus_species, PROV),
                                                            paste(phen_l1$genus_species, phen_l1$PROV))],
    n_observations_CA = sp_phenol_sum$n[match(genus_species, sp_phenol_sum$genus_species)],
    
    # Number of observations used for pheno estimates
    n_observations_used = case_when( # keep track of which scale the estimate came from
      (phen_scale == "JEPCODE") ~ n_observations_JEPCODE,
      (phen_scale == "REGION") ~ n_observations_REGION,
      (phen_scale == "PROV") ~ n_observations_PROV,
      (phen_scale %in% c("CALIFORNIA", "TYPE2", "TYPE3")) ~ n_observations_CA,
      TRUE ~ NA_real_
    )
  )

# Save as csv
write.csv(phen_best, "Data_Clean/Phenology/phenology_best_available.csv")

##### Remove Gymnosperms that made it through to here #####
phen_best <- read.csv("Data_Clean/Phenology/phenology_best_available.csv")
phen_simple <- read.csv("Data_Clean/Phenology/phenology_simple.csv")
phen_type2 <- read.csv("Data_Clean/Phenology/phenology_type2.csv")
phen_type3 <- read.csv("Data_Clean/Phenology/phenology_type3.csv")
phen_eco_jepcode <- read.csv("Data_Clean/Phenology/phenology_eco_jepcode.csv")
phen_eco_prov <- read.csv("Data_Clean/Phenology/phenology_eco_prov.csv")
phen_eco_region <- read.csv("Data_Clean/Phenology/phenology_eco_region.csv")

# checklists
checklists <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")


### Remove Gymnosperm flowering phenology
gymnosperm_orders <- c(
  "Ginkgoales", "Cycadales", "Araucariales",
  "Cupressales", "Pinales", "Ephedrales",
  "Welwitschiales", "Gnetales"
)

gymnosperma_taxa <- checklists %>%
  filter(order %in% gymnosperm_orders)

# Filter out gymnosperms
phen_best <- phen_best %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_simple <- phen_simple %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_type2 <- phen_type2 %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_type3 <- phen_type3 %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_eco_jepcode <- phen_eco_jepcode %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_eco_prov <- phen_eco_prov %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))
phen_eco_region <- phen_eco_region %>%
  filter(!(genus_species %in% gymnosperma_taxa$scientificName))


# Save files again
write.csv(phen_best, "Data_Clean/Phenology/phenology_best_available.csv", row.names=F)
write.csv(phen_simple, "Data_Clean/Phenology/phenology_simple.csv", row.names=F)
write.csv(phen_type2, "Data_Clean/Phenology/phenology_type2.csv", row.names=F)
write.csv(phen_type3, "Data_Clean/Phenology/phenology_type3.csv", row.names=F)
write.csv(phen_eco_jepcode, "Data_Clean/Phenology/phenology_eco_jepcode.csv", row.names=F)
write.csv(phen_eco_prov, "Data_Clean/Phenology/phenology_eco_prov.csv", row.names=F)
write.csv(phen_eco_region, "Data_Clean/Phenology/phenology_eco_region.csv", row.names=F)
