#### Load packages ####
library(brms)
library(future)
library(cmdstanr)
library(tidybayes)
library(furrr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(sf)
library(ggplot2)
library(patchwork)
library(nngeo)
library(bdc)
library(purrr)

#### Supporting scripts ####
source("Code/12_Analyses/12_00_LoadData.R")
source("Code/99_Supporting/name_harmonization.R")
source("Code/99_Supporting/genetic_algorithm.R")

gnverifier_path <- "/home/tb625/bin"

set.seed(123)

#### PART 1: Data Prep and plant selections ####
#### Load additional datasets ####
## Pollinator partnership lists
polpart_eco_plants <- read.csv("Data_Raw/External_NativePlantLists/P2/pollinatorPartnership_ecoregionPlantLists.csv") %>%
  mutate(list = "P2") %>%
  dplyr::select(Ecoregion, Botanical.Name, list) %>%
  rename(region = Ecoregion,
         scientificName = Botanical.Name)

# Bailey ecoregion shapefile
ecoregions_bailey <- read_sf("Data_Raw/Spatial/Bailey_Ecoregions/eco_us.shp") %>%
  group_by(PROVINCE) %>%
  summarise() %>%
  st_transform(st_crs(ecoregions)) %>%
  rename(P2_region = PROVINCE)

## Xerces
xerces_eco_plants <- read_excel("Data_Raw/External_NativePlantLists/Xerces/xerces_plantList.xlsx") %>%
  mutate(list = "Xerces")
xerces_regions <- read_sf("Data_Raw/External_NativePlantLists/Xerces/xerces_georeferenced_regions/xerces_regions.shp") %>%
  rename(xerces_region = region)

## Join together
external_lists <- bind_rows(polpart_eco_plants, xerces_eco_plants)

#### Harmonize names from external lists ####
external_list_names <- unique(external_lists$scientificName)
external_list_names_clean <- bdc_clean_names(sci_names = external_list_names, save_outputs = FALSE)

names_harm_external_plants <- harmonize_names(
  names = unique(external_list_names_clean$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_external_plants.tsv",
  names_harm_file = "Temp/names_harmonized_external_plants.tsv",
  gnv_path = gnverifier_path
)

names_harm_external_plants <- apply_jepson_taxonomy(names_harm_external_plants)
# Note: some plants not in CA will have names_harm == NA

# Join back to external lists
external_lists_cleaned <- external_lists %>%
  left_join(names_harm_external_plants, by = join_by("scientificName" == "ScientificName")) %>%
  dplyr::select(c(names(external_lists), names_harm)) %>%
  filter(!is.na(names_harm))

#### Combine geometries for each region ####
ecoregions_simp <- ecoregions %>%
  group_by(JEPCODE) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()
ecoregions_bailey_simp <- ecoregions_bailey %>%
  group_by(P2_region) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()
ecoregions_xerces_simp <- xerces_regions %>%
  st_make_valid() %>%
  group_by(xerces_region) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()

# Pairwise intersections of geometries
eco_12 <- st_intersection(
  ecoregions_simp %>% dplyr::select(JEPCODE),
  ecoregions_bailey_simp %>% dplyr::select(P2_region)
) %>% st_make_valid()

eco_123 <- st_intersection(
  eco_12,
  ecoregions_xerces_simp %>% dplyr::select(xerces_region)
) %>% st_make_valid()

eco_12_only <- st_difference(
  eco_12,
  st_union(ecoregions_xerces_simp)
) %>% st_make_valid() %>%
  mutate(xerces_region = NA_character_)

eco_13 <- st_intersection(
  ecoregions_simp %>% dplyr::select(JEPCODE),
  ecoregions_xerces_simp %>% dplyr::select(xerces_region)
) %>% st_make_valid()

eco_13_only <- st_difference(
  eco_13,
  st_union(ecoregions_bailey_simp)
) %>% st_make_valid() %>%
  mutate(P2_region = NA_character_)

eco_23 <- st_intersection(
  ecoregions_bailey_simp %>% dplyr::select(P2_region),
  ecoregions_xerces_simp %>% dplyr::select(xerces_region)
) %>% st_make_valid()

eco_23_only <- st_difference(
  eco_23,
  st_union(ecoregions_simp)
) %>% st_make_valid() %>%
  mutate(JEPCODE = NA_character_)

# Combine all where there is a match with Jepson region
eco_all_combos <- bind_rows(
  eco_123,
  eco_12_only,
  eco_13_only
) %>%
  st_make_valid()

#### List plants from external lists for each unique region combo ####
external_plants_P2 <- eco_all_combos %>%
  st_drop_geometry() %>%
  unique() %>%
  left_join(external_lists_cleaned[external_lists_cleaned$list == "P2", ], by = join_by("P2_region" == "region")) %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, names_harm)

external_plants_xerces <- eco_all_combos %>%
  st_drop_geometry() %>%
  unique() %>%
  left_join(external_lists_cleaned[external_lists_cleaned$list == "Xerces", ], by = join_by("xerces_region" == "region")) %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, names_harm)

external_plants <- bind_rows(external_plants_P2, external_plants_xerces) %>%
  unique()

external_plants_number <- external_plants %>%
  group_by(JEPCODE, P2_region, xerces_region) %>%
  summarise(n_plants = n()) %>%
  # filter regions with less than 10 plants from external lists
  filter(n_plants >= 10)

#### Some summary plots ####
eco_all_combos_type <- eco_all_combos %>%
  mutate(combo_type = case_when(
    !is.na(JEPCODE) & !is.na(P2_region) & !is.na(xerces_region) ~ "JEPSON + P2 + Xerves",
    !is.na(JEPCODE) & !is.na(P2_region) ~ "JEPSON + P2",
    !is.na(JEPCODE) & !is.na(xerces_region) ~ "JEPSON + Xerces",
    TRUE ~ "Single"
  ))

number_overlapping_lists <- ggplot(eco_all_combos_type) +
  geom_sf(aes(fill = combo_type)) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(fill = "Overlap Type")

ggsave("Figures/Supplementary/number_overlapping_lists.pdf", number_overlapping_lists, width = 10, height = 10)

eco_all_combos_numPlants <- eco_all_combos %>%
  left_join(external_plants_number)

numberPlants <- ggplot(eco_all_combos_numPlants) +
  geom_sf(aes(fill = n_plants)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Number Native Plants\nin External Lists")

ggsave("Figures/Supplementary/numberPlants.pdf", numberPlants, width = 10, height = 10)

#### Interaction data density per jepcode ####
ecoregions_density <- ecoregions_simp %>%
  mutate(
    n_points = lengths(st_intersects(., native_CA_interactions)),
    area_km2 = as.numeric(st_area(geometry)) / 1e6,
    density = n_points / area_km2
  )

ecoregions_density <- ecoregions_density %>%
  mutate(int_density_quartile = case_when(
    ntile(density, 4) == 1 ~ "low",
    ntile(density, 4) == 2 ~ "mid-low",
    ntile(density, 4) == 3 ~ "mid-high",
    ntile(density, 4) == 4 ~ "high"
  ))

#### Top plants per region, based on n from lists ####
predicted_interactions_filtered <- predicted_interactions %>%
  filter(interactionTypeName == "visitsFlowersOf",
         method_type %in% c("Level1_Simple", "Level2_GEOOS"))

pols_supported <- predicted_interactions_filtered %>%
  dplyr::select(JEPCODE, sourceTaxonType, sourceTaxonName_harm, targetTaxonName_harm) %>%
  group_by(targetTaxonName_harm, JEPCODE, sourceTaxonType) %>%
  summarise(n = n())

pols_supported_wide <- pols_supported %>%
  group_by(targetTaxonName_harm, JEPCODE) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sourceTaxonType,
    values_from = n,
    values_fill = 0
  )

top_plants_predicted <- pols_supported_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(total), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_predicted_bees <- pols_supported_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(bees), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_predicted_butterflies <- pols_supported_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(butterflies), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_predicted_hoverflies <- pols_supported_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(hoverflies), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_predicted_moths <- pols_supported_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(moths), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

#### Organize top plants by group ####
top_plants_by_group <- list(
  all_pollinators = top_plants_predicted,
  bees            = top_plants_predicted_bees,
  butterflies     = top_plants_predicted_butterflies,
  hoverflies      = top_plants_predicted_hoverflies,
  moths           = top_plants_predicted_moths
)

# Top 25 plants per region as supplementary table
top_25plants_predicted <- pols_supported_wide %>%
  group_by(JEPCODE) %>%
  arrange(desc(total), .by_group = TRUE) %>%
  filter(row_number() <= 25) %>%
  mutate(region_rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, region_rank, targetTaxonName_harm, total, bees, butterflies, hoverflies, moths) %>%
  rename(Jepson.Region.JEPCODE = JEPCODE,
         Region.Rank = region_rank,
         Plant.Name = targetTaxonName_harm,
         Pollinators.Supported = total,
         Bees.Supported = bees,
         Butterflies.Supported = butterflies,
         Hoverflies.Supported = hoverflies,
         Moths.Supported = moths)

write.csv(top_25plants_predicted, "Tables/Regional_HighDegreePlants.csv")

top_plants_external <- external_plants %>%
  rename(targetTaxonName_harm = names_harm)

plants_per_region <- sp_by_region %>%
  filter(genus_species %in% checklist_cleaned$genus_species[checklist_cleaned$taxon == "plants"]) %>%
  dplyr::select(genus_species, JEPCODE)

random_plants_prepList <- eco_all_combos_type %>%
  st_drop_geometry() %>%
  left_join(plants_per_region) %>%
  rename(targetTaxonName_harm = genus_species)

#### Raw Metaweb based plant optimization ####
interactions_nectar_observed <- interactions_nectar %>%
  distinct(targetTaxonName_harm, sourceTaxonName_harm) %>%
  filter(targetTaxonName_harm %in% checklist_cleaned$genus_species,
         sourceTaxonName_harm %in% checklist_cleaned$genus_species) %>%
  left_join(unique(checklist_cleaned[, c("genus_species", "taxon")]),
            by = join_by("sourceTaxonName_harm" == "genus_species")) %>%
  rename(sourceTaxonType = taxon)

raw_metaweb <- interactions_nectar_observed %>%
  merge(data.frame(JEPCODE = unique(ecoregions_simp$JEPCODE))) %>%
  semi_join(sp_by_region, by = c("targetTaxonName_harm" = "genus_species", "JEPCODE")) %>%
  semi_join(sp_by_region, by = c("sourceTaxonName_harm" = "genus_species", "JEPCODE"))

pols_supported_raw_metaweb <- raw_metaweb %>%
  dplyr::select(JEPCODE, sourceTaxonType, sourceTaxonName_harm, targetTaxonName_harm) %>%
  group_by(targetTaxonName_harm, JEPCODE, sourceTaxonType) %>%
  summarise(n = n())

pols_supported_raw_metaweb_wide <- pols_supported_raw_metaweb %>%
  group_by(targetTaxonName_harm, JEPCODE) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sourceTaxonType,
    values_from = n,
    values_fill = 0
  )

top_plants_raw_metaweb <- pols_supported_raw_metaweb_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(total), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_raw_metaweb_bees <- pols_supported_raw_metaweb_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(bees), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_raw_metaweb_butterflies <- pols_supported_raw_metaweb_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(butterflies), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_raw_metaweb_hoverflies <- pols_supported_raw_metaweb_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(hoverflies), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

top_plants_raw_metaweb_moths <- pols_supported_raw_metaweb_wide %>%
  inner_join(external_plants_number, by = "JEPCODE", relationship = "many-to-many") %>%
  group_by(JEPCODE, P2_region, xerces_region, n_plants) %>%
  arrange(desc(moths), .by_group = TRUE) %>%
  filter(row_number() <= n_plants) %>%
  ungroup() %>%
  dplyr::select(JEPCODE, P2_region, xerces_region, targetTaxonName_harm)

#### Organize top plants by group, raw metaweb based selection ####
top_plants_by_group_raw_metaweb <- list(
  all_pollinators = top_plants_raw_metaweb,
  bees            = top_plants_raw_metaweb_bees,
  butterflies     = top_plants_raw_metaweb_butterflies,
  hoverflies      = top_plants_raw_metaweb_hoverflies,
  moths           = top_plants_raw_metaweb_moths
)


#### Total pollinators per region ####
tax <- checklist_cleaned %>%
  dplyr::select(genus_species, taxon) %>%
  unique()

total_pollinators_by_region <- sp_by_region %>%
  dplyr::select(genus_species, JEPCODE) %>%
  left_join(tax) %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies")) %>%
  group_by(JEPCODE) %>%
  summarise(
    total_pollinators  = n_distinct(genus_species),
    total_bees         = n_distinct(genus_species[taxon == "bees"]),
    total_hoverflies   = n_distinct(genus_species[taxon == "hoverflies"]),
    total_moths        = n_distinct(genus_species[taxon == "moths"]),
    total_butterflies  = n_distinct(genus_species[taxon == "butterflies"])
  )

#### Function to sample plants and count unique pollinators ####
sample_and_count_unique_pollinators <- function(plant_list,
                                                interactions_df,
                                                dataset_name,
                                                totals_df,
                                                iter = 1,
                                                use_spatial = TRUE,
                                                sp_by_region_df = NULL) {
  sampled_plants <- plant_list %>%
    group_by(JEPCODE, P2_region, xerces_region) %>%
    slice_sample(n = 10, replace = FALSE) %>%
    ungroup()

  if (use_spatial) {
    plants_with_interactions <- sampled_plants %>%
      left_join(interactions_df,
                by = c("targetTaxonName_harm", "JEPCODE"),
                relationship = "many-to-many")
  } else {
    plants_with_interactions <- sampled_plants %>%
      left_join(interactions_df,
                by = "targetTaxonName_harm",
                relationship = "many-to-many")
    if (!is.null(sp_by_region_df)) {
      region_pols <- sp_by_region_df %>%
        dplyr::select(genus_species, JEPCODE) %>%
        distinct()
      plants_with_interactions <- plants_with_interactions %>%
        inner_join(region_pols,
                   by = c("JEPCODE", "sourceTaxonName_harm" = "genus_species"))
    }
  }

  pollinator_counts <- plants_with_interactions %>%
    group_by(JEPCODE, P2_region, xerces_region) %>%
    summarise(
      iteration          = iter,
      dataset            = dataset_name,
      n_unique_pollinators  = n_distinct(sourceTaxonName_harm, na.rm = TRUE),
      n_unique_bees         = n_distinct(sourceTaxonName_harm[sourceTaxonType == "bees"], na.rm = TRUE),
      n_unique_hoverflies   = n_distinct(sourceTaxonName_harm[sourceTaxonType == "hoverflies"], na.rm = TRUE),
      n_unique_moths        = n_distinct(sourceTaxonName_harm[sourceTaxonType == "moths"], na.rm = TRUE),
      n_unique_butterflies  = n_distinct(sourceTaxonName_harm[sourceTaxonType == "butterflies"], na.rm = TRUE),
      .groups = "drop"
    )

  pollinator_counts %>%
    left_join(totals_df, by = "JEPCODE") %>%
    mutate(
      pct_pollinators  = (n_unique_pollinators / total_pollinators) * 100,
      pct_bees         = (n_unique_bees / total_bees) * 100,
      pct_hoverflies   = (n_unique_hoverflies / total_hoverflies) * 100,
      pct_moths        = (n_unique_moths / total_moths) * 100,
      pct_butterflies  = (n_unique_butterflies / total_butterflies) * 100
    )
}

#### Calculate complementary plants for each region (JEPCODE) ####
message("Finding complementary plant combinations...")

unique_jepcodes <- ecoregions_simp %>%
  st_drop_geometry() %>%
  pull(JEPCODE) %>%
  unique()

complementary_results <- map_dfr(unique_jepcodes, function(jepcode) {
  message(paste("Processing region:", jepcode))
  results_list <- list()

  for (group in c("all", "bees", "butterflies", "hoverflies", "moths")) {
    message(paste("  Finding complementary plants for:", group))

    region_data <- data.frame(JEPCODE = jepcode)

    best_plants <- find_complementary_plants(
      region_data     = region_data,
      interactions_df = predicted_interactions_filtered,
      totals_df       = total_pollinators_by_region,
      target_group    = group,
      n_plants        = 10,
      n_generations   = 100,
      population_size = 200,
      mutation_rate   = 0.1
    )

    plants_df <- data.frame(JEPCODE = jepcode, targetTaxonName_harm = best_plants)

    plants_with_interactions <- plants_df %>%
      left_join(predicted_interactions_filtered,
                by = c("targetTaxonName_harm", "JEPCODE"),
                relationship = "many-to-many")

    result <- plants_with_interactions %>%
      summarise(
        JEPCODE              = first(JEPCODE),
        selection_method     = paste0("complementary_", group),
        n_unique_pollinators = n_distinct(sourceTaxonName_harm, na.rm = TRUE),
        n_unique_bees        = n_distinct(sourceTaxonName_harm[sourceTaxonType == "bees"], na.rm = TRUE),
        n_unique_hoverflies  = n_distinct(sourceTaxonName_harm[sourceTaxonType == "hoverflies"], na.rm = TRUE),
        n_unique_moths       = n_distinct(sourceTaxonName_harm[sourceTaxonType == "moths"], na.rm = TRUE),
        n_unique_butterflies = n_distinct(sourceTaxonName_harm[sourceTaxonType == "butterflies"], na.rm = TRUE)
      )

    results_list[[group]] <- result
  }

  bind_rows(results_list)
}) %>%
  left_join(total_pollinators_by_region, by = "JEPCODE") %>%
  mutate(
    pct_pollinators = (n_unique_pollinators / total_pollinators) * 100,
    pct_bees        = (n_unique_bees / total_bees) * 100,
    pct_hoverflies  = (n_unique_hoverflies / total_hoverflies) * 100,
    pct_moths       = (n_unique_moths / total_moths) * 100,
    pct_butterflies = (n_unique_butterflies / total_butterflies) * 100
  )

n_iterations <- 100

#### Run sampling iterations for each plant selection method ####
message("Running sampling iterations for each plant selection method...")

group_samples_list <- list()
for (group_name in names(top_plants_by_group)) {
  message(paste("Processing", group_name, "top plants..."))
  group_samples_list[[group_name]] <- map_dfr(
    1:n_iterations,
    ~sample_and_count_unique_pollinators(
      top_plants_by_group[[group_name]],
      predicted_interactions_filtered,
      paste0("top_", group_name),
      total_pollinators_by_region,
      iter = .x
    )
  )
}

group_samples_list_raw_metaweb <- list()
for (group_name in names(top_plants_by_group)) {
  message(paste("Processing", group_name, "top plants..."))
  group_samples_list_raw_metaweb[[group_name]] <- map_dfr(
    1:n_iterations,
    ~sample_and_count_unique_pollinators(
      top_plants_by_group_raw_metaweb[[group_name]],
      predicted_interactions_filtered,
      paste0("top_", group_name, "_raw_metaweb"),
      total_pollinators_by_region,
      iter = .x
    )
  )
}

message("Processing external lists...")
externalList_sample <- map_dfr(
  1:n_iterations,
  ~sample_and_count_unique_pollinators(
    top_plants_external,
    predicted_interactions_filtered,
    "external",
    total_pollinators_by_region,
    iter = .x
  )
)

message("Processing random selections...")
random_sample <- map_dfr(
  1:n_iterations,
  ~sample_and_count_unique_pollinators(
    random_plants_prepList,
    predicted_interactions_filtered,
    "random",
    total_pollinators_by_region,
    iter = .x
  )
)

#### Combine all plant selection results ####
message("Combining all results...")

plants_selected <- bind_rows(
  externalList_sample,
  group_samples_list$all_pollinators,
  group_samples_list$bees,
  group_samples_list$butterflies,
  group_samples_list$hoverflies,
  group_samples_list$moths,
  group_samples_list_raw_metaweb$all_pollinators,
  group_samples_list_raw_metaweb$bees,
  group_samples_list_raw_metaweb$butterflies,
  group_samples_list_raw_metaweb$hoverflies,
  group_samples_list_raw_metaweb$moths,
  random_sample
)

plants_selected <- bind_rows(
  plants_selected,
  complementary_results %>%
    mutate(iteration = 0) %>%
    rename(dataset = selection_method)
) %>%
  mutate(dummy_ext  = ifelse(dataset == "external", 1, 0),
         dummy_pred = ifelse(dataset != "external", 1, 0),
         P2_region     = if_else(is.na(P2_region), "None", P2_region),
         xerces_region = if_else(is.na(xerces_region), "None", xerces_region)) %>%
  rename(method = dataset)

saveRDS(plants_selected, "Temp/plants_selected.rds")

plants_selected <- readRDS("Temp/plants_selected.rds")

#### PART 2: Statistical models ####

# For a taxa, filter plant selections to 5 relevent methods for comparison
build_taxon_spec <- function(taxon, base_data, response_var, trials_var, output_path) {
  
  comp_taxon <- if (taxon == "all_pollinators") "all" else taxon
  
  keep_methods <- c(
    "random",
    "external",
    paste0("top_", taxon),
    paste0("complementary_", comp_taxon),
    paste0("top_", taxon, "_raw_metaweb")
  )
  
  data_taxon <- base_data %>%
    filter(method %in% keep_methods) %>%
    mutate(
      dummy_ext  = as.numeric(method == "external"),
      dummy_pred = as.numeric(method != "external")
    )
  
  formula_taxon <- as.formula(
    paste0(response_var, " | trials(", trials_var, ") ~ 0 + method + ",
           "(0 + dummy_pred | JEPCODE) + ",
           "(0 + dummy_ext | mm(P2_region, xerces_region))")
  )
  
  priors_taxon <- c(
    prior_string("normal(0,1)", class = "b", coef = "methodrandom"),
    prior_string("normal(0,1)", class = "b", coef = "methodexternal"),
    prior_string("normal(0,1)", class = "b", coef = paste0("methodtop_", taxon)),
    prior_string("normal(0,1)", class = "b", coef = paste0("methodcomplementary_", comp_taxon)),
    prior_string("normal(0,1)", class = "b", coef = paste0("methodtop_", taxon, "_raw_metaweb")),
    prior_string("student_t(3,0,2.5)", class = "sd", group = "JEPCODE"),
    prior_string("student_t(3,0,2.5)", class = "sd", group = "mmP2_regionxerces_region")
  )
  
  list(
    data    = data_taxon,
    formula = formula_taxon,
    prior   = priors_taxon,
    output  = output_path,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
}

# Same idea as above, but for the density-interaction models
build_density_spec <- function(taxon, base_data, response_var, trials_var, output_path) {
  
  comp_taxon <- if (taxon == "all_pollinators") "all" else taxon
  
  keep_methods <- c(
    "random",
    "external",
    paste0("top_", taxon),
    paste0("complementary_", comp_taxon),
    paste0("top_", taxon, "_raw_metaweb")
  )
  
  data_taxon <- base_data %>%
    filter(method %in% keep_methods) %>%
    mutate(
      dummy_ext  = as.numeric(method == "external"),
      dummy_pred = as.numeric(method != "external")
    )
  
  formula_taxon <- as.formula(
    paste0(response_var, " | trials(", trials_var, ") ~ 0 + method + method:log_density + ",
           "(0 + dummy_pred | JEPCODE) + ",
           "(0 + dummy_ext | mm(P2_region, xerces_region))")
  )
  
  method_coefs <- c(
    "methodrandom", "methodexternal",
    paste0("methodtop_", taxon),
    paste0("methodcomplementary_", comp_taxon),
    paste0("methodtop_", taxon, "_raw_metaweb")
  )
  interaction_coefs <- paste0(method_coefs, ":log_density")
  
  priors_taxon <- c(
    prior_string("normal(0,1)",   class = "b", coef = method_coefs[1]),
    prior_string("normal(0,1)",   class = "b", coef = method_coefs[2]),
    prior_string("normal(0,1)",   class = "b", coef = method_coefs[3]),
    prior_string("normal(0,1)",   class = "b", coef = method_coefs[4]),
    prior_string("normal(0,1)",   class = "b", coef = method_coefs[5]),
    prior_string("normal(0,0.5)", class = "b", coef = interaction_coefs[1]),
    prior_string("normal(0,0.5)", class = "b", coef = interaction_coefs[2]),
    prior_string("normal(0,0.5)", class = "b", coef = interaction_coefs[3]),
    prior_string("normal(0,0.5)", class = "b", coef = interaction_coefs[4]),
    prior_string("normal(0,0.5)", class = "b", coef = interaction_coefs[5]),
    prior_string("student_t(3,0,2.5)", class = "sd", group = "JEPCODE"),
    prior_string("student_t(3,0,2.5)", class = "sd", group = "mmP2_regionxerces_region")
  )
  
  list(
    data    = data_taxon,
    formula = formula_taxon,
    prior   = priors_taxon,
    output  = output_path,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    threads = 1
  )
}

taxa_response_map <- list(
  all_pollinators = c(response = "n_unique_pollinators", trials = "total_pollinators"),
  bees            = c(response = "n_unique_bees",        trials = "total_bees"),
  hoverflies      = c(response = "n_unique_hoverflies",  trials = "total_hoverflies"),
  butterflies     = c(response = "n_unique_butterflies", trials = "total_butterflies"),
  moths           = c(response = "n_unique_moths",       trials = "total_moths")
)

output_paths_predicted <- list(
  all_pollinators = "Data_Clean/Statistics/revision/brm_compare.rds",
  bees            = "Data_Clean/Statistics/revision/brm_compare_bees.rds",
  hoverflies      = "Data_Clean/Statistics/revision/brm_compare_hoverflies.rds",
  butterflies     = "Data_Clean/Statistics/revision/brm_compare_butterflies.rds",
  moths           = "Data_Clean/Statistics/revision/brm_compare_moths.rds"
)

model_specs_compare_predicted <- imap(taxa_response_map, function(rt, taxon) {
  build_taxon_spec(
    taxon        = taxon,
    base_data    = plants_selected,
    response_var = rt[["response"]],
    trials_var   = rt[["trials"]],
    output_path  = output_paths_predicted[[taxon]]
  )
})
names(model_specs_compare_predicted) <- paste0("predicted_", names(model_specs_compare_predicted))

#### Raw interaction metaweb comparison ####
unique_jepcodes <- ecoregions_simp %>%
  st_drop_geometry() %>%
  pull(JEPCODE) %>%
  unique()

message("Running raw metaweb comparison...")
interactions_nectar_observed <- interactions_nectar %>%
  distinct(targetTaxonName_harm, sourceTaxonName_harm) %>%
  filter(targetTaxonName_harm %in% checklist_cleaned$genus_species,
         sourceTaxonName_harm %in% checklist_cleaned$genus_species) %>%
  left_join(unique(checklist_cleaned[, c("genus_species", "taxon")]),
            by = join_by("sourceTaxonName_harm" == "genus_species")) %>%
  rename(sourceTaxonType = taxon)

complementary_results_rawCompare <- map_dfr(unique_jepcodes, function(jepcode) {
  message(paste("Processing region:", jepcode))
  results_list <- list()

  for (group in c("all", "bees", "butterflies", "hoverflies", "moths")) {
    message(paste("  Finding complementary plants for:", group))
    region_data <- data.frame(JEPCODE = jepcode)

    best_plants <- find_complementary_plants(
      region_data     = region_data,
      interactions_df = predicted_interactions_filtered,
      totals_df       = total_pollinators_by_region,
      target_group    = group,
      n_plants        = 10,
      n_generations   = 100,
      population_size = 200,
      mutation_rate   = 0.1
    )

    plants_df <- data.frame(JEPCODE = jepcode, targetTaxonName_harm = best_plants)

    plants_with_interactions <- plants_df %>%
      left_join(interactions_nectar_observed,
                by = "targetTaxonName_harm",
                relationship = "many-to-many")

    region_pols <- sp_by_region %>%
      filter(JEPCODE == jepcode) %>%
      dplyr::select(genus_species) %>%
      distinct()

    plants_with_interactions <- plants_with_interactions %>%
      inner_join(region_pols, by = c("sourceTaxonName_harm" = "genus_species"))

    result <- plants_with_interactions %>%
      summarise(
        JEPCODE              = first(JEPCODE),
        selection_method     = paste0("complementary_", group),
        n_unique_pollinators = n_distinct(sourceTaxonName_harm, na.rm = TRUE),
        n_unique_bees        = n_distinct(sourceTaxonName_harm[sourceTaxonType == "bees"], na.rm = TRUE),
        n_unique_hoverflies  = n_distinct(sourceTaxonName_harm[sourceTaxonType == "hoverflies"], na.rm = TRUE),
        n_unique_moths       = n_distinct(sourceTaxonName_harm[sourceTaxonType == "moths"], na.rm = TRUE),
        n_unique_butterflies = n_distinct(sourceTaxonName_harm[sourceTaxonType == "butterflies"], na.rm = TRUE)
      )

    results_list[[group]] <- result
  }

  bind_rows(results_list)
}) %>%
  left_join(total_pollinators_by_region, by = "JEPCODE") %>%
  mutate(
    pct_pollinators = (n_unique_pollinators / total_pollinators) * 100,
    pct_bees        = (n_unique_bees / total_bees) * 100,
    pct_hoverflies  = (n_unique_hoverflies / total_hoverflies) * 100,
    pct_moths       = (n_unique_moths / total_moths) * 100,
    pct_butterflies = (n_unique_butterflies / total_butterflies) * 100
  )

group_samples_list_rawCompare <- list()

for (group_name in names(top_plants_by_group)) {
  message(paste("Processing", group_name, "top plants..."))
  group_samples_list_rawCompare[[group_name]] <- map_dfr(
    1:n_iterations,
    ~sample_and_count_unique_pollinators(
      top_plants_by_group[[group_name]],
      interactions_nectar_observed,
      paste0("top_", group_name),
      total_pollinators_by_region,
      iter = .x,
      use_spatial = FALSE,
      sp_by_region_df = sp_by_region
    )
  )
}

group_samples_list_raw_metaweb_rawCompare <- list()

for (group_name in names(top_plants_by_group)) {
  message(paste("Processing", group_name, "top plants..."))
  group_samples_list_raw_metaweb_rawCompare[[group_name]] <- map_dfr(
    1:n_iterations,
    ~sample_and_count_unique_pollinators(
      top_plants_by_group_raw_metaweb[[group_name]],
      interactions_nectar_observed,
      paste0("top_", group_name, "_raw_metaweb"),
      total_pollinators_by_region,
      iter = .x,
      use_spatial = FALSE,
      sp_by_region_df = sp_by_region
    )
  )
}

message("Processing external lists (raw compare)...")
externalList_sample_rawCompare <- map_dfr(
  1:n_iterations,
  ~sample_and_count_unique_pollinators(
    top_plants_external,
    interactions_nectar_observed,
    "external",
    total_pollinators_by_region,
    iter = .x,
    use_spatial = FALSE,
    sp_by_region_df = sp_by_region
  )
)

message("Processing random selections (raw compare)...")
random_sample_rawCompare <- map_dfr(
  1:n_iterations,
  ~sample_and_count_unique_pollinators(
    random_plants_prepList,
    interactions_nectar_observed,
    "random",
    total_pollinators_by_region,
    iter = .x,
    use_spatial = FALSE,
    sp_by_region_df = sp_by_region
  )
)

plants_selected_rawCompare <- bind_rows(
  externalList_sample_rawCompare,
  group_samples_list_rawCompare$all_pollinators,
  group_samples_list_rawCompare$bees,
  group_samples_list_rawCompare$butterflies,
  group_samples_list_rawCompare$hoverflies,
  group_samples_list_rawCompare$moths,
  group_samples_list_raw_metaweb_rawCompare$all_pollinators,
  group_samples_list_raw_metaweb_rawCompare$bees,
  group_samples_list_raw_metaweb_rawCompare$butterflies,
  group_samples_list_raw_metaweb_rawCompare$hoverflies,
  group_samples_list_raw_metaweb_rawCompare$moths,
  random_sample_rawCompare
)

plants_selected_rawCompare <- bind_rows(
  plants_selected_rawCompare,
  complementary_results_rawCompare %>%
    mutate(iteration = 0) %>%
    rename(dataset = selection_method)
) %>%
  mutate(dummy_ext  = ifelse(dataset == "external", 1, 0),
         dummy_pred = ifelse(dataset != "external", 1, 0),
         P2_region     = if_else(is.na(P2_region), "None", P2_region),
         xerces_region = if_else(is.na(xerces_region), "None", xerces_region)) %>%
  rename(method = dataset)

saveRDS(plants_selected_rawCompare, "Temp/plants_selected_rawCompare.rds")

plants_selected_rawCompare <- readRDS("Temp/plants_selected_rawCompare.rds")


output_paths_raw <- list(
  all_pollinators = "Data_Clean/Statistics/revision/brm_compare_rawCompare.rds",
  bees            = "Data_Clean/Statistics/revision/brm_compare_bees_rawCompare.rds",
  hoverflies      = "Data_Clean/Statistics/revision/brm_compare_hoverflies_rawCompare.rds",
  butterflies     = "Data_Clean/Statistics/revision/brm_compare_butterflies_rawCompare.rds",
  moths           = "Data_Clean/Statistics/revision/brm_compare_moths_rawCompare.rds"
)

model_specs_compare_raw <- imap(taxa_response_map, function(rt, taxon) {
  build_taxon_spec(
    taxon        = taxon,
    base_data    = plants_selected_rawCompare,
    response_var = rt[["response"]],
    trials_var   = rt[["trials"]],
    output_path  = output_paths_raw[[taxon]]
  )
})
names(model_specs_compare_raw) <- paste0("raw_", names(model_specs_compare_raw))

#### Data Density Models #####
jepcode_areas <- ecoregions %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    total_area_km2 = sum(area_km2),
    .groups = "drop"
  )
jepcode_density <- native_CA_interactions %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  left_join(jepcode_areas, by = "JEPCODE") %>%
  mutate(
    density = n_interactions / total_area_km2
  )
plants_selected_with_density <- plants_selected %>%
  left_join(jepcode_density, by = "JEPCODE") %>%
  mutate(log_density = log10(density))
plants_selected_with_density_rawCompare <- plants_selected_rawCompare %>%
  left_join(jepcode_density, by = "JEPCODE") %>%
  mutate(log_density = log10(density))

model_specs_density <- list(
  density_predicted = build_density_spec(
    taxon        = "all_pollinators",
    base_data    = plants_selected_with_density,
    response_var = "n_unique_pollinators",
    trials_var   = "total_pollinators",
    output_path  = "Data_Clean/Statistics/revision/brm_density_interaction.rds"
  ),
  density_raw = build_density_spec(
    taxon        = "all_pollinators",
    base_data    = plants_selected_with_density_rawCompare,
    response_var = "n_unique_pollinators",
    trials_var   = "total_pollinators",
    output_path  = "Data_Clean/Statistics/revision/brm_density_interaction_rawCompare.rds"
  )
)

#### Fit everything at once ####
# Unified fit function used for every model, comparison or density
fit_model <- function(spec) {
  cat("Fitting model:", spec$output, "\n")
  
  model <- brm(
    formula = spec$formula,
    family  = beta_binomial(link = "logit"),
    data    = spec$data,
    prior   = spec$prior,
    chains  = 4,
    cores   = 4,
    iter    = 3000,
    warmup  = 1500,
    control = spec$control,
    backend = "cmdstanr",
    seed    = 123,
    refresh = 100,
    silent  = 0
  )
  saveRDS(model, spec$output)
  cat("Saved model:", spec$output, "\n")
  return(spec$output)
}

all_model_specs <- c(
  model_specs_compare_predicted,
  model_specs_compare_raw,
  model_specs_density
)

plan(multisession, workers = 12)
options(mc.cores = 4)  # 4 chains per model

results <- future.apply::future_lapply(
  all_model_specs,
  fit_model,
  future.seed = TRUE
)

# Reset to sequential processing when done
plan(sequential)

#### Part 3: Plotting ####
#### Support by taxa, NECTAR interaction networks ####
brm_compare <- readRDS("Data_Clean/Statistics/revision/brm_compare.rds")
brm_compare_bees <- readRDS("Data_Clean/Statistics/revision/brm_compare_bees.rds")
brm_compare_hoverflies <- readRDS("Data_Clean/Statistics/revision/brm_compare_hoverflies.rds")
brm_compare_butterflies <- readRDS("Data_Clean/Statistics/revision/brm_compare_butterflies.rds")
brm_compare_moths <- readRDS("Data_Clean/Statistics/revision/brm_compare_moths.rds")

# read in plants selected
plants_selected <- readRDS("Temp/plants_selected.rds")
plants_selected_rawCompare <- readRDS("Temp/plants_selected_rawCompare.rds")

## 1. Extract fixed effects by model/taxonomic group
extract_method_effects <- function(model, taxa_name, methods_to_extract) {
  # Get posterior draws
  post_draws <- as_draws_df(model)
  
  # Extract the specific method coefficients
  method_cols <- paste0("b_", methods_to_extract)
  
  method_draws <- post_draws %>%
    as.data.frame() %>%
    dplyr::select(all_of(method_cols))
  
  # Calculate mean and credible intervals in logit space, then transform
  fixed_effects <- method_draws %>%
    summarise(across(everything(), list(
      mean = ~mean(.x),
      lower = ~quantile(.x, 0.025),
      upper = ~quantile(.x, 0.975)
    ))) %>%
    # Transform to probability scale
    mutate(across(everything(), plogis)) %>%
    # Reshape to long format
    pivot_longer(
      everything(),
      names_to = c("method", ".value"),
      names_pattern = "b_(.*)_(mean|lower|upper)"
    ) %>%
    mutate(
      taxa = taxa_name,
      category = "Global"
    ) %>%
    rename(
      mean_response = mean,
      ymin = lower,
      ymax = upper
    )
  
  return(fixed_effects)
}

# Define which methods to extract for each taxonomic group
methods_by_taxa <- list(
  "All Pollinators" = c("methodrandom", "methodexternal",
                        "methodtop_all_pollinators", "methodcomplementary_all", "methodtop_all_pollinators_raw_metaweb"),
  "Bees" = c("methodrandom", "methodexternal",
             "methodtop_bees", "methodcomplementary_bees", "methodtop_bees_raw_metaweb"),
  "Butterflies" = c("methodrandom", "methodexternal",
                    "methodtop_butterflies", "methodcomplementary_butterflies", "methodtop_butterflies_raw_metaweb"),
  "Moths" = c("methodrandom", "methodexternal",
              "methodtop_moths", "methodcomplementary_moths", "methodtop_moths_raw_metaweb"),
  "Hoverflies" = c("methodrandom", "methodexternal",
                   "methodtop_hoverflies", "methodcomplementary_hoverflies", "methodtop_hoverflies_raw_metaweb")
)

# Extract fixed effects for all models
fixed_effects_all <- bind_rows(
  extract_method_effects(brm_compare, "All Pollinators",
                         methods_by_taxa[["All Pollinators"]]),
  extract_method_effects(brm_compare_bees, "Bees",
                         methods_by_taxa[["Bees"]]),
  extract_method_effects(brm_compare_butterflies, "Butterflies",
                         methods_by_taxa[["Butterflies"]]),
  extract_method_effects(brm_compare_moths, "Moths",
                         methods_by_taxa[["Moths"]]),
  extract_method_effects(brm_compare_hoverflies, "Hoverflies",
                         methods_by_taxa[["Hoverflies"]])
)

### Aside: What is absolute and relative increase compared to random for each method-taxa combo?
# Get random baseline per taxa
random_baseline <- fixed_effects_all %>%
  filter(grepl("random", method)) %>%
  dplyr::select(taxa, random_mean = mean_response)

# Join and compute increases
fixed_effects_all %>%
  filter(!grepl("random", method)) %>%
  left_join(random_baseline, by = "taxa") %>%
  mutate(
    absolute_increase = (mean_response - random_mean)*100,
    relative_increase = ((mean_response - random_mean) / random_mean)*100
  )
###

### Aside: What is absolute and relative increase compared to external list for each method-taxa combo?
# Get external list baseline per taxa
external_baseline <- fixed_effects_all %>%
  filter(grepl("external", method)) %>%
  dplyr::select(taxa, random_mean = mean_response)

# Join and compute increases
fixed_effects_all %>%
  filter(!grepl("external", method)) %>%
  filter(!grepl("random", method)) %>%
  left_join(external_baseline, by = "taxa") %>%
  mutate(
    absolute_increase = (mean_response - random_mean)*100,
    relative_increase = ((mean_response - random_mean) / random_mean)*100
  )
###

## 2. Prepare plotting data with positions

# Define order to display taxa groups and method for selecting plants
taxa_levels <- c("All Pollinators", "Bees", "Butterflies", "Moths", "Hoverflies")
method_display_order <- c("methodrandom", "methodexternal", "methodrawmetaweb",
                          "methodtop", "methodcomplementary")

# Clean up method names for display and create standardized method type
plot_data <- fixed_effects_all %>%
  mutate(
    # Extract the base method type (removing taxa-specific suffixes)
    method_type = case_when(
      grepl("^methodrandom$", method) ~ "methodrandom",
      grepl("^methodexternal$", method) ~ "methodexternal",
      grepl("raw_metaweb$", method) ~ "methodrawmetaweb",
      grepl("^methodtop_", method) ~ "methodtop",
      grepl("^methodcomplementary_", method) ~ "methodcomplementary",
      TRUE ~ method
    ),
    # Create display labels
    method_label = case_when(
      method_type == "methodrandom" ~ "Random",
      method_type == "methodexternal" ~ "External",
      method_type == "methodrawmetaweb" ~ "Raw Metaweb",
      method_type == "methodtop" ~ "NECTAR Top Species",
      method_type == "methodcomplementary" ~ "NECTAR Complementary",
      TRUE ~ method
    ),
    # Set levels of factor
    taxa = factor(taxa, levels = taxa_levels),
    method_type = factor(method_type, levels = method_display_order)
  ) %>%
  arrange(taxa, method_type)

## 3. set x positions and staggering manually

# Parameters for positioning
taxa_spacing <- 5        # Space between taxonomic groups
group_span <- (4-1) * 0.75
method_spacing <- group_span / (length(method_display_order) - 1)    # Space between methods within a group
stagger_offset <- 0.08   # Offset for each method point

plot_data <- plot_data %>%
  mutate(
    # Base position for each taxonomic group
    taxa_position = (as.numeric(taxa) - 1) * taxa_spacing + 1,
    # Offset within group based on method
    method_offset = (as.numeric(method_type) - 1) * method_spacing,
    # Additional small stagger
    stagger = (as.numeric(method_type) - (length(method_display_order) + 1) / 2) * stagger_offset,
    # Final x position
    x_position = taxa_position + method_offset + stagger
  )

## 4. Create x-axis breaks and labels
# Create breaks at the center of each taxonomic group
x_breaks <- seq(1, by = taxa_spacing, length.out = length(taxa_levels)) +
  (length(method_display_order) - 1) * method_spacing / 2

x_labels <- taxa_levels

## 5. Create the plot
improvement_predictedNetwork <- ggplot(plot_data,
                                       aes(x = x_position, y = mean_response,
                                           color = taxa,
                                           shape = method_label)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.1,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(add = c(1, 1))  # Increased from 0.5 to 1
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, by = 0.1)
  ) +
  scale_color_manual(
    values = c(
      "All Pollinators" = "#000000",
      "Bees" = "#FFCC00",
      "Butterflies" = "#E69190",
      "Moths" = "#BD9ECD",
      "Hoverflies" = "#F76700"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Random" = 16,      # Circle
      "External" = 17,    # Triangle
      "Raw Metaweb" = 8,  # Star
      "NECTAR Top Species" = 15, # Square
      "NECTAR Complementary" = 18 # Diamond
    )
  ) +
  labs(
    x = "",
    y = "% Pollinator Community Supported",
    color = "Taxonomic Group",
    shape = "Method"
  ) +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  ) +
  # Add vertical lines between taxonomic groups
  geom_vline(
    xintercept = seq(1, by = taxa_spacing, length.out = length(taxa_levels) - 1) +
      taxa_spacing / 2 + (length(method_display_order) - 1) * method_spacing / 2,
    linetype = "dashed",
    color = "gray70",
    linewidth = 0.3
  )

ggsave("Figures/P1_Comparison_ExternalLists/improvement_predictedNetwork.pdf",
       plot = improvement_predictedNetwork,
       width = 185,
       height = 60,
       units = "mm",
       dpi = 600)


#### Data Density vs Improvement Over Random, NECTAR interaction networks ####
brm_density_interaction <- readRDS("Data_Clean/Statistics/revision/brm_density_interaction.rds")

## 1. Calculate jepcode level data density
jepcode_areas <- ecoregions %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    total_area_km2 = sum(area_km2),
    .groups = "drop"
  )

# density of interaction data for this
jepcode_density <- native_CA_interactions %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  left_join(jepcode_areas, by = "JEPCODE") %>%
  mutate(
    density = n_interactions / total_area_km2,
    log_density = log10(density)
  )

## 2. Create prediction grid for marginal predictions
density_grid <- expand_grid(
  # create grid of interaction density values between min and max observed for predicted draws
  log_density = seq(min(jepcode_density$log_density),
                    max(jepcode_density$log_density),
                    length.out = 100),
  method = c("random", "external", "top_all_pollinators", "complementary_all", "top_all_pollinators_raw_metaweb")
) %>%
  mutate(
    density = 10^log_density,
    dummy_pred = ifelse(method == "external", 0, 1),
    dummy_ext = ifelse(method == "external", 1, 0),
    total_pollinators = 1
  )

## 3. Get marginal posterior predictions
grid_predictions <- add_epred_draws(
  density_grid,
  brm_density_interaction,
  re_formula = NA
)

## 4. Calculate expected improvement in pol community supported for each method relative to random
grid_improvements <- grid_predictions %>%
  ungroup() %>%
  dplyr::select(log_density, density, method, .draw, .epred) %>%
  pivot_wider(
    names_from = method,
    values_from = .epred
  ) %>%
  mutate(
    External_vs_Random = external - random,
    Top_vs_Random = top_all_pollinators - random,
    Complementary_vs_Random = complementary_all - random,
    RawMetaweb_vs_Random = top_all_pollinators_raw_metaweb - random
  ) %>%
  dplyr::select(log_density, density, .draw, ends_with("vs_Random")) %>%
  pivot_longer(
    cols = ends_with("vs_Random"),
    names_to = "comparison",
    values_to = "improvement"
  ) %>%
  mutate(
    comparison = case_when(
      comparison == "External_vs_Random" ~ "External vs Random",
      comparison == "Top_vs_Random" ~ "Top Species vs Random",
      comparison == "Complementary_vs_Random" ~ "Complementary vs Random",
      comparison == "RawMetaweb_vs_Random" ~ "Raw Metaweb vs Random"
    ),
    comparison = factor(comparison,
                        levels = c("External vs Random",
                                   "Top Species vs Random",
                                   "Complementary vs Random",
                                   "Raw Metaweb vs Random"))
  )

## 5. Summarize for smooth lines
grid_summary <- grid_improvements %>%
  group_by(comparison, density, log_density) %>%
  summarise(
    fit = mean(improvement),
    lower = quantile(improvement, 0.025),
    upper = quantile(improvement, 0.975),
    .groups = "drop"
  )

## 6. For each region, get observed expected improvement from original plants draws
observed_improvements <- plants_selected %>%
  filter(method %in% c("random", "external", "top_all_pollinators", "complementary_all", "top_all_pollinators_raw_metaweb")) %>%
  group_by(JEPCODE, method) %>%
  summarise(
    pct_pollinators = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = pct_pollinators
  ) %>%
  mutate(
    External_vs_Random = external - random,
    Top_vs_Random = top_all_pollinators - random,
    Complementary_vs_Random = complementary_all - random,
    RawMetaweb_vs_Random = top_all_pollinators_raw_metaweb - random
  ) %>%
  dplyr::select(JEPCODE, ends_with("vs_Random")) %>%
  pivot_longer(
    cols = ends_with("vs_Random"),
    names_to = "comparison",
    values_to = "improvement"
  ) %>%
  mutate(
    comparison = case_when(
      comparison == "External_vs_Random" ~ "External vs Random",
      comparison == "Top_vs_Random" ~ "Top Species vs Random",
      comparison == "Complementary_vs_Random" ~ "Complementary vs Random",
      comparison == "RawMetaweb_vs_Random" ~ "Raw Metaweb vs Random"
    ),
    comparison = factor(comparison,
                        levels = c("External vs Random",
                                   "Top Species vs Random",
                                   "Complementary vs Random",
                                   "Raw Metaweb vs Random"))
  ) %>%
  left_join(jepcode_density, by = "JEPCODE")

## 7. Create the plot
density_improvement_plot <- ggplot() +
  # "observed" - from actual regional draws
  geom_point(
    data = observed_improvements,
    aes(x = density, y = improvement, shape = comparison),
    size = 2,
    alpha = 0.5
  ) +
  # model predictions w.r.t interaction data density
  geom_line(
    data = grid_summary,
    aes(x = density, y = fit * 100, linetype = comparison),
    linewidth = 1
  ) +
  # Credible intervals
  geom_ribbon(
    data = grid_summary,
    aes(x = density, ymin = lower * 100, ymax = upper * 100,
        group = comparison),
    alpha = 0.15,
    fill = "gray50"
  ) +
  scale_x_log10(
    name = "Interaction Data Density (records/km²)",
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  scale_y_continuous(
    name = "Improvement Over Random (%)",
    breaks = seq(-20, 40, by = 10)
  ) +
  scale_linetype_manual(
    values = c(
      "External vs Random" = "solid",
      "Top Species vs Random" = "dashed",
      "Complementary vs Random" = "dotted",
      "Raw Metaweb vs Random" = "twodash"
    )
  ) +
  scale_shape_manual(
    values = c(
      "External vs Random" = 17,
      "Top Species vs Random" = 15,
      "Complementary vs Random" = 18,
      "Raw Metaweb vs Random" = 8
    )
  ) +
  labs(
    linetype = "Comparison",
    shape = "Comparison"
  ) +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave("Figures/P1_Comparison_ExternalLists/density_vs_improvement.pdf",
       plot = density_improvement_plot,
       width = 100,
       height = 60,
       units = "mm",
       dpi = 600)


#### Plot 1: Interaction data density ####
# Ecoregions boundary
ecoregions_boundary <- ecoregions %>%
  st_union() %>%
  nngeo::st_remove_holes()

# Create hex grid
hex_grid <- st_make_grid(
  ecoregions_boundary,
  cellsize = c(0.5, 0.5),
  square = FALSE,
  flat_topped = FALSE
) %>%
  st_sf(hex_id = 1:length(.))

hex_density <- hex_grid %>%
  st_join(native_CA_interactions) %>%
  group_by(hex_id) %>%
  summarise(
    count = sum(!is.na(JEPCODE)),
    .groups = "drop"
  ) %>%
  filter(count > 0)

dataDensity_plot <- ggplot() +
  geom_sf(
    data = hex_density,
    aes(fill = count),
    color = NA
  ) +
  geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_gradient2(
    low = "white",
    mid = "#EAA420",
    high = "#331E0F",
    midpoint = log10(100),
    name = "Interaction\nRecords (#)",
    na.value = "grey80",
    trans = "log10",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, family = "Arial"),
    plot.subtitle = element_text(size = 12, color = "grey40")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_sf()

ggsave("Figures/P1_Comparison_ExternalLists/data_density.svg",
       plot = dataDensity_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)

ggsave("Figures/P1_Comparison_ExternalLists/data_density.pdf",
       plot = dataDensity_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)

#### Plot 2: Area-weighted pct_improvement (Complementary vs Random), NECTAR interaction networks ####
## Get complementary results (JEPCODE level only)
complementary_by_jepcode <- plants_selected %>%
  filter(method == "complementary_all") %>%
  group_by(JEPCODE) %>%
  summarise(
    complementary_pct = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  )

## Get random selection results, mean per jepcode
random_by_jepcode <- plants_selected %>%
  filter(method == "random") %>%
  group_by(JEPCODE) %>%
  summarise(
    random_pct = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  )

## Calculate improvement at JEPCODE level
improvement_by_jepcode <- complementary_by_jepcode %>%
  left_join(random_by_jepcode, by = "JEPCODE") %>%
  mutate(
    pct_improvement = complementary_pct - random_pct,
    relative_improvement = ((complementary_pct - random_pct) / random_pct) * 100
  )

## Join with ALL region combinations for each JEPCODE
map_data <- eco_all_combos %>%
  left_join(improvement_by_jepcode, by = "JEPCODE")

# Ecoregions boundary
ecoregions_boundary <- ecoregions %>%
  st_union() %>%
  nngeo::st_remove_holes()

# Create hex grid
hex_grid <- st_make_grid(
  ecoregions_boundary,
  cellsize = c(0.5, 0.5),
  square = FALSE,
  flat_topped = FALSE
) %>%
  st_sf(hex_id = 1:length(.))

# Calculate area-weighted improvement per hex
hex_intersections <- st_intersection(hex_grid, map_data)

hex_pct_improvement <- hex_intersections %>%
  mutate(
    intersection_area = st_area(.)
  ) %>%
  st_drop_geometry() %>%
  group_by(hex_id) %>%
  summarise(
    weighted_pct_improvement = sum(pct_improvement * as.numeric(intersection_area), na.rm = TRUE) /
      sum(as.numeric(intersection_area), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(hex_grid, ., by = "hex_id") %>%
  filter(!is.na(weighted_pct_improvement))

pct_improvement_plot <- ggplot() +
  geom_sf(
    data = hex_pct_improvement,
    aes(fill = weighted_pct_improvement),
    color = NA
  ) +
  geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_gradient2(
    low = "#EAFFCC",
    mid = "#94C357",
    high = "#3D5E29",
    midpoint = 35,  # Adjust based on your data range
    name = "Complementary vs Random\n% Improvement",
    na.value = "grey80"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, family = "Arial"),
    plot.subtitle = element_text(size = 12, color = "grey40")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_sf()

ggsave("Figures/P1_Comparison_ExternalLists/pct_improvement_complementary_vs_random_hex.svg",
       plot = pct_improvement_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)

ggsave("Figures/P1_Comparison_ExternalLists/pct_improvement_complementary_vs_random_hex.pdf",
       plot = pct_improvement_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)



#### Supplementary, Improvement Plots ####
#### Plot: Support by taxa, raw metaweb ####
brm_compare <- readRDS("Data_Clean/Statistics/revision/brm_compare_rawCompare.rds")
brm_compare_bees <- readRDS("Data_Clean/Statistics/revision/brm_compare_bees_rawCompare.rds")
brm_compare_hoverflies <- readRDS("Data_Clean/Statistics/revision/brm_compare_hoverflies_rawCompare.rds")
brm_compare_butterflies <- readRDS("Data_Clean/Statistics/revision/brm_compare_butterflies_rawCompare.rds")
brm_compare_moths <- readRDS("Data_Clean/Statistics/revision/brm_compare_moths_rawCompare.rds")

## 1. Extract fixed effects by model/taxonomic group
extract_method_effects <- function(model, taxa_name, methods_to_extract) {
  # Get posterior draws
  post_draws <- as_draws_df(model)
  
  # Extract the specific method coefficients
  method_cols <- paste0("b_", methods_to_extract)
  
  method_draws <- post_draws %>%
    as.data.frame() %>%
    dplyr::select(all_of(method_cols))
  
  # Calculate mean and credible intervals in logit space, then transform
  fixed_effects <- method_draws %>%
    summarise(across(everything(), list(
      mean = ~mean(.x),
      lower = ~quantile(.x, 0.025),
      upper = ~quantile(.x, 0.975)
    ))) %>%
    # Transform to probability scale
    mutate(across(everything(), plogis)) %>%
    # Reshape to long format
    pivot_longer(
      everything(),
      names_to = c("method", ".value"),
      names_pattern = "b_(.*)_(mean|lower|upper)"
    ) %>%
    mutate(
      taxa = taxa_name,
      category = "Global"
    ) %>%
    rename(
      mean_response = mean,
      ymin = lower,
      ymax = upper
    )
  
  return(fixed_effects)
}

# Define which methods to extract for each taxonomic group
methods_by_taxa <- list(
  "All Pollinators" = c("methodrandom", "methodexternal",
                        "methodtop_all_pollinators", "methodcomplementary_all", "methodtop_all_pollinators_raw_metaweb"),
  "Bees" = c("methodrandom", "methodexternal",
             "methodtop_bees", "methodcomplementary_bees", "methodtop_bees_raw_metaweb"),
  "Butterflies" = c("methodrandom", "methodexternal",
                    "methodtop_butterflies", "methodcomplementary_butterflies", "methodtop_butterflies_raw_metaweb"),
  "Moths" = c("methodrandom", "methodexternal",
              "methodtop_moths", "methodcomplementary_moths", "methodtop_moths_raw_metaweb"),
  "Hoverflies" = c("methodrandom", "methodexternal",
                   "methodtop_hoverflies", "methodcomplementary_hoverflies", "methodtop_hoverflies_raw_metaweb")
)

# Extract fixed effects for all models
fixed_effects_all <- bind_rows(
  extract_method_effects(brm_compare, "All Pollinators",
                         methods_by_taxa[["All Pollinators"]]),
  extract_method_effects(brm_compare_bees, "Bees",
                         methods_by_taxa[["Bees"]]),
  extract_method_effects(brm_compare_butterflies, "Butterflies",
                         methods_by_taxa[["Butterflies"]]),
  extract_method_effects(brm_compare_moths, "Moths",
                         methods_by_taxa[["Moths"]]),
  extract_method_effects(brm_compare_hoverflies, "Hoverflies",
                         methods_by_taxa[["Hoverflies"]])
)

## 2. Prepare plotting data with positions

# Define order to display taxa groups and method for selecting plants
taxa_levels <- c("All Pollinators", "Bees", "Butterflies", "Moths", "Hoverflies")
method_display_order <- c("methodrandom", "methodexternal", "methodrawmetaweb",
                          "methodtop", "methodcomplementary")

# Clean up method names for display and create standardized method type
plot_data <- fixed_effects_all %>%
  mutate(
    # Extract the base method type (removing taxa-specific suffixes)
    method_type = case_when(
      grepl("^methodrandom$", method) ~ "methodrandom",
      grepl("^methodexternal$", method) ~ "methodexternal",
      grepl("raw_metaweb$", method) ~ "methodrawmetaweb",
      grepl("^methodtop_", method) ~ "methodtop",
      grepl("^methodcomplementary_", method) ~ "methodcomplementary",
      TRUE ~ method
    ),
    # Create display labels
    method_label = case_when(
      method_type == "methodrandom" ~ "Random",
      method_type == "methodexternal" ~ "External",
      method_type == "methodrawmetaweb" ~ "Raw Metaweb",
      method_type == "methodtop" ~ "Top Species",
      method_type == "methodcomplementary" ~ "Complementary",
      TRUE ~ method
    ),
    # Set levels of factor
    taxa = factor(taxa, levels = taxa_levels),
    method_type = factor(method_type, levels = method_display_order)
  ) %>%
  arrange(taxa, method_type)

## 3. set x positions and staggering manually

# Parameters for positioning
taxa_spacing <- 5        # Space between taxonomic groups
group_span <- (4 - 1) * 0.75
method_spacing <- group_span / (length(method_display_order) - 1) # Space between methods within a group
stagger_offset <- 0.08   # Offset for each method point

plot_data <- plot_data %>%
  mutate(
    # Base position for each taxonomic group
    taxa_position = (as.numeric(taxa) - 1) * taxa_spacing + 1,
    # Offset within group based on method
    method_offset = (as.numeric(method_type) - 1) * method_spacing,
    # Additional small stagger
    stagger = (as.numeric(method_type) - (length(method_display_order) + 1) / 2) * stagger_offset,
    # Final x position
    x_position = taxa_position + method_offset + stagger
  )

## 4. Create x-axis breaks and labels
# Create breaks at the center of each taxonomic group
x_breaks <- seq(1, by = taxa_spacing, length.out = length(taxa_levels)) +
  (length(method_display_order) - 1) * method_spacing / 2

x_labels <- taxa_levels

## 5. Create the plot
improvement_predictedNetwork <- ggplot(plot_data,
                                       aes(x = x_position, y = mean_response,
                                           color = taxa,
                                           shape = method_label)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.1,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(add = c(1, 1))  # Increased from 0.5 to 1
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 0.7),
    breaks = seq(0, 0.7, by = 0.1)
  ) +
  scale_color_manual(
    values = c(
      "All Pollinators" = "#000000",
      "Bees" = "#FFCC00",
      "Butterflies" = "#E69190",
      "Moths" = "#BD9ECD",
      "Hoverflies" = "#F76700"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Random" = 16,      # Circle
      "External" = 17,    # Triangle
      "Raw Metaweb" = 8,  # Star
      "Top Species" = 15, # Square
      "Complementary" = 18 # Diamond
    )
  ) +
  labs(
    x = "",
    y = "% Pollinator Community Supported",
    color = "Taxonomic Group",
    shape = "Method"
  ) +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  ) +
  # Add vertical lines between taxonomic groups
  geom_vline(
    xintercept = seq(1, by = taxa_spacing, length.out = length(taxa_levels) - 1) +
      taxa_spacing / 2 + (length(method_display_order) - 1) * method_spacing / 2,
    linetype = "dashed",
    color = "gray70",
    linewidth = 0.3
  )

ggsave("Figures/Supplementary/improvement_predictedNetwork_rawCompare.pdf",
       plot = improvement_predictedNetwork,
       width = 185,
       height = 60,
       units = "mm",
       dpi = 600)

#### Plot: Data Density vs Improvement Over Random, raw metaweb ####
brm_density_interaction <- readRDS("Data_Clean/Statistics/revision/brm_density_interaction_rawCompare.rds")

## 1. Calculate jepcode level data density
jepcode_areas <- ecoregions %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    total_area_km2 = sum(area_km2),
    .groups = "drop"
  )

# density of interaction data for this
jepcode_density <- native_CA_interactions %>%
  st_drop_geometry() %>%
  group_by(JEPCODE) %>%
  summarise(
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  left_join(jepcode_areas, by = "JEPCODE") %>%
  mutate(
    density = n_interactions / total_area_km2,
    log_density = log10(density)
  )

## 2. Create prediction grid for marginal predictions
density_grid <- expand_grid(
  # create grid of interaction density values between min and max observed for predicted draws
  log_density = seq(min(jepcode_density$log_density),
                    max(jepcode_density$log_density),
                    length.out = 100),
  method = c("random", "external", "top_all_pollinators", "complementary_all", "top_all_pollinators_raw_metaweb")
) %>%
  mutate(
    density = 10^log_density,
    dummy_pred = ifelse(method == "external", 0, 1),
    dummy_ext = ifelse(method == "external", 1, 0),
    total_pollinators = 1
  )

## 3. Get marginal posterior predictions
grid_predictions <- add_epred_draws(
  density_grid,
  brm_density_interaction,
  re_formula = NA
)

## 4. Calculate expected improvement in pol community supported for each method relative to random
grid_improvements <- grid_predictions %>%
  ungroup() %>%
  dplyr::select(log_density, density, method, .draw, .epred) %>%
  pivot_wider(
    names_from = method,
    values_from = .epred
  ) %>%
  mutate(
    External_vs_Random = external - random,
    Top_vs_Random = top_all_pollinators - random,
    Complementary_vs_Random = complementary_all - random,
    RawMetaweb_vs_Random = top_all_pollinators_raw_metaweb - random
  ) %>%
  dplyr::select(log_density, density, .draw, ends_with("vs_Random")) %>%
  pivot_longer(
    cols = ends_with("vs_Random"),
    names_to = "comparison",
    values_to = "improvement"
  ) %>%
  mutate(
    comparison = case_when(
      comparison == "External_vs_Random" ~ "External vs Random",
      comparison == "Top_vs_Random" ~ "Top Species vs Random",
      comparison == "Complementary_vs_Random" ~ "Complementary vs Random",
      comparison == "RawMetaweb_vs_Random" ~ "Raw Metaweb vs Random"
    ),
    comparison = factor(comparison,
                        levels = c("External vs Random",
                                   "Top Species vs Random",
                                   "Complementary vs Random",
                                   "Raw Metaweb vs Random"))
  )

## 5. Summarize for smooth lines
grid_summary <- grid_improvements %>%
  group_by(comparison, density, log_density) %>%
  summarise(
    fit = mean(improvement),
    lower = quantile(improvement, 0.025),
    upper = quantile(improvement, 0.975),
    .groups = "drop"
  )

## 6. For each region, get observed expected improvement from original plants draws
observed_improvements <- plants_selected_rawCompare %>%
  filter(method %in% c("random", "external", "top_all_pollinators", "complementary_all", "top_all_pollinators_raw_metaweb")) %>%
  group_by(JEPCODE, method) %>%
  summarise(
    pct_pollinators = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = pct_pollinators
  ) %>%
  mutate(
    External_vs_Random = external - random,
    Top_vs_Random = top_all_pollinators - random,
    Complementary_vs_Random = complementary_all - random,
    RawMetaweb_vs_Random = top_all_pollinators_raw_metaweb - random
  ) %>%
  dplyr::select(JEPCODE, ends_with("vs_Random")) %>%
  pivot_longer(
    cols = ends_with("vs_Random"),
    names_to = "comparison",
    values_to = "improvement"
  ) %>%
  mutate(
    comparison = case_when(
      comparison == "External_vs_Random" ~ "External vs Random",
      comparison == "Top_vs_Random" ~ "Top Species vs Random",
      comparison == "Complementary_vs_Random" ~ "Complementary vs Random",
      comparison == "RawMetaweb_vs_Random" ~ "Raw Metaweb vs Random"
    ),
    comparison = factor(comparison,
                        levels = c("External vs Random",
                                   "Top Species vs Random",
                                   "Complementary vs Random",
                                   "Raw Metaweb vs Random"))
  ) %>%
  left_join(jepcode_density, by = "JEPCODE")

## 7. Create the plot
density_improvement_plot <- ggplot() +
  # "observed" - from actual regional draws
  geom_point(
    data = observed_improvements,
    aes(x = density, y = improvement, shape = comparison),
    size = 2,
    alpha = 0.5
  ) +
  # model predictions w.r.t interaction data density
  geom_line(
    data = grid_summary,
    aes(x = density, y = fit * 100, linetype = comparison),
    linewidth = 1
  ) +
  # Credible intervals
  geom_ribbon(
    data = grid_summary,
    aes(x = density, ymin = lower * 100, ymax = upper * 100,
        group = comparison),
    alpha = 0.15,
    fill = "gray50"
  ) +
  scale_x_log10(
    name = "Interaction Data Density (records/km²)",
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  scale_y_continuous(
    name = "Improvement Over Random (%)",
    breaks = seq(-20, 40, by = 10)
  ) +
  scale_linetype_manual(
    values = c(
      "External vs Random" = "solid",
      "Top Species vs Random" = "dashed",
      "Complementary vs Random" = "dotted",
      "Raw Metaweb vs Random" = "twodash"
    )
  ) +
  scale_shape_manual(
    values = c(
      "External vs Random" = 17,
      "Top Species vs Random" = 15,
      "Complementary vs Random" = 18,
      "Raw Metaweb vs Random" = 8
    )
  ) +
  labs(
    linetype = "Comparison",
    shape = "Comparison"
  ) +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave("Figures/Supplementary/density_vs_improvement_rawCompare.pdf",
       plot = density_improvement_plot,
       width = 100,
       height = 60,
       units = "mm",
       dpi = 600)

### Area-weighted pct_improvement (Complementary vs Random), raw metaweb ####
## Get complementary results (JEPCODE level only)
complementary_by_jepcode <- plants_selected_rawCompare %>%
  filter(method == "complementary_all") %>%
  group_by(JEPCODE) %>%
  summarise(
    complementary_pct = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  )

## Get random selection results, mean per jepcode
random_by_jepcode <- plants_selected_rawCompare %>%
  filter(method == "random") %>%
  group_by(JEPCODE) %>%
  summarise(
    random_pct = mean(pct_pollinators, na.rm = TRUE),
    .groups = "drop"
  )

## Calculate improvement at JEPCODE level
improvement_by_jepcode <- complementary_by_jepcode %>%
  left_join(random_by_jepcode, by = "JEPCODE") %>%
  mutate(
    pct_improvement = complementary_pct - random_pct,
    relative_improvement = ((complementary_pct - random_pct) / random_pct) * 100
  )

## Join with ALL region combinations for each JEPCODE
map_data <- eco_all_combos %>%
  left_join(improvement_by_jepcode, by = "JEPCODE")

# Ecoregions boundary
ecoregions_boundary <- ecoregions %>%
  st_union() %>%
  nngeo::st_remove_holes()

# Create hex grid
hex_grid <- st_make_grid(
  ecoregions_boundary,
  cellsize = c(0.5, 0.5),
  square = FALSE,
  flat_topped = FALSE
) %>%
  st_sf(hex_id = 1:length(.))

# Calculate area-weighted improvement per hex
hex_intersections <- st_intersection(hex_grid, map_data)

hex_pct_improvement <- hex_intersections %>%
  mutate(
    intersection_area = st_area(.)
  ) %>%
  st_drop_geometry() %>%
  group_by(hex_id) %>%
  summarise(
    weighted_pct_improvement = sum(pct_improvement * as.numeric(intersection_area), na.rm = TRUE) /
      sum(as.numeric(intersection_area), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(hex_grid, ., by = "hex_id") %>%
  filter(!is.na(weighted_pct_improvement))

pct_improvement_plot <- ggplot() +
  geom_sf(
    data = hex_pct_improvement,
    aes(fill = weighted_pct_improvement),
    color = NA
  ) +
  geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_gradient2(
    low = "#EAFFCC",
    mid = "#94C357",
    high = "#3D5E29",
    midpoint = 10,  # Adjust based on your data range
    name = "Complementary vs Random\n% Improvement",
    na.value = "grey80"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, family = "Arial"),
    plot.subtitle = element_text(size = 12, color = "grey40")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_sf()

ggsave("Figures/Supplementary/pct_improvement_complementary_vs_random_hex_rawCompare.svg",
       plot = pct_improvement_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)

ggsave("Figures/Supplementary/pct_improvement_complementary_vs_random_hex_rawCompare.pdf",
       plot = pct_improvement_plot,
       width = 200,
       height = 100,
       units = "mm",
       dpi = 2000)
