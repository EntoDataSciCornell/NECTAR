#### Load packages ####
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(sf)
library(terra)
library(bdc)

#### Supporting scripts ####
source("Code/99_Supporting/name_harmonization.R")
gnverifier_path <- "/home/tb625/bin"

#### Read in data ####
# Core checklist
checklist <- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv") %>% 
  filter(taxon != "hummingbirds")

# Regional occurrence
sp_by_region <- read.csv("Data_Clean/Species_Checklists/ecoregion_checklist.csv")

# Phenometrics
phenology <- fread("Data_Clean/Phenology/phenology_best_available.csv") %>% distinct()

# Predicted interactions
predicted_interactions <- fread("Data_Clean/Interactions/Predicted_Interactions.csv")

# Calscape raw plant data (for drought tolerance and nursery availability)
calscape_plants_list <- read.csv("Data_Clean/Species_Traits_Attributes/Calscape_Plants_Harmonized.csv")
calscape_plant_data  <- read.csv("Data_Raw/Species_Traits_Attributes/plants_table_5-7.csv")
calscape_plant_data  <- left_join(calscape_plant_data, calscape_plants_list, by = "species")

# Plant community data (for chaparral flag)
calflora_communities <- read.csv("Data_Clean/Species_Traits_Attributes/Calflora_Communities.csv")

# Plant checklist (for CRPR, listing status)
plant_checklist_full <- read.csv("Data_Clean/Species_Checklists/plant_checklist.csv")

# Animal checklists (for listing status)
bee_checklist       <- read.csv("Data_Clean/Species_Checklists/bee_checklist.csv")
butterfly_checklist <- read.csv("Data_Clean/Species_Checklists/butterfly_checklist.csv")
moth_checklist      <- read.csv("Data_Clean/Species_Checklists/moth_checklist.csv")
hoverfly_checklist  <- read.csv("Data_Clean/Species_Checklists/hoverfly_checklist.csv")

# Jarrod Fowler oligolege list
jarrodFowler_oligoleges <- read.csv("Data_Raw/Interactions/fowler/jarrodFowler_scraped_interactions.csv") %>%
  rename(sourceTaxonName = higher,
         targetTaxonName = lower)

#### SDM flag ####
# Plants — record which tier each species has an SDM for
plant_sdm_names <- list.files("Data_Clean/SDMs/sdm_by_taxon/plants/continuous",
                              pattern = "\\.tif$", full.names = FALSE)
plant_sdm_names <- plant_sdm_names[!grepl("manifest", plant_sdm_names)] %>%
  gsub("\\_continuous.tif$", "", .) %>%
  gsub("_", " ", .)

# Pollinators
pollinator_folders <- c("bees", "butterflies", "moths", "hoverflies")
pollinator_sdm_names <- unlist(lapply(pollinator_folders, function(taxon) {
  files <- list.files(
    file.path("Data_Clean/SDMs/sdm_by_taxon", taxon, "continuous"),
    pattern = "\\.tif$", full.names = FALSE
  )
  files <- files[!grepl("manifest", files)]
  gsub("_", " ", gsub("\\_continuous.tif$", "", files))
}))

all_sdm_species <- c(plant_sdm_names, pollinator_sdm_names)

#### Phenometrics flag — record which method/tier was used ####
species_pheno_method <- phenology %>%
  filter(!is.na(phen_start) & !is.na(phen_end)) %>%
  dplyr::select(genus_species, phen_scale) %>%
  unique() %>%
  mutate(phenometrics_tier = case_when(
    phen_scale %in% c("PROV", "REGION", "JEPCODE", "CALIFORNIA") ~ "Tier1",
    phen_scale == "TYPE2" ~ "Tier2",
    phen_scale == "TYPE3" ~ "Tier3",
    TRUE ~ NA_character_
  )) %>%
  # If a species has multiple scales, keep the finest (Tier1 > Tier2 > Tier3)
  group_by(genus_species) %>%
  slice_min(order_by = phenometrics_tier, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(genus_species, phenometrics_tier)

#### Interaction network flag ####
species_in_network <- unique(c(
  predicted_interactions$sourceTaxonName_harm,
  predicted_interactions$targetTaxonName_harm
))

#### Regional occurrence — one column per JEPCODE ####
region_wide <- sp_by_region %>%
  dplyr::select(genus_species, JEPCODE) %>%
  distinct() %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = JEPCODE,
    values_from = present,
    values_fill = FALSE,
    names_prefix = "region_"
  )

#### Plant-specific attributes ####
plant_attrs <- plant_checklist_full %>%
  dplyr::select(scientificName, CRPR, State.Listing.Status, Federal.Listing.Status) %>%
  mutate(is_rare_plants = !is.na(CRPR) & CRPR != "") %>%
  rename(state_listing_plants  = State.Listing.Status,
         federal_listing_plants = Federal.Listing.Status) %>%
  unique()

#### Calscape attributes ####
calscape_attrs <- calscape_plant_data %>%
  filter(!is.na(genus_species_harm)) %>%
  mutate(
    # Available in at least one nursery
    available_in_nurseries = nurseries != "" & !is.na(nurseries),
    # Drought tolerant: low water need
    drought_tolerant = grepl("low", drought, ignore.case = TRUE)
  ) %>%
  dplyr::select(genus_species_harm, available_in_nurseries, drought_tolerant) %>%
  rename(genus_species = genus_species_harm) %>%
  # If species appears multiple times, take TRUE if any row is TRUE
  group_by(genus_species) %>%
  summarise(
    available_in_nurseries = any(available_in_nurseries, na.rm = TRUE),
    drought_tolerant       = any(drought_tolerant, na.rm = TRUE),
    .groups = "drop"
  )

#### Chaparral flag ####
chaparral_species <- calflora_communities %>%
  filter(community == "chaparral") %>%
  pull(genus_species_harm) %>%
  unique()

#### Oligolege flag for bees ####
# Harmonize Jarrod Fowler names to get CA oligolege list
oligolege_names   <- unique(jarrodFowler_oligoleges$sourceTaxonName)
names_olig_parsed <- bdc_clean_names(sci_names = oligolege_names, save_outputs = FALSE)

names_harm_olig <- harmonize_names(
  names           = unique(names_olig_parsed$names_clean),
  higher_tax      = "Hymenoptera",
  names_file      = "Temp/names_oligoleges.tsv",
  names_harm_file = "Temp/names_harmonized_oligoleges.tsv",
  gnv_path        = gnverifier_path
)

oligolege_species <- jarrodFowler_oligoleges %>%
  left_join(names_olig_parsed, by = join_by("sourceTaxonName" == "scientificName")) %>%
  left_join(names_harm_olig,   by = join_by("names_clean" == "ScientificName")) %>%
  filter(!is.na(names_harm)) %>%
  filter(names_harm %in% checklist$genus_species) %>%
  pull(names_harm) %>%
  unique()

#### Pollinator-specific attributes ####
pollinator_attrs <- bind_rows(
  bee_checklist       %>% dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA),
  butterfly_checklist %>% dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA),
  moth_checklist      %>% dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA),
  hoverfly_checklist  %>% dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA)
) %>%
  rename(genus_species = scientificName) %>%
  mutate(is_listed = (!is.na(ESA) & ESA != "") | (!is.na(CESA) & CESA != "")) %>%
  unique()

#### Assemble final supplementary table ####
# First filtering checklist to only one entry per genus_species (checklist has subspecies entries too)
checklist_filt <- checklist %>%
  # flag rows where scientificName matches genus_species (i.e. not a subspecies)
  mutate(is_species_level = scientificName == genus_species) %>%
  arrange(genus_species, desc(is_species_level)) %>%  # species-level rows first
  dplyr::select(genus_species, taxon, phylum, class, order, family, genus) %>%
  distinct(genus_species, .keep_all = TRUE) # Get first entry (preferred species level)

supplementary_checklist <- checklist_filt %>%
  dplyr::select(genus_species, taxon, phylum, class, order, family, genus) %>%
  unique() %>%
  ## Data availability flags
  mutate(
    has_SDM                = genus_species %in% all_sdm_species,
    in_predicted_interaction_network = genus_species %in% species_in_network
  ) %>%
  ## Phenometrics method
  left_join(species_pheno_method, by = "genus_species") %>%
  mutate(has_phenometrics = !is.na(phenometrics_tier)) %>%
  ## Regional occurrence
  left_join(region_wide, by = "genus_species") %>%
  ## Plant attributes
  left_join(plant_attrs, by = join_by("genus_species" == "scientificName")) %>%
  ## Calscape attributes
  left_join(calscape_attrs, by = "genus_species") %>%
  ## Chaparral flag
  mutate(in_chaparral_plants = genus_species %in% chaparral_species) %>%
  ## Pollinator attributes
  left_join(pollinator_attrs, by = "genus_species") %>%
  ## Oligolege flag (bees only)
  mutate(is_oligolege = genus_species %in% oligolege_species) %>%
  ## Set plant-only columns to NA for pollinators
  mutate(
    CRPR                   = ifelse(taxon != "plants", NA, CRPR),
    is_rare_plants         = ifelse(taxon != "plants", NA, is_rare_plants),
    state_listing_plants   = ifelse(taxon != "plants", NA, state_listing_plants),
    federal_listing_plants = ifelse(taxon != "plants", NA, federal_listing_plants),
    available_in_nurseries = ifelse(taxon != "plants", NA, available_in_nurseries),
    drought_tolerant       = ifelse(taxon != "plants", NA, drought_tolerant),
    in_chaparral_plants    = ifelse(taxon != "plants", NA, in_chaparral_plants)
  ) %>%
  ## Set pollinator-only columns to NA for plants
  mutate(
    Global.Rank  = ifelse(taxon == "plants", NA, Global.Rank),
    State.Rank   = ifelse(taxon == "plants", NA, State.Rank),
    ESA          = ifelse(taxon == "plants", NA, ESA),
    CESA         = ifelse(taxon == "plants", NA, CESA),
    is_listed    = ifelse(taxon == "plants", NA, is_listed),
    is_oligolege = ifelse(taxon != "bees",   NA, is_oligolege)
  ) %>%
  ## Sort
  arrange(taxon, genus_species)

#### Save ####
write.csv(supplementary_checklist, "Tables/Supplementary_Checklist.csv", row.names = FALSE)

