#### Load Packages ####
library(dplyr)
library(bdc)
library(sf)
library(stringr)

#### Read in data/supporting information ####
# name harmonization functions
source("Code/99_Supporting/name_harmonization.R")

# occurrence records
occurrences_clean <- read.csv("Temp/occurrences_clean.csv")
occurrences_clean_allCols <- read.csv("Temp/occurrences_clean_allColumns.csv")
occurrences_flagged <- read.csv("Temp/occurrences_flagged.csv")

# Invasive species
invasive_sp <- read.csv("Data_Raw/Species_Traits_Attributes/USRIISv2_MasterList.csv") %>%
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(scientificName, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean))) %>%
  # filter for only non-native taxa known to L48 states
  filter(locality == "L48")

# cisac invasives, CA (https://calinvasives.ucdavis.edu/)
all_invasives_CA <- read.csv("Data_Raw/Species_Traits_Attributes/cisac-species.csv") %>%
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(Scientific.Name, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean)))

# ipc invasive plants https://www.cal-ipc.org/plants/inventory/ 
plant_invasives_CA <- read.csv("Data_Raw/Species_Traits_Attributes/pafs.csv") %>%   
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(Latin.binomial, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean)))

# California listed threatened/endangered species
# I downloaded PDFs from 
# and extracted tables using https://tabula.technology
listed_sp_animals <- read.csv("Data_Raw/Species_Traits_Attributes/CNDDB_Special_Animals_List.csv") %>%
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(Scientific.Name, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean)))

# CNDDB list plants
listed_sp_plants <- read.csv("Data_Raw/Species_Traits_Attributes/CNDDB_Special_Plants_List.csv")  %>%
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(Scientific.Name, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean)))

# CNPS Rare Plant Inventory (RPI)
RPI_plants <- read.csv("Data_Raw/Species_Traits_Attributes/CNPS_RPI_2024-05-17.csv") %>%
  # remove extra bits in scientific name
  mutate(names_clean = stringr::str_remove_all(ScientificName, paste(c("ssp. ", "var. "), collapse = "|"))) %>%
  # and remove white space around name
  mutate(names_clean = gsub("\\s+", " ", trimws(names_clean)))

# Load Jepson names for use in checklist preparation below
jepson_names <- read.csv("Data_Raw/Species_Traits_Attributes/jepson_eflora_plants.csv") %>%
  dplyr::select(scientificNameClean, acceptedNameClean, status) %>%
  rename(names_harm_jepson = acceptedNameClean) %>%
  filter(!grepl(", in part", status),
         !status %in% c("Illegitimate name",
                        "Invalid name",
                        "Misapplied name",
                        "Noted name",
                        "Unabridged misapplied name")) %>%
  mutate(names_harm_jepson = ifelse(!is.na(status) & is.na(names_harm_jepson),
                                    scientificNameClean,
                                    names_harm_jepson))

## Prepare invasive species list
invasive_plants <- invasive_sp[invasive_sp$kingdom == "Plantae",]
invasive_animals <- invasive_sp[invasive_sp$kingdom == "Animalia",]

# CA only
CA_invasive_plant_names <- data.frame(names_clean = unique(c(all_invasives_CA[all_invasives_CA$Type == "plant",][["Scientific.Name"]], 
                                                             plant_invasives_CA$names_clean)))
CA_invasive_animal_names <- data.frame(names_clean = unique(all_invasives_CA[all_invasives_CA$Type %in% c("arthropod", "vertebrate", "invertebrate"),][["Scientific.Name"]]))


#### Harmonize names ####
# 1. Separate animals and plants due to possible duplication of names across Animalia and Plantae
names_plants <- unique(c(invasive_sp[invasive_sp$kingdom == "Plantae",][["scientificName"]],
                         all_invasives_CA[all_invasives_CA$Type == "plant",][["Scientific.Name"]],
                         plant_invasives_CA$Latin.binomial,
                         listed_sp_plants$Scientific.Name,
                         RPI_plants$ScientificName
))

names_animals <- unique(c(invasive_sp[invasive_sp$kingdom == "Animalia",][["scientificName"]],
                          all_invasives_CA[all_invasives_CA$Type %in% c("arthropod", "vertebrate", "invertebrate"),][["Scientific.Name"]],
                          listed_sp_animals$Scientific.Name
))

names_plants_parsed <- bdc_clean_names(sci_names = names_plants, save_outputs = FALSE)
names_animals_parsed <- bdc_clean_names(sci_names = names_animals, save_outputs = FALSE)

# 2. Harmonize names
names_harm_plants <- harmonize_names(
  names = unique(names_plants_parsed$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_supporting_plants.tsv",
  names_harm_file = "Temp/names_harmonized_supporting_plants.tsv"
)

names_harm_animals <- harmonize_names(
  names = unique(names_animals_parsed$names_clean),
  higher_tax = "Animalia",
  names_file = "Temp/names_supporting_animals.tsv",
  names_harm_file = "Temp/names_harmonized_supporting_animals.tsv"
)

# 3. Apply Jepson taxonomy to plants
names_harm_plants <- apply_jepson_taxonomy(names_harm_plants)
# Note: some plants not appearing in CA will have names_harm == NA. This is okay,
# since our occurrence dataset only includes names in Jepson anyway.

# Merge names with species lists
invasive_plants_harmonized <- left_join(invasive_plants, names_harm_plants, by = c("names_clean" = "ScientificName"))
invasive_animals_harmonized <- left_join(invasive_animals, names_harm_animals, by = c("names_clean" = "ScientificName"))
CA_invasive_plants_harmonized <- left_join(CA_invasive_plant_names, names_harm_plants, by = c("names_clean" = "ScientificName"))
CA_invasive_animals_harmonized <- left_join(CA_invasive_animal_names, names_harm_animals, by = c("names_clean" = "ScientificName"))
listed_plants_harmonized <- left_join(listed_sp_plants, names_harm_plants, by = c("names_clean" = "ScientificName"))
listed_animals_harmonized <- left_join(listed_sp_animals, names_harm_animals, by = c("names_clean" = "ScientificName"))
RPI_plants_harmonized <- left_join(RPI_plants, names_harm_plants, by = c("names_clean" = "ScientificName"))

## For invasive species, resort to species level (no ssp)
CA_invasive_plants_harmonized$names_harm <- word(CA_invasive_plants_harmonized$names_harm, 1, 2)
CA_invasive_animals_harmonized$names_harm <- word(CA_invasive_animals_harmonized$names_harm, 1, 2)

## select only relevant columns in supporting datasets 
invasive_plants_harmonized <- dplyr::select(invasive_plants_harmonized, names_harm, degreeOfEstablishment)
invasive_animals_harmonized <- dplyr::select(invasive_animals_harmonized, names_harm, degreeOfEstablishment)
listed_plants_harmonized <- dplyr::select(listed_plants_harmonized, names_harm, State.Listing.Status, Federal.Listing.Status, Heritage.Rank)
listed_animals_harmonized <- dplyr::select(listed_animals_harmonized, names_harm, Global.Rank, State.Rank, ESA, CESA)
RPI_plants_harmonized <- dplyr::select(RPI_plants_harmonized, names_harm, CRPR)

#### Prepare Checklists ####
#### Bees ####
# filter for only the chesshire data to create checklist
bee_data <- occurrences_clean[grepl("chesshire", occurrences_clean$databaseID, ignore.case = TRUE), ]

## Get unique species names
# Here, we check statuses at level of identification
bee_species_anyLevel <- bee_data %>%
  dplyr::select(scientificName, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Here, we downscale our observations to species level, so subspecies accrue any species level listings
bee_species_speciesLevel <- bee_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName, genus_species) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Coalesce the two data frames
bee_species <- bee_species_anyLevel %>%
  mutate(across(everything(), ~ coalesce(., bee_species_speciesLevel[[cur_column()]])))

## Filter out invasive species
# Note: this does not get rid of all non-native species, just the problematic invasives listed at https://calinvasives.ucdavis.edu/  and https://www.cal-ipc.org/plants/inventory/ 
bee_species = bee_species %>% 
  filter(is.na(CA_invasive)) %>%
  dplyr::select(-CA_invasive)

# Save as checklist
write.csv(bee_species, "Data_Clean/Species_Checklists/bee_checklist.csv")

### Plants ###
# filter for only the calflora data to create checklist
plant_data <- occurrences_clean[grepl("plants", occurrences_clean$taxon, ignore.case = TRUE), ]

## Get unique species names
# Here, we check statuses at level of identification
plant_species_anyLevel <- plant_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # get jepson status
  left_join(dplyr::select(jepson_names, c("scientificNameClean", "status")), by = join_by("scientificName" == "scientificNameClean"), multiple = "any") %>%
  # Get threatened status in CA and US
  left_join(listed_plants_harmonized, by = c("scientificName" = "names_harm"), multiple = "any") %>% # what species are listed in CA?
  # Get rare plant rank from CNPS RPI
  left_join(RPI_plants_harmonized, by = c("scientificName" = "names_harm"), multiple='any') %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_plants_harmonized, by = c("scientificName" = "names_harm"), multiple = "any") %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_plants_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, status, State.Listing.Status, Federal.Listing.Status, CRPR, Heritage.Rank, degreeOfEstablishment, CA_invasive,  kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Here, we downscale our observations to species level, so subspecies accrue any species level listings
plant_species_speciesLevel <- plant_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName, genus_species) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # get jepson status
  left_join(dplyr::select(jepson_names, c("scientificNameClean", "status")), by = join_by("genus_species" == "scientificNameClean"), multiple = "any") %>%
  # Get threatened status in CA and US
  left_join(listed_plants_harmonized, by = c("genus_species" = "names_harm"), multiple = "any") %>% # what species are listed in CA?
  # Get rare plant rank from CNPS RPI
  left_join(RPI_plants_harmonized, by = c("genus_species" = "names_harm"), multiple='any') %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_plants_harmonized, by = c("genus_species" = "names_harm"), multiple = "any") %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_plants_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, status, State.Listing.Status, Federal.Listing.Status, CRPR, Heritage.Rank, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Coalesce the two data frames
plant_species <- plant_species_anyLevel %>%
  mutate(across(everything(), ~ coalesce(., plant_species_speciesLevel[[cur_column()]])))

## Filter out invasive species
# Note: this does not get rid of all non-native species, just the problematic invasives listed at https://calinvasives.ucdavis.edu/  and https://www.cal-ipc.org/plants/inventory/ 
plant_species = plant_species %>% 
  filter(is.na(CA_invasive)) %>%
  dplyr::select(-CA_invasive)

## Filter to only species marked as "Native" in JEPSON, or appear in CNPS current plant list
# read cnps plant data, get native plants, clean names
CNPS_native_plants <- read.csv("Data_Raw/Species_Traits_Attributes/plants_table_5-7.csv") %>%
  filter(native_status %in% c("Native", "Native - Rare")) %>%
  pull(species) %>%
  bdc_clean_names(sci_names = ., save_outputs = FALSE)

plant_species = plant_species %>% 
  filter(status == "Native" | scientificName %in% CNPS_native_plants$names_clean) %>%
  dplyr::select(-status)

# Save as checklist
write.csv(plant_species, "Data_Clean/Species_Checklists/plant_checklist.csv")

#### Hoverflies ####
# filter for only the calflora data to create checklist
hoverfly_data <- occurrences_clean[grepl("hoverflies", occurrences_clean$taxon, ignore.case = TRUE), ]

## Get unique species names
# Here, we check statuses at level of identification
hoverfly_species_anyLevel <- hoverfly_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Here, we downscale our observations to species level, so subspecies accrue any species level listings
hoverfly_species_speciesLevel <- hoverfly_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName, genus_species) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Coalesce the two data frames
hoverfly_species <- hoverfly_species_anyLevel %>%
  mutate(across(everything(), ~ coalesce(., hoverfly_species_speciesLevel[[cur_column()]])))

## Filter out invasive species
# Note: this does not get rid of all non-native species, just the problematic invasives listed at https://calinvasives.ucdavis.edu/  and https://www.cal-ipc.org/plants/inventory/ 
hoverfly_species = hoverfly_species %>% 
  filter(is.na(CA_invasive)) %>%
  dplyr::select(-CA_invasive)

# Save as checklist
write.csv(hoverfly_species, "Data_Clean/Species_Checklists/hoverfly_checklist.csv")

#### Hummingbirds ####
# filter for only  data to create checklist
hummingbird_data <- occurrences_clean[grepl("hummingbirds", occurrences_clean$taxon, ignore.case = TRUE), ]

## Get unique species names
# Here, we check statuses at level of identification
hummingbird_species_anyLevel <- hummingbird_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Here, we downscale our observations to species level, so subspecies accrue any species level listings
hummingbird_species_speciesLevel <- hummingbird_data %>%
  dplyr::select(scientificName, genus_species, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(scientificName, genus_species) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(scientificName)) %>%
  arrange(scientificName) %>% # Order A-Z for output
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>% # keep only relevant columns
  distinct() # make sure no repeats after left join

# Coalesce the two data frames
hummingbird_species <- hummingbird_species_anyLevel %>%
  mutate(across(everything(), ~ coalesce(., hummingbird_species_speciesLevel[[cur_column()]])))

## Filter out invasive species
# Note: this does not get rid of all non-native species, just the problematic invasives listed at https://calinvasives.ucdavis.edu/  and https://www.cal-ipc.org/plants/inventory/ 
hummingbird_species = hummingbird_species %>% 
  filter(is.na(CA_invasive)) %>%
  dplyr::select(-CA_invasive)

# Save as checklist
write.csv(hummingbird_species, "Data_Clean/Species_Checklists/hummingbird_checklist.csv")

#### Leps ####
## Read in checklists from MPG and BAMONA
# MPG
MPG_checklist <- read.csv("Data_Raw/Species_Traits_Attributes/MPG_Checklist_CA.csv")$Genus_Species

# BAMONA
BAMONA_lines <- readLines("Data_Raw/Species_Traits_Attributes/BAMONA_CA.txt")

# Initialize an empty vector to store the processed strings
BAMONA_checklist <- character(length = length(BAMONA_lines))

# Iterate over each line of file
for (i in seq_along(BAMONA_lines)) {
  # Split the line into words
  words <- strsplit(BAMONA_lines[i], "\\s+")[[1]]
  
  # Get the first two words and join them with a space. this represents the genus and species name
  first_two_words <- trimws(paste(words[1], words[2], sep = " "))
  
  # Store scientific name in checklist vector
  BAMONA_checklist[i] <- first_two_words
}

# 1. get unique scientific names in checklist
lep_names <- unique(c(BAMONA_checklist, MPG_checklist))

# 2. Harmonize lep names
names_harm_leps <- harmonize_names(
  names = lep_names,
  higher_tax = "Animalia",
  names_file = "Temp/names_leps.tsv",
  names_harm_file = "Temp/names_harmonized_leps.tsv"
)

# 3. Update taxonomy
names_harm_leps <- update_taxonomy(names_harm_leps, tax = "animal")

## Manually fix issues in supporting datasets
listed_animals_harmonized <- listed_animals_harmonized %>%
  # monarch butterfly - in CNDDB list, listed as Danaus plexippus plexippus, but all observations at sp. level. need to match taxonomic resolution here
  mutate(names_harm = ifelse(names_harm == "Danaus plexippus plexippus", "Danaus plexippus", names_harm))

### Get unique species names
# Here, we check statuses at level of identification
lep_species_anyLevel <- names_harm_leps %>%
  dplyr::select(names_harm, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(names_harm) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(names_harm)) %>%
  rename(scientificName = names_harm) %>%
  # Order A-Z for output
  arrange(scientificName) %>%
  # how many parts of scientific name?
  mutate(num_parts_name = sapply(strsplit(trimws(scientificName), "\\s+"), length)) %>%
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("scientificName" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  # filter out if only 1 (genus/family)
  filter(num_parts_name > 1) %>%
  # remove previously created column
  dplyr::select(-num_parts_name) %>%
  # keep only relevant columns
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  distinct() # make sure no repeats after left join

# Here, we downscale our observations to species level, so subspecies accrue any species level listings
lep_species_speciesLevel <- names_harm_leps %>%
  dplyr::select(names_harm, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  group_by(names_harm) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(names_harm)) %>%
  rename(scientificName = names_harm) %>%
  # Order A-Z for output
  arrange(scientificName) %>%
  # how many parts of scientific name?
  mutate(num_parts_name = sapply(strsplit(trimws(scientificName), "\\s+"), length)) %>%
  # Get only genus, species parts of name
  mutate(genus_species = sapply(strsplit(scientificName, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  # Get threatened status in CA and US
  left_join(listed_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # Is the species invasive/introduced?
  left_join(invasive_animals_harmonized, by = c("genus_species" = "names_harm")) %>% # what species are listed in CA?
  # invasive using CA lists?
  mutate(CA_invasive = if_else(scientificName %in% CA_invasive_animals_harmonized$names_harm, 1, NA_real_)) %>%
  # filter out if only 1 (genus/family)
  filter(num_parts_name > 1) %>%
  # remove previously created column
  dplyr::select(-num_parts_name) %>%
  # keep only relevant columns
  dplyr::select(scientificName, Global.Rank, State.Rank, ESA, CESA, degreeOfEstablishment, CA_invasive, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  distinct() # make sure no repeats after left join

# Coalesce the two data frames
lep_species <- lep_species_anyLevel %>%
  mutate(across(everything(), ~ coalesce(., lep_species_speciesLevel[[cur_column()]])))

## Filter out invasive species
# Note: this does not get rid of all non-native species, just the problematic invasives listed at https://calinvasives.ucdavis.edu/  and https://www.cal-ipc.org/plants/inventory/ 
lep_species = lep_species %>% 
  filter(is.na(CA_invasive)) %>%
  dplyr::select(-CA_invasive)

### Seperate butterflies and moths
butterfly_species <- filter(lep_species, family %in% c("Papilionidae", "Nymphalidae", "Pieridae", "Lycaenidae", "Riodinidae", "Hesperiidae", "Heliconiidae"))

moth_species <- filter(lep_species, !family %in% c("Papilionidae", "Nymphalidae", "Pieridae", "Lycaenidae", "Riodinidae", "Hesperiidae", "Heliconiidae"))

# Save as checklist
write.csv(butterfly_species, "Data_Clean/Species_Checklists/butterfly_checklist.csv")
write.csv(moth_species, "Data_Clean/Species_Checklists/moth_checklist.csv")

#### Now, let's save all checklists in one!####
### Load species checklists
# Bees
a = read.csv("Data_Clean/Species_Checklists/bee_checklist.csv")
a$taxon = "bees"

# Hoverflies
b = read.csv("Data_Clean/Species_Checklists/hoverfly_checklist.csv")
b$taxon = "hoverflies"

# Leps
c1 = read.csv("Data_Clean/Species_Checklists/butterfly_checklist.csv")
c1$taxon = "butterflies"
c2 = read.csv("Data_Clean/Species_Checklists/moth_checklist.csv")
c2$taxon = "moths"

# Hummingbirds
d = read.csv("Data_Clean/Species_Checklists/hummingbird_checklist.csv")
d$taxon = "hummingbirds"

# Plants
e = read.csv("Data_Clean/Species_Checklists/plant_checklist.csv")
e$taxon = "plants"

#compile into 1
a = a[c(2,19,9:18)]
b = b[c(2,19,9:18)]
c1 = c1[c(2,19,9:18)]
c2 = c2[c(2,19,9:18)]
d = d[c(2,19,9:18)]
e = e[c(2,19,9:18)]
checklist = rbind(a,b,c1,c2,d,e)
checklist = unique(checklist)

# add genus_species column
checklist$genus_species = word(checklist$scientificName, 1, 2)

#Save output
write.csv(checklist, "Data_Clean/Species_Checklists/checklist_cleaned.csv", row.names = F)
write.csv(checklist, "Data_Calscape/Species_Checklists/checklist_cleaned.csv", row.names = F)

#### Update occurrences dataset, filtered and taxon group update (for leps) ####
# filter out non checklist species
not_in_checklist <- occurrences_clean %>%
  filter(!genus_species %in% checklist$genus_species)

# Add removed records to flagged file
occurrences_flagged <- bind_rows(occurrences_flagged, not_in_checklist)

# Keep only checklist species
occurrences_clean <- occurrences_clean %>%
  filter(genus_species %in% checklist$genus_species)

occurrences_clean_allCols <- occurrences_clean_allCols %>%
  filter(genus_species %in% checklist$genus_species)

# Save outputs
write.csv(occurrences_clean, "Data_Clean/Occurrences/occurrences_clean.csv", row.names = F)
write.csv(occurrences_clean_allCols, "Data_Clean/Occurrences/occurrences_clean_allColumns.csv", row.names = F)
write.csv(occurrences_flagged, "Data_Clean/Occurrences/occurrences_flagged.csv", row.names = F)


#### Create species x region dataset -- ecoregion checklist ####
# Read jepcode delineations shapefile
ecoregions <- st_read("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid() %>%
  st_transform(4326)

occurrences_sf <- st_as_sf(occurrences_clean, 
                           coords = c("decimalLongitude", "decimalLatitude"), 
                           crs = 4326,
                           remove = FALSE)  # Keep original lat/lon columns

# now get unique combinations of species+region
sp_by_region <- occurrences_sf %>%
  st_join(ecoregions) %>%
  st_drop_geometry() %>%
  filter(!is.na(JEPCODE)) %>%
  dplyr::select(scientificName, genus_species, JEPCODE, REGION, PROV) %>%
  unique() %>%
  arrange(scientificName, genus_species, PROV, REGION, JEPCODE)

write.csv(sp_by_region, "Data_Clean/Species_Checklists/ecoregion_checklist.csv", row.names = F)

# # For exploring which genera have duplicate species names - can help identify synonymous species
# # If species identified as synonyms, go back to taxonomy_corrections.csv, add, and re-run from step 02_02
# df <- checklist %>%
#   select(genus, specificEpithet, family) %>%
#   distinct()
# 
# # Self-join to find genus pairs sharing species names
# genera_with_dupes <- df %>%
#   inner_join(df, by = "specificEpithet", relationship = "many-to-many") %>%
#   filter(genus.x < genus.y) %>%  # Avoid duplicates and self-matches
#   filter(family.x == family.y) %>%  # Only keep pairs in the same family
#   group_by(genus.x, genus.y) %>%
#   summarise(
#     family = first(family.x),  # Include the family name
#     n_shared_species = n(),
#     shared_species = paste(sort(unique(specificEpithet)), collapse = ", "),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_shared_species))
# 
# # Rename columns for clarity
# genera_with_dupes <- genera_with_dupes %>%
#   rename(
#     genus_1 = genus.x,
#     genus_2 = genus.y
#   )
