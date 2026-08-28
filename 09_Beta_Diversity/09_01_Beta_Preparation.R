#### Load Packages ####
library(dplyr)
library(bdc)
library(stringr)

gnverifier_path <- ""

#### Read in calflora community data ####
community_files <- list.files(path = "Data_Raw/Species_Traits_Attributes/calflora_plants_by_community", 
                              pattern = "*.csv", full.names = TRUE)

# Combine all community files
communities <- list()
for (file in community_files) {
  temp_data <- read.csv(file)
  temp_data$community <- tools::file_path_sans_ext(basename(file))
  communities[[length(communities) + 1]] <- temp_data
}
combined_communities <- bind_rows(communities)

#### Harmonize names ####
source("Code/99_Supporting/name_harmonization.R")

# Parse names
parse_names <- bdc_clean_names(sci_names = combined_communities$Taxon, save_outputs = FALSE)

# Examine errors
wrong_names <- filter(parse_names, is.na(names_clean) | quality == 0)

# Append cleaned names
parse_names <- parse_names %>%
  dplyr::select(.uncer_terms, names_clean)

calflora_clean <- dplyr::bind_cols(combined_communities, parse_names)

# Get unique plant names
names_plants <- calflora_clean %>% pull(names_clean) %>% unique()

# Harmonize names
names_harm_plants <- harmonize_names(
  names = names_plants,
  higher_tax = "Plantae",
  names_file = "Temp/names_calflora_communities.tsv",
  names_harm_file = "Temp/names_harmonized_calflora_communities.tsv",
  gnv_path = gnverifier_path
)

# Apply Jepson taxonomy
names_harm_plants <- apply_jepson_taxonomy(names_harm_plants)

# Join back with data
calflora_harmonized <- left_join(calflora_clean, names_harm_plants, by = c("names_clean" = "ScientificName"))

## Final dataset
Calflora_Communities <- calflora_harmonized %>%
  select(names_harm, community) %>%
  mutate(genus_species_harm = str_squish(word(names_harm, 1, 2))) %>%
  select(genus_species_harm, names_harm, community) %>%
  distinct()

write.csv(Calflora_Communities, "Data_Clean/Species_Traits_Attributes/Calflora_Communities.csv", row.names = FALSE)
