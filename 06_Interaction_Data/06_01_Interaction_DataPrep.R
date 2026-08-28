#### 6: Interaction Data ####

# Clear workspace
rm(list = ls(all = TRUE))

#libraries
library(dplyr)
library(data.table)
library(bdc)
library(BeeBDC)
library(tidyr)

# Gnverifier for taxonomic harm
source("Code/99_Supporting/name_harmonization.R")
# Define gnverifier path
gnverifier_path <- ""

##### GloBI #####

# download "interactions.tsv.gz" from https://www.globalbioticinteractions.org/data

# Last downloaded Nov 18, 2025

# Extract the desired data
script_lines <- c(
  "#!/bin/bash",
  "",
  "DIR=\"Data_Raw/Interactions/globi\"",
  "cd \"$DIR\" || { echo \"Failed to cd to $DIR\"; exit 1; }",
  "",
  "echo \"Processing large file with parallel decompression...\"",
  "",
  "# Use pigz for parallel decompression (much faster than gzcat)",
  "# Word boundaries ensure exact matches",
  "pattern='\\<(Andrenidae|Apidae|Colletidae|Halictidae|Megachilidae|Melittidae|Stenotritidae|Syrphidae|Lepidoptera|Trochilidae)\\>'",
  "",
  "outfile=\"globi_filtered_unique.tsv\"",
  "",
  "echo \"Step 1: Extracting header...\"",
  "pigz -dc -p 4 interactions.tsv.gz | head -n 1 > \"$outfile\"",
  "",
  "echo \"Step 2: Filtering data (this will take several minutes)...\"",
  "pigz -dc -p 4 interactions.tsv.gz | tail -n +2 | grep -E \"$pattern\" | sort -u -S 2G --parallel=4 >> \"$outfile\"",
  "",
  "echo \"Done!\"",
  "echo \"Total lines (including header):\" $(wc -l < \"$outfile\")",
  "echo \"File size:\" $(ls -lh \"$outfile\" | awk '{print $5}')"
)

script_path <- tempfile(fileext = ".sh")
writeLines(script_lines, script_path)
Sys.chmod(script_path, mode = "755")
system(paste("bash", shQuote(script_path)))

#read in the GBIF data
globi <- fread(
  'Data_Raw/Interactions/globi/globi_filtered_unique.tsv',
  sep = "\t",
  quote = ""
)

#filter just for Plantae in the target taxon 
globi = filter(globi, targetTaxonKingdomName == "Plantae")

#fix common issue with species names
globi$sourceTaxonName <- ifelse(grepl("^BOLD", globi$sourceTaxonName), paste(globi$sourceTaxonGenusName, globi$sourceTaxonSpeciesName, sep = " "), globi$sourceTaxonName)

globi$sourceTaxonName = trimws(globi$sourceTaxonName, which = "both")

globi$targetTaxonName <- ifelse(grepl("^BOLD", globi$targetTaxonName), paste(globi$targetTaxonGenusName, globi$targetTaxonSpeciesName, sep = " "), globi$targetTaxonName)

globi$targetTaxonName = trimws(globi$targetTaxonName, which = "both")

# Create tally for each taxon (checking appropriate taxonomic rank)
taxa_tally <- tibble(
  Taxon = c("Andrenidae", "Apidae", "Colletidae", "Halictidae", 
            "Megachilidae", "Melittidae", "Stenotritidae", 
            "Syrphidae", "Trochilidae", "Lepidoptera"),
  Count = c(
    # Bee families
    sum(globi$sourceTaxonFamilyName == "Andrenidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Apidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Colletidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Halictidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Megachilidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Melittidae", na.rm = TRUE),
    sum(globi$sourceTaxonFamilyName == "Stenotritidae", na.rm = TRUE),
    # Syrphidae (hover flies)
    sum(globi$sourceTaxonFamilyName == "Syrphidae", na.rm = TRUE),
    # Trochilidae (hummingbirds)
    sum(globi$sourceTaxonFamilyName == "Trochilidae", na.rm = TRUE),
    # Lepidoptera (butterflies & moths - order level)
    sum(globi$sourceTaxonOrderName == "Lepidoptera", na.rm = TRUE)
  )
)

#take only needed columns and remove duplicates
globi = globi %>%
  select(
    sourceTaxonName,
    sourceTaxonRank,
    sourceTaxonSpeciesName,
    sourceTaxonSubgenusName,
    sourceTaxonGenusName,
    sourceTaxonFamilyName,
    sourceTaxonOrderName,
    sourceTaxonClassName,
    sourceTaxonPhylumName,
    sourceTaxonKingdomName,
    sourceLifeStageName,
    interactionTypeName,
    targetTaxonName,
    targetTaxonRank,
    targetTaxonSpeciesName,
    targetTaxonSubgenusName,
    targetTaxonGenusName,
    targetTaxonFamilyName,
    targetTaxonOrderName,
    targetTaxonClassName,
    targetTaxonPhylumName,
    targetTaxonKingdomName,
    targetLifeStageName,
    decimalLatitude,
    decimalLongitude,
    localityName
  ) %>%
  distinct()

#### Katja's Bee Interaction Data ####

#Redundant with the full globi data above, uses same process to get bee data

##### HOSTS #####

# Download from https://data.nhm.ac.uk/dataset/hosts 
# Cite HOSTS: Gaden S. Robinson; Phillip R. Ackery; Ian Kitching; George W Beccaloni; Luis M. Hernández (2023). HOSTS - a Database of the World's Lepidopteran Hostplants [Data set]. Natural History Museum. Downloaded November 18, 2025. https://doi.org/10.5519/qd.bsucrxdz. 

# Last downloaded Nov 18, 2025

# Read in
hosts = read.csv('Data_Raw/Interactions/hosts/resource.csv')

#make new name columns with any information we have
hosts <- hosts %>%
  mutate(
    targetTaxonName = case_when(
      !is.na(Hostplant.Species) & Hostplant.Species != "" &
        !is.na(Hostplant.Subspecies.var) & Hostplant.Subspecies.var != "" ~
        paste(Hostplant.Genus, Hostplant.Species, Hostplant.Subspecies.var),
      
      !is.na(Hostplant.Species) & Hostplant.Species != "" ~
        paste(Hostplant.Genus, Hostplant.Species),
      
      !is.na(Hostplant.Genus) & Hostplant.Genus != "" ~
        Hostplant.Genus,
      
      TRUE ~ Hostplant.Family
    ),
    
    sourceTaxonName = case_when(
      !is.na(Insect.Species) & Insect.Species != "" &
        !is.na(Insect.Subspecies) & Insect.Subspecies != "" ~
        paste(Insect.Genus, Insect.Species, Insect.Subspecies),
      
      !is.na(Insect.Species) & Insect.Species != "" ~
        paste(Insect.Genus, Insect.Species),
      
      !is.na(Insect.Genus) & Insect.Genus != "" ~
        Insect.Genus,
      
      TRUE ~ Insect.Family
    )
  )

#Add interaction type column
hosts$interactionTypeName = "hasHost"

#Add authorship to insect names
hosts$sourceTaxonName <- ifelse(
  is.na(hosts$Insect.Author) | hosts$Insect.Author == "",
  hosts$sourceTaxonName,
  paste0(hosts$sourceTaxonName, " (", hosts$Insect.Author, ")")
)

#rename to match GloBI, take needed column, de-duplicate
hosts = hosts %>%
  rename(
    sourceTaxonFamilyName = Insect.Family,
    sourceTaxonGenusName = Insect.Genus,
    targetTaxonFamilyName = Hostplant.Family,
    targetTaxonGenusName = Hostplant.Genus, 
    localityName = Location
  ) %>%
  select(
    sourceTaxonName,
    sourceTaxonFamilyName,
    sourceTaxonGenusName,
    targetTaxonName,
    targetTaxonFamilyName,
    targetTaxonGenusName,
    localityName,
    interactionTypeName
  ) %>%
  distinct()

##### Caldwell #####

# Raw data from Caldwell 2021 "California Plants as Resources for Lepidoptera"

# Digitized and cleaned by Chris Cosma

# Read in
caldwell = read.csv('Data_Raw/Interactions/caldwell/int_manual.csv')

#rename to match globi, de-duplicate
caldwell <- caldwell %>%
  rename(
    targetTaxonName  = lower,
    sourceTaxonName  = higher,
    interactionTypeName = interaction
  ) %>%
  mutate(
    interactionTypeName = case_when(
      interactionTypeName == "herbivory"   ~ "hasHost",
      interactionTypeName == "pollination" ~ "visitsFlowersOf",
      TRUE ~ interactionTypeName
    )
  ) %>%
  select(-type) %>%
  distinct()

#add locality name
caldwell$localityName = "California, USA"

##### Discover Life #####

#Here is the raw data from discover life (in species_records_raw). There is an rds file for each species, and in each rds file is a list of lists. Each of the higher level list items is an occurrence record with an observed interaction from our checklists, and the lower level list are just like columns of info for that record. All that needs to be done is the records merged into a normal, dataframe type of format

# and here is the script to get to this point in case
# https://github.com/MorphoNativePlants/interaction_scraping

# Each records started with a search of species from our pollinators checklist

# List all RDS files in the folder
rds_files <- list.files(path = "Data_Raw/Interactions/discover_life/species_records_raw", pattern = "\\.rds$", full.names = T)

# Initialize an empty list to store data from each RDS file
records_list <- list()

# Loop through each RDS file
for (file in rds_files) {
  
  # Read the RDS file
  species_data <- readRDS(file)
  
  # Flatten the list of lists into a data frame
  species_df <- do.call(bind_rows, lapply(species_data, as.data.frame))
  
  #prepare and add the record title (species) as a column
  file_name <- basename(file)
  file_name <- gsub("\\+", " ", file_name)
  file_name <- gsub('"|\\.rds', "", file_name)
  species_df$recordScientificName <- file_name
  
  # Add the data frame to the list
  records_list[[file]] <- species_df
}

# Combine all data frames into one
combined_records <- do.call(bind_rows, records_list)

#Clean records

dl_cleaning = combined_records %>%
  select(
    
    recordScientificName, #pollinator species
    
    family, Family, # pollinator family
    
    order, Order, ORDER, #pollinator Order
    
    plant.host.prey, host.prey, associatedTaxa, Host.name, host_plant, Plant.host, FloralSource, HOST.PLANT, host.name, host.prey.associate, # plant/host species
    
    Host.Family, Family.of.plant.host, host.family, #plant/host family
    
    latitude, decimalLatitude, Latitude, Digital.latitude, LATITUDE, #latitude 
    
    longitude, decimalLongitude, Longitude, Digital.latitude, Digital.longitude, LONGITUDE, #longitude
    
    country, Country, COUNTRY, CountryName,  #country
    
    StateName, state, STATE, State, state.province, stateProvince, State.or.province, State.Province,  #state/province
    
    county, COUNTYNAME, County, County.other  #county
    
  ) %>%
  
  mutate(
    higherClean = recordScientificName,
    higherFamilyClean = coalesce(family, Family),
    higherOrderClean  = coalesce(order, Order, ORDER),
    lowerClean = coalesce(plant.host.prey, host.prey, associatedTaxa, Host.name, host_plant, Plant.host, FloralSource, HOST.PLANT, host.name, host.prey.associate),
    lowerFamilyClean = coalesce(Host.Family, Family.of.plant.host, host.family),
    latitudeClean = coalesce(latitude, decimalLatitude, Latitude, Digital.latitude, LATITUDE),
    longitudeClean = coalesce(longitude, decimalLongitude, Longitude, Digital.latitude, Digital.longitude, LONGITUDE),
    countryClean = coalesce(country, Country, COUNTRY, CountryName),
    stateClean = coalesce(StateName, state, STATE, State, state.province, stateProvince, State.or.province, State.Province),
    countyClean = coalesce(county, COUNTYNAME, County, County.other)
  ) %>%
  
  select(higherOrderClean, higherFamilyClean, higherClean, lowerFamilyClean, lowerClean, latitudeClean, longitudeClean, countryClean, stateClean, countyClean)

#combine country, state, county into localityName column
dl_cleaning <- dl_cleaning %>%
  unite(
    localityName,
    countyClean, stateClean, countryClean,
    sep = ", ",
    na.rm = TRUE
  ) 

names(dl_cleaning)

# rename to match globu and take unique
dl = dl_cleaning %>%
  rename(
    sourceTaxonOrderName = higherOrderClean,
    sourceTaxonFamilyName = higherFamilyClean,
    sourceTaxonName = higherClean,
    targetTaxonFamilyName = lowerFamilyClean,
    targetTaxonName = lowerClean,
    decimalLatitude = latitudeClean,
    decimalLongitude = longitudeClean
  ) %>%
  distinct()

# Degrees.latitude, Minutes.latitude, Seconds.latitude,

# Degrees.longitude, Minutes.longitude, Seconds.longitude,

#date.yyyymmdd.hr.mn, date1.yyyymmdd, year, day, month, verbatimEventDate, date.yyyymmdd, CollectionYear, CollectionMonth, CollectionDate, date, MONTH, DATE, YEAR, DATE.DD.MM.YYYY # date information

#record_id, id, basisOfRecord, collectionCode, collectionID, type, datasetName, ownerInstitutionCode, institutionCode, #record identifiers

#### Big Book of Hymenoptera ####
bigbook = read.csv('Data_Raw/Interactions/bigbook/bigBook_scraped_interactions.csv')

bigbook = bigbook %>%
  rename(
    sourceTaxonName = higher,
    targetTaxonName = lower
  ) %>%
  select(-X) %>%
  distinct()

#for now leaving collectsPollenOf as separate interaction type

#### Jared Fowler Pollen Specialist Bees ####
fowler = read.csv('Data_Raw/Interactions/fowler/jarrodFowler_scraped_interactions.csv')

fowler = fowler %>%
  rename(
    sourceTaxonName = higher,
    targetTaxonName = lower
  ) %>%
  distinct()

#### Simon Doneski Butterfly Hostplants ####
simon = read_excel('Data_Raw/Interactions/simon/UNR_UNM hostDatabase.xlsx')

## Select only columns we need
simon = simon[,c("butterflyFam", "SpScott1986", "SpNABA2.6", "spBOA", "spPelhamCatalog", "plantFamUpdated", "plantGenus", "plantSp", "subSpOrVar1", "subSpOrVar2", "Scientific.Name.Plant", "Location", "Scientific.Name.with.Author")]

#convert to NA for coalesce
simon[simon == ""] <- NA

#coalesce for new targetTaxonName column, rename, de-duplicate
simon <- simon %>%
  mutate(
    sourceTaxonName = coalesce(SpScott1986, SpNABA2.6, spBOA, spPelhamCatalog)
  ) %>%
  rename(
    sourceTaxonFamilyName = butterflyFam,
    targetTaxonFamilyName = plantFamUpdated,
    targetTaxonGenusName = plantGenus,
    targetTaxonName = Scientific.Name.Plant,
    localityName = Location
  ) %>%
  select (
    sourceTaxonName,
    sourceTaxonFamilyName,
    targetTaxonFamilyName,
    targetTaxonGenusName,
    targetTaxonName,
    localityName
  ) %>%
  distinct()

#add interactionTypeName
simon$interactionTypeName = "hasHost"

#### Shropshire + Tallamy ####

### 1. Read in raw data
# Read in Shropshire-Tallamy data from https://zenodo.org/records/17747343
ShropshireTallamy <- read_excel('Data_Raw/Interactions/ShropshireTallamy2025/oo_1478226.xlsx', 
                                sheet = "UDELep May 2025", skip = 1)

### 2. Parse and clean Shropshire-Tallamy lepidopteran names
# Parse lep names using bdc
parse_names <- bdc_clean_names(
  sci_names = ShropshireTallamy$SPECIES_COMPLETE_BINOMIAL_NOMENCLATURE, 
  save_outputs = FALSE
) %>% unique()

# Join parsed names back
ShropshireTallamy <- left_join(
  ShropshireTallamy, parse_names,
  by = join_by("SPECIES_COMPLETE_BINOMIAL_NOMENCLATURE" == "scientificName")
)

### 3. Function to expand a column split on whitespace (as in host families and genera cols)
expand_host_col <- function(df, col, new_col) {
  df %>%
    mutate("{new_col}" := str_split(.data[[col]], "\\s+")) %>%
    unnest(!!sym(new_col)) %>%
    filter(!is.na(.data[[new_col]]), 
           .data[[new_col]] != "", 
           .data[[new_col]] != "NA",
           .data[[new_col]] != ".none_specified.")
}

### 4. Expand and clean host genera
# Split host genera col, clean bracket annotations and synonyms
ShropshireTallamy_long_genera <- ShropshireTallamy %>%
  expand_host_col("HOST_GENERA", "HOST_GENERA_split") %>%
  mutate(
    # extract notes from brackets e.g. [plant], [=Synonym]
    HOST_GENUS_NOTE  = str_extract(HOST_GENERA_split, "(?<=\\[).*?(?=\\])"),
    # only keep note if it doesn't start with = or &
    HOST_GENUS_NOTE  = if_else(str_detect(HOST_GENUS_NOTE, "^[=&]"), NA_character_, HOST_GENUS_NOTE),
    # remove bracket content, synonyms, and stray brackets
    HOST_GENUS_CLEAN = str_remove_all(HOST_GENERA_split, "\\[.*?\\]|=\\w+|[\\[\\]]"),
    NON_PLANT_FLAG   = !is.na(HOST_GENUS_NOTE)
  ) %>%
  # remove non-plant hosts and non-butterfly records
  filter(HOST_GENUS_CLEAN != "",
         HOST_GENUS_NOTE %in% c(NA, "plant")) %>%
  rename(sourceTaxonName   = names_clean,
         targetTaxonName   = HOST_GENUS_CLEAN) %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, targetTaxonName, HOST_FAMILIES)

### 5. Harmonize plant genera via GNVerifier
# Get unique plant names for harmonization
names_plants <- ShropshireTallamy_long_genera %>% pull(targetTaxonName) %>% unique()

# Harmonize names
names_harm_plants <- harmonize_names(
  names           = names_plants,
  higher_tax      = "Plantae",
  names_file      = "Temp/names_interactions_plants.tsv",
  names_harm_file = "Temp/names_harmonized_interactions_plants.tsv",
  gnv_path = gnverifier_path
)
names_harm_plants_withHigher <- update_taxonomy(names_harm_plants, tax = "plant")

# Join harmonized names back to get inferred family per genus
ShropshireTallamy_long_genera <- left_join(
  ShropshireTallamy_long_genera, names_harm_plants_withHigher,
  by = join_by("targetTaxonName" == "ScientificName")
) %>%
  rename(plantFamInferred  = family,
         originalNamePlant = targetTaxonName) %>%
  mutate(plantGenus = str_extract(names_harm, "^\\S+")) %>%
  dplyr::select(-names_harm)

### 6. Expand host families and find those not covered by genera
# Split host families col, clean bracket annotations
ShropshireTallamy_long_fams <- ShropshireTallamy %>%
  expand_host_col("HOST_FAMILIES", "HOST_FAMILY_split") %>%
  mutate(
    HOST_FAMILY_CLEAN = str_remove_all(HOST_FAMILY_split, "\\[.*?\\]|=\\w+|[\\[\\]]")
  ) %>%
  # remove non-butterfly records
  filter(HOST_FAMILY_CLEAN != "") %>%
  rename(sourceTaxonName = names_clean) %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, HOST_FAMILY_CLEAN)

# Get families already represented by genus-level data
fams_from_genera <- ShropshireTallamy_long_genera %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, plantFamInferred) %>%
  filter(!is.na(plantFamInferred)) %>%
  distinct()

# Keep only host family entries not already covered by a genus row
fam_only_rows <- ShropshireTallamy_long_fams %>%
  anti_join(fams_from_genera, 
            by = join_by("sourceTaxonName", 
                         "LEP_FAMILY", 
                         "HOST_FAMILY_CLEAN" == "plantFamInferred")) %>%
  mutate(plantGenus        = NA_character_,
         plantFamInferred  = HOST_FAMILY_CLEAN,
         originalNamePlant = NA_character_) %>%
  dplyr::select(-HOST_FAMILY_CLEAN)

# Clean a single plant name from references
clean_ref_name <- function(x) {
  x %>%
    # remove bracket content entirely (synonyms and notes)
    str_remove_all("\\[.*?\\]") %>%
    # remove "sp." to reduce to genus level
    str_remove_all("\\s+sp\\.?$") %>%
    str_squish()
}

# Parse HOST_references into long format with one row per butterfly x citation x plant name
refs_parsed <- ShropshireTallamy %>%
  filter(!is.na(HOST_references)) %>%
  rename(sourceTaxonName = names_clean) %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, HOST_references) %>%
  mutate(sourceTaxonName = str_extract(sourceTaxonName, "^\\S+\\s+\\S+")) %>%
  rename(HOST_references_full = HOST_references) %>%
  mutate(citation_chunk = str_split(HOST_references_full, ";\\s*")) %>%
  unnest(citation_chunk) %>%
  filter(citation_chunk != "") %>%
  mutate(
    # remove bracket notes about the citation itself e.g. [for Achyra similalis Guenee]
    # and captive reared notes e.g. & [captive reared] (...)
    citation_chunk = str_remove_all(citation_chunk, "\\[for[^\\]]*\\]"),
    citation_chunk = str_remove_all(citation_chunk, "&\\s*\\[[^\\]]*\\]\\s*\\([^)]*\\)"),
    citation       = str_trim(str_extract(citation_chunk, "^[^(]+")),
    # extract only the FIRST parenthetical to avoid greedy matching issues
    plants_raw     = str_extract(citation_chunk, "(?<=\\()[^)]+(?=\\))")
  ) %>%
  filter(!is.na(plants_raw)) %>%
  mutate(
    # remove qualifying prefixes before colon e.g. "Rosaceae host plants: " or "Eurasian records: "
    plants_raw = str_remove(plants_raw, "^[^,)]+:\\s*")
  ) %>%
  mutate(plantName_raw = str_split(plants_raw, ",\\s*|\\s*&\\s*")) %>%
  unnest(plantName_raw) %>%
  mutate(
    plantName_raw   = str_remove(plantName_raw, "^\\w+\\s+(?=[A-Z])"),
    plantName_clean = clean_ref_name(plantName_raw),    n_words          = str_count(plantName_clean, "\\S+"),
    plantGenus_ref   = str_extract(plantName_clean, "^\\S+"),
    plantSp_ref      = if_else(n_words >= 2, word(plantName_clean, 2), NA_character_),
    plantName_norank = str_remove_all(plantName_clean, "\\b(var|ssp|subsp|f)\\b\\.?\\s*"),
    n_words_norank   = str_count(str_squish(plantName_norank), "\\S+"),
    subSpOrVar1      = if_else(n_words_norank >= 3, word(str_squish(plantName_norank), 3), NA_character_),
    subSpOrVar2      = if_else(n_words_norank >= 4, word(str_squish(plantName_norank), 4), NA_character_),
    is_family_ref    = str_detect(plantName_clean, "aceae$|idae$"),
    is_species_ref   = n_words >= 2
  ) %>%
  filter(
    plantName_clean != "",
    # remove individual names that are family-level
    !is_family_ref,
    # genus must start with a capital letter
    str_detect(plantGenus_ref, "^[A-Z]")
  ) %>%
  dplyr::select(-n_words, -plantName_norank, -n_words_norank) %>%
  distinct()

### 6. Combine genera rows and family-only rows
ShropshireTallamy_cleanPlants <- bind_rows(
  ShropshireTallamy_long_genera %>% dplyr::select(-HOST_FAMILIES),
  fam_only_rows
) %>%
  rename(plantFamUpdated = plantFamInferred)

### Join species-level reference data
# For each genus-level row, find matching species from refs_parsed
# A match is when plantGenus_ref == plantGenus for the same butterfly species
species_from_refs <- ShropshireTallamy_cleanPlants %>%
  inner_join(
    refs_parsed %>% filter(is_species_ref) %>%
      dplyr::select(sourceTaxonName, LEP_FAMILY, plantGenus_ref,
                    plantName_clean, plantSp_ref, subSpOrVar1, subSpOrVar2,
                    citation, HOST_references_full),
    by = join_by("sourceTaxonName", "LEP_FAMILY", "plantGenus" == "plantGenus_ref")
  ) %>%
  rename(plantSp = plantSp_ref,
         citation_ref = citation) %>%
  # drop full binomial since we now have it split into cols
  dplyr::select(-plantName_clean)

# genus-level rows that have at least one species match get replaced by species rows
genera_with_species <- species_from_refs %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, plantGenus) %>%
  distinct()

# keep genus rows only when no species were found for that genus
genus_only_rows <- ShropshireTallamy_cleanPlants %>%
  anti_join(genera_with_species,
            by = join_by("sourceTaxonName", "LEP_FAMILY", "plantGenus")) %>%
  mutate(plantSp     = NA_character_,
         subSpOrVar1  = NA_character_,
         subSpOrVar2  = NA_character_,
         citation_ref      = NA_character_,
         HOST_references_full = NA_character_)

# assign citation to genus-only and family-only rows when only one citation 
# exists for that butterfly AND that citation mentions the same genus
single_citation <- refs_parsed %>%
  dplyr::select(sourceTaxonName, LEP_FAMILY, citation, HOST_references_full, plantGenus_ref) %>%
  distinct() %>%
  group_by(sourceTaxonName, LEP_FAMILY, plantGenus_ref) %>%
  filter(n_distinct(citation) == 1) %>%
  slice(1) %>%
  ungroup() %>%
  rename(citation_ref = citation)

genus_only_rows <- genus_only_rows %>%
  dplyr::select(-citation_ref, -HOST_references_full) %>%
  left_join(single_citation %>% 
              dplyr::select(sourceTaxonName, LEP_FAMILY, 
                            plantGenus_ref, citation_ref, HOST_references_full),
            by = join_by("sourceTaxonName", "LEP_FAMILY", "plantGenus" == "plantGenus_ref"))

### 8. Combine species rows and genus-only rows
ShropshireTallamy_cleanPlants_withSp <- bind_rows(
  species_from_refs,
  genus_only_rows
)

### 9. Formatting
shropTall = ShropshireTallamy_cleanPlants_withSp %>%
  mutate(targetTaxonName = if_else(is.na(plantSp), plantGenus, paste(plantGenus, plantSp)),
         interactionTypeName = "hasHost") %>%
  rename(sourceTaxonFamilyName = LEP_FAMILY,
         targetTaxonFamilyName = plantFamUpdated) %>%
  dplyr::select(sourceTaxonName, sourceTaxonFamilyName, targetTaxonName, targetTaxonFamilyName, interactionTypeName) %>%
  # Remove family level records for now
  filter(!is.na(sourceTaxonName), 
         !is.na(targetTaxonName)) %>%
  unique()

##### CropPoll #####

# Download from https://github.com/ibartomeus/OBservData/tree/master/Final_Data

#CropPoll data
crops = read.csv("Data_Raw/Interactions/croppol/CropPol_field_level_data.csv")
pollinators = read.csv("Data_Raw/Interactions/croppol/CropPol_sampling_data.csv")

#Cleaning CropPoll data
# select joining columns (`study_id` and `site_id`) and column of interest and
# get distinct rows
crops_slim <- crops %>%
  select(study_id, site_id, crop, latitude, longitude) %>%
  distinct()

pollinators_slim <- pollinators %>%
  select(study_id, site_id, pollinator) %>%
  distinct()

# join columns of interest 
interactions <- full_join(crops_slim, pollinators_slim, by = join_by(study_id, site_id))

# Filter for only complete rows in crop and pollinator
interactions[interactions == ""] <- NA

interactions <- interactions[complete.cases(interactions[, c(3, 6)]), ]

#take only needed columns and rename
crop_poll = interactions[,c(3,6,4,5)]

#add interaction type
crop_poll$interactionTypeName = "visitsFlowersOf"

#rename columns
crop_poll = crop_poll %>%
  rename(
    targetTaxonName = crop,
    sourceTaxonName = pollinator,
    decimalLatitude = latitude,
    decimalLongitude = longitude
  ) %>%
  distinct()

##### Crop Common Names #####

# Crop common names: partly from manual research and partly from
# https://www.fao.org/fileadmin/templates/ess/documents/world_census_of_agriculture/appendix4_r7.pdf

# Read in list of crop latin and common names
crop_names <- read.csv('Data_Raw/Species_Traits_Attributes/cropNames.csv')

# 1. Parse names
crop_names_parsed <- bdc_clean_names(sci_names = unique(crop_names$latin), save_outputs = FALSE)

# 2. Harmonize names
names_harm_crops <- harmonize_names(
  names = unique(crop_names_parsed$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_crops.tsv",
  names_harm_file = "Temp/names_harmonized_crops.tsv",
  gnv_path = gnverifier_path
)

# Merge with original crop names
crop_names <- crop_names %>%
  left_join(crop_names_parsed, by = join_by("latin" == "scientificName")) %>%
  left_join(names_harm_crops, by = c("names_clean" = "ScientificName")) %>%
  rename(
    original = latin,
    harmonized = names_harm
  ) %>%
  # Apply manual fixes where provided
  mutate(harmonized = ifelse(!is.na(latin_new) & latin_new != "", latin_new, harmonized)) %>%
  select(original, harmonized, common) %>%
  distinct()

#### Join Interactions ####

crop_poll$isCrop = NULL

#add source column to each
globi$source    <- "globi"
hosts$source    <- "hosts"
caldwell$source <- "caldwell"
simon$source    <- "simon"
dl$source       <- "dl"
bigbook$source  <- "bigbook"
fowler$source   <- "fowler"
crop_poll$source <- "cropPoll"
shropTall$source <- "ShropTall"

dfs <- list(globi, hosts, caldwell, simon, dl, bigbook, fowler, crop_poll, shropTall)

dfs_chr <- lapply(dfs, function(x) {
  x %>%
    as_tibble() %>%
    mutate(across(everything(), as.character))
})

interactions <- bind_rows(dfs_chr)

interactions = interactions %>%
  distinct()

# setwd("~/Desktop")
# 
# write.csv(interactions, "interactions_raw.csv", row.names = F)
# 
# write.csv(cropNames_for_harm, "crop_names.csv", row.names = F)

#### Harmonize Names ####

# Clean names with bdc
pollinator_names <- bdc_clean_names(
  sci_names = unique(interactions$sourceTaxonName),
  save_outputs = FALSE
)

plant_names <- interactions$targetTaxonName %>%
  as.character() %>%
  unique() %>%
  bdc_clean_names(save_outputs = FALSE)

# Join cleaned names back
interactions <- interactions %>%
  left_join(pollinator_names, by = join_by("sourceTaxonName" == "scientificName")) %>%
  dplyr::select(-c(.uncer_terms, .infraesp_names, quality)) %>%
  rename(sourceTaxonName_clean = names_clean) %>%
  left_join(plant_names, by = join_by("targetTaxonName" == "scientificName")) %>%
  dplyr::select(-c(.uncer_terms, .infraesp_names, quality)) %>%
  rename(targetTaxonName_clean = names_clean)

# Get unique names for harmonization
names_plants <- plant_names %>% pull(names_clean) %>% unique()
names_animals <- pollinator_names %>% pull(names_clean) %>% unique()

# Harmonize names
names_harm_plants <- harmonize_names(
  names = names_plants,
  higher_tax = "Plantae",
  names_file = "Temp/names_interactions_plants.tsv",
  names_harm_file = "Temp/names_harmonized_interactions_plants.tsv",
  gnv_path = gnverifier_path
)
# save all harmonized names
names_harm_plants_all <- names_harm_plants

names_harm_animals <- harmonize_names(
  names = names_animals,
  higher_tax = "Animalia",
  names_file = "Temp/names_interactions_animals.tsv",
  names_harm_file = "Temp/names_harmonized_interactions_animals.tsv",
  gnv_path = gnverifier_path
)

# Apply Jepson taxonomy to plants
names_harm_plants <- apply_jepson_taxonomy(names_harm_plants)

# Additional genus-level Jepson matching for interactions
# (plants may only be identified to genus level in interaction databases)
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

current_genera <- jepson_names %>%
  filter(status == "Native") %>%
  mutate(accepted_genus = word(names_harm_jepson, 1)) %>%
  pull(accepted_genus) %>%
  unique() %>%
  na.omit()

# Pass 5: genus-level match for plants not matched in apply_jepson_taxonomy()
names_harm_plants <- names_harm_plants %>%
  mutate(
    ScientificName_genus = word(ScientificName, 1),
    names_harm_jepson_genus = case_when(
      ScientificName_genus %in% current_genera ~ ScientificName_genus,
      TRUE ~ NA_character_
    ),
    names_harm = coalesce(names_harm, names_harm_jepson_genus)
  ) %>%
  dplyr::select(-c(ScientificName_genus, names_harm_jepson_genus)) %>%
  unique()

# Join harmonized names back with interactions
interactions_harm <- interactions %>%
  left_join(names_harm_animals,
            by = c("sourceTaxonName_clean" = "ScientificName")) %>%
  rename(sourceTaxonName_harm = names_harm) %>%
  dplyr::select(-c(Kind, SortScore, MatchType, EditDistance,
                   MatchedName, MatchedCanonical, TaxonId, CurrentName,
                   DataSourceId, DataSourceTitle, ClassificationPath, Error)) %>%
  left_join(names_harm_plants,
            by = c("targetTaxonName_clean" = "ScientificName")) %>%
  rename(targetTaxonName_harm = names_harm) %>%
  dplyr::select(-c(Kind, SortScore, MatchType, EditDistance,
                   MatchedName, MatchedCanonical, TaxonId, CurrentName,
                   DataSourceId, DataSourceTitle, ClassificationPath, Error))

# When plant is a crop, use crop harmonized name
interactions_harm <- interactions_harm %>%
  left_join(crop_names %>% select(original, harmonized), by = c("targetTaxonName_clean" = "original")) %>%
  mutate(targetTaxonName_harm = if_else(!is.na(harmonized), harmonized, targetTaxonName_harm)) %>%
  select(-harmonized)

# Remove rows where harmonization did not yield a match
interactions_harm <- interactions_harm %>%
  filter(!is.na(sourceTaxonName_harm),
         !is.na(targetTaxonName_harm)) %>%
  distinct()

#### ADD CROP CORRECTIONS HERE ####

## Remove rows where taxonomic harmonization did not yield match
# Mostly removes interactions with plants that we will not consider in our models (genus not in CA)
interactions_harm <- interactions_harm %>%
  filter(!is.na(sourceTaxonName_harm),
         !is.na(targetTaxonName_harm)) %>%
  distinct()

# Save
# write.csv(interactions_harm, "Data_Clean/Interactions/BigBook_ints.csv")

###### Identify Crops #######

crop_names = crop_names[,c(2,3)]
names(crop_names) = c("targetTaxonName_harm", "cropCommonName")
crop_names$isCrop = "yes"

interactions_harm = left_join(interactions_harm, crop_names)

#Set "yes" in 'isCrop' for rows where source is "cropPol"
interactions_harm <- interactions_harm %>%
  mutate(isCrop = ifelse(source == "cropPol", "yes", isCrop))

# Now for rows where plant matches any where isCrop == "yes", set 'isCrop' == "yes

# Get the unique 'lower' values for datasetName == "croppol"
crop_vals <- interactions_harm %>%
  filter(isCrop == "yes") %>%
  pull(targetTaxonName_harm) %>%
  unique()

# Update 'isCrop' for matching values
interactions_harm <- interactions_harm %>%
  mutate(isCrop = ifelse(targetTaxonName_harm %in% crop_vals, "yes", isCrop))

#take unique (assuming we don't want exact duplicates)
interactions_harm = interactions_harm %>%
  distinct()

###### De-duplicate ######

#Did not do this yet

#Probably only works for globi records: 

#Concatenate institution code, collection code, and catalog number to search for duplicates if it is museum records. Source citation will also work for most

#If observation record, use the source citation column

#### FOR OCCUPANCY MODEL, FIGURE OUT DUPLICATED INTERACTIONS ####

### Remove Gymnosperm flower visitation interactions
gymnosperm_orders <- c(
  "Ginkgoales", "Cycadales", "Araucariales",
  "Cupressales", "Pinales", "Ephedrales",
  "Welwitschiales", "Gnetales"
)

gymnosperm_genera <- names_harm_plants_all %>%
  filter(
    !is.na(names_harm),
    str_detect(ClassificationPath, paste(gymnosperm_orders, collapse = "|"))
  ) %>%
  mutate(genus = word(names_harm, 1)) %>%
  distinct(genus)

interactions_harm <- interactions_harm %>%
  filter(!(targetTaxonGenusName %in% gymnosperm_genera$genus & interactionTypeName == "visitsFlowersOf"),
         !(targetTaxonOrderName %in% gymnosperm_orders & interactionTypeName == "visitsFlowersOf")
  )

#Save outputs

write.csv(interactions_harm, "Temp/ints_all_clean.csv", row.names = F)

###### Final Interaction Data ######

# Clear workspace
rm(list = ls(all = TRUE))
# Load data

#List of all species that we have SDMs for
sdm_data = read.csv('Data_Raw/SDM_input/sdm_data.csv')

# All interactions
ints = read.csv("Temp/ints_all_clean.csv")

# Remove TaxonomicStatus.x and TaxonomicStatus.y
ints <- ints %>%
  select(-c(TaxonomicStatus.x, TaxonomicStatus.y))

# Clean decimalLatitude and decimalLongitude
bad_lat <- grepl("^[NSEWnsew]", ints$decimalLatitude)
bad_lon <- grepl("^[NSEWnsew]", ints$decimalLongitude)

ints$decimalLatitude[bad_lat] <- NA
ints$decimalLongitude[bad_lon] <- NA

ints$decimalLatitude <- as.numeric(ints$decimalLatitude)
ints$decimalLongitude <- as.numeric(ints$decimalLongitude)

# Gather lists of genera in the sdm pollinators and plants
sdm_lower = filter(sdm_data, taxon == "plants")
sdm_lower_genera  <- unique(tstrsplit(sdm_lower$species, " ", fixed = TRUE, keep = 1)[[1]])

sdm_higher = filter(sdm_data, taxon != "plants")
sdm_higher_genera <- unique(tstrsplit(sdm_higher$species, " ", fixed = TRUE, keep = 1)[[1]])

# Make genus and species columns in interaction dataset (don't use the genus columns that currently exist, don't have those for all records)
ints$lowerGenus  <- sub(" .*", "", ints$targetTaxonName_harm)
ints$higherGenus <- sub(" .*", "", ints$sourceTaxonName_harm)

# Filter interaction dataset only for the genera in the sdm data
#### Fix 6/30/2025: This is where the issue was of losing many of the crop interactions. Have to exclude crops from the lower filter #####
ints <- ints %>%
  filter(
    (isCrop == "yes") | 
      (is.na(isCrop) & lowerGenus %in% sdm_lower_genera)
  )

ints = filter(ints, higherGenus %in% sdm_higher_genera)

# Get rid of all undesirable interaction types (don't remove NAs for now)
unique(ints$interactionTypeName)
ints = filter(ints, !(interactionTypeName %in% c("parasiteOf", "hemiparasiteOf",  "pathogenOf", "livesInsideOf", "livesUnder")))

##### Separate Lep Host vs Nectar plant interactions #####

#identify Lepidoptera genera
leps = filter(sdm_data, taxon %in% c("butterflies", "moths"))
lep_genera <- unique(sub(" .*", "", leps$species))

## *** CHECK THIS STEP: may not be the best way to do it if there si overlap in genus names between taxa
# subset for non-Lep data
ints_non_lep = filter(ints, !(higherGenus %in% lep_genera))

#For now we are assuming that anything that is not a Lep that has made it this far into the cleaning is a flower-visitation interaction

#change everything that is not a lep to "visitsFlowersOf" (retain collectsPollenOf)

ints_non_lep$interactionTypeName <- as.character(ints_non_lep$interactionTypeName)

ints_non_lep$interactionTypeName <-
  ifelse(!is.na(ints_non_lep$interactionTypeName) &
           ints_non_lep$interactionTypeName == "collectsPollenOf",
         "collectsPollenOf",
         "visitsFlowersOf")

#filter for only lep interactions
ints_leps = filter(ints, higherGenus %in% lep_genera)

#for lep interactions, first filter out known host, flower visitation interactions
unique(ints_leps$interactionTypeName)
ints_leps_confirmed = filter(ints_leps, interactionTypeName %in% c("hasHost", "visitsFlowersOf", "pollinates", "mutualistOf", "laysEggsOn", "laysEggsIn"))

# Modify the interactionTypeName column
ints_leps_confirmed <- ints_leps_confirmed %>%
  mutate(interactionTypeName = case_when(
    interactionTypeName %in% c("pollinates", "mutualistOf") ~ "visitsFlowersOf",
    interactionTypeName %in% c("laysEggsOn", "laysEggsIn") ~ "hasHost",
    TRUE ~ interactionTypeName  # Keep other values unchanged
  ))

#eats, interactsWith, adjacentTo, visits, ecologicallyRelatedTo, coOccursWith, preysOn are all ambiguous, could be pollination or herbivory

#filter for non-confirmed interactions
ints_leps_non_confirmed = filter(ints_leps, !(interactionTypeName %in% c("hasHost", "visitsFlowersOf", "pollinates", "mutualistOf", "laysEggsOn", "laysEggsIn")))

# we only have interactions remaining from globi and discover life
unique(ints_leps_non_confirmed$source)
unique(ints_leps_non_confirmed$sourceLifeStageName)


# Harmonize the sourceLifeStageName column to the right categories
ints_leps_non_confirmed <- ints_leps_non_confirmed %>%
  mutate(sourceLifeStageName = case_when(
    # Convert "A" to "adult"
    grepl("^A$", sourceLifeStageName, ignore.case = TRUE) ~ "adult",
    
    # Convert "U" to NA
    grepl("^U$", sourceLifeStageName, ignore.case = TRUE) ~ NA_character_,
    
    # Convert "head capsules" to NA
    grepl("head capsules", sourceLifeStageName, ignore.case = TRUE) ~ NA_character_,
    
    # Convert "lifee" to NA
    grepl("lifee", sourceLifeStageName, ignore.case = TRUE) ~ NA_character_,
    
    # Convert "pupae" to "pupa"
    grepl("pupae", sourceLifeStageName, ignore.case = TRUE) ~ "pupa",
    
    # Convert "adullt" to "adult"
    grepl("adullt", sourceLifeStageName, ignore.case = TRUE) ~ "adult",
    
    # Replace "dult" with "adult"
    grepl("\\bdult\\b", sourceLifeStageName, ignore.case = TRUE) ~ "adult",
    
    # Replace "immature; juvenile" with "larva"
    grepl("immature; juvenile", sourceLifeStageName, ignore.case = TRUE) ~ "larva",
    
    # Match variations of "adult" (precedence over other stages if combined with "pupa" or "larva")
    grepl("\\badult\\b", sourceLifeStageName, ignore.case = TRUE) & !grepl("\\b(larva|larvae|caterpillar|larval|pupa|pupal)\\b", sourceLifeStageName, ignore.case = TRUE) ~ "adult",
    
    # Match variations of "egg"
    grepl("\\b(egg|ova|ovum)", sourceLifeStageName, ignore.case = TRUE) ~ "egg",
    
    # Match variations of "larva", "larvae", "caterpillar", "larval stage", "Larval", "immature
    grepl("\\b(larva|larvae|caterpillar|larval stage|larval|immature)\\b", sourceLifeStageName, ignore.case = TRUE) ~ "larva",
    
    # Match variations of "pupa", "pupal", or "cocoon" only if "adult" is not present
    grepl("\\b(pupa|pupal|cocoon)\\b", sourceLifeStageName, ignore.case = TRUE) &
      !grepl("\\badult\\b", sourceLifeStageName, ignore.case = TRUE) ~ "pupa",
    
    # If both "adult" and larval/pupal terms are present, classify accordingly
    grepl("\\badult\\b", sourceLifeStageName, ignore.case = TRUE) &
      grepl("\\b(larva|larvae|caterpillar|larval|immature|pupa|pupal|cocoon)\\b",
            sourceLifeStageName, ignore.case = TRUE) ~ case_when(
              
              # larva group
              grepl("\\b(larva|larvae|caterpillar|larval|immature)\\b",
                    sourceLifeStageName, ignore.case = TRUE) ~ "larva",
              
              # pupa group
              grepl("\\b(pupa|pupal|cocoon)\\b",
                    sourceLifeStageName, ignore.case = TRUE) ~ "pupa"
            ),
    
    # Assign NA or blank values to NA
    is.na(sourceLifeStageName) | sourceLifeStageName == "" ~ NA_character_,
    
    # Keep other values unchanged (optional, depends on your goals)
    TRUE ~ sourceLifeStageName
  )) %>%
  
  # Manual fixes after the case_when logic
  mutate(sourceLifeStageName = gsub("\\bdult\\b", "adult", sourceLifeStageName, ignore.case = TRUE)) %>%
  mutate(sourceLifeStageName = gsub("immature; juvenile", "larva", sourceLifeStageName, ignore.case = TRUE))

# Check that it worked
unique(ints_leps_non_confirmed$sourceLifeStageName)

#for globi, anything that lists larva, egg as the lifeStage is assumed to be hasHost, FIX 6/3/2025, we were previously forcing anything with "adult" lifestage to visitsFlowersOf, but not doing this anymore because these could equally likely be ovipositing interactions
ints_leps_non_confirmed <- ints_leps_non_confirmed %>%
  mutate(interactionTypeName = case_when(
    # If sourceLifeStageName is "larva" or "egg", set interactionTypeName to "hasHost" ("pupa" also an option, but not always safe to assume pupa on plant is the hostplant)
    sourceLifeStageName %in% c("larva", "egg") ~ "hasHost",
    
    # # If sourceLifeStageName is "adult", set interactionTypeName to "visitsFlowersOf"
    # sourceLifeStageName == "adult" ~ "visitsFlowersOf",
    
    # Otherwise, keep the current value of interactionTypeName
    TRUE ~ interactionTypeName
  ))

#separate out new confirmed, add to df, then add to rest of interactons
ints_leps_confirmed2 = filter(ints_leps_non_confirmed, interactionTypeName %in% c("hasHost", "visitsFlowersOf"))

ints_leps_confirmed = rbind(ints_leps_confirmed, ints_leps_confirmed2)

ints_clean = rbind(ints_non_lep, ints_leps_confirmed)

#take unique
ints_clean = ints_clean %>%
  distinct()

#### remove extra hummingbirds ####
hum_keep = c("Calypte anna",
             "Calypte costae",
             "Archilochus alexandri",
             "Selasphorus calliope",
             "Selasphorus Platycercus",
             "Selasphorus rufus",
             "Selasphorus sasin",
             "Selasphorus sasin sedentarius")

hum = sdm_data %>%
  filter(taxon == "hummingbirds") %>%
  select(species, species_full) %>%
  distinct()

# species not in hum_keep
sp1 <- setdiff(hum$species, hum_keep)

# species_full not in hum_keep
sp2 <- setdiff(hum$species_full, hum_keep)

# final unique list
hum_remove <- unique(c(sp1, sp2))

ints_clean = ints_clean %>%
  filter(!(sourceTaxonName_harm %in% hum_remove))

#save version with crops
ints_clean_withCrops = ints_clean

#separate out for only sdm plants
ints_clean = filter(ints_clean, lowerGenus %in% sdm_lower_genera)

#separate out still non-confirmed. Leaving out for now. 
ints_leps_non_confirmed = filter(ints_leps_non_confirmed, !(interactionTypeName %in% c("hasHost", "visitsFlowersOf")))

#take unique
ints_leps_non_confirmed = ints_leps_non_confirmed %>%
  distinct()


#Save full versions 
write.csv(ints_clean, "Data_Clean/Interactions/ints_final_clean_full.csv", row.names = F)
write.csv(ints_clean_withCrops, "Data_Clean/Interactions/ints_final_clean_withCrops_full.csv", row.names = F)
write.csv(ints_leps_non_confirmed, "Data_Clean/Interactions/int_leps_unconfirmed.csv", row.names = F)

#Create and save abbreviated versions
ints_clean_slim <- ints_clean %>%
  dplyr::select(
    higherGenus, sourceTaxonName_harm, interactionTypeName,
    lowerGenus, targetTaxonName_harm,
    decimalLatitude, decimalLongitude, localityName
  ) %>%
  rename(
    sourceTaxonGenusName = higherGenus,
    targetTaxonGenusName = lowerGenus
  ) %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%   # convert "" → NA
  distinct()

ints_clean_withCrops_slim <- ints_clean_withCrops %>%
  select(
    higherGenus, sourceTaxonName_harm, interactionTypeName,
    lowerGenus, targetTaxonName_harm,
    decimalLatitude, decimalLongitude, localityName, isCrop, cropCommonName
  ) %>%
  rename(
    sourceTaxonGenusName = higherGenus,
    targetTaxonGenusName = lowerGenus
  ) %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%   # convert "" → NA
  distinct()

# save slim versions
write.csv(ints_clean_slim, "Data_Clean/Interactions/ints_final_clean_slim.csv", row.names = F)
write.csv(ints_clean_withCrops_slim, "Data_Clean/Interactions/ints_final_clean_withCrops_slim.csv", row.names = F)


##### Check how many new interactions from other sources than globi #####
# int datasets filtered
globi_ints <- ints_clean_withCrops %>%
  filter(source == "globi") %>%
  select(higherGenus, sourceTaxonName_harm, interactionTypeName,
         lowerGenus, targetTaxonName_harm) %>%
  unique()

# species level
nonglobi_ints <- ints_clean_withCrops %>%
  filter(source != "globi") %>%
  select(higherGenus, sourceTaxonName_harm, interactionTypeName,
         lowerGenus, targetTaxonName_harm) %>%
  unique()

# compute shared/unique counts and percentages
compare_interactions <- function(df1, df2, cols, label = "") {
  
  set1 <- df1 %>% select(all_of(cols)) %>% distinct()
  set2 <- df2 %>% select(all_of(cols)) %>% distinct()
  
  shared     <- dplyr::intersect(set1, set2)
  unique_1   <- dplyr::setdiff(set1, set2)
  unique_2   <- dplyr::setdiff(set2, set1)
  total_union <- dplyr::union(set1, set2)
  
  n_shared   <- nrow(shared)
  n_unique_1 <- nrow(unique_1)
  n_unique_2 <- nrow(unique_2)
  n_total    <- nrow(total_union)
  
  tibble::tibble(
    level              = label,
    n_source1           = nrow(set1),
    n_source2           = nrow(set2),
    n_union            = n_total,
    n_shared           = n_shared,
    n_unique_source1   = n_unique_1,
    n_unique_source2   = n_unique_2,
    pct_shared_of_union      = n_shared / n_total * 100,
    pct_unique_source1_of_union = n_unique_1 / n_total * 100,
    pct_unique_source2_of_union = n_unique_2 / n_total * 100,
    pct_of_source1_shared = n_shared / nrow(set1) * 100,
    pct_of_source2_shared = n_shared / nrow(set2) * 100
  )
}

# species level (source-target-interactionType combo)
species_cols <- c("sourceTaxonName_harm", "interactionTypeName", "targetTaxonName_harm")

species_comparison <- compare_interactions(
  globi_ints, nonglobi_ints,
  cols = species_cols,
  label = "species"
)

# genus level
genus_globi_ints <- ints_clean_withCrops %>%
  filter(source == "globi") %>%
  select(higherGenus, interactionTypeName, lowerGenus) %>%
  unique()

genus_nonglobi_ints <- ints_clean_withCrops %>%
  filter(source != "globi") %>%
  select(higherGenus, interactionTypeName, lowerGenus) %>%
  unique()

genus_cols <- c("higherGenus", "interactionTypeName", "lowerGenus")

genus_comparison <- compare_interactions(
  genus_globi_ints, genus_nonglobi_ints,
  cols = genus_cols,
  label = "genus"
)

# combine into one summary table
overlap_summary <- bind_rows(species_comparison, genus_comparison)
overlap_summary


# genus only level