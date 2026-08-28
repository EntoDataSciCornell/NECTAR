#### Load packages ####
library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggalluvial)
library(tidyterra)
library(bipartite)
library(nngeo)

#### Supporting scripts ####
source("Code/12_Analyses/12_00_LoadData.R")

#### General summary information about datasets ####
## How many species did we aggregate data for?
plants <- checklist_cleaned %>%
  filter(taxon == "plants") %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(plants), "native plant species"))

pols <- checklist_cleaned %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies")) %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(pols), "native pollinator species"))

bees <- checklist_cleaned %>%
  filter(taxon == "bees") %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(bees), "native bee species"))

hoverflies <- checklist_cleaned %>%
  filter(taxon == "hoverflies") %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(hoverflies), "native hoverfly species"))

bfly <- checklist_cleaned %>%
  filter(taxon == "butterflies") %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(bfly), "native butterfly species"))

moths <- checklist_cleaned %>%
  filter(taxon == "moths") %>%
  pull(genus_species) %>%
  unique()

message(paste("Aggregated data for ", length(moths), "native butterfly species"))


#### Summarise Raw interaction datasets ####
interactions_full <- read.csv("Data_Clean/Interactions/ints_final_clean_full.csv") %>%
  mutate(sourceTaxonGenusName = higherGenus,
         targetTaxonGenusName = lowerGenus)

# which pollinator taxa does each genus belong to?
pol_genus_tax <- checklist_cleaned %>%
  dplyr::select(genus, taxon) %>%
  unique()

interactions_full_tax <- interactions_full %>% 
  left_join(pol_genus_tax, by = join_by("sourceTaxonGenusName" == "genus"))

# Remove hummingbirds (not included in our analysis)
interactions_full_tax <- interactions_full_tax %>%
  filter(taxon != "hummingbirds")

# unique sp-sp level ints
interactions_full_tax %>% 
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, taxon) %>% 
  filter(sourceTaxonName_harm %in% checklist_cleaned$genus_species & 
           targetTaxonName_harm %in% checklist_cleaned$genus_species) %>% 
  unique() %>%
  group_by(taxon) %>% 
  summarise(n=n())


interactions_full_tax_nectar <- interactions_full_tax %>%
  filter(interactionTypeName == "visitsFlowersOf")
interactions_full_tax_host <- interactions_full_tax %>%
  filter(interactionTypeName == "hasHost")


alluvial_plot <- function(dataset, genus_threshold = 10, panel_title = "", source_levels = NULL) {
  # Calculate total interactions (using all rows, no collapsing)
  total_interactions <- nrow(dataset)
  
  # Find genera that represent >= threshold% of data
  genus_freq <- dataset %>%
    group_by(sourceTaxonGenusName) %>%
    summarise(total_freq = n(), .groups = 'drop') %>%
    mutate(percentage = total_freq / total_interactions * 100)
  
  # Get genera that meet the threshold
  major_genera <- genus_freq$sourceTaxonGenusName[genus_freq$percentage >= genus_threshold]
  
  # Create display labels and count combinations
  plot_data <- dataset %>%
    mutate(
      sourceTaxonGenusName_display = ifelse(sourceTaxonGenusName %in% major_genera, 
                                            sourceTaxonGenusName, 
                                            "Other"),
      # Set factor levels if provided
      source = if(!is.null(source_levels)) factor(source, levels = source_levels) else source
    ) %>%
    # Now count the combinations for plotting efficiency
    count(source, taxon, sourceTaxonGenusName_display, name = "freq")
  
  # Create the alluvial diagram
  p <- ggplot(plot_data,
              aes(y = freq, axis1 = source, axis2 = taxon, axis3 = sourceTaxonGenusName_display)) +
    geom_alluvium(aes(fill = source), width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/12, fill = "lightgray", color = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Taxon", "Genus"), 
                     expand = c(.05, .05)) +
    scale_fill_viridis_d(name = "Source", drop = FALSE) +
    labs(title = panel_title,
         y = "Number of Interactions") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom")
  
  return(p)
}

# Get all unique sources from both datasets
all_sources <- unique(c(interactions_full_tax_nectar$source, 
                        interactions_full_tax_host$source))

# Create both plots with consistent color mapping
p1 <- alluvial_plot(interactions_full_tax_nectar, 
                    genus_threshold = 5, 
                    panel_title = "Flower Visitation",
                    source_levels = all_sources)

p2 <- alluvial_plot(interactions_full_tax_host, 
                    genus_threshold = 5, 
                    panel_title = "Host Interactions",
                    source_levels = all_sources)

p1 <- p1 + theme(legend.position = "bottom")
p2 <- p2 + theme(legend.position = "bottom")
combined_plot <- (p1 + theme(legend.position = "none")) + p2 +
  plot_layout(guides = "collect")

# Save the combined plot
ggsave("Figures/Supplementary/interactionNetworkRaw_combined_alluvial.png", 
       combined_plot, 
       width = 20, 
       height = 10, 
       dpi = 300)

#### # species from SDMs and phenometrics ####
## Plant sdms
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

# Filter pollinator SDMs to checklist species
pollinator_species_in_checklist <- checklist_cleaned %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies")) %>%
  pull(genus_species)

# Subset the raster stack by name
pollinator_sdms_filtered <- pollinator_sdms[[pollinator_names %in% pollinator_species_in_checklist]]

# Filter plant SDMs to checklist species
plant_species_in_checklist <- checklist_cleaned %>%
  filter(taxon == "plants") %>%
  pull(genus_species)

plant_sdms_filtered <- plant_sdms[[plant_names %in% plant_species_in_checklist]]

### how many species did we run SDMs for?
message(paste(dim(plant_sdms_filtered)[3], "plant species ran SDMs"))
message(paste(dim(pollinator_sdms_filtered)[3], "pollinator species ran SDMs"))

### How about phenometrics?
phenology_filtered <- phenology %>%
  filter(!is.na(phen_start) & 
           !is.na(phen_end))

message(paste(length(names(plant_sdms_filtered)[names(plant_sdms_filtered) %in% phenology_filtered$genus_species]), "plant species ran SDMs+Pheno"))
message(paste(length(names(pollinator_sdms_filtered)[names(pollinator_sdms_filtered) %in% phenology_filtered$genus_species]), "pollinator species ran SDMs+Pheno"))

checklist_cleaned %>% filter(genus_species %in% names(pollinator_sdms_filtered)[names(pollinator_sdms_filtered) %in% phenology_filtered$genus_species]) %>% group_by(taxon) %>% summarise(n=n())

### Richness
## NOTE: Maybe update to only show richness for species with both sdms+pheno? Or maybe even only richness of sp. in our predicted interaction network?? Or both
## sum continuous probabilities to get estimated richness
plant_richness <- sum(plant_sdms_filtered[[names(plant_sdms_filtered) %in% phenology_filtered$genus_species]], na.rm = TRUE)
pollinator_richness <- sum(pollinator_sdms_filtered[[names(pollinator_sdms_filtered) %in% phenology_filtered$genus_species]], na.rm = TRUE)

ecoregions_boundary <- ecoregions %>%
  st_union() %>%
  nngeo::st_remove_holes()  # Removes all internal holes

plant_richness_plot <- ggplot() +
  geom_spatraster(data = plant_richness) +
  geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_gradient2(
    low = "white",
    mid = "#BAFF76",
    high = "#3B7004",
    name = "Expected\nPlant\nRichness",
    midpoint = 250,
    na.value = "white"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL) +
  coord_sf()

pollinator_richness_plot <- ggplot() +
  geom_spatraster(data = pollinator_richness) +
  geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_gradient2(
    low = "white",
    mid = "#F9B753",
    high = "#331E0F",
    name = "Expected\nPollinator\nRichness",
    midpoint = 250,
    na.value = "white"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL) +
  coord_sf()

ggsave("Figures/richness/plant_richness.pdf", 
       plot = plant_richness_plot,
       width = 200, 
       height = 100, 
       units = "mm",
       dpi = 2000)
ggsave("Figures/richness/plant_richness.svg", 
       plot = plant_richness_plot,
       width = 200, 
       height = 100, 
       units = "mm",
       dpi = 2000)


# Define taxa groups and their midpoints/colors for the gradient
pollinator_taxa <- c("bees", "hoverflies", "moths", "butterflies")

taxa_palette <- list(
  bees        = list(mid = "#FFCC00", high = "#89700C", midpoint = 200),
  hoverflies  = list(mid = "#F76700", high = "#974A13", midpoint = 50),
  moths       = list(mid = "#BD9ECD", high = "#4A235A", midpoint = 600),
  butterflies = list(mid = "#E69190", high = "#9B4543", midpoint = 70)
)

# Compute per-taxa richness rasters
taxa_richness_rasters <- map(pollinator_taxa, function(tx) {
  # Get species names for this taxon that are in phenology
  tx_species <- checklist_cleaned %>%
    filter(taxon == tx) %>%
    pull(genus_species)
  
  tx_species_filtered <- tx_species[tx_species %in% phenology_filtered$genus_species]
  
  # Sum SDM layers for this taxon
  sum(pollinator_sdms_filtered[[names(pollinator_sdms_filtered) %in% tx_species_filtered]], 
      na.rm = TRUE)
})
names(taxa_richness_rasters) <- pollinator_taxa

# Build sub-panel maps
taxa_richness_plots <- imap(taxa_richness_rasters, function(rast, tx) {
  pal <- taxa_palette[[tx]]
  
  ggplot() +
    geom_spatraster(data = rast) +
    geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 0.5) +
    scale_fill_gradient2(
      low      = "white",
      mid      = pal$mid,
      high     = pal$high,
      midpoint = pal$midpoint,
      name     = paste0(str_to_title(tx), "\nRichness"),
      na.value = "white"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 9)
    ) +
    labs(title = str_to_title(tx), x = NULL, y = NULL) +
    coord_sf()
})

# Combine with patchwork
bottom_row <- wrap_plots(taxa_richness_plots, nrow = 1)

pollinator_richness_figure <- pollinator_richness_plot / bottom_row +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

ggsave("Figures/richness/pollinator_richness_by_taxon.pdf",
       plot   = pollinator_richness_figure,
       width  = 200,
       height = 160,
       units  = "mm",
       dpi    = 2000)

ggsave("Figures/richness/pollinator_richness_by_taxon.svg",
       plot   = pollinator_richness_figure,
       width  = 200,
       height = 160,
       units  = "mm",
       dpi    = 2000)

#### % of checklist with each type of data ####
data_completeness_by_taxon <- checklist_cleaned %>%
  group_by(taxon) %>% 
  dplyr::select(taxon, genus_species) %>%
  summarise(n=n())

occurrences_completeness_by_taxon <- occurrences %>%
  filter(genus_species %in% checklist_cleaned$genus_species) %>%
  group_by(taxon) %>%
  dplyr::select(taxon, genus_species) %>%
  summarise(n_sp_occurrence = length(unique(genus_species)))

sdms_completeness_by_taxon <- checklist_cleaned %>%
  filter(genus_species %in% c(names(plant_sdms_filtered), names(pollinator_sdms_filtered))) %>%
  group_by(taxon) %>%
  dplyr::select(taxon, genus_species) %>%
  summarise(n_sp_sdms = length(unique(genus_species)))

pheno_completeness_by_taxon <- checklist_cleaned %>%
  filter(genus_species %in% phenology_filtered$genus_species) %>%
  group_by(taxon) %>%
  dplyr::select(taxon, genus_species) %>%
  summarise(n_sp_pheno = length(unique(genus_species)))

sdmandpheno_completeness_by_taxon <- checklist_cleaned %>%
  filter(genus_species %in% phenology_filtered$genus_species &
           genus_species %in% c(names(plant_sdms_filtered), names(pollinator_sdms_filtered))) %>%
  group_by(taxon) %>%
  dplyr::select(taxon, genus_species) %>%
  summarise(n_sp_both = length(unique(genus_species)))

sdmandpheno_completeness_species <- checklist_cleaned %>%
  filter(genus_species %in% phenology_filtered$genus_species &
           genus_species %in% c(names(plant_sdms_filtered), names(pollinator_sdms_filtered))) %>%
  group_by(taxon) %>%
  dplyr::select(genus_species, taxon) %>%
  distinct()

write.csv(sdmandpheno_completeness_species, "Data_Clean/Analyses/checklist_SDMandPheno.csv")



data_completeness <- data_completeness_by_taxon %>%
  left_join(occurrences_completeness_by_taxon) %>%
  left_join(sdms_completeness_by_taxon) %>%
  left_join(pheno_completeness_by_taxon) %>%
  left_join(sdmandpheno_completeness_by_taxon) %>%
  filter(taxon != "hummingbirds")

data_completeness_long <- data_completeness %>%
  pivot_longer(cols = -taxon, 
               names_to = "category", 
               values_to = "count")

category_labels <- c(
  "n" = "# Species in Checklist",
  "n_sp_occurrence" = "# Species with Occurrence Records",
  "n_sp_sdms" = "# Species with SDMs Run",
  "n_sp_pheno" = "# Species with Phenometrics",
  "n_sp_both" = "# Species with Both SDMs and Phenometrics"
)

# grouped bar chart
p_data_completeness <- ggplot(data_completeness_long, aes(x = taxon, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_discrete(labels = category_labels) +
  scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
  theme_cowplot() +
  labs(x = "Taxon", y = "Number of Species", fill = "Category")

ggsave("Figures/Supplementary/data_completeness_by_taxon.png", p_data_completeness, width = 10, height = 6, dpi = 300)

#### Summarize occurrence cleaning filtering ####
occurrences_flagged <- read.csv("Data_Clean/Occurrences/occurrences_flagged.csv")

# Count totals
n_flagged <- nrow(occurrences_flagged)
n_clean <- nrow(occurrences)

# Panel A: Records kept vs removed
records_summary <- data.frame(
  status = c("Kept", "Removed"),
  count = c(n_clean, n_flagged)
)

panel_a <- ggplot(records_summary, aes(x = status, y = count)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = "", y = "Number of Occurrence Records")

# Panel B: Reasons for removal (pie chart)
flag_cols <- names(occurrences_flagged)[grepl("^\\.", names(occurrences_flagged))]

# Create flag labels dictionary
flag_labels <- c(
  "overlaps_California" = "Outside focal\nregion",
  "coordinates_Precise" = "Imprecise\ncoordinates",
  "valid_date" = "Invalid date",
  "correct_taxon" = "Invalid taxon",
  "name_missing" = "Missing Name",
  "inst" = "Other",
  "dupExact" = "Other",
  "scientificName_empty" = "Other",
  "basisOfRecords_notStandard" = "Other",
  "uncer_terms" = "Other",
  "presence_record" = "Other")

flag_summary <- occurrences_flagged %>%
  select(all_of(flag_cols)) %>%
  summarise(across(everything(), ~sum(.x == FALSE, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "flag", values_to = "n_records") %>%
  filter(n_records > 0) %>%
  arrange(desc(n_records)) %>%
  mutate(
    flag_clean = gsub("^\\.", "", flag),
    flag_label = ifelse(flag_clean %in% names(flag_labels), 
                        flag_labels[flag_clean], 
                        gsub("_", " ", flag_clean)),
    pct = round(n_records / n_flagged * 100, 1)
  )

panel_b <- ggplot(flag_summary, aes(x = "", y = n_records, fill = flag_label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  theme(legend.position = "top") +
  labs(fill = "Reason for\nRemoving")

# Combine panels
combined_plot_occurrences <- plot_grid(panel_a, panel_b, labels = c("A", "B"), ncol = 2)

ggsave("Figures/Supplementary/occurrence_data_cleaning.png", combined_plot_occurrences, width = 10, height = 6, dpi = 300)


#### Distribution of occurrence records throughout CA ####
occurrences_CA <- occurrences_spatial %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(st_crs(ecoregions)) %>%
  filter(lengths(st_intersects(., ecoregions)) > 0) %>% 
  filter(taxon != "hummingbirds")

# Overall map with legend
map_all <- ggplot() +
  geom_sf(data = ecoregions, fill = "grey95", color = "black", linewidth = 0.3) +
  stat_bin_hex(data = st_coordinates(occurrences_CA) %>% as.data.frame() %>% 
                 setNames(c("X", "Y")),
               aes(x = X, y = Y), bins = 30, alpha = 0.7) +
  scale_fill_viridis_c(name = "# Occurrence Records", trans = "log10", 
                       limits = c(1, 300000),
                       breaks = c(1, 10, 100, 1000, 10000, 100000),
                       labels = c("1", "10", "100", "1,000", "10,000", "100,000")) +
  ggtitle("All Taxa") +
  theme_void() +
  theme(legend.position = "right")

# Maps by taxon - same scale, no legend
taxa <- unique(occurrences_CA$taxon)
taxon_maps <- lapply(taxa, function(tx) {
  occ_subset <- occurrences_CA[occurrences_CA$taxon == tx, ]
  
  ggplot() +
    geom_sf(data = ecoregions, fill = "grey95", color = "black", linewidth = 0.2) +
    stat_bin_hex(data = st_coordinates(occ_subset) %>% as.data.frame() %>% 
                   setNames(c("X", "Y")),
                 aes(x = X, y = Y), bins = 20, alpha = 0.7) +
    scale_fill_viridis_c(name = "Count", trans = "log10",
                         limits = c(1, 300000),
                         breaks = c(1, 10, 100, 1000, 10000, 100000),
                         labels = c("1", "10", "100", "1,000", "10,000", "100,000")) +
    ggtitle(tools::toTitleCase(tx)) +
    theme_void() +
    theme(legend.position = "none")
})

# Combine
top_row <- plot_grid(map_all, ncol = 1)
bottom_rows <- plot_grid(plotlist = taxon_maps, ncol = 5)
combined_map <- plot_grid(top_row, bottom_rows, ncol = 1, rel_heights = c(1, 0.5))

ggsave("Figures/Supplementary/occurrence_density_CA.png", combined_map, 
       width = 14, height = 16, dpi = 300, bg = "white")

#### Breakdown of SDM method used ####
#Visualize the number/proportion of species in each SDM type per taxa
#Pre-processing 

#list all CSV files recursively
csv_files <- list.files(
  path = "Data_Clean/SDMs/sdm_by_taxon",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

sdm_manifests <- list.files(
  path = "Data_Clean/SDMs/sdm_by_taxon",
  pattern = "\\continuous_manifest.csv$",
  full.names = TRUE,
  recursive = TRUE
)

sdm_types <- do.call(rbind, lapply(sdm_manifests, read.csv))
setDT(sdm_types)

# # keep only the 2nd occurrence per species (if duplicates exist)
# sdm_types <- sdm_types[
#   , .SD[.N >= 2][2],
#   by = species
# ]

# save for later in script if needed
sdm_types_allCols <- sdm_types

#take needed columns
sdm_types = sdm_types[,c(2,3)]

#write
fwrite(
  sdm_types,
  file = "Temp/sdm_types.csv"
)

sdm_types <- read.csv("Temp/sdm_types.csv")

# clean + recode
sdm_types <- sdm_types %>%
  filter(taxon != "hummingbirds") %>%
  mutate(
    model_type = recode(
      model_type,
      "SDM" = "Standard Ensemble Model",
      "ESM" = "Ensemble of Small Models",
      "ultraRare_alphahull" = "Simple Buffer"
    ),
    model_type = factor(
      model_type,
      levels = c(
        "Simple Buffer",
        "Ensemble of Small Models",
        "Standard Ensemble Model"
      )
    )
  )

# calculate raw counts
count_df <- sdm_types %>%
  count(taxon, model_type)

# plot raw stacked counts
p <- ggplot(count_df, aes(x = taxon, y = n, fill = model_type)) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(
    x = "Taxon",
    y = "Number of Species",
    fill = "Model Type"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

# calculate proportions within each taxon
prop_df <- sdm_types %>%
  count(taxon, model_type) %>%
  group_by(taxon) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# plot stacked proportions 
p2 <- ggplot(prop_df, aes(x = taxon, y = prop, fill = model_type)) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Taxon",
    y = "Proportion of Species",
    fill = "Model Type"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2

#Save 
ggsave("Figures/Supplementary/sdm_model_type_counts_by_taxon.pdf", p, width = 8, height = 6)
ggsave("Figures/Supplementary/sdm_model_type_proportion_by_taxon.pdf", p2, width = 8, height = 6)

#### Breakdown of Phenometric method used ####
phenology_with_taxon <- phenology %>%
  left_join(unique(checklist_cleaned[, c("genus_species", "taxon")])) %>%
  filter(taxon != "hummingbirds")

# Panel A: Count unique species by estimation type
phenology_summary <- phenology_with_taxon %>%
  filter(!is.na(phen_scale),
         !is.na(phen_start),
         !is.na(phen_end)) %>%
  mutate(
    estimation_type = case_when(
      phen_scale %in% c("PROV", "REGION", "JEPCODE", "CALIFORNIA") ~ "Type_1",
      phen_scale == "TYPE2" ~ "Type_2",
      phen_scale == "TYPE3" ~ "Type_3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(estimation_type)) %>%
  group_by(taxon, estimation_type) %>%
  summarise(n_species = n_distinct(genus_species), .groups = "drop")

panel_a <- ggplot(phenology_summary, aes(x = taxon, y = n_species, fill = estimation_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
  theme_cowplot() +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  labs(x = "Taxon", y = "# Species", fill = "Estimation Type")

# Panel B: Count species-regions by scale (Type_1 only)
scale_summary <- phenology_with_taxon %>%
  filter(phen_scale %in% c("PROV", "REGION", "JEPCODE", "CALIFORNIA")) %>%
  mutate(phen_scale = factor(phen_scale, levels = c("CALIFORNIA", "PROV", "REGION", "JEPCODE"))) %>%
  group_by(taxon, phen_scale) %>%
  summarise(n_sp_regions = n(), .groups = "drop")

panel_b <- ggplot(scale_summary, aes(x = taxon, y = n_sp_regions, fill = phen_scale)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
  theme_cowplot() +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  labs(x = "Taxon", y = "# Species-Regions", fill = "Phenology Scale")

# Combine
combined_plot <- plot_grid(panel_a, panel_b, labels = c("A", "B"), nrow = 2)

ggsave("Figures/Supplementary/phenology_summary.png", combined_plot, 
       width = 10, height = 10, dpi = 300)

#### Summarise # species by method used ####
supplementary_checklist <- read.csv("Tables/Supplementary_Checklist.csv")

# Total
supplementary_checklist %>% 
  filter(has_SDM == T & has_phenometrics == T) %>%
  unique() %>%
  group_by(taxon) %>%
  summarise(n=n())


# How many removed if we did not use lower tier pheno methods and/or hull approach for sdm
sp_topMethods <- supplementary_checklist %>% 
  filter(has_SDM == TRUE & has_phenometrics == T) %>%
  group_by(taxon) %>%
  filter(phenometrics_tier == "Tier1") %>%
  filter(genus_species %in% sdm_types_allCols[sdm_types_allCols$model_type %in% c("ESM","SDM"),]$species) 

sp_topMethods %>% 
  group_by(taxon) %>%
  summarise(n=n())

# num ints 
predicted_interactions %>%
  filter(method_type != "Level3_GEOOSIS",
         sourceTaxonType != "hummingbirds") %>%
  select(targetTaxonName_harm, sourceTaxonName_harm, interactionTypeName) %>%
  group_by(interactionTypeName) %>%
  summarise(n= n())

# num ints filtered for top methods
predicted_interactions %>%
  filter(method_type != "Level3_GEOOSIS",
         sourceTaxonType != "hummingbirds",
         targetTaxonName_harm %in% sp_topMethods$genus_species,
         sourceTaxonName_harm %in% sp_topMethods$genus_species) %>%
  select(targetTaxonName_harm, sourceTaxonName_harm, interactionTypeName) %>%
  group_by(interactionTypeName) %>%
  summarise(n= n())

# how many pairs?
predicted_interactions %>%
  filter(method_type != "Level3_GEOOSIS",
         sourceTaxonType != "hummingbirds") %>%
  select(targetTaxonName_harm, sourceTaxonName_harm, interactionTypeName) %>%
  unique() %>%
  group_by(interactionTypeName) %>%
  summarise(n= n())

# how many pairs, top method only?
predicted_interactions %>%
  filter(method_type != "Level3_GEOOSIS",
         sourceTaxonType != "hummingbirds",
         targetTaxonName_harm %in% sp_topMethods$genus_species,
         sourceTaxonName_harm %in% sp_topMethods$genus_species) %>%
  select(targetTaxonName_harm, sourceTaxonName_harm, interactionTypeName) %>%
  unique() %>%
  group_by(interactionTypeName) %>%
  summarise(n= n())

#### Size of interaction data ####
ints_pols_CAsp <- interactions %>%
  filter(interactionTypeName == "visitsFlowersOf") %>%
  mutate(sourceTaxonName_harm_gsp = sapply(strsplit(sourceTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         targetTaxonName_harm_gsp = sapply(strsplit(targetTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  filter(sourceTaxonName_harm_gsp %in% checklist_cleaned$genus_species)

message(paste(length(unique(ints_pols_CAsp$sourceTaxonName_harm_gsp)), "pollinator sp in interaction data"))

ints_plants_CAsp <- interactions %>%
  filter(interactionTypeName == "visitsFlowersOf") %>%
  mutate(sourceTaxonName_harm_gsp = sapply(strsplit(sourceTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         targetTaxonName_harm_gsp = sapply(strsplit(targetTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  filter(targetTaxonName_harm_gsp %in% checklist_cleaned$genus_species)

message(paste(length(unique(ints_plants_CAsp$targetTaxonName_harm_gsp)), "plant sp in interaction data"))
message(paste(length(unique(interactions$targetTaxonGenusName[interactions$targetTaxonGenusName %in% checklist_cleaned$genus])), "plant genera in interaction data"))

ints_bflys_CAsp_host <- interactions %>%
  filter(interactionTypeName == "hasHost") %>%
  mutate(sourceTaxonName_harm_gsp = sapply(strsplit(sourceTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         targetTaxonName_harm_gsp = sapply(strsplit(targetTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  filter(sourceTaxonName_harm_gsp %in% checklist_cleaned$genus_species[checklist_cleaned$taxon == "butterflies"])
ints_moths_CAsp_host <- interactions %>%
  filter(interactionTypeName == "hasHost") %>%
  mutate(sourceTaxonName_harm_gsp = sapply(strsplit(sourceTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         targetTaxonName_harm_gsp = sapply(strsplit(targetTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  filter(sourceTaxonName_harm_gsp %in% checklist_cleaned$genus_species[checklist_cleaned$taxon == "moths"])

message(paste(length(unique(ints_bflys_CAsp_host$sourceTaxonName_harm_gsp)), "butterflies with host in interaction data"))
message(paste(length(unique(ints_moths_CAsp_host$sourceTaxonName_harm_gsp)), "Lepidoptera with host in interaction data"))

ints_plants_CAsp_host <- interactions %>%
  filter(interactionTypeName == "hasHost") %>%
  mutate(sourceTaxonName_harm_gsp = sapply(strsplit(sourceTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         targetTaxonName_harm_gsp = sapply(strsplit(targetTaxonName_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
  filter(targetTaxonName_harm_gsp %in% checklist_cleaned$genus_species)

message(paste(length(unique(ints_plants_CAsp_host$targetTaxonName_harm_gsp)), "Lepidoptera host plant sp in interaction data"))
message(paste(length(unique(ints_plants_CAsp_host$targetTaxonGenusName[ints_plants_CAsp_host$targetTaxonGenusName %in% checklist_cleaned$genus])), "plant genera used as hosts in interaction data"))


# Total interaction records?
ints_full <- read.csv("Data_Clean/Interactions/ints_final_clean_full.csv")

message(paste("Compiled", nrow(ints_full), "interaction records"))

length(unique(c(ints_full$sourceTaxonName_harm, ints_full$targetTaxonName_harm)))

# How many georeferenced?
ints_with_geolocation <- ints_full %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

message(paste(nrow(ints_with_geolocation), "interaction records georeferenced"))

ints_CAsp_with_geolocation_inCA <- interactions_spatial %>%
  st_drop_geometry() %>%
  filter(!is.na(JEPCODE)) %>%
  filter(sourceTaxonName_harm %in% pollinator_checklists$genus_species,
         targetTaxonName_harm %in% plant_checklist$genus_species)

message(paste(nrow(ints_CAsp_with_geolocation_inCA), "interaction records georeferenced"))



#### Summarize predicted interaction network
## How many sp in the flower visitation and host data for Level1+Level2?
pred_L12_regions <- predicted_interactions %>%
  filter(method_type %in% c("Level1_Simple", "Level2_GEOOS"))

message(paste(length(unique(pred_L12_regions$sourceTaxonName_harm)), 
              "pollinators;",
              length(unique(pred_L12_regions$targetTaxonName_harm)),
              "plants"))

pred_simple_flowervisitation_regions <- predicted_interactions %>%
  filter(interactionTypeName =="visitsFlowersOf",
         method_type == "Level1_Simple") %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1)) %>%
  group_by(sourceTaxonName_harm, targetTaxonName_harm, targetTaxonGenusName) %>%
  summarise(n_regions = n())

message(paste("Predicted", 
              nrow(pred_simple_flowervisitation_regions), 
              "unique plant-pol flower vistation interactions between",
              length(unique(pred_simple_flowervisitation_regions$sourceTaxonName_harm)),
              "pollinators and",
              length(unique(pred_simple_flowervisitation_regions$targetTaxonName_harm)),
              "plants (",
              length(unique(pred_simple_flowervisitation_regions$targetTaxonGenusName)),
              "genera), occurring on average in", 
              mean(pred_simple_flowervisitation_regions$n_regions), 
              "regions"))

pred_L2_flowervisitation_regions <- predicted_interactions %>%
  filter(interactionTypeName =="visitsFlowersOf",
         method_type == "Level2_GEOOS") %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1)) %>%
  group_by(sourceTaxonName_harm, targetTaxonName_harm, targetTaxonGenusName) %>%
  summarise(n_regions = n())

message(paste("Predicted", 
              nrow(pred_L2_flowervisitation_regions), 
              "unique plant-pol flower vistation interactions between",
              length(unique(pred_L2_flowervisitation_regions$sourceTaxonName_harm)),
              "pollinators and",
              length(unique(pred_L2_flowervisitation_regions$targetTaxonName_harm)),
              "plants (",
              length(unique(pred_L2_flowervisitation_regions$targetTaxonGenusName)),
              "genera)"))

pred_L3_flowervisitation_regions <- predicted_interactions %>%
  filter(interactionTypeName =="visitsFlowersOf",
         method_type == "Level3_GEOOSIS") %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1)) %>%
  group_by(sourceTaxonName_harm, targetTaxonName_harm, targetTaxonGenusName) %>%
  summarise(n_regions = n())

message(paste("Predicted", 
              nrow(pred_L3_flowervisitation_regions), 
              "unique plant-pol flower vistation interactions between",
              length(unique(pred_L3_flowervisitation_regions$sourceTaxonName_harm)),
              "pollinators and",
              length(unique(pred_L3_flowervisitation_regions$targetTaxonName_harm)),
              "plants (",
              length(unique(pred_L3_flowervisitation_regions$targetTaxonGenusName)),
              "genera)"))

#### % of species and genera in checklist represented in interaction data ####
sp_with_interactions <- checklist_cleaned %>%
  dplyr::select(taxon, genus, genus_species) %>%
  unique() %>% 
  filter(genus_species %in% 
           c(word(interactions$sourceTaxonName_harm, 1, 2), 
             word(interactions$targetTaxonName_harm, 1, 2))) %>% 
  group_by(taxon) %>%
  summarise(n_sp_interactions = length(unique(genus_species)))

genera_with_interactions <- checklist_cleaned %>%
  dplyr::select(taxon, genus, genus_species) %>%
  unique() %>% 
  filter(genus %in% 
           c(interactions$sourceTaxonGenusName, 
             interactions$targetTaxonGenusName)) %>% 
  group_by(taxon) %>%
  summarise(n_sp_interactions = length(unique(genus)))

## compared to total # species
sp_by_taxon <- checklist_cleaned %>%
  group_by(taxon) %>% 
  dplyr::select(taxon, genus_species) %>%
  unique() %>%
  filter(taxon != "hummingbirds") %>%
  summarise(n_sp=n())

genera_by_taxon <- checklist_cleaned %>%
  group_by(taxon) %>% 
  dplyr::select(taxon, genus) %>%
  unique() %>%
  filter(taxon != "hummingbirds") %>%
  summarise(n_sp=n())

# Panel A: Species data
species_data <- sp_by_taxon %>%
  left_join(sp_with_interactions, by = "taxon") %>%
  rename(total = n_sp, with_interactions = n_sp_interactions) %>%
  pivot_longer(cols = c(total, with_interactions), 
               names_to = "category", 
               values_to = "count") %>%
  mutate(category = factor(category, levels = c("total", "with_interactions")))

# Panel B: Genera data
genera_data <- genera_by_taxon %>%
  left_join(genera_with_interactions, by = "taxon") %>%
  rename(total = n_sp, with_interactions = n_sp_interactions) %>%
  pivot_longer(cols = c(total, with_interactions), 
               names_to = "category", 
               values_to = "count") %>%
  mutate(category = factor(category, levels = c("total", "with_interactions")))

# Use same legend labels for both
category_labels <- c("total" = "Total", 
                     "with_interactions" = "With Interaction Data")

panel_a <- ggplot(species_data, aes(x = taxon, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_discrete(labels = category_labels) +
  scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = "Taxon", y = "Number of Species", fill = "")

panel_b <- ggplot(genera_data, aes(x = taxon, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_discrete(labels = category_labels) +
  scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  labs(x = "Taxon", y = "Number of Genera", fill = "")

legend <- get_legend(
  panel_b + theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0))
)

# Combine with shared legend at bottom
combined_plot <- plot_grid(
  plot_grid(panel_a, panel_b, labels = c("A", "B"), nrow = 2),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

ggsave("Figures/Supplementary/interactionRaw_summary.png", combined_plot, width = 10, height = 10, dpi = 300)

#### Summarise interaction data - % of data with spatial coordinates in CA, with coordinates but outside CA, and % without coordinates ####
# Check which points are in CA
interactions_spatial$in_CA <- lengths(st_intersects(interactions_spatial, ecoregions)) > 0

# Summarize
coord_summary <- data.frame(
  category = c("Coordinates in CA", "Coordinates outside CA", "No coordinates"),
  count = c(
    sum(interactions_spatial$in_CA),
    sum(!interactions_spatial$in_CA),
    sum(is.na(interactions$decimalLatitude) | is.na(interactions$decimalLongitude))
  )
) %>%
  mutate(percentage = round(count / nrow(interactions) * 100, 1))

ggplot(coord_summary, aes(x = "", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(aes(label = paste0(percentage, "%")), 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  labs(fill = "Location Category")

ggsave("Figures/Supplementary/interaction_spatial_summary.png", width = 10, height = 6, dpi = 300)



#### Summarize interaction data - % of data with both plant and pol species in checklist, % with only plant species in checklist, % with only pollinator species in checklist, and % with neither species in checklist (which means only genus is in checklist) ####
# Extract genus_species from interactions
interactions <- interactions %>%
  mutate(
    source_genus_species = word(sourceTaxonName_harm, 1, 2),
    target_genus_species = word(targetTaxonName_harm, 1, 2)
  )

# Check if species are in checklist
checklist_species <- unique(checklist_cleaned$genus_species)

interactions_checklist <- interactions %>%
  mutate(
    source_in_checklist = source_genus_species %in% checklist_species,
    target_in_checklist = target_genus_species %in% checklist_species,
    category = case_when(
      source_in_checklist & target_in_checklist ~ "Plant and Pollinator Species in Checklist",
      source_in_checklist & !target_in_checklist ~ "Pollinator Species in Checklist",
      !source_in_checklist & target_in_checklist ~ "Plant Species in Checklist",
      TRUE ~ "Only Plant/Pollinator Genus in Checklist"
    )
  )

# Summarize
checklist_summary <- interactions_checklist %>%
  group_by(category) %>%
  summarise(count = n()) %>%
  mutate(percentage = round(count / nrow(interactions) * 100, 1))

print(checklist_summary)

# Pie chart
ggplot(checklist_summary, aes(x = "", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(aes(label = paste0(percentage, "%")), 
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  labs(fill = "Checklist Status")

ggsave("Figures/Supplementary/interaction_taxonomicScale_summary.png", width = 10, height = 6, dpi = 300)

int_plantPol_unique <- interactions_full_tax %>% 
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm) %>% 
  unique()
message(paste("Total of", dim(int_plantPol_unique)[1], "unique plant-pol pairs"))


int_plantPol_speciesLevel_CA <- interactions_full_tax %>% 
  mutate(sourceTaxonName_harm_sp = word(sourceTaxonName_harm,1,2),
         targetTaxonName_harm_sp = word(targetTaxonName_harm,1,2)) %>%
  filter(sourceTaxonName_harm_sp %in% checklist_cleaned$genus_species & 
           targetTaxonName_harm_sp %in% checklist_cleaned$genus_species) %>%
  dplyr::select(sourceTaxonName_harm_sp, targetTaxonName_harm_sp) %>% 
  unique()
message(paste("Total of", dim(int_plantPol_speciesLevel_CA)[1], "unqiue plant-pol pairs native to CA at sp level in raw int data"))

int_plantPol_speciesLevel_CA %>%
  left_join(checklist_cleaned, by = join_by(sourceTaxonName_harm_sp == "genus_species")) %>%
  group_by(taxon) %>%
  summarise(n=n())

#### Predicted Interactions Summary ####
# how many distinct plant-pol-region combos?
predicted_interactions_plantpolregionCombos <- predicted_interactions %>% 
  filter(interactionTypeName == "visitsFlowersOf") %>% 
  filter(method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>% 
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, JEPCODE) %>% 
  distinct()

# how many unique plant-pol pairs from geoos?
predicted_interactions %>% 
  filter(interactionTypeName == "visitsFlowersOf") %>% 
  filter(method_type %in% c("Level2_GEOOS")) %>% 
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm) %>% 
  distinct() %>% 
  nrow()

nrow(predicted_interactions_plantpolregionCombos)

predicted_interactions_plantpolregionCombos %>% 
  left_join(checklist_cleaned, by = join_by("sourceTaxonName_harm" == "genus_species")) %>% 
  group_by(taxon) %>% 
  summarise(n=n())

# And for host
# how many distinct plant-pol-region combos?
predicted_interactions_plantpolregionCombos_hasHost <- predicted_interactions %>% 
  filter(interactionTypeName == "hasHost") %>% 
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, JEPCODE) %>% 
  distinct()

nrow(predicted_interactions_plantpolregionCombos_hasHost)

predicted_interactions_plantpolregionCombos_hasHost %>% 
  left_join(checklist_cleaned, by = join_by("sourceTaxonName_harm" == "genus_species")) %>% 
  group_by(taxon) %>% 
  summarise(n=n())

#### Metaweb degree vs. range size, raw occurrences, and raw interactions (flower visitation) ####
metaweb_flowervisitation <- predicted_interactions %>%
  filter(interactionTypeName == "visitsFlowersOf",
         method_type %in% c("Level1_Simple", "Level2_GEOOS")) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm) %>%
  distinct()

message(paste("Metaweb contains", nrow(metaweb_flowervisitation), "unique plant-pollinator pairs"))

## Total predicted degree per species (flower visitation)
pollinator_degree <- metaweb_flowervisitation %>%
  group_by(sourceTaxonName_harm) %>%
  summarise(degree = n_distinct(targetTaxonName_harm), .groups = "drop") %>%
  rename(genus_species = sourceTaxonName_harm)

plant_degree <- metaweb_flowervisitation %>%
  group_by(targetTaxonName_harm) %>%
  summarise(degree = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
  rename(genus_species = targetTaxonName_harm)

## Range size: sum of presence cells from binary SDMs
plant_binary_files <- list.files("Data_Clean/SDMs/sdm_by_taxon/plants/binary",
                                 pattern = "\\.tif$",
                                 full.names = TRUE)

plant_binary_files <- plant_binary_files[!grepl("manifest", plant_binary_files)]

plant_sdms_binary <- rast(plant_binary_files)

names(plant_sdms_binary) <- basename(plant_binary_files) %>%
  gsub("\\_binary.tif$", "", .) %>%
  gsub("_", " ", .)

plant_sdms_binary_filtered <- plant_sdms_binary[[names(plant_sdms_binary) %in% plant_species_in_checklist]]
all_pollinator_binary_files <- c()

for (taxon in pollinator_folders) {
  taxon_files <- list.files(
    file.path("Data_Clean/SDMs/sdm_by_taxon", taxon, "binary"),
    pattern = "\\.tif$",
    full.names = TRUE
  )
  taxon_files <- taxon_files[!grepl("manifest", taxon_files)]
  all_pollinator_binary_files <- c(all_pollinator_binary_files, taxon_files)
}

pollinator_sdms_binary <- rast(all_pollinator_binary_files)

names(pollinator_sdms_binary) <- basename(all_pollinator_binary_files) %>%
  gsub("\\_binary.tif$", "", .) %>%
  gsub("_", " ", .)

pollinator_sdms_binary_filtered <- pollinator_sdms_binary[[names(pollinator_sdms_binary) %in% pollinator_species_in_checklist]]

plant_range_size <- terra::global(plant_sdms_binary_filtered, "sum", na.rm = TRUE) %>%
  tibble::rownames_to_column("genus_species") %>%
  rename(range_size = sum)

pollinator_range_size <- terra::global(pollinator_sdms_binary_filtered, "sum", na.rm = TRUE) %>%
  tibble::rownames_to_column("genus_species") %>%
  rename(range_size = sum)

## Total # raw occurrences per species
occurrences_n <- occurrences %>%
  filter(genus_species %in% checklist_cleaned$genus_species) %>%
  count(genus_species, name = "n_occurrences")

## Total # raw interactions per species (flower visitation)
interactions_full_tax_nectar <- interactions_full_tax_nectar %>%
  mutate(sourceTaxonName_harm_gsp = word(sourceTaxonName_harm, 1, 2),
         targetTaxonName_harm_gsp = word(targetTaxonName_harm, 1, 2))

pollinator_n_interactions <- interactions_full_tax_nectar %>%
  count(sourceTaxonName_harm_gsp, name = "n_raw_interactions") %>%
  rename(genus_species = sourceTaxonName_harm_gsp)

plant_n_interactions <- interactions_full_tax_nectar %>%
  count(targetTaxonName_harm_gsp, name = "n_raw_interactions") %>%
  rename(genus_species = targetTaxonName_harm_gsp)

## Combine into one table per taxon group
pollinator_degree_data <- pollinator_degree %>%
  left_join(pollinator_range_size, by = "genus_species") %>%
  left_join(occurrences_n, by = "genus_species") %>%
  left_join(pollinator_n_interactions, by = "genus_species") %>%
  mutate(across(c(n_occurrences, n_raw_interactions), ~replace_na(., 0)))

plant_degree_data <- plant_degree %>%
  left_join(plant_range_size, by = "genus_species") %>%
  left_join(occurrences_n, by = "genus_species") %>%
  left_join(plant_n_interactions, by = "genus_species") %>%
  mutate(across(c(n_occurrences, n_raw_interactions), ~replace_na(., 0)))

pollinator_degree_data <- pollinator_degree_data %>% filter(!is.na(range_size))

plant_degree_data <- plant_degree_data %>% filter(!is.na(range_size))

## Scatterplot + intercept-forced-through-0 regression + correlation, for one x-y pair
degree_scatter_plot <- function(data, x_var, y_var, x_lab, y_lab) {
  plot_data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  cor_val <- cor(plot_data[[x_var]], plot_data[[y_var]], method = "pearson")
  ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("r = ", round(cor_val, 2)), size = 3.5) +
    labs(x = x_lab, y = y_lab) +
    theme_cowplot(font_size = 10)
}
## Pollinator degree
pollinator_degree_panels <- list(
  degree_scatter_plot(pollinator_degree_data, "range_size", "degree",
                      "Total Range Size (# cells)", "Predicted Pollinator Degree"),
  degree_scatter_plot(pollinator_degree_data, "n_occurrences", "degree",
                      "Total # Raw Occurrences", "Predicted Pollinator Degree"),
  degree_scatter_plot(pollinator_degree_data, "n_raw_interactions", "degree",
                      "Total # Raw Interactions", "Predicted Pollinator Degree")
)

pollinator_degree_plot <- wrap_plots(pollinator_degree_panels, nrow = 1) +
  plot_annotation(tag_levels = "a")

ggsave("Figures/network_metrics/pollinator_metaweb_degree_vs_datavolume.png",
       pollinator_degree_plot, width = 15, height = 5, dpi = 300)

## Plant degree: 3-panel plot
plant_degree_panels <- list(
  degree_scatter_plot(plant_degree_data, "range_size", "degree",
                      "Total Range Size (# cells)", "Predicted Plant Degree"),
  degree_scatter_plot(plant_degree_data, "n_occurrences", "degree",
                      "Total # Raw Occurrences", "Predicted Plant Degree"),
  degree_scatter_plot(plant_degree_data, "n_raw_interactions", "degree",
                      "Total # Raw Interactions", "Predicted Plant Degree")
)

plant_degree_plot <- wrap_plots(plant_degree_panels, nrow = 1)  +
  plot_annotation(tag_levels = "a")

ggsave("Figures/network_metrics/plant_metaweb_degree_vs_datavolume.png",
       plant_degree_plot, width = 15, height = 5, dpi = 300)

#### Metaweb degree vs. range size, raw occurrences, and raw interactions (Lepidopteran host interactions) ####
metaweb_hostinteractions <- predicted_interactions %>%
  filter(interactionTypeName == "hasHost",
         sourceTaxonType %in% c("moths", "butterflies")) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm) %>%
  distinct()
message(paste("Host-interaction metaweb contains", nrow(metaweb_hostinteractions), "unique Lepidoptera-plant pairs"))

## Total predicted degree per species (host interactions)
lepidoptera_degree <- metaweb_hostinteractions %>%
  group_by(sourceTaxonName_harm) %>%
  summarise(degree = n_distinct(targetTaxonName_harm), .groups = "drop") %>%
  rename(genus_species = sourceTaxonName_harm)

hostplant_degree <- metaweb_hostinteractions %>%
  group_by(targetTaxonName_harm) %>%
  summarise(degree = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
  rename(genus_species = targetTaxonName_harm)

## Total # raw interactions per species (host interactions)
interactions_full_tax_host <- interactions_full_tax_host %>%
  mutate(sourceTaxonName_harm_gsp = word(sourceTaxonName_harm, 1, 2),
         targetTaxonName_harm_gsp = word(targetTaxonName_harm, 1, 2))

lepidoptera_n_interactions <- interactions_full_tax_host %>%
  count(sourceTaxonName_harm_gsp, name = "n_raw_interactions") %>%
  rename(genus_species = sourceTaxonName_harm_gsp)

hostplant_n_interactions <- interactions_full_tax_host %>%
  count(targetTaxonName_harm_gsp, name = "n_raw_interactions") %>%
  rename(genus_species = targetTaxonName_harm_gsp)

## Combine with range size (from binary SDMs) and raw occurrences
lepidoptera_degree_data <- lepidoptera_degree %>%
  left_join(pollinator_range_size, by = "genus_species") %>%
  left_join(occurrences_n, by = "genus_species") %>%
  left_join(lepidoptera_n_interactions, by = "genus_species") %>%
  mutate(across(c(n_occurrences, n_raw_interactions), ~replace_na(., 0)))

hostplant_degree_data <- hostplant_degree %>%
  left_join(plant_range_size, by = "genus_species") %>%
  left_join(occurrences_n, by = "genus_species") %>%
  left_join(hostplant_n_interactions, by = "genus_species") %>%
  mutate(across(c(n_occurrences, n_raw_interactions), ~replace_na(., 0)))

lepidoptera_degree_data <- lepidoptera_degree_data %>% filter(!is.na(range_size))
hostplant_degree_data <- hostplant_degree_data %>% filter(!is.na(range_size))

## Lepidoptera degree
lepidoptera_degree_panels <- list(
  degree_scatter_plot(lepidoptera_degree_data, "range_size", "degree",
                      "Total Range Size (# cells)", "Predicted Lepidoptera Degree (Host)"),
  degree_scatter_plot(lepidoptera_degree_data, "n_occurrences", "degree",
                      "Total # Raw Occurrences", "Predicted Lepidoptera Degree (Host)"),
  degree_scatter_plot(lepidoptera_degree_data, "n_raw_interactions", "degree",
                      "Total # Raw Interactions", "Predicted Lepidoptera Degree (Host)")
)

lepidoptera_degree_plot <- wrap_plots(lepidoptera_degree_panels, nrow = 1)  +
  plot_annotation(tag_levels = "a")

ggsave("Figures/network_metrics/lepidoptera_metaweb_degree_vs_datavolume_host.png",
       lepidoptera_degree_plot, width = 15, height = 5, dpi = 300)

## Host plant degree
hostplant_degree_panels <- list(
  degree_scatter_plot(hostplant_degree_data, "range_size", "degree",
                      "Total Range Size (# cells)", "Predicted Host Plant Degree"),
  degree_scatter_plot(hostplant_degree_data, "n_occurrences", "degree",
                      "Total # Raw Occurrences", "Predicted Host Plant Degree"),
  degree_scatter_plot(hostplant_degree_data, "n_raw_interactions", "degree",
                      "Total # Raw Interactions", "Predicted Host Plant Degree")
)

hostplant_degree_plot <- wrap_plots(hostplant_degree_panels, nrow = 1)  +
  plot_annotation(tag_levels = "a")

ggsave("Figures/network_metrics/plant_metaweb_degree_vs_datavolume_host.png",
       hostplant_degree_plot, width = 15, height = 5, dpi = 300)

#### Calculate network metrics by ecoregion (all taxa combined) ####
compute_network_metrics_by_ecoregion <- function(interaction_data, taxa_groups, all_label = "All") {
  
  ecoregions_list <- unique(interaction_data$JEPCODE)
  network_metrics_results <- data.frame()
  
  # function for compute metrics from filtered dataset
  compute_metrics <- function(df, ecoregion_name, taxa_label) {
    if(nrow(df) < 3) return(NULL)
    
    tryCatch({
      mat <- df %>%
        distinct(targetTaxonName_harm, sourceTaxonName_harm) %>%
        mutate(n = 1) %>%
        pivot_wider(names_from = targetTaxonName_harm, values_from = n, values_fill = 0) %>%
        column_to_rownames("sourceTaxonName_harm") %>%
        as.matrix()
      
      if(nrow(mat) < 2 || ncol(mat) < 2) return(NULL)
      
      connectance     <- sum(mat > 0) / (nrow(mat) * ncol(mat))
      nodf_val        <- nested(mat, method = "NODF")
      
      # Mean degree (avg number of plant partners per pollinator)
      mean_degree     <- mean(rowSums(mat > 0))
      
      data.frame(
        ecoregion       = ecoregion_name,
        taxa            = taxa_label,
        n_pollinators   = nrow(mat),
        n_plants        = ncol(mat),
        n_interactions  = sum(mat),
        connectance     = connectance,
        nestedness_nodf = nodf_val,
        mean_degree     = mean_degree
      )
    }, error = function(e) {
      message(paste("Error in", ecoregion_name, taxa_label, ":", e$message))
      NULL
    })
  }
  
  print("Network properties")
  for(eco in ecoregions_list) {
    message(paste("Processing", eco))
    
    eco_data <- interaction_data %>% filter(JEPCODE == eco)
    
    if(nrow(eco_data) < 3) next
    
    ## Whole community (all taxa passed to this function)
    metrics_all <- compute_metrics(eco_data, eco, all_label)
    
    ## By taxa group
    metrics_taxa <- taxa_groups %>%
      purrr::map(function(tx) {
        tx_data <- eco_data %>% filter(sourceTaxonType == tx)
        compute_metrics(tx_data, eco, tx)
      })
    
    network_metrics_results <- bind_rows(
      network_metrics_results,
      metrics_all,
      metrics_taxa
    )
  }
  
  return(network_metrics_results)
}

## Flower visitation, all pollinator taxa 
network_metrics_results <- compute_network_metrics_by_ecoregion(
  interaction_data = predicted_interactions %>%
    filter(method_type %in% c("Level1_Simple", "Level2_GEOOS"),
           interactionTypeName == "visitsFlowersOf"),
  taxa_groups = c("bees", "hoverflies", "moths", "butterflies")
)
saveRDS(network_metrics_results, "Data_Clean/Analyses/network_metrics_results.rds")

## Lepidopteran host interactions only (moths + butterflies)
network_metrics_results_hostLepidoptera <- compute_network_metrics_by_ecoregion(
  interaction_data = predicted_interactions %>%
    filter(interactionTypeName == "hasHost",
           sourceTaxonType %in% c("moths", "butterflies")),
  taxa_groups = c("moths", "butterflies"),
  all_label = "Lepidoptera Host Use"
)

saveRDS(network_metrics_results_hostLepidoptera, "Data_Clean/Analyses/network_metrics_results_hostLepidoptera.rds")

##### Analyze network metrics // specialization across space #####
network_metrics_results <- readRDS("Data_Clean/Analyses/network_metrics_results.rds")
network_metrics_results_hostLepidoptera <- readRDS("Data_Clean/Analyses/network_metrics_results_hostLepidoptera.rds")

## Define metrics to plot
metrics_to_plot <- list(
  connectance     = "Connectance",
  nestedness_nodf = "Nestedness\n(NODF)",
  mean_degree     = "Mean\nDegree"
)

# Map function
make_eco_map <- function(data_filtered, metric_col, metric_name, panel_title,
                         shared_limits = NULL, show_legend = TRUE, use_log = FALSE) {
  
  plot_data <- ecoregions %>%
    left_join(data_filtered, by = c("JEPCODE" = "ecoregion"))
  
  if(use_log) {
    # Transform the fill values to log scale
    plot_data <- plot_data %>% mutate(fill_val = log(.data[[metric_col]]))
    fill_var  <- "fill_val"
    
    # Compute log-transformed limits
    scale_limits <- if(!is.null(shared_limits)) log(shared_limits) else NULL
    
    # Generate pretty breaks on original scale, then log-transform for positioning
    raw_breaks   <- pretty(shared_limits, n = 5)
    raw_breaks   <- raw_breaks[raw_breaks > 0]  # log needs positive values
    log_breaks   <- log(raw_breaks)
    
    fill_scale <- scale_fill_viridis_c(
      name   = metric_name,
      limits = scale_limits,
      breaks = log_breaks,
      labels = round(raw_breaks, 2),  # show original values as labels
      na.value = "grey90"
    )
  } else {
    fill_var   <- metric_col
    fill_scale <- scale_fill_viridis_c(
      name     = metric_name,
      limits   = shared_limits,
      na.value = "grey90"
    )
  }
  
  p <- ggplot(plot_data) +
    geom_sf(aes(fill = .data[[fill_var]]), color = "white", linewidth = 0.1) +
    fill_scale +
    labs(title = panel_title) +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.position = if(show_legend) "right" else "none"
    )
  
  return(p)
}

## Produce one output file per network metric, for a given results table/taxa set
plot_network_metric_maps <- function(network_metrics_results, taxa_groups, all_label = "All",
                                     output_suffix = "") {
  
  iwalk(metrics_to_plot, function(metric_name, metric_col) {
    
    message("Processing: ", metric_col, output_suffix)
    
    use_log <- metric_col == "nestedness_nodf"
    
    all_vals <- network_metrics_results %>%
      filter(!is.na(.data[[metric_col]]), .data[[metric_col]] > 0) %>%  # guard against log(0)
      pull(.data[[metric_col]])
    
    shared_limits <- range(all_vals, na.rm = TRUE)
    
    df_all <- network_metrics_results %>% filter(taxa == all_label)
    
    p_all <- make_eco_map(
      data_filtered = df_all,
      metric_col    = metric_col,
      metric_name   = metric_name,
      panel_title   = all_label,
      shared_limits = shared_limits,
      show_legend   = TRUE,
      use_log       = use_log
    )
    
    taxa_maps <- map(taxa_groups, function(tx) {
      df_tx <- network_metrics_results %>% filter(taxa == tx)
      
      make_eco_map(
        data_filtered = df_tx,
        metric_col    = metric_col,
        metric_name   = metric_name,
        panel_title   = tx,
        shared_limits = shared_limits,
        show_legend   = FALSE,
        use_log       = use_log
      )
    })
    
    bottom_row <- wrap_plots(taxa_maps, nrow = 1) +
      plot_layout(guides = "collect")
    
    final_plot <- p_all / bottom_row +
      plot_layout(heights = c(2, 1)) +
      plot_annotation(
        theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
      )
    
    ggsave(
      filename = paste0("Figures/network_metrics/ecoregion_", metric_col, "_by_taxon", output_suffix, ".png"),
      plot     = final_plot,
      width    = 10,
      height   = 8,
      dpi      = 300
    )
    
    cat("Saved:", metric_col, output_suffix, "\n")
  })
}

## Flower visitation maps
plot_network_metric_maps(
  network_metrics_results,
  taxa_groups = c("bees", "hoverflies", "moths", "butterflies"),
  all_label = "All",
  output_suffix = ""
)

## Lepidopteran host interaction maps
plot_network_metric_maps(
  network_metrics_results_hostLepidoptera,
  taxa_groups = c("moths", "butterflies"),
  all_label = "Lepidopteran Host Use",
  output_suffix = "_hostLepidoptera"
)

#### Mean pairwise species and interaction dissimilarity across regions ####
compute_beta_diversity_by_ecoregion <- function(int_data_long, output_suffix = "",
                                                midpoint_WN = 0.75, midpoint_ST = 0.75,
                                                midpoint_OS = 0.1) {
  
  # jepcodes
  regions <- unique(int_data_long$region)
  
  ### Convert to network matrices
  network_list <- list()
  for (r in regions) {
    
    cat("Processing region:", r, "\n")
    
    region_data <- int_data_long %>%
      filter(region == r)
    
    plants <- sort(unique(region_data$plant))
    pollinators <- sort(unique(region_data$pollinator))
    
    mat <- matrix(0,
                  nrow = length(plants),
                  ncol = length(pollinators),
                  dimnames = list(plants, pollinators))
    
    for (i in 1:nrow(region_data)) {
      plant_name <- region_data$plant[i]
      pollinator_name <- region_data$pollinator[i]
      mat[plant_name, pollinator_name] <- 1
    }
    
    network_list[[r]] <- mat
  }
  
  ### append regional matrices to make array
  all_plants <- sort(unique(unlist(lapply(network_list, rownames))))
  all_pollinators <- sort(unique(unlist(lapply(network_list, colnames))))
  cat("Total unique plants across all regions:", length(all_plants), "\n")
  cat("Total unique pollinators across all regions:", length(all_pollinators), "\n")
  
  n_regions <- length(network_list)
  webarray <- array(0,
                    dim = c(length(all_plants), length(all_pollinators), n_regions),
                    dimnames = list(all_plants, all_pollinators, names(network_list)))
  
  for (i in 1:n_regions) {
    region_name <- names(network_list)[i]
    mat <- network_list[[i]]
    webarray[rownames(mat), colnames(mat), i] <- mat
  }
  cat("Array dimensions:", dim(webarray), "(plants x pollinators x regions)\n")
  
  ### Calculate network dissimilarity
  results_commondenom <- betalinkr_multi(
    webarray = webarray,
    partitioning = "commondenom",
    binary = TRUE,
    index = "sorensen",
    distofempty = "zero"
  )
  saveRDS(results_commondenom, paste0("Data_Clean/Analyses/results_commondenom", output_suffix, ".rds"))
  
  ### average results across regions (mean per region)
  region_averages <- results_commondenom %>%
    pivot_longer(cols = c(i, j),
                 names_to = "position",
                 values_to = "region") %>%
    group_by(region) %>%
    summarise(
      mean_beta_WN = mean(WN),
      mean_beta_OS = mean(OS),
      mean_beta_ST = mean(ST),
      mean_prop_rewiring = mean(OS / WN),
      mean_prop_turnover = mean(ST / WN),
      n_comparisons = n()
    )
  
  ecoregions_beta <- ecoregions %>%
    left_join(region_averages, by = join_by("JEPCODE" == "region"))
  
  ### WN - network dissimilarity
  beta_WN_plot <- ggplot() +
    geom_sf(data = ecoregions_beta, aes(fill = mean_beta_WN), color = NA) +
    geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
    scale_fill_gradient2(
      low = "white",
      mid = "#F9B753",
      high = "#331E0F",
      midpoint = midpoint_WN,
      name = "Mean\nNetwork\nDissimilarity",
      na.value = "grey90"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(x = NULL, y = NULL) +
    coord_sf()
  ggsave(paste0("Figures/networks/beta_WN", output_suffix, ".pdf"),
         plot = beta_WN_plot,
         width = 200,
         height = 100,
         units = "mm",
         dpi = 2000)
  
  ### beta ST - Species turnover
  beta_ST_plot <- ggplot() +
    geom_sf(data = ecoregions_beta, aes(fill = mean_beta_ST), color = NA) +
    geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
    scale_fill_gradient2(
      low = "white",
      mid = "#F9B753",
      high = "#331E0F",
      midpoint = midpoint_ST,
      name = "Contribution of\nSpecies Turnover",
      na.value = "grey90"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(x = NULL, y = NULL) +
    coord_sf()
  ggsave(paste0("Figures/networks/beta_ST", output_suffix, ".pdf"),
         plot = beta_ST_plot,
         width = 200,
         height = 100,
         units = "mm",
         dpi = 2000)
  
  ### beta OS - interaction rewiring
  beta_OS_plot <- ggplot() +
    geom_sf(data = ecoregions_beta, aes(fill = mean_beta_OS), color = NA) +
    geom_sf(data = ecoregions_boundary, fill = NA, color = "black", linewidth = 1) +
    scale_fill_gradient2(
      low = "white",
      mid = "#F9B753",
      high = "#331E0F",
      midpoint = midpoint_OS,
      name = "Contribution of\nInteraction Rewiring",
      na.value = "grey90"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(x = NULL, y = NULL) +
    coord_sf()
  ggsave(paste0("Figures/networks/beta_OS", output_suffix, ".pdf"),
         plot = beta_OS_plot,
         width = 200,
         height = 100,
         units = "mm",
         dpi = 2000)
  
  return(list(results_commondenom = results_commondenom,
              region_averages = region_averages,
              plots = list(WN = beta_WN_plot, ST = beta_ST_plot, OS = beta_OS_plot)))
}

## Flower visitation
int_data_long <- predicted_interactions %>%
  filter(method_type %in% c("Level1_Simple", "Level2_GEOOS"),
         interactionTypeName == "visitsFlowersOf") %>%
  dplyr::select(
    plant = targetTaxonName_harm,
    pollinator = sourceTaxonName_harm,
    region = JEPCODE
  )

beta_diversity_flowervisitation <- compute_beta_diversity_by_ecoregion(
  int_data_long,
  output_suffix = ""
)

## Lepidopteran host interactions  
int_data_long_host <- predicted_interactions %>%
  filter(interactionTypeName == "hasHost",
         sourceTaxonType %in% c("moths", "butterflies")) %>%
  dplyr::select(
    plant = targetTaxonName_harm,
    pollinator = sourceTaxonName_harm,
    region = JEPCODE
  )

beta_diversity_hostLepidoptera <- compute_beta_diversity_by_ecoregion(
  int_data_long_host,
  output_suffix = "_hostLepidoptera"
)

##### Plot Xerces digitized boundaries ####
xerces_regions <- xerces_regions %>%
  mutate(xerces_region = str_squish(str_replace_all(xerces_region, "(?<=[a-z])(?=[A-Z])", " ")))

xerces_map <- ggplot() +
  geom_sf(data = xerces_regions, aes(fill = xerces_region),
          color = "black", linewidth = 0.4, alpha = 0.8) +
  geom_sf(data = CA_boundary, fill = NA, color = "black", linewidth = 1) +
  scale_fill_brewer(palette = "Set2", name = "Xerces Region") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("Figures/Supplementary/xerces_regions_map.pdf",
       plot  = xerces_map,
       width = 150, height = 150, units = "mm", dpi = 300)
