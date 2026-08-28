#### Load packages ####
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(data.table)

#### Supporting scripts ####
source("Code/12_Analyses/12_00_LoadData.R")
source("Code/99_Supporting/genetic_algorithm.R")

#### Aggregate PredictedInteractions, save for network visualization ####
tax <- checklist_cleaned %>%
  dplyr::select(genus_species, taxon) %>%
  unique()

# Filter predicted_interactions for flower visitation with Level1_Simple or Level2_GEOOS
predicted_interactions_filtered <- predicted_interactions[
  interactionTypeName == "visitsFlowersOf" &
    method_type %in% c("Level1_Simple", "Level2_GEOOS")
][, `:=`(targetTaxonGenusName = sub(" .*", "", targetTaxonName_harm),
         sourceTaxonGenusName = sub(" .*", "", sourceTaxonName_harm))]

# Raw nectar interactions
interactions_nectar <- interactions_nectar %>%
  left_join(tax, by = join_by(sourceTaxonGenusName == genus_species)) %>%
  rename(sourceTaxonType = taxon)

# Filter raw interactions for CA-only nectar interactions
interactions_nectar_CA <- interactions_nectar %>%
  filter(sourceTaxonName_harm %in% checklist_cleaned[checklist_cleaned$taxon != "plants", ]$genus_species |
           targetTaxonName_harm %in% checklist_cleaned[checklist_cleaned$taxon == "plants", ]$genus_species) %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm) %>%
  unique()

#### Create node file (genus level) ####
plant_nodes <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm) %>%
  group_by(targetTaxonGenusName) %>%
  summarise(
    n_species = n_distinct(targetTaxonName_harm),
    taxon = "plants",
    .groups = "drop"
  ) %>%
  rename(ID = targetTaxonGenusName)

pollinator_nodes <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm, sourceTaxonType) %>%
  group_by(sourceTaxonGenusName, sourceTaxonType) %>%
  summarise(n_species = n_distinct(sourceTaxonName_harm), .groups = "drop") %>%
  rename(ID = sourceTaxonGenusName, taxon = sourceTaxonType)

node_file <- bind_rows(plant_nodes, pollinator_nodes) %>%
  arrange(taxon, ID)

### Degree
pred_degree_pollinators <- predicted_interactions_filtered %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName) %>%
  unique() %>%
  group_by(sourceTaxonGenusName) %>%
  summarise(pred_degree = n()) %>%
  rename(ID = sourceTaxonGenusName)

pred_degree_pollinators_region <- predicted_interactions_filtered %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName, JEPCODE) %>%
  unique() %>%
  group_by(sourceTaxonGenusName, JEPCODE) %>%
  summarise(pred_degree = n()) %>%
  rename(ID = sourceTaxonGenusName) %>%
  pivot_wider(names_from = JEPCODE, values_from = pred_degree, names_prefix = "pred_", values_fill = 0)

pred_degree_plants <- predicted_interactions_filtered %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName) %>%
  unique() %>%
  group_by(targetTaxonGenusName) %>%
  summarise(pred_degree = n()) %>%
  rename(ID = targetTaxonGenusName)

pred_degree_plants_region <- predicted_interactions_filtered %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName, JEPCODE) %>%
  unique() %>%
  group_by(targetTaxonGenusName, JEPCODE) %>%
  summarise(pred_degree = n()) %>%
  rename(ID = targetTaxonGenusName) %>%
  pivot_wider(names_from = JEPCODE, values_from = pred_degree, names_prefix = "pred_", values_fill = 0)

obs_degree_pollinators <- interactions_nectar_CA %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName) %>%
  unique() %>%
  group_by(sourceTaxonGenusName) %>%
  summarise(obs_degree = n()) %>%
  rename(ID = sourceTaxonGenusName)

obs_degree_plants <- interactions_nectar_CA %>%
  dplyr::select(targetTaxonName_harm, sourceTaxonName_harm, targetTaxonGenusName, sourceTaxonGenusName) %>%
  unique() %>%
  group_by(targetTaxonGenusName) %>%
  summarise(obs_degree = n()) %>%
  rename(ID = targetTaxonGenusName)

node_file <- node_file %>%
  left_join(bind_rows(pred_degree_pollinators, pred_degree_plants)) %>%
  left_join(bind_rows(obs_degree_pollinators, obs_degree_plants)) %>%
  left_join(bind_rows(pred_degree_pollinators_region, pred_degree_plants_region)) %>%
  mutate(obs_degree = replace_na(obs_degree, 0))

#### Create edge file (genus level) ####
edge_CA <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm) %>%
  unique() %>%
  group_by(sourceTaxonGenusName, targetTaxonGenusName) %>%
  summarise(
    interactionDensityCA = n_distinct(paste(sourceTaxonName_harm, targetTaxonName_harm, sep = "_")),
    .groups = "drop"
  )

edge_raw <- interactions_nectar_CA %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm) %>%
  unique() %>%
  group_by(sourceTaxonGenusName, targetTaxonGenusName) %>%
  summarise(
    interactionDensityRaw = n_distinct(paste(sourceTaxonName_harm, targetTaxonName_harm, sep = "_")),
    .groups = "drop"
  )

edge_by_ecoregion_predicted <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm, JEPCODE) %>%
  unique() %>%
  group_by(sourceTaxonGenusName, targetTaxonGenusName, JEPCODE) %>%
  summarise(
    interaction_density = n_distinct(paste(sourceTaxonName_harm, targetTaxonName_harm, sep = "_")),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = JEPCODE, values_from = interaction_density, names_prefix = "pred_", values_fill = 0)

edge_by_ecoregion_confirmed <- predicted_interactions_filtered %>%
  filter(categoricalInteraction == "Confirmed") %>%
  dplyr::select(sourceTaxonGenusName, sourceTaxonName_harm, targetTaxonGenusName, targetTaxonName_harm, JEPCODE) %>%
  unique() %>%
  group_by(sourceTaxonGenusName, targetTaxonGenusName, JEPCODE) %>%
  summarise(
    interaction_density = n_distinct(paste(sourceTaxonName_harm, targetTaxonName_harm, sep = "_")),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = JEPCODE, values_from = interaction_density, names_prefix = "conf_", values_fill = 0)

edge_file <- edge_CA %>%
  left_join(edge_raw, by = c("sourceTaxonGenusName", "targetTaxonGenusName")) %>%
  left_join(edge_by_ecoregion_predicted, by = c("sourceTaxonGenusName", "targetTaxonGenusName")) %>%
  left_join(edge_by_ecoregion_confirmed, by = c("sourceTaxonGenusName", "targetTaxonGenusName")) %>%
  replace_na(list(interactionDensityRaw = 0)) %>%
  mutate(across(starts_with("pred_") | starts_with("conf_"), ~replace_na(.x, 0)))

edge_file <- edge_file %>%
  mutate(
    n_ecoregions_predicted = rowSums(dplyr::select(., starts_with("pred_")) > 0),
    n_ecoregions_confirmed = rowSums(dplyr::select(., starts_with("conf_")) > 0),
    has_confirmation       = interactionDensityRaw > 0
  )

edge_file <- edge_file %>%
  rename(Source = sourceTaxonGenusName, Target = targetTaxonGenusName) %>%
  dplyr::select(
    Source, Target,
    interactionDensityRaw, interactionDensityCA,
    n_ecoregions_predicted, n_ecoregions_confirmed, has_confirmation,
    starts_with("pred_"),
    starts_with("conf_")
  )

node_file_reduced <- node_file %>%
  dplyr::select(ID, n_species, taxon, pred_degree, obs_degree, pred_SCo, pred_CaRF)
edge_file_reduced <- edge_file %>%
  dplyr::select(Source, Target, interactionDensityRaw, interactionDensityCA,
                n_ecoregions_predicted, n_ecoregions_confirmed, pred_SCo, pred_CaRF, conf_SCo, conf_CaRF)

#### Export ####
write.csv(node_file_reduced, "Data_Clean/Analyses/gephi_nodes_genus.csv", row.names = FALSE)
write.csv(edge_file_reduced, "Data_Clean/Analyses/gephi_edges_genus.csv", row.names = FALSE)

#### Taxon pie chart for predicted interactions ####
interaction_taxon_summary <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, sourceTaxonType) %>%
  distinct() %>%
  group_by(sourceTaxonType) %>%
  summarise(n = n()) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

# [Pie chart plot and find_topMix() call for SCo/CaRF continue here unchanged,
#  except find_topMix() is now sourced from genetic_algorithm.R]

focalPollinators <- checklist_cleaned %>%
  filter(taxon %in% c("bees", "butterflies", "moths", "hoverflies"))

topComplimentary <- lapply(c("SCo", "CaRF"), function(ecoregion) {
  find_top_mix(
    ecoregion,
    interaction_data_subset = predicted_interactions_filtered,
    candidate_plants        = unique(plant_checklist$genus_species),
    focal_pollinators       = unique(focalPollinators$genus_species),
    checklist               = checklist_cleaned,
    sp_by_region            = sp_by_region,
    n_plants                = 6,
    n_generations           = 200,
    population_size         = 200,
    mutation_rate           = 0.1
  )
})

topComplimentary_taxon_summary <- do.call(rbind, lapply(topComplimentary, function(x) x$taxon_summary))

region_plant_rank <- predicted_interactions_filtered %>%
  group_by(targetTaxonName_harm, JEPCODE) %>%
  summarise(n_supported_region=n())
# 
# topComplimentary[[1]]$top_mix %>%
#   left_join(region_plant_rank, 
#             by = join_by(ecoregion == JEPCODE, 
#                          scientificName == targetTaxonName_harm)) %>%
#   group_by(ecoregion) %>%
#   arrange(desc(n_supported_region)) %>%
#   View()

prepare_pollinator_plot_data <- function(ecoregion_code, 
                                         optimized_results) {
  
  # Just prepare optimized results
  pollinator_data <- optimized_results %>%
    filter(ecoregion == ecoregion_code) %>%
    mutate(
      taxon_group = tools::toTitleCase(taxon),
      proportion = total_in_region / sum(total_in_region),
      end_angle = cumsum(proportion) * 2 * pi,
      start_angle = lag(end_angle, default = 0),
      mid_angle = (start_angle + end_angle) / 2
    )
  
  return(pollinator_data)
}


plot_pollinator_support <- function(pollinator_data, ecoregion_name) {
  
  tick_positions <- seq(0, 100, by = 25)
  tick_labels <- c(0, 25, 50, 75, 100)
  
  colors_optimized <- c(
    "Bees" = "#FFCC00",
    "Butterflies" = "#E69190",
    "Moths" = "#BD9ECD",
    "Hoverflies" = "#F76700",
    "Hummingbirds" = "#52A273"
  )
  
  ggplot(pollinator_data) +
    
    # Radial grid lines
    geom_segment(
      data = data.frame(x = tick_positions / 100),
      aes(x = x, xend = x, y = 0, yend = 2 * pi),
      color = "grey70",
      linewidth = 0.2,
      alpha = 0.5
    ) +
    
    # Optimized scenario (solid colors)
    geom_rect(
      aes(xmin = 0, xmax = percent_supported/100, 
          ymin = start_angle, ymax = end_angle,
          fill = taxon_group),
      alpha = 1,
      color = "white",
      linewidth = 0.3
    ) +
    
    # Percentage labels
    geom_text(
      data = data.frame(
        x = tick_labels / 100,
        label = paste0(tick_labels, "%")
      ),
      aes(x = x, y = 0, label = label),
      size = 4 / .pt,
      color = "black",
      vjust = -0.5,
      family = "Arial"
    ) +
    
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors_optimized) +
    
    theme_void() +
    theme(
      legend.position = "none",
      text = element_text(family = "Arial", size = 4),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    xlim(0, 1.2)
}

# Make plots
SCo_data <- prepare_pollinator_plot_data("SCo",
                                         topComplimentary_taxon_summary)
plot_SCo <- plot_pollinator_support(SCo_data, "South Coast (SCo)")

CaRF_data <- prepare_pollinator_plot_data("CaRF",
                                          topComplimentary_taxon_summary)
plot_CaRF <- plot_pollinator_support(CaRF_data, "Central Coast Ranges (CaRF)")

# Save
ggsave("Figures/networks/SCo_complimentary.png", 
       plot = plot_SCo,
       width = 50, 
       height = 50, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/networks/SCo_complimentary.pdf", 
       plot = plot_SCo,
       width = 50, 
       height = 50, 
       units = "mm",
       device = cairo_pdf)

ggsave("Figures/networks/CaRF_complimentary.png", 
       plot = plot_CaRF,
       width = 50, 
       height = 50, 
       units = "mm",
       dpi = 4000)

ggsave("Figures/networks/CaRF_complimentary.pdf", 
       plot = plot_CaRF,
       width = 50, 
       height = 50, 
       units = "mm",
       device = cairo_pdf)



#### Taxon pie charts: % of interactions by pollinator taxon ####

# Raw (observed) flower visitation interactions, by taxon
raw_taxon_summary <- interactions_nectar_CA %>%
  left_join(tax %>% mutate(genus = word(genus_species, 1, 1)), by = c("sourceTaxonGenusName" = "genus")) %>%
  distinct(sourceTaxonName_harm, targetTaxonName_harm, taxon) %>%
  count(taxon, name = "n") %>%
  mutate(
    taxon_label = tools::toTitleCase(taxon),
    percentage = round(n / sum(n) * 100, 1)
  )

# Predicted flower visitation interactions, by taxon
predicted_taxon_summary <- predicted_interactions_filtered %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, sourceTaxonType) %>%
  distinct() %>%
  count(sourceTaxonType, name = "n") %>%
  mutate(
    taxon_label = tools::toTitleCase(sourceTaxonType),
    percentage = round(n / sum(n) * 100, 1)
  )

taxon_colors <- c(
  "Bees" = "#FFCC00",
  "Butterflies" = "#E69190",
  "Moths" = "#BD9ECD",
  "Hoverflies" = "#F76700",
  "Hummingbirds" = "#52A273"
)

plot_taxon_pie <- function(data) {
  ggplot(data, aes(x = "", y = n, fill = taxon_label)) +
    geom_col(color = "white", linewidth = 0.5) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(name = "Taxon", values = taxon_colors) +
    theme_void() +
    theme(
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom"
    )
}

raw_taxon_pie <- plot_taxon_pie(raw_taxon_summary)
predicted_taxon_pie <- plot_taxon_pie(predicted_taxon_summary)

ggsave("Figures/taxon_pie/raw_taxon_pie.png", raw_taxon_pie, width = 6, height = 6, dpi = 300)
ggsave("Figures/taxon_pie/predicted_taxon_pie.png", predicted_taxon_pie, width = 6, height = 6, dpi = 300)

ggsave("Figures/taxon_pie/raw_taxon_pie.pdf", raw_taxon_pie, width = 6, height = 6, device = cairo_pdf)
ggsave("Figures/taxon_pie/predicted_taxon_pie.pdf", predicted_taxon_pie, width = 6, height = 6, device = cairo_pdf)
