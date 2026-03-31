#### Load Packages ####
library(tidyverse)

#Note: since this is just the alpha summaries at the community level, for now we are keeping pollen and nectar interactions separate for bees
# Note: In version 1 of the calscape data, we counted "Confirmed" as only interactions that are confirmed in that eocregion. Now we are separating confirmed in the ecoregion from confirmed anywhwere, and Calscape can display those levels now and group them in any way they want
# Note: For now, still decided to count host and flower visitation for same species as two separate events, so int he tallies that would get 2 points. If the plant is serving two functions, should be prioritized more, is the reasoning. 

#### Read in Data ####
# Predicted interactions dataset
predicted_interactions <- read.csv("Data_Clean/Interactions/Predicted_Interactions.csv")

#For now, filter for only simple interaction predictions (more conservative approach)
predicted_interactions <- predicted_interactions %>%
  filter(
    (interactionTypeName == "hasHost" & method_type == "host") |
      (interactionTypeName == "visitsFlowersOf" & method_type == "Level1_Simple")
  )

#Add a column designating whether the interaction has been observed anywhere (Right now, our confirmed column is really for "confirmed in that ecoregion")

#read raw interactions
raw_ints = read.csv("Data_Clean/Interactions/ints_final_clean_withCrops_slim.csv")

raw_ints = raw_ints %>%
  select(sourceTaxonName_harm, interactionTypeName, targetTaxonName_harm) %>%
  mutate(interactionTypeName = ifelse(interactionTypeName == "collectsPollenOf",
                                      "visitsFlowersOf",
                                      interactionTypeName)) %>%
  distinct()

raw_ints$observedInteractionAnywhere = 1

predicted_interactions = left_join(predicted_interactions, raw_ints) 

predicted_interactions$observedInteractionAnywhere = ifelse(is.na(predicted_interactions$observedInteractionAnywhere), 0, 1)

#Now change categories in categoricalInteraction column
predicted_interactions <- predicted_interactions %>%
  mutate(
    categoricalInteraction = case_when(
      observedInteraction == 1L & observedInteractionAnywhere == 1L ~ "Confirmed_Ecoregion",
      observedInteraction == 0L & observedInteractionAnywhere == 1L ~ "Confirmed_General",
      TRUE ~ "Potential"
    )
  )

#### Create calscape summary table ####
# Get number of interactions per plant+ecoregion, broken up by confirmed vs potential interactions
## Spatialized by Jepson ecoregion, keep taxon seperate
alpha_spatial_byType <- predicted_interactions %>%
  group_by(JEPCODE, sourceTaxonType, targetTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  rename(region = JEPCODE,
         pollinatorGroup = sourceTaxonType,
         plantSpecies = targetTaxonName_harm)

## Spatialized by Jepson ecoregion, all taxa together
alpha_spatial_all <- predicted_interactions %>%
  group_by(JEPCODE, targetTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(pollinatorGroup = "all") %>%
  rename(region = JEPCODE,
         plantSpecies = targetTaxonName_harm)

## Statewide, keep taxon seperate
alpha_state_byType <- predicted_interactions %>%
  # First aggregate interactions to state scale
  group_by(targetTaxonName_harm, sourceTaxonName_harm, sourceTaxonType, interactionTypeName) %>%
  summarize(
    categoricalInteraction = ifelse(
      any(categoricalInteraction %in% c("Confirmed_Ecoregion", "Confirmed_General")), # modify to only "Confirmed_Ecoregion" if more conservative
      "Confirmed_General", 
      "Potential"
    ),
    .groups = "drop"
  ) %>%
  # now summarise number of interactions per pollinator taxa
  group_by(sourceTaxonType, targetTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(region = "California") %>%
  rename(pollinatorGroup = sourceTaxonType,
         plantSpecies = targetTaxonName_harm)

## Statewide, all taxa together
alpha_state_all <- predicted_interactions %>%
  # First aggregate interactions to state scale
  group_by(targetTaxonName_harm, sourceTaxonName_harm, sourceTaxonType, interactionTypeName) %>%
  summarize(
    categoricalInteraction = ifelse(
      any(categoricalInteraction %in% c("Confirmed_Ecoregion", "Confirmed_General")), # modify to only "Confirmed_Ecoregion" if more conservative
      "Confirmed_General", 
      "Potential"
    ),
    .groups = "drop"
  ) %>%
  # now summarise number of interactions across pollinator taxa
  group_by(targetTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(region = "California",
         pollinatorGroup = "all") %>%
  rename(plantSpecies = targetTaxonName_harm)

## Combine all scales, reformat
alpha_all <- bind_rows(alpha_spatial_all, alpha_spatial_byType, alpha_state_all, alpha_state_byType)

# Rearrange columns
alpha_all <- alpha_all %>%
  select(region, pollinatorGroup, plantSpecies, interactionTypeName, Confirmed_Ecoregion, Confirmed_General, Potential)

# Sort values
alpha_all <- alpha_all[order(alpha_all$region, alpha_all$pollinatorGroup, -(alpha_all$Confirmed_Ecoregion + alpha_all$Confirmed_General + alpha_all$Potential)), ]

## Add total confirmed column
alpha_all$Confirmed_Total <- rowSums(
  alpha_all[, c("Confirmed_Ecoregion",
                "Confirmed_General")],
  na.rm = TRUE
)

# Add total Sum column
alpha_all$Sum <- rowSums(
  alpha_all[, c("Confirmed_Ecoregion",
                "Confirmed_General",
                "Potential")],
  na.rm = TRUE
)


### Same, but for pollinators now

# Get number of interactions per plant+ecoregion, broken up by confirmed vs potential interactions
## Spatialized by Jepson ecoregion, keep taxon seperate
alpha_spatial_byType_poll <- predicted_interactions %>%
  group_by(JEPCODE, sourceTaxonType, sourceTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  rename(region = JEPCODE,
         pollinatorGroup = sourceTaxonType,
         pollinatorSpecies = sourceTaxonName_harm)

## Spatialized by Jepson ecoregion, all taxa together
alpha_spatial_all_poll <- predicted_interactions %>%
  group_by(JEPCODE, sourceTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(pollinatorGroup = "all") %>%
  rename(region = JEPCODE,
         pollinatorSpecies = sourceTaxonName_harm)

## Statewide, keep taxon seperate
alpha_state_byType_poll <- predicted_interactions %>%
  # First aggregate interactions to state scale
  group_by(targetTaxonName_harm, sourceTaxonName_harm, sourceTaxonType, interactionTypeName) %>%
  summarize(
    categoricalInteraction = ifelse(
      any(categoricalInteraction %in% c("Confirmed_Ecoregion", "Confirmed_General")), # modify to only "Confirmed_Ecoregion" if more conservative
      "Confirmed_General", 
      "Potential"
    ),
    .groups = "drop"
  ) %>%
  # now summarise number of interactions per pollinator taxa
  group_by(sourceTaxonType, sourceTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(region = "California") %>%
  rename(pollinatorGroup = sourceTaxonType,
         pollinatorSpecies = sourceTaxonName_harm)

## Statewide, all taxa together
alpha_state_all_poll <- predicted_interactions %>%
  # First aggregate interactions to state scale
  group_by(targetTaxonName_harm, sourceTaxonName_harm, sourceTaxonType, interactionTypeName) %>%
  summarize(
    categoricalInteraction = ifelse(
      any(categoricalInteraction %in% c("Confirmed_Ecoregion", "Confirmed_General")), # modify to only "Confirmed_Ecoregion" if more conservative
      "Confirmed_General", 
      "Potential"
    ),
    .groups = "drop"
  ) %>%
  # now summarise number of interactions across pollinator taxa
  group_by(sourceTaxonName_harm, interactionTypeName, categoricalInteraction) %>%
  summarize(count = n(), .groups = 'drop') %>%
  # Pivot wider to have columns for each interaction type and confirmation status combination
  pivot_wider(
    names_from = categoricalInteraction,
    values_from = count,
    values_fill = 0
  ) %>% 
  mutate(region = "California",
         pollinatorGroup = "all") %>%
  rename(pollinatorSpecies = sourceTaxonName_harm)

## Combine all scales, reformat
alpha_all_poll <- bind_rows(alpha_spatial_all_poll, alpha_spatial_byType_poll, alpha_state_all_poll, alpha_state_byType_poll)

# Rearrange columns
alpha_all_poll <- alpha_all_poll %>%
  select(region, pollinatorGroup, pollinatorSpecies, interactionTypeName, Confirmed_Ecoregion, Confirmed_General, Potential)

# Sort values
alpha_all_poll <- alpha_all_poll[order(alpha_all_poll$region, alpha_all_poll$pollinatorGroup, -(alpha_all_poll$Confirmed_Ecoregion + alpha_all_poll$Confirmed_General + alpha_all_poll$Potential)), ]

## Add total confirmed column
alpha_all_poll$Confirmed_Total <- rowSums(
  alpha_all_poll[, c("Confirmed_Ecoregion",
                "Confirmed_General")],
  na.rm = TRUE
)

# Add total Sum column
alpha_all_poll$Sum <- rowSums(
  alpha_all_poll[, c("Confirmed_Ecoregion",
                "Confirmed_General",
                "Potential")],
  na.rm = TRUE
)

#remove all column
alpha_all_poll = alpha_all_poll %>%
  filter(pollinatorGroup != "all")


## Save
#Summarized interactions
write.csv(alpha_all, "Data_Calscape/Interactionscalscape_interactions_allPredicted_summarised.csv")
write.csv(alpha_all_poll, "Data_Calscape/Interactions/calscape_interactions_allPredicted_pollinators_summarised.csv")

## Metadata

# The only thing that has changed since last time is that now instead of just 1 Confirmed columns, we have split this into two levels: Confirmed_Exoregion, and Confirmed_General:

# Confirmed_Ecoregion: We have a georeferenced data point indicating that this exact interaction occurs in this exact ecoregion (note that will always be NA when region = California since all the statewide summaries are inherently not considering the ecoregion)
# Confirmed_General: When the region is one of the 35 Jepson ecoregions (not California), this column refers to interactions where the two interacting species co-occur spatially and temporally in the ecoregion, and we have a data point for the interaction, but that data point is not in the focal ecoregion. So in other words, the interaction itself is confirmed, we just don't know for sure that it occurs in this ecoregion. When the region = California, this column refers to all of the confirmed interactions that we have data points for (agnostic of ecoregion, since its at the state level)
# Potential: All the predicted interactions from our model, that we don't have any raw data points for, but have high probability of occurring based on our modeling. 

#Reformat the raw predicted interactions
predicted_interactions = predicted_interactions %>%
  select(JEPCODE, sourceTaxonType, interactionTypeName, sourceTaxonName_harm, targetTaxonName_harm, categoricalInteraction) %>%
  rename(
    taxon = sourceTaxonType,
    higher = sourceTaxonName_harm,
    lower = targetTaxonName_harm
  )

#Flag pollen interactions
pollen = read.csv('Data_Clean/Interactions/ints_final_clean_full.csv')

pollen = pollen %>%
 filter(source %in% c("fowler", "bigbook"), interactionTypeName == "collectsPollenOf") %>%
  select(
    source, higherGenus, sourceTaxonName_harm, interactionTypeName,
    lowerGenus, targetTaxonName_harm
  ) %>%
  rename(
    sourceTaxonGenusName = higherGenus,
    targetTaxonGenusName = lowerGenus
  ) %>%
  mutate(across(everything(), ~ na_if(., ""))) %>%   # convert "" → NA
  distinct()

# If both interaction types exist for the same sourceTaxonName_harm × targetTaxonName_harm combination, retain only the specialist record and remove the generalist duplicate.
pollen <- pollen %>%
  mutate(
    interactionTypeName = ifelse(
      source == "fowler",
      "collectsPollenOf_specialist",
      "collectsPollenOf_general"
    )
  ) %>%
  group_by(sourceTaxonName_harm, targetTaxonName_harm) %>%
  filter(
    !(interactionTypeName == "collectsPollenOf_general" &
        any(interactionTypeName == "collectsPollenOf_specialist"))
  ) %>%
  ungroup() %>%
  select(-source) %>%
  distinct()

#First populate all known species-level associations
pollen1 = pollen %>%
  select(sourceTaxonName_harm, targetTaxonName_harm, interactionTypeName) %>% distinct() %>%
  rename(
    higher = sourceTaxonName_harm,
    lower = targetTaxonName_harm
  )

#replace interactionTypeName with the value in pollen1 for matching rows
predicted_interactions <- predicted_interactions %>%
  rows_update(
    pollen1,
    by = c("higher", "lower"),
    unmatched = "ignore"
  )

#Now populate the specialist interactions at the genus level
pollen2 = pollen %>%
  filter(interactionTypeName == "collectsPollenOf_specialist") %>%
  select(sourceTaxonName_harm, targetTaxonGenusName, interactionTypeName) %>% distinct() %>%
  rename(
    higher = sourceTaxonName_harm,
    lower_genus = targetTaxonGenusName
  )

predicted_interactions <- predicted_interactions %>%
  
  # Extract genus from lower (first word)
  mutate(lower_genus = word(lower, 1)) %>%
  
  # Join on higher + genus
  left_join(
    pollen2,
    by = c("higher", "lower_genus"),
    suffix = c("", "_pollen")
  ) %>%
  
  # Overwrite interaction type when pollen2 provides one
  mutate(
    interactionTypeName = coalesce(interactionTypeName_pollen,
                                   interactionTypeName)
  ) %>%
  
  # Clean up helper columns
  select(-lower_genus, -interactionTypeName_pollen)

#Save all predicted Interactions with pollen categories separated
write.csv(predicted_interactions, "Data_Calscape/Interactions/calscape_interactions_allPredicted_pollenSeparated.csv")

#now remove pollen separation
predicted_interactions = predicted_interactions %>%
  mutate(interactionTypeName = ifelse(interactionTypeName %in% c("collectsPollenOf_general", "collectsPollenOf_specialist"),
                                      "collectPollenOf",
                                      interactionTypeName)) %>%
  distinct()


#Add estimates of how common the pollinators are seen visiting plants, for Calscape to show common representatives of each taxon
pollRank = alpha_all_poll %>% 
  select(region, pollinatorGroup, pollinatorSpecies, interactionTypeName, Sum) %>%
  rename(JEPCODE = region,
         taxon = pollinatorGroup,
         higher = pollinatorSpecies) %>%
  distinct() %>%
  group_by(JEPCODE, taxon) %>%
  mutate(pollRank = dense_rank(desc(Sum))) %>%
  select(-Sum) %>%
  ungroup()

predicted_interactions = left_join(predicted_interactions, pollRank)

#Save all predicted Interactions
write.csv(predicted_interactions, "Data_Calscape/Interactions/calscape_interactions_allPredicted.csv")


#Metadata
#The only things that has changed since last time are: 

#(1) The categoricalInteraction column, which now has 3 levels instead of the previous 2 (Confirmed and Potential): 

# Confirmed_Ecoregion: We have a georeferenced data point indicating that this exact interaction occurs in this exact ecoregion (note that will always be NA when region = California since all the statewide summaries are inherently not considering the ecoregion)
# Confirmed_General: When the region is one of the 35 Jepson ecoregions (not California), this column refers to interactions where the two interacting species co-occur spatially and temporally in the ecoregion, and we have a data point for the interaction, but that data point is not in the focal ecoregion. So in other words, the interaction itself is confirmed, we just don't know for sure that it occurs in this ecoregion. When the region = California, this column refers to all of the confirmed interactions that we have data points for (agnostic of ecoregion, since its at the state level)
# Potential: All the predicted interactions from our model, that we don't have any raw data points for, but have high probability of occurring based on our modeling. 

#(2) The interactionTypeName column, which used to have 2 levels (hasHost and visitsFlowersOf) now has 3 levels:

#hasHost: Lepidoptera (butterflies and moths) larvel host plant interactions
#visitsFlowersOF: General flower visitation interactions across all taxa
#collectsPollenOf: Pollen collection interactions

# The "calscape_interactions_allPredicted_pollenSeparated.csv" file also has:
#collectsPollenOf_general: General pollen collection interactions, which only applies to bees on our dataset. Note that some of the pollen collection interactions in this column may represent specialist (oligolectic) bee species, but we do not currently have confirmation for that in our dataset, and it is just as likely that they are non-specialized interactions 
#collectsPollenOf_specialist: Pollen collection interactions involving specialist (oligolectic) bee species. These are the interactions we know, based on our data, are specific to these specialized bees. Note that pollen specialist bee species may still have visitsFlowersOf interactions, which may represent resource use outside of their pollen specialist plants (e.g., nectaring)

# prepare separate dataset of just pollen specialist bees
pollen = read.csv("Data_Calscape/Interactions/calscape_interactions_allPredicted_pollenSeparated.csv") %>%
  filter(interactionTypeName == "collectsPollenOf_specialist") %>%
  select(higher) %>%
  rename(pollinatorSpecies = higher) %>%
  distinct() %>%
  arrange(pollinatorSpecies)

write.csv(pollen, "Data_Calscape/Interactions/pollenSpecialistBees.csv", row.names = F)
