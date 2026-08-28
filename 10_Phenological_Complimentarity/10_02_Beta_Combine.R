# Combine the companion plant datasets for Calscape

#### Load Packages ####
library(tidyverse)

#### Load Data ####
phen = read.csv("Data_Calscape/Companion_Plants/community_constrained_phenological_complimentarity.csv")
nect = read.csv("Data_Calscape/Companion_Plants/community_constrained_companions_visitsFlowersOf.csv")
host = read.csv("Data_Calscape/Companion_Plants/community_constrained_companions_hasHost.csv")

#Load interactions file to make total_pollinators the total number of pollinators supported by the focal plant + the companion plant for each pair
ints = read.csv('Data_Calscape/Interactions/calscape_interactions_allPredicted.csv')

#### Combine Complimentary Plant Methods into One File ####
#create dataframes from original interaction data to identify number of pollinator species supported by each plant in each region

#nectar first 

unique(ints$interactionTypeName)

#jepcode scale
ints_nect_jepcode = ints %>%
  filter(interactionTypeName != "hasHost") %>%
  group_by(JEPCODE, lower) %>%
  summarise(total_polls = n_distinct(higher), .groups = "drop")

#state scale and join
ints_nect = ints %>%
  filter(interactionTypeName != "hasHost") %>%
  group_by(lower) %>%
  summarise(total_polls = n_distinct(higher), .groups = "drop") %>%
  mutate(JEPCODE = "California") %>%
  bind_rows(ints_nect_jepcode)

#now host 

#jepcode scale
ints_host_jepcode = ints %>%
  filter(interactionTypeName == "hasHost") %>%
  group_by(JEPCODE, lower) %>%
  summarise(total_polls = n_distinct(higher), .groups = "drop")

#state scale and join
ints_host = ints %>%
  filter(interactionTypeName == "hasHost") %>%
  group_by(lower) %>%
  summarise(total_polls = n_distinct(higher), .groups = "drop") %>%
  mutate(JEPCODE = "California") %>%
  bind_rows(ints_host_jepcode)

#join back to beta diversity dataframes

#nectar

#first for focal plant
nect1 = left_join(
  nect,
  ints_nect,
  by = c("focal_plant" = "lower", "JEPCODE" = "JEPCODE")
)

#now for companion plant
nect1 = left_join(
  nect1,
  ints_nect,
  by = c("companion_plant" = "lower", "JEPCODE" = "JEPCODE")
)

#To check that the numbers are adding up
#reorder and rename
# nect1 = nect1 %>%
#   mutate(category = "diversity_nectar") %>%
#   select(JEPCODE, category, focal_plant, companion_plant, total_polls.x, total_polls.y, added_pollinators, total_pollinators, rank) %>%
#   rename(region = JEPCODE,
#          total_focal_nectar = total_polls.x,
#          total_companion_nectar = total_polls.y,
#          added_by_companion_nectar = added_pollinators,
#          sum_focal_companion_nectar = total_pollinators,
#          companion_rank_nectar = rank)

nect1 = nect1 %>%
  mutate(category = "diversity_nectar") %>%
  select(JEPCODE, category, focal_plant, companion_plant, added_pollinators, total_pollinators, rank) %>%
  rename(region = JEPCODE,
         added_polls_nect = added_pollinators,
         total_polls_nect = total_pollinators,
         rank_nect = rank)

# host

#first for focal plant
host1 = left_join(
  host,
  ints_host,
  by = c("focal_plant" = "lower", "JEPCODE" = "JEPCODE")
)

#now for companion plant
host1 = left_join(
  host1,
  ints_host,
  by = c("companion_plant" = "lower", "JEPCODE" = "JEPCODE")
)

#reorder and rename
# host1 = host1 %>%
#   mutate(category = "diversity_host") %>%
#   select(JEPCODE, category, focal_plant, companion_plant, total_polls.x, total_polls.y, added_pollinators, total_pollinators, rank) %>%
#   rename(region = JEPCODE,
#          total_focal_host = total_polls.x,
#          total_companion_host = total_polls.y,
#          added_by_companion_host = added_pollinators,
#          sum_focal_companion_host = total_pollinators,
#          companion_rank_host = rank)

host1 = host1 %>%
  mutate(category = "diversity_host") %>%
  select(JEPCODE, category, focal_plant, companion_plant, added_pollinators, total_pollinators, rank) %>%
  rename(region = JEPCODE,
         added_polls_host = added_pollinators,
         total_polls_host = total_pollinators,
         rank_host = rank)

### Now, get number of unique pollinators added (irregardless of interaction type) by new plant
ints_jepcode = ints %>%
  select(JEPCODE, lower, higher) %>%
  unique()

# Based on whole state
ints_jepcode_state = ints %>%
  select(lower, higher) %>%
  unique() %>%
  mutate(JEPCODE = "California")

#join together ecoregion level and state level total pol data
ints_jepcode_full = bind_rows(ints_jepcode, ints_jepcode_state)

# Split into "focal" and "companion" versions for joining
focal_polls = ints_jepcode_full %>%
  rename(focal_plant = lower, higher_focal = higher)

companion_polls = ints_jepcode_full %>%
  rename(companion_plant = lower, higher_companion = higher)

# All focal-companion-region combos we need to compute this for
pairs = bind_rows(
  nect %>% select(JEPCODE, focal_plant, companion_plant),
  host %>% select(JEPCODE, focal_plant, companion_plant)
) %>% unique()

# For each pair, get every pollinator the companion supports in that region,
# flag whether the focal plant already supports it too, then summarize
added_total = pairs %>%
  # add pollinators supported by companion plant
  left_join(companion_polls, by = c("JEPCODE", "companion_plant")) %>%
  # flag which of the companion's pollinators are also supported by the focal plant
  left_join(
    focal_polls %>% mutate(in_focal = TRUE),
    by = c("JEPCODE", "focal_plant", "higher_companion" = "higher_focal")
  ) %>%
  # group and summarize: total pollinators supported by companion, and how many are new (not already covered by focal)
  group_by(JEPCODE, focal_plant, companion_plant) %>%
  summarise(
    total_pollinators = n_distinct(higher_companion),
    added_pollinators = n_distinct(higher_companion[is.na(in_focal)]),
    .groups = "drop"
  ) %>%
  # rank companions within each focal plant + region by pollinators added
  group_by(JEPCODE, focal_plant) %>%
  mutate(rank = dense_rank(desc(added_pollinators))) %>%
  ungroup()


#phenology 

phen1 = phen %>%
  mutate(category = "phenology") %>%
  rename(rank_pheno = rank,
         region = JEPCODE) %>%
  select(region, category, focal_plant, companion_plant, new_days, total_days, rank_pheno)


#Join them all

# Full join all three by the shared keys
combined = nect1 %>%
  full_join(host1, by = c("region", "focal_plant", "companion_plant")) %>%
  full_join(phen1, by = c("region", "focal_plant", "companion_plant"))

# Rename for consistency with the rest of the method
overall = added_total %>%
  rename(region = JEPCODE,
         added_polls_overall = added_pollinators,
         total_polls_overall = total_pollinators,
         rank_overall_diversity = rank)

combined = combined %>%
  left_join(overall, by = c("region", "focal_plant", "companion_plant"))


# Collapse the category columns into one comma-separated field
combined = combined %>%
  unite("category", category.x, category.y, category, sep = ", ", na.rm = TRUE) %>%
  # clean up any empty strings from all-NA rows
  mutate(category = na_if(category, ""))

combined %>%
  count(category) %>%
  mutate(proportion = n / sum(n))

demo = combined %>%
  group_by(category) %>%
  slice_sample(n = 3) %>%
  ungroup()

write.csv(demo, "Data_Calscape/Companion_Plants/companion_demo.csv", row.names = F)
write.csv(combined, "Data_Calscape/Companion_Plants/calscape_companions.csv", row.names = F)


# Metadata for combined companion plant file
# -------------------------------------------------
# region:            Jepson ecoregion code, or "California" for state-level
# category:          Comma-separated list of companion criteria met (diversity_nectar, diversity_host, phenology)
# focal_plant:       Focal plant species
# companion_plant:   Recommended companion plant
# added_polls_nect:  Additional nectar/pollen visitor species gained by adding the companion to the focal plant
# total_polls_nect:  Total nectar/pollen visitor species supported by the focal + companion pair
# rank_nect:         Companion rank for nectar/pollen visitor diversity (1 = best)
# added_polls_host:  Additional caterpillar species gained by adding the companion to the focal plant
# total_polls_host:  Total caterpillar species supported by the focal + companion pair
# rank_host:         Companion rank for host diversity (1 = best)
# new_days:          Additional days of floral resource availability added by the companion
# total_days:        Total days of floral resource availability for the focal + companion pair
# rank_pheno:        Companion rank for phenological complementarity (1 = best)
# added_polls_overall: Additional pollinator species (any interaction type) gained by adding the companion
# total_polls_overall: Total pollinator species (any interaction type) supported by the companion in that region
# rank_overall_diversity: Companion rank for overall pollinator diversity gained (among choices already identified), regardless of interaction type
# -------------------------------------------------
# Columns are NA when the companion was not identified under that criterion.
