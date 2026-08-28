# Morpho
# Chris Cosma
# 3/27/2025
# Final crop datasets

#### Part 1: Generate dataset of crop interactions at ecoregion and state level ####

#NOTE: the crop interactions did not go through the same interaction prediction framework as the native plants because we do not have reliable distribution or phenology data for the crops yet. So our approach right now is to assume 1) that anyone can choose to plant any crop in any ecoregion and 2) interactions between a crop and pollinator species is assumed to occur in the ecoregion if a) the pollinator species occurs there and b) we have a raw datapoint for the interaction between that pollinator and the crop 

#packages
library(tidyverse)

#load all final predicted interactions (exluding crops)
allInts = read.csv("Data_Calscape/Interactions/calscape_interactions_allPredicted.csv")

#load intermediate dataset that includes crop interactions (note that it does not matter that the plant names have not been standardized back to the Calscape names, because we are not using native plants here.)
cropInts = read.csv("Data_Clean/Interactions/ints_final_clean_withCrops_slim.csv")

#isolate all crop interactions, and keep only flower visitation
cropInts = filter(cropInts, isCrop == "yes" & interactionTypeName %in% c("visitsFlowersOf", "collectsPollenOf"))

#take only needed columns and remove duplicates
cropInts = cropInts %>%
  select(targetTaxonName_harm, sourceTaxonName_harm, interactionTypeName, cropCommonName) %>%
  rename(
    lower = targetTaxonName_harm,
    higher =sourceTaxonName_harm
  ) %>%
  distinct()

#filter out non animal-pollinated crops (these were the only two that Claude identified as having absolutely no evidence of benefit from insects for fruit/seed production)
cropInts = cropInts %>%
  filter((!cropCommonName %in% c("Corn", "Barley")))

#create a list of unique pollinators per ecoregion
ecoPoll = allInts %>%
  select(JEPCODE, taxon, higher) %>%
  distinct()

#join to crop interactions based on only shared values of "higher"
cropInts = inner_join(cropInts, ecoPoll)

#reorder
cropInts = cropInts %>%
  select(JEPCODE, taxon, higher, interactionTypeName, lower, cropCommonName)

#rename JEPCODE
names(cropInts)[1] = "region"

#add in all interactions at state level
ints = unique(cropInts[,-1])
ints$region = "California"

#join ecoregion and state level
cropInts = rbind(cropInts, ints)

write.csv(cropInts, "Data_Calscape/Interactions/calscape_interactions_crops.csv")

#### Part 2: Generate crop support plant dataset ####
##create ecoregion by crop lists of native support plants

#packages
library(tidyverse)

#load crop interactions
cropInts = read.csv("Data_Calscape/Interactions/calscape_interactions_crops.csv")

#load native plant interactions (harmonized to Calscape names)
nativeInts = read.csv('Data_Calscape/Interactions/calscape_interactions_allPredicted.csv')

#rename ecoregion column
names(nativeInts)[2] = "region"

#rename all pollen collection to visitsFlowersOf (note that this means that a plant that serves as both visitsFlowersOf and collectsPollenOf for a species will not get extra points. For now too complicated to think about those details, but preserve distinction for host and nectar plant for leps)
nativeInts <- nativeInts %>%
  mutate(
    interactionTypeName = ifelse(
      interactionTypeName %in% c("collectsPollenOf_specialist",
                                 "collectsPollenOf_general",
                                 "collectsPollenOf"),
      "visitsFlowersOf",
      interactionTypeName
    )
  )

#add state level interaction to nativeInts
ints = nativeInts %>% 
  select(-X, -region) %>%
  distinct()

#add region = California
ints$region = "California"

#join
nativeInts = bind_rows(nativeInts, ints)

#remove categoricalInteraction and interactionTypeName columns (for now we don't care how the native plant supports the crop pollinators, acting as a host plant for a crop pollinator is equally valid) but don't distinct after removing interactionTypeName, because should still get extra points for acting as host and nectar 
nativeInts = nativeInts %>%
  select(-categoricalInteraction) %>%
  distinct() %>%
  select(-interactionTypeName)

#NOTE: right now there is an issue with the var. and ssp. because we assumed they have all the same interactions and occur wherever the regular species occurr. Opting for just removing them right now because they are flooding the lists

#for now, also remove var. and ssp.
nativeInts <- nativeInts[sapply(strsplit(nativeInts$lower, "\\s+"), length) <= 2, ]

#unique crop latin and common names
crop_names = cropInts %>%
  select(lower, cropCommonName) %>%
  distinct()

#create list of unique regions
eco = unique(nativeInts$region)

#blank dataframe to fill
crop_support_plants = data.frame()

#loop through to create the full dataset
for (i in 1:length(eco)) {
  
  #filter for desired ecoregion
  dat = filter(cropInts, region == eco[[i]])
  
  #make a list of unique crop plants in that ecoregion (for now every crop is in every ecoregion)
  crops = unique(dat$lower)
  
  #loop through each crop plant
  for (j in 1:length(crops)) {
    
    #make a list of pollinators of that crop in that ecoregion
    crop_polls = unique(filter(dat, lower == crops[[j]])$higher)
    
    #filter the native plant dataset for interaction involving just these pollinator species in that ecoregion
    dat2 = filter(nativeInts, region == eco[[i]] & higher %in% crop_polls)
    
    #calculate alpha diversity by pollinator taxon
    alpha = dat2 %>%
      group_by(lower, taxon) %>%
      summarise(alpha = n())
    
    #calculate alpha diversity across all taxa
    alpha_all = dat2 %>%
      group_by(lower) %>%
      summarise(alpha = n())
    
    #add descriptor for all
    alpha_all$taxon = "all"
    
    #join the datasets
    alpha = rbind(alpha, alpha_all)
    
    #add region information
    alpha$region = eco[[i]]
    
    #add crop latin name
    alpha$cropLatinName = crops[[j]]
    
    #add crop common names
    alpha = left_join(alpha, crop_names, by = c("cropLatinName" = "lower"))
    
    #reorder and arrange
    alpha = alpha %>%
      select(region, cropLatinName, cropCommonName, taxon, lower, alpha) %>%
      arrange(region, cropLatinName, taxon, -alpha)
    
    #compile into final dataset
    crop_support_plants = rbind(crop_support_plants, alpha)
    
  }
  
}

#create dataset of total pollinator species supported per ecoregion for native plants
scores = nativeInts %>%
  group_by(region,taxon,lower) %>%
  summarise(sum = n())

#create category for all taxa
scores2 = nativeInts %>%
  group_by(region, lower) %>%
  summarise(sum = n())

#add taxon = all
scores2$taxon = "all"

#join
scores = rbind(scores,scores2)

#add to crop dataframe
crop_support_plants = left_join(crop_support_plants, scores)

#filter for just top 20 for each ecoregion + crop for a more manageable dataset (can change this to whatever number we want). In the case of ties, break the tie with the pollinator score (sum) column
calscape_cropSupportPlants20 <- crop_support_plants %>%
  group_by(region, cropLatinName, cropCommonName, taxon) %>%
  arrange(desc(alpha), desc(sum)) %>% 
  mutate(rank = row_number()) %>%
  slice_head(n = 20) %>%
  ungroup()

#remove unneeded column
calscape_cropSupportPlants20$sum = NULL

names(calscape_cropSupportPlants20)[c(5,6)] = c("nativeSupportPlant", "cropPollinatorsSupported")

#save final dataset
write.csv(calscape_cropSupportPlants20, "Data_Calscape/Crop_Support/calscape_cropSupportPlants20.csv")



