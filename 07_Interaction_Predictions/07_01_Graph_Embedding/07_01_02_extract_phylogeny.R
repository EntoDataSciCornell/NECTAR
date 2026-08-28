library(data.table)
library(rotl)
library(tidyverse)
library(ape)

#read in the interaction data
interaction_df <- fread("Data_Clean/Interactions/ints_final_clean_slim.csv")%>%
  filter(interactionTypeName %in% c("collectsPollenOf", "visitsFlowersOf"))

#read in the checklist
checklist_df<- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")


###After harmonization, we read in the harmonized list and compare it with the species name on the interaction data and the checklist

group_name <- "butterflies" #change names for other taxon groups

if (group_name == "bees" | group_name == "hoverflies"){
harmonized_names <- fread(paste0("Data_Raw/Phylogeny/", group_name,"_names_tree_harm.csv"))%>%
select(name_from_tree, names_harm)}else{harmonized_names<-fread("Data_Raw/Phylogeny/leps_names_tree_harm.csv")%>%
 select(name_from_tree, names_harm)}

group_genera<-checklist_df%>%
  filter(taxon== group_name)%>%
  select(genus)%>%
  unique()%>%
  pull()

checklist_sp<-checklist_df%>%
  filter(taxon==group_name)%>%
  select(genus_species)%>%
  unique()%>%
  arrange()%>%
  pull()

interaction_sp<-interaction_df%>%
  filter(sourceTaxonGenusName %in% group_genera)%>%
  mutate(genus_species = word(sourceTaxonName_harm, 1,2))%>%
  select(genus_species)%>%unique()%>%arrange()%>%filter(!is.na(genus_species))%>%
  pull()

####Subset the phylogentic tree from Open Tree of Life####

##create a dataframe with species names and source to input into the tnrs_match_names() function
sp_names_to_input <- data.frame(
  species = c(checklist_sp, interaction_sp), 
  source = c(rep("checklist", length(checklist_sp)), 
             rep("interaction_data", length(interaction_sp)))) %>%
  arrange(species, source) %>%  # This puts "checklist" before "interaction_data" alphabetically
  distinct(species, .keep_all = TRUE)%>%
  left_join(harmonized_names, join_by(species == names_harm))

#double check to make sure any names that overlap between the interaction data and the checklist has the source listed as "checklist"
sp_names_to_input%>%
  filter(species == intersect(interaction_sp, checklist_sp)[1])

#final species list to be input into tnrs_match_names()
poll_species <- sp_names_to_input%>%
  filter(!is.na(name_from_tree))%>%
  select(name_from_tree)%>%
  pull()

##Now we run tnrs_match_names() to get the species id for extracting the phylogentic tree
#we split into chunks of 1000 names to avoid overloading
chunk_size <- 1000
poll_matches_list <- list()

for(i in seq(1, length(poll_species), by = chunk_size)) {
  chunk <- poll_species[i:min(i + chunk_size - 1, length(poll_species))]
  poll_matches_list[[length(poll_matches_list) + 1]] <- tnrs_match_names(names = chunk)
  Sys.sleep(5)
}

poll_matches <- do.call(rbind, poll_matches_list)


#check to see if any genus has no matching id in the system
sum(is.na(poll_matches$ott_id))

#remove rows that has no match (id = NAs) in the system
good_matches <- poll_matches[!is.na(poll_matches$ott_id),]


#make sure there is no duplicated id 
sum(duplicated(good_matches$ott_id))

##Check which OTT IDs are actually in the tree (this step filters out problematic IDs)
ott_ids_to_check <- good_matches$ott_id

# Test IDs in small batches to find problematic ones
valid_id_df <- data.frame() #tracking the valid IDs
invalid_id_df<-data.frame() #tracking the problematic IDs

for(i in 1:length(ott_ids_to_check)) {
  
  test_id <- ott_ids_to_check[i]
  # Try to get a tree with just this ID and one known good ID
  if (group_name == "bees"){
    ref_id = 3272412
  }else if (group_name == "hoverflies"){
    ref_id = 956387
  }else if (group_name == "butterflies"){
    ref_id = 532411
  }else if (group_name == "moths"){
    ref_id = 7142984}
  
  try_tree <- try(tol_induced_subtree(ott_ids = c(ref_id, test_id)), silent = TRUE)
  
  if(!inherits(try_tree, "try-error")) {

    valid_df<-data.frame(species = good_matches$unique_name[good_matches$ott_id == test_id], 
      id= test_id, 
      score = good_matches$score[good_matches$ott_id == test_id],
      name_from_tree = good_matches$search_string[good_matches$ott_id == test_id])

    valid_id_df <- rbind(valid_id_df, valid_df)
  } else {

    cat("Saving problematic ID:", test_id, "-", good_matches$unique_name[good_matches$ott_id == test_id], "\n")
    
    invalid_df<-data.frame(species = good_matches$unique_name[good_matches$ott_id == test_id], 
      id= test_id, 
      name_from_tree = good_matches$search_string[good_matches$ott_id == test_id])
    
    invalid_id_df <- rbind(invalid_id_df, invalid_df)
    
  }
}


#now create a dataframe that link the id back to the names of checklist
data_to_match<-valid_id_df%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(name_from_tree)), join_by(name_from_tree == name_from_tree))%>%
  select(id, score, species.y)%>%#the harmonized names
  group_by(species.y)%>%
  slice_max(score, n =1, with_ties= FALSE) #remove duplicates, save the one with the higher matching scores


#then we build the tree with only valid IDs 
poll_tree <- tol_induced_subtree(ott_ids = data_to_match$id, label_format ="id")


# Now we match the tip labels with the name_from_tree using the ott_id
for(i in 1:length(poll_tree$tip.label)) {
  
  # Extract OTT ID from tip label
  ott_id_str <- sub(".*ott", "", poll_tree$tip.label[i])  
  
  if(ott_id_str %in% data_to_match$id) {
    
    # Match it with the original name
    original_name <- data_to_match$species.y[data_to_match$id == ott_id_str]
    
    poll_tree$tip.label[i]<-original_name
    
    print(poll_tree$tip.label[i])
  }
}


#Add branch lengths
final_tree <- compute.brlen(poll_tree, method = "Grafen", power = 0.5)

#test<-read.tree(paste0("Data_Clean/", group_name, "_tree.nwk"))

#save the tree
write.tree(final_tree, paste0("Data_Clean/Graph_Embedding/Phylogeny/", group_name, "_tree.nwk"))


## Create a summary table to summarize species that were not on the tree

#names that are only at genus level is not included
missing_genera<-interaction_df%>%
  filter(sourceTaxonGenusName %in% group_genera)%>%
  mutate(genus_species = word(sourceTaxonName_harm, 1,2))%>%
  filter(is.na(genus_species))%>%
  select(sourceTaxonName_harm)%>%
  unique()%>%
  rename(species = sourceTaxonName_harm)%>%
  mutate(source = "interaction_data",
         reason = "name at genus level")

#names that has no match with the harmonized taxonomy
sp_with_no_match<-sp_names_to_input%>%
  filter(is.na(name_from_tree))%>%
  select(species, source)%>%
  mutate(reason = "no match in the harmonized taxonomy")
  

#names that are not avaliable on the tree
names_not_on_the_tree<-invalid_id_df%>%
  select(name_from_tree)%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(name_from_tree)))%>%
  select(species, source)%>%
  mutate(reason = "names not avaliable on the tree")%>%
  distinct()

  
#the summary dataframe
missing_names_summary<-rbind(missing_genera, sp_with_no_match, names_not_on_the_tree)%>%
  distinct()

#Save the summary
write.csv(missing_names_summary, paste0("Data_Clean/Graph_Embedding/Phylogeny/", group_name, "_sp_without_phylogeny_summary.csv"), row.names = FALSE)







