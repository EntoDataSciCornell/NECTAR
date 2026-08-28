#### Load Packages ####
library(dplyr)

#### For Calscape Pheno Display ####
## read in data
# Phenometrics
phen_best <- read.csv("Data_Clean/Phenology/phenology_best_available.csv")
phen_simple <- read.csv("Data_Clean/Phenology/phenology_simple.csv")
phen_type2 <- read.csv("Data_Clean/Phenology/phenology_type2.csv")
phen_type3 <- read.csv("Data_Clean/Phenology/phenology_type3.csv")

# checklists
checklists <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# # read in predicted interactions, use to remove pollinators not modelled
# predicted_interactions <- read.csv("Data_Clean/Interactions/Predicted_Interactions.csv")

## supporting vectors
month_names <- c("January", "Febuary", "March", "April", "May", "June", 
                 "July", "August", "September", "October", "November", "December")
month_breaks <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

## Add confidence flag
# For phen_best(regional phenology), flag if pheno estimated at scale other then JEPCODE
# For phen_simple, flag if estimated using Type2 or Type3 method
phen_best$confidence_flag <- case_when(
  phen_best$phen_scale == "JEPCODE" ~ NA,
  phen_best$phen_scale == "REGION" ~ "dataPooled",
  phen_best$phen_scale == "PROV" ~ "dataPooled",
  phen_best$phen_scale == "CALIFORNIA" ~ NA,
  phen_best$phen_scale == "TYPE2" ~ "methodSimple",
  phen_best$phen_scale == "TYPE3" ~ "methodSimple",
  .default = "?"
)

phen_simple$confidence_flag <- NA
phen_type2$confidence_flag <- "methodSimple"
phen_type3$confidence_flag <- "methodSimple"

## Reformat phenology
# remove NA
phen_best <- phen_best[, c("genus_species", "JEPCODE", "phen_start", "phen_end", "confidence_flag")]
phen_best <- phen_best[!is.na(phen_best$genus_species) & 
                         !is.na(phen_best$phen_start) & 
                         !is.na(phen_best$phen_end), ]

phen_simple <- phen_simple[, c("genus_species", "phen_start", "phen_end", "confidence_flag")]
phen_simple <- phen_simple[!is.na(phen_simple$genus_species) & 
                             !is.na(phen_simple$phen_start) & 
                             !is.na(phen_simple$phen_end), ]

phen_type2 <- phen_type2[, c("genus_species", "phen_start", "phen_end", "confidence_flag")]
phen_type2 <- phen_type2[!is.na(phen_type2$genus_species) & 
                           !is.na(phen_type2$phen_start) & 
                           !is.na(phen_type2$phen_end), ]

phen_type3 <- phen_type3[, c("genus_species", "phen_start", "phen_end", "confidence_flag")]
phen_type3 <- phen_type3[!is.na(phen_type3$genus_species) & 
                           !is.na(phen_type3$phen_start) & 
                           !is.na(phen_type3$phen_end), ]

# Merge simple, type2, and type3 (all are at same stat-level resolution)
phen_simple <- rbind(phen_simple, phen_type2, phen_type3)

# constrain to [0, 364]
phen_simple$phen_start[phen_simple$phen_start > 364] <- 364
phen_simple$phen_end[phen_simple$phen_end > 364] <- 364
phen_best$phen_start[phen_best$phen_start > 364] <- 364
phen_best$phen_end[phen_best$phen_end > 364] <- 364

# add taxa column
phen_simple <- inner_join(phen_simple,
                          unique(select(checklists, genus_species, taxon)), 
                          by = "genus_species")
phen_best <- inner_join(phen_best, 
                        unique(select(checklists, genus_species, taxon)), 
                        by = "genus_species")

## Pheno reformatting
get_month <- function(days) {
  findInterval(days, month_breaks)
}

phenomatrix <- function(phenoframe_onset, phenoframe_offset, month_names){
  # Create onset matrix (1s from start month to December)
  onset_matrix <- matrix(0, nrow = length(phenoframe_onset), ncol = 12)
  start_months <- get_month(phenoframe_onset)
  for(i in 1:nrow(onset_matrix)) {
    onset_matrix[i, start_months[i]:12] <- 1
  }
  
  # Create end matrix (1s from January to end month)
  end_matrix <- matrix(0, nrow = length(phenoframe_offset), ncol = 12)
  end_months <- get_month(phenoframe_offset)
  for(i in 1:nrow(end_matrix)) {
    end_matrix[i, 1:end_months[i]] <- 1
  }
  
  # Multiply matrices element-wise to get overlap
  phenomatrix <- onset_matrix * end_matrix
  colnames(phenomatrix) <- month_names
  
  return(phenomatrix)
}

# Calscape phenoframes
phen_simple_mat <- phenomatrix(phen_simple$phen_start, phen_simple$phen_end, month_names)
phen_simple_calscape <- cbind(phen_simple, phen_simple_mat)

phen_best_mat <- phenomatrix(phen_best$phen_start, phen_best$phen_end, month_names)
phen_best_calscape <- cbind(phen_best, phen_best_mat)

## Add/Rename JEPCODE to region
phen_best_calscape <- phen_best_calscape %>%
  rename(region = JEPCODE)
phen_simple_calscape$region <- "CALIFORNIA"

## Merge
phen_all_calscape <- bind_rows(phen_best_calscape, phen_simple_calscape)

# Save
write.csv(phen_all_calscape, "Data_Calscape/Phenology/calscape_phenology.csv")
