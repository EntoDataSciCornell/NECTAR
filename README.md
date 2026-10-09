# Purpose
This repository hosts an example implementation of the **Network-Enhanced Conservation Tool for Analysis and Recommendation (NECTAR)** framework, as descibed in "Local interaction networks reconstructed from global biodiversity data improve pollinator restoration decision making" by Baiotto and Cosma et al. (2026). 

NECTAR is a reproducible framework for predicting large-scale plant-pollinator interaction metawebs using existing interaction data and phylogenetic information, and downscaling these large-scale metawebs to regional metawebs using spatial co-occurrence and the overlap of key life-stage events (e.g., for flower visitation interactions, regional flowering phenology and pollinator adult flight period). While the code below is our application of NECTAR to California plant-pollinator communities, it can be extended to other regions and pollinator taxa with relative ease by re-configuring the extent of input occurrence, interaction, and environmental data, the spatial boundaries used to downscale to regional metawebs, and the "last say" authority for name harmonization in the region of interest (e.g., for plants in California, we use Jepson as the ultimate taxonomic authority).

This repository includes the code required to reproduce the outputs of Baiotto and Cosma et al. (2026), organized so that the framework can be adapted to other source data and applications with moderate effort (see Organization below).

# Organization

NECTAR's code is organized into 12 steps that are intended to be completed in order. The results presented in Baiotto and Cosma generally rely on the code in steps 01-07 and 12 and skip steps 08-11, which are various data summary/aggregation steps to support the implementation of NECTAR on the California Native Plant Societies' Calscape Pollinator Companion tool (beta version available at: https://pollinators.calscape.org/pollinator-companion/welcome).

NECTAR is designed to be executed under a general project directory with the code placed in a `Code` subfolder. The code reads input data from `Data_Raw` and places cleaned datasets in `Data_Clean` and `Data_Calscape`, intermediate files in `Temp` and `logs`, and other outputs in `Figures` and `Tables`.

📁 01_Data_Download - Automated data download (e.g., for plant occurrence and supporting environmental/ecological data from Calflora) and descriptions for manual data download and very basic cleaning/formatting.

📁 02_Occurrence_Cleaning - Harmonization of occurrence records from all sources and spatial, temporal, taxonomic, and column cleaning generally following established occurrence cleaning workflows (e.g., biodiversity data cleaning workflow - `bdc`)

📁 03_Checklist_Generation - Generation of taxon-specific and regionally-explicit checklists of plant and pollinator species using external sources/authorities (e.g., using Jepson to constrain to only California native plants) or occurrence data (when there is not a good external authority, e.g., hoverflies).

📁 04_SDMs - Ensemble species distribution models for all species in our checklists using a three-tiered approach based on data availability.

📁 05_Phenology - Estimation of region-level phenometrics for key life stage events using a three-tiered approach based on data availability.

📁 06_Interaction_Data - Downloading, scraping, and aggregation/cleaning of plant-pollinator interaction data from multiple sources.

📁 07_Interaction_Predictions - Generation of large-scale metawebs using simple phylogenetic constraints and graph embedding with transfer learning, and downscaling of large-scale metawebs to regional metawebs using spatial and phenological co-occurrence for potentially interacting pairs identified in the large-scale metaweb.

📁 08_Alpha_Diversity - Summarization of pollinator support for all native plants included, broken up by ecoregion (this is the foundation of the ranked draw selection of pollinator-supporting native plants)

📁 09_Beta_Diversity - Identification of complimentary pollinator diversity companion plants (plants that support the most additional pollinator species) for each plant in each ecoregion using the predicted intercation networks and constraining companion plants to only those found in the same ecosystem(s) as the focal plant.

📁 10_Phenological_Diversity - Identification of complimentary phenological companion plants (plants that add the most additional contiguous flowering time) for each plant in each ecoregion using the estimated phenometrics and constraining companion plants to only those found in the same ecosystem(s) as the focal plant.

📁 11_Crop_Support - Interaction predictions for crops grown in California to aid in the identification of native plants that would best serve as companion plants for each crop (in each region), based on the diversity of polliators that we expect to visit the crops that also would be likely to be supported by each native plant.

📁 12_Analyses - Analyses conducted specifically for Baiotto and Cosma et al. (2026), including comparison of pollinator support native plants identified with NECTAR to random selection of plants and external lists (e.g., Xerces Society regional plant lists and Pollinator Partnership ecoregional planting guides), ecosystem restoration scenarios, and network/data summarisation. 

📁 99_Supporting - Supporting script files contain functions re-used many times (or for methods that could be easily used in other applications - e.g., phenometric estimated), including the single objective genetic algorithm we use for identifying complimentary sets of native plants in the optimized planting scenarios, taxonomic name harmonization, and phenometrics.

## Input data requirements

All raw inputs are in `Data_Raw/` in the project folder.

**Core inputs** (required to build regional metawebs):

| Input | Role in NECTAR | California implementation | Location (`Data_Raw/`) |
|---|---|---|---|
| Occurrence records | Species distribution models, phenology estimates, and regional species pools | Plants: Calflora, CCH2. Pollinators: GBIF, iDigBio, Ecdysis (formerly SCAN), Chesshire et al. bee records, Moth Photographers Group. Cleaned with bdc, CoordinateCleaner, and BeeBDC | `Occurrences/` (raw), `Formatted/` (formatted) |
| Flowering / life-stage annotations | Restrict plant phenology to flowering records and pollinators to adult life stages | CCH2 phenology records and GBIF records annotated for flowering (Phenobase); larval/juvenile records excluded for pollinators | `Occurrences/cch2_pheno/`, `Phenology/annotations_all2.csv` |
| Species pool and status lists | Define the native species pool and remove invasive species | Jepson eFlora and CNPS native plant lists; CNDDB Special Plants and Animals lists; Cal-IPC, CISAC, and US-RIIS invasive lists; MPG and BAMONA Lepidoptera checklists | `Species_Traits_Attributes/` |
| Interaction records | Observed interactions used to predict the metaweb and as a validation baseline | GloBI, Discover Life, Big Book of Hymenoptera (Apoidea), Fowler specialist bee host lists, NHM HOSTS database, Caldwell (2021) California plants as Lepidoptera resources, UNR/UNM butterfly host plant database, Shropshire & Tallamy (2025), CropPol | `Interactions/<source>/` |
| Phylogenies | Phylogenetic constraint and graph-embedding transfer to species without interaction data | Open Tree of Life subtrees for each pollinator group | `Phylogeny/*_names_tree_harm.csv` |
| Environmental predictors | Covariates for species distribution models | 5 arcmin stack over California + 100 km buffer: WorldClim 2.1 bioclimatic variables, SRTM elevation and derived slope/aspect/TRI, ESA WorldCover land cover, SoilGrids soil properties, human footprint, population density, and distance to coast and streams (Natural Earth) | `SDM_input/env_stack_final.tif` |
| Regional boundaries | Spatial units for downscaling the metaweb and SDM extent | Jepson ecoregions (province, region, and subregion levels); combined Jepson/EPA ecoregions and GADM state boundaries for SDMs | `Spatial/Jepson_Ecoregions/`, `SDM_input/` |
| Taxonomic authority | Final reference for name harmonization | Jepson eFlora for California plants; names parsed and verified with gnparser and gnverifier, plus a manual corrections table | `Species_Traits_Attributes/jepson_eflora_plants.csv`, `taxonomy_corrections.csv` |

**Application-specific inputs**:

| Input | Use | California implementation | Location (`Data_Raw/`) |
|---|---|---|---|
| External native plant lists and their regions | Benchmarking NECTAR plant recommendations | Pollinator Partnership ecoregion lists (mapped with Bailey ecoregions); Xerces Society regional lists with georeferenced regions | `External_NativePlantLists/`, `Spatial/Bailey_Ecoregions/` |
| Plant attribute data | Constraining candidate plants in restoration and companion-planting scenarios | Calflora plant communities (habitat), CNPS rarity rank, Calscape plant table (drought tolerance, nursery availability) | `Species_Traits_Attributes/` |
| Specialist and focal species lists | Scenario objectives | Fowler oligolectic bee host plants; declining Southwest butterflies from Edwards et al. (2025) | `Interactions/fowler/`, `Edwards2025/` |
| Crop list | Crop pollinator support analysis | Crop common and scientific names, linked to CropPol records | `Species_Traits_Attributes/cropNames.csv` |
| Independent interaction surveys | External validation of predicted interactions | Hung et al. (2019; South Coast), Hazlehurst et al. (2023; Southern Channel Islands), BBNA and BUGS datasets | `External_Validation/` |

# Software and dependencies

### R (v4.2.3)

Most analyses were run in R.

| Purpose | Packages |
|---|---|
| Data handling and manipulation | tidyverse (v2.0.0), data.table (v1.16.2) |
| Spatial data handling | sf (v1.0.21), terra (v1.8.70), raster (v3.6.32), tigris (v2.2.1), maps (v3.4.3), rmapshaper (v0.5.0), nngeo (v0.4.8) |
| Occurrence data download and access | rgbif (v3.8.3), ridigbio (v0.4.1), rvest (v1.0.4), httr (v1.4.7) |
| Occurrence data cleaning | CoordinateCleaner (v3.0.1), bdc (v1.1.5), BeeBDC (v1.3.1) |
| Species distribution modelling | flexsdm (v1.3.9), biomod2 (v4.2.6.2), ecospat (v4.1.2), blockCV (v3.2.0), CovSel (v0.1.0), dismo (v1.3.16), adehabitatHR (v0.4.22), mda (v0.5.5), earth (v5.3.4), geodata (v0.6.6), rnaturalearth (v1.1.0), nhdplusTools (v1.3.2), FedData (v4.3.0) |
| Phenology modelling | circular (v0.5.2), boot (v1.3.28.1) |
| Interaction data scraping and preparation | rvest (v1.0.4), rotl (v3.1.0) |
| Phylogenetic data handling | ape (v5.8.1) |
| Interaction prediction and network analysis | bipartite (v2.23.0), nngeo (v0.4.8) |
| Genetic algorithm plant selection | futureverse (v1.34.0) |
| Bayesian modelling and analysis | brms (v2.23.0), cmdstanr (v0.9.0), tidybayes (v3.0.7) |
| Visualization | patchwork (v1.3.0), cowplot (v1.1.3), ggdist (v3.3.3), ggalluvial (v0.12.5), tidyterra (v0.7.2), units (v0.8.5) |
| Reading external data files | readxl (v1.4.3) |

### Julia (v1.10.10)

Graph embedding and transfer learning were implemented in Julia.

| Package | Version |
|---|---|
| CSV.jl | v0.10.15 |
| DataFrames.jl | v1.7.0 |
| Distributions.jl | v0.25.120 |
| JLD2.jl | v0.6.2 |
| PhyloNetworks.jl | v0.16.4 |
| Pipe.jl | v1.3.0 |
| Random.jl | standard library |
| SpeciesInteractionNetworks.jl | v0.0.1 |
| StatsBase.jl | v0.34.6 |

### External command-line tools

Taxonomic name parsing and harmonization require two command-line tools from Global Names. Both must be installed and on your system `PATH`, or have their location passed to the harmonization functions (e.g. the `gnv_path` argument).

| Tool | Purpose | Source |
|---|---|---|
| gnparser | Parsing scientific names into their components | [github.com/gnames/gnparser](https://github.com/gnames/gnparser) |
| gnverifier | Verifying and harmonizing names against taxonomic databases | [github.com/gnames/gnverifier](https://github.com/gnames/gnverifier) |
