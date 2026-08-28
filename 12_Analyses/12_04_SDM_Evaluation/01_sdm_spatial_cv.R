#!/usr/bin/env Rscript

#### Spatial cross-validation for tiered SDMs ####

# USAGE
#
# LOCAL COMPUTER:
#   Rscript 01_sdm_spatial_cv.R
#       -> automatically runs a small end-to-end test subset locally
#
# CLUSTER (from terminal):
#   1. cd /local/storage/NECTAR/Code/12_Analyses/12_04_SDM_Evaluation
#   2. Rscript 01_sdm_spatial_cv.R prep
#   3. sbatch  sdm_spatial_cv-submit.sh
#   4. Rscript 01_sdm_spatial_cv.R compile
#   5. Rscript 02_sdm_eval_final_results.R
#
# The Slurm .sh passes only SLURM_ARRAY_TASK_ID

# ---- Simple settings you may actually change ---------------------------------

LOCAL_TEST_N_PER_TAXON <- 2L

# Total number of Slurm array slots. This MUST match --array=1-200 in the .sh.
# Preparation allocates these slots across taxa according to estimated workload,
# while keeping each individual task restricted to one taxon.
CLUSTER_ARRAY_TASKS <- 200L

# Local and cluster project roots. The script chooses whichever exists.
LOCAL_ROOT   <- "path/to/local/root"
CLUSTER_ROOT <- "/local/storage/NECTAR"

# ---- Determine how the script was called -------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)

SCRIPT_FILE <- {
  x <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(x)) {
    stop("Run this file with Rscript rather than source().")
  }
  normalizePath(sub("^--file=", "", x[1]), mustWork = TRUE)
}

if (!length(ARGS)) {
  RUN_MODE <- "local"
} else if (tolower(ARGS[1]) == "prep") {
  RUN_MODE <- "prepare"
} else if (tolower(ARGS[1]) == "prep-test") {
  RUN_MODE <- "prepare_test"
} else if (tolower(ARGS[1]) == "compile") {
  RUN_MODE <- "compile"
} else {
  TASK_ID <- suppressWarnings(as.integer(ARGS[1]))
  if (!is.finite(TASK_ID) || TASK_ID < 1L) {
    stop("Argument must be `prep`, `compile`, or a positive array task ID.")
  }
  RUN_MODE <- "worker"
}

IS_LOCAL_TEST <- RUN_MODE %in% c("local", "prepare_test") ||
  identical(Sys.getenv("SDM_LOCAL_TEST", unset = ""), "1")

# ---- Configuration -----------------------------------------------------------

SCRIPT_VERSION <- "2026-08-08-cluster-optimized-v3"

ROOT <- Sys.getenv("SDM_ROOT", unset = "")
if (!nzchar(ROOT)) {
  ROOT <- if (dir.exists(CLUSTER_ROOT)) CLUSTER_ROOT else LOCAL_ROOT
}
ROOT <- normalizePath(ROOT, mustWork = TRUE)

OCC_PATH <- file.path(
  ROOT,
  "Data_Raw", "SDM_input", "sdm_data.csv"
)

MANIFEST_DIR <- file.path(
  ROOT,
  "Data_Clean", "SDMs", "sdm_by_taxon"
)

OUTPUT_ROOT <- Sys.getenv("SDM_CV_OUTPUT_ROOT", unset = "")
if (!nzchar(OUTPUT_ROOT)) {
  OUTPUT_ROOT <- if (IS_LOCAL_TEST) {
    file.path(ROOT, "SDM_Eval", "local_test")
  } else {
    file.path(ROOT, "SDM_Eval")
  }
}
CACHE_DIR    <- file.path(OUTPUT_ROOT, "cache")
ARRAY_DIR    <- file.path(OUTPUT_ROOT, "array")
SPECIES_DIR  <- file.path(OUTPUT_ROOT, "species")
WORK_ROOT    <- file.path(OUTPUT_ROOT, "work")
TEMP_ROOT    <- file.path(OUTPUT_ROOT, "terra_tmp")
COMBINED_DIR <- file.path(OUTPUT_ROOT, "combined")

TARGET_TAXA <- c(
  "plants",
  "bees",
  "butterflies",
  "hoverflies",
  "moths",
  "hummingbirds"
)

TARGET_SPECIES <- NULL
BATCH_COUNT <- 1L
BATCH_ID <- 1L
COMPILE_ONLY <- FALSE

RESUME_COMPLETED <- TRUE
RETRY_FAILED     <- TRUE
FORCE_RERUN      <- FALSE

KEEP_SUCCESSFUL_MODEL_WORK <- FALSE
KEEP_FAILED_MODEL_WORK     <- FALSE
KEEP_SUCCESSFUL_CONSOLE_LOGS <- FALSE
KEEP_FAILED_CONSOLE_LOGS     <- TRUE

WRITE_COMBINED_FOLD_ASSIGNMENTS <- TRUE
WRITE_GLOBAL_FOLD_ASSIGNMENTS   <- FALSE

SEED <- 123L

N_BG_ESM <- 5000L
N_BG_SDM <- 10000L

MAX_FOLDS <- 5L
MIN_FOLDS <- 3L

MIN_TEST_PRES_ESM  <- 1L
MIN_TEST_PRES_SDM  <- 2L
MIN_TEST_BG        <- 5L
MIN_TRAIN_PRES_ESM <- 4L
MIN_TRAIN_PRES_SDM <- 10L
MIN_ESM_ALGORITHMS <- 2L

ALLOW_ENV_FILTER_FALLBACK <- TRUE

BLOCK_ITERATIONS          <- 200L
AUTOCOR_RANGE_MULTIPLIER  <- 1.25
MAX_BLOCK_EXTENT_FRACTION <- 0.25
MIN_BLOCK_CELLS           <- 5L
AUTOCOR_SAMPLE_SIZE       <- 3000L
AUTOCOR_DOMAIN_BUFFER_M   <- 125000

# Slurm supplies parallelism across array tasks. Each task still processes
# species sequentially, while BIOMOD2 may use two independent workers for the
# fold-by-algorithm fits within the current species. Local smoke tests remain
# single-core to keep them lightweight and maximally portable.
PHYSICAL_CORES <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(PHYSICAL_CORES) || PHYSICAL_CORES < 1L) PHYSICAL_CORES <- 1L

SPECIES_WORKERS <- 1L
BIOMOD_CORES <- if (IS_LOCAL_TEST) 1L else 2L
PROGRESS_CHUNK_SIZE <- 1L

if (RUN_MODE == "worker") {
  slurm_cpus_raw <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (nzchar(slurm_cpus_raw)) {
    slurm_cpus <- suppressWarnings(as.integer(slurm_cpus_raw))
    if (!is.finite(slurm_cpus) || slurm_cpus < BIOMOD_CORES) {
      stop(
        "This worker requires at least ", BIOMOD_CORES,
        " CPUs per Slurm task, but SLURM_CPUS_PER_TASK=", slurm_cpus_raw,
        ". Update --cpus-per-task in the submit script."
      )
    }
  }
}

TERRA_MEMFRAC <- suppressWarnings(as.numeric(
  Sys.getenv("SDM_TERRA_MEMFRAC", unset = "0.50")
))
if (!is.finite(TERRA_MEMFRAC) || TERRA_MEMFRAC <= 0 || TERRA_MEMFRAC > 0.9) {
  TERRA_MEMFRAC <- 0.50
}
TERRA_THREADS <- 1L

ENV_CANDIDATES <- c(
  file.path(
    ROOT,
    "Data_Raw", "SDM_input", "env_stack_final.tif"
  )
)

ECO_CANDIDATES <- c(
  file.path(
    ROOT,
    "Data_Raw", "SDM_input",
    "jepson_epa_ecoregions_combined.gpkg"
  )
)

CA_ALBERS <- "EPSG:3310"

options(device = function(...) grDevices::pdf(NULL))

# ---- Thread limits before model packages are loaded --------------------------

# Prevent nested BLAS/OpenMP thread storms inside BIOMOD2 parallel workers.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "data.table", "terra", "sf", "blockCV", "biomod2", "ecospat",
  "flexsdm", "pROC", "foreach", "mgcv", "randomForest", "maxnet",
  "R.utils"
)
if (BIOMOD_CORES > 1L) {
  required_packages <- c(required_packages, "doParallel")
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following packages before running:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
  library(sf)
  library(blockCV)
  library(biomod2)
  library(ecospat)
})

foreach::registerDoSEQ()

package_versions <- vapply(
  required_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
)

# ---- Directories and terra options -------------------------------------------

for (d in c(
  OUTPUT_ROOT, CACHE_DIR, ARRAY_DIR, SPECIES_DIR, WORK_ROOT, TEMP_ROOT,
  COMBINED_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

terra::terraOptions(
  memfrac = TERRA_MEMFRAC,
  tempdir = TEMP_ROOT,
  progress = 0,
  parallel = FALSE,
  threads = TERRA_THREADS
)

# ---- General helpers ---------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

sanitize <- function(x) {
  out <- gsub("[^[:alnum:]_-]", "_", enc2utf8(x))
  out <- gsub("_+", "_", out)
  out <- sub("^_+", "", out)
  out <- sub("_+$", "", out)
  ifelse(nzchar(out), out, "unnamed")
}

first_existing <- function(paths, label) {
  out <- paths[file.exists(paths)]
  if (!length(out)) {
    stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  }
  normalizePath(out[1], mustWork = TRUE)
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(is.finite(x)) > 1L) stats::sd(x, na.rm = TRUE) else NA_real_
}

normalize_biomod_score <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite_x <- x[is.finite(x)]
  if (length(finite_x) && max(abs(finite_x), na.rm = TRUE) > 1.5) {
    x <- x / 1000
  }
  x
}

standardize_metric <- function(x) {
  x <- toupper(as.character(x))
  x[x %in% c("ROC", "AUCROC", "AUC_ROC")] <- "AUC"
  x[x %in% c("MAXTSS", "MAX_TSS")] <- "TSS"
  x
}

file_signature <- function(path) {
  info <- file.info(path)
  paste(
    normalizePath(path, mustWork = TRUE),
    info$size,
    format(info$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
    sep = "|"
  )
}

md5_text <- function(x) {
  tf <- tempfile(pattern = "signature_", fileext = ".txt")
  on.exit(unlink(tf, force = TRUE), add = TRUE)
  writeLines(enc2utf8(x), tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(
    dirname(path),
    paste0(".", basename(path), ".", Sys.getpid(), ".tmp")
  )
  saveRDS(object, tmp, compress = compress)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp, force = TRUE)
    stop("Could not atomically write: ", path)
  }
  invisible(path)
}

atomic_fwrite <- function(x, path, compress = "none") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suffix <- if (identical(compress, "gzip")) ".gz" else ""
  tmp <- file.path(
    dirname(path),
    paste0(".", sub("\\.gz$", "", basename(path)), ".", Sys.getpid(), ".tmp", suffix)
  )
  data.table::fwrite(x, tmp, compress = compress)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp, force = TRUE)
    stop("Could not atomically write: ", path)
  }
  invisible(path)
}

clean_warning_text <- function(x) {
  x <- unique(as.character(x))
  x[
    !is.na(x) &
      nzchar(trimws(x)) &
      !trimws(x) %chin% c("NaN", "NULL")
  ]
}

stable_hash_int <- function(x) {
  codes <- utf8ToInt(enc2utf8(x))
  hash <- 0
  
  for (code in codes) {
    hash <- (hash * 131 + code) %% 2147483647
  }
  
  as.integer(hash)
}

species_seed <- function(taxon, species) {
  key <- paste(taxon, species, sep = "||")
  as.integer((SEED + stable_hash_int(key)) %% 2147483647)
}

species_folder_name <- function(species_key, species) {
  paste0(species_key, "_", sanitize(species))
}

species_result_dir <- function(row) {
  file.path(
    SPECIES_DIR,
    row$taxon,
    species_folder_name(row$species_key, row$species)
  )
}

species_work_dir <- function(row) {
  file.path(
    WORK_ROOT,
    row$taxon,
    species_folder_name(row$species_key, row$species)
  )
}

species_temp_dir <- function(row) {
  file.path(
    TEMP_ROOT,
    paste0(
      row$taxon, "_",
      species_folder_name(row$species_key, row$species),
      "_", Sys.getpid()
    )
  )
}

status_path_for <- function(row) {
  file.path(species_result_dir(row), "status.rds")
}

read_status <- function(row) {
  path <- status_path_for(row)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

species_checkpoint_files <- function(row) {
  out_dir <- species_result_dir(row)
  file.path(
    out_dir,
    c(
      "species_metrics.csv",
      "component_metrics.csv",
      "fold_summary.csv",
      "fold_assignments.csv.gz",
      "errors.csv",
      "warnings.csv",
      "predictor_mapping.csv",
      "status.rds"
    )
  )
}

status_is_current_complete <- function(row, run_signature) {
  status <- read_status(row)
  !is.null(status) &&
    identical(status$state, "complete") &&
    identical(status$run_signature, run_signature) &&
    all(file.exists(species_checkpoint_files(row)))
}

status_is_current_failed <- function(row, run_signature) {
  status <- read_status(row)
  !is.null(status) &&
    identical(status$state, "failed") &&
    identical(status$run_signature, run_signature) &&
    file.exists(file.path(species_result_dir(row), "errors.csv"))
}

# ---- Output schemas -----------------------------------------------------------

species_metrics_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  metric = character(),
  mean = numeric(),
  sd = numeric(),
  n_folds_valid = integer(),
  n_folds_expected = integer(),
  calibration_mean = numeric(),
  optimism_gap = numeric(),
  evaluation_model = character(),
  n_presences = integer(),
  n_background = integer(),
  n_spatial_folds = integer(),
  block_size_km = numeric(),
  autocor_range_km = numeric(),
  predictor_count = integer(),
  background_source = character(),
  background_candidates = integer(),
  n_algorithms_success = integer()
)

component_metrics_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  component_type = character(),
  component = character(),
  metric = character(),
  value = numeric(),
  n_successful_bivariate_models = integer(),
  n_bivariate_min = numeric(),
  n_bivariate_median = numeric(),
  n_bivariate_max = numeric(),
  fold_bivariate_counts = character()
)

fold_assignments_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  longitude = numeric(),
  latitude = numeric(),
  pr_ab = integer(),
  fold = integer()
)

fold_summary_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  fold = integer(),
  train_presences = integer(),
  train_background = integer(),
  test_presences = integer(),
  test_background = integer(),
  block_size_km = numeric(),
  autocor_range_km = numeric()
)

errors_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  stage = character(),
  error = character()
)

warnings_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  stage = character(),
  warning = character()
)

predictor_mapping_prototype <- data.table(
  species = character(),
  taxon = character(),
  model_type = character(),
  requested_predictor = character(),
  raster_layer = character(),
  match_type = character()
)

status_prototype <- data.table(
  species_id = integer(),
  species_key = character(),
  species = character(),
  taxon = character(),
  model_type = character(),
  state = character(),
  stage = character(),
  elapsed_seconds = numeric(),
  completed_at = character(),
  run_signature = character()
)

# ---- Manifest and predictor matching -----------------------------------------

read_manifests <- function(path) {
  files <- list.files(
    path,
    pattern = "_continuous_manifest\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(files)) stop("No continuous model manifests found under: ", path)
  
  tables <- lapply(files, function(z) {
    x <- tryCatch(
      fread(z, showProgress = FALSE),
      error = function(e) {
        stop("Could not read manifest ", z, ": ", conditionMessage(e))
      }
    )
    x[, manifest_file := normalizePath(z, mustWork = TRUE)]
    x
  })
  
  out <- rbindlist(tables, fill = TRUE, use.names = TRUE)
  
  needed <- c("species", "taxon", "model_type", "selected_predictors")
  missing <- setdiff(needed, names(out))
  if (length(missing)) {
    stop("Manifest is missing columns: ", paste(missing, collapse = ", "))
  }
  
  duplicate_check <- out[
    ,
    .(
      n_rows = .N,
      n_model_types = uniqueN(model_type),
      n_predictor_strings = uniqueN(selected_predictors)
    ),
    by = .(species, taxon)
  ][n_rows > 1L]
  
  conflicting <- duplicate_check[
    n_model_types > 1L | n_predictor_strings > 1L
  ]
  if (nrow(conflicting)) {
    stop(
      "Conflicting duplicate manifest rows were found for: ",
      paste(
        paste(conflicting$taxon, conflicting$species, sep = "::"),
        collapse = ", "
      )
    )
  }
  
  unique(out, by = c("species", "taxon"))
}

clean_predictor_names <- function(x) {
  make.unique(gsub("[^A-Za-z0-9_]", "_", x), sep = "_")
}

resolve_predictors <- function(x, available) {
  requested <- trimws(unlist(strsplit(x %||% "", ",", fixed = TRUE)))
  requested <- requested[nzchar(requested)]
  
  if (length(requested) < 2L) {
    stop("Manifest contains fewer than two selected predictors.")
  }
  if (anyDuplicated(requested)) {
    stop(
      "Manifest contains duplicated selected predictor names: ",
      paste(unique(requested[duplicated(requested)]), collapse = ", ")
    )
  }
  
  aliases <- clean_predictor_names(available)
  
  mapping <- rbindlist(lapply(requested, function(req) {
    exact <- which(available == req)
    
    if (length(exact) == 1L) {
      return(data.table(
        requested_predictor = req,
        raster_layer = available[exact],
        match_type = "exact"
      ))
    }
    if (length(exact) > 1L) {
      stop("Environmental stack contains duplicated layer name: ", req)
    }
    
    alias <- which(aliases == req)
    if (length(alias) == 1L) {
      return(data.table(
        requested_predictor = req,
        raster_layer = available[alias],
        match_type = "cleaned_alias"
      ))
    }
    if (length(alias) > 1L) {
      stop(
        "Predictor alias is ambiguous for ", req, ": ",
        paste(available[alias], collapse = ", ")
      )
    }
    
    stop("Manifest predictor did not match any environmental layer: ", req)
  }))
  
  if (anyDuplicated(mapping$raster_layer)) {
    stop(
      "Multiple manifest predictors resolved to the same raster layer: ",
      paste(mapping$raster_layer[duplicated(mapping$raster_layer)], collapse = ", ")
    )
  }
  
  list(
    raster_names = mapping$raster_layer,
    model_names = mapping$requested_predictor,
    mapping = mapping
  )
}

# ---- Occurrence access and cached target-group pools -------------------------

get_target_occurrences <- function(occ, target_name) {
  a <- occ[
    .(target_name),
    on = .(species),
    nomatch = 0L,
    .(species = target_name, longitude, latitude)
  ]
  
  b <- occ[
    .(target_name),
    on = .(species_full),
    nomatch = 0L,
    .(species = target_name, longitude, latitude)
  ]
  
  unique(rbindlist(list(a, b), use.names = TRUE))
}

one_per_cell <- function(x, template) {
  if (!nrow(x)) return(x)
  
  cells <- terra::cellFromXY(
    template,
    as.matrix(x[, .(longitude, latitude)])
  )
  
  keep <- !is.na(cells) & !duplicated(cells)
  x[keep]
}

background_pool_cache_paths <- function(taxon) {
  stem <- paste0("target_group_occurrences_", sanitize(taxon))
  list(
    data = file.path(CACHE_DIR, paste0(stem, ".rds")),
    metadata = file.path(CACHE_DIR, paste0(stem, "_metadata.rds"))
  )
}

background_pool_is_current <- function(taxon, pool_signature) {
  paths <- background_pool_cache_paths(taxon)
  
  if (!file.exists(paths$data) || !file.exists(paths$metadata)) {
    return(FALSE)
  }
  
  metadata <- tryCatch(
    readRDS(paths$metadata),
    error = function(e) NULL
  )
  
  !is.null(metadata) &&
    identical(metadata$signature, pool_signature)
}

prepare_background_pool <- function(
    taxon, occ, env_template, pool_signature, force = FALSE) {
  
  paths <- background_pool_cache_paths(taxon)
  
  if (!force && file.exists(paths$data) && file.exists(paths$metadata)) {
    metadata <- tryCatch(
      readRDS(paths$metadata),
      error = function(e) NULL
    )
    
    if (!is.null(metadata) &&
        identical(metadata$signature, pool_signature)) {
      return(paths$data)
    }
  }
  
  message("Building labeled target-group occurrence cache for ", taxon, "...")
  
  if (is.null(occ)) {
    stop("Occurrence data were not loaded for a stale target-group cache.")
  }
  
  taxon_value <- taxon
  x <- occ[
    taxon == taxon_value,
    .(species, species_full, longitude, latitude)
  ]
  
  if (!nrow(x)) {
    stop("No occurrence records were available for taxon: ", taxon)
  }
  
  # Exact duplicate records can be removed once without changing target-group
  # sampling. Species identities are retained so that all records belonging to
  # the focal species can be excluded before its background is sampled.
  x <- unique(x)
  
  cells <- terra::cellFromXY(
    env_template,
    as.matrix(x[, .(longitude, latitude)])
  )
  
  keep <- !is.na(cells)
  x <- x[keep]
  x[, cell := as.numeric(cells[keep])]
  
  atomic_save_rds(x, paths$data, compress = FALSE)
  atomic_save_rds(
    list(
      signature = pool_signature,
      taxon = taxon,
      n_records = nrow(x),
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    paths$metadata,
    compress = FALSE
  )
  
  paths$data
}

load_background_pool <- function(path) {
  x <- readRDS(path)
  setDT(x)
  
  required <- c(
    "species", "species_full", "longitude", "latitude", "cell"
  )
  missing <- setdiff(required, names(x))
  
  if (length(missing)) {
    stop(
      "Target-group cache is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  setindexv(x, "species")
  setindexv(x, "species_full")
  setindexv(x, "cell")
  x
}

# ---- Ecoregions and accessible areas -----------------------------------------

prepare_ecoregions <- function(path) {
  ecoregions <- st_read(path, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(CA_ALBERS)
  
  if (!"ecoregion_code" %in% names(ecoregions)) {
    stop("Ecoregion layer must contain an `ecoregion_code` column.")
  }
  
  ecoregions <- ecoregions[
    !is.na(ecoregions$ecoregion_code) &
      nzchar(as.character(ecoregions$ecoregion_code)),
  ]
  
  groups <- split(
    seq_len(nrow(ecoregions)),
    as.character(ecoregions$ecoregion_code)
  )
  
  unioned <- lapply(names(groups), function(code) {
    geom <- st_union(st_geometry(ecoregions[groups[[code]], ]))
    st_sf(
      ecoregion_code = code,
      geometry = geom,
      crs = st_crs(ecoregions)
    )
  })
  
  do.call(rbind, unioned) |>
    st_make_valid()
}

prepare_ecoregions_cached <- function(path, force = FALSE) {
  data_path <- file.path(CACHE_DIR, "ecoregions_by_code.rds")
  metadata_path <- file.path(
    CACHE_DIR,
    "ecoregions_by_code_metadata.rds"
  )
  signature <- md5_text(c(
    "ecoregions_by_code_v1",
    file_signature(path),
    CA_ALBERS,
    paste0("sf=", package_versions[["sf"]])
  ))
  
  if (!force && file.exists(data_path) && file.exists(metadata_path)) {
    metadata <- tryCatch(
      readRDS(metadata_path),
      error = function(e) NULL
    )
    
    if (!is.null(metadata) &&
        identical(metadata$signature, signature)) {
      cached <- tryCatch(readRDS(data_path), error = function(e) NULL)
      if (!is.null(cached) && inherits(cached, "sf")) {
        return(cached)
      }
    }
  }
  
  output <- prepare_ecoregions(path)
  atomic_save_rds(output, data_path, compress = FALSE)
  atomic_save_rds(
    list(
      signature = signature,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    metadata_path,
    compress = FALSE
  )
  output
}

make_accessible_area <- function(sp_clean, ecoregions_by_code) {
  points_albers <- st_as_sf(
    sp_clean,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    st_transform(CA_ALBERS)
  
  intersections <- st_intersects(points_albers, ecoregions_by_code)
  eco_index <- unique(unlist(intersections, use.names = FALSE))
  
  if (!length(eco_index)) {
    stop("No occurrence points intersected the ecoregion layer.")
  }
  
  codes <- unique(as.character(
    ecoregions_by_code$ecoregion_code[eco_index]
  ))
  codes <- codes[!is.na(codes) & nzchar(codes)]
  
  if (!length(codes)) {
    stop("Occurrence intersections did not return valid ecoregion codes.")
  }
  
  eco_geom <- st_union(st_geometry(
    ecoregions_by_code[
      as.character(ecoregions_by_code$ecoregion_code) %in% codes,
    ]
  ))
  
  n <- nrow(sp_clean)
  buffer_m <- if (n < 7L) {
    25000
  } else if (n < 30L) {
    50000
  } else {
    100000
  }
  
  accessible_albers <- st_buffer(eco_geom, dist = buffer_m) |>
    st_make_valid()
  
  accessible_albers_sf <- st_sf(
    geometry = accessible_albers,
    crs = st_crs(ecoregions_by_code)
  )
  
  accessible_wgs84_sf <- st_transform(accessible_albers_sf, 4326)
  
  list(
    vector = terra::vect(accessible_wgs84_sf),
    sf_albers = accessible_albers_sf,
    sf_wgs84 = accessible_wgs84_sf,
    ecoregion_codes = codes,
    buffer_m = buffer_m
  )
}

# ---- Cached projected environmental stack for fold design -------------------

autocor_cache_paths <- function() {
  list(
    raster = file.path(CACHE_DIR, "env_autocor_albers.tif"),
    metadata = file.path(CACHE_DIR, "env_autocor_albers_metadata.rds")
  )
}

prepare_autocor_cache <- function(
    env_stack, ecoregions_by_code, cache_signature, force = FALSE) {
  
  paths <- autocor_cache_paths()
  
  if (!force && file.exists(paths$raster) && file.exists(paths$metadata)) {
    meta <- tryCatch(readRDS(paths$metadata), error = function(e) NULL)
    cached <- tryCatch(terra::rast(paths$raster), error = function(e) NULL)
    
    if (!is.null(meta) &&
        !is.null(cached) &&
        identical(meta$signature, cache_signature) &&
        identical(meta$layer_names, names(cached))) {
      return(paths$raster)
    }
  }
  
  message("Creating one-time Albers environmental cache for spatial folds...")
  
  factor_layers <- tryCatch(
    as.logical(terra::is.factor(env_stack)),
    error = function(e) rep(FALSE, terra::nlyr(env_stack))
  )
  if (length(factor_layers) != terra::nlyr(env_stack)) {
    factor_layers <- rep(FALSE, terra::nlyr(env_stack))
  }
  
  continuous_names <- names(env_stack)[!factor_layers]
  if (!length(continuous_names)) {
    continuous_names <- names(env_stack)[1]
    projection_method <- "near"
  } else {
    projection_method <- "bilinear"
  }
  
  study_domain_albers <- st_union(st_geometry(ecoregions_by_code)) |>
    st_buffer(AUTOCOR_DOMAIN_BUFFER_M) |>
    st_make_valid()
  
  study_domain_wgs84 <- st_transform(
    st_sf(
      geometry = study_domain_albers,
      crs = st_crs(ecoregions_by_code)
    ),
    4326
  )
  
  domain_vect <- terra::vect(study_domain_wgs84)
  
  env_crop <- terra::crop(
    env_stack[[continuous_names]],
    domain_vect,
    snap = "out"
  )
  env_crop <- terra::mask(env_crop, domain_vect)
  
  tmp_raster <- file.path(
    CACHE_DIR,
    paste0(".env_autocor_albers_", Sys.getpid(), ".tif")
  )
  
  old_options <- terra::terraOptions(print = FALSE)
  on.exit(
    terra::terraOptions(
      memfrac = old_options$memfrac,
      tempdir = old_options$tempdir,
      progress = old_options$progress,
      parallel = old_options$parallel,
      threads = old_options$threads
    ),
    add = TRUE
  )
  
  cache_threads <- 1L
  terra::terraOptions(
    memfrac = max(TERRA_MEMFRAC, 0.5),
    progress = 0,
    parallel = FALSE,
    threads = 1L
  )
  
  projected <- terra::project(
    env_crop,
    CA_ALBERS,
    method = projection_method,
    filename = tmp_raster,
    overwrite = TRUE,
    threads = FALSE
  )
  
  # GeoTIFF band names are normally preserved. If the driver changes them,
  # rewrite once to a second temporary file with the intended layer names.
  if (!identical(names(projected), continuous_names)) {
    tmp_named <- file.path(
      CACHE_DIR,
      paste0(".env_autocor_albers_named_", Sys.getpid(), ".tif")
    )
    names(projected) <- continuous_names
    terra::writeRaster(
      projected,
      tmp_named,
      overwrite = TRUE
    )
    unlink(tmp_raster, force = TRUE)
    tmp_raster <- tmp_named
  }
  
  rm(projected)
  gc(verbose = FALSE)
  
  if (file.exists(paths$raster)) unlink(paths$raster, force = TRUE)
  if (!file.rename(tmp_raster, paths$raster)) {
    unlink(tmp_raster, force = TRUE)
    stop("Could not finalize projected environmental cache.")
  }
  
  atomic_save_rds(
    list(
      signature = cache_signature,
      layer_names = continuous_names,
      crs = CA_ALBERS,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    paths$metadata,
    compress = FALSE
  )
  
  paths$raster
}

# ---- Occurrence preprocessing ------------------------------------------------

prepare_species_occurrences <- function(
    target_name, occ, env_stack, expected_model_type) {
  
  sp <- get_target_occurrences(occ, target_name)
  if (!nrow(sp)) stop("No occurrence records matched this target.")
  
  n_original <- nrow(sp)
  sp <- one_per_cell(sp, env_stack[[1]])
  n_thinned <- nrow(sp)
  
  n_raw <- if (n_thinned >= 6L) n_original else n_thinned
  rare_or_ultra <- n_raw <= 24L
  
  if (!rare_or_ultra) {
    sp_env <- terra::extract(
      env_stack,
      as.matrix(sp[, .(longitude, latitude)])
    )
    
    if (sum(complete.cases(sp_env)) >= 3L) {
      sp <- tryCatch(
        as.data.table(flexsdm::occfilt_env(
          data = as.data.frame(sp),
          x = "longitude",
          y = "latitude",
          id = "species",
          env_layer = env_stack,
          nbins = 5
        )),
        error = function(e) {
          if (!ALLOW_ENV_FILTER_FALLBACK) {
            stop(
              target_name,
              ": environmental occurrence filtering failed. ",
              conditionMessage(e)
            )
          }
          
          warning(
            target_name,
            ": environmental occurrence filtering failed; using thinned records. ",
            conditionMessage(e)
          )
          sp
        }
      )
    }
  }
  
  sp <- sp[
    is.finite(longitude) & is.finite(latitude),
    .(species = target_name, longitude, latitude)
  ]
  sp <- unique(sp)
  
  if (!nrow(sp)) stop("No occurrences remained after preprocessing.")
  
  inferred <- if (nrow(sp) < 7L) {
    "ultraRare_alphahull"
  } else if (nrow(sp) <= 24L) {
    "ESM"
  } else {
    "SDM"
  }
  
  if (!identical(inferred, expected_model_type)) {
    warning(sprintf(
      "%s: reconstructed tier is %s but manifest says %s; using manifest tier.",
      target_name, inferred, expected_model_type
    ))
  }
  
  sp
}

# ---- Fast target-group background generation --------------------------------

generate_target_group_background <- function(
    env_selected, accessible, background_pool, target_name, presence_xy,
    n_target, local_seed) {
  
  bbox <- st_bbox(accessible$sf_wgs84)
  
  candidates <- background_pool[
    species != target_name &
      (is.na(species_full) | species_full != target_name) &
      longitude >= bbox["xmin"] &
      longitude <= bbox["xmax"] &
      latitude >= bbox["ymin"] &
      latitude <= bbox["ymax"]
  ]
  
  source <- "target_group"
  candidate_env <- NULL
  
  if (nrow(candidates)) {
    candidate_vect <- terra::vect(
      candidates[, .(longitude, latitude)],
      geom = c("longitude", "latitude"),
      crs = "EPSG:4326",
      keepgeom = FALSE
    )
    
    inside <- terra::is.related(
      candidate_vect,
      accessible$vector,
      relation = "within"
    )
    inside <- !is.na(inside) & as.logical(inside)
    
    candidates <- candidates[inside]
  }
  
  if (nrow(candidates)) {
    # Environmental values are constant within a raster cell. Keep the extracted
    # values aligned with candidate rows so the sampled backgrounds do not need a
    # second raster extraction later in process_species().
    candidates <- candidates[!duplicated(cell)]
    candidate_env <- as.data.frame(
      terra::extract(env_selected, candidates$cell)
    )
    keep <- complete.cases(candidate_env)
    candidates <- candidates[keep]
    candidate_env <- candidate_env[keep, , drop = FALSE]
  }
  
  presence_cells <- terra::cellFromXY(
    env_selected[[1]],
    as.matrix(presence_xy[, .(longitude, latitude)])
  )
  presence_cells <- unique(presence_cells[!is.na(presence_cells)])
  
  if (nrow(candidates)) {
    keep <- !(candidates$cell %in% presence_cells)
    candidates <- candidates[keep]
    candidate_env <- candidate_env[keep, , drop = FALSE]
  }
  
  n_candidates <- nrow(candidates)
  
  if (nrow(candidates)) {
    set.seed(local_seed)
    sampled <- sample.int(
      nrow(candidates),
      min(n_target, nrow(candidates)),
      replace = FALSE
    )
    bg <- candidates[sampled, .(longitude, latitude)]
    bg_env <- candidate_env[sampled, , drop = FALSE]
  } else {
    source <- "random_fallback"
    
    masked_template <- terra::mask(
      terra::crop(
        env_selected[[1]],
        accessible$vector,
        snap = "out"
      ),
      accessible$vector
    )
    
    n_valid <- as.integer(terra::global(
      !is.na(masked_template),
      "sum",
      na.rm = TRUE
    )[1, 1])
    
    if (!is.finite(n_valid) || n_valid < 1L) {
      stop("Accessible environmental raster contains no valid cells.")
    }
    
    set.seed(local_seed)
    bg <- as.data.table(terra::spatSample(
      masked_template,
      size = min(n_target * 2L, n_valid),
      method = "random",
      na.rm = TRUE,
      xy = TRUE,
      values = FALSE,
      warn = FALSE
    ))
    setnames(bg, c("x", "y"), c("longitude", "latitude"))
    
    bg_cells <- terra::cellFromXY(
      env_selected[[1]],
      as.matrix(bg[, .(longitude, latitude)])
    )
    
    keep <- !is.na(bg_cells) &
      !(bg_cells %in% presence_cells) &
      !duplicated(bg_cells)
    bg <- bg[keep]
    bg_cells <- bg_cells[keep]
    
    bg_env <- as.data.frame(
      terra::extract(env_selected, as.numeric(bg_cells))
    )
    keep <- complete.cases(bg_env)
    bg <- bg[keep]
    bg_env <- bg_env[keep, , drop = FALSE]
    
    if (nrow(bg) > n_target) {
      keep <- seq_len(n_target)
      bg <- bg[keep]
      bg_env <- bg_env[keep, , drop = FALSE]
    }
    n_candidates <- nrow(bg)
  }
  
  if (nrow(bg) < MIN_TEST_BG * MIN_FOLDS) {
    stop(
      "Too few valid background points for spatial cross-validation: ",
      nrow(bg)
    )
  }
  
  if (nrow(bg_env) != nrow(bg) || ncol(bg_env) != terra::nlyr(env_selected)) {
    stop("Background environmental values were not aligned with sampled points.")
  }
  
  list(
    data = bg,
    env = bg_env,
    source = source,
    n_candidates = n_candidates
  )
}

# ---- Spatial folds -----------------------------------------------------------

make_spatial_folds <- function(
    model_data, env_autocor, block_template_cache,
    accessible, model_type, local_seed) {
  
  n_pres <- sum(model_data$pr_ab == 1L)
  k_start <- min(MAX_FOLDS, max(MIN_FOLDS, floor(n_pres / 3L)))
  k_values <- unique(seq(k_start, MIN_FOLDS, by = -1L))
  
  points_albers <- st_as_sf(
    model_data,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    st_transform(CA_ALBERS)
  
  accessible_vect_albers <- terra::vect(accessible$sf_albers)
  
  block_template <- terra::mask(
    terra::crop(
      block_template_cache,
      accessible_vect_albers,
      snap = "out"
    ),
    accessible_vect_albers
  )
  
  if (!is.null(env_autocor) && terra::nlyr(env_autocor) > 0L) {
    env_autocor <- terra::mask(
      terra::crop(
        env_autocor,
        accessible_vect_albers,
        snap = "out"
      ),
      accessible_vect_albers
    )
  } else {
    env_autocor <- NULL
  }
  
  if (terra::ncell(block_template) < 1L) {
    stop("Projected block template was empty for the accessible area.")
  }
  
  extent_albers <- terra::ext(block_template)
  width <- extent_albers$xmax - extent_albers$xmin
  height <- extent_albers$ymax - extent_albers$ymin
  cell_resolution <- mean(terra::res(block_template))
  
  autocor_args <- list(
    x = points_albers,
    column = "pr_ab",
    num_sample = min(AUTOCOR_SAMPLE_SIZE, terra::ncell(block_template)),
    progress = FALSE,
    plot = FALSE
  )
  if (!is.null(env_autocor)) autocor_args$r <- env_autocor
  
  available_autocor_args <- names(formals(blockCV::cv_spatial_autocor))
  autocor_args <- autocor_args[
    names(autocor_args) %in% available_autocor_args
  ]
  
  autocor <- tryCatch(
    do.call(blockCV::cv_spatial_autocor, autocor_args),
    error = function(e) {
      warning(
        "Spatial autocorrelation estimation failed; using the minimum ",
        "five-cell block scale. ", conditionMessage(e)
      )
      NULL
    }
  )
  
  ranges <- if (!is.null(autocor)) {
    suppressWarnings(as.numeric(unlist(autocor$range)))
  } else {
    numeric()
  }
  ranges <- ranges[is.finite(ranges) & ranges > 0]
  
  raw_range <- if (length(ranges)) {
    max(ranges)
  } else {
    MIN_BLOCK_CELLS * cell_resolution
  }
  
  minimum_block_size <- MIN_BLOCK_CELLS * cell_resolution
  desired_block_size <- max(
    minimum_block_size,
    AUTOCOR_RANGE_MULTIPLIER * raw_range
  )
  maximum_feasible_size <-
    MAX_BLOCK_EXTENT_FRACTION * max(width, height)
  
  block_size <- max(
    minimum_block_size,
    min(maximum_feasible_size, desired_block_size)
  )
  
  if (block_size + sqrt(.Machine$double.eps) < desired_block_size) {
    warning(sprintf(
      paste0(
        "Desired spatial block size (%.1f km) exceeded the configured ",
        "accessible-extent cap and was reduced to %.1f km."
      ),
      desired_block_size / 1000,
      block_size / 1000
    ))
  }
  
  min_test_pres <- if (model_type == "ESM") {
    MIN_TEST_PRES_ESM
  } else {
    MIN_TEST_PRES_SDM
  }
  
  min_train_pres <- if (model_type == "ESM") {
    MIN_TRAIN_PRES_ESM
  } else {
    MIN_TRAIN_PRES_SDM
  }
  
  for (k in k_values) {
    set.seed(local_seed + k)
    
    cv_args <- list(
      x = points_albers,
      column = "pr_ab",
      r = block_template,
      k = k,
      size = block_size,
      selection = "random",
      iteration = BLOCK_ITERATIONS,
      biomod2 = TRUE,
      presence_bg = TRUE,
      seed = local_seed + k,
      progress = FALSE,
      report = FALSE,
      plot = FALSE
    )
    
    available_cv_args <- names(formals(blockCV::cv_spatial))
    cv_args <- cv_args[names(cv_args) %in% available_cv_args]
    
    folds <- tryCatch(
      do.call(blockCV::cv_spatial, cv_args),
      error = function(e) NULL
    )
    
    if (is.null(folds) || is.null(folds$folds_ids)) next
    
    fold_id <- as.integer(folds$folds_ids)
    if (length(fold_id) != nrow(model_data) || anyNA(fold_id)) next
    
    fold_levels <- sort(unique(fold_id))
    
    balance <- rbindlist(lapply(fold_levels, function(fold) {
      test <- fold_id == fold
      data.table(
        fold = fold,
        train_presences = sum(model_data$pr_ab[!test] == 1L),
        train_background = sum(model_data$pr_ab[!test] == 0L),
        test_presences = sum(model_data$pr_ab[test] == 1L),
        test_background = sum(model_data$pr_ab[test] == 0L)
      )
    }))
    
    valid <- all(balance$test_presences >= min_test_pres) &&
      all(balance$test_background >= MIN_TEST_BG) &&
      all(balance$train_presences >= min_train_pres) &&
      all(balance$train_background >= MIN_TEST_BG)
    
    if (!valid) next
    
    cv_biomod <- vapply(
      fold_levels,
      function(fold) fold_id != fold,
      logical(nrow(model_data))
    )
    cv_biomod <- as.matrix(cv_biomod)
    colnames(cv_biomod) <- paste0(
      "_allData_RUN",
      seq_len(ncol(cv_biomod))
    )
    storage.mode(cv_biomod) <- "logical"
    
    if (any(rowSums(!cv_biomod) != 1L)) next
    
    return(list(
      cv_biomod = cv_biomod,
      fold_id = fold_id,
      k = ncol(cv_biomod),
      block_size_m = block_size,
      autocor_range_m = raw_range,
      balance = balance
    ))
  }
  
  stop(
    "Could not construct ", MIN_FOLDS, "-", MAX_FOLDS,
    " spatial folds meeting minimum train/test presence and background counts."
  )
}

# ---- Metrics -----------------------------------------------------------------

calculate_oof_metrics <- function(response, prediction) {
  keep <- is.finite(response) & is.finite(prediction)
  response <- as.numeric(response[keep])
  prediction <- as.numeric(prediction[keep])
  
  if (!length(response) || length(unique(response)) < 2L) {
    stop("Out-of-fold predictions do not contain both response classes.")
  }
  
  roc_object <- pROC::roc(
    response = response,
    predictor = prediction,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc <- as.numeric(pROC::auc(roc_object))
  
  best <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = c("sensitivity", "specificity"),
    transpose = FALSE
  )
  best <- as.data.frame(best)
  tss <- max(best$sensitivity + best$specificity - 1, na.rm = TRUE)
  
  presence_predictions <- prediction[response == 1]
  
  boyce <- if (length(presence_predictions) >= 3L) {
    tryCatch(
      as.numeric(ecospat::ecospat.boyce(
        fit = prediction,
        obs = presence_predictions,
        nclass = 0,
        window.w = "default",
        res = 100,
        PEplot = FALSE
      )$cor),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  data.table(
    metric = c("AUC", "TSS", "BOYCE", "SOMERSD"),
    value = c(auc, tss, boyce, 2 * auc - 1)
  )
}

# ---- ESM evaluation: direct BIOMOD2 equal-weight OOF implementation ----------
#
# ecospat 4.1.4 implements ecospat.ESM.Modeling() by looping over all bivariate
# predictor combinations and calling BIOMOD_Modeling() with:
#   * user-defined CV folds
#   * OPT.strategy = "bigboss"
#   * prevalence = 0.5
#   * var.import = 0
#   * AUCroc internally when weighting.score = "SomersD"
#
# ecospat also appends an all-data/allRun calibration column before those calls.
# That full-data model is not used anywhere in this strict OOF evaluator. The
# direct implementation below reproduces the bivariate BIOMOD2 calls on the
# exact same spatial folds but passes CV.do.full.models = FALSE and no allRun
# column, eliminating scientifically irrelevant full-data fits. BIOMOD_CORES
# parallelizes independent fold-by-algorithm fits inside each bivariate model.

extract_esm_prediction_table <- function(model, n_expected) {
  if (!inherits(model, "BIOMOD.models.out")) return(NULL)
  
  pred <- tryCatch(
    as.data.frame(
      biomod2::get_predictions(model, model.as.col = TRUE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(pred) || !nrow(pred) || nrow(pred) != n_expected) {
    return(NULL)
  }
  
  pred
}

fit_esm_bivariates_direct <- function(
    bm_data, cv_biomod, algorithms, target_name, nb_cpu) {
  
  expected_names <- paste0(
    "_allData_RUN",
    seq_len(ncol(cv_biomod))
  )
  
  if (!identical(colnames(cv_biomod), expected_names)) {
    stop(
      "ESM CV table names were not the expected BIOMOD2 run names. Expected: ",
      paste(expected_names, collapse = ", "),
      "; received: ",
      paste(colnames(cv_biomod), collapse = ", ")
    )
  }
  
  if (any(rowSums(!cv_biomod) != 1L)) {
    stop("Each ESM record must be held out in exactly one spatial fold.")
  }
  
  storage.mode(cv_biomod) <- "logical"
  
  predictor_names <- names(bm_data@data.env.var)
  combinations <- utils::combn(seq_along(predictor_names), 2L)
  prediction_tables <- vector("list", ncol(combinations))
  pair_warnings <- character()
  cv_verified <- FALSE
  
  modeling_id <- paste0("SpatialCV_", sanitize(target_name))
  esm_workdir <- file.path(
    getwd(),
    "ecospat_models",
    paste0("ESM.BIOMOD.output_", sanitize(target_name))
  )
  dir.create(esm_workdir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair in seq_len(ncol(combinations))) {
    mydata <- bm_data
    mydata@dir.name <- esm_workdir
    mydata@data.env.var <- bm_data@data.env.var[
      ,
      combinations[, pair],
      drop = FALSE
    ]
    mydata@sp.name <- paste("ESM.BIOMOD", pair, sep = ".")
    
    bivariate_model <- tryCatch(
      biomod2::BIOMOD_Modeling(
        bm.format = mydata,
        models = algorithms,
        CV.strategy = "user.defined",
        CV.nb.rep = ncol(cv_biomod),
        CV.user.table = cv_biomod,
        CV.do.full.models = FALSE,
        OPT.strategy = "bigboss",
        metric.eval = "AUCroc",
        prevalence = 0.5,
        var.import = 0,
        modeling.id = modeling_id,
        nb.cpu = nb_cpu,
        do.progress = FALSE
      ),
      error = function(e) {
        pair_warnings <<- c(
          pair_warnings,
          sprintf(
            "ESM bivariate pair %d (%s + %s) failed: %s",
            pair,
            predictor_names[combinations[1, pair]],
            predictor_names[combinations[2, pair]],
            conditionMessage(e)
          )
        )
        NULL
      }
    )
    
    if (is.null(bivariate_model)) next
    
    if (!cv_verified) {
      # get_calib_lines(as.data.frame = TRUE) returns BIOMOD2's long-format
      # representation (PA/run/points/calib.lines), not the original wide
      # calibration matrix. The default FALSE returns the stored wide table,
      # which can be compared directly with the user-supplied spatial CV matrix.
      fitted_cv <- biomod2::get_calib_lines(
        bivariate_model,
        as.data.frame = FALSE
      )
      fitted_cv <- as.matrix(fitted_cv)
      storage.mode(fitted_cv) <- "logical"
      
      requested_cv <- as.matrix(cv_biomod)
      storage.mode(requested_cv) <- "logical"
      
      if (!identical(dim(fitted_cv), dim(requested_cv)) ||
          !identical(colnames(fitted_cv), expected_names) ||
          !identical(
            unname(fitted_cv),
            unname(requested_cv)
          )) {
        stop(
          "BIOMOD2 did not preserve the requested ESM spatial CV table exactly."
        )
      }
      cv_verified <- TRUE
    }
    
    pred <- extract_esm_prediction_table(
      bivariate_model,
      n_expected = nrow(bm_data@data.env.var)
    )
    
    if (is.null(pred)) {
      pair_warnings <- c(
        pair_warnings,
        sprintf(
          "ESM bivariate pair %d (%s + %s) returned no complete prediction table.",
          pair,
          predictor_names[combinations[1, pair]],
          predictor_names[combinations[2, pair]]
        )
      )
      next
    }
    
    allrun_columns <- grep(
      "allRun|Full",
      names(pred),
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(allrun_columns)) {
      stop(
        "Unexpected full-data ESM prediction column(s) were fitted despite ",
        "CV.do.full.models = FALSE: ",
        paste(allrun_columns, collapse = ", ")
      )
    }
    
    prediction_tables[[pair]] <- pred
  }
  
  prediction_tables <- Filter(Negate(is.null), prediction_tables)
  
  if (!length(prediction_tables)) {
    stop("No ESM bivariate BIOMOD models were successfully fitted.")
  }
  
  list(
    prediction_tables = prediction_tables,
    warnings = clean_warning_text(pair_warnings),
    n_successful_pairs = length(prediction_tables),
    n_expected_pairs = ncol(combinations)
  )
}

build_esm_equal_weight_oof <- function(
    prediction_tables, cv_biomod, algorithm, n_rows) {
  
  k <- ncol(cv_biomod)
  oof_prediction <- rep(NA_real_, n_rows)
  counts_by_fold <- integer(k)
  minimum_row_counts <- integer(k)
  local_warnings <- character()
  
  for (fold in seq_len(k)) {
    pattern <- paste0("_RUN", fold, "_", algorithm, "$")
    
    fold_vectors <- lapply(prediction_tables, function(pred) {
      matched <- grep(
        pattern,
        names(pred),
        value = TRUE,
        ignore.case = TRUE
      )
      
      if (length(matched) > 1L) {
        stop(
          "Multiple prediction columns matched ", pattern,
          " in one bivariate model: ",
          paste(matched, collapse = ", ")
        )
      }
      if (!length(matched)) return(NULL)
      
      values <- normalize_biomod_score(pred[[matched]])
      if (length(values) != n_rows) {
        stop("Bivariate prediction length did not match the model table.")
      }
      
      values
    })
    
    fold_vectors <- Filter(Negate(is.null), fold_vectors)
    counts_by_fold[fold] <- length(fold_vectors)
    
    if (length(fold_vectors) < 2L) {
      local_warnings <- c(
        local_warnings,
        sprintf(
          "%s ESM fold %d had only %d successfully fitted bivariate model(s).",
          algorithm, fold, length(fold_vectors)
        )
      )
      next
    }
    
    prediction_matrix <- do.call(cbind, fold_vectors)
    held_out <- !cv_biomod[, fold]
    
    finite_counts <- rowSums(
      is.finite(prediction_matrix[held_out, , drop = FALSE])
    )
    minimum_row_counts[fold] <- min(finite_counts)
    
    if (any(finite_counts < 2L)) {
      local_warnings <- c(
        local_warnings,
        sprintf(
          "%s ESM fold %d had held-out records with fewer than two finite bivariate predictions.",
          algorithm,
          fold
        )
      )
      next
    }
    
    fold_prediction <- rowMeans(
      prediction_matrix[held_out, , drop = FALSE],
      na.rm = TRUE
    )
    fold_prediction[!is.finite(fold_prediction)] <- NA_real_
    oof_prediction[held_out] <- fold_prediction
  }
  
  list(
    prediction = oof_prediction,
    counts_by_fold = counts_by_fold,
    minimum_row_counts = minimum_row_counts,
    complete = all(is.finite(oof_prediction)) &&
      all(counts_by_fold >= 2L),
    warnings = local_warnings
  )
}

evaluate_esm <- function(
    model_data, env_data, cv_biomod,
    target_name, local_seed, nb_cpu) {
  
  set.seed(local_seed)
  
  bm_data <- biomod2::BIOMOD_FormatingData(
    resp.var = model_data$pr_ab,
    expl.var = as.data.frame(env_data),
    resp.xy = as.data.frame(
      model_data[, .(longitude, latitude)]
    ),
    resp.name = sanitize(target_name)
  )
  
  algorithms <- c("GLM", "GAM", "MAXNET")
  
  direct_fit <- fit_esm_bivariates_direct(
    bm_data = bm_data,
    cv_biomod = cv_biomod,
    algorithms = algorithms,
    target_name = target_name,
    nb_cpu = nb_cpu
  )
  
  prediction_tables <- direct_fit$prediction_tables
  algorithm_predictions <- list()
  raw_metrics <- list()
  local_warnings <- direct_fit$warnings
  
  if (direct_fit$n_successful_pairs < direct_fit$n_expected_pairs) {
    local_warnings <- c(
      local_warnings,
      sprintf(
        "ESM fitted %d of %d expected bivariate predictor pairs.",
        direct_fit$n_successful_pairs,
        direct_fit$n_expected_pairs
      )
    )
  }
  
  for (algorithm in algorithms) {
    algorithm_result <- build_esm_equal_weight_oof(
      prediction_tables = prediction_tables,
      cv_biomod = cv_biomod,
      algorithm = algorithm,
      n_rows = nrow(model_data)
    )
    
    local_warnings <- c(
      local_warnings,
      algorithm_result$warnings
    )
    
    if (!algorithm_result$complete) {
      local_warnings <- c(
        local_warnings,
        sprintf(
          "%s ESM excluded because complete out-of-fold predictions were not available; bivariate counts by fold: %s.",
          algorithm,
          paste(
            paste0(
              "RUN",
              seq_along(algorithm_result$counts_by_fold),
              "=",
              algorithm_result$counts_by_fold
            ),
            collapse = ", "
          )
        )
      )
      next
    }
    
    prediction <- algorithm_result$prediction
    algorithm_predictions[[algorithm]] <- prediction
    
    metrics <- calculate_oof_metrics(
      response = model_data$pr_ab,
      prediction = prediction
    )
    
    metrics[, `:=`(
      component_type = "ESM_algorithm",
      component = algorithm,
      n_successful_bivariate_models =
        min(algorithm_result$counts_by_fold),
      n_bivariate_min =
        min(algorithm_result$counts_by_fold),
      n_bivariate_median =
        median(algorithm_result$counts_by_fold),
      n_bivariate_max =
        max(algorithm_result$counts_by_fold),
      fold_bivariate_counts = paste(
        paste0(
          "RUN",
          seq_along(algorithm_result$counts_by_fold),
          "=",
          algorithm_result$counts_by_fold
        ),
        collapse = ";"
      )
    )]
    
    raw_metrics[[length(raw_metrics) + 1L]] <- metrics
  }
  
  n_algorithms <- length(algorithm_predictions)
  
  if (n_algorithms < MIN_ESM_ALGORITHMS) {
    stop(
      "Only ", n_algorithms,
      " ESM algorithms produced complete equal-weight out-of-fold ",
      "predictions; at least ", MIN_ESM_ALGORITHMS,
      " are required. Details: ",
      paste(unique(local_warnings), collapse = " | ")
    )
  }
  
  prediction_matrix <- do.call(cbind, algorithm_predictions)
  finite_algorithm_counts <- rowSums(is.finite(prediction_matrix))
  
  if (any(finite_algorithm_counts < MIN_ESM_ALGORITHMS)) {
    stop(
      "Some records had predictions from fewer than ",
      MIN_ESM_ALGORITHMS,
      " successful ESM algorithms."
    )
  }
  
  ensemble_prediction <- rowMeans(
    prediction_matrix,
    na.rm = TRUE
  )
  ensemble_prediction[!is.finite(ensemble_prediction)] <- NA_real_
  
  if (any(!is.finite(ensemble_prediction))) {
    stop(
      "The equal-weight ESM ensemble contained non-finite ",
      "out-of-fold predictions."
    )
  }
  
  if (any(rowSums(!cv_biomod) != 1L)) {
    stop("Each record must be held out in exactly one spatial fold.")
  }
  
  heldout_fold <- max.col(
    !cv_biomod,
    ties.method = "first"
  )
  
  ensemble_metrics <- calculate_oof_metrics(
    response = model_data$pr_ab,
    prediction = ensemble_prediction
  )
  
  ensemble_metrics[, `:=`(
    component_type = "ESM_ensemble",
    component = "equal_mean_bivariates_and_algorithms",
    n_successful_bivariate_models = NA_integer_,
    n_bivariate_min = NA_real_,
    n_bivariate_median = NA_real_,
    n_bivariate_max = NA_real_,
    fold_bivariate_counts = NA_character_
  )]
  
  raw <- rbindlist(
    c(raw_metrics, list(ensemble_metrics)),
    fill = TRUE,
    use.names = TRUE
  )
  
  summary <- ensemble_metrics[
    metric %chin% c("AUC", "TSS", "BOYCE"),
    .(
      metric,
      mean = value,
      sd = NA_real_,
      n_folds_valid = uniqueN(
        heldout_fold[is.finite(ensemble_prediction)]
      ),
      n_folds_expected = ncol(cv_biomod),
      calibration_mean = NA_real_,
      optimism_gap = NA_real_,
      evaluation_model =
        "ESM_equal_mean_bivariates_algorithms"
    )
  ]
  
  if (summary$n_folds_valid[1] <
      summary$n_folds_expected[1]) {
    stop("Pooled ESM predictions did not cover every spatial fold.")
  }
  
  if (!is.finite(
    summary[metric == "BOYCE", mean][1]
  )) {
    local_warnings <- c(
      local_warnings,
      "Pooled spatial-CV Boyce was not finite for this ESM."
    )
  }
  
  list(
    summary = summary,
    raw = raw,
    warnings = clean_warning_text(local_warnings),
    n_algorithms_success = n_algorithms
  )
}

# ---- Common-tier BIOMOD2 evaluation -----------------------------------------

evaluate_sdm <- function(
    model_data, env_data, cv_biomod,
    target_name, local_seed, nb_cpu) {
  
  set.seed(local_seed)
  
  bm_data <- biomod2::BIOMOD_FormatingData(
    resp.var = model_data$pr_ab,
    expl.var = as.data.frame(env_data),
    resp.xy = as.data.frame(
      model_data[, .(longitude, latitude)]
    ),
    resp.name = sanitize(target_name),
    PA.strategy = NULL
  )
  
  metrics <- c("TSS", "AUCroc", "BOYCE")
  algorithms <- c("GLM", "GAM", "RF", "MAXNET")
  predictor_count <- ncol(env_data)
  
  user_values <- list(
    RF.binary.randomForest.randomForest = list(
      for_all_datasets = list(
        do.classif = TRUE,
        ntree = 800,
        mtry = max(
          1,
          floor(sqrt(predictor_count))
        ),
        nodesize = 5
      )
    ),
    MAXNET.binary.maxnet.maxnet = list(
      for_all_datasets = list(
        classes = "lq",
        regmult = 1.0
      )
    ),
    GAM.binary.mgcv.bam = list(
      for_all_datasets = list(
        select = TRUE,
        gamma = 1.4
      )
    )
  )
  
  model_options <- biomod2::bm_ModelingOptions(
    data.type = "binary",
    models = algorithms,
    strategy = "user.defined",
    user.base = "bigboss",
    user.val = user_values,
    bm.format = bm_data,
    calib.lines = cv_biomod
  )
  
  single_models <- biomod2::BIOMOD_Modeling(
    bm.format = bm_data,
    modeling.id = paste0(
      "SpatialCV_",
      sanitize(target_name)
    ),
    models = algorithms,
    OPT.user = model_options,
    CV.strategy = "user.defined",
    CV.user.table = cv_biomod,
    metric.eval = metrics,
    CV.do.full.models = FALSE,
    var.import = 0,
    nb.cpu = nb_cpu,
    seed.val = local_seed
  )
  
  built_models <- biomod2::get_built_models(single_models)
  
  if (!length(built_models)) {
    stop("BIOMOD2 did not build any common-tier models.")
  }
  
  models_per_run <- vapply(
    seq_len(ncol(cv_biomod)),
    function(run) {
      sum(grepl(
        paste0(
          "_RUN", run,
          "_(GLM|GAM|RF|MAXNET)$"
        ),
        built_models,
        ignore.case = TRUE
      ))
    },
    integer(1)
  )
  
  if (any(models_per_run < 2L)) {
    stop(
      "Fewer than two common-tier algorithms were fitted in at least ",
      "one fold: ",
      paste(
        paste0(
          "RUN",
          seq_along(models_per_run),
          "=",
          models_per_run
        ),
        collapse = ", "
      )
    )
  }
  
  ensemble_models <- biomod2::BIOMOD_EnsembleModeling(
    bm.mod = single_models,
    models.chosen = built_models,
    em.by = "PA+run",
    em.algo = "EMmean",
    # No validation-metric filtering: all successfully fitted algorithms
    # contribute equally to the fold-specific ensemble.
    metric.select = NULL,
    metric.eval = metrics,
    var.import = 0,
    nb.cpu = nb_cpu,
    seed.val = local_seed,
    do.progress = FALSE
  )
  
  single_evaluations <- as.data.table(
    biomod2::get_evaluations(single_models)
  )
  
  ensemble_evaluations <- as.data.table(
    biomod2::get_evaluations(ensemble_models)
  )
  
  if (!nrow(ensemble_evaluations)) {
    stop("BIOMOD2 ensemble evaluation table was empty.")
  }
  
  if ("algo" %in% names(ensemble_evaluations)) {
    keep <- grepl(
      "^EMmean$",
      ensemble_evaluations$algo,
      ignore.case = TRUE
    )
    if (any(keep)) {
      ensemble_evaluations <- ensemble_evaluations[keep]
    }
  }
  
  metric_column <- intersect(
    c("metric.eval", "metric"),
    names(ensemble_evaluations)
  )[1]
  
  if (is.na(metric_column)) {
    stop(
      "Could not identify the metric column in BIOMOD2 evaluations."
    )
  }
  
  if (!"validation" %in% names(ensemble_evaluations)) {
    stop(
      "BIOMOD2 ensemble evaluations lacked a validation column."
    )
  }
  
  ensemble_evaluations[
    ,
    value := normalize_biomod_score(validation)
  ]
  
  ensemble_evaluations[
    ,
    calibration_value := if (
      "calibration" %in% names(ensemble_evaluations)
    ) {
      normalize_biomod_score(calibration)
    } else {
      NA_real_
    }
  ]
  
  ensemble_evaluations[
    ,
    metric_std := standardize_metric(get(metric_column))
  ]
  
  run_column <- intersect(
    c("merged.by.run", "run"),
    names(ensemble_evaluations)
  )[1]
  
  if (is.na(run_column)) {
    stop(
      "Could not identify the fold/run column in BIOMOD2 ",
      "ensemble evaluations."
    )
  }
  
  ensemble_evaluations[
    ,
    fold_key := as.character(get(run_column))
  ]
  
  ensemble_evaluations <- ensemble_evaluations[
    !grepl(
      "allRun|Full",
      fold_key,
      ignore.case = TRUE
    )
  ]
  
  duplicate_check <- ensemble_evaluations[
    metric_std %chin% c("AUC", "TSS", "BOYCE"),
    .(
      n_rows = .N,
      n_validation_values =
        uniqueN(round(value, 12), na.rm = TRUE),
      n_calibration_values =
        uniqueN(round(calibration_value, 12), na.rm = TRUE)
    ),
    by = .(fold_key, metric_std)
  ]
  
  if (any(duplicate_check$n_validation_values > 1L) ||
      any(duplicate_check$n_calibration_values > 1L)) {
    stop(
      "BIOMOD2 returned multiple distinct EMmean variants for the ",
      "same fold and evaluation metric."
    )
  }
  
  setorder(ensemble_evaluations, fold_key, metric_std)
  
  ensemble_evaluations <- unique(
    ensemble_evaluations,
    by = c("fold_key", "metric_std")
  )
  
  summary <- ensemble_evaluations[
    metric_std %chin% c("AUC", "TSS", "BOYCE") &
      is.finite(value),
    .(
      mean = safe_mean(value),
      sd = safe_sd(value),
      n_folds_valid =
        uniqueN(fold_key[is.finite(value)]),
      n_folds_expected = ncol(cv_biomod),
      calibration_mean =
        safe_mean(calibration_value),
      optimism_gap =
        safe_mean(calibration_value) -
        safe_mean(value),
      evaluation_model =
        "BIOMOD2_EMmean_equal_weight"
    ),
    by = .(metric = metric_std)
  ]
  
  if (!nrow(summary)) {
    stop(
      "No finite AUC, TSS, or Boyce validation metrics were ",
      "returned for EMmean."
    )
  }
  
  missing_primary <- setdiff(
    c("AUC", "TSS"),
    summary$metric
  )
  
  if (length(missing_primary)) {
    stop(
      "Missing required common-tier validation metrics: ",
      paste(missing_primary, collapse = ", ")
    )
  }
  
  incomplete_primary <- summary[
    metric %chin% c("AUC", "TSS") &
      n_folds_valid < n_folds_expected
  ]
  
  if (nrow(incomplete_primary)) {
    stop(
      "AUC/TSS were not returned for every spatial fold: ",
      paste(
        paste0(
          incomplete_primary$metric,
          "=",
          incomplete_primary$n_folds_valid,
          "/",
          incomplete_primary$n_folds_expected
        ),
        collapse = ", "
      )
    )
  }
  
  local_warnings <- character()
  
  if (!"BOYCE" %in% summary$metric ||
      !is.finite(
        summary[metric == "BOYCE", mean][1]
      )) {
    local_warnings <- c(
      local_warnings,
      "BIOMOD2 did not return a finite spatial-CV Boyce value ",
      "for the ensemble."
    )
  }
  
  single_evaluations[, `:=`(
    component_type = "BIOMOD_single_algorithm",
    component = if ("algo" %in% names(single_evaluations)) {
      as.character(algo)
    } else {
      NA_character_
    },
    metric = if (
      "metric.eval" %in% names(single_evaluations)
    ) {
      standardize_metric(metric.eval)
    } else {
      NA_character_
    },
    value = if (
      "validation" %in% names(single_evaluations)
    ) {
      normalize_biomod_score(validation)
    } else {
      NA_real_
    }
  )]
  
  ensemble_evaluations[, `:=`(
    component_type = "BIOMOD_ensemble",
    component = "EMmean",
    metric = metric_std
  )]
  
  raw <- rbindlist(
    list(single_evaluations, ensemble_evaluations),
    fill = TRUE,
    use.names = TRUE
  )
  
  n_algorithms_success <- if (
    all(c("algo", "validation") %in% names(single_evaluations))
  ) {
    uniqueN(single_evaluations[
      is.finite(normalize_biomod_score(validation)),
      algo
    ])
  } else {
    NA_integer_
  }
  
  list(
    summary = summary,
    raw = raw,
    warnings = clean_warning_text(local_warnings),
    n_algorithms_success = n_algorithms_success
  )
}

# ---- Per-species checkpoint writer ------------------------------------------

write_species_outputs <- function(
    row, run_signature, state, stage, elapsed_seconds,
    metrics = NULL, components = NULL,
    fold_assignments = NULL, fold_summary = NULL,
    errors = NULL, warnings = NULL,
    predictor_mapping = NULL) {
  
  out_dir <- species_result_dir(row)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(metrics)) metrics <- species_metrics_prototype
  if (is.null(components)) components <- component_metrics_prototype
  if (is.null(fold_assignments)) {
    fold_assignments <- fold_assignments_prototype
  }
  if (is.null(fold_summary)) fold_summary <- fold_summary_prototype
  if (is.null(errors)) errors <- errors_prototype
  if (is.null(warnings)) warnings <- warnings_prototype
  if (is.null(predictor_mapping)) {
    predictor_mapping <- predictor_mapping_prototype
  }
  
  atomic_fwrite(
    metrics,
    file.path(out_dir, "species_metrics.csv")
  )
  
  atomic_fwrite(
    components,
    file.path(out_dir, "component_metrics.csv")
  )
  
  atomic_fwrite(
    fold_summary,
    file.path(out_dir, "fold_summary.csv")
  )
  
  atomic_fwrite(
    fold_assignments,
    file.path(out_dir, "fold_assignments.csv.gz"),
    compress = "gzip"
  )
  
  atomic_fwrite(
    errors,
    file.path(out_dir, "errors.csv")
  )
  
  atomic_fwrite(
    warnings,
    file.path(out_dir, "warnings.csv")
  )
  
  atomic_fwrite(
    predictor_mapping,
    file.path(out_dir, "predictor_mapping.csv")
  )
  
  status <- list(
    species_id = row$species_id,
    species_key = row$species_key,
    species = row$species,
    taxon = row$taxon,
    model_type = row$model_type,
    state = state,
    stage = stage,
    elapsed_seconds = elapsed_seconds,
    completed_at = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    ),
    run_signature = run_signature,
    script_version = SCRIPT_VERSION,
    package_versions = package_versions
  )
  
  # Status is written last; its presence marks an atomic completed checkpoint.
  atomic_save_rds(
    status,
    status_path_for(row),
    compress = FALSE
  )
  
  status
}

# ---- Per-species execution ---------------------------------------------------

process_species <- function(row_index) {
  row <- manifest[row_index]
  target <- row$species
  taxon <- row$taxon
  model_type <- row$model_type
  local_seed <- species_seed(row$taxon, row$species)
  
  out_dir <- species_result_dir(row)
  work_dir <- species_work_dir(row)
  temp_dir <- species_temp_dir(row)
  
  started <- Sys.time()
  stage <- "initialization"
  captured_warnings <- list()
  fold_assignments <- NULL
  fold_summary <- NULL
  
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (dir.exists(work_dir)) {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  console_log <- file.path(out_dir, "console.log")
  
  old_wd <- getwd()
  old_terra <- terra::terraOptions(print = FALSE)
  
  # The message stream requires an already open connection. Passing the path
  # directly to sink(type = "message") causes every species to fail before
  # occurrence preprocessing or modeling begins.
  log_con <- file(console_log, open = "wt")
  output_sink_active <- FALSE
  message_sink_active <- FALSE
  cleaned <- FALSE
  
  cleanup <- function(success) {
    if (cleaned) return(invisible(NULL))
    cleaned <<- TRUE
    
    if (message_sink_active) {
      try(sink(type = "message"), silent = TRUE)
      message_sink_active <<- FALSE
    }
    if (output_sink_active) {
      try(sink(type = "output"), silent = TRUE)
      output_sink_active <<- FALSE
    }
    try(close(log_con), silent = TRUE)
    try(setwd(old_wd), silent = TRUE)
    
    try(
      terra::tmpFiles(
        current = TRUE,
        orphan = TRUE,
        remove = TRUE
      ),
      silent = TRUE
    )
    
    try(
      terra::terraOptions(
        memfrac = old_terra$memfrac,
        tempdir = old_terra$tempdir,
        progress = old_terra$progress,
        parallel = old_terra$parallel,
        threads = old_terra$threads
      ),
      silent = TRUE
    )
    
    unlink(temp_dir, recursive = TRUE, force = TRUE)
    
    keep_work <- if (success) {
      KEEP_SUCCESSFUL_MODEL_WORK
    } else {
      KEEP_FAILED_MODEL_WORK
    }
    
    if (!keep_work) {
      unlink(work_dir, recursive = TRUE, force = TRUE)
    }
    
    keep_console <- if (success) {
      KEEP_SUCCESSFUL_CONSOLE_LOGS
    } else {
      KEEP_FAILED_CONSOLE_LOGS
    }
    
    if (!keep_console) {
      unlink(console_log, force = TRUE)
    }
    
    gc(verbose = FALSE)
    invisible(NULL)
  }
  
  on.exit(cleanup(FALSE), add = TRUE)
  
  sink(log_con, type = "output", split = FALSE)
  output_sink_active <- TRUE
  sink(log_con, type = "message")
  message_sink_active <- TRUE
  
  result <- withCallingHandlers(
    tryCatch({
      setwd(work_dir)
      
      terra::terraOptions(
        memfrac = TERRA_MEMFRAC,
        tempdir = temp_dir,
        progress = 0,
        parallel = FALSE,
        threads = TERRA_THREADS
      )
      
      set.seed(local_seed)
      
      # Reopen file-backed rasters for each species. This keeps terra state
      # isolated and allows aggressive removal of temporary files after fitting.
      env_stack_local <- terra::rast(ENV_PATH)
      env_stack_local <- env_stack_local[[valid_layer_names]]
      
      autocor_stack_local <- terra::rast(AUTOCOR_CACHE_PATH)
      block_template_cache <- autocor_stack_local[[1]]
      
      stage <- "occurrence preprocessing"
      
      sp_clean <- prepare_species_occurrences(
        target_name = target,
        occ = occ,
        env_stack = env_stack_local,
        expected_model_type = model_type
      )
      
      stage <- "accessible area"
      
      accessible <- make_accessible_area(
        sp_clean,
        ecoregions_by_code
      )
      
      stage <- "predictor matching"
      
      predictor_resolution <- predictor_lookup[[
        paste(taxon, target, sep = "||")
      ]]
      
      if (is.null(predictor_resolution)) {
        stop(
          "No preflight predictor mapping was available for this species."
        )
      }
      
      env_selected <- env_stack_local[[
        predictor_resolution$raster_names
      ]]
      names(env_selected) <- predictor_resolution$model_names
      
      autocor_names <- intersect(
        predictor_resolution$raster_names,
        names(autocor_stack_local)
      )
      
      env_autocor <- if (length(autocor_names)) {
        autocor_stack_local[[autocor_names]]
      } else {
        NULL
      }
      
      stage <- "presence environmental extraction"
      
      presence_values <- terra::extract(
        env_selected,
        as.matrix(sp_clean[, .(longitude, latitude)])
      )
      
      presence_keep <- complete.cases(presence_values)
      sp_clean <- sp_clean[presence_keep]
      presence_values <- as.data.frame(
        presence_values[presence_keep, , drop = FALSE]
      )
      
      if (!nrow(sp_clean)) {
        stop(
          "No presences had complete selected environmental predictors."
        )
      }
      
      stage <- "background generation"
      
      n_background_target <- if (model_type == "ESM") {
        N_BG_ESM
      } else {
        N_BG_SDM
      }
      
      background <- generate_target_group_background(
        env_selected = env_selected,
        accessible = accessible,
        background_pool = background_pool_current,
        target_name = target,
        presence_xy = sp_clean,
        n_target = n_background_target,
        local_seed = local_seed + 101L
      )
      
      presences <- copy(sp_clean)
      presences[, pr_ab := 1L]
      
      background_data <- copy(background$data)
      background_data[, `:=`(
        species = target,
        pr_ab = 0L
      )]
      
      model_data <- rbindlist(
        list(presences, background_data),
        fill = TRUE,
        use.names = TRUE
      )
      
      stage <- "model-table environmental assembly"
      
      environmental_data <- rbind(
        presence_values,
        background$env
      )
      environmental_data <- as.data.frame(environmental_data)
      
      if (nrow(environmental_data) != nrow(model_data) ||
          ncol(environmental_data) != terra::nlyr(env_selected)) {
        stop(
          "Environmental model table was not aligned with presence/background rows."
        )
      }
      
      # Both components were already filtered to complete cases before sampling.
      # Keep this defensive check so a future upstream change cannot silently pass
      # missing predictor values into BIOMOD2.
      complete_rows <- complete.cases(environmental_data)
      if (!all(complete_rows)) {
        model_data <- model_data[complete_rows]
        environmental_data <- environmental_data[
          complete_rows,
          ,
          drop = FALSE
        ]
      }
      
      if (sum(model_data$pr_ab == 1L) < 7L) {
        stop(
          "Fewer than seven valid presences remained; ",
          "ESM/SDM evaluation is not appropriate."
        )
      }
      
      if (sum(model_data$pr_ab == 0L) <
          MIN_TEST_BG * MIN_FOLDS) {
        stop(
          "Too few valid background records remained for spatial CV."
        )
      }
      
      stage <- "spatial fold construction"
      
      folds <- make_spatial_folds(
        model_data = model_data,
        env_autocor = env_autocor,
        block_template_cache = block_template_cache,
        accessible = accessible,
        model_type = model_type,
        local_seed = local_seed + 202L
      )
      
      fold_assignments <- data.table(
        species = target,
        taxon = taxon,
        model_type = model_type,
        longitude = model_data$longitude,
        latitude = model_data$latitude,
        pr_ab = as.integer(model_data$pr_ab),
        fold = as.integer(folds$fold_id)
      )
      
      fold_summary <- copy(folds$balance)
      fold_summary[, `:=`(
        species = target,
        taxon = taxon,
        model_type = model_type,
        block_size_km = folds$block_size_m / 1000,
        autocor_range_km =
          folds$autocor_range_m / 1000
      )]
      setcolorder(
        fold_summary,
        names(fold_summary_prototype)
      )
      
      stage <- if (model_type == "ESM") {
        "ESM fitting and pooled evaluation"
      } else {
        "BIOMOD2 fitting and ensemble evaluation"
      }
      
      fit <- if (model_type == "ESM") {
        evaluate_esm(
          model_data = model_data,
          env_data = environmental_data,
          cv_biomod = folds$cv_biomod,
          target_name = target,
          local_seed = local_seed + 303L,
          nb_cpu = BIOMOD_CORES
        )
      } else {
        evaluate_sdm(
          model_data = model_data,
          env_data = environmental_data,
          cv_biomod = folds$cv_biomod,
          target_name = target,
          local_seed = local_seed + 303L,
          nb_cpu = BIOMOD_CORES
        )
      }
      
      fit_warnings <- clean_warning_text(fit$warnings)
      
      if (length(fit_warnings)) {
        captured_warnings[[
          length(captured_warnings) + 1L
        ]] <- data.table(
          stage = stage,
          warning = fit_warnings
        )
      }
      
      metrics <- copy(fit$summary)
      metrics[, `:=`(
        species = target,
        taxon = taxon,
        model_type = model_type,
        n_presences = sum(model_data$pr_ab == 1L),
        n_background = sum(model_data$pr_ab == 0L),
        n_spatial_folds = folds$k,
        block_size_km = folds$block_size_m / 1000,
        autocor_range_km =
          folds$autocor_range_m / 1000,
        predictor_count =
          length(predictor_resolution$model_names),
        background_source = background$source,
        background_candidates =
          background$n_candidates,
        n_algorithms_success =
          fit$n_algorithms_success
      )]
      setcolorder(metrics, names(species_metrics_prototype))
      
      components <- copy(fit$raw)
      components[, `:=`(
        species = target,
        taxon = taxon,
        model_type = model_type
      )]
      
      mapping <- copy(predictor_resolution$mapping)
      mapping[, `:=`(
        species = target,
        taxon = taxon,
        model_type = model_type
      )]
      setcolorder(mapping, names(predictor_mapping_prototype))
      
      list(
        metrics = metrics,
        components = components,
        mapping = mapping
      )
    }, error = function(e) {
      structure(
        list(
          error = conditionMessage(e),
          stage = stage
        ),
        class = "species_error"
      )
    }),
    warning = function(w) {
      warning_text <- conditionMessage(w)
      
      if (!is.na(warning_text) &&
          nzchar(trimws(warning_text)) &&
          trimws(warning_text) != "NaN") {
        captured_warnings[[
          length(captured_warnings) + 1L
        ]] <<- data.table(
          stage = stage,
          warning = warning_text
        )
      }
      
      invokeRestart("muffleWarning")
    }
  )
  
  elapsed <- as.numeric(
    difftime(Sys.time(), started, units = "secs")
  )
  
  warnings_output <- if (length(captured_warnings)) {
    x <- unique(rbindlist(
      captured_warnings,
      fill = TRUE,
      use.names = TRUE
    ))
    x[, `:=`(
      species = target,
      taxon = taxon,
      model_type = model_type
    )]
    setcolorder(x, names(warnings_prototype))
    x
  } else {
    warnings_prototype
  }
  
  if (inherits(result, "species_error")) {
    errors <- data.table(
      species = target,
      taxon = taxon,
      model_type = model_type,
      stage = result$stage,
      error = result$error
    )
    
    status <- write_species_outputs(
      row = row,
      run_signature = RUN_SIGNATURE,
      state = "failed",
      stage = result$stage,
      elapsed_seconds = elapsed,
      fold_assignments = fold_assignments,
      fold_summary = fold_summary,
      errors = errors,
      warnings = warnings_output,
      predictor_mapping = if (
        exists("predictor_resolution", inherits = FALSE) &&
        !is.null(predictor_resolution)
      ) {
        mapping <- copy(predictor_resolution$mapping)
        mapping[, `:=`(
          species = target,
          taxon = taxon,
          model_type = model_type
        )]
        setcolorder(
          mapping,
          names(predictor_mapping_prototype)
        )
        mapping
      } else {
        predictor_mapping_prototype
      }
    )
    
    cleanup(FALSE)
    return(status)
  }
  
  status <- write_species_outputs(
    row = row,
    run_signature = RUN_SIGNATURE,
    state = "complete",
    stage = "complete",
    elapsed_seconds = elapsed,
    metrics = result$metrics,
    components = result$components,
    fold_assignments = fold_assignments,
    fold_summary = fold_summary,
    warnings = warnings_output,
    predictor_mapping = result$mapping
  )
  
  cleanup(TRUE)
  status
}

process_species_safe <- function(row_index) {
  tryCatch(
    process_species(row_index),
    error = function(e) {
      row <- manifest[row_index]
      elapsed <- NA_real_
      fallback_error <- data.table(
        species = row$species,
        taxon = row$taxon,
        model_type = row$model_type,
        stage = "unhandled species execution error",
        error = conditionMessage(e)
      )
      
      status <- tryCatch(
        write_species_outputs(
          row = row,
          run_signature = RUN_SIGNATURE,
          state = "failed",
          stage = "unhandled species execution error",
          elapsed_seconds = elapsed,
          errors = fallback_error
        ),
        error = function(write_error) {
          list(
            species_id = row$species_id,
            species_key = row$species_key,
            species = row$species,
            taxon = row$taxon,
            model_type = row$model_type,
            state = "failed",
            stage = "checkpoint writing failed",
            elapsed_seconds = elapsed,
            completed_at = format(
              Sys.time(),
              "%Y-%m-%dT%H:%M:%S%z"
            ),
            run_signature = RUN_SIGNATURE,
            error = paste(
              conditionMessage(e),
              conditionMessage(write_error),
              sep = " | "
            )
          )
        }
      )
      
      status
    }
  )
}

# ---- Combined-output streaming -----------------------------------------------

csv_header_names <- function(path) {
  names(fread(path, nrows = 0L, showProgress = FALSE))
}

combine_csv_files <- function(files, output_path, prototype) {
  files <- files[file.exists(files)]
  
  all_names <- names(prototype)
  
  if (length(files)) {
    headers <- lapply(files, function(path) {
      tryCatch(
        csv_header_names(path),
        error = function(e) character()
      )
    })
    all_names <- unique(c(all_names, unlist(headers)))
  }
  
  dir.create(
    dirname(output_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  tmp <- file.path(
    dirname(output_path),
    paste0(".", basename(output_path), ".", Sys.getpid(), ".tmp")
  )
  
  # Write an explicit header so empty combined files remain valid CSV files.
  writeLines(
    paste(all_names, collapse = ","),
    tmp,
    useBytes = TRUE
  )
  
  for (path in files) {
    x <- fread(path, showProgress = FALSE)
    
    missing <- setdiff(all_names, names(x))
    if (length(missing)) {
      for (nm in missing) x[, (nm) := NA]
    }
    
    setcolorder(x, all_names)
    
    if (nrow(x)) {
      fwrite(
        x,
        tmp,
        append = TRUE,
        col.names = FALSE
      )
    }
  }
  
  if (file.exists(output_path)) unlink(output_path, force = TRUE)
  if (!file.rename(tmp, output_path)) {
    unlink(tmp, force = TRUE)
    stop("Could not finalize combined CSV: ", output_path)
  }
  
  invisible(output_path)
}

combine_gzip_csv_files <- function(
    files, output_path, header_names) {
  
  files <- files[file.exists(files)]
  
  dir.create(
    dirname(output_path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  tmp <- file.path(
    dirname(output_path),
    paste0(
      ".",
      sub("\\.gz$", "", basename(output_path)),
      ".",
      Sys.getpid(),
      ".tmp.gz"
    )
  )
  
  out_con <- gzfile(tmp, open = "wt")
  on.exit(try(close(out_con), silent = TRUE), add = TRUE)
  
  writeLines(
    paste(header_names, collapse = ","),
    out_con,
    useBytes = TRUE
  )
  
  for (path in files) {
    in_con <- gzfile(path, open = "rt")
    
    tryCatch({
      # Discard each per-species header.
      readLines(in_con, n = 1L, warn = FALSE)
      
      repeat {
        lines <- readLines(
          in_con,
          n = 50000L,
          warn = FALSE
        )
        if (!length(lines)) break
        writeLines(lines, out_con, useBytes = TRUE)
      }
    }, finally = {
      try(close(in_con), silent = TRUE)
    })
  }
  
  close(out_con)
  
  if (file.exists(output_path)) unlink(output_path, force = TRUE)
  if (!file.rename(tmp, output_path)) {
    unlink(tmp, force = TRUE)
    stop(
      "Could not finalize combined compressed CSV: ",
      output_path
    )
  }
  
  invisible(output_path)
}

status_to_row <- function(status) {
  data.table(
    species_id = as.integer(status$species_id),
    species_key = as.character(status$species_key %||% NA_character_),
    species = as.character(status$species),
    taxon = as.character(status$taxon),
    model_type = as.character(status$model_type),
    state = as.character(status$state),
    stage = as.character(status$stage),
    elapsed_seconds = as.numeric(status$elapsed_seconds),
    completed_at = as.character(status$completed_at),
    run_signature = as.character(status$run_signature)
  )
}

collect_current_species_files <- function(rows) {
  dirs <- vapply(
    seq_len(nrow(rows)),
    function(i) species_result_dir(rows[i]),
    character(1)
  )
  
  statuses <- lapply(
    seq_len(nrow(rows)),
    function(i) {
      row <- rows[i]
      status <- read_status(row)
      
      if (is.null(status) ||
          !identical(status$run_signature, RUN_SIGNATURE)) {
        return(NULL)
      }
      
      valid_files <- if (identical(status$state, "complete")) {
        all(file.exists(species_checkpoint_files(row)))
      } else if (identical(status$state, "failed")) {
        file.exists(file.path(species_result_dir(row), "errors.csv"))
      } else {
        FALSE
      }
      
      if (!valid_files) return(NULL)
      status
    }
  )
  
  valid <- !vapply(statuses, is.null, logical(1))
  
  list(
    rows = rows[valid],
    dirs = dirs[valid],
    statuses = statuses[valid]
  )
}

compile_scope_outputs <- function(rows, scope_name) {
  scope_dir <- file.path(COMBINED_DIR, scope_name)
  dir.create(scope_dir, recursive = TRUE, showWarnings = FALSE)
  
  collected <- collect_current_species_files(rows)
  dirs <- collected$dirs
  statuses <- collected$statuses
  
  file_vector <- function(name) file.path(dirs, name)
  
  combine_csv_files(
    file_vector("species_metrics.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_species_metrics.csv"
    ),
    species_metrics_prototype
  )
  
  combine_csv_files(
    file_vector("component_metrics.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_component_metrics.csv"
    ),
    component_metrics_prototype
  )
  
  combine_csv_files(
    file_vector("fold_summary.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_fold_summary.csv"
    ),
    fold_summary_prototype
  )
  
  combine_csv_files(
    file_vector("errors.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_errors.csv"
    ),
    errors_prototype
  )
  
  combine_csv_files(
    file_vector("warnings.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_warnings.csv"
    ),
    warnings_prototype
  )
  
  combine_csv_files(
    file_vector("predictor_mapping.csv"),
    file.path(
      scope_dir,
      "sdm_spatial_cv_predictor_mapping.csv"
    ),
    predictor_mapping_prototype
  )
  
  status_table <- if (length(statuses)) {
    rbindlist(
      lapply(statuses, status_to_row),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    status_prototype
  }
  
  atomic_fwrite(
    status_table,
    file.path(
      scope_dir,
      "sdm_spatial_cv_status.csv"
    )
  )
  
  expected_coverage <- rows[
    , .(n_expected = .N),
    by = .(taxon, model_type)
  ]
  
  observed_coverage <- if (nrow(status_table)) {
    dcast(
      status_table[
        , .N,
        by = .(taxon, model_type, state)
      ],
      taxon + model_type ~ state,
      value.var = "N",
      fill = 0
    )
  } else {
    data.table(
      taxon = character(),
      model_type = character()
    )
  }
  
  coverage <- merge(
    expected_coverage,
    observed_coverage,
    by = c("taxon", "model_type"),
    all.x = TRUE,
    sort = FALSE
  )
  
  for (column in c("complete", "failed")) {
    if (!column %in% names(coverage)) coverage[, (column) := 0L]
    set(coverage, which(is.na(coverage[[column]])), column, 0L)
  }
  
  coverage[, `:=`(
    n_complete = as.integer(complete),
    n_failed = as.integer(failed)
  )]
  coverage[, n_missing := pmax(
    0L,
    n_expected - n_complete - n_failed
  )]
  coverage[, completion_fraction := ifelse(
    n_expected > 0L,
    n_complete / n_expected,
    NA_real_
  )]
  coverage[, c("complete", "failed") := NULL]
  
  atomic_fwrite(
    coverage,
    file.path(
      scope_dir,
      "sdm_spatial_cv_evaluation_coverage.csv"
    )
  )
  
  if (WRITE_COMBINED_FOLD_ASSIGNMENTS) {
    combine_gzip_csv_files(
      file_vector("fold_assignments.csv.gz"),
      file.path(
        scope_dir,
        "sdm_spatial_cv_fold_assignments.csv.gz"
      ),
      names(fold_assignments_prototype)
    )
  }
  
  species_metrics_path <- file.path(
    scope_dir,
    "sdm_spatial_cv_species_metrics.csv"
  )
  
  species_metrics <- fread(
    species_metrics_path,
    showProgress = FALSE
  )
  
  publication_summary <- if (nrow(species_metrics)) {
    species_metrics[
      is.finite(mean),
      .(
        n_species = uniqueN(species),
        median = stats::median(mean, na.rm = TRUE),
        q25 = as.numeric(
          stats::quantile(mean, 0.25, na.rm = TRUE)
        ),
        q75 = as.numeric(
          stats::quantile(mean, 0.75, na.rm = TRUE)
        ),
        mean = mean(mean, na.rm = TRUE),
        sd_among_species = safe_sd(mean)
      ),
      by = .(
        taxon,
        model_type,
        evaluation_model,
        metric
      )
    ]
  } else {
    data.table(
      taxon = character(),
      model_type = character(),
      evaluation_model = character(),
      metric = character(),
      n_species = integer(),
      median = numeric(),
      q25 = numeric(),
      q75 = numeric(),
      mean = numeric(),
      sd_among_species = numeric()
    )
  }
  
  atomic_fwrite(
    publication_summary,
    file.path(
      scope_dir,
      "sdm_spatial_cv_publication_summary.csv"
    )
  )
  
  completed <- sum(status_table$state == "complete")
  failed <- sum(status_table$state == "failed")
  
  metadata <- c(
    paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Script version:", SCRIPT_VERSION),
    paste("Run signature:", RUN_SIGNATURE),
    paste("Scope:", scope_name),
    paste("Root:", ROOT),
    paste("Occurrence file:", OCC_PATH),
    paste("Environmental stack:", ENV_PATH),
    paste("Ecoregion layer:", ECO_PATH),
    paste("Manifest directory:", MANIFEST_DIR),
    paste("Projected autocorrelation cache:", AUTOCOR_CACHE_PATH),
    paste("Target taxa:", paste(TARGET_TAXA, collapse = ", ")),
    paste(
      "Target species:",
      if (is.null(TARGET_SPECIES)) {
        "ALL"
      } else {
        paste(TARGET_SPECIES, collapse = ", ")
      }
    ),
    paste("Batch:", BATCH_ID, "of", BATCH_COUNT),
    paste("Species workers:", SPECIES_WORKERS),
    paste("BIOMOD2 cores per species:", BIOMOD_CORES),
    paste("terra memfrac per process:", TERRA_MEMFRAC),
    paste("Species in current manifest scope:", nrow(rows)),
    paste("Species with current checkpoints:", nrow(status_table)),
    paste("Species completed:", completed),
    paste("Species failed:", failed),
    paste(
      "Package versions:",
      paste(
        names(package_versions),
        package_versions,
        sep = "=",
        collapse = "; "
      )
    ),
    "",
    "Interpretation:",
    "- Boyce is the primary metric.",
    "- AUC and maximum TSS are secondary and conditional on the sampled background domain.",
    "- Predictors are fixed from the final manifests rather than reselected inside folds.",
    "- Common-tier performance uses equal-weight, fold-specific EMmean ensembles.",
    "- ESM performance pools strict out-of-fold predictions with fixed equal weights.",
    "- ESM bivariates use BIOMOD2 bigboss settings matching ecospat 4.1.4, with no unused full-data fit.",
    "- Ultra-rare alpha-hull species are intentionally excluded."
  )
  
  writeLines(
    metadata,
    file.path(
      scope_dir,
      "sdm_spatial_cv_run_metadata.txt"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(
      scope_dir,
      "sdm_spatial_cv_sessionInfo.txt"
    )
  )
  
  invisible(list(
    scope = scope_name,
    completed = completed,
    failed = failed,
    output_dir = scope_dir
  ))
}


# ---- Run requested mode -------------------------------------------------------

if (RUN_MODE == "local") {
  message(
    "LOCAL TEST: ", LOCAL_TEST_N_PER_TAXON,
    " species per taxon; output = ",
    OUTPUT_ROOT
  )
  
  Sys.setenv(
    SDM_LOCAL_TEST = "1",
    SDM_ROOT = ROOT,
    SDM_CV_OUTPUT_ROOT = OUTPUT_ROOT
  )
  
  rscript <- file.path(R.home("bin"), "Rscript")
  
  run_child <- function(args) {
    status <- system2(
      rscript,
      args = c(shQuote(SCRIPT_FILE), args)
    )
    if (!identical(status, 0L)) {
      stop("Local test subprocess failed: ", paste(args, collapse = " "))
    }
  }
  
  run_child("prep-test")
  
  task_count_file <- file.path(
    OUTPUT_ROOT,
    "array",
    "task_count.txt"
  )
  n_tasks <- as.integer(readLines(task_count_file, n = 1L))
  
  for (task_id in seq_len(n_tasks)) {
    message("\n===== LOCAL TASK ", task_id, " OF ", n_tasks, " =====")
    run_child(as.character(task_id))
  }
  
  run_child("compile")
  
  message(
    "\nLocal test complete. Results: ",
    file.path(OUTPUT_ROOT, "combined", "all_taxa")
  )
  
} else if (RUN_MODE %in% c("prepare", "prepare_test")) {
  
  TEST_N_PER_TAXON <- if (RUN_MODE == "prepare_test") {
    LOCAL_TEST_N_PER_TAXON
  } else {
    NULL
  }
  
  ENV_PATH <- first_existing(ENV_CANDIDATES, "Environmental stack")
  ECO_PATH <- first_existing(ECO_CANDIDATES, "Ecoregion layer")
  
  if (!file.exists(OCC_PATH)) stop("Occurrence file not found: ", OCC_PATH)
  if (!dir.exists(MANIFEST_DIR)) {
    stop("Manifest directory not found: ", MANIFEST_DIR)
  }
  
  message("Preparing Slurm array inputs under: ", OUTPUT_ROOT)
  
  manifest_all <- read_manifests(MANIFEST_DIR)
  manifest_all <- manifest_all[
    taxon %chin% TARGET_TAXA &
      model_type %chin% c("ESM", "SDM")
  ]
  
  if (!nrow(manifest_all)) {
    stop("No ESM or SDM manifest rows matched TARGET_TAXA.")
  }
  
  setorder(manifest_all, taxon, model_type, species)
  manifest_all[, species_id := seq_len(.N)]
  manifest_all[, species_key := sprintf(
    "%010d",
    vapply(
      paste(taxon, species, sep = "||"),
      stable_hash_int,
      integer(1)
    )
  )]
  
  key_collisions <- manifest_all[
    , .(n_species = uniqueN(paste(taxon, species, sep = "||"))),
    by = species_key
  ][n_species > 1L]
  if (nrow(key_collisions)) {
    stop("Stable species-key collision detected; revise stable_hash_int().")
  }
  
  # Optional deterministic smoke-test subset. IDs, keys, seeds, and checkpoint
  # paths remain identical to the later full run.
  full_manifest <- copy(manifest_all)
  
  if (!is.null(TEST_N_PER_TAXON)) {
    set.seed(SEED)
    
    pick_rows <- function(idx, n) {
      n <- min(as.integer(n), length(idx))
      if (n < 1L) return(integer())
      idx[sample.int(length(idx), n)]
    }
    
    full_manifest <- full_manifest[, {
      n <- min(TEST_N_PER_TAXON, .N)
      keep <- c(
        pick_rows(which(model_type == "ESM"), ceiling(n / 2)),
        pick_rows(which(model_type == "SDM"), floor(n / 2))
      )
      if (length(keep) < n) {
        remaining <- setdiff(seq_len(.N), keep)
        keep <- c(keep, pick_rows(remaining, n - length(keep)))
      }
      .SD[sort(keep)]
    }, by = taxon]
    
    setorder(full_manifest, taxon, model_type, species)
  }
  
  if (!nrow(full_manifest)) stop("No species remained after filtering.")
  
  # Build a run signature that is independent of task/chunk assignment, so species
  # checkpoints remain reusable if chunks are later rebalanced.
  signature_parts <- c(
    paste("script", SCRIPT_VERSION, sep = "="),
    file_signature(OCC_PATH),
    file_signature(ENV_PATH),
    file_signature(ECO_PATH),
    paste(
      sort(unique(manifest_all$manifest_file)),
      vapply(
        sort(unique(manifest_all$manifest_file)),
        file_signature,
        character(1)
      ),
      sep = "="
    ),
    paste(names(package_versions), package_versions, sep = "="),
    paste("seed", SEED, sep = "="),
    paste("n_bg_esm", N_BG_ESM, sep = "="),
    paste("n_bg_sdm", N_BG_SDM, sep = "="),
    paste("folds", MIN_FOLDS, MAX_FOLDS, sep = "="),
    paste(
      "fold_minimums",
      MIN_TEST_PRES_ESM,
      MIN_TEST_PRES_SDM,
      MIN_TEST_BG,
      MIN_TRAIN_PRES_ESM,
      MIN_TRAIN_PRES_SDM,
      sep = "="
    ),
    paste(
      "block",
      BLOCK_ITERATIONS,
      AUTOCOR_RANGE_MULTIPLIER,
      MAX_BLOCK_EXTENT_FRACTION,
      MIN_BLOCK_CELLS,
      AUTOCOR_SAMPLE_SIZE,
      sep = "="
    ),
    paste(
      "model_design",
      "ESM_direct_BIOMOD2_no_full_models_equal_mean_bivariates_algorithms",
      "BIOMOD2_EMmean_equal_weight_unfiltered",
      sep = "="
    ),
    paste(
      "allow_env_filter_fallback",
      ALLOW_ENV_FILTER_FALLBACK,
      sep = "="
    )
  )
  RUN_SIGNATURE <- md5_text(signature_parts)
  message("Run signature: ", RUN_SIGNATURE)
  
  # Valid environmental layers.
  env_stack_master <- terra::rast(ENV_PATH)
  valid_layer_cache <- file.path(CACHE_DIR, "valid_environment_layers.rds")
  valid_layer_signature <- file_signature(ENV_PATH)
  
  valid_layer_metadata <- if (!FORCE_RERUN && file.exists(valid_layer_cache)) {
    tryCatch(readRDS(valid_layer_cache), error = function(e) NULL)
  } else {
    NULL
  }
  
  if (!is.null(valid_layer_metadata) &&
      identical(valid_layer_metadata$signature, valid_layer_signature) &&
      all(valid_layer_metadata$layer_names %in% names(env_stack_master))) {
    valid_layer_names <- valid_layer_metadata$layer_names
  } else {
    valid_layer_indices <- which(
      terra::global(env_stack_master, "max", na.rm = TRUE)[, 1] != -Inf
    )
    valid_layer_names <- names(env_stack_master)[valid_layer_indices]
    
    atomic_save_rds(
      list(
        signature = valid_layer_signature,
        layer_names = valid_layer_names,
        created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ),
      valid_layer_cache,
      compress = FALSE
    )
  }
  
  if (length(valid_layer_names) < 2L) {
    stop("Environmental stack contains fewer than two valid layers.")
  }
  env_stack_master <- env_stack_master[[valid_layer_names]]
  
  # Shared spatial caches, built once to avoid array-task races.
  ecoregions_by_code <- prepare_ecoregions_cached(
    ECO_PATH,
    force = FORCE_RERUN
  )
  ecoregion_cache_path <- file.path(CACHE_DIR, "ecoregions_by_code.rds")
  
  autocor_signature <- md5_text(c(
    "autocor_cache_v1",
    file_signature(ENV_PATH),
    file_signature(ECO_PATH),
    paste(valid_layer_names, collapse = "|"),
    CA_ALBERS,
    AUTOCOR_DOMAIN_BUFFER_M,
    paste0("terra=", package_versions[["terra"]])
  ))
  
  AUTOCOR_CACHE_PATH <- prepare_autocor_cache(
    env_stack = env_stack_master,
    ecoregions_by_code = ecoregions_by_code,
    cache_signature = autocor_signature,
    force = FORCE_RERUN
  )
  
  autocor_check <- terra::rast(AUTOCOR_CACHE_PATH)
  if (terra::nlyr(autocor_check) < 1L) {
    stop("Projected autocorrelation cache contains no layers.")
  }
  rm(autocor_check)
  gc(verbose = FALSE)
  
  # Predictor preflight and auditable mapping.
  predictor_lookup <- list()
  predictor_mapping_list <- list()
  preflight_errors <- list()
  
  for (i in seq_len(nrow(full_manifest))) {
    row <- full_manifest[i]
    key <- paste(row$taxon, row$species, sep = "||")
    
    resolved <- tryCatch(
      resolve_predictors(row$selected_predictors, valid_layer_names),
      error = function(e) {
        preflight_errors[[length(preflight_errors) + 1L]] <<- data.table(
          species = row$species,
          taxon = row$taxon,
          model_type = row$model_type,
          stage = "predictor preflight",
          error = conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(resolved)) {
      predictor_lookup[[key]] <- resolved
      mapping <- copy(resolved$mapping)
      mapping[, `:=`(
        species = row$species,
        taxon = row$taxon,
        model_type = row$model_type
      )]
      setcolorder(mapping, names(predictor_mapping_prototype))
      predictor_mapping_list[[length(predictor_mapping_list) + 1L]] <- mapping
    }
  }
  
  predictor_mapping_table <- if (length(predictor_mapping_list)) {
    rbindlist(predictor_mapping_list, use.names = TRUE, fill = TRUE)
  } else {
    predictor_mapping_prototype
  }
  
  preflight_error_table <- if (length(preflight_errors)) {
    rbindlist(preflight_errors, use.names = TRUE, fill = TRUE)
  } else {
    errors_prototype
  }
  
  atomic_fwrite(
    predictor_mapping_table,
    file.path(ARRAY_DIR, "predictor_mapping.csv")
  )
  atomic_fwrite(
    preflight_error_table,
    file.path(OUTPUT_ROOT, "sdm_spatial_cv_preflight_errors.csv")
  )
  
  # Preflight failures receive normal failed checkpoints and are excluded only
  # from runnable array tasks; they remain in full_manifest and final coverage.
  if (nrow(preflight_error_table)) {
    failed_rows <- merge(
      full_manifest,
      preflight_error_table,
      by = c("species", "taxon", "model_type"),
      all = FALSE,
      sort = FALSE
    )
    
    for (i in seq_len(nrow(failed_rows))) {
      row <- failed_rows[i]
      write_species_outputs(
        row = row,
        run_signature = RUN_SIGNATURE,
        state = "failed",
        stage = "predictor preflight",
        elapsed_seconds = 0,
        errors = row[, .(species, taxon, model_type, stage, error)]
      )
    }
  }
  
  runnable_manifest <- full_manifest[
    !paste(taxon, species, sep = "||") %chin%
      paste(
        preflight_error_table$taxon,
        preflight_error_table$species,
        sep = "||"
      )
  ]
  
  if (!nrow(runnable_manifest)) {
    stop("Every selected species failed predictor preflight.")
  }
  
  # Build one complete labeled occurrence cache per target taxon. These caches
  # intentionally retain records from species that will not themselves be
  # evaluated, because those records are valid target-group sampling locations.
  required_occurrence_columns <- c(
    "taxon", "species", "species_full", "longitude", "latitude"
  )
  
  pool_signature_base <- md5_text(c(
    "target_group_labeled_occurrence_pool_v1",
    file_signature(OCC_PATH),
    file_signature(ENV_PATH),
    paste(
      nrow(env_stack_master),
      ncol(env_stack_master),
      terra::res(env_stack_master),
      terra::ext(env_stack_master),
      terra::crs(env_stack_master),
      sep = "|"
    )
  ))
  
  taxon_names <- unique(full_manifest$taxon)
  background_pool_paths <- setNames(character(length(taxon_names)), taxon_names)
  pool_signatures <- setNames(
    vapply(
      taxon_names,
      function(taxon_name) md5_text(c(pool_signature_base, taxon_name)),
      character(1)
    ),
    taxon_names
  )
  
  cache_current <- vapply(
    taxon_names,
    function(taxon_name) {
      !FORCE_RERUN &&
        background_pool_is_current(taxon_name, pool_signatures[[taxon_name]])
    },
    logical(1)
  )
  
  occ_all <- NULL
  if (any(!cache_current)) {
    message(
      "Reading occurrence CSV once to build ",
      sum(!cache_current),
      " missing/stale taxon cache(s)..."
    )
    occ_all <- fread(
      OCC_PATH,
      select = required_occurrence_columns,
      showProgress = TRUE
    )
    
    missing_columns <- setdiff(required_occurrence_columns, names(occ_all))
    if (length(missing_columns)) {
      stop("Occurrence file is missing: ", paste(missing_columns, collapse = ", "))
    }
    
    stale_taxa <- taxon_names[!cache_current]
    occ_all <- occ_all[
      taxon %chin% stale_taxa &
        is.finite(longitude) &
        is.finite(latitude) &
        !is.na(species) &
        nzchar(species)
    ]
  }
  
  for (taxon_name in taxon_names) {
    background_pool_paths[[taxon_name]] <- prepare_background_pool(
      taxon = taxon_name,
      occ = occ_all,
      env_template = env_stack_master[[1]],
      pool_signature = pool_signatures[[taxon_name]],
      force = FORCE_RERUN
    )
    gc(verbose = FALSE)
  }
  
  if (!is.null(occ_all)) rm(occ_all)
  
  # Estimate thinned occurrence counts and relative computational costs. These
  # estimates are used only for load balancing; they never affect model fitting.
  estimate_taxon_counts <- function(pool, target_names) {
    a <- pool[
      species %chin% target_names,
      .(target = species, cell)
    ]
    b <- pool[
      !is.na(species_full) & species_full %chin% target_names,
      .(target = species_full, cell)
    ]
    
    matched <- unique(
      rbindlist(list(a, b), use.names = TRUE),
      by = c("target", "cell")
    )
    
    counts <- matched[, .(n_occ_est = .N), by = target]
    out <- data.table(species = target_names)
    out <- merge(
      out,
      counts,
      by.x = "species",
      by.y = "target",
      all.x = TRUE,
      sort = FALSE
    )
    out[is.na(n_occ_est), n_occ_est := 0L]
    out
  }
  
  cost_tables <- list()
  
  for (taxon_name in unique(runnable_manifest$taxon)) {
    pool <- load_background_pool(background_pool_paths[[taxon_name]])
    targets <- runnable_manifest[taxon == taxon_name, species]
    counts <- estimate_taxon_counts(pool, targets)
    
    x <- merge(
      runnable_manifest[taxon == taxon_name],
      counts,
      by = "species",
      all.x = TRUE,
      sort = FALSE
    )
    
    x[, predictor_count := lengths(strsplit(
      selected_predictors %||% "",
      ",",
      fixed = TRUE
    ))]
    x[predictor_count < 2L, predictor_count := 2L]
    
    x[, estimated_folds := pmin(
      MAX_FOLDS,
      pmax(MIN_FOLDS, floor(pmax(n_occ_est, 7L) / 3L))
    )]
    
    # ESM cost grows roughly with the number of predictor pairs; SDM cost grows
    # with folds, algorithms, rows, and predictor count. Only relative scale
    # matters for greedy balancing.
    x[, estimated_weight := ifelse(
      model_type == "ESM",
      estimated_folds *
        3 *
        pmax(1, choose(predictor_count, 2)) *
        (N_BG_ESM + pmin(n_occ_est, 1000L)) / N_BG_ESM,
      estimated_folds *
        4 *
        (1 + predictor_count / 10) *
        (N_BG_SDM + pmin(n_occ_est, 5000L)) / N_BG_SDM
    )]
    
    cost_tables[[length(cost_tables) + 1L]] <- x
    rm(pool)
    gc(verbose = FALSE)
  }
  
  runnable_manifest <- rbindlist(cost_tables, use.names = TRUE, fill = TRUE)
  
  assign_balanced_chunks <- function(x, n_chunks) {
    n_chunks <- min(as.integer(n_chunks), nrow(x))
    if (n_chunks < 1L) return(x[0])
    
    setorder(x, -estimated_weight, species)
    loads <- numeric(n_chunks)
    assignments <- integer(nrow(x))
    
    for (i in seq_len(nrow(x))) {
      chunk <- which.min(loads)
      assignments[i] <- chunk
      loads[chunk] <- loads[chunk] + x$estimated_weight[i]
    }
    
    x[, taxon_chunk := assignments]
    x
  }
  
  # Allocate the fixed task budget across taxa in proportion to estimated workload.
  # Every task still contains exactly one taxon.
  taxon_load <- runnable_manifest[
    ,
    .(
      n_species = .N,
      total_weight = sum(estimated_weight)
    ),
    by = taxon
  ]
  taxon_load[, taxon_order := match(taxon, TARGET_TAXA)]
  setorder(taxon_load, taxon_order)
  
  TASK_BUDGET <- if (RUN_MODE == "prepare_test") {
    nrow(taxon_load)  # one task per represented taxon for a local smoke test
  } else {
    min(CLUSTER_ARRAY_TASKS, nrow(runnable_manifest))
  }
  
  taxon_load[, n_chunks := 1L]
  remaining <- TASK_BUDGET - sum(taxon_load$n_chunks)
  
  while (remaining > 0L) {
    eligible <- which(taxon_load$n_chunks < taxon_load$n_species)
    if (!length(eligible)) break
    
    # Split the taxon whose current expected per-task workload is largest.
    score <- taxon_load$total_weight[eligible] / taxon_load$n_chunks[eligible]
    j <- eligible[which.max(score)]
    taxon_load$n_chunks[j] <- taxon_load$n_chunks[j] + 1L
    remaining <- remaining - 1L
  }
  
  split_tables <- lapply(
    seq_len(nrow(taxon_load)),
    function(i) {
      taxon_name <- taxon_load$taxon[i]
      assign_balanced_chunks(
        copy(runnable_manifest[taxon == taxon_name]),
        taxon_load$n_chunks[i]
      )
    }
  )
  
  species_tasks <- rbindlist(split_tables, use.names = TRUE, fill = TRUE)
  species_tasks[, taxon_order := match(taxon, TARGET_TAXA)]
  
  # Interleave taxa across task IDs to spread taxon-cache I/O.
  task_keys <- unique(
    species_tasks[, .(taxon, taxon_order, taxon_chunk)]
  )
  setorder(task_keys, taxon_chunk, taxon_order)
  task_keys[, task_id := seq_len(.N)]
  task_keys[, taxon_order := NULL]
  
  species_tasks <- merge(
    species_tasks,
    task_keys,
    by = c("taxon", "taxon_chunk"),
    all.x = TRUE,
    sort = FALSE
  )
  setorder(species_tasks, task_id, -estimated_weight, species)
  species_tasks[, taxon_order := NULL]
  
  task_manifest <- species_tasks[
    ,
    .(
      n_species = .N,
      n_esm = sum(model_type == "ESM"),
      n_sdm = sum(model_type == "SDM"),
      estimated_weight = sum(estimated_weight),
      minimum_species_weight = min(estimated_weight),
      maximum_species_weight = max(estimated_weight)
    ),
    by = .(task_id, taxon, taxon_chunk)
  ]
  setorder(task_manifest, task_id)
  
  full_manifest_path <- file.path(ARRAY_DIR, "full_manifest.csv")
  species_tasks_path <- file.path(ARRAY_DIR, "species_tasks.csv")
  task_manifest_path <- file.path(ARRAY_DIR, "task_manifest.csv")
  run_config_path <- file.path(ARRAY_DIR, "run_config.rds")
  
  atomic_fwrite(full_manifest, full_manifest_path)
  atomic_fwrite(species_tasks, species_tasks_path)
  atomic_fwrite(task_manifest, task_manifest_path)
  writeLines(as.character(nrow(task_manifest)), file.path(ARRAY_DIR, "task_count.txt"))
  
  run_config <- list(
    script_version = SCRIPT_VERSION,
    run_signature = RUN_SIGNATURE,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    root = ROOT,
    output_root = OUTPUT_ROOT,
    occurrence_path = OCC_PATH,
    environmental_path = ENV_PATH,
    ecoregion_path = ECO_PATH,
    manifest_dir = MANIFEST_DIR,
    full_manifest_path = full_manifest_path,
    species_tasks_path = species_tasks_path,
    task_manifest_path = task_manifest_path,
    predictor_mapping_path = file.path(ARRAY_DIR, "predictor_mapping.csv"),
    ecoregion_cache_path = ecoregion_cache_path,
    autocor_cache_path = AUTOCOR_CACHE_PATH,
    valid_layer_names = valid_layer_names,
    background_pool_paths = background_pool_paths,
    target_taxa = TARGET_TAXA,
    cluster_array_slots = CLUSTER_ARRAY_TASKS,
    test_n_per_taxon = TEST_N_PER_TAXON,
    n_tasks = nrow(task_manifest),
    n_species_expected = nrow(full_manifest),
    n_species_runnable = nrow(species_tasks),
    package_versions = package_versions
  )
  atomic_save_rds(run_config, run_config_path, compress = FALSE)
  
  message("\nArray preparation complete.")
  message("Expected species: ", nrow(full_manifest))
  message("Runnable species: ", nrow(species_tasks))
  message("Array tasks: ", nrow(task_manifest))
  message("Task table: ", task_manifest_path)
  message("Run configuration: ", run_config_path)
  
} else if (RUN_MODE == "worker") {
  
  task_id <- TASK_ID
  
  config_path <- file.path(ARRAY_DIR, "run_config.rds")
  if (!file.exists(config_path)) {
    stop(
      "Prepared array configuration not found: ", config_path,
      "\nRun `Rscript sdm_spatial_cv.R prep` first."
    )
  }
  
  config <- readRDS(config_path)
  
  if (!identical(config$script_version, SCRIPT_VERSION)) {
    stop(
      "Prepared configuration uses script version ", config$script_version,
      " but this worker uses ", SCRIPT_VERSION, ". Rerun preparation."
    )
  }
  
  if (!identical(
    as.character(config$package_versions[names(package_versions)]),
    as.character(package_versions)
  )) {
    stop(
      "Worker package versions differ from the preparation job. Use the same ",
      "R module and package library for preparation and array tasks."
    )
  }
  
  RUN_SIGNATURE <- config$run_signature
  ENV_PATH <- config$environmental_path
  ECO_PATH <- config$ecoregion_path
  AUTOCOR_CACHE_PATH <- config$autocor_cache_path
  valid_layer_names <- config$valid_layer_names
  
  full_manifest <- fread(config$full_manifest_path, showProgress = FALSE)
  species_tasks <- fread(config$species_tasks_path, showProgress = FALSE)
  task_manifest <- fread(config$task_manifest_path, showProgress = FALSE)
  
  this_task_id <- as.integer(task_id)
  task_row <- task_manifest[task_id == this_task_id]
  manifest <- species_tasks[task_id == this_task_id]
  
  if (!nrow(manifest)) {
    message("Task ", task_id, " has no assigned species; exiting successfully.")
    quit(save = "no", status = 0)
  }
  
  if (nrow(task_row) != 1L) {
    stop("Task ID ", task_id, " was not represented exactly once in task_manifest.csv.")
  }
  
  taxon_name <- unique(manifest$taxon)
  if (length(taxon_name) != 1L) {
    stop("An array task must contain exactly one taxon.")
  }
  
  # Node-local scratch substantially reduces shared-filesystem traffic from
  # BIOMOD2/ecospat model objects and terra temporary files.
  scratch_root <- Sys.getenv("SLURM_TMPDIR", unset = "")
  if (!nzchar(scratch_root)) {
    job_tag <- paste0(
      Sys.getenv("SLURM_ARRAY_JOB_ID", unset = "local"),
      "_",
      task_id
    )
    scratch_root <- file.path(OUTPUT_ROOT, "slurm_tmp", job_tag)
  }
  
  WORK_ROOT <- file.path(scratch_root, "work")
  TEMP_ROOT <- file.path(scratch_root, "terra_tmp")
  dir.create(WORK_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
  
  terra::terraOptions(
    memfrac = TERRA_MEMFRAC,
    tempdir = TEMP_ROOT,
    progress = 0,
    parallel = FALSE,
    threads = TERRA_THREADS
  )
  
  on.exit({
    try(
      terra::tmpFiles(current = TRUE, orphan = TRUE, remove = TRUE),
      silent = TRUE
    )
    # SLURM_TMPDIR is automatically removed by many clusters; the fallback
    # directory is ours and can be safely removed here.
    if (!nzchar(Sys.getenv("SLURM_TMPDIR", unset = ""))) {
      unlink(scratch_root, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  
  ecoregions_by_code <- readRDS(config$ecoregion_cache_path)
  background_pool_current <- load_background_pool(
    config$background_pool_paths[[taxon_name]]
  )
  
  # The full taxon cache supplies both focal occurrences and target-group
  # background candidates, including occurrences from non-modeled species.
  occ <- background_pool_current
  
  mapping_table <- fread(config$predictor_mapping_path, showProgress = FALSE)
  mapping_table <- mapping_table[
    paste(taxon, species, sep = "||") %chin%
      paste(manifest$taxon, manifest$species, sep = "||")
  ]
  
  predictor_lookup <- lapply(
    split(
      mapping_table,
      paste(mapping_table$taxon, mapping_table$species, sep = "||")
    ),
    function(x) {
      list(
        raster_names = x$raster_layer,
        model_names = x$requested_predictor,
        mapping = x[
          ,
          .(requested_predictor, raster_layer, match_type)
        ]
      )
    }
  )
  
  task_status_dir <- file.path(ARRAY_DIR, "task_status")
  dir.create(task_status_dir, recursive = TRUE, showWarnings = FALSE)
  task_status_path <- file.path(
    task_status_dir,
    sprintf("task_%04d.csv", task_id)
  )
  
  should_run <- vapply(
    seq_len(nrow(manifest)),
    function(i) {
      row <- manifest[i]
      
      if (FORCE_RERUN) return(TRUE)
      if (RESUME_COMPLETED &&
          status_is_current_complete(row, RUN_SIGNATURE)) {
        return(FALSE)
      }
      if (!RETRY_FAILED &&
          status_is_current_failed(row, RUN_SIGNATURE)) {
        return(FALSE)
      }
      TRUE
    },
    logical(1)
  )
  
  pending <- which(should_run)
  message(
    "Array task ", task_id,
    ": taxon=", taxon_name,
    "; assigned=", nrow(manifest),
    "; pending=", length(pending),
    "; skipped current checkpoints=", sum(!should_run)
  )
  
  status_rows <- list()
  
  # Record already-current checkpoints in task status.
  if (any(!should_run)) {
    for (i in which(!should_run)) {
      existing <- read_status(manifest[i])
      if (!is.null(existing)) {
        status_rows[[length(status_rows) + 1L]] <- status_to_row(existing)
      }
    }
  }
  
  write_task_status <- function() {
    x <- if (length(status_rows)) {
      rbindlist(status_rows, use.names = TRUE, fill = TRUE)
    } else {
      status_prototype
    }
    x[, `:=`(
      task_id = task_id,
      task_taxon = taxon_name
    )]
    atomic_fwrite(x, task_status_path)
  }
  
  write_task_status()
  
  if (length(pending)) {
    for (position in seq_along(pending)) {
      row_index <- pending[position]
      row <- manifest[row_index]
      
      message(
        sprintf(
          "[task %d] species %d/%d: %s (%s, %s)",
          task_id,
          position,
          length(pending),
          row$species,
          row$taxon,
          row$model_type
        )
      )
      
      status <- process_species_safe(row_index)
      status_rows[[length(status_rows) + 1L]] <- status_to_row(status)
      write_task_status()
      
      message(
        "[task ", task_id, "] ",
        row$species, ": ", status$state,
        " at stage ", status$stage
      )
    }
  }
  
  final_status <- fread(task_status_path, showProgress = FALSE)
  message(
    "Array task ", task_id, " finished: ",
    sum(final_status$state == "complete"), " complete, ",
    sum(final_status$state == "failed"), " failed."
  )
  
  # Species-level scientific failures are checkpointed and reported but do not
  # make the Slurm task itself fail. Unhandled worker/script failures still exit
  # nonzero through stop(), which Slurm records correctly.
  
} else if (RUN_MODE == "compile") {
  
  config_path <- file.path(ARRAY_DIR, "run_config.rds")
  if (!file.exists(config_path)) {
    stop("Prepared array configuration not found: ", config_path)
  }
  
  config <- readRDS(config_path)
  if (!identical(config$script_version, SCRIPT_VERSION)) {
    stop(
      "Prepared configuration uses script version ", config$script_version,
      " but this compiler uses ", SCRIPT_VERSION, "."
    )
  }
  
  if (!identical(
    as.character(config$package_versions[names(package_versions)]),
    as.character(package_versions)
  )) {
    stop(
      "Compiler package versions differ from the preparation job. Use the same ",
      "R module and package library throughout the pipeline."
    )
  }
  
  RUN_SIGNATURE <- config$run_signature
  ENV_PATH <- config$environmental_path
  ECO_PATH <- config$ecoregion_path
  AUTOCOR_CACHE_PATH <- config$autocor_cache_path
  valid_layer_names <- config$valid_layer_names
  full_manifest <- fread(config$full_manifest_path, showProgress = FALSE)
  
  message(
    "Compiling current checkpoints for ",
    nrow(full_manifest),
    " expected ESM/SDM species..."
  )
  
  compile_results <- list()
  
  for (taxon_name in unique(full_manifest$taxon)) {
    compile_results[[taxon_name]] <- compile_scope_outputs(
      rows = full_manifest[taxon == taxon_name],
      scope_name = taxon_name
    )
  }
  
  old_write_combined <- WRITE_COMBINED_FOLD_ASSIGNMENTS
  WRITE_COMBINED_FOLD_ASSIGNMENTS <- WRITE_GLOBAL_FOLD_ASSIGNMENTS
  
  compile_results[["all_taxa"]] <- compile_scope_outputs(
    rows = full_manifest,
    scope_name = "all_taxa"
  )
  
  WRITE_COMBINED_FOLD_ASSIGNMENTS <- old_write_combined
  
  # Consolidate task-level status files for operational auditing.
  task_status_files <- list.files(
    file.path(ARRAY_DIR, "task_status"),
    pattern = "^task_[0-9]+\\.csv$",
    full.names = TRUE
  )
  
  task_status <- if (length(task_status_files)) {
    rbindlist(
      lapply(task_status_files, fread, showProgress = FALSE),
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table(
      task_id = integer(),
      task_taxon = character()
    )
  }
  
  atomic_fwrite(
    task_status,
    file.path(ARRAY_DIR, "array_task_species_status.csv")
  )
  
  coverage <- fread(
    file.path(
      COMBINED_DIR,
      "all_taxa",
      "sdm_spatial_cv_evaluation_coverage.csv"
    ),
    showProgress = FALSE
  )
  
  message("\nCompilation complete.")
  print(coverage)
  message(
    "Review: ",
    file.path(COMBINED_DIR, "all_taxa", "sdm_spatial_cv_errors.csv")
  )
  
}
