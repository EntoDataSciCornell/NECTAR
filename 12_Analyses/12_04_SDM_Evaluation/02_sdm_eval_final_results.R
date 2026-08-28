library(data.table)
library(ggplot2)

# ---- Paths ----

d <- "/local/storage/NECTAR/SDM_Eval/combined/all_taxa"

# ---- Read compiled evaluation outputs ----

m <- fread(file.path(d, "sdm_spatial_cv_species_metrics.csv"))
cov <- fread(file.path(d, "sdm_spatial_cv_evaluation_coverage.csv"))

taxa <- c("plants", "bees", "butterflies", "hoverflies", "moths")
labs <- c(
  plants = "Plants",
  bees = "Bees",
  butterflies = "Butterflies",
  hoverflies = "Hoverflies",
  moths = "Moths"
)

# ---- Species-level evaluation metrics ----

m <- m[
  taxon %chin% taxa &
    metric %chin% c("BOYCE", "AUC", "TSS")
]

s <- dcast(
  unique(m[, .(species, taxon, model_type, metric, mean)]),
  species + taxon + model_type ~ metric,
  value.var = "mean"
)

fmt <- function(x) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf("%.3f [%.3f-%.3f]", q[2], q[1], q[3])
}

pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)

# ---- Supplementary table ----

tab <- s[, .(
  `Boyce median [IQR]` = fmt(BOYCE),
  `Boyce >=0.25 (%)`   = pct(BOYCE >= 0.25),
  `AUC median [IQR]`   = fmt(AUC),
  `AUC >=0.70 (%)`     = pct(AUC >= 0.70),
  `TSS median [IQR]`   = fmt(TSS),
  `TSS >=0.40 (%)`     = pct(TSS >= 0.40)
), by = .(taxon, model_type)]

tab <- merge(
  cov[taxon %chin% taxa,
      .(taxon, model_type, N = n_complete)],
  tab,
  by = c("taxon", "model_type")
)

tab[, Taxon := labs[taxon]]
tab[, Model := model_type]

tab <- tab[
  order(match(taxon, taxa), match(model_type, c("ESM", "SDM"))),
  .(
    Taxon,
    Model,
    N,
    `Boyce median [IQR]`,
    `Boyce >=0.25 (%)`,
    `AUC median [IQR]`,
    `AUC >=0.70 (%)`,
    `TSS median [IQR]`,
    `TSS >=0.40 (%)`
  )
]

fwrite(tab, file.path(d, "Table_S1_SDM_evaluation.csv"))
print(tab)

# ---- Supplementary figure: Boyce distributions ----

s[, taxon := factor(taxon, levels = taxa, labels = labs[taxa])]
s[, model_type := factor(model_type, levels = c("ESM", "SDM"))]

p <- ggplot(s[is.finite(BOYCE)], aes(taxon, BOYCE)) +
  geom_boxplot(width = 0.65, outlier.size = 0.4) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0.25, linetype = 2) +
  facet_wrap(~model_type, nrow = 1) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, y = "Continuous Boyce index") +
  theme_classic(base_size = 11)

ggsave(
  file.path(d, "Figure_S1_Boyce_distribution.png"),
  p, width = 8, height = 4
)

#### List of all species with Boyce < 0 ####

m <- fread(file.path(d, "sdm_spatial_cv_species_metrics.csv"))

bad_boyce <- m[
  taxon %in% c("plants", "bees", "butterflies", "hoverflies", "moths") &
    metric == "BOYCE" &
    is.finite(mean) &
    mean <= 0.25,
  .(taxon, species)
]

setorder(bad_boyce, taxon, species)

fwrite(
  bad_boyce,
  file.path(d, "species_boyce_0.25_or_lower.csv")
)

#### List of all the ultra-rare species ####

manifest_files <- list.files(
  "/local/storage/NECTAR/Data_Clean/SDMs/sdm_by_taxon",
  pattern = "_continuous_manifest\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

ultra_rare <- rbindlist(lapply(manifest_files, fread), fill = TRUE)[
  model_type == "ultraRare_alphahull",
  .(taxon, species)
] |> unique()

setorder(ultra_rare, taxon, species)

fwrite(
  ultra_rare,
  file.path(d, "ultra_rare_species.csv")
)
