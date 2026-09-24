# ============================================================
# Script: 01_wlc_literature.R
#
# Purpose:
#   Build dry- and wet-season Ebola virus disease (EVD) risk
#   maps using a weighted linear combination multi-criteria
#   decision analysis (WLC-MCDA) framework.
#
#   Risk factors and their weights were derived from the
#   scientific literature reporting statistically significant
#   associations between candidate risk factors and EVD
#   occurrence.
#
# Inputs:
#   Preprocessed and standardised raster layers representing
#   EVD risk factors for the dry and wet seasons.
#
# Outputs:
#   - Literature-based EVD risk raster for the dry season
#   - Literature-based EVD risk raster for the wet season
#
# Notes:
#   - All input layers are assumed to have been previously
#     standardised to a common 0-1 scale.
#   - Raster preparation and upstream geoprocessing are not
#     performed in this script.
#   - All layers are aligned to a common reference raster before
#     implementation of the WLC model.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(here)


# ------------------------------------------------------------
# 2. Define input and output directories
# ------------------------------------------------------------

input_dir_dry <- here(
  "inputs",
  "standardised_rasters",
  "literature",
  "dry"
)

input_dir_wet <- here(
  "inputs",
  "standardised_rasters",
  "literature",
  "wet"
)

output_dir <- here("outputs", "wlc_literature")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 3. Define literature-derived weights
# ------------------------------------------------------------

# The same literature-derived weighting scheme is applied to
# both seasons. Seasonal variation is represented through the
# season-specific raster layers where applicable.

weights_literature <- c(
  rainfall                 = 0.090909091,
  bushmeat_hunting         = 0.068181818,
  temperature              = 0.227272727,
  forests                  = 0.227272727,
  NHP_density              = 0.045454546,
  bat_density              = 0.045454546,
  human_population_density = 0.204545500,
  elevation                = 0.090909091
)

# Check that weights sum approximately to 1
if (!isTRUE(all.equal(sum(weights_literature), 1, tolerance = 1e-6))) {
  stop("Literature-derived weights do not sum to 1.")
}


# ------------------------------------------------------------
# 4. Helper function for raster alignment
# ------------------------------------------------------------

align_to_reference <- function(raster, reference) {
  project(
    raster,
    reference,
    method = "bilinear"
  )
}


# ============================================================
# DRY SEASON
# ============================================================


# ------------------------------------------------------------
# 5. Load dry-season reference raster
# ------------------------------------------------------------

rainfall_dry <- rast(
  here(
    input_dir_dry,
    "Standardized_dry_season_BIO13_precipitation_wettest_month.tif"
  )
)

# Ensure that the reference raster has a defined CRS
if (is.na(crs(rainfall_dry))) {
  stop("The dry-season reference raster has no defined CRS.")
}


# ------------------------------------------------------------
# 6. Load and align dry-season risk-factor rasters
# ------------------------------------------------------------

bushmeat_hunting_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "Relative_hunting_pressure_index_1km_Cameroon.tif"
  )),
  rainfall_dry
)

temperature_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "MDTR_dry_utm33_1km_std01_minmax.tif"
  )),
  rainfall_dry
)

forests_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "Standardized_prop_foret_1km_2020_CMR.tif"
  )),
  rainfall_dry
)

NHP_density_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "great_ape_suitability_APES_1km_UTM33N.tif"
  )),
  rainfall_dry
)

bat_density_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
  )),
  rainfall_dry
)

human_population_density_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
  )),
  rainfall_dry
)

elevation_dry <- align_to_reference(
  rast(here(
    input_dir_dry,
    "Standardized_Elevation_low_altitude_index_1km_UTM33N_new.tif"
  )),
  rainfall_dry
)


# ------------------------------------------------------------
# 7. Verify raster geometry
# ------------------------------------------------------------

dry_alignment_check <- all(
  compareGeom(rainfall_dry, bushmeat_hunting_dry),
  compareGeom(rainfall_dry, temperature_dry),
  compareGeom(rainfall_dry, forests_dry),
  compareGeom(rainfall_dry, NHP_density_dry),
  compareGeom(rainfall_dry, bat_density_dry),
  compareGeom(rainfall_dry, human_population_density_dry),
  compareGeom(rainfall_dry, elevation_dry)
)

if (!dry_alignment_check) {
  stop("Dry-season rasters are not correctly aligned.")
}

message("Dry-season rasters successfully aligned.")


# ------------------------------------------------------------
# 8. Calculate literature-based dry-season EVD risk index
# ------------------------------------------------------------

risk_index_literature_dry <-
  rainfall_dry *
  weights_literature["rainfall"] +
  bushmeat_hunting_dry *
  weights_literature["bushmeat_hunting"] +
  temperature_dry *
  weights_literature["temperature"] +
  forests_dry *
  weights_literature["forests"] +
  NHP_density_dry *
  weights_literature["NHP_density"] +
  bat_density_dry *
  weights_literature["bat_density"] +
  human_population_density_dry *
  weights_literature["human_population_density"] +
  elevation_dry *
  weights_literature["elevation"]


# ------------------------------------------------------------
# 9. Export dry-season risk raster
# ------------------------------------------------------------

writeRaster(
  risk_index_literature_dry,
  here(
    output_dir,
    "risk_index_EVD_literature_dry.tif"
  ),
  overwrite = TRUE
)

message("Dry-season literature-based EVD risk map exported.")


# ============================================================
# WET SEASON
# ============================================================


# ------------------------------------------------------------
# 10. Load wet-season reference raster
# ------------------------------------------------------------

rainfall_wet <- rast(
  here(
    input_dir_wet,
    "Standardized_wet_season_BIO13_precipitation_wettest_month.tif"
  )
)

if (is.na(crs(rainfall_wet))) {
  stop("The wet-season reference raster has no defined CRS.")
}


# ------------------------------------------------------------
# 11. Load and align wet-season risk-factor rasters
# ------------------------------------------------------------

bushmeat_hunting_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "Relative_hunting_pressure_index_1km_Cameroon.tif"
  )),
  rainfall_wet
)

temperature_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "MDTR_wet_utm33_1km_std01_minmax.tif"
  )),
  rainfall_wet
)

forests_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "Standardized_prop_foret_1km_2020_CMR.tif"
  )),
  rainfall_wet
)

NHP_density_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "great_ape_suitability_APES_1km_UTM33N.tif"
  )),
  rainfall_wet
)

bat_density_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
  )),
  rainfall_wet
)

human_population_density_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
  )),
  rainfall_wet
)

elevation_wet <- align_to_reference(
  rast(here(
    input_dir_wet,
    "Standardized_Elevation_low_altitude_index_1km_UTM33N_new.tif"
  )),
  rainfall_wet
)


# ------------------------------------------------------------
# 12. Verify raster geometry
# ------------------------------------------------------------

wet_alignment_check <- all(
  compareGeom(rainfall_wet, bushmeat_hunting_wet),
  compareGeom(rainfall_wet, temperature_wet),
  compareGeom(rainfall_wet, forests_wet),
  compareGeom(rainfall_wet, NHP_density_wet),
  compareGeom(rainfall_wet, bat_density_wet),
  compareGeom(rainfall_wet, human_population_density_wet),
  compareGeom(rainfall_wet, elevation_wet)
)

if (!wet_alignment_check) {
  stop("Wet-season rasters are not correctly aligned.")
}

message("Wet-season rasters successfully aligned.")


# ------------------------------------------------------------
# 13. Calculate literature-based wet-season EVD risk index
# ------------------------------------------------------------

risk_index_literature_wet <-
  rainfall_wet *
  weights_literature["rainfall"] +
  bushmeat_hunting_wet *
  weights_literature["bushmeat_hunting"] +
  temperature_wet *
  weights_literature["temperature"] +
  forests_wet *
  weights_literature["forests"] +
  NHP_density_wet *
  weights_literature["NHP_density"] +
  bat_density_wet *
  weights_literature["bat_density"] +
  human_population_density_wet *
  weights_literature["human_population_density"] +
  elevation_wet *
  weights_literature["elevation"]


# ------------------------------------------------------------
# 14. Export wet-season risk raster
# ------------------------------------------------------------

writeRaster(
  risk_index_literature_wet,
  here(
    output_dir,
    "risk_index_EVD_literature_wet.tif"
  ),
  overwrite = TRUE
)

message("Wet-season literature-based EVD risk map exported.")


# ------------------------------------------------------------
# 15. Completion message
# ------------------------------------------------------------

message(
  "Literature-based WLC-MCDA modelling completed successfully ",
  "for both dry and wet seasons."
)