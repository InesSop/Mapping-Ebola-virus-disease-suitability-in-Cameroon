# ============================================================
# Script: 02_wlc_experts.R
#
# Purpose:
#   Build dry- and wet-season Ebola virus disease (EVD) risk
#   maps using a weighted linear combination multi-criteria
#   decision analysis (WLC-MCDA) framework.
#
#   Risk factors and their relative weights were derived from
#   field-expert elicitation. Separate weighting schemes are
#   applied to the dry and wet seasons to reflect seasonal
#   differences in the perceived importance of EVD risk factors.
#
# Inputs:
#   Preprocessed and standardised raster layers representing
#   EVD risk factors for the dry and wet seasons.
#
# Outputs:
#   - Expert-based EVD risk raster for the dry season
#   - Expert-based EVD risk raster for the wet season
#
# Notes:
#   - All input layers are assumed to have been previously
#     standardised to a common 0-1 scale.
#   - Upstream geoprocessing and standardisation of the original
#     spatial datasets are not performed in this script.
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
  "experts",
  "dry"
)

input_dir_wet <- here(
  "inputs",
  "standardised_rasters",
  "experts",
  "wet"
)

output_dir <- here("outputs", "wlc_experts")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 3. Define expert-derived seasonal weights
# ------------------------------------------------------------

# Dry-season expert-derived weights
weights_expert_dry <- c(
  rainfall                 = 0.007833333,
  bushmeat_hunting         = 0.307000000,
  bushmeat_trade           = 0.099777778,
  deforestation            = 0.106166667,
  forests                  = 0.059500000,
  NHP_density              = 0.051444444,
  bat_density              = 0.033388889,
  poverty                  = 0.028111111,
  human_population_density = 0.306777800
)

# Wet-season expert-derived weights
weights_expert_wet <- c(
  rainfall                 = 0.025541667,
  bushmeat_hunting         = 0.301687500,
  bushmeat_trade           = 0.098458333,
  deforestation            = 0.088666667,
  forests                  = 0.067625000,
  NHP_density              = 0.051375000,
  bat_density              = 0.028250000,
  poverty                  = 0.024916667,
  human_population_density = 0.313479200
)


# ------------------------------------------------------------
# 4. Check that weights sum to 1
# ------------------------------------------------------------

if (!isTRUE(all.equal(
  sum(weights_expert_dry),
  1,
  tolerance = 1e-6
))) {
  stop("Dry-season expert-derived weights do not sum to 1.")
}

if (!isTRUE(all.equal(
  sum(weights_expert_wet),
  1,
  tolerance = 1e-6
))) {
  stop("Wet-season expert-derived weights do not sum to 1.")
}


# ------------------------------------------------------------
# 5. Helper function for raster alignment
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
# 6. Load dry-season reference raster
# ------------------------------------------------------------

rainfall_dry <- rast(
  here(
    input_dir_dry,
    "Standardized_dry_season_BIO13_precipitation_wettest_month.tif"
  )
)

if (is.na(crs(rainfall_dry))) {
  stop("The dry-season reference raster has no defined CRS.")
}


# ------------------------------------------------------------
# 7. Load and align dry-season risk-factor rasters
# ------------------------------------------------------------

bushmeat_hunting_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "Relative_hunting_pressure_index_1km_Cameroon.tif"
    )
  ),
  rainfall_dry
)

bushmeat_trade_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "Standardized_CMR_bushmeat_trade_distance_based_0_1_1km_UTM33N.tif"
    )
  ),
  rainfall_dry
)

deforestation_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "DeforestationRisk_1km.tif"
    )
  ),
  rainfall_dry
)

forests_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "Standardized_prop_foret_1km_2020_CMR.tif"
    )
  ),
  rainfall_dry
)

NHP_density_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "great_ape_suitability_APES_1km_UTM33N.tif"
    )
  ),
  rainfall_dry
)

bat_density_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
    )
  ),
  rainfall_dry
)

poverty_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "Risk_index_RWI_1km_UTM33N_CMR.tif"
    )
  ),
  rainfall_dry
)

human_population_density_dry <- align_to_reference(
  rast(
    here(
      input_dir_dry,
      "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
    )
  ),
  rainfall_dry
)


# ------------------------------------------------------------
# 8. Verify dry-season raster geometry
# ------------------------------------------------------------

dry_alignment_check <- all(
  compareGeom(rainfall_dry, bushmeat_hunting_dry),
  compareGeom(rainfall_dry, bushmeat_trade_dry),
  compareGeom(rainfall_dry, deforestation_dry),
  compareGeom(rainfall_dry, forests_dry),
  compareGeom(rainfall_dry, NHP_density_dry),
  compareGeom(rainfall_dry, bat_density_dry),
  compareGeom(rainfall_dry, poverty_dry),
  compareGeom(rainfall_dry, human_population_density_dry)
)

if (!dry_alignment_check) {
  stop("Dry-season expert-model rasters are not correctly aligned.")
}

message("Dry-season expert-model rasters successfully aligned.")


# ------------------------------------------------------------
# 9. Calculate expert-based dry-season EVD risk index
# ------------------------------------------------------------

risk_index_expert_dry <-
  rainfall_dry *
  weights_expert_dry["rainfall"] +
  bushmeat_hunting_dry *
  weights_expert_dry["bushmeat_hunting"] +
  bushmeat_trade_dry *
  weights_expert_dry["bushmeat_trade"] +
  deforestation_dry *
  weights_expert_dry["deforestation"] +
  forests_dry *
  weights_expert_dry["forests"] +
  NHP_density_dry *
  weights_expert_dry["NHP_density"] +
  bat_density_dry *
  weights_expert_dry["bat_density"] +
  poverty_dry *
  weights_expert_dry["poverty"] +
  human_population_density_dry *
  weights_expert_dry["human_population_density"]


# ------------------------------------------------------------
# 10. Export dry-season expert-based risk raster
# ------------------------------------------------------------

writeRaster(
  risk_index_expert_dry,
  here(
    output_dir,
    "risk_index_EVD_expert_dry.tif"
  ),
  overwrite = TRUE
)

message("Dry-season expert-based EVD risk map exported.")


# ============================================================
# WET SEASON
# ============================================================


# ------------------------------------------------------------
# 11. Load wet-season reference raster
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
# 12. Load and align wet-season risk-factor rasters
# ------------------------------------------------------------

bushmeat_hunting_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "Relative_hunting_pressure_index_1km_Cameroon.tif"
    )
  ),
  rainfall_wet
)

bushmeat_trade_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "Standardized_CMR_bushmeat_trade_distance_based_0_1_1km_UTM33N.tif"
    )
  ),
  rainfall_wet
)

deforestation_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "DeforestationRisk_1km.tif"
    )
  ),
  rainfall_wet
)

forests_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "Standardized_prop_foret_1km_2020_CMR.tif"
    )
  ),
  rainfall_wet
)

NHP_density_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "great_ape_suitability_APES_1km_UTM33N.tif"
    )
  ),
  rainfall_wet
)

bat_density_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
    )
  ),
  rainfall_wet
)

poverty_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "Risk_index_RWI_1km_UTM33N_CMR.tif"
    )
  ),
  rainfall_wet
)

human_population_density_wet <- align_to_reference(
  rast(
    here(
      input_dir_wet,
      "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
    )
  ),
  rainfall_wet
)


# ------------------------------------------------------------
# 13. Verify wet-season raster geometry
# ------------------------------------------------------------

wet_alignment_check <- all(
  compareGeom(rainfall_wet, bushmeat_hunting_wet),
  compareGeom(rainfall_wet, bushmeat_trade_wet),
  compareGeom(rainfall_wet, deforestation_wet),
  compareGeom(rainfall_wet, forests_wet),
  compareGeom(rainfall_wet, NHP_density_wet),
  compareGeom(rainfall_wet, bat_density_wet),
  compareGeom(rainfall_wet, poverty_wet),
  compareGeom(rainfall_wet, human_population_density_wet)
)

if (!wet_alignment_check) {
  stop("Wet-season expert-model rasters are not correctly aligned.")
}

message("Wet-season expert-model rasters successfully aligned.")


# ------------------------------------------------------------
# 14. Calculate expert-based wet-season EVD risk index
# ------------------------------------------------------------

risk_index_expert_wet <-
  rainfall_wet *
  weights_expert_wet["rainfall"] +
  bushmeat_hunting_wet *
  weights_expert_wet["bushmeat_hunting"] +
  bushmeat_trade_wet *
  weights_expert_wet["bushmeat_trade"] +
  deforestation_wet *
  weights_expert_wet["deforestation"] +
  forests_wet *
  weights_expert_wet["forests"] +
  NHP_density_wet *
  weights_expert_wet["NHP_density"] +
  bat_density_wet *
  weights_expert_wet["bat_density"] +
  poverty_wet *
  weights_expert_wet["poverty"] +
  human_population_density_wet *
  weights_expert_wet["human_population_density"]


# ------------------------------------------------------------
# 15. Export wet-season expert-based risk raster
# ------------------------------------------------------------

writeRaster(
  risk_index_expert_wet,
  here(
    output_dir,
    "risk_index_EVD_expert_wet.tif"
  ),
  overwrite = TRUE
)

message("Wet-season expert-based EVD risk map exported.")


# ------------------------------------------------------------
# 16. Completion message
# ------------------------------------------------------------

message(
  "Expert-based WLC-MCDA modelling completed successfully ",
  "for both dry and wet seasons."
)