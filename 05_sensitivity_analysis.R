# ============================================================
# Script: 05_sensitivity_analysis.R
#
# Purpose:
#   Assess the sensitivity of the seasonal Ebola virus disease
#   (EVD) WLC-MCDA models to uncertainty in risk-factor weights.
#
#   For each risk factor, its baseline weight is varied from
#   -25% to +25% in 5% increments. Following each perturbation,
#   the remaining factor weights are proportionally rescaled so
#   that the total weight remains equal to 1.
#
#   The resulting perturbed risk surfaces are used to calculate:
#
#   1. A pixel-wise uncertainty surface, expressed as the
#      standard deviation across all perturbed risk maps.
#
#   2. An Average Relative Change (ARC) surface, describing the
#      mean relative deviation of perturbed maps from the
#      baseline model.
#
# Inputs:
#   Preprocessed and standardised risk-factor rasters used in the
#   seasonal combined-weight sensitivity models.
#
# Outputs:
#   - Dry-season baseline sensitivity model
#   - Wet-season baseline sensitivity model
#   - Dry-season uncertainty surface
#   - Wet-season uncertainty surface
#   - Dry-season Average Relative Change surface
#   - Wet-season Average Relative Change surface
#
# Notes:
#   - All input rasters are assumed to have previously been
#     standardised to a common 0-1 scale.
#   - The sensitivity weights used here correspond to the
#     seasonal combined-weight models employed for the sensitivity
#     analysis.
#   - The sensitivity model is distinct from the final annual
#     combined surface obtained by averaging the annual expert-
#     and literature-based models.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(here)


# ------------------------------------------------------------
# 2. Define directories
# ------------------------------------------------------------

input_dir_dry <- here(
  "inputs",
  "standardised_rasters",
  "sensitivity",
  "dry"
)

input_dir_wet <- here(
  "inputs",
  "standardised_rasters",
  "sensitivity",
  "wet"
)

output_dir <- here(
  "outputs",
  "sensitivity_analysis"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Define perturbation range
# ------------------------------------------------------------

# Weight perturbations from -25% to +25% in 5% increments

deltas <- seq(
  -0.25,
  0.25,
  by = 0.05
)


# ------------------------------------------------------------
# 4. Define seasonal baseline weights
# ------------------------------------------------------------

# Dry-season combined weights

weights_dry <- c(
  rainfall                 = 0.049371212,
  bushmeat_hunting         = 0.187590909,
  bushmeat_trade           = 0.049888889,
  deforestation            = 0.053083333,
  forests                  = 0.143386364,
  NHP_density              = 0.048449495,
  bat_density              = 0.039421717,
  poverty                  = 0.014055556,
  human_population_density = 0.255661600,
  temperature              = 0.113636364,
  elevation                = 0.045454545
)


# Wet-season combined weights

weights_wet <- c(
  rainfall                 = 0.058225379,
  bushmeat_hunting         = 0.184934659,
  bushmeat_trade           = 0.049229167,
  deforestation            = 0.044333333,
  forests                  = 0.147448864,
  NHP_density              = 0.048414773,
  bat_density              = 0.036852273,
  poverty                  = 0.012458333,
  human_population_density = 0.259012300,
  temperature              = 0.113636364,
  elevation                = 0.045454545
)


# ------------------------------------------------------------
# 5. Check weight sums
# ------------------------------------------------------------

if (!isTRUE(
  all.equal(
    sum(weights_dry),
    1,
    tolerance = 1e-6
  )
)) {
  stop(
    "Dry-season sensitivity weights do not sum to 1."
  )
}


if (!isTRUE(
  all.equal(
    sum(weights_wet),
    1,
    tolerance = 1e-6
  )
)) {
  stop(
    "Wet-season sensitivity weights do not sum to 1."
  )
}


# ------------------------------------------------------------
# 6. Helper function: raster alignment
# ------------------------------------------------------------

align_to_reference <- function(raster, reference) {
  
  project(
    raster,
    reference,
    method = "bilinear"
  )
  
}


# ------------------------------------------------------------
# 7. Helper function: geometry check
# ------------------------------------------------------------

check_geometry <- function(reference, raster_list) {
  
  checks <- vapply(
    raster_list,
    function(x) {
      
      compareGeom(
        reference,
        x,
        stopOnError = FALSE
      )
      
    },
    logical(1)
  )
  
  all(checks)
}


# ------------------------------------------------------------
# 8. Helper function: calculate WLC risk surface
# ------------------------------------------------------------

calculate_wlc <- function(rasters, weights) {
  
  if (!all(names(weights) %in% names(rasters))) {
    stop(
      "Not all weighted factors are available in the raster list."
    )
  }
  
  risk <-
    rasters[[names(weights)[1]]] *
    weights[1]
  
  for (i in 2:length(weights)) {
    
    factor_name <- names(weights)[i]
    
    risk <-
      risk +
      rasters[[factor_name]] *
      weights[i]
    
  }
  
  return(risk)
}


# ------------------------------------------------------------
# 9. Helper function: weight perturbation analysis
# ------------------------------------------------------------

run_weight_perturbation <- function(
    rasters,
    baseline_weights,
    deltas
) {
  
  factors <- names(baseline_weights)
  
  risk_maps <- list()
  
  for (factor_name in factors) {
    
    for (delta in deltas) {
      
      # Perturb the selected factor weight
      new_weight <-
        baseline_weights[factor_name] *
        (1 + delta)
      
      # Identify all remaining factors
      other_factors <-
        setdiff(
          factors,
          factor_name
        )
      
      # Remaining weight to be distributed
      remaining_weight <-
        1 - new_weight
      
      # Preserve the proportional contribution of the
      # non-perturbed factors
      proportional_weights <-
        baseline_weights[other_factors] /
        sum(
          baseline_weights[other_factors]
        )
      
      adjusted_weights <-
        setNames(
          numeric(
            length(baseline_weights)
          ),
          factors
        )
      
      adjusted_weights[factor_name] <-
        new_weight
      
      adjusted_weights[other_factors] <-
        proportional_weights *
        remaining_weight
      
      # Verify that adjusted weights sum to 1
      if (
        abs(
          sum(adjusted_weights) - 1
        ) > 1e-8
      ) {
        
        stop(
          paste(
            "Adjusted weights do not sum to 1 for",
            factor_name,
            "with delta",
            delta
          )
        )
      }
      
      # Recalculate risk surface
      perturbed_risk <-
        calculate_wlc(
          rasters,
          adjusted_weights
        )
      
      map_name <- paste0(
        factor_name,
        "_",
        sprintf("%+.2f", delta)
      )
      
      risk_maps[[map_name]] <-
        perturbed_risk
      
    }
  }
  
  return(risk_maps)
}


# ------------------------------------------------------------
# 10. Helper function: calculate ARC
# ------------------------------------------------------------

calculate_arc <- function(
    baseline_risk,
    perturbed_maps,
    epsilon = 1e-8
) {
  
  perturbed_stack <-
    rast(
      perturbed_maps
    )
  
  arc_stack <-
    c(
      baseline_risk,
      perturbed_stack
    )
  
  arc_function <- function(x) {
    
    reference_value <- x[1]
    perturbed_values <- x[-1]
    
    # Prevent division by values equal or extremely close to zero
    if (
      is.na(reference_value) ||
      abs(reference_value) <= epsilon
    ) {
      
      return(NA_real_)
      
    }
    
    mean(
      abs(
        (
          perturbed_values -
            reference_value
        ) /
          reference_value
      ),
      na.rm = TRUE
    )
  }
  
  arc_map <-
    app(
      arc_stack,
      fun = arc_function
    )
  
  return(arc_map)
}


# ============================================================
# DRY SEASON
# ============================================================


# ------------------------------------------------------------
# 11. Load dry-season reference raster
# ------------------------------------------------------------

rainfall_dry <- rast(
  here(
    input_dir_dry,
    "rainfall_dry_std.tif"
  )
)


if (is.na(crs(rainfall_dry))) {
  
  stop(
    "The dry-season reference raster has no defined CRS."
  )
}


# ------------------------------------------------------------
# 12. Load and align dry-season factor rasters
# ------------------------------------------------------------

rasters_dry <- list(
  
  rainfall =
    rainfall_dry,
  
  bushmeat_hunting =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "Relative_hunting_pressure_index_1km_Cameroon.tif"
        )
      ),
      rainfall_dry
    ),
  
  bushmeat_trade =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "Standardized_CMR_bushmeat_trade_distance_based_0_1_1km_UTM33N.tif"
        )
      ),
      rainfall_dry
    ),
  
  deforestation =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "DeforestationRisk_1km.tif"
        )
      ),
      rainfall_dry
    ),
  
  forests =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "Standardized_prop_foret_1km_2020_CMR.tif"
        )
      ),
      rainfall_dry
    ),
  
  NHP_density =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "NHP_UTM33N_1km.tif"
        )
      ),
      rainfall_dry
    ),
  
  bat_density =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
        )
      ),
      rainfall_dry
    ),
  
  poverty =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "Risk_index_RWI_1km_UTM33N_CMR.tif"
        )
      ),
      rainfall_dry
    ),
  
  human_population_density =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
        )
      ),
      rainfall_dry
    ),
  
  temperature =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "MDTR_dry_utm33_1km_std01_minmax.tif"
        )
      ),
      rainfall_dry
    ),
  
  elevation =
    align_to_reference(
      rast(
        here(
          input_dir_dry,
          "elevation_std.tif"
        )
      ),
      rainfall_dry
    )
)


# ------------------------------------------------------------
# 13. Check dry-season raster geometry
# ------------------------------------------------------------

if (
  !check_geometry(
    rainfall_dry,
    rasters_dry
  )
) {
  
  stop(
    "Dry-season sensitivity rasters are not correctly aligned."
  )
}


message(
  "Dry-season sensitivity rasters successfully aligned."
)


# ------------------------------------------------------------
# 14. Calculate dry-season baseline sensitivity model
# ------------------------------------------------------------

baseline_risk_dry <-
  calculate_wlc(
    rasters_dry,
    weights_dry
  )


names(baseline_risk_dry) <-
  "baseline_sensitivity_dry"


writeRaster(
  baseline_risk_dry,
  here(
    output_dir,
    "baseline_sensitivity_model_dry.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 15. Run dry-season weight perturbation analysis
# ------------------------------------------------------------

perturbed_maps_dry <-
  run_weight_perturbation(
    rasters = rasters_dry,
    baseline_weights = weights_dry,
    deltas = deltas
  )


# ------------------------------------------------------------
# 16. Calculate dry-season uncertainty surface
# ------------------------------------------------------------

perturbed_stack_dry <-
  rast(
    perturbed_maps_dry
  )


uncertainty_dry <-
  app(
    perturbed_stack_dry,
    fun = sd,
    na.rm = TRUE
  )


names(uncertainty_dry) <-
  "uncertainty_sd_dry"


writeRaster(
  uncertainty_dry,
  here(
    output_dir,
    "uncertainty_surface_dry.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 17. Calculate dry-season Average Relative Change
# ------------------------------------------------------------

arc_dry <-
  calculate_arc(
    baseline_risk = baseline_risk_dry,
    perturbed_maps = perturbed_maps_dry
  )


names(arc_dry) <-
  "average_relative_change_dry"


writeRaster(
  arc_dry,
  here(
    output_dir,
    "average_relative_change_dry.tif"
  ),
  overwrite = TRUE
)


message(
  "Dry-season sensitivity analysis completed."
)


# ============================================================
# WET SEASON
# ============================================================


# ------------------------------------------------------------
# 18. Load wet-season reference raster
# ------------------------------------------------------------

rainfall_wet <- rast(
  here(
    input_dir_wet,
    "rainfall_wet_std.tif"
  )
)


if (is.na(crs(rainfall_wet))) {
  
  stop(
    "The wet-season reference raster has no defined CRS."
  )
}


# ------------------------------------------------------------
# 19. Load and align wet-season factor rasters
# ------------------------------------------------------------

rasters_wet <- list(
  
  rainfall =
    rainfall_wet,
  
  bushmeat_hunting =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "Relative_hunting_pressure_index_1km_Cameroon.tif"
        )
      ),
      rainfall_wet
    ),
  
  bushmeat_trade =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "Standardized_CMR_bushmeat_trade_distance_based_0_1_1km_UTM33N.tif"
        )
      ),
      rainfall_wet
    ),
  
  deforestation =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "DeforestationRisk_1km.tif"
        )
      ),
      rainfall_wet
    ),
  
  forests =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "Standardized_prop_foret_1km_2020_CMR.tif"
        )
      ),
      rainfall_wet
    ),
  
  NHP_density =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "NHP_UTM33N_1km.tif"
        )
      ),
      rainfall_wet
    ),
  
  bat_density =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "African_Bat_occurrence_proximity_index_1km_UTM33N_01.tif"
        )
      ),
      rainfall_wet
    ),
  
  poverty =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "Risk_index_RWI_1km_UTM33N_CMR.tif"
        )
      ),
      rainfall_wet
    ),
  
  human_population_density =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "Standardized_CMR_pop_density_mean_2010_2020_scaled_0_1_1km_UTM33N.tif"
        )
      ),
      rainfall_wet
    ),
  
  temperature =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "MDTR_wet_utm33_1km_std01_minmax.tif"
        )
      ),
      rainfall_wet
    ),
  
  elevation =
    align_to_reference(
      rast(
        here(
          input_dir_wet,
          "elevation_std.tif"
        )
      ),
      rainfall_wet
    )
)


# ------------------------------------------------------------
# 20. Check wet-season raster geometry
# ------------------------------------------------------------

if (
  !check_geometry(
    rainfall_wet,
    rasters_wet
  )
) {
  
  stop(
    "Wet-season sensitivity rasters are not correctly aligned."
  )
}


message(
  "Wet-season sensitivity rasters successfully aligned."
)


# ------------------------------------------------------------
# 21. Calculate wet-season baseline sensitivity model
# ------------------------------------------------------------

baseline_risk_wet <-
  calculate_wlc(
    rasters_wet,
    weights_wet
  )


names(baseline_risk_wet) <-
  "baseline_sensitivity_wet"


writeRaster(
  baseline_risk_wet,
  here(
    output_dir,
    "baseline_sensitivity_model_wet.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 22. Run wet-season weight perturbation analysis
# ------------------------------------------------------------

perturbed_maps_wet <-
  run_weight_perturbation(
    rasters = rasters_wet,
    baseline_weights = weights_wet,
    deltas = deltas
  )


# ------------------------------------------------------------
# 23. Calculate wet-season uncertainty surface
# ------------------------------------------------------------

perturbed_stack_wet <-
  rast(
    perturbed_maps_wet
  )


uncertainty_wet <-
  app(
    perturbed_stack_wet,
    fun = sd,
    na.rm = TRUE
  )


names(uncertainty_wet) <-
  "uncertainty_sd_wet"


writeRaster(
  uncertainty_wet,
  here(
    output_dir,
    "uncertainty_surface_wet.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 24. Calculate wet-season Average Relative Change
# ------------------------------------------------------------

arc_wet <-
  calculate_arc(
    baseline_risk = baseline_risk_wet,
    perturbed_maps = perturbed_maps_wet
  )


names(arc_wet) <-
  "average_relative_change_wet"


writeRaster(
  arc_wet,
  here(
    output_dir,
    "average_relative_change_wet.tif"
  ),
  overwrite = TRUE
)


message(
  "Wet-season sensitivity analysis completed."
)


# ------------------------------------------------------------
# 25. Completion message
# ------------------------------------------------------------

message(
  paste(
    "Seasonal EVD sensitivity analyses completed successfully.",
    "Uncertainty and Average Relative Change surfaces were",
    "generated for both seasons."
  )
)