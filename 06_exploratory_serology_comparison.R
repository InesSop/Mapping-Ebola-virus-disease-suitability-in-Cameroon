# ============================================================
# Script: 06_exploratory_serology_comparison.R
#
# Purpose:
#   Conduct an exploratory spatial comparison between modelled
#   Ebola virus disease (EVD) spillover suitability and published
#   serological observations from Cameroon.
#
# Background:
#   Cameroon has not reported documented human EVD outbreaks or
#   confirmed clinical cases. Consequently, direct model validation
#   against observed human EVD occurrence was not possible.
#
#   The model outputs are therefore interpreted as spatial estimates
#   of relative EVD spillover suitability rather than predictions of
#   confirmed spillover occurrence.
#
#   Published serological observations are used as independent
#   spatial evidence to examine whether sampling sites showing
#   evidence of Ebola virus exposure occur in locations assigned
#   relatively high modelled suitability.
#
# Serological evidence:
#   - Human serological observations are used as the primary evidence
#     of potential human exposure.
#   - A strict human scenario includes confirmed seropositive sites.
#   - A sensitive human scenario additionally includes sites with
#     unconfirmed serological reactivity.
#   - Bat serological observations are treated as complementary
#     ecological evidence because seroreactivity may reflect
#     cross-reactivity and was not accompanied by viral RNA detection.
#
# Spatial comparison:
#   Because published sampling periods were not season-specific,
#   two non-season-specific composite suitability surfaces are used:
#
#     1. Mean suitability across dry and wet seasons
#     2. Maximum suitability across dry and wet seasons
#
#   Suitability is extracted:
#
#     1. At the raster cell containing each sampling location
#     2. As the mean suitability within a 10 km buffer
#
#   The 10 km buffer accounts for uncertainty in published sampling
#   locations and potential movements of humans and wildlife hosts.
#
# Outputs:
#   - Mean and maximum seasonal composite suitability rasters
#   - Site-level suitability values
#   - 10 km buffer mean suitability values
#   - Empirical percentiles relative to the national suitability
#     distribution
#   - Scenario-specific summaries
#   - Exploratory randomisation results
#   - Publication-ready tables
#
# Important:
#   This analysis is exploratory and should not be interpreted as
#   formal external validation of the EVD suitability model.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(sf)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)
library(here)


# ------------------------------------------------------------
# 2. Define input and output directories
# ------------------------------------------------------------

seasonal_model_dir <- here(
  "outputs",
  "combined_models"
)

serology_file <- here(
  "inputs",
  "serology",
  "EVD_serological_sites.xlsx"
)

output_dir <- here(
  "outputs",
  "exploratory_serology_comparison"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Load seasonal suitability surfaces
# ------------------------------------------------------------

suitability_dry <- rast(
  here(
    seasonal_model_dir,
    "combined_risk_index_EVD_dry.tif"
  )
)

suitability_wet <- rast(
  here(
    seasonal_model_dir,
    "combined_risk_index_EVD_wet.tif"
  )
)


names(suitability_dry) <- "dry_suitability"
names(suitability_wet) <- "wet_suitability"


# ------------------------------------------------------------
# 4. Verify seasonal raster geometry
# ------------------------------------------------------------

if (!compareGeom(
  suitability_dry,
  suitability_wet,
  stopOnError = FALSE
)) {
  
  stop(
    paste(
      "Dry- and wet-season suitability rasters",
      "do not have identical geometry."
    )
  )
}


if (is.na(crs(suitability_dry))) {
  stop(
    "The dry-season suitability raster has no defined CRS."
  )
}

if (is.na(crs(suitability_wet))) {
  stop(
    "The wet-season suitability raster has no defined CRS."
  )
}


# ------------------------------------------------------------
# 5. Create non-season-specific composite surfaces
# ------------------------------------------------------------

# Mean suitability represents the average spatial suitability
# across the two seasonal models.

suitability_mean <- (
  suitability_dry +
    suitability_wet
) / 2

names(suitability_mean) <-
  "mean_suitability"


# Maximum suitability represents the highest modelled suitability
# observed in either season.

suitability_max <- max(
  c(
    suitability_dry,
    suitability_wet
  ),
  na.rm = TRUE
)

names(suitability_max) <-
  "maximum_suitability"


# Export composite rasters

writeRaster(
  suitability_mean,
  here(
    output_dir,
    "EVD_mean_seasonal_suitability.tif"
  ),
  overwrite = TRUE
)

writeRaster(
  suitability_max,
  here(
    output_dir,
    "EVD_maximum_seasonal_suitability.tif"
  ),
  overwrite = TRUE
)


# Store surfaces to be evaluated

suitability_surfaces <- list(
  mean_suitability = suitability_mean,
  maximum_suitability = suitability_max
)


# ------------------------------------------------------------
# 6. Load published serological sampling sites
# ------------------------------------------------------------

serology <- read_excel(
  serology_file
)


# Expected minimum columns:
#
# site
# longitude
# latitude
# host_type
# validation_type
# scenario
#
# Recommended additional columns:
#
# reference
# sample_size
# confirmed_seropositive
# seroreactive
# pcr_positive
# evidence_level
# notes


required_columns <- c(
  "site",
  "longitude",
  "latitude",
  "host_type",
  "validation_type",
  "scenario"
)


missing_columns <- setdiff(
  required_columns,
  names(serology)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "The serological dataset is missing the following columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 7. Clean serological site information
# ------------------------------------------------------------

serology <- serology %>%
  mutate(
    site = as.character(site),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    host_type = as.character(host_type),
    validation_type = as.character(validation_type),
    scenario = as.character(scenario)
  ) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude)
  )


# ------------------------------------------------------------
# 8. Create spatial point layer
# ------------------------------------------------------------

# Published coordinates are assumed to be supplied in WGS84.

serology_sf_wgs84 <- st_as_sf(
  serology,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
)


# Reproject sampling sites to the suitability raster CRS.

serology_sf_raster <- st_transform(
  serology_sf_wgs84,
  crs(suitability_mean)
)

serology_points_vect <- vect(
  serology_sf_raster
)


# ------------------------------------------------------------
# 9. Create 10 km buffers
# ------------------------------------------------------------

# A projected CRS expressed in metres is required for buffering.
# EPSG:32633 corresponds to WGS 84 / UTM zone 33N.

buffer_crs <- 32633


serology_sf_metric <- st_transform(
  serology_sf_wgs84,
  buffer_crs
)


serology_buffer_10km_metric <- st_buffer(
  serology_sf_metric,
  dist = 10000
)


# Reproject buffers to raster CRS before extraction.

serology_buffer_10km <- st_transform(
  serology_buffer_10km_metric,
  crs(suitability_mean)
)


serology_buffer_10km_vect <- vect(
  serology_buffer_10km
)


# ------------------------------------------------------------
# 10. Function to extract suitability values
# ------------------------------------------------------------

extract_serology_values <- function(
    raster_layer,
    surface_name
) {
  
  # Exact raster-cell value
  pixel_values <- terra::extract(
    raster_layer,
    serology_points_vect
  )
  
  
  # Mean suitability within 10 km
  buffer_values <- terra::extract(
    raster_layer,
    serology_buffer_10km_vect,
    fun = mean,
    na.rm = TRUE
  )
  
  
  data.frame(
    serology,
    surface = surface_name,
    suitability_pixel = pixel_values[[2]],
    suitability_10km_mean = buffer_values[[2]]
  )
}


# ------------------------------------------------------------
# 11. Extract suitability at serological sampling locations
# ------------------------------------------------------------

serology_extractions <- map_dfr(
  names(suitability_surfaces),
  function(surface_name) {
    
    extract_serology_values(
      raster_layer =
        suitability_surfaces[[surface_name]],
      surface_name =
        surface_name
    )
  }
)


# ------------------------------------------------------------
# 12. Generate background suitability sample
# ------------------------------------------------------------

# Background pixels describe the national suitability distribution.
# They are not interpreted as EVD absences.

set.seed(123)

n_background <- 10000


sample_background <- function(
    raster_layer,
    surface_name
) {
  
  sampled_points <- spatSample(
    raster_layer,
    size = n_background,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE,
    values = TRUE
  )
  
  
  sampled_df <- as.data.frame(
    sampled_points
  )
  
  
  value_column <- names(raster_layer)[1]
  
  names(sampled_df)[
    names(sampled_df) == value_column
  ] <- "suitability"
  
  
  sampled_df %>%
    mutate(
      surface = surface_name
    )
}


background_values <- map_dfr(
  names(suitability_surfaces),
  function(surface_name) {
    
    sample_background(
      raster_layer =
        suitability_surfaces[[surface_name]],
      surface_name =
        surface_name
    )
  }
)


# ------------------------------------------------------------
# 13. Summarise background suitability
# ------------------------------------------------------------

background_summary <- background_values %>%
  group_by(surface) %>%
  summarise(
    n_background = n(),
    mean_suitability =
      mean(
        suitability,
        na.rm = TRUE
      ),
    median_suitability =
      median(
        suitability,
        na.rm = TRUE
      ),
    sd_suitability =
      sd(
        suitability,
        na.rm = TRUE
      ),
    q25 =
      quantile(
        suitability,
        0.25,
        na.rm = TRUE
      ),
    q75 =
      quantile(
        suitability,
        0.75,
        na.rm = TRUE
      ),
    q90 =
      quantile(
        suitability,
        0.90,
        na.rm = TRUE
      ),
    q95 =
      quantile(
        suitability,
        0.95,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 14. Function for empirical percentile calculation
# ------------------------------------------------------------

calculate_percentile <- function(
    value,
    background
) {
  
  if (is.na(value)) {
    return(NA_real_)
  }
  
  mean(
    background <= value,
    na.rm = TRUE
  ) * 100
}


# ------------------------------------------------------------
# 15. Calculate empirical percentiles
# ------------------------------------------------------------

serology_comparison <- serology_extractions %>%
  rowwise() %>%
  mutate(
    
    percentile_pixel =
      calculate_percentile(
        suitability_pixel,
        background_values$suitability[
          background_values$surface ==
            surface
        ]
      ),
    
    percentile_10km =
      calculate_percentile(
        suitability_10km_mean,
        background_values$suitability[
          background_values$surface ==
            surface
        ]
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 16. Add background context
# ------------------------------------------------------------

serology_comparison <- serology_comparison %>%
  left_join(
    background_summary,
    by = "surface"
  ) %>%
  mutate(
    
    pixel_minus_background_mean =
      suitability_pixel -
      mean_suitability,
    
    buffer_10km_minus_background_mean =
      suitability_10km_mean -
      mean_suitability,
    
    pixel_above_background_mean =
      suitability_pixel >
      mean_suitability,
    
    buffer_10km_above_background_mean =
      suitability_10km_mean >
      mean_suitability
  )


# ------------------------------------------------------------
# 17. Create long-format table
# ------------------------------------------------------------

comparison_long <- serology_comparison %>%
  select(
    site,
    longitude,
    latitude,
    host_type,
    validation_type,
    scenario,
    surface,
    suitability_pixel,
    suitability_10km_mean,
    percentile_pixel,
    percentile_10km
  ) %>%
  pivot_longer(
    cols = c(
      suitability_pixel,
      suitability_10km_mean
    ),
    names_to = "extraction_scale",
    values_to = "suitability"
  ) %>%
  mutate(
    
    percentile = case_when(
      
      extraction_scale ==
        "suitability_pixel" ~
        percentile_pixel,
      
      extraction_scale ==
        "suitability_10km_mean" ~
        percentile_10km
      
    ),
    
    extraction_scale = recode(
      extraction_scale,
      "suitability_pixel" =
        "Pixel",
      "suitability_10km_mean" =
        "10 km buffer mean"
    ),
    
    surface = recode(
      surface,
      "mean_suitability" =
        "Mean suitability",
      "maximum_suitability" =
        "Maximum suitability"
    )
  )


# ------------------------------------------------------------
# 18. Define evidence scenarios
# ------------------------------------------------------------

# strict_human:
#   confirmed human seropositive sites only
#
# sensitive_human:
#   confirmed human sites plus sites with unconfirmed human
#   serological reactivity
#
# ecological:
#   wildlife serological evidence considered separately


comparison_long <- comparison_long %>%
  mutate(
    
    scenario_group = case_when(
      
      scenario == "strict" &
        host_type == "Human" ~
        "strict_human",
      
      scenario %in% c(
        "strict",
        "sensitive"
      ) &
        host_type == "Human" ~
        "sensitive_human",
      
      scenario == "ecological" |
        validation_type == "ecological" ~
        "ecological",
      
      TRUE ~
        "other"
    )
  )


# ------------------------------------------------------------
# 19. Summarise comparison by evidence scenario
# ------------------------------------------------------------

scenario_summary <- comparison_long %>%
  filter(
    scenario_group != "other"
  ) %>%
  group_by(
    surface,
    scenario_group,
    extraction_scale
  ) %>%
  summarise(
    
    n_sites = n(),
    
    mean_suitability =
      mean(
        suitability,
        na.rm = TRUE
      ),
    
    median_suitability =
      median(
        suitability,
        na.rm = TRUE
      ),
    
    mean_percentile =
      mean(
        percentile,
        na.rm = TRUE
      ),
    
    median_percentile =
      median(
        percentile,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 20. Descriptive empirical percentile thresholds
# ------------------------------------------------------------

# These values are descriptive reference points only.
# They do not define suitability classes.

percentile_summary <- comparison_long %>%
  filter(
    scenario_group != "other"
  ) %>%
  mutate(
    
    above_50 =
      percentile >= 50,
    
    above_75 =
      percentile >= 75,
    
    above_90 =
      percentile >= 90
    
  ) %>%
  group_by(
    surface,
    scenario_group,
    extraction_scale
  ) %>%
  summarise(
    
    n_sites = n(),
    
    n_above_50 =
      sum(
        above_50,
        na.rm = TRUE
      ),
    
    n_above_75 =
      sum(
        above_75,
        na.rm = TRUE
      ),
    
    n_above_90 =
      sum(
        above_90,
        na.rm = TRUE
      ),
    
    proportion_above_50 =
      n_above_50 / n_sites,
    
    proportion_above_75 =
      n_above_75 / n_sites,
    
    proportion_above_90 =
      n_above_90 / n_sites,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 21. Exploratory randomisation function
# ------------------------------------------------------------

randomisation_test <- function(
    observed_values,
    background_values,
    n_iter = 10000
) {
  
  observed_values <-
    observed_values[
      is.finite(
        observed_values
      )
    ]
  
  background_values <-
    background_values[
      is.finite(
        background_values
      )
    ]
  
  
  n_sites <-
    length(
      observed_values
    )
  
  
  if (n_sites == 0) {
    
    return(
      data.frame(
        n_sites = 0,
        observed_mean = NA_real_,
        background_mean = NA_real_,
        empirical_p_value = NA_real_,
        observed_percentile_in_null = NA_real_
      )
    )
  }
  
  
  observed_mean <-
    mean(
      observed_values
    )
  
  
  null_means <- replicate(
    n_iter,
    {
      
      mean(
        sample(
          background_values,
          size = n_sites,
          replace = TRUE
        )
      )
      
    }
  )
  
  
  empirical_p_value <-
    (
      sum(
        null_means >= observed_mean
      ) + 1
    ) /
    (
      n_iter + 1
    )
  
  
  observed_percentile <-
    mean(
      null_means <= observed_mean
    ) * 100
  
  
  data.frame(
    
    n_sites =
      n_sites,
    
    observed_mean =
      observed_mean,
    
    background_mean =
      mean(
        background_values
      ),
    
    empirical_p_value =
      empirical_p_value,
    
    observed_percentile_in_null =
      observed_percentile
  )
}


# ------------------------------------------------------------
# 22. Run exploratory randomisation comparisons
# ------------------------------------------------------------

set.seed(123)


randomisation_input <- comparison_long %>%
  filter(
    extraction_scale ==
      "Pixel",
    scenario_group !=
      "other"
  )


randomisation_results <- expand_grid(
  
  surface =
    unique(
      randomisation_input$surface
    ),
  
  scenario_group =
    unique(
      randomisation_input$scenario_group
    )
  
) %>%
  
  pmap_dfr(
    function(
    surface,
    scenario_group
    ) {
      
      observed_values <-
        randomisation_input %>%
        filter(
          .data$surface ==
            surface,
          .data$scenario_group ==
            scenario_group
        ) %>%
        pull(
          suitability
        )
      
      
      internal_surface <- recode(
        surface,
        "Mean suitability" =
          "mean_suitability",
        "Maximum suitability" =
          "maximum_suitability"
      )
      
      
      national_background <-
        background_values %>%
        filter(
          .data$surface ==
            internal_surface
        ) %>%
        pull(
          suitability
        )
      
      
      result <-
        randomisation_test(
          observed_values,
          national_background,
          n_iter = 10000
        )
      
      
      result %>%
        mutate(
          surface =
            surface,
          scenario_group =
            scenario_group
        )
    }
  )


# ------------------------------------------------------------
# 23. Create manuscript-oriented wide table
# ------------------------------------------------------------

manuscript_table <- serology_comparison %>%
  select(
    site,
    host_type,
    validation_type,
    scenario,
    surface,
    suitability_pixel,
    suitability_10km_mean,
    percentile_pixel,
    percentile_10km
  ) %>%
  
  mutate(
    surface = recode(
      surface,
      "mean_suitability" =
        "mean",
      "maximum_suitability" =
        "maximum"
    )
  ) %>%
  
  pivot_wider(
    names_from =
      surface,
    
    values_from = c(
      suitability_pixel,
      suitability_10km_mean,
      percentile_pixel,
      percentile_10km
    ),
    
    names_glue =
      "{surface}_{.value}"
  )


# ------------------------------------------------------------
# 24. Export results
# ------------------------------------------------------------

write_xlsx(
  
  list(
    
    site_values =
      serology_comparison,
    
    comparison_long =
      comparison_long,
    
    background_summary =
      background_summary,
    
    scenario_summary =
      scenario_summary,
    
    percentile_summary =
      percentile_summary,
    
    randomisation =
      randomisation_results,
    
    manuscript_table =
      manuscript_table
  ),
  
  path = here(
    output_dir,
    "EVD_exploratory_serology_comparison.xlsx"
  )
)


write.csv(
  
  manuscript_table,
  
  here(
    output_dir,
    "EVD_exploratory_serology_manuscript_table.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Completion message
# ------------------------------------------------------------

message(
  paste(
    "Exploratory comparison between EVD spillover suitability",
    "and published serological observations completed successfully."
  )
)

message(
  paste(
    "These results should be interpreted as an exploratory",
    "assessment of spatial consistency rather than as formal",
    "external validation of EVD occurrence."
  )
)