# ── Dependencies ──────────────────────────────────────────────────────────────
install.packages(c("tidyverse", "gbm", "caret", "jsonlite", "pROC", "terra"))

library(tidyverse)
library(gbm)
library(caret)
library(jsonlite)
library(pROC)
library(terra)
library(tmap)

# ── 1. Load occurrences ───────────────────────────────────────────────────────
load_occurrences <- function(path) {
  df <- read_csv(path) |>
    filter(year <= 2025,
           lat  >= 25, lat  <= 50,
           lon  >= -130, lon <= -65)
  cat(sprintf("Loaded %d records (%d-%d)\n", nrow(df), min(df$year), max(df$year)))
  df
}

occurrences <- load_occurrences("./Data/insects_clean.csv")

# ── 2. Load and stack WorldClim rasters ───────────────────────────────────────
cat("Loading WorldClim rasters...\n")

bio1  <- rast("worldclim/bio1_mean_temp.tif")
bio4  <- rast("worldclim/bio4_temp_seasonality.tif")
bio6  <- rast("worldclim/bio6_min_temp.tif")
bio12 <- rast("worldclim/bio12_annual_precip.tif")

# Stack into one object and name the layers
clim <- c(bio1, bio4, bio6, bio12)
names(clim) <- c("bio1", "bio4", "bio6", "bio12")

# Crop to study extent
study_ext <- ext(-130, -65, 25, 50)
clim <- crop(clim, study_ext)
cat("Climate rasters loaded and cropped.\n")

# ── 3. Load ash basal area data ───────────────────────────────────────────────
cat("Loading ash basal area data...\n")

ash_raw <- read_csv("./Data/ash_ba.csv") |>
  rename_with(tolower) |>                        # lowercase all column names
  rename(lat = lat, lon = lon, baa = baa) |>
  filter(lat >= 25, lat <= 50,
         lon >= -130, lon <= -65,
         baa >= 0)

# Aggregate to 0.5 degree grid (same resolution as prediction grid)
resolution <- 0.5
ash_grid <- ash_raw |>
  mutate(
    lat_bin = floor(lat / resolution) * resolution + resolution / 2,
    lon_bin = floor(lon / resolution) * resolution + resolution / 2
  ) |>
  group_by(lat_bin, lon_bin) |>
  summarise(mean_baa = mean(baa, na.rm = TRUE), .groups = "drop") |>
  rename(lat = lat_bin, lon = lon_bin)

cat(sprintf("Ash grid: %d cells with BAA data\n", nrow(ash_grid)))

# ── 4. Extract climate values at occurrence + background points ───────────────

# Helper: extract all 4 climate vars at a set of points
extract_climate <- function(lats, lons) {
  pts <- cbind(lons, lats)          # terra expects (x=lon, y=lat)
  vals <- terra::extract(clim, pts)
  tibble(
    bio1  = vals$bio1  / 10,        # WorldClim temps stored as C x 10
    bio4  = vals$bio4  / 10,
    bio6  = vals$bio6  / 10,
    bio12 = vals$bio12              # precipitation already in mm
  )
}

# Helper: look up ash BAA for a set of points via grid join
extract_ash <- function(lats, lons) {
  pts <- tibble(lat = lats, lon = lons) |>
    mutate(
      lat_bin = floor(lat / resolution) * resolution + resolution / 2,
      lon_bin = floor(lon / resolution) * resolution + resolution / 2
    ) |>
    left_join(ash_grid, by = c("lat_bin" = "lat", "lon_bin" = "lon")) |>
    mutate(mean_baa = replace_na(mean_baa, 0))   # no FIA plot = assume no ash
  pts$mean_baa
}

# ── 5. Background points ──────────────────────────────────────────────────────
set.seed(42)
n_bg <- 8000

background <- tibble(
  lat  = runif(n_bg, 25.0, 50.0),
  lon  = runif(n_bg, -130.0, -65.0),
  year = sample(2001:2025, n_bg, replace = TRUE)
)

# ── 6. Build feature matrix ───────────────────────────────────────────────────
cat("Extracting climate values at occurrence points...\n")
occ_clim <- extract_climate(occurrences$lat, occurrences$lon)
occ_ash  <- extract_ash(occurrences$lat, occurrences$lon)

cat("Extracting climate values at background points...\n")
bg_clim  <- extract_climate(background$lat, background$lon)
bg_ash   <- extract_ash(background$lat, background$lon)

presence_df <- occurrences |>
  bind_cols(occ_clim) |>
  mutate(ash_baa = occ_ash, presence = 1)

background_df <- background |>
  bind_cols(bg_clim) |>
  mutate(ash_baa = bg_ash, presence = 0)

combined <- bind_rows(presence_df, background_df) |>
  drop_na()                                        # drop rows where raster had no data

y <- combined$presence
X <- combined |> select(lat, lon, year, bio1, bio4, bio6, bio12, ash_baa)

# Scale features
scaler_means <- colMeans(X)
scaler_sds   <- apply(X, 2, sd)
X_scaled <- scale(X, center = scaler_means, scale = scaler_sds) |> as_tibble()
train_df <- X_scaled |> mutate(presence = y)

cat(sprintf("Training set: %d presences | %d background\n", sum(y), sum(y == 0)))

# ── 7. Train GBM model ────────────────────────────────────────────────────────
cat("Running 5-fold cross-validation...\n")
folds <- createFolds(y, k = 5, returnTrain = TRUE)

auc_scores <- map_dbl(folds, function(train_idx) {
  trn <- train_df[train_idx, ]
  val <- train_df[-train_idx, ]

  m <- gbm(
    presence ~ lat + lon + year + bio1 + bio4 + bio6 + bio12 + ash_baa,
    data              = trn,
    distribution      = "bernoulli",
    n.trees           = 300,
    interaction.depth = 4,
    shrinkage         = 0.05,
    bag.fraction      = 0.8,
    verbose           = FALSE
  )
  preds <- predict(m, val, n.trees = 300, type = "response")
  as.numeric(auc(roc(val$presence, preds, quiet = TRUE)))
})

cat(sprintf("Cross-validated AUC: %.3f +/- %.3f\n", mean(auc_scores), sd(auc_scores)))

# Final model on full data
model <- gbm(
  presence ~ lat + lon + year + bio1 + bio4 + bio6 + bio12 + ash_baa,
  data              = train_df,
  distribution      = "bernoulli",
  n.trees           = 300,
  interaction.depth = 4,
  shrinkage         = 0.05,
  bag.fraction      = 0.8,
  verbose           = FALSE
)
cat("Model trained.\n")

# Variable importance
cat("\nVariable importance:\n")
importance <- summary(model, plotit = FALSE)
print(importance)

# ── 8. Prediction grid ────────────────────────────────────────────────────────
resolution  <- 0.5
spread_rate <- 0.45

grid_base <- expand.grid(
  lat = seq(25.0, 50.0, by = resolution),
  lon = seq(-130.0, -65.0, by = resolution)
) |> as_tibble()

# Extract climate + ash for every grid cell (done once, reused for all years)
cat("Extracting climate values for prediction grid...\n")
grid_clim <- extract_climate(grid_base$lat, grid_base$lon)
grid_ash  <- extract_ash(grid_base$lat, grid_base$lon)

grid_base <- grid_base |>
  bind_cols(grid_clim) |>
  mutate(ash_baa = grid_ash)

# Scale grid using training parameters
scale_grid <- function(df) {
  feat <- df |> select(lat, lon, year, bio1, bio4, bio6, bio12, ash_baa)
  scaled <- scale(feat, center = scaler_means, scale = scaler_sds) |> as_tibble()
  scaled
}

# Min distance to known occurrences (for dispersal kernel)
cat("Computing dispersal distances...\n")
occ_lats <- occurrences$lat
occ_lons <- occurrences$lon

min_dists <- map2_dbl(grid_base$lat, grid_base$lon, function(la, lo) {
  min(sqrt((la - occ_lats)^2 + (lo - occ_lons)^2))
})

# ── 9. Project future years ───────────────────────────────────────────────────
predict_year <- function(yr) {
  grid <- grid_base |> mutate(year = yr)
  grid_scaled <- scale_grid(grid)

  base_suit <- predict(model, grid_scaled, n.trees = 300, type = "response")

  # Dispersal kernel
  extra_years <- max(0, yr - 2025)
  max_reach   <- extra_years * spread_rate
  dispersal_w <- ifelse(
    min_dists <= 0.5,
    1.0,
    1.0 / (1.0 + exp(4.0 * (min_dists - 0.5 - max_reach)))
  )

  # Ash constraint — no ash = no EAB
  ash_w <- ifelse(grid$ash_baa > 0, 1.0, 0.0)

  grid |>
    mutate(suitability = base_suit * dispersal_w * ash_w) |>
    filter(suitability > 0.05)
}

target_years <- c(2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065)
predictions  <- map(target_years, predict_year)
names(predictions) <- target_years

walk2(target_years, predictions, function(yr, p) {
  high <- sum(p$suitability > 0.7)
  cat(sprintf("%d: %d cells, %d high-risk\n", yr, nrow(p), high))
})

# ── 10. Export JSON for Leaflet ───────────────────────────────────────────────
leaflet_data <- list(
  years = map2(target_years, predictions, function(yr, grid) {
    list(
      year  = yr,
      cells = map(seq_len(nrow(grid)), function(i) {
        list(
          lat = round(grid$lat[i], 4),
          lon = round(grid$lon[i], 4),
          s   = round(grid$suitability[i], 3)
        )
      })
    )
  }),
  occurrences = map(seq_len(nrow(occurrences)), function(i) {
    list(
      lat  = round(occurrences$lat[i], 4),
      lon  = round(occurrences$lon[i], 4),
      year = occurrences$year[i]
    )
  })
)

write_json(leaflet_data, "eab_leaflet_data.js", auto_unbox = TRUE, pretty = FALSE)
cat(sprintf("Exported: eab_leaflet_data.js (%.1f KB)\n",
            file.size("eab_leaflet_data.js") / 1024))


#---------------------------------------


