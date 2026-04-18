# =============================================================================
# Emerald Ash Borer (EAB) — Species Distribution Model + Spread Simulation
# =============================================================================
#
# Pipeline overview:
#   1. Load & prepare occurrence, background, and climate raster data
#   2. Fit SDM using maxnet (MaxEnt implementation in R)
#   3. Project suitability surface across CONUS
#   4. Run spread simulation in two phases:
#        PHASE 1 (2001–2025): step through observed data year by year,
#          locking in real detections as they accumulate — no simulation needed,
#          this reconstructs the observed invasion history
#        PHASE 2 (2026–2040): stratified diffusion forward projection from the
#          confirmed 2025 baseline
#          - Local diffusion:  Gaussian dispersal kernel (~20 km/yr)
#          - Long-distance:    Human-mediated jump dispersal
#   5. Export annual spread maps as PNG + animated GIF (2001–2040)
#
# Inputs (place all in `data/` subfolder relative to this script):
#   - insects_clean.csv          EAB occurrence points  (lat, lon, year)
#   - ash_ba.csv                 FIA ash basal area     (LAT, LON, BAA, ...)
#   - bio6_min_temp.tif          BIO6  Min temp coldest month   (°C)
#   - bio4_temp_seasonality.tif  BIO4  Temperature seasonality
#   - bio12_annual_precip.tif    BIO12 Annual precipitation     (mm)
#
# Outputs (written to `outputs/`):
#   - eab_suitability.tif        Projected suitability raster (0–1)
#   - spread_YYYY.png            Annual spread maps (one per year)
#   - eab_spread_animation.gif   Animated GIF of spread 2001–2040
#   - model_diagnostics.png      ROC, variable response, importance
#
# Required packages (install once, then comment out):
#   install.packages(c("terra", "maxnet", "sf", "dplyr", "ggplot2",
#                      "pROC", "gifski", "scales", "viridis", "maps"))
#
# Author:  Jack Ross
# R version tested: 4.3+
# =============================================================================


# ── 0. Setup ------------------------------------------------------------------

library(terra)       # raster operations (replaces raster package)
library(maxnet)      # MaxEnt via glmnet (Philips et al.)
library(sf)          # vector / points
library(dplyr)       # data wrangling
library(ggplot2)     # plotting
library(pROC)        # AUC / ROC
library(gifski)      # GIF animation
library(scales)      # color scales
library(viridis)     # perceptually uniform palettes

# ── SET THIS TO YOUR FOLDER LOCATION ─────────────────────────────────────────
# Point this at the eab_model folder you unzipped. Examples:
#   Windows:  SCRIPT_DIR <- "C:/Users/Jack Ross/Documents/eab_model"
#   Mac/Linux: SCRIPT_DIR <- "~/Documents/eab_model"
#
# TIP: In RStudio you can also use Session > Set Working Directory >
#      To Source File Location, then set SCRIPT_DIR <- getwd()

SCRIPT_DIR <- getwd()   # <-- change this if needed

DATA_DIR   <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Confirm paths look right before proceeding
cat("Looking for data in:", DATA_DIR, "\n")
stopifnot("data/ folder not found — update SCRIPT_DIR above" =
            dir.exists(DATA_DIR))

# Reproducibility
set.seed(42)

# Helper: convert a terra SpatRaster layer to a ggplot2-ready data frame
# Replaces tidyterra::geom_spatraster() — no extra package needed
rast_to_df <- function(r, name = "value") {
  df        <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", name)
  df[!is.na(df[[name]]), ]
}


# =============================================================================
# 1. LOAD DATA
# =============================================================================

# ── 1a. EAB occurrences -------------------------------------------------------
eab_all <- read.csv(file.path(DATA_DIR, "insects_clean.csv")) |>
  filter(!is.na(lat), !is.na(lon)) |>
  mutate(presence = 1L)

# Full record used for: invasion history (2001–2025) and map overlays
# SDM subset (2013+) used for: model fitting only
# Pre-2013 has only 47 records clustered near Detroit — too sparse and biased
# to reliably characterise the climate niche; including them would pull the
# model toward that one location. They are still shown on all maps.
SDM_YEAR_MIN <- 2013
eab_sdm <- eab_all |> filter(year >= SDM_YEAR_MIN, year <= 2025)
eab     <- eab_all |> filter(year <= 2025)   # full confirmed record

cat(sprintf("EAB full record: %d points (years %d–%d)\n",
            nrow(eab), min(eab$year), max(eab$year)))
cat(sprintf("EAB SDM subset: %d points (years %d–2025)\n",
            nrow(eab_sdm), SDM_YEAR_MIN))

# ── 1b. Ash basal area (host availability) ------------------------------------
ash <- read.csv(file.path(DATA_DIR, "ash_ba.csv")) |>
  filter(!is.na(LAT), !is.na(LON), BAA > 0) |>
  rename(lat = LAT, lon = LON, baa = BAA)

cat(sprintf("Ash FIA plots with BAA > 0: %d\n", nrow(ash)))

# ── 1c. Climate rasters -------------------------------------------------------
# Stack the 3 selected variables (VIF-filtered from diagnostic step)
clim <- c(
  rast(file.path(DATA_DIR, "bio6_min_temp.tif")),
  rast(file.path(DATA_DIR, "bio4_temp_seasonality.tif")),
  rast(file.path(DATA_DIR, "bio12_annual_precip.tif"))
)
names(clim) <- c("bio6_min_temp", "bio4_temp_seasonality", "bio12_annual_precip")

cat(sprintf("Climate rasters: %d layers | CRS: %s | res: %.4f°\n",
            nlyr(clim), crs(clim, describe = TRUE)$code, res(clim)[1]))


# =============================================================================
# 2. PREPARE PRESENCE / BACKGROUND DATA
# =============================================================================

# ── 2a. Convert SDM occurrences to SpatVector --------------------------------
eab_sv <- vect(eab_sdm, geom = c("lon", "lat"), crs = "EPSG:4326")

# ── 2b. Spatially thin SDM occurrences (1 per raster cell) -------------------
# Thinning applied only to the SDM subset to reduce spatial autocorrelation
eab_cells  <- cellFromXY(clim[[1]], geom(eab_sv)[, c("x","y")])
eab_unique <- eab_sdm[!duplicated(eab_cells), ]
eab_sv     <- vect(eab_unique, geom = c("lon", "lat"), crs = "EPSG:4326")
cat(sprintf("SDM occurrences after spatial thinning: %d cells\n", nrow(eab_unique)))

# ── 2c. Generate ash-weighted background points ------------------------------
# Weight sampling by ash basal area — EAB can only establish where ash exists.
# This is a key improvement over purely random background.
n_bg       <- 10000
bg_weights <- ash$baa / sum(ash$baa)
bg_idx     <- sample(nrow(ash), size = n_bg, replace = TRUE, prob = bg_weights)
bg         <- ash[bg_idx, c("lat", "lon")] |>
  mutate(presence = 0L, year = NA_integer_)

# Combine presence + background
occ_bg <- bind_rows(
  eab_unique |> select(lat, lon, presence, year),
  bg
) |> filter(!is.na(lat), !is.na(lon))

# ── 2d. Extract climate at all points ----------------------------------------
occ_sv   <- vect(occ_bg, geom = c("lon", "lat"), crs = "EPSG:4326")
clim_ext <- extract(clim, occ_sv, ID = FALSE)

model_df <- bind_cols(occ_bg, clim_ext) |>
  filter(complete.cases(pick(names(clim))))

cat(sprintf("Model data: %d presences + %d background\n",
            sum(model_df$presence), sum(model_df$presence == 0)))


# =============================================================================
# 3. FIT SDM — maxnet (MaxEnt)
# =============================================================================

CLIM_VARS <- c("bio6_min_temp", "bio4_temp_seasonality", "bio12_annual_precip")

X <- as.matrix(model_df[, CLIM_VARS])
y <- model_df$presence

# maxnet uses feature classes: l=linear, q=quadratic, p=product, h=hinge
# 'lqh' is the standard MaxEnt feature set for moderate sample sizes
sdm <- maxnet(
  p    = y,
  data = as.data.frame(X),
  f    = maxnet.formula(y, as.data.frame(X), classes = "lqh"),
  regmult = 1.5   # regularization — increase if model looks overfit
)

cat("MaxEnt model fitted.\n")
cat("Non-zero coefficients:", sum(sdm$betas != 0), "\n")


# ── 3a. 5-fold cross-validation AUC ------------------------------------------
cv_auc <- numeric(5)
folds  <- sample(rep(1:5, length.out = nrow(model_df)))

for (k in 1:5) {
  train <- model_df[folds != k, ]
  test  <- model_df[folds == k, ]
  Xtr   <- as.matrix(train[, CLIM_VARS])
  Xte   <- as.matrix(test[,  CLIM_VARS])
  ytr   <- train$presence
  yte   <- test$presence

  m_k <- maxnet(
    p    = ytr,
    data = as.data.frame(Xtr),
    f    = maxnet.formula(ytr, as.data.frame(Xtr), classes = "lqh"),
    regmult = 1.5
  )
  pred_k  <- predict(m_k, as.data.frame(Xte), type = "cloglog")
  cv_auc[k] <- as.numeric(auc(roc(yte, as.numeric(pred_k), quiet = TRUE)))
}

cat(sprintf("5-fold CV AUC: %.3f ± %.3f\n", mean(cv_auc), sd(cv_auc)))


# =============================================================================
# 4. PROJECT SUITABILITY ACROSS CONUS — masked to ash host footprint
# =============================================================================
# We project the fitted SDM across the full raster extent, then mask to the
# ash basal area footprint from the FIA data. Cells with no recorded ash are
# set to NA — suitability is only meaningful where the host tree exists.
# The ash mask also serves as a natural western boundary, since ash is largely
# absent from the arid interior West.

suit_rast <- predict(clim, sdm, type = "cloglog", na.rm = TRUE)
names(suit_rast) <- "suitability"

ash_sv   <- vect(ash, geom = c("lon", "lat"), crs = "EPSG:4326")
ash_rast <- rasterize(ash_sv, clim[[1]], field = "baa", fun = "mean")
ash_mask <- !is.na(ash_rast)
suit_rast_masked <- mask(suit_rast, ash_mask, maskvalue = FALSE)

writeRaster(suit_rast_masked,
            file.path(OUTPUT_DIR, "eab_suitability.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")
cat("Suitability raster saved (ash-masked).\n")


# =============================================================================
# 5. MODEL DIAGNOSTICS PLOT
# =============================================================================

# ── 5a. Variable response curves (marginal) -----------------------------------
response_data <- lapply(CLIM_VARS, function(var) {
  x_range <- seq(
    quantile(model_df[[var]], 0.02, na.rm = TRUE),
    quantile(model_df[[var]], 0.98, na.rm = TRUE),
    length.out = 200
  )
  # Hold other vars at their mean
  newdat <- as.data.frame(
    lapply(CLIM_VARS, function(v) {
      if (v == var) x_range else rep(mean(model_df[[v]], na.rm = TRUE), 200)
    })
  )
  names(newdat) <- CLIM_VARS
  data.frame(
    variable = var,
    x        = x_range,
    suit     = as.numeric(predict(sdm, newdat, type = "cloglog"))
  )
}) |> bind_rows()

var_labels <- c(
  bio6_min_temp        = "BIO6: Min Temp Coldest Month (°C)",
  bio4_temp_seasonality = "BIO4: Temp Seasonality",
  bio12_annual_precip  = "BIO12: Annual Precipitation (mm)"
)

p_response <- ggplot(response_data, aes(x, suit)) +
  geom_line(color = "#00d4aa", linewidth = 1.2) +
  geom_ribbon(aes(ymin = 0, ymax = suit), fill = "#00d4aa", alpha = 0.15) +
  facet_wrap(~variable, scales = "free_x",
             labeller = labeller(variable = var_labels)) +
  labs(title    = "MaxEnt Response Curves",
       subtitle = sprintf("5-fold CV AUC = %.3f ± %.3f", mean(cv_auc), sd(cv_auc)),
       x = "Predictor value", y = "Predicted suitability") +
  theme_dark() +
  theme(
    plot.background  = element_rect(fill = "#0d1117", color = NA),
    panel.background = element_rect(fill = "#161b22"),
    panel.grid       = element_line(color = "#30363d"),
    strip.background = element_rect(fill = "#21262d"),
    strip.text       = element_text(color = "#e6edf3", face = "bold"),
    axis.text        = element_text(color = "#8b949e"),
    axis.title       = element_text(color = "#8b949e"),
    plot.title       = element_text(color = "#e6edf3", face = "bold", size = 14),
    plot.subtitle    = element_text(color = "#8b949e")
  ) +
  ylim(0, 1)

ggsave(file.path(OUTPUT_DIR, "model_diagnostics.png"),
       p_response, width = 12, height = 4, dpi = 150, bg = "#0d1117")
cat("Diagnostics plot saved.\n")


# =============================================================================
# 6. STRATIFIED DIFFUSION SPREAD SIMULATION
# =============================================================================
#
# Model: at each annual time step —
#   A) LOCAL diffusion:     Gaussian kernel convolution, σ = local_sd km
#   B) LONG-DISTANCE jumps: random jumps drawn from exponential tail,
#                           scaled by road density proxy (population density
#                           correlate; uses lat/lon grid noise as stand-in
#                           until road raster is available)
#   C) Establishment:       a cell becomes "invaded" only if
#                           suitability > SUIT_THRESHOLD
#
# The result is a binary invaded/not-invaded raster updated each year.
# =============================================================================

# ── 6a. Parameters -----------------------------------------------------------
OBS_START      <- 2001      # first year of real detections
OBS_END        <- 2025      # last year of real detection data
SIM_START      <- 2026      # first purely projected year
SIM_END        <- 2040      # forecast horizon
SUIT_THRESHOLD <- 0.30      # minimum suitability to allow new establishment
LOCAL_SIGMA_KM <- 10        # local diffusion σ in km/yr
JUMP_RATE      <- 0.008     # fraction of invaded cells generating a long-dist jump
JUMP_MEAN_KM   <- 93       # mean jump distance (km) — firewood transport scale
N_JUMP_MAX     <- 500       # cap on jumps per year

# ── 6b. Resample suitability to coarser grid for speed (~0.125° ≈ ~14 km) ----
suit_sim <- aggregate(suit_rast_masked, fact = 3, fun = "mean", na.rm = TRUE)
cat(sprintf("Simulation grid: %d × %d cells (%.4f° res)\n",
            nrow(suit_sim), ncol(suit_sim), res(suit_sim)[1]))

suit_vals  <- values(suit_sim)[, 1]
n_cells    <- ncell(suit_sim)
all_xy     <- xyFromCell(suit_sim, 1:n_cells)
cell_res_deg <- res(suit_sim)[1]
km_per_deg   <- 111.0

local_sigma_d <- LOCAL_SIGMA_KM / km_per_deg
jump_mean_d   <- JUMP_MEAN_KM   / km_per_deg

# ── 6c. Build Gaussian kernel -------------------------------------------------
build_gaussian_kernel <- function(sigma_cells, max_radius_cells = NULL) {
  if (is.null(max_radius_cells)) max_radius_cells <- ceiling(3 * sigma_cells)
  r    <- max_radius_cells
  size <- 2 * r + 1
  mat  <- matrix(0, size, size)
  for (i in seq_len(size)) for (j in seq_len(size)) {
    d2 <- (i - r - 1)^2 + (j - r - 1)^2
    mat[i, j] <- exp(-d2 / (2 * sigma_cells^2))
  }
  mat / sum(mat)
}

sigma_cells  <- local_sigma_d / cell_res_deg
gauss_kernel <- build_gaussian_kernel(sigma_cells)

convolve_invaded <- function(invaded_vec, template_rast, kernel) {
  r <- template_rast
  values(r) <- as.numeric(invaded_vec)
  foc <- focal(r, w = kernel, fun = "sum", na.policy = "omit", na.rm = TRUE)
  values(foc)[, 1]
}

valid_cells    <- which(!is.na(suit_vals))
suitable_cells <- which(!is.na(suit_vals) & suit_vals >= SUIT_THRESHOLD)

# Helper: snap detection points to simulation grid cells
pts_to_cells <- function(df) {
  cells <- cellFromXY(suit_sim, cbind(df$lon, df$lat))
  cells[!is.na(cells)]
}

# =============================================================================
# PHASE 1: Observed invasion history (2001–2025)
# =============================================================================
# For each year, the invaded state is simply the set of all raster cells that
# contain a confirmed EAB detection up to and including that year.
# No simulation runs here — we are reconstructing history from real data.
# This gives the webmap a year-by-year replay of the actual invasion.

cat("Phase 1: building observed invasion history (2001–2025)…\n")

annual_invaded <- list()
invaded        <- rep(FALSE, n_cells)

obs_years <- OBS_START:OBS_END

for (yr in obs_years) {
  # Add any new detections in this year
  new_pts   <- eab |> filter(year == yr)
  new_cells <- pts_to_cells(new_pts)
  if (length(new_cells) > 0) invaded[new_cells] <- TRUE
  annual_invaded[[as.character(yr)]] <- invaded
}

# Lock the full 2025 observed state — simulation can never remove these cells
observed_invaded <- invaded

cat(sprintf("  2001: %d invaded cells\n", sum(annual_invaded[["2001"]])))
cat(sprintf("  2025: %d invaded cells (confirmed baseline)\n", sum(invaded)))

# =============================================================================
# PHASE 2: Forward projection (2026–2040)
# =============================================================================
# Starting from the confirmed 2025 footprint, run stratified diffusion:
#   A) Gaussian local spread
#   B) Long-distance jump dispersal
#   C) Re-lock observed cells every step

cat("Phase 2: projecting spread 2026–2040…\n")
pb <- txtProgressBar(min = SIM_START, max = SIM_END, style = 3)

for (yr in SIM_START:SIM_END) {

  # ── A) Local diffusion ──────────────────────────────────────────────────────
  spread_pressure <- convolve_invaded(invaded, suit_sim, gauss_kernel)

  new_local <- !invaded &
               !is.na(suit_vals) &
               suit_vals >= SUIT_THRESHOLD &
               spread_pressure > 0.05

  # ── B) Long-distance jump dispersal ─────────────────────────────────────────
  invaded_idx <- which(invaded)
  n_jumps     <- min(round(length(invaded_idx) * JUMP_RATE), N_JUMP_MAX)

  new_jump <- rep(FALSE, n_cells)
  if (n_jumps > 0 && length(invaded_idx) > 0) {
    src_cells <- sample(invaded_idx, n_jumps, replace = TRUE)
    src_xy    <- all_xy[src_cells, , drop = FALSE]

    jump_dist  <- rexp(n_jumps, rate = 1 / jump_mean_d)
    jump_angle <- runif(n_jumps, 0, 2 * pi)
    dest_x     <- src_xy[, 1] + jump_dist * cos(jump_angle)
    dest_y     <- src_xy[, 2] + jump_dist * sin(jump_angle)

    dest_cells <- cellFromXY(suit_sim, cbind(dest_x, dest_y))
    dest_cells <- dest_cells[!is.na(dest_cells)]
    dest_cells <- dest_cells[
      !is.na(suit_vals[dest_cells]) &
      suit_vals[dest_cells] >= SUIT_THRESHOLD
    ]
    if (length(dest_cells) > 0) new_jump[dest_cells] <- TRUE
  }

  # ── C) Update + re-lock observed baseline ───────────────────────────────────
  invaded <- invaded | new_local | new_jump | observed_invaded

  annual_invaded[[as.character(yr)]] <- invaded
  setTxtProgressBar(pb, yr)
}
close(pb)

n_inv_final <- sum(annual_invaded[[as.character(SIM_END)]], na.rm = TRUE)
area_km2    <- n_inv_final * (cell_res_deg * km_per_deg)^2
cat(sprintf("\nSimulation complete.\nProjected invaded cells by %d: %d (~%.0f km²)\n",
            SIM_END, n_inv_final, area_km2))

# =============================================================================
# 7. PLOT ANNUAL SPREAD MAPS + ANIMATED GIF
# =============================================================================

# ── 7a. Years to render -------------------------------------------------------
# Skip years with fewer than 10 cumulative detections (sparse early years
# produce blocky isolated cells that misrepresent the invasion state).
# The full annual_invaded list is still complete for the summary table and
# webmap export — this filter only affects the PNG/GIF output.
all_sim_years <- as.integer(names(annual_invaded))

det_by_year <- sapply(all_sim_years, function(yr)
  nrow(eab |> filter(year <= min(yr, OBS_END)))
)

sim_years <- all_sim_years[
  all_sim_years >= SIM_START |          # always include all projected years
  det_by_year   >= 10                   # observed years only if ≥ 10 detections
]

cat(sprintf("Rendering %d years (skipping %d sparse early years before 2013)\n",
            length(sim_years),
            sum(all_sim_years < SIM_START & det_by_year < 10)))

# ── 7b. Color palette ---------------------------------------------------------
suit_pal   <- colorRampPalette(c("#0d1117","#1a3a4a","#0d6e8a",
                                  "#00b4a0","#7dde82","#ffd60a",
                                  "#ff6b35","#c9184a"))(100)
obs_col  <- "#3fb950"   # observed spread colour (green)
proj_col <- "#e63946"   # projected spread colour (red)

# ── 7c. Generate one PNG per year then stitch into GIF -----------------------
png_paths <- character(length(sim_years))

# US state + country boundaries
if (requireNamespace("maps", quietly = TRUE)) {
  us_states <- sf::st_as_sf(maps::map("state",   plot = FALSE, fill = TRUE))
  us_nation <- sf::st_as_sf(maps::map("usa",     plot = FALSE, fill = TRUE))
  canada    <- sf::st_as_sf(maps::map("world",   regions = "Canada",
                                       plot = FALSE, fill = TRUE))
} else {
  us_states <- us_nation <- canada <- NULL
  message("Install the 'maps' package for boundary overlays.")
}

cat("Generating annual maps…\n")
for (i in seq_along(sim_years)) {
  yr <- sim_years[i]

  # Invaded raster for this year
  inv_r         <- suit_sim
  values(inv_r) <- as.numeric(annual_invaded[[as.character(yr)]])
  inv_r[inv_r == 0] <- NA

  # Full CONUS extent
  ext_plot  <- ext(-130, -65, 25, 50)
  suit_crop <- crop(suit_rast_masked, ext_plot)
  inv_crop  <- crop(inv_r,  ext_plot)

  is_observed <- yr <= OBS_END

  # Per-year label and subtitle
  year_label <- if (is_observed)
    sprintf("EAB Observed Invasion — %d", yr)
  else
    sprintf("EAB Projected Spread — %d", yr)

  n_det_yr    <- nrow(eab |> filter(year <= yr))
  area_km2_yr <- sum(annual_invaded[[as.character(yr)]], na.rm = TRUE) *
                   (cell_res_deg * km_per_deg)^2

  subtitle_txt <- if (is_observed)
    sprintf("Cumulative detections: %d  |  Invaded area: %s km²",
            n_det_yr,
            formatC(round(area_km2_yr), big.mark = ",", format = "d"))
  else
    sprintf("Projected invaded area: %s km²  |  %d yr forecast",
            formatC(round(area_km2_yr), big.mark = ",", format = "d"),
            yr - OBS_END)

  # Raster → data frame
  suit_df <- rast_to_df(suit_crop, "suitability")
  inv_df  <- rast_to_df(inv_crop,  "invaded")
  inv_df  <- inv_df[!is.na(inv_df$invaded) & inv_df$invaded == 1, ]

  # Detection points up to current year (for observed phase only)
  det_pts <- eab |> filter(year <= min(yr, OBS_END))

  fill_col <- if (is_observed) alpha(obs_col, 0.60) else alpha(proj_col, 0.65)

  # Build ggplot
  p <- ggplot() +
    # Suitability background
    geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradientn(
      colors   = suit_pal,
      na.value = "#0d1117",
      name     = "Suitability",
      limits   = c(0, 1),
      guide    = guide_colorbar(barheight = 8, barwidth = 0.8)
    ) +
    # Invaded / spread overlay — green for observed years, red for projected
    { if (nrow(inv_df) > 0)
        geom_raster(data = inv_df, aes(x = x, y = y),
                    fill = fill_col, inherit.aes = FALSE) } +
    # Country / state outlines
    { if (!is.null(canada))
        geom_sf(data = canada, fill = "#161b22", color = "#30363d",
                linewidth = 0.2, inherit.aes = FALSE) } +
    { if (!is.null(us_states))
        geom_sf(data = us_states, fill = NA, color = "#30363d",
                linewidth = 0.25, inherit.aes = FALSE) } +
    # Detection points — show cumulative up to current year (observed phase)
    # In projection phase show all confirmed 2001–2025 points
    { if (nrow(det_pts) > 0)
        geom_point(data = det_pts, aes(x = lon, y = lat),
                   color = "#00d4aa", size = 0.7, alpha = 0.55,
                   inherit.aes = FALSE) } +
    coord_sf(xlim = c(-130, -65), ylim = c(25, 50), expand = FALSE) +
    labs(
      title    = year_label,
      subtitle = subtitle_txt,
      caption  = paste0(
        "SDM: MaxEnt (maxnet, 2013–2025) | Spread: stratified diffusion + jump dispersal\n",
        "Background: ash-weighted FIA plots | Climate: BIO6, BIO4, BIO12 | Full CONUS extent"
      ),
      x = NULL, y = NULL
    ) +
    theme_void() +
    theme(
      plot.background   = element_rect(fill = "#0d1117", color = NA),
      panel.background  = element_rect(fill = "#0d1117"),
      plot.title        = element_text(color = "#e6edf3", face = "bold",
                                       size = 16, hjust = 0.5, margin = margin(t=8)),
      plot.subtitle     = element_text(color = "#8b949e", size = 9,
                                       hjust = 0.5, margin = margin(b=4)),
      plot.caption      = element_text(color = "#4a5568", size = 7,
                                       hjust = 0.5, margin = margin(t=6, b=4)),
      legend.text       = element_text(color = "#8b949e", size = 8),
      legend.title      = element_text(color = "#e6edf3", size = 9, face = "bold"),
      legend.background = element_rect(fill = "#161b22", color = NA),
      legend.position   = "right"
    )

  out_path <- file.path(OUTPUT_DIR, sprintf("spread_%d.png", yr))
  ggsave(out_path, p, width = 10, height = 6, dpi = 120, bg = "#0d1117")
  png_paths[i] <- out_path

  if (i %% 5 == 0) cat(sprintf("  Saved %d / %d maps\n", i, length(sim_years)))
}

# ── 7d. Stitch into animated GIF ---------------------------------------------
cat("Rendering GIF…\n")

# Frame durations — find positions by value since sim_years may skip early years
frame_dur <- rep(0.18, length(sim_years))
frame_dur[sim_years <= OBS_END]                    <- 0.25   # observed phase
frame_dur[1]                                        <- 0.8    # hold on first frame
obs_end_pos <- match(OBS_END, sim_years)
if (!is.na(obs_end_pos)) frame_dur[obs_end_pos]    <- 1.2    # pause at 2025 boundary
frame_dur[length(sim_years)]                        <- 2.0    # hold on 2040

gifski::gifski(
  png_files  = png_paths,
  gif_file   = file.path(OUTPUT_DIR, "eab_spread_animation.gif"),
  width      = 1200,
  height     = 720,
  delay      = frame_dur,
  loop       = TRUE,
  progress   = TRUE
)
cat("GIF saved to:", file.path(OUTPUT_DIR, "eab_spread_animation.gif"), "\n")


# =============================================================================
# 8. SUMMARY STATISTICS TABLE
# =============================================================================

summary_tbl <- data.frame(
  year = all_sim_years,
  invaded_cells = sapply(annual_invaded, sum, na.rm = TRUE),
  invaded_km2   = sapply(annual_invaded, sum, na.rm = TRUE) *
                    (cell_res_deg * km_per_deg)^2
) |>
  mutate(
    status        = ifelse(year <= OBS_END, "Observed", "Projected"),
    new_cells_yr  = c(0, diff(invaded_cells)),
    pct_of_conus  = round(invaded_cells / length(suitable_cells) * 100, 1)
  )

write.csv(summary_tbl,
          file.path(OUTPUT_DIR, "spread_summary.csv"),
          row.names = FALSE)

cat("\n── Spread Summary (selected years) ──\n")
print(summary_tbl[summary_tbl$year %in% c(2001, 2010, 2020, OBS_END, 2030, 2035, 2040), ],
      row.names = FALSE)

cat("\n── Session info ──\n")
cat(sprintf("R %s | terra %s | maxnet %s\n",
            R.Version()$version.string,
            packageVersion("terra"),
            packageVersion("maxnet")))

cat("\nDone. All outputs written to:", OUTPUT_DIR, "\n")

# =============================================================================
# END OF SCRIPT
# =============================================================================
