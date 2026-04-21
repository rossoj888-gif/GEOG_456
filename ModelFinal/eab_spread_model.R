# =============================================================================
# Emerald Ash Borer (EAB) — Species Distribution Model + Spread Simulation
# Author: Jack Ross | R 4.3+
#
# Inputs  (data/):  insects_clean.csv | ash_ba.csv | bio6/bio4/bio12 .tif
# Outputs (outputs/): eab_suitability.tif | spread_YYYY.png | animation.gif
#
# Run order: this script → export_webmap.R (same session)
# =============================================================================

# ── 0. Setup ------------------------------------------------------------------

library(terra)
library(maxnet)
library(sf)
library(dplyr)
library(ggplot2)
library(pROC)
library(gifski)
library(scales)
library(viridis)

# Set SCRIPT_DIR to the eab_model folder, e.g.:
#   Windows:  SCRIPT_DIR <- "C:/Users/Jack/Documents/eab_model"
#   Mac/Linux: SCRIPT_DIR <- "~/Documents/eab_model"
SCRIPT_DIR <- getwd()

DATA_DIR   <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("Data directory:", DATA_DIR, "\n")
stopifnot("data/ folder not found — update SCRIPT_DIR" = dir.exists(DATA_DIR))

set.seed(42)

# Convert SpatRaster to ggplot2-ready data frame
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

# SDM trained on 2013+ only — pre-2013 has 47 points clustered near Detroit,
# too sparse and geographically biased to characterise the climate niche.
SDM_YEAR_MIN <- 2013
eab_sdm <- eab_all |> filter(year >= SDM_YEAR_MIN, year <= 2025)
eab     <- eab_all |> filter(year <= 2025)

cat(sprintf("EAB full record: %d points (%d–%d)\n",
            nrow(eab), min(eab$year), max(eab$year)))
cat(sprintf("EAB SDM subset:  %d points (%d–2025)\n",
            nrow(eab_sdm), SDM_YEAR_MIN))

# ── 1b. Ash basal area — used for background weighting and ash mask -----------
ash <- read.csv(file.path(DATA_DIR, "ash_ba.csv")) |>
  filter(!is.na(LAT), !is.na(LON), BAA > 0) |>
  rename(lat = LAT, lon = LON, baa = BAA)

cat(sprintf("Ash FIA plots (BAA > 0): %d\n", nrow(ash)))

# ── 1c. Climate rasters — VIF-filtered to BIO6, BIO4, BIO12 ------------------
clim <- c(
  rast(file.path(DATA_DIR, "bio6_min_temp.tif")),
  rast(file.path(DATA_DIR, "bio4_temp_seasonality.tif")),
  rast(file.path(DATA_DIR, "bio12_annual_precip.tif"))
)
names(clim) <- c("bio6_min_temp", "bio4_temp_seasonality", "bio12_annual_precip")

cat(sprintf("Climate: %d layers | %s | %.4f°\n",
            nlyr(clim), crs(clim, describe = TRUE)$code, res(clim)[1]))


# =============================================================================
# 2. PREPARE PRESENCE / BACKGROUND DATA
# =============================================================================

# ── 2a. SDM occurrences to SpatVector ----------------------------------------
eab_sv <- vect(eab_sdm, geom = c("lon", "lat"), crs = "EPSG:4326")

# ── 2b. Thin to 1 point per raster cell to reduce spatial autocorrelation ----
eab_cells  <- cellFromXY(clim[[1]], geom(eab_sv)[, c("x","y")])
eab_unique <- eab_sdm[!duplicated(eab_cells), ]
eab_sv     <- vect(eab_unique, geom = c("lon", "lat"), crs = "EPSG:4326")
cat(sprintf("After thinning: %d SDM presence cells\n", nrow(eab_unique)))

# ── 2c. Ash-weighted background (n = 10,000) ----------------------------------
# Sampling proportional to basal area restricts the background to ash-occupied
# environments — avoids conflating unsuitable climate with absence of host.
n_bg       <- 10000
bg_weights <- ash$baa / sum(ash$baa)
bg_idx     <- sample(nrow(ash), size = n_bg, replace = TRUE, prob = bg_weights)
bg         <- ash[bg_idx, c("lat", "lon")] |>
  mutate(presence = 0L, year = NA_integer_)

occ_bg <- bind_rows(
  eab_unique |> select(lat, lon, presence, year),
  bg
) |> filter(!is.na(lat), !is.na(lon))

# ── 2d. Extract climate values at all points ----------------------------------
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

# lqh: linear + quadratic + hinge features — standard for moderate sample sizes
# regmult 1.5: slightly above default to reduce overfitting from clustered detections
sdm <- maxnet(
  p       = y,
  data    = as.data.frame(X),
  f       = maxnet.formula(y, as.data.frame(X), classes = "lqh"),
  regmult = 1.5
)

cat("MaxEnt fitted | non-zero coefficients:", sum(sdm$betas != 0), "\n")


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
  pred_k    <- predict(m_k, as.data.frame(Xte), type = "cloglog")
  cv_auc[k] <- as.numeric(auc(roc(yte, as.numeric(pred_k), quiet = TRUE)))
}

cat(sprintf("5-fold CV AUC: %.3f ± %.3f\n", mean(cv_auc), sd(cv_auc)))


# =============================================================================
# 4. PROJECT SUITABILITY — masked to ash host footprint
# =============================================================================
# Cells outside the FIA ash footprint are set to NA (no host = no suitability).
# The ash mask also defines the natural western study area boundary.

suit_rast <- predict(clim, sdm, type = "cloglog", na.rm = TRUE)
names(suit_rast) <- "suitability"

ash_sv   <- vect(ash, geom = c("lon", "lat"), crs = "EPSG:4326")
ash_rast <- rasterize(ash_sv, clim[[1]], field = "baa", fun = "mean")
ash_mask <- !is.na(ash_rast)
suit_rast_masked <- mask(suit_rast, ash_mask, maskvalue = FALSE)

writeRaster(suit_rast_masked,
            file.path(OUTPUT_DIR, "eab_suitability.tif"),
            overwrite = TRUE, datatype = "FLT4S")
cat("Suitability raster saved.\n")


# =============================================================================
# 5. MODEL DIAGNOSTICS — marginal response curves
# =============================================================================

# Each variable varied across its 2nd–98th percentile, others held at mean
response_data <- lapply(CLIM_VARS, function(var) {
  x_range <- seq(
    quantile(model_df[[var]], 0.02, na.rm = TRUE),
    quantile(model_df[[var]], 0.98, na.rm = TRUE),
    length.out = 200
  )
  newdat <- as.data.frame(
    lapply(CLIM_VARS, function(v) {
      if (v == var) x_range else rep(mean(model_df[[v]], na.rm = TRUE), 200)
    })
  )
  names(newdat) <- CLIM_VARS
  data.frame(variable = var, x = x_range,
             suit = as.numeric(predict(sdm, newdat, type = "cloglog")))
}) |> bind_rows()

var_labels <- c(
  bio6_min_temp         = "BIO6: Min Temp Coldest Month (°C)",
  bio4_temp_seasonality = "BIO4: Temp Seasonality",
  bio12_annual_precip   = "BIO12: Annual Precipitation (mm)"
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
cat("Diagnostics saved.\n")


# =============================================================================
# 6. STRATIFIED DIFFUSION SPREAD SIMULATION
# =============================================================================
# Each annual time step applies:
#   A) Gaussian local diffusion  — natural beetle flight
#   B) Long-distance jump dispersal — human-mediated firewood/nursery transport
#   C) Suitability filter — cells below SUIT_THRESHOLD cannot establish

# ── 6a. Parameters -----------------------------------------------------------
OBS_START      <- 2001   # first year of detections
OBS_END        <- 2025   # last confirmed detection year
SIM_START      <- 2026   # first projected year
SIM_END        <- 2040   # forecast horizon
SUIT_THRESHOLD <- 0.30   # min suitability for establishment (Valladares et al. 2012)
LOCAL_SIGMA_KM <- 10     # Gaussian diffusion σ — calibrated to observed front advance
JUMP_RATE      <- 0.008  # fraction of invaded cells generating a jump per year
JUMP_MEAN_KM   <- 93     # mean jump distance — Ward et al. (2020) 93 ± 7 SE km
N_JUMP_MAX     <- 200    # annual jump cap — calibrated to observed footprint

# ── 6b. Resample to ~0.125° (~14 km) for computational efficiency ------------
suit_sim <- aggregate(suit_rast_masked, fact = 3, fun = "mean", na.rm = TRUE)
cat(sprintf("Simulation grid: %d × %d cells (%.4f°)\n",
            nrow(suit_sim), ncol(suit_sim), res(suit_sim)[1]))

suit_vals    <- values(suit_sim)[, 1]
n_cells      <- ncell(suit_sim)
all_xy       <- xyFromCell(suit_sim, 1:n_cells)
cell_res_deg <- res(suit_sim)[1]
km_per_deg   <- 111.0

local_sigma_d <- LOCAL_SIGMA_KM / km_per_deg
jump_mean_d   <- JUMP_MEAN_KM   / km_per_deg

# ── 6c. Gaussian dispersal kernel --------------------------------------------
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

# Snap lat/lon points to simulation grid cells
pts_to_cells <- function(df) {
  cells <- cellFromXY(suit_sim, cbind(df$lon, df$lat))
  cells[!is.na(cells)]
}

# =============================================================================
# PHASE 1: Observed invasion history (2001–2025)
# =============================================================================
# Directly reconstruct from detection records — no simulation.
# The 2025 state is locked as the baseline for Phase 2.

cat("Phase 1: observed history 2001–2025…\n")

annual_invaded <- list()
invaded        <- rep(FALSE, n_cells)

for (yr in OBS_START:OBS_END) {
  new_cells <- pts_to_cells(eab |> filter(year == yr))
  if (length(new_cells) > 0) invaded[new_cells] <- TRUE
  annual_invaded[[as.character(yr)]] <- invaded
}

observed_invaded <- invaded   # locked baseline — never removed in Phase 2

cat(sprintf("  2001: %d cells | 2025: %d cells\n",
            sum(annual_invaded[["2001"]]), sum(invaded)))

# =============================================================================
# PHASE 2: Forward projection (2026–2040)
# =============================================================================

cat("Phase 2: projecting 2026–2040…\n")
pb <- txtProgressBar(min = SIM_START, max = SIM_END, style = 3)

for (yr in SIM_START:SIM_END) {

  # ── A) Local diffusion ──────────────────────────────────────────────────────
  spread_pressure <- convolve_invaded(invaded, suit_sim, gauss_kernel)

  new_local <- !invaded &
               !is.na(suit_vals) &
               suit_vals >= SUIT_THRESHOLD &
               spread_pressure > 0.05

  # ── B) Long-distance jumps (exponential distances, random directions) ───────
  invaded_idx <- which(invaded)
  n_jumps     <- min(round(length(invaded_idx) * JUMP_RATE), N_JUMP_MAX)

  new_jump <- rep(FALSE, n_cells)
  if (n_jumps > 0 && length(invaded_idx) > 0) {
    src_cells  <- sample(invaded_idx, n_jumps, replace = TRUE)
    src_xy     <- all_xy[src_cells, , drop = FALSE]
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

  # ── C) Update — re-lock observed baseline each step ─────────────────────────
  invaded <- invaded | new_local | new_jump | observed_invaded
  annual_invaded[[as.character(yr)]] <- invaded
  setTxtProgressBar(pb, yr)
}
close(pb)

n_inv_final <- sum(annual_invaded[[as.character(SIM_END)]], na.rm = TRUE)
area_km2    <- n_inv_final * (cell_res_deg * km_per_deg)^2
cat(sprintf("\nProjected invaded cells by %d: %d (~%.0f km²)\n",
            SIM_END, n_inv_final, area_km2))


# =============================================================================
# 7. ANNUAL SPREAD MAPS + ANIMATED GIF
# =============================================================================

# ── 7a. Filter years — skip observed years with < 10 cumulative detections ----
all_sim_years <- as.integer(names(annual_invaded))

det_by_year <- sapply(all_sim_years, function(yr)
  nrow(eab |> filter(year <= min(yr, OBS_END)))
)

sim_years <- all_sim_years[
  all_sim_years >= SIM_START |
  det_by_year   >= 10
]

cat(sprintf("Rendering %d years\n", length(sim_years)))

# ── 7b. Colour palette --------------------------------------------------------
suit_pal <- colorRampPalette(c("#0d1117","#1a3a4a","#0d6e8a",
                                "#00b4a0","#7dde82","#ffd60a",
                                "#ff6b35","#c9184a"))(100)
obs_col  <- "#3fb950"
proj_col <- "#e63946"

# ── 7c. Generate one PNG per year ---------------------------------------------
png_paths <- character(length(sim_years))

if (requireNamespace("maps", quietly = TRUE)) {
  us_states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
  us_nation <- sf::st_as_sf(maps::map("usa",   plot = FALSE, fill = TRUE))
  canada    <- sf::st_as_sf(maps::map("world", regions = "Canada",
                                       plot = FALSE, fill = TRUE))
} else {
  us_states <- us_nation <- canada <- NULL
  message("Install the 'maps' package for boundary overlays.")
}

cat("Generating maps…\n")
for (i in seq_along(sim_years)) {
  yr <- sim_years[i]

  inv_r         <- suit_sim
  values(inv_r) <- as.numeric(annual_invaded[[as.character(yr)]])
  inv_r[inv_r == 0] <- NA

  ext_plot  <- ext(-130, -65, 25, 50)
  suit_crop <- crop(suit_rast_masked, ext_plot)
  inv_crop  <- crop(inv_r, ext_plot)

  is_observed  <- yr <= OBS_END
  year_label   <- if (is_observed) sprintf("EAB Observed Invasion — %d", yr) else
                                   sprintf("EAB Projected Spread — %d", yr)
  n_det_yr     <- nrow(eab |> filter(year <= yr))
  area_km2_yr  <- sum(annual_invaded[[as.character(yr)]], na.rm = TRUE) *
                    (cell_res_deg * km_per_deg)^2
  subtitle_txt <- if (is_observed)
    sprintf("Cumulative detections: %d  |  Invaded area: %s km²",
            n_det_yr, formatC(round(area_km2_yr), big.mark = ",", format = "d"))
  else
    sprintf("Projected invaded area: %s km²  |  %d yr forecast",
            formatC(round(area_km2_yr), big.mark = ",", format = "d"), yr - OBS_END)

  suit_df  <- rast_to_df(suit_crop, "suitability")
  inv_df   <- rast_to_df(inv_crop, "invaded")
  inv_df   <- inv_df[!is.na(inv_df$invaded) & inv_df$invaded == 1, ]
  det_pts  <- eab |> filter(year <= min(yr, OBS_END))
  fill_col <- if (is_observed) alpha(obs_col, 0.60) else alpha(proj_col, 0.65)

  p <- ggplot() +
    geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradientn(colors = suit_pal, na.value = "#0d1117",
                         name = "Suitability", limits = c(0, 1),
                         guide = guide_colorbar(barheight = 8, barwidth = 0.8)) +
    { if (nrow(inv_df) > 0)
        geom_raster(data = inv_df, aes(x = x, y = y),
                    fill = fill_col, inherit.aes = FALSE) } +
    { if (!is.null(canada))
        geom_sf(data = canada, fill = "#161b22", color = "#30363d",
                linewidth = 0.2, inherit.aes = FALSE) } +
    { if (!is.null(us_states))
        geom_sf(data = us_states, fill = NA, color = "#30363d",
                linewidth = 0.25, inherit.aes = FALSE) } +
    { if (nrow(det_pts) > 0)
        geom_point(data = det_pts, aes(x = lon, y = lat),
                   color = "#00d4aa", size = 0.7, alpha = 0.55,
                   inherit.aes = FALSE) } +
    coord_sf(xlim = c(-130, -65), ylim = c(25, 50), expand = FALSE) +
    labs(title = year_label, subtitle = subtitle_txt,
         caption = paste0(
           "SDM: MaxEnt (maxnet, 2013–2025) | Spread: stratified diffusion + jump dispersal\n",
           "Background: ash-weighted FIA | Climate: BIO6, BIO4, BIO12 | CONUS extent"),
         x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "#0d1117", color = NA),
      panel.background = element_rect(fill = "#0d1117"),
      plot.title       = element_text(color = "#e6edf3", face = "bold",
                                      size = 16, hjust = 0.5, margin = margin(t=8)),
      plot.subtitle    = element_text(color = "#8b949e", size = 9,
                                      hjust = 0.5, margin = margin(b=4)),
      plot.caption     = element_text(color = "#4a5568", size = 7,
                                      hjust = 0.5, margin = margin(t=6, b=4)),
      legend.text      = element_text(color = "#8b949e", size = 8),
      legend.title     = element_text(color = "#e6edf3", size = 9, face = "bold"),
      legend.background = element_rect(fill = "#161b22", color = NA),
      legend.position  = "right"
    )

  out_path    <- file.path(OUTPUT_DIR, sprintf("spread_%d.png", yr))
  ggsave(out_path, p, width = 10, height = 6, dpi = 120, bg = "#0d1117")
  png_paths[i] <- out_path

  if (i %% 5 == 0) cat(sprintf("  %d / %d\n", i, length(sim_years)))
}

# ── 7d. Stitch into GIF — variable frame timing ------------------------------
cat("Rendering GIF…\n")

frame_dur <- rep(0.18, length(sim_years))
frame_dur[sim_years <= OBS_END] <- 0.25       # slower through observed phase
frame_dur[1]                    <- 0.8         # hold first frame
obs_end_pos <- match(OBS_END, sim_years)
if (!is.na(obs_end_pos)) frame_dur[obs_end_pos] <- 1.2   # pause at 2025
frame_dur[length(sim_years)]    <- 2.0         # hold final frame

gifski::gifski(
  png_files = png_paths,
  gif_file  = file.path(OUTPUT_DIR, "eab_spread_animation.gif"),
  width     = 1200, height = 720,
  delay     = frame_dur, loop = TRUE, progress = TRUE
)
cat("GIF saved.\n")


# =============================================================================
# 8. SUMMARY STATISTICS TABLE
# =============================================================================

summary_tbl <- data.frame(
  year          = all_sim_years,
  invaded_cells = sapply(annual_invaded, sum, na.rm = TRUE),
  invaded_km2   = sapply(annual_invaded, sum, na.rm = TRUE) *
                    (cell_res_deg * km_per_deg)^2
) |>
  mutate(
    status       = ifelse(year <= OBS_END, "Observed", "Projected"),
    new_cells_yr = c(0, diff(invaded_cells)),
    pct_of_conus = round(invaded_cells / length(suitable_cells) * 100, 1)
  )

write.csv(summary_tbl, file.path(OUTPUT_DIR, "spread_summary.csv"), row.names = FALSE)

cat("\n── Spread summary (selected years) ──\n")
print(summary_tbl[summary_tbl$year %in% c(2001, 2010, 2020, OBS_END, 2030, 2035, 2040), ],
      row.names = FALSE)

cat(sprintf("\nR %s | terra %s | maxnet %s\n",
            R.Version()$version.string,
            packageVersion("terra"),
            packageVersion("maxnet")))

cat("\nDone. Outputs:", OUTPUT_DIR, "\n")
