# =============================================================================
# EAB Economic Loss Export — Low / Mid / High Scenarios
# =============================================================================
# Run AFTER eab_spread_model.R IN THE SAME R SESSION.
#
# Expected folder structure:
#
#   eab_model/
#   ├── eab_spread_model.R
#   ├── export_webmap.R
#   ├── export_economic_loss.R      ← this script
#   ├── data/
#   │   ├── insects_clean.csv
#   │   ├── ash_ba.csv
#   │   ├── bio6_min_temp.tif
#   │   ├── bio4_temp_seasonality.tif
#   │   ├── bio12_annual_precip.tif
#   │   └── fia/                    ← FIA downloads cached here (auto-created)
#   │       ├── MI_CSV.zip
#   │       ├── OH_CSV.zip
#   │       └── ...
#   ├── outputs/
#   │   ├── eab_suitability.tif
#   │   ├── spread_YYYY.png
#   │   └── eab_spread_animation.gif
#   └── webmap/
#       ├── index.html
#       └── data/
#           ├── manifest.js
#           ├── suitability.js
#           ├── detections.js
#           ├── spread_YYYY.js
#           └── economic_loss.js    ← output of this script
#
# Price source: Appalachian Hardwood Sawtimber Pricing by Species
#   NW Georgia, E Tennessee, W West Virginia (Doyle scale)
#   Low = $165/MBF  |  Mid = $400/MBF  |  High = $700/MBF
#
# install.packages(c("rFIA", "terra", "dplyr", "tidyr", "jsonlite"))
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(jsonlite)

# =============================================================================
# CONFIGURATION — edit these before running
# =============================================================================

# ── Stumpage price scenarios ($/MBF Doyle) ────────────────────────────────────
SCENARIOS <- list(
  low  = list(label = "Low",  price_mbf = 165),
  mid  = list(label = "Mid",  price_mbf = 400),
  high = list(label = "High", price_mbf = 700)
)

# Doyle scale underestimates actual cubic volume, especially on smaller logs.
# ~0.65 is a reasonable recovery factor for mixed ash stands.
DOYLE_RECOVERY <- 0.65

# ── States to pull FIA data for ───────────────────────────────────────────────
FIA_STATES <- c("MI", "OH", "IN", "IL", "WI", "MN", "IA", "MO",
                "PA", "NY", "WV", "VA", "KY", "TN", "NC", "GA",
                "AL", "MS", "AR", "LA", "CT", "MA", "VT", "NH",
                "ME", "NJ", "DE", "MD", "RI")

# ── Directories ───────────────────────────────────────────────────────────────
FIA_DIR    <- file.path(SCRIPT_DIR, "data", "fia")
WEBMAP_DIR <- file.path(SCRIPT_DIR, "webmap", "data")
dir.create(FIA_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(WEBMAP_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. DOWNLOAD FIA DATA
# =============================================================================
# FIA state zips are large (100MB–600MB each). The default R timeout of 60s
# is too short. We set 1 hour and download one state at a time so a single
# failure doesn't abort everything. Already-downloaded states are skipped.
#
# If a state keeps failing, download the zip manually from:
#   https://apps.fs.usda.gov/fia/datamart/CSV/{ST}_CSV.zip
# and place it in:  data/fia/{ST}_CSV.zip
# getFIA() will detect it and skip the download.

cat("1/4  Downloading FIA data...\n")
cat("     Large files — this may take 10–30 minutes on first run.\n")
cat("     Already-downloaded states will be skipped.\n\n")

library(rFIA)

# Increase timeout to 1 hour per file
options(timeout = 3600)

# Download one state at a time — skips states already in FIA_DIR
failed_states <- character(0)

for (st in FIA_STATES) {
  zip_path <- file.path(FIA_DIR, paste0(st, "_CSV.zip"))

  # Check if a complete zip already exists (> 1MB as a sanity check)
  if (file.exists(zip_path) && file.size(zip_path) > 1e6) {
    cat(sprintf("  %-3s  already downloaded — skipping\n", st))
    next
  }

  # Remove any incomplete/corrupt partial download before retrying
  if (file.exists(zip_path)) {
    cat(sprintf("  %-3s  incomplete download found — removing and retrying\n", st))
    file.remove(zip_path)
  }

  cat(sprintf("  %-3s  downloading...\n", st))
  tryCatch({
    getFIA(states = st, dir = FIA_DIR, load = FALSE)
    cat(sprintf("  %-3s  done (%.0f MB)\n", st,
                file.size(zip_path) / 1e6))
  }, error = function(e) {
    cat(sprintf("  %-3s  FAILED: %s\n", st, conditionMessage(e)))
    failed_states <<- c(failed_states, st)
  })
}

if (length(failed_states) > 0) {
  cat(sprintf(
    "\nWarning: %d state(s) failed to download: %s\n",
    length(failed_states), paste(failed_states, collapse = ", ")
  ))
  cat("Download these manually from:\n")
  cat("  https://apps.fs.usda.gov/fia/datamart/CSV/{ST}_CSV.zip\n")
  cat(sprintf("and place in: %s\n\n", FIA_DIR))
}

# Load all successfully downloaded states
cat("\nLoading FIA data into memory...\n")
fia <- readFIA(FIA_DIR, common = TRUE)
cat(sprintf("  Loaded FIA data for states in: %s\n", FIA_DIR))

# =============================================================================
# 2. CALCULATE PLOT-LEVEL ASH VOLUME
# =============================================================================
# Ash species codes:
#   541 = white ash  (Fraxinus americana)
#   543 = black ash  (Fraxinus nigra)
#   544 = green ash  (Fraxinus pennsylvanica)

#
# VOLCFNET  = net merchantable cubic foot volume per tree (taper-corrected)
# TPA_UNADJ = trees per acre expansion factor
# STATUSCD  = 1 for live trees only

cat("2/4  Calculating plot-level ash volume...\n")

ASH_SPCD <- c(541, 543, 544)

ash_vol_plot <- fia$TREE |>
  filter(
    SPCD     %in% ASH_SPCD,
    STATUSCD == 1,
    !is.na(VOLCFNET),
    VOLCFNET > 0
  ) |>
  group_by(PLT_CN) |>
  summarise(
    vol_ft3_per_acre = sum(VOLCFNET * TPA_UNADJ, na.rm = TRUE),
    n_trees          = n(),
    .groups          = "drop"
  ) |>
  # ft³/acre → m³/ha
  # 1 ft³ = 0.0283168 m³
  # 1 acre = 0.404686 ha  →  1 ft³/acre = 0.0283168 / 0.404686 m³/ha
  mutate(vol_m3_per_ha = vol_ft3_per_acre * (0.0283168 / 0.404686)) |>
  left_join(
    fia$PLOT |>
      select(PLT_CN, LAT, LON, INVYR) |>
      filter(!is.na(LAT), !is.na(LON)),
    by = "PLT_CN"
  ) |>
  filter(!is.na(LAT), !is.na(LON)) |>
  # Keep most recent plot visit only
  group_by(PLT_CN) |>
  slice_max(INVYR, n = 1, with_ties = FALSE) |>
  ungroup()

cat(sprintf("  Ash plots with volume data: %d\n",   nrow(ash_vol_plot)))
cat(sprintf("  Volume range: %.1f – %.1f m³/ha\n",
            min(ash_vol_plot$vol_m3_per_ha),
            max(ash_vol_plot$vol_m3_per_ha)))

# =============================================================================
# 3. RASTERIZE ONTO SIMULATION GRID
# =============================================================================
# suit_sim, cell_res_deg, km_per_deg all inherited from eab_spread_model.R

cat("3/4  Rasterizing ash volume onto simulation grid...\n")

ash_sv   <- vect(ash_vol_plot, geom = c("LON", "LAT"), crs = "EPSG:4326")
vol_rast <- rasterize(ash_sv, suit_sim, field = "vol_m3_per_ha",
                      fun = "mean", na.rm = TRUE)
vol_vals <- values(vol_rast)[, 1]
vol_vals[is.na(vol_vals)] <- 0

# Cell area in hectares
cell_area_ha <- (cell_res_deg * km_per_deg * 100)^2 / 10000

cat(sprintf("  Cell area:              %.1f ha\n",          cell_area_ha))
cat(sprintf("  Cells with ash volume:  %d\n",               sum(vol_vals > 0)))
cat(sprintf("  Total ash volume:       %.1f million m³\n",
            sum(vol_vals * cell_area_ha, na.rm = TRUE) / 1e6))

# =============================================================================
# 4. CALCULATE ANNUAL LOSS — ALL THREE SCENARIOS
# =============================================================================
# Price conversion:
#   $/MBF Doyle → $/ft³:  divide by 1000
#   $/ft³ → $/m³:         divide by 0.0283168
#   Doyle recovery:        multiply by DOYLE_RECOVERY
#
# Loss is calculated only on NEWLY invaded cells each year to avoid
# double-counting timber that was already counted in a previous year.

cat("4/4  Calculating annual economic loss across scenarios...\n\n")

# Export years — same filter used in export_webmap.R
all_export_years <- as.integer(names(annual_invaded))
det_cumul <- sapply(all_export_years, function(yr)
  nrow(eab |> filter(year <= min(yr, OBS_END))))
export_years <- all_export_years[
  all_export_years >= SIM_START | det_cumul >= 10]

all_scenarios <- lapply(names(SCENARIOS), function(s) {

  price_mbf    <- SCENARIOS[[s]]$price_mbf
  price_usd_m3 <- (price_mbf / 1000) / 0.0283168 * DOYLE_RECOVERY

  value_per_cell <- vol_vals * cell_area_ha * price_usd_m3

  cat(sprintf("  %s scenario: $%d/MBF → $%.2f/m³ | total at risk: $%.1fB\n",
              SCENARIOS[[s]]$label, price_mbf, price_usd_m3,
              sum(value_per_cell, na.rm = TRUE) / 1e9))

  loss_rows <- lapply(seq_along(export_years), function(i) {
    yr       <- export_years[i]
    inv_curr <- annual_invaded[[as.character(yr)]]

    if (i == 1) {
      new_cells <- inv_curr
    } else {
      inv_prev  <- annual_invaded[[as.character(export_years[i - 1])]]
      new_cells <- inv_curr & !inv_prev
    }

    loss <- sum(value_per_cell[new_cells], na.rm = TRUE)

    list(
      year      = yr,
      scenario  = s,
      label     = SCENARIOS[[s]]$label,
      price_mbf = price_mbf,
      new_cells = sum(new_cells),
      loss_usd  = round(loss)
    )
  })

  bind_rows(loss_rows) |>
    mutate(cumulative_loss_usd = cumsum(loss_usd))
})

loss_df <- bind_rows(all_scenarios)

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n── Cumulative Loss by Year ──\n")
summary_years <- c(2003, 2010, 2015, 2020, 2025, 2030, 2035, 2040)

loss_wide <- loss_df |>
  filter(year %in% summary_years) |>
  mutate(cumul_B = paste0("$", round(cumulative_loss_usd / 1e9, 2), "B")) |>
  select(year, scenario, cumul_B) |>
  pivot_wider(names_from = scenario, values_from = cumul_B) |>
  select(year, low, mid, high)

print(loss_wide, row.names = FALSE)

cat("\n── Total Projected Loss 2003–2040 ──\n")
for (s in names(SCENARIOS)) {
  total <- sum(loss_df$loss_usd[loss_df$scenario == s], na.rm = TRUE) / 1e9
  cat(sprintf("  %-5s ($%d/MBF Doyle):  $%.2f billion\n",
              SCENARIOS[[s]]$label, SCENARIOS[[s]]$price_mbf, total))
}

# ── Write JS ──────────────────────────────────────────────────────────────────
out <- file.path(WEBMAP_DIR, "economic_loss.js")
writeLines(
  paste0("var ECONOMIC_LOSS = ", toJSON(loss_df, auto_unbox = TRUE), ";"),
  out
)
cat(sprintf("\n  economic_loss.js written → %s (%.0f KB)\n",
            out, file.size(out) / 1024))

cat("\nDone. Refresh the webmap to see the cumulative loss stat.\n")

# =============================================================================
# END
# =============================================================================
