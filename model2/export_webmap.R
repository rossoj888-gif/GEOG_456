# =============================================================================
# EAB Webmap Export
# =============================================================================
# Run AFTER eab_spread_model.R IN THE SAME R SESSION.
#
# Exports model results as .js files (GeoJSON wrapped in a var assignment)
# so index.html can load them with <script src> — no server needed.
#
# Output: webmap/
#   index.html
#   data/
#     manifest.js          var MANIFEST = {...};
#     suitability.js       var SUIT     = {FeatureCollection};
#     detections.js        var DETECTIONS = {FeatureCollection};
#     spread_YYYY.js       var SPREAD_YYYY = {FeatureCollection};
#
# install.packages(c("terra", "sf", "dplyr", "jsonlite"))
# =============================================================================

library(terra)
library(dplyr)
library(jsonlite)

WEBMAP_DIR <- file.path(SCRIPT_DIR, "webmap", "data")
dir.create(WEBMAP_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Exporting to:", WEBMAP_DIR, "\n\n")

# Helper: write a GeoJSON FeatureCollection as a JS variable file
write_js_var <- function(var_name, sf_obj, path) {
  geojson <- geojsonsf::sf_geojson(sf_obj)   # fast conversion
  writeLines(paste0("var ", var_name, " = ", geojson, ";"), path)
  cat(sprintf("  %-20s  %d features  %.1f MB\n",
              basename(path), nrow(sf_obj),
              file.size(path) / 1e6))
}

# =============================================================================
# 1. SUITABILITY
# =============================================================================
# 1. SUITABILITY
# =============================================================================
cat("1/3  Suitability...\n")

SUIT_STEP   <- 4       # finer grid — matches GIF resolution
SUIT_THRESH <- 0.30

suit_down <- aggregate(suit_rast_masked, fact = SUIT_STEP, fun = "mean", na.rm = TRUE)
suit_down[suit_down < SUIT_THRESH] <- NA

# Extract cell centroids and values as a plain data frame — never use sfg
suit_df   <- as.data.frame(suit_down, xy = TRUE, na.rm = TRUE)
names(suit_df) <- c("cx", "cy", "suit")
suit_df$suit   <- round(suit_df$suit, 3)
half_s         <- res(suit_down)[1] / 2   # half cell width in degrees

# Build GeoJSON manually from centroid + half-width — no sf, no sfg
suit_features <- lapply(seq_len(nrow(suit_df)), function(i) {
  cx <- suit_df$cx[i]; cy <- suit_df$cy[i]
  list(
    type       = "Feature",
    properties = list(suit = suit_df$suit[i]),
    geometry   = list(
      type        = "Polygon",
      coordinates = list(list(
        list(cx - half_s, cy - half_s),
        list(cx + half_s, cy - half_s),
        list(cx + half_s, cy + half_s),
        list(cx - half_s, cy + half_s),
        list(cx - half_s, cy - half_s)
      ))
    )
  )
})

suit_list <- list(type = "FeatureCollection", features = suit_features)

out <- file.path(WEBMAP_DIR, "suitability.js")
writeLines(paste0("var SUIT = ", toJSON(suit_list, auto_unbox = TRUE), ";"), out)
cat(sprintf("  suitability.js  %d cells  %.1f MB\n",
            nrow(suit_df), file.size(out) / 1e6))

# =============================================================================
# 2. DETECTIONS
# =============================================================================
cat("2/3  Detections...\n")

det_list <- list(
  type = "FeatureCollection",
  features = lapply(seq_len(nrow(eab)), function(i)
    list(type = "Feature",
         properties = list(year = eab$year[i]),
         geometry = list(type = "Point",
                         coordinates = list(eab$lon[i], eab$lat[i]))))
)

out <- file.path(WEBMAP_DIR, "detections.js")
writeLines(paste0("var DETECTIONS = ", toJSON(det_list, auto_unbox = TRUE), ";"), out)
cat(sprintf("  detections.js  %d points  %.1f MB\n",
            nrow(eab), file.size(out)/1e6))

# =============================================================================
# 3. PER-YEAR SPREAD
# =============================================================================
cat("3/3  Per-year spread...\n")

all_export_years <- as.integer(names(annual_invaded))
det_cumul <- sapply(all_export_years, function(yr)
  nrow(eab |> filter(year <= min(yr, OBS_END))))
export_years <- all_export_years[
  all_export_years >= SIM_START | det_cumul >= 10]

for (yr in export_years) {
  inv_vec <- annual_invaded[[as.character(yr)]]
  idx     <- which(inv_vec)
  if (length(idx) == 0) next

  xy <- xyFromCell(suit_sim, idx)

  # Each invaded cell → a GeoJSON polygon
  half <- cell_res_deg / 2
  feat_list <- list(
    type = "FeatureCollection",
    features = lapply(seq_along(idx), function(i) {
      cx <- xy[i, 1]; cy <- xy[i, 2]
      list(
        type = "Feature",
        properties = list(year = yr),
        geometry = list(
          type = "Polygon",
          coordinates = list(list(
            list(cx-half, cy-half), list(cx+half, cy-half),
            list(cx+half, cy+half), list(cx-half, cy+half),
            list(cx-half, cy-half)
          ))
        )
      )
    })
  )

  out <- file.path(WEBMAP_DIR, sprintf("spread_%d.js", yr))
  writeLines(
    paste0("var SPREAD_", yr, " = ", toJSON(feat_list, auto_unbox = TRUE), ";"),
    out
  )
  cat(sprintf("  spread_%d.js  %d cells  %.0f KB\n",
              yr, length(idx), file.size(out)/1024))
}

# =============================================================================
# 4. MANIFEST
# =============================================================================
manifest <- list(
  obs_end      = OBS_END,
  sim_start    = SIM_START,
  sim_end      = SIM_END,
  export_years = export_years
)
out <- file.path(WEBMAP_DIR, "manifest.js")
writeLines(paste0("var MANIFEST = ", toJSON(manifest, auto_unbox = TRUE), ";"), out)
cat("  manifest.js written\n")

# =============================================================================
# 5. COPY index.html
# =============================================================================
html_src <- file.path(SCRIPT_DIR, "index.html")
html_dst <- file.path(SCRIPT_DIR, "webmap", "index.html")

if (file.exists(html_src)) {
  # Inject <script> tags for every spread year
  script_tags <- paste(
    sapply(export_years, function(yr)
      sprintf('    <script src="./data/spread_%d.js"></script>', yr)),
    collapse = "\n"
  )
  html <- readLines(html_src, warn = FALSE)
  html <- gsub("<!-- SPREAD_SCRIPTS -->", script_tags, html, fixed = TRUE)
  writeLines(html, html_dst)
  cat("  index.html written\n")
} else {
  message("index.html not found at ", html_src)
}

cat("\n------------------------------------------------------------\n")
cat("Done. Open directly in browser (no server needed):\n")
cat("  ", html_dst, "\n")
cat("------------------------------------------------------------\n")
