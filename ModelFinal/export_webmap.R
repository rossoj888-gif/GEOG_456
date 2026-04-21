# =============================================================================
# EAB Webmap Export
# Run AFTER eab_spread_model.R IN THE SAME R SESSION.
#
# Wraps GeoJSON in JS variable assignments so index.html can load them via
# <script src> without a backend — only a local HTTP server is needed.
#
# Output: webmap/data/
#   manifest.js   suitability.js   detections.js
#   spread_YYYY.js (one per year)   isolines.js
# =============================================================================

library(terra)
library(dplyr)
library(jsonlite)
library(smoothr)

WEBMAP_DIR <- file.path(SCRIPT_DIR, "webmap", "data")
dir.create(WEBMAP_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Exporting to:", WEBMAP_DIR, "\n\n")


# =============================================================================
# 1. SUITABILITY
# =============================================================================
cat("1/4  Suitability…\n")

SUIT_STEP   <- 4      # downsample factor — reduces file size, minimal visual loss
SUIT_THRESH <- 0.30   # drop cells below threshold (not part of potential range)

suit_down <- aggregate(suit_rast_masked, fact = SUIT_STEP, fun = "mean", na.rm = TRUE)
suit_down[suit_down < SUIT_THRESH] <- NA

suit_df        <- as.data.frame(suit_down, xy = TRUE, na.rm = TRUE)
names(suit_df) <- c("cx", "cy", "suit")
suit_df$suit   <- round(suit_df$suit, 3)
half_s         <- res(suit_down)[1] / 2

# Build GeoJSON polygons manually from cell centroids
suit_features <- lapply(seq_len(nrow(suit_df)), function(i) {
  cx <- suit_df$cx[i]; cy <- suit_df$cy[i]
  list(
    type       = "Feature",
    properties = list(suit = suit_df$suit[i]),
    geometry   = list(
      type        = "Polygon",
      coordinates = list(list(
        list(cx - half_s, cy - half_s), list(cx + half_s, cy - half_s),
        list(cx + half_s, cy + half_s), list(cx - half_s, cy + half_s),
        list(cx - half_s, cy - half_s)
      ))
    )
  )
})

out <- file.path(WEBMAP_DIR, "suitability.js")
writeLines(paste0("var SUIT = ", toJSON(
  list(type = "FeatureCollection", features = suit_features),
  auto_unbox = TRUE), ";"), out)
cat(sprintf("  suitability.js  %d cells  %.1f MB\n", nrow(suit_df), file.size(out)/1e6))


# =============================================================================
# 2. DETECTIONS
# =============================================================================
cat("2/4  Detections…\n")

# Year stored as a property so the browser can filter by slider position
det_list <- list(
  type     = "FeatureCollection",
  features = lapply(seq_len(nrow(eab)), function(i)
    list(type       = "Feature",
         properties = list(year = eab$year[i]),
         geometry   = list(type        = "Point",
                           coordinates = list(eab$lon[i], eab$lat[i]))))
)

out <- file.path(WEBMAP_DIR, "detections.js")
writeLines(paste0("var DETECTIONS = ", toJSON(det_list, auto_unbox = TRUE), ";"), out)
cat(sprintf("  detections.js  %d points  %.1f MB\n", nrow(eab), file.size(out)/1e6))


# =============================================================================
# 3. PER-YEAR SPREAD
# =============================================================================
cat("3/4  Spread layers…\n")

# Same year filter as the GIF — skip sparse early observed years
all_export_years <- as.integer(names(annual_invaded))
det_cumul <- sapply(all_export_years, function(yr)
  nrow(eab |> filter(year <= min(yr, OBS_END))))
export_years <- all_export_years[all_export_years >= SIM_START | det_cumul >= 10]

half <- cell_res_deg / 2

for (yr in export_years) {
  inv_vec <- annual_invaded[[as.character(yr)]]
  idx     <- which(inv_vec)
  if (length(idx) == 0) next

  xy <- xyFromCell(suit_sim, idx)

  feat_list <- list(
    type     = "FeatureCollection",
    features = lapply(seq_along(idx), function(i) {
      cx <- xy[i, 1]; cy <- xy[i, 2]
      list(type       = "Feature",
           properties = list(year = yr),
           geometry   = list(
             type        = "Polygon",
             coordinates = list(list(
               list(cx-half, cy-half), list(cx+half, cy-half),
               list(cx+half, cy+half), list(cx-half, cy+half),
               list(cx-half, cy-half)
             ))
           ))
    })
  )

  out <- file.path(WEBMAP_DIR, sprintf("spread_%d.js", yr))
  writeLines(paste0("var SPREAD_", yr, " = ",
                    toJSON(feat_list, auto_unbox = TRUE), ";"), out)
  cat(sprintf("  spread_%d.js  %d cells  %.0f KB\n",
              yr, length(idx), file.size(out)/1024))
}


# =============================================================================
# 3b. TEMPORAL ISOLINES
# =============================================================================
# Invasion front boundaries for 4 key years — dissolved, boundary-extracted,
# and Chaikin-smoothed to remove jagged grid edges.
cat("3b   Isolines…\n")

library(sf)

ISOLINE_YEARS <- c(2025, 2030, 2035, 2040)
half_i        <- cell_res_deg / 2

isoline_features <- lapply(ISOLINE_YEARS, function(yr) {
  inv_vec <- annual_invaded[[as.character(yr)]]
  if (is.null(inv_vec)) return(NULL)
  idx <- which(inv_vec)
  if (length(idx) == 0) return(NULL)

  xy <- xyFromCell(suit_sim, idx)

  polys <- lapply(seq_along(idx), function(i) {
    cx <- xy[i, 1]; cy <- xy[i, 2]
    st_polygon(list(matrix(c(
      cx-half_i, cy-half_i, cx+half_i, cy-half_i,
      cx+half_i, cy+half_i, cx-half_i, cy+half_i,
      cx-half_i, cy-half_i
    ), ncol = 2, byrow = TRUE)))
  })

  dissolved <- st_union(st_sfc(polys, crs = 4326))
  boundary  <- st_boundary(dissolved)
  boundary  <- smooth(boundary, method = "chaikin", refinements = 5)

  st_sf(year = yr, observed = yr <= OBS_END, geometry = boundary)
})

isolines_sf <- do.call(rbind, Filter(Negate(is.null), isoline_features))

out <- file.path(WEBMAP_DIR, "isolines.js")
writeLines(paste0("var ISOLINES = ", geojsonsf::sf_geojson(isolines_sf), ";"), out)
cat(sprintf("  isolines.js  %d boundaries  %.1f MB\n",
            nrow(isolines_sf), file.size(out)/1e6))


# =============================================================================
# 4. MANIFEST
# =============================================================================
# Tells the webmap which years have data, and where the observed/projected
# boundary sits — avoids hardcoding these values in index.html.

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
# 5. COPY index.html — inject spread <script> tags
# =============================================================================
# The <!-- SPREAD_SCRIPTS --> placeholder is replaced with a <script src> tag
# for every exported year so the browser can load each year's data on demand.

html_src <- file.path(SCRIPT_DIR, "index.html")
html_dst <- file.path(SCRIPT_DIR, "webmap", "index.html")

if (file.exists(html_src)) {
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
cat("Done. Serve the map with:\n")
cat("  cd", file.path(SCRIPT_DIR, "webmap"), "\n")
cat("  python3 -m http.server 8000\n")
cat("------------------------------------------------------------\n")
