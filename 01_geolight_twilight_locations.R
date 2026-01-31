# 01_geolight_twilight_locations.R
# 1) Read .lig (light data)
# 2) Detect twilights using twilightCalc (threshold)
# 3) Loess-filter twilight outliers (loessFilter)
# 4) Estimate sun elevation from a known stationary window (getElevation)
# 5) Convert twilights to coarse locations (coord), with equinox tolerance (tol)
# 6) Infer residency/stopovers (changeLight)
# 7) Export intermediate files (CSV + RDS)

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GeoLight)
  library(maps)
})

threshold <- cfg$geolight$threshold %||% 32
k_loess   <- cfg$geolight$loess_k %||% 2
tol_coord <- cfg$geolight$coord_tol %||% 0.13
known_lon <- cfg$geolight$known_coord$lon
known_lat <- cfg$geolight$known_coord$lat
inter_dir <- cfg$paths$intermediate_dir

for (lig in cfg$paths$lig_files) {
  tag_id <- tools::file_path_sans_ext(basename(lig))
  message("=== GeoLight: ", tag_id, " ===")

  d_lux <- read_lig_basic(lig, tz = "GMT")
  if (nrow(d_lux) < 10) stop("Too few rows in .lig file: ", lig)

  # --- Twilight detection ---
  twl <- twilightCalc(datetime = d_lux$Date,
                      light    = d_lux$Light,
                      LightThreshold = threshold,
                      ask = TRUE)

  # --- Loess filter (automatic outlier removal) ---
  twl_loess <- loessFilter(tFirst = twl[,1],
                           tSecond = twl[,2],
                           type = 2,
                           twl = twl,
                           k = k_loess,
                           plot = TRUE)

  # --- Sun elevation calibration on known stationary dates (example indices) ---
  si <- cfg$geolight$sun_elev_calib_window$start_index %||% 33
  ei <- cfg$geolight$sun_elev_calib_window$end_index %||% 63
  idx <- si:ei
  idx <- idx[idx >= 1 & idx <= nrow(twl)]

  sun_elev <- getElevation(tFirst = twl[idx,1],
                           tSecond = twl[idx,2],
                           type = twl[idx,3],
                           known.coord = c(known_lon, known_lat),
                           plot = TRUE)

  # --- Coarse locations (will be NA near equinox per tol) ---
  locs <- coord(tFirst = twl[,1],
                tSecond = twl[,2],
                type = twl[,3],
                degElevation = sun_elev,
                tol = tol_coord)

  # --- Quick map ---
  plot(locs, pch="*", col="red", xlab="Longitude", ylab="Latitude")
  map("world", add = TRUE)

  # --- Residency inference ---
  stop_tbl <- changeLight(tFirst = twl[,1],
                          tSecond = twl[,2],
                          type = 2,
                          twl = twl,
                          quantile = cfg$changelight$quantile %||% 0.9,
                          rise.prob = NA,
                          set.prob = NA,
                          days = cfg$changelight$days %||% 2,
                          plot = TRUE,
                          summary = TRUE)

  # --- Export ---
  write_csv_safe(as.data.frame(twl), file.path(inter_dir, paste0(tag_id, "_twilights.csv")))
  write_csv_safe(as.data.frame(twl_loess), file.path(inter_dir, paste0(tag_id, "_twilights_loess.csv")))
  write_csv_safe(as.data.frame(locs), file.path(inter_dir, paste0(tag_id, "_locations.csv")))
  if (!is.null(stop_tbl$schedule)) {
    write_csv_safe(as.data.frame(stop_tbl$schedule), file.path(inter_dir, paste0(tag_id, "_changelight_schedule.csv")))
  }

  saveRDS(list(lig = lig, d_lux = d_lux, twl = twl, twl_loess = twl_loess, locs = locs, stop = stop_tbl),
          file.path(cfg$outputs$results_dir, paste0(tag_id, "_geolight_objects.rds")))

  message("Saved outputs for ", tag_id)
}
