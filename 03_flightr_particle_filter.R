# 03_flightr_particle_filter.R
# Runs the full FLightR reconstruction and saves Result objects.
#
# Inputs:
#  - *_TAGS_twilights.csv
# Outputs:
#  - *_FLightR_result.rds

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(FLightR) })

inter_dir <- cfg$paths$intermediate_dir
res_dir <- cfg$outputs$results_dir

start_date <- as.POSIXct(cfg$flightr$start_date, tz = cfg$flightr$tz %||% "GMT")
end_date   <- as.POSIXct(cfg$flightr$end_date,   tz = cfg$flightr$tz %||% "GMT")

calib_start <- as.POSIXct(cfg$flightr$calibration_start, tz = cfg$flightr$tz %||% "GMT")
calib_stop  <- as.POSIXct(cfg$flightr$calibration_stop,  tz = cfg$flightr$tz %||% "GMT")

known_lon <- cfg$geolight$known_coord$lon
known_lat <- cfg$geolight$known_coord$lat

grid_cfg <- cfg$flightr$grid

for (lig in cfg$paths$lig_files) {
  tag_id <- tools::file_path_sans_ext(basename(lig))
  message("=== FLightR: ", tag_id, " ===")

  tags_csv <- file.path(inter_dir, paste0(tag_id, "_TAGS_twilights.csv"))
  if (!file.exists(tags_csv)) stop("Missing TAGS CSV (run script 02): ", tags_csv)

  Proc.data <- get.tags.data(
    tags_csv,
    start_date,
    end_date,
    log.light.borders = "auto",
    log.irrad.borders = "auto",
    saves = c("auto","mean","max"),
    measurement.period = 120,
    impute.on.boundaries = FALSE
  )

  Calibration.periods <- data.frame(
    calibration.start = calib_start,
    calibration.stop  = calib_stop,
    lon = known_lon,
    lat = known_lat
  )

  Calib <- make.calibration(
    Proc.data,
    Calibration.periods,
    model.ageing = FALSE,
    plot.each = FALSE,
    plot.final = FALSE,
    likelihood.correction = "auto",
    fixed.logSlope = c(NA, NA),
    suggest.irrad.borders = FALSE,
    return.slopes = FALSE
  )

  Grid <- make.grid(
    left   = grid_cfg$left,
    bottom = grid_cfg$bottom,
    right  = grid_cfg$right,
    top    = grid_cfg$top,
    distance.from.land.allowed.to.use  = grid_cfg$distance_from_land_allowed_to_use,
    distance.from.land.allowed.to.stay = grid_cfg$distance_from_land_allowed_to_stay
  )

  prerun <- make.prerun.object(
    Proc.data, Grid,
    start = c(known_lon, known_lat),
    end = NA,
    Calibration = Calib
  )

  Result <- run.particle.filter(
    prerun,
    threads = cfg$flightr$threads %||% -1,
    nParticles = cfg$flightr$n_particles %||% 10000,
    known.last = cfg$flightr$known_last %||% FALSE,
    precision.sd = cfg$flightr$precision_sd %||% 25,
    check.outliers = cfg$flightr$check_outliers %||% TRUE,
    b = cfg$flightr$b %||% 1500
  )

  out_rds <- file.path(res_dir, paste0(tag_id, "_FLightR_result.rds"))
  saveRDS(Result, out_rds)
  message("Saved: ", out_rds)
}
