# 02_flightr_prepare_tags.R
# Convert GeoLight twilight events into TAGS format for FLightR.
#
# Inputs:
#  - .lig file
#  - *_twilights.csv from script 01
# Output:
#  - *_TAGS_twilights.csv

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GeoLight)
  library(FLightR)
})

threshold <- cfg$geolight$threshold %||% 32
inter_dir <- cfg$paths$intermediate_dir

for (lig in cfg$paths$lig_files) {
  tag_id <- tools::file_path_sans_ext(basename(lig))
  message("=== TAGS conversion: ", tag_id, " ===")

  d_lux <- read_lig_basic(lig, tz = cfg$flightr$tz %||% "GMT")
  GL_light <- subset(d_lux, select = c("Date", "Light"))

  twl_path <- file.path(inter_dir, paste0(tag_id, "_twilights.csv"))
  if (!file.exists(twl_path)) stop("Missing twilights CSV (run script 01): ", twl_path)
  twl <- read.csv(twl_path, stringsAsFactors = FALSE)

  # Coerce to expected types
  twl[,1] <- as.POSIXct(twl[,1], tz = cfg$flightr$tz %||% "GMT")
  twl[,2] <- as.POSIXct(twl[,2], tz = cfg$flightr$tz %||% "GMT")
  twl[,3] <- as.numeric(twl[,3])
  twl_mat <- as.matrix(twl)

  TAGS_twilights <- GeoLight2TAGS(GL_light, twl_mat, threshold = threshold)
  TAGS_twilights <- na.omit(TAGS_twilights)

  out_csv <- file.path(inter_dir, paste0(tag_id, "_TAGS_twilights.csv"))
  write.csv(TAGS_twilights, out_csv, quote = FALSE, row.names = FALSE)
  message("Wrote: ", out_csv)
}
