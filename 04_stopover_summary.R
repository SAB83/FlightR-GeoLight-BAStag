# 04_stopover_summary.R
# Extract tidy stopover tables from FLightR results.
#
# Inputs: outputs/results/*_FLightR_result.rds
# Outputs: outputs/results/*_stopovers.csv + *_stationary_summary.rds

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(FLightR)
  library(dplyr)
})

res_dir <- cfg$outputs$results_dir
rds_files <- list.files(res_dir, pattern = "_FLightR_result\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No *_FLightR_result.rds found. Run script 03 first.")

for (rf in rds_files) {
  tag_id <- sub("_FLightR_result\.rds$", "", basename(rf))
  message("=== Stopovers: ", tag_id, " ===")

  Result <- readRDS(rf)

  Summary <- stationary.migration.summary(Result, prob.cutoff = 0.1, min.stay = 1)
  saveRDS(Summary, file.path(res_dir, paste0(tag_id, "_stationary_summary.rds")))

  st <- Summary$Stationary.periods
  if (is.null(st) || nrow(st) == 0) {
    message("No stationary periods for ", tag_id)
    next
  }

  st <- st %>%
    mutate(stopover_duration_days = as.numeric(difftime(Departure.Q.50, Arrival.Q.50, units = "days"))) %>%
    filter(!is.na(stopover_duration_days)) %>%
    select(Meanlat, Meanlon, Distance2, Distance2cumulative, Arrival.Q.50, Departure.Q.50, stopover_duration_days)

  out_csv <- file.path(res_dir, paste0(tag_id, "_stopovers.csv"))
  write.csv(st, out_csv, quote = FALSE, row.names = FALSE)
  message("Wrote: ", out_csv)
}
