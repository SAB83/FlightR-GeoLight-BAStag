suppressPackageStartupMessages({
  library(yaml)
  source("R/helpers.R")
})

cfg <- yaml::read_yaml("config/config.yml")

dir_create(
  cfg$outputs$dir,
  cfg$outputs$figures_dir,
  cfg$outputs$results_dir,
  cfg$paths$intermediate_dir
)

message("Setup ok.")
