suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(stringr)
  library(lubridate)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

dir_create <- function(...) {
  for (p in c(...)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# Read a .lig file with BAStag-style structure (Valid, Date, Julian, Light)
read_lig_basic <- function(file, skip = 0, tz = "GMT") {
  d <- read.csv(
    file, header = FALSE, skip = skip,
    col.names = c("Valid", "Date", "Julian", "Light"),
    colClasses = c("character", "character", "numeric", "integer"),
    stringsAsFactors = FALSE
  )
  d$Date <- as.POSIXct(strptime(d$Date, "%d/%m/%y %H:%M:%S", tz = tz))
  d
}

write_csv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
}
