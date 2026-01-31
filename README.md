# Purple Martin geolocator migration pipeline (GeoLight + FLightR)

A **GitHub-ready, data-free** workflow to reconstruct Purple Martin migration routes and stopovers from **light-level geolocator**
`.lig` files.

This repo implements two complementary analyses:

1) **GeoLight**: twilight detection → loess outlier filtering → coarse positions → residency inference (`changeLight`)  
2) **FLightR**: TAGS conversion → calibration → spatial grid → particle filter → stopover summary

> **No data are included**. Point the pipeline to your local `.lig` files via `config/config.yml` (git-ignored).

---

## Repository layout

- `scripts/01_geolight_twilight_locations.R`  
  Detect twilights, filter outliers, compute coarse coordinates, run `changeLight`, and export intermediate CSVs.

- `scripts/02_flightr_prepare_tags.R`  
  Convert `.lig` + GeoLight twilights to **TAGS** format for FLightR.

- `scripts/03_flightr_particle_filter.R`  
  Run the full FLightR workflow (Proc.data → calibration → grid → particle filter) and save results.

- `scripts/04_stopover_summary.R`  
  Extract stopover table from FLightR results and export tidy CSVs (duration, arrival, departure, distance).

- `scripts/05_mapping_leaflet_optional.R`  
  Optional mapping template (leaflet). API keys stay in config and out of git.

---

## Install packages

```r
install.packages(c("maps","GeoLight","raster","RCurl","yaml","dplyr","stringr","lubridate"))
install.packages("FLightR")
```

Optional:
```r
install.packages(c("leaflet","ggplot2"))
```

---

## Configure

```bash
cp config/config_example.yml config/config.yml
```

Edit `config/config.yml`:
- `paths.lig_files`: your local `.lig` files
- `geolight.threshold`, `geolight.offset`, `geolight.coord_tol`
- calibration coordinates + dates
- grid bounding box
- `flightr.n_particles` (trial 1e4; final 1e6)

---

## Run

```r
source("scripts/01_geolight_twilight_locations.R")
source("scripts/02_flightr_prepare_tags.R")
source("scripts/03_flightr_particle_filter.R")
source("scripts/04_stopover_summary.R")
# optional:
source("scripts/05_mapping_leaflet_optional.R")
```

---

## Outputs

All outputs go under `outputs/`:

- `outputs/intermediate/`  
  `*_twilights.csv`, `*_twilights_loess.csv`, `*_locations.csv`, `*_TAGS_twilights.csv`

- `outputs/results/`  
  `*_geolight_objects.rds`, `*_FLightR_result.rds`, `*_stationary_summary.rds`, `*_stopovers.csv`

- `outputs/figures/`  
  diagnostic plots (as you add them)

---

## Suggested repo name

**purple-martin-geolocator-flightr-pipeline**
