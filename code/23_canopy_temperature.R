#!/usr/bin/env Rscript

# Canopy temperature: how far it departs from air temperature, and an
# independent check on the sensor.
#
# Why this exists. The manuscript stated that canopy-top leaf temperatures at
# this site "typically run several degrees C below air temperature". That was
# uncited, and it is wrong in sign: the canopy runs WARMER than air throughout
# daylight. This script produces the numbers that replace it, so every value
# quoted in the text traces to a script and a machine-readable output.
#
# Two independent estimates are compared:
#   T_canopy  measured, the mean of two infrared radiometers (SI-111 and a
#             second IR sensor), 27 Oct 2022 - 18 Oct 2024, exported by
#             'from Matlab/Export_data.m'
#   T_lw      derived from outgoing longwave via Stefan-Boltzmann, using the
#             NR01 four-component radiometer. Available for the whole record,
#             read straight from the biomet database because LW_IN/LW_OUT are
#             deliberately not exported.
#
# The derived series is used here only to validate the sensor. It is not a
# reported product.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

# ─────────────────────────────────────────────────────────────────────────────
# config

tz_local <- "America/Lima"

# broadband emissivity of a closed tropical canopy. the sign of the canopy-air
# offset is insensitive to this over any defensible range; the sensitivity is
# reported in the output so the assumption is visible.
emissivity      <- 0.98
emissivity_test <- c(0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
stefan          <- 5.670374419e-8

# day/night follow the project convention: night = Rg < 10 W m-2,
# day = Rg > 44 W m-2 (equivalent to PAR > 100). high light is a stricter cut
# used to characterise the offset when the radiative load is largest.
night_sw_max <- 10
day_sw_min   <- 44
high_sw_min  <- 600

# Root of the raw Biomet.net database, used ONLY for the independent longwave
# validation below. It is not part of this repository. Set PETNR_DB_ROOT to
# point at your own copy; when it is unset or absent the script still writes
# table_canopy_temperature.csv from data/dataset_from_matlab.csv, but WITHOUT
# the Stefan-Boltzmann validation rows -- so the table will be shorter than the
# published one. See README ("Reproducing the tables").
db_root  <- Sys.getenv("PETNR_DB_ROOT", unset = "")
if (!nzchar(db_root)) {
  message("PETNR_DB_ROOT is not set; skipping the independent longwave validation.")
}
site     <- "PETNR"
lw_years <- 2022:2024          # only the sensor overlap is needed for validation

flux_csv <- file.path(output_path, "dataset_from_matlab.csv")
out_csv  <- file.path(out_base, "table_canopy_temperature.csv")

# ─────────────────────────────────────────────────────────────────────────────
# load the exported record

if (!file.exists(flux_csv)) {
  stop("exported flux table not found: ", flux_csv,
       "\nRun 'from Matlab/Export_data.m' and copy the result into data/.")
}

required_cols <- c("tv_dt", "T_canopy", "TA_1_1_1", "SW_IN_1_1_1")
df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nT_canopy is written by 'from Matlab/Export_data.m'. If it is absent, the",
       " sensor CSV has not been ingested by PETNR_convert_to_database.m,",
       " or the cleaning stages have not been re-run since.")
}

d <- df_raw %>%
  transmute(
    ts   = dmy_hms(tv_dt, tz = tz_local, quiet = TRUE),
    Tcan = as.numeric(T_canopy),
    Ta   = as.numeric(TA_1_1_1),
    SW   = as.numeric(SW_IN_1_1_1)
  )

if (!any(is.finite(d$Tcan))) {
  stop("T_canopy is present but entirely NaN. The trace exists in the export but",
       " was never populated upstream — check the ingest block in",
       " PETNR_convert_to_database.m and that cleaning stages 1, 2 and 7 have run.")
}

cov <- d %>% filter(is.finite(Tcan))
message(sprintf("T_canopy: %d finite half-hours, %s to %s, range %.1f to %.1f C",
                nrow(cov), format(min(cov$ts), "%d %b %Y"),
                format(max(cov$ts), "%d %b %Y"), min(cov$Tcan), max(cov$Tcan)))

# ─────────────────────────────────────────────────────────────────────────────
# canopy - air offset from the measured sensor

k   <- is.finite(d$Tcan) & is.finite(d$Ta)
off <- d$Tcan - d$Ta

subset_mask <- list(
  all        = k,
  daytime    = k & is.finite(d$SW) & d$SW > day_sw_min,
  high_light = k & is.finite(d$SW) & d$SW > high_sw_min,
  nighttime  = k & is.finite(d$SW) & d$SW < night_sw_max
)

offset_rows <- bind_rows(lapply(names(subset_mask), function(nm) {
  m <- subset_mask[[nm]]
  tibble(
    quantity        = "canopy minus air temperature",
    source          = "measured sensor",
    subset          = nm,
    n               = sum(m),
    mean_C          = mean(off[m]),
    median_C        = median(off[m]),
    sd_C            = sd(off[m]),
    p05_C           = unname(quantile(off[m], 0.05)),
    p95_C           = unname(quantile(off[m], 0.95)),
    frac_canopy_warmer = mean(off[m] > 0)
  )
}))

# ─────────────────────────────────────────────────────────────────────────────
# independent estimate from outgoing longwave

readtrace <- function(p, size = 4) {
  if (!file.exists(p)) return(NULL)
  readBin(p, "numeric", n = file.info(p)$size / size, size = size, endian = "little")
}
readtv <- function(p) readBin(p, "numeric", n = file.info(p)$size / 8, size = 8, endian = "little")

lw <- bind_rows(lapply(lw_years, function(y) {
  cn <- file.path(db_root, y, site, "Clean", "ThirdStage")
  tvf <- file.path(cn, "clean_tv")
  if (!file.exists(tvf)) return(NULL)
  tv <- readtv(tvf)
  g <- function(v) { x <- readtrace(file.path(cn, v)); if (is.null(x)) rep(NA_real_, length(tv)) else x }
  # matlab datenum is a float, so snap to the exact half-hour grid before joining
  tibble(
    ts    = as.POSIXct(round(((tv - 719529) * 86400) / 1800) * 1800,
                       origin = "1970-01-01", tz = "UTC"),
    LWout = g("LW_OUT_1_1_1"),
    LWin  = g("LW_IN_1_1_1")
  )
}))

if (!nrow(lw) || !any(is.finite(lw$LWout))) {
  message("longwave traces unavailable; skipping the independent validation")
  validation_rows <- tibble()
} else {
  # outgoing longwave contains a reflected component, so LW_IN is subtracted out
  # before inverting Stefan-Boltzmann
  canopy_from_lw <- function(LWout, LWin, eps) {
    ((LWout - (1 - eps) * LWin) / (eps * stefan))^0.25 - 273.15
  }
  lw$Tlw <- canopy_from_lw(lw$LWout, lw$LWin, emissivity)

  key <- function(x) format(x, "%Y-%m-%d %H:%M")
  j <- inner_join(d %>% mutate(k = key(ts)) %>% select(k, Tcan, Ta, SW),
                  lw %>% mutate(k = key(ts)) %>% select(k, Tlw, LWout, LWin),
                  by = "k") %>%
       filter(is.finite(Tcan), is.finite(Tlw))

  if (nrow(j) < 100) {
    message("too few paired half-hours for validation; skipping")
    validation_rows <- tibble()
  } else {
    fit <- lm(Tcan ~ Tlw, data = j)
    validation_rows <- tibble(
      quantity = "measured vs longwave-derived canopy temperature",
      source   = sprintf("Stefan-Boltzmann, emissivity = %.2f", emissivity),
      subset   = "paired half-hours",
      n        = nrow(j),
      slope     = unname(coef(fit)[2]),
      intercept = unname(coef(fit)[1]),
      r         = cor(j$Tcan, j$Tlw),
      R2        = summary(fit)$r.squared,
      RMSE_C    = sqrt(mean((j$Tcan - j$Tlw)^2)),
      bias_C    = mean(j$Tcan - j$Tlw)
    )
    message(sprintf("validation: n = %d, slope %.3f, r %.3f, RMSE %.2f C, bias %+.2f C",
                    nrow(j), coef(fit)[2], cor(j$Tcan, j$Tlw),
                    sqrt(mean((j$Tcan - j$Tlw)^2)), mean(j$Tcan - j$Tlw)))

    # how much of the agreement depends on the emissivity choice
    sens <- bind_rows(lapply(emissivity_test, function(e) {
      tlw <- canopy_from_lw(j$LWout, j$LWin, e)
      tibble(quantity = "emissivity sensitivity",
             source   = sprintf("emissivity = %.2f", e),
             subset   = "paired half-hours",
             n        = nrow(j),
             bias_C   = mean(j$Tcan - tlw),
             RMSE_C   = sqrt(mean((j$Tcan - tlw)^2)))
    }))
    validation_rows <- bind_rows(validation_rows, sens)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# save

res <- bind_rows(offset_rows, validation_rows) %>%
  mutate(coverage_start = format(min(cov$ts), "%Y-%m-%d"),
         coverage_end   = format(max(cov$ts), "%Y-%m-%d"),
         n_finite_T_canopy = nrow(cov))

readr::write_csv(res, out_csv)
cat("saved:", out_csv, "\n\n")

print(as.data.frame(
  offset_rows %>%
    transmute(subset, n,
              mean_C = round(mean_C, 2), median_C = round(median_C, 2),
              p05_C = round(p05_C, 2), p95_C = round(p95_C, 2),
              pct_canopy_warmer = round(100 * frac_canopy_warmer))
), row.names = FALSE)

dm <- offset_rows$mean_C[offset_rows$subset == "daytime"]
hm <- offset_rows$mean_C[offset_rows$subset == "high_light"]
nm <- offset_rows$mean_C[offset_rows$subset == "nighttime"]
cat(sprintf(
  "\ncanopy exceeded air by %.1f C during daylight and %.1f C under high light,\nand fell %.1f C below air at night.\n",
  dm, hm, abs(nm)))
