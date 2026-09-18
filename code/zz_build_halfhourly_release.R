#!/usr/bin/env Rscript
# =============================================================================
# zz_build_halfhourly_release.R
# -----------------------------------------------------------------------------
# Assemble the half-hourly Zenodo DATASET record: the CO2, energy and
# meteorological record actually analysed by the paper, at native 30-min
# resolution.
#
# Scope decision: this is a companion to the aggregated deposit built by
# zz_build_zenodo_deposit.R. The aggregates back the figures and tables; this
# file is the series they were computed from, so that a third party can redo the
# processing rather than only re-plot the summaries.
#
# Selection is an ALLOWLIST, not an exclusion list. Only the 24 columns named in
# `keep_cols` are ever written, and the script hard-stops if the output contains
# anything else. Non-CO2 trace columns present in the source export are
# therefore never carried into the release, and cannot be reintroduced by an
# upstream change to the export.
#
# The record is truncated at 2024-10-15 23:30 local, the end of the processed
# record reported in Table S2. The source export runs a container past that date
# and carries provisional values that are not part of this study.
#
# Run:     Rscript code/zz_build_halfhourly_release.R
# Output:  <repo>/halfhourly_release/   (gitignored; upload to Zenodo)
# =============================================================================

## ---- locate repo root (standalone; does not source paths.R) ----------------
get_script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(normalizePath(f))
  normalizePath(file.path(getwd(), "code", "zz_build_halfhourly_release.R"))
}
script_path <- tryCatch(get_script_path(), error = function(e) normalizePath("."))
repo_root <- dirname(dirname(script_path))
data_dir  <- file.path(repo_root, "data")
out_dir   <- file.path(repo_root, "halfhourly_release")

cat("repo_root :", repo_root, "\n")
cat("out_dir   :", out_dir, "\n\n")

## ---- the allowlist ----------------------------------------------------------
# tv_dt is consumed and replaced by DateTime; the other 23 pass through.
time_col  <- "tv_dt"
keep_cols <- c(
  "NEE", "NEE_f", "FC",
  "H", "LE",
  "TA_1_1_1", "RH_1_1_1", "VPD_1_1_1",
  "WS_1_1_1", "WD_1_1_1", "USTAR",
  "TS_3", "SWC_1_1_1",
  "PPFD_IN_1_1_1", "SW_IN_1_1_1", "NETRAD_1_1_1",
  "G_1_1_1", "G_1_1_2", "G",
  "SC_profiler", "SC_model", "SC",
  "T_canopy"
)
out_cols <- c("DateTime", keep_cols)   # 24 columns

record_end <- "2024-10-15 23:30:00"
tz_local   <- "America/Lima"

## ---- units and descriptions -------------------------------------------------
dict <- read.csv(text = "column,unit,description
DateTime,,Half-hourly timestamp in local time (America/Lima); labels the END of the averaging interval
NEE,umol CO2 m-2 s-1,Net ecosystem exchange = FC + SC; negative = uptake. Not u*-filtered here
NEE_f,umol CO2 m-2 s-1,Gap-filled NEE from REddyProc; empty for years where REddyProc did not converge
FC,umol CO2 m-2 s-1,Turbulent vertical CO2 flux above the canopy
H,W m-2,Sensible heat flux
LE,W m-2,Latent heat flux
TA_1_1_1,degC,Air temperature at 47 m (HMP155)
RH_1_1_1,%,Relative humidity at 47 m (HMP155)
VPD_1_1_1,hPa,Vapour pressure deficit as exported; divide by 10 for kPa
WS_1_1_1,m s-1,Horizontal wind speed
WD_1_1_1,degrees,Wind direction (0 = north, clockwise)
USTAR,m s-1,Friction velocity
TS_3,degC,Soil temperature, depth-weighted mean over 0-100 cm (SoilVue; 2022 onward)
SWC_1_1_1,cm,Soil water content as integrated 0-100 cm storage in centimetres, NOT volumetric (2022 onward)
PPFD_IN_1_1_1,umol photons m-2 s-1,Incoming PAR. DERIVED as SW_IN_1_1_1*2.3-1.22, not measured; not truncated at zero
SW_IN_1_1_1,W m-2,Incoming shortwave radiation (NR01)
NETRAD_1_1_1,W m-2,Net radiation (NR01)
G_1_1_1,W m-2,Soil heat flux, plate 1 at 10 cm depth (from 2020)
G_1_1_2,W m-2,Soil heat flux, plate 2 at 10 cm depth (from 2023)
G,W m-2,Soil heat flux, mean of the available plates
SC_profiler,umol CO2 m-2 s-1,Storage flux from the 8-level CO2 profile (profiler years only)
SC_model,umol CO2 m-2 s-1,Storage flux reconstructed by the DNN (v31). The FirstStage [-100 100] bound is NOT applied here
SC,umol CO2 m-2 s-1,Storage flux actually used for NEE: SC_model before 2024, SC_profiler from 2024 onward
T_canopy,degC,Canopy temperature, mean of the SI-111 and IR radiometer sensors (2022-10-27 onward)
", stringsAsFactors = FALSE, colClasses = "character")

## ---- read, subset, truncate --------------------------------------------------
src <- file.path(data_dir, "dataset_from_matlab.csv")
if (!file.exists(src)) {
  stop("source export not found: ", src,
       "\nRun the upstream MATLAB export and copy the result into data/.")
}

raw <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
cat("source:", nrow(raw), "rows x", ncol(raw), "columns\n")

missing_cols <- setdiff(c(time_col, keep_cols), names(raw))
if (length(missing_cols)) {
  stop("the source export is missing expected columns: ",
       paste(missing_cols, collapse = ", "))
}

# parse "01-Jan-2017 00:30:00" locale-independently
old_lc <- Sys.getlocale("LC_TIME"); invisible(Sys.setlocale("LC_TIME", "C"))
dt <- as.POSIXct(raw[[time_col]], format = "%d-%b-%Y %H:%M:%S", tz = tz_local)
invisible(Sys.setlocale("LC_TIME", old_lc))
if (all(is.na(dt))) stop("could not parse ", time_col, " as '%d-%b-%Y %H:%M:%S'")

tmax <- as.POSIXct(record_end, tz = tz_local)
sel  <- !is.na(dt) & dt <= tmax

out <- raw[sel, keep_cols, drop = FALSE]
out <- cbind(DateTime = format(dt[sel], "%Y-%m-%d %H:%M:%S"), out,
             stringsAsFactors = FALSE)

for (cc in keep_cols) {
  out[[cc]] <- suppressWarnings(as.numeric(out[[cc]]))
  out[[cc]][is.nan(out[[cc]])] <- NA
}

cat("released:", nrow(out), "rows x", ncol(out), "columns",
    "( truncated", nrow(raw), "->", nrow(out), "rows at", record_end, ")\n")

## ---- hard stop: the output must be exactly the allowlist ---------------------
if (!identical(names(out), out_cols)) {
  stop("output columns do not match the allowlist.\n  got     : ",
       paste(names(out), collapse = ", "), "\n  expected: ",
       paste(out_cols, collapse = ", "))
}
if (ncol(out) != 24L) stop("expected 24 columns, got ", ncol(out))

## ---- write -------------------------------------------------------------------
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, showWarnings = FALSE)

csv_name <- "tambopata_PE-TNR_halfhourly_2017_2024.csv"
csv_path <- file.path(out_dir, csv_name)
write.csv(out, csv_path, row.names = FALSE, na = "NA")

write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(lines, con)
}

# data dictionary, restricted to the columns actually written
dd <- dict[match(names(out), dict$column), ]
write.csv(dd, file.path(out_dir, "data_dictionary.csv"), row.names = FALSE, na = "")

write_utf8(c(
  "Tambopata (PE-TNR) half-hourly record",
  "Copyright (c) 2026 Rafael Stern and the Tambopata flux-tower authors.",
  "",
  "This dataset is licensed under the Creative Commons Attribution 4.0",
  "International License (CC BY 4.0).",
  "You are free to share and adapt the material for any purpose, provided you",
  "give appropriate credit (cite the associated paper and this deposit).",
  "",
  "Full license text: https://creativecommons.org/licenses/by/4.0/legalcode"
), file.path(out_dir, "LICENSE.txt"))

n_meas <- sum(is.finite(out$NEE))
write_utf8(c(
  "# Tambopata (PE-TNR) half-hourly eddy-covariance record, 2017-2024",
  "",
  "Half-hourly CO2, energy and meteorological measurements underlying:",
  "",
  "> Stern, R., Cruz, R., Negron-Juarez, R., Ramos, E., Knox, S. H., Hoyt, A. M.,",
  "> Jackson, R. B., Li, F., Dayalu, A., Pernak, R., Grace, J., Gloor, E.,",
  "> Salinas, N., & Cosio, E. G. *Long-Term Flux Measurements in Southwestern",
  "> Amazonia Reveal a Near-Neutral CO2 Balance and Unique Effects of El",
  "> Nino-Southern Oscillation.* Global Change Biology (submitted).",
  "",
  "**Site:** Tambopata National Reserve, Madre de Dios, Peru.",
  "AmeriFlux **PE-TNR** (AndesFlux/PUCP), 12.8314 S, 69.2836 W, 214 m a.s.l.",
  "",
  sprintf("**Coverage:** %s to %s local time (America/Lima), %s rows, %d columns.",
          format(min(dt[sel]), "%Y-%m-%d"), format(max(dt[sel]), "%Y-%m-%d"),
          format(nrow(out), big.mark = ","), ncol(out)),
  sprintf("%s half-hours carry a finite NEE.", format(n_meas, big.mark = ",")),
  "",
  "## Files",
  "",
  sprintf("- `%s` - the record", csv_name),
  "- `data_dictionary.csv` - column, unit, description",
  "- `LICENSE.txt` - CC BY 4.0",
  "",
  "## Conventions",
  "",
  "- Timestamps label the **end** of each 30-minute averaging interval.",
  "- NEE and FC are in umol CO2 m-2 s-1; **negative = uptake**, positive = release.",
  "- `NEE = FC + SC`, summed with omitnan, so a half-hour with valid FC and",
  "  missing SC yields NEE = FC.",
  "- `SC` is the storage term actually used: the DNN reconstruction (`SC_model`)",
  "  before 2024, the profile measurement (`SC_profiler`) from 2024 onward.",
  "  Roughly 77% of the record uses the modelled term.",
  "- `PPFD_IN_1_1_1` is **derived**, not measured: PAR = SW_IN * 2.3 - 1.22,",
  "  applied to the whole record because the original PAR sensor had a",
  "  calibration fault. Values are not truncated at zero.",
  "- `SWC_1_1_1` is integrated 0-100 cm soil-water **storage in centimetres**,",
  "  not volumetric water content. Soil sensors were installed in September 2022.",
  "- `VPD_1_1_1` is in hPa as exported; divide by 10 for kPa.",
  "",
  "## What this file is not",
  "",
  "- **Not u*-filtered.** The paper applies a per-year u* threshold",
  "  (0.076-0.157 m s-1, pooled fallback 0.115) before analysis. See the",
  "  analysis code.",
  "- **Not the derived analysis series.** Season and ENSO labels, the screened",
  "  NEE, and the Reco/GEP partitioning are produced by the analysis code from",
  "  this file.",
  "- **Not the full instrument record.** Only the CO2, energy and meteorological",
  "  variables used by this paper are included.",
  "",
  "## Related",
  "",
  "- Analysis code: https://github.com/rafastern/tambopata-co2-gcb-2026",
  "- Aggregated products behind the figures and tables: separate Zenodo record",
  "- Curated long-term record: AmeriFlux site PE-TNR,",
  "  https://ameriflux.lbl.gov/sites/siteinfo/PE-TNR",
  "",
  "Generated by `code/zz_build_halfhourly_release.R`."
), file.path(out_dir, "README.md"))

cat("\nHalf-hourly release assembled at:", out_dir, "\n")
cat("  ", csv_name, sprintf("(%.1f MB)", file.size(csv_path) / 1024^2), "\n")
cat("   columns  :", ncol(out), "\n")
cat("   rows     :", format(nrow(out), big.mark = ","), "\n")
cat("   finite NEE:", format(n_meas, big.mark = ","), "\n\n")
cat("Next: upload the contents of halfhourly_release/ to Zenodo as a dataset record.\n")
