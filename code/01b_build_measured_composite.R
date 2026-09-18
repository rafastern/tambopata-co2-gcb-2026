# =============================================================================
# 01b_build_measured_composite.R
# -----------------------------------------------------------------------------
# Build the measured-only (non-gap-filled) 48-point monthly diel composite.
#
# Purpose: give the driver/relationship analyses a NEE response built from the
#   raw non-gap-filled trace (column NEE) instead of the gap-filled NEE_ok.
#
# Design: mirrors the export_diel_month grouping and mean aggregation in
#   code/04, but filters on finite raw NEE (not NEE_ok). This recovers the 15
#   whole months where the REddyProc gap-filled trace NEE_f is empty but raw
#   NEE exists (mostly El Nino).
#
# Output: data/tambopata_48points_per_month_measured.csv -- a SEPARATE file so
#   the gap-filled composite (used by the carbon balance and surface
#   conductance) is left untouched.
#
# Note: GEP/Reco/Gs are NOT included here. They are partitioned/modelled from
#   the gap-filled series and remain the responsibility of code/04.
#
# Provenance: this script previously lived outside code/ in a separate working
#   directory. It is the sole producer of the composite that Figures 4, 5 and 6
#   and Figures S14-S16 and S18 are built from, so it was consolidated here to
#   make those figures reproducible from the repository alone.
#   Run it after code/01 and before code/05, 09, 11, 26 and 27.
# =============================================================================

suppressWarnings(suppressMessages({
  library(dplyr); library(readr); library(lubridate)
}))

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

enriched_fp <- file.path(input_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
if (!file.exists(enriched_fp))
  stop("missing input: ", enriched_fp,
       "\nRun code/01_prepare_flux_timeseries.R first.")

df <- read_csv(enriched_fp, show_col_types = FALSE)
stopifnot(all(c("NEE", "PAR", "season", "ENSO", "DateTime") %in% names(df)))

halfhour_levels <- format(
  seq(as.POSIXct("2000-01-01 00:00:00", tz = "UTC"), by = "30 min", length.out = 48), "%H:%M:%S")

# optional columns carried through if present (env predictors used by 05/09/11)
opt <- intersect(c("TA_1_1_1", "VPD_kPa", "SWC_1_1_1", "TS_3", "USTAR", "WS_1_1_1",
                   "PAR_corrected_SWin"), names(df))

measured <- df %>%
  filter(!is.na(DateTime), is.finite(NEE)) %>%
  mutate(ym   = format(DateTime, "%Y-%m"),
         hour = factor(format(floor_date(DateTime, "30 minutes"), "%H:%M:%S"), levels = halfhour_levels)) %>%
  group_by(ym, season, ENSO, hour) %>%
  summarise(
    NEE_meas = mean(NEE, na.rm = TRUE),
    PAR      = mean(PAR, na.rm = TRUE),
    across(all_of(opt), ~ mean(.x, na.rm = TRUE)),
    n_meas   = dplyr::n(),
    .groups  = "drop"
  ) %>%
  arrange(ym, hour)

out_fp <- file.path(output_path, "tambopata_48points_per_month_measured.csv")
write_csv(measured, out_fp)

cat("saved:", out_fp, "\n")
cat("rows:", nrow(measured), "| months:", dplyr::n_distinct(measured$ym),
    "| columns:", paste(names(measured), collapse = ", "), "\n")
cat("groups (season x ENSO) present:\n")
print(measured %>% distinct(season, ENSO) %>% arrange(season, ENSO) %>% as.data.frame())
