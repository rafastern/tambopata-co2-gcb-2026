# -----------------------------------------------------------------------------
# u* (friction velocity) filtering for the measured (non-gap-filled) NEE trace.
#
# Why this exists:
#   The MATLAB export (from Matlab/Export_data.m) builds the `NEE` column with a
#   priority + fallback:
#       1) NEE_PI_SC_JSZ_MAD_RP_uStar_orig  -- u*-filtered, non-gap-filled
#       2) fallback NEE_PI_SC_JSZ_MAD        -- NOT u*-filtered
#   REddyProc only produced (1) for 2018, 2021, 2022 and 2024. For 2017, 2019,
#   2020, 2023 (and 2025) it crashed / returned NA, so the export fell back to
#   (2). The raw `NEE` column is therefore u*-filtered for some years and not
#   others -- and the un-filtered years are concentrated in El Nino, which
#   biases the ENSO comparison (low-u* calm nights lower nighttime NEE).
#
# What this does:
#   Re-imposes a consistent u* threshold in R, using each REddyProc-valid year's
#   own estimate and a pooled fallback for the failed years. For the already-
#   filtered years this is a no-op (their low-u* points are already NA); the
#   only NEW filtering happens on the failed years via the pooled fallback.
#
# Scope: NEE only. Any other measured flux sharing the same fallback would need
#   the same treatment if it is ever used.
# -----------------------------------------------------------------------------

# Per-year u* thresholds (m s-1) = REddyProc `uStar` estimate ("Ustar filtering
# Th_1"), read from from Matlab/PETNR_ThirdStageCleaning<year>.log.
USTAR_THRESHOLDS_BY_YEAR <- c(
  "2018" = 0.0947,
  "2021" = 0.1565,
  "2022" = 0.1343,
  "2024" = 0.0764
)

# Pooled fallback for the REddyProc-failed years (2017, 2019, 2020, 2023, 2025):
# median of the four valid per-year estimates (mean ~ 0.116). Change here if a
# genuine pooled bootstrap estimate is preferred.
USTAR_POOLED_FALLBACK <- 0.115

# Threshold for each element of `years` (falls back where the year is unlisted).
ustar_threshold_for_year <- function(years) {
  thr <- unname(USTAR_THRESHOLDS_BY_YEAR[as.character(years)])
  thr[is.na(thr)] <- USTAR_POOLED_FALLBACK
  thr
}

# Apply the u* filter: set `nee_col` to NA where u* is below the year's
# threshold and (if drop_ustar_na) where u* is missing. Returns the data frame
# with `nee_col` modified. Emits a per-year removal report when verbose.
apply_ustar_filter <- function(df,
                               nee_col = "NEE",
                               ustar_col = "USTAR",
                               datetime_col = "tv_dt",
                               year_vec = NULL,
                               drop_ustar_na = TRUE,
                               verbose = TRUE) {
  if (!nee_col %in% names(df))
    stop("apply_ustar_filter: missing NEE column '", nee_col, "'")
  if (!ustar_col %in% names(df))
    stop("apply_ustar_filter: missing u* column '", ustar_col, "'")

  # derive the year for each row if not supplied
  if (is.null(year_vec)) {
    if (!datetime_col %in% names(df))
      stop("apply_ustar_filter: need year_vec or datetime_col '", datetime_col, "'")
    dt <- suppressWarnings(lubridate::dmy_hms(df[[datetime_col]], quiet = TRUE))
    if (all(is.na(dt))) dt <- suppressWarnings(lubridate::ymd_hms(df[[datetime_col]], quiet = TRUE))
    if (all(is.na(dt))) dt <- suppressWarnings(as.POSIXct(df[[datetime_col]]))
    year_vec <- lubridate::year(dt)
  }

  nee   <- suppressWarnings(as.numeric(df[[nee_col]]))
  ustar <- suppressWarnings(as.numeric(df[[ustar_col]]))
  thr   <- ustar_threshold_for_year(year_vec)

  below   <- is.finite(nee) & is.finite(ustar) & (ustar < thr)
  ustarna <- is.finite(nee) & !is.finite(ustar)
  remove  <- below | (drop_ustar_na & ustarna)

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\napply_ustar_filter on '%s' (pooled fallback = %.4f m s-1, drop_ustar_na = %s):\n",
      nee_col, USTAR_POOLED_FALLBACK, drop_ustar_na))
    for (y in sort(unique(year_vec[!is.na(year_vec)]))) {
      idx <- !is.na(year_vec) & year_vec == y
      fn  <- sum(is.finite(nee[idx]))
      nb  <- sum(below[idx])
      na_ <- sum(ustarna[idx])
      rm_ <- nb + (if (drop_ustar_na) na_ else 0L)
      cat(sprintf(
        "  %s: thr=%.4f  finiteNEE=%5d  below=%5d  ustarNA=%5d  removed=%5d  kept=%5d\n",
        y, ustar_threshold_for_year(y), fn, nb, na_, rm_, fn - rm_))
    }
    cat(sprintf("  TOTAL: finiteNEE=%d  removed=%d  kept=%d\n",
                sum(is.finite(nee)), sum(remove), sum(is.finite(nee)) - sum(remove)))
  }

  nee[remove] <- NA_real_
  df[[nee_col]] <- nee
  df
}
