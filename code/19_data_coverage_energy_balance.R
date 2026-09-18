#!/usr/bin/env Rscript

# =============================================================================
# 19_data_coverage_energy_balance.R
# -----------------------------------------------------------------------------
# Per-year data availability (after QA/QC and u* filtering) and energy-balance
# closure. Requested by a co-author as a likely reviewer question.
#
# Two blocks per year:
#
#  (1) DATA AVAILABILITY
#      - measured (non-gap-filled) NEE after QA/QC + u* filtering  <- the honest
#        "how much real data" number; day/night split reported because the u*
#        filter removes calm nights preferentially.
#      - gap-filled NEE_ok coverage (what most analyses actually use).
#      - the u* threshold applied that year, and whether it is the year's own
#        REddyProc estimate or the pooled fallback (see code/ustar_filter.R).
#
#  (2) ENERGY-BALANCE CLOSURE, reported BOTH ways so the effect of the soil
#      term is explicit rather than left to the reader:
#      - against Rn:      EBR = sum(H + LE) / sum(Rn)
#      - against Rn - G:  EBR = sum(H + LE) / sum(Rn - G)
#      plus the OLS slope/intercept/R2 of (H + LE) on each, daytime and all-data.
#
#      G is the mean of two co-located replicate plates (G_1_1_1, G_1_1_2) at
#      10 cm depth, ~20 m west of the tower and ~1 m apart. The mean is taken
#      HERE rather than using the exported `G` column, because that column is
#      masked to half-hours where BOTH plates are valid and therefore drops
#      2020-2022 entirely (16k vs 42k half-hours). The plates agree closely
#      (r ~ 0.98, mean difference ~0.15 W m-2), so plate 1 alone represents the
#      pair well for the years before plate 2 was installed; n_plates records
#      how many plates contributed to each half-hour.
#
#      CAVEAT: the plates measure the flux at 10 cm, not at the surface. No
#      0-10 cm heat-storage correction is applied, so sub-daily G is damped and
#      lagged relative to the surface flux; over a year sum(dS) ~ 0, so annual
#      sums are essentially unaffected. In practice G is only ~0.2% of Rn at
#      this closed-canopy site, so including it moves the closure by only ~0.005.
#
# Run: Rscript code/19_data_coverage_energy_balance.R
# Outputs: output/table_data_coverage_energy_balance.csv   (full precision, all columns)
#          output/table_data_coverage_energy_balance.md    (markdown + captions)
#          output/table_data_coverage_energy_balance.html  (open in a browser, select
#            all, copy, paste into Google Docs -> real tables with the italics,
#            subscripts and superscripts already correct)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

in_fp <- file.path(input_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
if (!file.exists(in_fp)) stop("missing input: ", in_fp)

df <- readr::read_csv(in_fp, show_col_types = FALSE, progress = FALSE)

## ---- required columns -------------------------------------------------------
need <- c("DateTime", "NEE", "NEE_f", "NEE_ok", "H", "LE", "Rn", "SW_IN", "USTAR")
miss <- setdiff(need, names(df))
if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))

## ---- soil heat flux: two-plate mean, computed here (see header) -------------
g_cols <- intersect(c("G_1_1_1", "G_1_1_2"), names(df))
have_G <- length(g_cols) > 0
if (have_G) {
  gm <- as.matrix(df[, g_cols, drop = FALSE])
  gm[!is.finite(gm)] <- NA_real_
  df$n_plates <- rowSums(!is.na(gm))
  df$G_mean   <- ifelse(df$n_plates > 0, rowMeans(gm, na.rm = TRUE), NA_real_)
  cat("\nsoil heat flux: using", paste(g_cols, collapse = " + "),
      "->", sum(is.finite(df$G_mean)), "half-hours with G\n")
  if (length(g_cols) == 2) {
    ok <- is.finite(gm[, 1]) & is.finite(gm[, 2])
    cat(sprintf("  plate agreement over %d shared half-hours: r = %.3f, mean(p1 - p2) = %.2f W m-2\n",
                sum(ok), stats::cor(gm[ok, 1], gm[ok, 2]), mean(gm[ok, 1] - gm[ok, 2])))
  }
} else {
  df$n_plates <- 0L
  df$G_mean <- NA_real_
  cat("\nsoil heat flux: no G_1_1_* columns found; reporting closure against Rn only\n")
}

df <- df %>%
  mutate(
    DateTime = as.POSIXct(DateTime, tz = "America/Lima"),
    year = lubridate::year(DateTime)
  ) %>%
  filter(!is.na(DateTime))

## ---- re-apply the shared u* filter (idempotent for already-filtered years) --
# Guarantees the counts below use exactly the paper's filter, whatever state the
# stored NEE column is in.
df <- apply_ustar_filter(
  df, nee_col = "NEE", ustar_col = "USTAR",
  datetime_col = "DateTime", year_vec = df$year,
  drop_ustar_na = TRUE, verbose = TRUE
)

## ---- day / night classification (the paper's radiation criteria) -----------
# day   = PAR > 100 umol m-2 s-1 (equivalently Rg > 44 W m-2) together with Ta > 17 C
# night = Rg < 10 W m-2
# the intervening 10-44 W m-2 twilight band belongs to neither and is reported
# separately.
#
# Coverage is expressed as the COMPOSITION of each year's measured NEE (day,
# night and twilight as a share of the measured half-hours, summing to 100%).
# Using the criterion-meeting half-hours as the denominator instead would be
# misleading, because that denominator itself needs PAR/Ta/Rg to exist and so
# shrinks in poorly covered years: 2017 would read "68% daytime" for a year that
# is only 3.6% covered overall.
df <- df %>%
  mutate(
    is_day_p   = is.finite(PAR) & PAR > 100 & is.finite(Ta) & Ta > 17,
    is_night_p = is.finite(SW_IN) & SW_IN < 10,
    # used only to select well-lit half-hours for the energy-balance regression
    is_day_rg  = is.finite(SW_IN) & SW_IN > 44
  )

## ---- helper: OLS of (H+LE) on Rn -------------------------------------------
eb_fit <- function(rn, turb) {
  ok <- is.finite(rn) & is.finite(turb)
  if (sum(ok) < 30) return(list(slope = NA_real_, intercept = NA_real_, r2 = NA_real_, n = sum(ok)))
  m <- stats::lm(turb[ok] ~ rn[ok])
  list(
    slope     = unname(coef(m)[2]),
    intercept = unname(coef(m)[1]),
    r2        = summary(m)$r.squared,
    n         = sum(ok)
  )
}

## ---- per-year summary -------------------------------------------------------
years <- sort(unique(df$year))

rows <- lapply(years, function(y) {
  d <- df[df$year == y, , drop = FALSE]

  # half-hours actually spanned by the record within this year (the record ends
  # mid-October 2024, so a plain 17 520 denominator would understate 2024)
  n_rows_in_record <- nrow(d)
  n_hh_full_year <- if (lubridate::leap_year(y)) 366 * 48 else 365 * 48

  meas_ok  <- is.finite(d$NEE)
  gf_ok    <- is.finite(d$NEE_ok)
  turb     <- d$H + d$LE
  eb_ok    <- is.finite(turb) & is.finite(d$Rn)

  f_all <- eb_fit(d$Rn, turb)
  f_day <- eb_fit(d$Rn[d$is_day_rg], turb[d$is_day_rg])

  # same, but against available energy (Rn - G)
  avail    <- d$Rn - d$G_mean
  eb_ok_G  <- is.finite(turb) & is.finite(avail)
  fG_day   <- eb_fit(avail[d$is_day_rg], turb[d$is_day_rg])
  # matched baseline: Rn-only closure restricted to the SAME half-hours that have
  # G, so the Rn vs Rn-G comparison isolates the soil term instead of also mixing
  # in a different sample of half-hours.
  dayG     <- d$is_day_rg & is.finite(d$G_mean)
  fRn_dayG <- eb_fit(d$Rn[dayG], turb[dayG])

  # ---- the row actually reported in Table B -------------------------------
  # G counts as usable for a year only if it covers a meaningful share of that
  # year's energy-balance half-hours. 2020 fails this (129 of 2584, ~2.7 days),
  # so it is reported against Rn like the pre-installation years rather than on
  # a 2.7-day sample. Every number in a row then comes from one sample.
  g_usable <- sum(eb_ok_G) >= 500 && sum(eb_ok_G) >= 0.25 * max(1, sum(eb_ok))
  use      <- if (g_usable) eb_ok_G else eb_ok
  ref      <- if (g_usable) avail else d$Rn

  thr <- ustar_threshold_for_year(y)
  thr_src <- if (as.character(y) %in% names(USTAR_THRESHOLDS_BY_YEAR)) "per-year (REddyProc)" else "pooled fallback"

  tibble(
    year = y,
    ustar_threshold = thr,
    ustar_source = thr_src,
    n_hh_record = n_rows_in_record,
    n_hh_full_year = n_hh_full_year,
    # --- measured NEE after QA/QC + u* ---
    n_meas = sum(meas_ok),
    pct_meas_of_record = 100 * sum(meas_ok) / n_rows_in_record,
    pct_meas_of_year   = 100 * sum(meas_ok) / n_hh_full_year,
    n_meas_day      = sum(meas_ok & d$is_day_p),
    n_meas_night    = sum(meas_ok & d$is_night_p),
    n_meas_twilight = sum(meas_ok & !d$is_day_p & !d$is_night_p),
    pct_meas_day      = 100 * sum(meas_ok & d$is_day_p) / max(1, sum(meas_ok)),
    pct_meas_night    = 100 * sum(meas_ok & d$is_night_p) / max(1, sum(meas_ok)),
    pct_meas_twilight = 100 * sum(meas_ok & !d$is_day_p & !d$is_night_p) / max(1, sum(meas_ok)),
    # --- gap-filled coverage: as produced by REddyProc, and the subset that
    #     survives re-invalidation of runs of >24 consecutive missing half-hours ---
    n_gapfill_produced = sum(is.finite(d$NEE_f)),
    pct_gapfill_produced = 100 * sum(is.finite(d$NEE_f)) / n_rows_in_record,
    n_gapfilled = sum(gf_ok),
    pct_gapfilled_of_record = 100 * sum(gf_ok) / n_rows_in_record,
    # --- energy balance vs Rn ---
    n_eb = sum(eb_ok),
    EBR = sum(turb[eb_ok]) / sum(d$Rn[eb_ok]),
    slope_all = f_all$slope, intercept_all = f_all$intercept, r2_all = f_all$r2,
    slope_day = f_day$slope, intercept_day = f_day$intercept, r2_day = f_day$r2,
    n_eb_day = f_day$n,
    # --- soil heat flux + energy balance vs (Rn - G) ---
    n_G = sum(is.finite(d$G_mean)),
    n_G_two_plates = sum(d$n_plates == 2L),
    G_mean_day = if (any(d$is_day_rg & is.finite(d$G_mean))) mean(d$G_mean[d$is_day_rg], na.rm = TRUE) else NA_real_,
    n_eb_G = sum(eb_ok_G),
    EBR_Rn_Gsubset = if (sum(eb_ok_G)) sum(turb[eb_ok_G]) / sum(d$Rn[eb_ok_G]) else NA_real_,
    EBR_RnG = if (sum(eb_ok_G)) sum(turb[eb_ok_G]) / sum(avail[eb_ok_G]) else NA_real_,
    slope_day_Rn_Gsubset = fRn_dayG$slope, slope_day_RnG = fG_day$slope,
    r2_day_RnG = fG_day$r2, n_eb_day_RnG = fG_day$n,
    # --- the reported row: one sample, all means from it ---
    g_usable = g_usable,
    n_report = sum(use),
    turb_bar = if (sum(use)) mean(turb[use]) else NA_real_,
    Rn_bar   = if (sum(use)) mean(d$Rn[use]) else NA_real_,
    G_bar    = if (g_usable && sum(use)) mean(d$G_mean[use]) else NA_real_,
    ratio    = if (sum(use)) sum(turb[use]) / sum(ref[use]) else NA_real_
  )
})

tab <- bind_rows(rows)

## ---- all-years row ----------------------------------------------------------
turb_all <- df$H + df$LE
eb_ok_all <- is.finite(turb_all) & is.finite(df$Rn)
f_all_pooled <- eb_fit(df$Rn, turb_all)
f_day_pooled <- eb_fit(df$Rn[df$is_day_rg], turb_all[df$is_day_rg])
avail_all     <- df$Rn - df$G_mean
eb_ok_G_all   <- is.finite(turb_all) & is.finite(avail_all)
fG_day_pooled <- eb_fit(avail_all[df$is_day_rg], turb_all[df$is_day_rg])
dayG_all        <- df$is_day_rg & is.finite(df$G_mean)
fRn_dayG_pooled <- eb_fit(df$Rn[dayG_all], turb_all[dayG_all])

overall <- tibble(
  year = NA_integer_,
  ustar_threshold = NA_real_,
  ustar_source = "all years pooled",
  n_hh_record = nrow(df),
  n_hh_full_year = sum(tab$n_hh_full_year),
  n_meas = sum(is.finite(df$NEE)),
  pct_meas_of_record = 100 * sum(is.finite(df$NEE)) / nrow(df),
  pct_meas_of_year = 100 * sum(is.finite(df$NEE)) / sum(tab$n_hh_full_year),
  n_meas_day      = sum(is.finite(df$NEE) & df$is_day_p),
  n_meas_night    = sum(is.finite(df$NEE) & df$is_night_p),
  n_meas_twilight = sum(is.finite(df$NEE) & !df$is_day_p & !df$is_night_p),
  pct_meas_day      = 100 * sum(is.finite(df$NEE) & df$is_day_p) / max(1, sum(is.finite(df$NEE))),
  pct_meas_night    = 100 * sum(is.finite(df$NEE) & df$is_night_p) / max(1, sum(is.finite(df$NEE))),
  pct_meas_twilight = 100 * sum(is.finite(df$NEE) & !df$is_day_p & !df$is_night_p) / max(1, sum(is.finite(df$NEE))),
  n_gapfill_produced = sum(is.finite(df$NEE_f)),
  pct_gapfill_produced = 100 * sum(is.finite(df$NEE_f)) / nrow(df),
  n_gapfilled = sum(is.finite(df$NEE_ok)),
  pct_gapfilled_of_record = 100 * sum(is.finite(df$NEE_ok)) / nrow(df),
  n_eb = sum(eb_ok_all),
  EBR = sum(turb_all[eb_ok_all]) / sum(df$Rn[eb_ok_all]),
  slope_all = f_all_pooled$slope, intercept_all = f_all_pooled$intercept, r2_all = f_all_pooled$r2,
  slope_day = f_day_pooled$slope, intercept_day = f_day_pooled$intercept, r2_day = f_day_pooled$r2,
  n_eb_day = f_day_pooled$n,
  n_G = sum(is.finite(df$G_mean)),
  n_G_two_plates = sum(df$n_plates == 2L),
  G_mean_day = mean(df$G_mean[df$is_day_rg], na.rm = TRUE),
  n_eb_G = sum(eb_ok_G_all),
  EBR_Rn_Gsubset = if (sum(eb_ok_G_all)) sum(turb_all[eb_ok_G_all]) / sum(df$Rn[eb_ok_G_all]) else NA_real_,
  EBR_RnG = if (sum(eb_ok_G_all)) sum(turb_all[eb_ok_G_all]) / sum(avail_all[eb_ok_G_all]) else NA_real_,
  slope_day_Rn_Gsubset = fRn_dayG_pooled$slope, slope_day_RnG = fG_day_pooled$slope,
  r2_day_RnG = fG_day_pooled$r2, n_eb_day_RnG = fG_day_pooled$n,
  g_usable = FALSE,
  n_report = sum(eb_ok_all),
  turb_bar = mean(turb_all[eb_ok_all]),
  Rn_bar   = mean(df$Rn[eb_ok_all]),
  G_bar    = NA_real_,
  ratio    = sum(turb_all[eb_ok_all]) / sum(df$Rn[eb_ok_all])
)

tab_out <- bind_rows(tab, overall)

csv_fp <- file.path(out_base, "table_data_coverage_energy_balance.csv")
readr::write_csv(tab_out, csv_fp)

## ---- manuscript-ready markdown ---------------------------------------------
fmt <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
fmt0 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "d", big.mark = " "))

## ---- manuscript-ready markdown: two tables ---------------------------------
fmt <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
fmt0 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "d", big.mark = " "))
yr_label <- function(y) if (is.na(y)) "**All years**" else as.character(y)
ustar_cell <- function(r) ifelse(is.na(r$ustar_threshold), "--",
  paste0(fmt(r$ustar_threshold, 4), ifelse(grepl("pooled f", r$ustar_source), "*", "")))

## Table A: data availability
mdA <- c(
  "### Table A. Data availability per year, after QA/QC and u* filtering",
  "",
  "| Year | u* threshold (m s^-1) | Half-hours in record | Measured NEE after QA/QC + u* (half-hours) | % of record | Daytime (%) | Nighttime (%) | Twilight (%) | Gap-filled NEE produced (%) | Usable after long-gap QC (%) |",
  "|---|---|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(tab_out))) {
  r <- tab_out[i, ]
  mdA <- c(mdA, paste0(
    "| ", yr_label(r$year),
    " | ", ustar_cell(r),
    " | ", fmt0(r$n_hh_record),
    " | ", fmt0(r$n_meas),
    " | ", fmt(r$pct_meas_of_record, 1),
    " | ", fmt(r$pct_meas_day, 1),
    " | ", fmt(r$pct_meas_night, 1),
    " | ", fmt(r$pct_meas_twilight, 1),
    " | ", fmt(r$pct_gapfill_produced, 1),
    " | ", fmt(r$pct_gapfilled_of_record, 1),
    " |"))
}
mdA <- c(mdA, "",
  "*Pooled fallback u* threshold (REddyProc did not converge for that year); see code/ustar_filter.R.",
  "",
  "Percentages are relative to the half-hours actually spanned by the processed record in each year (the",
  "record ends 15 October 2024, so 2024 is a partial year). The daytime, nighttime and twilight columns",
  "give the composition of that year's measured NEE using the same radiation criteria as the rest of the",
  "study: daytime is PAR > 100 umol m-2 s-1 (equivalently Rg > 44 W m-2) together with Ta > 17 C,",
  "nighttime is Rg < 10 W m-2, and twilight is the intervening 10-44 W m-2 band, which is excluded from",
  "both the daytime and nighttime analyses. The three columns sum to 100 %. Nighttime is a much smaller",
  "share of the retained data than daytime because the friction-velocity filter preferentially removes",
  "calm nights. This imbalance does not propagate into the flux analyses: rather than using the raw",
  "half-hourly observations, those analyses are based on mean diel composites of 48 half-hourly bins, in",
  "which each half-hour of the day contributes one mean value regardless of how many days were sampled,",
  "subject to minimum-coverage thresholds (Sections 2.6 and 2.9).",
  "",
  "The two gap-filling columns separate two distinct limitations. *Gap-filled NEE produced* is the",
  "coverage of the REddyProc gap-filled series as generated: it is 0 % for 2017, 2019, 2020 and 2023,",
  "years for which REddyProc did not converge and no gap-filled series exists at all (the same years",
  "that use the pooled fallback u* threshold). *Usable after long-gap QC* is the subset retained after",
  "re-invalidating any run of more than 24 consecutive missing half-hours (more than 12 h) in the",
  "measured series, since gap-filling across such a gap is not defensible. The two together explain why",
  "year-integrated carbon budgets are not computed from the gap-filled series: in 2021, for example,",
  "gap-filling covered the whole year but only 42.9 % of it survives that criterion."
)

## Table B: energy-balance closure
## Deliberately framed as the physical comparison H + LE vs Rn - G rather than as an
## abstract ratio: showing the fluxes in W m-2 lets the reader verify the ratio and
## makes it obvious why G barely matters here (a few W m-2 against an Rn of ~250).
## One sample per row, so the half-hour count applies to every number in it.
mdB <- c(
  "### Table B. Energy-balance closure per year",
  "",
  "| Year | Half-hours | H + LE (W m-2) | Rn (W m-2) | G (W m-2) | (H + LE) / (Rn - G) |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(tab_out))) {
  r <- tab_out[i, ]
  mdB <- c(mdB, paste0(
    "| ", yr_label(r$year),
    " | ", fmt0(r$n_report),
    " | ", fmt(r$turb_bar, 1),
    " | ", fmt(r$Rn_bar, 1),
    " | ", if (is.na(r$G_bar)) "--" else fmt(r$G_bar, 2),
    " | ", fmt(r$ratio, 3),
    " |"))
}
mdB <- c(mdB, "",
  "Each row compares the turbulent fluxes (sensible plus latent heat, H + LE) with the energy available to",
  "drive them (net radiation minus the soil heat flux, Rn - G), over the half-hours in which every term in",
   "that row is present; the half-hour count therefore applies to all values in the row. H + LE, Rn and G",
  "are means over those half-hours, so the final column can be checked directly against them. That final",
  "ratio is conventionally called the energy-balance ratio (EBR); a value below 1 means the turbulent",
  "fluxes account for less than the available energy, as is typical of eddy-covariance sites.",
  "",
  "G is the mean of two co-located replicate plates at 10 cm depth, ~20 m west of the tower and ~1 m apart,",
  "which agree closely (r ~ 0.98, mean difference ~0.16 W m-2); 2021-2022 rests on plate 1 alone, as plate",
  "2 was installed in 2023. A dash means no usable soil heat flux for that year and the ratio is then",
  "H + LE over Rn alone: the plates were installed in 2020, and 2020 itself has G for only 129 half-hours",
  "(about 2.7 days), too few to represent the year, so it is reported against Rn as well.",
  "",
  "Including G changes the closure very little, because at this closed-canopy site G is only a few W m-2",
  "against an Rn of order 250 W m-2: over the 16 619 half-hours for which G is available the ratio is 0.733",
  "with G included and 0.728 without it. Because the plates sit at 10 cm rather than at the surface and no",
  "0-10 cm heat-storage correction is applied, sub-daily G is damped and lagged relative to the surface",
  "flux, though the storage term integrates to about zero over a year. An ordinary least-squares fit of",
  "(H + LE) on Rn over daytime half-hours (Rg > 44 W m-2) gives a slope of 0.68 and R2 = 0.76 for the whole",
  "record; per-year slopes, R2 values and two-plate counts are in the accompanying CSV."
)


## ---- HTML version -----------------------------------------------------------
## Markdown pasted into Google Docs arrives as literal text (asterisks, "m-2"),
## so also emit HTML: opened in a browser and copied, it pastes into Docs as a
## real table with the italics, subscripts and superscripts already correct.
sup  <- function(x) paste0("<sup>", x, "</sup>")
sub_ <- function(x) paste0("<sub>", x, "</sub>")
ital <- function(x) paste0("<i>", x, "</i>")
MINUS <- "−"          # true minus sign, not a hyphen
EMDASH <- "—"
# typography helpers for the values
nice <- function(x) ifelse(is.na(x) | x == "--", EMDASH, sub("^-", MINUS, x))

W_M2   <- paste0("W m", sup(paste0(MINUS, "2")))
MS1    <- paste0("m s", sup(paste0(MINUS, "1")))
Rn_h   <- paste0(ital("R"), sub_("n"))
G_h    <- ital("G")
ustar_h <- paste0(ital("u"), "*")

th <- function(x) paste0('<th style="border:1px solid #999;padding:5px 8px;
  text-align:center;vertical-align:bottom;background:#f2f2f2">', x, "</th>")
td <- function(x, align = "right") paste0('<td style="border:1px solid #999;',
  'padding:4px 8px;text-align:', align, '">', x, "</td>")

html_table <- function(title, headers, rows) {
  c(paste0("<p><b>", title, "</b></p>"),
    '<table style="border-collapse:collapse;font-family:Calibri,Arial,sans-serif;font-size:10pt">',
    "<tr>", vapply(headers, th, ""), "</tr>",
    unlist(lapply(rows, function(r) c("<tr>",
      td(r[1], "left"), vapply(r[-1], td, "", align = "right"), "</tr>"))),
    "</table>")
}

## Table A rows
rowsA <- lapply(seq_len(nrow(tab_out)), function(i) {
  r <- tab_out[i, ]
  yr <- if (is.na(r$year)) "<b>All years</b>" else as.character(r$year)
  c(yr, nice(ustar_cell(r)), nice(fmt0(r$n_hh_record)), nice(fmt0(r$n_meas)),
    nice(fmt(r$pct_meas_of_record, 1)), nice(fmt(r$pct_meas_day, 1)),
    nice(fmt(r$pct_meas_night, 1)), nice(fmt(r$pct_meas_twilight, 1)),
    nice(fmt(r$pct_gapfill_produced, 1)), nice(fmt(r$pct_gapfilled_of_record, 1)))
})
headA <- c("Year", paste0(ustar_h, " threshold (", MS1, ")"), "Half-hours in record",
           paste0("Measured NEE after QA/QC + ", ustar_h, " (half-hours)"),
           "% of record", "Daytime (%)", "Nighttime (%)", "Twilight (%)",
           "Gap-filled NEE produced (%)", "Usable after long-gap QC (%)")

## Table B rows
rowsB <- lapply(seq_len(nrow(tab_out)), function(i) {
  r <- tab_out[i, ]
  yr <- if (is.na(r$year)) "<b>All years</b>" else as.character(r$year)
  c(yr, nice(fmt0(r$n_report)), nice(fmt(r$turb_bar, 1)), nice(fmt(r$Rn_bar, 1)),
    if (is.na(r$G_bar)) EMDASH else nice(fmt(r$G_bar, 2)), nice(fmt(r$ratio, 3)))
})
headB <- c("Year", "Half-hours", paste0(ital("H"), " + ", ital("LE"), " (", W_M2, ")"),
           paste0(Rn_h, " (", W_M2, ")"), paste0(G_h, " (", W_M2, ")"),
           paste0("(", ital("H"), " + ", ital("LE"), ") / (", Rn_h, " ", MINUS, " ", G_h, ")"))

html <- c("<!DOCTYPE html>", '<html><head><meta charset="utf-8"></head><body>',
  html_table("Table S2. Data availability per year, after QA/QC and u* filtering.", headA, rowsA),
  "<p><br></p>",
  html_table("Table S3. Energy-balance closure per year.", headB, rowsB),
  "<p style='font-size:9pt;color:#444'>Select all and copy this page, then paste into",
  "Google Docs: the tables arrive with formatting intact. Captions are kept separately",
  "in the .md file.</p>",
  "</body></html>")
html_fp <- file.path(out_base, "table_data_coverage_energy_balance.html")
writeLines(html, html_fp, useBytes = TRUE)

md <- c(mdA, "", mdB)
md_fp <- file.path(out_base, "table_data_coverage_energy_balance.md")
writeLines(md, md_fp, useBytes = TRUE)

## ---- console report ---------------------------------------------------------
cat("\n================ per-year data coverage & energy balance ================\n")
print(as.data.frame(tab_out %>% transmute(
  year, ustar = round(ustar_threshold, 4), src = substr(ustar_source, 1, 12),
  n_meas, pct_meas = round(pct_meas_of_record, 1),
  pct_day = round(pct_meas_day, 1), pct_night = round(pct_meas_night, 1),
  pct_gf = round(pct_gapfilled_of_record, 1),
  EBR = round(EBR, 3), slope_day = round(slope_day, 3), r2_day = round(r2_day, 3)
)), row.names = FALSE)

cat("\nwrote:\n  ", csv_fp, "\n  ", md_fp, "\n")
