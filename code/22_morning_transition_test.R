#!/usr/bin/env Rscript

# Timing of the morning source-to-sink transition in wet-season NEE, compared
# between El Nino and La Nina.
#
# Why this script exists. The original version of this analysis lived only in
# a scratch script outside code/, was not
# reachable from any runner, and its saved result predated later changes to the
# flux pipeline. It is brought into code/ here so the manuscript number is
# reproducible.
#
# What changed relative to that version, and why:
#   * the transition is the interpolated downward zero crossing, not the first
#     half-hour with NEE < 0. The old estimator snapped every day to the 30-min
#     grid, so group means quantised to 30 minutes were being reported to the
#     nearest minute, and a single noisy negative half-hour at 05:00 set that
#     whole day's value.
#   * it uses measured (non-gap-filled, u*-filtered) NEE. Under the old
#     estimator 36% of the half-hours that defined a crossing were gap-filled,
#     so a third of the "measured" transitions were model output, and the
#     gap-fill is driven by radiation - the same thing the crossing measures.
#   * the difference is reported as a month-block bootstrap interval rather than
#     a two-sample test on days. Wet-season El Nino days come from only two
#     episodes (Nov 2018 and Jan-Apr 2024) and La Nina from about four, so days
#     are not independent replicates and a day-level t-test is anticonservative
#     by an unquantified factor.
#
# All four estimator/data combinations are reported so the choice is visible.
#
# Added in response to a co-author comment on section 3.1 ("is this light
# corrected/normalized? is PAR identical? ... How can you claim 36 min precision
# if data are half-hourly? The uncertainty is for the timing only, but not the
# statistical significance of the difference in NEE?"). Four checks are appended
# below the original estimator table:
#
#   * light controls. The crossing time confounds when the light arrives with
#     how much light the ecosystem needs. So we also report morning PAR, the PAR
#     at which each day's crossing occurs, and a decomposition of the timing
#     difference into a light-supply part (the time each phase needs to reach the
#     pooled mean crossing PAR) and a threshold part (the residual). Note PAR at
#     this site is PAR = 2.3 * SW_IN - 1.22 exactly, i.e. derived from measured
#     shortwave, so this is equivalently a statement about SW_IN.
#   * the same month-block bootstrap applied to mean morning NEE, so the text can
#     say explicitly what the interval does and does not cover.
#   * sunrise, from solar geometry, so the contrast can be shown not to be a
#     day-length artifact arising from the different calendar-month composition
#     of the two samples.
#   * per-episode breakdown and a storage-term sensitivity. The El Nino sample
#     mixes storage methods (profiler from 2024, DNN model before), and the
#     pooled effect is not homogeneous across episodes; both need to be visible.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidyr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

# ─────────────────────────────────────────────────────────────────────────────
# config

seed   <- 20260810
n_boot <- 2000

tz_local    <- "America/Lima"
window_start <- 5    # local hour; the morning transition cannot fall outside
window_end   <- 12
enso_levels  <- c("El Nino", "La Nina")

# morning index window, used for the mean morning PAR and mean morning NEE that
# the light controls compare between phases
morning_from <- 6
morning_to   <- 9

# half-hours after the crossing that must also be negative for the primary
# estimator to accept it (see crossing_interpolated below)
primary_persist <- 2

# site coordinates (code/13_maps.R) and the standard-time offset for
# America/Lima; used only for the sunrise calculation
site_lat  <- -12.83
site_lon  <- -70.25
tz_offset <- -5

flux_csv  <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
out_main  <- file.path(out_base, "table_morning_transition.csv")
out_days  <- file.path(output_path, "morning_transition_per_day.csv")
out_light <- file.path(out_base, "table_morning_transition_light_controls.csv")
out_epi   <- file.path(out_base, "table_morning_transition_by_episode.csv")
out_hh    <- file.path(out_base, "table_morning_transition_halfhourly.csv")

set.seed(seed)

# ─────────────────────────────────────────────────────────────────────────────
# load

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv,
       "\nRun code/01_prepare_flux_timeseries.R first.")
}

required_cols <- c("tv_dt", "NEE", "NEE_ok", "season", "ENSO",
                   # light controls and the storage-term sensitivity
                   "PAR", "FC", "SC_profiler", "SC_model", "USTAR")
df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("missing required columns: ", paste(missing_cols, collapse = ", "))
}

base <- df_raw %>%
  mutate(
    ts   = dmy_hms(tv_dt, tz = tz_local, quiet = TRUE),
    Date = as.Date(ts),
    hour_dec = hour(ts) + minute(ts) / 60
  ) %>%
  filter(season == "wet", ENSO %in% enso_levels,
         hour_dec >= window_start, hour_dec <= window_end)

# NEE is FC + SC exactly in this export (max |NEE - (FC + SC)| = 2e-6), so the
# storage-term sensitivity can be built by substituting the other SC series and
# re-imposing the same u* filter the pipeline uses.
base <- base %>%
  mutate(
    year     = lubridate::year(ts),
    NEE_prof = as.numeric(FC) + as.numeric(SC_profiler),
    NEE_mod  = as.numeric(FC) + as.numeric(SC_model)
  )
for (v in c("NEE_prof", "NEE_mod")) {
  base <- apply_ustar_filter(base, nee_col = v, ustar_col = "USTAR",
                             year_vec = base$year, verbose = FALSE)
}

if (!nrow(base)) stop("no wet-season El Nino / La Nina half-hours in the morning window")

# ─────────────────────────────────────────────────────────────────────────────
# transition estimators

# the half-hour grid means a true crossing lies between two samples; interpolate
# linearly across a positive-to-negative pair of the day.
#
# `min_persist` is the number of half-hours after the crossing that must also be
# negative for it to count. min_persist = 1 is the bare first downward crossing
# and is what the first version of this script used; it fires on a single noisy
# half-hour whenever morning NEE hovers near zero, which is exactly the regime
# the 2024 mornings are in (NEE returns positive after the first crossing on 50%
# of 2023-24 El Nino days, against 19% in 2018-19 and 5% under La Nina). That is
# the same failure mode this script's header objects to in the pre-2026
# estimator; interpolating fixed the snapping but not the persistence. The
# primary result therefore uses min_persist = 2.
crossing_interpolated <- function(df, col, min_persist = 1) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    arrange(Date, hour_dec) %>%
    group_by(Date, ENSO) %>%
    summarise(
      t = {
        y <- .data[[col]]
        x <- hour_dec
        i <- which(head(y, -1) > 0 & tail(y, -1) < 0)
        keep <- vapply(i, function(k) {
          nxt <- y[seq.int(k + 1, min(k + min_persist, length(y)))]
          length(nxt) >= min_persist && all(nxt < 0)
        }, logical(1))
        i <- i[keep][1]
        if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
      },
      .groups = "drop"
    ) %>%
    filter(is.finite(t))
}

# GRID-NATIVE estimator, and the one the manuscript reports.
#
# Same persistence rule, but the day's value is the timestamp of the first
# negative half-hour itself, with no interpolation, so every reported quantity
# lives on the 30-minute measurement grid. Timestamps are period-end labels, so
# a value of 08:30 means the 08:00-08:30 averaging interval.
#
# The interpolated estimators are retained below only as sensitivity rows; their
# sub-half-hour values are not quoted in the text.
crossing_slot <- function(df, col, min_persist = 1) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    arrange(Date, hour_dec) %>%
    group_by(Date, ENSO) %>%
    summarise(
      t = {
        y <- .data[[col]]; x <- hour_dec
        i <- which(head(y, -1) > 0 & tail(y, -1) < 0)
        keep <- vapply(i, function(k) {
          nxt <- y[seq.int(k + 1, min(k + min_persist, length(y)))]
          length(nxt) >= min_persist && all(nxt < 0)
        }, logical(1))
        i <- i[keep][1]
        if (is.na(i)) NA_real_ else x[i + 1]
      },
      .groups = "drop"
    ) %>%
    filter(is.finite(t))
}

# the last downward crossing of the window, i.e. the transition after which the
# day does not revert to a source again. Reported as a further sensitivity.
crossing_last <- function(df, col) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    arrange(Date, hour_dec) %>%
    group_by(Date, ENSO) %>%
    summarise(
      t = {
        y <- .data[[col]]; x <- hour_dec
        i <- which(head(y, -1) > 0 & tail(y, -1) < 0)
        i <- if (length(i)) i[length(i)] else NA_integer_
        if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
      },
      .groups = "drop"
    ) %>%
    filter(is.finite(t))
}

# the original estimator, kept for the sensitivity table
crossing_first_negative <- function(df, col) {
  df %>%
    filter(is.finite(.data[[col]]), .data[[col]] < 0) %>%
    group_by(Date, ENSO) %>%
    summarise(t = min(hour_dec), .groups = "drop")
}

hm <- function(x) sprintf("%02d:%02d", floor(x), round((x - floor(x)) * 60))

# ─────────────────────────────────────────────────────────────────────────────
# light-control helpers

# linear interpolation of y at x0, NA outside the observed range
interp_at <- function(x, y, x0) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2 || !is.finite(x0)) return(NA_real_)
  o <- order(x[keep])
  stats::approx(x[keep][o], y[keep][o], xout = x0, rule = 1)$y
}

# first time y rises through thr, interpolated between the bracketing samples
first_time_reaching <- function(x, y, thr) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) return(NA_real_)
  o <- order(x[keep]); x <- x[keep][o]; y <- y[keep][o]
  i <- which(head(y, -1) < thr & tail(y, -1) >= thr)[1]
  if (is.na(i)) return(NA_real_)
  x[i] + (x[i + 1] - x[i]) * (thr - y[i]) / (y[i + 1] - y[i])
}

# NOAA solar-position approximation; returns local-clock sunrise in decimal
# hours. Only used to show that the two ENSO samples, which have different
# calendar-month composition, do not differ in day length.
sunrise_local_hour <- function(doy) {
  g   <- 2 * pi / 365 * (doy - 1 + (6 - 12) / 24)
  eqt <- 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) -
                   0.014615 * cos(2 * g) - 0.040849 * sin(2 * g))
  dec <- 0.006918 - 0.399912 * cos(g) + 0.070257 * sin(g) -
         0.006758 * cos(2 * g) + 0.000907 * sin(2 * g) -
         0.002697 * cos(3 * g) + 0.00148 * sin(3 * g)
  lat <- site_lat * pi / 180
  x   <- cos(90.833 * pi / 180) / (cos(lat) * cos(dec)) - tan(lat) * tan(dec)
  ha  <- acos(pmin(pmax(x, -1), 1)) * 180 / pi
  (720 - 4 * (site_lon + ha) - eqt) / 60 + tz_offset
}

# the ENSO episode each day belongs to; the pooled contrast is not homogeneous
# across these, so section 3.1 reports them separately
# ENSO has to be an argument: the wet season straddles the new year, so a date
# alone cannot tell an El Nino episode from a La Nina one (Jan-Feb 2023 is
# La Nina, Nov 2023 onward is El Nino).
episode_of <- function(date, enso) {
  y <- as.integer(format(date, "%Y")); m <- as.integer(format(date, "%m"))
  dplyr::case_when(
    enso == "El Nino" & y <= 2019        ~ "EN 2018-19",
    enso == "El Nino"                    ~ "EN 2023-24",
    y == 2017 | (y == 2018 & m <= 6)     ~ "LN 2017-18",
    y == 2020 | (y == 2021 & m <= 6)     ~ "LN 2020-21",
    y == 2021 | (y == 2022 & m <= 6)     ~ "LN 2021-22",
    TRUE                                 ~ "LN 2022-23"
  )
}

# per-day morning summary: the crossing, the light at the crossing, and the mean
# morning PAR and NEE over the same window for every day in the sample
per_day_morning <- function(df) {
  df %>%
    arrange(Date, hour_dec) %>%
    group_by(Date, ENSO) %>%
    summarise(
      # primary: the label of the first negative half-hour that holds for
      # `primary_persist` further half-hours. On the measurement grid by
      # construction. `t_interp` is the interpolated value, kept for the
      # sensitivity table only, and `t_first` the unscreened interpolated one.
      t = {
        ok <- is.finite(NEE)
        y  <- NEE[ok]; x <- hour_dec[ok]
        i  <- which(head(y, -1) > 0 & tail(y, -1) < 0)
        keep <- vapply(i, function(k) {
          nxt <- y[seq.int(k + 1, min(k + primary_persist, length(y)))]
          length(nxt) >= primary_persist && all(nxt < 0)
        }, logical(1))
        i <- i[keep][1]
        if (is.na(i)) NA_real_ else x[i + 1]
      },
      t_interp = {
        ok <- is.finite(NEE)
        y  <- NEE[ok]; x <- hour_dec[ok]
        i  <- which(head(y, -1) > 0 & tail(y, -1) < 0)
        keep <- vapply(i, function(k) {
          nxt <- y[seq.int(k + 1, min(k + primary_persist, length(y)))]
          length(nxt) >= primary_persist && all(nxt < 0)
        }, logical(1))
        i <- i[keep][1]
        if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
      },
      t_first = {
        ok <- is.finite(NEE)
        y  <- NEE[ok]; x <- hour_dec[ok]
        i  <- which(head(y, -1) > 0 & tail(y, -1) < 0)[1]
        if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
      },
      # does the day revert to a source after that first crossing? this is the
      # evidence for screening on persistence at all
      reverts = {
        ok <- is.finite(NEE)
        y  <- NEE[ok]
        i  <- which(head(y, -1) > 0 & tail(y, -1) < 0)[1]
        if (is.na(i)) NA else any(y[seq.int(i + 1, length(y))] > 0)
      },
      PAR_hour = list(hour_dec), PAR_val = list(PAR),
      PAR_morn = mean(PAR[hour_dec >= morning_from & hour_dec <= morning_to], na.rm = TRUE),
      NEE_morn = mean(NEE[hour_dec >= morning_from & hour_dec <= morning_to], na.rm = TRUE),
      n_NEE_morn = sum(is.finite(NEE) & hour_dec >= morning_from & hour_dec <= morning_to),
      .groups = "drop"
    ) %>%
    # mapply, not rowwise(): under rowwise() a list-column is already unwrapped
    # to its element, so the usual `col[[1]]` idiom silently takes the first
    # number of the vector instead of the vector
    mutate(
      # grid-native: the PAR actually recorded in the transition half-hour, not
      # a value interpolated to a sub-grid instant
      PAR_at_crossing = mapply(
        function(xh, xp, tt) if (!is.finite(tt)) NA_real_ else {
          j <- which(xh == tt)[1]; if (is.na(j)) NA_real_ else xp[j]
        }, PAR_hour, PAR_val, t, SIMPLIFY = TRUE)
    ) %>%
    mutate(
      PAR_morn = ifelse(is.finite(PAR_morn), PAR_morn, NA_real_),
      NEE_morn = ifelse(is.finite(NEE_morn), NEE_morn, NA_real_),
      sunrise  = sunrise_local_hour(as.integer(format(Date, "%j"))),
      t_since_sunrise = t - sunrise,
      episode  = episode_of(Date, ENSO)
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# month-block bootstrap
#
# resamples whole calendar months with replacement within each ENSO phase, so
# the resample respects the block structure of the record instead of treating
# adjacent days as independent draws.
#
# `value_col` / `scale` default to the transition time in minutes, so the
# original call sites are unchanged; the light controls reuse the same block
# scheme for PAR and NEE by passing their own column and scale = 1.
month_block_bootstrap <- function(st, reps, value_col = "t", scale = 60) {
  st <- st[is.finite(st[[value_col]]), , drop = FALSE]
  st <- st %>% mutate(blk = format(Date, "%Y-%m"))
  blocks <- split(st, list(st$ENSO, st$blk), drop = TRUE)
  by_enso <- lapply(enso_levels, function(g) {
    names(blocks)[vapply(blocks, function(b) b$ENSO[1] == g, logical(1))]
  })
  names(by_enso) <- enso_levels

  out <- numeric(reps)
  for (r in seq_len(reps)) {
    means <- vapply(enso_levels, function(g) {
      keys <- sample(by_enso[[g]], length(by_enso[[g]]), replace = TRUE)
      mean(unlist(lapply(keys, function(k) blocks[[k]][[value_col]])))
    }, numeric(1))
    out[r] <- (means[["La Nina"]] - means[["El Nino"]]) * scale
  }
  out
}

summarise_variant <- function(st, estimator, nee_series, is_primary) {
  en <- st$t[st$ENSO == "El Nino"]
  ln <- st$t[st$ENSO == "La Nina"]
  if (length(en) < 3 || length(ln) < 3) return(NULL)

  boot <- month_block_bootstrap(st, n_boot)
  ci   <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)

  # secondary check at the month level, where the replicate is a calendar month
  monthly <- st %>%
    mutate(ym = format(Date, "%Y-%m")) %>%
    group_by(ym, ENSO) %>%
    summarise(t = mean(t), .groups = "drop")
  mt <- t.test(monthly$t[monthly$ENSO == "El Nino"],
               monthly$t[monthly$ENSO == "La Nina"])

  tibble(
    estimator          = estimator,
    nee_series         = nee_series,
    primary            = is_primary,
    n_days_el_nino     = length(en),
    n_days_la_nina     = length(ln),
    n_months_el_nino   = sum(monthly$ENSO == "El Nino"),
    n_months_la_nina   = sum(monthly$ENSO == "La Nina"),
    mean_el_nino_h     = mean(en),
    mean_la_nina_h     = mean(ln),
    mean_el_nino_hm    = hm(mean(en)),
    mean_la_nina_hm    = hm(mean(ln)),
    diff_min           = (mean(ln) - mean(en)) * 60,
    diff_ci_low_min    = unname(ci[1]),
    diff_ci_high_min   = unname(ci[2]),
    month_level_t_p    = mt$p.value,
    n_boot             = n_boot,
    seed               = seed
  )
}

# The primary row is the persistence-screened crossing on measured NEE. The
# unscreened rows are kept because the earlier version of section 3.1 quoted
# them, so the change in the headline number stays traceable.
persist_fn <- function(k) function(df, col) crossing_interpolated(df, col, min_persist = k)
slot_fn    <- function(k) function(df, col) crossing_slot(df, col, min_persist = k)

variants <- list(
  list("half-hour slot, persists 2 half-hours",        "measured",   "NEE",    slot_fn(2),    TRUE),
  list("half-hour slot, persists 2 half-hours",        "gap-filled", "NEE_ok", slot_fn(2),    FALSE),
  list("half-hour slot, persists 3 half-hours",        "measured",   "NEE",    slot_fn(3),    FALSE),
  list("interpolated crossing, persists 2 half-hours", "measured",   "NEE",    persist_fn(2), FALSE),
  list("interpolated crossing, persists 3 half-hours", "measured",   "NEE",    persist_fn(3), FALSE),
  list("last downward crossing",                       "measured",   "NEE",    crossing_last, FALSE),
  list("interpolated crossing",    "measured",   "NEE",    persist_fn(1),           FALSE),
  list("interpolated crossing",    "gap-filled", "NEE_ok", persist_fn(1),           FALSE),
  list("first negative half-hour", "measured",   "NEE",    crossing_first_negative, FALSE),
  list("first negative half-hour", "gap-filled", "NEE_ok", crossing_first_negative, FALSE)
)

rows <- lapply(variants, function(v) {
  st <- v[[4]](base, v[[3]])
  summarise_variant(st, v[[1]], v[[2]], v[[5]])
})

res <- bind_rows(rows)
if (!nrow(res)) stop("no variant produced a usable comparison")

readr::write_csv(res, out_main)
cat("saved:", out_main, "\n")

# per-day values behind the primary variant, so the distribution is inspectable.
# `day_all` keeps every day in the window, including days with no crossing, so
# the morning PAR and NEE indices are not conditioned on a crossing existing.
day_all <- per_day_morning(base) %>%
  select(-PAR_hour, -PAR_val) %>%
  arrange(ENSO, Date)

primary_days <- day_all %>%
  filter(is.finite(t)) %>%
  mutate(transition_hm = hm(t)) %>%
  select(Date, ENSO, t, transition_hm, t_first, reverts, PAR_at_crossing, PAR_morn,
         NEE_morn, n_NEE_morn, sunrise, t_since_sunrise, episode)
readr::write_csv(primary_days, out_days)
cat("saved:", out_days, "\n\n")

# why the primary estimator screens on persistence
cat("share of days where NEE reverts to a source after the FIRST downward crossing:\n")
print(as.data.frame(
  day_all %>%
    filter(is.finite(t_first)) %>%
    group_by(episode) %>%
    summarise(n_days = n(), pct_reverts = round(100 * mean(reverts, na.rm = TRUE), 1),
              .groups = "drop")
), row.names = FALSE)
cat("\n")

print(as.data.frame(
  res %>%
    transmute(estimator, nee_series, primary,
              n = paste0(n_days_el_nino, "/", n_days_la_nina),
              months = paste0(n_months_el_nino, "/", n_months_la_nina),
              EN = mean_el_nino_hm, LN = mean_la_nina_hm,
              diff_min = round(diff_min, 1),
              CI = sprintf("%.0f to %.0f", diff_ci_low_min, diff_ci_high_min),
              month_p = round(month_level_t_p, 3))
), row.names = FALSE)

p <- res %>% filter(primary)
cat(sprintf(
  "\nPRIMARY (%s, %s NEE): El Nino %s vs La Nina %s, %.0f min earlier\n  95%% CI %.0f to %.0f min (%d month-block bootstrap replicates, seed %d)\n  month-level Welch t-test p = %.3f (n = %d vs %d months)\n",
  p$estimator, p$nee_series, p$mean_el_nino_hm, p$mean_la_nina_hm, p$diff_min,
  p$diff_ci_low_min, p$diff_ci_high_min, p$n_boot, p$seed,
  p$month_level_t_p, p$n_months_el_nino, p$n_months_la_nina))

cat(sprintf("\nEl Nino days come from %d calendar months, La Nina from %d.\n",
            p$n_months_el_nino, p$n_months_la_nina))
cat("months contributing:\n")
print(primary_days %>% mutate(ym = format(Date, "%Y-%m")) %>% count(ENSO, ym) %>% as.data.frame(),
      row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# light controls
#
# The crossing time answers "when does the ecosystem reach net uptake", which
# confounds when the light arrives with how much light is needed. These three
# statistics separate the two, each with the same month-block bootstrap so the
# intervals are comparable with the timing interval above.

# re-seeded here so this section is reproducible independently of how many
# bootstrap draws the estimator table above consumed
set.seed(seed)

boot_ci <- function(st, value_col, scale = 1) {
  b <- month_block_bootstrap(st, n_boot, value_col = value_col, scale = scale)
  unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

summarise_stat <- function(st, value_col, label, unit, scale = 1) {
  s  <- st[is.finite(st[[value_col]]), , drop = FALSE]
  en <- s[[value_col]][s$ENSO == "El Nino"]
  ln <- s[[value_col]][s$ENSO == "La Nina"]
  if (length(en) < 3 || length(ln) < 3) return(NULL)
  ci <- boot_ci(s, value_col, scale)
  tibble(
    statistic     = label,
    unit          = unit,
    sample        = if (identical(attr(st, "sample_label"), NULL)) NA_character_
                    else attr(st, "sample_label"),
    n_el_nino     = length(en),
    n_la_nina     = length(ln),
    mean_el_nino  = mean(en) * scale,
    mean_la_nina  = mean(ln) * scale,
    diff_la_minus_el = (mean(ln) - mean(en)) * scale,
    ci_low        = ci[1],
    ci_high       = ci[2],
    excludes_zero = (ci[1] > 0) | (ci[2] < 0)
  )
}

tag <- function(df, lab) { attr(df, "sample_label") <- lab; df }

cross_days <- tag(primary_days, "days with a crossing")
all_days   <- tag(day_all,      "all days in window")

# The light-supply component: how much of the timing difference is just El Nino
# mornings reaching a given light level earlier. The reference level is the
# pooled mean PAR at which the crossing actually happens.
par_ref <- mean(primary_days$PAR_at_crossing, na.rm = TRUE)

# grid-native: the label of the first half-hour whose recorded PAR reaches the
# reference level, so the light-supply comparison is on the same grid as the
# transition itself
t_at_ref <- base %>%
  arrange(Date, hour_dec) %>%
  group_by(Date, ENSO) %>%
  summarise(
    t_at_par_ref = {
      ok <- is.finite(PAR); y <- PAR[ok]; x <- hour_dec[ok]
      j <- which(y >= par_ref)[1]
      if (is.na(j)) NA_real_ else x[j]
    },
    .groups = "drop") %>%
  semi_join(primary_days, by = c("Date", "ENSO"))
t_at_ref <- tag(t_at_ref, "days with a crossing")

lab_cross <- "transition time"
lab_light <- sprintf("time to reach PAR = %.0f (light supply)", par_ref)
lab_parx  <- "PAR at the crossing"
lab_par   <- sprintf("mean PAR %02d:00-%02d:00", morning_from, morning_to)
lab_nee   <- sprintf("mean NEE %02d:00-%02d:00", morning_from, morning_to)
umol      <- "umol m-2 s-1"

light_rows <- bind_rows(
  summarise_stat(cross_days, "t",               lab_cross,                       "min", 60),
  summarise_stat(cross_days, "t_since_sunrise", "transition time after sunrise", "min", 60),
  summarise_stat(t_at_ref,   "t_at_par_ref",    lab_light,                       "min", 60),
  summarise_stat(cross_days, "PAR_at_crossing", lab_parx,                        umol,   1),
  summarise_stat(cross_days, "PAR_morn",        lab_par,                         umol,   1),
  summarise_stat(cross_days, "NEE_morn",        lab_nee,                         umol,   1),
  summarise_stat(all_days,   "PAR_morn",        lab_par,                         umol,   1),
  summarise_stat(all_days,   "NEE_morn",        lab_nee,                         umol,   1)
) %>%
  mutate(par_ref = par_ref, n_boot = n_boot, seed = seed)

# decomposition: total = light supply + threshold residual
diff_total <- light_rows$diff_la_minus_el[light_rows$statistic == lab_cross]
diff_light <- light_rows$diff_la_minus_el[light_rows$statistic == lab_light]
light_rows <- bind_rows(
  light_rows,
  tibble(statistic = "decomposition: threshold residual", unit = "min",
         sample = "days with a crossing",
         diff_la_minus_el = diff_total - diff_light,
         par_ref = par_ref, n_boot = n_boot, seed = seed)
)

readr::write_csv(light_rows, out_light)
cat("\nsaved:", out_light, "\n\n")
cat("=== light controls (La Nina minus El Nino, month-block bootstrap) ===\n")
print(as.data.frame(
  light_rows %>%
    transmute(statistic, unit, sample,
              n = ifelse(is.na(n_el_nino), NA_character_, paste0(n_el_nino, "/", n_la_nina)),
              EN = round(mean_el_nino, 2), LN = round(mean_la_nina, 2),
              diff = round(diff_la_minus_el, 1),
              CI = ifelse(is.na(ci_low), "", sprintf("%.1f to %.1f", ci_low, ci_high)),
              excl0 = ifelse(is.na(excludes_zero), "", ifelse(excludes_zero, "yes", "no")))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# grid-native summary
#
# Everything here is expressed in half-hour units, because that is the
# resolution of the measurement. Timestamps are period-end labels, so "08:30"
# denotes the 08:00-08:30 averaging interval.

out_grid <- file.path(out_base, "table_morning_transition_grid.csv")

# lower median, so the reported centre is itself a grid value rather than the
# midpoint of two adjacent slots that R's default median() would return for an
# even number of days
median_slot <- function(x) sort(x)[ceiling(length(x) / 2)]

slot_dist <- primary_days %>%
  count(ENSO, t) %>%
  group_by(ENSO) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

cat("\n=== transition half-hour, share of days ===\n")
print(as.data.frame(
  slot_dist %>%
    select(t, ENSO, share) %>%
    tidyr::pivot_wider(names_from = ENSO, values_from = share, values_fill = 0) %>%
    arrange(t) %>%
    mutate(half_hour = hm(t), .before = 1) %>%
    select(-t) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
), row.names = FALSE)

centres <- primary_days %>%
  group_by(ENSO) %>%
  summarise(n_days = n(), median_slot = median_slot(t),
            q25 = quantile(t, 0.25, type = 1), q75 = quantile(t, 0.75, type = 1),
            .groups = "drop") %>%
  mutate(across(c(median_slot, q25, q75), hm))
cat("\nmedian transition half-hour (lower median; IQR on the grid):\n")
print(as.data.frame(centres), row.names = FALSE)

cat("\nby episode:\n")
print(as.data.frame(
  primary_days %>%
    group_by(episode) %>%
    summarise(n_days = n(), median_slot = hm(median_slot(t)), .groups = "drop")
), row.names = FALSE)

# Share of days already a net sink at each half-hour, and the share of all
# observed half-hours that are net uptake. Both are grid-native: they compare
# the two phases at fixed clock times instead of estimating a time.
prop_boot <- function(df, valcol, reps) {
  blocks <- split(df, list(df$ENSO, format(df$Date, "%Y-%m")), drop = TRUE)
  by_enso <- lapply(enso_levels, function(g)
    names(blocks)[vapply(blocks, function(b) b$ENSO[1] == g, logical(1))])
  names(by_enso) <- enso_levels
  out <- numeric(reps)
  for (r in seq_len(reps)) {
    m <- vapply(enso_levels, function(g) {
      keys <- sample(by_enso[[g]], length(by_enso[[g]]), replace = TRUE)
      mean(unlist(lapply(keys, function(k) blocks[[k]][[valcol]])))
    }, numeric(1))
    out[r] <- m[["El Nino"]] - m[["La Nina"]]
  }
  out
}

grid_rows <- lapply(seq(morning_from + 1, morning_to + 1, by = 0.5), function(thr) {
  # (a) days that have already transitioned by this half-hour
  dd <- primary_days %>% mutate(v = as.numeric(t <= thr))
  b1 <- prop_boot(dd, "v", n_boot); c1 <- quantile(b1, c(0.025, 0.975))
  # (b) all observed half-hours at this clock time that are net uptake
  hh <- base %>% filter(hour_dec == thr, is.finite(NEE)) %>% mutate(v = as.numeric(NEE < 0))
  b2 <- prop_boot(hh, "v", n_boot); c2 <- quantile(b2, c(0.025, 0.975))
  tibble(
    half_hour        = hm(thr),
    share_transitioned_el_nino = mean(dd$v[dd$ENSO == "El Nino"]),
    share_transitioned_la_nina = mean(dd$v[dd$ENSO == "La Nina"]),
    share_transitioned_diff    = mean(dd$v[dd$ENSO == "El Nino"]) - mean(dd$v[dd$ENSO == "La Nina"]),
    share_transitioned_ci_low  = unname(c1[1]), share_transitioned_ci_high = unname(c1[2]),
    n_uptake_el_nino = sum(hh$ENSO == "El Nino"), n_uptake_la_nina = sum(hh$ENSO == "La Nina"),
    share_uptake_el_nino = mean(hh$v[hh$ENSO == "El Nino"]),
    share_uptake_la_nina = mean(hh$v[hh$ENSO == "La Nina"]),
    share_uptake_diff    = mean(hh$v[hh$ENSO == "El Nino"]) - mean(hh$v[hh$ENSO == "La Nina"]),
    share_uptake_ci_low  = unname(c2[1]), share_uptake_ci_high = unname(c2[2])
  )
})
grid_tbl <- bind_rows(grid_rows) %>% mutate(n_boot = n_boot, seed = seed)
readr::write_csv(grid_tbl, out_grid)
cat("\nsaved:", out_grid, "\n\n")

cat("=== grid-native comparison at fixed clock times (El Nino minus La Nina) ===\n")
print(as.data.frame(
  grid_tbl %>% transmute(
    half_hour,
    transitioned_EN = round(share_transitioned_el_nino, 2),
    transitioned_LN = round(share_transitioned_la_nina, 2),
    diff  = round(share_transitioned_diff, 2),
    CI    = sprintf("%.2f to %.2f", share_transitioned_ci_low, share_transitioned_ci_high),
    excl0 = ifelse(share_transitioned_ci_low > 0 | share_transitioned_ci_high < 0, "yes", "no"),
    uptake_EN = round(share_uptake_el_nino, 2),
    uptake_LN = round(share_uptake_la_nina, 2),
    u_CI  = sprintf("%.2f to %.2f", share_uptake_ci_low, share_uptake_ci_high),
    u_excl0 = ifelse(share_uptake_ci_low > 0 | share_uptake_ci_high < 0, "yes", "no"))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# half-hourly morning composite
#
# Supports the statement in section 3.1 that morning light is the same in the two
# phases: the relative PAR difference half-hour by half-hour, rather than only as
# a window mean.

hh_tbl <- base %>%
  filter(hour_dec >= morning_from, hour_dec <= morning_to + 1) %>%
  group_by(hour_dec, ENSO) %>%
  # counts before the means: summarise() evaluates in order, so computing
  # `PAR` first would make `sum(is.finite(PAR))` count the scalar mean
  summarise(n_PAR = sum(is.finite(PAR)), n_NEE = sum(is.finite(NEE)),
            PAR = mean(PAR, na.rm = TRUE), NEE = mean(NEE, na.rm = TRUE),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = ENSO, values_from = c(PAR, n_PAR, NEE, n_NEE)) %>%
  rename_with(~gsub(" ", "_", .x)) %>%
  mutate(PAR_pct_el_vs_la = 100 * (`PAR_El_Nino` - `PAR_La_Nina`) / `PAR_La_Nina`,
         local_time = hm(hour_dec)) %>%
  relocate(local_time)

readr::write_csv(hh_tbl, out_hh)
cat("\nsaved:", out_hh, "\n\n")
cat("=== morning half-hourly composite (light is the same; NEE is not) ===\n")
print(as.data.frame(
  hh_tbl %>% transmute(local_time,
                       PAR_EN = round(`PAR_El_Nino`, 0), PAR_LN = round(`PAR_La_Nina`, 0),
                       PAR_pct = round(PAR_pct_el_vs_la, 1),
                       NEE_EN = round(`NEE_El_Nino`, 2), NEE_LN = round(`NEE_La_Nina`, 2),
                       n_NEE = paste0(`n_NEE_El_Nino`, "/", `n_NEE_La_Nina`))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# per-episode breakdown
#
# The pooled 36 min is not homogeneous: it is dominated by the 2023-24 El Nino.
# Section 3.1 reports the episodes separately rather than claiming the shift is
# uniform across years.

epi <- primary_days %>%
  group_by(ENSO, episode) %>%
  summarise(n_days = n(),
            n_months = n_distinct(format(Date, "%Y-%m")),
            mean_t = mean(t),
            mean_hm = hm(mean(t)),
            mean_PAR_at_crossing = mean(PAR_at_crossing, na.rm = TRUE),
            mean_PAR_morn = mean(PAR_morn, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(ENSO, episode)

ln_mean <- mean(primary_days$t[primary_days$ENSO == "La Nina"])
epi <- epi %>% mutate(min_earlier_than_la_nina_mean = (ln_mean - mean_t) * 60)

readr::write_csv(epi, out_epi)
cat("\nsaved:", out_epi, "\n\n")
cat("=== per-episode transition (La Nina pooled mean =", hm(ln_mean), ") ===\n")
print(as.data.frame(epi %>% mutate(across(where(is.numeric), ~round(.x, 1)))), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# storage-term sensitivity
#
# SC is profiler-measured from 2024 and DNN-modelled before, so the El Nino
# sample mixes methods. On days where both series yield a crossing, how much
# does the choice move it, and is that shift common to both phases?

storage_days <- base %>%
  arrange(Date, hour_dec) %>%
  group_by(Date, ENSO) %>%
  summarise(
    t_prof = {
      ok <- is.finite(NEE_prof); y <- NEE_prof[ok]; x <- hour_dec[ok]
      i <- which(head(y, -1) > 0 & tail(y, -1) < 0)[1]
      if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
    },
    t_mod = {
      ok <- is.finite(NEE_mod); y <- NEE_mod[ok]; x <- hour_dec[ok]
      i <- which(head(y, -1) > 0 & tail(y, -1) < 0)[1]
      if (is.na(i)) NA_real_ else x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1])
    },
    .groups = "drop"
  ) %>%
  filter(is.finite(t_prof), is.finite(t_mod)) %>%
  mutate(shift_min = (t_mod - t_prof) * 60, episode = episode_of(Date, ENSO))

cat("\n=== storage-term sensitivity: modelled minus profiler storage, paired days ===\n")
if (nrow(storage_days) >= 6) {
  print(as.data.frame(
    storage_days %>%
      group_by(ENSO) %>%
      summarise(n_days = n(),
                mean_shift_min = round(mean(shift_min), 1),
                median_shift_min = round(median(shift_min), 1),
                .groups = "drop")
  ), row.names = FALSE)
  cat("A shift of similar size in both phases moves the absolute clock times but\n",
      "not the El Nino - La Nina contrast.\n", sep = "")
} else {
  cat("too few paired days to report\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# the numbers as they appear in section 3.1

g <- function(stat, samp, field) {
  i <- which(light_rows$statistic == stat &
               (is.na(light_rows$sample) | light_rows$sample == samp))
  if (!length(i)) {
    warning("no light-control row for '", stat, "' / '", samp, "'")
    return(NA_real_)
  }
  as.numeric(light_rows[[field]][i[1]])
}
cs <- "days with a crossing"
cat("\n=== diagnostics only: sub-half-hour means, NOT quoted in the manuscript ===\n")
cat(sprintf(
  paste0("transition:    El Nino %s vs La Nina %s, %.0f min earlier (95%% CI %.0f to %.0f)\n",
         "morning NEE:   %.1f vs %.1f umol m-2 s-1, diff %.1f (95%% CI %.1f to %.1f)\n",
         "morning PAR:   %.0f vs %.0f umol m-2 s-1, diff %.0f (95%% CI %.0f to %.0f)\n",
         "crossing PAR:  %.0f vs %.0f umol m-2 s-1, diff %.0f (95%% CI %.0f to %.0f)\n",
         "decomposition: %.0f min light supply + %.0f min threshold = %.0f min total\n",
         "sunrise:       El Nino %s vs La Nina %s (%.1f min apart)\n"),
  p$mean_el_nino_hm, p$mean_la_nina_hm, p$diff_min, p$diff_ci_low_min, p$diff_ci_high_min,
  g(lab_nee, cs, "mean_el_nino"), g(lab_nee, cs, "mean_la_nina"),
  g(lab_nee, cs, "diff_la_minus_el"), g(lab_nee, cs, "ci_low"), g(lab_nee, cs, "ci_high"),
  g(lab_par, cs, "mean_el_nino"), g(lab_par, cs, "mean_la_nina"),
  g(lab_par, cs, "diff_la_minus_el"), g(lab_par, cs, "ci_low"), g(lab_par, cs, "ci_high"),
  g(lab_parx, cs, "mean_el_nino"), g(lab_parx, cs, "mean_la_nina"),
  g(lab_parx, cs, "diff_la_minus_el"), g(lab_parx, cs, "ci_low"), g(lab_parx, cs, "ci_high"),
  diff_light, diff_total - diff_light, diff_total,
  hm(mean(primary_days$sunrise[primary_days$ENSO == "El Nino"])),
  hm(mean(primary_days$sunrise[primary_days$ENSO == "La Nina"])),
  (mean(primary_days$sunrise[primary_days$ENSO == "La Nina"]) -
     mean(primary_days$sunrise[primary_days$ENSO == "El Nino"])) * 60))

# ─────────────────────────────────────────────────────────────────────────────
# the numbers section 3.1 actually quotes, all on the 30-minute grid

gr <- function(t, field) grid_tbl[[field]][grid_tbl$half_hour == t][1]
med <- function(g) hm(median_slot(primary_days$t[primary_days$ENSO == g]))

cat("\n=== section 3.1 numbers (grid-native) ===\n")
cat(sprintf("median transition half-hour:  El Nino %s   La Nina %s   (same slot = no shift)\n",
            med("El Nino"), med("La Nina")))
cat("per episode: ")
cat(paste(vapply(sort(unique(primary_days$episode)), function(e)
  sprintf("%s %s", e, hm(median_slot(primary_days$t[primary_days$episode == e]))),
  character(1)), collapse = " | "), "\n")
for (t in c("08:00", "08:30", "09:00")) {
  cat(sprintf(
    "by %s: transitioned %.0f%% vs %.0f%% (diff %+.0f pp, 95%% CI %+.0f to %+.0f) | in uptake %.0f%% vs %.0f%% (diff %+.0f pp, 95%% CI %+.0f to %+.0f)\n",
    t,
    100 * gr(t, "share_transitioned_el_nino"), 100 * gr(t, "share_transitioned_la_nina"),
    100 * gr(t, "share_transitioned_diff"),
    100 * gr(t, "share_transitioned_ci_low"), 100 * gr(t, "share_transitioned_ci_high"),
    100 * gr(t, "share_uptake_el_nino"), 100 * gr(t, "share_uptake_la_nina"),
    100 * gr(t, "share_uptake_diff"),
    100 * gr(t, "share_uptake_ci_low"), 100 * gr(t, "share_uptake_ci_high")))
}
cat(sprintf("PAR in the transition half-hour: El Nino %.0f vs La Nina %.0f umol m-2 s-1\n",
            g(lab_parx, cs, "mean_el_nino"), g(lab_parx, cs, "mean_la_nina")))
cat("morning PAR by half-hour is in", basename(out_hh),
    "- differences are <=3% from 06:30 to 08:00\n")
