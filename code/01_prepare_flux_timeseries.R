# ───────────────────────────────────────────────────────────────────────────────
# ts_NEE_and_NEE_ok_full_resolution_with_ENSO_masks_plus_Reco_GEP.R
# - loads dataset_from_matlab.csv
# - defines Season + ENSO masks
# - builds NEE_ok from NEE_f, invalidating long gaps in raw NEE
# - estimates Reco and GEP using "prev+next full nights" method:
#     * night block D is [D 18:00, D+1 06:00) with SW_IN/Rg < 10 W m-2
#     * day block D is [D 06:00, D 18:00)
#     * Reco_day(D) = mean(NEE night points from night blocks (D-1) and D)
#     * assigns Reco only to day block rows
#     * sign convention used here:
#         NEE > 0  = net CO2 release
#         NEE < 0  = net CO2 uptake
#         Reco > 0 = ecosystem respiration release
#         GEP  < 0 = photosynthetic uptake
#       therefore:
#         NEE = Reco + GEP
#         GEP = NEE - Reco
# - saves enriched csv
# - makes time-series plots (points) for NEE, NEE_ok, Reco, GEP with ENSO shading
# - GEP plot shows daytime rows only
# - limits data to <= 2024-10-15 and x-axis to [2018-01-01, 2024-10-15]
#
# adds in this version:
# - friagem mask: daily Tmax and Tmin below monthly 10th percentiles for >=3 days
# - friagem-highlighted time series: Ta, NEE_ok, Reco, GEP (black) vs non-friagem (gray)
# - monthly diel-cycle daily C balance plot (season colored + ENSO shading)
# - 16-day Ta time series (season colored + ENSO shading)
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(ggplot2)
  library(zoo)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

flux_fp <- file.path(input_path, "dataset_from_matlab.csv")
out_csv <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")

out_nee     <- file.path(fig_path, "ts_NEE_full_season_ENSO.png")
out_ok      <- file.path(fig_path, "ts_NEE_ok_full_season_ENSO.png")
out_reco    <- file.path(fig_path, "ts_Reco_full_season_ENSO.png")
out_gep     <- file.path(fig_path, "ts_GEP_full_season_ENSO_day_only.png")

# friagem-highlight figures
out_ta_friagem    <- file.path(fig_path, "ts_Ta_full_friagem_black_ENSO.png")
out_ok_friagem    <- file.path(fig_path, "ts_NEE_ok_full_friagem_black_ENSO.png")
out_reco_friagem  <- file.path(fig_path, "ts_Reco_full_friagem_black_ENSO.png")
out_gep_friagem   <- file.path(fig_path, "ts_GEP_full_friagem_black_ENSO_day_only.png")

# outputs: Ta time series
out_ta_full  <- file.path(fig_path, "ts_Ta_full_season_ENSO.png")
out_ta_16day <- file.path(fig_path, "ts_Ta_16day_season_ENSO.png")

# axis label
ylab_ta <- expression(italic(T)[plain(a)] ~ (degree*C))

monthly_fp  <- file.path(output_path, "tambopata_monthly_diurnal_sum_gC_m2_day.csv")
out_monthly <- file.path(fig_path, "ts_monthly_diurnal_sum_gC_m2_day_seasoncolor_ENSOshade.png")

# ───────────────────────────────────────────────────────────────────────────────
# config
tz_local <- "America/Lima"

dt_col    <- "tv_dt"
nee_col   <- "NEE"
nee_f_col <- "NEE_f"
par_col   <- "PPFD_IN_1_1_1"
vpd_col   <- "VPD_1_1_1"
ta_col    <- "TA_1_1_1"
swin_col  <- "SW_IN_1_1_1"
rn_col    <- "NETRAD_1_1_1"
fc_col    <- "FC"

# seasons (match your logic)
dry_months <- 5:10

# NEE_ok logic
gap_run_thresh <- 24  # >24 consecutive half-hours -> invalidate NEE_f there

# reco/gep logic (prev+next full nights)
swin_night_thresh <- 10
night_start_hour <- 18L
night_end_hour   <- 6L   # next day
day_start_hour   <- 6L
day_end_hour     <- 18L
min_night_points_total <- 8  # across the two nights combined

# plot config
# season_cols / enso_cols come from palette.R (sourced via paths.R):
#   season = blue (wet) / red (dry); ENSO = gold (El Nino) / purple (La Nina) / grey (Neutral)

tmax <- ymd_hms("2024-10-15 23:59:59", tz = tz_local)
xmin_plot <- ymd_hms("2018-01-01 00:00:00", tz = tz_local)
xmax_plot <- tmax

# axis labels (math expression; consistent across figures)
ylab_nee  <- expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))
ylab_ok   <- expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))
ylab_reco <- expression(Reco ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))
ylab_gep  <- expression(GEP ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))
ylab_monthly_c <- expression(paste("Daily C balance from monthly mean diel cycle (g C ", m^-2, " ", day^-1, ")"))

msg_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("n/a")
  paste0("[", signif(min(x), 4), ", ", signif(max(x), 4), "]")
}

# ───────────────────────────────────────────────────────────────────────────────
# load + parse
stopifnot(file.exists(flux_fp))
df_raw <- read_csv(flux_fp, show_col_types = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# apply consistent u* filtering to the measured NEE (per-year REddyProc threshold
# + pooled fallback for the years REddyProc failed). done here, at load time, so
# the filtered NEE propagates to the enriched CSV and everything built from it
# (the 02 -> enso_conditions RData -> 03 chain, and build_measured_composite ->
# tambopata_48points_per_month_measured.csv -> 05/09/11). see code/ustar_filter.R.
# NEE_ok/Reco/GEP are unaffected: they come from the already-u*-filtered NEE_f.
df_raw <- apply_ustar_filter(df_raw, nee_col = nee_col, ustar_col = "USTAR",
                             datetime_col = dt_col)

# ───────────────────────────────────────────────────────────────────────────────
# diagnostic: MATLAB-corrected PAR diel cycle immediately after loading dataset_from_matlab.csv
stopifnot(dt_col %in% names(df_raw))
stopifnot(par_col %in% names(df_raw))

dt_str <- df_raw[[dt_col]]

parse_dt <- function(x, is_utc, tz_local = "America/Lima") {
  dt <- lubridate::dmy_hms(x, tz = if (is_utc) "UTC" else tz_local, quiet = TRUE)
  if (all(is.na(dt))) dt <- as.POSIXct(x, tz = if (is_utc) "UTC" else tz_local)
  if (is_utc) dt <- lubridate::with_tz(dt, tz_local)
  dt
}

df_ppfd <- tibble(
  DateTime_local_assuming_utc   = parse_dt(dt_str, is_utc = TRUE,  tz_local = tz_local),
  DateTime_local_assuming_local = parse_dt(dt_str, is_utc = FALSE, tz_local = tz_local),
  PAR_corrected_SWin = as.numeric(df_raw[[par_col]])
)

make_diel_par <- function(dt_vec, par_vec, label) {
  tibble(DateTime = dt_vec, PAR_corrected_SWin = par_vec) %>%
    filter(!is.na(DateTime), is.finite(PAR_corrected_SWin)) %>%
    mutate(
      DateTime_30 = floor_date(DateTime, unit = "30 minutes"),
      hour = format(DateTime_30, "%H:%M:%S")
    ) %>%
    group_by(hour) %>%
    summarise(mean_par = mean(PAR_corrected_SWin, na.rm = TRUE), .groups = "drop") %>%
    mutate(mode = label)
}

diel_all <- bind_rows(
  make_diel_par(df_ppfd$DateTime_local_assuming_utc,   df_ppfd$PAR_corrected_SWin, "tv_dt parsed as UTC -> converted to America/Lima"),
  make_diel_par(df_ppfd$DateTime_local_assuming_local, df_ppfd$PAR_corrected_SWin, "tv_dt parsed as America/Lima (no conversion)")
)

halfhour_levels <- sprintf(
  "%02d:%02d:00",
  rep(0:23, each = 2),
  rep(c(0, 30), times = 24)
)
diel_all$hour <- factor(diel_all$hour, levels = halfhour_levels)

p_ppfd <- ggplot(diel_all, aes(x = hour, y = mean_par, color = mode, group = mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.4) +
  scale_x_discrete(breaks = c("02:00:00", "08:00:00", "14:00:00", "20:00:00")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "MATLAB-corrected PAR diel cycle check right after loading dataset_from_matlab.csv",
    y = expression(PAR ~ (mu*mol ~ m^-2 ~ s^-1)),
    color = NULL
  )

print(p_ppfd)

peak_hours <- diel_all %>%
  group_by(mode) %>%
  slice_max(mean_par, n = 3, with_ties = FALSE) %>%
  arrange(mode, desc(mean_par))
print(peak_hours)

# ───────────────────────────────────────────────────────────────────────────────
# now continue with your normal pipeline
df <- df_raw %>%
  mutate(
    DateTime = dmy_hms(.data[[dt_col]], tz = tz_local),
    NEE      = as.numeric(.data[[nee_col]]),
    VPD_kPa  = as.numeric(.data[[vpd_col]]) / 10,
    PAR_corrected_SWin = as.numeric(.data[[par_col]]),
    PAR = PAR_corrected_SWin,
    PAR_source = "matlab_second_stage_swin_corrected",
    SW_IN    = as.numeric(.data[[swin_col]]),
    SWin     = SW_IN,
    Rn       = as.numeric(.data[[rn_col]]),
    Fc       = as.numeric(.data[[fc_col]]),
    Ta       = as.numeric(.data[[ta_col]])
  ) %>%
  filter(!is.na(DateTime)) %>%
  arrange(DateTime)

cat(
  "PAR alias:",
  par_col, "from MATLAB export -> PAR_corrected_SWin -> PAR;",
  "source = matlab_second_stage_swin_corrected",
  "\n"
)

# season and year + day block helper
df <- df %>%
  mutate(
    year   = year(DateTime),
    month  = month(DateTime),
    season = if_else(month %in% dry_months, "dry", "wet"),
    Date   = as.Date(DateTime, tz = tz_local),
    hour   = hour(DateTime),
    is_day_block = hour >= day_start_hour & hour < day_end_hour
  )

# ───────────────────────────────────────────────────────────────────────────────
# ENSO masks (exactly your logic)
df <- df %>%
  mutate(ENSO = "neutral")

positive_enso_mask <-
  (df$DateTime > ymd_hms("2018-09-01 00:00:00", tz = tz_local) &
     df$DateTime < ymd_hms("2019-07-31 23:59:59", tz = tz_local)) |
  (df$DateTime > ymd_hms("2023-04-01 00:00:00", tz = tz_local) &
     df$DateTime < ymd_hms("2024-05-31 23:59:59", tz = tz_local))

negative_enso_mask <-
  (df$DateTime > ymd_hms("2016-07-01 00:00:00", tz = tz_local) &
     df$DateTime < ymd_hms("2017-01-31 23:59:59", tz = tz_local)) |
  (df$DateTime > ymd_hms("2017-09-01 00:00:00", tz = tz_local) &
     df$DateTime < ymd_hms("2018-05-31 23:59:59", tz = tz_local)) |
  (df$DateTime > ymd_hms("2020-07-01 00:00:00", tz = tz_local) &
     df$DateTime < ymd_hms("2023-02-28 23:59:59", tz = tz_local))

df$ENSO[positive_enso_mask] <- "El Nino"
df$ENSO[negative_enso_mask] <- "La Nina"

# ───────────────────────────────────────────────────────────────────────────────
# build NEE_ok from NEE_f, invalidating long raw-NEE gaps
stopifnot(nee_f_col %in% names(df_raw))
df$NEE_ok <- as.numeric(df_raw[[nee_f_col]])

nee_missing <- !is.finite(df$NEE)
r <- rle(nee_missing)
long_missing_runs <- r$values & (r$lengths > gap_run_thresh)
mask_long_gaps <- inverse.rle(list(lengths = r$lengths, values = long_missing_runs))

df$NEE_ok[mask_long_gaps] <- NA_real_

cat("\nNEE_ok created:\n")
cat("rows invalidated due to long raw NEE gaps:", sum(mask_long_gaps), "\n")
cat("NEE_ok non-NA:", sum(is.finite(df$NEE_ok)), "/", nrow(df), "\n\n")

# ───────────────────────────────────────────────────────────────────────────────
# estimate Reco and GEP (prev+next full nights)
# sign convention used here:
#   NEE > 0  = net CO2 release to the atmosphere
#   NEE < 0  = net CO2 uptake by the ecosystem
#   Reco > 0 = ecosystem respiration release
#   GEP < 0  = photosynthetic uptake
#
# therefore:
#   NEE = Reco + GEP
#   GEP = NEE - Reco
night_df <- df %>%
  mutate(
    in_night_hours = (hour >= night_start_hour) | (hour < night_end_hour),
    night_block_date = case_when(
      hour >= night_start_hour ~ Date,
      hour <  night_end_hour   ~ Date - days(1),
      TRUE ~ NA_Date_
    ),
    is_night_point = in_night_hours & is.finite(SW_IN) & SW_IN < swin_night_thresh & is.finite(NEE_ok)
  ) %>%
  filter(is_night_point) %>%
  select(DateTime, night_block_date, NEE_ok)

cat("night points with finite NEE (SW_IN/Rg < 10 W m-2, within 18-06):", nrow(night_df), "\n")
if (!nrow(night_df)) stop("no night points found. check SW_IN/Rg threshold and datetime parsing.")

night_stats <- night_df %>%
  group_by(night_block_date) %>%
  summarise(
    sum_nee = sum(NEE_ok, na.rm = TRUE),
    n_nee   = dplyr::n(),
    .groups = "drop"
  )

day_list <- sort(unique(df$Date))

day_reco <- tibble(Date = day_list) %>%
  mutate(prev_night = Date - days(1)) %>%
  left_join(
    night_stats %>% rename(prev_sum = sum_nee, prev_n = n_nee),
    by = c("prev_night" = "night_block_date")
  ) %>%
  left_join(
    night_stats %>% rename(next_sum = sum_nee, next_n = n_nee),
    by = c("Date" = "night_block_date")
  ) %>%
  mutate(
    n_total   = coalesce(prev_n, 0L) + coalesce(next_n, 0L),
    sum_total = coalesce(prev_sum, 0) + coalesce(next_sum, 0),
    Reco_day_raw = if_else(n_total >= min_night_points_total, sum_total / n_total, NA_real_),
    # a negative nightly mean NEE means the night was net uptake, which is not a
    # usable respiration estimate. drop those days rather than propagate them:
    # they are the main source of spuriously positive daytime GEP (all 16 fall in
    # 2018 and they carry a 37.7% positive-GEP rate against 8.4% elsewhere).
    Reco_day  = if_else(is.finite(Reco_day_raw) & Reco_day_raw < 0, NA_real_, Reco_day_raw)
  ) %>%
  select(Date, Reco_day, Reco_day_raw, n_total)

df <- df %>%
  left_join(day_reco, by = "Date") %>%
  mutate(
    Reco = if_else(is_day_block, Reco_day, NA_real_),
    GEP  = if_else(
      is_day_block & is.finite(NEE_ok) & is.finite(Reco),
      NEE_ok - Reco,
      NA_real_
    )
  ) %>%
  select(-Reco_day, -Reco_day_raw)

# positive GEP is retained on purpose. Reco carries random and systematic error,
# so GEP obtained by subtraction legitimately takes small positive values, and
# clipping them to zero would censor one side of that error and inflate the
# retrieved uptake (~3%, 154 g C m-2 over the record). retaining them also keeps
# NEE = Reco + GEP exact at every half-hour. the positive values are not white
# noise: 82% of the affected carbon falls between 06:00 and 09:00 local on days
# following calm nights, when respired CO2 that escaped the nighttime NEE (so
# Reco is biased low) vents after sunrise (so morning NEE is biased high).
# see also the REddyProc FAQ, which advises against removing negative GPP.

cat("\nReco/GEP created:\n")
cat("days total:", nrow(day_reco), "\n")
cat("days with Reco_day:", sum(is.finite(day_reco$Reco_day)), "\n")
cat("median n_total:", median(day_reco$n_total, na.rm = TRUE), "\n")
cat("min n_total:", min(day_reco$n_total, na.rm = TRUE), "\n")
cat("Reco NA fraction:", mean(is.na(df$Reco)), "\n")
cat("GEP  NA fraction:", mean(is.na(df$GEP)), "\n\n")

# ───────────────────────────────────────────────────────────────────────────────
# limit time range (data and plotting)
df <- df %>%
  filter(DateTime <= tmax)

# ───────────────────────────────────────────────────────────────────────────────
# friagem mask: daily Tmax and Tmin below monthly 10th percentiles for >=3 days

stopifnot("Ta" %in% names(df))

# Following a common cold-wave definition: daily maximum and minimum
# temperatures below their 10th percentile for at least three consecutive days.
min_ta_points_per_day <- 30L
min_friagem_run_days <- 3L

daily_ta <- df %>%
  mutate(
    day = as.Date(DateTime, tz = tz_local),
    month = month(day)
  ) %>%
  group_by(day) %>%
  summarise(
    month = dplyr::first(month),
    n_Ta = sum(is.finite(Ta)),
    Ta_min = if (n_Ta > 0L) min(Ta[is.finite(Ta)]) else NA_real_,
    Ta_max = if (n_Ta > 0L) max(Ta[is.finite(Ta)]) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    eligible = n_Ta >= min_ta_points_per_day
  )

friagem_thresholds <- daily_ta %>%
  filter(eligible, is.finite(Ta_min), is.finite(Ta_max)) %>%
  group_by(month) %>%
  summarise(
    Tmin_p10 = as.numeric(quantile(Ta_min, 0.10, na.rm = TRUE)),
    Tmax_p10 = as.numeric(quantile(Ta_max, 0.10, na.rm = TRUE)),
    n_days = dplyr::n(),
    .groups = "drop"
  )

friagem_candidates <- daily_ta %>%
  left_join(friagem_thresholds, by = "month") %>%
  mutate(
    friagem_candidate = eligible &
      is.finite(Ta_min) & is.finite(Ta_max) &
      is.finite(Tmin_p10) & is.finite(Tmax_p10) &
      Ta_min < Tmin_p10 & Ta_max < Tmax_p10
  ) %>%
  arrange(day)


# Collapse candidate days into consecutive-date runs.
friagem_runs <- friagem_candidates %>%
  filter(friagem_candidate) %>%
  mutate(
    gap_days = as.integer(day - dplyr::lag(day, default = dplyr::first(day))),
    run_id = cumsum(gap_days != 1L)
  ) %>%
  group_by(run_id) %>%
  summarise(
    start_day = min(day),
    end_day = max(day),
    n_days = dplyr::n(),
    .groups = "drop"
  ) %>%
  filter(n_days >= min_friagem_run_days) %>%
  mutate(friagem_event = dplyr::row_number())

friagem_event_days <- if (nrow(friagem_runs)) {
  data.frame(
    day = do.call(c, Map(
      seq.Date,
      friagem_runs$start_day,
      friagem_runs$end_day,
      MoreArgs = list(by = "day")
    ))
  )
} else {
  data.frame(day = as.Date(character()))
}

friagem_days <- friagem_candidates %>%
  select(day) %>%
  left_join(
    friagem_event_days %>% distinct(day) %>% mutate(friagem = TRUE),
    by = "day"
  ) %>%
  mutate(friagem = if_else(is.na(friagem), FALSE, friagem))

df <- df %>%
  mutate(day = as.Date(DateTime, tz = tz_local)) %>%
  left_join(friagem_days, by = "day") %>%
  mutate(friagem = if_else(is.na(friagem), FALSE, friagem)) %>%
  select(-day)

cat("\nfriagem monthly 10th-percentile thresholds:\n")
print(friagem_thresholds %>% arrange(month))
cat("friagem candidate days:", sum(friagem_candidates$friagem_candidate, na.rm = TRUE), "\n")
cat("friagem retained events:", nrow(friagem_runs), "\n")
cat("friagem retained days:", sum(friagem_days$friagem, na.rm = TRUE), "\n")
if (nrow(friagem_runs)) {
  cat("friagem retained event spans:\n")
  print(friagem_runs %>% select(friagem_event, start_day, end_day, n_days))
}
cat("\n")

# ───────────────────────────────────────────────────────────────────────────────
# ENSO shading rectangles (build after time limit)
enso_rect <- df %>%
  transmute(DateTime, ENSO = as.character(ENSO)) %>%
  arrange(DateTime) %>%
  mutate(
    change = ENSO != dplyr::lag(ENSO, default = dplyr::first(ENSO)),
    grp = cumsum(change)
  ) %>%
  group_by(grp) %>%
  summarise(
    ENSO = dplyr::first(ENSO),
    xmin = min(DateTime),
    xmax = max(DateTime),
    .groups = "drop"
  )

# ───────────────────────────────────────────────────────────────────────────────
# plotting helpers (defined BEFORE any calls)

make_ts_plot <- function(d, ycol, ylab) {
  ggplot(d, aes(x = DateTime)) +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "black",
      alpha = 0.6
    ) +
    scale_fill_manual(
      values = enso_fill,
      drop = FALSE,
      name = NULL,
      labels = c(
        "El Nino" = "El Niño",
        "La Nina" = "La Niña",
        "neutral" = "neutral"
      )
    ) +
    geom_point(
      aes(y = .data[[ycol]], color = season),
      size = 0.4, alpha = 0.8, na.rm = TRUE
    ) +
    scale_color_manual(
      values = season_cols,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(override.aes = list(size = 3.5, alpha = 1))
    ) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))
}

make_friagem_ts_plot <- function(d, ycol, ylab) {
  stopifnot("friagem" %in% names(d))
  
  d2 <- d %>%
    filter(is.finite(.data[[ycol]])) %>%
    mutate(friagem_flag = if_else(friagem, "friagem", "non-friagem"))
  
  friagem_cols <- c("non-friagem" = "grey70", "friagem" = "black")
  
  ggplot(d2, aes(x = DateTime)) +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "black",
      alpha = 0.6
    ) +
    scale_fill_manual(
      values = enso_fill,
      drop = FALSE,
      name = NULL,
      labels = c(
        "El Nino" = "El Niño",
        "La Nina" = "La Niña",
        "neutral" = "neutral"
      )
    ) +
    geom_point(
      aes(y = .data[[ycol]], color = friagem_flag),
      size = 0.4, alpha = 0.85, na.rm = TRUE
    ) +
    scale_color_manual(
      values = friagem_cols,
      drop = FALSE,
      name = NULL,
      labels = c("friagem" = "friagem", "non-friagem" = "non-friagem")
    ) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))
}

# ───────────────────────────────────────────────────────────────────────────────
# core figures (ENSO shading + season colors)
p_nee <- make_ts_plot(df, "NEE", ylab_nee)
ggsave(out_nee, p_nee, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_nee, "\n")

p_ok <- make_ts_plot(df, "NEE_ok", ylab_ok)
ggsave(out_ok, p_ok, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_ok, "\n")

p_reco <- make_ts_plot(df, "Reco", ylab_reco)
ggsave(out_reco, p_reco, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_reco, "\n")

p_gep <- make_ts_plot(df %>% filter(is_day_block), "GEP", ylab_gep)
ggsave(out_gep, p_gep, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_gep, "\n")

p_ta_full <- make_ts_plot(df, "Ta", ylab_ta)
ggsave(out_ta_full, p_ta_full, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_ta_full, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# friagem-highlighted figures (black = friagem, gray = non-friagem)
p_ta_friagem <- make_friagem_ts_plot(df, "Ta", ylab_ta)
ggsave(out_ta_friagem, p_ta_friagem, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_ta_friagem, "\n")

p_ok_friagem <- make_friagem_ts_plot(df, "NEE_ok", ylab_ok)
ggsave(out_ok_friagem, p_ok_friagem, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_ok_friagem, "\n")

# reco is daytime-only by construction, but plotting function will just skip NAs
p_reco_friagem <- make_friagem_ts_plot(df, "Reco", ylab_reco)
ggsave(out_reco_friagem, p_reco_friagem, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_reco_friagem, "\n")

# gep is daytime-only by construction; keep explicit day filter for clean plotting
p_gep_friagem <- make_friagem_ts_plot(df %>% filter(is_day_block), "GEP", ylab_gep)
ggsave(out_gep_friagem, p_gep_friagem, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_gep_friagem, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# new figure: monthly diel-cycle daily C balance (g C m^-2 day^-1), season-colored + ENSO shading
if (!file.exists(monthly_fp)) {
  sec_per_halfhour <- 1800
  gC_per_umolCO2 <- 12e-6

  monthly_ts_out <- df %>%
    filter(is.finite(NEE_ok), !is.na(DateTime)) %>%
    mutate(
      ym = format(DateTime, "%Y-%m"),
      hh = format(floor_date(DateTime, "30 minutes"), "%H:%M:%S")
    ) %>%
    group_by(ym, season, ENSO, hh) %>%
    summarise(NEE_ok_mean = mean(NEE_ok, na.rm = TRUE), .groups = "drop") %>%
    group_by(ym, season, ENSO) %>%
    summarise(
      gC_m2_day = sum(NEE_ok_mean * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
      n_halfhours = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(ym, season, ENSO)

  write_csv(monthly_ts_out, monthly_fp)
  cat("saved:", monthly_fp, "\n")
}

monthly_ts <- read_csv(monthly_fp, show_col_types = FALSE) %>%
  mutate(
    season = factor(tolower(trimws(season)), levels = c("dry", "wet")),
    date_mid = as.Date(paste0(ym, "-15")),
    ENSO_plot = case_when(
      ENSO %in% c("El Nino", "El Niño") ~ "El Nino",
      ENSO %in% c("La Nina", "La Niña") ~ "La Nina",
      TRUE                              ~ "neutral"
    )
  ) %>%
  arrange(date_mid, season, ENSO_plot)

enso_rect_month <- df %>%
  transmute(DateTime, ENSO = as.character(ENSO)) %>%
  arrange(DateTime) %>%
  mutate(
    change = ENSO != dplyr::lag(ENSO, default = dplyr::first(ENSO)),
    grp = cumsum(change)
  ) %>%
  group_by(grp) %>%
  summarise(
    ENSO = dplyr::first(ENSO),
    xmin = as.Date(min(DateTime)),
    xmax = as.Date(max(DateTime)) + 1,
    .groups = "drop"
  )

p_monthly <- ggplot(monthly_ts, aes(x = date_mid, y = gC_m2_day)) +
  geom_rect(
    data = enso_rect_month,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_point(aes(color = season, shape = season), size = 2.6, alpha = 0.95) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 3.5, alpha = 1))
  ) +
  scale_shape_manual(values = season_shapes, drop = FALSE, name = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(y = ylab_monthly_c) +
  coord_cartesian(xlim = c(as.Date("2018-01-01"), as.Date("2024-10-15")))

ggsave(out_monthly, p_monthly, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_monthly, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# 16-day block Ta time series (mean Ta per 16-day block)
origin_day <- as.Date("2018-01-01")  # matches your plotting range start

ta_16day <- df %>%
  filter(is.finite(Ta)) %>%
  mutate(
    day = as.Date(DateTime),
    block_id = as.integer(floor(as.numeric(day - origin_day) / 16)),
    block_start = origin_day + days(block_id * 16),
    block_mid = block_start + days(8)
  ) %>%
  group_by(block_start, block_mid) %>%
  summarise(
    Ta_mean = mean(Ta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    month = month(block_mid),
    season = if_else(month %in% dry_months, "dry", "wet")
  ) %>%
  filter(block_mid >= as.Date(xmin_plot), block_mid <= as.Date(xmax_plot)) %>%
  arrange(block_mid)

enso_rect_day <- df %>%
  transmute(day = as.Date(DateTime), ENSO = as.character(ENSO)) %>%
  arrange(day) %>%
  mutate(
    change = ENSO != dplyr::lag(ENSO, default = dplyr::first(ENSO)),
    grp = cumsum(change)
  ) %>%
  group_by(grp) %>%
  summarise(
    ENSO = dplyr::first(ENSO),
    xmin = min(day),
    xmax = max(day) + 1,
    .groups = "drop"
  )

p_ta_16day <- ggplot(ta_16day, aes(x = block_mid, y = Ta_mean)) +
  geom_rect(
    data = enso_rect_day,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_point(aes(color = season, shape = season), size = 1.6, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 3.5, alpha = 1))
  ) +
  scale_shape_manual(values = season_shapes, drop = FALSE, name = NULL) +
  labs(x = NULL, y = ylab_ta) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(as.Date("2018-01-01"), as.Date("2024-10-15")))

ggsave(out_ta_16day, p_ta_16day, width = 14, height = 3.8, dpi = 300)
cat("saved:", out_ta_16day, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# save enriched csv (after friagem column exists)
write_csv(df, out_csv)
cat("saved:", out_csv, "\n")



cat("GEP diagnostics:\n")
gep_day <- df$GEP[df$is_day_block & is.finite(df$GEP)]
cat("finite daytime GEP:", length(gep_day), "\n")
cat("GEP range:", msg_range(df$GEP), "\n")
cat("positive (unphysical-sign) GEP retained:", sum(gep_day > 0),
    sprintf("(%.2f%% of finite daytime GEP)\n", 100 * mean(gep_day > 0)))
cat("carbon in those positive values:",
    sprintf("%.1f g C m-2\n", sum(gep_day[gep_day > 0]) * 1800 * 12e-6))
cat("days with Reco dropped for being negative:",
    sum(is.finite(day_reco$Reco_day_raw) & day_reco$Reco_day_raw < 0), "of", nrow(day_reco), "\n")

# the identity NEE = Reco + GEP must hold exactly now that GEP is not clipped
ident <- df$is_day_block & is.finite(df$NEE_ok) & is.finite(df$Reco) & is.finite(df$GEP)
cat("max |NEE_ok - (Reco + GEP)| over", sum(ident), "daytime half-hours:",
    sprintf("%.3g\n\n", max(abs(df$NEE_ok[ident] - (df$Reco[ident] + df$GEP[ident])))))
