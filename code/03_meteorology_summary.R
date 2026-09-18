# ───────────────────────────────────────────────────────────────────────────────
# env + light-response + precipitation + ENSO figures
# + NEW: 3pm VPD 98.4% quantile (overall + optional 2012–2017 window)
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────
# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(minpack.lm)
  library(rlang)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(data.table)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
graphs_path <- file.path(graphs_path, "meteorology")

# ───────────────────────────────────────────────────────────────────────────────
# config
MIN_LIGHT_RESPONSE_POINTS <- 100

# nee_source for the per-year light-response FITS: "measured" = non-gap-filled raw NEE
#   (primary; matches the non-gap-filled contract), "gapfilled" = NEE_ok (sensitivity).
#   note: only the curve fits use this; yearly environmental means and the meteorology ANOVA
#   keep the full daytime record (filter_for_fit on NEE_ok) so site climate is characterised on all data.
nee_source <- "measured"
nee_col    <- if (identical(nee_source, "measured")) "NEE" else "NEE_ok"

pretty_names <- c(
  wet_el_nino   = "El Niño – wet",
  dry_el_nino   = "El Niño – dry",
  wet_la_nina   = "La Niña – wet",
  dry_la_nina   = "La Niña – dry",
  wet_neutral   = "Neutral – wet",
  dry_neutral   = "Neutral – dry"
)

# Two figures here place the individual-year points with position_jitter() and
# nothing seeded the RNG, so Fig1_env_ENSO_faceted_SWC_Ts.png and
# Fig1_env_ENSO_faceted_MEAN_SD_with_years.png (= supplement Fig. S3) came out
# byte-different on every run and drifted away from whatever was pasted into the
# manuscript. Only the jitter consumes randomness in this script, so seeding here
# makes both reproducible without touching any number.
set.seed(20260917)

# vpd-quantile settings
tz_local <- "America/Lima"
vpd_quantile_prob <- 0.984
vpd_use_exact_1500 <- FALSE  # if TRUE: requires timestamps exactly at 15:00:00
vpd_restrict_2012_2017 <- FALSE  # if TRUE: restrict quantile calculation to 2012-01-01..2017-12-31 (local time)

# ───────────────────────────────────────────────────────────────────────────────
# helpers

# add Year if absent
extract_year <- function(df) {
  if (!"Year" %in% names(df) && "DateTime" %in% names(df)) df$Year <- lubridate::year(df$DateTime)
  df
}

# normalize common legacy names
#
# This is deliberately a no-op. It used to hold
#   if ("SWC" %in% nm && !"SWC" %in% nm) d$SWC <- d$SWC
# whose condition is `X && !X`, i.e. never true -- so it never renamed anything.
# Nothing depends on it: SWC_1_1_1 -> SWC is already handled by the rename below
# (see the rename_with/all_of block), and reinstating a second rename here would
# risk producing a duplicate column. Kept as a named passthrough so the six call
# sites stay readable.
normalize_cols <- function(d) d

# safe mean helpers
mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# data filtering for light-curve fits
filter_for_fit <- function(df) {
  if (!all(c("PAR", "NEE_ok") %in% names(df))) return(df[0, , drop = FALSE])
  out <- df[is.finite(df$PAR) & is.finite(df$NEE_ok) & df$PAR > 20, , drop = FALSE]
  if (!nrow(out)) return(out)
  out$FC_pos <- -out$NEE_ok
  out
}

# light-response (non-rectangular hyperbola)
fit_light_curve <- function(df) {
  # fit on the selected NEE response (nee_col); env means / counts keep NEE_ok separately
  if (!all(c("PAR", nee_col) %in% names(df))) return(NULL)
  keep <- is.finite(df$PAR) & is.finite(df[[nee_col]]) & df$PAR > 20
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  df$FC_pos <- -df[[nee_col]]
  if (nrow(df) < MIN_LIGHT_RESPONSE_POINTS) return(NULL)
  if (length(unique(df$PAR)) < 3 || length(unique(df$FC_pos)) < 3) return(NULL)
  tryCatch({
    fit <- nlsLM(
      FC_pos ~ ((phi0 * Pmax * PAR) / (phi0 * PAR + Pmax)) - Rd,
      data = df,
      start = list(phi0 = 0.03, Pmax = 30, Rd = 5),
      lower = c(phi0 = 0.005, Pmax = 5, Rd = 0),
      upper = c(phi0 = 0.30,  Pmax = 80, Rd = 20),
      control = nls.lm.control(maxiter = 500)
    )
    list(fit = fit, N = nrow(df))
  }, error = function(e) NULL)
}

# plotting helper (year-colored scatter)
plot_var_by_ENSO_faceted <- function(df, var, ylab, season_means_df, year_colors) {
  v <- rlang::ensym(var)
  ggplot(df, aes(x = ENSO, y = !!v, color = factor(Year))) +
    geom_point(size = 2.6, alpha = 0.9,
               position = position_jitter(width = 0.15, height = 0)) +
    geom_hline(data = season_means_df,
               aes(yintercept = mean_val),
               inherit.aes = FALSE, linetype = "dashed") +
    facet_wrap(~ Season, nrow = 1) +
    scale_color_manual(values = year_colors, name = "Year") +
    labs(x = NULL, y = ylab, color = "Year") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank()
    )
}

# significance helpers
p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "ns",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}
skip_test_row <- function(variable, variable_label, test, comparison, note, Season = NA_character_) {
  out <- tibble(
    variable = variable,
    variable_label = variable_label,
    test = test,
    comparison = comparison,
    p_value = NA_real_,
    note = note
  )
  if (!is.na(Season)) out <- mutate(out, Season = Season, .before = test)
  out
}
sig_tests <- function(df, var) {
  v <- rlang::ensym(var)
  dat <- df %>% dplyr::select(Season, !!v) %>% dplyr::filter(!is.na(!!v))
  if (!nrow(dat) || length(unique(dat$Season)) < 2) return(list(stars = "ns", ymax = NA_real_))
  tt <- try(t.test(dat[[2]] ~ dat[[1]]), silent = TRUE)
  pval <- if (inherits(tt, "try-error")) NA_real_ else tt$p.value
  ymax <- suppressWarnings(max(dat[[2]], na.rm = TRUE))
  if (!is.finite(ymax)) ymax <- NA_real_
  list(stars = p_to_stars(pval), ymax = ymax)
}
add_sig <- function(p, sig_info) {
  if (is.null(sig_info) || is.na(sig_info$ymax) || !is.finite(sig_info$ymax)) return(p)
  p +
    expand_limits(y = sig_info$ymax * 1.10) +
    annotate("text", x = 1.5, y = sig_info$ymax * 1.06,
             label = sig_info$stars, size = 5, fontface = "bold")
}

# normalize ENSO/Season labels (handles case and accents)
normalize_enso <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x0 %in% c("el niño", "el nino") ~ "El Niño",
    x0 %in% c("la niña", "la nina") ~ "La Niña",
    x0 %in% c("neutral")            ~ "Neutral",
    TRUE                            ~ "Neutral"
  )
}
normalize_season <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(startsWith(x0, "w"), "Wet", "Dry")
}

# robust local-time conversion for a DateTime column
ensure_datetime_local <- function(d, tz_local) {
  if (!"DateTime" %in% names(d)) return(d)
  if (!inherits(d$DateTime, "POSIXct")) {
    # try flexible parsing (works for most ymd_hms strings)
    d$DateTime <- suppressWarnings(lubridate::ymd_hms(d$DateTime, quiet = TRUE, tz = tz_local))
    if (all(is.na(d$DateTime))) {
      # last resort
      d$DateTime <- as.POSIXct(d$DateTime, tz = tz_local)
    }
  } else {
    d$DateTime <- lubridate::with_tz(d$DateTime, tz_local)
  }
  d
}

# ───────────────────────────────────────────────────────────────────────────────
# data
load(file.path(output_path, "enso_conditions_data.RData"))

df_list <- list(
  wet_el_nino   = extract_year(normalize_cols(df_wet_el_nino)),
  dry_el_nino   = extract_year(normalize_cols(df_dry_el_nino)),
  wet_la_nina   = extract_year(normalize_cols(df_wet_la_nina)),
  dry_la_nina   = extract_year(normalize_cols(df_dry_la_nina)),
  wet_neutral   = extract_year(normalize_cols(df_wet_neutral)),
  dry_neutral   = extract_year(normalize_cols(df_dry_neutral))
)

# ───────────────────────────────────────────────────────────────────────────────
# cutoff: limit all datasets to <= Oct 15, 2024 (end of day) in UTC
# note: this is used for the yearly tables/figures as you had it.
cutoff_datetime <- as.POSIXct("2024-10-15 23:59:59", tz = "UTC")

df_list <- lapply(df_list, function(d) {
  
  # ensure DateTime is POSIXct in UTC (robust to character inputs)
  if ("DateTime" %in% names(d)) {
    if (!inherits(d$DateTime, "POSIXct")) {
      d$DateTime <- lubridate::ymd_hms(d$DateTime, quiet = TRUE, tz = "UTC")
      if (all(is.na(d$DateTime))) {
        d$DateTime <- as.POSIXct(d$DateTime, tz = "UTC")
      }
    } else {
      attr(d$DateTime, "tzone") <- "UTC"
    }
    
    # apply cutoff
    d <- d %>% dplyr::filter(!is.na(DateTime), DateTime <= cutoff_datetime)
  }
  
  # recompute Year from the filtered DateTime
  d <- extract_year(d)
  
  d
})

# ───────────────────────────────────────────────────────────────────────────────
# rename columns inside each df in df_list
rename_map <- c(
  "TA_1_1_1"       = "Ta",
  "WS_1_1_1"       = "Wind_speed",
  "USTAR"          = "u_star",
  "VPD_kPa"        = "VPD",
  "SWC_1_1_1"      = "SWC",
  "TS_3"           = "Ts"
)

df_list <- lapply(df_list, function(d) {
  if (!"PAR_corrected_SWin" %in% names(d) && "PPFD_IN_1_1_1" %in% names(d)) {
    d$PAR_corrected_SWin <- suppressWarnings(as.numeric(d$PPFD_IN_1_1_1))
  }
  if (!"PAR" %in% names(d) && "PAR_corrected_SWin" %in% names(d)) {
    d$PAR <- suppressWarnings(as.numeric(d$PAR_corrected_SWin))
  }
  if (!"PAR" %in% names(d)) {
    stop("MATLAB-corrected PAR column is missing from enso_conditions_data.RData")
  }
  par_mismatch <- is.finite(d$PAR) & is.finite(d$PAR_corrected_SWin) & abs(d$PAR - d$PAR_corrected_SWin) > 1e-8
  if (any(par_mismatch, na.rm = TRUE)) {
    stop("canonical PAR does not match PAR_corrected_SWin in enso_conditions_data.RData")
  }
  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(d) && !new %in% names(d)) {
      d <- dplyr::rename(d, !!new := !!rlang::sym(old))
    }
  }
  d
})

# ───────────────────────────────────────────────────────────────────────────────
# NEW: compute 3pm VPD 98.4% quantile on the full dataset (local time)

df_all_for_vpd <- bind_rows(df_list) %>%
  ensure_datetime_local(tz_local = tz_local) %>%
  dplyr::filter(!is.na(DateTime), !is.na(VPD))

if (vpd_restrict_2012_2017) {
  df_all_for_vpd <- df_all_for_vpd %>%
    dplyr::filter(
      DateTime >= as.POSIXct("2012-01-01 00:00:00", tz = tz_local),
      DateTime <= as.POSIXct("2017-12-31 23:59:59", tz = tz_local)
    )
}

df_3pm <- df_all_for_vpd %>%
  dplyr::filter(lubridate::hour(DateTime) == 15)

if (vpd_use_exact_1500) {
  df_3pm <- df_3pm %>%
    dplyr::filter(lubridate::minute(DateTime) == 0, lubridate::second(DateTime) == 0)
}

vpd_q984 <- as.numeric(quantile(df_3pm$VPD, probs = vpd_quantile_prob, na.rm = TRUE))

vpd_quantile_out_csv <- file.path(output_path, "vpd_3pm_quantile_0p984.csv")
write.csv(
  data.frame(
    tz_local = tz_local,
    prob = vpd_quantile_prob,
    restrict_2012_2017 = vpd_restrict_2012_2017,
    exact_1500 = vpd_use_exact_1500,
    n_points = nrow(df_3pm),
    vpd_quantile_kPa = vpd_q984
  ),
  vpd_quantile_out_csv,
  row.names = FALSE
)

p_vpd_hist <- ggplot(df_3pm, aes(x = VPD)) +
  geom_histogram(bins = 60, color = "white", fill = "steelblue", alpha = 0.75) +
  geom_vline(xintercept = vpd_q984, color = "red", linewidth = 1.1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "3 pm VPD distribution",
    subtitle = paste0("Red line = ", 100 * vpd_quantile_prob, "th percentile; n = ", nrow(df_3pm),
                      if (vpd_restrict_2012_2017) " (restricted to 2012–2017)" else ""),
    x = "VPD (kPa)",
    y = "Count"
  )

ggsave(file.path(paper_figures_path, "Fig_vpd_3pm_quantile_hist.png"),
       p_vpd_hist, width = 7.2, height = 4.6, dpi = 300)

cat(
  "\ncomputed 3pm VPD quantile:\n",
  "  - tz_local:", tz_local, "\n",
  "  - prob:", vpd_quantile_prob, "\n",
  "  - restrict_2012_2017:", vpd_restrict_2012_2017, "\n",
  "  - exact_1500:", vpd_use_exact_1500, "\n",
  "  - n_points:", nrow(df_3pm), "\n",
  "  - q:", round(vpd_q984, 4), "kPa\n",
  "wrote:\n  -", vpd_quantile_out_csv, "\n",
  "saved:\n  -", file.path(paper_figures_path, "Fig_vpd_3pm_quantile_hist.(png|pdf)"), "\n"
)

# ───────────────────────────────────────────────────────────────────────────────
# 1) yearly environmental means from the light-response fit rows
env_means_by_year <- bind_rows(lapply(names(df_list), function(gn) {
  df <- df_list[[gn]]
  if (!"Year" %in% names(df) && "DateTime" %in% names(df)) df$Year <- year(df$DateTime)
  df <- filter_for_fit(df)
  if (!nrow(df)) return(NULL)
  df %>%
    group_by(Year) %>%
    summarise(
      Ta      = mean_or_na(Ta),
      VPD     = mean_or_na(VPD),
      SWC     = mean_or_na(SWC),
      PAR     = mean_or_na(PAR),
      u_star  = mean_or_na(u_star),
      Ts      = mean_or_na(Ts),
      .groups = "drop"
    ) %>%
    mutate(GroupKey = gn, Group = pretty_names[[gn]])
}))

# 2) yearly light-curve parameters
params_all_years <- bind_rows(lapply(names(df_list), function(gn) {
  d0 <- df_list[[gn]]; if (!nrow(d0)) return(NULL)
  yrs <- sort(unique(d0$Year))
  bind_rows(lapply(yrs, function(y) {
    dy <- subset(d0, Year == y)
    res <- fit_light_curve(dy)
    if (is.null(res)) return(NULL)
    co <- coef(res$fit)
    phi0 <- unname(co["phi0"]); Pmax <- unname(co["Pmax"]); Rd <- unname(co["Rd"])
    LCP <- if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0)
      Rd * Pmax / (phi0 * (Pmax - Rd)) else NA_real_
    P2000 <- if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd))
      ((phi0 * Pmax * 2000) / (phi0 * 2000 + Pmax)) - Rd else NA_real_
    tibble(Group = pretty_names[[gn]], Year = y,
           phi0 = round(phi0, 4), Pmax = round(Pmax, 2),
           P2000 = round(P2000, 2), Rd = round(Rd, 2),
           LCP = round(LCP, 2), N = res$N)
  }))
}))

# 3) combined table
summary_table_yearly <- env_means_by_year %>%
  select(-GroupKey) %>%
  left_join(params_all_years, by = c("Group","Year")) %>%
  mutate(
    Season = ifelse(grepl("wet", tolower(Group)), "Wet", "Dry"),
    ENSO = dplyr::case_when(
      grepl("El Niño", Group) ~ "El Niño",
      grepl("La Niña", Group) ~ "La Niña",
      grepl("Neutral", Group) ~ "Neutral",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = c("Wet","Dry")),
    ENSO   = factor(ENSO,   levels = c("El Niño","La Niña","Neutral"))
  )

# 4) filter combos with too few raw points
counts_by_year <- bind_rows(lapply(names(df_list), function(gn) {
  d0 <- df_list[[gn]]
  if (!"Year" %in% names(d0) && "DateTime" %in% names(d0)) d0$Year <- year(d0$DateTime)
  d0 <- filter_for_fit(d0); if (!nrow(d0)) return(NULL)
  summarise(group_by(d0, Year), N_points = n(), .groups = "drop") %>%
    mutate(Group = pretty_names[[gn]])
}))
cat("\nharmonized light-response rows by year and ENSO-season:\n")
print(counts_by_year %>% arrange(Group, Year))

valid_keys <- counts_by_year %>%
  filter(N_points >= MIN_LIGHT_RESPONSE_POINTS) %>%
  select(Group, Year)

summary_table_yearly <- summary_table_yearly %>%
  inner_join(valid_keys, by = c("Group","Year")) %>%
  mutate(
    Season = factor(Season, levels = c("Wet","Dry")),
    ENSO   = factor(ENSO,   levels = c("El Niño","La Niña","Neutral"))
  )

# ───────────────────────────────────────────────────────────────────────────────
# precipitation (IMERG 30-min) → mean mm/day by Year × ENSO × Season

summary_table_yearly_filt <- summary_table_yearly

# guard against the historical fan-out duplication (C1): exactly one row per Group x Year
stopifnot(!any(duplicated(summary_table_yearly[, c("Group", "Year")])))
stopifnot(!any(duplicated(summary_table_yearly[, c("Season", "ENSO", "Year")])))

out_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year.csv")
out_csv_with_ts <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
# NOTE: these by-year tables are written AFTER the precipitation join below (see
# "write the by-year parameter tables" further down), so that precip_mm_day is
# available as a candidate predictor for the light-response LMM (code/06).

# find and read the most recent IMERG 30-min CSV from input_folder
imerg_csv <- list.files(input_folder, pattern = "(?i)imerg.*30min.*\\.csv$", full.names = TRUE)
if (!length(imerg_csv)) stop("no IMERG 30-min CSV found in input_folder")
imerg_csv <- imerg_csv[order(file.info(imerg_csv)$mtime, decreasing = TRUE)][1]
imerg <- data.table::fread(imerg_csv)

# parse timestamp (prefer explicit local column; otherwise convert from UTC)
utc_col   <- intersect(names(imerg), c("time_utc_iso","time_utc","time","system:time_start"))[1]
local_col <- intersect(names(imerg), c("time_local_minus05","time_local"))[1]
if (!is.na(local_col)) {
  DateTime <- lubridate::parse_date_time(
    imerg[[local_col]], orders = c("Y-m-d H:M:S","Y-m-d H:M","Ymd HMS","Ymd HM"), tz = tz_local
  )
} else if (!is.na(utc_col)) {
  t_utc <- lubridate::parse_date_time(
    imerg[[utc_col]], orders = c("Y-m-d H:M:S","Y-m-d H:M","Ymd HMS","Ymd HM","YmdHMS","YmdHM"), tz = "UTC"
  )
  DateTime <- lubridate::with_tz(t_utc, tz_local)
} else stop("could not find a time column (time_utc*/time_local*).")

# precipitation per 30 min (mm)
if ("precip_mm_30min" %in% names(imerg)) {
  pr30 <- as.numeric(imerg$precip_mm_30min)
} else {
  rate_col <- intersect(names(imerg), c("precip_mm_per_hr","precipitation","precipitationCal"))[1]
  if (is.na(rate_col)) stop("no precipitation column found (precip_mm_30min or rate in mm/hr)")
  pr30 <- as.numeric(imerg[[rate_col]]) * 0.5
}

imerg30 <- data.frame(DateTime = DateTime, precip_mm_30min = pr30) %>%
  dplyr::filter(DateTime >= as.POSIXct("2017-01-01 00:00:00", tz = tz_local),
                DateTime <= as.POSIXct("2024-12-31 23:59:59", tz = tz_local))

# daily totals in local time
pr_daily <- imerg30 %>%
  mutate(Date = as.Date(DateTime)) %>%
  group_by(Date) %>%
  summarise(precip_day_mm = sum(precip_mm_30min, na.rm = TRUE), .groups = "drop")

# attach ENSO/Season from flux csv if available, else fallback with rules
flux_fp <- file.path(input_folder, "dataset_from_matlab.csv")
if (file.exists(flux_fp)) {
  df_flux <- read.csv(flux_fp)
  if (all(c("DateTime","ENSO","season") %in% names(df_flux))) {
    df_flux$DateTime <- lubridate::ymd_hms(df_flux$DateTime, quiet = TRUE, tz = tz_local)
    enso_map <- df_flux %>%
      mutate(
        Date   = as.Date(DateTime),
        ENSO   = normalize_enso(ENSO),
        Season = normalize_season(season)
      ) %>%
      select(Date, ENSO, Season) %>%
      distinct()
  } else enso_map <- NULL
} else enso_map <- NULL

if (!is.null(enso_map)) {
  pr_labeled <- pr_daily %>%
    left_join(enso_map, by = "Date") %>%
    mutate(
      Season = ifelse(is.na(Season), ifelse(month(Date) %in% 5:10, "Dry", "Wet"), Season),
      ENSO   = ifelse(is.na(ENSO), "Neutral", ENSO),
      ENSO   = normalize_enso(ENSO),
      Season = normalize_season(Season)
    )
} else {
  pr_labeled <- pr_daily %>%
    mutate(
      DateTime0 = as.POSIXct(paste0(Date, " 00:00:00"), tz = tz_local),
      Season    = ifelse(month(Date) %in% 5:10, "Dry", "Wet"),
      ENSO      = "Neutral",
      ENSO      = ifelse(
        (DateTime0 > ymd_hms("2018-09-01 00:00:00") & DateTime0 < ymd_hms("2019-07-31 23:59:59")) |
          (DateTime0 > ymd_hms("2023-04-01 00:00:00") & DateTime0 < ymd_hms("2024-05-31 23:59:59")),
        "El Niño", ENSO),
      ENSO      = ifelse(
        (DateTime0 > ymd_hms("2016-07-01 00:00:00") & DateTime0 < ymd_hms("2017-01-31 23:59:59")) |
          (DateTime0 > ymd_hms("2017-09-01 00:00:00") & DateTime0 < ymd_hms("2018-05-31 23:59:59")) |
          (DateTime0 > ymd_hms("2020-07-01 00:00:00") & DateTime0 < ymd_hms("2023-02-28 23:59:59")),
        "La Niña", ENSO),
      ENSO      = normalize_enso(ENSO),
      Season    = normalize_season(Season)
    ) %>%
    select(-DateTime0)
}

# yearly precipitation by Season × ENSO with canonical labels
precip_yearly <- pr_labeled %>%
  mutate(Year = year(Date)) %>%
  group_by(Year, Season, ENSO) %>%
  summarise(precip_mm_day = mean(precip_day_mm, na.rm = TRUE), .groups = "drop") %>%
  mutate(Group = paste(ENSO, "–", tolower(Season))) %>%
  select(Group, Year, precip_mm_day)

# join into the master table
summary_table_yearly_filt <- summary_table_yearly_filt %>%
  left_join(precip_yearly, by = c("Group","Year"))

# guard: the precip join must not fan out rows (one row per Group x Year)
stopifnot(!any(duplicated(summary_table_yearly_filt[, c("Group", "Year")])))

# write the by-year parameter tables WITH precipitation merged in, so precip_mm_day
# reaches the LMM input (code/06). (paths defined above, before the precip block.)
write.csv(summary_table_yearly_filt, file = out_csv,         row.names = FALSE)
write.csv(summary_table_yearly_filt, file = out_csv_with_ts, row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# season-wide means (for dashed lines) and significance stars
season_means <- summary_table_yearly_filt %>%
  group_by(Season) %>%
  summarise(
    mean_Ta     = mean(Ta, na.rm = TRUE),
    mean_VPD    = mean(VPD, na.rm = TRUE),
    mean_SWC    = mean(SWC, na.rm = TRUE),
    mean_PAR    = mean(PAR, na.rm = TRUE),
    mean_Ts     = mean(Ts, na.rm = TRUE),
    mean_Precip = mean(precip_mm_day, na.rm = TRUE),
    .groups     = "drop"
  )

season_means_Ta     <- season_means %>% select(Season, mean_val = mean_Ta)
season_means_VPD    <- season_means %>% select(Season, mean_val = mean_VPD)
season_means_SWC    <- season_means %>% select(Season, mean_val = mean_SWC)
season_means_PAR    <- season_means %>% select(Season, mean_val = mean_PAR)
season_means_Ts     <- season_means %>% select(Season, mean_val = mean_Ts)
season_means_Precip <- season_means %>% select(Season, mean_val = mean_Precip)

sig_Ta     <- sig_tests(summary_table_yearly_filt, Ta)
sig_VPD    <- sig_tests(summary_table_yearly_filt, VPD)
sig_SWC    <- sig_tests(summary_table_yearly_filt, SWC)
sig_PAR    <- sig_tests(summary_table_yearly_filt, PAR)
sig_Ts     <- sig_tests(summary_table_yearly_filt, Ts)
sig_Precip <- sig_tests(summary_table_yearly_filt, precip_mm_day)

# ───────────────────────────────────────────────────────────────────────────────
# statistical tests for meteorological variables
# 1) wet vs dry season
# 2) ENSO differences within each season
# ───────────────────────────────────────────────────────────────────────────────

met_vars <- c(
  "Ta",
  "VPD",
  "SWC",
  "PAR",
  "Ts",
  "precip_mm_day"
)

met_var_labels <- c(
  Ta = "Air temperature",
  VPD = "VPD",
  SWC = "Soil water content",
  PAR = "PAR",
  Ts = "Soil temperature",
  precip_mm_day = "Precipitation"
)

# wet vs dry season tests
season_tests <- bind_rows(lapply(met_vars, function(v) {
  
  dat <- summary_table_yearly_filt %>%
    select(Season, value = all_of(v)) %>%
    filter(!is.na(value), !is.na(Season))
  
  if (length(unique(dat$Season)) < 2) {
    return(skip_test_row(v, met_var_labels[[v]], "Welch two-sample t-test", "Wet vs Dry", "fewer than two seasons"))
  }

  group_n <- table(dat$Season)
  wet_n <- unname(group_n[["Wet"]]); if (is.na(wet_n)) wet_n <- 0L
  dry_n <- unname(group_n[["Dry"]]); if (is.na(dry_n)) dry_n <- 0L

  if (wet_n < 2 || dry_n < 2) {
    return(skip_test_row(
      v,
      met_var_labels[[v]],
      "Welch two-sample t-test",
      "Wet vs Dry",
      paste0("not enough observations: Wet n=", wet_n, ", Dry n=", dry_n)
    ) %>% mutate(wet_n = wet_n, dry_n = dry_n))
  }
  
  tt <- tryCatch(t.test(value ~ Season, data = dat), error = function(e) e)
  if (inherits(tt, "error")) {
    return(skip_test_row(
      v,
      met_var_labels[[v]],
      "Welch two-sample t-test",
      "Wet vs Dry",
      tt$message
    ) %>% mutate(wet_n = wet_n, dry_n = dry_n))
  }
  
  tibble(
    variable = v,
    variable_label = met_var_labels[[v]],
    test = "Welch two-sample t-test",
    comparison = "Wet vs Dry",
    wet_mean = mean(dat$value[dat$Season == "Wet"], na.rm = TRUE),
    dry_mean = mean(dat$value[dat$Season == "Dry"], na.rm = TRUE),
    wet_n = sum(dat$Season == "Wet"),
    dry_n = sum(dat$Season == "Dry"),
    statistic = unname(tt$statistic),
    df = unname(tt$parameter),
    p_value = tt$p.value,
    note = NA_character_
  )
}))

# adjust p-values across variables for the wet vs dry tests
season_tests <- season_tests %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    signif = p_to_stars(p_value),
    signif_adj = p_to_stars(p_adj_BH)
  )

# ENSO tests within each season
enso_anova_tests <- bind_rows(lapply(met_vars, function(v) {
  
  bind_rows(lapply(levels(summary_table_yearly_filt$Season), function(ss) {
    
    dat <- summary_table_yearly_filt %>%
      filter(Season == ss) %>%
      select(ENSO, value = all_of(v)) %>%
      filter(!is.na(value), !is.na(ENSO))
    
    group_n <- table(dat$ENSO)
    usable_groups <- names(group_n[group_n > 0])

    if (length(usable_groups) < 2 || nrow(dat) <= length(usable_groups)) {
      return(skip_test_row(
        v,
        met_var_labels[[v]],
        "one-way ANOVA",
        "ENSO phases within season",
        paste0("not enough observations by ENSO: ", paste(names(group_n), group_n, sep = " n=", collapse = "; ")),
        Season = ss
      ))
    }
    
    fit <- tryCatch(aov(value ~ ENSO, data = dat), error = function(e) e)
    if (inherits(fit, "error")) {
      return(skip_test_row(v, met_var_labels[[v]], "one-way ANOVA", "ENSO phases within season", fit$message, Season = ss))
    }
    aov_tab <- summary(fit)[[1]]
    
    tibble(
      variable = v,
      variable_label = met_var_labels[[v]],
      Season = ss,
      test = "one-way ANOVA",
      comparison = "ENSO phases within season",
      df_enso = aov_tab["ENSO", "Df"],
      df_residual = aov_tab["Residuals", "Df"],
      statistic = aov_tab["ENSO", "F value"],
      p_value = aov_tab["ENSO", "Pr(>F)"],
      n = nrow(dat),
      note = NA_character_
    )
  }))
}))

# adjust p-values across all variable × season ENSO tests
enso_anova_tests <- enso_anova_tests %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    signif = p_to_stars(p_value),
    signif_adj = p_to_stars(p_adj_BH)
  )

# post-hoc pairwise ENSO comparisons within each season
enso_tukey_tests <- bind_rows(lapply(met_vars, function(v) {
  
  bind_rows(lapply(levels(summary_table_yearly_filt$Season), function(ss) {
    
    dat <- summary_table_yearly_filt %>%
      filter(Season == ss) %>%
      select(ENSO, value = all_of(v)) %>%
      filter(!is.na(value), !is.na(ENSO))
    
    group_n <- table(dat$ENSO)
    usable_groups <- names(group_n[group_n > 0])

    if (length(usable_groups) < 3 || nrow(dat) <= length(usable_groups)) {
      return(skip_test_row(
        v,
        met_var_labels[[v]],
        "Tukey HSD",
        "ENSO pairwise comparisons",
        paste0("not enough observations by ENSO: ", paste(names(group_n), group_n, sep = " n=", collapse = "; ")),
        Season = ss
      ))
    }
    
    fit <- tryCatch(aov(value ~ ENSO, data = dat), error = function(e) e)
    if (inherits(fit, "error")) {
      return(skip_test_row(v, met_var_labels[[v]], "Tukey HSD", "ENSO pairwise comparisons", fit$message, Season = ss))
    }
    tk <- tryCatch(TukeyHSD(fit, "ENSO")$ENSO, error = function(e) e)
    if (inherits(tk, "error")) {
      return(skip_test_row(v, met_var_labels[[v]], "Tukey HSD", "ENSO pairwise comparisons", tk$message, Season = ss))
    }
    
    as.data.frame(tk) %>%
      tibble::rownames_to_column("comparison") %>%
      as_tibble() %>%
      mutate(
        variable = v,
        variable_label = met_var_labels[[v]],
        Season = ss,
        test = "Tukey HSD",
        note = NA_character_,
        .before = 1
      ) %>%
      rename(
        diff = diff,
        lower = lwr,
        upper = upr,
        p_value = `p adj`
      )
  }))
}))

enso_tukey_tests <- enso_tukey_tests %>%
  mutate(
    signif = p_to_stars(p_value)
  )

# save statistical results
write.csv(
  season_tests,
  file.path(output_path, "meteorological_tests_wet_vs_dry.csv"),
  row.names = FALSE
)

write.csv(
  enso_anova_tests,
  file.path(output_path, "meteorological_tests_ENSO_within_season_ANOVA.csv"),
  row.names = FALSE
)

write.csv(
  enso_tukey_tests,
  file.path(output_path, "meteorological_tests_ENSO_within_season_TukeyHSD.csv"),
  row.names = FALSE
)

cat(
  "\nsaved meteorological statistical tests:\n",
  "  -", file.path(output_path, "meteorological_tests_wet_vs_dry.csv"), "\n",
  "  -", file.path(output_path, "meteorological_tests_ENSO_within_season_ANOVA.csv"), "\n",
  "  -", file.path(output_path, "meteorological_tests_ENSO_within_season_TukeyHSD.csv"), "\n"
)

# ───────────────────────────────────────────────────────────────────────────────
# colors by year (fixed palette)
years <- sort(unique(summary_table_yearly_filt$Year))
custom_colors <- c("black","gray40","blue","red","purple","darkgreen","orange","brown","cyan","magenta")
year_colors <- setNames(rep(custom_colors, length.out = length(years)), years)

# plots (a–f): Ta, VPD, SWC, PAR, Precip, Ts
p_Ta <- plot_var_by_ENSO_faceted(summary_table_yearly_filt, Ta, expression(italic(T)[plain(a)] ~ "(°C)"), season_means_Ta, year_colors) +
  labs(tag = "a)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_Ta <- add_sig(p_Ta, sig_Ta)

p_VPD <- plot_var_by_ENSO_faceted(summary_table_yearly_filt, VPD, "VPD (kPa)", season_means_VPD, year_colors) +
  labs(tag = "b)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_VPD <- add_sig(p_VPD, sig_VPD)

p_SWC <- plot_var_by_ENSO_faceted(summary_table_yearly_filt, SWC, "SWC (0–1 m, cm)", season_means_SWC, year_colors) +
  labs(tag = "c)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_SWC <- add_sig(p_SWC, sig_SWC)

p_PAR <- plot_var_by_ENSO_faceted(summary_table_yearly_filt, PAR, "PAR (µmol m⁻² s⁻¹)", season_means_PAR, year_colors) +
  labs(tag = "d)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_PAR <- add_sig(p_PAR, sig_PAR)

p_Precip <- plot_var_by_ENSO_faceted(summary_table_yearly_filt, precip_mm_day,
                                     "Precipitation (mm d⁻¹)", season_means_Precip, year_colors) +
  labs(tag = "e)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_Precip <- add_sig(p_Precip, sig_Precip)

p_Ts <- plot_var_by_ENSO_faceted(
  summary_table_yearly_filt, Ts,
  expression(italic(T)[plain(s)] ~ "(°C)"), season_means_Ts, year_colors
) +
  labs(tag = "f)") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.01, 0.99))
p_Ts <- add_sig(p_Ts, sig_Ts)

# ───────────────────────────────────────────────────────────────────────────────
# combine and save (year-colored figure with sig)
combined_fig_8 <- (p_Ta | p_VPD) /
  (p_SWC | p_PAR) /
  (p_Precip | p_Ts) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_fig_8 <- combined_fig_8 +
  plot_annotation(
    title = "Meteorological and soil variables at the Tambopata station\n(*, **, *** indicate Wet vs Dry differences by two-sample t-test)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  )

ggsave(file.path(paper_figures_path, "Fig1_env_ENSO_faceted_SWC_Ts.png"),
       combined_fig_8, width = 12, height = 14, dpi = 300)

cat("\nwrote tables:\n  - yearly:", out_csv,
    "\n  - joined precipitation source: IMERG (30-min → daily mean by Year × Season × ENSO, with canonical labels)",
    "\nsaved combined figure (with significance):\n  -",
    file.path(paper_figures_path, "Fig1_env_ENSO_faceted_SWC_Ts.(png|pdf)"),
    "\n")

# ───────────────────────────────────────────────────────────────────────────────
# FIGURE 2: ENSO means ± SD (collapsed across years)
# each point = mean of yearly values within ENSO × Season
# error bars = ± 1 SD across years
# ───────────────────────────────────────────────────────────────────────────────

enso_summary <- summary_table_yearly_filt %>%
  group_by(Season, ENSO) %>%
  summarise(
    Ta_mean      = mean(Ta, na.rm = TRUE),
    Ta_sd        = sd(Ta, na.rm = TRUE),
    
    VPD_mean     = mean(VPD, na.rm = TRUE),
    VPD_sd       = sd(VPD, na.rm = TRUE),
    
    SWC_mean     = mean(SWC, na.rm = TRUE),
    SWC_sd       = sd(SWC, na.rm = TRUE),
    
    PAR_mean     = mean(PAR, na.rm = TRUE),
    PAR_sd       = sd(PAR, na.rm = TRUE),
    
    Ts_mean      = mean(Ts, na.rm = TRUE),
    Ts_sd        = sd(Ts, na.rm = TRUE),
    
    Precip_mean  = mean(precip_mm_day, na.rm = TRUE),
    Precip_sd    = sd(precip_mm_day, na.rm = TRUE),
    
    n_years      = sum(!is.na(Ta)),
    .groups      = "drop"
  ) %>%
  mutate(
    Season = factor(Season, levels = c("Wet", "Dry")),
    ENSO   = factor(ENSO, levels = c("El Niño", "La Niña", "Neutral"))
  )

season_means_enso <- enso_summary %>%
  group_by(Season) %>%
  summarise(
    Ta      = mean(Ta_mean, na.rm = TRUE),
    VPD     = mean(VPD_mean, na.rm = TRUE),
    SWC     = mean(SWC_mean, na.rm = TRUE),
    PAR     = mean(PAR_mean, na.rm = TRUE),
    Ts      = mean(Ts_mean, na.rm = TRUE),
    Precip  = mean(Precip_mean, na.rm = TRUE),
    .groups = "drop"
  )

plot_mean_sd_by_ENSO <- function(df, mean_col, sd_col, ylab, season_means_df) {
  
  m <- rlang::ensym(mean_col)
  s <- rlang::ensym(sd_col)
  
  ggplot(df, aes(x = ENSO, y = !!m)) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(
      aes(ymin = !!m - !!s, ymax = !!m + !!s),
      width = 0.18,
      linewidth = 0.6
    ) +
    geom_hline(
      data = season_means_df,
      aes(yintercept = mean_val),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Season, nrow = 1) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0)
    )
}

season_means_Ta_enso     <- season_means_enso %>% select(Season, mean_val = Ta)
season_means_VPD_enso    <- season_means_enso %>% select(Season, mean_val = VPD)
season_means_SWC_enso    <- season_means_enso %>% select(Season, mean_val = SWC)
season_means_PAR_enso    <- season_means_enso %>% select(Season, mean_val = PAR)
season_means_Precip_enso <- season_means_enso %>% select(Season, mean_val = Precip)
season_means_Ts_enso     <- season_means_enso %>% select(Season, mean_val = Ts)

p2_Ta <- plot_mean_sd_by_ENSO(enso_summary, Ta_mean, Ta_sd, expression(italic(T)[plain(a)] ~ "(°C)"), season_means_Ta_enso) +
  labs(tag = "a)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

p2_VPD <- plot_mean_sd_by_ENSO(enso_summary, VPD_mean, VPD_sd, "VPD (kPa)", season_means_VPD_enso) +
  labs(tag = "b)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

p2_SWC <- plot_mean_sd_by_ENSO(enso_summary, SWC_mean, SWC_sd, "SWC (0–1 m, cm)", season_means_SWC_enso) +
  labs(tag = "c)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

p2_PAR <- plot_mean_sd_by_ENSO(enso_summary, PAR_mean, PAR_sd, "PAR (µmol m⁻² s⁻¹)", season_means_PAR_enso) +
  labs(tag = "d)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

p2_Precip <- plot_mean_sd_by_ENSO(enso_summary, Precip_mean, Precip_sd, "Precipitation (mm d⁻¹)", season_means_Precip_enso) +
  labs(tag = "e)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

p2_Ts <- plot_mean_sd_by_ENSO(enso_summary, Ts_mean, Ts_sd, expression(italic(T)[plain(s)] ~ "(°C)"), season_means_Ts_enso) +
  labs(tag = "f)") +
  theme(plot.tag = element_text(face = "bold", size = 14))

combined_fig_mean_sd <- (p2_Ta | p2_VPD) /
  (p2_SWC | p2_PAR) /
  (p2_Precip | p2_Ts)

combined_fig_mean_sd <- combined_fig_mean_sd +
  plot_annotation(
    title = "Meteorological and soil variables at the Tambopata station\n(Means ± SD across years for each ENSO phase)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  )

ggsave(file.path(paper_figures_path, "Fig1_env_ENSO_faceted_MEAN_SD.png"),
       combined_fig_mean_sd, width = 12, height = 14, dpi = 300)

cat("\nsaved ENSO-mean figure (mean ± SD across years):\n  -",
    file.path(paper_figures_path, "Fig1_env_ENSO_faceted_MEAN_SD.(png|pdf)"),
    "\n")

# ───────────────────────────────────────────────────────────────────────────────
# FIGURE 3: ENSO means ± SD + individual yearly points (transparent)
# ───────────────────────────────────────────────────────────────────────────────

yearly_for_overlay <- summary_table_yearly_filt %>%
  mutate(
    Season = factor(Season, levels = c("Wet", "Dry")),
    ENSO   = factor(ENSO, levels = c("El Niño", "La Niña", "Neutral"))
  )

plot_mean_sd_with_years_by_ENSO <- function(mean_df, yearly_df, var_yearly, mean_col, sd_col, ylab, season_means_df) {
  
  v  <- rlang::ensym(var_yearly)
  m  <- rlang::ensym(mean_col)
  s  <- rlang::ensym(sd_col)
  
  ggplot() +
    geom_point(
      data = yearly_df,
      aes(x = ENSO, y = !!v),
      position = position_jitter(width = 0.14, height = 0),
      alpha = 0.35,
      size = 2.2,
      color = "gray30"
    ) +
    geom_point(
      data = mean_df,
      aes(x = ENSO, y = !!m),
      size = 3,
      color = "black"
    ) +
    geom_errorbar(
      data = mean_df,
      aes(x = ENSO, ymin = !!m - !!s, ymax = !!m + !!s),
      width = 0.18,
      linewidth = 0.6,
      color = "black"
    ) +
    geom_hline(
      data = season_means_df,
      aes(yintercept = mean_val),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Season, nrow = 1) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      # 30 deg: at 14.4 pt the three ENSO categories fill ~90% of each narrow
      # facet and ran together horizontally
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

p3_Ta <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = Ta, mean_col = Ta_mean, sd_col = Ta_sd,
  ylab = expression(italic(T)[plain(a)] ~ "(°C)"), season_means_df = season_means_Ta_enso
) +
  labs(tag = "a)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

p3_VPD <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = VPD, mean_col = VPD_mean, sd_col = VPD_sd,
  ylab = "VPD (kPa)", season_means_df = season_means_VPD_enso
) +
  labs(tag = "b)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

p3_SWC <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = SWC, mean_col = SWC_mean, sd_col = SWC_sd,
  ylab = "SWC (0–1 m, cm)", season_means_df = season_means_SWC_enso
) +
  labs(tag = "c)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

p3_PAR <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = PAR, mean_col = PAR_mean, sd_col = PAR_sd,
  ylab = "PAR (µmol m⁻² s⁻¹)", season_means_df = season_means_PAR_enso
) +
  labs(tag = "d)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

p3_Precip <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = precip_mm_day, mean_col = Precip_mean, sd_col = Precip_sd,
  ylab = "Precipitation (mm d⁻¹)", season_means_df = season_means_Precip_enso
) +
  labs(tag = "e)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

p3_Ts <- plot_mean_sd_with_years_by_ENSO(
  mean_df = enso_summary, yearly_df = yearly_for_overlay,
  var_yearly = Ts, mean_col = Ts_mean, sd_col = Ts_sd,
  ylab = expression(italic(T)[plain(s)] ~ "(°C)"), season_means_df = season_means_Ts_enso
) +
  labs(tag = "f)") +
  theme(plot.tag = element_text(face = "bold", size = 21))

combined_fig_mean_sd_years <- (p3_Ta | p3_VPD) /
  (p3_SWC | p3_PAR) /
  (p3_Precip | p3_Ts)

combined_fig_mean_sd_years <- combined_fig_mean_sd_years +
  plot_annotation(
    # the parenthetical is wrapped onto its own two lines: at 22 pt the old
    # single 87-character line was ~4050 px wide against a 3600 px canvas and
    # ran off both edges.
    title = "Meteorological and soil variables at the Tambopata station\n(Means ± SD across years for each ENSO phase;\ntransparent points show individual years)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  )

ggsave(
  file.path(paper_figures_path, "Fig1_env_ENSO_faceted_MEAN_SD_with_years.png"),
  combined_fig_mean_sd_years,
  width = 12, height = 14, dpi = 300
)

cat(
  "\nsaved ENSO-mean figure (mean ± SD + transparent yearly points):\n  -",
  file.path(paper_figures_path, "Fig1_env_ENSO_faceted_MEAN_SD_with_years.(png|pdf)"),
  "\n"
)

# ───────────────────────────────────────────────────────────────────────────────
# Supplement Figure S2: monthly precipitation with ENSO-phase shading
#
# This block is purely ADDITIVE: it reuses pr_daily (built above from the IMERG
# 30-min file) and writes one new PNG plus one CSV. Nothing above is touched.
#
# The ENSO bands deliberately come from the labelled flux dataset rather than
# from the fallback date rules at the pr_labeled step further up, so that the
# bands land in exactly the same places as every other ENSO-shaded figure (the
# contiguous-run collapse below is the same one used at code/04:1302-1315).
# ───────────────────────────────────────────────────────────────────────────────

# Daily totals are rebuilt here from imerg30 rather than reusing pr_daily above,
# because pr_daily assigns each half-hour with as.Date(DateTime), whose default
# is tz = "UTC": at UTC-5 that pushes the 19:00-24:00 local block into the NEXT
# day, moving up to five hours of rain across every day and month boundary
# (e.g. Jan 2017 loses its last 10 half-hours to February). Pre-existing, and out
# of scope to fix here because pr_daily feeds precip_mm_day -> Fig. S3 panel (e)
# and the light-response LMM input; see the note in the revision guide. S2 needs
# the local-date assignment its caption claims, so it does its own aggregation.
precip_daily_local <- imerg30 %>%
  mutate(Date = as.Date(DateTime, tz = tz_local)) %>%
  group_by(Date) %>%
  summarise(
    precip_day_mm = sum(precip_mm_30min, na.rm = TRUE),
    n_halfhours   = dplyr::n(),
    .groups = "drop"
  )

precip_monthly <- precip_daily_local %>%
  mutate(Month = lubridate::floor_date(Date, "month")) %>%
  group_by(Month) %>%
  summarise(
    precip_mm    = sum(precip_day_mm, na.rm = TRUE),
    n_days       = dplyr::n(),
    n_halfhours  = sum(n_halfhours),
    .groups = "drop"
  ) %>%
  mutate(
    days_in_month = lubridate::days_in_month(Month),
    coverage      = n_halfhours / (48 * as.numeric(days_in_month)),
    # May-Oct / Nov-Apr, the split used throughout the paper (cf. code/04:1165)
    season = factor(
      if_else(lubridate::month(Month) %in% 5:10, "dry", "wet"),
      levels = c("dry", "wet")
    )
  )

stopifnot(nrow(precip_monthly) == 96L)

short_months <- precip_monthly %>% filter(coverage < 0.99)
if (nrow(short_months)) {
  cat("\nS2 coverage note - months below 99% of their half-hours:\n")
  for (i in seq_len(nrow(short_months))) {
    cat(sprintf(
      "  %s  %d/%d half-hours (%.1f%%)\n",
      format(short_months$Month[i], "%Y-%m"),
      short_months$n_halfhours[i], 48L * as.integer(short_months$days_in_month[i]),
      100 * short_months$coverage[i]
    ))
  }
} else {
  cat("\nS2 coverage note: all 96 months complete\n")
}

# ENSO shading rectangles: contiguous runs of the canonical flux-site labels.
# These end with the flux record (15 Oct 2024); the remaining weeks of 2024 are
# left unshaded rather than assigned a phase the dataset does not label.
enso_src <- file.path(input_folder, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
if (!file.exists(enso_src)) stop("S2: labelled ENSO dataset not found: ", enso_src)

# NOTE: parse tv_dt, not the DateTime column. tv_dt is local ("01-Jan-2017
# 00:30:00") while DateTime is UTC ("2017-01-01T05:30:00Z"); code/04:1132 uses
# dmy_hms(tv_dt, tz = tz_local), and the bands must match it exactly.
enso_rect_s2 <- data.table::fread(enso_src, select = c("tv_dt", "ENSO")) %>%
  as_tibble() %>%
  mutate(DateTime = lubridate::dmy_hms(as.character(tv_dt), tz = tz_local)) %>%
  filter(!is.na(DateTime), !is.na(ENSO)) %>%
  arrange(DateTime) %>%
  mutate(
    ENSO   = as.character(ENSO),
    change = (ENSO != dplyr::lag(ENSO, default = dplyr::first(ENSO))),
    grp    = cumsum(change)
  ) %>%
  group_by(grp) %>%
  summarise(
    ENSO = dplyr::first(ENSO),
    xmin = as.Date(min(DateTime)),
    xmax = as.Date(max(DateTime)),
    .groups = "drop"
  ) %>%
  mutate(ENSO = factor(ENSO, levels = c("El Nino", "La Nina", "neutral")))

cat("\nS2 ENSO bands (from the flux-site labels):\n")
for (i in seq_len(nrow(enso_rect_s2))) {
  cat(sprintf("  %-8s %s -> %s\n", as.character(enso_rect_s2$ENSO[i]),
              enso_rect_s2$xmin[i], enso_rect_s2$xmax[i]))
}

p_precip_s2 <- ggplot() +
  geom_rect(
    data = enso_rect_s2,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE, alpha = 0.9
  ) +
  scale_fill_manual(
    values = enso_fill, drop = FALSE, name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  # second fill scale so the bars can carry season without colliding with the
  # ENSO bands (same idiom as code/14:298)
  ggnewscale::new_scale_fill() +
  geom_col(
    data = precip_monthly,
    aes(x = Month + 15, y = precip_mm, fill = season),
    width = 26, colour = NA
  ) +
  scale_fill_manual(values = season_cols, drop = FALSE, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = expression(Precipitation ~ (mm ~ month^{-1}))) +
  theme_bw(base_size = 9) +
  theme(
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    # the 1% x-expansion above puts the final break (2025) only ~20 px inside the
    # panel border, and the default 5.5 pt margin left the "2025" label 3 px short
    # of fitting, so its last digit was clipped at the canvas edge. Widen the
    # right margin rather than the expansion: more expansion would open a gap
    # inside the panel next to the already-unshaded post-Oct-2024 stretch, which
    # would read as missing data.
    plot.margin       = margin(t = 5.5, r = 14, b = 5.5, l = 5.5, unit = "pt")
  )

out_fig_s2 <- file.path(paper_figures_path, "FigS_monthly_precipitation_ENSO.png")
ggsave(out_fig_s2, p_precip_s2, width = 7.10, height = 3.20, dpi = 300, bg = "white")
cat("saved:", out_fig_s2, "\n")

out_csv_s2 <- file.path(out_base, "table_monthly_precipitation.csv")
precip_monthly %>%
  transmute(
    ym = format(Month, "%Y-%m"),
    season = as.character(season),
    precip_mm = round(precip_mm, 2),
    n_days, days_in_month, n_halfhours
  ) %>%
  write.csv(out_csv_s2, row.names = FALSE)
cat("saved:", out_csv_s2, "\n")

precip_annual_s2 <- precip_monthly %>%
  mutate(Year = lubridate::year(Month)) %>%
  group_by(Year) %>%
  summarise(precip_mm_year = sum(precip_mm), .groups = "drop")
cat("\nS2 annual precipitation totals (mm):\n")
print(as.data.frame(precip_annual_s2), row.names = FALSE)
cat(sprintf(
  "  mean %.0f mm/yr; monthly max %.1f (%s), min %.1f (%s)\n",
  mean(precip_annual_s2$precip_mm_year),
  max(precip_monthly$precip_mm),
  format(precip_monthly$Month[which.max(precip_monthly$precip_mm)], "%Y-%m"),
  min(precip_monthly$precip_mm),
  format(precip_monthly$Month[which.min(precip_monthly$precip_mm)], "%Y-%m")
))
