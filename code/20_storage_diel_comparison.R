#!/usr/bin/env Rscript

# Supplementary figure: does FC + modelled storage reproduce the diel cycle of
# measured NEE (FC + profile-measured storage)?
#
# Background. NEE = FC + SC. The storage term comes from two independent sources:
#   SC_profiler  profile measurement, computed from the 8-level CO2 profile
#   SC_model     DNN v31 reconstruction (raw trace SC_MULTIPLE_DNN_v31)
# PETNR_SecondStage.ini selects SC_model before 2024 and SC_profiler from 2024
# onward, so the two are never used together in the published NEE. Their raw
# traces do overlap, which is what makes this comparison possible.
#
# Sample. 2024 only. It is the single year that is both well covered by the
# profiler and the year whose published NEE actually uses profile storage, so
# here FC + SC_profiler IS the published measured NEE and FC + SC_model is the
# direct counterfactual. 2021 and 2022 support the same conclusion and are
# reported in the stats table. 2023 is reported but flagged: its profiler trace
# is corrupt (see below) and it is excluded from the figure.
#
# THIS COMPARISON IS IN-SAMPLE. The DNN was trained on the profile record itself
# (Manuscript_v09.md sec 2.5: "we used the period with overlapping single-point
# storage estimates from EddyPro and multi-level storage estimates from the
# profile measurements to train a deep neural network model"), and the adopted
# batch-wise 60/20/20 split interleaves training samples through the whole
# 2021-2024 target record. The profiler only ever ran 2021-2025, so no
# out-of-sample year exists. Read the figure as a consistency check on the
# reconstruction -- it shows the model reproduces the diel shape it was fitted
# to -- and NOT as independent validation. Per-year fit quality is flat across
# 2021-2024 (R2 0.69/0.68/0.73/0.76), as expected when every year is in-sample.
#
# Note "measured" here refers to the storage term, not to gap filling: a
# non-gap-filled NEE value can still carry a model-derived storage term.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(lubridate)
  library(readr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# config

tz_local <- "America/Lima"

# the year the figure is built from
focus_year <- 2024L

# years carried in the stats table for context
context_years <- c(2021L, 2022L, 2023L)

# the 2023 profiler campaign (2023-09-24 to 2023-10-24) returned a median of
# +32 umol m-2 s-1 with more than half its values beyond |50|, which is not
# physically plausible for a storage flux. FirstStage never ingested it (it only
# reads SC_profiler from 2024-01-01), so it never reached the published NEE.
corrupt_profiler_years <- 2023L

# FirstStage minMax on SC_model (PETNR_FirstStage.ini). the export deliberately
# writes the raw trace, so the bound is imposed here where it is visible.
sc_model_abs_max <- 100

# day/night convention used throughout the project: night = Rg < 10 W m-2,
# day = PAR > 100 umol m-2 s-1 and Ta > 17 C. the 10-44 W m-2 twilight band is
# in neither class.
night_sw_max  <- 10
day_par_min   <- 100
day_ta_min    <- 17

# series colours. season (blue/red) and ENSO (gold/purple/grey) are already
# taken in palette.R, so the two storage sources get a disjoint pair, with
# linetype as the redundant non-colour channel.
series_levels    <- c("FC + profile SC (measured)", "FC + modelled SC")
series_cols      <- setNames(c("black", "#009E73"), series_levels)
series_linetypes <- setNames(c("solid", "22"), series_levels)

sc_series_levels <- c("SC profile (measured)", "SC modelled (DNN v31)")
sc_series_cols   <- setNames(c("black", "#009E73"), sc_series_levels)
sc_series_lty    <- setNames(c("solid", "22"), sc_series_levels)

# publication figure text settings, matching code/04_diel_cycles_carbon_balance.R
fig_base_size <- 16

theme_pub <- function(base_size = fig_base_size, legend_position = "right") {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "plain", hjust = 0.5),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1),
      strip.text = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      legend.position = legend_position
    )
}

# x-axis ticks for half-hour binned diel plots
x_breaks <- c("02:00:00", "08:00:00", "14:00:00", "20:00:00")
x_labels <- c("02:00",    "08:00",    "14:00",    "20:00")

halfhour_levels <- sprintf(
  "%02d:%02d:00",
  rep(0:23, each = 2),
  rep(c(0, 30), times = 24)
)

flux_csv  <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
out_diel  <- file.path(output_path, "storage_diel_profiler_vs_model_2024.csv")
out_stats <- file.path(output_path, "storage_profiler_vs_model_stats_2024.csv")
out_fig   <- file.path(paper_figures_path, "FigS_storage_diel_profiler_vs_model.png")

# ─────────────────────────────────────────────────────────────────────────────
# load and check

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv,
       "\nRun code/01_prepare_flux_timeseries.R first.")
}

required_cols <- c("tv_dt", "FC", "SC_profiler", "SC_model", "USTAR",
                   "SW_IN_1_1_1", "PPFD_IN_1_1_1", "TA_1_1_1")

df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("missing required columns in the enriched flux table: ",
       paste(missing_cols, collapse = ", "),
       "\nSC_profiler and SC_model are written by 'from Matlab/Export_data.m'",
       " and carried through by code/01_prepare_flux_timeseries.R.",
       "\nIf they are absent, re-run the MATLAB export and then script 01.")
}

df <- df_raw %>%
  mutate(
    DateTime = dmy_hms(.data[["tv_dt"]], tz = tz_local),
    year     = year(DateTime),
    hour     = factor(format(DateTime, "%H:%M:%S"), levels = halfhour_levels),
    FC          = as.numeric(FC),
    SC_profiler = as.numeric(SC_profiler),
    USTAR       = as.numeric(USTAR),
    # impose the FirstStage bound on the modelled storage only
    SC_model = ifelse(is.finite(SC_model) & abs(SC_model) <= sc_model_abs_max,
                      as.numeric(SC_model), NA_real_),
    # u* screening, using each year's own REddyProc threshold (ustar_filter.R)
    ustar_ok = is.finite(USTAR) & USTAR >= ustar_threshold_for_year(year),
    daynight = case_when(
      is.finite(SW_IN_1_1_1) & SW_IN_1_1_1 < night_sw_max ~ "night",
      is.finite(PPFD_IN_1_1_1) & is.finite(TA_1_1_1) &
        PPFD_IN_1_1_1 > day_par_min & TA_1_1_1 > day_ta_min ~ "day",
      TRUE ~ "twilight"
    )
  )

# paired sample: both storage sources and FC present on the same half-hour, so
# the two NEE series differ only in their storage term
paired <- df %>%
  filter(is.finite(FC), is.finite(SC_profiler), is.finite(SC_model)) %>%
  mutate(
    NEE_profiler = FC + SC_profiler,
    NEE_model    = FC + SC_model,
    # storage provenance, per the AGENTS.md storage_source vocabulary; both
    # sources are present on every row of this table by construction
    storage_source = "measured_profile+modeled_profile"
  )

if (!nrow(paired)) {
  stop("no half-hour has FC, SC_profiler and SC_model all finite; nothing to compare.")
}

focus <- paired %>% filter(year == focus_year)

if (!nrow(focus)) {
  stop("no paired half-hours in ", focus_year, "; check the SC_profiler and ",
       "SC_model columns in ", flux_csv)
}

message(sprintf("paired half-hours: %d total, %d in %d, across %d days",
                nrow(paired), nrow(focus), focus_year,
                dplyr::n_distinct(as.Date(focus$DateTime))))

# the pairing is limited by the DNN, not by the profiler: report both so the
# sample size is not mistaken for a profiler coverage limit
focus_all <- df %>% filter(year == focus_year)
message(sprintf(
  "  %d finite in %d: SC_profiler = %d, SC_model (after the +/-%d bound) = %d",
  focus_year, focus_year, sum(is.finite(focus_all$SC_profiler)),
  sc_model_abs_max, sum(is.finite(focus_all$SC_model))))

# ─────────────────────────────────────────────────────────────────────────────
# diel means

# mean, standard error and n per half-hour for one variable
diel_summarise <- function(d, ustar_filtered) {
  if (isTRUE(ustar_filtered)) d <- d %>% filter(ustar_ok)
  d %>%
    group_by(hour, .drop = FALSE) %>%
    summarise(
      across(
        c(FC, SC_profiler, SC_model, NEE_profiler, NEE_model),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          se   = ~ sd(.x, na.rm = TRUE) / sqrt(sum(is.finite(.x)))
        ),
        .names = "{.col}_{.fn}"
      ),
      n = sum(is.finite(NEE_profiler) & is.finite(NEE_model)),
      .groups = "drop"
    ) %>%
    mutate(ustar_filtered = ustar_filtered, .before = 1)
}

diel <- bind_rows(
  diel_summarise(focus, FALSE),
  diel_summarise(focus, TRUE)
) %>%
  mutate(year = focus_year, .before = 1)

# ─────────────────────────────────────────────────────────────────────────────
# agreement statistics

# OLS of the measured storage on the modelled storage, plus the error metrics
# that do not depend on which variable is regressed on which
agreement_stats <- function(x_measured, y_modelled) {
  k <- is.finite(x_measured) & is.finite(y_modelled)
  x <- x_measured[k]
  y <- y_modelled[k]
  if (length(x) < 5) {
    return(tibble(n = length(x), slope = NA_real_, intercept = NA_real_,
                  r = NA_real_, R2 = NA_real_, RMSE = NA_real_, bias = NA_real_))
  }
  fit <- lm(x ~ y)
  tibble(
    n         = length(x),
    slope     = unname(coef(fit)[2]),
    intercept = unname(coef(fit)[1]),
    r         = cor(x, y),
    R2        = summary(fit)$r.squared,
    RMSE      = sqrt(mean((x - y)^2)),
    bias      = mean(x - y)
  )
}

# one stats row per (year, subset, u* variant)
stats_row <- function(d, yr, subset_label, ustar_filtered, note) {
  if (isTRUE(ustar_filtered)) d <- d %>% filter(ustar_ok)
  d <- switch(
    subset_label,
    all   = d,
    day   = d %>% filter(daynight == "day"),
    night = d %>% filter(daynight == "night")
  )
  bind_cols(
    tibble(year = yr, subset = subset_label, ustar_filtered = ustar_filtered),
    agreement_stats(d$SC_profiler, d$SC_model),
    tibble(note = note)
  )
}

# every row of this table is an in-sample statistic. the note travels with the
# numbers so they cannot be lifted out of the CSV and quoted as validation.
in_sample_note <- paste(
  "IN-SAMPLE: the DNN was trained on the profile record itself, so this measures",
  "how well the reconstruction reproduces the data it was fitted to, NOT",
  "out-of-sample skill. no out-of-sample year exists (profiler ran 2021-2025 and",
  "2021-2024 all fed training)."
)

focus_stats <- bind_rows(lapply(
  list(
    list("all", FALSE), list("day", FALSE), list("night", FALSE),
    list("all", TRUE),  list("day", TRUE),  list("night", TRUE)
  ),
  function(a) stats_row(focus, focus_year, a[[1]], a[[2]],
                        paste("figure sample.", in_sample_note))
))

context_stats <- bind_rows(lapply(context_years, function(yr) {
  d <- paired %>% filter(year == yr)
  note <- if (yr %in% corrupt_profiler_years) {
    "excluded_corrupt: profiler trace not physically plausible; never ingested by FirstStage"
  } else {
    paste("context only: profiler predates the 2024-01-01 FirstStage ingestion window.",
          in_sample_note)
  }
  if (!nrow(d)) {
    return(bind_cols(tibble(year = yr, subset = "all", ustar_filtered = FALSE),
                     agreement_stats(numeric(0), numeric(0)), tibble(note = note)))
  }
  stats_row(d, yr, "all", FALSE, note)
}))

stats_tbl <- bind_rows(focus_stats, context_stats)

# the quantity the co-author actually asked about: how far apart are the two
# diel NEE curves?
diel_gap <- diel %>%
  filter(!ustar_filtered) %>%
  mutate(gap = NEE_profiler_mean - NEE_model_mean) %>%
  summarise(mean_abs = mean(abs(gap), na.rm = TRUE),
            max_abs  = max(abs(gap), na.rm = TRUE))

stats_tbl <- bind_rows(
  stats_tbl,
  tibble(year = focus_year, subset = "diel_NEE_curve_gap", ustar_filtered = FALSE,
         n = sum(diel$n[!diel$ustar_filtered]),
         slope = NA_real_, intercept = NA_real_, r = NA_real_, R2 = NA_real_,
         RMSE = diel_gap$max_abs, bias = diel_gap$mean_abs,
         note = paste("RMSE column holds max |mean NEE_profiler - mean NEE_model|",
                      "over the 48 half-hour bins; bias column holds the mean of",
                      "that absolute difference (umol m-2 s-1)"))
)

readr::write_csv(diel, out_diel)
readr::write_csv(stats_tbl, out_stats)
cat("saved:", out_diel, "\n")
cat("saved:", out_stats, "\n")

print(as.data.frame(stats_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))

# ─────────────────────────────────────────────────────────────────────────────
# figure

# the DNN output carries no value at the 00:00 half-hour in any year, so that
# bin is empty for both series (the paired sample needs both). drop it rather
# than let ggplot warn about it, and say so on the panel.
empty_bins <- diel %>% filter(!ustar_filtered, n == 0) %>% pull(hour) %>% as.character()
if (length(empty_bins)) {
  message("half-hour bins with no paired data (dropped from the figure): ",
          paste(empty_bins, collapse = ", "))
}

diel_unfiltered <- diel %>% filter(!ustar_filtered, n > 0)

nee_long <- bind_rows(
  diel_unfiltered %>% transmute(hour, series = series_levels[1],
                                mean = NEE_profiler_mean, se = NEE_profiler_se),
  diel_unfiltered %>% transmute(hour, series = series_levels[2],
                                mean = NEE_model_mean, se = NEE_model_se)
) %>% mutate(series = factor(series, levels = series_levels))

sc_long <- bind_rows(
  diel_unfiltered %>% transmute(hour, series = sc_series_levels[1],
                                mean = SC_profiler_mean, se = SC_profiler_se),
  diel_unfiltered %>% transmute(hour, series = sc_series_levels[2],
                                mean = SC_model_mean, se = SC_model_se)
) %>% mutate(series = factor(series, levels = sc_series_levels))

diel_panel <- function(d, cols, ltys, y_label, subtitle) {
  ggplot(d, aes(x = hour, y = mean, colour = series, linetype = series, group = series)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = series),
                alpha = 0.18, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.5, show.legend = FALSE) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    labs(x = "Local time", y = y_label, subtitle = subtitle) +
    theme_pub(legend_position = "top") +
    theme(plot.subtitle = element_text(size = fig_base_size - 2, hjust = 0),
          legend.margin = margin(0, 0, 0, 0))
}

p_nee <- diel_panel(
  nee_long, series_cols, series_linetypes,
  expression(NEE ~ (mu * mol ~ m^-2 ~ s^-1)),
  sprintf("(a) mean diel NEE, %d paired half-hours in %d (%d days); shading is ± 1 SE%s",
          nrow(focus), focus_year, dplyr::n_distinct(as.Date(focus$DateTime)),
          if (length(empty_bins))
            sprintf("; no modelled SC at %s", paste(substr(empty_bins, 1, 5), collapse = ", "))
          else "")
)

p_sc <- diel_panel(
  sc_long, sc_series_cols, sc_series_lty,
  expression(SC ~ (mu * mol ~ m^-2 ~ s^-1)),
  "(b) the storage term alone, where the two series actually differ"
) + theme(plot.margin = margin(5.5, 16, 5.5, 5.5))

s_all <- stats_tbl %>% filter(year == focus_year, subset == "all", !ustar_filtered)
stats_label <- sprintf(
  "n = %d\nslope = %.2f\nintercept = %.2f\nr = %.3f\nR2 = %.3f\nRMSE = %.2f\nbias = %.2f",
  s_all$n, s_all$slope, s_all$intercept, s_all$r, s_all$R2, s_all$RMSE, s_all$bias
)

sc_lim <- range(c(focus$SC_profiler, focus$SC_model), na.rm = TRUE)

p_scatter <- ggplot(focus, aes(x = SC_model, y = SC_profiler)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linetype = "22", linewidth = 0.6) +
  geom_point(alpha = 0.25, size = 1.1, colour = "#009E73") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.12,
           label = stats_label, size = 4.2, lineheight = 1.05) +
  coord_equal(xlim = sc_lim, ylim = sc_lim) +
  labs(x = expression(SC[modelled] ~ (mu * mol ~ m^-2 ~ s^-1)),
       y = expression(SC[profile] ~ (mu * mol ~ m^-2 ~ s^-1)),
       subtitle = "(c) half-hourly agreement") +
  theme_pub(legend_position = "none") +
  theme(plot.subtitle = element_text(size = fig_base_size - 2, hjust = 0))

r2_ustar <- stats_tbl$R2[stats_tbl$year == focus_year & stats_tbl$subset == "all" &
                           stats_tbl$ustar_filtered]

# one entry per paragraph; each is wrapped independently so the caption never
# runs past the figure edge
caption_paragraphs <- c(
  paste0("Panel (c): grey dashed line is 1:1, black line is the OLS fit of SC_profile on ",
         "SC_modelled. FC is wind-sector and precipitation filtered but not storage-corrected ",
         "and not u*-filtered; u*-filtering the pair moves R2 from ",
         sprintf("%.3f to %.3f.", s_all$R2, r2_ustar)),
  paste0("IN-SAMPLE COMPARISON: the DNN was trained on this same profile record (Manuscript ",
         "sec. 2.5), and the profiler never ran outside 2021-2025, so no out-of-sample year ",
         "exists. This shows the reconstruction reproduces the diel shape it was fitted to; it ",
         "is not independent validation."),
  paste0("Sample size is set by the model, not the profiler: ", focus_year, " has ",
         sum(is.finite(focus_all$SC_model)), " modelled vs ",
         sum(is.finite(focus_all$SC_profiler)), " profile half-hours. 2021 and 2022 give ",
         "R2 0.68 on the same comparison; 2023 is excluded (corrupt profiler trace).")
)

fig_caption <- paste(
  vapply(caption_paragraphs,
         function(p) paste(strwrap(p, width = 165), collapse = "\n"),
         character(1)),
  collapse = "\n"
)

final_fig <- p_nee / (p_sc | p_scatter) +
  plot_layout(heights = c(1, 1), widths = c(1)) +
  plot_annotation(
    caption = fig_caption,
    theme = theme(plot.caption = element_text(size = fig_base_size - 4, hjust = 0,
                                              lineheight = 1.15,
                                              margin = margin(t = 10)))
  )

ggsave(out_fig, plot = final_fig, width = 14, height = 11.6, dpi = 300, bg = "white")
cat("saved:", out_fig, "\n")

cat(sprintf(
  "\ndiel NEE curves differ by %.2f umol m-2 s-1 on average (max %.2f) over the 48 bins\n",
  diel_gap$mean_abs, diel_gap$max_abs))
