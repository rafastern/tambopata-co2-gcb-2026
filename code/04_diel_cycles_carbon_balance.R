# ───────────────────────────────────────────────────────────────────────────────
# ts + enso shading + nee_ok/reco/gep + photosynthetic capacity (pc) + surface conductance (Gs)
#
# update (this version):
#   - computes surface conductance (Gs, m s^-1) from Lee et al. (2021) using the enriched half-hourly file
#   - creates Gs time series with:
#       * ENSO-phase background shading
#       * wet/dry season point colors
#       * daily and 16-day window aggregation (matches the existing pattern)
#   - computes photosynthetic capacity (Pc) from GEP (no sign flip), day-only, met-constraint subset
#       * 16-day Pc (existing)
#       * daily Pc (Option A)
#   - keeps all existing outputs (Reco/GEP/NEE_ok, diel export, means table)
#   - NEW: friagem monthly diel cycles for NEE_ok, GEP, Reco (friagem vs no-friagem)
#
# inputs:
#   - dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv  (half-hourly, enriched)
#   - expects VPD_kPa already in kPa (no conversion)
#
# outputs added/kept:
#   - ts_daily_Gs_mps_points_ENSOshade.png
#   - ts_Gs_16day_mps_points_ENSOshade.png
#   - gs_16day_timeseries.csv
#   - ts_photosynthetic_capacity_pc_16day_ENSOshade.png
#   - ts_daily_photosynthetic_capacity_pc_ENSOshade.png
#   - photosynthetic_capacity_pc_16day_strict.csv
#   - photosynthetic_capacity_pc_16day_fallback.csv
#   - photosynthetic_capacity_pc_daily_strict.csv
#   - photosynthetic_capacity_pc_daily_fallback.csv
#   - FigS4_diel_by_year_season_ENSO_VPD_kPa.png
#   - Fig2_diel_monthly_background_plus_ENSO_means_<VAR>.png
#   - FigS4_diel_by_year_season_ENSO_<VAR>.png
#   - FigS5_diel_mean_by_ENSO_season_<VAR>.png
#   - FigX_diel_energy_partition_Rn_H_LE.png
#   - FigX_diel_FC_vs_NEE_mean_all_ENSO.png
#   - FigX_monthly_diurnal_sum_gC_m2_day_yearshape_keep_ENSO.png
#   - tambopata_48points_per_month.csv
#   - flux_means_16day_by_ENSO_season.csv
#   - NEW: diel_monthly_friagem_NEEok.png + .csv
#   - NEW: diel_monthly_friagem_GEP.png   + .csv
#   - NEW: diel_monthly_friagem_Reco.png  + .csv
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(grid)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# read enriched file as input (already contains ENSO/season/NEE_ok/Reco/GEP and VPD_kPa)
flux_fp <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")

out_csv_enriched    <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")

# pc outputs (kept + daily)
out_csv_pc_strict         <- file.path(output_path, "photosynthetic_capacity_pc_16day_strict.csv")
out_csv_pc_fallback       <- file.path(output_path, "photosynthetic_capacity_pc_16day_fallback.csv")
out_fig_pc_16day          <- file.path(fig_path,    "ts_photosynthetic_capacity_pc_16day_ENSOshade.png")

out_csv_pc_daily_strict   <- file.path(output_path, "photosynthetic_capacity_pc_daily_strict.csv")
out_csv_pc_daily_fallback <- file.path(output_path, "photosynthetic_capacity_pc_daily_fallback.csv")
out_fig_pc_daily          <- file.path(fig_path,    "ts_daily_photosynthetic_capacity_pc_ENSOshade.png")

# yearly diel output stubs (S4-style)
out_fig_s4_vpd <- file.path(fig_path, "FigS4_diel_by_year_season_ENSO_VPD_kPa.png")


# gs time series outputs
out_fig_gs_daily <- file.path(fig_path, "ts_daily_Gs_mps_points_ENSOshade.png")
out_fig_gs_16day <- file.path(fig_path, "ts_Gs_16day_mps_points_ENSOshade.png")
out_csv_gs_16day <- file.path(output_path, "gs_16day_timeseries.csv")

# NEW: friagem monthly diel outputs
out_fig_diel_monthly_friagem_nee  <- file.path(fig_path,    "diel_monthly_friagem_NEEok.png")
out_fig_diel_monthly_friagem_gep  <- file.path(fig_path,    "diel_monthly_friagem_GEP.png")
out_fig_diel_monthly_friagem_reco <- file.path(fig_path,    "diel_monthly_friagem_Reco.png")

out_csv_diel_monthly_friagem_nee  <- file.path(output_path, "diel_monthly_friagem_NEEok.csv")
out_csv_diel_monthly_friagem_gep  <- file.path(output_path, "diel_monthly_friagem_GEP.csv")
out_csv_diel_monthly_friagem_reco <- file.path(output_path, "diel_monthly_friagem_Reco.csv")

out_restored_diel_patterns <- c(
  "Fig2_diel_monthly_background_plus_ENSO_means_<VAR>.png",
  "FigS4_diel_by_year_season_ENSO_<VAR>.png",
  "FigS5_diel_mean_by_ENSO_season_<VAR>.png",
  "FigX_diel_energy_partition_Rn_H_LE.png",
  "FigX_diel_FC_vs_NEE_mean_all_ENSO.png",
  "FigX_monthly_diurnal_sum_gC_m2_day_yearshape_keep_ENSO.png"
)

# ───────────────────────────────────────────────────────────────────────────────
# config
tz_local <- "America/Lima"

dt_col     <- "tv_dt"
nee_col    <- "NEE"
nee_ok_col <- "NEE_ok"

par_col   <- "PAR"
tair_col  <- "TA_1_1_1"
vpd_col   <- "VPD_kPa"        # this is the kPa column in the enriched file

reco_col   <- "Reco"
gep_col    <- "GEP"
enso_col   <- "ENSO"
season_col <- "season"

dry_months <- 5:10

day_start_hour <- 6L
day_end_hour   <- 18L

tmax <- ymd_hms("2024-10-15 23:59:59", tz = tz_local)
xmin_plot <- ymd_hms("2018-01-01 00:00:00", tz = tz_local)
xmax_plot <- tmax

# season_cols / enso_cols now come from palette.R (sourced via paths.R):
#   season = blue (wet) / red (dry);
#   ENSO   = gold (El Nino) / purple (La Nina) / grey (Neutral)

# publication figure text settings
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

sec_per_halfhour <- 1800
gC_per_umolCO2   <- 12 / 1e6

# strict targets from excerpt (CI ignored) - VPD here is kPa
par_target  <- 836
par_halfwin <- 200
par_min <- par_target - par_halfwin
par_max <- par_target + par_halfwin

tair_target <- 27.22
tair_sd     <- 2.04
tair_min <- tair_target - tair_sd
tair_max <- tair_target + tair_sd

vpd_target <- 1.02
vpd_sd     <- 0.45
vpd_min <- vpd_target - vpd_sd
vpd_max <- vpd_target + vpd_sd

# coverage requirements (16-day; used for plot selection)
min_days_pc_strict   <- 8
min_points_pc_strict <- 200
min_days_pc_loose    <- 2
min_points_pc_loose  <- 30

# daily Pc coverage (so we don't plot days with a handful of half-hours)
min_points_pc_daily_strict <- 30
min_points_pc_daily_loose  <- 10

# selection thresholds (fallback)
fallback_use_if_total_points_lt <- 200
fallback_use_if_windows_lt       <- 10

# s4-style yearly diel config
min_total_points_year_s4 <- 200
min_points_month_diel <- 80

# x-axis ticks for half-hour binned diel plots
x_breaks <- c("02:00:00", "08:00:00", "14:00:00", "20:00:00")
x_labels <- c("02:00",    "08:00",    "14:00",    "20:00")

halfhour_levels <- sprintf(
  "%02d:%02d:00",
  rep(0:23, each = 2),
  rep(c(0, 30), times = 24)
)

# friagem column candidates
friagem_candidates <- c("friagem", "Friagem", "FRIAGEM", "friagem_event", "friagem_flag")

# ───────────────────────────────────────────────────────────────────────────────
# helpers
safe_stopifnot_cols <- function(d, cols) {
  miss <- cols[!cols %in% names(d)]
  if (length(miss)) stop(paste0("missing columns: ", paste(miss, collapse = ", ")))
}

msg_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("n/a")
  paste0("[", signif(min(x), 4), ", ", signif(max(x), 4), "]")
}

apply_coverage_16day <- function(pc_tbl, min_days, min_points) {
  pc_tbl %>% filter(n_days >= min_days, n_points >= min_points)
}

apply_coverage_daily <- function(pc_tbl, min_points) {
  pc_tbl %>% filter(n_points >= min_points)
}

# core selector: build the filtered point set used to define Pc
select_pc_points <- function(d, par_min, par_max, tair_min, tair_max, vpd_min, vpd_max) {
  d %>%
    mutate(gC_halfhour = GEP * sec_per_halfhour * gC_per_umolCO2) %>%
    filter(
      is_day_block,
      is.finite(gC_halfhour),
      is.finite(PAR),
      is.finite(Tair),
      is.finite(VPD),
      PAR  >= par_min  & PAR  <= par_max,
      Tair >= tair_min & Tair <= tair_max,
      VPD  >= vpd_min  & VPD  <= vpd_max
    )
}

# 16-day Pc from selected points
compute_pc_16day_from_points <- function(pts) {
  pts %>%
    mutate(
      window_start = as.Date(floor_date(Date, unit = "16 days")),
      window_mid   = window_start + days(8)
    ) %>%
    group_by(window_start, window_mid, season, ENSO) %>%
    summarise(
      Pc_gC_m2_day = mean(gC_halfhour, na.rm = TRUE) * 48,
      n_points     = dplyr::n(),
      n_days       = n_distinct(Date),
      .groups = "drop"
    ) %>%
    arrange(window_start)
}

# daily Pc from selected points
compute_pc_daily_from_points <- function(pts) {
  pts %>%
    group_by(Date, season, ENSO) %>%
    summarise(
      Pc_gC_m2_day = mean(gC_halfhour, na.rm = TRUE) * 48,
      n_points     = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      DateTime = as.POSIXct(Date, tz = tz_local),
      season = factor(as.character(season), levels = c("dry", "wet")),
      ENSO   = factor(as.character(ENSO),   levels = c("El Nino", "La Nina", "neutral"))
    ) %>%
    arrange(DateTime)
}

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx)) g$grobs[[idx[1]]] else NULL
}

legend_grob_panel <- function(lg) {
  if (is.null(lg)) {
    plot_spacer()
  } else {
    patchwork::wrap_elements(full = lg)
  }
}

row_strip <- function(lbl) {
  patchwork::wrap_elements(
    full = grid::textGrob(
      lbl,
      rot = 270,
      gp = grid::gpar(fontsize = 12)
    )
  )
}

make_yearly_lines_by_year_season_enso <- function(d, var, y_label, out_png, use_lines = TRUE) {
  if (!var %in% names(d)) stop(paste0("missing column: ", var))
  
  d2 <- d %>%
    filter(
      !is.na(.data[[var]]),
      !is.na(NEE_ok),
      !is.na(DateTime),
      !is.na(season),
      !is.na(ENSO),
      !is.na(year)
    ) %>%
    mutate(
      DateTime_30 = floor_date(DateTime, unit = "30 minutes"),
      hour_str    = format(DateTime_30, "%H:%M:%S"),
      hour        = factor(hour_str, levels = halfhour_levels),
      season      = factor(as.character(season), levels = c("dry", "wet")),
      ENSO        = factor(as.character(ENSO), levels = c("El Nino", "La Nina", "neutral")),
      year        = as.integer(year)
    )
  
  group_totals <- d2 %>%
    count(ENSO, season, year, name = "total_points") %>%
    filter(total_points >= min_total_points_year_s4)
  
  d2 <- d2 %>% inner_join(group_totals, by = c("ENSO", "season", "year"))
  
  diel_df <- d2 %>%
    group_by(ENSO, season, year, hour) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  if (!nrow(diel_df)) stop("no data after filters (min_total_points_year_s4 too strict?)")
  
  years_chr <- sort(unique(as.character(diel_df$year)))
  diel_df <- diel_df %>% mutate(year_chr = factor(as.character(year), levels = years_chr))
  
  custom_colors <- c("black","gray40","blue","red","purple","darkgreen","orange","brown","cyan","magenta")
  year_colors <- setNames(rep(custom_colors, length.out = length(years_chr)), years_chr)
  
  make_one_panel <- function(season_val, enso_val, show_legend = FALSE, show_y = TRUE, title_txt = NULL) {
    dd <- diel_df %>% filter(season == season_val, ENSO == enso_val)
    
    p <- ggplot(dd, aes(x = hour, y = mean_value, color = year_chr, group = year_chr)) +
      { if (use_lines) geom_line(linewidth = 0.8) } +
      geom_point(size = 1.4) +
      geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
      scale_x_discrete(breaks = x_breaks, labels = x_labels) +
      scale_color_manual(
        values = year_colors,
        breaks = years_chr,
        limits = years_chr,
        drop   = FALSE
      ) +
      theme_bw() +
      theme(
        axis.text.x  = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = if (show_y) element_text() else element_blank(),
        axis.text.y  = if (show_y) element_text() else element_blank(),
        axis.ticks.y = if (show_y) element_line() else element_blank(),
        legend.position = if (show_legend) "right" else "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      labs(x = "Hour of Day", y = y_label, color = NULL)
    
    if (!is.null(title_txt)) p <- p + labs(title = title_txt)
    p
  }
  
  p_dry_el  <- make_one_panel("dry", "El Nino",  show_legend = FALSE, show_y = TRUE,  title_txt = "El Niño")
  p_dry_la  <- make_one_panel("dry", "La Nina",  show_legend = FALSE, show_y = FALSE, title_txt = "La Niña")
  p_dry_neu <- make_one_panel("dry", "neutral",  show_legend = FALSE, show_y = FALSE, title_txt = "Neutral")
  
  p_wet_el  <- make_one_panel("wet", "El Nino",  show_legend = FALSE, show_y = TRUE)
  p_wet_la  <- make_one_panel("wet", "La Nina",  show_legend = FALSE, show_y = FALSE)
  
  legend_df <- data.frame(
    year_chr = factor(years_chr, levels = years_chr),
    x = 1, y = 1
  )
  
  p_leg_src <- ggplot(legend_df, aes(x = x, y = y, color = year_chr)) +
    geom_point(size = 3) +
    scale_color_manual(
      values = year_colors,
      breaks = years_chr,
      limits = years_chr,
      drop   = FALSE
    ) +
    theme_void() +
    theme(
      legend.position   = "right",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 11),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key        = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(title = NULL))
  
  lg <- extract_legend_grob(p_leg_src)
  
  p_wet_neu_as_legend <- legend_grob_panel(lg)
  
  top_row <- (p_dry_el | p_dry_la | p_dry_neu | row_strip("dry")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  bot_row <- (p_wet_el | p_wet_la | p_wet_neu_as_legend | row_strip("wet")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  
  full <- (top_row / bot_row) + plot_layout(heights = c(1, 1))
  
  ggsave(out_png, full, width = 12, height = 6, dpi = 300)
  invisible(full)
}

normalise_enso_for_diel_display <- function(x) {
  dplyr::recode(
    as.character(x),
    "El Nino" = "El Ni\u00f1o",
    "El Ni\u00f1o" = "El Ni\u00f1o",
    "La Nina" = "La Ni\u00f1a",
    "La Ni\u00f1a" = "La Ni\u00f1a",
    "neutral" = "Neutral",
    "Neutral" = "Neutral",
    .default = as.character(x)
  )
}

prepare_restored_diel_df <- function(d) {
  d %>%
    mutate(
      DateTime_30 = floor_date(DateTime, unit = "30 minutes"),
      hour = factor(format(DateTime_30, "%H:%M:%S"), levels = halfhour_levels),
      season = factor(tolower(trimws(as.character(season))), levels = c("dry", "wet")),
      ENSO_plot = factor(
        normalise_enso_for_diel_display(ENSO),
        levels = c("El Ni\u00f1o", "La Ni\u00f1a", "Neutral")
      ),
      year = as.integer(year),
      Rn = dplyr::coalesce(
        if ("Rn" %in% names(.)) as.numeric(Rn) else NA_real_,
        if ("NETRAD_1_1_1" %in% names(.)) as.numeric(NETRAD_1_1_1) else NA_real_
      )
    )
}

save_restored_diel_plot <- function(plot, file_stub, width, height) {
  out_png <- file.path(paper_figures_path, paste0(file_stub, ".png"))
  ggsave(out_png, plot = plot, width = width, height = height, dpi = 300)
  cat("saved:", out_png, "\n")
  invisible(out_png)
}

# `base_size`/`width`/`height` default to the standalone figure geometry, so every
# existing call is unchanged. They exist so the same panel can also be rendered
# small enough to drop into the two-column Figure 2 composite at 1:1 scale, where
# rescaling a 12-inch raster down to a 3.7-inch slot would leave the type at
# about 4.6 pt. See code/15_combine_figure2_env_stack_and_nee_diel.R.
# `save = FALSE` returns the panel without writing a PNG, so the supplement
# composites below can assemble panels re-rendered at a small base_size without
# also emitting a standalone file for each.
make_monthly_bg_plus_enso_plot <- function(d, var, var_pretty, y_label, file_stub,
                                           base_size = 16, width = 12, height = 6.5,
                                           legend_inside = TRUE, save = TRUE) {
  if (!var %in% names(d)) {
    message("skipping ", var, " (missing column)")
    return(invisible(NULL))
  }

  d2 <- d %>%
    filter(!is.na(.data[[var]]), !is.na(DateTime), !is.na(season), !is.na(hour)) %>%
    mutate(
      ym = floor_date(DateTime, unit = "month"),
      ym_lab = format(ym, "%Y-%m")
    )

  month_totals <- d2 %>%
    count(ENSO_plot, season, ym_lab, name = "n_pts") %>%
    filter(n_pts >= min_points_month_diel)

  d2 <- d2 %>% inner_join(month_totals, by = c("ENSO_plot", "season", "ym_lab"))

  if (!nrow(d2)) {
    message("skipping ", var, " (no monthly groups after min_points_month_diel filter)")
    return(invisible(NULL))
  }

  diel_month <- d2 %>%
    group_by(ENSO_plot, season, ym_lab, hour) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  diel_enso_mean <- d2 %>%
    group_by(ENSO_plot, season, hour) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  # month-to-month spread: SD across the individual monthly diel means -> +/-1 SD envelope
  diel_spread <- diel_month %>%
    group_by(ENSO_plot, season, hour) %>%
    summarise(sd_value = sd(mean_value, na.rm = TRUE), .groups = "drop")
  diel_enso_mean <- diel_enso_mean %>%
    left_join(diel_spread, by = c("ENSO_plot", "season", "hour")) %>%
    mutate(ymin = mean_value - sd_value, ymax = mean_value + sd_value)

  background_colors <- enso_cols
  background_linetypes <- enso_linetypes
  background_shapes <- c("El Ni\u00f1o" = 16, "La Ni\u00f1a" = 17)
  mean_shapes <- c("El Ni\u00f1o" = 16, "La Ni\u00f1a" = 17, "Neutral" = NA)

  p <- ggplot() +
    # +/-1 SD envelope of month-to-month variation per ENSO phase (replaces the
    # faint per-month "spaghetti" lines, which were nearly impossible to see)
    geom_ribbon(
      data = diel_enso_mean,
      aes(x = hour, ymin = ymin, ymax = ymax, group = ENSO_plot, fill = ENSO_plot),
      alpha = 0.18,
      color = NA,
      na.rm = TRUE
    ) +
    geom_line(
      data = diel_enso_mean,
      aes(x = hour, y = mean_value, color = ENSO_plot, linetype = ENSO_plot, group = ENSO_plot),
      linewidth = 1.5 * base_size / 16
    ) +
    geom_point(
      data = diel_enso_mean,
      aes(x = hour, y = mean_value, color = ENSO_plot, shape = ENSO_plot),
      size = 2.1 * base_size / 16,
      show.legend = FALSE
    ) +
    facet_wrap(~ season, ncol = 1) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = background_colors, breaks = c("El Ni\u00f1o", "La Ni\u00f1a", "Neutral")) +
    scale_linetype_manual(values = background_linetypes, breaks = c("El Ni\u00f1o", "La Ni\u00f1a", "Neutral")) +
    scale_shape_manual(values = mean_shapes, breaks = c("El Ni\u00f1o", "La Ni\u00f1a", "Neutral")) +
    scale_fill_manual(values = background_colors, breaks = c("El Ni\u00f1o", "La Ni\u00f1a", "Neutral"), guide = "none") +
    guides(color = guide_legend(order = 1), linetype = "none", shape = "none") +
    theme_bw(base_size = base_size) +
    theme(
      # sizes are ratios of base_size chosen to reproduce the original absolute
      # values (15/18/15) exactly at the default base_size of 16
      axis.text.x = element_text(size = base_size * 15 / 16, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = base_size * 15 / 16),
      axis.title.y = element_text(size = base_size * 18 / 16),
      axis.title.x = element_blank(),
      strip.text = element_text(size = base_size * 18 / 16),
      strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
      legend.position = if (legend_inside) c(0.98, 0.05) else "none",
      legend.justification = c(1, 0),
      legend.text = element_text(size = base_size * 15 / 16),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.3),
      # 8 pt at the default base_size of 16, so the standalone figures are byte-
      # identical; scales down with the panel version
      plot.margin = margin(rep(base_size / 2, 4))
    ) +
    labs(x = "Hour of Day", y = y_label, color = NULL)

  if (save) save_restored_diel_plot(p, file_stub, width = width, height = height)
  invisible(p)
}

make_yearly_lines_by_year_season_enso_restored <- function(d, var, var_pretty, y_label, file_stub, use_lines = TRUE) {
  if (!var %in% names(d)) {
    message("skipping ", var, " (missing column)")
    return(invisible(NULL))
  }

  d2 <- d %>%
    filter(!is.na(.data[[var]]), !is.na(DateTime), !is.na(season), !is.na(ENSO_plot), !is.na(year), !is.na(hour))

  group_totals <- d2 %>%
    count(ENSO_plot, season, year, name = "total_points") %>%
    filter(total_points >= min_total_points_year_s4)

  d2 <- d2 %>% inner_join(group_totals, by = c("ENSO_plot", "season", "year"))

  diel_df <- d2 %>%
    group_by(ENSO_plot, season, year, hour) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  if (!nrow(diel_df)) {
    message("skipping ", var, " (no data after min_total_points_year_s4 filter)")
    return(invisible(NULL))
  }

  years_chr <- sort(unique(as.character(diel_df$year)))
  diel_df <- diel_df %>% mutate(year_chr = factor(as.character(year), levels = years_chr))

  custom_colors <- c("black", "gray40", "blue", "red", "purple", "darkgreen", "orange", "brown", "cyan", "magenta")
  year_colors <- setNames(rep(custom_colors, length.out = length(years_chr)), years_chr)

  make_one_panel <- function(season_val, enso_val, show_y = TRUE) {
    dd <- diel_df %>% filter(season == season_val, ENSO_plot == enso_val)

    ggplot(dd, aes(x = hour, y = mean_value, color = year_chr, group = year_chr)) +
      { if (use_lines) geom_line(linewidth = 0.8) } +
      geom_point(size = 1.4) +
      geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
      scale_x_discrete(breaks = x_breaks, labels = x_labels) +
      scale_color_manual(values = year_colors, breaks = years_chr, limits = years_chr, drop = FALSE) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = if (show_y) element_text() else element_blank(),
        axis.text.y = if (show_y) element_text() else element_blank(),
        axis.ticks.y = if (show_y) element_line() else element_blank(),
        legend.position = "none",
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      labs(x = "Hour of Day", y = y_label, color = NULL)
  }

  p_dry_el <- make_one_panel("dry", "El Ni\u00f1o", show_y = TRUE) +
    labs(title = "El Ni\u00f1o") + theme(plot.title = element_text(hjust = 0.5))
  p_dry_la <- make_one_panel("dry", "La Ni\u00f1a", show_y = FALSE) +
    labs(title = "La Ni\u00f1a") + theme(plot.title = element_text(hjust = 0.5))
  p_dry_neu <- make_one_panel("dry", "Neutral", show_y = FALSE) +
    labs(title = "Neutral") + theme(plot.title = element_text(hjust = 0.5))
  p_wet_el <- make_one_panel("wet", "El Ni\u00f1o", show_y = TRUE)
  p_wet_la <- make_one_panel("wet", "La Ni\u00f1a", show_y = FALSE)

  legend_df <- data.frame(year_chr = factor(years_chr, levels = years_chr), x = 1, y = 1)
  p_leg_src <- ggplot(legend_df, aes(x = x, y = y, color = year_chr)) +
    geom_point(size = 3) +
    scale_color_manual(values = year_colors, breaks = years_chr, limits = years_chr, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(title = NULL))

  lg <- extract_legend_grob(p_leg_src)
  p_wet_neu_as_legend <- legend_grob_panel(lg)

  top_row <- (p_dry_el | p_dry_la | p_dry_neu | row_strip("dry")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  bot_row <- (p_wet_el | p_wet_la | p_wet_neu_as_legend | row_strip("wet")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  full <- (top_row / bot_row) + plot_layout(heights = c(1, 1))

  save_restored_diel_plot(full, file_stub, width = 12, height = 6)
  invisible(full)
}

# `base_size`/`width`/`height`/`legend_inside`/`save` default to the standalone
# geometry, so every existing call and its PNG are unchanged. They exist so this
# panel can also be rebuilt small for the supplement composite below; the sizes
# are expressed as ratios of base_size for exactly that reason (theme_bw()'s own
# default base_size is 11).
make_mean_lines_by_enso_with_seasons <- function(d, var, var_pretty, y_label, file_stub,
                                                 base_size = 11, width = 12, height = 5,
                                                 legend_inside = TRUE, save = TRUE) {
  if (!var %in% names(d)) {
    message("skipping ", var, " (missing column)")
    return(invisible(NULL))
  }

  d2 <- d %>%
    filter(!is.na(.data[[var]]), !is.na(DateTime), !is.na(season), !is.na(ENSO_plot), !is.na(year), !is.na(hour))

  group_totals <- d2 %>%
    count(ENSO_plot, season, year, name = "total_points") %>%
    filter(total_points >= min_total_points_year_s4)

  d2 <- d2 %>% inner_join(group_totals, by = c("ENSO_plot", "season", "year"))

  diel_mean <- d2 %>%
    group_by(ENSO_plot, season, hour) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  if (!nrow(diel_mean)) {
    message("skipping ", var, " (no data after min_total_points_year_s4 filter)")
    return(invisible(NULL))
  }

  p <- ggplot(diel_mean, aes(x = hour, y = mean_value, color = season, group = season)) +
    geom_line(linewidth = 1.0 * base_size / 11) +
    geom_point(size = 1.6 * base_size / 11) +
    facet_wrap(~ ENSO_plot, ncol = 3, drop = TRUE) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("dry" = "red", "wet" = "blue"), guide = guide_legend(title = NULL)) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = if (legend_inside) c(0.99, 0.02) else "right",
      legend.justification = if (legend_inside) c(1, 0) else "center",
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA),
      plot.margin = margin(5.5, 20, 5.5, 5.5) * base_size / 11
    ) +
    labs(x = "Hour of Day", y = y_label, color = NULL)

  if (save) save_restored_diel_plot(p, file_stub, width = width, height = height)
  invisible(p)
}

make_energy_partition_figure <- function(d, y_label, file_stub) {
  needed_vars <- c("Rn", "H", "LE")
  missing_vars <- needed_vars[!needed_vars %in% names(d)]
  if (length(missing_vars) > 0) {
    message("skipping energy partition figure; missing columns: ", paste(missing_vars, collapse = ", "))
    return(invisible(NULL))
  }

  d2 <- d %>%
    filter(!is.na(DateTime), !is.na(season), !is.na(ENSO_plot), !is.na(year), !is.na(hour))

  group_totals <- d2 %>%
    count(ENSO_plot, season, year, name = "total_points") %>%
    filter(total_points >= min_total_points_year_s4)

  d2 <- d2 %>% inner_join(group_totals, by = c("ENSO_plot", "season", "year"))

  diel_mean <- d2 %>%
    group_by(ENSO_plot, season, hour) %>%
    summarise(
      Rn = mean(Rn, na.rm = TRUE),
      H = mean(H, na.rm = TRUE),
      LE = mean(LE, na.rm = TRUE),
      .groups = "drop"
    )

  if (!nrow(diel_mean)) {
    message("skipping energy partition figure (no data after filtering)")
    return(invisible(NULL))
  }

  diel_long <- bind_rows(
    diel_mean %>% transmute(ENSO_plot, season, hour, variable = "Rn", mean_value = Rn),
    diel_mean %>% transmute(ENSO_plot, season, hour, variable = "H", mean_value = H),
    diel_mean %>% transmute(ENSO_plot, season, hour, variable = "LE", mean_value = LE)
  ) %>%
    mutate(variable = factor(variable, levels = c("Rn", "H", "LE")))

  energy_colors <- c("Rn" = "black", "H" = "red", "LE" = "blue")
  y_limits <- range(diel_long$mean_value, na.rm = TRUE)

  make_one_panel <- function(season_val, enso_val, show_y = TRUE) {
    dd <- diel_long %>% filter(season == season_val, ENSO_plot == enso_val)

    ggplot(dd, aes(x = hour, y = mean_value, color = variable, group = variable)) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.6) +
      geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
      scale_x_discrete(breaks = x_breaks, labels = x_labels) +
      scale_color_manual(values = energy_colors, drop = FALSE) +
      coord_cartesian(ylim = y_limits) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = if (show_y) element_text() else element_blank(),
        axis.text.y = if (show_y) element_text() else element_blank(),
        axis.ticks.y = if (show_y) element_line() else element_blank(),
        legend.position = "none",
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      labs(x = "Hour of Day", y = y_label, color = NULL)
  }

  p_dry_el <- make_one_panel("dry", "El Ni\u00f1o", show_y = TRUE) +
    labs(title = "El Ni\u00f1o") + theme(plot.title = element_text(hjust = 0.5))
  p_dry_la <- make_one_panel("dry", "La Ni\u00f1a", show_y = FALSE) +
    labs(title = "La Ni\u00f1a") + theme(plot.title = element_text(hjust = 0.5))
  p_dry_neu <- make_one_panel("dry", "Neutral", show_y = FALSE) +
    labs(title = "Neutral") + theme(plot.title = element_text(hjust = 0.5))
  p_wet_el <- make_one_panel("wet", "El Ni\u00f1o", show_y = TRUE)
  p_wet_la <- make_one_panel("wet", "La Ni\u00f1a", show_y = FALSE)

  legend_df <- data.frame(variable = factor(c("Rn", "H", "LE"), levels = c("Rn", "H", "LE")), x = 1, y = 1)
  p_leg_src <- ggplot(legend_df, aes(x = x, y = y, color = variable, group = variable)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.4) +
    scale_color_manual(values = energy_colors, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(title = NULL))

  lg <- extract_legend_grob(p_leg_src)
  p_wet_neu_as_legend <- legend_grob_panel(lg)

  top_row <- (p_dry_el | p_dry_la | p_dry_neu | row_strip("dry")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  bot_row <- (p_wet_el | p_wet_la | p_wet_neu_as_legend | row_strip("wet")) +
    plot_layout(widths = c(1, 1, 1, 0.06))
  full <- (top_row / bot_row) + plot_layout(heights = c(1, 1))

  save_restored_diel_plot(full, file_stub, width = 12, height = 6)
  invisible(full)
}

make_fc_nee_comparison_figure <- function(d, y_label, file_stub) {
  needed_vars <- c("FC", "NEE_ok")
  missing_vars <- needed_vars[!needed_vars %in% names(d)]
  if (length(missing_vars) > 0) {
    message("skipping FC vs NEE figure; missing columns: ", paste(missing_vars, collapse = ", "))
    return(invisible(NULL))
  }

  d2 <- d %>%
    filter(!is.na(DateTime), !is.na(season), !is.na(hour), !is.na(FC), !is.na(NEE_ok))

  diel_mean <- d2 %>%
    group_by(season, hour) %>%
    summarise(FC = mean(FC, na.rm = TRUE), NEE_ok = mean(NEE_ok, na.rm = TRUE), .groups = "drop")

  diel_long <- bind_rows(
    diel_mean %>% transmute(season, hour, variable = "NEE", mean_value = NEE_ok),
    diel_mean %>% transmute(season, hour, variable = "FC", mean_value = FC)
  ) %>%
    mutate(variable = factor(variable, levels = c("NEE", "FC")))

  p <- ggplot(diel_long, aes(x = hour, y = mean_value, color = variable, group = variable)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.6) +
    facet_wrap(~ season, ncol = 1) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("NEE" = "black", "FC" = "red"), guide = guide_legend(title = NULL)) +
    # 18, not ggplot's default 11: this is supplement Fig. S9, authored 12 in wide
    # and placed at 180 mm, so default type reached the page at 5-6 pt. Matches
    # the base_size used for Figs S3 and S5, which share the same 12 in width.
    theme_bw(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = c(0.98, 0.08),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Hour of Day", y = y_label, color = NULL)

  save_restored_diel_plot(p, file_stub, width = 12, height = 6)
  invisible(p)
}

run_restored_diel_figures <- function(d) {
  restored_df <- prepare_restored_diel_df(d)

  vars_to_plot <- list(
    list(var = "NEE_ok", pretty = "NEE", y = expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))),
    list(var = "FC", pretty = "FC", y = expression(FC ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))),
    list(var = "GEP", pretty = "GEP", y = expression(GEP ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))),
    list(var = "TA_1_1_1", pretty = "Ta", y = expression(italic(T)[plain(a)] ~ (degree*C))),
    list(var = "VPD_kPa", pretty = "VPD", y = expression(VPD ~ (kPa))),
    list(var = "SWC_1_1_1", pretty = "SWC", y = expression(SWC~"(cm)")),
    list(var = "PAR", pretty = "PAR", y = expression(PAR ~ (mu*mol ~ m^-2 ~ s^-1))),
    list(var = "LE", pretty = "LE", y = expression(LE ~ (W ~ m^-2))),
    list(var = "TS_3", pretty = "Ts", y = expression(italic(T)[plain(s)] ~ (degree*C))),
    list(var = "H", pretty = "H", y = expression(H ~ (W ~ m^-2))),
    list(var = "Rn", pretty = "Rn", y = expression(R[n] ~ (W ~ m^-2)))
  )

  for (v in vars_to_plot) {
    make_monthly_bg_plus_enso_plot(
      d = restored_df,
      var = v$var,
      var_pretty = v$pretty,
      y_label = v$y,
      file_stub = paste0("Fig2_diel_monthly_background_plus_ENSO_means_", v$var)
    )

    make_yearly_lines_by_year_season_enso_restored(
      d = restored_df,
      var = v$var,
      var_pretty = v$pretty,
      y_label = v$y,
      file_stub = paste0("FigS4_diel_by_year_season_ENSO_", v$var),
      use_lines = TRUE
    )

    make_mean_lines_by_enso_with_seasons(
      d = restored_df,
      var = v$var,
      var_pretty = v$pretty,
      y_label = v$y,
      file_stub = paste0("FigS5_diel_mean_by_ENSO_season_", v$var)
    )
  }

  # small-format copy of the NEE panel for the two-column Figure 2 composite.
  # Authored at its final size so code/15 can place it 1:1; the legend is off
  # because the ENSO key is carried by the environmental-stack panel.
  make_monthly_bg_plus_enso_plot(
    d = restored_df,
    var = "NEE_ok",
    var_pretty = "NEE",
    y_label = expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1)),
    file_stub = "Fig2a_diel_NEE_panel",
    base_size = 7, width = 3.70, height = 3.10, legend_inside = FALSE
  )

  # ─────────────────────────────────────────────────────────────────────────────
  # Supplement composite: the five meteorological diel cycles in one figure.
  #
  # These were supplementary figures S4-S8, one full page each. Reviewers asked
  # for a shorter supplement, so they are combined. The panels are REBUILT at
  # base_size 8 rather than the standalone rasters being rescaled: dropping a
  # 12-inch panel into a 3.5-inch slot would leave the type near 4.6 pt (the
  # failure documented in code/15). The standalone PNGs are still written above,
  # so nothing else that references them breaks.
  #
  # Legends: each panel places its key INSIDE the panel by default, which ggplot2
  # assembles into the panel gtable where patchwork's guide collection cannot see
  # it. legend_inside = FALSE moves it out so `guides = "collect"` works, and the
  # single collected key lands in the empty sixth slot via guide_area().
  met_panel_specs <- list(
    list(var = "PAR",       pretty = "PAR", y = expression(PAR ~ (mu*mol ~ m^-2 ~ s^-1))),
    list(var = "SWC_1_1_1", pretty = "SWC", y = expression(SWC~"(cm)")),
    list(var = "TA_1_1_1",  pretty = "Ta",  y = expression(italic(T)[plain(a)] ~ (degree*C))),
    list(var = "TS_3",      pretty = "Ts",  y = expression(italic(T)[plain(s)] ~ (degree*C))),
    list(var = "VPD_kPa",   pretty = "VPD", y = expression(VPD ~ (kPa)))
  )

  met_panels <- lapply(met_panel_specs, function(v) {
    make_monthly_bg_plus_enso_plot(
      d = restored_df, var = v$var, var_pretty = v$pretty, y_label = v$y,
      file_stub = paste0("unused_", v$var),
      base_size = 8, legend_inside = FALSE, save = FALSE
    )
  })
  met_panels <- Filter(Negate(is.null), met_panels)

  # The builder blanks axis.title.x (see the theme above), which is fine for a
  # standalone diel figure where the tick labels are obviously clock times. In a
  # 3x2 composite it is worth naming the axis once per column, so restore it on
  # the bottom panel of each: with the legend occupying slot 6, those are (e)
  # bottom-left and (d) bottom-right.
  if (length(met_panels) == length(met_panel_specs)) {
    for (i in c(4L, 5L)) {
      met_panels[[i]] <- met_panels[[i]] +
        theme(axis.title.x = element_text(size = 8 * 18 / 16))
    }
  }

  if (length(met_panels) == length(met_panel_specs)) {
    met_composite <- patchwork::wrap_plots(
      c(met_panels, list(patchwork::guide_area())), ncol = 2
    ) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
      theme(legend.position = "right", legend.direction = "vertical",
            plot.tag = element_text(size = 10, face = "bold"))

    out_met <- file.path(paper_figures_path, "FigS_diel_met_composite.png")
    ggsave(out_met, met_composite, width = 7.10, height = 8.60, dpi = 300, bg = "white")
    cat("saved:", out_met, "(replaces the five standalone met diel figures)\n")
  } else {
    message("met composite skipped: only ", length(met_panels), " of ",
            length(met_panel_specs), " panels available")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # Supplement composite: diel carbon fluxes.
  #
  # (a) mean diel NEE by ENSO phase (was S20) over (b) the diel GEP cycle (was
  # S21a). The per-year diel NEE figure (S19) is deliberately NOT included: at
  # 180 mm all three together would leave each sub-panel row about 0.65 in tall,
  # and S19 was never cited on its own. The full per-year set remains in the
  # Zenodo deposit.
  #
  # The two legends are different variables (season for (a), ENSO for (b)), so
  # `guides = "collect"` de-duplicates within each but keeps both keys.
  nee_mean_panel <- make_mean_lines_by_enso_with_seasons(
    d = restored_df, var = "NEE_ok", var_pretty = "NEE",
    y_label = expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1)),
    file_stub = "unused_composite_NEE", base_size = 8,
    legend_inside = FALSE, save = FALSE
  )
  gep_diel_panel <- make_monthly_bg_plus_enso_plot(
    d = restored_df, var = "GEP", var_pretty = "GEP",
    y_label = expression(GEP ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1)),
    file_stub = "unused_composite_GEP", base_size = 8,
    legend_inside = FALSE, save = FALSE
  )

  if (!is.null(nee_mean_panel) && !is.null(gep_diel_panel)) {
    gep_diel_panel <- gep_diel_panel +
      theme(axis.title.x = element_text(size = 8 * 18 / 16))

    carbon_composite <- (nee_mean_panel / gep_diel_panel) +
      patchwork::plot_layout(heights = c(1, 1.7), guides = "collect") +
      patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
      theme(legend.position = "right", plot.tag = element_text(size = 10, face = "bold"))

    out_carbon <- file.path(paper_figures_path, "FigS_diel_carbon_composite.png")
    ggsave(out_carbon, carbon_composite, width = 7.10, height = 7.00, dpi = 300, bg = "white")
    cat("saved:", out_carbon, "(replaces the mean diel NEE and diel GEP figures)\n")
  } else {
    message("diel carbon composite skipped: a panel was unavailable")
  }

  make_energy_partition_figure(
    d = restored_df,
    y_label = expression(Flux ~ (W ~ m^-2)),
    file_stub = "FigX_diel_energy_partition_Rn_H_LE"
  )

  make_fc_nee_comparison_figure(
    d = restored_df,
    y_label = expression(Flux ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1)),
    file_stub = "FigX_diel_FC_vs_NEE_mean_all_ENSO"
  )

  invisible(restored_df)
}

round_smart <- function(x, digits = 2) {
  ifelse(is.na(x), NA_real_, round(x, digits))
}

# monthly diel by friagem (48 half-hours)
make_monthly_diel_by_friagem <- function(d, var, ylab, out_png, out_csv) {
  if (!var %in% names(d)) stop(paste0("missing column: ", var))
  if (!"friagem" %in% names(d)) stop("missing column: friagem")
  
  dd <- d %>%
    filter(!is.na(DateTime), !is.na(friagem), is.finite(.data[[var]])) %>%
    mutate(
      DateTime_30 = lubridate::floor_date(DateTime, unit = "30 minutes"),
      hour_str    = format(DateTime_30, "%H:%M:%S"),
      hour        = factor(hour_str, levels = halfhour_levels),
      month_i     = lubridate::month(DateTime),
      month_lbl   = factor(month_i, levels = 1:12, labels = month.abb),
      friagem_lbl = factor(if_else(friagem, "friagem", "no-friagem"), levels = c("no-friagem", "friagem"))
    )
  
  diel <- dd %>%
    group_by(month_lbl, friagem_lbl, hour) %>%
    summarise(
      mean_value = mean(.data[[var]], na.rm = TRUE),
      n_obs      = dplyr::n(),
      .groups = "drop"
    )
  
  if (!nrow(diel)) stop("no data for monthly diel by friagem (after filters).")
  
  write_csv(diel, out_csv)
  cat("saved:", out_csv, "\n")
  
  p <- ggplot(diel, aes(x = hour, y = mean_value, group = friagem_lbl, linetype = friagem_lbl)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.2) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    facet_wrap(~ month_lbl, ncol = 4) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = c(0.09, 0.79),   # top-right inside plotting area
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = alpha("white", 0.85), color = "grey60"),
      legend.key = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70")
    ) +
    labs(y = ylab, linetype = NULL)
  
  ggsave(out_png, p, width = 14, height = 8, dpi = 300)
  cat("saved:", out_png, "\n")
  
  invisible(p)
}

# ───────────────────────────────────────────────────────────────────────────────
# gs helpers (Lee et al. 2021 style; returns m s^-1)
esat_kPa <- function(Tair_C) {
  0.6108 * exp(17.27 * Tair_C / (Tair_C + 237.3))
}

slope_svp <- function(Tair_C) {
  es <- esat_kPa(Tair_C)
  4098 * es / (Tair_C + 237.3)^2
}

Gaero_from_u <- function(u_star, u) {
  u_safe <- ifelse(is.finite(u) & u > 0, u, NA_real_)
  
  ge <- (u_star^2) / u_safe
  re <- 1 / ge
  rb <- 6.2 * u_star^(-2/3)
  
  inv_Gaero <- re + rb
  Gaero <- 1 / inv_Gaero
  Gaero[!is.finite(Gaero)] <- NA_real_
  Gaero
}

Gs_from_fluxes <- function(LE, H, VPD, Tair_C, u_star, u,
                           rho_air = 1.2,
                           cp_air  = 1007,
                           gamma   = 0.066) {
  beta  <- H / LE
  s     <- slope_svp(Tair_C)
  Gaero <- Gaero_from_u(u_star, u)
  
  term1 <- (rho_air * cp_air * VPD) / (gamma * LE)
  term2 <- (s * beta - gamma) / (gamma * Gaero)
  
  inv_Gs <- term1 + term2
  Gs <- 1 / inv_Gs
  Gs[!is.finite(Gs)] <- NA_real_
  Gs
}

# pick the first existing name from a candidate list
find_col <- function(candidates, names_vec) {
  hits <- candidates[candidates %in% names_vec]
  if (length(hits) == 0) NA_character_ else hits[1]
}

# ───────────────────────────────────────────────────────────────────────────────
# load + parse (from enriched file)
stopifnot(file.exists(flux_fp))
df_raw <- read_csv(flux_fp, show_col_types = FALSE)

safe_stopifnot_cols(
  df_raw,
  c(dt_col, nee_col, nee_ok_col, par_col, tair_col, vpd_col, reco_col, gep_col, enso_col, season_col)
)

df <- df_raw %>%
  mutate(
    DateTime = dmy_hms(.data[[dt_col]], tz = tz_local),
    NEE      = as.numeric(.data[[nee_col]]),
    NEE_ok   = as.numeric(.data[[nee_ok_col]]),
    PAR      = as.numeric(.data[[par_col]]),
    PAR_corrected_SWin = if ("PAR_corrected_SWin" %in% names(df_raw)) {
      as.numeric(.data[["PAR_corrected_SWin"]])
    } else {
      if ("PPFD_IN_1_1_1" %in% names(df_raw)) as.numeric(.data[["PPFD_IN_1_1_1"]]) else PAR
    },
    PAR_source = if ("PAR_source" %in% names(df_raw)) {
      as.character(.data[["PAR_source"]])
    } else {
      "matlab_second_stage_swin_corrected"
    },
    SW_IN = if ("SW_IN" %in% names(df_raw)) as.numeric(.data[["SW_IN"]]) else {
      if ("SWin" %in% names(df_raw)) as.numeric(.data[["SWin"]]) else NA_real_
    },
    SWin = SW_IN,
    Tair     = as.numeric(.data[[tair_col]]),
    VPD      = as.numeric(.data[[vpd_col]]),
    Reco     = as.numeric(.data[[reco_col]]),
    GEP      = as.numeric(.data[[gep_col]]),
    ENSO     = as.character(.data[[enso_col]]),
    season   = as.character(.data[[season_col]])
  ) %>%
  filter(!is.na(DateTime)) %>%
  arrange(DateTime) %>%
  mutate(
    year   = year(DateTime),
    month  = month(DateTime),
    Date   = as.Date(DateTime, tz = tz_local),
    hour   = hour(DateTime),
    is_day_block = hour >= day_start_hour & hour < day_end_hour,
    season = if_else(season %in% c("dry", "wet"), season, if_else(month %in% dry_months, "dry", "wet")),
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = if_else(is.na(ENSO) | ENSO == "", "neutral", ENSO)
  )

par_mismatch <- is.finite(df$PAR) & is.finite(df$PAR_corrected_SWin) &
  abs(df$PAR - df$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin in enriched input")
}

# positive GEP is retained, matching code/01_prepare_flux_timeseries.R. clipping
# it to zero censors one side of the Reco error and breaks NEE = Reco + GEP; see
# the note at the GEP calculation in 01 for the morning-transition mechanism.

# ───────────────────────────────────────────────────────────────────────────────
# friagem parsing (logical)
friagem_col <- find_col(friagem_candidates, names(df_raw))

if (is.na(friagem_col)) {
  cat("warning: no friagem column found. skipping friagem diel plots.\n")
  df <- df %>% mutate(friagem = NA)
} else {
  cat("friagem column mapping:\n")
  print(c(friagem = friagem_col))
  
  df <- df %>%
    mutate(
      friagem_raw = .data[[friagem_col]],
      friagem = dplyr::case_when(
        is.logical(friagem_raw) ~ friagem_raw,
        is.numeric(friagem_raw) ~ (friagem_raw > 0),
        is.character(friagem_raw) ~ tolower(friagem_raw) %in% c("true", "t", "1", "yes", "y"),
        is.factor(friagem_raw) ~ tolower(as.character(friagem_raw)) %in% c("true", "t", "1", "yes", "y"),
        TRUE ~ NA
      )
    )
  
  cat("friagem diagnostics:\n")
  cat("non-na:", sum(!is.na(df$friagem)), "\n")
  cat("true :", sum(df$friagem %in% TRUE, na.rm = TRUE), "\n")
  cat("false:", sum(df$friagem %in% FALSE, na.rm = TRUE), "\n\n")
}

cat("\ninput diagnostics:\n")
cat("read:", flux_fp, "\n")
cat("rows:", nrow(df), "\n")
cat("datetime range:", as.character(min(df$DateTime)), "to", as.character(max(df$DateTime)), "\n")
cat("vpd (kPa) range:", msg_range(df$VPD), "\n")
cat("par range:", msg_range(df$PAR), "\n")
cat("tair range:", msg_range(df$Tair), "\n")
cat("finite gep daytime:", sum(df$is_day_block & is.finite(df$GEP)), "\n\n")

# ───────────────────────────────────────────────────────────────────────────────
# limit time range
df <- df %>% filter(DateTime <= tmax)



# ───────────────────────────────────────────────────────────────────────────────
# monthly GPP climatology: average daily GPP for each calendar month

out_csv_gpp_monthly_year <- file.path(
  output_path,
  "gpp_daily_mean_by_year_month.csv"
)

out_csv_gpp_monthly_climatology <- file.path(
  output_path,
  "gpp_daily_mean_climatology_by_calendar_month.csv"
)

min_hh_gep_per_day <- 18  # require at least 75% of daytime half-hours, 18 of 24

daily_gpp <- df %>%
  filter(
    is_day_block,
    !is.na(Date),
    is.finite(GEP)
  ) %>%
  group_by(Date) %>%
  summarise(
    n_hh_gep = dplyr::n(),
    GEP_gC_m2_day = sum(GEP * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_hh_gep >= min_hh_gep_per_day) %>%
  mutate(
    GPP_gC_m2_day = -GEP_gC_m2_day,
    year = lubridate::year(Date),
    month = lubridate::month(Date),
    month_name = factor(month.abb[month], levels = month.abb)
  ) %>%
  filter(is.finite(GPP_gC_m2_day))

# a handful of complete days integrate to a negative daily GPP once positive
# half-hourly GEP is no longer clipped. dropping them would reimpose the same
# one-sided censoring at daily scale, so they are kept and reported instead.
neg_gpp_days <- daily_gpp %>% filter(GPP_gC_m2_day < 0) %>% arrange(GPP_gC_m2_day)
if (nrow(neg_gpp_days) > 0) {
  cat("\ndays with negative daily GPP (retained, not dropped):", nrow(neg_gpp_days),
      "of", nrow(daily_gpp), "\n")
  print(neg_gpp_days %>% dplyr::select(Date, n_hh_gep, GPP_gC_m2_day))
}

gpp_by_year_month <- daily_gpp %>%
  group_by(year, month, month_name) %>%
  summarise(
    GPP_mean_gC_m2_day = mean(GPP_gC_m2_day, na.rm = TRUE),
    GPP_sd_gC_m2_day   = sd(GPP_gC_m2_day, na.rm = TRUE),
    n_days             = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(year, month)

gpp_monthly_climatology <- gpp_by_year_month %>%
  group_by(month, month_name) %>%
  summarise(
    GPP_mean_climatology_gC_m2_day = mean(GPP_mean_gC_m2_day, na.rm = TRUE),
    GPP_sd_across_year_months_gC_m2_day = sd(GPP_mean_gC_m2_day, na.rm = TRUE),
    n_year_months = dplyr::n(),
    mean_valid_days_per_year_month = mean(n_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month)

write_csv(gpp_by_year_month, out_csv_gpp_monthly_year)
cat("saved:", out_csv_gpp_monthly_year, "\n")

write_csv(gpp_monthly_climatology, out_csv_gpp_monthly_climatology)
cat("saved:", out_csv_gpp_monthly_climatology, "\n")

cat("\nmonthly GPP climatology:\n")
print(gpp_monthly_climatology, n = 12)

# ───────────────────────────────────────────────────────────────────────────────
# enso shading rectangles across contiguous periods
enso_rect <- df %>%
  transmute(DateTime, ENSO = as.character(ENSO)) %>%
  arrange(DateTime) %>%
  mutate(
    change = (ENSO != dplyr::lag(ENSO, default = dplyr::first(ENSO))),
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
# compute Gs (m s^-1) at half-hourly resolution using columns available in df_raw
nms <- names(df_raw)

le_candidates    <- c("LE", "LE_f", "LE_fqcOK", "LE_qcOK", "LE_CORR", "LE_orig")
h_candidates     <- c("H", "H_f", "H_fqcOK", "H_qcOK", "H_CORR", "H_orig")
ustar_candidates <- c("USTAR", "u_star", "Ustar")
wind_candidates  <- c("WIND", "WS_1_1_1", "WS", "U", "WS_ms")

col_LE    <- find_col(le_candidates,    nms)
col_H     <- find_col(h_candidates,     nms)
col_USTAR <- find_col(ustar_candidates, nms)
col_WIND  <- find_col(wind_candidates,  nms)

cat("gs column mapping:\n")
print(c(LE = col_LE, H = col_H, USTAR = col_USTAR, WIND = col_WIND))

if (any(is.na(c(col_LE, col_H, col_USTAR, col_WIND)))) {
  cat("\nwarning: cannot compute Gs because at least one of LE/H/USTAR/WIND is missing in the enriched file.\n")
  cat("warning: skipping Gs time series.\n\n")
  df <- df %>% mutate(Gs_mps = NA_real_)
} else {
  df <- df %>%
    mutate(
      LE_gs    = as.numeric(.data[[col_LE]]),
      H_gs     = as.numeric(.data[[col_H]]),
      USTAR_gs = as.numeric(.data[[col_USTAR]]),
      WIND_gs  = as.numeric(.data[[col_WIND]])
    ) %>%
    mutate(
      LE_gs    = if_else(is.finite(LE_gs)    & LE_gs > 0,      LE_gs,    NA_real_),
      H_gs     = if_else(is.finite(H_gs),                    H_gs,     NA_real_),
      USTAR_gs = if_else(is.finite(USTAR_gs) & USTAR_gs > 0,  USTAR_gs, NA_real_),
      WIND_gs  = if_else(is.finite(WIND_gs)  & WIND_gs > 0,   WIND_gs,  NA_real_),
      VPD      = if_else(is.finite(VPD)      & VPD >= 0,      VPD,      NA_real_),
      Tair     = if_else(is.finite(Tair),                   Tair,     NA_real_),
      Gs_mps = Gs_from_fluxes(
        LE     = LE_gs,
        H      = H_gs,
        VPD    = VPD,
        Tair_C = Tair,
        u_star = USTAR_gs,
        u      = WIND_gs
      )
    )
  
  cat("\nGs diagnostics:\n")
  cat("finite gs:", sum(is.finite(df$Gs_mps)), "\n")
  cat("gs range:", msg_range(df$Gs_mps), "\n\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# pc: strict (points -> daily + 16-day)
pc_points_strict <- select_pc_points(df, par_min, par_max, tair_min, tair_max, vpd_min, vpd_max)
pc_16day_strict  <- compute_pc_16day_from_points(pc_points_strict)
pc_daily_strict  <- compute_pc_daily_from_points(pc_points_strict)

cat("pc strict diagnostics:\n")
cat("strict points:", nrow(pc_points_strict), "\n")
cat("strict 16-day windows:", nrow(pc_16day_strict), "\n")
cat("strict daily rows:", nrow(pc_daily_strict), "\n\n")

# ───────────────────────────────────────────────────────────────────────────────
# fallback pc (dataset-centered mean ± 1 sd), if strict too sparse
fallback_used <- FALSE
fallback_box <- NULL

pc_points_fallback <- tibble()
pc_16day_fallback  <- tibble()
pc_daily_fallback  <- tibble()

if (nrow(pc_points_strict) < fallback_use_if_total_points_lt || nrow(pc_16day_strict) < fallback_use_if_windows_lt) {
  fallback_used <- TRUE
  
  base <- df %>%
    filter(is_day_block, is.finite(GEP), is.finite(PAR), is.finite(Tair), is.finite(VPD))
  
  if (nrow(base) > 0) {
    m_par  <- mean(base$PAR,  na.rm = TRUE)
    sd_par <- sd(base$PAR,    na.rm = TRUE)
    
    m_t    <- mean(base$Tair, na.rm = TRUE)
    sd_t   <- sd(base$Tair,   na.rm = TRUE)
    
    m_v    <- mean(base$VPD,  na.rm = TRUE)
    sd_v   <- sd(base$VPD,    na.rm = TRUE)
    
    par_min_fb  <- max(0, m_par - sd_par)
    par_max_fb  <- m_par + sd_par
    
    tair_min_fb <- m_t - sd_t
    tair_max_fb <- m_t + sd_t
    
    vpd_min_fb  <- max(0, m_v - sd_v)
    vpd_max_fb  <- m_v + sd_v
    
    fallback_box <- list(
      par_min = par_min_fb, par_max = par_max_fb,
      tair_min = tair_min_fb, tair_max = tair_max_fb,
      vpd_min = vpd_min_fb, vpd_max = vpd_max_fb
    )
    
    cat("\nfallback pc box (dataset-centered, mean ± 1 sd):\n")
    cat("par :", signif(par_min_fb,4), "to", signif(par_max_fb,4), "\n")
    cat("tair:", signif(tair_min_fb,4), "to", signif(tair_max_fb,4), "\n")
    cat("vpd :", signif(vpd_min_fb,4), "to", signif(vpd_max_fb,4), " (kPa)\n\n")
    
    pc_points_fallback <- select_pc_points(df, par_min_fb, par_max_fb, tair_min_fb, tair_max_fb, vpd_min_fb, vpd_max_fb)
    pc_16day_fallback  <- compute_pc_16day_from_points(pc_points_fallback)
    pc_daily_fallback  <- compute_pc_daily_from_points(pc_points_fallback)
    
    cat("fallback points:", nrow(pc_points_fallback), "\n")
    cat("fallback 16-day windows:", nrow(pc_16day_fallback), "\n")
    cat("fallback daily rows:", nrow(pc_daily_fallback), "\n\n")
  } else {
    cat("\nwarning: no rows with finite gep + met vars; cannot build fallback pc.\n")
  }
}

# ───────────────────────────────────────────────────────────────────────────────
# choose which Pc to use for plotting (16-day) and (daily)
use_fallback <- fallback_used && nrow(pc_16day_fallback) > 0

pc_16day_plot_df <- tibble()
pc_daily_plot_df <- tibble()

pc_plot_note_16day <- ""
pc_plot_note_daily <- ""

if (use_fallback) {
  a16 <- apply_coverage_16day(pc_16day_fallback, min_days_pc_strict, min_points_pc_strict)
  b16 <- apply_coverage_16day(pc_16day_fallback, min_days_pc_loose,  min_points_pc_loose)
  
  if (nrow(a16) > 0) {
    pc_16day_plot_df <- a16
    pc_plot_note_16day <- paste0("pc 16-day (fallback) + coverage: n_days>=", min_days_pc_strict, ", n_points>=", min_points_pc_strict)
  } else if (nrow(b16) > 0) {
    pc_16day_plot_df <- b16
    pc_plot_note_16day <- paste0("pc 16-day (fallback) + loose coverage: n_days>=", min_days_pc_loose, ", n_points>=", min_points_pc_loose)
  } else {
    pc_16day_plot_df <- pc_16day_fallback
    pc_plot_note_16day <- "pc 16-day (fallback; no coverage filter passed)"
  }
  
  ad <- apply_coverage_daily(pc_daily_fallback, min_points_pc_daily_strict)
  bd <- apply_coverage_daily(pc_daily_fallback, min_points_pc_daily_loose)
  
  if (nrow(ad) > 0) {
    pc_daily_plot_df <- ad
    pc_plot_note_daily <- paste0("pc daily (fallback) + coverage: n_points>=", min_points_pc_daily_strict)
  } else if (nrow(bd) > 0) {
    pc_daily_plot_df <- bd
    pc_plot_note_daily <- paste0("pc daily (fallback) + loose coverage: n_points>=", min_points_pc_daily_loose)
  } else {
    pc_daily_plot_df <- pc_daily_fallback
    pc_plot_note_daily <- "pc daily (fallback; no coverage filter passed)"
  }
} else if (nrow(pc_16day_strict) > 0) {
  a16 <- apply_coverage_16day(pc_16day_strict, min_days_pc_strict, min_points_pc_strict)
  b16 <- apply_coverage_16day(pc_16day_strict, min_days_pc_loose,  min_points_pc_loose)
  
  if (nrow(a16) > 0) {
    pc_16day_plot_df <- a16
    pc_plot_note_16day <- paste0("pc 16-day (strict) + coverage: n_days>=", min_days_pc_strict, ", n_points>=", min_points_pc_strict)
  } else if (nrow(b16) > 0) {
    pc_16day_plot_df <- b16
    pc_plot_note_16day <- paste0("pc 16-day (strict) + loose coverage: n_days>=", min_days_pc_loose, ", n_points>=", min_points_pc_loose)
  } else {
    pc_16day_plot_df <- pc_16day_strict
    pc_plot_note_16day <- "pc 16-day (strict; no coverage filter passed)"
  }
  
  ad <- apply_coverage_daily(pc_daily_strict, min_points_pc_daily_strict)
  bd <- apply_coverage_daily(pc_daily_strict, min_points_pc_daily_loose)
  
  if (nrow(ad) > 0) {
    pc_daily_plot_df <- ad
    pc_plot_note_daily <- paste0("pc daily (strict) + coverage: n_points>=", min_points_pc_daily_strict)
  } else if (nrow(bd) > 0) {
    pc_daily_plot_df <- bd
    pc_plot_note_daily <- paste0("pc daily (strict) + loose coverage: n_points>=", min_points_pc_daily_loose)
  } else {
    pc_daily_plot_df <- pc_daily_strict
    pc_plot_note_daily <- "pc daily (strict; no coverage filter passed)"
  }
} else {
  pc_plot_note_16day <- "no pc data available (strict+fallback empty)"
  pc_plot_note_daily <- "no pc data available (strict+fallback empty)"
  cat("error: no pc data to plot.\n")
}

cat("pc plotting selection:\n")
cat(pc_plot_note_16day, "\n")
cat("pc_16day_plot_df windows:", nrow(pc_16day_plot_df), "\n")
cat(pc_plot_note_daily, "\n")
cat("pc_daily_plot_df days:", nrow(pc_daily_plot_df), "\n\n")

# ───────────────────────────────────────────────────────────────────────────────
# NEW: diel cycles: monthly, split by friagem vs no-friagem
if (any(!is.na(df$friagem))) {
  make_monthly_diel_by_friagem(
    d = df,
    var = "NEE_ok",
    ylab = expression(NEE ~ "(" * mu * "mol" ~ m^-2 ~ s^-1 * ")"),
    out_png = out_fig_diel_monthly_friagem_nee,
    out_csv = out_csv_diel_monthly_friagem_nee
  )
  
  make_monthly_diel_by_friagem(
    d = df %>% filter(is_day_block),
    var = "GEP",
    ylab = expression(GEP ~ "(" * mu * "mol" ~ m^-2 ~ s^-1 * ")"),
    out_png = out_fig_diel_monthly_friagem_gep,
    out_csv = out_csv_diel_monthly_friagem_gep
  )
  
  make_monthly_diel_by_friagem(
    d = df,
    var = "Reco",
    ylab = expression(R[eco] ~ "(" * mu * "mol" ~ m^-2 ~ s^-1 * ")"),
    out_png = out_fig_diel_monthly_friagem_reco,
    out_csv = out_csv_diel_monthly_friagem_reco
  )
} else {
  cat("skipped friagem diel plots (friagem column missing or all NA).\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# daily Gs time series (m s^-1)
if (any(is.finite(df$Gs_mps))) {
  gs_daily <- df %>%
    filter(!is.na(Date), is.finite(Gs_mps), Gs_mps >= 0) %>%
    group_by(Date) %>%
    summarise(
      Gs_mps = median(Gs_mps, na.rm = TRUE),
      season   = dplyr::first(na.omit(as.character(season))),
      ENSO     = dplyr::first(na.omit(as.character(ENSO))),
      n_hh_obs = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      DateTime = as.POSIXct(Date, tz = tz_local),
      season = factor(season, levels = c("dry", "wet")),
      ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
    ) %>%
    arrange(DateTime)
  
  ylab_gs <- expression(italic(G)[plain(s)] ~ "(" * m ~ s^-1 * ")")
  
  p_gs_daily <- ggplot() +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
    scale_fill_manual(
      values = enso_fill,
      drop = FALSE,
      name = NULL,
      labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
    ) +
    geom_point(
      data = gs_daily,
      aes(x = DateTime, y = Gs_mps, color = season),
      size = 1.6,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = season_cols,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = ylab_gs) +
    theme_pub(base_size = 16, legend_position = "right") +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))
  
  ggsave(out_fig_gs_daily, p_gs_daily, width = 14, height = 4.2, dpi = 300)
  cat("saved:", out_fig_gs_daily, "\n")
} else {
  cat("skipped daily gs plot (no finite gs).\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# 16-day window Gs time series (m s^-1)
if (any(is.finite(df$Gs_mps))) {
  gs_16day <- df %>%
    filter(!is.na(Date), is.finite(Gs_mps), Gs_mps >= 0) %>%
    mutate(
      window_start = as.Date(floor_date(Date, unit = "16 days")),
      window_mid   = window_start + days(8)
    ) %>%
    group_by(window_start, window_mid) %>%
    summarise(
      Gs_mps = median(Gs_mps, na.rm = TRUE),
      season   = dplyr::first(na.omit(as.character(season))),
      ENSO     = dplyr::first(na.omit(as.character(ENSO))),
      n_hh_obs = dplyr::n(),
      n_days   = n_distinct(Date),
      .groups = "drop"
    ) %>%
    mutate(
      DateTime = as.POSIXct(window_mid, tz = tz_local),
      season = factor(season, levels = c("dry", "wet")),
      ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
    ) %>%
    arrange(DateTime)
  
  write_csv(gs_16day, out_csv_gs_16day)
  cat("saved:", out_csv_gs_16day, "\n")
  
  ylab_gs <- expression(italic(G)[plain(s)] ~ "(" * m ~ s^-1 * ")")
  
  p_gs_16day <- ggplot() +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
    scale_fill_manual(
      values = enso_fill,
      drop = FALSE,
      name = NULL,
      labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
    ) +
    geom_point(
      data = gs_16day,
      aes(x = DateTime, y = Gs_mps, color = season),
      size = 2.0,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = season_cols,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = ylab_gs) +
    theme_pub(base_size = 16, legend_position = "right") +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))
  
  ggsave(out_fig_gs_16day, p_gs_16day, width = 14, height = 4.2, dpi = 300)
  cat("saved:", out_fig_gs_16day, "\n")
} else {
  cat("skipped 16-day gs plot (no finite gs).\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# daily Reco time series (g C m-2 d-1)
out_fig_reco_daily <- file.path(fig_path, "ts_daily_Reco_gC_m2_d_points_ENSOshade.png")

daily_reco <- df %>%
  filter(!is.na(Date), is.finite(Reco)) %>%
  group_by(Date) %>%
  summarise(
    Reco_umol_m2_s = median(Reco, na.rm = TRUE),
    season   = dplyr::first(na.omit(as.character(season))),
    ENSO     = dplyr::first(na.omit(as.character(ENSO))),
    n_hh_obs = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    Reco_gC_m2_day = Reco_umol_m2_s * 86400 * gC_per_umolCO2,
    DateTime = as.POSIXct(Date, tz = tz_local),
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  filter(is.finite(Reco_gC_m2_day), Reco_gC_m2_day >= 0) %>%
  arrange(DateTime)

ylab_reco_day <- expression(R[eco] ~ "(" * g ~ C ~ m^-2 ~ d^-1 * ")")

p_reco_day <- ggplot() +
  geom_rect(
    data = enso_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_point(
    data = daily_reco,
    aes(x = DateTime, y = Reco_gC_m2_day, color = season),
    size = 1.6,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_reco_day) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_reco_daily, p_reco_day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_reco_daily, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# 16-day window Reco time series (g C m-2 d-1)
out_fig_reco_16day <- file.path(fig_path, "ts_Reco_16day_gC_m2_d_points_ENSOshade.png")

# 16-day window Reco integrated diurnal sum (g C m^-2 d^-1)
reco_16day <- df %>%
  filter(!is.na(Date), is.finite(Reco)) %>%
  mutate(
    window_start = as.Date(floor_date(Date, unit = "16 days")),
    window_mid   = window_start + days(8),
    hh_slot      = format(DateTime, "%H:%M:%S")
  ) %>%
  # Group by window and specific half-hour slot to build the representative diel cycle
  group_by(window_start, window_mid, hh_slot) %>%
  summarise(
    Reco_hh_mean = mean(Reco, na.rm = TRUE),
    season       = dplyr::first(na.omit(as.character(season))),
    ENSO         = dplyr::first(na.omit(as.character(ENSO))),
    .groups = "drop_last"
  ) %>%
  # Sum the 48 slots (each representing 1800 seconds) to get the total daily mass flux
  summarise(
    # Reco is continuous over 24 h but is finite mostly in daytime half-hours
    # (nighttime is u*-filtered), so summing only the populated slots under-integrates
    # by ~2x. Extrapolate the mean respiration rate to the full day (matches daily Reco).
    Reco_gC_m2_day = mean(Reco_hh_mean, na.rm = TRUE) * 86400 * gC_per_umolCO2,
    season         = dplyr::first(season),
    ENSO           = dplyr::first(ENSO),
    .groups = "drop"
  ) %>%
  mutate(
    DateTime = as.POSIXct(window_mid, tz = tz_local),
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  filter(is.finite(Reco_gC_m2_day), Reco_gC_m2_day >= 0) %>%
  arrange(DateTime)

ylab_reco_16day <- expression(R[eco] ~ "(" * g ~ C ~ m^-2 ~ d^-1 * ")")

p_reco_16day <- ggplot() +
  geom_rect(
    data = enso_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_point(
    data = reco_16day,
    aes(x = DateTime, y = Reco_gC_m2_day, color = season),
    size = 2.0,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_reco_16day) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_reco_16day, p_reco_16day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_reco_16day, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# 16-day window GEP time series (g C m-2 d-1)
out_fig_gep_16day <- file.path(fig_path, "ts_GEP_16day_gC_m2_d_points_ENSOshade.png")

# 16-day window GEP integrated diurnal sum (g C m^-2 d^-1)
gep_16day <- df %>%
  filter(!is.na(Date), is.finite(GEP)) %>%
  mutate(
    window_start = as.Date(floor_date(Date, unit = "16 days")),
    window_mid   = window_start + days(8),
    hh_slot      = format(DateTime, "%H:%M:%S")
  ) %>%
  # Group by window and specific half-hour slot to build the representative diel cycle
  group_by(window_start, window_mid, hh_slot) %>%
  summarise(
    GEP_hh_mean = mean(GEP, na.rm = TRUE),
    season      = dplyr::first(na.omit(as.character(season))),
    ENSO        = dplyr::first(na.omit(as.character(ENSO))),
    .groups = "drop_last"
  ) %>%
  # Sum the 48 slots (each representing 1800 seconds) to get the total daily mass flux
  summarise(
    GEP_gC_m2_day = sum(GEP_hh_mean * 1800, na.rm = TRUE) * gC_per_umolCO2,
    season        = dplyr::first(season),
    ENSO          = dplyr::first(ENSO),
    .groups = "drop"
  ) %>%
  mutate(
    DateTime = as.POSIXct(window_mid, tz = tz_local),
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  filter(is.finite(GEP_gC_m2_day)) %>%
  arrange(DateTime)

deadband_gC <- 0.1
gep_16day <- gep_16day %>% filter(abs(GEP_gC_m2_day) >= deadband_gC)

deadband_gC <- 0.1
gep_16day <- gep_16day %>% filter(abs(GEP_gC_m2_day) >= deadband_gC)

ylab_gep_16day <- expression(GEP ~ "(" * g ~ C ~ m^-2 ~ d^-1 * ")")

p_gep_16day <- ggplot() +
  geom_rect(
    data = enso_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_point(
    data = gep_16day,
    aes(x = DateTime, y = GEP_gC_m2_day, color = season),
    size = 2.0,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_gep_16day) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_gep_16day, p_gep_16day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_gep_16day, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Supplement composite: the three 16-day time series in one figure.
#
# These were supplementary figures for GEP, Reco and Gs, one full-width strip
# each. Reviewers asked for a shorter supplement, so they are combined.
#
# This block is deliberately ADDITIVE: it rebuilds the three panels at a small
# base_size rather than reusing p_gep_16day / p_reco_16day / p_gs_16day, whose
# type is authored at base_size 16 for a 14-inch canvas and would be far too
# large once the panel is halved. The three standalone ggsave() calls above are
# untouched, so those PNGs stay byte-identical.
make_16day_panel <- function(df, yvar, ylab, base_size = 9) {
  ggplot() +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE, alpha = 0.9
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
    scale_fill_manual(
      values = enso_fill, drop = FALSE, name = NULL,
      labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
    ) +
    geom_point(
      data = df, aes(x = DateTime, y = .data[[yvar]], color = season),
      size = 2.0 * base_size / 16, alpha = 0.85
    ) +
    scale_color_manual(
      values = season_cols, drop = FALSE, name = NULL,
      guide = guide_legend(override.aes = list(size = 4.0 * base_size / 16, alpha = 1))
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = ylab) +
    theme_pub(base_size = base_size, legend_position = "bottom") +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))
}

ts16_parts <- list()
if (exists("gep_16day"))  ts16_parts$GEP  <- make_16day_panel(gep_16day,  "GEP_gC_m2_day",  ylab_gep_16day)
if (exists("reco_16day")) ts16_parts$Reco <- make_16day_panel(reco_16day, "Reco_gC_m2_day", ylab_reco_16day)
# gs_16day and ylab_gs are created inside a conditional further up, so guard both
if (exists("gs_16day") && exists("ylab_gs")) ts16_parts$Gs <- make_16day_panel(gs_16day, "Gs_mps", ylab_gs)

if (length(ts16_parts) == 3L) {
  ts16_composite <- patchwork::wrap_plots(ts16_parts, ncol = 1) +
    patchwork::plot_layout(guides = "collect", axes = "collect_x") +
    patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(legend.position = "bottom", plot.tag = element_text(size = 10, face = "bold"))

  out_ts16 <- file.path(fig_path, "FigS_16day_GEP_Reco_Gs_composite.png")
  ggsave(out_ts16, ts16_composite, width = 7.10, height = 6.40, dpi = 300, bg = "white")
  cat("saved:", out_ts16, "(replaces the three standalone 16-day series)\n")
} else {
  message("16-day composite skipped: only ", length(ts16_parts), " of 3 panels available")
}

# ───────────────────────────────────────────────────────────────────────────────
# 16-day window NEE_ok time series (g C m-2 d-1)
out_fig_nee_16day <- file.path(fig_path, "ts_NEEok_16day_gC_m2_d_points_ENSOshade.png")

nee_deadband_gC_m2_day <- 0.0

# 16-day window NEE integrated diurnal sum (g C m^-2 d^-1)
nee_16day <- df %>%
  filter(!is.na(Date), is.finite(NEE_ok)) %>%
  mutate(
    window_start = as.Date(floor_date(Date, unit = "16 days")),
    window_mid   = window_start + days(8),
    hh_slot      = format(DateTime, "%H:%M:%S")
  ) %>%
  # Group by window and specific half-hour slot to build the representative diel cycle
  group_by(window_start, window_mid, hh_slot) %>%
  summarise(
    NEE_hh_mean = mean(NEE_ok, na.rm = TRUE),
    season      = dplyr::first(na.omit(as.character(season))),
    ENSO        = dplyr::first(na.omit(as.character(ENSO))),
    .groups = "drop_last"
  ) %>%
  # Sum the 48 slots (each representing 1800 seconds) to get the total daily mass flux
  summarise(
    NEE_gC_m2_day = sum(NEE_hh_mean * 1800, na.rm = TRUE) * gC_per_umolCO2,
    season        = dplyr::first(season),
    ENSO          = dplyr::first(ENSO),
    .groups = "drop"
  ) %>%
  mutate(
    DateTime = as.POSIXct(window_mid, tz = tz_local),
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = factor(ENSO,   levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  filter(is.finite(NEE_gC_m2_day)) %>%
  arrange(DateTime)

ylab_nee_16day <- expression(NEE[ok] ~ "(" * g ~ C ~ m^-2 ~ d^-1 * ")")

p_nee_16day <- ggplot() +
  geom_rect(
    data = enso_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_point(
    data = nee_16day,
    aes(x = DateTime, y = NEE_gC_m2_day, color = season),
    size = 2.0,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_nee_16day) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_nee_16day, p_nee_16day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_nee_16day, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# 16-day daily C balance, with diel coverage checks

out_fig_dailyC_16day <- file.path(
  fig_path,
  "ts_16day_dailyCbalance_meanSD_seasoncolor_ENSOshade_checked.png"
)

out_csv_dailyC_16day <- file.path(
  output_path,
  "check_16day_dailyCbalance_checked.csv"
)

out_csv_coverage_16day <- file.path(
  output_path,
  "check_16day_dailyC_diel_coverage.csv"
)

out_csv_biomet_16day_fig3 <- file.path(
  output_path,
  "biomet_16day_timeseries_for_Fig3.csv"
)

out_fig_dailyC_biomet_16day <- file.path(
  fig_path,
  "Fig3_dailyCbalance_plus_biomet_16day_ENSOshade.png"
)

min_hh_per_day      <- 40
min_days_per_window <- 6
min_diel_coverage   <- 0.90

daily_cbal <- df %>%
  filter(!is.na(Date), is.finite(NEE_ok)) %>%
  group_by(Date) %>%
  summarise(
    n_hh = dplyr::n(),
    cbal_gC_m2_day = sum(NEE_ok * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
    season = dplyr::first(na.omit(as.character(season))),
    ENSO   = dplyr::first(na.omit(as.character(ENSO))),
    .groups = "drop"
  ) %>%
  filter(n_hh >= min_hh_per_day) %>%
  mutate(
    window_start = as.Date(lubridate::floor_date(Date, unit = "16 days")),
    window_mid   = window_start + lubridate::days(8)
  )

diel_16day <- df %>%
  filter(!is.na(DateTime), is.finite(NEE_ok)) %>%
  mutate(
    DateTime_30 = lubridate::floor_date(DateTime, unit = "30 minutes"),
    hour_str    = format(DateTime_30, "%H:%M:%S"),
    hour        = factor(hour_str, levels = halfhour_levels),
    window_start = as.Date(lubridate::floor_date(Date, unit = "16 days")),
    window_mid   = window_start + lubridate::days(8),
    season = factor(as.character(season), levels = c("dry", "wet")),
    ENSO   = as.character(ENSO)
  ) %>%
  group_by(window_start, window_mid, season, ENSO, hour) %>%
  summarise(
    NEE_ok_mean = mean(NEE_ok, na.rm = TRUE),
    n_obs_bin   = dplyr::n(),
    .groups = "drop"
  )

coverage_16day <- diel_16day %>%
  group_by(window_start, window_mid, season, ENSO) %>%
  summarise(
    n_halfhour_bins   = n_distinct(hour),
    frac_diel_covered = n_halfhour_bins / 48,
    missing_bins      = 48 - n_halfhour_bins,
    .groups = "drop"
  ) %>%
  arrange(frac_diel_covered)

write_csv(coverage_16day, out_csv_coverage_16day)
cat("saved:", out_csv_coverage_16day, "\n")

cat("\n16-day diel coverage summary:\n")
coverage_16day %>%
  summarise(
    n_windows = n(),
    min_coverage = min(frac_diel_covered, na.rm = TRUE),
    median_coverage = median(frac_diel_covered, na.rm = TRUE),
    n_below_90 = sum(frac_diel_covered < 0.90, na.rm = TRUE),
    n_below_75 = sum(frac_diel_covered < 0.75, na.rm = TRUE)
  ) %>%
  print()

cat("\nlow-coverage windows:\n")
coverage_16day %>%
  filter(frac_diel_covered < min_diel_coverage) %>%
  arrange(frac_diel_covered) %>%
  print(n = 50)

# Data volume behind each window mean.
#
# A co-author read the "mean valid days" column of Table 1 (~10.4) as the sample
# size and concluded the record was too sparse. It is not the sample size:
# dailyC_mean below is the integral of a MEAN DIEL CYCLE over 48 bins, and every
# finite half-hour in the window feeds those bins -- including half-hours on days
# that n_days excludes, since n_days counts only days with >= min_hh_per_day (40)
# of 48 half-hours and exists to support dailyC_sd and a conservative screen.
# n_hh_window and obs_per_bin are the quantities that actually govern whether the
# mean is well sampled, so they are carried through to Table 1 alongside n_days.
cbal_16day_raw <- diel_16day %>%
  group_by(window_start, window_mid, season, ENSO) %>%
  summarise(
    dailyC_mean = sum(NEE_ok_mean * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
    n_hh_window = sum(n_obs_bin),          # half-hourly observations behind the window
    obs_per_bin = mean(n_obs_bin),         # ... averaged over the 48 diel bins
    .groups = "drop"
  ) %>%
  left_join(
    # distinct days contributing ANY half-hour, as opposed to n_days below
    df %>%
      filter(!is.na(DateTime), is.finite(NEE_ok)) %>%
      mutate(
        window_start = as.Date(lubridate::floor_date(Date, unit = "16 days")),
        window_mid   = window_start + lubridate::days(8),
        season = factor(as.character(season), levels = c("dry", "wet")),
        ENSO   = as.character(ENSO)
      ) %>%
      group_by(window_start, window_mid, season, ENSO) %>%
      summarise(n_days_any = dplyr::n_distinct(Date), .groups = "drop"),
    by = c("window_start", "window_mid", "season", "ENSO")
  ) %>%
  left_join(
    daily_cbal %>%
      group_by(window_start, window_mid, season, ENSO) %>%
      summarise(
        dailyC_sd = sd(cbal_gC_m2_day, na.rm = TRUE),
        n_days    = dplyr::n(),
        .groups = "drop"
      ),
    by = c("window_start", "window_mid", "season", "ENSO")
  ) %>%
  left_join(
    coverage_16day,
    by = c("window_start", "window_mid", "season", "ENSO")
  ) %>%
  mutate(
    DateTime = as.POSIXct(window_mid, tz = tz_local),
    season = factor(as.character(season), levels = c("dry", "wet")),
    ENSO   = factor(as.character(ENSO), levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  arrange(DateTime)

cat("\nchecking duplicate 16-day carbon-balance windows:\n")
cbal_16day_raw %>%
  count(window_start, window_mid, season, ENSO) %>%
  filter(n > 1) %>%
  print(n = 50)

cbal_16day_checked <- cbal_16day_raw %>%
  filter(
    n_days >= min_days_per_window,
    frac_diel_covered >= min_diel_coverage
  )

write_csv(cbal_16day_checked, out_csv_dailyC_16day)
cat("saved:", out_csv_dailyC_16day, "\n")

cat("\n16-day C balance sign summary after filtering:\n")
cbal_16day_checked %>%
  summarise(
    n_windows = n(),
    n_source = sum(dailyC_mean > 0, na.rm = TRUE),
    n_sink   = sum(dailyC_mean < 0, na.rm = TRUE),
    frac_source = mean(dailyC_mean > 0, na.rm = TRUE),
    mean_dailyC = mean(dailyC_mean, na.rm = TRUE),
    median_dailyC = median(dailyC_mean, na.rm = TRUE)
  ) %>%
  print()

# Same tally broken out by ENSO, by ENSO x season, and by El Nino episode.
# Section 3.1 reports the *consistency* of the sign across windows, which is what
# Fig. 3 uniquely shows and which the ENSO x season means in Table 1 do not
# convey: every La Nina window is a source, whereas El Nino reverses between the
# 2018 and 2024 episodes.
window_sign_summary <- function(d, ...) {
  d %>%
    group_by(...) %>%
    summarise(
      n_windows   = n(),
      n_source    = sum(dailyC_mean > 0, na.rm = TRUE),
      n_sink      = sum(dailyC_mean < 0, na.rm = TRUE),
      mean_dailyC = mean(dailyC_mean, na.rm = TRUE),
      min_dailyC  = min(dailyC_mean, na.rm = TRUE),
      max_dailyC  = max(dailyC_mean, na.rm = TRUE),
      .groups = "drop"
    )
}

cbal_sign_by_enso <- bind_rows(
  window_sign_summary(cbal_16day_checked, ENSO) %>% mutate(grouping = "ENSO", .before = 1),
  window_sign_summary(cbal_16day_checked, ENSO, season) %>% mutate(grouping = "ENSO x season", .before = 1),
  window_sign_summary(
    cbal_16day_checked %>%
      filter(ENSO == "El Nino") %>%
      mutate(ENSO = if_else(lubridate::year(window_mid) >= 2023, "El Nino 2023-24", "El Nino 2018-19")),
    ENSO
  ) %>% mutate(grouping = "El Nino episode", .before = 1)
)

write_csv(cbal_sign_by_enso, file.path(output_path, "table_16day_window_sign_by_ENSO.csv"))
cat("\n16-day window sign tally (screened windows, the ones Fig 3 plots):\n")
print(as.data.frame(cbal_sign_by_enso %>% mutate(across(where(is.numeric), ~round(.x, 2)))),
      row.names = FALSE)

ylab_cbal <- "Daily net C flux (g C m\u207B\u00B2 d\u207B\u00B9)"

p_cbal_16day <- ggplot() +
  geom_rect(
    data = enso_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black", alpha = 0.6) +
  scale_fill_manual(
    values = enso_fill,
    drop = FALSE,
    name = NULL,
    labels = c("El Nino" = "El Niño", "La Nina" = "La Niña", "neutral" = "neutral")
  ) +
  geom_errorbar(
    data = cbal_16day_checked,
    aes(
      x = DateTime,
      ymin = dailyC_mean - dailyC_sd,
      ymax = dailyC_mean + dailyC_sd,
      color = season
    ),
    width = 0,
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    data = cbal_16day_checked,
    aes(x = DateTime, y = dailyC_mean, color = season, shape = season),
    size = 3.2,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_shape_manual(values = season_shapes, drop = FALSE, name = NULL) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_cbal) +
  # sign-convention arrows, placed in the empty 2020 data gap so they overlap no points
  annotate("segment", x = ymd_hms("2020-01-01 00:00:00", tz = tz_local),
           xend = ymd_hms("2020-01-01 00:00:00", tz = tz_local), y = 0.5, yend = 3.4,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
           color = "grey25", linewidth = 0.7) +
  annotate("text", x = ymd_hms("2020-01-01 00:00:00", tz = tz_local), y = 4.2,
           label = "net release", size = 4.2, color = "grey25") +
  annotate("segment", x = ymd_hms("2020-01-01 00:00:00", tz = tz_local),
           xend = ymd_hms("2020-01-01 00:00:00", tz = tz_local), y = -0.5, yend = -2.4,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
           color = "grey25", linewidth = 0.7) +
  annotate("text", x = ymd_hms("2020-01-01 00:00:00", tz = tz_local), y = -3.0,
           label = "net uptake", size = 4.2, color = "grey25") +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_dailyC_16day, p_cbal_16day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_dailyC_16day, "\n")

biomet_specs_fig3 <- tibble::tribble(
  ~variable, ~source_column, ~label, ~value_type,
  "Ta", "TA_1_1_1", "atop(italic(T)[plain(a)], (degree*C))", "16-day mean",
  "VPD", "VPD_kPa", "atop(VPD, (kPa))", "16-day mean",
  "PAR", "PAR", "atop(PAR, (mu*mol~m^{-2}~s^{-1}))", "16-day mean",
  "SWC", "SWC_1_1_1", "atop(SWC, (cm))", "16-day mean",
  "Ts", "TS_3", "atop(italic(T)[plain(s)], (degree*C))", "16-day mean"
)

summarise_biomet_16day <- function(d, specs) {
  available_specs <- specs %>%
    filter(source_column %in% names(d))

  missing_specs <- specs %>%
    filter(!source_column %in% names(d))

  if (nrow(missing_specs) > 0) {
    cat(
      "warning: skipping missing Fig3 biomet columns:",
      paste(missing_specs$source_column, collapse = ", "),
      "\n"
    )
  }

  bind_rows(lapply(seq_len(nrow(available_specs)), function(i) {
    spec <- available_specs[i, ]

    d %>%
      transmute(
        Date,
        value = suppressWarnings(as.numeric(.data[[spec$source_column]]))
      ) %>%
      filter(!is.na(Date), is.finite(value)) %>%
      mutate(
        window_start = as.Date(lubridate::floor_date(Date, unit = "16 days")),
        window_mid = window_start + lubridate::days(8)
      ) %>%
      group_by(window_start, window_mid) %>%
      summarise(
        value = mean(value, na.rm = TRUE),
        n_obs = dplyr::n(),
        .groups = "drop"
      ) %>%
      mutate(
        variable = spec$variable,
        source_column = spec$source_column,
        label = spec$label,
        value_type = spec$value_type
      )
  }))
}

parse_imerg_datetime_local <- function(imerg) {
  local_col <- find_col(c("time_local_minus05", "time_local"), names(imerg))
  utc_col <- find_col(c("time_utc_iso", "time_utc", "time", "system:time_start"), names(imerg))

  if (!is.na(local_col)) {
    return(suppressWarnings(lubridate::parse_date_time(
      imerg[[local_col]],
      orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM"),
      tz = tz_local
    )))
  }

  if (!is.na(utc_col)) {
    t_utc <- suppressWarnings(lubridate::parse_date_time(
      imerg[[utc_col]],
      orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM"),
      tz = "UTC"
    ))
    return(lubridate::with_tz(t_utc, tz_local))
  }

  stop("could not find a recognized IMERG time column")
}

summarise_precip_16day <- function(input_folder) {
  imerg_csv <- list.files(input_folder, pattern = "(?i)imerg.*30min.*\\.csv$", full.names = TRUE)
  if (!length(imerg_csv)) {
    cat("warning: skipping Fig3 precipitation: no IMERG 30-minute CSV found\n")
    return(tibble())
  }

  imerg_csv <- imerg_csv[order(file.info(imerg_csv)$mtime, decreasing = TRUE)][1]
  imerg <- readr::read_csv(imerg_csv, show_col_types = FALSE)
  DateTime <- parse_imerg_datetime_local(imerg)

  if ("precip_mm_30min" %in% names(imerg)) {
    precip_30min <- suppressWarnings(as.numeric(imerg$precip_mm_30min))
    precip_source <- "precip_mm_30min"
  } else {
    rate_col <- find_col(c("precip_mm_per_hr", "precipitation", "precipitationCal"), names(imerg))
    if (is.na(rate_col)) {
      cat("warning: skipping Fig3 precipitation: no recognized precipitation column\n")
      return(tibble())
    }
    precip_30min <- suppressWarnings(as.numeric(imerg[[rate_col]])) * 0.5
    precip_source <- paste0(rate_col, " converted from mm hr-1")
  }

  tibble(DateTime = DateTime, precip_mm_30min = precip_30min) %>%
    filter(
      !is.na(DateTime),
      DateTime <= tmax,
      is.finite(precip_mm_30min)
    ) %>%
    mutate(
      Date = as.Date(DateTime, tz = tz_local),
      window_start = as.Date(lubridate::floor_date(Date, unit = "16 days")),
      window_mid = window_start + lubridate::days(8)
    ) %>%
    group_by(window_start, window_mid) %>%
    summarise(
      value = sum(precip_mm_30min, na.rm = TRUE),
      n_obs = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      variable = "Precipitation",
      source_column = precip_source,
      label = "atop(\"Precip.\", \"mm/16 d\")",
      value_type = "16-day total"
    )
}

cbal_windows_fig3 <- cbal_16day_checked %>%
  transmute(
    window_start,
    window_mid,
    DateTime,
    season = factor(as.character(season), levels = c("dry", "wet")),
    ENSO = factor(as.character(ENSO), levels = c("El Nino", "La Nina", "neutral"))
  ) %>%
  distinct()

biomet_16day_flux <- summarise_biomet_16day(df, biomet_specs_fig3)
precip_16day_fig3 <- summarise_precip_16day(input_folder)

biomet_panel_order_fig3 <- c(
  "atop(italic(T)[plain(a)], (degree*C))",
  "atop(VPD, (kPa))",
  "atop(\"Precip.\", \"mm/16 d\")",
  "atop(PAR, (mu*mol~m^{-2}~s^{-1}))",
  "atop(SWC, (cm))",
  "atop(italic(T)[plain(s)], (degree*C))"
)

biomet_16day_fig3 <- bind_rows(biomet_16day_flux, precip_16day_fig3) %>%
  inner_join(cbal_windows_fig3, by = c("window_start", "window_mid")) %>%
  mutate(
    label = factor(label, levels = biomet_panel_order_fig3),
    variable = factor(variable, levels = c("Ta", "VPD", "Precipitation", "PAR", "SWC", "Ts"))
  ) %>%
  arrange(variable, DateTime)

write_csv(biomet_16day_fig3, out_csv_biomet_16day_fig3)
cat("saved:", out_csv_biomet_16day_fig3, "\n")

if (nrow(biomet_16day_fig3) > 0) {
  p_cbal_16day_combined <- p_cbal_16day +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  p_biomet_16day <- ggplot() +
    geom_rect(
      data = enso_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ENSO),
      inherit.aes = FALSE,
      alpha = 0.9
    ) +
    # no ENSO key here: the labelled strip above the panels names the phases, so
    # a colour-only legend for them would be redundant
    scale_fill_manual(values = enso_fill, drop = FALSE, guide = "none") +
    geom_point(
      data = biomet_16day_fig3,
      aes(x = DateTime, y = value, color = season, shape = season),
      size = 2.2,
      alpha = 0.85
    ) +
    facet_grid(label ~ ., scales = "free_y", switch = "y", labeller = label_parsed) +
    scale_color_manual(
      values = season_cols,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(order = 2, override.aes = list(size = 3.0, alpha = 1))
    ) +
    scale_shape_manual(
      values = season_shapes,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(order = 2)
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = NULL) +
    theme_pub(base_size = 12, legend_position = c(0.04, 0.08)) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 11, lineheight = 0.95),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 10),
      legend.justification = c(0, 0),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 12),
      legend.background = element_rect(fill = "white", color = "grey60"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.spacing.y = unit(0.25, "lines"),
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    ) +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))

  # ENSO phase strip above the panels.
  #
  # The background bands encode ENSO phase by fill colour alone. Even with the
  # separated tints in palette.R the three fills differ by only ~9 percentage
  # points of greyscale luminance, so the phase is unreadable in black and white
  # and marginal for a reader with a colour-vision deficiency (a co-author
  # flagged exactly this). Naming the phases in a strip makes the encoding
  # colour-independent; the fills then just reinforce it.
  enso_strip <- enso_rect %>%
    mutate(
      xmin = pmax(xmin, xmin_plot),
      xmax = pmin(xmax, xmax_plot),
      width_days = as.numeric(difftime(xmax, xmin, units = "days")),
      xmid = xmin + (xmax - xmin) / 2,
      label = dplyr::recode(ENSO, "El Nino" = "El Niño", "La Nina" = "La Niña",
                            "neutral" = "neutral", .default = ENSO)
    ) %>%
    filter(xmax > xmin)

  # a label needs roughly this much room; narrower bands stay unlabelled rather
  # than overplotting their neighbours
  min_days_for_label <- 150
  dropped <- enso_strip %>% filter(width_days < min_days_for_label)
  if (nrow(dropped)) {
    cat("Fig3 ENSO strip: ", nrow(dropped), " band(s) too narrow to label (<",
        min_days_for_label, " days): ",
        paste(sprintf("%s %s", dropped$label, format(dropped$xmin, "%Y-%m")), collapse = ", "),
        "\n", sep = "")
  }

  p_enso_strip <- ggplot(enso_strip) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = ENSO),
              colour = "grey70", linewidth = 0.2, alpha = 0.9) +
    geom_text(data = enso_strip %>% filter(width_days >= min_days_for_label),
              aes(x = xmid, y = 0.5, label = label),
              size = 3.9, colour = "grey15", vjust = 0.5) +
    scale_fill_manual(values = enso_fill, drop = FALSE, guide = "none") +
    # same scale + coord as the panels below, so patchwork lines the strip up
    # with them; coord_cartesian keeps ggplot's default 5% expansion on both
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(2, 5.5, 1, 5.5)) +
    coord_cartesian(xlim = c(xmin_plot, xmax_plot))

  p_dailyC_biomet_16day <- p_enso_strip / p_cbal_16day_combined / p_biomet_16day +
    plot_layout(heights = c(0.16, 1.25, 2.8))

  ggsave(out_fig_dailyC_biomet_16day, p_dailyC_biomet_16day, width = 14, height = 10.5, dpi = 300)
  cat("saved:", out_fig_dailyC_biomet_16day, "\n")
} else {
  cat("warning: no Fig3 biomet data available; skipping combined C balance + biomet figure\n")
}

# ───────────────────────────────────────────────────────────────────────────────
# save outputs
write_csv(df, out_csv_enriched)
cat("saved (enriched passthrough):", out_csv_enriched, "\n")

write_csv(pc_16day_strict, out_csv_pc_strict)
cat("saved:", out_csv_pc_strict, "\n")

write_csv(pc_16day_fallback, out_csv_pc_fallback)
cat("saved:", out_csv_pc_fallback, "\n")

write_csv(pc_daily_strict, out_csv_pc_daily_strict)
cat("saved:", out_csv_pc_daily_strict, "\n")

write_csv(pc_daily_fallback, out_csv_pc_daily_fallback)
cat("saved:", out_csv_pc_daily_fallback, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# plot pc 16-day with enso shading + season-colored points
ylab_pc <- expression(P[c] ~ "(" * g ~ C ~ m^-2 ~ d^-1 * ")")

p_pc_16day <- ggplot() +
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
    data = pc_16day_plot_df,
    aes(
      x = as.POSIXct(window_mid, tz = tz_local),
      y = Pc_gC_m2_day,
      color = season
    ),
    size = 2.4,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_pc) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_pc_16day, p_pc_16day, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_pc_16day, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# plot pc daily with enso shading + season-colored points (Option A)
p_pc_daily <- ggplot() +
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
    data = pc_daily_plot_df,
    aes(
      x = DateTime,
      y = Pc_gC_m2_day,
      color = season
    ),
    size = 1.6,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = season_cols,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(override.aes = list(size = 4.0, alpha = 1))
  ) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = ylab_pc) +
  theme_pub(base_size = 16, legend_position = "right") +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot))

ggsave(out_fig_pc_daily, p_pc_daily, width = 14, height = 4.2, dpi = 300)
cat("saved:", out_fig_pc_daily, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# restored old-style diel figure families
run_restored_diel_figures(df)

# ───────────────────────────────────────────────────────────────────────────────
# export 48-pt monthly diel means (by month × season × ENSO)
out_csv_48pt <- file.path(output_path, "tambopata_48points_per_month.csv")

export_diel_month <- df %>%
  filter(
    !is.na(DateTime),
    !is.na(NEE_ok)
  ) %>%
  mutate(
    ym = format(DateTime, "%Y-%m"),
    DateTime_30 = floor_date(DateTime, unit = "30 minutes"),
    hour_str = format(DateTime_30, "%H:%M:%S"),
    hour = factor(hour_str, levels = halfhour_levels)
  ) %>%
  group_by(ym, season, ENSO, hour) %>%
  summarise(
    NEE_ok = mean(NEE_ok, na.rm = TRUE),
    GEP    = mean(GEP,  na.rm = TRUE),
    Reco   = mean(Reco, na.rm = TRUE),
    TA_1_1_1      = mean(Tair, na.rm = TRUE),
    VPD_kPa       = mean(VPD,  na.rm = TRUE),
    PAR = mean(PAR, na.rm = TRUE),
    PAR_corrected_SWin = mean(PAR_corrected_SWin, na.rm = TRUE),
    PAR_source = paste(sort(unique(na.omit(PAR_source))), collapse = ";"),
    PPFD_IN_1_1_1 = if ("PPFD_IN_1_1_1" %in% names(df)) {
      mean(PPFD_IN_1_1_1, na.rm = TRUE)
    } else {
      NA_real_
    },
    SW_IN = mean(SW_IN, na.rm = TRUE),
    SWin = mean(SWin, na.rm = TRUE),
    SWC_1_1_1 = if ("SWC_1_1_1" %in% names(df)) mean(SWC_1_1_1, na.rm = TRUE) else NA_real_,
    TS_3      = if ("TS_3"      %in% names(df)) mean(TS_3,      na.rm = TRUE) else NA_real_,
    Rn    = if ("Rn"    %in% names(df)) mean(Rn,    na.rm = TRUE) else NA_real_,
    FC    = if ("FC"    %in% names(df)) mean(FC,    na.rm = TRUE) else NA_real_,
    H     = if ("H"     %in% names(df)) mean(H,     na.rm = TRUE) else NA_real_,
    LE    = if ("LE"    %in% names(df)) mean(LE,    na.rm = TRUE) else NA_real_,
    USTAR = if ("USTAR" %in% names(df)) mean(USTAR, na.rm = TRUE) else NA_real_,
    WS    = if ("WS_1_1_1" %in% names(df)) mean(WS_1_1_1, na.rm = TRUE) else NA_real_,
    Gs_mps = if ("Gs_mps" %in% names(df)) mean(Gs_mps, na.rm = TRUE) else NA_real_,
    n_obs = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(ym, hour)

write_csv(export_diel_month, out_csv_48pt)
cat("saved:", out_csv_48pt, "\n")

# monthly diurnal sum from 48 half-hour monthly mean diel cycles
out_csv_monthly_diurnal_sum <- file.path(output_path, "tambopata_monthly_diurnal_sum_gC_m2_day.csv")
out_fig_monthly_diurnal_sum <- file.path(fig_path, "FigX_monthly_diurnal_sum_gC_m2_day_yearshape_keep_ENSO.png")

monthly_diurnal_sum <- export_diel_month %>%
  filter(!is.na(NEE_ok)) %>%
  group_by(ym, season, ENSO) %>%
  summarise(
    gC_m2_day = sum(NEE_ok * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
    n_bins = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    year = as.integer(substr(ym, 1, 4)),
    month = as.integer(substr(ym, 6, 7)),
    carbon_balance = case_when(
      gC_m2_day < 0 ~ "sink",
      gC_m2_day > 0 ~ "source",
      TRUE ~ "neutral"
    )
  ) %>%
  arrange(year, month, season, ENSO)

write_csv(monthly_diurnal_sum, out_csv_monthly_diurnal_sum)
cat("saved:", out_csv_monthly_diurnal_sum, "\n")

monthly_ts <- monthly_diurnal_sum %>%
  mutate(
    date_midmonth = as.Date(paste0(ym, "-15")),
    year_f = factor(year),
    season = factor(as.character(season), levels = c("dry", "wet"))
  ) %>%
  arrange(date_midmonth, season, ENSO)

shape_pool <- c(16, 15, 17, 18, 3, 4, 8, 7, 9, 10, 11, 12, 13, 14)
year_levels <- levels(monthly_ts$year_f)
year_shapes <- setNames(rep(shape_pool, length.out = length(year_levels)), year_levels)

p_monthly_c <- ggplot(
  monthly_ts,
  aes(x = date_midmonth, y = gC_m2_day, color = season, shape = year_f)
) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.95) +
  scale_color_manual(values = c("dry" = "red", "wet" = "blue"), breaks = c("dry", "wet")) +
  scale_shape_manual(values = year_shapes) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    legend.box = "vertical"
  ) +
  labs(
    y = expression(paste("Daily C balance from monthly mean diel cycle (g C ", m^-2, " ", day^-1, ")")),
    color = NULL,
    shape = NULL
  )

ggsave(out_fig_monthly_diurnal_sum, p_monthly_c, width = 12, height = 4.5, dpi = 300)
cat("saved:", out_fig_monthly_diurnal_sum, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# table-1 style means by ENSO × season (using 16-day window products)
reco_tbl <- reco_16day %>%
  transmute(
    window_mid = as.Date(DateTime, tz = tz_local),
    ENSO = as.character(ENSO),
    season = as.character(season),
    Reco_gC_m2_day = Reco_gC_m2_day
  )

gep_tbl <- gep_16day %>%
  transmute(
    window_mid = as.Date(DateTime, tz = tz_local),
    ENSO = as.character(ENSO),
    season = as.character(season),
    GEP_gC_m2_day = GEP_gC_m2_day
  )

nee_tbl <- cbal_16day_checked %>%
  transmute(
    window_mid = as.Date(DateTime, tz = tz_local),
    ENSO = as.character(ENSO),
    season = as.character(season),
    NEE_gC_m2_day = dailyC_mean,
    NEE_sd_within_window = dailyC_sd,
    n_days_NEE = n_days,
    frac_diel_covered = frac_diel_covered
  )

pc_tbl_16day <- pc_16day_plot_df %>%
  transmute(
    window_mid = as.Date(window_mid),
    ENSO = as.character(ENSO),
    season = as.character(season),
    Pc_gC_m2_day = Pc_gC_m2_day
  )

gs_tbl <- gs_16day %>%
  transmute(
    window_mid = as.Date(DateTime, tz = tz_local),
    ENSO = as.character(ENSO),
    season = as.character(season),
    Gs_mps = Gs_mps
  )

# ───────────────────────────────────────────────────────────────────────────────
# table-1 style means by ENSO × season using matched 16-day windows only

nee_summary_table <- cbal_16day_checked %>%
  group_by(ENSO, season) %>%
  summarise(
    NEE_mean = mean(dailyC_mean, na.rm = TRUE),
    NEE_sd   = sd(dailyC_mean, na.rm = TRUE),
    mean_n_days_NEE = mean(n_days, na.rm = TRUE),
    min_n_days_NEE = min(n_days, na.rm = TRUE),
    mean_diel_coverage_NEE = mean(frac_diel_covered, na.rm = TRUE),
    # the data volume actually behind each mean (see note at cbal_16day_raw)
    mean_n_hh_NEE       = mean(n_hh_window, na.rm = TRUE),
    mean_obs_per_bin    = mean(obs_per_bin, na.rm = TRUE),
    mean_n_days_any_NEE = mean(n_days_any, na.rm = TRUE),
    total_n_hh_NEE      = sum(n_hh_window, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  arrange(season, ENSO)

nee_table_presentation <- nee_summary_table %>%
  transmute(
    Season = season,
    ENSO = ENSO,
    `NEE daily C balance (g C m-2 d-1)` = paste0(
      round_smart(NEE_mean, 2), " (", round_smart(NEE_sd, 2), ")"
    ),
    `mean valid days` = round_smart(mean_n_days_NEE, 1),
    `min valid days` = min_n_days_NEE,
    `mean diel coverage` = round_smart(mean_diel_coverage_NEE, 2),
    # printed next to "mean valid days" in Table 1 so the conservative screening
    # statistic and the actual data volume are never seen in isolation
    `mean obs per diel bin` = round_smart(mean_obs_per_bin, 1),
    `mean half-hours per window` = round(mean_n_hh_NEE),
    n_windows = n_windows
  )

cat("\n=== data volume behind Table 1 (answers the 'only 10.4 valid days' reading) ===\n")
print(as.data.frame(nee_summary_table %>% transmute(
  season, ENSO, n_windows,
  `mean n_days`   = round(mean_n_days_NEE, 1),
  `mean days any` = round(mean_n_days_any_NEE, 1),
  `half-hours`    = round(mean_n_hh_NEE),
  `obs/bin`       = round(mean_obs_per_bin, 1),
  `diel cov`      = round(mean_diel_coverage_NEE, 3))), row.names = FALSE)
cat(sprintf(
  paste0("  pooled: %d windows | %d half-hours | %.0f per window | %.1f per diel bin\n",
         "          n_days sum %d (mean %.2f) vs distinct days contributing %d\n",
         "  n_days counts only days with >= %d of 48 half-hours; it is a screen, not the sample size.\n"),
  nrow(cbal_16day_checked), sum(cbal_16day_checked$n_hh_window),
  mean(cbal_16day_checked$n_hh_window), mean(cbal_16day_checked$obs_per_bin),
  sum(cbal_16day_checked$n_days), mean(cbal_16day_checked$n_days),
  sum(cbal_16day_checked$n_days_any), min_hh_per_day))

write_csv(
  nee_summary_table,
  file.path(output_path, "nee_dailyC_16day_by_ENSO_season.csv")
)

write_csv(
  nee_table_presentation,
  file.path(output_path, "nee_dailyC_16day_by_ENSO_season_presentation.csv")
)


component_16day_joined <- inner_join(
  reco_tbl,
  gep_tbl,
  by = c("window_mid", "ENSO", "season")
) %>%
  left_join(
    gs_tbl,
    by = c("window_mid", "ENSO", "season")
  ) %>%
  mutate(
    season = factor(season, levels = c("dry", "wet")),
    ENSO   = factor(ENSO, levels = c("El Nino", "La Nina", "neutral"))
  )

component_summary_table <- component_16day_joined %>%
  group_by(ENSO, season) %>%
  summarise(
    GEP_mean  = mean(GEP_gC_m2_day, na.rm = TRUE),
    Reco_mean = mean(Reco_gC_m2_day, na.rm = TRUE),
    Gs_mean   = mean(Gs_mps, na.rm = TRUE),
    
    GEP_sd  = sd(GEP_gC_m2_day, na.rm = TRUE),
    Reco_sd = sd(Reco_gC_m2_day, na.rm = TRUE),
    Gs_sd   = sd(Gs_mps, na.rm = TRUE),
    
    n_windows = n(),
    .groups = "drop"
  ) %>%
  arrange(season, ENSO)

component_table_presentation <- component_summary_table %>%
  transmute(
    Season = season,
    ENSO = ENSO,
    `GEP (g C m-2 d-1)` = paste0(
      round_smart(GEP_mean, 2), " (", round_smart(GEP_sd, 2), ")"
    ),
    `Reco (g C m-2 d-1)` = paste0(
      round_smart(Reco_mean, 2), " (", round_smart(Reco_sd, 2), ")"
    ),
    `Gs (m s-1)` = paste0(
      round_smart(Gs_mean, 4), " (", round_smart(Gs_sd, 4), ")"
    ),
    n_windows = n_windows
  )

write_csv(
  component_summary_table,
  file.path(output_path, "gep_reco_gs_16day_by_ENSO_season.csv")
)

write_csv(
  component_table_presentation,
  file.path(output_path, "gep_reco_gs_16day_by_ENSO_season_presentation.csv")
)

cat("\n=== NEE-ONLY TABLE PREVIEW ===\n")
print(as.data.frame(nee_table_presentation), row.names = FALSE)

cat("\n=== COMPONENT TABLE PREVIEW ===\n")
print(as.data.frame(component_table_presentation), row.names = FALSE)

# ── Extra diagnostic checks ────────────────────────────────────────────────────


coverage_16day %>%
  summarise(
    n_windows = n(),
    min_coverage = min(frac_diel_covered, na.rm = TRUE),
    median_coverage = median(frac_diel_covered, na.rm = TRUE),
    n_below_90 = sum(frac_diel_covered < 0.90, na.rm = TRUE),
    n_below_75 = sum(frac_diel_covered < 0.75, na.rm = TRUE)
  )

coverage_16day %>%
  filter(frac_diel_covered < 0.90) %>%
  arrange(frac_diel_covered) %>%
  print(n = 50)

cbal_16day_checked %>%
  summarise(
    n_windows = n(),
    n_source = sum(dailyC_mean > 0, na.rm = TRUE),
    n_sink = sum(dailyC_mean < 0, na.rm = TRUE),
    frac_source = mean(dailyC_mean > 0, na.rm = TRUE),
    mean_dailyC = mean(dailyC_mean, na.rm = TRUE),
    median_dailyC = median(dailyC_mean, na.rm = TRUE)
  )

coverage_16day %>%
  summarise(
    n_windows = n(),
    min_coverage = min(frac_diel_covered, na.rm = TRUE),
    median_coverage = median(frac_diel_covered, na.rm = TRUE),
    n_below_90 = sum(frac_diel_covered < 0.90, na.rm = TRUE),
    n_below_75 = sum(frac_diel_covered < 0.75, na.rm = TRUE)
  )
