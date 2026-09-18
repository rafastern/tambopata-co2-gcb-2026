# ───────────────────────────────────────────────────────────────────────────────
# merged: Fig6-style binned boxplots for multiple response variables (48 points per month)
# - reads tambopata_48points_per_month.csv once
# - for each response (GEP, NEE_ok, LE, Gs_mps) generates + saves:
#     * 8 single panels (Ta day; VPD day/night; SWC day/night; Ts day/night; Ta night)
#     * 1 composite figure:
#         - GEP / LE / Gs_mps: daytime-only 2×2 composite with shared y-grob
#         - NEE_ok: 2×2 composite (day panels + night panels) with collected y title
# updates in this version:
#   - legend moved to top-left in all SINGLE panels and in the NEE composite
#   - adds "Daytime"/"Nighttime" label INSIDE each panel, top-center (true panel center via npc)
#   - caps DISPLAYED y-range for Gs only (coord_cartesian: keeps stats/outliers)
#       * Gs: ylim = c(-0.2, 0.2)
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(patchwork)
  library(grid)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# data (48 points per month)
flux_fp <- file.path(input_path, "tambopata_48points_per_month.csv")
df0_raw <- read_csv(flux_fp, show_col_types = FALSE)

# nee_source: swap the NEE response to non-gap-filled measured NEE (primary), joined from the
#   measured 48-point composite by (ym, season, ENSO, hour); GEP/LE/Gs stay gap-filled.
#   "gapfilled" leaves NEE_ok as-is (sensitivity).
nee_source <- "measured"
if (identical(nee_source, "measured")) {
  meas_fp <- file.path(input_path, "tambopata_48points_per_month_measured.csv")
  meas <- readr::read_csv(meas_fp, show_col_types = FALSE) %>%
    dplyr::mutate(hour = as.character(hour)) %>%
    dplyr::select(ym, season, ENSO, hour, NEE_meas)
  # Keep the gap-filled NEE_ok AND join the measured NEE_meas as a SEPARATE column.
  # Daytime panels use measured NEE (assigned after df_day_ta below); nighttime panels
  # keep the gap-filled NEE_ok, because measured nighttime NEE is dominated by storage/
  # advection artifacts (raw range ~ -90..+120 umol m-2 s-1) and is unreliable for respiration.
  df0_raw <- df0_raw %>%
    dplyr::mutate(hour = as.character(hour)) %>%
    dplyr::left_join(meas, by = c("ym", "season", "ENSO", "hour"))
}

# ───────────────────────────────────────────────────────────────────────────────
# shared config
tz_local <- "America/Lima"

axis_title_size <- 15
axis_text_size <- 12
annotation_text_size <- 12
legend_text_size <- 12
plot_tag_size <- 15
shared_axis_title_size <- 15

ta_axis_label <- expression(italic(T)[plain(a)] ~ "(°C)")
ts_axis_label <- expression(italic(T)[plain(s)] ~ "(°C)")

ta_col  <- "TA_1_1_1"        # °C
par_col <- "PAR"             # µmol m^-2 s^-1
swin_col <- "SW_IN"          # incoming shortwave, W m^-2 (night criterion: Rg < 10)
vpd_col <- "VPD_kPa"         # kPa
swc_col <- "SWC_1_1_1"
ts_col  <- "TS_3"

if (!"PAR_corrected_SWin" %in% names(df0_raw) && "PPFD_IN_1_1_1" %in% names(df0_raw)) {
  df0_raw$PAR_corrected_SWin <- suppressWarnings(as.numeric(df0_raw$PPFD_IN_1_1_1))
}
if (!par_col %in% names(df0_raw) && "PAR_corrected_SWin" %in% names(df0_raw)) {
  df0_raw$PAR <- suppressWarnings(as.numeric(df0_raw$PAR_corrected_SWin))
}
if (!par_col %in% names(df0_raw)) {
  stop("missing MATLAB-corrected PAR column in 48-point monthly input")
}
par_mismatch <- is.finite(df0_raw$PAR) & is.finite(df0_raw$PAR_corrected_SWin) &
  abs(df0_raw$PAR - df0_raw$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin in 48-point monthly input")
}
if (!swin_col %in% names(df0_raw)) {
  stop("missing SW_IN (incoming shortwave) column required for the Rg<10 night criterion")
}

par_day_threshold  <- 100     # daytime: PAR > 100 (= Rg > 44 W m^-2) together with Ta > 17
rg_night_threshold <- 10      # nighttime: Rg (SW_IN) < 10 W m^-2 (matches the Section 2.6 partitioning cut)

# ───────────────────────────────────────────────────────────────────────────────
# bins
ta_breaks_day  <- seq(5, 39, by = 2)
ta_centers_day <- seq(5, 37, by = 2)

ta_breaks_night  <- seq(5, 39, by = 2)
ta_centers_night <- head(ta_breaks_night, -1) + diff(ta_breaks_night) / 2

ts_breaks  <- seq(10, 40, by = 2)
ts_centers <- head(ts_breaks, -1) + diff(ts_breaks) / 2

# swc_1_1_1 is integrated soil-water storage over 0-100 cm, in cm (~13-64); bin on the cm scale
swc_breaks  <- seq(0, 100, by = 5)
swc_centers <- head(swc_breaks, -1) + diff(swc_breaks) / 2

# dynamic vpd bins from full dataset, with kPa-friendly step
vpd_vals <- as.numeric(df0_raw[[vpd_col]])
vpd_vals <- vpd_vals[is.finite(vpd_vals)]

if (length(vpd_vals) < 2) {
  stop("vpd column has <2 finite values. check vpd_col name and input file.")
}

vpd_min <- min(vpd_vals)
vpd_max <- max(vpd_vals)
vpd_rng <- vpd_max - vpd_min

vpd_step <- dplyr::case_when(
  vpd_rng <= 1.5 ~ 0.2,
  vpd_rng <= 3.0 ~ 0.5,
  vpd_rng <= 8.0 ~ 1.0,
  TRUE           ~ 2.0
)

vpd_lo <- floor(vpd_min / vpd_step) * vpd_step
vpd_hi <- ceiling(vpd_max / vpd_step) * vpd_step

vpd_breaks  <- seq(vpd_lo, vpd_hi, by = vpd_step)
vpd_centers <- head(vpd_breaks, -1) + diff(vpd_breaks) / 2

# base theme used in NEE script
base_theme <- theme_bw(base_size = axis_text_size) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = legend_text_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    plot.tag = element_text(size = plot_tag_size, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

# ───────────────────────────────────────────────────────────────────────────────
# 1) parse datetime from ym + hour, define Season
df0 <- df0_raw %>%
  mutate(
    ym = as.character(ym),
    season = tolower(as.character(season)),
    hour = as.character(hour),
    DateTime = ymd_hms(paste0(ym, "-01 ", hour), tz = tz_local),
    month = month(DateTime),
    Season = factor(season, levels = c("dry", "wet")),
    # swc_1_1_1 is manaus-calibrated integrated soil-water storage over 0-100 cm, reported in cm as-is
    "{swc_col}" := as.numeric(.data[[swc_col]])
  )

cat("\nseason counts BEFORE filters:\n")
print(table(df0$Season, useNA = "ifany"))

# 2) filtered subsets (shared)
df_day_ta <- df0 %>%
  filter(
    is.finite(.data[[par_col]]), .data[[par_col]] > par_day_threshold,
    is.finite(.data[[ta_col]]),  .data[[ta_col]]  > 17
  )

# Daytime uses measured NEE; nighttime (df_night, below) keeps gap-filled NEE_ok.
# Measured nighttime NEE is corrupted by storage/advection artifacts, so respiration
# panels (SWC/Ts nighttime) rest on the u*-filtered gap-filled series instead.
if (identical(nee_source, "measured") && "NEE_meas" %in% names(df_day_ta)) {
  df_day_ta$NEE_ok <- df_day_ta$NEE_meas
}

df_night <- df0 %>%
  filter(
    is.finite(.data[[swin_col]]), .data[[swin_col]] < rg_night_threshold
  )

cat(sprintf("\nseason counts DAY (par>%s & ta>17):\n", par_day_threshold))
print(table(df_day_ta$Season, useNA = "ifany"))

cat(sprintf("\nseason counts NIGHT (Rg<%s W m-2):\n", rg_night_threshold))
print(table(df_night$Season, useNA = "ifany"))

# ───────────────────────────────────────────────────────────────────────────────
# helper: build binned dataframe for a given var + response
build_binned <- function(df, x_col, y_col, breaks, centers, mid_name) {
  df %>%
    transmute(
      DateTime = DateTime,
      X = as.numeric(.data[[x_col]]),
      Y = as.numeric(.data[[y_col]]),
      Season = Season
    ) %>%
    filter(is.finite(X), is.finite(Y)) %>%
    mutate(
      X_bin = cut(X, breaks = breaks, right = FALSE, include.lowest = TRUE),
      X_mid = centers[pmax(1, pmin(length(centers), as.integer(X_bin)))]
    ) %>%
    filter(!is.na(X_bin)) %>%
    rename(!!mid_name := X_mid) %>%
    select(DateTime, Season, Y, all_of(mid_name))
}

# helper: true top-center label (panel coordinates via npc)
add_daynight_label <- function(label) {
  annotation_custom(
    grob = grid::textGrob(
      label,
      x = unit(0.5, "npc"),
      y = unit(0.92, "npc"),
      hjust = 0.5,
      vjust = 1,
      gp = grid::gpar(fontsize = annotation_text_size)
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
}

# helper: plot binned boxplot (legend top-left by default)
# note: if ylim is provided, uses coord_cartesian so outliers/stats are preserved
make_boxplot <- function(df_binned, mid_col, y_label, x_label,
                         legend_pos = c(0.02, 0.98),
                         ylim = NULL) {
  
  p <- ggplot(df_binned, aes(x = factor(.data[[mid_col]]), y = Y, fill = Season)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_boxplot(width = 0.6, outlier.size = 0.7, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = season_box_fills) +
    labs(x = x_label, y = y_label) +
    base_theme +
    theme(
      legend.position = legend_pos,
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.key.size = unit(0.9, "lines")
    )
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}

# ───────────────────────────────────────────────────────────────────────────────
# response config: label, units, output prefixes, composite style, y-cap (optional)
get_response_spec <- function(response_col) {
  if (response_col == "GEP") {
    list(
      prefix = "Fig6_Tambopata_GEP",
      y_grob_label = "GEP (µmol m⁻² s⁻¹)",
      ylab = expression(GEP*"("*mu*mol~m^{-2}~s^{-1}*")"),
      composite_type = "day_only",
      ycap = NULL
    )
  } else if (response_col == "NEE_ok") {
    list(
      prefix = "Fig6_Tambopata_NEE",
      y_grob_label = "NEE (µmol m⁻² s⁻¹)",
      ylab = expression(NEE~"("*mu*mol~m^{-2}~s^{-1}*")"),
      composite_type = "day_night",
      ycap = NULL
    )
  } else if (response_col == "LE") {
    list(
      prefix = "Fig6_Tambopata_LE",
      y_grob_label = "LE (W m⁻²)",
      ylab = expression(LE~"(W"~m^{-2}*")"),
      composite_type = "day_only",
      ycap = NULL
    )
  } else if (response_col == "Gs_mps") {
    list(
      prefix = "Fig6_Tambopata_Gs",
      y_grob_label = "Gs (m s⁻¹)",
      ylab = expression(italic(G)[plain(s)]~"(m"~s^{-1}*")"),
      composite_type = "day_only",
      ycap = c(-0.2, 0.2)
    )
  } else {
    stop("unknown response_col: ", response_col)
  }
}

# ───────────────────────────────────────────────────────────────────────────────
# runner: create + save all panels for one response variable
run_set <- function(response_col) {
  spec <- get_response_spec(response_col)
  
  # build binned dataframes
  df_ta_day  <- build_binned(df_day_ta, ta_col,  response_col, ta_breaks_day,   ta_centers_day,   "Ta_mid")
  df_vpd_day <- build_binned(df_day_ta, vpd_col, response_col, vpd_breaks,      vpd_centers,      "VPD_mid")
  df_vpd_ngt <- build_binned(df_night,  vpd_col, response_col, vpd_breaks,      vpd_centers,      "VPD_mid")
  df_swc_day <- build_binned(df_day_ta, swc_col, response_col, swc_breaks,      swc_centers,      "SWC_mid")
  df_swc_ngt <- build_binned(df_night,  swc_col, response_col, swc_breaks,      swc_centers,      "SWC_mid")
  df_ts_day  <- build_binned(df_day_ta, ts_col,  response_col, ts_breaks,       ts_centers,       "Ts_mid")
  df_ts_ngt  <- build_binned(df_night,  ts_col,  response_col, ts_breaks,       ts_centers,       "Ts_mid")
  df_ta_ngt  <- build_binned(df_night,  ta_col,  response_col, ta_breaks_night, ta_centers_night, "Ta_mid")
  
  # plots (single panels)
  p_ta_day  <- make_boxplot(df_ta_day,  "Ta_mid",  spec$ylab, ta_axis_label, ylim = spec$ycap) +
    add_daynight_label("Daytime")
  
  p_vpd_day <- make_boxplot(df_vpd_day, "VPD_mid", spec$ylab, "VPD (kPa)", ylim = spec$ycap) +
    add_daynight_label("Daytime")
  
  p_vpd_ngt <- make_boxplot(df_vpd_ngt, "VPD_mid", spec$ylab, "VPD (kPa)", ylim = spec$ycap) +
    add_daynight_label("Nighttime")
  
  p_swc_day <- make_boxplot(df_swc_day, "SWC_mid", spec$ylab, "SWC (cm)", ylim = spec$ycap) +
    add_daynight_label("Daytime")
  
  p_swc_ngt <- make_boxplot(df_swc_ngt, "SWC_mid", spec$ylab, "SWC (cm)", ylim = spec$ycap) +
    add_daynight_label("Nighttime")
  
  p_ts_day  <- make_boxplot(df_ts_day,  "Ts_mid",  spec$ylab, ts_axis_label, ylim = spec$ycap) +
    add_daynight_label("Daytime")
  
  p_ts_ngt  <- make_boxplot(df_ts_ngt,  "Ts_mid",  spec$ylab, ts_axis_label, ylim = spec$ycap) +
    add_daynight_label("Nighttime")
  
  p_ta_ngt  <- make_boxplot(df_ta_ngt,  "Ta_mid",  spec$ylab, ta_axis_label, ylim = spec$ycap) +
    add_daynight_label("Nighttime")
  
  # save 8 single panels
  prefix <- spec$prefix
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_Tair_bins_daytime_Ta17_48ppmonth.png")),
         p_ta_day, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_VPD_bins_daytime_Ta17_48ppmonth.png")),
         p_vpd_day, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_VPD_bins_nighttime_Rglt10_48ppmonth.png")),
         p_vpd_ngt, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_SWC_bins_daytime_Ta17_48ppmonth.png")),
         p_swc_day, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_SWC_bins_nighttime_Rglt10_48ppmonth.png")),
         p_swc_ngt, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_Ts_bins_daytime_Ta17_48ppmonth.png")),
         p_ts_day, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_Ts_bins_nighttime_Rglt10_48ppmonth.png")),
         p_ts_ngt, width = 8.5, height = 5.0, dpi = 300)
  
  ggsave(file.path(paper_figures_path, paste0(prefix, "_vs_Tair_bins_nighttime_Rglt10_48ppmonth.png")),
         p_ta_ngt, width = 8.5, height = 5.0, dpi = 300)
  
  # composite figure
  if (spec$composite_type == "day_only") {
    p_a <- p_ta_day + theme(legend.position = "none") + labs(tag = "(a)", y = NULL)
    p_b <- p_vpd_day + theme(legend.position = "none") + labs(tag = "(b)", y = NULL)
    p_c <- p_swc_day + theme(legend.position = "none") + labs(tag = "(c)", y = NULL)
    p_d <- p_ts_day  + theme(legend.position = "none") + labs(tag = "(d)", y = NULL)
    
    panel_grid <- (p_a | p_b) / (p_c | p_d)
    
    y_grob <- wrap_elements(
      full = textGrob(
        label = spec$y_grob_label,
        rot = 90,
        gp = gpar(fontsize = shared_axis_title_size)
      )
    )
    
    main_fig <- (y_grob | panel_grid) + plot_layout(widths = c(0.06, 1))
    print(main_fig)
    
    out_fp <- file.path(paper_figures_path, paste0(prefix, "_MainControls_DayOnly_48ppmonth.png"))
    ggsave(out_fp, main_fig, width = 11, height = 9, dpi = 300)
    message("saved: ", out_fp)
    
  } else if (spec$composite_type == "day_night") {
    
    # ---------------------------------------------------------
    # EXISTING MAIN FIGURE 6 CODE
    # ---------------------------------------------------------
    p_a <- p_ta_day  + theme(legend.position = "none") + labs(tag = "(a)")
    p_b <- p_vpd_day + theme(legend.position = "none") + labs(tag = "(b)")
    p_c <- p_swc_ngt + theme(legend.position = "none") + labs(tag = "(c)")
    p_d <- p_ts_ngt  + labs(tag = "(d)")
    
    main_fig <- (p_a | p_b) / (p_c | p_d) +
      plot_layout(axis_titles = "collect_y") +
      plot_annotation(theme = theme(legend.position = c(0.02, 0.98)))
    
    main_fig <- main_fig & theme(
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.key.size = unit(0.9, "lines")
    )
    
    print(main_fig)
    
    out_fp <- file.path(paper_figures_path, "Fig6_Tambopata_MainControls_DayNight_48ppmonth.png")
    ggsave(out_fp, main_fig, width = 11, height = 9, dpi = 300)
    message("saved: ", out_fp)
    
    # ---------------------------------------------------------
    # NEW SUPPLEMENTARY FIGURE CODE
    # ---------------------------------------------------------
    # Assign the requested nighttime and daytime panels
    p_a_supp <- p_ta_ngt  + theme(legend.position = "none") + labs(tag = "(a)")
    p_b_supp <- p_vpd_ngt + theme(legend.position = "none") + labs(tag = "(b)")
    p_c_supp <- p_swc_day + theme(legend.position = "none") + labs(tag = "(c)")
    p_d_supp <- p_ts_day  + labs(tag = "(d)") # Legend stays here
    
    # Combine using patchwork
    supp_fig <- (p_a_supp | p_b_supp) / (p_c_supp | p_d_supp) +
      plot_layout(axis_titles = "collect_y") +
      plot_annotation(theme = theme(legend.position = c(0.02, 0.98)))
    
    # Apply identical theme adjustments to match main figure
    supp_fig <- supp_fig & theme(
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.key.size = unit(0.9, "lines")
    )
    
    print(supp_fig)
    
    # Save the supplementary figure using dynamic prefix
    out_fp_supp <- file.path(paper_figures_path, paste0(prefix, "_SuppControls_NightDay_48ppmonth.png"))
    ggsave(out_fp_supp, supp_fig, width = 11, height = 9, dpi = 300)
    message("saved: ", out_fp_supp)
    
  } else {
    stop("unknown composite_type: ", spec$composite_type)
  }
  
} # <-- THIS correctly closes the run_set() function

# ───────────────────────────────────────────────────────────────────────────────
# run all figure sets
run_set("GEP")
run_set("NEE_ok")
run_set("LE")
run_set("Gs_mps")
