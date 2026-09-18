# ───────────────────────────────────────────────────────────────────────────────
# nighttime NEE_ok vs soil moisture and soil temperature (48 points per month)
# - uses ym + hour (no DateTime in this aggregated file)
# - defines nighttime hours (20:00–05:59)
# - filters to rows with needed vars (NEE_ok, SWC_1_1_1, TS_3, VPD_kPa)
# - fits simple linear models (NEE_ok ~ SWC, NEE_ok ~ Ts)
# - plots NEE_ok vs SWC and NEE_ok vs Ts, colored by VPD in grayscale
# - legend inside the NEE_ok vs Ts panel (top-left), smaller
# ───────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)
library(grid)

# ───────────────────────────────────────────────────────────────────────────────
# paths
# ───────────────────────────────────────────────────────────────────────────────
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# load data (CSV)
# ───────────────────────────────────────────────────────────────────────────────
# nee_source: "measured" = non-gap-filled raw NEE (primary; recovers the 15 RP-empty months),
#   "gapfilled" = NEE_ok (= NEE_f), sensitivity.
nee_source <- "measured"
dat <- read.csv(
  file.path(input_folder,
            if (identical(nee_source, "measured")) "tambopata_48points_per_month_measured.csv"
            else "tambopata_48points_per_month.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
# internal response column NEE_ok holds the non-gap-filled measured NEE in measured mode
if (identical(nee_source, "measured")) dat$NEE_ok <- dat$NEE_meas

# fix common mojibake in ENSO strings (e.g., "La NiÃ±a" -> "La Niña")
if ("ENSO" %in% names(dat)) {
  dat$ENSO <- iconv(dat$ENSO, from = "UTF-8", to = "UTF-8")
  dat$ENSO <- gsub("NiÃ±a", "Niña", dat$ENSO, fixed = TRUE)
}

# ───────────────────────────────────────────────────────────────────────────────
# prepare time fields (no DateTime in this file)
# ───────────────────────────────────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    year = as.integer(substr(ym, 1, 4)),
    # parse hour like "3:30:00" into seconds, then to decimal hours
    hour_decimal = as.numeric(hms(hour)) / 3600
  )

# nighttime: 20:00–05:59
dat_night <- dat %>%
  filter(hour_decimal >= 20 | hour_decimal < 6) %>%
  filter(
    is.finite(NEE_ok),
    is.finite(SWC_1_1_1),
    is.finite(TS_3),
    is.finite(VPD_kPa)
  ) %>%
  mutate(
    VPD_kPa = as.numeric(VPD_kPa)
  ) %>%
  filter(is.finite(VPD_kPa))

# ───────────────────────────────────────────────────────────────────────────────
# split day/night using PAR
# ───────────────────────────────────────────────────────────────────────────────

# use corrected PAR as the canonical light variable
if (!"PAR_corrected_SWin" %in% names(dat) && "PPFD_IN_1_1_1" %in% names(dat)) {
  dat$PAR_corrected_SWin <- suppressWarnings(as.numeric(dat$PPFD_IN_1_1_1))
}
if (!"PAR" %in% names(dat) && "PAR_corrected_SWin" %in% names(dat)) {
  dat$PAR <- suppressWarnings(as.numeric(dat$PAR_corrected_SWin))
}
if (!"PAR" %in% names(dat)) {
  stop("missing MATLAB-corrected PAR column")
}
par_mismatch <- is.finite(dat$PAR) & is.finite(dat$PAR_corrected_SWin) &
  abs(dat$PAR - dat$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin")
}

par_day_threshold <- 20

dat_common <- dat %>%
  mutate(
    NEE_ok    = as.numeric(NEE_ok),
    # swc_1_1_1 is manaus-calibrated integrated soil-water storage over 0-100 cm, reported in cm as-is
    SWC_1_1_1 = as.numeric(SWC_1_1_1),
    TS_3      = as.numeric(TS_3),
    TA_1_1_1  = as.numeric(TA_1_1_1),
    VPD_kPa   = as.numeric(VPD_kPa),
    PAR       = as.numeric(PAR)
  ) %>%
  filter(
    is.finite(NEE_ok),
    is.finite(SWC_1_1_1),
    is.finite(TS_3),
    is.finite(TA_1_1_1),
    is.finite(VPD_kPa),
    is.finite(PAR)
  )

dat_day <- dat_common %>%
  filter(PAR >= par_day_threshold)

dat_night <- dat_common %>%
  filter(PAR < par_day_threshold)

lm_label_stats <- function(df, xvar) {
  if (
    nrow(df) < 3 ||
      dplyr::n_distinct(df[[xvar]]) < 2 ||
      dplyr::n_distinct(df$NEE_ok) < 2
  ) {
    return(data.frame(r2 = NA_real_, p_value = NA_real_))
  }

  fit <- lm(stats::as.formula(paste("NEE_ok ~", xvar)), data = df)
  fit_summary <- summary(fit)

  data.frame(
    r2 = fit_summary$r.squared,
    p_value = coef(fit_summary)[2, "Pr(>|t|)"]
  )
}

format_lm_p <- function(p_value) {
  ifelse(
    is.na(p_value),
    "NA",
    ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )
}

nee_x_labels <- list(
  SWC_1_1_1 = expression(SWC~"(cm)"),
  TS_3 = expression(italic(T)[plain(s)]~"(" * degree*C * ")"),
  VPD_kPa = "VPD (kPa)",
  TA_1_1_1 = expression(italic(T)[plain(a)]~"(" * degree*C * ")")
)

decile_units <- c(
  SWC_1_1_1 = "cm",
  TS_3 = "\u00b0C",
  VPD_kPa = "kPa",
  TA_1_1_1 = "\u00b0C"
)

make_nee_decile_plot <- function(df, xvar, decile_var, title_label) {
  decile_col <- "plot_decile"
  plot_df <- df %>%
    filter(
      is.finite(NEE_ok),
      is.finite(.data[[xvar]]),
      is.finite(.data[[decile_var]])
    ) %>%
    mutate("{decile_col}" := dplyr::ntile(.data[[decile_var]], 10L))

  if (!nrow(plot_df)) {
    stop("No finite rows available for ", title_label)
  }

  decile_labels <- plot_df %>%
    group_by(.data[[decile_col]]) %>%
    summarise(
      decile_min = min(.data[[decile_var]], na.rm = TRUE),
      decile_max = max(.data[[decile_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      facet_label = sprintf(
        "%d-%d\n(%.2f-%.2f %s)",
        (.data[[decile_col]] - 1L) * 10L,
        .data[[decile_col]] * 10L,
        decile_min,
        decile_max,
        decile_units[[decile_var]]
      )
    )

  plot_df <- plot_df %>%
    left_join(decile_labels, by = decile_col) %>%
    mutate(facet_label = factor(facet_label, levels = decile_labels$facet_label))

  label_df <- plot_df %>%
    group_by(facet_label) %>%
    group_modify(~ lm_label_stats(.x, xvar = xvar)) %>%
    ungroup() %>%
    mutate(
      label = sprintf("R\u00b2=%.2f\np=%s", r2, format_lm_p(p_value))
    )

  lm_df <- plot_df %>%
    group_by(facet_label) %>%
    filter(
      n() >= 3,
      dplyr::n_distinct(.data[[xvar]]) >= 2,
      dplyr::n_distinct(NEE_ok) >= 2
    ) %>%
    ungroup()

  ggplot(plot_df, aes(x = .data[[xvar]], y = NEE_ok)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.35) +
    geom_point(color = "grey40", alpha = 0.45, size = 1.6) +
    geom_smooth(
      data = lm_df,
      method = "lm",
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_text(
      data = label_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.15,
      size = 4.2
    ) +
    facet_wrap(~ facet_label, ncol = 5) +
    labs(
      title = title_label,
      x = nee_x_labels[[xvar]],
      y = expression(NEE~"(" * mu * "mol CO"[2]~m^-2~s^-1 * ")")
    ) +
    # Kept deliberately identical to base_theme_facets in
    # code/10_surface_conductance.R:625-635. This figure (supplement Fig. S16) and
    # the Gs decile figure (S17) are the same analysis on two response variables
    # and sit next to each other in the supplement, so they share one style. The
    # ggsave geometry below mirrors code/10:93-95 for the same reason -- change
    # both files together or they drift apart again.
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 13, lineheight = 0.95),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      plot.title = element_text(size = 18)
    )
}

nee_decile_specs <- data.frame(
  xvar = c(
    "SWC_1_1_1", "SWC_1_1_1",
    "SWC_1_1_1", "SWC_1_1_1",
    "TS_3", "TS_3",
    "TS_3", "TS_3",
    "VPD_kPa", "VPD_kPa",
    "TA_1_1_1", "TA_1_1_1"
  ),
  decile_var = c(
    "VPD_kPa", "VPD_kPa",
    "TA_1_1_1", "TA_1_1_1",
    "TA_1_1_1", "TA_1_1_1",
    "VPD_kPa", "VPD_kPa",
    "TA_1_1_1", "TA_1_1_1",
    "VPD_kPa", "VPD_kPa"
  ),
  period = rep(c("day", "night"), 6),
  x_label_text = c(
    "SWC", "SWC",
    "SWC", "SWC",
    "Ts", "Ts",
    "Ts", "Ts",
    "VPD", "VPD",
    "Ta", "Ta"
  ),
  decile_label_text = c(
    "VPD", "VPD",
    "TA", "TA",
    "TA", "TA",
    "VPD", "VPD",
    "TA", "TA",
    "VPD", "VPD"
  ),
  out_path = c(
    file.path(paper_figures_path, "FigX_day_NEEok_vs_SWC_VPDdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_SWC_VPDdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_day_NEEok_vs_SWC_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_SWC_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_day_NEEok_vs_Ts_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_Ts_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_day_NEEok_vs_Ts_VPDdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_Ts_VPDdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_day_NEEok_vs_VPD_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_VPD_TAdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_day_NEEok_vs_Ta_VPDdeciles_LM_48points_per_month.png"),
    file.path(paper_figures_path, "FigX_night_NEEok_vs_Ta_VPDdeciles_LM_48points_per_month.png")
  ),
  stringsAsFactors = FALSE
)

# Supplement reduction (2026-09-08). These twelve decile panels were figures
# S48-S59 and existed to show that Ta and VPD cannot be separated in this record.
# Twelve figures is far more than that claim needs, so the supplement keeps only
# the one the manuscript sentence actually rests on -- daytime NEE against VPD,
# conditioned on Ta deciles -- and points readers to the Zenodo deposit for the
# rest. The full spec table is retained above so the others can be regenerated by
# commenting out this filter.
nee_decile_keep <- with(nee_decile_specs,
                        period == "day" & xvar == "VPD_kPa" & decile_var == "TA_1_1_1")
message(sprintf("decile panels: keeping %d of %d (supplement reduction)",
                sum(nee_decile_keep), nrow(nee_decile_specs)))
nee_decile_specs <- nee_decile_specs[nee_decile_keep, , drop = FALSE]
stopifnot(nrow(nee_decile_specs) == 1L)

for (i in seq_len(nrow(nee_decile_specs))) {
  spec <- nee_decile_specs[i, ]
  period_df <- if (identical(spec$period, "day")) dat_day else dat_night
  period_title <- if (identical(spec$period, "day")) {
    sprintf("Day (PAR \u2265 %s)", par_day_threshold)
  } else {
    sprintf("Night (PAR < %s)", par_day_threshold)
  }
  plot_title <- sprintf(
    "%s: NEE vs %s (panels = %s deciles)",
    period_title,
    spec$x_label_text,
    spec$decile_label_text
  )
  plot_obj <- make_nee_decile_plot(
    df = period_df,
    xvar = spec$xvar,
    decile_var = spec$decile_var,
    title_label = plot_title
  )
  ggsave(
    filename = spec$out_path,
    plot = plot_obj,
    width = 230,
    height = 170,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

# ───────────────────────────────────────────────────────────────────────────────
# optional statistical summaries
# ───────────────────────────────────────────────────────────────────────────────

cat("\n──────── day NEE_ok ~ SWC_1_1_1 ────────\n")
print(summary(lm(NEE_ok ~ SWC_1_1_1, data = dat_day)))

cat("\n──────── day NEE_ok ~ TS_3 ─────────────\n")
print(summary(lm(NEE_ok ~ TS_3, data = dat_day)))

cat("\n──────── night NEE_ok ~ SWC_1_1_1 ────────\n")
print(summary(lm(NEE_ok ~ SWC_1_1_1, data = dat_night)))

cat("\n──────── night NEE_ok ~ TS_3 ─────────────\n")
print(summary(lm(NEE_ok ~ TS_3, data = dat_night)))

# ───────────────────────────────────────────────────────────────────────────────
# Figure 8 removed (2026-07-07)
# ───────────────────────────────────────────────────────────────────────────────
# The half-hourly NEE-vs-SWC/Ts scatter (with LOESS), formerly saved here as
# "Fig8_NEEok_vs_SWC_TS_PARday_PARnight_4panel_48points_per_month.{png,pdf}", was
# dropped from the manuscript following co-author review: it was largely redundant
# with the Fig. 7 boxplots (night SWC/Ts) and the SI day-side boxplot companion
# (day SWC/Ts, from code/09), and its LOESS smooths implied a daytime soil-temperature
# "optimum" that a later reanalysis does not support. The day/night
# NEE-vs-SWC/Ts relationships are retained as the decile-faceted scatter figures and
# lm summaries above. The panel builder and its helpers were removed with the figure.
