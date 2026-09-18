#!/usr/bin/env Rscript

# Supplementary GEP-vs-environment figure using current daily integrated GEP
# and the same leave-one-out slope-stability rule as the parameter figures.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridisLite)
  library(grid)
  library(readr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)

flux_csv <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
env_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
plot_csv <- file.path(output_path, "GEP_mean_env_by_year_season_ENSO.csv")

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv, "\nRun code/01_prepare_flux_timeseries.R first.")
}
if (!file.exists(env_csv)) {
  stop("yearly environmental table not found: ", env_csv, "\nRun code/03_meteorology_summary.R first.")
}

flux_required <- c("Date", "year", "season", "ENSO", "GEP")
env_required <- c("Year", "Season", "ENSO", "Ta", "VPD", "SWC", "PAR", "Ts")

flux <- readr::read_csv(flux_csv, show_col_types = FALSE)
env <- readr::read_csv(env_csv, show_col_types = FALSE)

missing_flux <- setdiff(flux_required, names(flux))
missing_env <- setdiff(env_required, names(env))
if (length(missing_flux)) {
  stop("missing required columns in enriched flux table: ", paste(missing_flux, collapse = ", "))
}
if (length(missing_env)) {
  stop("missing required columns in yearly environmental table: ", paste(missing_env, collapse = ", "))
}

normalise_season <- function(x) {
  dplyr::case_when(
    tolower(as.character(x)) == "wet" ~ "Wet",
    tolower(as.character(x)) == "dry" ~ "Dry",
    TRUE ~ as.character(x)
  )
}

normalise_enso <- function(x) {
  dplyr::case_when(
    grepl("El", x, ignore.case = TRUE) ~ "El Nino",
    grepl("La", x, ignore.case = TRUE) ~ "La Nina",
    grepl("Neutral", x, ignore.case = TRUE) | tolower(as.character(x)) == "neutral" ~ "Neutral",
    TRUE ~ as.character(x)
  )
}

gC_per_umolCO2 <- 12.0107e-6
sec_per_halfhour <- 1800

daily_gep <- flux %>%
  mutate(
    Date = as.Date(.data$Date),
    Year = as.integer(.data$year),
    Season = normalise_season(.data$season),
    ENSO = normalise_enso(.data$ENSO),
    GEP = as.numeric(.data$GEP)
  ) %>%
  filter(!is.na(Date), is.finite(GEP), !is.na(Year), !is.na(Season), !is.na(ENSO)) %>%
  group_by(Year, Season, ENSO, Date) %>%
  summarise(
    GEP_gC_m2_day = sum(GEP * sec_per_halfhour * gC_per_umolCO2, na.rm = TRUE),
    n_halfhours = dplyr::n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(GEP_gC_m2_day))

seasonal_gep <- daily_gep %>%
  group_by(Year, Season, ENSO) %>%
  summarise(
    GEP_gC_m2_day = mean(GEP_gC_m2_day, na.rm = TRUE),
    GEP_sd_gC_m2_day = stats::sd(GEP_gC_m2_day, na.rm = TRUE),
    n_days = dplyr::n(),
    n_halfhours = sum(n_halfhours, na.rm = TRUE),
    .groups = "drop"
  )

env_summary <- env %>%
  mutate(
    Year = as.integer(.data$Year),
    Season = normalise_season(.data$Season),
    ENSO = normalise_enso(.data$ENSO)
  ) %>%
  select(Year, Season, ENSO, Ta, VPD, SWC, PAR, Ts)

plot_data <- env_summary %>%
  left_join(seasonal_gep, by = c("Year", "Season", "ENSO")) %>%
  arrange(Year, Season, ENSO)

readr::write_csv(plot_data, plot_csv)
message("wrote plotting table: ", plot_csv)

plot_data <- plot_data %>%
  mutate(
    Season = factor(Season, levels = c("Wet", "Dry")),
    ENSO = factor(ENSO, levels = c("El Nino", "La Nina", "Neutral"))
  )

season_fills <- c("Wet" = "gray70", "Dry" = "white")
enso_shapes <- c("El Nino" = 21, "La Nina" = 24, "Neutral" = 22)
xvars_order <- c("Ta", "VPD", "SWC", "PAR", "Ts")
xvars_available <- intersect(xvars_order, names(plot_data))
years_chr <- sort(unique(as.character(plot_data$Year)))
year_colors <- setNames(viridisLite::turbo(length(years_chr)), years_chr)
years_with_data <- plot_data %>%
  tidyr::drop_na(GEP_gC_m2_day) %>%
  pull(Year) %>%
  unique() %>%
  as.character() %>%
  sort()

xlab_map <- list(
  Ta = expression(italic(T)[plain(a)] ~ "(°C)"),
  VPD = "VPD (kPa)",
  SWC = "SWC (cm)",
  PAR = "PAR (µmol m⁻² s⁻¹)",
  Ts = expression(italic(T)[plain(s)] ~ "(°C)")
)

ylab_gep <- expression(GEP ~ "(" * g ~ C ~ m^{-2} ~ d^{-1} * ")")
gep_limits <- range(plot_data$GEP_gC_m2_day, na.rm = TRUE)

lm_label <- function(d, xvar, yvar) {
  if (!nrow(d) || dplyr::n_distinct(d[[xvar]]) < 2 || dplyr::n_distinct(d[[yvar]]) < 2) {
    return(NULL)
  }

  fit <- stats::lm(stats::as.formula(paste(yvar, "~", xvar)), data = d)
  s <- summary(fit)
  sprintf(
    "y = %.3g + %.3g*x\nR2 = %.3f; p = %.3g",
    coef(fit)[1],
    coef(fit)[2],
    s$r.squared,
    coef(s)[2, "Pr(>|t|)"]
  )
}

has_stable_loo_slope <- function(d, xvar, yvar) {
  if (nrow(d) < 4 || dplyr::n_distinct(d[[xvar]]) < 2 || dplyr::n_distinct(d[[yvar]]) < 2) {
    return(FALSE)
  }

  fit <- tryCatch(
    stats::lm(stats::as.formula(paste(yvar, "~", xvar)), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(FALSE)

  all_slope <- unname(stats::coef(fit)[[xvar]])
  if (!is.finite(all_slope) || all_slope == 0) return(FALSE)

  loo_slopes <- vapply(seq_len(nrow(d)), function(i) {
    dd <- d[-i, , drop = FALSE]
    if (nrow(dd) < 3 || dplyr::n_distinct(dd[[xvar]]) < 2 || dplyr::n_distinct(dd[[yvar]]) < 2) {
      return(NA_real_)
    }

    loo_fit <- tryCatch(
      stats::lm(stats::as.formula(paste(yvar, "~", xvar)), data = dd),
      error = function(e) NULL
    )
    if (is.null(loo_fit)) return(NA_real_)

    unname(stats::coef(loo_fit)[[xvar]])
  }, numeric(1))

  all(is.finite(loo_slopes)) && !any(sign(loo_slopes) != sign(all_slope))
}

make_gep_plot <- function(df, xvar, show_legend = FALSE, stable_lm_only = TRUE) {
  yvar <- "GEP_gC_m2_day"
  d <- df %>%
    tidyr::drop_na(all_of(c(xvar, yvar, "Year", "ENSO", "Season"))) %>%
    mutate(Year = factor(Year, levels = years_chr))

  if (!nrow(d)) return(ggplot() + theme_void())

  draw_lm <- !isTRUE(stable_lm_only) || has_stable_loo_slope(d, xvar, yvar)

  x_range <- range(d[[xvar]], na.rm = TRUE)
  jitter_width <- if (all(is.finite(x_range)) && diff(x_range) > 0) 0.01 * diff(x_range) else 0

  p <- ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(
      aes(color = Year, fill = Season, shape = ENSO),
      size = 3.2,
      stroke = 1.6,
      alpha = 0.95,
      position = position_jitter(width = jitter_width, height = 0, seed = 123)
    ) +
    scale_color_manual(
      name = "Year",
      values = year_colors,
      breaks = years_with_data,
      limits = years_with_data,
      drop = FALSE
    ) +
    scale_fill_manual(name = "Season", values = season_fills) +
    scale_shape_manual(name = "ENSO", values = enso_shapes) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(order = 2, override.aes = list(fill = "white", color = "black")),
      color = guide_legend(
        order = 3,
        ncol = 2,
        override.aes = list(
          shape = 21,
          fill = unname(year_colors[years_with_data]),
          size = 3.2,
          stroke = 1.2
        )
      )
    ) +
    labs(x = xlab_map[[xvar]], y = NULL) +
    theme_minimal() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_blank(),
      axis.title.y = element_blank()
    )

  if (draw_lm) {
    p <- p +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", formula = y ~ x)
  }

  if (all(is.finite(gep_limits))) {
    p <- p + coord_cartesian(ylim = gep_limits)
  }

  lbl <- if (draw_lm) lm_label(d, xvar, yvar) else NULL
  if (!is.null(lbl)) {
    yr <- gep_limits
    if (all(is.finite(x_range)) && diff(x_range) > 0 && all(is.finite(yr)) && diff(yr) > 0) {
      p <- p + annotate(
        "label",
        x = x_range[1] + 0.02 * diff(x_range),
        y = yr[2] - 0.10 * diff(yr),
        label = lbl,
        size = 3,
        label.padding = unit(0.2, "lines"),
        alpha = 0.9,
        hjust = 0
      )
    }
  }

  p
}

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx)) g$grobs[[idx[1]]] else NULL
}

message("Generating GEP-vs-environment supplement figure...")

plot_list <- list()
for (xvar in xvars_available) {
  plot_list[[xvar]] <- make_gep_plot(plot_data, xvar, stable_lm_only = TRUE)
}

legend_dummy <- make_gep_plot(plot_data, xvars_available[1], show_legend = TRUE, stable_lm_only = TRUE)
legend_grob <- extract_legend_grob(legend_dummy)
plot_list[["legend"]] <- if (is.null(legend_grob)) plot_spacer() else wrap_elements(legend_grob)

supp_grid <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

y_axis_grob <- grid::textGrob(ylab_gep, rot = 90, gp = grid::gpar(fontsize = 14))
final_fig <- wrap_elements(y_axis_grob) | supp_grid
final_fig <- final_fig + plot_layout(widths = c(0.04, 1))

canonical_path <- file.path(paper_figures_path, "GEP_vs_environment_daily_yearcolor.png")
legacy_path <- file.path(paper_figures_path, "FigS1_yearcolor_GEP_vs_env.png")

ggplot2::ggsave(canonical_path, plot = final_fig, width = 14, height = 9, dpi = 300, bg = "white")
invisible(file.copy(canonical_path, legacy_path, overwrite = TRUE))

message("wrote: ", canonical_path)
message("wrote legacy alias: ", legacy_path)

line_status <- lapply(xvars_available, function(xvar) {
  d <- plot_data %>%
    tidyr::drop_na(all_of(c(xvar, "GEP_gC_m2_day", "Year", "ENSO", "Season")))
  tibble::tibble(
    xvar = xvar,
    n = nrow(d),
    line_shown = has_stable_loo_slope(d, xvar, "GEP_gC_m2_day")
  )
}) %>%
  dplyr::bind_rows()

message("leave-one-out line status:")
print(line_status)
