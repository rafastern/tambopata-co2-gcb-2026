#!/usr/bin/env Rscript

# Direct parameter-vs-environment highlights for yearly light-response fits.
# This preserves the plotting logic from Inspiration/Fig4_v083_combined_P2000.R
# while using repository-relative paths and analysis-based output names.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridisLite)
  library(grid)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)

yearly_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
if (!file.exists(yearly_csv)) {
  stop(
    "yearly parameter table not found: ", yearly_csv,
    "\nRun code/03_meteorology_summary.R first to create it."
  )
}

summary_table_yearly <- read.csv(yearly_csv, check.names = FALSE)
# guard: one row per year x group (fails loudly if the C1 fan-out duplication recurs)
if (all(c("Year", "Season", "ENSO") %in% names(summary_table_yearly)))
  stopifnot(!any(duplicated(summary_table_yearly[, c("Year", "Season", "ENSO")])))

required_cols <- c("Year", "Season", "ENSO", "Ta", "Ts", "SWC", "VPD", "P2000", "Rd", "phi0", "LCP")
missing_cols <- setdiff(required_cols, names(summary_table_yearly))
if (length(missing_cols)) {
  stop(
    "missing required columns in yearly parameter table: ",
    paste(missing_cols, collapse = ", "),
    "\nFile used: ", yearly_csv
  )
}

summary_table_yearly <- summary_table_yearly %>%
  mutate(
    Season = factor(Season, levels = c("Wet", "Dry")),
    ENSO = case_when(
      grepl("El", ENSO, ignore.case = TRUE) ~ "El Nino",
      grepl("La", ENSO, ignore.case = TRUE) ~ "La Nina",
      grepl("Neutral", ENSO, ignore.case = TRUE) ~ "Neutral",
      TRUE ~ as.character(ENSO)
    ),
    ENSO = factor(ENSO, levels = c("El Nino", "La Nina", "Neutral"))
  )

season_fills <- c("Wet" = "gray70", "Dry" = "white")
enso_shapes <- c("El Nino" = 21, "La Nina" = 24, "Neutral" = 22)
xvars_order <- c("Ta", "VPD", "SWC", "PAR", "Ts")
param_vars <- c("P2000", "Rd", "phi0", "LCP")

ylab_shared_map <- list(
  phi0 = expression(italic(phi)[0] ~ "(" * mu * "mol " * mu * "mol"^{-1} * ")"),
  P2000 = expression(italic(P)[2000] ~ "(" * mu * "mol CO"[2] ~ "m"^{-2} ~ "s"^{-1} * ")"),
  Rd = expression(italic(R)[plain(d)] ~ "(" * mu * "mol CO"[2] ~ "m"^{-2} ~ "s"^{-1} * ")"),
  LCP = expression("LCP (" * mu * "mol m"^{-2} ~ "s"^{-1} * ")")
)

xlab_map <- list(
  Ta = expression(italic(T)[plain(a)] ~ "(°C)"),
  VPD = "VPD (kPa)",
  SWC = "SWC (cm)",
  PAR = "PAR (µmol m⁻² s⁻¹)",
  Ts = expression(italic(T)[plain(s)] ~ "(°C)")
)

xvars_available <- intersect(xvars_order, names(summary_table_yearly))
years_chr <- sort(unique(as.character(summary_table_yearly$Year)))
year_colors <- setNames(viridisLite::turbo(length(years_chr)), years_chr)
years_with_data <- summary_table_yearly %>%
  tidyr::drop_na(any_of(param_vars)) %>%
  pull(Year) %>%
  unique() %>%
  as.character() %>%
  sort()
ylims_map <- setNames(lapply(param_vars, function(yv) {
  range(summary_table_yearly[[yv]], na.rm = TRUE)
}), param_vars)

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

make_xy_plot_year <- function(df, xvar, yvar, y_limits = NULL, show_legend = FALSE,
                              show_lm = TRUE, sensitive_label = NULL,
                              stable_lm_only = FALSE) {
  d <- df %>%
    tidyr::drop_na(all_of(c(xvar, yvar, "Year", "ENSO", "Season"))) %>%
    mutate(Year = factor(Year, levels = years_chr))

  if (!nrow(d)) return(ggplot() + theme_void())

  draw_lm <- isTRUE(show_lm) && (!isTRUE(stable_lm_only) || has_stable_loo_slope(d, xvar, yvar))

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

  if (!is.null(y_limits) && all(is.finite(y_limits))) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  lbl <- if (draw_lm) lm_label(d, xvar, yvar) else sensitive_label
  if (!is.null(lbl)) {
    yr <- if (!is.null(y_limits) && all(is.finite(y_limits))) y_limits else range(d[[yvar]], na.rm = TRUE)
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

save_png_with_alias <- function(plot_obj, canonical_name, legacy_name = NULL, width, height) {
  canonical_path <- file.path(paper_figures_path, canonical_name)
  ggplot2::ggsave(canonical_path, plot = plot_obj, width = width, height = height, dpi = 300, bg = "white")
  message("wrote: ", canonical_path)

  if (!is.null(legacy_name)) {
    legacy_path <- file.path(paper_figures_path, legacy_name)
    file.copy(canonical_path, legacy_path, overwrite = TRUE)
    message("wrote legacy alias: ", legacy_path)
  }

  invisible(canonical_path)
}

message("Generating direct light-response parameter highlight figure...")

p_phi0 <- make_xy_plot_year(summary_table_yearly, "Ta", "phi0", ylims_map[["phi0"]]) +
  labs(y = ylab_shared_map[["phi0"]]) +
  labs(tag = "(a)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_lcp <- make_xy_plot_year(summary_table_yearly, "Ts", "LCP", ylims_map[["LCP"]]) +
  labs(y = ylab_shared_map[["LCP"]]) +
  labs(tag = "(b)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_p2000 <- make_xy_plot_year(summary_table_yearly, "SWC", "P2000", ylims_map[["P2000"]]) +
  labs(y = ylab_shared_map[["P2000"]]) +
  labs(tag = "(c)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_rd <- make_xy_plot_year(summary_table_yearly, "VPD", "Rd", ylims_map[["Rd"]]) +
  labs(y = ylab_shared_map[["Rd"]]) +
  labs(tag = "(d)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )

legend_grob <- extract_legend_grob(
  make_xy_plot_year(summary_table_yearly, "Ta", "phi0", ylims_map[["phi0"]], show_legend = TRUE)
)

combined_grid <- (p_phi0 | p_lcp) / (p_p2000 | p_rd)

combined_fig <- combined_grid |
  if (is.null(legend_grob)) plot_spacer() else wrap_elements(legend_grob)
combined_fig <- combined_fig +
  plot_layout(widths = c(1, 0.25))

fig_parameter_highlight_name <- "light_response_parameter_highlights_P2000.png"
fig_parameter_highlight_legacy_name <- "FigS_Combined_Highlights_P2000.png"

save_png_with_alias(
  combined_fig,
  canonical_name = fig_parameter_highlight_name,
  legacy_name = fig_parameter_highlight_legacy_name,
  width = 10,
  height = 8
)

message("Generating stable-regression Figure 4 variant...")
message("OLS lines are shown only where leave-one-out slope direction was stable; all observations are retained.")

p_phi0_stable <- make_xy_plot_year(
  summary_table_yearly,
  "Ta",
  "phi0",
  ylims_map[["phi0"]],
  show_lm = FALSE,
  sensitive_label = "slope sensitive"
) +
  labs(y = ylab_shared_map[["phi0"]]) +
  labs(tag = "(a)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_lcp_stable <- make_xy_plot_year(
  summary_table_yearly,
  "Ts",
  "LCP",
  ylims_map[["LCP"]],
  show_lm = FALSE,
  sensitive_label = "slope sensitive"
) +
  labs(y = ylab_shared_map[["LCP"]]) +
  labs(tag = "(b)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_p2000_stable <- make_xy_plot_year(summary_table_yearly, "SWC", "P2000", ylims_map[["P2000"]]) +
  labs(y = ylab_shared_map[["P2000"]]) +
  labs(tag = "(c)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )
p_rd_stable <- make_xy_plot_year(summary_table_yearly, "VPD", "Rd", ylims_map[["Rd"]]) +
  labs(y = ylab_shared_map[["Rd"]]) +
  labs(tag = "(d)") +
  theme(
    axis.title.y = element_text(angle = 90, margin = margin(r = 10), size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )

stable_combined_grid <- (p_phi0_stable | p_lcp_stable) / (p_p2000_stable | p_rd_stable)
stable_combined_fig <- stable_combined_grid |
  if (is.null(legend_grob)) plot_spacer() else wrap_elements(legend_grob)
stable_combined_fig <- stable_combined_fig +
  plot_layout(widths = c(1, 0.25))

stable_fig_png <- file.path(paper_figures_path, "FigS_Combined_Highlights_stable_regression_lines.png")
ggplot2::ggsave(stable_fig_png, plot = stable_combined_fig, width = 10, height = 8, dpi = 300, bg = "white")
message("wrote stable-regression Figure 4 variant: ", stable_fig_png)

message("Generating supplementary parameter-vs-environment figures...")

fig_param_vs_env_phi0_name <- "light_response_parameter_vs_environment_phi0.png"
fig_param_vs_env_P2000_name <- "light_response_parameter_vs_environment_P2000.png"
fig_param_vs_env_Rd_name <- "light_response_parameter_vs_environment_Rd.png"
fig_param_vs_env_LCP_name <- "light_response_parameter_vs_environment_LCP.png"

fig_supplementary_names <- c(
  phi0 = fig_param_vs_env_phi0_name,
  P2000 = fig_param_vs_env_P2000_name,
  Rd = fig_param_vs_env_Rd_name,
  LCP = fig_param_vs_env_LCP_name
)

for (param in param_vars) {
  plot_list <- list()

  panel_tags <- paste0("(", letters[seq_along(xvars_available)], ")")
  for (xvar in xvars_available) {
    plot_list[[xvar]] <- make_xy_plot_year(
      summary_table_yearly,
      xvar,
      param,
      ylims_map[[param]],
      stable_lm_only = TRUE
    ) +
      labs(tag = panel_tags[[match(xvar, xvars_available)]]) +
      theme(plot.tag = element_text(face = "bold", size = 14))
  }

  legend_dummy <- make_xy_plot_year(
    summary_table_yearly,
    xvars_available[1],
    param,
    ylims_map[[param]],
    show_legend = TRUE,
    stable_lm_only = TRUE
  )
  legend_grob <- extract_legend_grob(legend_dummy)
  plot_list[["legend"]] <- if (is.null(legend_grob)) plot_spacer() else wrap_elements(legend_grob)

  supp_grid <- wrap_plots(plot_list, ncol = 3)

  y_axis_grob <- grid::textGrob(ylab_shared_map[[param]], rot = 90, gp = grid::gpar(fontsize = 14))
  final_fig <- wrap_elements(y_axis_grob) | supp_grid
  final_fig <- final_fig + plot_layout(widths = c(0.04, 1))

  save_png_with_alias(
    final_fig,
    canonical_name = fig_supplementary_names[[param]],
    width = 14,
    height = 9
  )
}

message("All direct parameter highlight figures generated in: ", paper_figures_path)
