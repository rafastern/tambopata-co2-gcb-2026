# ───────────────────────────────────────────────────────────────────────────────
# residual panels of light-response parameters vs environment (current analysis)
# - update: use plotmath expressions in figure titles (Pmax -> P[max], Rd -> R[d], phi0 -> phi[0])
# - update: remove panel titles
# - update: option b for label box placement (place by a fraction of x/y ranges)
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(rlang)
  library(patchwork)
  library(viridisLite)
  library(grid)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths (update if needed)
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)

old_lp <- list.files(paper_figures_path, pattern = "_legendpanel\\.(png|pdf)$", full.names = TRUE)
if (length(old_lp)) file.remove(old_lp)

# ───────────────────────────────────────────────────────────────────────────────
# config
MIN_LIGHT_RESPONSE_POINTS <- 100
param_vars  <- c("P2000", "Rd", "phi0", "LCP")

xvars_order <- c("Ts", "SWC", "Ta", "VPD")
base_xs     <- c("Ts", "SWC", "Ta", "VPD")

season_fills       <- c("Wet" = "gray70", "Dry" = "white")
enso_shapes_filled <- c("El Niño" = 21, "La Niña" = 24, "Neutral" = 22)

xlab_map <- list(
  Ta  = expression(italic(T)[plain(a)] ~ "(°C)"),
  VPD = "VPD (kPa)",
  SWC = "SWC (cm)",
  Ts  = expression(italic(T)[plain(s)] ~ "(°C)")
)

# plotmath labels for figure titles
yvar_expr <- list(
  P2000 = quote(italic(P)[2000]),
  Rd    = quote(italic(R)[plain(d)]),
  phi0  = quote(italic(phi)[0]),
  LCP   = quote(LCP)
)

xvar_expr <- list(
  Ta  = quote(italic(T)[plain(a)]),
  Ts  = quote(italic(T)[plain(s)]),
  VPD = quote(VPD),
  SWC = quote(SWC)
)

# label-box placement (option b)
# - these are fractions of the x/y ranges used to "pull" the label in from top-right
LABEL_DX_FRAC <- 0.02
LABEL_DY_FRAC <- 0.01

# ───────────────────────────────────────────────────────────────────────────────
# load yearly parameter table produced by code/03_meteorology_summary.R
joined_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
if (!file.exists(joined_csv)) {
  stop("yearly parameter table not found: ", joined_csv, "\nRun code/03_meteorology_summary.R first to create it.")
}
summary_table_yearly <- read.csv(joined_csv, check.names = FALSE)
# guard: one row per year x group (fails loudly if the C1 fan-out duplication recurs)
if (all(c("Year", "Season", "ENSO") %in% names(summary_table_yearly)))
  stopifnot(!any(duplicated(summary_table_yearly[, c("Year", "Season", "ENSO")])))

if (!"P2000" %in% names(summary_table_yearly)) {
  stop(
    "P2000 column not found in light_response_mean_env_and_fit_params_by_year_with_Ts.csv. ",
    "Regenerate the yearly parameter table with code/03_meteorology_summary.R."
  )
}

# standardize key columns
if (!"Season" %in% names(summary_table_yearly)) {
  summary_table_yearly$Season <- ifelse(grepl("wet", tolower(summary_table_yearly$Group)), "Wet", "Dry")
}
if (!"ENSO" %in% names(summary_table_yearly)) {
  summary_table_yearly$ENSO <- dplyr::case_when(
    grepl("El Niño", summary_table_yearly$Group) ~ "El Niño",
    grepl("La Niña", summary_table_yearly$Group) ~ "La Niña",
    grepl("Neutral", summary_table_yearly$Group) ~ "Neutral",
    TRUE ~ NA_character_
  )
}

summary_table_yearly$Season <- factor(summary_table_yearly$Season, levels = c("Wet", "Dry"))
summary_table_yearly$ENSO   <- factor(summary_table_yearly$ENSO,   levels = c("El Niño", "La Niña", "Neutral"))

xvars_available <- intersect(xvars_order, names(summary_table_yearly))
if (length(xvars_available) < 2) stop("not enough x variables found in csv (need at least 2).")

if ("N" %in% names(summary_table_yearly)) {
  n_too_low <- sum(is.finite(summary_table_yearly$N) & summary_table_yearly$N < MIN_LIGHT_RESPONSE_POINTS)
  if (n_too_low > 0) {
    stop(
      "yearly parameter table contains ", n_too_low,
      " rows with N < ", MIN_LIGHT_RESPONSE_POINTS,
      ". Regenerate it with code/03_meteorology_summary.R."
    )
  }
  message("confirmed yearly light-response N >= ", MIN_LIGHT_RESPONSE_POINTS, " for retained rows.")
}

years_all_chr   <- sort(unique(as.character(summary_table_yearly$Year)))
year_colors_all <- setNames(viridisLite::turbo(length(years_all_chr)), years_all_chr)

# ───────────────────────────────────────────────────────────────────────────────
# helpers
lm_label <- function(d, xvar, yvar) {
  if (!nrow(d) || length(unique(d[[xvar]])) < 2 || length(unique(d[[yvar]])) < 2) return(NULL)
  fit <- lm(stats::as.formula(paste(yvar, "~", xvar)), data = d)
  s   <- summary(fit)
  a   <- unname(coef(fit)[1])
  b   <- unname(coef(fit)[2])
  r2  <- s$r.squared
  p   <- coef(s)[2, "Pr(>|t|)"]
  sprintf("resid = %.3g %s %.3g·x\nR²=%.2f; p=%.3g", a, ifelse(b >= 0, "+", "−"), abs(b), r2, p)
}

residual_lm_summary <- function(d, xvar, yvar) {
  if (nrow(d) < 4 || length(unique(d[[xvar]])) < 2 || length(unique(d[[yvar]])) < 2) {
    return(list(n = nrow(d), r2 = NA_real_, p = NA_real_, line_shown = FALSE))
  }

  fit <- tryCatch(
    lm(stats::as.formula(paste(yvar, "~", xvar)), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(n = nrow(d), r2 = NA_real_, p = NA_real_, line_shown = FALSE))
  }

  sm <- summary(fit)
  all_slope <- unname(coef(fit)[[xvar]])
  p_value <- unname(coef(sm)[xvar, "Pr(>|t|)"])
  if (!is.finite(all_slope) || all_slope == 0) {
    return(list(n = nrow(d), r2 = sm$r.squared, p = p_value, line_shown = FALSE))
  }

  loo_slopes <- vapply(seq_len(nrow(d)), function(i) {
    dd <- d[-i, , drop = FALSE]
    if (nrow(dd) < 3 || length(unique(dd[[xvar]])) < 2 || length(unique(dd[[yvar]])) < 2) {
      return(NA_real_)
    }

    loo_fit <- tryCatch(
      lm(stats::as.formula(paste(yvar, "~", xvar)), data = dd),
      error = function(e) NULL
    )
    if (is.null(loo_fit)) return(NA_real_)

    unname(coef(loo_fit)[[xvar]])
  }, numeric(1))

  list(
    n = nrow(d),
    r2 = sm$r.squared,
    p = p_value,
    line_shown = all(is.finite(loo_slopes)) && !any(sign(loo_slopes) != sign(all_slope))
  )
}

compute_residuals <- function(df, fit_x, yvar, xvars_available) {
  need <- unique(c(fit_x, yvar, "Season", "ENSO", "Year", xvars_available))
  d <- df %>%
    dplyr::select(any_of(need)) %>%
    dplyr::filter(
      is.finite(.data[[fit_x]]),
      is.finite(.data[[yvar]]),
      !is.na(Season),
      !is.na(ENSO),
      !is.na(Year)
    )
  
  if (!nrow(d) || length(unique(d[[fit_x]])) < 2) return(NULL)
  
  fml <- stats::as.formula(sprintf("`%s` ~ `%s`", yvar, fit_x))
  fit <- stats::lm(fml, data = d)
  
  d$residual <- stats::resid(fit)
  d$fitted   <- stats::fitted(fit)
  d
}

make_resid_panel <- function(dres, xv, years_use, year_colors_use,
                             show_lm = TRUE, show_lm_label = TRUE,
                             stable_lm_only = TRUE) {
  d <- dres %>%
    dplyr::filter(is.finite(.data[[xv]]), is.finite(residual)) %>%
    dplyr::mutate(
      x_value = .data[[xv]],
      Year = factor(as.character(Year), levels = years_use)
    )
  
  if (!nrow(d)) {
    return(ggplot() + theme_void() + labs(x = xlab_map[[xv]] %||% xv))
  }

  lm_data <- dplyr::rename(d, x = x_value, y = residual)
  lm_status <- residual_lm_summary(lm_data, "x", "y")
  draw_lm <- isTRUE(show_lm) && (!isTRUE(stable_lm_only) || isTRUE(lm_status$line_shown))
  
  xr <- range(d$x_value, na.rm = TRUE)
  w  <- 0.01 * diff(xr)
  if (!is.finite(w) || w <= 0) w <- 0
  set.seed(123)
  
  p <- ggplot(d, aes(x = x_value, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(
      aes(color = Year, fill = Season, shape = ENSO),
      size = 3.2,
      stroke = 1.6,
      alpha = 0.95,
      position = position_jitter(width = w, height = 0)
    ) +
    scale_color_manual(
      name   = "Year",
      values = year_colors_use,
      breaks = years_use,
      limits = years_use,
      drop   = FALSE
    ) +
    scale_fill_manual(
      name   = "Season",
      values = season_fills,
      breaks = c("Wet", "Dry"),
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "ENSO",
      values = enso_shapes_filled,
      breaks = c("El Niño", "La Niña", "Neutral"),
      drop   = FALSE
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(order = 2, override.aes = list(fill = "white", color = "black")),
      color = guide_legend(
        order = 3,
        ncol  = 2,
        byrow = TRUE,
        override.aes = list(
          shape  = 21,
          fill   = unname(year_colors_use[years_use]),
          size   = 3.2,
          stroke = 1.2
        )
      )
    ) +
    labs(
      x = xlab_map[[xv]] %||% xv,
      y = NULL
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  if (draw_lm) {
    p <- p +
      geom_smooth(
        aes(group = 1),
        method = "lm",
        se = FALSE,
        color = "black",
        linewidth = 1
      )
  }
  
  lbl <- if (isTRUE(show_lm_label) && draw_lm) lm_label(lm_data, "x", "y") else NULL
  if (!is.null(lbl)) {
    # option b: place by a fraction of axis ranges (pulled in from top-right)
    xr2 <- range(d$x_value, na.rm = TRUE)
    yr2 <- range(d$residual, na.rm = TRUE)
    
    dx <- LABEL_DX_FRAC * diff(xr2)
    dy <- LABEL_DY_FRAC * diff(yr2)
    if (!is.finite(dx) || dx <= 0) dx <- 0
    if (!is.finite(dy) || dy <= 0) dy <- 0
    
    p <- p + annotate(
      "label",
      x = Inf, y = Inf,
      label = lbl,
      hjust = 1.05, vjust = 1.05,  # increase these to push the box further up/right
      size = 3,
      alpha = 0.9
    )
    
  }
  
  p
}

make_legend_plot <- function(dres, years_use, year_colors_use) {
  # make a stable "legend-only" plot that always contains the full key set
  d <- dres %>%
    dplyr::filter(is.finite(residual), !is.na(Year), !is.na(Season), !is.na(ENSO)) %>%
    dplyr::mutate(
      Year = factor(as.character(Year), levels = years_use),
      x_dummy = 0,
      y_dummy = 0
    )
  
  ggplot(d, aes(x = x_dummy, y = y_dummy)) +
    geom_point(
      aes(color = Year, fill = Season, shape = ENSO),
      size = 3.2,
      stroke = 1.6,
      alpha = 0.95,
      position = position_jitter(width = 0.1, height = 0.1)
    ) +
    scale_color_manual(
      name   = "Year",
      values = year_colors_use,
      breaks = years_use,
      limits = years_use,
      drop   = FALSE
    ) +
    scale_fill_manual(
      name   = "Season",
      values = season_fills,
      breaks = c("Wet", "Dry"),
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "ENSO",
      values = enso_shapes_filled,
      breaks = c("El Niño", "La Niña", "Neutral"),
      drop   = FALSE
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(order = 2, override.aes = list(fill = "white", color = "black")),
      color = guide_legend(
        order = 3,
        ncol  = 2,
        byrow = TRUE,
        override.aes = list(
          shape  = 21,
          fill   = unname(year_colors_use[years_use]),
          size   = 3.2,
          stroke = 1.2
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11),
      legend.key.size = unit(4.5, "mm"),
      legend.spacing.y = unit(2, "mm")
    )
}

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx)) g$grobs[[idx[1]]] else NULL
}

save_plot_safely <- function(plot_obj, png_path, width, height, dpi = 300) {
  dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

residual_line_status <- list()

make_and_save_residual_panels <- function(yvar, base_fit_x, file_stub_prefix,
                                          summary_table_yearly, xvars_available,
                                          year_colors_all, paper_figures_path) {
  if (!base_fit_x %in% xvars_available) return(invisible(NULL))
  if (!yvar %in% names(summary_table_yearly)) return(invisible(NULL))
  
  dres <- compute_residuals(summary_table_yearly, base_fit_x, yvar, xvars_available)
  if (is.null(dres)) return(invisible(NULL))
  
  years_use       <- sort(unique(as.character(dres$Year)))
  year_colors_use <- year_colors_all[years_use]
  
  others <- setdiff(xvars_available, base_fit_x)
  if (!length(others)) return(invisible(NULL))

  status_rows <- lapply(others, function(xv) {
    d <- dres %>%
      dplyr::filter(is.finite(.data[[xv]]), is.finite(residual)) %>%
      dplyr::mutate(x_value = .data[[xv]])
    st <- residual_lm_summary(dplyr::rename(d, x = x_value, y = residual), "x", "y")
    data.frame(
      parameter = yvar,
      base_predictor = base_fit_x,
      xvar = xv,
      n = st$n,
      r_squared = st$r2,
      p_value = st$p,
      line_shown = st$line_shown,
      stringsAsFactors = FALSE
    )
  })
  residual_line_status[[paste(yvar, base_fit_x, sep = "_")]] <<- dplyr::bind_rows(status_rows)
  
  plots <- lapply(others, function(xv) {
    make_resid_panel(
      dres = dres,
      xv = xv,
      years_use = years_use,
      year_colors_use = year_colors_use
    )
  })
  
  # one legend (extracted)
  legend_plot  <- make_legend_plot(dres, years_use, year_colors_use)
  legend_grob  <- extract_legend_grob(legend_plot)
  legend_patch <- if (is.null(legend_grob)) plot_spacer() else wrap_elements(legend_grob)
  
  out_h <- 9.2
  
  # fixed 2x2: 3 panels + legend in bottom-right
  panel_tag_theme <- theme(
    plot.margin = margin(t = 24, r = 5, b = 5, l = 34),
    plot.tag = element_text(face = "bold", size = 16, margin = margin(r = 6, b = 4)),
    plot.tag.position = c(0, 1.02)
  )
  p1 <- plots[[1]] + labs(tag = "(a)") + panel_tag_theme
  p2 <- if (length(plots) >= 2) plots[[2]] + labs(tag = "(b)") + panel_tag_theme else plot_spacer()
  p3 <- if (length(plots) >= 3) plots[[3]] + labs(tag = "(c)") + panel_tag_theme else plot_spacer()
  
  panel_grid_2x2 <- (p1 | p2) / (p3 | legend_patch) +
    plot_layout(widths = c(1, 1), heights = c(1, 1))
  
  ylab_grob <- grid::textGrob("Residual", rot = 90, gp = grid::gpar(fontsize = 14))
  
  # build a plotmath title (with literal "~" printed)
  y_expr <- yvar_expr[[yvar]] %||% as.name(yvar)
  x_expr <- xvar_expr[[base_fit_x]] %||% as.name(base_fit_x)
  title_expr <- bquote(Residuals: ~ .(y_expr) ~ paste("~") ~ .(x_expr))
  
  fig <- (wrap_elements(ylab_grob) + panel_grid_2x2 + plot_layout(widths = c(0.06, 0.94))) +
    plot_annotation(title = title_expr) &
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0),
      plot.margin = margin(t = 8, r = 6, b = 6, l = 6)
    )
  
  outfile_png <- file.path(
    paper_figures_path,
    sprintf("%s_base%s_%s_panels.png", file_stub_prefix, base_fit_x, yvar)
  )
  if (grepl("_legendpanel", outfile_png, fixed = TRUE)) return(invisible(NULL))
  
  save_plot_safely(fig, outfile_png, width = 12.5, height = out_h, dpi = 300)
  message("wrote: ", outfile_png)

  legacy_png <- file.path(
    paper_figures_path,
    sprintf("Fig5_residuals_base%s_%s_residuals_panels.png", base_fit_x, yvar)
  )
  invisible(file.copy(outfile_png, legacy_png, overwrite = TRUE))
  message("wrote legacy alias: ", legacy_png)
  invisible(outfile_png)
}

# ───────────────────────────────────────────────────────────────────────────────
# run
for (yv in param_vars) {
  for (bx in base_xs) {
    make_and_save_residual_panels(
      yvar = yv,
      base_fit_x = bx,
      file_stub_prefix = "light_response_residuals",
      summary_table_yearly = summary_table_yearly,
      xvars_available = xvars_available,
      year_colors_all = year_colors_all,
      paper_figures_path = paper_figures_path
    )
  }
}

message("wrote residual grid figures in: ", paper_figures_path)

residual_line_status_table <- dplyr::bind_rows(residual_line_status) %>%
  dplyr::arrange(parameter, base_predictor, xvar)
status_csv <- file.path(output_path, "light_response_residual_line_status.csv")
write.csv(residual_line_status_table, status_csv, row.names = FALSE)
message("wrote residual line-status table: ", status_csv)
print(residual_line_status_table)

# ───────────────────────────────────────────────────────────────────────────────
# Generate MAIN TEXT Figure 5: Highlighted Residual Panels
# ───────────────────────────────────────────────────────────────────────────────
message("Generating Main Text Figure 5 (Residual Highlights)...")

# 1. Determine universal years to keep the legend consistent
years_use <- sort(unique(as.character(summary_table_yearly$Year[!is.na(summary_table_yearly$Year)])))
year_colors_use <- year_colors_all[years_use]

# 2. Compute the specific residuals needed for the 4 highlight panels
dres_phi0_Ta  <- compute_residuals(summary_table_yearly, "Ta", "phi0", xvars_available)
dres_LCP_Ts   <- compute_residuals(summary_table_yearly, "Ts", "LCP", xvars_available)
dres_P2000_SWC <- compute_residuals(summary_table_yearly, "SWC", "P2000", xvars_available)

# 3. Create the 4 individual panels with specific Y-axis math expressions
y_theme <- theme(
  axis.title.y = element_text(margin = margin(r = 10), size = 12),
  plot.margin = margin(t = 24, r = 5, b = 5, l = 34)
)
highlight_tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 16, margin = margin(r = 6, b = 4)),
  plot.tag.position = c(0, 1.02)
)

# Panel a: phi0 (base Ta) vs SWC
p_a <- make_resid_panel(dres_phi0_Ta, "SWC", years_use, year_colors_use) +
  labs(y = expression("Residual" ~ italic(phi)[0] ~ "(base" ~ italic(T)[plain(a)] * ")"), tag = "(a)") +
  y_theme + highlight_tag_theme

# Panel b: phi0 (base Ta) vs Ts
p_b <- make_resid_panel(dres_phi0_Ta, "Ts", years_use, year_colors_use) +
  labs(y = expression("Residual" ~ italic(phi)[0] ~ "(base" ~ italic(T)[plain(a)] * ")"), tag = "(b)") +
  y_theme + highlight_tag_theme

# Panel c: LCP (base Ts) vs Ta
p_c <- make_resid_panel(dres_LCP_Ts, "Ta", years_use, year_colors_use) +
  labs(y = expression("Residual LCP (base" ~ italic(T)[plain(s)] * ")"), tag = "(c)") +
  y_theme + highlight_tag_theme

# Panel d: P2000 (base SWC) vs VPD
p_d <- make_resid_panel(dres_P2000_SWC, "VPD", years_use, year_colors_use) +
  labs(y = expression("Residual" ~ italic(P)[2000] ~ "(base SWC)"), tag = "(d)") +
  y_theme + highlight_tag_theme

# 4. Generate the universal legend
legend_plot <- make_legend_plot(dres_phi0_Ta, years_use, year_colors_use)
legend_grob <- extract_legend_grob(legend_plot)

# 5. & 6. Assemble the 2x2 grid and attach legend
# We combine the panels and legend first, THEN apply the tags to the whole layout
data_grid <- (p_a + p_b) / (p_c + p_d)

legend_col <- if (is.null(legend_grob)) {
  plot_spacer()
} else {
  wrap_elements(legend_grob)
}

fig5_main <- data_grid | legend_col

fig5_main <- fig5_main +
  plot_layout(widths = c(1, 0.2))

# 7. Save the final Highlight Figure
fig5_main_path <- file.path(paper_figures_path, "light_response_residual_highlights_P2000.png")
save_plot_safely(fig5_main, fig5_main_path, width = 12, height = 9, dpi = 300)

fig5_legacy_path <- file.path(paper_figures_path, "Fig5_MainText_Residual_Highlights_P2000.png")
invisible(file.copy(fig5_main_path, fig5_legacy_path, overwrite = TRUE))

message("Saved residual highlight figure successfully at: ", fig5_main_path)
message("Saved legacy residual highlight alias at: ", fig5_legacy_path)

message("Generating exploratory Figure 5 residual variant...")
message("Regression lines and statistics are omitted; all residual points and zero-residual reference lines are retained.")

p_a_exploratory <- make_resid_panel(
  dres_phi0_Ta,
  "SWC",
  years_use,
  year_colors_use,
  show_lm = FALSE,
  show_lm_label = FALSE
) +
  labs(y = expression("Residual" ~ italic(phi)[0] ~ "(base" ~ italic(T)[plain(a)] * ")"), tag = "(a)") +
  y_theme + highlight_tag_theme

p_b_exploratory <- make_resid_panel(
  dres_phi0_Ta,
  "Ts",
  years_use,
  year_colors_use,
  show_lm = FALSE,
  show_lm_label = FALSE
) +
  labs(y = expression("Residual" ~ italic(phi)[0] ~ "(base" ~ italic(T)[plain(a)] * ")"), tag = "(b)") +
  y_theme + highlight_tag_theme

p_c_exploratory <- make_resid_panel(
  dres_LCP_Ts,
  "Ta",
  years_use,
  year_colors_use,
  show_lm = FALSE,
  show_lm_label = FALSE
) +
  labs(y = expression("Residual LCP (base" ~ italic(T)[plain(s)] * ")"), tag = "(c)") +
  y_theme + highlight_tag_theme

p_d_exploratory <- make_resid_panel(
  dres_P2000_SWC,
  "VPD",
  years_use,
  year_colors_use,
  show_lm = FALSE,
  show_lm_label = FALSE
) +
  labs(y = expression("Residual" ~ italic(P)[2000] ~ "(base SWC)"), tag = "(d)") +
  y_theme + highlight_tag_theme

data_grid_exploratory <- (p_a_exploratory + p_b_exploratory) / (p_c_exploratory + p_d_exploratory)

data_grid_exploratory <- data_grid_exploratory +
  plot_annotation(
    caption = "n = 4 complete cases per panel; exploratory residual patterns"
  ) &
  theme(
    plot.caption = element_text(size = 10, hjust = 0)
  )

fig5_exploratory <- data_grid_exploratory | legend_col

fig5_exploratory <- fig5_exploratory +
  plot_layout(widths = c(1, 0.2))

fig5_exploratory_path <- file.path(paper_figures_path, "Fig5_MainText_Residual_Highlights_exploratory.png")
save_plot_safely(fig5_exploratory, fig5_exploratory_path, width = 12, height = 9, dpi = 300)

message("Saved exploratory residual highlight figure at: ", fig5_exploratory_path)
# ───────────────────────────────────────────────────────────────────────────────
