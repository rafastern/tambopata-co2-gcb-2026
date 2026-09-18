#!/usr/bin/env Rscript

# influence and sensitivity diagnostics for Figure 4 and Figure 5 regressions.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(tibble)
  library(MASS)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

MIN_LIGHT_RESPONSE_POINTS <- 100

input_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
out_dir <- file.path(out_base, "diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_csv)) {
  stop("yearly parameter table not found: ", input_csv, "\nRun code/03_meteorology_summary.R first.")
}

required_cols <- c("Year", "Season", "ENSO", "Ta", "Ts", "SWC", "VPD", "P2000", "Rd", "phi0", "LCP")
raw <- readr::read_csv(input_csv, show_col_types = FALSE)
# guard: one row per year x group (fails loudly if the C1 fan-out duplication recurs)
if (all(c("Year", "Season", "ENSO") %in% names(raw)))
  stopifnot(!any(duplicated(raw[, c("Year", "Season", "ENSO")])))
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols)) {
  stop("missing required columns in ", input_csv, ": ", paste(missing_cols, collapse = ", "))
}

dat <- raw %>%
  dplyr::mutate(
    observation_id = dplyr::row_number(),
    Year = as.character(Year),
    Season = as.character(Season),
    ENSO = as.character(ENSO)
  ) %>%
  dplyr::filter(
    if ("N" %in% names(.)) is.na(N) | !is.finite(N) | N >= MIN_LIGHT_RESPONSE_POINTS else TRUE
  )

if ("N" %in% names(dat) && any(is.finite(dat$N) & dat$N < MIN_LIGHT_RESPONSE_POINTS)) {
  stop("rows with N < ", MIN_LIGHT_RESPONSE_POINTS, " remain after filtering.")
}

fig4_specs <- tibble::tribble(
  ~panel, ~response, ~predictor,
  "a", "phi0", "Ta",
  "b", "LCP", "Ts",
  "c", "P2000", "SWC",
  "d", "Rd", "VPD"
)

fig5_specs <- tibble::tribble(
  ~panel, ~response, ~base_predictor, ~secondary_predictor,
  "a", "phi0", "Ta", "SWC",
  "b", "phi0", "Ta", "Ts",
  "c", "LCP", "Ts", "Ta",
  "d", "P2000", "SWC", "VPD"
)

xlab_map <- list(
  Ta = expression(italic(T)[plain(a)] ~ "(°C)"),
  Ts = expression(italic(T)[plain(s)] ~ "(°C)"),
  SWC = "SWC (cm)",
  VPD = "VPD (kPa)"
)

ylab_map <- c(
  phi0 = "phi0",
  LCP = "LCP",
  P2000 = "P2000",
  Rd = "Rd"
)

safe_formula <- function(response, predictor) {
  stats::as.formula(sprintf("`%s` ~ `%s`", response, predictor))
}

finite_complete <- function(df, vars) {
  df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(vars), is.finite))
}

identity_cols <- function(df) {
  keep <- intersect(c("observation_id", "Year", "Season", "ENSO", "Group", "N"), names(df))
  df %>% dplyr::select(dplyr::all_of(keep))
}

identity_label <- function(df) {
  paste0("Year=", df$Year, "; Season=", df$Season, "; ENSO=", df$ENSO)
}

fit_lm_checked <- function(df, response, predictor) {
  if (nrow(df) < 3 || dplyr::n_distinct(df[[predictor]]) < 2 || dplyr::n_distinct(df[[response]]) < 2) {
    return(NULL)
  }
  stats::lm(safe_formula(response, predictor), data = df)
}

ols_summary_row <- function(fit, panel, response, predictor, model_label = "OLS") {
  sm <- summary(fit)
  co <- coef(sm)
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  slope_ci_lower <- if (!is.null(ci) && predictor %in% rownames(ci)) ci[predictor, 1] else NA_real_
  slope_ci_upper <- if (!is.null(ci) && predictor %in% rownames(ci)) ci[predictor, 2] else NA_real_

  tibble::tibble(
    panel = panel,
    model = model_label,
    response = response,
    predictor = predictor,
    n = stats::nobs(fit),
    residual_df = stats::df.residual(fit),
    intercept = unname(stats::coef(fit)[["(Intercept)"]]),
    slope = unname(stats::coef(fit)[[predictor]]),
    r_squared = unname(sm$r.squared),
    adjusted_r_squared = unname(sm$adj.r.squared),
    p_value = unname(co[predictor, "Pr(>|t|)"]),
    slope_ci_lower = slope_ci_lower,
    slope_ci_upper = slope_ci_upper
  )
}

robust_summary_row <- function(df, panel, response, predictor, model_label = "robust rlm") {
  out <- tryCatch({
    fit <- MASS::rlm(safe_formula(response, predictor), data = df, maxit = 100)
    tibble::tibble(
      panel = panel,
      model = model_label,
      response = response,
      predictor = predictor,
      n = nrow(stats::model.frame(fit)),
      robust_intercept = unname(stats::coef(fit)[["(Intercept)"]]),
      robust_slope = unname(stats::coef(fit)[[predictor]])
    )
  }, error = function(e) {
    tibble::tibble(
      panel = panel,
      model = model_label,
      response = response,
      predictor = predictor,
      n = nrow(df),
      robust_intercept = NA_real_,
      robust_slope = NA_real_
    )
  })
  out
}

influence_rows <- function(df, fit, panel, response, predictor, model_label = "OLS") {
  p <- length(stats::coef(fit))
  n <- stats::nobs(fit)
  lev <- stats::hatvalues(fit)
  cooks <- stats::cooks.distance(fit)
  stud <- stats::rstudent(fit)
  dff <- stats::dffits(fit)
  flags <- mapply(function(h, ck, st, dfv) {
    hit <- c(
      if (is.finite(h) && h > 2 * p / n) "leverage_gt_2p_over_n" else character(),
      if (is.finite(h) && h > 3 * p / n) "leverage_gt_3p_over_n" else character(),
      if (is.finite(ck) && ck > 4 / n) "cooks_gt_4_over_n" else character(),
      if (is.finite(st) && abs(st) > 3) "abs_studentized_residual_gt_3" else character(),
      if (is.finite(dfv) && abs(dfv) > 2 * sqrt(p / n)) "abs_dffits_gt_2sqrt_p_over_n" else character()
    )
    if (length(hit)) paste(hit, collapse = ";") else NA_character_
  }, lev, cooks, stud, dff, USE.NAMES = FALSE)

  dplyr::bind_cols(
    tibble::tibble(
      panel = panel,
      model = model_label,
      response = response,
      predictor = predictor
    ),
    identity_cols(df),
    tibble::tibble(
      predictor_value = df[[predictor]],
      response_value = df[[response]],
      fitted_value = stats::fitted(fit),
      residual = stats::resid(fit),
      leverage = lev,
      cooks_distance = cooks,
      studentized_residual = stud,
      dffits = dff,
      leverage_threshold_2p_over_n = 2 * p / n,
      leverage_threshold_3p_over_n = 3 * p / n,
      cooks_threshold_4_over_n = 4 / n,
      dffits_threshold_2sqrt_p_over_n = 2 * sqrt(p / n),
      flags_triggered = flags
    )
  )
}

loo_rows <- function(df, panel, response, predictor, all_fit, model_label = "OLS") {
  all_slope <- unname(stats::coef(all_fit)[[predictor]])
  rows <- lapply(seq_len(nrow(df)), function(i) {
    dd <- df[-i, , drop = FALSE]
    fit <- fit_lm_checked(dd, response, predictor)
    if (is.null(fit)) {
      return(dplyr::bind_cols(
        tibble::tibble(
          panel = panel,
          model = model_label,
          response = response,
          predictor = predictor
        ),
        identity_cols(df[i, , drop = FALSE]),
        tibble::tibble(
          omitted_predictor_value = df[[predictor]][i],
          omitted_response_value = df[[response]][i],
          n = nrow(dd),
          intercept = NA_real_,
          slope = NA_real_,
          r_squared = NA_real_,
          p_value = NA_real_,
          slope_sign_changes = NA,
          slope_abs_change = NA_real_,
          slope_relative_change = NA_real_
        )
      ))
    }
    sm <- summary(fit)
    slope <- unname(stats::coef(fit)[[predictor]])
    dplyr::bind_cols(
      tibble::tibble(
        panel = panel,
        model = model_label,
        response = response,
        predictor = predictor
      ),
      identity_cols(df[i, , drop = FALSE]),
      tibble::tibble(
        omitted_predictor_value = df[[predictor]][i],
        omitted_response_value = df[[response]][i],
        n = stats::nobs(fit),
        intercept = unname(stats::coef(fit)[["(Intercept)"]]),
        slope = slope,
        r_squared = unname(sm$r.squared),
        p_value = unname(coef(sm)[predictor, "Pr(>|t|)"]),
        slope_sign_changes = sign(slope) != sign(all_slope),
        slope_abs_change = abs(slope - all_slope),
        slope_relative_change = ifelse(is.finite(all_slope) && all_slope != 0, abs(slope - all_slope) / abs(all_slope), NA_real_)
      )
    )
  })
  dplyr::bind_rows(rows)
}

panel_data <- function(df, response, predictor) {
  df %>%
    dplyr::filter(
      is.finite(.data[[response]]),
      is.finite(.data[[predictor]]),
      !is.na(Year),
      !is.na(Season),
      !is.na(ENSO)
    )
}

fig4_ols <- list()
fig4_influence <- list()
fig4_loo <- list()
fig4_robust <- list()
fig4_plot_data <- list()

for (i in seq_len(nrow(fig4_specs))) {
  spec <- fig4_specs[i, ]
  d <- panel_data(dat, spec$response, spec$predictor)
  fit <- fit_lm_checked(d, spec$response, spec$predictor)
  if (is.null(fit)) next

  fig4_ols[[spec$panel]] <- ols_summary_row(fit, spec$panel, spec$response, spec$predictor)
  fig4_influence[[spec$panel]] <- influence_rows(d, fit, spec$panel, spec$response, spec$predictor)
  fig4_loo[[spec$panel]] <- loo_rows(d, spec$panel, spec$response, spec$predictor, fit)
  fig4_robust[[spec$panel]] <- robust_summary_row(d, spec$panel, spec$response, spec$predictor)
  fig4_plot_data[[spec$panel]] <- list(data = d, fit = fit, spec = spec)
}

fig4_ols_summary <- dplyr::bind_rows(fig4_ols)
fig4_influence_diagnostics <- dplyr::bind_rows(fig4_influence)
fig4_leave_one_out <- dplyr::bind_rows(fig4_loo)
fig4_robust_summary <- dplyr::bind_rows(fig4_robust)

make_sensitivity_summary <- function(ols_tbl, infl_tbl, loo_tbl, robust_tbl) {
  ols_tbl %>%
    dplyr::left_join(
      robust_tbl %>% dplyr::select(panel, robust_slope),
      by = "panel"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      leave_one_out_slope_min = min(loo_tbl$slope[loo_tbl$panel == panel], na.rm = TRUE),
      leave_one_out_slope_max = max(loo_tbl$slope[loo_tbl$panel == panel], na.rm = TRUE),
      leave_one_out_sign_changes = sum(loo_tbl$slope_sign_changes[loo_tbl$panel == panel], na.rm = TRUE),
      largest_slope_change_observation = {
        d <- loo_tbl %>% dplyr::filter(panel == .env$panel) %>% dplyr::arrange(dplyr::desc(slope_abs_change)) %>% dplyr::slice(1)
        if (nrow(d)) identity_label(d) else NA_character_
      },
      largest_slope_change_abs = {
        d <- loo_tbl %>% dplyr::filter(panel == .env$panel) %>% dplyr::arrange(dplyr::desc(slope_abs_change)) %>% dplyr::slice(1)
        if (nrow(d)) d$slope_abs_change else NA_real_
      },
      max_cooks_observation = {
        d <- infl_tbl %>% dplyr::filter(panel == .env$panel) %>% dplyr::arrange(dplyr::desc(cooks_distance)) %>% dplyr::slice(1)
        if (nrow(d)) identity_label(d) else NA_character_
      },
      max_cooks_distance = {
        d <- infl_tbl %>% dplyr::filter(panel == .env$panel) %>% dplyr::arrange(dplyr::desc(cooks_distance)) %>% dplyr::slice(1)
        if (nrow(d)) d$cooks_distance else NA_real_
      },
      robust_slope_sign_agrees = ifelse(is.finite(robust_slope), sign(robust_slope) == sign(slope), NA)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      panel, response, predictor, n,
      all_points_slope = slope,
      all_points_r_squared = r_squared,
      all_points_p_value = p_value,
      leave_one_out_slope_min,
      leave_one_out_slope_max,
      leave_one_out_sign_changes,
      largest_slope_change_observation,
      largest_slope_change_abs,
      max_cooks_observation,
      max_cooks_distance,
      robust_slope,
      robust_slope_sign_agrees
    )
}

fig4_sensitivity_summary <- make_sensitivity_summary(
  fig4_ols_summary,
  fig4_influence_diagnostics,
  fig4_leave_one_out,
  fig4_robust_summary
)

readr::write_csv(fig4_ols_summary, file.path(out_dir, "fig4_ols_summary.csv"))
readr::write_csv(fig4_influence_diagnostics, file.path(out_dir, "fig4_influence_diagnostics.csv"))
readr::write_csv(fig4_leave_one_out, file.path(out_dir, "fig4_leave_one_out.csv"))
readr::write_csv(fig4_sensitivity_summary, file.path(out_dir, "fig4_sensitivity_summary.csv"))

prediction_lines <- function(info) {
  d <- info$data
  spec <- info$spec
  x <- spec$predictor
  y <- spec$response
  xr <- range(d[[x]], na.rm = TRUE)
  grid <- tibble::tibble(x_value = seq(xr[1], xr[2], length.out = 100))
  names(grid)[1] <- x
  all_line <- grid %>%
    dplyr::mutate(
      y_value = stats::predict(info$fit, newdata = grid),
      line_type = "all points",
      omitted_label = NA_character_
    )

  loo_lines <- lapply(seq_len(nrow(d)), function(i) {
    fit <- fit_lm_checked(d[-i, , drop = FALSE], y, x)
    if (is.null(fit)) return(NULL)
    grid %>%
      dplyr::mutate(
        y_value = stats::predict(fit, newdata = grid),
        line_type = "leave one out",
        omitted_label = identity_label(d[i, , drop = FALSE])
      )
  })

  robust_fit <- tryCatch(MASS::rlm(safe_formula(y, x), data = d, maxit = 100), error = function(e) NULL)
  robust_line <- if (is.null(robust_fit)) {
    NULL
  } else {
    grid %>%
      dplyr::mutate(
        y_value = stats::predict(robust_fit, newdata = grid),
        line_type = "robust rlm",
        omitted_label = NA_character_
      )
  }

  dplyr::bind_rows(all_line, dplyr::bind_rows(loo_lines), robust_line) %>%
    dplyr::mutate(panel = spec$panel, response = y, predictor = x)
}

make_fig4_lines_panel <- function(info) {
  d <- info$data
  spec <- info$spec
  x <- spec$predictor
  y <- spec$response
  lines <- prediction_lines(info)
  infl <- fig4_influence_diagnostics %>% dplyr::filter(panel == spec$panel)
  loo <- fig4_leave_one_out %>% dplyr::filter(panel == spec$panel)
  max_cook <- infl %>% dplyr::arrange(dplyr::desc(cooks_distance)) %>% dplyr::slice(1)
  max_change <- loo %>% dplyr::arrange(dplyr::desc(slope_abs_change)) %>% dplyr::slice(1)
  annotate_ids <- unique(c(max_cook$observation_id, max_change$observation_id))
  ann <- d %>%
    dplyr::filter(observation_id %in% annotate_ids) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        observation_id == max_cook$observation_id & observation_id == max_change$observation_id ~ "max Cook / slope change",
        observation_id == max_cook$observation_id ~ "max Cook",
        observation_id == max_change$observation_id ~ "max slope change",
        TRUE ~ NA_character_
      )
    )

  ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
    geom_line(
      data = lines %>% dplyr::filter(line_type == "leave one out"),
      aes(x = .data[[x]], y = y_value, group = omitted_label),
      inherit.aes = FALSE,
      color = "gray70",
      linewidth = 0.35,
      alpha = 0.75
    ) +
    geom_line(
      data = lines %>% dplyr::filter(line_type == "all points"),
      aes(x = .data[[x]], y = y_value),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    geom_line(
      data = lines %>% dplyr::filter(line_type == "robust rlm"),
      aes(x = .data[[x]], y = y_value),
      inherit.aes = FALSE,
      color = "#0072B2",
      linewidth = 0.9
    ) +
    geom_point(aes(fill = Season), shape = 21, size = 3.0, stroke = 1.0, color = "black") +
    geom_text(
      data = ann,
      aes(label = label),
      nudge_y = 0.04 * diff(range(d[[y]], na.rm = TRUE)),
      size = 3,
      check_overlap = TRUE
    ) +
    scale_fill_manual(values = c(Wet = "gray70", Dry = "white"), na.value = "gray85") +
    labs(
      title = paste0("(", spec$panel, ") ", ylab_map[[y]], " ~ ", x),
      x = xlab_map[[x]] %||% x,
      y = ylab_map[[y]] %||% y
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
}

line_panels <- lapply(fig4_plot_data, make_fig4_lines_panel)
fig4_lines_plot <- wrap_plots(line_panels, ncol = 2) +
  plot_annotation(caption = "black dashed = all-points OLS; gray = leave-one-out OLS; blue = robust rlm")
ggplot2::ggsave(
  file.path(out_dir, "fig4_leave_one_out_lines.png"),
  plot = fig4_lines_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

fig4_cooks_data <- fig4_influence_diagnostics %>%
  dplyr::mutate(observation_label = identity_label(.))

fig4_cooks_plot <- ggplot(fig4_cooks_data, aes(x = reorder(observation_label, cooks_distance), y = cooks_distance)) +
  geom_col(fill = "gray55") +
  geom_hline(aes(yintercept = cooks_threshold_4_over_n), linetype = "dashed", color = "red") +
  facet_wrap(~ panel + response + predictor, scales = "free_x") +
  coord_flip() +
  labs(x = NULL, y = "Cook's distance", caption = "red dashed line = 4/n threshold") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))
ggplot2::ggsave(
  file.path(out_dir, "fig4_cooks_distance.png"),
  plot = fig4_cooks_plot,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

compute_residual_dataset <- function(df, response, base_predictor) {
  base_df <- panel_data(df, response, base_predictor)
  fit <- fit_lm_checked(base_df, response, base_predictor)
  if (is.null(fit)) return(NULL)
  base_df %>%
    dplyr::mutate(
      primary_fitted = stats::fitted(fit),
      residual = stats::resid(fit)
    )
}

fig5_ols <- list()
fig5_influence <- list()
fig5_loo <- list()
fig5_robust <- list()

for (i in seq_len(nrow(fig5_specs))) {
  spec <- fig5_specs[i, ]
  dres <- compute_residual_dataset(dat, spec$response, spec$base_predictor)
  if (is.null(dres)) next

  rd <- dres %>%
    dplyr::filter(is.finite(.data[[spec$secondary_predictor]]), is.finite(residual)) %>%
    dplyr::mutate(residual_response = residual)
  fit <- fit_lm_checked(rd, "residual_response", spec$secondary_predictor)
  if (is.null(fit)) next

  fig5_ols[[spec$panel]] <- ols_summary_row(
    fit,
    spec$panel,
    paste0("residual(", spec$response, " ~ ", spec$base_predictor, ")"),
    spec$secondary_predictor,
    "residual OLS"
  ) %>%
    dplyr::mutate(base_predictor = spec$base_predictor, original_response = spec$response)

  fig5_influence[[spec$panel]] <- influence_rows(
    rd,
    fit,
    spec$panel,
    paste0("residual(", spec$response, " ~ ", spec$base_predictor, ")"),
    spec$secondary_predictor,
    "residual OLS"
  ) %>%
    dplyr::mutate(base_predictor = spec$base_predictor, original_response = spec$response)

  all_slope <- unname(stats::coef(fit)[[spec$secondary_predictor]])
  all_p <- unname(coef(summary(fit))[spec$secondary_predictor, "Pr(>|t|)"])
  loo_panel <- lapply(seq_len(nrow(dres)), function(j) {
    primary_minus <- dres[-j, , drop = FALSE] %>% dplyr::select(-primary_fitted, -residual)
    primary_fit <- fit_lm_checked(primary_minus, spec$response, spec$base_predictor)
    omitted <- dres[j, , drop = FALSE]
    if (is.null(primary_fit)) {
      return(dplyr::bind_cols(
        tibble::tibble(
          panel = spec$panel,
          model = "residual OLS",
          original_response = spec$response,
          base_predictor = spec$base_predictor,
          response = paste0("residual(", spec$response, " ~ ", spec$base_predictor, ")"),
          predictor = spec$secondary_predictor
        ),
        identity_cols(omitted),
        tibble::tibble(
          omitted_base_predictor_value = omitted[[spec$base_predictor]],
          omitted_secondary_predictor_value = omitted[[spec$secondary_predictor]],
          omitted_response_value = omitted[[spec$response]],
          omitted_in_secondary_regression = is.finite(omitted[[spec$secondary_predictor]]),
          n = NA_integer_,
          intercept = NA_real_,
          slope = NA_real_,
          r_squared = NA_real_,
          p_value = NA_real_,
          slope_sign_changes = NA,
          slope_abs_change = NA_real_,
          slope_relative_change = NA_real_,
          p_lt_0.05 = NA,
          conclusion_changed_p05 = NA
        )
      ))
    }

    recomputed <- primary_minus %>%
      dplyr::mutate(
        residual_response = stats::resid(primary_fit)
      ) %>%
      dplyr::filter(is.finite(.data[[spec$secondary_predictor]]), is.finite(residual_response))
    resid_fit <- fit_lm_checked(recomputed, "residual_response", spec$secondary_predictor)
    if (is.null(resid_fit)) {
      slope <- NA_real_
      sm <- NULL
      p_val <- NA_real_
      r2 <- NA_real_
      intercept <- NA_real_
      n_fit <- nrow(recomputed)
    } else {
      sm <- summary(resid_fit)
      slope <- unname(stats::coef(resid_fit)[[spec$secondary_predictor]])
      p_val <- unname(coef(sm)[spec$secondary_predictor, "Pr(>|t|)"])
      r2 <- unname(sm$r.squared)
      intercept <- unname(stats::coef(resid_fit)[["(Intercept)"]])
      n_fit <- stats::nobs(resid_fit)
    }

    dplyr::bind_cols(
      tibble::tibble(
        panel = spec$panel,
        model = "residual OLS",
        original_response = spec$response,
        base_predictor = spec$base_predictor,
        response = paste0("residual(", spec$response, " ~ ", spec$base_predictor, ")"),
        predictor = spec$secondary_predictor
      ),
      identity_cols(omitted),
      tibble::tibble(
        omitted_base_predictor_value = omitted[[spec$base_predictor]],
        omitted_secondary_predictor_value = omitted[[spec$secondary_predictor]],
        omitted_response_value = omitted[[spec$response]],
        omitted_in_secondary_regression = is.finite(omitted[[spec$secondary_predictor]]),
        n = n_fit,
        intercept = intercept,
        slope = slope,
        r_squared = r2,
        p_value = p_val,
        slope_sign_changes = ifelse(is.finite(slope), sign(slope) != sign(all_slope), NA),
        slope_abs_change = ifelse(is.finite(slope), abs(slope - all_slope), NA_real_),
        slope_relative_change = ifelse(is.finite(slope) && is.finite(all_slope) && all_slope != 0, abs(slope - all_slope) / abs(all_slope), NA_real_),
        p_lt_0.05 = ifelse(is.finite(p_val), p_val < 0.05, NA),
        conclusion_changed_p05 = ifelse(is.finite(p_val), (p_val < 0.05) != (all_p < 0.05), NA)
      )
    )
  })
  fig5_loo[[spec$panel]] <- dplyr::bind_rows(loo_panel)
  fig5_robust[[spec$panel]] <- robust_summary_row(rd, spec$panel, "residual_response", spec$secondary_predictor, "residual robust rlm") %>%
    dplyr::mutate(
      response = paste0("residual(", spec$response, " ~ ", spec$base_predictor, ")"),
      base_predictor = spec$base_predictor,
      original_response = spec$response
    )
}

fig5_residual_ols_summary <- dplyr::bind_rows(fig5_ols) %>%
  dplyr::select(panel, model, original_response, base_predictor, response, predictor, dplyr::everything())
fig5_residual_influence_diagnostics <- dplyr::bind_rows(fig5_influence) %>%
  dplyr::select(panel, model, original_response, base_predictor, response, predictor, dplyr::everything())
fig5_residual_leave_one_out <- dplyr::bind_rows(fig5_loo)
fig5_residual_robust_summary <- dplyr::bind_rows(fig5_robust)

fig5_residual_sensitivity_summary <- make_sensitivity_summary(
  fig5_residual_ols_summary,
  fig5_residual_influence_diagnostics,
  fig5_residual_leave_one_out,
  fig5_residual_robust_summary
) %>%
  dplyr::left_join(
    fig5_residual_ols_summary %>% dplyr::select(panel, original_response, base_predictor),
    by = "panel"
  ) %>%
  dplyr::select(panel, original_response, base_predictor, dplyr::everything())

readr::write_csv(fig5_residual_ols_summary, file.path(out_dir, "fig5_residual_ols_summary.csv"))
readr::write_csv(fig5_residual_influence_diagnostics, file.path(out_dir, "fig5_residual_influence_diagnostics.csv"))
readr::write_csv(fig5_residual_leave_one_out, file.path(out_dir, "fig5_residual_leave_one_out.csv"))
readr::write_csv(fig5_residual_sensitivity_summary, file.path(out_dir, "fig5_residual_sensitivity_summary.csv"))

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "g"))
}

summarize_panel_sentence <- function(row, figure_label) {
  sensitivity <- if (isTRUE(row$leave_one_out_sign_changes > 0)) {
    "leave-one-out fits changed the slope sign"
  } else {
    "leave-one-out fits did not change the slope sign"
  }
  robust_txt <- if (isTRUE(row$robust_slope_sign_agrees)) {
    "the robust slope had the same sign as OLS"
  } else if (identical(row$robust_slope_sign_agrees, FALSE)) {
    "the robust slope changed sign relative to OLS"
  } else {
    "the robust slope could not be estimated"
  }
  paste0(
    "- ", figure_label, " panel ", row$panel, " (`", row$response, " ~ ", row$predictor, "`): ",
    "OLS slope = ", fmt(row$all_points_slope, 6),
    ", R2 = ", fmt(row$all_points_r_squared, 3),
    ", p = ", fmt_p(row$all_points_p_value),
    "; LOO slope range = [", fmt(row$leave_one_out_slope_min, 6), ", ", fmt(row$leave_one_out_slope_max, 6), "]; ",
    sensitivity, "; ", robust_txt, "."
  )
}

fig4_a <- fig4_sensitivity_summary %>% dplyr::filter(panel == "a") %>% dplyr::slice(1)
fig4_b <- fig4_sensitivity_summary %>% dplyr::filter(panel == "b") %>% dplyr::slice(1)
fig4_high_flags <- fig4_influence_diagnostics %>% dplyr::filter(!is.na(flags_triggered))
fig5_high_flags <- fig5_residual_influence_diagnostics %>% dplyr::filter(!is.na(flags_triggered))
fig5_conclusion_changes <- fig5_residual_leave_one_out %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(any_p05_change = any(conclusion_changed_p05, na.rm = TRUE), .groups = "drop")

report <- c(
  "# Influence and Sensitivity Report",
  "",
  paste0("Input CSV: `", input_csv, "`"),
  paste0("Rows after existing `N >= ", MIN_LIGHT_RESPONSE_POINTS, "` light-response QC: ", nrow(dat)),
  "",
  "## Figure 4",
  "",
  vapply(seq_len(nrow(fig4_sensitivity_summary)), function(i) summarize_panel_sentence(fig4_sensitivity_summary[i, ], "Figure 4"), character(1)),
  "",
  "Figure 4 panels a and b should be interpreted cautiously if their leave-one-out slope ranges are wide or if their maximum Cook's distance rows also trigger leverage or DFFITS criteria. The point was retained because there was no predefined QC basis for exclusion.",
  paste0(
    "Panel a sensitivity summary: largest slope-change observation = ",
    fig4_a$largest_slope_change_observation,
    "; maximum Cook's distance observation = ",
    fig4_a$max_cooks_observation,
    "."
  ),
  paste0(
    "Panel b sensitivity summary: largest slope-change observation = ",
    fig4_b$largest_slope_change_observation,
    "; maximum Cook's distance observation = ",
    fig4_b$max_cooks_observation,
    "."
  ),
  "",
  "## Formal Influence Flags",
  "",
  if (nrow(fig4_high_flags)) {
    paste0(
      "- Figure 4 panel ", fig4_high_flags$panel,
      " observation Year=", fig4_high_flags$Year,
      ", Season=", fig4_high_flags$Season,
      ", ENSO=", fig4_high_flags$ENSO,
      " triggered: ", fig4_high_flags$flags_triggered,
      "."
    )
  } else {
    "- No Figure 4 observations triggered the specified influence thresholds."
  },
  "",
  "## Figure 5",
  "",
  vapply(seq_len(nrow(fig5_residual_sensitivity_summary)), function(i) summarize_panel_sentence(fig5_residual_sensitivity_summary[i, ], "Figure 5"), character(1)),
  "",
  "The Figure 5 residual relationships are based on four secondary-regression observations per main panel in the current data. They are useful as exploratory diagnostics, but manuscript language should avoid strong mechanistic claims from these regressions alone.",
  if (any(fig5_conclusion_changes$any_p05_change, na.rm = TRUE)) {
    paste0(
      "At least one Figure 5 leave-one-out fit changed the p < 0.05 conclusion in panels: ",
      paste(fig5_conclusion_changes$panel[fig5_conclusion_changes$any_p05_change], collapse = ", "),
      "."
    )
  } else {
    "No Figure 5 leave-one-out fit changed a p < 0.05 conclusion; however, all main residual panels have very small n."
  },
  "",
  "Figure 5 residual-regression influence flags:",
  "",
  if (nrow(fig5_high_flags)) {
    paste0(
      "- Figure 5 panel ", fig5_high_flags$panel,
      " observation Year=", fig5_high_flags$Year,
      ", Season=", fig5_high_flags$Season,
      ", ENSO=", fig5_high_flags$ENSO,
      " triggered: ", fig5_high_flags$flags_triggered,
      "."
    )
  } else {
    "- No Figure 5 residual-regression observations triggered the specified influence thresholds."
  },
  "",
  "## Manuscript Guidance",
  "",
  "- Keep all observations in the main analysis unless a predefined QC rule fails.",
  "- Use language such as: the apparent bivariate slope is sensitive to individual observations, the relationship should be interpreted cautiously, and the observation is influential in the regression sense.",
  "- Avoid labeling retained observations as quality-control failures unless separate QC evidence supports that label.",
  "- For Figure 5, emphasize exploratory residual associations and the small number of complete SWC/Ts rows.",
  "",
  "## Output Files",
  "",
  "- `fig4_ols_summary.csv`",
  "- `fig4_influence_diagnostics.csv`",
  "- `fig4_leave_one_out.csv`",
  "- `fig4_sensitivity_summary.csv`",
  "- `fig4_leave_one_out_lines.png`",
  "- `fig4_cooks_distance.png`",
  "- `fig5_residual_ols_summary.csv`",
  "- `fig5_residual_influence_diagnostics.csv`",
  "- `fig5_residual_leave_one_out.csv`",
  "- `fig5_residual_sensitivity_summary.csv`"
)

writeLines(report, file.path(out_dir, "influence_sensitivity_report.md"))

message("wrote influence and sensitivity diagnostics to: ", out_dir)
