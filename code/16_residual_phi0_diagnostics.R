#!/usr/bin/env Rscript

# Diagnostic companion to Figure 5 residual analysis.
# This script does not overwrite figures or alter the current Figure 5 workflow.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

MIN_LIGHT_RESPONSE_POINTS <- 100
input_csv <- file.path(output_path, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")
out_dir <- file.path(out_base, "residual_phi0_diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv)
}

needed <- c("Year", "Group", "Season", "ENSO", "N", "phi0", "Ta", "SWC", "Ts", "VPD")
raw <- readr::read_csv(input_csv, show_col_types = FALSE)
# guard: one row per year x group (fails loudly if the C1 fan-out duplication recurs)
if (all(c("Year", "Season", "ENSO") %in% names(raw)))
  stopifnot(!any(duplicated(raw[, c("Year", "Season", "ENSO")])))
missing <- setdiff(c("phi0", "Ta"), names(raw))
if (length(missing)) {
  stop("Missing required columns in ", input_csv, ": ", paste(missing, collapse = ", "))
}

if (!"Season" %in% names(raw)) {
  raw$Season <- ifelse(grepl("wet", tolower(raw$Group)), "Wet", "Dry")
}
if (!"ENSO" %in% names(raw)) {
  raw$ENSO <- NA_character_
}

dat <- raw %>%
  dplyr::select(dplyr::any_of(needed)) %>%
  dplyr::filter(
    if ("N" %in% names(.)) is.na(N) | !is.finite(N) | N >= MIN_LIGHT_RESPONSE_POINTS else TRUE
  )

finite_complete <- function(df, vars) {
  df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(vars), is.finite))
}

tidy_lm_terms <- function(fit, response, predictor_label, controls = character(),
                          model_type, formula_label, partial = FALSE) {
  sm <- summary(fit)
  co <- as.data.frame(coef(sm))
  co$term <- rownames(co)
  rownames(co) <- NULL

  tibble::tibble(
    response = response,
    predictor = co$term,
    control_variables = if (length(controls)) paste(controls, collapse = " + ") else "none",
    model_formula = formula_label,
    model_type = model_type,
    n = stats::nobs(fit),
    residual_df = stats::df.residual(fit),
    slope = co$Estimate,
    standard_error = co$`Std. Error`,
    statistic = co$`t value`,
    p_value = co$`Pr(>|t|)`,
    r_squared = unname(sm$r.squared),
    adjusted_r_squared = unname(sm$adj.r.squared),
    partial_r_squared = ifelse(partial & co$term != "(Intercept)",
                               co$`t value`^2 / (co$`t value`^2 + stats::df.residual(fit)),
                               NA_real_),
    p_value_type = ifelse(model_type == "FWL residual-on-residual diagnostic",
                          "diagnostic two-sided OLS t-test; use multiple-regression term p-value for formal inference",
                          "two-sided OLS t-test")
  ) %>%
    dplyr::filter(predictor != "(Intercept)") %>%
    dplyr::mutate(predictor = ifelse(predictor == predictor_label, predictor_label, predictor))
}

fit_lm_safe <- function(df, formula) {
  tryCatch(stats::lm(formula, data = df), error = function(e) e)
}

model_rows <- list()

# A. Current-style residual analysis: phi0 ~ Ta, then residual(phi0|Ta) ~ x.
base_df <- finite_complete(dat, c("phi0", "Ta"))
base_fit <- stats::lm(phi0 ~ Ta, data = base_df)
model_rows[["current_base"]] <- tidy_lm_terms(
  base_fit, "phi0", "Ta", character(), "current base residual model",
  "phi0 ~ Ta"
)

resid_df <- base_df %>%
  dplyr::mutate(phi0_resid_Ta = stats::resid(base_fit))

for (x in c("SWC", "Ts", "VPD")) {
  if (!x %in% names(resid_df)) next
  d <- finite_complete(resid_df, c("phi0_resid_Ta", x))
  if (nrow(d) < 3 || length(unique(d[[x]])) < 2) next
  fit <- stats::lm(stats::as.formula(paste("phi0_resid_Ta ~", x)), data = d)
  model_rows[[paste0("current_resid_", x)]] <- tidy_lm_terms(
    fit, "residual(phi0 | Ta)", x, "Ta", "current residual-on-predictor model",
    paste0("residual(phi0 ~ Ta) ~ ", x)
  )
}

# B. Multiple regression models.
multiple_specs <- list(
  c("Ta", "SWC"),
  c("Ta", "Ts"),
  c("Ta", "SWC", "Ts"),
  c("Ta", "SWC", "Ts", "VPD")
)

for (xs in multiple_specs) {
  xs <- xs[xs %in% names(dat)]
  if (length(xs) < 2) next
  d <- finite_complete(dat, c("phi0", xs))
  if (nrow(d) <= length(xs)) next
  formula_label <- paste("phi0 ~", paste(xs, collapse = " + "))
  fit <- fit_lm_safe(d, stats::as.formula(formula_label))
  if (inherits(fit, "error")) next
  for (x in xs) {
    model_rows[[paste0("multiple_", paste(xs, collapse = "_"), "_", x)]] <-
      tidy_lm_terms(
        fit, "phi0", x, setdiff(xs, x), "multiple regression",
        formula_label, partial = TRUE
      ) %>%
      dplyr::filter(predictor == x)
  }
}

# C. Frisch-Waugh-Lovell residualization: residualize y and x against Ta.
for (x in c("SWC", "Ts", "VPD")) {
  if (!x %in% names(dat)) next
  d <- finite_complete(dat, c("phi0", "Ta", x))
  if (nrow(d) < 4 || length(unique(d[[x]])) < 2) next
  y_fit <- stats::lm(phi0 ~ Ta, data = d)
  x_fit <- stats::lm(stats::as.formula(paste(x, "~ Ta")), data = d)
  fwl_df <- tibble::tibble(
    phi0_resid_Ta = stats::resid(y_fit),
    x_resid_Ta = stats::resid(x_fit)
  )
  fit <- stats::lm(phi0_resid_Ta ~ x_resid_Ta, data = fwl_df)
  model_rows[[paste0("fwl_", x)]] <- tidy_lm_terms(
    fit,
    "residual(phi0 | Ta)",
    "x_resid_Ta",
    paste0("residualized ", x, " against Ta"),
    "FWL residual-on-residual diagnostic",
    paste0("residual(phi0 ~ Ta) ~ residual(", x, " ~ Ta)")
  ) %>%
    dplyr::mutate(predictor = paste0("residual(", x, " | Ta)"))
}

summary_table <- dplyr::bind_rows(model_rows) %>%
  dplyr::mutate(
    interpretation = dplyr::case_when(
      model_type == "current base residual model" ~
        "Base linear association removed before current-style residual panels.",
      grepl("SWC", predictor) & response %in% c("phi0", "residual(phi0 | Ta)") ~
        "Positive association in the checked-in data, but inference is limited by small n.",
      grepl("Ts", predictor) & response %in% c("phi0", "residual(phi0 | Ta)") ~
        "Association is sensitive to model formulation and limited by small n.",
      TRUE ~ "Exploratory diagnostic; interpret as association, not causation."
    )
  )

readr::write_csv(summary_table, file.path(out_dir, "phi0_residual_model_summary.csv"))

# Pairwise Pearson correlation table among environmental covariates.
env_vars <- intersect(c("Ta", "SWC", "Ts", "VPD"), names(dat))
cor_rows <- list()
for (i in seq_along(env_vars)) {
  for (j in seq_along(env_vars)) {
    if (j <= i) next
    a <- env_vars[[i]]
    b <- env_vars[[j]]
    d <- finite_complete(dat, c(a, b))
    if (nrow(d) < 3) next
    ct <- suppressWarnings(stats::cor.test(d[[a]], d[[b]], method = "pearson"))
    cor_rows[[paste(a, b, sep = "_")]] <- tibble::tibble(
      variable_1 = a,
      variable_2 = b,
      n = nrow(d),
      pearson_r = unname(ct$estimate),
      p_value = ct$p.value
    )
  }
}
cor_table <- dplyr::bind_rows(cor_rows)
readr::write_csv(cor_table, file.path(out_dir, "environment_pairwise_correlations.csv"))

# VIF diagnostics for each candidate multiple-regression predictor set.
vif_for_predictors <- function(df, predictors) {
  predictors <- predictors[predictors %in% names(df)]
  d <- finite_complete(df, predictors)
  if (nrow(d) <= length(predictors) || length(predictors) < 2) return(NULL)

  rows <- lapply(predictors, function(x) {
    others <- setdiff(predictors, x)
    fit <- fit_lm_safe(d, stats::as.formula(paste(x, "~", paste(others, collapse = " + "))))
    if (inherits(fit, "error")) {
      return(tibble::tibble(
        predictor_set = paste(predictors, collapse = " + "),
        predictor = x,
        n = nrow(d),
        r_squared_predictor_model = NA_real_,
        vif = NA_real_,
        note = fit$message
      ))
    }
    r2 <- summary(fit)$r.squared
    tibble::tibble(
      predictor_set = paste(predictors, collapse = " + "),
      predictor = x,
      n = nrow(d),
      r_squared_predictor_model = r2,
      vif = ifelse(is.finite(r2) && r2 < 1, 1 / (1 - r2), NA_real_),
      note = ifelse(nrow(d) <= length(predictors) + 1,
                    "very low residual df; VIF is unstable",
                    NA_character_)
    )
  })
  dplyr::bind_rows(rows)
}

vif_table <- dplyr::bind_rows(lapply(multiple_specs, function(xs) vif_for_predictors(dat, xs)))
readr::write_csv(vif_table, file.path(out_dir, "environment_vif_diagnostics.csv"))

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

extract_row <- function(model_type_value, predictor_pattern) {
  summary_table %>%
    dplyr::filter(model_type == model_type_value, grepl(predictor_pattern, predictor)) %>%
    dplyr::slice(1)
}

current_swc <- extract_row("current residual-on-predictor model", "SWC")
current_ts <- extract_row("current residual-on-predictor model", "Ts")
base_ta <- extract_row("current base residual model", "Ta")

report <- c(
  "# phi0 Residual Diagnostic Report",
  "",
  paste0("Input CSV: `", input_csv, "`"),
  paste0("Rows after shared `N >= ", MIN_LIGHT_RESPONSE_POINTS, "` filter: ", nrow(dat)),
  paste0("Complete rows for phi0, Ta, SWC, Ts, VPD: ", nrow(finite_complete(dat, c("phi0", "Ta", "SWC", "Ts", "VPD")))),
  "",
  "## Existing Figure 5 Logic",
  "",
  "- `phi0` is first regressed against `Ta` with ordinary least squares.",
  "- Residuals from `phi0 ~ Ta` are then regressed separately against `SWC`, `Ts`, and other displayed covariates.",
  "- R2 and p-values are from two-sided OLS slope tests reported by `summary(lm(...))`.",
  "- No transformations, weights, robust standard errors, outlier removal, or bootstrap are used in the residual-panel statistics.",
  "- FWL residualization rows are included to verify slope equivalence with multiple regression; formal p-values should be taken from the corresponding multiple-regression term.",
  "",
  "## Current Reproduction",
  "",
  paste0("- `phi0 ~ Ta`: n = ", base_ta$n, ", slope = ", fmt(base_ta$slope, 6),
         ", R2 = ", fmt(base_ta$r_squared, 3), ", p = ", fmt(base_ta$p_value, 3), "."),
  paste0("- `residual(phi0 ~ Ta) ~ SWC`: n = ", current_swc$n, ", slope = ", fmt(current_swc$slope, 6),
         ", R2 = ", fmt(current_swc$r_squared, 3), ", p = ", fmt(current_swc$p_value, 3), "."),
  paste0("- `residual(phi0 ~ Ta) ~ Ts`: n = ", current_ts$n, ", slope = ", fmt(current_ts$slope, 6),
         ", R2 = ", fmt(current_ts$r_squared, 3), ", p = ", fmt(current_ts$p_value, 3), "."),
  "",
  "These checked-in data do not reproduce the manuscript values `R2 = 0.528, p = 0.011` for SWC or `R2 = 0.589, p = 0.006` for Ts.",
  "",
  "## Recommended Manuscript Language",
  "",
  "Methods/statistics:",
  "",
  "Year-specific light-response parameters were estimated for each ENSO x season group by fitting a rectangular-hyperbola light-response model to flux observations. For Figure 5, we examined residual associations among fitted parameters and environmental covariates. For phi0, we first fit an ordinary least-squares model, phi0 ~ Ta, using the Year x ENSO x season parameter table. We then regressed the resulting residuals separately against SWC and Ts using ordinary least squares. Reported R2 and p-values are from two-sided OLS slope tests. Because SWC and Ts were available for only a small subset of fitted parameter rows, these residual analyses were interpreted as exploratory.",
  "",
  "Results:",
  "",
  "Variation in phi0 not explained by the linear association with Ta was positively associated with SWC in the current residual analysis, but this relationship was based on only four SWC-complete fitted parameter rows and should be interpreted cautiously. The corresponding association with Ts was weaker in the current code/data. These results suggest possible secondary soil-related structure in phi0 beyond Ta, but they do not by themselves demonstrate co-limitation or causal control.",
  "",
  "Figure 5 caption:",
  "",
  "Figure 5. Residual relationships among fitted light-response parameters and environmental covariates. Panel a shows residuals from an OLS regression of phi0 against Ta plotted against SWC; panel b shows the same phi0 residuals plotted against Ts. Lines are ordinary least-squares fits to the displayed residuals. Points represent Year x ENSO x season fitted light-response parameter estimates, colored by year and marked by season/ENSO phase. These panels are exploratory because SWC and Ts availability limits the sample size.",
  "",
  "## Output Files",
  "",
  "- `phi0_residual_model_summary.csv`",
  "- `environment_pairwise_correlations.csv`",
  "- `environment_vif_diagnostics.csv`",
  "- `phi0_residual_diagnostic_report.md`"
)

writeLines(report, file.path(out_dir, "phi0_residual_diagnostic_report.md"))

message("Wrote diagnostics to: ", out_dir)
