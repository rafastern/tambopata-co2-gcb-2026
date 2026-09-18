# ───────────────────────────────────────────────────────────────────────────────
# yearly LMM/LM analysis + model selection + diagnostics for light-response params
# Adds: nested ANOVA (minimal vs adjusted) using ML
# covariate pool (z-scored): {Ta, PAR, VPD, Precip}
# factor add-ons tested per subset: {none, +Season, +ENSO, +Season+ENSO, +Season*ENSO}
# precipitation preferred column: precip_mm_day_mean (with fallbacks)
# automatically switches to lm() (no random effect) when (1|Year) is not identifiable
# uses AICc when MuMIn is available
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(rlang)
})

# optional namespaces (script works without them, but enables extra outputs)
required_model_pkgs <- c("broom.mixed", "lmerTest")
missing_model_pkgs <- required_model_pkgs[
  !vapply(required_model_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_model_pkgs) > 0) {
  paths_file <- c("paths.R", file.path("code", "paths.R"))
  source(paths_file[file.exists(paths_file)][1])
  skip_dir <- file.path(out_base, "lmm_yearly_params_P2000_skipped")
  dir.create(skip_dir, showWarnings = FALSE, recursive = TRUE)
  skip_note <- c(
    "code/06_lmm_light_response_params.R was skipped.",
    "",
    paste0("Missing required modeling package(s): ", paste(missing_model_pkgs, collapse = ", ")),
    "This script writes LMM/LM diagnostics under output/ and does not generate canonical figures/ files.",
    "Install the missing package(s) and rerun this script to regenerate those diagnostics."
  )
  writeLines(skip_note, file.path(skip_dir, "README.txt"))
  message(paste(skip_note, collapse = "\n"))
  quit(status = 0)
}

dharma_available <- requireNamespace("DHARMa", quietly = TRUE)
if (!dharma_available) {
  message("optional: install.packages('DHARMa') for simulated residual diagnostics")
}
performance_available <- requireNamespace("performance", quietly = TRUE)
if (!performance_available) {
  message("optional: install.packages('performance') for r2_nakagawa(), r2(), and check_model()")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("optional: install.packages('ggplot2')")
}
if (!requireNamespace("see", quietly = TRUE)) {
  message("optional: install.packages('see') for prettier check_model plots")
}
if (!requireNamespace("MuMIn", quietly = TRUE)) {
  message("optional: install.packages('MuMIn') for AICc")
}
if (!requireNamespace("car", quietly = TRUE)) {
  message("optional: install.packages('car') for VIF and type-III Anova on lm()")
}
if (!requireNamespace("influence.ME", quietly = TRUE)) {
  message("optional: install.packages('influence.ME') for LMM influence")
}

# ───────────────────────────────────────────────────────────────────────────────
# paths and inputs
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
ts       <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
out_dir  <- file.path(out_base, paste0("lmm_yearly_params_P2000_", ts))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# use yearly light-response parameter table
yearly_csv <- file.path(data_dir, "light_response_mean_env_and_fit_params_by_year_with_Ts.csv")

if (!file.exists(yearly_csv)) {
  stop("yearly input file not found: ", yearly_csv)
}

yearly_raw <- readr::read_csv(yearly_csv, show_col_types = FALSE)
# guard: one row per year x group (fails loudly if the C1 fan-out duplication recurs)
if (all(c("Year", "Season", "ENSO") %in% names(yearly_raw)))
  stopifnot(!any(duplicated(yearly_raw[, c("Year", "Season", "ENSO")])))

required_cols <- c("Year", "Ta", "PAR", "P2000", "Rd", "phi0", "LCP")

missing_cols <- required_cols[!required_cols %in% names(yearly_raw)]

if (length(missing_cols) > 0) {
  stop(
    "missing required columns in yearly input file: ",
    paste(missing_cols, collapse = ", "),
    "\nFile used: ", yearly_csv,
    "\nRerun the yearly light-response script after adding P2000."
  )
}

# ───────────────────────────────────────────────────────────────────────────────
# load and prep data

if (!"P2000" %in% names(yearly_raw)) {
  stop(
    "P2000 column not found in the yearly input file: ", yearly_csv,
    "\nUpdate the yearly light-response parameter file to include P2000 before running this LMM script."
  )
}

# expected base columns: Year, Ta, PAR, targets {P2000,Rd,phi0,LCP}; optional: N
MIN_LIGHT_RESPONSE_POINTS <- 100

yearly <- yearly_raw %>%
  mutate(Year = factor(Year)) %>%
  { if ("N" %in% names(.)) dplyr::filter(., is.finite(N), N >= MIN_LIGHT_RESPONSE_POINTS) else . } %>%
  tidyr::drop_na(Year, Ta, PAR)

# find VPD and Precip columns (robust to common variants)
.pick_first_present <- function(nms, dat) {
  cand <- nms[nms %in% names(dat)]
  if (length(cand)) cand[1] else NA_character_
}

# vpd: try common naming patterns used in flux / met products
vpd_var <- .pick_first_present(
  c("VPD","vpd","VPD_mean","vpd_mean","VPD_kPa","vpd_kPa","vpd_kpa","VPDkPa","vpdKPa"),
  yearly
)

pr_var <- .pick_first_present(
  c("precip_mm_day_mean","precip_mm_day","Precip","precip","precip_mm_day_sum"),
  yearly
)

# z-standardize helpers
z <- function(x) as.numeric(scale(x))
add_if_exists <- function(dat, src, dst, fun) {
  if (!is.na(src)) dat[[dst]] <- fun(dat[[src]])
  dat
}

# create z_ covariates (always for Ta,PAR; optionally for VPD, Precip)
yearly <- yearly %>%
  mutate(z_Ta  = z(Ta),
         z_PAR = z(PAR))
yearly <- add_if_exists(yearly, vpd_var, "z_VPD", z)
yearly <- add_if_exists(yearly, pr_var,  "z_Precip", z)

# report which covariates are present
present_covs <- c("z_Ta","z_PAR",
                  if ("z_VPD" %in% names(yearly)) "z_VPD",
                  if ("z_Precip" %in% names(yearly)) "z_Precip")
message("covariates available: ", paste(present_covs, collapse = ", "))

# targets to model (keep only those that exist)
targets <- c("P2000","Rd","phi0","LCP")
targets <- targets[targets %in% names(yearly)]
stopifnot(length(targets) > 0)

# optional factor column names — change here if your names differ
season_var <- if ("Season" %in% names(yearly)) "Season" else if ("season" %in% names(yearly)) "season" else NA_character_
enso_var   <- if ("ENSO"   %in% names(yearly)) "ENSO"   else if ("enso_phase" %in% names(yearly)) "enso_phase" else NA_character_
if (!is.na(season_var)) yearly[[season_var]] <- factor(yearly[[season_var]])
if (!is.na(enso_var))   yearly[[enso_var]]   <- factor(yearly[[enso_var]])

# quick sanity check of z-scales
chk <- tibble(
  var  = present_covs,
  mean = vapply(present_covs, function(v) mean(yearly[[v]], na.rm = TRUE), numeric(1)),
  sd   = vapply(present_covs, function(v) sd(  yearly[[v]], na.rm = TRUE), numeric(1))
)
print(chk)

# ───────────────────────────────────────────────────────────────────────────────
# plotting helpers and diagnostics
.save_png <- function(path, expr, width = 1200, height = 900) {
  png(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# dharma residual diagnostics (save main grid + uniformity + dispersion; silence pop-ups)
save_dharma <- function(model, path_png, path_txt, n_sims = 1000) {
  if (!dharma_available) {
    writeLines(
      "DHARMa diagnostics skipped because the optional 'DHARMa' package is not installed.",
      path_txt
    )
    return(invisible(NULL))
  }

  sim <- try(DHARMa::simulateResiduals(model, n = n_sims), silent = TRUE)
  if (inherits(sim, "try-error")) return(invisible(NULL))
  
  base_main <- path_png
  base_noext <- tools::file_path_sans_ext(path_png)
  path_uniform <- paste0(base_noext, "_uniformity.png")
  path_disp    <- paste0(base_noext, "_dispersion.png")
  
  .save_png(base_main, { plot(sim, quantreg = FALSE) })
  .save_png(path_uniform, { DHARMa::testUniformity(sim, plot = TRUE) })
  .save_png(path_disp,    { DHARMa::testDispersion(sim,  plot = TRUE) })
  
  out <- capture.output({
    cat("DHARMa tests\n============\n\n")
    cat("testUniformity (plot saved to ", basename(path_uniform), "):\n", sep = "")
    print(DHARMa::testUniformity(sim, plot = FALSE)); cat("\n")
    cat("testDispersion (plot saved to ", basename(path_disp), "):\n", sep = "")
    print(DHARMa::testDispersion(sim, plot = FALSE)); cat("\n")
  })
  writeLines(out, path_txt)
  invisible(NULL)
}

# performance::check_model() plot (grid)
save_check_model <- function(model, path_png) {
  if (!performance_available) return(invisible(NULL))
  cm <- performance::check_model(model)  # ggplot/see object
  .save_png(path_png, { print(cm) })
}

# type-III anova saver for both lm and lmm
save_typeIII <- function(model, path_csv) {
  if (inherits(model, "lmerMod")) {
    res <- try(stats::anova(model, type = 3), silent = TRUE)
    if (!inherits(res, "try-error")) {
      readr::write_csv(tibble::rownames_to_column(as.data.frame(res), "term"), path_csv)
    }
  } else if (inherits(model, "lm")) {
    if (requireNamespace("car", quietly = TRUE)) {
      res <- try(car::Anova(model, type = 3), silent = TRUE)
      if (!inherits(res, "try-error")) {
        readr::write_csv(tibble::rownames_to_column(as.data.frame(res), "term"), path_csv)
      }
    } else {
      res <- try(stats::anova(model), silent = TRUE)  # type-I fallback
      if (!inherits(res, "try-error")) {
        readr::write_csv(tibble::rownames_to_column(as.data.frame(res), "term"), path_csv)
      }
    }
  }
}

# base diagnostics bundle
run_all_diagnostics <- function(model, dat, resp, out_dir, model_tag = "model") {
  md <- file.path(out_dir, paste0("diagnostics_", resp, "_", model_tag))
  dir.create(md, showWarnings = FALSE, recursive = TRUE)
  is_lmm <- inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")
  
  .save_png(file.path(md, "resid_vs_fitted.png"), {
    plot(fitted(model), resid(model),
         xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
    abline(h = 0, lty = 2)
  })
  .save_png(file.path(md, "plot_model_base.png"), { plot(model) })
  
  .save_png(file.path(md, "qqplot_resid.png"), { qqnorm(resid(model)); qqline(resid(model)) })
  .save_png(file.path(md, "hist_resid.png"), { hist(resid(model), breaks = 30, main = "Histogram of residuals", xlab = "Residuals") })
  .save_png(file.path(md, "acf_resid.png"), { acf(resid(model), main = "ACF of residuals") })
  .save_png(file.path(md, "residuals_vs_index.png"), {
    plot(seq_len(nrow(dat)), resid(model), xlab = "Row index", ylab = "Residuals")
    abline(h = 0, lty = 2)
  })
  
  save_dharma(model, file.path(md, "DHARMa_plots.png"), file.path(md, "DHARMa_tests.txt"))
  save_check_model(model, file.path(md, "check_model_grid.png"))
  
  if (is_lmm) {
    out_re <- capture.output({
      cat("Random effects (ranef):\n"); print(lme4::ranef(model)); cat("\n")
      cat("VarCorr:\n"); print(lme4::VarCorr(model)); cat("\n")
    })
    writeLines(out_re, file.path(md, "random_effects.txt"))
  }
  
  if (requireNamespace("car", quietly = TRUE) && inherits(model, "lm")) {
    tl <- attr(terms(model), "term.labels")
    has_interactions <- length(tl) > 0 && (any(grepl(":", tl, fixed = TRUE)) || any(grepl("\\*", tl)))
    if (!has_interactions && length(tl) >= 2) {
      capture.output(car::vif(model, type = "predictor"),
                     file = file.path(md, "VIF_fixed_effects.txt"))
    } else {
      writeLines("VIF skipped (interactions present or <2 predictors).",
                 file.path(md, "VIF_fixed_effects.txt"))
    }
  }
  
  if (inherits(model, "lm")) {
    .save_png(file.path(md, "cooks_distance.png"), {
      plot(cooks.distance(model), ylab = "Cook's distance", main = "Cook's distance (lm)")
    })
  } else if (inherits(model, "lmerMod") && requireNamespace("influence.ME", quietly = TRUE)) {
    infl <- try(influence.ME::influence(model, obs = TRUE), silent = TRUE)
    if (!inherits(infl, "try-error")) {
      .save_png(file.path(md, "influenceME.png"), { plot(infl) })
      capture.output(summary(infl), file = file.path(md, "influenceME_summary.txt"))
    }
  }
  
  if (inherits(model, "lmerMod") && requireNamespace("MuMIn", quietly = TRUE)) {
    r2 <- MuMIn::r.squaredGLMM(model)
    capture.output(r2, file = file.path(md, "R2_marginal_conditional.txt"))
  } else if (performance_available) {
    r2_lm <- try(performance::r2(model), silent = TRUE)
    if (!inherits(r2_lm, "try-error")) {
      capture.output(r2_lm, file = file.path(md, "R2_lm.txt"))
    }
  } else {
    writeLines(
      "R2 output skipped because the optional 'performance' package is not installed.",
      file.path(md, "R2_skipped.txt")
    )
  }
  
  invisible(md)
}

# ───────────────────────────────────────────────────────────────────────────────
# cap on fixed-effect terms. with only ~12 aggregated observations (ENSO x season x
# year) a Season*ENSO interaction plus several covariates is badly overparameterised
# (conditional R^2 ~ 0.99). per reviewer request ("reduce parameters"), the candidate
# set is restricted to main-effects models with at most this many fixed terms (each
# covariate and each factor counts as one term); no interactions are considered.
MAX_FIXED_TERMS <- 2

# ───────────────────────────────────────────────────────────────────────────────
# formula generator — main-effects subsets (size 0..MAX_FIXED_TERMS) of the unified
# predictor pool {covariates + Season + ENSO}; adds (1|Year) only if estimable
make_all_formulas <- function(resp, has_season, has_enso, cov_pool, season_var, enso_var, has_replication_for_year) {
  add_re <- function(rhs) {
    if (has_replication_for_year) paste(rhs, "+ (1|Year)") else rhs
  }
  # unified main-effect pool: continuous covariates plus categorical factors,
  # each treated as ONE selectable term. no interactions.
  pred_pool <- cov_pool
  if (has_season) pred_pool <- c(pred_pool, season_var)
  if (has_enso)   pred_pool <- c(pred_pool, enso_var)

  term_key <- function(vec) if (length(vec) == 0) "none" else paste(sub("^z_", "", vec), collapse = "+")

  # all subsets of size 0..MAX_FIXED_TERMS (intercept-only through the capped set)
  max_k <- min(MAX_FIXED_TERMS, length(pred_pool))
  subsets <- list(character(0))
  if (length(pred_pool) && max_k >= 1) {
    for (k in seq_len(max_k)) subsets <- c(subsets, combn(pred_pool, k, simplify = FALSE))
  }

  forms <- list()
  for (s in subsets) {
    rhs <- if (length(s)) paste(s, collapse = " + ") else "1"
    nm  <- paste0("TERMS_", term_key(s))
    forms[[nm]] <- as.formula(paste(resp, "~", add_re(rhs)))
  }
  forms
}

# fit list of formulas with ML (or OLS) and return AIC/BIC/AICc table
fit_forms_ml <- function(forms, dat, use_lmm) {
  fits <- lapply(forms, function(f) {
    if (use_lmm) lmerTest::lmer(f, data = dat, REML = FALSE) else lm(f, data = dat)
  })
  aic  <- vapply(fits, AIC,  numeric(1))
  bic  <- vapply(fits, BIC,  numeric(1))
  aicc <- if (requireNamespace("MuMIn", quietly = TRUE)) {
    vapply(fits, MuMIn::AICc, numeric(1))
  } else rep(NA_real_, length(fits))
  tab <- tibble::tibble(
    model   = names(fits),
    formula = vapply(forms, function(x) paste(deparse(x), collapse = " "), character(1)),
    AIC = aic, BIC = bic, AICc = aicc
  ) %>% arrange(dplyr::coalesce(AICc, AIC))
  tab$deltaAICc <- tab$AICc - min(tab$AICc, na.rm = TRUE)
  list(fits = fits, table = tab)
}

# ───────────────────────────────────────────────────────────────────────────────
# NEW: nested ANOVA helper (minimal vs adjusted), using ML, writes CSV + TXT
save_nested_anova_min_vs_adj <- function(resp_dir, resp, dat, use_lmm, rhs_adj) {
  if (rhs_adj %in% c("", "1")) {
    note <- c("nested ANOVA skipped: adjusted model has no fixed covariates (intercept-only).")
    writeLines(note, file.path(resp_dir, "nested_ANOVA_min_vs_adj_note.txt"))
    return(invisible(NULL))
  }
  if (use_lmm) {
    m_min_ml <- lmerTest::lmer(as.formula(paste0(resp, " ~ (1|Year)")), data = dat, REML = FALSE)
    m_cov_ml <- lmerTest::lmer(as.formula(paste0(resp, " ~ ", rhs_adj, " + (1|Year)")), data = dat, REML = FALSE)
  } else {
    m_min_ml <- lm(as.formula(paste0(resp, " ~ 1")), data = dat)
    m_cov_ml <- lm(as.formula(paste0(resp, " ~ ", rhs_adj)), data = dat)
  }
  cmp <- anova(m_min_ml, m_cov_ml)
  cmp_df <- tibble::rownames_to_column(as.data.frame(cmp), "model")
  readr::write_csv(cmp_df, file.path(resp_dir, "nested_ANOVA_min_vs_adj_ml.csv"))
  summary_txt <- capture.output({
    cat("Nested model comparison (minimal vs adjusted) using ",
        if (use_lmm) "LRT (ML, lmer)" else "F-test (lm)", "\n", sep = "")
    print(cmp)
  })
  writeLines(summary_txt, file.path(resp_dir, "nested_ANOVA_min_vs_adj_ml.txt"))
  invisible(NULL)
}

# ───────────────────────────────────────────────────────────────────────────────
# main worker for each response variable
fit_select_and_save_for <- function(resp) {
  resp_dir <- file.path(out_dir, resp)
  dir.create(resp_dir, showWarnings = FALSE, recursive = TRUE)
  
  dat <- yearly %>% dplyr::filter(!is.na(.data[[resp]]))
  if (!nrow(dat)) return(invisible(NULL))
  
  # detect if (1|Year) is estimable: needs >=2 obs in at least some groups
  has_replication_for_year <- ("Year" %in% names(dat)) && any(table(dat$Year) >= 2)
  use_lmm <- has_replication_for_year
  
  # define covariate pool dynamically (Ta + PAR + VPD + Precip)
  cov_pool <- intersect(c("z_Ta","z_PAR","z_VPD","z_Precip"), names(dat))
  
  # baselines: minimal and adjusted that uses *all present* covariates
  rhs_adj <- if (length(cov_pool)) paste(cov_pool, collapse = " + ") else "1"
  if (use_lmm) {
    f_min <- as.formula(paste0(resp, " ~ (1|Year)"))
    f_cov <- as.formula(paste0(resp, " ~ ", rhs_adj, " + (1|Year)"))
    m_min <- lmerTest::lmer(f_min, data = dat, REML = TRUE)
    m_cov <- lmerTest::lmer(f_cov, data = dat, REML = TRUE)
    r2_min <- if (performance_available) try(performance::r2_nakagawa(m_min), silent = TRUE) else "skipped: optional package 'performance' is not installed"
    r2_cov <- if (performance_available) try(performance::r2_nakagawa(m_cov), silent = TRUE) else "skipped: optional package 'performance' is not installed"
  } else {
    f_min <- as.formula(paste0(resp, " ~ 1"))
    f_cov <- as.formula(paste0(resp, " ~ ", rhs_adj))
    m_min <- lm(f_min, data = dat)
    m_cov <- lm(f_cov, data = dat)
    r2_min <- if (performance_available) try(performance::r2(m_min), silent = TRUE) else "skipped: optional package 'performance' is not installed"
    r2_cov <- if (performance_available) try(performance::r2(m_cov), silent = TRUE) else "skipped: optional package 'performance' is not installed"
  }
  
  save_typeIII(m_min, file.path(resp_dir, "typeIII_minimal.csv"))
  save_typeIII(m_cov, file.path(resp_dir, "typeIII_adjusted.csv"))
  writeLines(capture.output(summary(m_min)), file.path(resp_dir, "summary_minimal.txt"))
  writeLines(capture.output(summary(m_cov)), file.path(resp_dir, "summary_adjusted.txt"))
  
  if (inherits(m_cov, "lmerMod")) {
    fe_cov <- broom.mixed::tidy(m_cov, effects = "fixed", conf.int = TRUE)
  } else {
    cf <- broom.mixed::tidy(m_cov, conf.int = TRUE)
    fe_cov <- cf[cf$effect %in% c("fixed","") | is.na(cf$effect), , drop = FALSE]
  }
  readr::write_csv(fe_cov, file.path(resp_dir, "fixed_effects_adjusted.csv"))
  
  out_fit <- capture.output({
    cat("R2 (minimal):\n"); print(r2_min); cat("\n")
    cat("R2 (adjusted):\n"); print(r2_cov); cat("\n")
    cat("Adjusted RHS used:\n", rhs_adj, "\n")
    cat("Estimator used for baseline summaries: ", if (use_lmm) "LMM (REML)" else "LM (OLS)", "\n")
  })
  writeLines(out_fit, file.path(resp_dir, "model_fit_stats_baselines.txt"))
  
  run_all_diagnostics(model = m_min, dat = dat, resp = resp, out_dir = resp_dir, model_tag = "minimal")
  run_all_diagnostics(model = m_cov, dat = dat, resp = resp, out_dir = resp_dir, model_tag = "adjusted")
  
  save_nested_anova_min_vs_adj(resp_dir, resp, dat, use_lmm, rhs_adj)
  
  has_season <- !is.na(season_var) && season_var %in% names(dat) && is.factor(dat[[season_var]])
  has_enso   <- !is.na(enso_var)   && enso_var   %in% names(dat) && is.factor(dat[[enso_var]])
  all_forms  <- make_all_formulas(resp, has_season, has_enso, cov_pool, season_var, enso_var, use_lmm)
  cmp        <- fit_forms_ml(all_forms, dat, use_lmm)
  
  readr::write_csv(cmp$table, file.path(resp_dir, "model_selection_AIC_BIC_ML_allCovSets.csv"))
  
  best_row <- cmp$table %>% arrange(dplyr::coalesce(AICc, AIC), BIC) %>% slice(1)
  best_formula_chr <- best_row$formula
  best_formula <- formula(best_formula_chr)
  
  if (use_lmm) {
    best_reml <- lmerTest::lmer(best_formula, data = dat, REML = TRUE)
  } else {
    best_reml <- lm(best_formula, data = dat)
  }
  writeLines(capture.output(summary(best_reml)), file.path(resp_dir, "summary_best_REML.txt"))
  
  if (inherits(best_reml, "lmerMod")) {
    fe_best <- broom.mixed::tidy(best_reml, effects = "fixed", conf.int = TRUE)
  } else {
    cf <- broom.mixed::tidy(best_reml, conf.int = TRUE)
    fe_best <- cf[cf$effect %in% c("fixed","") | is.na(cf$effect), , drop = FALSE]
  }
  readr::write_csv(fe_best, file.path(resp_dir, "fixed_effects_best_REML.csv"))
  
  save_typeIII(best_reml, file.path(resp_dir, "typeIII_best_REML.csv"))
  
  if (inherits(best_reml, "lmerMod") && requireNamespace("MuMIn", quietly = TRUE)) {
    capture.output(MuMIn::r.squaredGLMM(best_reml),
                   file = file.path(resp_dir, "R2_best_REML.txt"))
  } else if (performance_available) {
    r2_best <- try(performance::r2(best_reml), silent = TRUE)
    if (!inherits(r2_best, "try-error")) {
      capture.output(r2_best, file = file.path(resp_dir, "R2_best_lm.txt"))
    }
  } else {
    writeLines(
      "R2 output skipped because the optional 'performance' package is not installed.",
      file.path(resp_dir, "R2_best_skipped.txt")
    )
  }
  
  run_all_diagnostics(model = best_reml, dat = dat, resp = resp, out_dir = resp_dir, model_tag = "best_REML")
  
  note_lines <- c(
    paste0("best model by ", if (all(!is.na(cmp$table$AICc))) "AICc" else "AIC", " among covariate subsets & factor add-ons:"),
    paste0("Formula: ", best_formula_chr),
    paste0("Refit with ", if (use_lmm) "REML (LMM)" else "OLS (LM)"),
    "See model_selection_AIC_BIC_ML_allCovSets.csv for the full comparison.",
    "Nested ANOVA (minimal vs adjusted) results: nested_ANOVA_min_vs_adj_ml.*"
  )
  writeLines(note_lines, file.path(resp_dir, "MODEL_CHOICE.txt"))
  
  invisible(NULL)
}

# ───────────────────────────────────────────────────────────────────────────────
# run all targets
for (resp in targets) fit_select_and_save_for(resp)

# ───────────────────────────────────────────────────────────────────────────────
# write a simple readme
readme_lines <- c(
  paste0("results saved to: ", out_dir),
  "",
  "per-parameter folders contain:",
  "- baselines: minimal and adjusted (uses all available of {z_Ta, z_PAR, z_VPD, z_Precip})",
  "- nested ANOVA (minimal vs adjusted) using ML: nested_ANOVA_min_vs_adj_ml.csv/.txt",
  "- model selection across covariate subsets and factor add-ons {none, +Season, +ENSO, +Season+ENSO, +Season*ENSO}",
  "- full table: model_selection_AIC_BIC_ML_allCovSets.csv (sorted by AICc if available)",
  "- best model refit (REML for LMM / OLS for LM): summary_best_REML.txt, fixed_effects_best_REML.csv, typeIII_best_REML.csv, R2_*",
  "- diagnostics: residual plots, optional DHARMa, optional performance::check_model grid, and LMM/LM-specific extras",
  "",
  paste0("yearly file used: ", yearly_csv),
  if (!is.na(vpd_var)) paste0("VPD column used: ", vpd_var) else "VPD column not found.",
  if (!is.na(pr_var))  paste0("Precip column used: ", pr_var) else "Precip column not found."
)
writeLines(readme_lines, file.path(out_dir, "README.txt"))

message("done. results saved to: ", out_dir)
