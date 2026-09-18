#!/usr/bin/env Rscript

# Sensitivity check for the light-response analysis: refit the same rectangular
# hyperbola directly to the raw half-hourly observations, instead of to the
# monthly 48-bin mean diel composites used for the main result.
#
# Why this exists. Section 2.9 fits the pooled ENSO x season light-response
# curves to monthly diel composites so that every half-hour of the day carries
# equal weight regardless of how many days were observed. Averaging before
# fitting a saturating function is not free: by Jensen's inequality the mean of
# a concave function exceeds the function of the mean, so the composite fit is
# expected to report HIGHER high-light uptake than a fit to individual
# half-hours (Falge et al. 2001 make exactly this point about mean-diurnal-
# variation averaging and saturating light response). This script quantifies
# that difference so the manuscript can report it rather than assert it away.
#
# The model, the fitting routine, the bootstrap and the screening thresholds are
# identical to code/05_light_response_curves.R; only the input data differ.
# The functions below are copied verbatim from that script (lines 174-262) and
# must be kept in sync with it.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

# bootstrap resampling is stochastic; seed reported in the output for reproducibility
seed <- 20260808
set.seed(seed)

n_boot <- 500          # matches code/05_light_response_curves.R
q_ref  <- 2000         # high-light reference PAR used instead of the asymptotic Pmax
par_min <- 20          # low-light exclusion stated in section 2.9

flux_csv <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
comp_csv <- file.path(output_path, "season_enso_light_response_P2000_bootstrap_CI.csv")
out_csv  <- file.path(output_path, "light_response_halfhourly_sensitivity_by_enso_season.csv")

# ─────────────────────────────────────────────────────────────────────────────
# functions copied verbatim from code/05_light_response_curves.R

light_response <- function(Q, phi0, Pmax, Rd) {
  ((phi0 * Pmax * Q) / (phi0 * Q + Pmax)) - Rd
}

compute_lcp <- function(phi0, Pmax, Rd) {
  if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0) {
    Rd * Pmax / (phi0 * (Pmax - Rd))
  } else {
    NA_real_
  }
}

compute_fit_summary <- function(fit, dff, q_ref = 2000) {
  cf <- coef(fit)
  phi0 <- cf["phi0"]
  Pmax <- cf["Pmax"]
  Rd   <- cf["Rd"]

  data.frame(
    phi0  = as.numeric(phi0),
    Pmax  = as.numeric(Pmax),
    Rd    = as.numeric(Rd),
    LCP   = as.numeric(compute_lcp(phi0, Pmax, Rd)),
    P2000 = as.numeric(light_response(q_ref, phi0, Pmax, Rd)),
    N     = nrow(dff)
  )
}

filter_for_fit <- function(df) {
  if (!all(c("PAR", "NEE_ok") %in% names(df))) return(df[0, , drop = FALSE])
  out <- df[df$PAR > par_min & is.finite(df$PAR) & is.finite(df$NEE_ok), , drop = FALSE]
  out$FC_pos <- -out$NEE_ok
  out
}

fit_light_curve <- function(df) {
  df <- filter_for_fit(df)

  # basic sanity so nlsLM doesn't crash immediately
  if (is.null(df) || nrow(df) < 5) return(NULL)
  if (length(unique(df$PAR)) < 3 || length(unique(df$FC_pos)) < 3) return(NULL)

  tryCatch({
    fit <- nlsLM(
      FC_pos ~ ((phi0 * Pmax * PAR) / (phi0 * PAR + Pmax)) - Rd,
      data = df,
      start = list(phi0 = 0.03, Pmax = 30, Rd = 5),
      lower = c(phi0 = 0.005, Pmax = 5, Rd = 0),
      upper = c(phi0 = 0.30,  Pmax = 100, Rd = 20),
      control = nls.lm.control(maxiter = 500)
    )
    list(fit = fit, df = df)
  }, error = function(e) NULL)
}

bootstrap_light_ci <- function(df, n_boot = 500, q_ref = 2000) {
  dff <- filter_for_fit(df)
  if (nrow(dff) < 10) return(NULL)

  boot_rows <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    sampled <- dff[sample(seq_len(nrow(dff)), replace = TRUE), , drop = FALSE]
    res <- fit_light_curve(sampled)

    if (!is.null(res)) {
      boot_rows[[i]] <- compute_fit_summary(res$fit, res$df, q_ref = q_ref)
    }
  }

  boot_tbl <- dplyr::bind_rows(boot_rows)
  if (!nrow(boot_tbl)) return(NULL)

  boot_tbl %>%
    summarise(
      phi0_low  = quantile(phi0,  0.025, na.rm = TRUE),
      phi0_high = quantile(phi0,  0.975, na.rm = TRUE),
      Rd_low    = quantile(Rd,    0.025, na.rm = TRUE),
      Rd_high   = quantile(Rd,    0.975, na.rm = TRUE),
      LCP_low   = quantile(LCP,   0.025, na.rm = TRUE),
      LCP_high  = quantile(LCP,   0.975, na.rm = TRUE),
      P2000_low = quantile(P2000, 0.025, na.rm = TRUE),
      P2000_high = quantile(P2000, 0.975, na.rm = TRUE),
      n_boot_success = n()
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# load half-hourly data

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv,
       "\nRun code/01_prepare_flux_timeseries.R first.")
}

required_cols <- c("NEE", "PAR", "season", "ENSO")
df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nExpected them in ", flux_csv)
}

# NEE is the measured, non-gap-filled, u*-filtered series. it is renamed to
# NEE_ok only because filter_for_fit() expects that name -- code/05 does the
# same when nee_source = "measured" (it renames NEE_meas the same way).
hh <- df_raw %>%
  transmute(
    PAR    = as.numeric(PAR),
    NEE_ok = as.numeric(NEE),
    season = tolower(trimws(as.character(season))),
    ENSO   = trimws(as.character(ENSO))
  ) %>%
  filter(is.finite(PAR), is.finite(NEE_ok), PAR > par_min,
         season %in% c("wet", "dry"), !is.na(ENSO))

# label the five regimes exactly as they appear in the supplement table
enso_label <- function(x) {
  xx <- tolower(iconv(x, from = "", to = "ASCII//TRANSLIT"))
  dplyr::case_when(
    grepl("nina", xx)    ~ "La Niña",
    grepl("nino", xx)    ~ "El Niño",
    grepl("neutral", xx) ~ "Neutral",
    TRUE                 ~ NA_character_
  )
}

hh <- hh %>%
  mutate(ENSO_lab = enso_label(ENSO),
         Season   = ifelse(season == "wet", "Wet", "Dry"),
         regime   = paste0(ENSO_lab, "-", tolower(Season))) %>%
  filter(!is.na(ENSO_lab))

message(sprintf("half-hourly observations available for fitting (PAR > %d): %d",
                par_min, nrow(hh)))
print(table(hh$ENSO_lab, hh$Season))

# ─────────────────────────────────────────────────────────────────────────────
# fit each ENSO x season regime

regimes <- hh %>% distinct(Season, ENSO_lab, regime) %>% arrange(Season, ENSO_lab)

rows <- list()
for (i in seq_len(nrow(regimes))) {
  r  <- regimes[i, ]
  sub <- hh %>% filter(Season == r$Season, ENSO_lab == r$ENSO_lab)

  res <- fit_light_curve(sub)
  if (is.null(res)) {
    message(sprintf("  %s-%s: FIT FAILED (n = %d)", r$ENSO_lab, tolower(r$Season), nrow(sub)))
    next
  }

  est <- compute_fit_summary(res$fit, res$df, q_ref = q_ref)
  ci  <- bootstrap_light_ci(sub, n_boot = n_boot, q_ref = q_ref)

  fmt <- function(v, lo, hi, digits) {
    if (is.null(ci) || !is.finite(lo) || !is.finite(hi)) return(sprintf(paste0("%.", digits, "f"), v))
    sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), v, lo, hi)
  }

  rows[[length(rows) + 1]] <- data.frame(
    Season   = r$Season,
    ENSO     = r$ENSO_lab,
    P2000_CI = fmt(est$P2000, ci$P2000_low, ci$P2000_high, 1),
    Rd_CI    = fmt(est$Rd,    ci$Rd_low,    ci$Rd_high,    2),
    phi0_CI  = fmt(est$phi0,  ci$phi0_low,  ci$phi0_high,  3),
    LCP_CI   = fmt(est$LCP,   ci$LCP_low,   ci$LCP_high,   0),
    N        = est$N,
    n_boot_success = if (is.null(ci)) NA_integer_ else ci$n_boot_success,
    P2000_point = round(est$P2000, 3),
    Pmax_point  = round(est$Pmax, 3),
    stringsAsFactors = FALSE
  )
  message(sprintf("  %s-%s: P2000 = %.2f, Pmax = %.2f, n = %d",
                  r$ENSO_lab, tolower(r$Season), est$P2000, est$Pmax, est$N))
}

if (!length(rows)) stop("no regime produced a usable fit")

sens <- bind_rows(rows) %>% mutate(seed = seed, n_boot = n_boot)
readr::write_csv(sens, out_csv)
cat("saved:", out_csv, "\n\n")

print(as.data.frame(sens %>% select(Season, ENSO, P2000_CI, Rd_CI, phi0_CI, LCP_CI, N)),
      row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# compare against the composite fit that produced the main table

if (file.exists(comp_csv)) {
  comp <- read.csv(comp_csv, check.names = FALSE, fileEncoding = "UTF-8")
  num1 <- function(x) suppressWarnings(as.numeric(sub(" .*", "", x)))
  cmp <- comp %>%
    transmute(Season, ENSO, P2000_composite = num1(P2000_CI), N_composite = N) %>%
    left_join(sens %>% transmute(Season, ENSO,
                                 P2000_halfhourly = P2000_point, N_halfhourly = N),
              by = c("Season", "ENSO")) %>%
    mutate(difference = round(P2000_halfhourly - P2000_composite, 2),
           pct = round(100 * (P2000_halfhourly - P2000_composite) / P2000_composite, 1))

  cat("\n=== P2000: composite fit (main result) vs raw half-hourly fit ===\n")
  print(as.data.frame(cmp), row.names = FALSE)

  n_lower <- sum(cmp$difference < 0, na.rm = TRUE)
  cat(sprintf(
    "\nthe raw half-hourly fit gives lower P2000 in %d of %d regimes (mean difference %.2f umol m-2 s-1, %.1f%%)\n",
    n_lower, sum(is.finite(cmp$difference)),
    mean(cmp$difference, na.rm = TRUE), mean(cmp$pct, na.rm = TRUE)))
  cat("Jensen's inequality predicts the composite fit to be the higher of the two.\n")
} else {
  message("composite comparison table not found; run code/05_light_response_curves.R for the side-by-side")
}
