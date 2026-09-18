#!/usr/bin/env Rscript

# Splitting the light-response result into its gross and respiratory parts.
#
# Why this exists. A co-author, on section 3.2's statement that high daytime
# uptake capacity did not translate into net daily uptake, wrote: "True,
# connection of autotrophic R and GEP. Can you construct NEP? This may be
# revealing...."
#
# Taken literally NEP = -NEE, which the paper already reports; a sign flip
# reveals nothing. The substance is the first clause. The fitted rectangular
# hyperbola gives, by construction,
#
#     P2000 = Ag(2000) - Rd,     Ag(2000) = (phi0 * Pmax * 2000)/(phi0 * 2000 + Pmax)
#
# so the reported P2000 is already a NET quantity: gross light-saturated
# assimilation minus the fitted respiration. Reporting Ag separately, and the
# ratio Rd/Ag, separates the part of respiration that scales with photosynthetic
# capacity (the autotrophic link the comment names) from what is left over.
#
# THE CATCH, and the reason for the bootstrap. Pmax and Rd trade off inside a
# rectangular hyperbola: raising Rd pushes the whole curve down and is
# compensated by a higher Pmax. So an across-regime correlation between Ag and Rd
# can be manufactured by the fit itself, exactly as GEP-Reco correlation is
# manufactured by the partitioning identity (see code/25). This script therefore
# refits with the bootstrap draws RETAINED, so the within-regime parameter
# correlation can be compared against the across-regime one, and so Rd/Ag can be
# reported with an interval rather than as a bare point estimate.
#
# The model, the screening thresholds and the fitting routine below are copied
# verbatim from code/05_light_response_curves.R and must be kept in sync with it;
# code/21 does the same for the same reason. The script asserts that its point
# estimates reproduce code/05's saved parameters before reporting anything.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
  library(tidyr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

seed   <- 20260902
n_boot <- 500          # matches code/05
q_ref  <- 2000
set.seed(seed)

in_csv     <- file.path(input_folder, "tambopata_48points_per_month_measured.csv")
params_csv <- file.path(output_path, "season_enso_means_phi0_Pmax_Rd_LCP.csv")
out_csv    <- file.path(out_base, "table_gross_vs_respiration.csv")

# ─────────────────────────────────────────────────────────────────────────────
# functions copied verbatim from code/05_light_response_curves.R

light_response <- function(Q, phi0, Pmax, Rd) {
  ((phi0 * Pmax * Q) / (phi0 * Q + Pmax)) - Rd
}

compute_lcp <- function(phi0, Pmax, Rd) {
  if (is.finite(phi0) && is.finite(Pmax) && is.finite(Rd) && (Pmax - Rd) > 0) {
    Rd * Pmax / (phi0 * (Pmax - Rd))
  } else NA_real_
}

filter_for_fit <- function(df) {
  if (!all(c("PAR", "NEE_ok") %in% names(df))) return(df[0, , drop = FALSE])
  out <- df[df$PAR > 20 & is.finite(df$PAR) & is.finite(df$NEE_ok), , drop = FALSE]
  out$FC_pos <- -out$NEE_ok
  out
}

fit_light_curve <- function(df) {
  df <- filter_for_fit(df)
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

# the decomposition itself: gross assimilation at the reference light level
gross_at <- function(phi0, Pmax, Q = q_ref) (phi0 * Pmax * Q) / (phi0 * Q + Pmax)

# ─────────────────────────────────────────────────────────────────────────────
# load and build the five ENSO x season regimes

if (!file.exists(in_csv)) stop("input not found: ", in_csv, "\nRun code/05 first.")

df_all <- readr::read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(NEE_ok = NEE_meas,
         season = tolower(trimws(as.character(season))),
         ENSO   = trimws(as.character(ENSO)))

enso_label <- function(x) {
  xx <- tolower(iconv(x, from = "", to = "ASCII//TRANSLIT"))
  dplyr::case_when(grepl("nina", xx) ~ "La Niña",
                   grepl("nino", xx) ~ "El Niño",
                   grepl("neutral", xx) ~ "Neutral",
                   TRUE ~ NA_character_)
}

dat <- df_all %>%
  mutate(ENSO_lab = enso_label(ENSO),
         Season   = ifelse(season == "wet", "Wet", "Dry")) %>%
  filter(!is.na(ENSO_lab), season %in% c("wet", "dry"))

regimes <- dat %>% distinct(Season, ENSO_lab) %>% arrange(Season, ENSO_lab)

# ─────────────────────────────────────────────────────────────────────────────
# fit each regime, keeping the joint bootstrap draws

fit_one <- function(Season, ENSO_lab) {
  sub <- dat %>% filter(Season == !!Season, ENSO_lab == !!ENSO_lab)
  res <- fit_light_curve(sub)
  if (is.null(res)) return(NULL)
  cf <- coef(res$fit)
  point <- tibble(
    Season = Season, ENSO = ENSO_lab,
    phi0 = unname(cf["phi0"]), Pmax = unname(cf["Pmax"]), Rd = unname(cf["Rd"]),
    Ag2000 = gross_at(cf["phi0"], cf["Pmax"]),
    P2000 = light_response(q_ref, cf["phi0"], cf["Pmax"], cf["Rd"]),
    LCP = compute_lcp(cf["phi0"], cf["Pmax"], cf["Rd"]),
    N = nrow(res$df)
  ) %>% mutate(Rd_over_Ag = Rd / Ag2000)

  dff <- filter_for_fit(sub)
  draws <- vector("list", n_boot)
  for (i in seq_len(n_boot)) {
    r <- fit_light_curve(dff[sample(seq_len(nrow(dff)), replace = TRUE), , drop = FALSE])
    if (is.null(r)) next
    cb <- coef(r$fit)
    draws[[i]] <- tibble(phi0 = unname(cb["phi0"]), Pmax = unname(cb["Pmax"]),
                         Rd = unname(cb["Rd"]),
                         Ag2000 = gross_at(cb["phi0"], cb["Pmax"]))
  }
  draws <- bind_rows(draws) %>%
    mutate(P2000 = Ag2000 - Rd, Rd_over_Ag = Rd / Ag2000,
           Season = Season, ENSO = ENSO_lab)
  list(point = point, draws = draws)
}

fits <- lapply(seq_len(nrow(regimes)), function(i)
  fit_one(regimes$Season[i], regimes$ENSO_lab[i]))
fits <- Filter(Negate(is.null), fits)

point <- bind_rows(lapply(fits, `[[`, "point"))
draws <- bind_rows(lapply(fits, `[[`, "draws"))

# ─────────────────────────────────────────────────────────────────────────────
# assert the point estimates reproduce code/05's saved parameters

if (file.exists(params_csv)) {
  ref <- readr::read_csv(params_csv, show_col_types = FALSE) %>%
    transmute(Season, ENSO, phi0_ref = phi0, Pmax_ref = Pmax, Rd_ref = Rd)
  chk <- point %>% inner_join(ref, by = c("Season", "ENSO")) %>%
    mutate(dphi = abs(phi0 - phi0_ref), dP = abs(Pmax - Pmax_ref), dR = abs(Rd - Rd_ref))
  bad <- chk %>% filter(dphi > 0.001 | dP > 0.1 | dR > 0.02)
  if (nrow(bad)) {
    print(as.data.frame(bad %>% select(Season, ENSO, phi0, phi0_ref, Pmax, Pmax_ref, Rd, Rd_ref)),
          row.names = FALSE)
    stop("refit does not reproduce ", basename(params_csv),
         "; the decomposition below would not describe the published curves")
  }
  cat("refit reproduces", basename(params_csv), "for all", nrow(chk), "regimes\n")
} else {
  warning("code/05 parameter file not found; cannot verify the refit", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. the decomposition

summary_tbl <- point %>%
  left_join(
    draws %>% group_by(Season, ENSO) %>%
      summarise(Ag_lo = quantile(Ag2000, 0.025), Ag_hi = quantile(Ag2000, 0.975),
                ratio_lo = quantile(Rd_over_Ag, 0.025),
                ratio_hi = quantile(Rd_over_Ag, 0.975),
                within_corr_Ag_Rd = cor(Ag2000, Rd),
                n_draws = n(), .groups = "drop"),
    by = c("Season", "ENSO")
  ) %>%
  arrange(Rd_over_Ag)

cat("\n1) P2000 = Ag(2000) - Rd, i.e. the reported value is already net of respiration\n")
print(as.data.frame(
  summary_tbl %>% transmute(
    Season, ENSO,
    Ag2000 = sprintf("%.1f [%.1f, %.1f]", Ag2000, Ag_lo, Ag_hi),
    Rd = round(Rd, 2), P2000 = round(P2000, 1),
    `Rd/Ag` = sprintf("%.3f [%.3f, %.3f]", Rd_over_Ag, ratio_lo, ratio_hi))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2. is the Ag-Rd association ecology, or the fit's own parameter covariance?

across <- cor(point$Ag2000, point$Rd)
cat(sprintf("\n2) association between gross capacity and respiration\n"))
cat(sprintf("   across the %d regimes : corr(Ag2000, Rd) = %+.2f\n", nrow(point), across))
cat("   within each regime, from the bootstrap draws (the fit's own parameter covariance):\n")
print(as.data.frame(
  summary_tbl %>% transmute(Season, ENSO, within_corr = round(within_corr_Ag_Rd, 2), n_draws)
), row.names = FALSE)
wmax <- max(abs(summary_tbl$within_corr_Ag_Rd))
cat(sprintf("   largest |within-regime correlation| = %.2f\n", wmax))
if (wmax >= 0.8) {
  cat("   WARNING: the fit generates nearly as much Ag-Rd correlation on its own as is seen\n",
      "  across regimes. Report Rd/Ag with its interval and do not claim the ranking.\n", sep = "")
} else {
  cat("   the across-regime association is substantially stronger than the fit's internal\n",
      "  covariance, so it is not simply a parameter artefact.\n", sep = "")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. does the Rd/Ag ranking separate the ENSO phases?
#
# Pairwise: fraction of bootstrap draws in which one regime's ratio exceeds the
# other's. Draws are independent across regimes, so this is a simple comparison
# of the two bootstrap distributions.

pairs <- expand.grid(a = seq_len(nrow(point)), b = seq_len(nrow(point))) %>%
  filter(a < b) %>%
  rowwise() %>%
  mutate(
    A = paste(point$Season[a], point$ENSO[a]),
    B = paste(point$Season[b], point$ENSO[b]),
    p_A_gt_B = {
      da <- draws %>% filter(Season == point$Season[a], ENSO == point$ENSO[a]) %>% pull(Rd_over_Ag)
      db <- draws %>% filter(Season == point$Season[b], ENSO == point$ENSO[b]) %>% pull(Rd_over_Ag)
      n <- min(length(da), length(db))
      mean(sample(da, n) > sample(db, n))
    }
  ) %>% ungroup() %>% select(A, B, p_A_gt_B)

cat("\n3) Rd/Ag, pairwise probability that the first regime exceeds the second\n")
print(as.data.frame(pairs %>% mutate(p_A_gt_B = round(p_A_gt_B, 2))), row.names = FALSE)

readr::write_csv(
  summary_tbl %>% mutate(across_regime_corr_Ag_Rd = across, n_boot = n_boot, seed = seed),
  out_csv)
cat("\nsaved:", out_csv, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# the numbers section 3.2 quotes

cat("\n=== section 3.2 numbers ===\n")
for (i in seq_len(nrow(summary_tbl))) with(summary_tbl[i, ], cat(sprintf(
  "%-3s %-8s gross Ag(2000) %.1f [%.1f, %.1f], Rd %.2f, net P2000 %.1f, Rd/Ag %.3f [%.3f, %.3f]\n",
  Season, ENSO, Ag2000, Ag_lo, Ag_hi, Rd, P2000, Rd_over_Ag, ratio_lo, ratio_hi)))
cat(sprintf("across-regime corr(Ag, Rd) = %+.2f ; largest within-regime = %+.2f\n", across, wmax))
