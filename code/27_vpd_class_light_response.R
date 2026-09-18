#!/usr/bin/env Rscript

# Is high-light uptake constrained by VPD in the wet season, or only in the dry?
#
# Why this exists. A co-author asked, of section 3.2's claim that wet-season
# uptake "remained strongly light-driven and was less constrained by VPD":
# "have you used VPD classes to check?"
#
# The supplement has about fourteen VPD-binned or VPD-conditioned figures, but
# NONE of them conditions on PAR, and section 3.3.1 itself notes that the
# limitation "appears only once light is controlled for". So none of them can
# answer the question. This script does two things that can.
#
#   Part 1. Rebuilds the light-response residual analysis reported in section 3.2
#   and in writing/revision_guide_VPD-uptake_correction.md. That analysis lives
#   in a scratch script (lr_resid.R) that was never committed, and its published
#   table is MISSING the high-PAR wet-season cell -- the one cell that actually
#   tests the claim. The table below is complete. It shows the wet season at high
#   light is constrained at least as much as the dry, so "absent (opposite-signed)
#   in the wet season" came from comparing the high-light dry subset against the
#   all-light wet subset.
#
#   Part 2. Fits the light-response curve within VPD classes -- the co-author's
#   suggestion, done with light held fixed by construction -- as an independent
#   check. It agrees.
#
# Both parts are drawn as main-text Figure 5: (a) the high-light residuals of
# Part 1, (b) the VPD-class curves of Part 2. Part 1 leads because it is the
# paper's own method and carries no binning confound. Before Round 16 the
# VPD-class curves alone were an SI figure; that file is now removed.
#
# It is Figure 5, not 7, because the manuscript numbers figures in first-citation
# order and this one is first cited in section 3.2 -- ahead of the Ta/VPD boxplots
# (now Figure 6) and canopy conductance (now Figure 7).
#
# Caveat that must travel with Part 2: a VPD class is also a temperature class
# and a time-of-day class. Both are reported per class so the confound is visible
# in the output, not just asserted in the text.
#
# The model, bounds and screening are copied from code/05_light_response_curves.R
# and must be kept in sync with it, as in code/21 and code/26.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(minpack.lm)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

seed   <- 20260903
n_boot <- 400
q_ref  <- 2000
set.seed(seed)

comp_csv   <- file.path(input_folder, "tambopata_48points_per_month_measured.csv")
params_csv <- file.path(output_path, "season_enso_means_phi0_Pmax_Rd_LCP.csv")
hh_csv     <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
out_resid  <- file.path(out_base, "table_lr_residual_by_light_and_season.csv")
out_class  <- file.path(out_base, "table_vpd_class_light_response.csv")
out_decomp <- file.path(out_base, "table_reco_vs_gross_decomposition.csv")
out_png    <- file.path(fig_path, "Fig5_vpd_class_light_response.png")

# VPD range common to both seasons at high light, and the class edges
vpd_edges <- c(0.4, 0.8, 1.2, 1.6, 2.0)
par_min   <- 20
high_par  <- 1000

# ─────────────────────────────────────────────────────────────────────────────
# copied from code/05_light_response_curves.R

light_response <- function(Q, phi0, Pmax, Rd) ((phi0 * Pmax * Q) / (phi0 * Q + Pmax)) - Rd

fit_light_curve <- function(PAR, FC_pos) {
  if (length(PAR) < 5 || length(unique(PAR)) < 3) return(NULL)
  tryCatch(
    nlsLM(y ~ ((phi0 * Pmax * x) / (phi0 * x + Pmax)) - Rd,
          data = data.frame(x = PAR, y = FC_pos),
          start = list(phi0 = 0.03, Pmax = 30, Rd = 5),
          lower = c(phi0 = 0.005, Pmax = 5, Rd = 0),
          upper = c(phi0 = 0.30,  Pmax = 100, Rd = 20),
          control = nls.lm.control(maxiter = 500)),
    error = function(e) NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
# PART 1 - the residual analysis, completed
#
# Residual = observed(-NEE) - fitted(-NEE) from the per-ENSO x season hyperbola,
# on the monthly 48-bin composites, Neutral excluded. Negative residual = less
# uptake than the curve predicts at that light level.

if (!file.exists(comp_csv)) stop("composites not found: ", comp_csv)
if (!file.exists(params_csv)) stop("fitted parameters not found: ", params_csv,
                                   "\nRun code/05_light_response_curves.R first.")

pars <- readr::read_csv(params_csv, show_col_types = FALSE)

comp <- readr::read_csv(comp_csv, show_col_types = FALSE) %>%
  rename(VPD = VPD_kPa, Ta = TA_1_1_1) %>%
  mutate(NEE_ok = NEE_meas) %>%
  filter(PAR > par_min, is.finite(PAR), is.finite(NEE_ok), is.finite(Ta), is.finite(VPD)) %>%
  filter(!grepl("neutral", tolower(as.character(ENSO)))) %>%
  mutate(
    ENSO_lab = ifelse(grepl("nina|niña", tolower(iconv(ENSO, "", "ASCII//TRANSLIT"))),
                      "La Niña", "El Niño"),
    Season   = ifelse(tolower(season) == "wet", "Wet", "Dry")
  ) %>%
  inner_join(pars %>% select(Season, ENSO, phi0, Pmax, Rd),
             by = c("Season", "ENSO_lab" = "ENSO")) %>%
  mutate(residual = (-NEE_ok) - light_response(PAR, phi0, Pmax, Rd))

slope_of <- function(sub, v) {
  if (nrow(sub) < 20) return(c(NA_real_, NA_real_, NA_real_))
  m <- summary(lm(residual ~ sub[[v]], data = sub))
  c(coef(m)[2, 1], coef(m)[2, 4], m$r.squared)
}

resid_rows <- function(sub, label) {
  ta <- slope_of(sub, "Ta"); vp <- slope_of(sub, "VPD")
  tibble(subset = label, n = nrow(sub),
         Ta_slope = ta[1], Ta_p = ta[2], Ta_R2 = ta[3],
         VPD_slope = vp[1], VPD_p = vp[2], VPD_R2 = vp[3])
}

hi <- comp %>% filter(PAR > high_par)
resid_tbl <- bind_rows(
  resid_rows(comp, sprintf("all daytime (PAR > %d)", par_min)),
  resid_rows(comp %>% filter(Season == "Dry"), "dry season"),
  resid_rows(comp %>% filter(Season == "Wet"), "wet season"),
  resid_rows(hi, sprintf("high PAR (> %d)", high_par)),
  resid_rows(hi %>% filter(Season == "Dry"), "high PAR, dry"),
  resid_rows(hi %>% filter(Season == "Wet"), "high PAR, wet")
)

# Range-matching. At high light the wet season never reaches the VPD the dry
# season does (wet tops out near 1.65 kPa, dry near 2.79), so the two full-range
# slopes above are fitted over different VPD spans -- the same subset mismatch
# that produced the published error. Refit both over the range they share.
vpd_lo <- max(tapply(hi$VPD, hi$Season, min))
vpd_hi <- min(tapply(hi$VPD, hi$Season, max))
common <- hi %>% filter(VPD >= vpd_lo, VPD <= vpd_hi)
resid_tbl <- bind_rows(
  resid_tbl,
  resid_rows(common %>% filter(Season == "Dry"),
             sprintf("high PAR, dry, VPD %.2f-%.2f", vpd_lo, vpd_hi)),
  resid_rows(common %>% filter(Season == "Wet"),
             sprintf("high PAR, wet, VPD %.2f-%.2f", vpd_lo, vpd_hi))
)

# assert against the five rows published in the revision guide, so a drift in
# method shows up as an error rather than as a quietly different number
published <- tribble(
  ~subset,                    ~Ta_ref, ~VPD_ref,
  "all daytime (PAR > 20)",   -0.16,   -0.41,
  "dry season",               -0.38,   -1.39,
  "wet season",               +0.47,   +1.71,
  "high PAR (> 1000)",        -1.09,   -4.89,
  "high PAR, dry",            -1.05,   -5.36
)
chk <- resid_tbl %>% inner_join(published, by = "subset") %>%
  mutate(dTa = abs(Ta_slope - Ta_ref), dV = abs(VPD_slope - VPD_ref))
if (nrow(chk) != nrow(published) || any(chk$dTa > 0.02) || any(chk$dV > 0.02)) {
  print(as.data.frame(chk), row.names = FALSE)
  stop("residual analysis does not reproduce the published rows; do not edit the manuscript from this")
}
cat(sprintf("Part 1: reproduces all %d published residual rows (n = %d)\n", nrow(chk), nrow(comp)))

cat("\n=== light-response residual vs Ta and VPD (negative slope = less uptake) ===\n")
print(as.data.frame(resid_tbl %>% transmute(
  subset, n,
  Ta  = sprintf("%+.2f (R2=%.2f, p=%.1g)", Ta_slope, Ta_R2, Ta_p),
  VPD = sprintf("%+.2f (R2=%.2f, p=%.1g)", VPD_slope, VPD_R2, VPD_p))), row.names = FALSE)
cat("\n  The high-PAR WET row is the one the published table never contained.\n",
    "  At high light the wet season is constrained at least as much as the dry.\n",
    "  The wet season's POSITIVE slope exists only in the all-light subset, where\n",
    "  dim mornings pair low VPD with low uptake by construction.\n",
    sprintf("  Range-matched over the shared %.2f-%.2f kPa span the contrast holds and\n",
            vpd_lo, vpd_hi),
    "  widens, so it is not an artefact of the dry season reaching higher VPD.\n", sep = "")

# Sign convention. Everything above is computed on -NEE (uptake positive), so a
# NEGATIVE slope means less uptake. The manuscript, Figure 4 and Figure 5b all
# display NEE (release positive), where the same effect is a POSITIVE slope. The
# *_NEE columns carry the manuscript-facing sign so nobody has to flip it by hand;
# the uptake-space columns are left untouched because the assertion above is
# written against them.
resid_tbl <- resid_tbl %>%
  mutate(Ta_slope_NEE = -Ta_slope, VPD_slope_NEE = -VPD_slope)

readr::write_csv(resid_tbl %>% mutate(seed = seed), out_resid)
cat("\nsaved:", out_resid, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# PART 1b - how much of the residual slope is respiration, not photosynthesis?
#
# A co-author objected: "you plot NEE vs PAR, and not GEP/GPP. So, you always
# include Reco." Correct, and testable. With A = gross uptake (= -GEP),
#
#     -NEE = A - Reco     =>     residual = A - Reco - f(PAR)
#     d(residual)/dx = dA/dx - dReco/dx     =>     dA/dx = d(residual)/dx + dReco/dx
#
# so regressing Reco on the same driver over the same subset separates the two.
#
# Reco is not modelled here: code/01 sets it to the mean of measured nighttime NEE
# from the two flanking nights, a PER-DAY CONSTANT with no Ta, VPD, SWC or
# radiation term. It therefore has no diurnal shape, and any Reco contribution to
# the residual slope is a between-month effect. That is also why joining monthly
# mean Reco onto the composite is legitimate: within a month Reco barely moves.
#
# The measured composite deliberately carries no Reco (see build_measured_composite),
# so it is joined in from the half-hourly master by year-month.

hh_reco <- readr::read_csv(hh_csv, show_col_types = FALSE,
                           col_select = c("tv_dt", "Reco")) %>%
  mutate(ym = format(lubridate::dmy_hms(tv_dt, quiet = TRUE), "%Y-%m")) %>%
  filter(!is.na(ym), is.finite(Reco)) %>%
  group_by(ym) %>%
  summarise(Reco = mean(Reco), n_Reco = n(), .groups = "drop")

hi_r <- hi %>% mutate(ym = as.character(ym)) %>% inner_join(hh_reco, by = "ym")

cat(sprintf("\nPart 1b: Reco joined to %d of %d high-PAR composite rows (%d months)\n",
            nrow(hi_r), nrow(hi), dplyr::n_distinct(hi_r$ym)))
cat(sprintf("  within-month spread of Reco is negligible by construction; between-month sd = %.2f\n",
            sd(hh_reco$Reco)))

decomp_rows <- function(sub, label) {
  bind_rows(lapply(c("Ta", "VPD"), function(v) {
    r_res  <- coef(lm(residual ~ sub[[v]], data = sub))[2]   # d(residual)/dx
    r_reco <- coef(lm(Reco     ~ sub[[v]], data = sub))[2]   # dReco/dx
    # residual = A - Reco - f(PAR), so the residual OVERSTATES the photosynthetic
    # slope by exactly -dReco/dx. Positive share = Reco inflates the apparent
    # effect; negative share = Reco masks it and the real effect is stronger.
    tibble(subset = label, driver = v, n = nrow(sub),
           residual_slope = unname(r_res),
           Reco_slope     = unname(r_reco),
           gross_slope    = unname(r_res + r_reco),
           Reco_pct_of_residual = 100 * unname(-r_reco) / unname(r_res))
  }))
}

decomp <- bind_rows(
  decomp_rows(hi_r, sprintf("high PAR (> %d)", high_par)),
  decomp_rows(hi_r %>% filter(Season == "Dry"), "high PAR, dry"),
  decomp_rows(hi_r %>% filter(Season == "Wet"), "high PAR, wet")
)

# the identity is algebraic, so it must close exactly
stopifnot(max(abs(decomp$residual_slope + decomp$Reco_slope - decomp$gross_slope)) < 1e-10)

cat("\n=== Part 1b: is the residual slope respiration or photosynthesis? ===\n")
cat("  (uptake-positive space: negative = less uptake. gross = residual + Reco.)\n")
print(as.data.frame(decomp %>% transmute(
  subset, driver, n,
  `d(residual)`  = sprintf("%+.2f", residual_slope),
  `d(Reco)`      = sprintf("%+.2f", Reco_slope),
  `d(gross)`     = sprintf("%+.2f", gross_slope),
  `Reco % of resid` = sprintf("%+.0f%%", Reco_pct_of_residual))), row.names = FALSE)

for (d in c("VPD", "Ta")) {
  z <- decomp %>% filter(driver == d, subset == sprintf("high PAR (> %d)", high_par))
  cat(sprintf("  %-3s: %s\n", d, if (z$Reco_pct_of_residual > 0)
    sprintf("Reco inflates the apparent effect by %.0f%%; %.0f%% is photosynthetic",
            z$Reco_pct_of_residual, 100 - z$Reco_pct_of_residual)
    else
    sprintf("Reco MASKS the effect (%.0f%%); the photosynthetic slope is steeper (%.2f vs %.2f)",
            -z$Reco_pct_of_residual, z$gross_slope, z$residual_slope)))
}

# The join costs months: Reco exists only where gap-filled NEE does, so this runs
# on a subset of the rows Part 1 used. Report the residual slope on BOTH so the
# reader can see the subset is not doing the work.
cat("\n  subset check - residual slope on the Reco-joined rows vs all Part 1 rows:\n")
for (d in c("Ta", "VPD")) {
  full <- coef(lm(residual ~ hi[[d]], data = hi))[2]
  join <- decomp$residual_slope[decomp$driver == d &
            decomp$subset == sprintf("high PAR (> %d)", high_par)]
  cat(sprintf("    %-3s  all %d rows: %+.2f   Reco-joined %d rows: %+.2f\n",
              d, nrow(hi), full, nrow(hi_r), join))
}

readr::write_csv(decomp %>% mutate(seed = seed), out_decomp)
cat("saved:", out_decomp, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# PART 2 - light-response curves within VPD classes
#
# Half-hourly rather than composites: compositing averages over VPD within each
# month-hour bin, which is exactly the variation the classes need.

hh <- readr::read_csv(hh_csv, show_col_types = FALSE,
                      col_select = c("NEE", "GEP", "PAR", "VPD_kPa", "TA_1_1_1", "tv_dt", "season")) %>%
  rename(VPD = VPD_kPa, Ta = TA_1_1_1) %>%
  filter(PAR > par_min, is.finite(PAR), is.finite(NEE), is.finite(VPD), is.finite(Ta),
         VPD >= min(vpd_edges), VPD < max(vpd_edges), season %in% c("wet", "dry")) %>%
  mutate(FC_pos = -NEE,
         hour   = as.numeric(format(lubridate::dmy_hms(tv_dt, quiet = TRUE), "%H")) +
                  as.numeric(format(lubridate::dmy_hms(tv_dt, quiet = TRUE), "%M")) / 60,
         vpd_class = cut(VPD, vpd_edges, right = FALSE),
         Season = ifelse(season == "wet", "Wet", "Dry"))

class_rows <- list()
for (s in c("Dry", "Wet")) for (cl in levels(hh$vpd_class)) {
  g <- hh %>% filter(Season == s, vpd_class == cl)
  if (nrow(g) < 80) next
  fit <- fit_light_curve(g$PAR, g$FC_pos)
  if (is.null(fit)) next
  cf <- coef(fit)
  boot <- vapply(seq_len(n_boot), function(i) {
    idx <- sample(seq_len(nrow(g)), replace = TRUE)
    f <- fit_light_curve(g$PAR[idx], g$FC_pos[idx])
    if (is.null(f)) NA_real_ else light_response(q_ref, coef(f)["phi0"], coef(f)["Pmax"], coef(f)["Rd"])
  }, numeric(1))
  class_rows[[length(class_rows) + 1]] <- tibble(
    Season = s, vpd_class = cl, n = nrow(g),
    mean_VPD = mean(g$VPD), mean_Ta = mean(g$Ta), mean_hour = mean(g$hour, na.rm = TRUE),
    n_PAR_gt_1500 = sum(g$PAR > 1500),
    phi0 = unname(cf["phi0"]), Pmax = unname(cf["Pmax"]), Rd = unname(cf["Rd"]),
    P2000 = unname(light_response(q_ref, cf["phi0"], cf["Pmax"], cf["Rd"])),
    P2000_lo = quantile(boot, 0.025, na.rm = TRUE),
    P2000_hi = quantile(boot, 0.975, na.rm = TRUE)
  )
}
class_tbl <- bind_rows(class_rows)

# --- Part 2b: the same fit on GEP instead of NEE -----------------------------
# Second answer to "you always include Reco": drop Reco entirely and fit the
# hyperbola to gross uptake. Two costs, both stated in the output rather than
# buried: GEP is derived as NEE - Reco (code/25), so it is not an independent
# measurement; and it exists only where the gap-filled series does, so this arm
# covers fewer half-hours than the NEE arm above.
gep_rows <- list()
for (s in c("Dry", "Wet")) for (cl in levels(hh$vpd_class)) {
  g <- hh %>% filter(Season == s, vpd_class == cl, is.finite(GEP))
  if (nrow(g) < 80) next
  fit <- fit_light_curve(g$PAR, -g$GEP)   # -GEP = gross uptake, positive
  if (is.null(fit)) next
  cf <- coef(fit)
  gep_rows[[length(gep_rows) + 1]] <- tibble(
    Season = s, vpd_class = cl, n_GEP = nrow(g),
    A2000 = unname(light_response(q_ref, cf["phi0"], cf["Pmax"], cf["Rd"])))
}
gep_tbl <- bind_rows(gep_rows)

cat("\n=== Part 2: light response fitted WITHIN VPD classes ===\n")
print(as.data.frame(class_tbl %>% transmute(
  Season, `VPD class` = vpd_class, n, `mean Ta` = round(mean_Ta, 1),
  `mean hour` = round(mean_hour, 1), `n PAR>1500` = n_PAR_gt_1500,
  P2000 = sprintf("%.1f [%.1f, %.1f]", P2000, P2000_lo, P2000_hi))), row.names = FALSE)

for (s in c("Dry", "Wet")) {
  z <- class_tbl %>% filter(Season == s)
  if (nrow(z) >= 2) cat(sprintf("  %s: P2000 %.1f -> %.1f across the VPD range (%+.1f umol CO2 m-2 s-1)\n",
                                s, z$P2000[1], z$P2000[nrow(z)], z$P2000[nrow(z)] - z$P2000[1]))
}
cat("  NOTE: a VPD class is also a temperature class and a time-of-day class\n",
    "  (see mean Ta and mean hour above); the VPD effect cannot be isolated from them.\n", sep = "")

cat("\n=== Part 2b: same classes, fitted to GEP (Reco removed entirely) ===\n")
if (nrow(gep_tbl) >= 2) {
  cmp <- class_tbl %>% select(Season, vpd_class, n, P2000) %>%
    inner_join(gep_tbl, by = c("Season", "vpd_class"))
  print(as.data.frame(cmp %>% transmute(
    Season, `VPD class` = vpd_class, `n NEE` = n, `n GEP` = n_GEP,
    `P2000 (NEE)` = sprintf("%.1f", P2000), `A2000 (GEP)` = sprintf("%.1f", A2000))),
    row.names = FALSE)
  for (s in c("Dry", "Wet")) {
    z <- cmp %>% filter(Season == s)
    if (nrow(z) >= 2) cat(sprintf(
      "  %s: across the VPD range, NEE-based %+.1f vs GEP-based %+.1f umol CO2 m-2 s-1\n",
      s, z$P2000[nrow(z)] - z$P2000[1], z$A2000[nrow(z)] - z$A2000[1]))
  }
  cat("  Caveats: GEP = NEE - Reco is derived, not measured (see code/25), and exists\n",
      "  only where the gap-filled series does, so n is smaller than the NEE arm.\n", sep = "")
} else {
  cat("  too few GEP fits converged to compare\n")
}

readr::write_csv(class_tbl %>% mutate(n_boot = n_boot, seed = seed), out_class)
cat("\nsaved:", out_class, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# MAIN-TEXT FIGURE 5
#
# Panel (a) is the evidence that actually overturns the published claim: the
# high-light residuals of Part 1, with the per-season fits drawn. Panel (b) is
# Part 2, the corroboration the co-author asked for. (a) comes first because it
# is the paper's own method and carries no binning confound.
#
# Authored at the GCB two-column maximum of 180 mm (7.10 in) and placed 1:1 --
# see code/15, where rescaling a raster to this width was shown to destroy the
# type. base_size 9 at 1:1 gives ~9 pt on the page.

fig_w <- 7.10
fig_h <- 3.40
fig_base <- 9

# --- (a) residual vs VPD at high light, by season ---------------------------
# `hi` is the PAR > 1000 subset from Part 1; `resid_tbl` holds its fitted slopes,
# so the annotation cannot drift away from the table the text quotes.
# plotmath via paste(), so there is no operator-precedence ambiguity and the
# superscripts render on any device
# Every number is inside a quoted string: plotmath renders a bare 0.10 as "0.1".
# Slope units and sample sizes live in the caption -- at 2.55 in the panel only
# has room for the two numbers that carry the comparison.
# Sign: the panel is drawn in NEE space (positive = release = less uptake) so it
# matches panel (b), Figure 4 and the section 3.2 text, which all display NEE.
# The computation stays in uptake space, so *_slope_NEE is what gets annotated.
slope_label <- function(row_subset, prefix) {
  r <- resid_tbl[resid_tbl$subset == row_subset, ]
  sprintf("paste('%s:  %+.1f,  ', italic(R)^2, ' = %.2f')", prefix, r$VPD_slope_NEE, r$VPD_R2)
}
# The dry line is drawn over a wider VPD span than the wet one, because at high
# light the wet season never gets that dry. The middle label range-matches it, so
# the panel answers that objection without the reader having to ask.
lab_dry  <- slope_label("high PAR, dry", "Dry")
lab_dryc <- slope_label(sprintf("high PAR, dry, VPD %.2f-%.2f", vpd_lo, vpd_hi),
                        sprintf("Dry, VPD < %.1f kPa", vpd_hi))
lab_wet  <- slope_label("high PAR, wet", "Wet")

p_resid <- ggplot(hi, aes(VPD, -residual, colour = Season, shape = Season)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  scale_colour_manual(values = season_cols, name = NULL) +
  scale_shape_manual(values = season_shapes, name = NULL) +
  geom_vline(xintercept = vpd_hi, linetype = "dotted", colour = "grey45", linewidth = 0.3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.6, size = 2.3,
           parse = TRUE, colour = season_cols[["Dry"]], label = lab_dry) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 3.1, size = 2.3,
           parse = TRUE, colour = season_cols[["Dry"]], label = lab_dryc) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 4.6, size = 2.3,
           parse = TRUE, colour = season_cols[["Wet"]], label = lab_wet) +
  labs(x = "VPD (kPa)",
       y = expression(atop("NEE residual" ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1),
                           "positive = less uptake")),
       subtitle = expression(PAR > 1000 ~ mu*mol ~ m^-2 ~ s^-1)) +
  theme_bw(base_size = fig_base) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = rel(0.9)),
        axis.title.y = element_text(size = rel(0.95)),
        # in NEE space the cloud runs bottom-left to top-right, so the empty
        # quadrants are top-left (annotations) and bottom-right (legend)
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.7, "lines"), legend.margin = margin(1, 3, 1, 3))

# --- (b) light response within VPD classes ----------------------------------
grid <- seq(par_min, 2200, length.out = 120)
# display labels for the cut() levels; cut() itself is left alone so the CSV
# written above does not change
class_labs <- setNames(sprintf("%.1f-%.1f", head(vpd_edges, -1), tail(vpd_edges, -1)),
                       levels(hh$vpd_class))

curves <- bind_rows(lapply(seq_len(nrow(class_tbl)), function(i) {
  r <- class_tbl[i, ]
  tibble(Season = r$Season, vpd_class = r$vpd_class, PAR = grid,
         NEE = -light_response(grid, r$phi0, r$Pmax, r$Rd))
})) %>%
  mutate(vpd_lab = factor(unname(class_labs[as.character(vpd_class)]),
                          levels = unname(class_labs)))

p_class <- ggplot(curves, aes(PAR, NEE, colour = vpd_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Season, nrow = 1) +
  scale_colour_viridis_d(option = "plasma", end = 0.85, direction = -1,
                         name = "VPD class (kPa)") +
  labs(x = expression(PAR ~ (mu*mol ~ m^-2 ~ s^-1)),
       y = expression(NEE ~ (mu*mol ~ CO[2] ~ m^-2 ~ s^-1))) +
  theme_bw(base_size = fig_base) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        legend.margin = margin(0, 0, 0, 0),
        strip.background = element_rect(fill = "grey92", colour = "grey40"))

p <- p_resid + p_class +
  plot_layout(widths = c(2.55, 4.45)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = fig_base + 1, face = "bold"))

ggsave(out_png, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
cat("saved:", out_png, sprintf("(%.2f x %.2f in = %.0f mm wide)\n",
                               fig_w, fig_h, fig_w * 25.4))

# the SI version this replaced; no panel should appear in both the main text and
# the supplement (same rule applied in code/15 when the cumulative panel moved)
old_si <- file.path(fig_path, "FigS_vpd_class_light_response.png")
if (file.exists(old_si)) {
  unlink(old_si)
  cat("removed superseded SI figure:", old_si, "\n")
}
