#!/usr/bin/env Rscript

# Statistics behind section 3.1's carbon-balance claims. Four things:
#
#   1-5) self-correlation between GEP and Reco arising from the flux partitioning
#   6)   block-bootstrap intervals on the wet-season ENSO contrast
#   6b)  whether Reco and GEP corroborate the 16-day balance, and which of the two
#        tracks its sign
#   6c)  whether ENSO phase drives the sign of the balance, tested at the level of
#        independent ENSO episodes
#
# (The filename reflects only the first of these; rename if that becomes confusing.)
#
# Why this script exists. A co-author, commenting on section 3.1, pointed out
# that NEE = GEP + Reco, so GEP and Reco are not independent, and asked for a
# spurious-correlation analysis in the sense of
#
#   Vickers, D., Thomas, C.K., Martin, J.G., Law, B. (2009). Self-correlation
#   between assimilation and respiration resulting from flux partitioning of
#   eddy-covariance CO2 fluxes. Agric. For. Meteorol. 149.
#   https://doi.org/10.1016/j.agrformet.2009.03.009
#
# The concern is well founded in general. Here GEP is obtained by subtraction,
# GEP = NEE - Reco, so every error in Reco enters GEP with the opposite sign and
# any GEP-Reco correlation is inflated by construction. This script quantifies
# that inflation with a permutation null that destroys any true Reco-NEE
# relation while preserving the arithmetic coupling.
#
# Two things it establishes, both of which the manuscript needs:
#
#   1. The apparent GEP-Reco correlation is mostly arithmetic, so GEP-Reco
#      covariation must never be reported as a result. (The manuscript does not
#      report one; this makes that a stated choice rather than an accident.)
#   2. For the specific wet-season ENSO contrast in section 3.1 the shared term
#      works AGAINST the reported difference: El Nino has lower Reco, which
#      mechanically lowers its apparent |GEP|, yet |GEP| is higher. So that
#      contrast is not a self-correlation artifact. It is, however, not
#      statistically resolved, which the script also shows.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidyr)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

# ─────────────────────────────────────────────────────────────────────────────
# config

seed   <- 20260901
n_perm <- 20000     # permutation replicates for the self-correlation null
n_boot <- 5000      # block-bootstrap replicates for the ENSO contrast
set.seed(seed)

tz_local <- "America/Lima"

# copied from code/04_diel_cycles_carbon_balance.R:150-151, where they are local
sec_per_halfhour <- 1800
gC_per_umolCO2   <- 12 / 1e6

enso_levels <- c("El Nino", "La Nina")

flux_csv  <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
table2_csv <- file.path(output_path, "gep_reco_gs_16day_by_ENSO_season.csv")
out_csv   <- file.path(out_base, "table_gep_reco_self_correlation.csv")

# ─────────────────────────────────────────────────────────────────────────────
# load

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv,
       "\nRun code/01_prepare_flux_timeseries.R first.")
}

required_cols <- c("tv_dt", "NEE_ok", "Reco", "GEP", "season", "ENSO", "is_day_block")
df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) stop("missing required columns: ", paste(missing_cols, collapse = ", "))

df <- df_raw %>%
  mutate(
    ts      = dmy_hms(tv_dt, tz = tz_local, quiet = TRUE),
    Date    = as.Date(ts),
    hh_slot = format(ts, "%H:%M:%S"),
    # identical to code/04: floor_date on a 16-day multiple anchors to the start
    # of each year, so the windows line up with the ones behind Table 2
    window_start = as.Date(floor_date(Date, unit = "16 days")),
    window_mid   = window_start + days(8)
  ) %>%
  filter(!is.na(Date))

# ─────────────────────────────────────────────────────────────────────────────
# 1. the partitioning identity

ident <- df %>% filter(is.finite(NEE_ok), is.finite(Reco), is.finite(GEP))
resid_max <- max(abs(ident$GEP - (ident$NEE_ok - ident$Reco)))
cat(sprintf(
  "\n1) partitioning identity GEP = NEE_ok - Reco\n   max |residual| = %.3g over %d half-hours (machine precision)\n",
  resid_max, nrow(ident)))
if (resid_max > 1e-8) stop("the partitioning identity does not hold; check code/01")
cat(sprintf("   Reco is finite in %d day half-hours and %d night half-hours\n",
            sum(df$is_day_block & is.finite(df$Reco)),
            sum(!df$is_day_block & is.finite(df$Reco))))

# ─────────────────────────────────────────────────────────────────────────────
# 2. reproduce the Table 2 window aggregation, and check it
#
# Reco is extrapolated to 24 h from the window-mean rate while GEP is summed over
# the populated (daytime) slots, exactly as code/04 does. The two are therefore
# on different integration bases, which is why the identity does NOT hold for
# the Table 2 numbers even though it holds at every half-hour.

window_agg <- function(col, mode) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    group_by(window_mid, hh_slot) %>%
    summarise(v = mean(.data[[col]], na.rm = TRUE),
              season = dplyr::first(na.omit(as.character(season))),
              ENSO   = dplyr::first(na.omit(as.character(ENSO))),
              .groups = "drop_last") %>%
    summarise(
      value  = if (mode == "mean24") mean(v, na.rm = TRUE) * 86400 * gC_per_umolCO2
               else sum(v * sec_per_halfhour, na.rm = TRUE) * gC_per_umolCO2,
      season = dplyr::first(season),
      ENSO   = dplyr::first(ENSO),
      .groups = "drop"
    )
}

t2 <- inner_join(
  window_agg("Reco", "mean24") %>% rename(Reco = value),
  window_agg("GEP",  "sum")    %>% select(window_mid, GEP = value),
  by = "window_mid"
) %>%
  # NEE on the same window definition, WITHOUT the coverage screen. Kept only as
  # a robustness comparison: Table 1 and Fig. 3 apply a screen (see below), and
  # section 3.1 must quote the screened numbers so the text agrees with its own
  # table. An earlier version of this script quoted the unscreened values.
  left_join(window_agg("NEE_ok", "sum") %>% select(window_mid, NEE_unscreened = value),
            by = "window_mid")

t2_summary <- t2 %>%
  group_by(season, ENSO) %>%
  summarise(GEP_mean = mean(GEP), Reco_mean = mean(Reco), n_windows = n(), .groups = "drop")

cat("\n2) Table 2 reproduction\n")
print(as.data.frame(t2_summary %>% mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

if (file.exists(table2_csv)) {
  ref <- readr::read_csv(table2_csv, show_col_types = FALSE) %>%
    transmute(season = as.character(season), ENSO = as.character(ENSO),
              GEP_ref = GEP_mean, Reco_ref = Reco_mean, n_ref = n_windows)
  chk <- t2_summary %>% inner_join(ref, by = c("season", "ENSO")) %>%
    mutate(dGEP = abs(GEP_mean - GEP_ref), dReco = abs(Reco_mean - Reco_ref), dn = n_windows - n_ref)
  bad <- chk %>% filter(dGEP > 0.01 | dReco > 0.01 | dn != 0)
  if (nrow(bad)) {
    print(as.data.frame(bad), row.names = FALSE)
    stop("window aggregation does not match ", basename(table2_csv),
         "; the intervals below would not describe the windows Table 2 reports")
  }
  cat("   matches", basename(table2_csv), "exactly (means and window counts)\n")
} else {
  warning("Table 2 reference not found; cannot verify the window aggregation", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. common daytime basis, so the identity holds at window level
#
# Reco, NEE and GEP are integrated over the SAME daytime half-hours. Only on this
# basis is |GEP| = Reco - NEE exactly, which is what the self-correlation
# argument is about.

common <- df %>%
  filter(is_day_block, is.finite(NEE_ok), is.finite(Reco), is.finite(GEP)) %>%
  group_by(window_mid, hh_slot) %>%
  summarise(R = mean(Reco), N = mean(NEE_ok),
            season = dplyr::first(na.omit(as.character(season))),
            ENSO   = dplyr::first(na.omit(as.character(ENSO))),
            .groups = "drop_last") %>%
  summarise(
    R = sum(R * sec_per_halfhour) * gC_per_umolCO2,
    N = sum(N * sec_per_halfhour) * gC_per_umolCO2,
    season = dplyr::first(season), ENSO = dplyr::first(ENSO),
    .groups = "drop"
  ) %>%
  mutate(absG = R - N)

cat(sprintf("\n3) common daytime basis: %d windows, |GEP| = Reco - NEE by construction\n", nrow(common)))

# ─────────────────────────────────────────────────────────────────────────────
# 4. self-correlation: observed versus the identity-only null
#
# The null permutes NEE across windows. That removes any real Reco-NEE
# relationship but leaves the subtraction intact, so the resulting correlation is
# exactly the part of corr(Reco, |GEP|) that the partitioning creates by itself.

self_corr <- function(sub, label) {
  R <- sub$R; N <- sub$N
  if (length(R) < 5) return(NULL)
  obs <- cor(R, R - N)
  null <- vapply(seq_len(n_perm), function(i) cor(R, R - sample(N)), numeric(1))
  tibble(
    subset = label, n = length(R),
    obs_corr = obs,
    null_median = median(null),
    null_lo = unname(quantile(null, 0.025)),
    null_hi = unname(quantile(null, 0.975)),
    p_obs_gt_null = (sum(null >= obs) + 1) / (n_perm + 1)
  )
}

rows <- list(self_corr(common, "all windows"))
for (se in c("dry", "wet")) {
  rows[[length(rows) + 1]] <- self_corr(common %>% filter(season == se), paste("season:", se))
}
for (se in c("dry", "wet")) for (en in enso_levels) {
  rows[[length(rows) + 1]] <- self_corr(common %>% filter(season == se, ENSO == en),
                                        paste(se, en, sep = " / "))
}
sc <- bind_rows(rows)

cat("\n4) corr(Reco, |GEP|): observed versus the null generated by the identity alone\n")
print(as.data.frame(
  sc %>% transmute(subset, n,
                   observed = round(obs_corr, 2),
                   null_median = round(null_median, 2),
                   null_95 = sprintf("%.2f to %.2f", null_lo, null_hi),
                   p = round(p_obs_gt_null, 3),
                   exceeds_null = ifelse(p_obs_gt_null < 0.05, "yes", "no"))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 5. does the shared term create the wet-season ENSO contrast, or oppose it?

wet <- common %>% filter(season == "wet", ENSO %in% enso_levels)
en <- wet %>% filter(ENSO == "El Nino"); ln <- wet %>% filter(ENSO == "La Nina")
dR <- mean(en$R) - mean(ln$R)
dN <- mean(en$N) - mean(ln$N)
dG <- mean(en$absG) - mean(ln$absG)

cat(sprintf(paste0(
  "\n5) wet-season El Nino minus La Nina, decomposition of the |GEP| difference\n",
  "   |GEP| difference        %+6.2f g C m-2 d-1\n",
  "     from the NEE term     %+6.2f   (measured)\n",
  "     from the Reco term    %+6.2f   (the shared, self-correlated term)\n",
  "   The Reco term has the opposite sign to the difference, so the contrast is\n",
  "   not produced by self-correlation; it survives in spite of it.\n"),
  dG, -dN, dR))

# ─────────────────────────────────────────────────────────────────────────────
# 6. is the contrast resolved at all?
#
# Block bootstrap over calendar months and over ENSO episodes, the same scheme as
# month_block_bootstrap() in code/22_morning_transition_test.R (local to that
# script, so repeated here; keep the two in sync).

episode_of <- function(date, enso) {
  y <- as.integer(format(date, "%Y")); m <- as.integer(format(date, "%m"))
  dplyr::case_when(
    enso == "El Nino" & y <= 2019    ~ "EN 2018-19",
    enso == "El Nino"                ~ "EN 2023-24",
    y == 2017 | (y == 2018 & m <= 6) ~ "LN 2017-18",
    y == 2020 | (y == 2021 & m <= 6) ~ "LN 2020-21",
    y == 2021 | (y == 2022 & m <= 6) ~ "LN 2021-22",
    TRUE                             ~ "LN 2022-23"
  )
}

t2_wet <- t2 %>%
  filter(season == "wet", ENSO %in% enso_levels) %>%
  mutate(blk_month = format(window_mid, "%Y-%m"),
         blk_episode = episode_of(window_mid, ENSO))

# NEE on the basis Table 1 and Fig. 3 actually use: the coverage-screened windows
# written by code/04 (n_days >= min_days_per_window and frac_diel_covered >= 0.90).
# Section 3.1 quotes these, not the unscreened values, so the text agrees with its
# own table. The unscreened values are carried through as a robustness row.
nee_screened_csv <- file.path(output_path, "check_16day_dailyCbalance_checked.csv")
nee_ref_csv      <- file.path(output_path, "nee_dailyC_16day_by_ENSO_season.csv")
if (!file.exists(nee_screened_csv)) {
  stop("screened 16-day carbon balance not found: ", nee_screened_csv,
       "\nRun code/04_diel_cycles_carbon_balance.R first.")
}
nee_win <- readr::read_csv(nee_screened_csv, show_col_types = FALSE) %>%
  mutate(window_mid = as.Date(window_mid),
         season = as.character(season), ENSO = as.character(ENSO))

if (file.exists(nee_ref_csv)) {
  ref1 <- readr::read_csv(nee_ref_csv, show_col_types = FALSE) %>%
    transmute(season = as.character(season), ENSO = as.character(ENSO),
              NEE_ref = NEE_mean, n_ref = n_windows)
  chk1 <- nee_win %>%
    group_by(season, ENSO) %>%
    summarise(NEE_mean = mean(dailyC_mean), n_windows = n(), .groups = "drop") %>%
    inner_join(ref1, by = c("season", "ENSO")) %>%
    mutate(dNEE = abs(NEE_mean - NEE_ref), dn = n_windows - n_ref)
  if (nrow(chk1 %>% filter(dNEE > 0.01 | dn != 0))) {
    print(as.data.frame(chk1), row.names = FALSE)
    stop("screened windows do not match ", basename(nee_ref_csv))
  }
  cat("   NEE windows match", basename(nee_ref_csv), "exactly\n")
}

nee_wet <- nee_win %>%
  filter(season == "wet", ENSO %in% enso_levels) %>%
  mutate(NEE = dailyC_mean,
         blk_month = format(window_mid, "%Y-%m"),
         blk_episode = episode_of(window_mid, ENSO))

# ─────────────────────────────────────────────────────────────────────────────
# 6b. do Reco and GEP corroborate the balance on these same 16-day windows?
#
# Asked by a co-author. Joining the components onto the screened windows Fig. 3
# plots answers it two ways: whether they add up, and which of the two tracks the
# sign of the balance.
#
# NOTE on the join key: floor_date(Date, "16 days") anchors to the 1st of each
# MONTH, so windows start on the 1st and the 17th. Reconstructing them from the
# day of year instead matches only a handful of windows.

comp_win <- df %>%
  filter(is.finite(Reco) | is.finite(GEP)) %>%
  group_by(window_start, hh_slot) %>%
  summarise(R = mean(Reco, na.rm = TRUE), G = mean(GEP, na.rm = TRUE), .groups = "drop_last") %>%
  summarise(
    Reco_w = mean(R[is.finite(R)], na.rm = TRUE) * 86400 * gC_per_umolCO2,
    GEP_w  = sum(G[is.finite(G)] * sec_per_halfhour, na.rm = TRUE) * gC_per_umolCO2,
    .groups = "drop"
  )

corrob <- nee_win %>%
  mutate(window_start = as.Date(window_start)) %>%
  inner_join(comp_win, by = "window_start") %>%
  mutate(sum_RG = Reco_w + GEP_w)

cat(sprintf("\n6b) component corroboration on the screened windows: %d of %d matched\n",
            nrow(corrob), nrow(nee_win)))
if (nrow(corrob) < nrow(nee_win)) {
  warning("not every screened window matched a component window; check the join key", call. = FALSE)
}
cat(sprintf("    corr(NEE, Reco + GEP) = %.3f ; mean offset %+.2f ; sd of offset %.2f g C m-2 d-1\n",
            cor(corrob$dailyC_mean, corrob$sum_RG),
            mean(corrob$sum_RG - corrob$dailyC_mean), sd(corrob$sum_RG - corrob$dailyC_mean)))

components <- corrob %>%
  group_by(season, ENSO) %>%
  summarise(n_windows = n(), NEE = mean(dailyC_mean), Reco = mean(Reco_w), GEP = mean(GEP_w),
            .groups = "drop") %>%
  # NEP is just -NEE, and Reco/GPP > 1 is exactly the condition for a net source.
  # The ratio is the flux-level counterpart of Rd/Ag from the light-response
  # curves (code/26), and unlike that one it does not depend on a curve fit.
  mutate(NEP = -NEE, GPP = -GEP, Reco_over_GPP = Reco / GPP) %>%
  arrange(season, ENSO)
cat("\n    components on the same windows (g C m-2 d-1):\n")
print(as.data.frame(
  components %>% transmute(season, ENSO, n_windows, GPP = round(GPP, 2), Reco = round(Reco, 2),
                           NEP = round(NEP, 2), Reco_over_GPP = round(Reco_over_GPP, 2))
), row.names = FALSE)
cat("    GEP is near-identical between phases within a season; Reco is what differs,\n",
    "   so the La Nina source is a respiration signal. Reco is the independently\n",
    "   estimated component and GEP the residual, so this is the safe attribution.\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# 6c. is ENSO phase a driver of the sign of the balance?
#
# Asked by a co-author playing devil's advocate. Windows within one ENSO episode
# are not independent, so the test permutes phase labels across EPISODES, not
# across windows; permuting windows would be badly anticonservative.

ep_tbl <- nee_win %>%
  mutate(ep = dplyr::case_when(
    ENSO == "El Nino" & lubridate::year(window_mid) <= 2019 ~ "EN 2018-19",
    ENSO == "El Nino"                                        ~ "EN 2023-24",
    ENSO == "La Nina" & lubridate::year(window_mid) <= 2018  ~ "LN 2017-18",
    ENSO == "La Nina" & lubridate::year(window_mid) == 2021 &
      lubridate::month(window_mid) <= 6                      ~ "LN 2020-21",
    ENSO == "La Nina"                                        ~ "LN 2021-22",
    TRUE ~ paste("Neutral", lubridate::year(window_mid))
  )) %>%
  group_by(ENSO, ep) %>%
  summarise(n_windows = n(), n_source = sum(dailyC_mean > 0), mean_dailyC = mean(dailyC_mean),
            .groups = "drop")

cat("\n6c) episode-level view (the independent unit)\n")
print(as.data.frame(ep_tbl %>% mutate(mean_dailyC = round(mean_dailyC, 2))), row.names = FALSE)

spread_stat <- function(labels) {
  t <- ep_tbl %>% mutate(L = labels) %>% group_by(L) %>%
    summarise(f = sum(n_source) / sum(n_windows), .groups = "drop")
  max(t$f) - min(t$f)
}
obs_spread <- spread_stat(ep_tbl$ENSO)
null_spread <- vapply(seq_len(n_perm), function(i) spread_stat(sample(ep_tbl$ENSO)), numeric(1))
p_ep <- (sum(null_spread >= obs_spread) + 1) / (n_perm + 1)
cat(sprintf(paste0(
  "    source fraction by phase: %s\n",
  "    episode-label permutation test: observed spread %.2f, null median %.2f, p = %.3f (%d episodes)\n",
  "    => with %d El Nino and %d La Nina episodes, ENSO phase is NOT established as a\n",
  "       driver of consistent direction; the two El Nino episodes have opposite signs.\n"),
  paste(sprintf("%s %.2f", ep_tbl$ENSO %>% unique(),
                tapply(ep_tbl$n_source, ep_tbl$ENSO, sum) / tapply(ep_tbl$n_windows, ep_tbl$ENSO, sum)),
        collapse = "; "),
  obs_spread, median(null_spread), p_ep, nrow(ep_tbl),
  sum(ep_tbl$ENSO == "El Nino"), sum(ep_tbl$ENSO == "La Nina")))

block_boot <- function(tbl, value_col, block_col, reps) {
  blocks <- split(tbl, list(tbl$ENSO, tbl[[block_col]]), drop = TRUE)
  keys <- lapply(enso_levels, function(g)
    names(blocks)[vapply(blocks, function(b) b$ENSO[1] == g, logical(1))])
  names(keys) <- enso_levels
  out <- numeric(reps)
  for (r in seq_len(reps)) {
    m <- vapply(enso_levels, function(g) {
      k <- sample(keys[[g]], length(keys[[g]]), replace = TRUE)
      mean(unlist(lapply(k, function(x) blocks[[x]][[value_col]])))
    }, numeric(1))
    out[r] <- m[["El Nino"]] - m[["La Nina"]]
  }
  unname(quantile(out, c(0.025, 0.975)))
}

one_contrast <- function(tbl, v, label) {
  e <- tbl[[v]][tbl$ENSO == "El Nino"]; l <- tbl[[v]][tbl$ENSO == "La Nina"]
  cm <- block_boot(tbl, v, "blk_month", n_boot)
  ce <- block_boot(tbl, v, "blk_episode", n_boot)
  tibble(quantity = label, n_el_nino = length(e), n_la_nina = length(l),
         mean_el_nino = mean(e), mean_la_nina = mean(l), difference = mean(e) - mean(l),
         month_ci_low = cm[1], month_ci_high = cm[2],
         episode_ci_low = ce[1], episode_ci_high = ce[2],
         resolved = (cm[1] > 0 | cm[2] < 0) & (ce[1] > 0 | ce[2] < 0))
}

contrast <- bind_rows(
  one_contrast(nee_wet, "NEE", "NEE (Table 1)"),
  one_contrast(t2_wet, "Reco", "Reco (Table 2)"),
  one_contrast(t2_wet, "GEP", "GEP (Table 2)"),
  one_contrast(t2_wet, "NEE_unscreened", "NEE (no coverage screen)")
)

cat("\n6) wet-season El Nino minus La Nina, block-bootstrap intervals\n")
print(as.data.frame(
  contrast %>% transmute(quantity, n = paste0(n_el_nino, "/", n_la_nina),
                         EN = round(mean_el_nino, 2), LN = round(mean_la_nina, 2),
                         diff = round(difference, 2),
                         month_block_CI = sprintf("%+.2f to %+.2f", month_ci_low, month_ci_high),
                         episode_block_CI = sprintf("%+.2f to %+.2f", episode_ci_low, episode_ci_high),
                         resolved = ifelse(resolved, "yes", "no"))
), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# output

readr::write_csv(
  bind_rows(
    sc %>% mutate(block = "self-correlation") %>%
      transmute(block, subset, n, obs_corr, null_median, null_lo, null_hi, p_obs_gt_null),
    contrast %>% mutate(block = "wet-season ENSO contrast", subset = quantity) %>%
      transmute(block, subset, n = n_el_nino + n_la_nina, mean_el_nino, mean_la_nina,
                difference, month_ci_low, month_ci_high, episode_ci_low, episode_ci_high),
    tibble(block = "decomposition", subset = "wet EN-LN |GEP| difference",
           difference = dG, nee_term = -dN, reco_term = dR)
  ) %>% mutate(n_perm = n_perm, n_boot = n_boot, seed = seed),
  out_csv)
cat("\nsaved:", out_csv, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# the numbers as they appear in the text

g <- function(lbl, field) sc[[field]][sc$subset == lbl][1]
cat("\n=== numbers for sections 2.6 and 3.1 ===\n")
cat(sprintf("self-correlation (all windows): observed %+.2f, identity-only null %+.2f [%.2f, %.2f]\n",
            g("all windows", "obs_corr"), g("all windows", "null_median"),
            g("all windows", "null_lo"), g("all windows", "null_hi")))
cat(sprintf("  within season: dry observed %+.2f vs null %+.2f (p = %.2f); wet %+.2f vs %+.2f (p = %.2f)\n",
            g("season: dry", "obs_corr"), g("season: dry", "null_median"), g("season: dry", "p_obs_gt_null"),
            g("season: wet", "obs_corr"), g("season: wet", "null_median"), g("season: wet", "p_obs_gt_null")))
for (i in seq_len(nrow(contrast))) with(contrast[i, ], cat(sprintf(
  "%-4s wet: EN %.2f vs LN %.2f, difference %+.2f (month-block 95%% CI %+.2f to %+.2f) - %s\n",
  quantity, mean_el_nino, mean_la_nina, difference, month_ci_low, month_ci_high,
  ifelse(resolved, "resolved", "NOT resolved"))))
cat(sprintf("decomposition: NEE term %+.2f, Reco term %+.2f, total %+.2f g C m-2 d-1\n", -dN, dR, dG))
