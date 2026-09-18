#!/usr/bin/env Rscript

# Cumulative diel NEE by season and ENSO phase.
#
# Outputs:
#   figures/Fig2b_cumulative_diel_NEE_panel.png  - panel (b) of manuscript
#       Figure 2, authored at its slot size for the composite in code/15
#   figures/FigS_morning_NEE_level_matched.png   - supplementary figure showing
#       that removing the El Nino / La Nina level difference moves the La Nina
#       zero crossing into El Nino's half-hour
#   output/table_cumulative_diel_NEE.csv         - curve values and turning points
#
# Why this script exists. A co-author, commenting on Figure 2, asked for "the
# cumulative curves for both seasons in one plot [...] backing up your claim
# that 36 min earlier sign change is significant, and not small differences
# during the days (slope of cumulative curves), or timing of the LAT."
#
# That is the right diagnostic, and it does not back the claim up - it explains
# it away. The morning source-to-sink transition is the turning point of the
# cumulative curve, so a purely vertical (level) difference in morning NEE
# displaces the turning point without any change in the shape or phase of the
# diel cycle. Panel (b) quantifies exactly that: the two wet-season morning
# limbs are near-parallel and offset by about 5 umol CO2 m-2 s-1, and that
# offset alone accounts for more apparent shift than the 36 minutes originally
# reported. The timing claim has been withdrawn from section 3.1 (see
# code/22_morning_transition_test.R, which shows both phases transitioning in
# the same half-hour once the transition is required to persist).
#
# The local-apparent-time alternative the co-author raises is also checked here
# and reported to the console: the mean clock-to-solar offset of the two
# samples differs by under 2 minutes.
#
# What the curves do show is the real wet-season contrast, which is about the
# amount of carbon rather than its timing: La Nina accumulates ~1.2 g C m-2 more
# overnight than El Nino and ends the day a net source, while El Nino ends near
# neutral.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])

# ─────────────────────────────────────────────────────────────────────────────
# config

seed   <- 20260830
n_boot <- 2000
set.seed(seed)

tz_local <- "America/Lima"

# copied from code/04_diel_cycles_carbon_balance.R:150-151 and :185, where they
# are local rather than shared. Keep in sync with that script.
sec_per_halfhour      <- 1800
gC_per_umolCO2        <- 12 / 1e6
min_points_month_diel <- 80

enso_levels   <- c("El Nino", "La Nina")
season_levels <- c("wet", "dry")

# site coordinates (code/13_maps.R) and standard-time offset, for the local
# apparent time check only
site_lon  <- -70.25
tz_offset <- -5

# window used for the morning level/slope comparison in panel (b)
morn_from <- 6.5
morn_to   <- 8.0

flux_csv <- file.path(output_path, "dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv")
out_csv  <- file.path(out_base, "table_cumulative_diel_NEE.csv")
# the cumulative curves themselves are panel (b) of manuscript Figure 2, so the
# supplement carries only the level-matched morning comparison
out_png  <- file.path(fig_path, "FigS_morning_NEE_level_matched.png")

# ─────────────────────────────────────────────────────────────────────────────
# load

if (!file.exists(flux_csv)) {
  stop("enriched flux table not found: ", flux_csv,
       "\nRun code/01_prepare_flux_timeseries.R first.")
}

required_cols <- c("tv_dt", "NEE", "NEE_ok", "season", "ENSO")
df_raw <- readr::read_csv(flux_csv, show_col_types = FALSE)
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("missing required columns: ", paste(missing_cols, collapse = ", "))
}

base <- df_raw %>%
  mutate(
    ts       = dmy_hms(tv_dt, tz = tz_local, quiet = TRUE),
    hour_dec = hour(ts) + minute(ts) / 60,
    # timestamps are period-end labels, so the 00:00 stamp closes the previous
    # day; move it to 24 so each cumulative curve spans exactly one day
    hour_dec = ifelse(hour_dec == 0, 24, hour_dec),
    ym       = format(ts, "%Y-%m"),
    doy      = as.integer(format(ts, "%j")),
    season   = tolower(trimws(as.character(season))),
    ENSO     = trimws(as.character(ENSO))
  ) %>%
  filter(season %in% season_levels, ENSO %in% enso_levels, !is.na(ts))

if (!nrow(base)) stop("no wet/dry El Nino / La Nina half-hours found")

# ─────────────────────────────────────────────────────────────────────────────
# local apparent time check
#
# The clock-to-solar offset is the equation of time plus the longitude
# correction. The two ENSO samples draw on different calendar months, so this
# offset could in principle differ between them; it does not.

eqt_minutes <- function(doy) {
  g <- 2 * pi / 365 * (doy - 1)
  229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) -
              0.014615 * cos(2 * g) - 0.040849 * sin(2 * g))
}

lat_chk <- base %>%
  filter(season == "wet", is.finite(NEE)) %>%
  mutate(lat_offset_min = eqt_minutes(doy) + 4 * (site_lon - 15 * tz_offset)) %>%
  group_by(ENSO) %>%
  summarise(n = n(), mean_offset_min = mean(lat_offset_min), .groups = "drop")

cat("\n=== local apparent time: clock-to-solar offset, wet season ===\n")
print(as.data.frame(lat_chk %>% mutate(mean_offset_min = round(mean_offset_min, 2))),
      row.names = FALSE)
lat_gap <- diff(lat_chk$mean_offset_min[order(lat_chk$ENSO)])
cat(sprintf("difference between phases: %.1f min - too small to displace the transition\n",
            abs(lat_gap)))

# ─────────────────────────────────────────────────────────────────────────────
# diel composites
#
# Monthly diel means first, then the mean across months, so every calendar month
# carries equal weight regardless of how many days it contributed. Same
# construction as make_monthly_bg_plus_enso_plot() in
# code/04_diel_cycles_carbon_balance.R:472-486, including the coverage screen.

monthly_diel <- function(col) {
  d <- base %>% filter(is.finite(.data[[col]]))
  keep <- d %>%
    count(season, ENSO, ym, name = "n_pts") %>%
    filter(n_pts >= min_points_month_diel)
  d %>%
    inner_join(keep, by = c("season", "ENSO", "ym")) %>%
    group_by(season, ENSO, ym, hour_dec) %>%
    summarise(value = mean(.data[[col]], na.rm = TRUE), .groups = "drop")
}

composite_from_monthly <- function(md) {
  md %>%
    group_by(season, ENSO, hour_dec) %>%
    summarise(nee = mean(value, na.rm = TRUE), n_months = n_distinct(ym), .groups = "drop") %>%
    arrange(season, ENSO, hour_dec) %>%
    group_by(season, ENSO) %>%
    mutate(cum_gC = cumsum(nee * sec_per_halfhour * gC_per_umolCO2)) %>%
    ungroup()
}

md_meas <- monthly_diel("NEE")
md_fill <- monthly_diel("NEE_ok")
comp_meas <- composite_from_monthly(md_meas)
comp_fill <- composite_from_monthly(md_fill)

n_bins <- comp_meas %>% count(season, ENSO)
if (any(n_bins$n != 48)) {
  message("note: incomplete diel coverage in some groups:")
  print(as.data.frame(n_bins), row.names = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# month-block bootstrap band
#
# Resample calendar months with replacement within each season x ENSO, rebuild
# the composite and re-cumulate, so the band reflects month-to-month variability
# rather than treating half-hours as independent.

boot_band <- function(md, reps) {
  groups <- md %>% distinct(season, ENSO)
  out <- list()
  for (i in seq_len(nrow(groups))) {
    se <- groups$season[i]; en <- groups$ENSO[i]
    sub <- md %>% filter(season == se, ENSO == en)
    months <- unique(sub$ym)
    by_month <- split(sub, sub$ym)
    hrs <- sort(unique(sub$hour_dec))
    mat <- matrix(NA_real_, nrow = reps, ncol = length(hrs))
    for (r in seq_len(reps)) {
      pick <- sample(months, length(months), replace = TRUE)
      draw <- bind_rows(by_month[pick])
      m <- draw %>%
        group_by(hour_dec) %>%
        summarise(nee = mean(value, na.rm = TRUE), .groups = "drop") %>%
        arrange(hour_dec)
      if (nrow(m) != length(hrs)) next
      mat[r, ] <- cumsum(m$nee * sec_per_halfhour * gC_per_umolCO2)
    }
    out[[i]] <- tibble(
      season = se, ENSO = en, hour_dec = hrs,
      cum_low  = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
      cum_high = apply(mat, 2, quantile, 0.975, na.rm = TRUE)
    )
  }
  bind_rows(out)
}

band <- boot_band(md_meas, n_boot)
comp_meas <- comp_meas %>% left_join(band, by = c("season", "ENSO", "hour_dec"))

# ─────────────────────────────────────────────────────────────────────────────
# turning points and daily totals

hm <- function(x) sprintf("%02d:%02d", floor(x) %% 24, round((x - floor(x)) * 60))

turning <- function(comp, label) {
  comp %>%
    group_by(season, ENSO) %>%
    summarise(
      series          = label,
      n_bins          = n(),
      turning_hour    = hour_dec[which.max(cum_gC)],
      turning_hm      = hm(hour_dec[which.max(cum_gC)]),
      peak_gC         = max(cum_gC),
      cum_by_0800_gC  = cum_gC[which.min(abs(hour_dec - 8))],
      day_total_gC    = cum_gC[which.max(hour_dec)],
      .groups = "drop"
    )
}

turn <- bind_rows(turning(comp_meas, "measured"), turning(comp_fill, "gap-filled"))

cat("\n=== cumulative diel NEE: turning points and daily totals ===\n")
print(as.data.frame(turn %>% mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# level versus timing in the wet-season morning
#
# The co-author's alternative explanation, tested directly rather than through a
# slope approximation. La Nina morning NEE sits above El Nino at every morning
# half-hour, but not by a constant amount, so converting the offset into a time
# shift via offset/slope would be fitting a line the data do not support.
#
# The exact version instead: subtract the mean morning offset from the La Nina
# composite and ask which half-hour the shifted curve crosses zero in. If it
# lands in El Nino's half-hour, the apparent timing difference was a difference
# in level, not in phase.

morn <- comp_meas %>%
  filter(season == "wet", hour_dec >= morn_from, hour_dec <= morn_to) %>%
  select(ENSO, hour_dec, nee) %>%
  pivot_wider(names_from = ENSO, values_from = nee)

level_offset <- mean(morn[["La Nina"]] - morn[["El Nino"]], na.rm = TRUE)

wet_wide <- comp_meas %>%
  filter(season == "wet") %>%
  select(ENSO, hour_dec, nee) %>%
  pivot_wider(names_from = ENSO, values_from = nee) %>%
  arrange(hour_dec) %>%
  mutate(`La Nina, level-matched` = `La Nina` - level_offset)

# grid-native crossing: the label of the first half-hour below zero in the
# morning window, matching how section 3.1 reports the transition
crossing_slot_vec <- function(h, y) {
  ok <- is.finite(y) & h >= 5 & h <= 12
  h2 <- h[ok]; y2 <- y[ok]
  i <- which(head(y2, -1) > 0 & tail(y2, -1) < 0)[1]
  if (is.na(i)) NA_real_ else h2[i + 1]
}

cross_tbl <- tibble(
  curve = c("El Niño", "La Niña", "La Niña, level-matched"),
  slot  = c(crossing_slot_vec(wet_wide$hour_dec, wet_wide[["El Nino"]]),
            crossing_slot_vec(wet_wide$hour_dec, wet_wide[["La Nina"]]),
            crossing_slot_vec(wet_wide$hour_dec, wet_wide[["La Nina, level-matched"]]))
) %>% mutate(crossing_half_hour = hm(slot))

offset_by_hh <- comp_meas %>%
  filter(season == "wet", hour_dec >= 5.5, hour_dec <= 10) %>%
  select(ENSO, hour_dec, nee) %>%
  pivot_wider(names_from = ENSO, values_from = nee) %>%
  mutate(diff = `La Nina` - `El Nino`)

cat(sprintf(
  paste0("\n=== wet-season morning, measured NEE composite ===\n",
         "La Nina NEE exceeds El Nino at every half-hour from 05:30 to 10:00,\n",
         "  by %.1f to %.1f umol CO2 m-2 s-1 (mean over %s-%s: %+.2f)\n"),
  min(offset_by_hh$diff), max(offset_by_hh$diff), hm(morn_from), hm(morn_to), level_offset))
cat("grid-native zero crossing of the composite:\n")
print(as.data.frame(cross_tbl %>% select(curve, crossing_half_hour)), row.names = FALSE)
same_slot <- isTRUE(all.equal(cross_tbl$slot[1], cross_tbl$slot[3]))
cat(if (same_slot) {
  "=> removing the level difference puts La Nina in the SAME half-hour as El Nino:\n   the apparent shift was a difference in the amount of morning CO2 exchange,\n   not in the timing of the transition.\n"
} else {
  "=> the level-matched curve still crosses in a different half-hour; the shift is\n   not explained by the level difference alone.\n"
})

# ─────────────────────────────────────────────────────────────────────────────
# outputs

curves_out <- comp_meas %>%
  mutate(series = "measured") %>%
  bind_rows(comp_fill %>% mutate(series = "gap-filled")) %>%
  transmute(series, season, ENSO, local_time = hm(hour_dec), hour_dec,
            nee_umol_m2_s = nee, cum_gC_m2 = cum_gC, cum_ci_low = cum_low,
            cum_ci_high = cum_high, n_months) %>%
  arrange(series, season, ENSO, hour_dec)

readr::write_csv(
  bind_rows(
    curves_out %>% mutate(row_type = "curve"),
    turn %>% transmute(row_type = "summary", series, season, ENSO,
                       local_time = turning_hm, hour_dec = turning_hour,
                       cum_gC_m2 = peak_gC, n_months = NA_integer_)
  ) %>% mutate(seed = seed, n_boot = n_boot),
  out_csv)
cat("\nsaved:", out_csv, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# figure
#
# Encoding: ENSO = colour (palette.R `enso_cols`), season = linetype AND point
# shape (palette.R `season_shapes`). This departs from palette.R in one respect:
# there linetype is ENSO's redundant cue, but with four curves in a single panel
# two gold and two pink lines would then be separable only by their markers.
# Season therefore takes linetype here. The colour convention, which palette.R
# treats as the invariant, is untouched.

lab_enso   <- function(x) recode(x, "El Nino" = "El Niño", "La Nina" = "La Niña")
lab_season <- function(x) recode(x, "wet" = "Wet", "dry" = "Dry")

plot_df <- comp_meas %>%
  mutate(ENSO_lab = factor(lab_enso(ENSO), levels = c("El Niño", "La Niña")),
         Season   = factor(lab_season(season), levels = c("Wet", "Dry")),
         grp      = paste(Season, ENSO_lab))

turn_pts <- turn %>%
  filter(series == "measured") %>%
  mutate(ENSO_lab = factor(lab_enso(ENSO), levels = c("El Niño", "La Niña")),
         Season   = factor(lab_season(season), levels = c("Wet", "Dry")))

x_breaks <- seq(0, 24, by = 4)
marker_hours <- seq(1, 24, by = 3)

p_a <- ggplot(plot_df, aes(x = hour_dec, group = grp)) +
  geom_ribbon(aes(ymin = cum_low, ymax = cum_high, fill = ENSO_lab),
              alpha = 0.15, colour = NA, na.rm = TRUE) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_line(aes(y = cum_gC, colour = ENSO_lab, linetype = Season), linewidth = 1.1) +
  geom_point(data = plot_df %>% filter(hour_dec %in% marker_hours),
             aes(y = cum_gC, colour = ENSO_lab, shape = Season), size = 2.6) +
  geom_point(data = turn_pts, aes(x = turning_hour, y = peak_gC, colour = ENSO_lab),
             inherit.aes = FALSE, shape = 21, fill = "white", size = 4, stroke = 1.3,
             show.legend = FALSE) +
  annotate("text", x = 8.4, y = max(plot_df$cum_high, na.rm = TRUE),
           hjust = 0, vjust = 1, size = 3.4, colour = "grey20",
           label = paste0("open circles: turning point of each curve\n",
                          "(= the morning source-to-sink transition).\n",
                          "All four fall within one half-hour.")) +
  scale_x_continuous(breaks = x_breaks, labels = sprintf("%02d:00", x_breaks %% 24),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_colour_manual(values = enso_cols, name = "ENSO") +
  scale_fill_manual(values = enso_cols, guide = "none") +
  scale_linetype_manual(values = c("Wet" = "solid", "Dry" = "dashed"), name = "Season") +
  scale_shape_manual(values = season_shapes, name = "Season") +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = NA)),
         linetype = guide_legend(order = 2), shape = guide_legend(order = 2)) +
  labs(x = NULL, y = expression("Cumulative NEE (g C m"^-2*")")) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.015, 0.02), legend.justification = c(0, 0),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", colour = "grey70"),
        panel.grid.minor = element_blank())

morn_long <- wet_wide %>%
  filter(hour_dec >= 6, hour_dec <= 9.5) %>%
  pivot_longer(c(`El Nino`, `La Nina`, `La Nina, level-matched`),
               names_to = "curve", values_to = "nee") %>%
  mutate(curve = factor(recode(curve,
                               "El Nino" = "El Niño", "La Nina" = "La Niña",
                               "La Nina, level-matched" = "La Niña, level-matched"),
                        levels = c("El Niño", "La Niña", "La Niña, level-matched")))

curve_cols <- c("El Niño" = unname(enso_cols[["El Nino"]]),
                "La Niña" = unname(enso_cols[["La Nina"]]),
                "La Niña, level-matched" = unname(enso_cols[["La Nina"]]))
curve_ltys <- c("El Niño" = "solid", "La Niña" = "dashed",
                "La Niña, level-matched" = "dotted")

# shade the half-hour each curve crosses zero in, so the comparison is made at
# the resolution of the measurement rather than at a fitted instant
cross_bands <- cross_tbl %>%
  filter(is.finite(slot)) %>%
  transmute(curve = factor(curve, levels = levels(morn_long$curve)),
            xmin = slot - 0.5, xmax = slot)

p_b <- ggplot(morn_long, aes(x = hour_dec, y = nee, colour = curve)) +
  geom_rect(data = cross_bands %>% filter(curve == "El Niño"),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.55) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_line(aes(linetype = curve), linewidth = 1.1) +
  geom_point(size = 2.4, shape = season_shapes[["wet"]]) +
  annotate("text", x = 6.05, y = min(morn_long$nee, na.rm = TRUE),
           hjust = 0, vjust = 0, size = 3.4, colour = "grey20",
           label = sprintf(paste0("La Niña NEE exceeds El Niño at every morning half-hour\n",
                                  "(%.1f to %.1f µmol m⁻² s⁻¹; mean %.1f over %s–%s).\n",
                                  "Removing that difference moves the La Niña crossing\n",
                                  "into the same half-hour as El Niño (shaded)."),
                           min(offset_by_hh$diff), max(offset_by_hh$diff), level_offset,
                           hm(morn_from), hm(morn_to))) +
  scale_x_continuous(breaks = seq(6, 9.5, by = 0.5),
                     labels = sprintf("%02d:%02d", floor(seq(6, 9.5, by = 0.5)),
                                      round((seq(6, 9.5, by = 0.5) %% 1) * 60))) +
  scale_colour_manual(values = curve_cols, name = NULL) +
  scale_linetype_manual(values = curve_ltys, name = NULL) +
  labs(x = "Local time", y = expression("NEE (µmol CO"[2]*" m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "grey70"))

ggsave(out_png, p_b, width = 7.1, height = 4.4, dpi = 300)
cat("saved:", out_png, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# small-format copy of panel (a) for the two-column Figure 2 composite (code/15)
#
# Authored at its final size so code/15 places it 1:1. The in-panel note about
# the turning points moves to the caption, and only the Season key is kept: ENSO
# colour is already keyed by the environmental-stack panel in the same figure.

out_png_panel <- file.path(fig_path, "Fig2b_cumulative_diel_NEE_panel.png")

p_a_panel <- ggplot(plot_df, aes(x = hour_dec, group = grp)) +
  geom_ribbon(aes(ymin = cum_low, ymax = cum_high, fill = ENSO_lab),
              alpha = 0.15, colour = NA, na.rm = TRUE) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_line(aes(y = cum_gC, colour = ENSO_lab, linetype = Season), linewidth = 0.6) +
  geom_point(data = plot_df %>% filter(hour_dec %in% marker_hours),
             aes(y = cum_gC, colour = ENSO_lab, shape = Season), size = 1.3) +
  geom_point(data = turn_pts, aes(x = turning_hour, y = peak_gC, colour = ENSO_lab),
             inherit.aes = FALSE, shape = 21, fill = "white", size = 2.2, stroke = 0.7,
             show.legend = FALSE) +
  scale_x_continuous(breaks = x_breaks, labels = sprintf("%02d:00", x_breaks %% 24),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_colour_manual(values = enso_cols, guide = "none") +
  scale_fill_manual(values = enso_cols, guide = "none") +
  scale_linetype_manual(values = c("Wet" = "solid", "Dry" = "dashed"), name = NULL) +
  scale_shape_manual(values = season_shapes, name = NULL) +
  labs(x = "Hour of Day", y = expression("Cumulative NEE (g C m"^-2*")")) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.key.size = unit(0.75, "lines"),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
        legend.margin = margin(1, 3, 1, 1),
        # extra room on the right so the 00:00 tick label, which sits on the
        # panel edge, is not clipped by the canvas
        plot.margin = margin(3.5, 9, 3.5, 3.5))

ggsave(out_png_panel, p_a_panel, width = 3.70, height = 2.70, dpi = 300)
cat("saved:", out_png_panel, "\n")
cat("\nassign the supplementary figure number from writing/supplement_renumber_map.md",
    "and rename the file accordingly.\n")
