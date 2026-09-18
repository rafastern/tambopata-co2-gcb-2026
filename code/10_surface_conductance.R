# ───────────────────────────────────────────────────────────────────────────────
# surface conductance (gs) from lee et al. (2021) + day/night + decile-panel figures
#
# what this script does:
#   1) loads tambopata_48points_per_month.csv
#   2) computes gs (m s^-1) from inverted penman–monteith
#   3) applies filters (day/night, optional pm-transition/unstable removal, high-gs removal)
#   4) saves:
#        a) old-style two-panel gs vs vpd and gs vs ta (cap + full-range)
#        b) NEW: faceted decile-panel figures (same style as your NEE_ok facet script)
#
# requested outputs (the new part):
#   split by day/night (separate figures) and parsed by bins (deciles) of vpd or ta
#   x-axis variables included:
#     - vpd
#     - ta
#     - ts
#     - swc
#
# panels (deciles) included:
#   - panels = vpd deciles : gs vs ta, gs vs ts, gs vs swc
#   - panels = ta  deciles : gs vs vpd, gs vs ts, gs vs swc
#
# day/night versions:
#   for each of the above, the script saves both:
#     - day   (par >= threshold, else hour fallback)
#     - night (par <  threshold, else hour fallback)
#
# each faceted panel includes:
#   - scatter + lm line
#   - y=0 horizontal line
#   - per-panel r^2, p-value (slope), N
#   - strip label includes the decile range (numeric min–max) of the binned variable
#
# outputs:
#   - png + pdf for each figure into paper_figures_path
#   - diagnostic tables + diagnostic scatter plots into output_path / diag_path
# ───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
  library(hms)
})

# ───────────────────────────────────────────────────────────────────────────────
# paths
paths_file <- c("paths.R", file.path("code", "paths.R"))
source(paths_file[file.exists(paths_file)][1])
graphs_root        <- graphs_path
graphs_path        <- file.path(graphs_root, "Gs_decile_conditioned")
diag_path          <- file.path(graphs_root, "Gs_outlier_diagnostics")

dir.create(paper_figures_path, showWarnings = FALSE, recursive = TRUE)
dir.create(graphs_path, showWarnings = FALSE, recursive = TRUE)
dir.create(diag_path, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# config
gs_thresh <- 0.15

le_transition_thresh <- 20
exclude_pm_transition_in_figures <- TRUE
exclude_pm_transition_in_saved_gs <- FALSE
exclude_pm_unstable_in_figures <- FALSE

par_day_threshold <- 20
day_start_hm <- "06:30"
day_end_hm   <- "18:30"

# old-style figure configs
gs_cap <- 0.10
median_bins <- 22
ta_line_min <- 20   # panel (b): draw the Gs-vs-Ta bin-median line only for Tair >= this (°C);
                    # low-Ta bins have too few points to summarize. all points stay plotted.
x_q_hi_vpd <- 0.99
x_q_lo_ta  <- 0.01
x_q_hi_ta  <- 0.99
vpd_color_q_hi <- 0.95
ta_color_q_lo  <- 0.01
ta_color_q_hi  <- 0.99

# facet style configs
point_alpha  <- 0.45
point_size   <- 1.6
lm_linewidth <- 0.8
facet_ncol   <- 5

fig_width_mm  <- 230
fig_height_mm <- 170
fig_dpi <- 600

p_digits <- 3
r2_digits <- 2
range_digits <- 2

# x-limits quantiles per variable (robust limits on the full filtered set)
x_q_lo <- 0.01
x_q_hi <- 0.99

# ───────────────────────────────────────────────────────────────────────────────
# load data
csv_fp <- file.path(input_folder, "tambopata_48points_per_month.csv")
dat <- read.csv(csv_fp, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# fix common mojibake for enso labels
if ("ENSO" %in% names(dat)) {
  dat$ENSO <- gsub("NiÃ±a", "Niña", dat$ENSO, fixed = TRUE)
  dat$ENSO <- gsub("NiÃ±o", "Niño", dat$ENSO, fixed = TRUE)
}

# ───────────────────────────────────────────────────────────────────────────────
# helper functions: saturation vapor pressure and slope
esat_kPa <- function(Tair_C) {
  0.6108 * exp(17.27 * Tair_C / (Tair_C + 237.3))
}

slope_svp <- function(Tair_C) {
  es <- esat_kPa(Tair_C)
  4098 * es / (Tair_C + 237.3)^2
}

# aerodynamic conductance (m s^-1)
Gaero_from_u <- function(u_star, u) {
  u_safe <- ifelse(is.finite(u) & u > 0, u, NA_real_)
  ge <- (u_star^2) / u_safe
  re <- 1 / ge
  rb <- 6.2 * u_star^(-2/3)
  inv_Gaero <- re + rb
  ga <- 1 / inv_Gaero
  ga[!is.finite(ga)] <- NA_real_
  ga
}

# surface conductance (m s^-1), inverted penman–monteith
Gs_from_fluxes <- function(LE, H, VPD, Tair_C,
                           u_star, u,
                           rho_air = 1.2,
                           cp_air  = 1007,
                           gamma   = 0.066) {
  
  beta  <- H / LE
  s     <- slope_svp(Tair_C)
  gaero <- Gaero_from_u(u_star, u)
  
  term1 <- (rho_air * cp_air * VPD) / (gamma * LE)
  term2 <- (s * beta - gamma) / (gamma * gaero)
  
  inv_gs <- term1 + term2
  gs <- 1 / inv_gs
  gs[!is.finite(gs)] <- NA_real_
  gs
}

# ───────────────────────────────────────────────────────────────────────────────
# column mapping helpers
find_col <- function(candidates, names_vec) {
  hits <- candidates[candidates %in% names_vec]
  if (length(hits) == 0) NA_character_ else hits[1]
}

nms <- names(dat)

le_candidates    <- c("LE", "LE_f", "LE_fqcOK", "LE_qcOK", "LE_CORR", "LE_orig")
h_candidates     <- c("H", "H_f", "H_fqcOK", "H_qcOK", "H_CORR", "H_orig")
vpd_candidates   <- c("VPD_kPa", "VPD_kpa", "VPD_kP", "VPD_kpA", "VPD")
tair_candidates  <- c("Tair", "TA_1_1_1", "TA", "Ta")
ustar_candidates <- c("USTAR", "u_star", "Ustar")
wind_candidates  <- c("WIND", "WS_1_1_1", "WS", "U", "WS_ms")

if (!"PAR_corrected_SWin" %in% names(dat) && "PPFD_IN_1_1_1" %in% names(dat)) {
  dat$PAR_corrected_SWin <- suppressWarnings(as.numeric(dat$PPFD_IN_1_1_1))
}
if (!"PAR" %in% names(dat) && "PAR_corrected_SWin" %in% names(dat)) {
  dat$PAR <- suppressWarnings(as.numeric(dat$PAR_corrected_SWin))
}
if (!"PAR" %in% names(dat)) {
  stop("MATLAB-corrected PAR column is missing from surface-conductance input")
}
par_mismatch <- is.finite(dat$PAR) & is.finite(dat$PAR_corrected_SWin) &
  abs(dat$PAR - dat$PAR_corrected_SWin) > 1e-8
if (any(par_mismatch, na.rm = TRUE)) {
  stop("canonical PAR does not match PAR_corrected_SWin in surface-conductance input")
}

par_candidates   <- c("PAR", "PAR_corrected_SWin")
hour_candidates  <- c("hour", "Hour", "HOUR", "time_of_day", "tod")

ts_candidates    <- c("TS_3", "TS", "Ts", "Tsoil", "SoilT")
swc_candidates   <- c("SWC_1_1_1", "SWC", "SWC_1", "SWC_vol", "SWC_vf")

col_LE    <- find_col(le_candidates,    nms)
col_H     <- find_col(h_candidates,     nms)
col_VPD   <- find_col(vpd_candidates,   nms)
col_Tair  <- find_col(tair_candidates,  nms)
col_USTAR <- find_col(ustar_candidates, nms)
col_WIND  <- find_col(wind_candidates,  nms)
col_PAR   <- find_col(par_candidates,   nms)
col_HOUR  <- find_col(hour_candidates,  nms)

col_TS    <- find_col(ts_candidates,    nms)
col_SWC   <- find_col(swc_candidates,   nms)

message("mapped columns:")
print(c(
  LE    = col_LE,
  H     = col_H,
  VPD   = col_VPD,
  Tair  = col_Tair,
  USTAR = col_USTAR,
  WIND  = col_WIND,
  PAR   = col_PAR,
  hour  = col_HOUR,
  TS    = col_TS,
  SWC   = col_SWC
))

required <- c(col_LE, col_H, col_VPD, col_Tair, col_USTAR, col_WIND)
if (any(is.na(required))) stop("missing at least one required column for gs computation. check csv header.")

# ───────────────────────────────────────────────────────────────────────────────
# rename to standard names and compute gs
rename_map <- c(
  LE    = col_LE,
  H     = col_H,
  VPD   = col_VPD,
  Tair  = col_Tair,
  USTAR = col_USTAR,
  WIND  = col_WIND
)
if (!is.na(col_PAR))  rename_map <- c(rename_map, PAR  = col_PAR)
if (!is.na(col_HOUR)) rename_map <- c(rename_map, hour = col_HOUR)
if (!is.na(col_TS))   rename_map <- c(rename_map, TS   = col_TS)
if (!is.na(col_SWC))  rename_map <- c(rename_map, SWC  = col_SWC)

dat <- dat %>%
  dplyr::rename(!!!rename_map) %>%
  mutate(
    LE    = ifelse(is.finite(LE) & LE > 0, LE, NA_real_),
    H     = ifelse(is.finite(H), H, NA_real_),
    VPD   = ifelse(is.finite(VPD) & VPD >= 0, VPD, NA_real_),
    Tair  = ifelse(is.finite(Tair), Tair, NA_real_),
    USTAR = ifelse(is.finite(USTAR) & USTAR > 0, USTAR, NA_real_),
    WIND  = ifelse(is.finite(WIND) & WIND > 0, WIND, NA_real_),
    PAR   = if ("PAR" %in% names(.)) ifelse(is.finite(PAR) & PAR >= 0, PAR, NA_real_) else NA_real_,
    TS    = if ("TS"  %in% names(.)) ifelse(is.finite(TS), TS, NA_real_) else NA_real_,
    # swc_1_1_1 is manaus-calibrated integrated soil-water storage over 0-100 cm, reported in cm as-is
    SWC   = if ("SWC" %in% names(.)) ifelse(is.finite(SWC), SWC, NA_real_) else NA_real_,
    Gs_mps = Gs_from_fluxes(
      LE     = LE,
      H      = H,
      VPD    = VPD,
      Tair_C = Tair,
      u_star = USTAR,
      u      = WIND
    )
  )

# ───────────────────────────────────────────────────────────────────────────────
# day/night flag (par if available, else hour fallback)
make_hour_hms <- function(x) {
  if (inherits(x, "hms")) return(x)
  if (inherits(x, "difftime")) return(as_hms(x))
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as_hms(format(x, "%H:%M:%S")))
  out <- suppressWarnings(hms::as_hms(x))
  if (all(is.na(out))) {
    out <- suppressWarnings(as_hms(parse_date_time(x, orders = c("HMS", "HM"))))
  }
  out
}

day_start <- hms::as_hms(paste0(day_start_hm, ":00"))
day_end   <- hms::as_hms(paste0(day_end_hm,   ":00"))

dat <- dat %>%
  mutate(
    hour_hms = if ("hour" %in% names(.)) make_hour_hms(hour) else as_hms(NA_character_),
    day_by_hour = is.finite(as.numeric(hour_hms)) & (hour_hms >= day_start) & (hour_hms < day_end),
    day_by_par  = if ("PAR" %in% names(.)) is.finite(PAR) & (PAR >= par_day_threshold) else NA,
    is_day = if ("PAR" %in% names(.)) ifelse(is.na(day_by_par), day_by_hour, day_by_par) else day_by_hour
  )

# optionally remove pm-transition rows from saved gs too
if (exclude_pm_transition_in_saved_gs) {
  dat <- dat %>%
    mutate(flag_pm_transition = is.finite(LE) & (LE < le_transition_thresh) & is.finite(H) & (H < 0)) %>%
    filter(!flag_pm_transition) %>%
    select(-flag_pm_transition)
  message(sprintf("note: excluded pm-transition rows from saved gs: le < %.1f & h < 0", le_transition_thresh))
}

# save updated data with gs
# keep only finite, non-negative gs in the saved product, matching methods 2.10
# ("retained only finite, non-negative Gs"); raw values remain in `dat` for the
# diagnostics below so the term1/term2 cancellation analysis is unaffected
out_csv <- file.path(output_path, "tambopata_48points_per_month_with_gs.csv")
dat_save <- dat %>%
  mutate(Gs_mps = ifelse(is.finite(Gs_mps) & Gs_mps >= 0, Gs_mps, NA_real_))
write.csv(dat_save, out_csv, row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# diagnostics terms for investigating large gs
dat_diag <- dat %>%
  mutate(
    beta  = H / LE,
    s     = slope_svp(Tair),
    Gaero = Gaero_from_u(USTAR, WIND),
    term1 = (1.2 * 1007 * VPD) / (0.066 * LE),
    term2 = (s * beta - 0.066) / (0.066 * Gaero),
    inv_Gs = term1 + term2,
    abs_inv_Gs = abs(inv_Gs),
    cancellation = abs(term1 + term2),
    flag_highGs = is.finite(Gs_mps) & (Gs_mps > gs_thresh),
    flag_pm_transition = is.finite(LE) & (LE < le_transition_thresh) & is.finite(H) & (H < 0),
    flag_pm_unstable = is.finite(beta) & (beta < 0) & is.finite(H) & (H < 0)
  ) %>%
  filter(is.finite(Gs_mps))

dat_highGs <- dat_diag %>% filter(flag_highGs)

message(sprintf("high-gs threshold = %.3f m s^-1", gs_thresh))
message(sprintf("n(high-gs rows) = %d", nrow(dat_highGs)))

high_csv <- file.path(
  output_path,
  paste0("tambopata_high_gs_rows_over_", gsub("\\.", "p", sprintf("%.2f", gs_thresh)), ".csv")
)
write.csv(dat_highGs, high_csv, row.names = FALSE)

vars_to_check <- c(
  "Gs_mps",
  "LE", "H", "beta", "VPD", "Tair", "PAR",
  "USTAR", "WIND", "Gaero", "s",
  "term1", "term2", "inv_Gs", "abs_inv_Gs",
  "flag_pm_unstable", "flag_pm_transition", "is_day"
)

summary_compare <- lapply(vars_to_check, function(v) {
  if (!v %in% names(dat_diag)) return(NULL)
  data.frame(
    variable = v,
    n_all = sum(!is.na(dat_diag[[v]])),
    n_high = sum(!is.na(dat_highGs[[v]])),
    median_all  = suppressWarnings(median(dat_diag[[v]], na.rm = TRUE)),
    median_high = suppressWarnings(median(dat_highGs[[v]], na.rm = TRUE)),
    q05_all     = suppressWarnings(unname(quantile(dat_diag[[v]], 0.05, na.rm = TRUE))),
    q95_all     = suppressWarnings(unname(quantile(dat_diag[[v]], 0.95, na.rm = TRUE))),
    q05_high    = suppressWarnings(unname(quantile(dat_highGs[[v]], 0.05, na.rm = TRUE))),
    q95_high    = suppressWarnings(unname(quantile(dat_highGs[[v]], 0.95, na.rm = TRUE)))
  )
}) %>% bind_rows()

sum_csv <- file.path(
  output_path,
  paste0("tambopata_high_gs_summary_compare_over_", gsub("\\.", "p", sprintf("%.2f", gs_thresh)), ".csv")
)
write.csv(summary_compare, sum_csv, row.names = FALSE)

rank_csv <- file.path(
  output_path,
  paste0("tambopata_gs_ranked_by_abs_invGs_over_", gsub("\\.", "p", sprintf("%.2f", gs_thresh)), ".csv")
)

dat_ranked <- dat_diag %>%
  filter(is.finite(inv_Gs)) %>%
  arrange(abs_inv_Gs) %>%
  mutate(row_id = row_number()) %>%
  select(
    row_id,
    Gs_mps,
    flag_highGs,
    flag_pm_unstable,
    flag_pm_transition,
    is_day,
    hour,
    PAR,
    inv_Gs,
    abs_inv_Gs,
    cancellation,
    term1,
    term2,
    LE, H, beta,
    VPD, Tair,
    USTAR, WIND, Gaero, s,
    everything()
  )

write.csv(dat_ranked, rank_csv, row.names = FALSE)

message("wrote diagnostics:")
message(paste0("  ", out_csv))
message(paste0("  ", high_csv))
message(paste0("  ", sum_csv))
message(paste0("  ", rank_csv))

# diagnostic plots
p_diag1 <- ggplot(dat_diag, aes(x = beta, y = Gs_mps, color = flag_highGs)) +
  geom_point(alpha = 0.45) +
  scale_color_manual(values = c("grey60", "red")) +
  labs(x = "beta = H/LE", y = "Gs (m s^-1)", color = paste0("Gs > ", gs_thresh)) +
  theme_minimal()

p_diag2 <- ggplot(dat_diag, aes(x = H, y = LE, color = flag_highGs)) +
  geom_point(alpha = 0.45) +
  scale_color_manual(values = c("grey60", "red")) +
  labs(x = "H", y = "LE", color = paste0("Gs > ", gs_thresh)) +
  theme_minimal()

p_diag3 <- ggplot(dat_diag, aes(x = term1, y = term2, color = flag_highGs)) +
  geom_point(alpha = 0.45) +
  scale_color_manual(values = c("grey60", "red")) +
  labs(x = "term1", y = "term2", color = paste0("Gs > ", gs_thresh)) +
  theme_minimal()

p_diag4 <- ggplot(dat_diag, aes(x = inv_Gs, y = Gs_mps, color = flag_highGs)) +
  geom_point(alpha = 0.45) +
  scale_color_manual(values = c("grey60", "red")) +
  labs(x = "inv_Gs = term1 + term2", y = "Gs (m s^-1)", color = paste0("Gs > ", gs_thresh)) +
  theme_minimal()

p_diag5 <- ggplot(dat_diag, aes(x = abs_inv_Gs, y = Gs_mps, color = flag_highGs)) +
  geom_point(alpha = 0.45) +
  scale_color_manual(values = c("grey60", "red")) +
  labs(x = "|inv_Gs|", y = "Gs (m s^-1)", color = paste0("Gs > ", gs_thresh)) +
  theme_minimal()

ggsave(file.path(diag_path, "diag_beta_vs_gs.png"),     p_diag1, width = 6.5, height = 4.2, dpi = 300)
ggsave(file.path(diag_path, "diag_H_vs_LE.png"),        p_diag2, width = 6.5, height = 4.2, dpi = 300)
ggsave(file.path(diag_path, "diag_term1_vs_term2.png"), p_diag3, width = 6.5, height = 4.2, dpi = 300)
ggsave(file.path(diag_path, "diag_invGs_vs_gs.png"),    p_diag4, width = 6.5, height = 4.2, dpi = 300)
ggsave(file.path(diag_path, "diag_absInvGs_vs_gs.png"), p_diag5, width = 6.5, height = 4.2, dpi = 300)

# ───────────────────────────────────────────────────────────────────────────────
# filtered data for figures
dat_gs <- dat_diag %>%
  filter(
    is.finite(Gs_mps),
    Gs_mps >= 0,
    is.finite(VPD),
    is.finite(Tair)
  )

if (exclude_pm_transition_in_figures) {
  dat_gs <- dat_gs %>% filter(!flag_pm_transition)
  message(sprintf("note: excluded pm-transition rows from figures: le < %.1f & h < 0", le_transition_thresh))
}

if (exclude_pm_unstable_in_figures) {
  dat_gs <- dat_gs %>% filter(!flag_pm_unstable)
  message("note: excluded pm-unstable rows from figures: beta < 0 & h < 0")
}

dat_gs <- dat_gs %>% filter(!flag_highGs)

# split day/night
dat_gs_day   <- dat_gs %>% filter(is_day)
dat_gs_night <- dat_gs %>% filter(!is_day)

message(sprintf("n day rows:   %d", nrow(dat_gs_day)))
message(sprintf("n night rows: %d", nrow(dat_gs_night)))

# robust plotting ranges for old-style
vpd_max <- unname(quantile(dat_gs$VPD,  x_q_hi_vpd, na.rm = TRUE))
ta_max  <- unname(quantile(dat_gs$Tair, x_q_hi_ta,  na.rm = TRUE))
ta_min  <- unname(quantile(dat_gs$Tair, x_q_lo_ta,  na.rm = TRUE))



# ───────────────────────────────────────────────────────────────────────────────
# old-style two-panel figure (kept)
base_theme_old <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA),
    legend.key = element_blank()
  )

tag_theme <- theme(
  plot.tag = element_text(size = 10, hjust = 0, vjust = 1),
  plot.tag.position = c(0.08, 0.98)
)

y_lab <- expression(italic(G)[plain(s)]~"(m s"^{-1}*")")
y_grob <- grid::textGrob(label = y_lab, rot = 90, gp = grid::gpar(fontsize = 12))

ta_col_lo <- unname(quantile(dat_gs$Tair, ta_color_q_lo, na.rm = TRUE))
ta_col_hi <- unname(quantile(dat_gs$Tair, ta_color_q_hi, na.rm = TRUE))

p1_base <- ggplot(dat_gs, aes(x = VPD, y = Gs_mps, color = Tair)) +
  geom_point(alpha = 0.35, size = 1.3) +
  stat_summary_bin(fun = median, bins = median_bins, geom = "line", linewidth = 0.9, color = "black") +
  scale_color_gradient(low = "grey85", high = "grey20", limits = c(ta_col_lo, ta_col_hi), oob = scales::squish) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.002)) +
  labs(x = "VPD (kPa)", y = NULL, color = expression(italic(T)[plain(a)] ~ "(°C)"), tag = "(a)") +
  base_theme_old +
  tag_theme

vpd_col_hi <- unname(quantile(dat_gs$VPD, vpd_color_q_hi, na.rm = TRUE))

p2_base <- ggplot(dat_gs, aes(x = Tair, y = Gs_mps, color = VPD)) +
  geom_point(alpha = 0.35, size = 1.3) +
  stat_summary_bin(fun = median, bins = median_bins, geom = "line", linewidth = 0.9, color = "black") +
  scale_color_gradient(low = "grey85", high = "grey20", limits = c(0, vpd_col_hi), oob = scales::squish) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.002)) +
  labs(x = expression(italic(T)[plain(a)] ~ "(°C)"), y = NULL, color = "VPD (kPa)", tag = "(b)") +
  base_theme_old +
  tag_theme

p1_cap <- p1_base + coord_cartesian(xlim = c(0, vpd_max), ylim = c(0, gs_cap))
p2_cap <- p2_base + coord_cartesian(xlim = c(ta_min, ta_max), ylim = c(0, gs_cap))
fig_gs_cap <- (patchwork::wrap_elements(full = y_grob) | (p1_cap | p2_cap)) +
  patchwork::plot_layout(widths = c(0.055, 1))

p1_full <- p1_base + coord_cartesian(xlim = c(0, vpd_max))
p2_full <- p2_base + coord_cartesian(xlim = c(ta_min, ta_max))
fig_gs_full <- (patchwork::wrap_elements(full = y_grob) | (p1_full | p2_full)) +
  patchwork::plot_layout(widths = c(0.055, 1))

# Superseded capped variant, kept as a diagnostic. FigX_ is this script's prefix
# for superseded capped output (12 others below). It was named Fig7_ under a dead
# figure numbering; the main-text Gs figure is now Figure 7 and is the UNCAPPED
# FigS_Gs_two_panel_DAY_PHYSICAL_FILTER.png written further down -- its caption
# states "no artificial upper cap was imposed", so this file must not read as
# Figure 7.
ggsave(file.path(paper_figures_path, "FigX_Gs_two_panel_CAP_0p10.png"),
       fig_gs_cap, width = 7.2, height = 3.4, dpi = 600, bg = "white")

ggsave(file.path(paper_figures_path, "FigS_Gs_two_panel_full_range.png"),
       fig_gs_full, width = 7.2, height = 3.4, dpi = 600, bg = "white")

# ───────────────────────────────────────────────────────────────────────────────
# decile-panel figures (day/night + ts/swc added)

make_pct_labels <- function(n_bins) {
  pct_lo <- seq(0, 90, by = 10)
  pct_hi <- seq(10, 100, by = 10)
  labs <- paste0(pct_lo, "–", pct_hi)
  if (n_bins != 10) labs <- labs[seq_len(n_bins)]
  labs
}

make_deciles <- function(x) {
  probs <- seq(0, 1, by = 0.1)
  br <- unname(quantile(x, probs = probs, na.rm = TRUE, type = 7))
  br <- sort(unique(br))
  if (length(br) < 3) stop("too few unique values to form deciles.")
  list(breaks = br, n_bins = length(br) - 1)
}

assign_decile_col <- function(df, value_col, breaks, labels, out_col) {
  df %>%
    mutate(
      "{out_col}" := cut(
        .data[[value_col]],
        breaks = breaks,
        include.lowest = TRUE,
        right = TRUE,
        labels = labels
      )
    ) %>%
    filter(!is.na(.data[[out_col]])) %>%
    mutate("{out_col}" := factor(.data[[out_col]], levels = labels, ordered = TRUE))
}

make_decile_range_label_map <- function(df, decile_col, value_col, unit = "", digits = 2) {
  df %>%
    group_by(.data[[decile_col]]) %>%
    summarise(
      lo = min(.data[[value_col]], na.rm = TRUE),
      hi = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dec = as.character(.data[[decile_col]]),
      label = paste0(
        dec,
        "\n(",
        formatC(lo, format = "f", digits = digits),
        "–",
        formatC(hi, format = "f", digits = digits),
        if (unit != "") paste0(" ", unit) else "",
        ")"
      )
    ) %>%
    { setNames(.$label, .$dec) }
}

lm_stats_by_panel <- function(df, facet_col, x_col, y_col = "Gs_mps") {
  df %>%
    group_by(.data[[facet_col]]) %>%
    summarise(
      n = sum(is.finite(.data[[x_col]]) & is.finite(.data[[y_col]])),
      R2 = {
        ok <- is.finite(.data[[x_col]]) & is.finite(.data[[y_col]])
        if (sum(ok) >= 3) summary(lm(.data[[y_col]][ok] ~ .data[[x_col]][ok]))$r.squared else NA_real_
      },
      p = {
        ok <- is.finite(.data[[x_col]]) & is.finite(.data[[y_col]])
        if (sum(ok) >= 3) summary(lm(.data[[y_col]][ok] ~ .data[[x_col]][ok]))$coefficients[2, 4] else NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_txt = ifelse(!is.finite(p), "p=NA",
                     ifelse(p < 0.001, "p<0.001",
                            paste0("p=", formatC(p, format = "f", digits = p_digits)))),
      r2_txt = ifelse(!is.finite(R2), "R²=NA",
                      paste0("R²=", formatC(R2, format = "f", digits = r2_digits))),
      label = paste0(r2_txt, "\n", p_txt, "\nN=", n)
    )
}

base_theme_facets <- theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 13, lineheight = 0.95),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

plot_gs_facets <- function(
    df,
    facet_col,
    x_col,
    x_lab,
    title_txt,
    xlim_vec = NULL,
    ylim_vec = NULL,
    decile_label_map = NULL
) {
  
  stats <- lm_stats_by_panel(df, facet_col = facet_col, x_col = x_col, y_col = "Gs_mps")
  
  facet_layer <- if (is.null(decile_label_map)) {
    facet_wrap(stats::as.formula(paste0("~", facet_col)), ncol = facet_ncol)
  } else {
    facet_wrap(
      stats::as.formula(paste0("~", facet_col)),
      ncol = facet_ncol,
      labeller = as_labeller(decile_label_map)
    )
  }
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = Gs_mps)) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    geom_point(alpha = point_alpha, size = point_size, color = "grey40") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = lm_linewidth) +
    facet_layer +
    geom_text(
      data = stats,
      aes(x = -Inf, y = Inf, label = label, group = .data[[facet_col]]),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.10,
      size = 4.2
    ) +
    labs(
      title = title_txt,
      x = x_lab,
      y = expression(italic(G)[plain(s)]~"(m s"^{-1}*")")
    ) +
    base_theme_facets
  
  if ((is.null(xlim_vec) || !all(is.finite(xlim_vec))) && (is.null(ylim_vec) || !all(is.finite(ylim_vec)))) {
    return(p)
  }
  
  p + coord_cartesian(
    xlim = if (!is.null(xlim_vec) && all(is.finite(xlim_vec))) xlim_vec else NULL,
    ylim = if (!is.null(ylim_vec) && all(is.finite(ylim_vec))) ylim_vec else NULL
  )
}

save_plot_mm <- function(plot_obj, stem) {
  ggsave(
    filename = file.path(paper_figures_path, paste0(stem, ".png")),
    plot     = plot_obj,
    width    = fig_width_mm,
    height   = fig_height_mm,
    units    = "mm",
    dpi      = fig_dpi,
    bg       = "white"
  )
}

# ensure ts and swc exist if you want those plots
if (!("TS" %in% names(dat_gs)))  dat_gs$TS  <- NA_real_
if (!("SWC" %in% names(dat_gs))) dat_gs$SWC <- NA_real_
dat_gs_day$TS  <- if ("TS" %in% names(dat_gs_day)) dat_gs_day$TS else NA_real_
dat_gs_day$SWC <- if ("SWC" %in% names(dat_gs_day)) dat_gs_day$SWC else NA_real_
dat_gs_night$TS  <- if ("TS" %in% names(dat_gs_night)) dat_gs_night$TS else NA_real_
dat_gs_night$SWC <- if ("SWC" %in% names(dat_gs_night)) dat_gs_night$SWC else NA_real_

# build deciles on the full filtered dataset (consistent across day/night)
vpd_dec <- make_deciles(dat_gs$VPD)
vpd_labels <- make_pct_labels(vpd_dec$n_bins)

ta_dec <- make_deciles(dat_gs$Tair)
ta_labels <- make_pct_labels(ta_dec$n_bins)

# assign deciles to each split
dat_day_vpddec   <- assign_decile_col(dat_gs_day,   "VPD",  vpd_dec$breaks, vpd_labels, "VPD_decile")
dat_night_vpddec <- assign_decile_col(dat_gs_night, "VPD",  vpd_dec$breaks, vpd_labels, "VPD_decile")

dat_day_tadec    <- assign_decile_col(dat_gs_day,   "Tair", ta_dec$breaks,  ta_labels,  "TA_decile")
dat_night_tadec  <- assign_decile_col(dat_gs_night, "Tair", ta_dec$breaks,  ta_labels,  "TA_decile")

# label maps (ranges computed on the full dataset’s decile assignments)
dat_all_vpddec <- assign_decile_col(dat_gs, "VPD",  vpd_dec$breaks, vpd_labels, "VPD_decile")
dat_all_tadec  <- assign_decile_col(dat_gs, "Tair", ta_dec$breaks,  ta_labels,  "TA_decile")

vpd_decile_label_map <- make_decile_range_label_map(
  df = dat_all_vpddec, decile_col = "VPD_decile", value_col = "VPD", unit = "kPa", digits = range_digits
)
ta_decile_label_map <- make_decile_range_label_map(
  df = dat_all_tadec,  decile_col = "TA_decile",  value_col = "Tair", unit = "°C",  digits = range_digits
)

# x-lims computed on the full filtered dataset
xlim_ta <- c(unname(quantile(dat_gs$Tair, x_q_lo, na.rm = TRUE)), unname(quantile(dat_gs$Tair, x_q_hi, na.rm = TRUE)))
xlim_vpd <- c(unname(quantile(dat_gs$VPD, x_q_lo, na.rm = TRUE)), unname(quantile(dat_gs$VPD, x_q_hi, na.rm = TRUE)))

xlim_ts <- c(NA_real_, NA_real_)
if (any(is.finite(dat_gs$TS))) {
  xlim_ts <- c(unname(quantile(dat_gs$TS, x_q_lo, na.rm = TRUE)), unname(quantile(dat_gs$TS, x_q_hi, na.rm = TRUE)))
}

xlim_swc <- c(NA_real_, NA_real_)
if (any(is.finite(dat_gs$SWC))) {
  xlim_swc <- c(unname(quantile(dat_gs$SWC, x_q_lo, na.rm = TRUE)), unname(quantile(dat_gs$SWC, x_q_hi, na.rm = TRUE)))
}

# y-lims (full vs capped)
gs_ylim_full <- c(0, max(dat_gs$Gs_mps, na.rm = TRUE))
gs_ylim_cap  <- c(0, min(gs_ylim_full[2], gs_cap))

# Supplement reduction (2026-09-08). The twelve save_full_and_cap() calls below
# emit 24 decile-panel PNGs at 600 dpi; five of them were supplement figures
# S60-S64, and their only job in the manuscript is to show that soil temperature
# exerts no independent control on Gs. One figure carries that claim, so the
# supplement keeps just Gs against VPD conditioned on Ta deciles.
#
# The gate lives here rather than at the twelve call sites so it is one edit and
# so it also catches the CAP twins. The plot objects are still built (cheap
# relative to rendering 5433 x 4016 px), so nothing downstream changes and the
# others can be restored by emptying the whitelist.
gs_supplement_keep <- c("FigS_day_Gs_vs_VPD_split_by_TAdeciles_full_range")

save_plot_mm_gated <- function(plot_obj, stem) {
  if (length(gs_supplement_keep) && !(stem %in% gs_supplement_keep)) {
    message("  skipped (supplement reduction): ", stem)
    return(invisible(NULL))
  }
  save_plot_mm(plot_obj, stem)
}

# helper: build and save a matched full+cap pair
save_full_and_cap <- function(p_full, p_cap, stem_full, stem_cap) {
  save_plot_mm_gated(p_full, stem_full)
  save_plot_mm_gated(p_cap,  stem_cap)
}

# ───────────────────────────────────────────────────────────────────────────────
# figures: panels = vpd deciles (x = ta/ts/swc), day & night
# day: gs vs ta, panels = vpd deciles
p_day_gs_ta_by_vpd_full <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "Tair",
  x_lab = expression(italic(T)[plain(a)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ta (panels = VPD deciles)",
  xlim_vec = xlim_ta, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_day_gs_ta_by_vpd_cap <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "Tair",
  x_lab = expression(italic(T)[plain(a)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ta (panels = VPD deciles) — capped",
  xlim_vec = xlim_ta, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

# night: gs vs ta, panels = vpd deciles
p_night_gs_ta_by_vpd_full <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "Tair",
  x_lab = expression(italic(T)[plain(a)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ta (panels = VPD deciles)",
  xlim_vec = xlim_ta, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_night_gs_ta_by_vpd_cap <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "Tair",
  x_lab = expression(italic(T)[plain(a)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ta (panels = VPD deciles) — capped",
  xlim_vec = xlim_ta, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

# day/night: gs vs ts, panels = vpd deciles
p_day_gs_ts_by_vpd_full <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ts (panels = VPD deciles)",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_day_gs_ts_by_vpd_cap <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ts (panels = VPD deciles) — capped",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

p_night_gs_ts_by_vpd_full <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ts (panels = VPD deciles)",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_night_gs_ts_by_vpd_cap <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ts (panels = VPD deciles) — capped",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

# day/night: gs vs swc, panels = vpd deciles
p_day_gs_swc_by_vpd_full <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "day: Gs vs SWC (panels = VPD deciles)",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_day_gs_swc_by_vpd_cap <- plot_gs_facets(
  df = dat_day_vpddec, facet_col = "VPD_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "day: Gs vs SWC (panels = VPD deciles) — capped",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

p_night_gs_swc_by_vpd_full <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "night: Gs vs SWC (panels = VPD deciles)",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_full,
  decile_label_map = vpd_decile_label_map
)
p_night_gs_swc_by_vpd_cap <- plot_gs_facets(
  df = dat_night_vpddec, facet_col = "VPD_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "night: Gs vs SWC (panels = VPD deciles) — capped",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_cap,
  decile_label_map = vpd_decile_label_map
)

# ───────────────────────────────────────────────────────────────────────────────
# figures: panels = ta deciles (x = vpd/ts/swc), day & night

# day/night: gs vs vpd, panels = ta deciles
p_day_gs_vpd_by_ta_full <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "VPD",
  x_lab = "VPD (kPa)",
  title_txt = "day: Gs vs VPD (panels = Ta deciles)",
  xlim_vec = xlim_vpd, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_day_gs_vpd_by_ta_cap <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "VPD",
  x_lab = "VPD (kPa)",
  title_txt = "day: Gs vs VPD (panels = Ta deciles) — capped",
  xlim_vec = xlim_vpd, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

p_night_gs_vpd_by_ta_full <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "VPD",
  x_lab = "VPD (kPa)",
  title_txt = "night: Gs vs VPD (panels = Ta deciles)",
  xlim_vec = xlim_vpd, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_night_gs_vpd_by_ta_cap <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "VPD",
  x_lab = "VPD (kPa)",
  title_txt = "night: Gs vs VPD (panels = Ta deciles) — capped",
  xlim_vec = xlim_vpd, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

# day/night: gs vs ts, panels = ta deciles
p_day_gs_ts_by_ta_full <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ts (panels = Ta deciles)",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_day_gs_ts_by_ta_cap <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "day: Gs vs Ts (panels = Ta deciles) — capped",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

p_night_gs_ts_by_ta_full <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ts (panels = Ta deciles)",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_night_gs_ts_by_ta_cap <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "TS",
  x_lab = expression(italic(T)[plain(s)]~"("*degree*C*")"),
  title_txt = "night: Gs vs Ts (panels = Ta deciles) — capped",
  xlim_vec = xlim_ts, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

# day/night: gs vs swc, panels = ta deciles
p_day_gs_swc_by_ta_full <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "day: Gs vs SWC (panels = Ta deciles)",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_day_gs_swc_by_ta_cap <- plot_gs_facets(
  df = dat_day_tadec, facet_col = "TA_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "day: Gs vs SWC (panels = Ta deciles) — capped",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

p_night_gs_swc_by_ta_full <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "night: Gs vs SWC (panels = Ta deciles)",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_full,
  decile_label_map = ta_decile_label_map
)
p_night_gs_swc_by_ta_cap <- plot_gs_facets(
  df = dat_night_tadec, facet_col = "TA_decile", x_col = "SWC",
  x_lab = expression(SWC~"(cm)"),
  title_txt = "night: Gs vs SWC (panels = Ta deciles) — capped",
  xlim_vec = xlim_swc, ylim_vec = gs_ylim_cap,
  decile_label_map = ta_decile_label_map
)

# ───────────────────────────────────────────────────────────────────────────────
# save all requested figures

# (A `gs_requested_day_output_manifest` vector used to sit here. It was defined
# once and referenced nowhere, and the supplement reduction of 2026-09-08
# replaced its intent with the live `gs_supplement_keep` whitelist enforced in
# save_plot_mm_gated() above. Removed to avoid two competing manifests.)

# panels = vpd deciles
save_full_and_cap(p_day_gs_ta_by_vpd_full,    p_day_gs_ta_by_vpd_cap,
                  "FigS_day_Gs_vs_Ta_split_by_VPDdeciles_full_range",
                  "FigX_day_Gs_vs_Ta_split_by_VPDdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_ta_by_vpd_full,  p_night_gs_ta_by_vpd_cap,
                  "FigS_night_Gs_vs_Ta_split_by_VPDdeciles_full_range",
                  "FigX_night_Gs_vs_Ta_split_by_VPDdeciles_CAP_0p10")

save_full_and_cap(p_day_gs_ts_by_vpd_full,    p_day_gs_ts_by_vpd_cap,
                  "FigS_day_Gs_vs_Ts_split_by_VPDdeciles_full_range",
                  "FigX_day_Gs_vs_Ts_split_by_VPDdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_ts_by_vpd_full,  p_night_gs_ts_by_vpd_cap,
                  "FigS_night_Gs_vs_Ts_split_by_VPDdeciles_full_range",
                  "FigX_night_Gs_vs_Ts_split_by_VPDdeciles_CAP_0p10")

save_full_and_cap(p_day_gs_swc_by_vpd_full,   p_day_gs_swc_by_vpd_cap,
                  "FigS_day_Gs_vs_SWC_split_by_VPDdeciles_full_range",
                  "FigX_day_Gs_vs_SWC_split_by_VPDdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_swc_by_vpd_full, p_night_gs_swc_by_vpd_cap,
                  "FigS_night_Gs_vs_SWC_split_by_VPDdeciles_full_range",
                  "FigX_night_Gs_vs_SWC_split_by_VPDdeciles_CAP_0p10")

# panels = ta deciles
save_full_and_cap(p_day_gs_vpd_by_ta_full,    p_day_gs_vpd_by_ta_cap,
                  "FigS_day_Gs_vs_VPD_split_by_TAdeciles_full_range",
                  "FigX_day_Gs_vs_VPD_split_by_TAdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_vpd_by_ta_full,  p_night_gs_vpd_by_ta_cap,
                  "FigS_night_Gs_vs_VPD_split_by_TAdeciles_full_range",
                  "FigX_night_Gs_vs_VPD_split_by_TAdeciles_CAP_0p10")

save_full_and_cap(p_day_gs_ts_by_ta_full,     p_day_gs_ts_by_ta_cap,
                  "FigS_day_Gs_vs_Ts_split_by_TAdeciles_full_range",
                  "FigX_day_Gs_vs_Ts_split_by_TAdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_ts_by_ta_full,   p_night_gs_ts_by_ta_cap,
                  "FigS_night_Gs_vs_Ts_split_by_TAdeciles_full_range",
                  "FigX_night_Gs_vs_Ts_split_by_TAdeciles_CAP_0p10")

save_full_and_cap(p_day_gs_swc_by_ta_full,    p_day_gs_swc_by_ta_cap,
                  "FigS_day_Gs_vs_SWC_split_by_TAdeciles_full_range",
                  "FigX_day_Gs_vs_SWC_split_by_TAdeciles_CAP_0p10")

save_full_and_cap(p_night_gs_swc_by_ta_full,  p_night_gs_swc_by_ta_cap,
                  "FigS_night_Gs_vs_SWC_split_by_TAdeciles_full_range",
                  "FigX_night_Gs_vs_SWC_split_by_TAdeciles_CAP_0p10")

message("done. wrote:")
message(paste0("  ", out_csv))
message(paste0("  ", high_csv))
message(paste0("  ", sum_csv))
message(paste0("  ", rank_csv))
message(paste0("  old two-panel: ", file.path(paper_figures_path, "FigX_Gs_two_panel_CAP_0p10.png")))
message(paste0("  old two-panel: ", file.path(paper_figures_path, "FigS_Gs_two_panel_full_range.png")))
message(paste0("  diag plots in: ", diag_path))

if (is.na(col_TS))  message('note: ts column not found; gs vs ts figures will be empty/NA unless you map TS_3 (or set ts_candidates).')
if (is.na(col_SWC)) message('note: swc column not found; gs vs swc figures will be empty/NA unless you map SWC_1_1_1 (or set swc_candidates).')



# ───────────────────────────────────────────────────────────────────────────────
# NEW FIGURE: Daytime only, physically filtered (VPD, USTAR, LE)
# ───────────────────────────────────────────────────────────────────────────────
dat_gs_day_phys_filtered <- dat_diag %>%
  filter(
    is_day == TRUE,      # Limit to daytime only
    is.finite(Gs_mps),   # Must have a valid computed number
    Gs_mps >= 0,         # Ignore physically impossible negative conductance
    is.finite(VPD),
    is.finite(Tair),
    # --- The New Physical Filters ---
    VPD > 0.1,           # Avoid near-zero VPD causing math artifacts
    USTAR > 0.2,         # Ensure sufficient boundary-layer mixing
    LE > 10              # Avoid near-zero latent heat flux
  )

# Create the VPD plot (physically filtered daytime)
p1_day_phys <- ggplot(dat_gs_day_phys_filtered, aes(x = VPD, y = Gs_mps, color = Tair)) +
  geom_point(alpha = 0.35, size = 1.3) +
  stat_summary_bin(fun = median, bins = median_bins, geom = "line", linewidth = 0.9, color = "black") +
  scale_color_viridis_c(option = "turbo", limits = c(ta_col_lo, ta_col_hi), oob = scales::squish) +
  labs(x = "VPD (kPa)", y = NULL, color = expression(italic(T)[plain(a)] ~ "(°C)"), tag = "(a)") +
  base_theme_old +
  tag_theme +
  theme(plot.tag.position = c(0.09, 1)) # <-- Nudge tag right (x) and keep it at the top (y)

# Create the Ta plot (physically filtered daytime)
p2_day_phys <- ggplot(dat_gs_day_phys_filtered, aes(x = Tair, y = Gs_mps, color = VPD)) +
  geom_point(alpha = 0.35, size = 1.3) +
  # median line only where Ta is well sampled (>= ta_line_min); points below stay shown
  stat_summary_bin(
    data = dplyr::filter(dat_gs_day_phys_filtered, Tair >= ta_line_min),
    fun = median, bins = median_bins, geom = "line", linewidth = 0.9, color = "black"
  ) +
  scale_color_viridis_c(option = "turbo", limits = c(0, vpd_col_hi), oob = scales::squish) +
  labs(x = expression(italic(T)[plain(a)] ~ "(°C)"), y = NULL, color = "VPD (kPa)", tag = "(b)") +
  base_theme_old +
  tag_theme +
  theme(plot.tag.position = c(0.09, 1)) # <-- Nudge tag right (x) and keep it at the top (y)

# Combine them using patchwork
fig_gs_day_phys <- (patchwork::wrap_elements(full = y_grob) | (p1_day_phys | p2_day_phys)) +
  patchwork::plot_layout(widths = c(0.055, 1))

# Save the plot
ggsave(file.path(paper_figures_path, "FigS_Gs_two_panel_DAY_PHYSICAL_FILTER.png"),
       fig_gs_day_phys, width = 7.2, height = 3.4, dpi = 600, bg = "white")
