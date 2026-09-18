#!/usr/bin/env Rscript
# =============================================================================
# zz_build_zenodo_deposit.R
# -----------------------------------------------------------------------------
# Assemble the Zenodo DATA deposit for the Global Change Biology submission.
#
# Scope decision (see writing/ and the project plan): deposit the AGGREGATED
# processed products that underlie the figures and tables ONLY. This satisfies
# GCB's "data underlying the figures and results" requirement while excluding
# the raw/processed half-hourly series and the calibration files (per co-author
# guidance). The complete half-hourly eddy-covariance record is released
# separately at AmeriFlux (site PE-TNR) upon publication.
#
# INCLUDED : 48-point diel composites, light-response parameter tables,
#            daily/16-day/monthly carbon-balance & meteorology summaries,
#            the per-year coverage and energy-balance closure tables
#            (Tables S2/S3), the storage-term comparison behind Fig S10, the
#            morning-transition, partitioning and VPD-constraint result tables,
#            and the flux-site metadata used for the map (Fig 1).
# EXCLUDED : dataset_from_matlab*.csv (half-hourly), IMERG (third-party),
#            WorldClim (third-party; redistribution is expressly forbidden by
#            its own licence), *calibration*, *preliminary*/*FConly*
#            (not-for-citation), resampled_data, half-hourly soil / ENSO
#            products, internal diagnostics.
#
# NOT REPRODUCIBLE FROM THIS DEPOSIT, and named as such in the paper's Data
# Availability Statement: Fig S1 (site locator) and Fig S11 (footprints), both
# produced outside this repository; Figs S7/S8 (storage DNN skill and
# permutation importance), produced externally by AER; Figs S6 and S10 and the
# Section 2.6 partitioning counts, which need the half-hourly series.
#
# Run:     "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" code/zz_build_zenodo_deposit.R
#   (or)   Rscript code/zz_build_zenodo_deposit.R
# Output:  <repo>/data_release/   (gitignored; upload its contents to Zenodo)
# =============================================================================

## ---- locate repo root (standalone; does not source paths.R) ----------------
get_script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(normalizePath(f))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
  normalizePath(file.path(getwd(), "code", "zz_build_zenodo_deposit.R"))
}
script_path <- tryCatch(get_script_path(), error = function(e) normalizePath("."))
code_dir    <- dirname(script_path)
repo_root   <- dirname(code_dir)
data_dir    <- file.path(repo_root, "data")
output_dir  <- file.path(repo_root, "output")
out_dir     <- file.path(repo_root, "data_release")

cat("repo_root :", repo_root, "\n")
cat("data_dir  :", data_dir, "\n")
cat("out_dir   :", out_dir, "\n\n")

## ---- manifest: source file -> subfolder + provenance/notes -----------------
# `src` is relative to data_dir unless it starts with "output/".
M <- function(src, subdir, supports) data.frame(src = src, subdir = subdir,
                                                 supports = supports, stringsAsFactors = FALSE)
manifest <- do.call(rbind, list(
  # --- diel composites (48 half-hour bins per month x season x ENSO) --------
  M("tambopata_48points_per_month.csv",         "diel_composites", "Fig 2a, Fig 4, Fig 6 (core 48-point diel composite)"),
  M("tambopata_48points_per_month_measured.csv", "diel_composites", "measured-NEE variant; Fig 4, Fig 5, Fig 6, Fig S17 inputs (and the archived decile panels)"),
  M("tambopata_48points_per_month_with_gs.csv",  "diel_composites", "Fig 7 (adds canopy conductance Gs)"),
  M("diel_monthly_friagem_NEEok.csv",            "diel_composites", "friagem (cold-spell) diel composite - NEE"),
  M("diel_monthly_friagem_GEP.csv",              "diel_composites", "friagem diel composite - GEP"),
  M("diel_monthly_friagem_Reco.csv",             "diel_composites", "friagem diel composite - Reco"),
  # --- light-response parameters --------------------------------------------
  # NOTE: the plain `_by_year.csv` is byte-identical to `_by_year_with_Ts.csv`
  # (md5 9f1030ed...), so only the _with_Ts copy is deposited.
  M("light_response_mean_env_and_fit_params_by_year_with_Ts.csv", "light_response", "per-year light-response params + env (with Ts); inputs to the Section 3.5 mixed-effects model selection"),
  M("light_response_mean_environmental_and_fit_params.csv",       "light_response", "pooled season x ENSO params + bootstrap CIs; Table S6"),
  M("season_enso_means_phi0_Pmax_Rd_LCP.csv",                     "light_response", "season x ENSO mean light-response params"),
  M("season_enso_light_response_P2000_bootstrap_CI.csv",          "light_response", "season x ENSO params with bootstrap CIs; Table S6"),
  M("light_response_halfhourly_sensitivity_by_enso_season.csv",   "light_response", "Table S7 (hyperbola refitted to raw half-hourly observations; Section 2.9 sensitivity check)"),
  M("light_response_residual_line_status.csv",                    "light_response", "residual-regression fit status (which lines are shown)"),
  # --- carbon balance: daily / 16-day / monthly -----------------------------
  M("biomet_16day_timeseries_for_Fig3.csv",              "carbon_balance", "Fig 3 biomet panels (16-day)"),
  M("nee_dailyC_16day_by_ENSO_season.csv",               "carbon_balance", "Table 1 (daily C balance by ENSO/season)"),
  M("nee_dailyC_16day_by_ENSO_season_presentation.csv",  "carbon_balance", "Table 1 (presentation form)"),
  M("gep_reco_gs_16day_by_ENSO_season.csv",              "carbon_balance", "Table 3 (GEP/Reco/Gs by ENSO/season)"),
  M("gep_reco_gs_16day_by_ENSO_season_presentation.csv", "carbon_balance", "Table 3 (presentation form)"),
  M("gs_16day_timeseries.csv",                           "carbon_balance", "Fig S13c (16-day canopy conductance)"),
  M("photosynthetic_capacity_pc_daily_strict.csv",       "carbon_balance", "daily photosynthetic capacity (strict coverage)"),
  M("photosynthetic_capacity_pc_16day_strict.csv",       "carbon_balance", "16-day photosynthetic capacity (strict coverage)"),
  M("check_16day_dailyCbalance_checked.csv",             "carbon_balance", "Fig 3 daily C-balance (16-day) with diel-coverage QC"),
  M("tambopata_monthly_diurnal_sum_gC_m2_day.csv",       "carbon_balance", "monthly carbon balance (g C m-2 d-1)"),
  M("gpp_daily_mean_by_year_month.csv",                  "carbon_balance", "monthly GPP means by year-month"),
  M("gpp_daily_mean_climatology_by_calendar_month.csv",  "carbon_balance", "GPP monthly climatology"),
  # --- meteorology summaries / statistical tests ----------------------------
  M("env_ENSO_monthly_composites.csv",                        "meteorology", "Fig 2b (monthly environmental composites by ENSO)"),
  M("vpd_3pm_quantile_0p984.csv",                             "meteorology", "98.4th-percentile 15:00 VPD (text result)"),
  M("meteorological_tests_wet_vs_dry.csv",                    "meteorology", "wet-vs-dry meteorology tests"),
  M("meteorological_tests_ENSO_within_season_ANOVA.csv",      "meteorology", "ENSO-within-season ANOVA"),
  M("meteorological_tests_ENSO_within_season_TukeyHSD.csv",   "meteorology", "ENSO-within-season Tukey HSD"),
  M("GEP_mean_env_by_year_season_ENSO.csv",                   "meteorology", "GEP + environment by year/season/ENSO"),
  M("table_16day_window_sign_by_ENSO.csv",               "carbon_balance", "Section 3.1 window-sign counts (16 La Nina windows all sources; El Nino 2018-19 vs 2023-24 split)"),
  M("check_16day_dailyC_diel_coverage.csv",              "carbon_balance", "the 72 window x season x ENSO groups screened down to the 35 retained (Section 2.6)"),
  M("output/table_cumulative_diel_NEE.csv",              "carbon_balance", "Fig 2b (cumulative mean diel NEE by season and ENSO)"),
  # --- per-year coverage and energy-balance closure -------------------------
  M("output/table_data_coverage_energy_balance.csv", "coverage_and_closure", "Tables S2 and S3 (measured NEE retained per year after QA/QC + u* filtering; annual energy-balance closure)"),
  # --- CO2 storage term: profile-measured vs DNN-reconstructed --------------
  M("storage_diel_profiler_vs_model_2024.csv",  "storage_term", "Fig S10 (mean diel NEE with measured vs modelled storage, 2024)"),
  M("storage_profiler_vs_model_stats_2024.csv", "storage_term", "Section 2.5 storage agreement statistics (n = 2468, slope, r by day/night, diel curve gap)"),
  # --- morning source-to-sink transition (Section 3.1) ----------------------
  M("output/table_morning_transition.csv",                "morning_transition", "Section 3.1 transition timing by ENSO phase, month-block bootstrap"),
  M("output/table_morning_transition_by_episode.csv",     "morning_transition", "transition timing per ENSO episode"),
  M("output/table_morning_transition_halfhourly.csv",     "morning_transition", "mean diel NEE around the morning transition, by half-hour"),
  M("output/table_morning_transition_light_controls.csv", "morning_transition", "morning PAR controls for the transition comparison"),
  M("output/table_morning_transition_grid.csv",           "morning_transition", "share of days transitioned and of half-hours in uptake at fixed clock times"),
  M("morning_transition_per_day.csv",                     "morning_transition", "per-day transition half-hour (the unit the medians are taken over)"),
  # --- flux partitioning diagnostics ----------------------------------------
  M("output/table_gep_reco_self_correlation.csv", "partitioning", "Section 2.6 GEP-Reco self-correlation against the partitioning-identity null (r, null median, 95% range, p)"),
  # --- VPD / light constraint on high-light uptake (Section 3.2, Fig 5) -----
  M("output/table_lr_residual_by_light_and_season.csv", "vpd_light_constraint", "Fig 5a light-response residual slopes vs Ta and VPD, by light class and season"),
  M("output/table_vpd_class_light_response.csv",        "vpd_light_constraint", "Fig 5b hyperbolae fitted within fixed VPD classes"),
  M("output/table_reco_vs_gross_decomposition.csv",     "vpd_light_constraint", "Section 3.2 decomposition of the residual against independently estimated Reco"),
  M("output/table_gross_vs_respiration.csv",            "vpd_light_constraint", "Section 3.2 split of P2000 into gross assimilation and Rd"),
  # --- meteorology summaries / statistical tests ----------------------------
  M("output/table_monthly_precipitation.csv", "meteorology", "monthly IMERG precipitation totals (Discussion 2024 dry-season context)"),
  M("output/table_canopy_temperature.csv",    "meteorology", "canopy-minus-air temperature by day/night/high-light, with the Stefan-Boltzmann cross-check (Discussion)"),
  # --- flux-site metadata (Fig 1 map) ---------------------------------------
  M("panamazon_flux_sites.csv",           "site_metadata", "Fig 1 pan-Amazon flux-site list (scraped from panamazonflux.org)"),
  M("all_known_flux_sites_map_export.csv", "site_metadata", "Fig 1 flux-site map export")
))

## ---- files that must NEVER appear in the deposit (safety net) ---------------
exclude_patterns <- c("dataset_from_matlab", "calibration", "preliminary",
                       "FConly", "resampled_data", "soil_water_content",
                       "enso_conditions", "imerg", "high_gs", "2024_only")

## ---- curated column dictionary ---------------------------------------------
# unit + description keyed by exact column name; unmatched columns fall through
# to rule-based defaults in describe_col().
dict <- read.csv(text = "column,unit,description
ym,,Year-month label (YYYY-MM)
season,,Season: wet or dry
ENSO,,ENSO phase category (El Nino / La Nina / Neutral; multi-year runs distinguished)
hour,half-hour-of-day,Half-hour-of-day bin (48 bins per day) in local time
NEE_ok,umol CO2 m-2 s-1,Net ecosystem exchange - quality-screened MEASURED NEE (primary); gaps re-invalidated
NEE_meas,umol CO2 m-2 s-1,Net ecosystem exchange - measured (u*-filtered)
GEP,umol CO2 m-2 s-1,Gross ecosystem photosynthesis (= GPP) from nighttime-respiration partitioning
Reco,umol CO2 m-2 s-1,Ecosystem respiration
FC,umol CO2 m-2 s-1,Turbulent CO2 flux (measured)
TA_1_1_1,degC,Air temperature
Tair,degC,Air temperature
Ta,degC,Air temperature (group mean)
VPD_kPa,kPa,Vapor pressure deficit
VPD,kPa,Vapor pressure deficit (group mean)
PAR,umol m-2 s-1,Photosynthetically active radiation (PPFD)
PAR_corrected_SWin,umol m-2 s-1,PAR derived/corrected from incoming shortwave radiation
PAR_source,,Provenance flag for the PAR value used
PPFD_IN_1_1_1,umol m-2 s-1,Incoming photosynthetic photon flux density
SW_IN,W m-2,Incoming shortwave radiation
SWin,W m-2,Incoming shortwave radiation
Rn,W m-2,Net radiation
H,W m-2,Sensible heat flux
LE,W m-2,Latent heat flux
SWC_1_1_1,cm,Soil water content - integrated 0-100 cm storage (cm; NOT volumetric)
SWC,cm,Soil water content - integrated 0-100 cm storage (cm; NOT volumetric)
TS_3,degC,Soil temperature (probe 3)
TS,degC,Soil temperature
Ts,degC,Soil temperature (group mean)
USTAR,m s-1,Friction velocity
u_star,m s-1,Friction velocity (group mean)
Ustar,m s-1,Friction velocity (group mean)
WS,m s-1,Wind speed
WS_1_1_1,m s-1,Wind speed
WIND,m s-1,Wind speed
WD_1_1_1,degrees,Wind direction (0-360; meteorological convention)
Gs_mps,m s-1,Canopy/surface conductance to water vapor (Penman-Monteith inversion)
Gs_mean,m s-1,Mean canopy conductance
Gs_sd,m s-1,SD of canopy conductance
phi0,mol CO2 (mol photons)-1,Apparent quantum yield (initial slope of light-response curve)
Pmax,umol CO2 m-2 s-1,Maximum gross photosynthesis (light-response asymptote)
P2000,umol CO2 m-2 s-1,Gross photosynthesis at PAR = 2000 umol m-2 s-1
Rd,umol CO2 m-2 s-1,Respiration intercept from light-response fit
LCP,umol m-2 s-1,Light compensation point (PAR)
Group,,Season x ENSO grouping label
Year,,Calendar year
year,,Calendar year
month,,Calendar month (1-12)
month_name,,Calendar month name
precip_mm_day,mm d-1,Precipitation
gC_m2_day,g C m-2 d-1,Carbon balance (positive = source)
carbon_balance,,Monthly carbon-balance label (source/sink)
dailyC_mean,g C m-2 d-1,Mean daily carbon balance in the window
dailyC_sd,g C m-2 d-1,SD of daily carbon balance in the window
Pc_gC_m2_day,g C m-2 d-1,Daily photosynthetic capacity
GPP_mean_gC_m2_day,g C m-2 d-1,Mean daily GPP
GPP_sd_gC_m2_day,g C m-2 d-1,SD of daily GPP
GPP_mean_climatology_gC_m2_day,g C m-2 d-1,Climatological mean GPP for the calendar month
GPP_sd_across_year_months_gC_m2_day,g C m-2 d-1,SD of GPP across year-months
GEP_gC_m2_day,g C m-2 d-1,Mean daily GEP
GEP_sd_gC_m2_day,g C m-2 d-1,SD of daily GEP
window_start,,16-day window start date
window_mid,,16-day window midpoint date
Date,,Date (YYYY-MM-DD)
DateTime,,Timestamp (local time)
Month,,Calendar month (1-12)
Month_label,,Calendar month label
frac_diel_covered,fraction,Fraction of the 48 diel bins with data
mean_diel_coverage_NEE,fraction,Mean diel coverage across windows
missing_bins,count,Number of missing diel bins
value,,Value (see variable/source_column)
mean_value,,Mean value
mean,,Mean
sd,,Standard deviation
variable,,Variable identifier
variable_label,,Human-readable variable label
source_column,,Source column name in the half-hourly dataset
panel_label,,Figure-panel label
value_type,,Type of value (e.g. mean/sum)
label,,Label
unit,,Unit string (as plotted)
vpd_quantile_kPa,kPa,VPD at the requested quantile
prob,probability,Quantile probability (e.g. 0.984)
tz_local,,Local timezone used
statistic,,Test statistic
p_value,,p-value
p_adj_BH,,Benjamini-Hochberg adjusted p-value
signif,,Significance flag (raw p)
signif_adj,,Significance flag (adjusted p)
n_boot_success,count,Number of successful bootstrap iterations
Season,,Season: wet or dry
month_lbl,,Month label
friagem_lbl,,Friagem (cold-spell) category label
parameter,,Light-response parameter name
base_predictor,,Base predictor variable in the residual regression
xvar,,Predictor variable on the x-axis
line_shown,,Whether a fitted regression line is displayed (TRUE/FALSE)
r_squared,,Coefficient of determination (R-squared)
test,,Statistical test name
comparison,,Groups being compared
note,,Free-text note
df,,Degrees of freedom
df_enso,,Degrees of freedom (ENSO factor)
df_residual,,Residual degrees of freedom
diff,,Tukey HSD mean difference between groups
lower,,Tukey HSD lower confidence bound
upper,,Tukey HSD upper confidence bound
wet_mean,,Wet-season mean of the variable
dry_mean,,Dry-season mean of the variable
exact_1500,,Flag: restrict to exactly 15:00 local time
restrict_2012_2017,,Flag: restrict to the 2012-2017 subset
is_day,,Daytime flag (TRUE/FALSE)
day_by_hour,,Day/night flag derived from hour-of-day
day_by_par,,Day/night flag derived from a PAR threshold
hour_hms,,Hour-of-day as HH:MM:SS
min_year,,First year contributing to the composite
max_year,,Last year contributing to the composite
mean_valid_days_per_year_month,count,Mean valid days per year-month
ymin,,Lower ribbon bound (mean - SD)
ymax,,Upper ribbon bound (mean + SD)
NEE_mean,g C m-2 d-1,Mean daily NEE carbon balance
NEE_sd,g C m-2 d-1,SD of daily NEE carbon balance
GEP_mean,g C m-2 d-1,Mean daily GEP
GEP_sd,g C m-2 d-1,SD of daily GEP
Reco_mean,g C m-2 d-1,Mean daily Reco
Reco_sd,g C m-2 d-1,SD of daily Reco
NEE daily C balance (g C m-2 d-1),g C m-2 d-1,Mean daily NEE carbon balance
GEP (g C m-2 d-1),g C m-2 d-1,Mean daily GEP
Reco (g C m-2 d-1),g C m-2 d-1,Mean daily Reco
Gs (m s-1),m s-1,Mean canopy conductance
mean diel coverage,fraction,Mean fraction of the 48 diel bins covered
mean valid days,count,Mean number of valid days per window
min valid days,count,Minimum number of valid days per window
id,,Flux-site identifier
name,,Flux-site name
acronym,,Flux-site acronym
site_id,,Flux-site identifier
site_name,,Flux-site name
country,,Country
state,,State or region
vegetation,,Vegetation type
status,,Site status
pi,,Site principal investigator
data_availability,,Flux-site data-availability status
seed,,Random seed used for the bootstrap or permutation
subset,,Sample the statistic was computed over (all / day / night or equivalent)
slope,,Ordinary least-squares slope
intercept,,Ordinary least-squares intercept
r,,Pearson correlation coefficient
R2,,Coefficient of determination
RMSE,,Root mean squared error (on the diel-curve row this holds the maximum absolute difference across the 48 bins)
bias,,Mean signed difference (on the diel-curve row this holds the mean absolute difference across the 48 bins)
mean obs per diel bin,count,Mean number of observations per half-hour-of-day bin
mean half-hours per window,count,Mean number of half-hourly observations per retained window
obs_per_bin,count,Mean number of observations per half-hour-of-day bin
ustar_threshold,m s-1,Friction-velocity threshold applied to that year
ustar_source,,Whether the u* threshold was estimated for that year or taken from the pooled fallback
pct_meas_of_record,%,Measured NEE as a percentage of the half-hours spanned by the processed record that year
pct_meas_of_year,%,Measured NEE as a percentage of all half-hours in the calendar year
pct_meas_day,%,Measured NEE as a percentage of daytime half-hours
pct_meas_night,%,Measured NEE as a percentage of nighttime half-hours
pct_meas_twilight,%,Measured NEE as a percentage of twilight half-hours (Rg between 10 and 44 W m-2)
pct_gapfill_produced,%,Gap-filled NEE produced by REddyProc as a percentage of the record
pct_gapfilled_of_record,%,Gap-filled NEE still usable after the 24-half-hour gap-length screen
EBR,,Energy-balance ratio (H + LE) / Rn over valid half-hours
EBR_RnG,,Energy-balance ratio closed against available energy (Rn - G)
EBR_Rn_Gsubset,,Energy-balance ratio against Rn alone over only the half-hours where G is available
turb_bar,W m-2,Mean turbulent flux (H + LE)
Rn_bar,W m-2,Mean net radiation
G_bar,W m-2,Mean soil heat flux
G_mean_day,W m-2,Mean daytime soil heat flux
g_usable,,Whether soil heat flux is available for that year
slope_all,,OLS slope of (H + LE) on Rn over all valid half-hours
slope_day,,OLS slope of (H + LE) on Rn over daytime half-hours
slope_day_RnG,,OLS slope of (H + LE) on (Rn - G) over daytime half-hours
slope_day_Rn_Gsubset,,OLS slope of (H + LE) on Rn over the G-available daytime subset
intercept_all,W m-2,OLS intercept over all valid half-hours
intercept_day,W m-2,OLS intercept over daytime half-hours
r2_all,,Coefficient of determination over all valid half-hours
r2_day,,Coefficient of determination over daytime half-hours
r2_day_RnG,,Coefficient of determination for the daytime (Rn - G) fit
ratio,,Ratio reported for the closure comparison on that row
FC_mean,umol CO2 m-2 s-1,Mean turbulent CO2 flux in the half-hour-of-day bin
FC_se,umol CO2 m-2 s-1,Standard error of mean FC in the bin
NEE_profiler_mean,umol CO2 m-2 s-1,Mean NEE using the profile-measured storage term
NEE_profiler_se,umol CO2 m-2 s-1,Standard error of mean profiler-based NEE
NEE_model_mean,umol CO2 m-2 s-1,Mean NEE using the DNN-reconstructed storage term
NEE_model_se,umol CO2 m-2 s-1,Standard error of mean model-based NEE
SC_profiler_mean,umol CO2 m-2 s-1,Mean profile-measured storage flux in the bin
SC_profiler_se,umol CO2 m-2 s-1,Standard error of mean profiler storage
SC_model_mean,umol CO2 m-2 s-1,Mean DNN-reconstructed storage flux in the bin
SC_model_se,umol CO2 m-2 s-1,Standard error of mean modelled storage
ustar_filtered,,Whether the u* filter was applied before the comparison
estimator,,Which morning-transition estimator the row reports
primary,,Whether this row is the estimator reported in the paper
nee_series,,Which NEE series the transition was detected on
mean_el_nino_h,h,Mean transition time under El Nino in decimal hours
mean_la_nina_h,h,Mean transition time under La Nina in decimal hours
mean_el_nino_hm,,Mean El Nino transition time as HH:MM
mean_la_nina_hm,,Mean La Nina transition time as HH:MM
diff_min,min,La Nina minus El Nino difference in transition time
diff_ci_low_min,min,Lower bound of the month-block bootstrap interval on the difference
diff_ci_high_min,min,Upper bound of the month-block bootstrap interval on the difference
month_level_t_p,,p-value of a month-level t-test on the transition time
episode,,ENSO episode label
mean_t,h,Mean transition time in decimal hours
mean_hm,,Mean transition time as HH:MM
mean_PAR_at_crossing,umol photons m-2 s-1,Mean PAR in the transition half-hour
mean_PAR_morn,umol photons m-2 s-1,Mean morning (06:00-09:00) PAR
min_earlier_than_la_nina_mean,min,Minutes by which this episode transitions earlier than the La Nina mean
half_hour,,Clock time of the fixed comparison point
share_transitioned_el_nino,fraction,Share of El Nino days that have transitioned by this time
share_transitioned_la_nina,fraction,Share of La Nina days that have transitioned by this time
share_transitioned_diff,fraction,El Nino minus La Nina difference in transitioned share
share_uptake_el_nino,fraction,Share of observed El Nino half-hours in net uptake at this time
share_uptake_la_nina,fraction,Share of observed La Nina half-hours in net uptake at this time
share_uptake_diff,fraction,El Nino minus La Nina difference in uptake share
NEE_El_Nino,umol CO2 m-2 s-1,Mean El Nino NEE at this half-hour of day
NEE_La_Nina,umol CO2 m-2 s-1,Mean La Nina NEE at this half-hour of day
PAR_El_Nino,umol photons m-2 s-1,Mean El Nino PAR at this half-hour of day
PAR_La_Nina,umol photons m-2 s-1,Mean La Nina PAR at this half-hour of day
PAR_pct_el_vs_la,%,El Nino PAR as a percentage of La Nina PAR at this half-hour
hour_dec,h,Half-hour of day in decimal hours
local_time,,Half-hour of day as HH:MM local time (America/Lima)
sample,,Sample the light-matched comparison was computed over
par_ref,umol photons m-2 s-1,Reference PAR level for the light-matched comparison
mean_el_nino,,Mean value under El Nino
mean_la_nina,,Mean value under La Nina
diff_la_minus_el,,La Nina minus El Nino difference
excludes_zero,,Whether the bootstrap interval excludes zero
t,h,Transition half-hour for that day in decimal hours
t_first,h,First half-hour below zero before the persistence rule is applied
transition_hm,,Transition half-hour as HH:MM
sunrise,,Sunrise time for that day
t_since_sunrise,h,Hours between sunrise and the transition
PAR_at_crossing,umol photons m-2 s-1,PAR in the transition half-hour
PAR_morn,umol photons m-2 s-1,Mean morning (06:00-09:00) PAR for that day
NEE_morn,umol CO2 m-2 s-1,Mean morning (06:00-09:00) NEE for that day
reverts,,Whether NEE returns above zero later that morning after the transition
nee_term,,Which NEE-derived term the correlation is computed on
reco_term,,Which Reco estimate the correlation is computed against
obs_corr,,Observed correlation between Reco and absolute GEP
null_median,,Median correlation under the partitioning-identity null
null_lo,,2.5th percentile of the null distribution
null_hi,,97.5th percentile of the null distribution
p_obs_gt_null,,Proportion of null replicates at or above the observed correlation
block,,Block used for the bootstrap resampling
difference,,Difference between the two groups contrasted on that row
Ta_slope,umol CO2 m-2 s-1 per degC,Slope of the light-response residual on air temperature
Ta_slope_NEE,umol CO2 m-2 s-1 per degC,Residual slope on air temperature expressed on the NEE sign convention
Ta_R2,,Coefficient of determination for the air-temperature residual fit
Ta_p,,p-value for the air-temperature residual slope
VPD_slope,umol CO2 m-2 s-1 per kPa,Slope of the light-response residual on VPD
VPD_slope_NEE,umol CO2 m-2 s-1 per kPa,Residual slope on VPD expressed on the NEE sign convention
VPD_R2,,Coefficient of determination for the VPD residual fit
VPD_p,,p-value for the VPD residual slope
vpd_class,kPa,VPD class within which the hyperbola was fitted
mean_VPD,kPa,Mean VPD of the observations in the class
mean_Ta,degC,Mean air temperature of the observations in the class
mean_hour,h,Mean half-hour of day of the observations in the class
P2000_lo,umol CO2 m-2 s-1,Lower bootstrap bound on P2000 for the class
P2000_hi,umol CO2 m-2 s-1,Upper bootstrap bound on P2000 for the class
driver,,Driver the light-response residual was decomposed against
residual_slope,,Slope of the light-response residual on the driver
Reco_slope,,Slope of the independently estimated Reco on the same driver
gross_slope,,Slope attributable to gross assimilation once Reco is removed
Reco_pct_of_residual,%,Share of the residual slope accounted for by Reco
Ag2000,umol CO2 m-2 s-1,Gross assimilation at PAR = 2000
Ag_lo,umol CO2 m-2 s-1,Lower bootstrap bound on Ag2000
Ag_hi,umol CO2 m-2 s-1,Upper bootstrap bound on Ag2000
Rd_over_Ag,,Ratio of dark respiration to gross assimilation at PAR = 2000
ratio_lo,,Lower bootstrap bound on the Rd/Ag ratio
ratio_hi,,Upper bootstrap bound on the Rd/Ag ratio
within_corr_Ag_Rd,,Bootstrap correlation between Ag2000 and Rd within a regime
across_regime_corr_Ag_Rd,,Correlation between Ag2000 and Rd across regimes
P2000_point,umol CO2 m-2 s-1,Point estimate of P2000 from the half-hourly refit
Pmax_point,umol CO2 m-2 s-1,Point estimate of the fitted asymptote from the half-hourly refit
series,,Which series the row belongs to
row_type,,Whether the row is a diel value or a summary marker
nee_umol_m2_s,umol CO2 m-2 s-1,Mean NEE in the half-hour-of-day bin
cum_gC_m2,g C m-2,Cumulative carbon since 00:30 local time
quantity,,Quantity the row reports
source,,Measurement or derivation the row is based on
mean_C,degC,Mean canopy-minus-air temperature difference
median_C,degC,Median canopy-minus-air temperature difference
sd_C,degC,Standard deviation of the canopy-minus-air difference
p05_C,degC,5th percentile of the canopy-minus-air difference
p95_C,degC,95th percentile of the canopy-minus-air difference
frac_canopy_warmer,fraction,Fraction of half-hours in which the canopy is warmer than the air
RMSE_C,degC,Root mean squared error against the longwave-derived canopy temperature
bias_C,degC,Mean signed difference against the longwave-derived canopy temperature
coverage_start,,First date of the paired record
coverage_end,,Last date of the paired record
precip_mm,mm,Monthly precipitation total from IMERG
days_in_month,count,Number of days in the month
grouping,,Grouping over which the window-sign counts are tabulated
mean_dailyC,g C m-2 d-1,Mean daily carbon balance across the windows in the group
min_dailyC,g C m-2 d-1,Most sink-like window in the group
max_dailyC,g C m-2 d-1,Most source-like window in the group
", stringsAsFactors = FALSE, colClasses = "character")

describe_col <- function(col) {
  key <- trimws(col)
  hit <- match(key, dict$column)
  if (!is.na(hit)) return(c(unit = dict$unit[hit], description = dict$description[hit]))
  kl <- tolower(key)
  # rule-based fallbacks
  if (grepl("_low$|_high$|_ci$|_CI$", key)) return(c(unit = "", description = "Bootstrap confidence-interval bound (see base parameter)"))
  if (grepl("gc_m2_day", kl))              return(c(unit = "g C m-2 d-1", description = "Carbon flux"))
  if (grepl("^n_|_n_|(^|_)n$|_n$|n_obs|n_meas|n_days|n_bins|n_points|n_windows|n_halfhour|n_hh|n_years|n_monthly", kl))
                                           return(c(unit = "count", description = "Sample count"))
  if (grepl("(^| )N$", key))               return(c(unit = "count", description = "Sample count"))
  if (grepl("lat|latitude", kl))           return(c(unit = "degrees_north", description = "Site latitude"))
  if (grepl("lon|longitude", kl))          return(c(unit = "degrees_east", description = "Site longitude"))
  if (grepl("gc m-2 d-1", kl))             return(c(unit = "g C m-2 d-1", description = "Carbon flux"))
  c(unit = "", description = "")  # unmatched -> fill in manually
}

## ---- helpers ---------------------------------------------------------------
write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(lines, con)
}
read_header <- function(path) {
  h <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  scan(text = h, what = character(), sep = ",", quote = "\"", quiet = TRUE, strip.white = TRUE)
}
count_rows <- function(path) max(0L, length(readLines(path, warn = FALSE)) - 1L)

## ---- (re)create staging tree -----------------------------------------------
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, showWarnings = FALSE)
for (sd in unique(manifest$subdir)) dir.create(file.path(out_dir, sd), showWarnings = FALSE)

## ---- copy + collect manifest/dictionary rows -------------------------------
missing <- character(0)
man_rows <- list(); dict_rows <- list()
for (i in seq_len(nrow(manifest))) {
  src <- manifest$src[i]
  src_path <- if (startsWith(src, "output/")) file.path(output_dir, sub("^output/", "", src)) else file.path(data_dir, src)
  base <- basename(src)
  dest <- file.path(out_dir, manifest$subdir[i], base)
  if (!file.exists(src_path)) { missing <- c(missing, src); next }
  file.copy(src_path, dest, overwrite = TRUE)
  cols <- read_header(dest)
  man_rows[[length(man_rows) + 1]] <- data.frame(
    file = file.path(manifest$subdir[i], base),
    md5 = unname(tools::md5sum(dest)),
    n_rows = count_rows(dest),
    n_cols = length(cols),
    supports = manifest$supports[i],
    stringsAsFactors = FALSE)
  for (cc in cols) {
    d <- describe_col(cc)
    dict_rows[[length(dict_rows) + 1]] <- data.frame(
      file = file.path(manifest$subdir[i], base), column = cc,
      unit = unname(d["unit"]), description = unname(d["description"]),
      stringsAsFactors = FALSE)
  }
}

## ---- half-hourly METEOROLOGY-ONLY extract (Ta + wind) ----------------------
# The only half-hourly file in the deposit: air temperature and wind, no fluxes
# and no calibration. Needed to reproduce Fig S5 (wind rose), which cannot be
# rebuilt from aggregates. It also carries the friagem flag: the friagem Ta
# figure was dropped from the supplement in the 2026-09 reduction, but the
# classification is retained here so the analysis stays reproducible. Derived
# from data/dataset_from_matlab.csv by column subsetting.
raw_fp <- file.path(data_dir, "dataset_from_matlab.csv")
if (!file.exists(raw_fp)) {
  missing <- c(missing, "dataset_from_matlab.csv (source for the half-hourly Ta+wind extract)")
} else {
  hh_subdir <- "halfhourly_meteorology"
  dir.create(file.path(out_dir, hh_subdir), showWarnings = FALSE)
  raw  <- read.csv(raw_fp, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- c("tv_dt", "TA_1_1_1", "WS_1_1_1", "WD_1_1_1")
  if (!all(keep %in% names(raw)))
    stop("dataset_from_matlab.csv is missing expected columns: ",
         paste(setdiff(keep, names(raw)), collapse = ", "))
  ex <- raw[, keep]
  # parse "01-Jan-2017 00:30:00" locale-independently -> ISO local time
  old_lc <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
  dt <- as.POSIXct(ex$tv_dt, format = "%d-%b-%Y %H:%M:%S", tz = "America/Lima")
  Sys.setlocale("LC_TIME", old_lc)
  ex$DateTime <- format(dt, "%Y-%m-%d %H:%M:%S")
  # Truncate to the record the paper actually analyses. The MATLAB export runs a
  # container to 2026-01-01 and carries real 2025 values, but the processed
  # record ends 15 Oct 2024 (Table S2) and the 2025 fluxes are explicitly
  # not-for-citation. Shipping the extra year would contradict the deposit's own
  # "Study period: 2017-2024".
  tmax <- as.POSIXct("2024-10-15 23:30:00", tz = "America/Lima")
  n_before <- nrow(ex)
  ex <- ex[!is.na(dt) & dt <= tmax, , drop = FALSE]
  cat("half-hourly extract truncated at", format(tmax), ":", n_before, "->", nrow(ex), "rows\n")
  ex <- ex[, c("DateTime", "TA_1_1_1", "WS_1_1_1", "WD_1_1_1")]
  for (cc in c("TA_1_1_1", "WS_1_1_1", "WD_1_1_1")) {
    ex[[cc]] <- suppressWarnings(as.numeric(ex[[cc]]))
    ex[[cc]][is.nan(ex[[cc]])] <- NA
  }
  hh_dest <- file.path(out_dir, hh_subdir, "tambopata_halfhourly_Ta_wind.csv")
  write.csv(ex, hh_dest, row.names = FALSE, na = "NA")
  man_rows[[length(man_rows) + 1]] <- data.frame(
    file = file.path(hh_subdir, basename(hh_dest)),
    md5 = unname(tools::md5sum(hh_dest)),
    n_rows = count_rows(hh_dest), n_cols = ncol(ex),
    supports = "Fig S5 (wind rose); half-hourly meteorology only (2017-01-01 to 2024-10-15 local) - no fluxes and no calibration files",
    stringsAsFactors = FALSE)
  for (cc in names(ex)) {
    d <- describe_col(cc)
    dict_rows[[length(dict_rows) + 1]] <- data.frame(
      file = file.path(hh_subdir, basename(hh_dest)), column = cc,
      unit = unname(d["unit"]), description = unname(d["description"]),
      stringsAsFactors = FALSE)
  }
  cat("half-hourly Ta+wind extract:", nrow(ex), "rows\n")
}

## ---- hard stops: missing includes or forbidden files present ---------------
if (length(missing)) {
  stop("Missing source files (not found under data/ or output/):\n  - ",
       paste(missing, collapse = "\n  - "))
}
staged <- list.files(out_dir, recursive = TRUE)
bad <- staged[Reduce(`|`, lapply(exclude_patterns, function(p) grepl(p, staged, ignore.case = TRUE)))]
if (length(bad)) {
  stop("EXCLUDED file(s) leaked into the deposit:\n  - ", paste(bad, collapse = "\n  - "))
}

## ---- write MANIFEST.csv and data_dictionary.csv ----------------------------
man <- do.call(rbind, man_rows)
dct <- do.call(rbind, dict_rows)
write.csv(man, file.path(out_dir, "MANIFEST.csv"), row.names = FALSE)
write.csv(dct, file.path(out_dir, "data_dictionary.csv"), row.names = FALSE)

## ---- LICENSE (CC BY 4.0 notice) --------------------------------------------
write_utf8(c(
  "Tambopata (PE-TNR) processed data deposit",
  "Copyright (c) 2026 Rafael Stern and the tambopata-co2 authors.",
  "",
  "This dataset is licensed under the Creative Commons Attribution 4.0",
  "International License (CC BY 4.0).",
  "You are free to share and adapt the material for any purpose, provided you",
  "give appropriate credit (cite the associated paper and this deposit).",
  "",
  "Full license text: https://creativecommons.org/licenses/by/4.0/legalcode",
  "Human-readable summary: https://creativecommons.org/licenses/by/4.0/"
), file.path(out_dir, "LICENSE.txt"))

## ---- README.md -------------------------------------------------------------
git_commit <- tryCatch(system2("git", c("-C", shQuote(repo_root), "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
readme <- c(
  "# Tambopata (PE-TNR) processed data underlying Stern et al.",
  "",
  "Processed, **aggregated** datasets underlying the figures and tables of:",
  "",
  "> Stern, R., Cruz, R., Negron-Juarez, R., Ramos, E., Knox, S. H., Hoyt, A. M.,",
  "> Jackson, R. B., Li, F., Dayalu, A., Pernak, R., Grace, J., Gloor, E.,",
  "> Salinas, N., & Cosio, E. G. *Long-Term Flux Measurements in Southwestern",
  "> Amazonia Reveal a Near-Neutral CO2 Balance and Unique Effects of El",
  "> Nino-Southern Oscillation.* Global Change Biology (submitted).",
  "",
  "**Site:** Tambopata National Reserve, Madre de Dios, Peru (AmeriFlux **PE-TNR**;",
  "AndesFlux/PUCP). **Study period:** 2017-2024.",
  "",
  "## Scope of this deposit",
  "",
  "This record contains the **aggregated processed products** that the paper's",
  "figures and tables are built from: monthly 48-point diel composites,",
  "light-response parameters, and daily / 16-day / monthly carbon-balance and",
  "meteorology summaries. The only half-hourly file is a **meteorology-only**",
  "extract (air temperature and wind, in `halfhourly_meteorology/`) needed to",
  "reproduce the wind-rose and half-hourly-temperature supplementary figures; it",
  "does **not** include the half-hourly eddy-covariance flux series or the sensor",
  "calibration files. The **complete half-hourly** eddy-covariance and meteorological record is",
  "released separately through the AmeriFlux network (site PE-TNR):",
  "https://ameriflux.lbl.gov/sites/siteinfo/PE-TNR",
  "",
  "## Files",
  "",
  "See `MANIFEST.csv` (file, md5, row/column counts, and which figure/table each",
  "supports) and `data_dictionary.csv` (per-column units and definitions).",
  "",
  paste0("Contents (", nrow(man), " data files):"),
  "")
for (sd in unique(dirname(man$file))) {
  readme <- c(readme, paste0("- **", sd, "/**"))
  sub <- man[startsWith(man$file, paste0(sd, "/")), ]
  for (j in seq_len(nrow(sub)))
    readme <- c(readme, paste0("  - `", basename(sub$file[j]), "` - ", sub$supports[j]))
}
readme <- c(readme, "",
  "## Key conventions (see data_dictionary.csv for all columns)",
  "",
  "- **NEE** sign: positive = net CO2 release to the atmosphere (source),",
  "  negative = net uptake (sink). Measured, u*-filtered NEE is the primary product.",
  "- **GEP** (gross ecosystem photosynthesis) is equivalent to GPP.",
  "- **SWC** is integrated 0-100 cm soil-water **storage in centimetres**, not a",
  "  volumetric fraction.",
  "- Half-hourly fluxes are expressed in umol CO2 m-2 s-1; daily/monthly",
  "  carbon-balance sums in g C m-2 d-1.",
  "- `hour` in the diel-composite files is the half-hour-of-day bin (48 per day),",
  "  in local time.",
  "",
  "## Provenance",
  "",
  "Generated by `code/zz_build_zenodo_deposit.R` in the analysis repository",
  "https://github.com/rafastern/tambopata-co2-gcb-2026",
  if (length(git_commit) && !is.na(git_commit[1])) paste0("(commit ", git_commit[1], ").") else "",
  "",
  "## License",
  "",
  "CC BY 4.0 - see `LICENSE.txt`.")
write_utf8(readme, file.path(out_dir, "README.md"))

## ---- summary ---------------------------------------------------------------
cat("\nDeposit assembled at:", out_dir, "\n")
cat("  data files :", nrow(man), "\n")
cat("  columns    :", nrow(dct), "(", sum(dct$unit == "" & dct$description == ""), "need manual description )\n")
cat("  total rows :", sum(man$n_rows), "\n")
cat("\nNext: fill any blank rows in data_dictionary.csv, then upload the contents\n")
cat("of data_release/ to Zenodo as a dataset record.\n")
