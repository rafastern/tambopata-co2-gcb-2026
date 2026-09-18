#!/usr/bin/env Rscript
# =============================================================================
# run_all.R -- cross-platform equivalent of scripts/run_all.ps1
#
# Runs the full analysis in dependency order. Each script is run in a fresh
# child process so that one script's globals cannot leak into the next, which
# matches how the PowerShell runner behaves.
#
#   Rscript scripts/run_all.R
#
# Order notes (see README):
#   01  -> the enriched half-hourly hub CSV
#   01b -> the measured-only composite behind Figs 4, 5, 6 and S14-S16, S18
#   04  -> rewrites the hub CSV in place, adding the Gs columns
#   24  -> Figure 2 panel b; MUST run before 15, which assembles Figure 2
# =============================================================================

# locate the repo root from --file= when run via Rscript, else from getwd()
find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  cand <- if (length(f)) dirname(dirname(normalizePath(f[1]))) else normalizePath(getwd())
  # walk up until we find a directory holding both code/ and data/
  repeat {
    if (dir.exists(file.path(cand, "code")) && dir.exists(file.path(cand, "data")))
      return(cand)
    parent <- dirname(cand)
    if (identical(parent, cand))
      stop("could not find the repository root (needs both code/ and data/). ",
           "Run as: Rscript scripts/run_all.R")
    cand <- parent
  }
}
repo_root <- find_root()

scripts <- c(
  "01_prepare_flux_timeseries.R",
  "01b_build_measured_composite.R",
  "02_resample_anova.R",
  "03_meteorology_summary.R",
  "07_light_response_parameter_highlights.R",
  "08_light_response_residual_highlights.R",
  "17_influence_sensitivity_diagnostics.R",
  "04_diel_cycles_carbon_balance.R",
  "24_cumulative_diel_nee.R",
  "14_env_enso_seasonal_stack.R",
  "15_combine_figure2_env_stack_and_nee_diel.R",
  "05_light_response_curves.R",
  "06_lmm_light_response_params.R",
  "09_temperature_vpd_swc_bins.R",
  "10_surface_conductance.R",
  "11_nee_vs_soil_moisture_temperature.R",
  "12_wind_rose.R",
  "13_maps.R",
  "18_gep_environment_supplement.R",
  "19_data_coverage_energy_balance.R",
  "20_storage_diel_comparison.R",
  "21_light_response_halfhourly_sensitivity.R",
  "22_morning_transition_test.R",
  "25_gep_reco_self_correlation.R",
  "26_gross_vs_respiration_decomposition.R",
  "27_vpd_class_light_response.R"
)

# needs the raw Biomet.net database (set PETNR_DB_ROOT); warn instead of failing
optional_scripts <- c("23_canopy_temperature.R")

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

run_one <- function(name, optional = FALSE) {
  path <- file.path(repo_root, "code", name)
  if (!file.exists(path)) stop("script not found: ", path)
  message("\n==> code/", name, if (optional) " (optional)" else "")
  status <- system2(rscript, shQuote(path), stdout = "", stderr = "")
  if (status != 0) {
    if (optional) {
      warning("skipped code/", name, " (exit ", status,
              "); it needs the raw Biomet.net database.", call. = FALSE)
    } else {
      stop("stopped after failure in code/", name, " with exit code ", status)
    }
  }
  invisible(status)
}

setwd(repo_root)
for (s in scripts) run_one(s)
for (s in optional_scripts) run_one(s, optional = TRUE)

message("\nAll analysis scripts completed.")
message("Build the Zenodo data deposit with: Rscript code/zz_build_zenodo_deposit.R")
