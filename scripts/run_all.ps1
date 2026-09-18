$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "_common.ps1")

# Order matters. Notable dependencies:
#   01  -> data/dataset_from_matlab_with_ENSO_season_NEEok_Reco_GEP.csv (the hub)
#   01b -> data/tambopata_48points_per_month_measured.csv, the input to Figs 4, 5,
#          6 and S14-S16, S18. Must run after 01 and before 05/09/11/26/27.
#   04  -> rewrites the hub CSV in place, adding the Gs columns, and produces
#          Fig 2a plus the 48-point composites and the 16-day tables.
#   24  -> Fig2b_cumulative_diel_NEE_panel.png, which 15 assembles into Figure 2.
#          It MUST run before 15, or Figure 2 is built from a stale/missing panel.
#   14  -> Fig2c_env_stack_panel.png, also consumed by 15.
$scripts = @(
    "code\01_prepare_flux_timeseries.R",
    "code\01b_build_measured_composite.R",
    "code\02_resample_anova.R",
    "code\03_meteorology_summary.R",
    "code\07_light_response_parameter_highlights.R",
    "code\08_light_response_residual_highlights.R",
    "code\17_influence_sensitivity_diagnostics.R",
    "code\04_diel_cycles_carbon_balance.R",
    "code\24_cumulative_diel_nee.R",
    "code\14_env_enso_seasonal_stack.R",
    "code\15_combine_figure2_env_stack_and_nee_diel.R",
    "code\05_light_response_curves.R",
    "code\06_lmm_light_response_params.R",
    "code\09_temperature_vpd_swc_bins.R",
    "code\10_surface_conductance.R",
    "code\11_nee_vs_soil_moisture_temperature.R",
    "code\12_wind_rose.R",
    "code\13_maps.R",
    "code\18_gep_environment_supplement.R",
    "code\19_data_coverage_energy_balance.R",
    "code\20_storage_diel_comparison.R",
    "code\21_light_response_halfhourly_sensitivity.R",
    "code\22_morning_transition_test.R",
    "code\25_gep_reco_self_correlation.R",
    "code\26_gross_vs_respiration_decomposition.R",
    "code\27_vpd_class_light_response.R"
)

# Needs the raw Biomet.net database (set PETNR_DB_ROOT) for LW_IN/LW_OUT, which
# is not part of the repository. Warn rather than stop when it is absent.
$optional_scripts = @(
    "code\23_canopy_temperature.R"
)

foreach ($script in $scripts) {
    Write-Host ""
    Write-Host "==> $script"

    Invoke-RepoRscript -ScriptPath $script
    $exitCode = $script:LastRscriptExitCode
    if ($exitCode -ne 0) {
        Write-Error "Stopped after failure in ${script} with exit code ${exitCode}."
        exit $exitCode
    }
}

foreach ($script in $optional_scripts) {
    Write-Host ""
    Write-Host "==> $script (optional)"

    Invoke-RepoRscript -ScriptPath $script
    $exitCode = $script:LastRscriptExitCode
    if ($exitCode -ne 0) {
        Write-Warning "Skipped ${script} (exit code ${exitCode}); it needs the raw Biomet.net database."
    }
}

Write-Host ""
Write-Host "All analysis scripts completed."
Write-Host "Build the Zenodo data deposit with: code\zz_build_zenodo_deposit.R"
