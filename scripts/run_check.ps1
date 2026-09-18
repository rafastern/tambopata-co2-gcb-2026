$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "_common.ps1")

Invoke-RepoRscript -ScriptPath "code\check_data.R"
exit $script:LastRscriptExitCode
