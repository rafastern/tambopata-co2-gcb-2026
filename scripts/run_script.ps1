param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$ScriptPath,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ScriptArgs
)

$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "_common.ps1")

Invoke-RepoRscript -ScriptPath $ScriptPath -ScriptArgs $ScriptArgs
exit $script:LastRscriptExitCode
