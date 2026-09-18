function Get-RepoRoot {
    param(
        [string]$StartPath = $PSScriptRoot
    )

    $current = (Resolve-Path -LiteralPath $StartPath).Path

    while ($true) {
        if (
            (Test-Path -LiteralPath (Join-Path $current "code\paths.R")) -and
            (Test-Path -LiteralPath (Join-Path $current "data"))
        ) {
            return $current
        }

        $parent = Split-Path -Parent $current
        if ($parent -eq $current -or [string]::IsNullOrWhiteSpace($parent)) {
            throw "Could not find repository root containing code\paths.R and data\ from: $StartPath"
        }

        $current = $parent
    }
}

function Get-RscriptPath {
    $onPath = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($onPath) {
        return $onPath.Source
    }

    $fallback = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
    if (Test-Path -LiteralPath $fallback) {
        return $fallback
    }

    throw "Could not find Rscript. Add Rscript to PATH or install R at: $fallback"
}

function Invoke-RepoRscript {
    param(
        [Parameter(Mandatory = $true)]
        [string]$ScriptPath,

        [string[]]$ScriptArgs = @()
    )

    $repoRoot = Get-RepoRoot
    $rscript = Get-RscriptPath
    $fullScriptPath = Join-Path $repoRoot $ScriptPath

    if (-not (Test-Path -LiteralPath $fullScriptPath)) {
        throw "R script not found: $ScriptPath"
    }

    Push-Location $repoRoot
    try {
        $repoLib = Join-Path $repoRoot ".r-lib"
        if (Test-Path -LiteralPath $repoLib) {
            $env:R_LIBS_USER = $repoLib
        }

        & $rscript $ScriptPath @ScriptArgs
        $script:LastRscriptExitCode = $LASTEXITCODE
    }
    finally {
        Pop-Location
    }
}
