$ErrorActionPreference = "Stop"

$workspace = Split-Path -Parent $PSScriptRoot
$localPython = Join-Path $workspace "work\test-venv\Scripts\python.exe"
$python = if ($env:MOLECULE2FBX_PYTHON) { $env:MOLECULE2FBX_PYTHON } elseif (Test-Path -LiteralPath $localPython) { $localPython } else { "python" }
$orca = $env:ORCA_EXECUTABLE
if (-not $orca) { throw "Set ORCA_EXECUTABLE before starting the Freq jobs." }
$comparison = Join-Path $workspace "outputs\Bz_vs_SB_preliminary_comparison"
$statusPath = Join-Path $comparison "missing_freq_status.json"

Set-Location -LiteralPath $workspace

$jobs = @(
    @{
        molecule = "1Bz-LSD_RR"
        conformer = "conf008"
        xyz = Join-Path $workspace "outputs\1Bz-LSD_RR\conformers\conformer_008\conformer_008.xyz"
        output = Join-Path $workspace "outputs\1Bz-LSD_RR\conformers\frequency_additions\conformer_008"
    },
    @{
        molecule = "1SB-LSD_RR"
        conformer = "conf006"
        xyz = Join-Path $workspace "outputs\1SB-LSD_RR\conformers\conformer_006\conformer_006.xyz"
        output = Join-Path $workspace "outputs\1SB-LSD_RR\conformers\frequency_additions\conformer_006"
    }
)

$state = [ordered]@{
    status = "RUNNING"
    started_at_utc = [DateTime]::UtcNow.ToString("o")
    finished_at_utc = $null
    current_job = $null
    jobs = @()
}
$state | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding utf8

foreach ($job in $jobs) {
    if (Test-Path -LiteralPath $job.output) {
        $existing = @(Get-ChildItem -LiteralPath $job.output -Force -ErrorAction Stop)
        if ($existing.Count -gt 0) {
            throw "Refusing to overwrite non-empty frequency directory: $($job.output)"
        }
    }
    New-Item -ItemType Directory -Path $job.output -Force | Out-Null
    $state.current_job = "$($job.molecule) $($job.conformer)"
    $state | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding utf8
    $stdout = Join-Path $job.output "frequency_only.stdout.log"
    $stderr = Join-Path $job.output "frequency_only.stderr.log"
    $started = [DateTime]::UtcNow
    # Windows PowerShell converts native stderr records into terminating
    # errors when ErrorActionPreference is Stop. molecule2fbx intentionally
    # writes progress messages to stderr, so keep native execution nonfatal
    # and decide from LASTEXITCODE instead.
    $savedErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & $python -m molecule2fbx `
        --frequency-only $job.xyz `
        --output-dir $job.output `
        --method dft `
        --functional B3LYP `
        --basis def2-SVP `
        --charge 0 `
        --multiplicity 1 `
        --nprocs 16 `
        --maxcore 1000 `
        --quantum-timeout 14400 `
        --orca $orca 1>> $stdout 2>> $stderr
    $exitCode = $LASTEXITCODE
    $ErrorActionPreference = $savedErrorActionPreference
    $finished = [DateTime]::UtcNow
    $state.jobs += [ordered]@{
        molecule = $job.molecule
        conformer = $job.conformer
        exit_code = $exitCode
        started_at_utc = $started.ToString("o")
        finished_at_utc = $finished.ToString("o")
        elapsed_minutes = [Math]::Round(($finished - $started).TotalMinutes, 2)
        output_directory = $job.output
    }
    $state | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding utf8
}

$state.current_job = $null
$state.finished_at_utc = [DateTime]::UtcNow.ToString("o")
$state.status = if (@($state.jobs | Where-Object { $_.exit_code -ne 0 }).Count -eq 0) { "SUCCESS" } else { "PARTIAL" }
$state | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding utf8
