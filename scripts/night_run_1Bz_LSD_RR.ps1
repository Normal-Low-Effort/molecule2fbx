$ErrorActionPreference = "Stop"

$workspace = Split-Path -Parent $PSScriptRoot
$localPython = Join-Path $workspace "work\test-venv\Scripts\python.exe"
$python = if ($env:MOLECULE2FBX_PYTHON) { $env:MOLECULE2FBX_PYTHON } elseif (Test-Path -LiteralPath $localPython) { $localPython } else { "python" }
$orca = $env:ORCA_EXECUTABLE
$blender = $env:BLENDER_EXECUTABLE
if (-not $orca) { throw "Set ORCA_EXECUTABLE before starting the night run." }
if (-not $blender) { throw "Set BLENDER_EXECUTABLE before starting the night run." }
$output = Join-Path $workspace "outputs\1Bz-LSD_RR"
$stdoutLog = Join-Path $output "night_run.stdout.log"
$stderrLog = Join-Path $output "night_run.stderr.log"
$statusPath = Join-Path $output "night_run_status.json"
$smiles = "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccccc4)c5cccc(C2=C1)c35"
$optionalReuse = Join-Path $workspace "outputs\1Bz-LSD_RR_redo\1Bz-LSD_RR_dft_calculations"

New-Item -ItemType Directory -Path $output -Force | Out-Null
Set-Location -LiteralPath $workspace

$arguments = @(
    "-m", "molecule2fbx",
    "--ensemble",
    "--smiles", $smiles,
    "--name", "1Bz-LSD_RR",
    "--method", "dft",
    "--functional", "B3LYP",
    "--basis", "def2-SVP",
    "--charge", "0",
    "--multiplicity", "1",
    "--conformer-pool", "200",
    "--conformers", "10",
    "--forcefield-energy-window-kj", "10",
    "--conformer-rmsd-threshold", "0.75",
    "--dft-rmsd-threshold", "0.75",
    "--frequency",
    "--frequency-window-kj", "5",
    "--frequency-max", "3",
    "--nprocs", "16",
    "--maxcore", "1000",
    "--output-dir", $output,
    "--orca", $orca,
    "--blender", $blender,
    "--expected-stereocenters", "2",
    "--stereochemistry-label", "(6aR,9R)",
    "--quantum-timeout", "14400"
)

$reuseSource = $null
if (Test-Path -LiteralPath $optionalReuse -PathType Container) {
    $arguments += @("--reuse-calculations", $optionalReuse)
    $reuseSource = $optionalReuse
}

$startedAt = [DateTime]::UtcNow.ToString("o")
@{
    status = "RUNNING"
    started_at_utc = $startedAt
    process_id = $PID
    command = "$python " + ($arguments -join " ")
    electronic_structure = "ORCA 6.1.1 B3LYP/def2-SVP charge=0 multiplicity=1"
    resources = @{ nprocs = 16; maxcore_mb_per_process = 1000 }
    source_reuse_directory = $reuseSource
    new_calculation_directory = (Join-Path $output "conformers")
} | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusPath -Encoding UTF8

Add-Content -LiteralPath $stdoutLog -Encoding UTF8 -Value "[$startedAt] Night run started"
$exitCode = 1
$failure = $null
try {
    $ErrorActionPreference = "Continue"
    & $python @arguments 1>> $stdoutLog 2>> $stderrLog
    $exitCode = $LASTEXITCODE
    $ErrorActionPreference = "Stop"
}
catch {
    $failure = $_.Exception.Message
    Add-Content -LiteralPath $stderrLog -Encoding UTF8 -Value $failure
}
finally {
    $finishedAt = [DateTime]::UtcNow.ToString("o")
    $reportedStatus = $(if ($exitCode -eq 0) { "SUCCESS" } else { "FAILED" })
    $ensemblePath = Join-Path $output "ensemble.json"
    if (Test-Path -LiteralPath $ensemblePath) {
        try {
            $ensembleStatus = (Get-Content -LiteralPath $ensemblePath -Raw -Encoding UTF8 | ConvertFrom-Json).calculation_status
            if ($ensembleStatus) { $reportedStatus = $ensembleStatus }
        }
        catch {
            Add-Content -LiteralPath $stderrLog -Encoding UTF8 -Value "Could not read final ensemble status: $($_.Exception.Message)"
        }
    }
    @{
        status = $reportedStatus
        exit_code = $exitCode
        started_at_utc = $startedAt
        finished_at_utc = $finishedAt
        process_id = $PID
        failure = $failure
        stdout_log = $stdoutLog
        stderr_log = $stderrLog
    } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusPath -Encoding UTF8
}

exit $exitCode
