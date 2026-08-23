param(
    [string]$OrcaExecutable = $(
        if ($env:ORCA_EXECUTABLE) { $env:ORCA_EXECUTABLE }
        else { "C:\ORCA_6.1.1\orca.exe" }
    ),
    [string]$BlenderExecutable = $(
        if ($env:BLENDER_EXECUTABLE) { $env:BLENDER_EXECUTABLE }
        else { "C:\Program Files\Blender Foundation\Blender 4.0\blender.exe" }
    ),
    [int]$NProcs = 16,
    [int]$MaxCore = 1000,
    [int]$QuantumTimeoutSeconds = 14400
)

$ErrorActionPreference = "Stop"

$workspace = Split-Path -Parent $PSScriptRoot
$python = Join-Path $workspace "work\test-venv\Scripts\python.exe"
$outputsRoot = Join-Path $workspace "outputs"
$queueStatusPath = Join-Path $outputsRoot "para_alkyl_controls_queue_status.json"
$queueLogPath = Join-Path $outputsRoot "para_alkyl_controls_queue.log"

foreach ($required in @($python, $OrcaExecutable, $BlenderExecutable)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Required executable was not found: $required"
    }
}
if ($NProcs -lt 1) { throw "NProcs must be at least 1" }
if ($MaxCore -lt 1) { throw "MaxCore must be at least 1" }
if ($QuantumTimeoutSeconds -le 0) { throw "QuantumTimeoutSeconds must be positive" }

$controls = @(
    [pscustomobject]@{
        Name = "1pMeBz-LSD_RR"
        Label = "1-(4-methylbenzoyl)-LSD (6aR,9R)"
        Smiles = "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccc(cc4)C)c5cccc(C2=C1)c35"
    },
    [pscustomobject]@{
        Name = "1p-iPrBz-LSD_RR"
        Label = "1-[4-(propan-2-yl)benzoyl]-LSD (6aR,9R)"
        Smiles = "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccc(cc4)C(C)C)c5cccc(C2=C1)c35"
    },
    [pscustomobject]@{
        Name = "1ptBuBz-LSD_RR"
        Label = "1-(4-tert-butylbenzoyl)-LSD (6aR,9R)"
        Smiles = "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccc(cc4)C(C)(C)C)c5cccc(C2=C1)c35"
    }
)

$startedAt = [DateTime]::UtcNow.ToString("o")
$records = foreach ($control in $controls) {
    [ordered]@{
        name = $control.Name
        label = $control.Label
        smiles = $control.Smiles
        status = "PENDING"
        started_at_utc = $null
        finished_at_utc = $null
        exit_code = $null
        output_directory = (Join-Path $outputsRoot $control.Name)
        ensemble_json = (Join-Path $outputsRoot "$($control.Name)\ensemble.json")
    }
}

function Write-QueueStatus([string]$Status) {
    $payload = [ordered]@{
        status = $Status
        run_type = "para_alkyl_control_ensemble_dft_screening"
        started_at_utc = $startedAt
        updated_at_utc = [DateTime]::UtcNow.ToString("o")
        process_id = $PID
        execution = "sequential"
        existing_reference_ensembles = @(
            (Join-Path $outputsRoot "1Bz-LSD_RR"),
            (Join-Path $outputsRoot "1SB-LSD_RR")
        )
        electronic_structure = "ORCA 6.1.1 B3LYP/def2-SVP gas phase, charge=0, multiplicity=1"
        ensemble_settings = [ordered]@{
            conformer_pool = 200
            forcefield_energy_window_kj_mol = 10.0
            pre_dft_heavy_atom_rmsd_angstrom = 0.75
            maximum_dft_conformers = 10
            post_dft_heavy_atom_rmsd_angstrom = 0.75
            frequency_window_kj_mol = 5.0
            frequency_maximum = 3
            random_seed = 61453
        }
        resources = [ordered]@{
            nprocs = $NProcs
            maxcore_mb_per_process = $MaxCore
            quantum_timeout_seconds_per_job = $QuantumTimeoutSeconds
        }
        controls = $records
    }
    $payload | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $queueStatusPath -Encoding UTF8
}

New-Item -ItemType Directory -Path $outputsRoot -Force | Out-Null
Set-Location -LiteralPath $workspace
Write-QueueStatus "RUNNING"
Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "[$startedAt] Para-alkyl control queue started (PID $PID)"

for ($position = 0; $position -lt $controls.Count; $position++) {
    $control = $controls[$position]
    $record = $records[$position]
    $output = [string]$record.output_directory
    $ensemblePath = [string]$record.ensemble_json
    $stdoutLog = Join-Path $output "night_run.stdout.log"
    $stderrLog = Join-Path $output "night_run.stderr.log"
    $statusPath = Join-Path $output "night_run_status.json"

    if (Test-Path -LiteralPath $ensemblePath -PathType Leaf) {
        try {
            $existingStatus = (Get-Content -LiteralPath $ensemblePath -Raw -Encoding UTF8 | ConvertFrom-Json).calculation_status
            if ($existingStatus -eq "SUCCESS") {
                $record.status = "REUSED_COMPLETED"
                $record.finished_at_utc = [DateTime]::UtcNow.ToString("o")
                Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "[$($record.finished_at_utc)] $($control.Name) already complete; skipped"
                Write-QueueStatus "RUNNING"
                continue
            }
        }
        catch {
            Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "Could not parse existing ensemble for $($control.Name): $($_.Exception.Message)"
        }
    }

    New-Item -ItemType Directory -Path $output -Force | Out-Null
    $arguments = @(
        "-m", "molecule2fbx",
        "--ensemble",
        "--smiles", $control.Smiles,
        "--name", $control.Name,
        "--method", "dft",
        "--functional", "B3LYP",
        "--basis", "def2-SVP",
        "--charge", "0",
        "--multiplicity", "1",
        "--conformer-pool", "200",
        "--conformers", "10",
        "--random-seed", "61453",
        "--forcefield-energy-window-kj", "10",
        "--conformer-rmsd-threshold", "0.75",
        "--dft-rmsd-threshold", "0.75",
        "--frequency",
        "--frequency-window-kj", "5",
        "--frequency-max", "3",
        "--nprocs", "$NProcs",
        "--maxcore", "$MaxCore",
        "--output-dir", $output,
        "--orca", $OrcaExecutable,
        "--blender", $BlenderExecutable,
        "--expected-stereocenters", "2",
        "--stereochemistry-label", "(6aR,9R)",
        "--quantum-timeout", "$QuantumTimeoutSeconds"
    )

    $record.status = "RUNNING"
    $record.started_at_utc = [DateTime]::UtcNow.ToString("o")
    Write-QueueStatus "RUNNING"
    @{
        status = "RUNNING"
        queue_process_id = $PID
        queue_position = $position + 1
        queue_total = $controls.Count
        started_at_utc = $record.started_at_utc
        command = "$python " + ($arguments -join " ")
        output_directory = $output
    } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusPath -Encoding UTF8
    Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "[$($record.started_at_utc)] Starting $($control.Name) ($($position + 1)/$($controls.Count))"

    $exitCode = 1
    try {
        $ErrorActionPreference = "Continue"
        & $python @arguments 1>> $stdoutLog 2>> $stderrLog
        $exitCode = $LASTEXITCODE
    }
    catch {
        Add-Content -LiteralPath $stderrLog -Encoding UTF8 -Value $_.Exception.Message
    }
    finally {
        $ErrorActionPreference = "Stop"
    }

    $record.exit_code = $exitCode
    $record.finished_at_utc = [DateTime]::UtcNow.ToString("o")
    $record.status = $(if ($exitCode -eq 0) { "SUCCESS" } else { "FAILED" })
    if (Test-Path -LiteralPath $ensemblePath -PathType Leaf) {
        try {
            $reported = (Get-Content -LiteralPath $ensemblePath -Raw -Encoding UTF8 | ConvertFrom-Json).calculation_status
            if ($reported) { $record.status = $reported }
        }
        catch {
            Add-Content -LiteralPath $stderrLog -Encoding UTF8 -Value "Could not read final ensemble status: $($_.Exception.Message)"
        }
    }
    @{
        status = $record.status
        exit_code = $exitCode
        queue_process_id = $PID
        started_at_utc = $record.started_at_utc
        finished_at_utc = $record.finished_at_utc
        stdout_log = $stdoutLog
        stderr_log = $stderrLog
        ensemble_json = $ensemblePath
    } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusPath -Encoding UTF8
    Write-QueueStatus "RUNNING"
    Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "[$($record.finished_at_utc)] Finished $($control.Name): $($record.status) (exit $exitCode)"
}

$successful = @($records | Where-Object { $_.status -in @("SUCCESS", "REUSED_COMPLETED") }).Count
$finalStatus = if ($successful -eq $records.Count) { "SUCCESS" } elseif ($successful -gt 0) { "PARTIAL" } else { "FAILED" }
Write-QueueStatus $finalStatus
Add-Content -LiteralPath $queueLogPath -Encoding UTF8 -Value "[$([DateTime]::UtcNow.ToString('o'))] Queue finished: $finalStatus"
exit $(if ($finalStatus -eq "SUCCESS") { 0 } else { 1 })
