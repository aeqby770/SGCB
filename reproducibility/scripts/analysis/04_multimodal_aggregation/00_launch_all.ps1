# =============================================================================
#       [-Tier <sgcb|other|all>] [-MaxJobs 12] [-MaxMemGB 46] [-ReserveGB 10]
#       [-SkipExisting]
# =============================================================================

param(
    [string]$Category    = "all",
    [string]$Tier        = "all",      # sgcb | other | all
    [int]$MaxJobs        = 12,
    [double]$MaxMemGB    = 46.0,
    [double]$ReserveGB   = 10.0,
    [double]$PressureGB  = 8.0,
    [double]$CriticalGB  = 5.5,
    [switch]$SerialTail   = $false,
    [switch]$SkipExisting = $true,
    [string[]]$IncludeMethods = @()
)

$RSCRIPT = if ($env:RSCRIPT) { $env:RSCRIPT } elseif (Get-Command Rscript -ErrorAction SilentlyContinue) { (Get-Command Rscript).Source } else { "Rscript" }
$BASE    = $PSScriptRoot   # Runtime-resolved: avoids PS 5.1 encoding issues with Chinese paths
$LOG_DIR = "$BASE\logs"
New-Item -ItemType Directory -Force -Path $LOG_DIR | Out-Null

$CPU_THREADS = [Environment]::ProcessorCount
$script:ThreadCount = if ($SerialTail) { $CPU_THREADS } else { 1 }
if ($SerialTail) {
    $MaxJobs = 1
    $MaxMemGB = [double]::MaxValue
    $ReserveGB = 0.0
    $PressureGB = 0.0
    $CriticalGB = 0.0
}

function Get-ThreadCountForJob([string]$Name) {
    if ($SerialTail) { return $CPU_THREADS }
    return 1
}

$HISTORY_CSV = [IO.Path]::Combine($BASE, "job_history.csv")

$script:histProfile = @{}
if (Test-Path $HISTORY_CSV) {
    $rows = Import-Csv $HISTORY_CSV -Encoding UTF8
    $grouped = @{}
    foreach ($row in $rows) {
        if ($row.Success -ne 'True') { continue }
        $sig = $row.Signature
        if (-not $grouped.ContainsKey($sig)) { $grouped[$sig] = [System.Collections.Generic.List[double]]::new() }
        $v = 0.0; if ([double]::TryParse($row.PeakPrivMB, [ref]$v) -and $v -gt 0) { $grouped[$sig].Add($v) }
    }
    foreach ($kv in $grouped.GetEnumerator()) {
        $vals = $kv.Value | Sort-Object
        $idx  = [math]::Min($vals.Count - 1, [math]::Floor($vals.Count * 0.9))
        $script:histProfile[$kv.Key] = [math]::Ceiling($vals[$idx] * 1.15)
    }
    Write-Host "  [HIST] Loaded $($script:histProfile.Count) task signatures from job_history.csv"
}

function Get-TaskSignature([string]$Name) {
    $n = $Name -replace '_null\d+_', '_null_'
    $n = $n -replace '_null\d+$', '_null'
    return $n
}

function Get-DefaultMemMB([string]$Name) {
    if ($Name -match '^prot_')                                     { return 1000 }
    if ($Name -match '^bulk_other_null')                           { return 1500 }
    if ($Name -match '^bulk_other_dv')                             { return 4000 }
    if ($Name -match '^bulk_sgcb_dv')                              { return 4000 }
    if ($Name -match '^bulk_sgcb_gtex')                            { return 3000 }
    if ($Name -match '^bulk_sgcb_tcga')                            { return 3500 }
    if ($Name -match '^bulk_sgcb_(pb|null)')                       { return 2500 }
    if ($Name -match '^bulk_sgcb')                                 { return 1500 }
    if ($Name -match 'NEBULA')                                     { return 7000 }
    if ($Name -match 'scrna_other.*BPSC_xin')                      { return 12000 }
    if ($Name -match 'scrna.*zilionis')                            { return 5000 }
    if ($Name -match 'scrna.*muscat_sim')                          { return 2500 }
    if ($Name -match 'scrna_other.*(MAST|scDD|DEsingle|BPSC)')    { return 4500 }
    if ($Name -match 'scrna.*kang18')                              { return 4000 }
    if ($Name -match '^scrna_sgcb')                                { return 3500 }
    if ($Name -match '^scrna_other')                               { return 3000 }
    if ($Name -match '^paper_sgcb_reviewer')                        { return 4000 }
    if ($Name -match '^paper_sgcb_ablation')                        { return 3000 }
    if ($Name -match '^paper_sgcb_gg_gof')                          { return 3500 }
    return 2000
}

function Get-MemMB([string]$Name) {
    $sig = Get-TaskSignature $Name
    if ($script:histProfile.ContainsKey($sig)) {
        $hist = $script:histProfile[$sig]
        $def  = Get-DefaultMemMB $Name
        return [math]::Max($def, $hist)
    }
    return Get-DefaultMemMB $Name
}

$script:_memCounter = [System.Diagnostics.PerformanceCounter]::new('Memory', 'Available MBytes')
function Get-FreeMemGB {
    [math]::Round($script:_memCounter.NextValue() / 1024, 2)
}

function Get-MethodCap([string]$Name) {
    if ($Name -match 'NEBULA')               { return 1 }
    if ($Name -match 'GAMLSS')               { return 1 }
    if ($Name -match 'BPSC')                 { return 1 }
    if ($Name -match 'MDSeq')                { return 2 }
    if ($Name -match '(scDD|DEsingle)')      { return 2 }
    if ($Name -match 'MAST')                 { return 2 }
    return 999
}

function Get-RunningCountLike([string]$pattern) {
    $n = 0
    foreach ($p in $running) {
        if (-not $p.Proc.HasExited -and $p.Name -match $pattern) { $n++ }
    }
    return $n
}

function Test-MethodCapOK([string]$Name) {
    if ($Name -match 'NEBULA')               { return (Get-RunningCountLike 'NEBULA') -lt 1 }
    if ($Name -match 'GAMLSS')               { return (Get-RunningCountLike 'GAMLSS') -lt 1 }
    if ($Name -match 'BPSC')                 { return (Get-RunningCountLike 'BPSC') -lt 1 }
    if ($Name -match 'MDSeq')                { return (Get-RunningCountLike 'MDSeq') -lt 2 }
    if ($Name -match '(scDD|DEsingle)')      { return (Get-RunningCountLike '(scDD|DEsingle)') -lt 2 }
    if ($Name -match 'MAST')                 { return (Get-RunningCountLike 'MAST') -lt 2 }
    return $true
}

function Get-TimeSec([string]$Name) {
    if ($Name -match '^prot_.*_ttest')                             { return 10 }
    if ($Name -match '^prot_.*_limma')                             { return 15 }
    if ($Name -match '^prot_.*_(DEqMS|ROTS|proDA)')                { return 30 }
    if ($Name -match '^prot_.*_msqrob2')                           { return 60 }
    if ($Name -match '^prot_.*_prolfqua')                          { return 60 }
    if ($Name -match '^prot_.*_MSstats')                           { return 120 }
    if ($Name -match '^prot_.*dv_diffVar')                         { return 20 }
    if ($Name -match '^prot_.*dv_GAMLSS')                          { return 600 }
    if ($Name -match '^prot_sgcb')                                 { return 120 }
    if ($Name -match '^prot_')                                     { return 60 }
    if ($Name -match '^bulk.*(limma|edgeR)')                       { return 30 }
    if ($Name -match '^bulk.*diffVar')                             { return 40 }
    if ($Name -match '^bulk.*clrDV')                               { return 15 }
    if ($Name -match '^bulk.*(glmGamPoi|samr|NOISeq)')             { return 120 }
    if ($Name -match '^bulk.*EBSeq')                               { return 300 }
    if ($Name -match '^bulk.*DESeq2')                              { return 900 }
    if ($Name -match '^bulk.*MDSeq')                               { return 3600 }
    if ($Name -match '^bulk.*GAMLSS')                              { return 7200 }
    if ($Name -match '^bulk_sgcb_dv')                              { return 1800 }
    if ($Name -match '^bulk_sgcb')                                 { return 600 }
    if ($Name -match '^bulk_')                                     { return 300 }
    if ($Name -match 'scrna.*(limma|edgeR|Seurat)')                { return 30 }
    if ($Name -match 'scrna.*glmGamPoi')                           { return 60 }
    if ($Name -match 'scrna.*aggregateBioVar')                     { return 60 }
    if ($Name -match 'scrna.*glmmTMB')                             { return 120 }
    if ($Name -match 'scrna.*DESeq2')                              { return 180 }
    if ($Name -match 'scrna.*diffVar')                             { return 60 }
    if ($Name -match 'scrna.*muscat')                              { return 120 }
    if ($Name -match 'scrna.*(MAST|distinct)')                     { return 600 }
    if ($Name -match 'scrna.*DEsingle')                            { return 1800 }
    if ($Name -match 'scrna.*(scDD|BPSC)')                         { return 3600 }
    if ($Name -match 'scrna.*NEBULA')                              { return 1200 }
    if ($Name -match 'scrna.*GAMLSS')                              { return 7200 }
    if ($Name -match 'scrna_sgcb')                                 { return 600 }
    if ($Name -match 'scrna_')                                     { return 300 }
    if ($Name -match '^paper_sgcb_reviewer')                       { return 3600 }
    if ($Name -match '^paper_sgcb_ablation')                       { return 1800 }
    if ($Name -match '^paper_sgcb_gg_gof')                         { return 600 }
    return 300
}

function Get-JobTimeout([string]$Name) {
    if ($SerialTail) { return [int]::MaxValue }
    if ($Name -match '(scDD|DEsingle|BPSC)') { return 86400 }
    return 14400
}

$jobs = [System.Collections.ArrayList]@()

# ==================== 1. Bulk SGCB ====================
if (($Category -eq "all" -or $Category -eq "bulk") -and ($Tier -eq "all" -or $Tier -eq "sgcb")) {
    $s = "$BASE\01_bulk_rnaseq\sgcb\run_bulk_sgcb.R"
    $r = "$BASE\01_bulk_rnaseq\sgcb\results"

    # Classic (6)
    @("bottomly","hammer","pasilla","airway","sultan","wang") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="classic $_"; Name="bulk_sgcb_classic_$_"; ResultDir=$r; ResultCSV="classic__$_.csv" }) | Out-Null
    }
    # Simulation (6)
    @("sim_n3_pDE10_ef1.5","sim_n3_pDE10_ef2.0","sim_n5_pDE10_ef1.5","sim_n5_pDE10_ef2.0","sim_n10_pDE10_ef1.5","sim_n10_pDE10_ef2.0") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="simulation $_"; Name="bulk_sgcb_sim_$_"; ResultDir=$r; ResultCSV="simulation__$_.csv" }) | Out-Null
    }
    # GTEx (5)
    @("liver_vs_kidney","brain_vs_heart","lung_vs_colon","adipose_tissue_vs_muscle","blood_vs_skin") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="gtex_pair $_"; Name="bulk_sgcb_gtex_$_"; ResultDir=$r; ResultCSV="gtex_pair__$_.csv" }) | Out-Null
    }
    # TCGA (7)
    @("brca","luad","kirc","lihc","coad","hnsc","stad") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="tcga $_"; Name="bulk_sgcb_tcga_$_"; ResultDir=$r; ResultCSV="tcga__$_.csv" }) | Out-Null
    }
    # Pseudobulk (4)
    @("kang18","segerstolpe","xin","zilionis") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="pseudobulk $_"; Name="bulk_sgcb_pb_$_"; ResultDir=$r; ResultCSV="pseudobulk__$_.csv" }) | Out-Null
    }
    # Null FDR (10)
    1..10 | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="null_fdr $_"; Name="bulk_sgcb_null_$_"; ResultDir=$r; ResultCSV="null_fdr__$_.csv" }) | Out-Null
    }
    # SEQC (1)
    $jobs.Add(@{ Script=$s; Args="seqc seqc"; Name="bulk_sgcb_seqc"; ResultDir=$r; ResultCSV="seqc__seqc.csv" }) | Out-Null
    # DV (3)
    @("cerebellum_vs_cortex","cerehemi_vs_cortex","hippo_vs_hypo") | ForEach-Object {
        $jobs.Add(@{ Script=$s; Args="dv $_"; Name="bulk_sgcb_dv_$_"; ResultDir=$r; ResultCSV="dv__$_.csv" }) | Out-Null
    }
}

# ==================== 1b. Bulk Other Methods ====================
if (($Category -eq "all" -or $Category -eq "bulk") -and ($Tier -eq "all" -or $Tier -eq "other")) {
    $s = "$BASE\01_bulk_rnaseq\other_methods\scripts\run_other_bulk.R"
    $r = "$BASE\01_bulk_rnaseq\other_methods\results"
    # DE methods (7) + DV methods (4) = 11
    $bulk_de_methods = @("DESeq2","edgeR","limma","glmGamPoi","EBSeq","NOISeq","samr")
    $bulk_dv_methods = @("dv_MDSeq","dv_clrDV","dv_diffVar","dv_GAMLSS")
    $bulk_methods = $bulk_de_methods + $bulk_dv_methods

    # Null FDR (11 x 10 = 110)
    foreach ($m in $bulk_methods) {
        1..10 | ForEach-Object {
            $jobs.Add(@{ Script=$s; Args="$m null_fdr $_"; Name="bulk_other_null${_}_$m"; ResultDir=$r; ResultCSV="${m}__null_fdr__$_.csv" }) | Out-Null
        }
    }
    # DV real data: GTEx brain pairs (4 DV x 3 = 12)
    $dv_datasets = @("cerebellum_vs_cortex","cerehemi_vs_cortex","hippo_vs_hypo")
    foreach ($m in $bulk_dv_methods) {
        foreach ($d in $dv_datasets) {
            $jobs.Add(@{ Script=$s; Args="$m dv $d"; Name="bulk_other_dv_${m}_$d"; ResultDir=$r; ResultCSV="${m}__dv__$d.csv" }) | Out-Null
        }
    }
}

# ==================== 2. scRNA SGCB ====================
if (($Category -eq "all" -or $Category -eq "scrna") -and ($Tier -eq "all" -or $Tier -eq "sgcb")) {
    $s = "$BASE\02_scrna_pseudobulk\sgcb\scripts\run_sgcb_scrna.R"
    $r = "$BASE\02_scrna_pseudobulk\sgcb\results"
    $datasets = @("kang18","segerstolpe","xin","zilionis","jakel","muscat_sim")

    # Normal DE (6)
    foreach ($d in $datasets) {
        $jobs.Add(@{ Script=$s; Args="$d"; Name="scrna_sgcb_$d"; ResultDir=$r; ResultCSV="SGCB__$d.csv" }) | Out-Null
    }
    # Null FDR (6 x 10 = 60)
    foreach ($d in $datasets) {
        1..10 | ForEach-Object {
            $jobs.Add(@{ Script=$s; Args="$d $_"; Name="scrna_sgcb_null${_}_$d"; ResultDir=$r; ResultCSV="null_${_}__SGCB__$d.csv" }) | Out-Null
        }
    }
}

# ==================== 3. scRNA Other Methods ====================
if (($Category -eq "all" -or $Category -eq "scrna") -and ($Tier -eq "all" -or $Tier -eq "other")) {
    $s = "$BASE\02_scrna_pseudobulk\other_methods\scripts\run_other_scrna.R"
    $r = "$BASE\02_scrna_pseudobulk\other_methods\results"
    # DE methods (14) + DV methods (2) = 16
    $de_methods = @("DESeq2","edgeR","limma","glmGamPoi","Seurat",
                    "MAST","scDD","DEsingle","distinct","BPSC",
                    "muscat","aggregateBioVar","glmmTMB","NEBULA")
    $dv_methods = @("dv_diffVar","dv_GAMLSS")
    $methods  = $de_methods + $dv_methods
    $datasets = @("kang18","segerstolpe","xin","zilionis","jakel","muscat_sim")
    # NEBULA cannot run on muscat_sim (no single-cell data)
    $nebula_datasets = @("kang18","segerstolpe","xin","zilionis","jakel")

    # Normal: non-NEBULA (15 x 6 = 90) + NEBULA (1 x 5 = 5) = 95
    foreach ($m in $methods) {
        $ds_list = switch($m) { "NEBULA" { $nebula_datasets } default { $datasets } }
        foreach ($d in $ds_list) {
            $jobs.Add(@{ Script=$s; Args="$m $d"; Name="scrna_other_${m}_$d"; ResultDir=$r; ResultCSV="${m}__$d.csv" }) | Out-Null
        }
    }
    # Null FDR: non-NEBULA (15 x 6 x 10 = 900) + NEBULA (1 x 5 x 10 = 50) = 950
    foreach ($m in $methods) {
        $ds_list = switch($m) { "NEBULA" { $nebula_datasets } default { $datasets } }
        foreach ($d in $ds_list) {
            1..10 | ForEach-Object {
                $jobs.Add(@{ Script=$s; Args="$m $d $_"; Name="scrna_other_null${_}_${m}_$d"; ResultDir=$r; ResultCSV="null_${_}__${m}__$d.csv" }) | Out-Null
            }
        }
    }
}

# ==================== 4. Proteomics SGCB ====================
if (($Category -eq "all" -or $Category -eq "proteomics") -and ($Tier -eq "all" -or $Tier -eq "sgcb")) {
    $s = "$BASE\03_proteomics\sgcb\scripts\run_sgcb_prot.R"
    $r = "$BASE\03_proteomics\sgcb\results"
    $datasets = @("cptac","sim_proteomics","ecoli_lfq","sgsds_ratio2","sgsds_ratio2.5","tmt_mir","cptac_ups1","obrien_3species")

    # Normal DE (8)
    foreach ($d in $datasets) {
        $jobs.Add(@{ Script=$s; Args="$d"; Name="prot_sgcb_$d"; ResultDir=$r; ResultCSV="SGCB__$d.csv" }) | Out-Null
    }
    # Null FDR (8 x 10 = 80)
    foreach ($d in $datasets) {
        1..10 | ForEach-Object {
            $jobs.Add(@{ Script=$s; Args="$d $_"; Name="prot_sgcb_null${_}_$d"; ResultDir=$r; ResultCSV="null_${_}__SGCB__$d.csv" }) | Out-Null
        }
    }
}

# ==================== 5. Proteomics Other Methods ====================
if (($Category -eq "all" -or $Category -eq "proteomics") -and ($Tier -eq "all" -or $Tier -eq "other")) {
    $s = "$BASE\03_proteomics\other_methods\scripts\run_other_prot.R"
    $r = "$BASE\03_proteomics\other_methods\results"
    # DE methods (7) + robust (1) + DV methods (2) = 10
    $de_methods = @("limma","DEqMS","ROTS","proDA","msqrob2","ttest","MSstats","prolfqua")
    $dv_methods = @("dv_diffVar","dv_GAMLSS")
    $methods  = $de_methods + $dv_methods
    $datasets = @("cptac","sim_proteomics","ecoli_lfq","sgsds_ratio2","sgsds_ratio2.5","tmt_mir","cptac_ups1","obrien_3species")

    # Normal (10 x 8 = 80)
    foreach ($m in $methods) {
        foreach ($d in $datasets) {
            $jobs.Add(@{ Script=$s; Args="$m $d"; Name="prot_other_${m}_$d"; ResultDir=$r; ResultCSV="${m}__$d.csv" }) | Out-Null
        }
    }
    # Null FDR (10 x 8 x 10 = 800)
    foreach ($m in $methods) {
        foreach ($d in $datasets) {
            1..10 | ForEach-Object {
                $jobs.Add(@{ Script=$s; Args="$m $d $_"; Name="prot_other_null${_}_${m}_$d"; ResultDir=$r; ResultCSV="null_${_}__${m}__$d.csv" }) | Out-Null
            }
        }
    }
}

# ==================== 6. Paper Extra SGCB (standalone scripts feeding paper) ====================
# These scripts call sgcbDE directly and produce outputs referenced in the paper.
# They are included here so that -Tier sgcb reruns everything SGCB in one shot.
if ($Tier -eq "all" -or $Tier -eq "sgcb") {
    $BENCH = (Resolve-Path "$BASE\..\benchmark").Path

    $jobs.Add(@{
        Script    = "$BENCH\35_reviewer_concern_tests.R"
        Args      = ""
        Name      = "paper_sgcb_reviewer_checks"
        ResultDir = "$BENCH\output\reviewer_checks"
        ResultCSV = "00_reviewer_checks_summary.csv"
    }) | Out-Null

    $jobs.Add(@{
        Script    = "$BENCH\39_extended_ablation_benchmark.R"
        Args      = ""
        Name      = "paper_sgcb_ablation"
        ResultDir = "$BENCH\output\methodology_checks"
        ResultCSV = "ablation_simulation_summary.csv"
    }) | Out-Null

    $jobs.Add(@{
        Script    = "$BENCH\36_gg_gof_figure.R"
        Args      = ""
        Name      = "paper_sgcb_gg_gof"
        ResultDir = "$BENCH\output\paper_figures_v2"
        ResultCSV = "FigNew_gg_goodness_of_fit.png"
    }) | Out-Null
}

if ($SkipExisting) {
    $before = $jobs.Count
    $jobs = $jobs | Where-Object {
        $csv = [IO.Path]::Combine($_.ResultDir, $_.ResultCSV)
        -not (Test-Path $csv)
    }
    Write-Host "SkipExisting: $before -> $($jobs.Count) jobs (skipped $($before - $jobs.Count))"
}

if ($IncludeMethods.Count -gt 0) {
    $before = $jobs.Count
    $pattern = [string]::Join("|", ($IncludeMethods | ForEach-Object { [regex]::Escape($_) }))
    $jobs = $jobs | Where-Object {
        $_.Name -match "(^|_)(?:$pattern)(_|$)"
    }
    Write-Host "IncludeMethods: $before -> $($jobs.Count) jobs (methods=$($IncludeMethods -join ', '))"
}

$TIMEOUT_SEC = if ($SerialTail) { [int]::MaxValue } else { 14400 }    # serial tail mode: unlimited timeout

foreach ($j in $jobs) { $j.MemMB = Get-MemMB $j.Name; $j.TimeSec = Get-TimeSec $j.Name; $j.TimeoutSec = Get-JobTimeout $j.Name }

$tierGroups = @{}
foreach ($j in $jobs) {
    $tier = "$([math]::Ceiling($j.MemMB / 1000)) GB"
    if (-not $tierGroups.ContainsKey($tier)) { $tierGroups[$tier] = 0 }
    $tierGroups[$tier]++
}
$totalEstGB = [math]::Round(($jobs | ForEach-Object { $_.MemMB } | Measure-Object -Sum).Sum / 1024, 1)

Write-Host "============================================"
Write-Host " Dynamic Memory Scheduler"
Write-Host " Total jobs  : $($jobs.Count)"
Write-Host " Memory budget: ${MaxMemGB} GB  (OS reserve: ${ReserveGB} GB)"
Write-Host " Max parallel : $MaxJobs (safety cap)"
Write-Host " Timeout      : $(if ($SerialTail) { 'unlimited' } else { "${TIMEOUT_SEC}s" })"
Write-Host " Memory tiers :"
foreach ($k in ($tierGroups.Keys | Sort-Object)) { Write-Host "   $k : $($tierGroups[$k]) jobs" }
Write-Host " Sequential total: ${totalEstGB} GB"
$timeTiers = @{}
foreach ($j in $jobs) {
    $m = [math]::Floor($j.TimeSec / 60)
    if     ($m -lt 1)  { $tt = "<1 min" }
    elseif ($m -lt 5)  { $tt = "1-5 min" }
    elseif ($m -lt 15) { $tt = "5-15 min" }
    elseif ($m -lt 60) { $tt = "15-60 min" }
    else               { $tt = ">60 min" }
    if (-not $timeTiers.ContainsKey($tt)) { $timeTiers[$tt] = 0 }
    $timeTiers[$tt]++
}
Write-Host " Time tiers :"
foreach ($k in @("<1 min","1-5 min","5-15 min","15-60 min",">60 min")) {
    if ($timeTiers.ContainsKey($k)) { Write-Host "   $k : $($timeTiers[$k]) jobs" }
}
$totalHrs = [math]::Round(($jobs | ForEach-Object { $_.TimeSec } | Measure-Object -Sum).Sum / 3600, 1)
Write-Host " Sequential time : ~${totalHrs} hours"
Write-Host "============================================`n"

$SHORT_THRESHOLD_SEC = if ($SerialTail) { -1 } else { 600 }
$shortQ = [System.Collections.Generic.List[hashtable]]::new()
$longQ  = [System.Collections.Generic.List[hashtable]]::new()
$sorted = if ($SerialTail) { $jobs | Sort-Object { $_.TimeSec } -Descending } else { $jobs | Sort-Object { $_.TimeSec } }
foreach ($j in $sorted) {
    if ($j.TimeSec -le $SHORT_THRESHOLD_SEC) { $shortQ.Add($j) | Out-Null }
    else { $longQ.Add($j) | Out-Null }
}
Write-Host "  Queue split: short=$($shortQ.Count) long=$($longQ.Count) (threshold=${SHORT_THRESHOLD_SEC}s)"

$running      = [System.Collections.Generic.List[hashtable]]::new()
$completed    = 0
$failed       = 0
$failedList   = @()
$deferredList = @()
$memEstGB     = 0.0
$launchIdx    = 0
$totalJobs    = $shortQ.Count + $longQ.Count
$t_global     = Get-Date
$lastHB       = [datetime]::MinValue
$longLaunchN  = 0

function Reap {
    $toRemove = @()
    foreach ($p in $running) {
        $elapsed = ((Get-Date) - $p.StartTime).TotalSeconds
        if (-not $p.Proc.HasExited -and $elapsed -gt $p.TimeoutSec) {
            try { $p.Proc.Kill(); $p.Proc.WaitForExit(5000) } catch {}
            Write-Host "  [TIMEOUT] $($p.Name)" -ForegroundColor Yellow
        }
        if ($p.Proc.HasExited) {
            $p.Proc.WaitForExit()
            $ec  = $p.Proc.ExitCode
            $dur = [math]::Round(((Get-Date) - $p.StartTime).TotalSeconds, 1)

            $stdout = $p.OutSB.ToString()
            $stderr = $p.ErrSB.ToString()
            try {
                [IO.File]::WriteAllText($p.LogFile,
                    ($stdout + "`n--- STDERR ---`n" + $stderr),
                    [System.Text.UTF8Encoding]::new($false))
            } catch {}

            try { Unregister-Event -SourceIdentifier $p.OutEvtId -ErrorAction SilentlyContinue } catch {}
            try { Unregister-Event -SourceIdentifier $p.ErrEvtId -ErrorAction SilentlyContinue } catch {}
            try { $p.Proc.Close() } catch {}

            $csvPath = [IO.Path]::Combine($p.ResultDir, $p.ResultCSV)
            $ok = (Test-Path $csvPath)
            if (-not $ok -and $ec -eq 0) { $ok = $true }

            $script:memEstGB -= ($p.MemMB / 1024.0)
            if ($script:memEstGB -lt 0) { $script:memEstGB = 0 }

            $peakWSMB   = [math]::Round($p.PeakWS / 1MB, 1)
            $peakPrivMB = [math]::Round($p.PeakPrivate / 1MB, 1)

            if ($ok) {
                $script:completed++
                Write-Host "  [OK]   $($p.Name) (${dur}s, peak=${peakPrivMB}MB) [run=$($running.Count-1) mem=$([math]::Round($script:memEstGB,1))GB]" -ForegroundColor Green
            } else {
                $script:failed++
                $script:failedList += $p.Name
                Write-Host "  [FAIL] $($p.Name) (exit=$ec, ${dur}s, peak=${peakPrivMB}MB) log=$($p.LogFile)" -ForegroundColor Red
            }

            $histLine = [pscustomobject]@{
                Timestamp   = (Get-Date).ToString("s")
                Name        = $p.Name
                Signature   = Get-TaskSignature $p.Name
                ExitCode    = $ec
                DurationSec = $dur
                PeakWSMB    = $peakWSMB
                PeakPrivMB  = $peakPrivMB
                EstMemMB    = $p.MemMB
                Success     = $ok
            }
            $histLine | Export-Csv -Path $HISTORY_CSV -Append -NoTypeInformation -Encoding UTF8

            $toRemove += $p
        }
    }
    foreach ($r in $toRemove) { $running.Remove($r) | Out-Null }
}

function Launch-Job([hashtable]$job) {
    $script:launchIdx++
    New-Item -ItemType Directory -Force -Path $job.ResultDir | Out-Null
    $logFile = [IO.Path]::Combine($LOG_DIR, "$($job.Name).log")
    $jobThreadCount = Get-ThreadCountForJob $job.Name

    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName  = $RSCRIPT
    $psi.Arguments = "--encoding=UTF-8 `"$($job.Script)`" $($job.Args)"
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError  = $true
    $psi.CreateNoWindow = $true
    $psi.EnvironmentVariables["BENCHMARK_THREADS"]        = "$jobThreadCount"
    $psi.EnvironmentVariables["OMP_NUM_THREADS"]          = "$jobThreadCount"
    $psi.EnvironmentVariables["OPENBLAS_NUM_THREADS"]     = "$jobThreadCount"
    $psi.EnvironmentVariables["MKL_NUM_THREADS"]          = "$jobThreadCount"
    $psi.EnvironmentVariables["VECLIB_MAXIMUM_THREADS"]   = "$jobThreadCount"
    $psi.EnvironmentVariables["NUMEXPR_NUM_THREADS"]      = "$jobThreadCount"
    $psi.EnvironmentVariables["BLIS_NUM_THREADS"]         = "$jobThreadCount"
    $psi.EnvironmentVariables["R_DATATABLE_NUM_THREADS"]  = "$jobThreadCount"
    if ($SerialTail) {
        $psi.EnvironmentVariables["R_FUTURE_PLAN"]        = "sequential"
        $psi.EnvironmentVariables["R_FUTURE_FORK_ENABLE"] = "false"
    }
    $proc = [System.Diagnostics.Process]::Start($psi)

    $outSB = [System.Text.StringBuilder]::new()
    $errSB = [System.Text.StringBuilder]::new()
    $outEvtId = "out_$($job.Name)"
    $errEvtId = "err_$($job.Name)"
    Register-ObjectEvent -InputObject $proc -EventName OutputDataReceived `
        -SourceIdentifier $outEvtId -MessageData $outSB -Action {
            if ($null -ne $Event.SourceEventArgs.Data) {
                $Event.MessageData.AppendLine($Event.SourceEventArgs.Data)
            }
        } | Out-Null
    Register-ObjectEvent -InputObject $proc -EventName ErrorDataReceived `
        -SourceIdentifier $errEvtId -MessageData $errSB -Action {
            if ($null -ne $Event.SourceEventArgs.Data) {
                $Event.MessageData.AppendLine($Event.SourceEventArgs.Data)
            }
        } | Out-Null
    $proc.BeginOutputReadLine()
    $proc.BeginErrorReadLine()

    $script:memEstGB += ($job.MemMB / 1024.0)

    $tEst = if ($job.TimeSec -ge 3600) { "$([math]::Round($job.TimeSec/3600,1))h" } elseif ($job.TimeSec -ge 60) { "$([math]::Round($job.TimeSec/60))min" } else { "$($job.TimeSec)s" }
    Write-Host "  [LAUNCH $($script:launchIdx)/$totalJobs] $($job.Name) (~$($job.MemMB)MB, ~$tEst, thr=$jobThreadCount) [run=$($running.Count+1) mem=$([math]::Round($script:memEstGB,1))GB]"

    $running.Add(@{
        Proc        = $proc
        Name        = $job.Name
        Threads     = $jobThreadCount
        MemMB       = $job.MemMB
        TimeSec     = $job.TimeSec
        Priority    = $(if ($job.Name -match 'sgcb') { 2 } else { 1 })
        PeakWS      = [long]0
        PeakPrivate = [long]0
        StartTime   = Get-Date
        LogFile     = $logFile
        OutSB       = $outSB
        ErrSB       = $errSB
        OutEvtId    = $outEvtId
        ErrEvtId    = $errEvtId
        ResultDir   = $job.ResultDir
        ResultCSV   = $job.ResultCSV
        TimeoutSec  = $job.TimeoutSec
    }) | Out-Null
}

function Invoke-EmergencyFuse([double]$freeGB) {
    if ($freeGB -ge $CriticalGB -or $running.Count -eq 0) { return }
    $best = $null; $bestScore = -1e18
    foreach ($p in $running) {
        if ($p.Proc.HasExited) { continue }
        $elapsed = ((Get-Date) - $p.StartTime).TotalSeconds
        $remain  = [math]::Max(0, $p.TimeSec - $elapsed)
        $wsGB    = 0
        try { $wsGB = $p.Proc.WorkingSet64 / 1GB } catch {}
        $score = 3 * $wsGB + 0.001 * $remain - 3 * $p.Priority
        if ($score -gt $bestScore) { $bestScore = $score; $best = $p }
    }
    if ($null -eq $best) { return }
    $memMB = [math]::Round($best.PeakWS / 1MB)
    $el = [math]::Round(((Get-Date) - $best.StartTime).TotalSeconds)
    Write-Host "  [FUSE] FREE=${freeGB}GB < ${CriticalGB}GB! Killing $($best.Name) (peak=${memMB}MB, elapsed=${el}s, score=$([math]::Round($bestScore,1)))" -ForegroundColor Magenta
    try {
        $best.Proc.Kill()
        $null = $best.Proc.WaitForExit(3000)
    } catch {}
}

function Update-PeakAndGetActualGB {
    $totalWS = [long]0; $totalPM = [long]0
    try {
        $allRs = Get-Process Rscript -ErrorAction SilentlyContinue
        foreach ($r in $allRs) { $totalWS += $r.WorkingSet64; $totalPM += $r.PrivateMemorySize64 }
    } catch {}

    $estTotal = 0.0
    foreach ($p in $running) { if (-not $p.Proc.HasExited) { $estTotal += $p.MemMB } }
    if ($estTotal -gt 0) {
        foreach ($p in $running) {
            if ($p.Proc.HasExited) { continue }
            $frac = $p.MemMB / $estTotal
            $ws = [long]($totalWS * $frac)
            $pm = [long]($totalPM * $frac)
            if ($ws -gt $p.PeakWS) { $p.PeakWS = $ws }
            if ($pm -gt $p.PeakPrivate) { $p.PeakPrivate = $pm }
        }
    }
    [math]::Round($totalPM / 1GB, 2)
}

function Try-DequeueJob([System.Collections.Generic.List[hashtable]]$queue, [double]$freeGB) {
    $i = 0
    while ($i -lt $queue.Count) {
        $j   = $queue[$i]
        $jGB = $j.MemMB / 1024.0
        $soloRelax = ($running.Count -eq 0 -and $script:memEstGB -le 0.01 -and (($freeGB - $jGB) -ge $CriticalGB))
        if (($script:memEstGB + $jGB) -gt $MaxMemGB) { $i++; continue }
        if ((($freeGB - $jGB) -lt $ReserveGB) -and (-not $soloRelax)) { $i++; continue }
        if (-not (Test-MethodCapOK $j.Name))           { $i++; continue }
        $queue.RemoveAt($i)
        return $j
    }
    return $null
}

$fuseCount = 0
$throttleMsg = 0
$shortLaunched = 0

$initFree = Get-FreeMemGB
Write-Host "  [INIT] free=${initFree}GB  pressure=${PressureGB}GB  critical=${CriticalGB}GB  reserve=${ReserveGB}GB" -ForegroundColor Cyan

$pendingTotal = { $shortQ.Count + $longQ.Count }

while ((& $pendingTotal) -gt 0 -or $running.Count -gt 0) {
    Reap

    $freeGB   = Get-FreeMemGB
    $actualGB = Update-PeakAndGetActualGB

    if ($freeGB -lt $CriticalGB -and $running.Count -gt 0) {
        Invoke-EmergencyFuse $freeGB
        $fuseCount++
        Start-Sleep -Seconds 5
        Reap
        $freeGB = Get-FreeMemGB
    }

    $canLaunch = ($freeGB -ge $PressureGB)
    if (-not $canLaunch -and $running.Count -eq 0 -and (& $pendingTotal) -gt 0) {
        $canLaunch = $true
    }
    if (-not $canLaunch -and (& $pendingTotal) -gt 0) {
        $throttleMsg++
        if ($throttleMsg -le 3 -or ($throttleMsg % 15) -eq 0) {
            Write-Host "  [THROTTLE] free=${freeGB}GB < pressure=${PressureGB}GB, waiting... (run=$($running.Count) short=$($shortQ.Count) long=$($longQ.Count))" -ForegroundColor Yellow
        }
    } else { $throttleMsg = 0 }

    $dynMaxJobs = $MaxJobs
    if ($running.Count -gt 2 -and $freeGB -lt ($ReserveGB + 4)) {
        $dynMaxJobs = [math]::Max(2, [math]::Floor($running.Count * 0.7))
    }

    while ($canLaunch -and (& $pendingTotal) -gt 0 -and $running.Count -lt $dynMaxJobs) {
        $job = $null
        $tryLong = ($shortLaunched -ge 3 -and $longQ.Count -gt 0)
        if ($tryLong) {
            $job = Try-DequeueJob $longQ $freeGB
            if ($null -ne $job) { $shortLaunched = 0 }
        }
        if ($null -eq $job -and $shortQ.Count -gt 0) {
            $job = Try-DequeueJob $shortQ $freeGB
            if ($null -ne $job) { $shortLaunched++ }
        }
        if ($null -eq $job -and $longQ.Count -gt 0) {
            $job = Try-DequeueJob $longQ $freeGB
        }
        if ($null -eq $job) { break }
        Launch-Job $job
        $freeGB -= ($job.MemMB / 1024.0)
    }

    if ((& $pendingTotal) -gt 0 -and $running.Count -eq 0 -and $canLaunch) {
        Write-Host "`n  [DEFER] $($shortQ.Count + $longQ.Count) tasks exceed budget or method cap:" -ForegroundColor Yellow
        foreach ($d in $shortQ) { Write-Host "    - $($d.Name) (~$($d.MemMB)MB)" -ForegroundColor Yellow; $deferredList += $d.Name }
        foreach ($d in $longQ)  { Write-Host "    - $($d.Name) (~$($d.MemMB)MB)" -ForegroundColor Yellow; $deferredList += $d.Name }
        $shortQ.Clear(); $longQ.Clear()
    }

    $now = Get-Date
    if (($now - $lastHB).TotalSeconds -ge 30 -and ($running.Count -gt 0 -or (& $pendingTotal) -gt 0)) {
        $script:lastHB = $now
        $el = [math]::Round(($now - $t_global).TotalMinutes, 1)
        Write-Host "  [STATUS ${el}min] done=$completed fail=$failed run=$($running.Count) short=$($shortQ.Count) long=$($longQ.Count) est=$([math]::Round($memEstGB,1))GB actual=$([math]::Round($actualGB,1))GB free=$([math]::Round($freeGB,1))GB dynMax=$dynMaxJobs fuse=$fuseCount" -ForegroundColor Cyan
    }

    if ($running.Count -gt 0 -or (& $pendingTotal) -gt 0) { Start-Sleep -Milliseconds 500 }
}

$elapsed_total = [math]::Round(((Get-Date) - $t_global).TotalMinutes, 1)
Write-Host "`n============================================"
Write-Host " DONE in ${elapsed_total} min"
Write-Host "   OK: $completed  FAIL: $failed  FUSE: $fuseCount"
if ($failedList.Count -gt 0) {
    Write-Host "   Failed:" -ForegroundColor Red
    $failedList | ForEach-Object { Write-Host "     - $_" -ForegroundColor Red }
}
if ($deferredList.Count -gt 0) {
    Write-Host "   Deferred (memory too large, run later):" -ForegroundColor Yellow
    $deferredList | ForEach-Object { Write-Host "     - $_" -ForegroundColor Yellow }
}
Write-Host " Logs: $LOG_DIR"
Write-Host "============================================"

