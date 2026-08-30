param(
    [string]$Rscript = "Rscript.exe",
    [string]$Quarto = "quarto",
    [string]$OutputDirectory = ""
)

$ErrorActionPreference = "Stop"
$env:LC_ALL = "C"
$env:LANG = "C"
$project = Split-Path -Parent $PSScriptRoot
$removeTemporary = [string]::IsNullOrWhiteSpace($OutputDirectory)
$temporary = if ($removeTemporary) {
    Join-Path ([System.IO.Path]::GetTempPath()) ("dnaprs-smoke-" + [guid]::NewGuid().ToString("N"))
}
else {
    [System.IO.Path]::GetFullPath($OutputDirectory)
}
New-Item -ItemType Directory -Path $temporary -Force | Out-Null

try {
    Get-ChildItem -LiteralPath (Join-Path $project "bin") -Filter "*.R" | ForEach-Object {
        $parsePath = $_.FullName.Replace("\", "/").Replace("'", "\\'")
        & $Rscript -e "parse(file='$parsePath')" | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "R parse failed for $($_.Name)." }
    }

    Push-Location $temporary
    try {
        & $Rscript (Join-Path $project "bin/validate_manifests.R") `
            --run-name "smoke" `
            --target-manifest (Join-Path $project "tests/data/target_manifest.tsv") `
            --gwas-manifest (Join-Path $project "tests/data/gwas_manifest.tsv") `
            --reference-manifest (Join-Path $project "tests/data/reference_manifest.tsv") `
            --phenotype-file "tests/data/phenotype.tsv" `
            --phenotype-models (Join-Path $project "tests/data/phenotype_models.tsv") `
            --methods "plink_ct,sbayesrc" `
            --genome-build "GRCh37" `
            --seed "20260829" `
            --report-enabled "true" `
            --launch-dir $project
        if ($LASTEXITCODE -ne 0) { throw "Manifest-preflight smoke test failed." }

        & $Rscript (Join-Path $project "bin/harmonise_gwas.R") `
            --input (Join-Path $project "tests/data/gwas/mdd.tsv") `
            --trait-id "MDD" `
            --prs-name "MDD_PRS" `
            --effect-type "log_or" `
            --sample-size "100000" `
            --snp-col "SNP" `
            --chr-col "CHR" `
            --bp-col "BP" `
            --effect-allele-col "A1" `
            --other-allele-col "A2" `
            --beta-col "BETA" `
            --se-col "SE" `
            --p-col "P" `
            --freq-col "EAF" `
            --n-col "N"
        if ($LASTEXITCODE -ne 0) { throw "GWAS-harmonisation smoke test failed." }

        & $Rscript (Join-Path $project "bin/build_ct_weights.R") `
            --cojo "MDD.cojo.ma" `
            --clumps (Join-Path $project "tests/data/gwas/mdd.clumps") `
            --trait-id "MDD" `
            --prs-name "MDD_PRS"
        if ($LASTEXITCODE -ne 0) { throw "C+T weight smoke test failed." }

        & $Rscript (Join-Path $project "bin/parse_plink_score.R") `
            --score (Join-Path $project "tests/data/scores/plink.sscore") `
            --used (Join-Path $project "tests/data/scores/plink.sscore.vars") `
            --cohort "TEST" `
            --role "target" `
            --trait-id "MDD" `
            --prs-name "MDD_PRS"
        if ($LASTEXITCODE -ne 0) { throw "PLINK-score parser smoke test failed." }

        $score1 = Join-Path $project "tests/data/scores/test_plink.score.tsv"
        $score2 = Join-Path $project "tests/data/scores/test_sbayesrc.score.tsv"
        & $Rscript (Join-Path $project "bin/combine_scores.R") --scores "$score1,$score2"
        if ($LASTEXITCODE -ne 0) { throw "Score-combination smoke test failed." }

        & $Rscript (Join-Path $project "bin/phenotype_association.R") `
            --scores "prs_scores_long.tsv" `
            --phenotype (Join-Path $project "tests/data/phenotype_glm.tsv") `
            --models (Join-Path $project "tests/data/phenotype_models_glm.tsv")
        if ($LASTEXITCODE -ne 0) { throw "glm2 phenotype-model smoke test failed." }
        $glmResult = Import-Csv "phenotype_associations.tsv" -Delimiter "`t"
        $glmModel = Import-Csv "phenotype_models_fitted.tsv" -Delimiter "`t"
        if (($glmResult | Where-Object status -ne "ESTIMATED").Count -gt 0 -or
            ($glmModel | Where-Object estimator -ne "glm2::glm2").Count -gt 0) {
            throw "glm2 phenotype models were not estimated with the declared engine."
        }

        & $Rscript (Join-Path $project "bin/phenotype_association.R") `
            --scores "prs_scores_long.tsv" `
            --phenotype (Join-Path $project "tests/data/phenotype_lme.tsv") `
            --models (Join-Path $project "tests/data/phenotype_models_lme.tsv")
        if ($LASTEXITCODE -ne 0) { throw "lme4 phenotype-model smoke test failed." }
        $lmeResult = Import-Csv "phenotype_associations.tsv" -Delimiter "`t"
        $lmeModel = Import-Csv "phenotype_models_fitted.tsv" -Delimiter "`t"
        if (($lmeResult | Where-Object status -ne "ESTIMATED").Count -gt 0 -or
            ($lmeModel | Where-Object estimator -ne "lme4::lmer").Count -gt 0) {
            throw "lme4 phenotype models were not estimated with the declared engine."
        }

        & $Rscript (Join-Path $project "bin/phenotype_association.R") `
            --scores "prs_scores_long.tsv" `
            --phenotype (Join-Path $project "tests/data/phenotype.tsv") `
            --models (Join-Path $project "tests/data/phenotype_models.tsv")
        if ($LASTEXITCODE -ne 0) { throw "Phenotype-model smoke test failed." }

        $generationQC = @(
            (Join-Path $project "tests/data/generation_qc/harmonisation_qc.tsv"),
            (Join-Path $project "tests/data/generation_qc/plink_qc.tsv"),
            (Join-Path $project "tests/data/generation_qc/tidy_qc.tsv"),
            (Join-Path $project "tests/data/generation_qc/impute_qc.tsv"),
            (Join-Path $project "tests/data/generation_qc/model_qc.tsv")
        ) -join ","
        & $Rscript (Join-Path $project "bin/summarise_generation_qc.R") `
            --qc-files $generationQC `
            --score-qc "score_qc.tsv"
        if ($LASTEXITCODE -ne 0) { throw "Variant-flow smoke test failed." }

        $reportInputs = Join-Path $temporary "report_inputs"
        $reportRoot = Join-Path $temporary "reports"
        $reportProject = Join-Path $temporary "report_project"
        @($reportInputs, $reportRoot, $reportProject) | ForEach-Object {
            if (Test-Path -LiteralPath $_) {
                Remove-Item -LiteralPath $_ -Recurse -Force
            }
        }
        New-Item -ItemType Directory -Path $reportInputs -Force | Out-Null
        "DNAPRS_SMOKE:`n  R: synthetic-test" |
            Set-Content -LiteralPath (Join-Path $temporary "software_versions.yml") -Encoding utf8

        $reportFiles = @(
            @{ Path = "target_manifest.validated.tsv"; Publish = "run/preflight" },
            @{ Path = "gwas_manifest.validated.tsv"; Publish = "run/preflight" },
            @{ Path = "reference_manifest.validated.tsv"; Publish = "run/preflight" },
            @{ Path = "phenotype_models.validated.tsv"; Publish = "run/preflight" },
            @{ Path = "resolved_params.yaml"; Publish = "run/preflight" },
            @{ Path = "input_checksums.tsv"; Publish = "run/preflight" },
            @{ Path = "preflight_qc.tsv"; Publish = "qc/preflight" },
            @{ Path = "MDD.cojo.ma"; Publish = "gwas/MDD" },
            @{ Path = "MDD.harmonisation_qc.tsv"; Publish = "qc/gwas/MDD" },
            @{ Path = "MDD.plink_ct.weights.tsv"; Publish = "plink_ct/MDD" },
            @{ Path = "MDD.plink_ct.weight_qc.tsv"; Publish = "qc/plink_ct/MDD" },
            @{ Path = "TEST.MDD.plink_ct.score.tsv"; Publish = "scores/TEST/MDD" },
            @{ Path = "TEST.MDD.plink_ct.score_qc.tsv"; Publish = "qc/scores/TEST/MDD" },
            @{ Path = "prs_scores_long.tsv"; Publish = "scores/combined" },
            @{ Path = "prs_scores_wide.tsv"; Publish = "scores/combined" },
            @{ Path = "score_qc.tsv"; Publish = "qc/scores" },
            @{ Path = "method_concordance.tsv"; Publish = "qc/scores" },
            @{ Path = "variant_flow.tsv"; Publish = "qc/generation" },
            @{ Path = "phenotype_associations.tsv"; Publish = "phenotype" },
            @{ Path = "phenotype_models_fitted.tsv"; Publish = "phenotype" },
            @{ Path = "phenotype_plot_data.tsv"; Publish = "phenotype" },
            @{ Path = "software_versions.yml"; Publish = "pipeline_info" }
        )
        foreach ($record in $reportFiles) {
            Copy-Item -LiteralPath (Join-Path $temporary $record.Path) `
                -Destination (Join-Path $reportInputs $record.Path) -Force
        }

        $sbayesrcFiles = @(
            @{ Source = "tests/data/generation_qc/tidy_qc.tsv"; Name = "MDD.tidy_qc.tsv" },
            @{ Source = "tests/data/generation_qc/impute_qc.tsv"; Name = "MDD.impute_qc.tsv" },
            @{ Source = "tests/data/generation_qc/model_qc.tsv"; Name = "MDD.sbayesrc.model_qc.tsv" }
        )
        foreach ($record in $sbayesrcFiles) {
            Copy-Item -LiteralPath (Join-Path $project $record.Source) `
                -Destination (Join-Path $reportInputs $record.Name) -Force
            $reportFiles += @{ Path = $record.Name; Publish = "qc/sbayesrc/MDD" }
        }

        $manifestLines = @("publish_path`tfile_name")
        $manifestLines += $reportFiles | ForEach-Object { "$($_.Publish)`t$($_.Path)" }
        $outputManifest = Join-Path $temporary "output_files.tsv"
        $manifestLines | Set-Content -LiteralPath $outputManifest -Encoding utf8

        Copy-Item -LiteralPath (Join-Path $project "assets/report") `
            -Destination $reportProject -Recurse
        $env:DNAPRS_REPORT_INPUTS = $reportInputs
        $env:DNAPRS_OUTPUT_MANIFEST = $outputManifest
        $env:QUARTO_VERSION = (& $Quarto --version | Select-Object -First 1)
        Push-Location $reportProject
        try {
            & $Rscript "prepare-report.R"
            if ($LASTEXITCODE -ne 0) { throw "Report-preparation smoke test failed." }
            $previousErrorAction = $ErrorActionPreference
            $ErrorActionPreference = "Continue"
            $null = @(& $Quarto render "." 2>&1)
            $quartoStatus = $LASTEXITCODE
            $ErrorActionPreference = $previousErrorAction
            if ($quartoStatus -ne 0) { throw "Report-render smoke test failed." }
            Copy-Item -LiteralPath "_site" -Destination $reportRoot -Recurse
        }
        finally {
            Pop-Location
        }

        @("execution_report.html", "execution_timeline.html", "pipeline_dag.html") |
            ForEach-Object {
                "<html><body><p>Synthetic Nextflow test artifact.</p></body></html>" |
                    Set-Content -LiteralPath (Join-Path $reportRoot "provenance/$_") -Encoding utf8
            }
        "task_id`tname`tstatus" |
            Set-Content -LiteralPath (Join-Path $reportRoot "provenance/execution_trace.tsv") -Encoding utf8

        $pageNames = @(
            "index.html",
            "inputs.html",
            "target-gwas.html",
            "plink.html",
            "sbayesrc.html",
            "prs.html",
            "phenotype.html",
            "dictionary.html",
            "logs.html"
        )
        $required = @(
            "prs_scores_long.tsv",
            "prs_scores_wide.tsv",
            "MDD.cojo.ma",
            "MDD.plink_ct.weights.tsv",
            "TEST.MDD.plink_ct.score.tsv",
            "score_qc.tsv",
            "method_concordance.tsv",
            "phenotype_associations.tsv",
            "phenotype_models_fitted.tsv",
            "phenotype_plot_data.tsv",
            "variant_flow.tsv",
            "reports/provenance/output_files.tsv",
            "reports/provenance/data_dictionary.tsv",
            "reports/provenance/figure_manifest.tsv",
            "reports/provenance/report_state.tsv",
            "reports/provenance/report_software_versions.tsv"
        )
        $required += $pageNames | ForEach-Object { "reports/$_" }
        $figureNames = @(
            "variant_generation_flow",
            "participant_prs_distributions",
            "scoring_variant_counts",
            "cross_score_correlations",
            "method_agreement",
            "participant_method_ranks",
            "phenotype_prs_effects",
            "phenotype_incremental_fit",
            "adjusted_phenotype_prs_patterns",
            "phenotype_observed_fitted",
            "phenotype_residual_diagnostics"
        )
        foreach ($name in $figureNames) {
            $required += "reports/figures/tiff/$name.tiff"
            $required += "reports/figures/png/$name.png"
            $required += "reports/figures/jpeg/$name.jpeg"
        }
        foreach ($path in $required) {
            if (-not (Test-Path -LiteralPath $path)) { throw "Smoke-test output is missing: $path" }
        }
        $navigation = @(
            "Overview", "Inputs", "Target and GWAS", "PLINK", "SBayesRC",
            "PRS", "Phenotype", "Dictionary", "Logs"
        )
        foreach ($pageName in $pageNames) {
            $pagePath = Join-Path $reportRoot $pageName
            $reportHTML = Get-Content -LiteralPath $pagePath -Raw
            foreach ($pattern in $navigation) {
                if ($reportHTML -notmatch [regex]::Escape($pattern)) {
                    throw "$pageName is missing navigation content: $pattern"
                }
            }
            foreach ($pattern in @(
                "dnaprs-report.css", "dnaprs-report.js",
                "nf-core-dnaprs_logo_light.svg", "nf-core-dnaprs_logo_dark.svg"
            )) {
                if ($reportHTML -notmatch [regex]::Escape($pattern)) {
                    throw "$pageName is missing its report asset: $pattern"
                }
            }
            if ($reportHTML -match "POLYGENIC SCORING") {
                throw "$pageName contains the removed logo subtitle."
            }

            $localLinks = [regex]::Matches($reportHTML, '(?:href|src)="([^"]+)"') |
                ForEach-Object { $_.Groups[1].Value } |
                Where-Object { $_ -notmatch '^(#|https?:|mailto:|data:|javascript:)' } |
                ForEach-Object { ($_ -split '[#?]')[0] } |
                Where-Object { -not [string]::IsNullOrWhiteSpace($_) } |
                Sort-Object -Unique
            foreach ($link in $localLinks) {
                $decoded = [System.Uri]::UnescapeDataString($link)
                if (-not (Test-Path -LiteralPath (Join-Path $reportRoot $decoded))) {
                    throw "$pageName contains a broken local link: $link"
                }
            }
        }

        $overviewHTML = Get-Content -LiteralPath (Join-Path $reportRoot "index.html") -Raw
        foreach ($pattern in @("Polygenic risk score report", "data-table-viewer")) {
            if ($overviewHTML -notmatch [regex]::Escape($pattern)) {
                throw "The overview is missing required content: $pattern"
            }
        }
        $reportScript = Get-Content -LiteralPath `
            (Join-Path $reportRoot "assets/dnaprs-report.js") -Raw
        foreach ($pattern in @(
            "data-dnaprs-theme", "data-table-viewer", "data-table-expand",
            "data-figure-gallery", "data-figure-zoom-in", "dnaprs:figure-zoom"
        )) {
            if ($reportScript -notmatch [regex]::Escape($pattern)) {
                throw "The report script is missing required behaviour: $pattern"
            }
        }
        foreach ($pageName in @("target-gwas.html", "prs.html", "phenotype.html")) {
            $galleryHTML = Get-Content -LiteralPath (Join-Path $reportRoot $pageName) -Raw
            foreach ($pattern in @(
                "data-figure-gallery", "data-figure-stage", "data-figure-zoom-reset",
                "Download TIFF", "Download PNG", "Download JPEG"
            )) {
                if ($galleryHTML -notmatch [regex]::Escape($pattern)) {
                    throw "$pageName is missing figure-gallery content: $pattern"
                }
            }
        }

        $figureDictionary = Import-Csv `
            (Join-Path $reportRoot "provenance/figure_manifest.tsv") -Delimiter "`t"
        if ($figureDictionary.Count -ne 11) {
            throw "The figure dictionary does not contain all 11 expected figures."
        }
        if (($figureDictionary | Where-Object {
            -not $_.inspection -or -not $_.tiff -or -not $_.png -or -not $_.jpeg
        }).Count -gt 0) {
            throw "The figure dictionary contains incomplete inspection or download fields."
        }
    }
    finally {
        Pop-Location
    }

    Write-Host "dnaprs R and report smoke tests passed."
    if (-not $removeTemporary) { Write-Host "Smoke-test outputs: $temporary" }
}
finally {
    if ($removeTemporary -and (Test-Path -LiteralPath $temporary)) {
        Remove-Item -LiteralPath $temporary -Recurse -Force
    }
}
