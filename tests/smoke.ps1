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
$working = Join-Path $temporary ".work"
if (Test-Path -LiteralPath $working) {
    Remove-Item -LiteralPath $working -Recurse -Force
}
New-Item -ItemType Directory -Path $working -Force | Out-Null

try {
    Get-ChildItem -LiteralPath (Join-Path $project "bin") -Filter "*.R" | ForEach-Object {
        $parsePath = $_.FullName.Replace("\", "/").Replace("'", "\\'")
        & $Rscript -e "parse(file='$parsePath')" | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "R parse failed for $($_.Name)." }
    }

    Push-Location $working
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
            --target-imputation "true" `
            --launch-dir $project
        if ($LASTEXITCODE -ne 0) { throw "Input-validation smoke test failed." }
        $preflightResult = Import-Csv "input_checks.tsv" -Delimiter "`t"
        $referenceChecksumCheck = $preflightResult | Where-Object check -eq "reference_checksums"
        if ($referenceChecksumCheck.status -ne "PASS" -or
            $referenceChecksumCheck.value -notmatch "supplied and verified") {
            throw "Reference-checksum input-validation status was not recorded."
        }

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

        & $Rscript (Join-Path $project "bin/harmonise_gwas.R") `
            --input (Join-Path $project "tests/data/gwas/mdd_raw.tsv") `
            --trait-id "MDDRAW" `
            --prs-name "MDD_PRS" `
            --effect-type "log_or" `
            --sample-size "2000702" `
            --snp-col "ID" `
            --chr-col "CHROM" `
            --bp-col "POS" `
            --effect-allele-col "EA" `
            --other-allele-col "NEA" `
            --beta-col "BETA" `
            --se-col "SE" `
            --p-col "PVAL" `
            --case-freq-col "FCAS" `
            --control-freq-col "FCON" `
            --case-n-col "NCAS" `
            --control-n-col "NCON" `
            --n-col "NEFF" `
            --info-col "IMPINFO" `
            --info-min "0.8" `
            --maf-min "0.01" `
            --source-format "pgc_mdd2025"
        if ($LASTEXITCODE -ne 0) { throw "Raw MDD GWAS-adapter smoke test failed." }
        $rawMdd = Import-Csv "MDDRAW.cojo.ma" -Delimiter "`t"
        if ([math]::Abs(([double]$rawMdd[0].freq) - 0.35) -gt 1e-10) {
            throw "Raw MDD case/control frequencies were not combined correctly."
        }

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

        & $Rscript (Join-Path $project "bin/compare_direct_score.R") `
            --primary-score "TEST.MDD.plink_ct.score.tsv" `
            --primary-used (Join-Path $project "tests/data/scores/plink.sscore.vars") `
            --direct-score "TEST.MDD.plink_ct.score.tsv" `
            --direct-used (Join-Path $project "tests/data/scores/plink.sscore.vars")
        if ($LASTEXITCODE -ne 0) { throw "Direct-genotype score-sensitivity smoke test failed." }
        $sensitivityResult = Import-Csv "TEST.MDD.plink_ct.sensitivity_qc.tsv" -Delimiter "`t"
        if ($sensitivityResult.status -ne "PASS" -or [double]$sensitivityResult.pearson_r -ne 1) {
            throw "Direct-genotype score sensitivity did not record exact agreement."
        }

        $score1 = Join-Path $project "tests/data/scores/test_plink.score.tsv"
        $score2 = Join-Path $project "tests/data/scores/test_sbayesrc.score.tsv"
        & $Rscript (Join-Path $project "bin/combine_scores.R") `
            --scores "$score1,$score2" `
            --participant-decisions (Join-Path $project "tests/data/generation_qc/participant_decisions.tsv")
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

        $reportInputs = Join-Path $working "report_inputs"
        $reportRoot = Join-Path $working "reports"
        $reportProject = Join-Path $working "report_project"
        @($reportInputs, $reportRoot, $reportProject) | ForEach-Object {
            if (Test-Path -LiteralPath $_) {
                Remove-Item -LiteralPath $_ -Recurse -Force
            }
        }
        New-Item -ItemType Directory -Path $reportInputs -Force | Out-Null
        "DNAPRS_SMOKE:`n  R: synthetic-test" |
            Set-Content -LiteralPath (Join-Path $working "software_versions.yml") -Encoding utf8
        Get-ChildItem -LiteralPath (Join-Path $project "tests/data/genotype_eda") -File |
            Copy-Item -Destination $working -Force
        $targetReportFixtures = @(
            @{ Source = "target_sample_decisions.tsv"; Name = "TEST.sample_decisions.tsv" },
            @{ Source = "target_variant_decisions.tsv"; Name = "TEST.variant_decisions.tsv" },
            @{ Source = "participant_decisions.tsv"; Name = "TEST.participant_decisions.tsv" },
            @{ Source = "target_ancestry.tsv"; Name = "TEST.target_ancestry.tsv" },
            @{ Source = "reference_projection.tsv"; Name = "TEST.reference_projection.tsv" },
            @{ Source = "ancestry_summary.tsv"; Name = "TEST.ancestry_summary.tsv" },
            @{ Source = "target_imputation_qc.tsv"; Name = "TEST.imputation_qc.tsv" },
            @{ Source = "target_imputation_manifest.tsv"; Name = "TEST.imputation_manifest.tsv" },
            @{ Source = "target_imputation_dr2.tsv"; Name = "TEST.imputation_dr2.tsv" },
            @{ Source = "plink_sensitivity.tsv"; Name = "TEST.MDD.plink_ct.sensitivity.tsv" },
            @{ Source = "plink_sensitivity_qc.tsv"; Name = "TEST.MDD.plink_ct.sensitivity_qc.tsv" }
        )
        foreach ($record in $targetReportFixtures) {
            Copy-Item -LiteralPath (Join-Path $project "tests/data/generation_qc/$($record.Source)") `
                -Destination (Join-Path $working $record.Name) -Force
        }

        $reportFiles = @(
            @{ Path = "targets.tsv"; Publish = "inputs" },
            @{ Path = "gwas.tsv"; Publish = "inputs" },
            @{ Path = "references.tsv"; Publish = "inputs" },
            @{ Path = "models.tsv"; Publish = "inputs" },
            @{ Path = "run_settings.yml"; Publish = "inputs" },
            @{ Path = "input_checksums.tsv"; Publish = "inputs" },
            @{ Path = "input_checks.tsv"; Publish = "inputs/checks" },
            @{ Path = "TEST.genotype_eda_summary.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.genotype_eda_checks.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.chromosome_counts.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.marker_density.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.identifier_classes.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.allele_states.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.sample_missingness.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.variant_missingness.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.allele_frequency.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.heterozygosity.tsv"; Publish = "genotype_eda/TEST" },
            @{ Path = "TEST.target_prep_summary.tsv"; Publish = "target_prep/TEST" },
            @{ Path = "TEST.target_qc.tsv"; Publish = "target_qc/TEST" },
            @{ Path = "TEST.sample_decisions.tsv"; Publish = "target_qc/TEST" },
            @{ Path = "TEST.variant_decisions.tsv"; Publish = "target_qc/TEST" },
            @{ Path = "TEST.participant_decisions.tsv"; Publish = "target_qc/TEST" },
            @{ Path = "TEST.target_ancestry.tsv"; Publish = "target_qc/TEST/ancestry" },
            @{ Path = "TEST.reference_projection.tsv"; Publish = "target_qc/TEST/ancestry" },
            @{ Path = "TEST.ancestry_summary.tsv"; Publish = "target_qc/TEST/ancestry" },
            @{ Path = "TEST.imputation_qc.tsv"; Publish = "target_imputation/TEST" },
            @{ Path = "TEST.imputation_manifest.tsv"; Publish = "target_imputation/TEST" },
            @{ Path = "TEST.imputation_dr2.tsv"; Publish = "target_imputation/TEST" },
            @{ Path = "TEST.MDD.plink_ct.sensitivity.tsv"; Publish = "scores/TEST/MDD/sensitivity" },
            @{ Path = "TEST.MDD.plink_ct.sensitivity_qc.tsv"; Publish = "qc/scores/TEST/MDD/sensitivity" },
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
            @{ Path = "phenotype_permutations.tsv"; Publish = "phenotype" },
            @{ Path = "phenotype_influence.tsv"; Publish = "phenotype" },
            @{ Path = "phenoPRS.csv"; Publish = "phenotype" },
            @{ Path = "software_versions.yml"; Publish = "pipeline_info" }
        )
        foreach ($record in $reportFiles) {
            Copy-Item -LiteralPath (Join-Path $working $record.Path) `
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
        Copy-Item -LiteralPath (Join-Path $project "tests/data/generation_qc/sbayesrc_weights.tsv") `
            -Destination (Join-Path $reportInputs "MDD.sbayesrc.txt") -Force
        $reportFiles += @{ Path = "MDD.sbayesrc.txt"; Publish = "sbayesrc/MDD/model" }

        $manifestLines = @("publish_path`tfile_name")
        $manifestLines += $reportFiles | ForEach-Object { "$($_.Publish)`t$($_.Path)" }
        $outputManifest = Join-Path $working "output_files.tsv"
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
            $quartoOutput = @(& $Quarto render "." 2>&1)
            $quartoStatus = $LASTEXITCODE
            $ErrorActionPreference = $previousErrorAction
            if ($quartoStatus -ne 0) {
                throw "Report-render smoke test failed.`n$($quartoOutput -join [Environment]::NewLine)"
            }
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
            Set-Content -LiteralPath (Join-Path $reportRoot "provenance/execution_trace.txt") -Encoding utf8

        $pageNames = @(
            "index.html",
            "genotype-eda.html",
            "target-prep.html",
            "target-qc.html",
            "target-imputation.html",
            "gwas-qc.html",
            "plink-prs.html",
            "sbayesrc-prs.html",
            "phenotype.html",
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
            "phenotype_permutations.tsv",
            "phenotype_influence.tsv",
            "phenoPRS.csv",
            "variant_flow.tsv",
            "reports/downloads/dnaprs_report_tables.xlsx",
            "reports/provenance/output_files.tsv",
            "reports/provenance/data_dictionary.tsv",
            "reports/provenance/figure_manifest.tsv",
            "reports/provenance/report_state.tsv",
            "reports/provenance/report_software_versions.tsv"
        )
        $required += $pageNames | ForEach-Object { "reports/$_" }
        $figureNames = @(
            "genotype_variants_by_chromosome",
            "genotype_marker_density",
            "genotype_identifier_classes",
            "genotype_allele_states",
            "genotype_sample_missingness",
            "genotype_variant_missingness",
            "genotype_allele_frequency",
            "genotype_heterozygosity_missingness",
            "target_preparation_summary",
            "target_qc_summary",
            "target_sample_decisions",
            "target_variant_decisions",
            "participant_analysis_eligibility",
            "reference_ancestry_projection",
            "reference_ancestry_distance",
            "target_imputation_counts",
            "target_imputation_dr2",
            "gwas_variant_retention",
            "gwas_qq",
            "gwas_effect_distribution",
            "gwas_minor_allele_frequency",
            "gwas_manhattan",
            "variant_generation_flow",
            "plink_scoring_variant_counts",
            "plink_participant_prs_distributions",
            "plink_cross_score_correlations",
            "plink_imputation_sensitivity",
            "plink_typed_imputed_coverage",
            "sbayesrc_scoring_variant_counts",
            "sbayesrc_participant_prs_distributions",
            "sbayesrc_cross_score_correlations",
            "sbayesrc_posterior_effect_distribution",
            "sbayesrc_pip_distribution",
            "method_agreement",
            "participant_method_ranks",
            "phenotype_outcome_distributions",
            "phenotype_outcome_correlations",
            "phenotype_prs_effects",
            "phenotype_incremental_fit",
            "adjusted_phenotype_prs_patterns",
            "phenotype_observed_fitted",
            "phenotype_residual_diagnostics",
            "phenotype_residual_permutations",
            "phenotype_leave_one_out"
        )
        foreach ($name in $figureNames) {
            $required += "reports/figures/svg/$name.svg"
            $required += "reports/figures/tiff/$name.tiff"
            $required += "reports/figures/png/$name.png"
            $required += "reports/figures/jpeg/$name.jpeg"
        }
        foreach ($path in $required) {
            if (-not (Test-Path -LiteralPath $path)) { throw "Smoke-test output is missing: $path" }
        }
        $navigation = @(
            "Overview", "Genotype EDA", "Target PREP", "Target QC", "Target Imputation", "GWAS QC",
            "PLINK PRS", "SBayesRC PRS", "Phenotype", "Logs"
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
            if ($reportHTML -notmatch 'class="dnaprs-github-footer"' -or
                $reportHTML -notmatch 'href="https://github.com/paulYRP/dnaprs"') {
                throw "$pageName is missing the standalone GitHub footer icon."
            }
            if ($reportHTML -match '<nav[\s\S]*?href="https://github.com/paulYRP/dnaprs"[\s\S]*?</nav>') {
                throw "$pageName incorrectly places GitHub inside the report navigation."
            }
            if ($reportHTML -match '>\s*GitHub\s*<') {
                throw "$pageName contains visible GitHub link text instead of an icon only."
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
        if ($overviewHTML -match [regex]::Escape("This website records the inputs, score-generation checks")) {
            throw "The Overview contains the removed hero sentence."
        }
        foreach ($pageName in $pageNames | Where-Object { $_ -ne "index.html" }) {
            $pageHTML = Get-Content -LiteralPath (Join-Path $reportRoot $pageName) -Raw
            if ($pageHTML -match "dnaprs-opening-card") {
                throw "$pageName contains an opening card that should appear only on Overview."
            }
        }
        $logsHTML = Get-Content -LiteralPath (Join-Path $reportRoot "logs.html") -Raw
        foreach ($pattern in @(
            "dnaprs-log-panel", "<details", "<iframe",
            "provenance/execution_report.html", "provenance/execution_timeline.html",
            "provenance/pipeline_dag.html", "provenance/execution_trace.txt", "data-log-source",
            "<pre>", "Report software"
        )) {
            if ($logsHTML -notmatch [regex]::Escape($pattern)) {
                throw "The Logs page is missing inline execution content: $pattern"
            }
        }
        $reportScript = Get-Content -LiteralPath `
            (Join-Path $reportRoot "assets/dnaprs-report.js") -Raw
        foreach ($pattern in @(
            "data-dnaprs-theme", "data-table-viewer", "data-table-expand",
            "data-figure-gallery", "data-figure-zoom-in", "dnaprs:figure-zoom",
            "data-log-source", "fetch(target.dataset.logSource"
        )) {
            if ($reportScript -notmatch [regex]::Escape($pattern)) {
                throw "The report script is missing required behaviour: $pattern"
            }
        }
        if ($reportScript -match "data-figure-index-button") {
            throw "The report script still contains the removed duplicate figure-name buttons."
        }
        foreach ($pageName in @(
            "genotype-eda.html", "target-prep.html", "target-qc.html", "target-imputation.html",
            "gwas-qc.html", "plink-prs.html", "sbayesrc-prs.html", "phenotype.html"
        )) {
            $galleryHTML = Get-Content -LiteralPath (Join-Path $reportRoot $pageName) -Raw
            foreach ($pattern in @(
                "data-figure-gallery", "data-figure-stage", "data-figure-zoom-reset",
                "Download SVG", "Download TIFF", "Download PNG", "Download JPEG"
            )) {
                if ($galleryHTML -notmatch [regex]::Escape($pattern)) {
                    throw "$pageName is missing figure-gallery content: $pattern"
                }
            }
            if ($galleryHTML -match "data-figure-index") {
                throw "$pageName still contains the removed duplicate figure-name row."
            }
        }

        $figureDictionary = Import-Csv `
            (Join-Path $reportRoot "provenance/figure_manifest.tsv") -Delimiter "`t"
        if ($figureDictionary.Count -ne $figureNames.Count) {
            throw "The figure dictionary does not contain all $($figureNames.Count) expected figures."
        }
        if (($figureDictionary | Where-Object {
            -not $_.inspection -or -not $_.svg -or -not $_.tiff -or -not $_.png -or -not $_.jpeg -or -not $_.source_table
        }).Count -gt 0) {
            throw "The figure dictionary contains incomplete inspection, image, or source-table fields."
        }

        if (-not $removeTemporary) {
            $finalReport = Join-Path $temporary "reports"
            if (Test-Path -LiteralPath $finalReport) {
                Remove-Item -LiteralPath $finalReport -Recurse -Force
            }
            Move-Item -LiteralPath $reportRoot -Destination $finalReport
        }
    }
    finally {
        Pop-Location
    }

    Write-Host "dnaprs R and report smoke tests passed."
    if (-not $removeTemporary) { Write-Host "Smoke-test report: $(Join-Path $temporary 'reports')" }
}
finally {
    if ($removeTemporary -and (Test-Path -LiteralPath $temporary)) {
        Remove-Item -LiteralPath $temporary -Recurse -Force
    }
    elseif (Test-Path -LiteralPath $working) {
        Remove-Item -LiteralPath $working -Recurse -Force
    }
}
