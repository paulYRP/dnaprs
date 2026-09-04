# dnaprs software environments

dnaprs assigns a versioned container to each process. Reference panels, GWAS files,
target data, and phenotype data are not stored in these images.

| Image                                                 | Processes                                                                        | Fixed software                                                                                            |
| ----------------------------------------------------- | -------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `dnaprs-analysis:1.0.0`                               | Input validation, GWAS harmonisation, score processing, QC, and phenotype models | R 4.4.1 and the R packages in `analysis/R-packages.tsv`                                                   |
| `quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0` | Target QC, reference frequency, clumping, and C+T scoring                        | PLINK 2.0 alpha 6.9                                                                                       |
| `dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21`           | Raw genotype EDA and identity by descent                                         | PLINK 1.90b6.21 and PLINK 2.0 alpha 6.12                                                                  |
| `dnaprs-imputation:1.1.0`                             | Target marker resolution, reference preparation, QC, and imputation              | Beagle/unbref3 27Feb25.75f, bcftools/tabix, Java 17, PLINK 2.0 alpha 6.12, R 4.4.1, and data.table 1.18.0 |
| `zhiliz/sbayesrc:0.2.6`                               | SBayesRC preparation, modelling, and scoring                                     | The authors' SBayesRC 0.2.6 environment                                                                   |
| `dnaprs-report:1.0.0`                                 | Portable report generation                                                       | The analysis environment, Quarto 1.7.32, and the packages in `report/R-packages.tsv`                      |

The PLINK 2-only processes use the matching versioned BioContainer and native
Singularity image. `GENOTYPE_EDA` keeps the combined dnaprs image because it runs
PLINK 2 and then PLINK 1.9 `--genome full` on the same intermediate files. SBayesRC
processes use the authors' versioned 0.2.6 image.

Nextflow process declarations use stable version tags and do not append a Docker
digest. This avoids the `tag@sha256` reference that the Singularity version on Aqua
cannot parse. The dnaprs Dockerfiles still use digest-pinned base images. Downloaded
installers are SHA-256 checked, and R dependencies resolve against the dated
`2026-08-01` Posit Package Manager CRAN snapshot.

The component audit retained the local processes because they produce
pipeline-specific decision tables, QC records, or checkpoints that the generic
nf-core modules do not provide. Their PLINK 2 commands use a standard BioContainer
where one tool is sufficient. The SBayesRC image is used only by SBayesRC processes.

## Local build

Build the images from the repository root in dependency order:

```bash
docker build \
  --file containers/analysis/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-analysis:1.0.0 \
  .

docker build \
  --file containers/plink2/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21 \
  .

docker build \
  --file containers/report/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-report:1.0.0 \
  .

docker build \
  --file containers/imputation/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-imputation:1.1.0 \
  .
```

Verify the environments before a pipeline test:

```bash
docker run --rm ghcr.io/paulyrp/dnaprs-analysis:1.0.0 \
  Rscript -e "stopifnot(as.character(getRversion()) == '4.4.1')"
docker run --rm quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0 \
  plink2 --version
docker run --rm ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21 plink2 --version
docker run --rm ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21 plink --version
docker run --rm --entrypoint Rscript docker.io/zhiliz/sbayesrc:0.2.6 \
  -e "stopifnot(as.character(packageVersion('SBayesRC')) == '0.2.6')"
docker run --rm ghcr.io/paulyrp/dnaprs-report:1.0.0 quarto --version
docker run --rm ghcr.io/paulyrp/dnaprs-imputation:1.1.0 \
  bash -c "bcftools --version | head -n 1; java -jar /opt/beagle/beagle.jar 2>&1 | head -n 2; java -jar /opt/beagle/unbref3.jar help 2>&1 | head -n 2; Rscript -e \"stopifnot(as.character(packageVersion('data.table')) == '1.18.0')\""
```

`.github/workflows/containers.yml` builds each dnaprs image once and verifies that
exact local image before any registry login. It also checks the external PLINK 2 and
SBayesRC images. Images are pushed only by a matching stable release or a manually
dispatched publication. Pull requests run the build and tests, but their publication
conditions are false, so they do not log in or push. Published dnaprs packages must be
publicly pullable before an HPC run is accepted.

The Docker profile clears container entry points before starting each task. This is
required for the authors' SBayesRC image, whose default entry point is its command-line
wrapper rather than a Nextflow task shell.

## HPC reuse

PLINK 2-only processes use the native Singularity image on Aqua. Nextflow converts and
caches the other versioned Docker images as needed. An HPC user does not install PLINK
or R packages inside a personal library. Reference PGEN companions, SBayesRC LD
directories, and annotation files are declared as Nextflow `path` inputs and staged as
read-only task links. A shared Apptainer cache can be configured in an external
institutional configuration.
