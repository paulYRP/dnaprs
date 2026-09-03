# dnaprs software environments

dnaprs assigns a versioned container to each process. Reference panels, GWAS files,
target data, and phenotype data are not stored in these images.

| Image                                       | Processes                                                                        | Fixed software                                                                                            |
| ------------------------------------------- | -------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `dnaprs-analysis:1.0.0`                     | Input validation, GWAS harmonisation, score processing, QC, and phenotype models | R 4.4.1 and the R packages in `analysis/R-packages.tsv`                                                   |
| `dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21` | Target QC, identity by descent, and PLINK C+T operations                         | PLINK 1.90b6.21 and PLINK 2.0 alpha 6.12                                                                  |
| `dnaprs-imputation:1.1.0`                   | Target marker resolution, reference preparation, QC, and imputation              | Beagle/unbref3 27Feb25.75f, bcftools/tabix, Java 17, PLINK 2.0 alpha 6.12, R 4.4.1, and data.table 1.18.0 |
| `zhiliz/sbayesrc:0.2.6`                     | SBayesRC preparation, modelling, and scoring                                     | The authors' SBayesRC 0.2.6 environment                                                                   |
| `dnaprs-report:1.0.0`                       | Portable report generation                                                       | The analysis environment, Quarto 1.7.32, and the packages in `report/R-packages.tsv`                      |

The SBayesRC image is pinned to the authors' published image digest in each SBayesRC
module. The four dnaprs images are built from digest-pinned base images. Downloaded
installers are SHA-256 checked, and R dependencies resolve against the dated
`2026-08-01` Posit Package Manager CRAN snapshot.

The component audit retained these local environments because the active dnaprs
processes produce pipeline-specific decision tables, QC records, or checkpoints that
the generic nf-core modules do not provide. The standard PLINK 2 BioContainer also
does not currently provide the exact alpha 6.12 build used and tested here. The
SBayesRC authors' image is used only by SBayesRC processes; it is no longer used as a
general R or PLINK environment.

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
docker run --rm ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21 plink2 --version
docker run --rm ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12-plink1.90b6.21 plink --version
docker run --rm --entrypoint Rscript docker.io/zhiliz/sbayesrc:0.2.6 \
  -e "stopifnot(as.character(packageVersion('SBayesRC')) == '0.2.6')"
docker run --rm ghcr.io/paulyrp/dnaprs-report:1.0.0 quarto --version
docker run --rm ghcr.io/paulyrp/dnaprs-imputation:1.1.0 \
  bash -c "bcftools --version | head -n 1; java -jar /opt/beagle/beagle.jar 2>&1 | head -n 2; java -jar /opt/beagle/unbref3.jar help 2>&1 | head -n 2; Rscript -e \"stopifnot(as.character(packageVersion('data.table')) == '1.18.0')\""
```

`.github/workflows/containers.yml` builds each image once and verifies that exact local
image before any registry login. Images are pushed only by a matching stable release
or a manually dispatched publication. Pull requests run the build and tests,
but their publication conditions are false, so they do not log in or push. After the
first publish, each package must be made public, tested by anonymous Docker and
Apptainer/Singularity pulls, and its remote digest recorded in the process definitions
before the HPC run is accepted.

The Docker profile clears container entry points before starting each task. This is
required for the authors' SBayesRC image, whose default entry point is its command-line
wrapper rather than a Nextflow task shell.

## HPC reuse

The same Docker image addresses are used by the Apptainer profile. Nextflow pulls and
caches them as Apptainer images; an HPC user does not install the programs or R packages
inside a personal library. Reference PGEN companions, SBayesRC LD directories, and
annotation files are declared as Nextflow `path` inputs and staged as read-only task
links; they are not embedded in an image or written by a process. A shared Apptainer
cache can be configured in an external institutional configuration.
