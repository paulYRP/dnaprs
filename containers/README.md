# dnaprs software environments

dnaprs assigns a versioned container to each process. Reference panels, GWAS files,
target data, and phenotype data are not stored in these images.

| Image                        | Processes                                                                        | Fixed software                                                                       |
| ---------------------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `dnaprs-analysis:1.0.0`      | Input validation, GWAS harmonisation, score processing, QC, and phenotype models | R 4.4.1 and the R packages in `analysis/R-packages.tsv`                              |
| `dnaprs-plink2:2.0.0-a.6.12` | Target conversion and PLINK C+T operations                                       | PLINK 2.0 alpha 6.12                                                                 |
| `dnaprs-imputation:1.1.0`     | Source-reference preparation and target-genotype imputation                       | Beagle/unbref3 27Feb25.75f, bcftools/tabix, Java 17, unzip, and PLINK 2.0 alpha 6.12 |
| `zhiliz/sbayesrc:0.2.6`      | SBayesRC preparation, modelling, and scoring                                     | The authors' SBayesRC 0.2.6 environment                                              |
| `dnaprs-report:1.0.0`        | Portable report generation                                                       | The analysis environment, Quarto 1.7.32, and the packages in `report/R-packages.tsv` |

The SBayesRC image is pinned to the authors' published image digest in each SBayesRC
module. The four dnaprs images are built from the Dockerfiles in this directory and
published to GitHub Container Registry only for a matching stable release.

## Local build

Build the images from the repository root in dependency order:

```bash
docker build \
  --file containers/analysis/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-analysis:1.0.0 \
  .

docker build \
  --file containers/plink2/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12 \
  .

docker build \
  --build-arg ANALYSIS_IMAGE=ghcr.io/paulyrp/dnaprs-analysis:1.0.0 \
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
docker run --rm ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12 plink2 --version
docker run --rm --entrypoint Rscript docker.io/zhiliz/sbayesrc:0.2.6 \
  -e "stopifnot(as.character(packageVersion('SBayesRC')) == '0.2.6')"
docker run --rm ghcr.io/paulyrp/dnaprs-report:1.0.0 quarto --version
docker run --rm ghcr.io/paulyrp/dnaprs-imputation:1.1.0 \
  bash -c "bcftools --version | head -n 1; java -jar /opt/beagle/beagle.jar 2>&1 | head -n 2; java -jar /opt/beagle/unbref3.jar help 2>&1 | head -n 2"
```

`.github/workflows/containers.yml` repeats these builds and checks. Images are pushed
only when the corresponding stable pipeline release is published, preventing a fixed
version tag from being silently replaced during development.

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
