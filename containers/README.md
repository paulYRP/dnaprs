# dnaprs container

The local Docker profile uses one image so a process never receives conflicting software
environments when it needs both PLINK and R. The image extends the authors' versioned
`zhiliz/sbayesrc:0.2.6` image, which includes SBayesRC, PLINK, and its compiled
dependencies, and adds the fixed Quarto release and R packages used for modelling and
reporting. It also installs the system libraries required to compile those R packages.
Direct R package versions are declared in `containers/R-packages.tsv`; the build stops
if the final image does not contain those versions.

The SBayesRC base image starts with its own `sbayesrc` entry point. dnaprs clears that
entry point so Nextflow can execute each process command normally. The image retains
`/bin/bash` as its default command.

Build and record the immutable image identifier:

```bash
docker build --pull --tag dnaprs:1.0.0 --file containers/Dockerfile .
docker image inspect dnaprs:1.0.0 \
  --format 'image={{.Id}} entrypoint={{json .Config.Entrypoint}} cmd={{json .Config.Cmd}}'
docker run --rm --volume "$PWD:/work" dnaprs:1.0.0 Rscript -e \
  "for (p in read.delim('containers/R-packages.tsv')\$package) cat(p, as.character(packageVersion(p)), '\n')"
docker run --rm dnaprs:1.0.0 quarto --version
docker run --rm dnaprs:1.0.0 plink2 --version
```

The Docker profile uses the local tag `dnaprs:1.0.0`. The image has not been published
to a registry. Record its image ID for local runs. Before a public production release,
publish the tested image and replace the profile tag with its immutable registry digest.
The pipeline checks that use this image are in the
[contributing guide](../docs/CONTRIBUTING.md).
