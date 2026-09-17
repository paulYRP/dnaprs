# Workflow metro map

The diagram is generated with [nf-metro](https://seqeralabs.github.io/nf-metro/)
from `assets/metro_map.mmd`. Edit this source and regenerate the images; do not
adjust the SVG geometry by hand. nf-metro is a documentation tool, not a pipeline
runtime dependency.

## Layout

The map uses the nf-core/rnaseq light style: grey section panels, numbered headings,
coloured routes, white stations and a logo beside the legend. Five sections fold
across two rows. Layout settings and labels are stored in the Mermaid source.

Keep labels short and check connections against `workflows/dnaprs.nf`:

- Raw genotype EDA is a diagnostic branch, not an input to target preparation.
- REF icons identify independent consumers of the verified reference cache.
- Beagle imputes target genotypes; SBayesRC imputes GWAS summary statistics.
- Eligibility applies to primary scoring and direct-score sensitivity.
- Direct scores use unimputed genotypes and do not enter the primary score table.
- Primary scores reach the outputs without requiring phenotype models.

The map summarises scientific stages, not every process or auxiliary input.
Reference frequency and annotation inputs support their downstream scoring and
model stages. The three hidden stations in the phenotype section are routing
anchors, not analysis steps. They keep the shared output route straight while
the optional phenotype branch passes through its models.

## Regenerate

Use the tested renderer version in a development Python environment. From the
repository root:

```bash
python -m pip install nf-metro==2.1.0

nf-metro validate assets/metro_map.mmd

nf-metro render assets/metro_map.mmd \
  --strict --validate \
  -o docs/images/dnaprs-workflow-static.svg

nf-metro render assets/metro_map.mmd \
  --strict --validate --animate \
  -o docs/images/dnaprs-workflow.svg

nf-metro render assets/metro_map.mmd \
  --strict --validate --mode light --raster-width 2265 \
  -o docs/images/dnaprs-workflow.png

prek run --files assets/metro_map.mmd docs/dev/metro_map.md README.md \
  docs/images/dnaprs-workflow.svg docs/images/dnaprs-workflow-static.svg \
  docs/images/dnaprs-workflow.png
```

The README displays the animated SVG and links to the static PNG. The static SVG
keeps selectable labels and supports enlargement without losing image quality.
All images use the same source; no intermediate renders belong in the repository.

## Review

The strict layout and render checks must pass without warnings. Inspect the map
at full size and at README width on light and dark backgrounds. Check for crossed
headings, overlapping labels, clipped text, misleading junctions and crowded
reference inputs. Confirm that the static map remains clear without animation.

The light theme has a transparent background around its grey panels. This area
follows the page background; it does not represent another analysis stage.
