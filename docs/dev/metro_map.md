# Workflow metro map

The README workflow figure is generated from `assets/metro_map.mmd` with
[`nf-metro`](https://github.com/pinin4fjords/nf-metro), following the same source-first
approach used by nf-core pipelines. Update the Mermaid source whenever stages or data
connections change, then regenerate the SVG:

```bash
pip install 'nf-metro>=1.1.0' cairosvg
nf-metro render assets/metro_map.mmd \
  -o docs/images/dnaprs-workflow.svg \
  --theme nfcore-light --x-spacing 90 --y-spacing 78 \
  --section-x-gap 120 --diamond-style symmetric --directional --responsive \
  --no-chrome-css --no-strict
```

Inspect the rendered SVG at desktop and narrow README widths. Labels must not overlap
stations, section headings, or other labels, and every conditional route must remain
distinguishable without relying on animation. `nf-metro` 1.1 may warn when a routed
line crosses a station marker; treat that warning as a manual-review prompt and do not
accept any label or station overlap in the rendered SVG.
