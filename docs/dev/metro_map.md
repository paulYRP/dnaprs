# Workflow metro map

`assets/metro_map.mmd` is the semantic source for the routes and dependencies shown in
the README. The nf-core/rnaseq reference is generated from the same type of Mermaid
source with `nf-metro`; it is not a separate nf-core pipeline-template image. The dnaprs
production SVG therefore reuses the reference's visual tokens: a white canvas,
`#ededed` section panels, `#333333` labels, narrow rounded station markers, four-pixel
rounded coloured routes, parallel route bundles, a translucent white logo-and-legend
panel, and animated white route markers. The complete six-colour rnaseq palette is used
for the dnaprs route families, including a separate orange direct-score sensitivity
branch and red phenotype-model branch.

`docs/images/dnaprs-workflow.svg` retains a hand-tuned 1795 by 930 layout because the
automatic renderer compressed the larger dnaprs workflow into an unreadable single row.
This preserves the nf-core/rnaseq appearance while making the PRS-specific branches and
their dependencies legible in the README.

The diagram includes the generated `nf-core-dnaprs` pipeline logo. Keep both repository
logo variants under `docs/images/`; the SVG embeds the light-canvas variant so that the
logo also remains visible when GitHub displays the SVG as an image.

Update the Mermaid source and SVG together whenever a stage or dependency changes. Keep
station labels short; put explanations in the smaller secondary labels. The production
SVG must retain its accessible title and description, the reference's neutral palette,
and the six route colours defined in the Mermaid source.

Validate source edits from the repository root with:

```bash
nf-metro validate assets/metro_map.mmd
```

The source must remain semantically valid and keep the same five section names as the
production SVG. Do not replace the production SVG with an automatic render without
reviewing it at README size: dnaprs has cross-method target, GWAS and reference
dependencies that require the hand-tuned routing to remain legible.

Inspect the SVG at its native size and inside the rendered README. No label may overlap a
station, route, section heading, or another label. Optional paths must remain identifiable
without animation, and the map must still make sense when SVG animation is disabled.
