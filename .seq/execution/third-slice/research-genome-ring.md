# Research: a genome ring like Rondo's (GR, 2026-11-19)

Lane GR's read-only report (`briefs/gr-genome-ring-research.md`), condensed, with the orchestrator's checks. Input to
the fifth slice's plan (Genomes) and the author's questions at the end.

## Rondo

CSIRO Data61 with the Garvan Institute (Sabir et al., VIZBI 2015 poster); "inspired by Circos"; chromosome
connectivity from Hi-C, ChIA-PET and 5C (default: 38,482 conserved interactions over five human cell types). Its pages
never say "edge bundling": connections are **clustered hierarchically, precomputed at several resolutions**, internal
connections of a cluster shown only zoomed in; so the trumpets are cluster stems and the ribbons cluster-to-cluster
links, thick by count (GR's inference). Gene detail from about 1 Mb; filters by region, Q-value, minimum connections;
export as SVG, PNG, PDF, BEDPE, BED, a matrix. Nothing public on the square navigator or the shapes.

## By layer

**The spec**

| Element | Status | Note |
|---|---|---|
| chromosomes as blocks on a ring, gaps | **built** | blocks of several sequences on one ring: `rings-fixtures.ts:126` `blocksSpec({ layout: 'rings' })`, tested (orchestrator; GR thought it waited for Genomes) |
| cytobands per block | built | `bands` layer, `giemsa`, `acen` triangles (`language.ts:790-800`) |
| links | built | `pairs` on one or two systems |
| a hue per block | close | no theme reads the block; a field on the rows or a block measure |
| bundling | missing | no reduce form; a layer param plus a derived bin-pair table, or a `pairs` level read by region (GR leans to the level) |
| pairs read by region | missing | grids-and-rows.md: waits for Genomes (EC-208) |
| chromosome 4 as its own ring | close | a second ring section over one system; a navigation step, not a gesture |
| outward links to other chromosomes | missing | an end not shown is a caret today; needs a stub colored by the other block |
| the square navigator | missing | navigators move a window on one axis; navigators.md refuses switching systems |
| an extension's layer in a spec | missing | the layer kind is a fixed enum (`language.ts:790-800`); `ext:` kinds designed (extensibility.md, EC-235), not built |

**Rendering**

| Element | Status | Note |
|---|---|---|
| a ring of 23 blocks, bands, labels | close | built for a few blocks; 23 untested |
| chords between blocks on a ring | close | chords on a ring and arcs between blocks each tested (`pairs-grids.test.ts:636, 655`); together untested |
| 12,000+ links | close | `pairs.maxArcs` 4,096 (`config.ts:332`), culled by device column plus a density strip: a screen cull, not a count by bin pair |
| ribbons wide by count, trumpets, bulbs | missing | a new layer type (or a `pairs` mode); the polar warp already bends paths |
| a pyramid of bin-pair counts | missing | grids have a 2D pyramid; pairs none |
| the zoom into one chromosome's ring | close | a second section; the transition missing |
| the square navigator | missing | a navigator kind over systems |

**A story**

| Step | Status | Note |
|---|---|---|
| the whole ring | yes | a state |
| one chromosome hovered | partly | hover is session emphasis; a step names an emphasis (a query "pairs touching block N") |
| zoom into chromosome 4's ring | a cut | changing a section's `of`/`range` is a cut (stories.md); a ring-zoom channel missing |
| back | yes | a reverse step |

## As an extension

Layer types are an open registry by design, and the renderer takes extra types (`renderer.ts:41, 359`), but a spec
cannot name one (no `ext:` kinds), the viewer's passing of `types` is unverified, layouts and navigators are internal,
and no pairs reduction is public. GR's estimate: a ribbon-and-fan layer 400–700 lines, a generator and bin-pair table
150, a builder 100. The navigator and the ring zoom cannot be an extension today.

## GR's recommendation

1. Smallest honest demo: 23 chromosomes as blocks on a ring, bands per block, a hue per block, generated links as
   today's chords (culled at 4,096), chromosome 4 as a second ring. Circos-like, not Rondo-like, and said so.
2. Rondo-like: a bundled-ribbons layer as an extension package, testing the layer contract at genome scale and
   needing `ext:` kinds first.
3. The fifth slice, with Genomes.

What it teaches the designs: a pairs level by bin pair; `ext:` layer kinds; an emphasis "everything touching this
block"; a navigator over systems; a ring-zoom story channel.

## Open, for the author

1. Circos-like (today's chords, culled) or Rondo-like (bundled, count-weighted)?
2. Bundling as a `pairs` level read by region (recommended) or a layer param over a derived table?
3. `ext:` layer kinds before the demo, so it is a true extension?
4. A navigator over several systems: against navigators.md's refusal, or the point of Genomes?
5. The ring zoom as an animated story channel, or a cut to a second section?
6. UCSC cytobands committed (licence to check) or synthetic bands over GRCh38 lengths?
