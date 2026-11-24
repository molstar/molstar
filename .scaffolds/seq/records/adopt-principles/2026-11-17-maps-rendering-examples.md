# adopt-principles: maps, rendering, examples and edge cases (2026-11-17)

## Asked

- Whether maps could be computed or compressed rather than spelled out position by position; then "yes, fix this".
- "few example datasets, but … more views", keeping "the richness, the editing experiences", and watching for edge
  cases "so the 2nd implmentation doesn't smooth over them".
- Keep "the render ports, the ability render to svg", since the JBrowse reply's author rendered Seq\* documents in
  JBrowse beside Seq\*'s own SVG.

## Decided

- Maps: written per kind in its shortest form (identity, offset, numbering, CIGAR, exons, gapped rows, segments,
  computed, or a source); resolved always to segment tables ([design/building-blocks.md](../../design/building-blocks.md)).
- Rendering: the spec, the resolved view and draw lists as public forms; backends as the render port; SVG and the
  recording backend in the core, so a spec renders headless; a text-metrics port; a render CLI
  ([design/rendering.md](../../design/rendering.md)).
- Examples: six datasets, many views each, editing views included, every view in Storybook, a test and as portable
  files ([design/examples.md](../../design/examples.md)).
- Edge cases: the previous implementation's 184 items catalogued; each is kept, changed or dropped by decision. The
  proposed changes and flaw fixes await the author ([design/edge-cases.md](../../design/edge-cases.md)).

## Delivered

Uncommitted on `packages-seq`: the designs above, the catalog in `sources/seqstar/edge-cases.md`, the core and viewer
specs (maps, rendering forms), the TODO, and the Storybook pages Building blocks, Examples and Rendering. Checked with
the three scaffold checks and `pnpm format:check`.

## Left open

The confirmation of the proposed edge-case outcomes; the root nouns; the datasets' exact files.
