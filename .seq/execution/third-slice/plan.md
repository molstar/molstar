# Plan: the third slice

**Rewritten 2026-11-18 as the second slice closed**, from what that slice taught (its record:
[records/second-slice](../../records/second-slice/2026-11-18-second-slice.md)). It is a plan for the author to agree at
Gate 0; no loop runs on it before then, and its loop file is written after.

Work id `third-slice`. Version 1's second milestone, **figures** ([version-1](../../design/version-1.md)): the axis
finished (direction, letters through a map, collapsed introns), the layouts beyond tracks, the remaining layers and
themes, the stories library and its animation port, navigators, PNG and the render CLI, and the other datasets with
their views ([examples](../../design/examples.md)).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first (closed) | the end-to-end slice | the tree core, the language, resolve, draw lists, Canvas2D and SVG, the viewer's models, data by region, measured |
| second (closed) | 1, core editing and addressing | patches, the history port, previews, addresses, refs and links, selection, security bounds, MessagePack, navigation |
| **third (this plan)** | 2, figures | the axis, the text, ring and block layouts, alignments, dot plots, arcs and logos, themes and legends, stories, navigators, PNG and the CLI, the datasets |
| fourth | 3, integration | structures (the port, structure maps, the Mol\* showcase), services, bigWig and range reads, the Nightingale reference, archives |
| later | 4, Genomes | layers across sections, synteny, contact maps, MAF, two assemblies on one axis, the genome formats |

## What the second slice changes in how this one runs

From the principles ("Folded back from the second slice"), written into this plan's gates and its loop file:
- **Each new design gets the self-check when it is agreed** (the axis, the layouts, the theme sublanguage, navigators,
  the stories library): its premise negated, what integrates planned, the rest kept as a tension in the design.
- **The interaction matrix** ([ux-foundations](../../design/ux-foundations.md), section 8): this slice fills the cells
  of a track or layer (select a row, focus a track), the selection's bookmark command, a bookmark's removal and label,
  rows and clades (collapse, pin), the dot plot's rectangle, the ring's window; Gate 0 names them, each gate checks
  them.
- **Precedents cover interactions** at Gate 0: what Jalview, ESPript, IGV, SnapGene-style plasmid viewers, R2DT/forna
  and the current Seq\* let a user do with these figures, each listed as built, planned or refused.
- **Use gates**: a short one after stage 1 (the axis and the first themes, on the p53 views), and a novice pass by a
  light agent before the author's use gate at stage 6. Each gate's findings are counted by kind (bug, missing design,
  look); the fold-back reads the count (the carried question: does an early use gate pay for figures?).
- **At every gate**: the size check (`pnpm -C packages/seq/core size`) beside the tests; the gaps raised since the last
  gate asked as yes or no; drawing checked against the data in the example views; heavy runs never two at once.

## Done when

- Every dataset has its views from examples.md running three ways (Storybook with canvas and SVG side by side, a
  quick-lane test that also checks drawing against the data, a spec file a reader without Seq\* could take):
  - **p53**: the protein and its gene on one axis (the protein's letters through the codon map), the gene following
    its protein reversed, introns collapsed; the text figure, secondary structure and variants on the tracks, the
    themes view;
  - **Globins**: the family alignment with logo and rows by tree, collapsed clades, the ESPript figure, hemoglobin with
    copies, the 10,000-row alignment with lazy rows and the keyboard cursor;
  - **Antibody**: Kabat against IMGT, the design against the wild type, CDR bookmarks with undo;
  - **Plasmid**: the ring with a linear detail following its window, a feature across the origin;
  - **tRNA**: arcs nested and crossing, the base-pair dot plot with crosshair hover and rectangle selection;
  - **Complex** (CDK2–cyclin A with an AlphaFill ligand): the two chains with their interface and the ligand's
    contacts, as examples.md describes.

  Each dataset has an editing view.
- The section layout interface carries four layouts (tracks, text, blocks, rings) and the dot plot as its own section,
  with no layout-specific code in the viewer.
- A section's direction (`genomic`, `reversed`, `follow`), letters through a one-to-one or codon path, collapsed
  introns: built on the axis's pieces, hit-testing, selection and export following.
- **Stories**: `stateOf`, `plan` and `sample` in `seqx-core`, the `AnimationPort` with its default and headless
  adapters, the p53 story played, scrubbed and exported ([stories](../../design/stories.md)).
- **Navigators**: designed (`navigators.md`, with its self-check) and built: a strip and a miniature that move the view.
- Themes: field palettes (binned, continuous, log), rules, the residue schemes, legends, shared domains.
- PNG renders headless with the same picture as SVG; the render CLI turns a spec file into SVG or PNG under the headless
  asset policy (EC-126).
- Accessibility beyond the chrome: titles in exported SVG, the cursor and hover announced, forced colors.
- Measured: the 10,000 × 300 alignment, a ring's frame, the text layout's reflow, PNG export, a story's frame, the
  bundles at every gate. The author has used the views at both use gates.

## Gate 0 (the author, before any lane)

- **The datasets' exact data**: which globins, antibody, plasmid, tRNA and p53 structure; the Complex's entries.
  Recommendation: small public entries with permissive licenses, fixtures under 1 MB each, the large variants generated
  in code.
- **PNG in Node**: rasterize the SVG (resvg as a dev or optional dependency) or a Node canvas. Recommendation: rasterize
  the SVG (one picture, no native canvas in the library).
- **The render CLI's home**: `cli/` as `seqx-render`. Recommendation: yes.
- **Designs drafted first, each with its self-check**: the axis (direction, letters through a map, collapsed introns),
  the theme sublanguage (the `layer.color` slot's grammar, legends' place), `navigators.md`, the layout interface's
  extension for text, blocks and rings. Each agreed before its lane.
- **The interaction matrix's cells this slice fills**, and the precedent capability list.
- **Glyph atlas and band copy**: not needed by the second slice's browser numbers (frames about 1 ms). Recommendation:
  out of this slice unless the alignment's measurement asks.
- **Size targets with room**: headless rendering has 32 kB left under 300 kB, and this slice adds layouts, themes and
  layers to it. Recommendation: measure at stage 1's gate and bring a target change, if any, then.

## Stages

0. **Carry-overs** (light): exporting `projectPoint`, the codon-split follow edges, a link followed before load kept
   waiting, the Features links story's double message, `track.add`, the history's toolbar commands, drawing checked
   against the data in the existing example views (`packages/seq/TODO.md`).
1. **The axis, then themes and the small layers** (strong + fast): the axis lane first (layouts build on it); the theme
   sublanguage, residue schemes, legends, shared domains, size channels; the ruler and numbering, secondary structure,
   variants, values' line and area, ranges' line styles and arrowheads. **Use gate**: the author tries the p53 views for
   ten minutes.
2. **Layouts** (two strong lanes): the text layout and alignment figures (groups of ten, lanes, inline labels, numbers
   per row, frames); blocks and rings (labels with leaders, a window a linear detail follows, EC-108 and EC-114).
3. **Pairs, grids and rows** (strong): arcs and chords; the grid and dot plot as a section of its own (crosshair hover,
   rectangle selection); the logo; the rows layer (tree or data order, collapsed clades, lazy rows for 10,000).
4. **Stories and navigators** (strong + fast): the stories library and the `AnimationPort`; navigators.
5. **Datasets and views** (fast lanes, one per dataset): fixtures, builders, views, stories and their headless twins,
   an editing view per dataset.
6. **Output and accessibility** (fast): PNG, the render CLI, SVG titles, announcements, forced colors. **Novice pass**
   (a light agent with a task script per view), then **the author's use gate**.
7. **Measure**, then finalize (reconcile, the fold-back as its own step, the record, the next plan).

## Known risks

- **The layout interface** has carried one layout (tracks); the text layout and rings are its real test, and a change
  to it touches the draw lists' key and the hover index.
- **The axis lane first**: direction and collapsed introns change what an offset means on a section; hit-testing,
  selection, follow and export all read it. A twin per consumer, and drawing checked against the data.
- **Lazy rows for 10,000 × 300**: the store's want sets and the layout's bounded allocation meet here.
- **Headless rendering's size**: 32 kB of room for layouts, themes and layers (Gate 0).
- **Seq\* inside a bigger app**: legends, navigators and the animation loop stay host-neutral (parts a host places, a
  port the host drives).

## Out of scope

Structures and services (the fourth slice), bigWig, Genomes, writers beyond SVG and PNG, settings layers.
