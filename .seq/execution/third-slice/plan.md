# Plan: the third slice (draft)

**A draft for direction, written 2026-11-18 while the first slice closes.** It is rewritten as the second slice closes,
from what that slice teaches. No loop runs on it until then.

Work id `third-slice`. Version 1's second milestone, **figures** ([version-1](../../design/version-1.md), Milestones
after the slice): the layouts beyond tracks, the remaining layers and themes, PNG and the render CLI, and the four
other datasets with their views ([examples](../../design/examples.md)).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first (closing) | the end-to-end slice | the tree core, the language, resolve, draw lists, Canvas2D and SVG, the viewer's models, data by region, measured |
| second (planned) | 1, core editing and addressing | patches, Seq\*'s history port, previews, snapshots, addresses, refs and links, selection, security bounds, MessagePack |
| **third (this draft)** | 2, figures | the text, ring and block layouts, alignments, dot plots, arcs and logos, themes and legends, PNG and the CLI, four datasets |
| fourth | 3, integration | structures (the port, structure maps, the Mol\* showcase), services, bigWig and range reads, the Nightingale reference, archives |
| later | 4, Genomes | layers across sections, synteny, contact maps, MAF, two assemblies on one axis, the genome formats |

Version 1 is the end of the fourth slice: every planned row in [features](../../design/features.md) built and shown by
a view.

## Done when

- Every dataset but Genomes has its views from examples.md running three ways (Storybook with canvas and SVG side by
  side, a quick-lane test, a spec file a reader without Seq\* could take):
  - **p53**: the text figure, secondary structure and variants on the tracks, the themes view;
  - **Globins**: the family alignment with logo and rows by tree, collapsed clades, the ESPript figure, hemoglobin with
    copies, the 10,000-row alignment with lazy rows and the keyboard cursor;
  - **Antibody**: Kabat against IMGT, the design against the wild type, CDR bookmarks with undo;
  - **Plasmid**: the ring with a linear detail following its window, a feature across the origin;
  - **tRNA**: arcs nested and crossing, the base-pair dot plot with crosshair hover and rectangle selection.

  Each dataset has an editing view (from the second slice's editing and history).
- The section layout interface carries four layouts (tracks, text, blocks, rings), and the dot plot as its own
  section, with no layout-specific code in the viewer.
- PNG renders headless with the same picture as SVG. The render CLI turns a spec file into SVG or PNG under the
  headless asset policy (EC-126).
- Themes: field palettes (binned, continuous, log), rules, the residue schemes, legends, shared domains.
- Accessibility beyond the chrome: titles in exported SVG, the cursor and hover announced, forced colors (version-1's
  open item for this milestone).
- Measured: the 10,000 × 300 alignment, a ring's frame, the text layout's reflow, PNG export, the bundles. The author
  has used the views (§used-before-done).

## Gate 0 (to decide when the slice is rewritten)

- **The datasets' exact data**: which globins, which antibody, which plasmid and tRNA, which p53 structure.
  examples.md leaves this provisional. Recommendation: small public entries with permissive licenses, fixtures under
  1 MB each, the large variants generated in code.
- **PNG in Node**: rasterize the SVG (resvg as a dev or optional dependency) or a Node canvas. Recommendation: rasterize
  the SVG, so there is one picture and no native canvas in the library.
- **The render CLI's home**: `cli/` (the repository's command-line packages) as `seqx-render`. Recommendation: yes.
- **The theme sublanguage**: the `layer.color` slot's grammar (fields, palettes, rules, schemes) and legends' place in
  the chrome. Its design is drafted first, in the slice's first stage.
- **Glyph atlas and band copy** in the Canvas2D backend: in this slice only if the second slice's browser measurements
  ask for it.

## Stages (a sketch)

1. **Themes and the small layers** (strong + fast):
   - the theme sublanguage, residue schemes, legends, shared domains, size channels;
   - the ruler and numbering, secondary structure, variants, values' line and area, ranges' line styles and arrowheads.
2. **Layouts** (two strong lanes):
   - the text layout and alignment figures (groups of ten, lanes, inline labels, numbers per row, frames);
   - blocks and rings (labels with leaders, a window a linear detail follows, EC-108 and EC-114), navigators.
3. **Pairs, grids and rows** (strong):
   - arcs and chords;
   - the grid and dot plot as a section of its own, with crosshair hover and rectangle selection;
   - the logo;
   - the rows layer (tree or data order, collapsed clades, lazy rows for 10,000).
4. **Datasets and views** (fast lanes, one per dataset): fixtures, builders, views, stories and their headless twins,
   and an editing view per dataset.
5. **Output and accessibility** (fast): PNG, the render CLI, SVG titles, announcements, forced colors.
6. **Measure**, then finalize.

## Known risks to watch from here

- **The layout interface** has carried one layout (tracks). The text layout and rings are its real test; a change to
  it touches the draw lists' key and the hover index.
- **Lazy rows for 10,000 × 300**: the first slice never laid out rows on demand. The store's want sets and the layout's
  bounded allocation meet here.
- **The config surface** grows again with each layout: the Configuration tables and the gaps queue's "unsure only"
  rule (the first slice's fold-back) are tested at this size.
- **Seq\* inside a bigger app**: history is already a port (the second slice). The third slice keeps every new piece
  host-neutral as well: legends and navigators as parts a host's chrome may place, not a frame Seq\* owns.

## Out of scope

Structures and services (the fourth slice), bigWig, Genomes, writers beyond SVG and PNG, settings layers.
