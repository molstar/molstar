# Status

Start here: what is built, what is open, what comes next, and where to read. The current state only; how each slice
went is in its record.

## Built

- **`@molstar/seqx-core`** (`packages/seq/core`): SeqViewSpec, the language at `0.11.0` (the `group` node, the
  `reference` line on values, `height: 'content'`), with its builder, the checks
  of a document and of a change, and the query language (by role, `filter`, `expand`, `linked`, set operations, Keys);
  the model (sections, axes with directions, numberings, collapsed introns and reversed pieces; units, molecules and
  atoms; assemblies, blocks and the genome axis) and resolve with its maps (offset, segments, structure, codon) and the liftover across sections;
  addresses, refs, links and story steps; the source protocol, the data store (tiles, levels, 2D tiles, read costs,
  cancellation, eviction, summaries) and the readers (FASTA, aligned FASTA, Stockholm, Newick, records, GFF3, GenBank,
  opt-in mmCIF, bigWig, BED, GFF/GTF, VCF runtimes); themes, scales, legends; the layouts (tracks, text, alignment,
  blocks, rings, navigators) and layers (letters, ranges, values, ruler, variants, pairs, grids, the logo, rows); draw
  lists, the SVG, PNG and recording backends, `renderSvg` and the `seqx-render` command line; structures (facts, units,
  colors by run).
- **`@molstar/seqx-viewer`** (`packages/seq/viewer`): `createSeqViewer` and a minimal viewer; the models (spec as
  commits with previews and history, view with scroll groups and follow, hover, selection, highlights, address, store);
  behaviors (hover with its contract, navigate, the keyboard cursor, headers with menus, the square, bridges across
  sections); the hover card and tooltips, the keys panel, accessibility; extensions (`SeqExtension`, `extend()`); the
  structure link; the Canvas2D backend and the default ports.
- **`@molstar/seqx-ext-rondo`** (`packages/seq/ext-rondo`): the genome ring (`section.center`): stems, ribbons, stubs,
  cytobands, the square as a window, its menus.
- **`@molstar/seqx-examples`** and **`@molstar/seqx-storybook`**: the example views as builder modules with tests (p53,
  molecules over 4HHB, 1EHZ, 2OCJ and others, the generated datasets, Genomes), the stories with their headless twins,
  the Features, Examples and Design pages, the showcases and the perf tier.
- **Shared, in `@molstar/core`**: `repr` (the tree core, its JSON and S-expression codecs, `repr/patch`), half-open
  range sets in `data/int`, CHAMP maps and `Vector` in `data/immutable`, PCG32, the palette module, the asset manager's
  options and range reads, standard components. MessagePack in `data/msgpack` and the tree's codec on it in `repr/codecs/msgpack`. **In `@molstar/io`**: BGZF and tabix
  (and the MessagePack paths re-exported).
  **In `@molstar/ui`**: the history adapters.
- **Built by the close-out** (2026-11-24): search (matchers in the core, a registry for an application's, `viewer.search`
  and its commands, two highlight channels, Mod+G and F3; the box and results list are the application's); track
  groups (the `group` node, folds as session state, headers, summary rows, `lanes: one`); values marks (points, log
  scale, reference lines with their values); the legend's rows with an eye per layer and a Layers menu with Show all
  and Hide all; the ring's gestures (Cmd/Ctrl+wheel cut, a double-click on a chromosome's label or arc, the stem
  menu); `track.height` and `height: 'content'`; patches on streamed nodes waiting, never dangling; drags, the
  frame's allocation and the stories' twins at their bounds. Rows moved: [design/features.md](design/features.md).
- **Measured** at the close-out's last gate: every frame target met but the 200-copy split patch and merge (carried,
  below); sizes 265.0 / 995.2 / 1,590.8 / 1,120.9 kB min (reading, headless, default, minimal) against 268 / 1,000 /
  1,592 / 1,124 ([performance](design/performance.md), Size).

Which features are built: [design/features.md](design/features.md); which edge cases have an outcome:
[design/edge-cases.md](design/edge-cases.md), Outcomes.

Slices, each closed with its record:

- The scaffold and its reviews: [records/adopt-principles/](records/adopt-principles/2026-11-17-scaffold.md).
- First, the foundations: [record](records/first-slice/2026-11-18-first-slice.md).
- Second, editing and addressing: [record](records/second-slice/2026-11-18-second-slice.md).
- Third, figures: [record](records/third-slice/2026-11-20-third-slice.md).
- Fourth, integration (molecules, maps, structures): [record](records/fourth-slice/2026-11-21-fourth-slice.md).
- Fifth, Genomes: [record](records/fifth-slice/2026-11-23-fifth-slice.md).
- Sixth, scaffolds (`.scaffolds/`, the specs reconciled with the designs): [record](records/sixth-slice/2026-11-23-sixth-slice.md).

## Open

1. **Carried**: the 200-copy split patch's first frame and task and the merge's task over their edit bound (the author,
   2026-11-24: "Carry it"; sharing lays between identical copies is a design question, the microtask tails profiled
   first); the split and merge bound is in the TODO. The old Seq\* viewer's label churn over its allocation baseline
   (outside this package). A folded group's summary layer is named by its raw id ("ranges#summary") in its legend row.
2. **Provisional**: the root nouns ([design/building-blocks.md](design/building-blocks.md)) and the UX foundations'
   numbers ([design/ux-foundations.md](design/ux-foundations.md)), confirmed by use where a story shows them.
3. **Known gaps**: [packages/seq/TODO.md](../../packages/seq/TODO.md), whose new lines are search's (rows and cells
   not searched, a source not open skipped), a group's own menu, a values y axis, `rowsValue` over a lazy alignment,
   a layer naming a dynamic node (refused), the tracks layout's allocation pass for the zoom strip; and past version 1
   by decision: the `pairwise` and `alignment` maps, split (ranges, values), keyed tables, keeping a streamed region,
   bigWig's own levels, the glyph atlas and band copy, the rest of the deferred rows with their triggers.
4. **Defaults no story shows** (store sizes, bounds) are never confirmed by use: they stand as documented and tested.
5. **To decide as the work reaches them**: worker boundaries, versioning under 0.x
   ([design/version-1.md](design/version-1.md), Open).

## Next

The close-out's fold-back with the author and its record, then the next plan, not yet written. Prism on the new Seq\*,
then removing `packages/sequences` and renaming `seqx` to `seq`, are separate work (the author, 2026-11-23). Design
before code ([AGENTS.md](AGENTS.md), How work happens here).

## Which design to read

- Any work: [design/version-1.md](design/version-1.md) (scope, order of work), then the area's own design
  ([design/README.md](design/README.md) lists them).
- Language or tree: [design/format.md](design/format.md), [design/language-description.md](design/language-description.md),
  [design/tree-core.md](design/tree-core.md).
- Queries and addresses: [design/locations-and-queries.md](design/locations-and-queries.md),
  [design/addressing.md](design/addressing.md).
- Data by region: [design/protocol.md](design/protocol.md), [design/data-store.md](design/data-store.md),
  [design/dynamic-nodes.md](design/dynamic-nodes.md), [design/region-reads.md](design/region-reads.md).
- Drawing: [design/layout-and-draw-lists.md](design/layout-and-draw-lists.md), [design/rendering.md](design/rendering.md),
  [design/ux-foundations.md](design/ux-foundations.md).
- The viewer: [design/viewer-composition.md](design/viewer-composition.md), [design/viewer-models.md](design/viewer-models.md),
  [design/tooltips.md](design/tooltips.md).
- Genomes: [design/genome-coordinates.md](design/genome-coordinates.md), [design/genome-formats.md](design/genome-formats.md),
  [design/across-sections.md](design/across-sections.md), [design/ring-extension.md](design/ring-extension.md);
  [design/examples.md](design/examples.md) (Genomes).
- What the code says today: each package's `_spec/README.md`.

## Language version

`0.11.0`, unreleased: not yet versioned for readers outside this repository; no migration between the 0.x versions.
