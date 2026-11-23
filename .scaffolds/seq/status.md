# Status

Start here: what is built, what is open, what comes next, and where to read. The current state only; how each slice
went is in its record.

## Built

- **`@molstar/seqx-core`** (`packages/seq/core`): SeqViewSpec, the language at `0.10.0`, with its builder, the checks
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
- **Measured** at the fifth slice's close: every frame target met but the 200-copy split patch (median 16.2 ms) and one
  Signal zoom frame (25.7 ms); sizes 260.2 / 974.3 / 1,539.7 / 1,086.9 kB min (reading, headless, default, minimal)
  against 265 / 980 / 1,545 / 1,090 ([performance](design/performance.md), Size).

Which features are built: [design/features.md](design/features.md); which edge cases have an outcome:
[design/edge-cases.md](design/edge-cases.md), Outcomes.

Slices, each closed with its record:

- The scaffold and its reviews: [records/adopt-principles/](records/adopt-principles/2026-11-17-scaffold.md).
- First, the foundations: [record](records/first-slice/2026-11-18-first-slice.md).
- Second, editing and addressing: [record](records/second-slice/2026-11-18-second-slice.md).
- Third, figures: [record](records/third-slice/2026-11-20-third-slice.md).
- Fourth, integration (molecules, maps, structures): [record](records/fourth-slice/2026-11-21-fourth-slice.md).
- Fifth, Genomes: [record](records/fifth-slice/2026-11-23-fifth-slice.md).
- Sixth, scaffolds (`.scaffolds/`, the specs reconciled with the designs): in progress, `execution/sixth-slice/`.

## Open

1. **For the author, carried to the next plan**: the use gate and the novice passes' 14 design questions; the minimal
   viewer's +150.7 kB (accept or trim); five feature rows out of `deferred`.
2. **Provisional**: the root nouns ([design/building-blocks.md](design/building-blocks.md)) and the UX foundations'
   numbers ([design/ux-foundations.md](design/ux-foundations.md)), confirmed by use where a story shows them.
3. **Known gaps**: the close-out plan (`execution/close-out/plan.md`: the drag rule, the stage 1 novice decisions,
   the `pairwise` and `alignment` maps, the query forms `search`, `onto`, `group`, track groups, the Manhattan family,
   `split`) and [packages/seq/TODO.md](../../packages/seq/TODO.md) (M5's leftovers and the rest).
4. **Defaults no story shows** (store sizes, bounds) are never confirmed by use: they stand as documented and tested.
5. **To decide as the work reaches them**: worker boundaries, versioning under 0.x
   ([design/version-1.md](design/version-1.md), Open).

## Next

The sixth slice's finalize. Then version 1's close-out (the planned rows no milestone took, the carried items, the use
gate). Prism on the new Seq\*, then removing `packages/sequences` and renaming `seqx` to `seq`, are separate work (the
author, 2026-11-23). Design before code ([AGENTS.md](AGENTS.md), How work happens here).

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

`0.10.0`, unreleased: not yet versioned for readers outside this repository; no migration between the 0.x versions.
