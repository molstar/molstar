# Status

Start here. What is built, what is open, what comes next, and where to read.

## Built: the first slice

- The scaffold: `.scaffolds/seq/` (principles manifest, sources, designs, records), the checks wired ([AGENTS.md](AGENTS.md),
  Checks).
- `@molstar/core/repr`: the tree core (nodes, typed columns, hashes, addresses, issues, schema, builder, JSON and
  S-expression codecs); half-open range sets in `@molstar/core/data/int`.
- `packages/seq/core`: the language's first kinds and builder (`createSeqSpec`); the query language's forms; the
  source protocol, the `seq` format, FASTA with its index, the data store (tiles, levels, cancellation, eviction,
  summaries); resolve (the static model, paths, codon and offset maps; the overlay of dynamic nodes; track templates);
  draw lists for `letters`, `ranges`, `values`; the SVG and recording backends; `renderSvg`.
- `packages/seq/viewer`: `createSeqViewer`, the models (spec, view with scroll groups and follow, hover with
  projection, store), the `hover` and `navigate` behaviors, the Canvas2D backend, the default ports.
- `packages/seq/examples` and `storybook`: the p53 views (Tracks, Protein on its gene, Gene by region) as builder
  modules, tests and stories, canvas and SVG side by side.
- Measured: frame, hover, export, resolve, allocation, bundle sizes (the [first slice's record](records/first-slice/2026-11-18-first-slice.md), Measurements).
  After the author's stage 7 decision and lane J: every pan and zoom frame under 2.3 ms headless (3M bars included),
  draw, hover, selection and cursor redraws 0 bytes warm, a layout's allocation bounded by what is in view.

## Built: the second slice

Version 1's first milestone, core editing and addressing (work folder `execution/second-slice/`, its record when it
closes). Started at `fc2bba6ce`; lanes S0, S0b, K, K2, K3, R, R+K, L, M, N, O, P, Q, S, T, U merged.

- `@molstar/core`: `repr/patch` (five ops by address, all or nothing, inverses, `diff`, `merge3`, patches as trees,
  stored patches whose refs are typed by the language, the reverse ref index derived across revisions);
  `data/immutable` (CHAMP maps and sets, ordered ones, `Vector`); the asset manager's `ResolveOptions` (credentials,
  redirect, signal, a streamed byte limit) and `checkAssetManager`. `@molstar/io`: the repr tree's MessagePack codec.
- `packages/seq/core`: the language at `0.3.0` (`refs` and `ref`, `snapshots` and `snapshot`, `layer.hidden`,
  `track.hidden`); the checks of a change (incremental, their own entry `lang/changes`); story steps; addresses, refs
  followed, links in Mol\*'s grammar and the way back; the spec's and answers' bounds; the headless asset policy;
  emphasis drawn (selection as one band a section, highlights, projected hover and selection, the cursor); every drawn
  range mark hit, `Hit.hull`.
- `packages/seq/viewer`: `SpecModel` as commits (patches, previews, story steps, a stale patch rebased on opt-in), the
  history port and `@molstar/ui/history` adapters; `SelectionModel` (several ranges, the anchor, the keyboard cursor),
  `HighlightModel`, `AddressModel` (links, copy link, the session following edits); followers that lead; the view's
  history of jumps (back, forward, the zoom keys, a double-click focusing an item's extent); the editing commands.
- `packages/seq/storybook`: the p53 Edit and undo, Bookmarks and links and Selection stories and the feature stories
  (Editing, Selection, Links, Bookmarks), each with a headless twin; the Features and Examples docs pages.
- Measured (stage 5): every time and allocation target passes (Node and Chrome at dpr 2); the read bundle 150.6 kB
  against the 160 kB target (the author), headless rendering 267.9 kB.

## Built: the third slice

Version 1's second milestone, figures (work folder `execution/third-slice/`, its record when it closes). Gate 0
decided 2026-11-18; over sixty lanes merged (C0, AX, SL, TH, TX, TA, TV, TE, RB, G2, SY, RW, CH, PG, RF, NV, GD, SO,
FS-A1, FS-A2, IZ, CR, HL, PR, QS, LB, RC, DP, DT, SW, AY, OP, DG, DC, NB, CY, TW, DA, QT to QT3, NP, HC, SF, RK, ET, AT,
CS, SC, MF, CF, MV, UN, PG (PCG32), LP, RX, DB, DS, PZ, PF, LX; the list in `execution/third-slice/state.md`). CE
(copy ends) was parked unmerged, then merged by the fourth slice's CU.

- The axis: a section's direction, another system's letters through a map, collapsed introns of zero width on every
  layout, reversed pieces; the ruler's numberings (insertions, `unnumbered`, `units`).
- Themes: the theme sublanguage (field palettes, rules as words, schemes with safe substitutes, `scale` nodes, size
  channels), legends as data and a part a host places, the reader's color-blind substitution; the lists-free palette
  module and PCG32 in `@molstar/core`.
- Layouts: tracks fitting their content (`heights`, the `trackSizing` port), text and alignment figures, blocks and
  chips, rings with a window and label groups, through one layout interface; navigators (the ruler rail, strips, an
  ideogram, a vertical overview, a miniature, the zoom strip).
- Layers: the small layers (ruler, secondary, variants, values as line and area, ranges' styles), pairs (arcs, chords,
  linked pairs), grids and dot plots with their navigation, the logo, the rows layer (columns, clades, pinned and lazy
  rows).
- Data: aligned FASTA, Stockholm, Newick, records (TSV, CSV, JSON, MessagePack, GFF3), a source's own numbering,
  lengths through resolve; the six datasets and p53's chromosome by region, generated where large.
- The viewer: stories (the player, a step's view), the keyboard cursor over rows and grids, the hover card (a port) and
  the keys panel, highlights that overlap, accessibility (an SVG's figure and items, announcements), a minimal viewer
  composed by the spec's registry.
- Export: PNG, and the `seqx-render` command line.
- Sizes (after RR, the render reduction): 175.4 / 711.8 / 1142.5 / 801.7 kB min (reading a spec, headless rendering,
  the default and the minimal viewer), against bounds of 230 / 780 / 1250 / 880 kB (reading a spec moved to 230 for language 0.5.0, the author, 2026-11-20) ([performance](design/performance.md), Size).

## Built: the fourth slice

Version 1's third milestone, integration (work folder `execution/fourth-slice/`, its record when it closes). Gate 0
decided 2026-11-20; lanes AR, NG, ML1, FX1, SZ, DP, BW, ML2, MC, FX2, MK, QE, ST1, AC, HP1, MR, CU, MX, HP2, ST2, FX3,
ST3, HD, HP3, FX4, EX, LE, CN, NP, SV, FX5, PC, FX6, HT merged (the findings by lane in its `state.md`).

- Molecules: language 0.5.0 (molecules, compounds, complexes, units, links, modified positions; `source.ref`,
  `values.field`), the model (unit, molecule and atoms systems, copies, links per unit, `maxUnits`), a section of a
  complex drawn as blocks with block forms, molecule letters, modified marks, links as arcs, a pair's ends on split
  copies by unit; the examples on the new kinds (4HHB, 1EHZ, p53 over 2OCJ, 3AVE, 1FIN, the tRNA seed, a T-DM1 ADC);
  `@molstar/core`'s standard components.
- Maps: numberings by `offset` and `segments`, `segments` and `structure` maps, two maps of one pair by `priority` and
  `via`; pairs, ranges and maps on an alignment's columns; a protein over its coding DNA with the codon connector.
- Queries by role: `elements` and `groups`, `filter` over typed properties, `expand`, `linked`, `merge`, `intersect`,
  `except`, Keys.
- Structures: `@molstar/seqx-core/structure` (`structureFacts`, `toUnits`/`fromUnits` with the size gate, colors by
  run), the viewer's wiring and `linkStructure`, a fake structure viewer and its stories; an opt-in mmCIF runtime.
- Data: seq data ports (`changes`, `withAdapter`, `io/structure`, the viewer reopening in place; the fake host's
  adapter), range reads of assets in `@molstar/core` and an opt-in bigWig runtime, archives over the assets port,
  services (UniProt, PDBe SIFTS, AlphaFold DB) from recorded fixtures in the Storybook.
- The viewer: track headers (core's header line, the `header-action` registry kind, the `headers` behavior with its
  menu and parts ports, flyout submenus; hidden tracks listed and restored by the host), extensions (`SeqExtension`,
  `<id>:<name>` types with fallbacks, `extend()` at run time), custom navigator pictures and a host-drawn navigator,
  legends' library surface.
- Rendering: a test-only second consumer of the resolved view; `render/tracks.ts`, `rows.ts`, `ranges.ts`,
  `viewer.ts` and `models/view.ts` split under §module-size (the size baseline empty).
- Sizes at stage 4's gate: 223.1 / 806.0 / 1296.4 / 929.5 kB min (reading a spec, headless rendering, the default and
  the minimal viewer) against bounds of 230 / 830 / 1330 / 940 kB, each bound moved by the author with the feature
  behind it ([performance](design/performance.md), Size). Not measured since FX5, SV, PC, FX6 and HT merged: the
  finalize's gate.

## Built: the fifth slice

Version 1's fourth milestone, Genomes (work folder `execution/fifth-slice/`, its record when it closes). Gate 0 decided
2026-11-20, the six designs promoted 2026-11-21; 56 lanes merged: BZ, GC1, SC, GC2, BE, GC3, GG, NS, BS, VC, GB, RG
(with RS), GA2, FR, GA3, RL, LA, LA2, LA3, LA4, VL, SX, VP, PF, HC, TT, GX, RX, RX2, RX3, DX and FX1 to FX24 (the
findings by lane in its `state.md`).

- Coordinates: assemblies, `section.of` a list, `section.include`, `block.direction` and `.scale`, `navigator.switch`,
  the locus and block forms, the names rule; the genome axis on rows and rings (runs of small blocks, `×N`); navigators
  over blocks and the square; sections over a source's complex reflowing on its changes.
- Data by region: the `pairs`, `matrix`, `grid` and `rows` shapes, `RegionRequest.y`, 2D tiles, `readCost` and
  `maxBytesPerLayer`, the pairs level and `binPairs`, grid summaries; data over an assembly by sequence, inline by name.
- Formats: BGZF and tabix in `@molstar/io`; `runtimes/bed` (BED, bedGraph, BEDPE, chrom.sizes), `gff` (GFF3, GTF as
  transcripts, names indexed), `vcf` (sites, genotype classes), GenBank; tabix levels with exact counts.
- Across sections: the bridge (ribbons drawn once, lock, align, select across), grids over two sections, categorical
  axes, the liftover.
- The viewer: the hover contract and its sweep, legends both ways, tooltips composed by providers, values as points and
  on a log scale, sections scrolling past a cap, track sizing, no track an empty well; a card over blocks lists the
  pointer's block's items once; the square's moves are history steps.
- The genome ring, `@molstar/seqx-ext-rondo` (`section.center`): stems, ribbons by count, stubs, hues, cytobands, the
  square a window, element and stem menus, kinds of links; `PairsReader` and `LayerType` hooks in core.
- The Genomes dataset generated byte-identically and its views; the builder showcase; a Nightingale showcase; the
  storybook's Docs (the engine, not the game).
- Measured (M5, [the record](records/fifth-slice/2026-11-23-fifth-slice.md)): every frame target met but the 200-copy split
  patch (median 16.2 ms) and one Signal zoom frame (25.7 ms); vsync without drops. Sizes 260.2 / 974.3 / 1,539.7 /
  1,086.9 kB (reading, headless, default, minimal) against 265 / 980 / 1,545 / 1,090; the ring 43.9 of 45; each runtime
  9.3 to 23.1 of 30; the bridge chunk 25.4 of 28.

Which features are built, partly or not at all: [design/features.md](design/features.md). Which edge cases have an
outcome: [design/edge-cases.md](design/edge-cases.md), Outcomes.

## Open

1. **For the author**: the fifth slice's fold-back (`execution/fifth-slice/foldback.md`) and record; the minimal
   viewer's growth of 150.7 kB ([the record](records/fifth-slice/2026-11-23-fifth-slice.md)); the novice passes' 14 design
   questions (`execution/fifth-slice/novice-findings.md`) and the use gate, carried.
2. **Provisional, still to confirm**: the root nouns ([design/building-blocks.md](design/building-blocks.md)) and the UX
   foundations' numbers ([design/ux-foundations.md](design/ux-foundations.md)), confirmed by use where a story shows
   them.
3. **Known gaps**: [packages/seq/TODO.md](../../packages/seq/TODO.md): the drag rule (a drag's work bounded by the
   screen), the stage 1 novice decisions F3, F9, F11 and F14 (carried), the `pairwise` and `alignment` map kinds, the
   query forms `search`, `onto`, `group`, drags snapped to codons, track groups and summaries, the Manhattan family,
   `split`, EC-224's join (the close-out).
4. **Tooling, optional at gates**: a script that reverts a fix's non-test files on a scratch worktree and expects
   its new test to fail (the fold-back's "a fix's test is shown failing without the fix",
   [principles/README.md](principles/README.md)); until then the report names the command and the reviewer reruns it.
5. **To decide as the work reaches them**: worker boundaries, versioning under 0.x
   ([design/version-1.md](design/version-1.md), Open).

## Next

The fifth slice closes: the fold-back with the author's decision, the gate once more, the record. Then the sixth slice
(`execution/sixth-slice/plan.md`, draft): `.scaffolds/seq/` reconciled into a common place for scaffolds. Then version 1's
close-out (the planned rows no milestone took, the fifth slice's carried items, the use gate). Prism on the new Seq\*
and the old package's removal are separate work (the author, 2026-11-23). Design before code, as in
[AGENTS.md](AGENTS.md), How work happens here.

## Which design to read

- Any work: [design/version-1.md](design/version-1.md) (scope, order of work), then the area's own design
  ([design/README.md](design/README.md) lists them).
- Language or tree: [design/format.md](design/format.md), [design/language-description.md](design/language-description.md),
  [design/tree-core.md](design/tree-core.md).
- Queries and addresses: [design/locations-and-queries.md](design/locations-and-queries.md),
  [design/addressing.md](design/addressing.md).
- Data by region: [design/protocol.md](design/protocol.md), [design/data-store.md](design/data-store.md),
  [design/dynamic-nodes.md](design/dynamic-nodes.md).
- Drawing: [design/layout-and-draw-lists.md](design/layout-and-draw-lists.md), [design/rendering.md](design/rendering.md),
  [design/ux-foundations.md](design/ux-foundations.md).
- The viewer: [design/viewer-composition.md](design/viewer-composition.md), [design/viewer-models.md](design/viewer-models.md).
- Genomes: [design/genome-coordinates.md](design/genome-coordinates.md), [design/region-reads.md](design/region-reads.md),
  [design/genome-formats.md](design/genome-formats.md), [design/across-sections.md](design/across-sections.md),
  [design/ring-extension.md](design/ring-extension.md), [design/tooltips.md](design/tooltips.md); [design/examples.md](design/examples.md) (Genomes).
- What the code says today: each package's `_spec/README.md`.

## Language version

`0.10.0` (unreleased, pre-release: the fifth slice's bumps, 0.6.0 assemblies and the declared params, 0.7.0 inline
data by name, 0.8.0 points, 0.9.0 tooltips, 0.10.0 a ring's `center`, with no migration); not yet versioned for readers
outside this repository.
