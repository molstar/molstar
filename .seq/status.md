# Status

Start here. What is built, what is open, what comes next, and where to read.

## Built: the first slice

- The scaffold: `.seq/` (principles manifest, sources, designs, records), the checks wired ([AGENTS.md](AGENTS.md),
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
(copy ends) is parked unmerged for the fourth slice.

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
  the default and the minimal viewer), against bounds of 195 / 780 / 1250 / 880 kB ([performance](design/performance.md), Size).

Which features are built, partly or not at all: [design/features.md](design/features.md). Which edge cases have an
outcome: [design/edge-cases.md](design/edge-cases.md), Outcomes.

## Open

1. **For the author**: the third slice's fold-back (`execution/third-slice/foldback.md`) and the questions its close
   left (the drift notes in `execution/third-slice/drift/`: configuration the designs name and the code does not have,
   EC-63's slices against glyphs, EC-138's outline against saturated fills).
2. **Range reads of assets**: wait for bigWig ([design/base-library.md](design/base-library.md)); the asset manager's
   options are built.
3. **Provisional, still to confirm**: the root nouns ([design/building-blocks.md](design/building-blocks.md)) and the UX
   foundations' numbers ([design/ux-foundations.md](design/ux-foundations.md)), confirmed by use where a story shows
   them. The source protocol's shape is agreed with the files design ([design/files.md](design/files.md)).
4. **Known gaps**: [packages/seq/TODO.md](../packages/seq/TODO.md) and the third slice's `followups.md`: copy ends with
   units (CE, parked), the query forms `search`, `onto`, `expand`, `group`, pairs, ranges and maps on an alignment's
   columns (0.5.0), drags snapped to codons, the cursor across block gaps, track groups and summaries, and what is
   deferred to the Genomes milestone.
5. **Tooling, optional at gates**: a script that reverts a fix's non-test files on a scratch worktree and expects
   its new test to fail (the fold-back's "a fix's test is shown failing without the fix",
   [principles/README.md](principles/README.md)); until then the report names the command and the reviewer reruns it.
6. **To decide as the work reaches them**: worker boundaries, structure and map runtimes, versioning under 0.x
   ([design/version-1.md](design/version-1.md), Open).

## Next

The third slice closes: stage 7's measurements (M3), the fold-back with the author's decision, the record. Then the
fourth slice (`execution/fourth-slice/plan.md`, its Gate 0 decided 2026-11-20): version 1's third milestone,
integration (molecules and units, structure syncing on a fake 3D viewer, seq data ports, extension openings, track
header parts, services, bigWig and range reads, archives). The fifth slice is the Genomes milestone, its genome ring an
extension ([design/version-1.md](design/version-1.md), Milestones). Design before code for each, as in
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
- What the code says today: each package's `_spec/README.md`.

## Language version

`0.4.0` (unreleased, pre-release: each kind added moved it, with no migration; `0.5.0` is the fourth slice's, with
molecules); not yet versioned for readers outside this repository.
