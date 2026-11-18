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

Which features are built, partly or not at all: [design/features.md](design/features.md). Which edge cases have an
outcome: [design/edge-cases.md](design/edge-cases.md), Outcomes.

## Open

1. **For the author**: the second slice's fold-back (`execution/second-slice/foldback.md`) and its "built but not
   agreed" queue (`execution/second-slice/gaps.md`: `SizeLimitExceeded`, the redirect option, MessagePack's extra
   extension codes, io's one-line fix, `ChildIndexThreshold` and the patch tree's layout, `FlatThreshold` 16 and
   `Vector.indexOf`, `track.hidden`, the selection's look). Still to confirm from before: the selection tint's order
   (over the data, against "under") and the density units (letters in CSS px, the rest in device px). The first slice's
   open decisions 1–9 were confirmed at the second slice's Gate 0.
2. **Range reads of assets**: wait for bigWig ([design/base-library.md](design/base-library.md)); the asset manager's
   options are built.
3. **Provisional, still to confirm**: the root nouns ([design/building-blocks.md](design/building-blocks.md)), the
   source protocol's shape beyond the slice ([design/sources.md](design/sources.md)), the UX foundations' numbers
   ([design/ux-foundations.md](design/ux-foundations.md)). `values.coarser` is 1 (decided at Gate 0).
4. **Known gaps**: [packages/seq/TODO.md](../packages/seq/TODO.md): one window for all leading sections, a ref's
   placement in data params and a section's `of`/`follow`, codon snapping of drags, the view history's chrome, the
   stories library and the animation port, navigators, every template re-expanding on each arrival, and what is
   deferred to the Genomes milestone.
5. **Tooling, optional at gates**: a script that reverts a fix's non-test files on a scratch worktree and expects
   its new test to fail (the fold-back's "a fix's test is shown failing without the fix",
   [principles/README.md](principles/README.md)); until then the report names the command and the reviewer reruns it.
6. **To decide as the work reaches them**: worker boundaries, structure and map runtimes, accessibility, versioning
   under 0.x ([design/version-1.md](design/version-1.md), Open).

## Next

The second slice closes: the fold-back with the author's decision, the record, then the third slice's plan rewritten
from its draft (`execution/third-slice/plan.md`): version 1's second milestone, figures, starting with the axis (a
section's direction, letters through a map, collapsed introns: the stage 4 gate's findings). Then integration and the
Genomes milestone ([design/version-1.md](design/version-1.md), Milestones). Design before code for each, as in
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

`0.3.0` (unreleased, pre-release: each kind added moved it, with no migration); not yet versioned for readers outside
this repository.
