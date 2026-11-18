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

Which features are built, partly or not at all: [design/features.md](design/features.md). Which edge cases have an
outcome: [design/edge-cases.md](design/edge-cases.md), Outcomes.

## Open

1. **For the author**: the choices the slice made and that await confirmation
   (the [first slice's record](records/first-slice/2026-11-18-first-slice.md), Open decisions): undeclared columns tagged
   in JSON; wrapping ranges on cyclic systems; the query decisions (`chain` as Nodes, quoted `.n` labels, boundary at a
   system's end, structural Location equality); the protocol's and store's details; the language's codes and alphabets;
   the configurable defaults of every lane; the theme source; the retention of shown nodes; the per-source summaries.
   Also the selection tint's order (over the data, against "under") and the density units (letters in CSS px, the rest
   in device px).
2. **Needs the author's call, a base-library change**: the asset manager's fetch (credentials, redirects, signal, byte
   limit) and range reads ([design/base-library.md](design/base-library.md)).
3. **Provisional, still to confirm**: the root nouns ([design/building-blocks.md](design/building-blocks.md)), the
   source protocol's shape beyond the slice ([design/sources.md](design/sources.md)), `values.coarser`'s default (2, or 1 now that the direct path makes level 1 cheap), the UX foundations' numbers
   ([design/ux-foundations.md](design/ux-foundations.md)).
4. **Known gaps**: [packages/seq/TODO.md](../packages/seq/TODO.md): projected hover not drawn, one window for all
   leading sections, a source-backed leader's length not fed to the view model, every template re-expands on each
   arrival, the TP53 test geometry unchecked against RefSeq, and what is deferred to the Genomes milestone.
5. **Tooling, optional at gates**: a script that reverts a fix's non-test files on a scratch worktree and expects
   its new test to fail (the fold-back's "a fix's test is shown failing without the fix",
   [principles/README.md](principles/README.md)); until then the report names the command and the reviewer reruns it.
6. **To decide as the work reaches them**: worker boundaries, structure and map runtimes, accessibility, versioning
   under 0.x ([design/version-1.md](design/version-1.md), Open).

## Next

The second slice, planned and waiting for its Gate 0: version 1's first milestone, core editing and addressing, with
selection, the security bounds and MessagePack; history is a Seq\* port, the host owning the undo stack. Then figures
(the third slice, drafted), integration, and the Genomes milestone ([design/version-1.md](design/version-1.md),
Milestones). Design before code for each, as in [AGENTS.md](AGENTS.md), How work happens here.

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

`0.0.0`: the first slice runs; the language is not yet versioned for readers outside this repository.
