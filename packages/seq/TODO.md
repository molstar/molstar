# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

The first slice runs: the tree core, the language's first kinds with the query language's first forms, resolve, draw
lists, the SVG and Canvas2D backends, the viewer's models, a source read by region (synthetic, FASTA with its index),
dynamic nodes with templates, the examples package and the measurements. In this order:

- **Selection, highlights, cursor**: the models, the click and drag-select bindings in the default preset, projected
  selection, and the selection tint's place against the data.
- **Patches and `diff`** in the tree core, then `SpecModel.applyPatch`, history, previews, undo.
- **bigWig** as a format runtime, with range reads of assets (a base-library change, below).
- **Structures**: the structure port and runtime, structure maps, the Mol\* showcase.
- **The remaining layers and layouts**, the other datasets, the UX foundation stories, the typed configuration.
- **MessagePack** (from io).

## Base library (`@molstar/core`, `@molstar/io`)

- **`repr`**: `checkRefs` (the resolve pass), `composeLanguage`, `patch` and `diff`, `descriptors`, issue dedupe, the
  catalog docs generation, a version-migration hook, `fallback` handling, a read-side `ParamValue` with typed columns
  only; a version-aware validate (an unknown kind in a newer minor version is a warning); staged `pending` unioned
  across a chain of resolved segments (only the last resolver's survives); union normalization re-validates slot
  subtrees (O(depth × size)); spelling checks of kinds and ids on trees built in code; integer columns with missing
  values (a `{ values, valid }` shape or a sentinel per type; refused until a column needs one).
- **`data/int` half-open ranges**: `combine` as a linear two-pointer merge into a presized array (it sorts the cuts and
  builds number arrays today); tests for `Empty`, non-integers and `overlaps` without a length.
- **The asset manager**: range reads (a page cache by asset and page, `Blob.slice` for files).
- **MessagePack** codec in io runs the codec suite (`checkCodec`).

## `@molstar/seqx-core`

- **Addressing**: an item's coordinate system (`<node>/(item :id k)`) as a prefix of location forms (the `item` form
  gives Elements today, so `/refs/t1/(item :id T1)/(range …)` fails); declared numberings (built-ins only, so the way
  back prints `ordinal` or `index`); rects projected for a focus; a whole Location's size bound.

- **Language**: the theme sublanguage slot for `layer.color` and a union of a ref and a query slot for `layer.data`
  (both loose today); the `snapshots` level; a strand column check; codon exons checked against the protein's length
  (stop codons); a section's `range` as a Location; refs followed in data nodes' params and a section's `of` and
  `follow` (a layer's `data` follows them); a `label` param on nodes (an instance's label is the template's plus its top
  ancestor's key); `seq.lang.ref.target` overlaps what resolve reports: keep one; the package index re-exports
  `SeqLanguageVersion` as `LanguageVersion`: one definition.
- **Queries**: dynamic nodes in `not`'s universe; open or streamed systems in `(all)` and `not` (a pending length on
  `QuerySystem` instead of `openQueryLength`); a Location limit (what it counts); a label lookup across a codon step (a
  unit's numbering does not cross one).
- **Resolve**: the isoform half of the closure rule; NCBI codes beyond 1, 2 and 11; projecting rects; reverse paths
  carried across models; an item system's length when letters show a stop (it counts as a residue now); mixed strands
  within one item (the first row's is taken); an allocation check of `projectRanges`. The TP53 exon coordinates in
  `resolve/_test/p53-geometry.ts` were written from memory: check them against RefSeq NM_000546.
- **Overlay and templates**: re-expand only the template an arrival touches, keep instance layers across overlay
  versions, key the overlay on what it reads of the model (a view-only edit rebuilds it); a static layer whose `data`
  names a dynamic node waits for it instead of reporting `seq.view.layer.data`; the selection and highlights hold their
  nodes (shown placements and `hold` already do); the overlay's sync issues in the viewer's `layoutIssues`.
- **Data**: `summaries` in `ViewerConfig`; summaries for systems of unknown length or past `maxLength`; byte dedupe
  across tiles; `SourceEntry` carries only `address`, `format`, `uri`, `value`.
- **Render, the axis** (the third slice, stage 1): letters through a one-to-one or codon path (a protein's letters on
  its gene, a residue a codon wide); a section's direction (`genomic`, `reversed`, `follow`) on reversed pieces, the
  reversed axis counting down, reversed DNA letters as their complement; collapsed introns.
- **Render**: a per-section viewport map in the frame input (one window serves every leading section now); an
  `Emphasis.projected` so the hover model's projection onto linked sections is drawn in the `hover` pass; highlight
  channels (the pass exists, empty); a summary pyramid for inline values zoomed out (inline data is level 1 only, so a
  large one zoomed out hits the tile cap); a legend's place in the chrome; the Nightingale reference adapter; a `types`
  option on `renderSvg`, then layer-type, theme and palette registry entries, atomic, counted and undoable at run time;
  minus-strand followers show the leader's hull in genomic order (reversed pieces and a config item later); `query`
  layers read the whole annotation (filtering comes with the query sublanguage); `values` line and area, points and
  split; the ranges layout's allocation an item in view (the `earlier` and `occurrences` sets, `StickyLanes.keep`'s map,
  keys of tiles without ids, `placeholderSpans` and `dataRegions` arrays: about 300 B an item); the values direct path
  for linear maps (offset paths) too; the renderer's rebuild sorts the reads and builds the key's text each layout
  (bounded by tiles in view).

## `@molstar/seqx-viewer`

- Addressing: the highlight channels drawn (the `highlight` pass is empty); the selection's gestures, anchor and cursor
  over `SelectionModel`; copy link from a context menu and a key; a link's hover projected onto linked sections.
- Tests: a frame tick asked for on fake timers a test then drops leaves `@molstar/ui/dom`'s frame lane waiting for it
  (`requested` stays true for the same host function), so a later file's viewer draws no frame; the addressing twins
  swap `DomHost.requestAnimationFrame` per test. Fix in the frame lane (forget a request across a timer reset) or in
  `vitest.setup.ts`.
- Wrapping pans on a cyclic axis (only clamped now); scroll groups across viewers and views (EC-113); ring windows
  (EC-114).
- Touch pinch (two pointers); the Mol\* structure port; the `@molstar/ui` keymap and command manager, replacing the
  shortcut match scoped to the container; `ParamDefinition`s for the config items.
- An allocation tier for the viewer (`test:perf`: the Canvas2D draw, `store.overlay`, `_viewport`, `hover.set`, the
  bindings' resolve); an allocation check by hand of the canvas draw and hover loops.
- Browser measurements: canvas frame times (draw included) and first paint; the Gene by region story rechecked in a
  visible browser (frames do not run in a hidden pane).

## `@molstar/seqx-examples`, `@molstar/seqx-storybook`

- The other datasets and their views (Antibody, Plasmid, tRNA, Globins).
- `vitest` in the Storybook package's dev dependencies; the Vite plugin's test moves under `src`.
- The Features and Examples pages.

## Deferred

- **The Genomes milestone**: layers across sections (ribbons, a dot plot with its own y section), synteny, contact maps,
  MAF over a genome, values for genomics (points, split, log scale, reference lines), variants by sample, matrices,
  keyed joins, gene hierarchies, a name index, two assemblies on one axis, and the formats they need (2bit, BigBed, VCF,
  PAF, GFF, chain files). The coordinate model must not rule them out.
- **Smaller layer features** (numeric filters, lane order, px heights): each with the example that needs it.
- **Definitions and uses** (`def`, arguments, templates per member of a set): when an example needs them.
- **Layouts and chrome parts as public extension points**: when an outside application needs one.
- **Settings in layers and a generated settings form**: when the typed configuration has users.
- **The compute port (workers)**: when a profile shows the main thread is the limit.
- **A JSON Schema** derived from the language description: when an outside reader or an agent tool needs to validate
  without the core.
