# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

The first slice built the tree core, the language's first kinds, resolve, draw lists, the SVG and Canvas2D backends, the
viewer's models and data by region; the second, editing and addressing (patches, history, previews, story steps, refs,
links, selection, the security bounds, MessagePack). Next, version 1's second milestone, figures, in this order:

- **Themes and the small layers**: the theme sublanguage, residue schemes, legends, shared domains; the ruler, secondary
  structure, variants, values' line and area.
- **Layouts**: text and alignment figures, blocks and rings, navigators; then pairs, grids and dot plots, the logo, the
  rows layer.
- **The other datasets** and their views, each with an editing view; PNG and the render CLI.
- Later: **bigWig** with range reads of assets (a base-library change, below); **structures** (the structure port and
  runtime, structure maps, the Mol\* showcase).

## Base library (`@molstar/core`, `@molstar/io`)

- **`repr`**: `checkRefs` (the resolve pass), `composeLanguage`, `descriptors`, issue dedupe, the catalog docs
  generation, a version-migration hook, `fallback` handling, a read-side `ParamValue` with typed columns only; a
  version-aware validate (an unknown kind in a newer minor version is a warning); staged `pending` unioned across a
  chain of resolved segments (only the last resolver's survives); union normalization re-validates slot subtrees
  (O(depth × size)); spelling checks of kinds and ids on trees built in code; integer columns with missing values (a
  `{ values, valid }` shape or a sentinel per type; refused until a column needs one).
- **`data/int` half-open ranges**: `combine` as a linear two-pointer merge into a presized array (it sorts the cuts and
  builds number arrays today); tests for `Empty`, non-integers and `overlaps` without a length.
- **`repr/patch`**:
  - derive the next revision's ref index at the end of each apply: a structural edit on a 20,000-ref spec costs about 10
    ms, mostly rebuilding the old root's index (0.7 ms cached), against the 16 ms edit target;
  - many moves under one parent are quadratic in `applyPatch` (reversing 16k children: 105 ms, 229 ms in a merge;
    inserts and diffs are linear): done when such a batch is linear, held by a scaling test;
  - `RefSite.path` holds child indices, so an insert in a wide parent re-creates its later siblings' sites, and nothing
    reads it: drop it or make it lazy;
  - a key typed as a ref with different types in different kinds is walked with each: a conformance case once a language
    has one.
- **`repr/tree`**: `tree/address.ts` (the child index) is reachable through core's `./*`; io imports `repr/tree/index`
  by path, so closing `./repr/tree/*` needs io moved to `./repr/tree` first. `ChildIndexThreshold` (8) is not measured.
- **Codecs**: `checkCodec`'s trees gain a leading-U+FEFF string, key and category entry; the MessagePack codec's issue
  codes (`syntax`, `duplicate_key`, `key_type`) move into `CodecIssues` beside JSON's and the S-expression's; the
  `Codec` documentation says `encode` throws on a tree no reader would give back (a lone surrogate).
- **MessagePack** (`common/msgpack`): a duplicate map key keeps the last silently (JSON reports it): an option; payload
  alignment is by position: encoder padding when zero-copy decoding matters more; a canonical-category fast path;
  key-order determinism shared with JSON. Each when a measured use asks.
- **The asset manager**: range reads (a page cache by asset and page, `Blob.slice` for files), with bigWig. `maxBytes`
  is ignored on `file://` reads in Node (`fetchData`'s file path ignores the four options): done when the file path
  honors them and the conformance suite covers it. The suite sees a manager through `fetch` only; a manager that
  downloads by other means is outside it.
- **`util/data-source`**: a `file:` URL's path is opened without percent-decoding, so a root or file with a space or a
  non-ASCII character is unreadable and Windows paths open as `/C:/…` (fails safe). When it moves to `fileURLToPath`,
  the headless policy's file check is reviewed again for Windows segments ending in a dot or space, device names and `:`
  streams.
- **Type checks of core's tests**: core's other `_test` folders have 73 type errors; `tsconfig.test.json` covers
  `data/immutable` only. Done when every `_test` folder is in it and clean.

## `@molstar/seqx-core`

- **Addressing**: an item's coordinate system (`<node>/(item :id k)`) as a prefix of location forms (the `item` form
  gives Elements today, so `/refs/t1/(item :id T1)/(range …)` fails); declared numberings (built-ins only, so the way
  back prints `ordinal` or `index`); rects projected for a focus; a whole Location's size bound.

- **Language**: the theme sublanguage slot for `layer.color` and a union of a ref and a query slot for `layer.data`
  (both loose today); a strand column check; codon exons checked against the protein's length (stop codons); a section's
  `range` as a Location; refs followed in data nodes' params and a section's `of` and `follow` (a layer's `data` follows
  them); a `label` param on nodes (an instance's label is the template's plus its top ancestor's key);
  `seq.lang.ref.target` overlaps what resolve reports: keep one; the package index re-exports `SeqLanguageVersion` as
  `LanguageVersion`: one definition.
- **Queries**: dynamic nodes in `not`'s universe; open or streamed systems in `(all)` and `not` (a pending length on
  `QuerySystem` instead of `openQueryLength`); a Location limit (what it counts); a label lookup across a codon step (a
  unit's numbering does not cross one).
- **Resolve**: the isoform half of the closure rule; NCBI codes beyond 1, 2 and 11; projecting rects; reverse paths
  carried across models; an item system's length when letters show a stop (it counts as a residue now); mixed strands
  within one item (the first row's is taken); an allocation check of `projectRanges`.
- **Story steps** (the Features page's Snapshots row): the stories library and the animation port, placed with the
  player; a step that no longer applies shows the ops that still apply, flagged (today: the root without the step); a
  step's `from`, view (EC-85) and timing need language fields; an op that no longer applies has no issue in
  `repr/patch`, so Seq reports it on the step; a host's move can leave a step's later refs dangling in the step's own
  coordinates (refused unless `allowDangling`): watched for in the editing UI.
- **Bounds**: regular expressions in search and archives get theirs with those features; nested-source breadth (a
  fan-out of W^depth) and joins allocating before `checkAnswer` (at most the file's size) are recorded, not bounded.
- **Overlay and templates**: re-expand only the template an arrival touches, keep instance layers across overlay
  versions, key the overlay on what it reads of the model (a view-only edit rebuilds it); a static layer whose `data`
  names a dynamic node waits for it instead of reporting `seq.view.layer.data`; highlights hold their nodes (the
  selection does, `selection.holdNodes`); the overlay's sync issues in the viewer's `layoutIssues`.
- **Data**: `summaries` in `ViewerConfig`; summaries for systems of unknown length or past `maxLength`; byte dedupe
  across tiles; `SourceEntry` carries only `address`, `format`, `uri`, `value`.
- **Render**: a summary pyramid for inline values zoomed out (inline data is level 1 only, so a large one zoomed out
  hits the tile cap); a legend's place in the chrome; the Nightingale reference adapter; a `types` option on
  `renderSvg`, then layer-type, theme and palette registry entries, atomic, counted and undoable at run time; `query`
  layers read the whole annotation (filtering comes with the query sublanguage); `values` line and area, points and
  split; the ranges layout's allocation an item in view (the `earlier` and `occurrences` sets, `StickyLanes.keep`'s map,
  keys of tiles without ids, `placeholderSpans` and `dataRegions` arrays: about 300 B an item); the values direct path
  for linear maps (offset paths) too; the renderer's rebuild sorts the reads and builds the key's text each layout
  (bounded by tiles in view); `render/export.ts` still requests data for hidden layers (done when an export of a spec
  with a hidden layer reads nothing for it).
- **Emphasis**: projected selection bands on a whole axis are rounded at intron edges, where the frame rule wants the
  ends a projection cut square (a projection-aware cut flag; on a collapsed axis an end at a break is square); the
  keyboard cursor reads as a full-height hairline below 1 px a cell (a minimum width or a caret).
- **Tests**: `render/_test/frame-time.test.ts`'s per-frame max bounds flake under a parallel heavy run (18 and 50 ms
  against 16): best of three, or a percentile instead of the max. The p53 fixtures round-tripped JSON ↔ MessagePack in a
  test here (io cannot import them).

## `@molstar/seqx-viewer`

- Addressing: copy link from a context menu (and its key).
- Selection: drag snapping to whole codons over nucleotides (EC-149); rects (dot plots) projected and drag-selected;
  context-menu entries for `selection.add` and `selection.subtract`; the `genome-browser` preset without the add,
  subtract and extend drags (a drag pans there): decided with the presets; a cursor moved off-screen on a follower
  revealing itself by moving the leader.
- Editing: a rebase's `invalid` discards every change of theirs from `partial`, an unrelated one included (by design):
  the editing UI says so.
- Tests: the frame lane's dropped-tick leak is fixed in `@molstar/ui/dom`, so the stories' twins can drop their per-test
  `DomHost.requestAnimationFrame` swap.
- Wrapping pans on a cyclic axis (only clamped now); scroll groups across viewers and views (EC-113); ring windows
  (EC-114); navigators (strip, miniature, ideogram, vertical overview; EC-110: the Features page's Navigators row, a
  design first).
- Touch pinch (two pointers); the Mol\* structure port; the `@molstar/ui` keymap and command manager, replacing the
  shortcut match scoped to the container; `ParamDefinition`s for the config items.
- An allocation tier for the viewer (`test:perf`: the Canvas2D draw, `store.overlay`, `_viewport`, `hover.set`, the
  bindings' resolve); an allocation check by hand of the canvas draw and hover loops.
- Browser measurements: canvas frame times (draw included) and first paint; the Gene by region story rechecked in a
  visible browser (frames do not run in a hidden pane).

## `@molstar/seqx-examples`, `@molstar/seqx-storybook`

- The other datasets and their views (Antibody, Plasmid, tRNA, Globins), each carrying a drawing-against-the-data test
  as the p53 views do (`drawn-vs-data.test.ts`); the Complex dataset (CDK2 with cyclin A, its interface, PAE, ligands
  with AlphaFill's transplants) with pairs across chains, dot plots and ligands, its entries chosen and checked against
  their sources then.
- `vitest` in the Storybook package's dev dependencies; the Vite plugin's test moves under `src`.
- The examples' root keeps its own copy of the SVG-follow block (`_stories/follow-svg.ts`).

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
