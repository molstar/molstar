# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

The first slice built the tree core, the language's first kinds, resolve, draw lists, the SVG and Canvas2D backends, the
viewer's models and data by region; the second, editing and addressing (patches, history, previews, story steps, refs,
links, selection, the security bounds, MessagePack). Next, version 1's second milestone, figures, in this order:

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

- **The small layers** (built: ruler, secondary, variants, values' line and area, ranges' line styles and arrowheads):
  - channels they do not read yet: `size` (a variant's radius, a line's thickness, a band's share), `opacity` (the
    area's 0.35, a highlight's, a variant's), `shape` (variants' glyph: circle only), `line` and `arrow` by a form (the
    slots take forms; the layers read a word); a values color by a field other than its value;
  - `secondary`'s `style: cartoon` (with the text figure, layouts); a declared numbering for the ruler (the language
    declares none: only `index` and `ordinal` can be asked for, the insertion labels are tested through a numbering
    handed to the layer);
  - a tiled source's items have no arrowhead (their other rows may lie outside what was read): a source that answers an
    item's extent would give them one; the three layers' `nextMark` and hover passes are still their own;
  - `arrow` points by an item's own start and end, not by its strand: a strand column (EC-232 asks arrows only when
    asked) would say which way an item runs;
  - a ruler's labels over a numbering that is not monotone (IMGT) in a piece of more than 4,096 positions are looked up
    by number and skipped when the numbering jumps over 16,000 labels a piece;
  - a line over a reversed collapsed axis in the device-column step does not break at the axis's breaks (the columns are
    contiguous).

- **Addressing**: an item's coordinate system (`<node>/(item :id k)`) as a prefix of location forms (the `item` form
  gives Elements today, so `/refs/t1/(item :id T1)/(range …)` fails); declared numberings (built-ins only, so the way
  back prints `ordinal` or `index`); rects projected for a focus; a whole Location's size bound.

- **Language**: a union of a ref and a query slot for `layer.data` (loose today); a strand column check; codon exons
  checked against the protein's length (stop codons); a section's `range` as a Location; refs followed in data nodes'
  params and a section's `of` and `follow` (a layer's `data` follows them); a `label` param on nodes (an instance's
  label is the template's plus its top ancestor's key); `seq.lang.ref.target` overlaps what resolve reports: keep one;
  the package index re-exports `SeqLanguageVersion` as `LanguageVersion`: one definition.
- **Queries**: dynamic nodes in `not`'s universe; open or streamed systems in `(all)` and `not` (a pending length on
  `QuerySystem` instead of `openQueryLength`); a Location limit (what it counts); a label lookup across a codon step (a
  unit's numbering does not cross one).
- **Resolve**: the isoform half of the closure rule; NCBI codes beyond 1, 2 and 11; projecting rects; reverse paths
  carried across models; an item system's length when letters show a stop (it counts as a residue now); mixed strands
  within one item (the first row's is taken); an allocation check of `projectRanges`.
- **Stories** (built: `stateOf`, `plan`, `sample`, the encoder, the animation port and its two adapters, the frame slot,
  the p53 story): `order` switches at the cut (each child moving to its new place needs a y offset a track the layout
  does not have); a change of the view shown is a cut at the start, not a crossfade; a step whose state changes an axis
  (a collapse, a reverse) is placed on the step's axis while the frame shows the start's, so a section's viewport goes
  through the update rule (EC-111) on the way; **the bands' seam for the rows layer**: a band's `bandTop` and `rowScale`
  animate in the core (the view channel's `bands.<track>.*` fields, a frame's viewport `bands`), but the viewer's
  `ViewModel` holds no bands (its `Viewport` and `setViewport` take none, `_modelViews` copies none), so a step's or a
  jump's bands reach the frame only from the step's view: the rows layer adds `bands` to the viewport record and
  `setViewport`, and `bandTop` as its top row's key (display rows today); the story's own UI (steps listed, edit this
  step, store the view as a step) and the render CLI's video and frame-sequence outputs; an op that no longer applies
  has no issue in `repr/patch`, so Seq reports it on the step; a host's move can leave a step's later refs dangling in
  the step's own coordinates (refused unless `allowDangling`): watched for in the editing UI.
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
  with a hidden layer reads nothing for it); the ranges layout walks every row of the entries it reads (a 150,000-row
  inline annotation: about 114 ms a frame); a letter through a codon path allocates its parts (`projectRanges` sets its
  buffer's length to 0, which lets the storage go: 161 KB a pan frame on the coding gene).
- **Rings and blocks** (built: the rings layout with labels, the title and the window; blocks in tracks and on rings):
  - blocks in the text layout (a block starting a line, its label inline on its first line; one-position blocks as chips
    flowing after the chains, EC-99): waits for the text layout on this branch;
  - the ring's track labels (no label column on a ring), its ruler's ticks across the ruler ring only;
  - a copy group is the same system and range here: the complex kind's molecules and units name it once they exist
    (EC-65), and a merged block's hover then lights every member's system;
  - a ring's frame where only its window moved draws every pass again (the design asks for the emphasis passes alone):
    the window moves with a viewport, and a viewport the ring's window shows is one a follower (or its leader) lays out
    and draws by, so the frame's data passes are drawn again for that section anyway; the ring's own data kept needs
    passes split by section (its lays are kept already);
  - a ring's `labels` sink: only `ranges` gives its labels (secondary, variants keep theirs inside);
  - `export.window` read by `renderSvg` alone: the PNG export with it;
  - a section of blocks on several systems (its axis on no system): a drag-select stays in the block it began on (across
    blocks: a location of several systems), the keyboard cursor does not step there (a step across blocks), a focus
    moves to the first block holding the positions (to their hull across blocks), a bookmark of the view names none (an
    address of several systems' ranges), and a column-shaped selection on it draws nothing (per block).
- **Emphasis**: projected selection bands on a whole axis are rounded at intron edges, where the frame rule wants the
  ends a projection cut square (a projection-aware cut flag; on a collapsed axis an end at a break is square); the
  keyboard cursor reads as a full-height hairline below 1 px a cell (a minimum width or a caret).
- **Tests**: the p53 fixtures round-tripped JSON ↔ MessagePack in a test here (io cannot import them).

- **Themes**: a filter by a value's range (steps, a ramp: `(item :<field>)` compares atoms only) and filters through a
  tiled source's server; the legend's menu entry for `layer.filter` (the command and the chip are built); the pLDDT
  view's data is a seeded stand-in (a model's confidence from a real fixture); a hover halo on the legend's entry in the
  canvas (EC-138: the legend highlights, the marks do not ring); a values layer's y on a `log` scale (its color reads
  one, its bars are linear); a `field` form's `:data` (a field of another annotation); forced colors
  (`forced-colors: active`: the system's palette); the column measures (`share`, `identity`, `matches-reference`)
  registered by an alignment's layers; `classes` (the default protein scheme) is not safe under the three deficiencies
  (its substitute is drawn under `theme.cvdSafe`): a safe default, a design decision; the dark band (L 0.30–0.62) reads
  dim against the dark page (fills 1.4–1.7:1 against it at the low end) and cannot rise without giving up 4.5:1 ink: an
  outline or a lighter page, a design decision; an inline table's column named `kind` is refused by the s-expression
  codec (a reserved key), so `secondary` reads `category` in the examples; the `theme/` entry is reached through the
  `./*` export (`theme/index`): its own `./theme` export.

- Layouts (text): chips (one-position blocks: with blocks; the line loop is their seam), the `figure` preset's spacing,
  inline group labels before a line's first group; the line heights are kept full until laid out (the text's `end` an
  estimate); the number gutter is measured from the total's digits, not the widest label of a numbering; a frame's title
  row stays on every line its region is on (a band collapsed is not laid out, so its outline would go: the design's "a
  title row only where a region starts" waits for an outline drawn in the bands it covers); the zoom's new breaking
  under `pieceBreak: 'line'` reads the pieces kept in the geometry, and collapsed lines' heights make its "same line" an
  estimate.
- Layouts (interface): `FrameInput.viewport` is still the default for sections without an entry (the core's tests and
  exports use it); the viewer passes every leading section in `viewports`. Its removal is a mechanical change of the
  tests (one commit, after the layout lanes merge).
- Layouts (interface): a section's relayout still allocates more than before the per-section cache (closures and
  contexts a section; p53 on its gene pan 75 KB a step against 70.5 KB at the base, its 64 KB bound already passed at
  the base on this Node): an engine object a section would take the closures away.

## `@molstar/seqx-viewer`

- Layouts: a leader's cyclic window wraps only where its layout `wraps` (a ring); a tracks leader over the origin (two
  runs, as a follower's) is not drawn, so its window stays within the axis; a click on a window is the ring's own (its
  hit must carry the position under it, and the selection's clicks take the `window` target); the window's slider is a
  mounted element (its arrows and Shift+arrows bound to the `window.*` commands, `aria-valuetext` from `value`); text
  with rows (↑↓ through a line block's rows) waits for the rows layer.

- Legends: the Legend part's place in the default chrome (the host mounts it today); a legend entry's color picker
  (`theme.value-color` has no UI); the scheme picker as a menu (the stories' buttons stand in).

- Addressing: copy link from a context menu (and its key).
- Selection: drag snapping to whole codons over nucleotides (EC-149); rects (dot plots) projected and drag-selected;
  context-menu entries for `selection.add` and `selection.subtract`; the `genome-browser` preset without the add,
  subtract and extend drags (a drag pans there): decided with the presets; a cursor moved off-screen on a follower
  revealing itself by moving the leader.
- Editing: a rebase's `invalid` discards every change of theirs from `partial`, an unrelated one included (by design):
  the editing UI says so.
- Tests: the frame lane's dropped-tick leak is fixed in `@molstar/ui/dom`, so the stories' twins can drop their per-test
  `DomHost.requestAnimationFrame` swap.
- Scroll groups across viewers and views (EC-113); navigators (strip, miniature, ideogram, vertical overview; EC-110:
  the Features page's Navigators row, a design first).
- Rings in the viewer: the matrix cells of a ring section and a block as binding tests; a drag-select longer than half
  the ring selects the other side (the shorter arc is taken: the drag's swept direction would say which); a reversed
  ring's arcs (`selection.ts` reads an arc on a forward axis only, `dir === 1`).
- Touch pinch (two pointers); the Mol\* structure port; the `@molstar/ui` keymap and command manager, replacing the
  shortcut match scoped to the container; `ParamDefinition`s for the config items.
- An allocation tier for the viewer (`test:perf`: the Canvas2D draw, `store.overlay`, `_viewport`, `hover.set`, the
  bindings' resolve); an allocation check by hand of the canvas draw and hover loops.
- Browser measurements: canvas frame times (draw included) and first paint; the Gene by region story rechecked in a
  visible browser (frames do not run in a hidden pane).

## `@molstar/seqx-examples`, `@molstar/seqx-storybook`

- The Plasmid and the Complex are stand-ins (generated letters at the entries' lengths, features at the published
  layout's places, named as such): pUC19c and 1FIN recorded with their dataset lanes replace them, the ring's record
  written again.
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
