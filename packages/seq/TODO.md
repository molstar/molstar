# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

The first slice built the tree core, the language's first kinds, resolve, draw lists, the SVG and Canvas2D backends, the
viewer's models and data by region; the second, editing and addressing; the third, figures (the axis, themes, the text,
block and ring layouts, navigators, pairs, grids, the logo, the rows layer, the files and their readers, the datasets,
PNG and the render CLI); the fourth, integration (molecules and units in language 0.5.0, structure maps and the
structure link on a fake 3D viewer, seq data ports, extension openings, track header parts, services, bigWig with range
reads of assets, archives, a test consumer of the resolved view). Next, the Genomes milestone (the fifth slice's plan;
built so far: BGZF and tabix in io, language 0.6.0 with assemblies and genome resolve loaded when used, sections over a
source's complex, the genome axis and layouts), on extension kinds.

- **GFF's follow-ups** (GG, 2026-11-21): a name hit focuses its range, not its item address (the item waits for its data
  after open); the GFF row is 25.8 kB (the design's 12 kB was low); by region is FR's.
- **Genome navigators' follow-ups** (GC3, 2026-11-21): bands on the ideogram navigator (GA2, done); the square's cells
  and grips (RG); the label menu's genome entries are static (about 0.6 kB: `extend` refuses late header entries by
  design); Shift+arrows extending a selection across blocks.
- **BED's follow-ups** (BE, 2026-11-21): BEDPE onto the protocol's pairs shape (RG, FR); `itemStart`/`itemEnd` kept
  through `sliceTable` and the store (RG); BED12 span rows by region; region line issues not yet said; the renderer's
  `label` against the design's `name` for a BED name (split convention: decide once in files.md).
- **A component query on a source-backed sequence** (`(elements :component X)`) is pending until the source's letters
  are read (EC-169); a format with a component index could answer it without the letters.

- **Region reads' follow-ups** (RG, 2026-11-22): sided pairs read on side `a` only (side `b` for a track on genome b,
  LA's); matrix state labels in the header, `order` and `hidden`; a pairs bin card's per-category breakdown;
  `MatrixItems` as a ranges list plus a cells column; `store.ts` (854 lines) and `viewer.ts` (818) past 800 with no
  natural split found; held pairs over blocks still binned by device column (Decided 7, partial); the whole-assembly
  level not cut into block pairs (Default 5); inline data over an assembly (GA3); BEDPE through the store (FR).

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
- **MessagePack** (`@molstar/core/data/msgpack`): a duplicate map key keeps the last silently (JSON reports it): an
  option; payload alignment is by position: encoder padding when zero-copy decoding matters more; a canonical-category
  fast path; key-order determinism shared with JSON. Each when a measured use asks.
- **The asset manager**: `maxBytes` is ignored on `file://` reads in Node (`fetchData`'s file path ignores the four
  options): done when the file path honors them and the conformance suite covers it. The suite sees a manager through
  `fetch` only; a manager that downloads by other means is outside it.
- **`util/data-source`**: a `file:` URL's path is opened without percent-decoding, so a root or file with a space or a
  non-ASCII character is unreadable and Windows paths open as `/C:/…` (fails safe). When it moves to `fileURLToPath`,
  the headless policy's file check is reviewed again for Windows segments ending in a dot or space, device names and `:`
  streams.
- **Type checks of core's tests**: core's other `_test` folders have 73 type errors; `tsconfig.test.json` covers
  `data/immutable` only. Done when every `_test` folder is in it and clean.

- `"sideEffects": false` in `@molstar/core` and `@molstar/ui`: about 14.4 kB off the viewer and 2.5 kB off reading a
  spec (immutable's unused barrel members, DOM watch and work). Not now (the author, 2026-11-20): a repository-wide
  change, after an audit of both packages for import-time effects (CSS, registrations) and a build of every app.

## `@molstar/seqx-core`

- **A dynamic (overlay) values node's `field`**: read its source's channel, as a static one does (DP, 2026-11-20).

- **Track heights and labels beside for a source's items** (`render/ranges-content.ts`): a source-backed layer (a GFF3,
  a table, a tiled genome) is not known whole before it is read, so its track keeps its fixed height and has no label
  lane (its labels inside or dropped). A content from what the store holds whole (a small file read in one tile) would
  give it the inline rule.
- **The label column over a template's instances**: measured over the view's own track labels, so an instance's name
  (`Domains T1-domains`) ends in an ellipsis; measuring the instances needs the viewer's data box (`dataX`) to follow
  the frame's column, not the view's.
- **Ring labels inside their arc** when they fit (the linear rule): today a ring places every label outside.
- **Records** (`data/records.ts`): the GFF3 preset and a row past a cyclic system's length wrapped by the row reader;
  the delimited tokenizer reads a string whole: a streaming form when a file past `maxBytes` read by region asks; a
  values file is bounded by `maxRows` positions from its first (a genome's values wait for region reads).
- **A text section's hit allocates** about 70 B past its caller's call (`render/_test/text.quick.ts` reports it; the
  placement's `at` measures 0 on its grouped axes): its other hit paths (the gutter's `hitChrome`, the layers' hits
  through a grouped axis) to read their numbers from arrays, then a bound in that case.
- **The small layers** (built: ruler, secondary, variants, values' line and area, ranges' line styles and arrowheads):
  - channels they do not read yet: `size` (a variant's radius, a line's thickness, a band's share), `opacity` (the
    area's 0.35, a highlight's, a variant's), `shape` (variants' glyph: circle only), `line` and `arrow` by a form (the
    slots take forms; the layers read a word); a values color by a field other than its value;
  - `secondary`'s `style: cartoon` (with the text figure, layouts); a text's numbers in a numbering whose labels are a
    source's (a section layout reads no store: they print ordinals there, the ruler the source's labels);
  - a tiled source's items have no arrowhead (their other rows may lie outside what was read): a source that answers an
    item's extent would give them one; the three layers' `nextMark` and hover passes are still their own;
  - `arrow` points by an item's own start and end, not by its strand: a strand column (EC-232 asks arrows only when
    asked) would say which way an item runs;
  - a ruler's labels over a numbering that is not monotone (IMGT) in a piece of more than 4,096 positions are looked up
    by number and skipped when the numbering jumps over 16,000 labels a piece;
  - a line over a reversed collapsed axis in the device-column step does not break at the axis's breaks (the columns are
    contiguous).

- **Addressing**: an item's coordinate system (`<node>/(item :id k)`) as a prefix of location forms (the `item` form
  gives Elements today, so `/refs/t1/(item :id T1)/(range …)` fails); rects projected for a focus; a whole Location's
  size bound.

- **Language**: a source's labels keyed by something other than a position (a numbering's `offset` and `segments` are
  built, MK); a structure map's `unpaired` status written (derived from the runs today); a `segments` or `structure` map
  from a source (inline only); `via` on the overlay's own item maps (a spec map's rivals only); a union of a ref and a
  query slot for `layer.data` (loose today); a strand column check; codon exons checked against the protein's length
  (stop codons); a section's `range` as a Location; refs followed in data nodes' params and a section's `of` and
  `follow` (a layer's `data` follows them); a `label` param on nodes (an instance's label is the template's plus its top
  ancestor's key); `seq.lang.ref.target` overlaps what resolve reports: keep one; the package index re-exports
  `SeqLanguageVersion` as `LanguageVersion`: one definition.
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
  through the update rule (EC-111) on the way; a jump's motion tweens no band (a step's does: the model's bands are its
  start); story steps that fold clades (a step's view holds no folds: grids-and-rows, "Opening and closing clades"); the
  story's own UI (steps listed, edit this step, store the view as a step) and the render CLI's video and frame-sequence
  outputs; an op that no longer applies has no issue in `repr/patch`, so Seq reports it on the step; a host's move can
  leave a step's later refs dangling in the step's own coordinates (refused unless `allowDangling`): watched for in the
  editing UI.
- **Bounds**: regular expressions in search and archives get theirs with those features; nested-source breadth (a
  fan-out of W^depth) and joins allocating before `checkAnswer` (at most the file's size) are recorded, not bounded.
- **Overlay and templates**: re-expand only the template an arrival touches, keep instance layers across overlay
  versions, key the overlay on what it reads of the model (a view-only edit rebuilds it); a static layer whose `data`
  names a dynamic node waits for it instead of reporting `seq.view.layer.data`; highlights hold their nodes (the
  selection does, `selection.holdNodes`); the overlay's sync issues in the viewer's `layoutIssues`.
- **Data**: `summaries` in `ViewerConfig`; summaries for systems of unknown length or past `maxLength`; byte dedupe
  across tiles; `SourceEntry` carries only `address`, `format`, `uri`, `value`.
- **Render**: a summary pyramid for inline values zoomed out (inline data is level 1 only, so a large one zoomed out
  hits the tile cap); a legend's place in the chrome; palette registry entries (layer types and schemes are an
  extension's, atomic and undoable at run time); `query` layers read the whole annotation (filtering comes with the
  query sublanguage); `values` line and area, points and split; the ranges layout's allocation an item in view (the
  `earlier` and `occurrences` sets, `StickyLanes.keep`'s map, keys of tiles without ids, `placeholderSpans` and
  `dataRegions` arrays: about 300 B an item); the values direct path for linear maps (offset paths) too; the renderer's
  rebuild sorts the reads and builds the key's text each layout (bounded by tiles in view); `render/export.ts` still
  requests data for hidden layers (done when an export of a spec with a hidden layer reads nothing for it); the ranges
  layout walks every row of the entries it reads (a 150,000-row inline annotation: about 114 ms a frame); a letter
  through a codon path allocates its parts (`projectRanges` sets its buffer's length to 0, which lets the storage go:
  161 KB a pan frame on the coding gene).
- **The rows layer** (built: the rows half of grids-and-rows; each item here with its design section):
  - the keyboard cursor across text line blocks (its rows, part and a grid's y are built, RK: "Selection");
  - block drag-select (rects) and Shift extension over rows by drag ("Selection");
  - a bound row's projection from its sequence; `(row :key k)` and `(clade :of [a b])` as query forms ("The rows axis",
    "Selection");
  - `filter` and `rows.filter` (by a field's values or range: "Order");
  - the logo of the rows shown (`rows: shown`, taken out of the language until built: "The logo"); a small logo as a
    summary row (`summary: logo`, taken out likewise: "Opening and closing clades");
  - a reduction by letter beside the theme's classes (the reducer is fixed, decided by the author 2026-11-20: an item
    `rows.reduce` only with a second use: "Drawing by density") and a lazy source's served row bins, coarse levels drawn
    for 10M rows at the fit (a neutral fill today: "Lazy rows");
  - a drag of the distance cut line across the tree (`clades.collapseTo { height }` is built); a reveal focusing a
    folded clade's summary row instead of opening it (a reveal always opens, fixed: decided by the author 2026-11-20; an
    item only with a second use) ("Opening and closing clades");
  - a toggle's motion (the rows' heights animated as a jump's: `view.motion`, `view.jumpMs`); the drag past the band's
    edge scrolling the rows ("Opening and closing clades", "Lazy rows");
  - accessibility: a tree node announced with its state through the canvas's accessible layer ("Opening and closing
    clades");
  - icons in the header's fixed slots (pin, a bound row's link) and the rows table's data columns ("The rows axis");
  - a rows track's header line as its title row (`rows.titleRow: auto`): until then a rows track draws no states, no
    slots and no label hit (core's "Track headers"); the header line in the text layout's first line, a grid's label and
    a ring's menu (the core lays it out in the tracks layout only, its blocks' first row included);
  - the track headers' `seq:pending` state (`loader-circle`, data pending in view: no per-track pending fact yet), a
    run's message said by the live region, a header entry's widget slot (`kind: 'widget'`, when a host asks), a template
    instance's hide (the template's);
  - the renderer's section cache keyed by the export's `folds`; the cell hit's 16 B (a number boxed from the point to
    its column through the axis, every layer's hit) measured by the allocation tier;
  - a Design page (the Pfam globin seed and its tree are recorded: the seed as published, the tree by neighbor joining,
    inline);
  - no form names an alignment's rows by position; `(row :key …)` on an alignment read from a source names a row only
    once its rows are held (a rows layer read them): a row form does not ask for them;
  - a heme per heavy atom (4HHB's hemes are one position each this slice) and the map from 4HHB's units onto the seed's
    columns (a `segments` map through the closest row, as the tRNA's chain; the `alignment` map kind, the gapped row
    itself with no runs written, is not built); pairs on an alignment's columns are drawn on its sections only (a
    chain's section says `seq.render.pairs.systems`: pairs are not projected);
  - the text layout's `chrome` throws when an alignment's source has not opened (`render/layout/text.ts`, its
    placement's `marks` undefined: a figure over a seed whose runtime is missing).
- **Rings and blocks** (built: the rings layout with labels, the title and the window; blocks in tracks and on rings):
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
- **Wells and padding**: a template's tracks take the template track's `padding` as written; whether each instance
  should have its own (a param of the dynamic node, a per-template default) is undecided. The column-shaped tints
  (`selectionShape: 'column'`, the default) span the tracks and the gaps between them and are not clipped per well; with
  `selectionPlace: 'under'` they show in the gaps and past the wells' round corners (per-well clips would split the one
  band into a band a track). A band's hover is square at the view's edge by its mark's place, not by a cut the layout
  records (a band there may end exactly at the edge).
- **Tests**: the p53 fixtures round-tripped JSON ↔ MessagePack in a test here (io cannot import them).

- **Themes**: a filter by a value's range (steps, a ramp: `(item :<field>)` compares atoms only) and filters through a
  tiled source's server; the legend's menu entry for `layer.filter` (the command and the chip are built); the pLDDT
  view's data is a seeded stand-in (a model's confidence from a real fixture); a hover halo on the legend's entry in the
  canvas (EC-138: the legend highlights, the marks do not ring); the hover outline's double outline (light and dark)
  over a saturated fill under 3:1 against it (EC-138, decided by the author 2026-11-20: 3:1 is against the track); a
  values layer's y on a `log` scale (its color reads one, its bars are linear); a `field` form's `:data` (a field of
  another annotation); forced colors (`forced-colors: active`: the system's palette); the column measures (`share`,
  `identity`, `matches-reference`) registered by an alignment's layers; `classes` (the default protein scheme) is not
  safe under the three deficiencies (its substitute is drawn under `theme.cvdSafe`): a safe default, a design decision;
  the dark band (L 0.30–0.62) reads dim against the dark page (fills 1.4–1.7:1 against it at the low end) and cannot
  rise without giving up 4.5:1 ink: an outline or a lighter page, a design decision; an inline table's column named
  `kind` is refused by the s-expression codec (a reserved key), so `secondary` reads `category` in the examples; the
  `theme/` entry is reached through the `./*` export (`theme/index`): its own `./theme` export.

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
- Pairs and grids: pair labels under an apex (the `label` column is read, not drawn); kinds at one cell as glyphs with
  "+" when they do not fit, if a figure asks (EC-63: slices decided by the author 2026-11-20; arcs do not share);
  contacts on copies folded by copy group (EC-62, EC-65); a grid's image tiles assume a one-piece axis (a grid over a
  section with breaks draws rects or stretches the tiles); `across_lines` in text, and arcs on a stack of tracks
  sections read as lines; pairs and grids are held in the model (`resolve/held.ts`), not through the seq format and the
  store (read by region waits for the Genomes milestone); `split` with the letters between its two sides (one band: the
  split arcs stand on its middle).

## `@molstar/seqx-viewer`

- **The minimal viewer, what still holds it** (viewer-composition "Injected, opt in"): the tracks layout reads the rows
  layer (`RowsLay`, the row columns of `rows-model`), the text layout (`TextLayout` by identity) and the grid's
  `GridPlacement` and `rectOver` itself, so those stay in a viewer that lists none (about 30 kB minified of the rows
  layer); a layout's name or a layer type's hook in its place drops them. The input's `RingPlacement` check keeps part
  of the rings layout. The minimal viewer's check assumes `sideEffects: false` in the core's and the viewer's
  `package.json` (set at the lane's merge): without it the bundler keeps module-level statements of unlisted modules
  (the navigator's scratch arrays, the story's issue codes, the rows' zoom bounds); the story config's defaults stay
  (the viewer's config reads them).
- **A records source's row issues** (`RecordsSource.issues`: rows left out, clipped, ragged) are not reported: the store
  model keeps one issue a source address, its opening's refusal; the opened source's warnings to join `issues`.
- Layouts: a leader's cyclic window wraps only where its layout `wraps` (a ring); a tracks leader over the origin (two
  runs, as a follower's) is not drawn, so its window stays within the axis; a click on a window is the ring's own (its
  hit must carry the position under it, and the selection's clicks take the `window` target); the window's slider is a
  mounted element (its arrows and Shift+arrows bound to the `window.*` commands, `aria-valuetext` from `value`); text
  with rows (↑↓ through a line block's rows) waits for the rows layer.

- Legends: the Legend part's place in the default chrome (the host mounts it today); a legend entry's color picker
  (`theme.value-color` has no UI); the scheme picker as a menu (the stories' buttons stand in).

- Addressing: copy link from a context menu (and its key).
- Grids: a grid's own y window in the history and sessions (kept by the view and drawn, not yet a step's); a drag on a
  grid's ruler panning that axis alone.
- Selection: drag snapping to whole codons over nucleotides (EC-149); a rows block selected as a rect; context-menu
  entries for `selection.add` and `selection.subtract`; the `genome-browser` preset without the add, subtract and extend
  drags (a drag pans there): decided with the presets; a cursor moved off-screen on a follower revealing itself by
  moving the leader.
- Editing: a rebase's `invalid` discards every change of theirs from `partial`, an unrelated one included (by design):
  the editing UI says so.
- Tests: the frame lane's dropped-tick leak is fixed in `@molstar/ui/dom`, so the stories' twins can drop their per-test
  `DomHost.requestAnimationFrame` swap.
- Scroll groups across viewers and views (EC-113).
- Navigators, after the first build: a row zoom dragged on an overview's ends (the rows axis's ends only pan); the rows'
  one-value-a-row fallback (`rowsValue`) and then the rows' structure where a source serves no row summary (a lazy
  alignment's rows are left out now); a navigator over a text section (the positions on its lines in view); the window's
  hover readout as an ARIA live region; a slider per covered axis (one per window, roving by PageDown and PageUp); a
  grid's own y window drawn on a navigator (the view's history keeps it; open until a Hi-C view); the rows rail as a
  slider for a screen reader (a canvas mark now: `]`, `[` and the wheel are its keys); `navigator.idle: hidden` (fixed
  quiet, decided by the author 2026-11-20; idle is quiet only, no item: until a host asks, navigators.md).
- Rings in the viewer: the matrix cells of a ring section and a block as binding tests; a drag-select longer than half
  the ring selects the other side (the shorter arc is taken: the drag's swept direction would say which); a reversed
  ring's arcs (`selection.ts` reads an arc on a forward axis only, `dir === 1`).
- Touch pinch (two pointers); the `@molstar/ui` keymap and command manager, replacing the shortcut match scoped to the
  container; `ParamDefinition`s for the config items.
- An allocation tier for the viewer (`test:perf`: the Canvas2D draw, `store.overlay`, `_viewport`, `hover.set`, the
  bindings' resolve); an allocation check by hand of the canvas draw and hover loops.
- Browser measurements: canvas frame times (draw included) and first paint; the Gene by region story rechecked in a
  visible browser (frames do not run in a hidden pane).

## `@molstar/seqx-examples`, `@molstar/seqx-storybook`

- **Watched at the gates' runs**: `@molstar/ui`'s `playback.test.tsx:112` failed once in a combined run (FX1's; 12 runs
  green alone); the VCF ClinVar story test (`vcf-stories.test.ts`) failed once under load during GB's run (2026-11-22),
  then passed alone. The services story test (fixed by FX8) and Prism's `commands.test.ts` and `collapse.test.tsx` (6 of
  6 green, FX1) are off the list. A second failure is a fix lane's item, with the load it ran at.

- The Plasmid and the Complex are stand-ins (generated letters at the entries' lengths, features at the published
  layout's places, named as such): pUC19c and 1FIN recorded with their dataset lanes replace them, the ring's record
  written again.
- The other datasets and their views (Antibody, Plasmid, tRNA), each carrying a drawing-against-the-data test as the p53
  views do (`drawn-vs-data.test.ts`); the Complex dataset (CDK2 with cyclin A, its interface, PAE, ligands with
  AlphaFill's transplants) with pairs across chains, dot plots and ligands, its entries chosen and checked against their
  sources then.
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
- **The ring's site labels fan out near a cluster** (the plasmid's restriction sites crowd the MCS, their leaders a
  narrow fan, the labels a line apart on one side): a placement that spreads a cluster over both sides or stacks it in
  columns; the novice pass found the names hard to follow.
- **Long legends do not wrap** (the tRNA and complex views' legends run off a line in a narrow view): wrap the entries
  by the available width, a group a line when it fits, and keep the pointer's hit boxes on the wrapped entries.
- **The ring's lanes carry no names of their own**: the strands are named by the legend's words and the detail's track
  labels; a lane's label drawn on the ring where no feature crosses it would name them in place.
- **The grid rows' highlight stripes allocate a few bytes a run** (`grid.ts` `emitRows`, the quick case in
  `pairs-grids.quick.ts` bounds it at 128 B): find the boxed value; the tracks' stripes allocate nothing.
- **Host commands cannot declare their change kind** (structure-syncing: a host's color command changes no Seq\* state):
  `CommandChangeKinds` has no hook for a host's ids.
- **A gap mark on the rows layer** (FX4): the rows layer skips an alignment's gap cells, so the Antibody's against view
  shows its gaps faintly at the fit; a `gap_marks` param (a gap list in `RowsLay`, canvas and SVG, the image-tile step
  below 3 px a cell).
- **The fourth slice's novice pass, small items** (novice-findings.md): the fake's log time-stamped (F2); a selection
  tint hiding a ruler tick label (F4); a tick at the focus at high zoom (F6); one of the two Color by this track stories
  goes (F10); short tracks keep the layer's colors and a mark's hit target grows (F16); the host entry's copy slot shows
  it ran (F19); a chip's count by width, a swatch outline on a light ground (F24); the keys panel scrolls or opens
  upward (F26).
