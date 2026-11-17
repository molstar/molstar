# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.0.0` until the first slice runs). Why: one language, one
  description, one number for readers to check.

## Language

`lang/`: `SeqLanguage`, the description of SeqViewSpec on `@molstar/core/repr` (kinds, parents, params, docs), and the
checks it cannot express.

- **Kinds** (first slice): `spec` (root: `version`, `title`); the containers `data` and `views` (reserved ids of the
  same name); under `data`: `source`, `sequence`, `ranges`, `values`, `map`; under `views`: `view`, then `section`,
  `track`, `layer` (ids optional on a layer only). Reference params (`target`, `of`, `from`, `to`, `follow`, `source`,
  `data`) are `ref` params naming the kinds they take; a reference is an address.
- **Positions** are 0-based, ranges half-open. A range with `start > end` wraps the origin where its sequence is
  `cyclic` and is an error elsewhere; a section's `range` follows the same rule. Strands are `+` or `-`. Exons of a
  codon map are in genomic order and do not overlap (`strand` says which way the protein runs). Why: one convention for
  every coordinate the language writes.
- **Gate order of `checkSeqSpec`**: (1) the version: a newer major is `seq.lang.version.unsupported` alone, since this
  version's description says nothing of a newer spec; (2) the description's validation with normalization (`repr.*`,
  columns typed, so a spec read from JSON is checked as one built); (3) the language's checks (`seq.lang.*`), only when
  step 2 found no error, so they read a valid, normalized tree. Why: each gate assumes what the one before settled.
- **What the checks cover**: one of two params (`letters` or `source`, `rows` or `source`, `values` or `source`, `uri`
  or `value`), an `index` only for a file, letters of the alphabet, row intervals, bounds and strands (the issue's path
  carries the row), what a map method needs, exons, a codon map's alphabets, section ranges, that a reference names a
  node of a kind its param takes, that a layer's data is what its type draws (reported once: a kind the param does not
  take is the reference's issue alone), that `follow` chains end. Extents are checked where the letters are inline.
- **Open**: the theme slot of a layer's `color` and a query in its `data` wait for the sublanguages (TODO).

## Language builder

`createSeqSpec(options)`: the way most specs are written; `build()` is the normalized tree, `check()` is `checkSeqSpec`
of it.

- **Handles are addresses, branded by kind**: `Address<'sequence'>` is the string `/data/p53` that says what it names.
  Every call that adds a node with an id returns its address; a `layer` without an id has none and returns `undefined`.
  Section, track and view builders carry their `address` and a `toString`, and stand where an address of their kind
  does.
- **A reference param takes** a handle (or builder) of a kind its description allows, or an address written out (any
  string); a handle of another kind is a type error. Which params are references is read from the description, in the
  types and when the builder converts a handle to its address. Why: no list of param names to keep in step with the
  language.
- **Misuse** that types do not reach (an id twice under one parent, a malformed id) throws where it is written; the rest
  is found by `check`. Calls are recorded and replayed on the engine of `@molstar/core/repr/builder` at `build`.
- **Type tests** (`_test/types.test.ts`) mark wrong kinds and params with `@ts-expect-error`; they are checked by
  `tsc -p tsconfig.test.json`, which the package build does not run.

## Issues

`seq.lang.*` codes (`LangIssues`), each with a severity, a description and a hint; a test fails on a code a source uses
and the catalog lacks, and on one cataloged and never used, and a misuse for each code is in `_test/lang.test.ts`. The
namespaces are `repr.*` (the tree core's), `seq.lang.*`, `seq.query.*`, and the later `seq.resolve.*`, `seq.source.*`,
`seq.store.*`, `seq.view.*`, `seq.render.*`, `seq.address.*`. An issue's `path` is the step path and its `address` the
node's, where it has one.

## Queries

`src/query/` (`@molstar/seqx-core/query`): locations, element keys, the three result types, and the query forms as
`repr` form resolvers over a model of systems, numberings and annotations. Forms and addresses are one language: a form
is a segment of an address, `/data/p53/(range :from 94 :to 292)`.

### Result types

- **Nodes** (exact addresses of spec nodes and of nodes data makes), **Elements** (keys of what data makes: an
  annotation's items) and **Location** (positions on systems). A form gives one of them; they travel as `repr`'s `Named`
  (`{ kind: 'nodes' }`, `{ kind: 'other', type: 'elements' | 'location', value }`). A `Named` whose value lacks its
  type's shape is no query result (`fromNamed` checks it). Why: a hover, a click and a layer each hold the type they act
  on, and no form guesses a conversion.
- **Conversions are forms**: `(extent :of <elements>)` gives a Location, `(items :overlapping <location>)` Elements. A
  form given another type is `seq.query.form.type`.

### Location

- `{ parts: { system, ranges: Int32Array }[], rects: { a: part, b: part }[] }`: one part per system (sorted by address),
  its ranges sorted, merged (adjacent ones too), half-open `[start, end)` index pairs, none empty; a cyclic system's
  wrapped range is two ranges. Rects pair two parts (two axes), sorted, each once.
- **A value**: built normalized (`location`), frozen, never mutated. Equality is structural (`locationsEqual`): equal
  parts and an equal list of rects. Rects are never merged, so two rects covering the cells of one are another value.
  `hashLocation` is equal for equal Locations, cached per value.
- **Positions are int32**: `normalizeRanges` throws on an end that is not an int32 integer (a programmer error, never a
  wrap); the forms check data first and report `seq.query.position.invalid`.
- **Set operations** per system: union, intersection, subtraction. Rects intersect side by side over the same two
  systems; a rect is not subtracted, and its complement is refused (`seq.query.form.unsupported`).
- **Overlap** of an item `[start, end)`: a boundary (zero length) at `p` overlaps `[s, e)` when `s ≤ p < e`; one at the
  system's end (`p = length`) overlaps a range ending there, so `(all)` holds it. O(log n).

### Element keys

- `{ owner, key }`: the owner's canonical address and the data's own key, a string or a number. A number never equals a
  string (`1` is not `"1"`), as in an address. Elements are sorted (owner, numbers before strings) and each once;
  compared and hashed as values.

### Numberings and labels

- **Indices** are 0-based half-open inside; **labels** are inclusive names at the edge, converted once where a form is
  read.
- **Built-ins**: `index` (0-based) and `ordinal` (1-based), on every system; without `:numbering` a form reads
  `ordinal`. A declared numbering is looked up on the system, then on its molecule (`parent`, shifted by the unit's
  start, clipped to the unit): a unit's own beats its molecule's of the same id. Why: the author's numbering of a chain
  wins over the entity's. Not across a codon step: that waits for paths.
- **A label** is an integer with an optional suffix, letters (`94A`) or `.n` (`"100.1"`), or else a name. One total
  order: by integer; no suffix, then letters (shorter first, then alphabetic ignoring case), then `.n` by n; names after
  every number; spelling breaks the last ties. A `.n` label is written quoted: bare, `100.1` is the number 100.1 (the
  atoms' rule) and is refused with a hint to quote it; quoted, `"100.10"` and `"100.1"` are distinct labels. An integer
  outside int32 is refused (`seq.query.label.int32`). A repeated label names its first position.
- **Insertion codes**: `(at :label 94 :ins A)` is `94A`; `""`, `?`, `.` mean none; a code on a range, or on a label with
  a suffix, is `seq.query.label.insertion`.
- **Resolving a range**: each end is its label's position; a missing end snaps inward (a start to the first position
  whose label is at or after it, an end to the last at or before it, `seq.query.label.snapped`, info). Then one rule in
  index order for exact and snapped ends alike, never assuming labels rise with the index (IMGT): `from ≤ to` is the
  range; ends snapped past each other over positions with no label are empty (`seq.query.range.empty`); otherwise a
  cyclic system wraps as two ranges short of a whole lap (`seq.query.range.lap`), and any other is
  `seq.query.range.backwards`. A missing `:from` or `:to` is the system's own end.
- **Bounded label work**: a numbering may answer `snap` and `nearest` itself (`index` and `ordinal` do, in O(1));
  otherwise the positions are scanned, each a visit (signal, deadline). A hint lists the labels nearest an unknown one,
  kept a few at a time (no list of every label), and none past `HintScanLimit` positions.

### Forms

| Form                                           | Prefix            | Result                                                                                                          |
| ---------------------------------------------- | ----------------- | --------------------------------------------------------------------------------------------------------------- |
| `(range :from l :to l :numbering n)`           | systems           | Location                                                                                                        |
| `(at :label l :ins c :numbering n)`            | systems           | Location (one position)                                                                                         |
| `(all)`                                        | systems           | Location                                                                                                        |
| `(rect :a <location> :b <location>)`           | as its sides need | Location (one rect; each side one system; sides never wrap, `seq.query.range.wrap`)                             |
| `(item :id k)`, `(item :<field> v …)`          | annotations       | Elements (every param matches; a list param matches any of its atoms; a field the annotations lack is an issue) |
| `(items :overlapping <location>)`              | annotations       | Elements (every item without `:overlapping`)                                                                    |
| `(node :key k)`, `(chain :id k)`               | nodes             | Nodes, addressed by the form itself (`/data/1tup/(chain :id A)`), so location forms may follow a chain          |
| `(extent :of <elements>)`                      | any               | Location                                                                                                        |
| `(and :of […])`, `(or :of […])`, `(not :of x)` | any               | the operands' one type                                                                                          |

- **Combinators take `:of`**: a form's params are keywords only (the address grammar), so `and`, `or` and `not` read a
  list or one operand under `:of`. Operands of different types are `seq.query.form.type`.
- **`not`** complements within what the prefix names: its systems' positions, its annotations' items, its nodes'
  children with ids. Dynamic nodes are not in that universe yet.
- **Operands**: a nested form is evaluated against the same prefix (its resolver, else `repr`'s child selection); an
  address from the root. The location of `items :overlapping` is evaluated against the annotations' targets; a part on a
  system no target is reports `seq.query.location.unreachable` (projection is not built), a rect is refused. An operand
  naming a ref is refused (`seq.query.ref.unresolved`): refs are resolved by a later pass.
- **The table**: `queryResolvers(model, extra?)`; an extra kind named as a query form throws (a form's meaning is never
  replaced silently). `evaluateQuery` is `selectAddress` with the table, typed.

### The model the forms read

- An interface (`QueryModel`): `system(address)` (length, cyclic, numberings, parent), `annotation(address)` (target,
  fields, `items(over?)` with what is pending), `node(owner, key)` and `chain(owner, id)` for nodes data makes. The
  static model implements it; tests use a hand-written one.

### Staged and bounded

- **Staged**: over data read by region an evaluation gives `{ result, complete, pending }`, `pending` a Location of the
  regions not loaded; a form unions its operands' pending regions. Across two resolved segments of one address the
  selection keeps the last segment's (the union is `repr`'s, to come).
- **Bounded**: the caller's `limit` bounds the nodes and elements a form gives (more is
  `repr.tree.address.limit_exceeded`, never a truncation); the signal and the deadline are read every `CancelEvery`
  visits, one counter across nested forms and nested addresses (`repr.tree.address.cancelled`, `deadline_exceeded`). A
  Location's size is not bounded yet.

### Issues

- `seq.query.*` (`QueryIssues`): `form.prefix`, `form.param`, `form.type`, `form.unsupported`, `numbering.unknown`,
  `label.syntax`, `label.int32`, `label.unknown` (with the nearest labels), `label.insertion`, `label.snapped` (info),
  `range.backwards`, `range.lap`, `range.wrap`, `range.empty` (warning), `position.invalid`, `location.unreachable`
  (warning), `ref.unresolved`, `node.missing` (warning). A resolver's errors fail the selection; its warnings and infos
  are carried. A test checks the catalog against the codes the sources use.

## Data

`src/data/` (`@molstar/seqx-core/data`): the source protocol, the `seq` format and its runtime, the data store.

### Protocol

- A **format runtime** (`FormatRuntime`: `name`, `formats`, `open(source, { assets, signal })`) opens a `SourceEntry`
  (`address`, `format`, `uri` or inline `value`) to an `OpenedSource`: its `coverage` and `read(request, signal)`.
  `assets` is a narrow port (`readText(uri, signal)`). Why: a runtime is tested with a fake, never the network. A file
  that cannot be read is a `SourceOpenError` with `seq.source.unreadable`.
- **Coverage**: systems (unique names, lengths ≥ 0, alphabet, cyclic, aliases), shapes per listed system (with a ranges
  table's fields), integer levels ≥ 1 as bin widths, ascending, at least one (`[1]` without levels), nodes known up
  front (`checkCoverage`).
- **An answer** is tables by shape and the dynamic nodes the region brings (`{ key, node }`, a spec node of its kind),
  or a failure: `failed` (retried with backoff), `malformed` (fails at once, with a cooldown), `not-served` and
  `too-large` (never retried; too large asks a coarser level). Cancellation rejects with an `AbortError`, never an error
  state; an `AbortError` the store did not ask for is the source's failure.

### Tables

- Typed columns, each table with a `form`: `letters` (`start`, a string or uint8 codes with an alphabet); `ranges` as a
  `list` (`start`, `end` int32, `strand` int8, `id` category, `row` int32, `fields` as `{ name, column }`) or `density`
  (`bin` int32, `count` float32, `category`); `values` as `positions` at level 1 (`value` float32 from `start`) or
  `bins` above it (`min`, `max`, `mean`, `count` float32 from `start`, aligned to the level). Why a `form`: one shape
  has two encodings by level. Why fields as a list: any name, `kind` and `$x` too, in memory and in a file.
- **A row's identity** is its `id`, else its `row`: a range over a tile's edge is in each tile it touches, and readers
  dedupe by it.
- `checkTables` checks forms, column types and equal lengths. `checkAnswer(answer, request)` checks an answer against
  its request: a record, a known failure with an issue, nodes with string keys, only the shapes asked, each table within
  the region, ranges with `start ≤ end` and strands −1, 0, 1, codes within their alphabet, at most `maxRows` rows.
- Cutting (`sliceTable`) keeps letters, positions and bins by offset, and a range or a density bin by its anchor (its
  first position in the region) within the part cut; `joinTables` joins consecutive parts, so each row is held once.

### The seq format

- The protocol written down: a coverage and `{ request, answer }` pairs as a tree value. A file is a `seq` node through
  the JSON codec (typed columns tagged, atoms lossless); inline, the same record. A failure's `kind` is written
  `failure` (a record with `kind` reads as a node). Reading reports the first misfit (`seq.source.malformed`, at its
  path) and never throws. Why a tree value: one lossless encoding, the codec's, and MessagePack with it later.
- **Its runtime** joins the pairs of the request's system (by name or alias) and level that overlap the region, in
  position order, each cut to its part; a gap no pair holds is `not-served`, naming its positions (the format is
  explicit: an empty region is a pair with empty tables). A failure among the pairs is the answer; an unlisted system or
  level is not served. A list without ids gets `row` from its place in the file. Inline data (`inlineSeq`) is one pair
  over the whole system at level 1, its length in the coverage.

### Store

- **Keys**: `<source>|<system>|<level>|<tile>` and `<source>#<key>`; a part holding `|`, `#` or `\` escapes it, so no
  two keys of other parts are equal. A tile is `tileBins` bins (1,024) at the level's width; a want is clipped to the
  system's length when known and snaps to whole tiles; an unlisted level falls back to the coarsest listed no wider
  (EC-170). A want over `maxTilesPerWant` tiles is `seq.store.too_many_tiles` and not asked: work is bounded before it
  starts. A tile is asked for every shape its system serves: the key has no shapes.
- **Sources** are checked on `attach`: a malformed coverage is an issue and the source unavailable. `detach` cancels its
  requests, drops its entries, releases its holds and clears its reports.
- **Entries** are `pending`, `loaded`, `empty` (no rows, no nodes), `failed` (kind and issue) or `evicted` (wanted,
  dropped for the bound), values of a `KeyedState`, so only what reads an entry is told. Each change has a new version
  from one counter; the store's `version` moves with any; an unchanged entry keeps its object. Pending is never empty
  (EC-169).
- **Want sets**: each consumer has one (`wanter()`: `want`, `whenSettled`, `dispose`; the live view's is `want`), the
  whole set it needs now. They merge by best priority, so an export requests and waits without touching the view's wants
  (EC-125). Tiles are deduped, ordered by best priority (view, margin, prefetch) then first want (EC-172), at most
  `maxInFlight` (6) in flight.
- **Leaving the set**: a request in flight is cancelled after `graceMs` (a jitter back keeps it), or at once, oldest
  stale first, when view work waits for a slot (EC-173); a tile not yet asked for or evicted is dropped; an answered one
  stays in the cache; one in backoff or cooldown keeps its timer and its attempts, so going away and back neither resets
  them nor asks a failing source at once.
- **Failures** (EC-174): `failed` is retried `attempts` times with backoff, the entry as it was, then `failed` with a
  cooldown; `malformed` (also a result that is no promise, or an answer that throws when read) goes to the cooldown at
  once; after it a wanted tile is asked again, the entry failed meanwhile (no flicker). `not-served` and `too-large` are
  never retried. Nothing stays in flight after an answer of any kind.
- **Reports**: each cause once per source, at most `maxIssues`, then `seq.store.limit_exceeded`.
- **Eviction**: the loaded bytes are bounded (`maxBytes`, 64 MB); least recently used first among what is not wanted;
  still over, the wanted prefetch then margin tiles become `evicted` and are not asked again until they leave the set
  and come back, so a planner that wants again on every change does not loop. Never the view, never a request on its
  way. Still over: `overBudget`, reported once each time (`seq.store.over_budget`); a coarser level is the planner's.
- **Settle**: no wanted entry pending (a failed one waiting to be asked again is settled), kept on every transition;
  `whenSettled(timeout?)` resolves `false` on a timeout.
- **Nodes** are held, not cached: by each data entry that carried them, by their source's coverage while attached, and
  by `hold`; with no holder they go quietly. A re-answer takes in its nodes before releasing the old ones. The same key
  with other content keeps the held node (`seq.store.node_conflict`).
- Time is a `StoreClock` (real timers by default): tests drive a fake one and never wait.
