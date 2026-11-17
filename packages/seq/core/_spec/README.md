# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.0.0` until the first slice runs). Why: one language, one
  description, one number for readers to check.

## Data

`src/data/` (`@molstar/seqx-core/data`): the source protocol, the `seq` format and its runtime, the data store.

### Protocol

- A **format runtime** (`FormatRuntime`: `name`, `formats`, `open(source, { assets, signal })`) opens a `SourceEntry`
  (`address`, `format`, `uri` or inline `value`) to an `OpenedSource`: its `coverage` and `read(request, signal)`.
  `assets` is a narrow port (`readText(uri, signal)`). Why: a runtime is tested with a fake, never the network.
- **Coverage**: systems (name, length, alphabet, cyclic, aliases), shapes per system (with a ranges table's fields),
  levels as bin widths ascending (`[1]` without levels), nodes known up front.
- **An answer** is tables by shape and the dynamic nodes the region brings (`{ key, node }`, a spec node of its kind),
  or a failure: `failed` (retried), `not-served` (never), `malformed` (a failure), `too-large` (a coarser level).
  Cancellation rejects with an `AbortError`, never an error state.

### Tables

- Typed columns, each table with a `form`: `letters` (`start`, a string or uint8 codes with an alphabet); `ranges` as
  `rows` (`start`, `end` int32, `strand` int8, `id` category, `fields` columns) or `density` (`bin` int32, `count`
  float32, `category`); `values` as `positions` (`value` float32 from `start`) or `bins` (`min`, `max`, `mean`, `count`
  float32 from `start`, at the level's width). Why a `form`: one shape has two encodings by level.
- `checkTables` checks forms, column types and equal lengths; the store treats a misfit as `malformed`.
- A table cut to a region (`sliceTable`): letters, positions and bins by offset, rows by overlap (a range over a tile's
  edge is in each tile it touches; readers dedupe by id), density by bin index.

### The seq format

- The protocol written down: a coverage and `{ request, answer }` pairs as a tree value. A file is a `seq` node through
  the JSON codec (typed columns tagged, atoms lossless); inline, the same record. A failure's `kind` is written
  `failure` (a record with `kind` reads as a node). Reading reports the first misfit (`seq.source.malformed`, at its
  path) and never throws. Why a tree value: one lossless encoding, the codec's, and MessagePack with it later.
- **Its runtime** answers from the pair of the request's system (by name or alias) and level that holds the region, else
  one that overlaps it, cut to the region and the shapes asked; a region no pair touches is empty; an unlisted system or
  level is not served. Inline data (`inlineSeq`) is one pair over the whole system at level 1.

### Store

- **Keys**: `<source>|<system>|<level>|<tile>` and `<source>#<key>`; a part holding `|`, `#` or `\` escapes it, so no
  two keys of other parts are equal. A tile is `tileBins` bins (1,024) at the level's width; a want snaps to whole
  tiles; an unlisted level falls back to the coarsest listed no wider (EC-170). A tile is asked for every shape its
  system serves: the key has no shapes.
- **Entries** are `pending`, `loaded`, `empty` (no rows, no nodes) or `failed` (kind and issue), values of a
  `KeyedState`, so only what reads an entry is told. Each change has a new version from one counter; the store's
  `version` moves with any. Pending is never empty (EC-169).
- **`want(wants)`** is the whole set needed now: tiles deduped, ordered by best priority (view, margin, prefetch) then
  first want (EC-172), at most `maxInFlight` (6) in flight. A tile that leaves the set is cancelled after `graceMs` if
  in flight (a jitter back keeps it; the oldest stale goes first, EC-173), dropped at once if not yet asked for, and
  kept in the cache if answered.
- **Failures** (EC-174): `failed` and `malformed` are retried `attempts` times with backoff, the entry pending, then
  `failed` with a cooldown, after which a wanted tile is asked again (the entry stays failed meanwhile: no flicker);
  `not-served` and `too-large` never. Each cause is reported once per source (`issues`).
- **Eviction**: the loaded bytes are bounded (`maxBytes`, 64 MB); least recently used first among what is not wanted,
  then prefetched, then in the margin; never the view or a request on its way. Still over: `overBudget`, reported once
  each time (`seq.store.over_budget`); asking a coarser level is the planner's.
- **Settle**: no wanted tile queued, in flight or in backoff; `whenSettled(timeout?)` resolves `false` on a timeout.
- **Nodes** are held, not cached: by each data entry that carried them, by their source's coverage while attached, and
  by `hold`; with no holder they go quietly. The same key with other content keeps the held node
  (`seq.store.node_conflict`).
- Time is a `StoreClock` (real timers by default): tests drive a fake one and never wait.
