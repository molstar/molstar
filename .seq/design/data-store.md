# The data store

## Intent

The review found the pure-pass story broken at the data: data read by region arrives asynchronously, is not in the
spec, and the draw lists depend on it. This design says where that data lives, how it is keyed and versioned, how the
viewer learns it arrived, and what "an open domain" means when data comes by region (accepted as recommended, the
author, 2026-11-17).

## Decision

- **One store per viewer** (or per headless render context), outside the compiled forms. Data entries are keyed by
  `(source, system, level, tile)` and node entries by `(source, key)`, encoded as strings ([protocol.md](protocol.md)).
- **An entry's state** is one of `pending`, `loaded`, `empty`, `failed` (with its issue). Pending is never drawn or hit
  as empty (EC-169); failed shows its issue on the track and is retried by policy (EC-174).
- **Versions**: each entry has a version that changes when its state or data changes; the store has a version that
  changes when any entry does. The draw lists' cache key includes the versions of the entries they read, so a frame
  after an arrival redraws only what read the new entry.
- **Notification**: the store is a reactive model (`@molstar/core/reactive`): `entry(key)` gives a keyed state
  (`KeyedState`), so only what reads an entry is told; the viewer schedules a frame on change.
- **Requests**: from the resolved view and the viewport, each visible layer asks for (source, region, level); the store
  dedupes, orders (view, margin, prefetch; EC-172), bounds what is in flight, cancels what left the view after a grace
  (EC-173), and passes the rest to the source's runtime through the source protocol ([sources.md](sources.md)).
- **Summaries are per attached source**, not per data entry: a coarse tile per entry would need unbounded level-1
  tiles. Factor 16, bounded by `maxLength` (2^24); over it, binned per tile from a bounded level-1 read (or not
  served). The pyramid counts in the byte budget.
- **A detached source's tiles leave the request queue.** Sources opened for the overlay are shared by count on one
  store: the last release closes them.
- **What two consumers share is held per holder**: anything the store or a registry shares between consumers (the live
  view and an export, two views) is counted by holder, so one consumer ending releases only its own hold. Its tests
  include two consumers interleaved, one ending while the other reads. Why: an export ending once detached the live
  view's nested sources, and a detached source left tiles pending; tests with one consumer cannot see either.
- **Eviction**: data entries byte-bounded (64 MB by default, a setting), least recently used, never an entry in view;
  node entries held by their holders, not cached ([protocol.md](protocol.md), Lifetime).
- **Answers carry nodes**: a region's answer holds tables and the data nodes it brings; the nodes are entries of their
  source, identified by (source, key), and form the resolved model's dynamic overlay ([dynamic-nodes.md](dynamic-nodes.md)).
- **The store's clock is readable** (`clock`): an export's deadline runs on it, so a shared store's clock (or a fake one in tests) times the export too.
- **Inline data** goes through the same store, as entries that are loaded at once; one path for drawing.
- **Want sets per consumer**: `store.wanter()` gives each consumer (the live view, an export) its own want set, merged
  by best priority; an export never touches the live view's wants (EC-125). Over the byte bound, wanted prefetch and
  margin tiles may become `evicted` (not asked again until they leave and re-enter the want set); view tiles and requests
  in flight never are.
- **Headless and export**: `renderSvg` is asynchronous: it requests what the export's size and viewport need and waits
  for those entries to settle (or a timeout, after which the export says it is unsettled, EC-125). The live viewer's
  viewport is never changed by an export.
- **Open domains**: a value domain left open is computed from the loaded entries in the view when the view settles (no
  pending requests, no gesture in progress), not every frame; during a pan or zoom the last settled domain holds, so
  colors and heights do not shimmer. An export computes it for its own viewport and records it in the SVG's metadata.
  A legend states the domain it shows.

- **Tooltip tables** ([tooltips](tooltips.md); TT): read by region beside the windows of the sections whose marks they
  decorate, at the view's priority, keyed by the tooltip's address; kept as those windows are; pending shows a
  `loading` row, failed its issue (EC-257).

## Consequences

- Tier 1 (the static model and the resolved view) never reads the store; node entries feed the overlay (tier 2); data
  entries feed the draw lists and legends (tier 3) ([composition.md](composition.md)).
- The store's conformance tests: pending vs empty, dedupe, cancellation, eviction never in view, versions, settle and
  open domains, with a synthetic asynchronous source and a fake clock.

## In the specs

| Decision | Where |
|---|---|
| One store per viewer or headless context; keys by source, system, level, tile and by node | viewer › StoreModel |
| An entry's state; pending never empty; failed shown and retried | core › Store |
| Versions per entry and store; the frame keyed by the versions it read | core › The key |
| Notification: entries a `KeyedState`, the viewer's frame on change | core › Store |
| Requests: deduped, ordered, bounded in flight, cancelled after a grace | core › Want sets |
| Summaries per attached source, factor 16, `maxLength`, binned per tile past it | core › Summaries |
| A detached source's tiles leave the queue; overlay sources shared by count | core › Store |
| What two consumers share held per holder, tested with two interleaved | core › Shared per holder |
| Eviction byte-bounded, LRU, never in view; nodes held, not cached | core › Eviction |
| Answers carry nodes, the overlay's dynamic nodes | core › Overlay |
| The store's clock readable; an export's deadline on it | core › One deadline |
| Inline data through the same store | core › Its runtime |
| Want sets per consumer; prefetch and margin may be `evicted`, view and in flight never | core › Want sets |
| Headless and export: `renderSvg` asynchronous, waits to settle or times out, the live viewport untouched | core › One deadline |
| Open domains computed at settle, held during a gesture, an export's recorded in the SVG | core › Scales |
| Tooltip tables read by region beside the sections they decorate | core › Tooltip tables wanted |
