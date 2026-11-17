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

## Consequences

- Tier 1 (the static model and the resolved view) never reads the store; node entries feed the overlay (tier 2); data
  entries feed the draw lists and legends (tier 3) ([composition.md](composition.md)).
- The store's conformance tests: pending vs empty, dedupe, cancellation, eviction never in view, versions, settle and
  open domains, with a synthetic asynchronous source and a fake clock.
