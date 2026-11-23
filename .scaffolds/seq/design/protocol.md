# The source protocol, the store's keys and request planning

## Intent

The second review found the source protocol without types, two designs disagreeing on what would settle it, the
store's keys unable to name a system or a nested source, no plan for requests that start on one system and load from
another, and node lifetime colliding with the cache's eviction. Inline data already goes through the protocol in the
first slice's step 3, so the types are needed first. This design gives them.

## Settled

Accepted as recommended (the author, 2026-11-17, at the first slice's Gate 0): tiles of 1,024 bins, the shapes' column
names and levels as bin widths, exercised by the slice's synthetic source and the `seq` format before FASTA and bigWig.
Range reads in `AssetManager` wait for bigWig; the slice reads its small FASTA whole.

## Decision

### Runtimes

A format runtime is a registry entry ([ports.md](ports.md)):

```ts
interface FormatRuntime {
  readonly name: string;                        // 'fasta', 'bigwig', 'seq'
  readonly formats: readonly string[];          // the `format` values it reads
  open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource>;
}
interface OpenContext {
  readonly assets: AssetReader;                 // the assets port's narrow shape: readText, readBytes (files.md)
  readonly signal: AbortSignal;
  readonly limits: { maxRows: number; maxBytes: number };  // the store's settings; a parse's bound before it allocates
}
interface OpenedSource {
  readonly coverage: Coverage;                  // known after open
  read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer>;
  dispose(): void;
}
```

### Coverage

```ts
interface Coverage {
  systems: { name: string; length?: number; alphabet?: string; cyclic?: boolean; aliases?: string[] }[];
  shapes: { system: string; shape: 'letters' | 'ranges' | 'values'; fields?: FieldInfo[] }[];
  levels: number[];                             // bin widths in positions, ascending; 1 is per position
  nodes?: DynamicNode[];                        // nodes known up front (dynamic-nodes.md)
}
```

- A system's name is the source's own (a chromosome, a transcript); aliases map names a spec may use (`chr1`, `1`).
- **Levels are bin widths** in positions of the system. A format's own levels (bigWig's zoom levels) are listed as they
  are; a source without levels serves `1` only, and the store asks it for per-position data and builds summaries.

### Requests and answers

```ts
interface RegionRequest { system: string; start: number; end: number; level: number; shapes: Shape[] }
type RegionAnswer =
  | { ok: true; tables: Partial<Record<Shape, Table>>; nodes?: DynamicNode[] }
  | { ok: false; kind: 'failed' | 'not-served' | 'malformed' | 'too-large'; issue: Issue };
```

- **Table shapes** (typed columns, [format.md](format.md)):
  - `letters`: `{ start, letters }` (a string, or codes with an alphabet);
  - every table carries a `form` (ranges: `list` or `density`; values: per position or per bin);
  - `ranges` as a `list`: `start`, `end` (int32), `strand` (int8: −1, 0, 1), `id` (category) or else a stable `row`
    (int32, the row's place in its source) so a range sliced across tiles is one range, and `fields` as a list of
    `{ name, column }`; as `density` at a coarse level: `bin` (int32), `count` (float32), by `category`;
  - `values`: at level 1, `value` (float32) per position from `start`; at a coarser level, per bin: `min`, `max`,
    `mean`, `count` (float32) from `start` with the level's bin width.
  - `pairs`, `matrix`, `grid`, `rows`: [region-reads](region-reads.md) (the fifth slice).
- **Nodes** in an answer refer to each other by `(node :key K)`, resolved against the answering source, and to spec
  nodes by absolute address ([dynamic-nodes.md](dynamic-nodes.md)); the language's `ref` has no relative form, so a
  runtime resolves `(node :key K)` against the answering source and writes the absolute address.
- **Errors**: `failed` is retried with backoff then a cooldown, reported once (EC-174); `not-served` is never retried;
  `malformed` is a failure with its issue, never retried (failed with a cooldown); `too-large` asks for a coarser level. Cancellation rejects with an
  `AbortError`, never an error state.
- **The `seq` format** is the protocol written down: a source with `format: 'seq'` holds a coverage and a list of
  `{ request, answer }` pairs, inline (`value`) or as a JSON or MessagePack file. It serves fixtures, tests and the
  synthetic source, and inline spec data is read through it (one region, the whole system).

### The store's keys

- **Data entries**: `(source address, system, level, tile)`, a tile being `TileBins` bins (1,024 by default, a setting)
  at the level's bin width, so a region at a level snaps to whole tiles and overlapping requests share them. Encoded as
  one string for `KeyedState`: `<source>|<system>|<level>|<tile>`.
- **Node entries**: `(source address, key)`, encoded `<source>#<key>`.
- **Nested sources** are addressed through their parent's node, `/data/genes/(node :key T1)/(node :key domains)`, and
  their entries use that address as their source.
- **Inline data** is a `seq` source of one tile per system.

### Request planning across maps

- Each visible layer gives (its section's system, the viewport's range on it, the layer's data). When the data is on
  another system, the viewport's range is projected back along the path (static, or through the overlay) to the data's
  system: the inverse image as runs, each run turned into tiles at the level that fits.
- **A path through a dynamic map that is not held**: the map's own source is requested first (its region on the
  section's system); the layer is pending until the map arrives.
- **The level that fits**: the coarsest listed bin width no wider than the positions a device pixel covers (EC-170).
- **Order and bounds**: view, margin, prefetch (EC-172); at most a set number in flight; cancel what left the view after
  a grace (EC-173).

### Lifetime and eviction

- **Data entries are a byte-bounded cache** (64 MB by default), least recently used, never an entry in view. When what
  is in view alone exceeds the bound, the store reports it once (`seq.store.over_budget`, info) and never asks a coarser
  level itself: a coarser level is the planner's, and the per-layer bounds (`maxTilesPerLayer`, `maxBytesPerLayer`,
  `maxTileBytes`, a region layer reading the first coarser level that fits) keep a view's reads within the cache before
  it is asked (amended 2026-11-23: the code's way, the author as recommended). Revisit if a view is seen over budget with its layers each within their bounds.
- **Node entries are held, not cached**: each counts its holders (data entries that carried it, placements shown, the
  selection, highlights, other held nodes naming it). A data entry's eviction releases its holds; a node with no holder
  is evicted, quietly, with what it alone named. A ref or bookmark in the spec does not hold a node: it waits for it.
- Summaries belong to their attached source, not a data entry, and go when the source is detached ([data-store.md](data-store.md)).

## Consequences

- [sources.md](sources.md)'s provisional protocol is settled here; [data-store.md](data-store.md) takes its keys and
  lifetime from here.
- Tests: the `seq` format round trip; a synthetic source with a fake clock covering levels, tiles, errors, cancellation,
  planning through a static and a dynamic map, and node lifetime against eviction.

## In the specs

| Decision | Where |
|---|---|
| Tiles of 1,024 bins, the shapes' columns, levels as bin widths (Settled) | core › Store |
| Range reads in `AssetManager` wait for bigWig; the slice's FASTA read whole | design-only: history (bigWig reads by range now: core › bigWig) |
| A format runtime, its `OpenContext` (assets, signal, limits) and `OpenedSource` | core › format runtime |
| Coverage: systems with aliases, shapes, levels as bin widths, nodes up front; `[1]` without levels | core › Coverage |
| Requests and answers; the failure kinds and their retries; cancellation an `AbortError` | core › An answer |
| Table shapes and their `form`s | core › Tables |
| The region shapes (`pairs`, `matrix`, `grid`, `rows`) | core › The region shapes |
| Nodes in an answer: refs absolute | core › Checked by the description |
| The `seq` format: the protocol written down, inline or a JSON or MessagePack file | core › The seq format |
| The store's keys (data, node), the tile of `tileBins` bins | core › Store |
| Nested sources addressed through their parent's node | core › Nested sources |
| Inline data a `seq` source of one tile per system | core › Its runtime |
| Planning: the view carried back along the layer's path, the level that fits a device pixel | core › The level rules by shape |
| A dynamic map not held: the layer pending until it arrives | core › Described ends |
| Order and bounds: view, margin, prefetch; in flight; cancel after a grace | core › Want sets |
| Data entries a byte-bounded LRU cache, never in view | core › Eviction |
| Over the bound the store evicts and reports; a coarser level is the planner's, by the per-layer bounds | core › Eviction |
| Node entries held, not cached; a ref does not hold a node, it waits | core › Store |
| Summaries belong to the attached source, go when it is detached | core › Summaries |
