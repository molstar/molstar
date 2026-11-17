# Composition: the spec's levels and the compiled forms

## Intent

The author, 2026-11-17: "expand on how the "spec" will be composed? the views etc .. similar to how meso composes many
concepts. and do we need intermediate representation the spec compiles/expands to? i think this should streamline many
things, including mapping the "streaming ports" etc?"; after the review: one version, and the architecture reconciled
("this is not just about scope, but about reconciling the architecture").

A spec is short and declarative. What the viewer, the source runtimes and other renderers need is explicit: refs
followed, defaults filled, maps resolved, requests known. Between them sit three compiled forms.

## Precedents

- meso\*: levels in one tree, slots embedding sublanguages, definitions and uses, IRs compiled by pure passes cached by
  hash and keyed by address, streamed data never written back ([sources/meso](../sources/meso/README.md)).
- The previous implementation: `data` and `view` under a root; `compile` and `resolveView` as two passes; streamed
  sources resolved against a manifest ([sources/seqstar](../sources/seqstar/README.md)).

## Decision

### One language, four levels

```
spec
├── data        sources, sequences (units, complexes, copies, structure chains), alignments, trees,
│               numberings, maps, annotations (the six shapes), computed annotations
├── views       views → sections → tracks → layers; themes; focus; highlights
├── refs        user-defined names for addresses (bookmarks, regions, sets)
└── snapshots   named patches over the base
```

- **One language, one version** (the author: "one"); the levels are container kinds with reserved ids, so addresses
  read `/data/p53`, `/views/main`, `/refs/dbd`, `/snapshots/step3`. A small spec writes `data` and `views` only.
- **Data kinds sit under `data`** (a source, a sequence, a map, an annotation) and refer to each other by address;
  their order is free.
- **Sublanguages in slots**, described in the same language description: queries (a layer's data, a filter, an
  address's forms), themes (a channel's value), patches (snapshots). Addresses are the query language's forms after a
  path ([addressing.md](addressing.md)).
- **Definitions and uses wait**: themes and styles may be named once in `views` and referred to by address (a plain
  ref); `def` with arguments and overrides come when an example needs them. Track templates over a source's dynamic
  nodes are in ([dynamic-nodes.md](dynamic-nodes.md)).
- **Extension kinds** sit in any level under their namespace, with a fallback ([extensibility.md](extensibility.md)).

### Session is not the spec

- The view shown, the viewport, the selection, the hover, the cursor, open panels: session state, transient, held by
  the viewer's models. Its shape and whether it is saved (beside the spec, never in it) is decided with the viewer's
  models; snapshots and refs are how a spec names a view state on purpose.

### Three compiled forms, in three tiers by what keys them

```
                 keyed by the spec's hashes                       keyed by the store's versions
spec ─read─▶ checked spec ─resolve─▶ static model ─▶ resolved view │ overlay ─▶ draw lists ─▶ backend
                                                                   │    ▲          ▲
                                                 data store ◀── runtimes ◀── requests (resolved view + viewport)
```

- **Tier 1, pure, keyed by the spec's hashes**:
  - **the static model**: every coordinate system the spec declares and its numberings; units and copy groups; maps
    the spec declares, resolved to their encodings, with paths precomposed ([maps.md](maps.md)); each source with the
    runtime bound to it (unavailable with an issue when none); each annotation's shape, fields and access plan.
  - **the resolved view** (public, renderer-independent, [rendering.md](rendering.md)): built over the static model
    only. Sections, tracks, layers with themes as values or value functions, fallback chains, and for each layer its
    data as an address or a query and the path from its data's system to its section's when that path is static. It
    holds no data, no dynamic node, no pixels, no fonts, no target.
- **Tier 2, keyed by the store's versions**: **the overlay**: the dynamic nodes held in the data store (systems, maps,
  annotations, nested sources), with the paths through dynamic maps composed when the overlay's version changes, for the
  systems it touched ([dynamic-nodes.md](dynamic-nodes.md)). It is never part of tier 1, so a stream never invalidates
  the resolved view.
- **Tier 3, per frame, keyed by both plus the size, viewport and target**: **the draw lists**: the resolved view laid
  out from the static model, the overlay and the store's data entries; a layer whose data or path is in the overlay
  resolves it here (a pending layer waits, never "unplaced"); a template expands here into a track per matching node;
  the target picks each layer's first drawable fallback. Built into reused buffers, pure by contract
  ([layout-and-draw-lists.md](layout-and-draw-lists.md)).
- **The data store** ([data-store.md](data-store.md)) sits beside the tiers: requests flow from the resolved view and the
  viewport to the runtimes; answers land in the store with versions. Data entries feed tier 3; node entries feed tier
  2. Nothing from the store goes into tier 1.
- **Summaries** (pyramids, indexes) belong to their store entry, built when it loads, evicted with it
  ([performance.md](performance.md)).

## Rejected alternatives

- A language and a version per level, and five IRs (data, model, annotation, view, render), as first drafted: more
  descriptions and boundaries than the work needs; data, model and annotation resolve together.

## Consequences

- Each consumer reads one form: the viewer and exports the draw lists, other viewers (Nightingale, JBrowse) the
  resolved view with the overlay and the store, runtimes their source entries in the static model, tests any of them.
- The core's folders follow the pipeline: `lang/` (the description, the builder), `resolve/`, `view/`, `data/` (the
  protocol, the store, the runtime contract), `render/` (layout, draw lists, the SVG and recording backends).
