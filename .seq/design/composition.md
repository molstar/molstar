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
  ref); `def` with arguments and overrides, and templates used per member of a set, come when an example needs them.
- **Extension kinds** sit in any level under their namespace, with a fallback ([extensibility.md](extensibility.md)).

### Session is not the spec

- The view shown, the viewport, the selection, the hover, the cursor, open panels: session state, transient, held by
  the viewer's models. Its shape and whether it is saved (beside the spec, never in it) is decided with the viewer's
  models; snapshots and refs are how a spec names a view state on purpose.

### Three compiled forms

```
spec ─read─▶ checked spec ─resolve─▶ resolved model ─▶ resolved view ─layout─▶ draw lists ─▶ backend
 (decode, migrate,          (systems, maps, paths;      (sections, tracks,      (for a size,
  validate, defaults)        sources bound; annotations  layers, themes,         a viewport and
                             and their access plans)     placement through maps) the data store)
                                       ▲                                              │
                                       └── data store ◀── runtimes ◀── requests ◀─────┘
```

- **Resolved model**: every coordinate system and its numberings; units and copy groups; maps resolved to their
  encodings with paths precomposed ([maps.md](maps.md)); each source with the runtime bound to it (unavailable with an
  issue when none); each annotation's shape, fields and access plan (inline columns, or regions and levels from a
  source). Pure, cached by the spec's hashes, keyed by address.
- **Resolved view** (public, renderer-independent, [rendering.md](rendering.md)): sections with their axes and layout
  maps, tracks, layers with themes as values or value functions, each layer's fallback chain, and how its data reaches
  its section. No data, no pixels, no fonts, no target: it is the same for every renderer. Pure, cached, keyed.
- **Draw lists**: the resolved view laid out for a size, a viewport, a target (which picks each layer's first drawable
  fallback) and the data store's current entries. Recomputed for what changed; built into reused buffers, so pure by
  contract (same inputs, same lists), not by allocation ([performance.md](performance.md)).
- **The data store** ([data-store.md](data-store.md)) sits outside the passes: requests flow from the resolved view and
  the viewport to the runtimes; answers land in the store with versions; the draw lists' cache key includes the
  versions of the entries they read. Nothing from the store goes into the spec, the resolved model or the resolved view.
- **Keys and versions**: every entry of a compiled form carries its address and a version from its inputs, so the
  viewer, an export or another renderer applies differences, not rebuilds.

## Rejected alternatives

- A language and a version per level, and five IRs (data, model, annotation, view, render), as first drafted: more
  descriptions and boundaries than the work needs; data, model and annotation resolve together.

## Consequences

- Each consumer reads one form: the viewer and exports the draw lists, other viewers (Nightingale, JBrowse) the
  resolved view, runtimes their source entries in the resolved model, tests any of them.
- The core's folders follow the pipeline: `lang/` (the description, the builder), `resolve/`, `view/`, `data/` (the
  protocol, the store, the runtime contract), `render/` (layout, draw lists, the SVG and recording backends).
