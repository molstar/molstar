# Composition: languages, levels and the IRs

## Intent

The author, 2026-11-17: "expand on how the "spec" will be composed? the views etc .. similar to how meso composes many
concepts. and do we need intermediate representation the spec compiles/expands to? i think this should streamline many
things, including mapping the "streaming ports" etc?"

A spec is written by people and agents, short and declarative. What the viewer, the source runtimes and other
renderers need is explicit: every ref followed, every default filled, every map resolved, every request known. Between
them sit intermediate representations (IRs), compiled by pure passes, so each consumer reads one compiled form and none
re-derives it.

## Precedents

- meso\*: languages by level (data, molecular, scene, time, view, session), one tree; slots embed one language in
  another; definitions in containers, used by reference with arguments and overrides; IRs (data, scene, time, view,
  render) compiled by pure passes cached by hash and keyed by address; streamed data declared in the spec, loaded by
  the runtime through ports, never written back; the render IR keyed and versioned so an executor commits differences
  ([sources/meso](../sources/meso/README.md): its packages, state and runtime specs).
- The previous implementation: one document grammar with `data` and `view` under a root; `compile` and `resolveView`
  as two passes; streamed sources resolved by host stream providers against a manifest ([sources/seqstar](../sources/seqstar/README.md)).

## Decision

### The spec is composed of levels, each a language

| Level | Language | Holds |
|---|---|---|
| `sources` | data | where data comes from: a source per file, service or inline value, in its own format ([sources.md](sources.md)) |
| `model` | model | coordinate systems: sequences, units, complexes and copies, structures' chains, alignments and their rows, trees; numberings; maps between systems ([maps.md](maps.md)) |
| `annotations` | annotation | data on systems in the six shapes, read from sources or written inline; computed nodes (a query's result, coverage, conservation) |
| `views` | view | views, sections (axes and pieces), tracks, layers, themes, focus, highlights; the layouts (tracks, text, rings, dot plots) |
| `refs` | refs | user-defined names for addresses: bookmarks, regions, sets ([addressing.md](addressing.md)) |
| `snapshots` | patch | named patches over the base: a story's steps |

- **One tree, one root**: `spec` holds the levels as its children, each with a reserved id, so `/views/main` and
  `/model/p53` are addresses. A level may be absent; a small spec writes `sources`, `model` and `views` only.
- **Each language has its own description and version**, as meso\*'s do; the spec's version is the set of them. A
  level's language is described once ([format.md](format.md)) and validated, documented and typed for the builder from
  that description.
- **Languages embed each other through slots**: a layer's `data` is a ref or a query (the query language), a filter
  is a query, a theme is a node of the theme sublanguage, a computed annotation holds a query or a provider's params.
  Sublanguages: queries, themes, patches, addresses.
- **Definitions and uses** (meso\*'s `def` and `use`): anything definable is written in place or defined once and used
  by address, with arguments and overrides. Definition containers: `views/themes` (named themes and palettes),
  `views/tracks` (track templates: a track with layers, used per sequence or per sample), `views/styles` (section
  styles and presets), `model/numberings` (schemes a lab defines). A template used per member of a set (every row,
  every sample) is how many tracks come from one definition.
- **Extensions add kinds to a level** under their namespace, with a fallback ([extensibility.md](extensibility.md)).

### Session is not the spec

- What is shown now (the view chosen, the viewport, the selection, the hover, the cursor, open panels) is session
  state: transient, reconstructible, saved beside the spec when saved at all, never in it. Snapshots and refs are how
  a spec names a view state on purpose.

### The IRs

Pure passes, each cached by its input's hash and keyed by address, so a patch recompiles only what it touched:

```
spec ─read─▶ checked spec ─compile─▶ data IR ──▶ model IR ──▶ annotation IR ──▶ view IR ─layout─▶ render IR ─▶ backend
 (decode, migrate,          (sources bound   (systems,      (shapes, fields,     (the resolved   (draw lists
  validate, defaults)        to runtimes)     maps, paths)   access plans)        view)           for a viewport)
                                  ▲                                                    │
                                  └──── requests: what each visible layer needs ◀──────┘  (session: viewport)
```

- **Data IR**: one entry per source: its address, its format and the runtime bound to it (by the application's
  composition), its location and index, what it declares it covers. This is where the streaming ports map: a runtime
  binds to a data IR entry, never to a spec node, and answers the one protocol for it. A source with no runtime is an
  entry marked unavailable with its issue.
- **Model IR**: the coordinate systems table (every system with its length, alphabet, cyclic flag, numberings), units
  and copy groups, the maps resolved to their encodings ([maps.md](maps.md)) with their sourcing (inline, a data IR
  entry by region, a computed job), and the map graph with its paths precomposed (never on hover). Nodes made from
  data (a structure's chains, an alignment's rows) are here, keyed by address.
- **Annotation IR**: per annotation, its shape, its system, its fields and their types, and its access plan: inline
  tables (typed columns), or regions and levels from a data IR entry, with the multi-resolution summaries it needs
  ([performance.md](performance.md)).
- **View IR** (the public resolved view, [rendering.md](rendering.md)): sections with their axes and layout maps,
  tracks, layers with themes resolved to values or value functions, fallbacks chosen for the target, and for each
  layer how its data reaches its section (direct, or through a path). Renderer-independent.
- **Render IR**: the view IR laid out for a size and a viewport, as draw lists; recomputed per frame for what changed,
  bounded by the pixels on screen ([performance.md](performance.md)).
- **Requests flow back**: from the view IR and the session's viewport, each visible layer says which regions at which
  level it needs; the data IR's runtimes answer through the protocol; answers land in a data store keyed by (source,
  region, level), never in the spec or the IRs, and the render IR draws from it.
- **Keys and versions**: every IR entry carries its address as key and a version from its inputs' hashes, so the
  viewer, an export or another renderer applies differences, not rebuilds.

## Consequences

- Each consumer reads one compiled form: the viewer the view and render IRs, a source runtime its data IR entry, an
  export the render IR, Nightingale and other viewers the view IR, tests any of them.
- The storybook's Format, Sources and Rendering pages describe the levels and the IRs; a Composition page shows a
  spec and each IR it compiles to, side by side.
- The core's folders follow the pipeline: `lang/` (one folder per level and sublanguage), `ir/`, `compile/`, `data/`
  (the protocol, the store, the runtimes' port), `render/`; the packages design lists them as they are built.
