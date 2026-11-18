# Designs

One file per topic. A design says what is true and why: intent, precedents, tensions, the decision, rejected
alternatives, consequences (§design-first). No status field: what is here is agreed. A section marked **Provisional**
is the agreed direction with its details still open; it says what would settle them.

| Design | Covers |
|---|---|
| [version-1.md](version-1.md) | **start here**: version 1's architecture as reconciled after the review, what it fixes, what it defers, and the order of work |
| [language-description.md](language-description.md) | the tree core's delta from meso\*, the description's shape (column and table params), the slice's kinds, the builder, issue codes |
| [locations-and-queries.md](locations-and-queries.md) | Nodes, Elements and Locations as result types; address forms by kind; numberings in forms; the first forms; staged evaluation |
| [protocol.md](protocol.md) | runtimes, coverage, requests and answers with table shapes, errors, the `seq` format, the store's keys, planning across maps, node lifetime |
| [layout-and-draw-lists.md](layout-and-draw-lists.md) | the layer-type contract, the layout interface, draw lists (units, passes, primitives, buffers), hit-testing, stable lanes, colors and fonts without CSS |
| [viewer-models.md](viewer-models.md) | the models and session state, viewports, scroll groups and follow, hover and selection, behaviors, gesture bindings, command change kinds, the slice's config items |
| [examples-and-stories.md](examples-and-stories.md) | the examples package, fixtures in Node and Storybook, one story per view |
| [security.md](security.md) | the asset policy, defaults by host, bounds on everything a spec sizes, escaping |
| [features.md](features.md) | the scope of this version: one row per feature with its status, design, edge cases and view (provisional: the first pass) |
| [vision.md](vision.md) | what Seq\* is, what it is for, and what is deliberately simpler this time |
| [building-blocks.md](building-blocks.md) | the root nouns (§name-the-roots) and how they fit: spec, data, coordinate systems, maps, views, layers, themes, queries, issues |
| [maps.md](maps.md) | maps: kind, sourcing (inline, a source, computed by a runtime) and encoding (implicit, runs, expanded); structure maps (SIFTS in abstract form) |
| [composition.md](composition.md) | one language with four levels (data, views, refs, snapshots); session apart; three compiled forms (resolved model, resolved view, draw lists) |
| [changes.md](changes.md) | streaming and patches as one incremental model; every kind of change with where it lives and whether it is undoable, saved, addressable |
| [patches.md](patches.md) | git's model for a spec: revisions as snapshots, commits, six patch operations with inverses, `diff` and `merge3`, the API of `repr/patch` |
| [merge.md](merge.md) | `merge3`: a three-way merge by node identity, conflicts as data, orders compared as orders, rebase from a base root |
| [immutable.md](immutable.md) | `@molstar/core/data/immutable`: CHAMP maps and sets, ordered maps and sets, a vector, transients; how `repr` uses them |
| [dynamic-nodes.md](dynamic-nodes.md) | data that arrives as nodes (sequences, maps, annotations, nested sources, item coordinate systems), the resolved model's overlay, track templates |
| [data-store.md](data-store.md) | where data read by region lives: keyed, versioned entries, pending vs empty, requests, eviction, settle and open domains |
| [base-library.md](base-library.md) | what Seq\* reuses from `@molstar/core`, `io` and `ui`, what to extend, promote and add there (provisional, each the author's call) |
| [tree-core.md](tree-core.md) | meso\*'s tree core moved into `@molstar/core/repr`, codec-independent; codecs at the edge |
| [format.md](format.md) | one tree, one language and version, typed columns, JSON by default and MessagePack for large data, S-expression printing, the language as data, the builder as the main way to write a spec, versions |
| [rendering.md](rendering.md) | the render port and backends (Canvas2D, SVG, PNG), headless rendering, and the forms other viewers read (the spec, the resolved view, draw lists) |
| [sources.md](sources.md) | where data comes from: one protocol, formats converted onto it by runtimes the application provides |
| [packages.md](packages.md) | the packages, their layers and names, and what each may import |
| [ux-foundations.md](ux-foundations.md) | the data's visual foundations on top of `@molstar/ui`'s (color, type, grid, density, when to squircle, emphasis, motion) and one interaction grammar for selection, toggling and navigation |
| [viewer-composition.md](viewer-composition.md) | the viewer composed from an explicit spec, as Mol\*'s plugins are: registry entries, behaviors, bindings, config, port factories with ownership; defaults as a spec; inside a Mol\* plugin |
| [ports.md](ports.md) | one per viewer is a port (assets, text metrics, structure), named providers are registry entries; defaults out of the box |
| [extensibility.md](extensibility.md) | what can be configurable is (the spec, the viewer configuration, user settings, layered); every registry open to extensions, passed at composition |
| [addressing.md](addressing.md) | everything in a view has an address (nodes, locations in labels, elements); commands and Mol\*-grammar links act on addresses |
| [performance.md](performance.md) | time, memory, size and scale targets, and how each is held by a test (provisional) |
| [testing.md](testing.md) | what holds each promise (properties, conformance suites, records, scaling ratios), in which lane, and how tests are written |
| [edge-cases.md](edge-cases.md) | the previous implementation's edge cases and rich behaviors: each kept, changed or dropped by decision, never by omission |
| [examples.md](examples.md) | six datasets, many views each (editing included); every view runs in Storybook, in a test and as a portable spec |
