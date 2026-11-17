# Designs

One file per topic. A design says what is true and why: intent, precedents, tensions, the decision, rejected
alternatives, consequences (§design-first). No status field: what is here is agreed. A section marked **Provisional**
is the agreed direction with its details still open; it says what would settle them.

| Design | Covers |
|---|---|
| [features.md](features.md) | the scope of this version: one row per feature with its status, design, edge cases and view (provisional: the first pass) |
| [vision.md](vision.md) | what Seq\* is, what it is for, and what is deliberately simpler this time |
| [building-blocks.md](building-blocks.md) | the root nouns (§name-the-roots) and how they fit: spec, data, coordinate systems, maps, views, layers, themes, queries, issues |
| [maps.md](maps.md) | maps: kind, sourcing (inline, a source, computed by a runtime) and encoding (implicit, runs, expanded); structure maps (SIFTS in abstract form) |
| [composition.md](composition.md) | the spec's levels (sources, model, annotations, views, refs, snapshots), each a language; definitions and uses; session apart; the IRs (data, model, annotation, view, render) and requests flowing back |
| [base-library.md](base-library.md) | what Seq\* reuses from `@molstar/core`, `io` and `ui`, what to extend, promote and add there (provisional, each the author's call) |
| [tree-core.md](tree-core.md) | meso\*'s tree core (tree, schema, patches, addresses, builder, descriptors) moved into `@molstar/core` |
| [format.md](format.md) | the spec's one tree: MessagePack (default) and JSON encodings, S-expression printing, the language as data, the builder as the main way to write a spec, versions |
| [rendering.md](rendering.md) | the render port and backends (Canvas2D, SVG, PNG), headless rendering, and the forms other viewers read (the spec, the resolved view, draw lists) |
| [sources.md](sources.md) | where data comes from: one protocol, formats converted onto it by runtimes the application provides |
| [packages.md](packages.md) | the packages, their layers and names, and what each may import |
| [ux-foundations.md](ux-foundations.md) | the data's visual foundations on top of `@molstar/ui`'s (color, type, grid, density, when to squircle, emphasis, motion) and one interaction grammar for selection, toggling and navigation |
| [viewer-composition.md](viewer-composition.md) | the viewer composed from an explicit spec, as Mol\*'s plugins are: registry entries, behaviors, bindings, config, port factories with ownership; defaults as a spec; inside a Mol\* plugin |
| [ports.md](ports.md) | every effect through a port bound at composition, with default adapters out of the box (fetch, archive, runtimes, metrics, backends, writers, host, compute, storage, clipboard, links, services) |
| [extensibility.md](extensibility.md) | what can be configurable is (the spec, the viewer configuration, user settings, layered); every registry open to extensions, passed at composition |
| [addressing.md](addressing.md) | everything in a view has an address (nodes, locations in labels, elements); commands and Mol\*-grammar links act on addresses |
| [performance.md](performance.md) | time, memory, size and scale targets, and how each is held by a test (provisional) |
| [testing.md](testing.md) | what holds each promise (properties, conformance suites, records, scaling ratios), in which lane, and how tests are written |
| [edge-cases.md](edge-cases.md) | the previous implementation's edge cases and rich behaviors: each kept, changed or dropped by decision, never by omission |
| [examples.md](examples.md) | six datasets, many views each (editing included); every view runs in Storybook, in a test and as a portable spec |
