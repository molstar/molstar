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
| [stories.md](stories.md) | a story as a sequence of full states stored as patches; transitions planned from diffs, sampled purely; the animation port by which the application owns the clock |
| [ux-foundations.md](ux-foundations.md) | the data's visual foundations on top of `@molstar/ui`'s (color, type, grid, density, when to squircle, emphasis, motion) and one interaction grammar for selection, toggling and navigation |
| [axis.md](axis.md) | a section's axis: its direction (`genomic`, `reversed`, `follow`, `{ along }`), another system's letters through a path (a residue a codon wide), collapsed introns as pieces with breaks, and what reads the axis (hits, selection, follow, focus, history, export) |
| [themes.md](themes.md) | the theme sublanguage (a slot per channel: uniform, field, rules, scheme, letters), palettes and `scale` nodes for shared domains, the residue schemes, tokens in light and dark, legends as data and a part a host places, the small layers |
| [layouts.md](layouts.md) | the section layout interface (strips, bands, warps), placement and hit-testing as its inverse, the one viewport record, text, blocks and rings, the ring's window and its follower |
| [grids-and-rows.md](grids-and-rows.md) | pairs as arcs, chords and dots; the grid section and dot plot (crosshair, rectangles); the rows layer (order, clades folded as session state, pins, lazy rows, the cursor); the logo; one reduction and one image cache |
| [navigators.md](navigators.md) | one part over an axis or two (the ruler rail, strip, ideogram, vertical overview, miniature), its window moving the view, idle and floored windows, a part a host places |
| [files.md](files.md) | files as sources: TSV, CSV, JSON, MessagePack (records and columns read the same) and GFF3 through one records runtime, the transform as params of the `source` node (`fields`, `records`, `base`, `ends`, `types`), one track a feature type by a leading item form, examples read from recorded files |
| [molecules.md](molecules.md) | small molecules, glycans, conjugates and modified residues as pointers and labels, never chemistry: the `molecule` kind, a sequence's modified positions, complexes and units, links, the atoms form; built in the fourth slice (0.5.0) |
| [viewer-composition.md](viewer-composition.md) | the viewer composed from an explicit spec, as Mol\*'s plugins are: registry entries, behaviors, bindings, config, port factories with ownership; defaults as a spec; inside a Mol\* plugin |
| [ports.md](ports.md) | one per viewer is a port (assets, text metrics; a structure viewer is not one), named providers are registry entries; defaults out of the box |
| [extensibility.md](extensibility.md) | what can be configurable is (the spec, the viewer configuration, user settings, layered); every registry open to extensions, passed at composition |
| [structure-syncing.md](structure-syncing.md) | a library, not a port: pure translation between Seq\*'s Locations on units and a structure's residues (`@molstar/seqx-core/structure`), `linkStructure` wiring a viewer's public events both ways, "color by this track" as runs, choosing between two maps of one pair (`priority`, `via`), the fake structure viewer |
| [track-header-parts.md](track-header-parts.md) | track and section headers: actions as `headerActions` registry entries in fixed slots, count chips, the menu behind its port, drawn on the canvas or a host's DOM parts over the published geometry, the secondary click a `headers` param |
| [seq-data-ports.md](seq-data-ports.md) | a host's data and Seq\*'s own reading on one path: data adapters by format (the core's runtimes the defaults, `withAdapter`), `StructureDescription` and `structureNodes` (`io/structure`), updates reopened in place by key, `source.ref`, `values.field`, annotations given back as the host's |
| [extension-openings.md](extension-openings.md) | the extension descriptor and `<id>:<name>` kinds with their fallbacks, an extension's params (plain on its types, namespaced keys elsewhere), the `extension` entry; custom navigator pictures and `NavigatorModel`; legends' library surface |
| [addressing.md](addressing.md) | everything in a view has an address (nodes, locations in labels, elements); commands and Mol\*-grammar links act on addresses |
| [performance.md](performance.md) | time, memory, size and scale targets, and how each is held by a test (provisional) |
| [testing.md](testing.md) | what holds each promise (properties, conformance suites, records, scaling ratios), in which lane, and how tests are written |
| [edge-cases.md](edge-cases.md) | the previous implementation's edge cases and rich behaviors: each kept, changed or dropped by decision, never by omission |
| [examples.md](examples.md) | six datasets, many views each (editing included); the exact data per dataset (entries, licences, sizes, what is generated), the fixture rules, what a user does in each view; every view runs in Storybook, in a test and as a portable spec |
