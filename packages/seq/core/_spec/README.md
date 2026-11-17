# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description, validation, the
builder, printing, patches, the coordinate model, the source protocol, queries, and the pure passes from a spec to draw
lists. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

## 1. The spec is one tree

- A spec is a tree with two lossless encodings: MessagePack (the default for files and transfer; typed columns of
  numbers as extension types `0x11`–`0x18`) and JSON (readable; typed columns as arrays, NaN as `null`). Both decode to
  an equal tree; a typed column means the array of its numbers, a float32 its shortest decimals. Nodes are
  `{ kind, id?, params?, children? }`. Why: JSON every tool reads, MessagePack for real data, one meaning.
- `id` is local (unique among siblings); an address is the path of ids from the root, hierarchical as an IPv6 address is
  a path of groups (`/data/p53/features`). Why: short, composable specs, and an address names rather than counts, so
  inserts and reorders keep it.
- The S-expression spelling is a printing for people, docs, tests and agents, with a reader that lowers to the same
  tree; it is never a file type. Why: compact and diffable to inspect, without a second canonical form.

## 2. The language is data

- One description declares the kinds, where each may sit, the id policy, params (type, default, bounds, docs) and refs.
  Validation, defaults, the builder's types, the docs and the agent skill are derived from it. Why: nothing is written
  twice, so nothing drifts.
- Validation reports every issue, `{ code, severity, path, message, hint? }`, with a stable namespaced code. Why:
  authors (people and agents) fix a spec in one pass, and callers branch on codes.
- An unknown kind from a newer minor version or a missing extension is kept, reported and skipped (drawn as its fallback
  when it has one). Why: a spec opens wherever most of it can be read.

## 3. Versions

- The language has one semver (`LanguageVersion`); a spec states its version. Older specs are migrated by explicit
  tree-to-tree functions; a newer major version is refused with an issue.
- No compatibility with earlier SeqViewSpec documents before 1.0 of this language.

## 4. Patches

- A spec changes only by patches: `set`, `unset`, `insert`, `remove`, `move`, each naming a node by address. A batch is
  all or nothing and returns its inverse; `diff(a, b)` is the patch from `a` to `b`. Why: one way to change a spec for
  the UI, commands and agents, with undo for free.

## 5. Positions and maps

- Positions are indices, 0-based; ranges are half-open. Labels (residue numbers, insertion codes, chromosome names) are
  converted at the edge, once. Why: one arithmetic everywhere.
- A map has a kind (what it means: `identity`, `offset`, `numbering`, `segments`, `pairwise`, `alignment`, `codon`,
  `structure`), a sourcing (inline in its kind's compact form, a source read by region, or computed by a runtime the
  application binds) and an encoding (implicit, runs, or expanded per position for the region asked). Computed is how a
  map is sourced, never its kind. Why: the spec stays small, and a map is held no larger than its meaning needs.
- Every resolved map answers one interface: the image of a range, forward and back, as runs with a direction (+1 or −1)
  and a step (1, or 3 for codons); the boundary rules; each element's status. Paths compose encodings without expanding
  them. Why: an inversion, a reverse-strand gene, a protein on its codons and a chain in a structure are one mechanism.
- A structure system is a chain instance of a structure, its residues keyed by the structure's numbering and flagged
  observed; a `structure` map relates a sequence to it, per segment with its provenance and per residue with its status
  (SIFTS in abstract form). The core reads no structure file; a structure runtime binds the system to atoms.

## 6. Sources

- A source names its data in the data's own terms (format, location, index, what the format does not say) and is read
  through one protocol: what it covers, and a region at a level of detail, answered as the model's tables. A runtime per
  format, provided by the application, converts the format onto the protocol. A source whose format has no runtime is
  unavailable, with an issue. Why: a spec opens without host code, and the model reads one thing.
- Nothing depends on a whole-data manifest: an open value domain means the data in view; an unlisted category takes a
  color by a stable hash of its value.

## 7. Levels and IRs

- The spec's levels are languages, each with its description and version: `sources`, `model`, `annotations`, `views`,
  `refs`, `snapshots`. Slots embed sublanguages (queries, themes, patches, addresses); anything definable is written in
  place or defined once and used by address with arguments and overrides. Session state (viewport, selection, hover) is
  not the spec. Why: each concern is small to learn and versions on its own.
- The spec compiles through pure passes to IRs: data (sources bound to runtimes), model (systems, resolved maps,
  precomposed paths), annotation (shapes, fields, access plans), view (the resolved view) and render (draw lists for a
  viewport). Each entry is keyed by address and versioned by its inputs' hashes; a patch recompiles what it touched.
  Requests for regions flow back from the view IR and the viewport to the data IR's runtimes; answers go to a data
  store, never into the spec or the IRs. Why: every consumer reads one compiled form, and runtimes bind to IR entries,
  not to spec nodes.
- A frame costs what is on screen, not the data's size: every layer reads a multi-resolution structure (summary
  pyramids, interval indexes, virtualized rows, 2D pyramids) by region and resolution, and draw lists are bounded by
  device pixels. Why: the same frame time for a protein and a genome.

## 8. Rendering forms

- Three forms are public contracts, each a place another renderer can attach: the spec (documented by the generated
  language reference); the resolved view (refs followed, defaults filled, maps resolved, themes resolved to values,
  fallbacks chosen for what the target can draw; no pixels, no fonts); and draw lists (the resolved view laid out for a
  size and viewport, as typed primitives in passes). Why: other viewers draw Seq\*'s meaning without reimplementing it,
  and every output draws the same picture.
- A backend draws draw lists. The core has the SVG backend (a string, no DOM) and a recording backend, so a spec renders
  headless in Node (`renderSvg`, and a CLI to SVG or PNG); text metrics are a port (a canvas in the browser, a bundled
  metrics table in Node).

## 9. Configuration and extensions

- What can be configurable, is. Settings are typed data with defaults, bounds and docs, in layers: built-in defaults,
  the viewer configuration, user settings; the spec holds what is shown. The effective value of a setting can be asked
  for with its origin. Why: applications differ, and a fixed choice is a fork waiting to happen.
- `createSeqContext(spec)` composes what rendering needs without a viewer (registry, config, ports: assets, text
  metrics, compute, backends), for the CLI, servers and tests.
- Every registry is an extension point (node kinds, custom fields, formats to read and write, layer types, themes and
  palettes, map runtimes, numbering schemes and codon tables, query functions, backends and text metrics). Custom fields
  are namespaced (`custom: { "<ns>": … }`, and extra table columns) and carried unchanged through parse, patch, diff and
  print by readers that do not own them. An extension is a descriptor passed at composition; nothing registers on
  import; built-ins use the same points; conflicts are errors at composition, never resolved by order. Extension kinds
  in a spec are namespaced and declare a built-in fallback.

## 10. Ports

- The core does no I/O and reaches no platform: fetching, archives, format and map runtimes, text metrics, backends,
  writers, time and frames, compute, storage, the clipboard, links and services are ports, each with a conformance
  suite. A viewer composed with nothing binds the default adapters (archives as zip among them); an application replaces
  any. A port with no adapter makes what needs it unavailable, with an issue. Why: the same core runs in a browser, in
  Node, in a worker and in tests, and what a viewer can reach is visible at its composition root.

## 11. Addresses and refs

- Addresses are `path/(query :key value)`: ids as segments, query forms (the query language's) as segments selecting
  among what the prefix names (child nodes by kind and params; positions, items, rows below a node that holds them),
  `/*` for a subtree; one word wherever it appears, printed canonically. Why: one generic grammar for nodes, locations
  and elements, shared with meso\*.
- Rects (two ranges over two axes: positions × positions, rows × columns, states × positions) are addressable,
  selectable and linkable like ranges.
- Authors name addresses as refs in the `refs` container (bookmarks and location sets are refs); a `ref` param holds an
  address word, all checked in one resolve pass; a move or rename rewrites what names it.
- Commands take addresses; links carry an action and an address in Mol\*'s link grammar; the view gives back the address
  of what is under the pointer or selected. An address that no longer resolves is reported with what it named. Why:
  text, links, agents, tests and other viewers point at the same things the same way.

## Folders

None yet beyond the entry; each folder added becomes an entry (`@molstar/seqx-core/<folder>`) and joins the folder graph
test.
