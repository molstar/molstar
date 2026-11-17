# Features

## Intent

The scope of this version: what it promises, one row per feature. The designs say how things work, the edge-case
catalog ([edge-cases.md](edge-cases.md)) which behaviors must survive, `packages/seq/TODO.md` what is next. This says
what the version contains, so scope never grows or shrinks silently (the author, 2026-11-17: "ok" to the proposal).

## Rules

- **Status**: `planned` (in this version, not built), `building`, `done`, `deferred` (not in this version, with its
  trigger). A feature is `done` when its example view runs ([examples.md](examples.md)) and every edge case it covers
  has an outcome.
- **Changes are decisions**: moving a row to or from `deferred` is the author's call, with the reason in its row.
- **Sources**: the previous implementation's features page, the designs here, the JBrowse reply. A feature the
  previous implementation had that is not listed is an omission to fix, not a decision.
- **Provisional**: the first pass below is a proposal. `planned` rows are the recommended scope of this version; the
  `deferred` rows are the ones to review first.
- **Before `.seq/` is deleted**, this list becomes the Storybook's Features page, each row linking its view; the Cases
  column goes with the catalog it cites.
- **Version 1 and later milestones**: what version 1 contains and the order of work is [version-1.md](version-1.md);
  rows deferred to the Genomes milestone come back together.

Columns: the feature; what it does; (in Streaming and patches) what undo does to it; its status; the design; the edge cases it covers (EC numbers, catalog); the view
that shows it (dataset: view).

## Language and spec

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| One tree, one node shape | `{ kind, id, params, children }`, addresses, refs | planned | [format](format.md) | 68–71, 77 | every view |
| Language as data | one description: validation, defaults, builder types, docs, agent skill | planned | [format](format.md) | 86–88 | — |
| Typed builder | the main way to write a spec; misuse is a type error | planned | [format](format.md) | — | every view |
| S-expression printing and reading | inspection, docs, tests, agent input; never stored | planned | [format](format.md) | — | every view (printed beside it) |
| Validation, every issue | stable codes, hints, counts, severity with one meaning | planned | [format](format.md) | 6, 76–78, 86–89 | Antibody: invalid labels |
| Versions and migrations | one semver, explicit migrations, lenient minors | planned | [format](format.md) | 72–74, 235 | — |
| Encodings | JSON by default; MessagePack for large inline data (typed columns as extension types); columns declared, hashed by bytes | planned | [format](format.md) | 53, 90 | every view (both encodings in its test) |
| Archives, through a port | a spec with its data files and an index, read and written by the archive port; a zip adapter out of the box | planned | [ports](ports.md) | 91 | every dataset, saved and opened as an archive |
| Published JSON Schema | derived from the description | deferred | [format](format.md) | — | trigger: an outside reader or agent tool needs it |

## Coordinates and maps

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Indices and labels | 0-based indices, labels converted at the edge | planned | [building-blocks](building-blocks.md) | 1–8 | Antibody: Kabat vs IMGT |
| Numbering schemes | offset, segments, labels, names; Kabat, IMGT | planned | [maps](maps.md) | 2–8 | Antibody: Kabat vs IMGT |
| Cyclic systems | origin, wrapped ranges, shortest arc | planned | [building-blocks](building-blocks.md) | 9–14 | Plasmid: across the origin |
| Map kinds | identity, offset, numbering, segments, pairwise, alignment, codon, structure | planned | [maps](maps.md) | 15–30 | p53: protein on its gene |
| Map sourcing | inline, from a source by region, computed by a runtime | planned | [maps](maps.md) | 212 | p53: structure map computed by alignment |
| Map encodings | implicit, runs, expanded by region; composed without expanding | planned | [maps](maps.md) | 27, 33, 34 | Globins: 10,000 rows |
| Direction on every run | inversions, reverse strands, reversed pairs rows | planned | [maps](maps.md) | 206 | p53: the minus-strand gene; Plasmid: features on both strands |
| Paths and projection | bounded shortest paths, the copy rule, ambiguity reported | planned | [maps](maps.md) | 19–28, 207 | p53: structure chains |
| Closure | a nucleotide lifts to its residue, one step | planned | [building-blocks](building-blocks.md) | 30, 144 | p53: protein on its gene |
| Structure systems and maps | chain instances, SIFTS in abstract form | planned | [maps](maps.md) | 65, 118 | p53: structure chains |
| Codon tables | NCBI tables, Sec, Pyl, stops | planned | [maps](maps.md) | 16–18 | p53: protein on its gene |
| Pieces and layout maps | ranges, reversed pieces, collapsed introns, letter maps | planned | [building-blocks](building-blocks.md) | 31–33, 233 | p53: introns collapsed |

## Data and sources

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Sequences, units, complexes, copies | molecules shared by units, copy groups | planned | [building-blocks](building-blocks.md) | 35, 36, 65–67 | Globins: hemoglobin |
| Ligands | one position or one per atom; chips in text | planned | [building-blocks](building-blocks.md) | 64, 99 | Globins: hemoglobin with heme |
| Alignments | gaps, insert columns, numbered rows, lazy rows | planned | [building-blocks](building-blocks.md) | 38–44 | Globins: family |
| Trees | given (Newick) or computed; clades collapse | planned | — | 43, 100, 103 | Globins: family by tree |
| Annotations, six shapes | ranges, values, pairs, matrix, grid, rows | planned | [building-blocks](building-blocks.md) | 45–63 | every dataset |
| Keyed tables and joins | rows keyed by id joined to items, rows, states | deferred | — | 106, 224 | trigger: the Genomes milestone ([version-1](version-1.md)) |
| Gene hierarchies | gene, transcript, exon (several levels) | deferred | — | 227 | trigger: the Genomes milestone ([version-1](version-1.md)) |
| Formats in the core | FASTA, aligned FASTA, A3M, Stockholm, Clustal, Newick, TSV, CSV, JSON | planned | [sources](sources.md) | 38–40, 49, 50, 92 | every dataset |
| More genome formats | 2bit, BigBed, VCF with tabix, PAF | deferred | [sources](sources.md) | 208 | trigger: the bigWig runtime done; each with the dataset that needs it |
| Services at the application layer | UniProt, PDBe (SIFTS), AlphaFold DB as runtimes an application composes, not in the core; recorded fixtures offline, live optional | planned | [sources](sources.md) | 51 | Storybook: an application composing services |
| Structure files at the application layer | mmCIF and BinaryCIF read (through Mol\*'s readers) for chains, `label` and `auth` numbering, insertion codes, observed residues | planned | [maps](maps.md) | 3–5, 118 | Storybook: predictions onto a structure |
| Name index | open on a gene by name from a file read by region | deferred | — | 228 | trigger: the Genomes milestone ([version-1](version-1.md)) |

## Streaming and patches

How state changes ([changes.md](changes.md)): patches edit the spec and are the only undoable changes; streams bring
data into the store and are never undone or saved. Both are increments by address with versions. The Undo column says
what undo does to each.

### Patches

| Feature | What | Undo | Status | Design | Cases | View |
|---|---|---|---|---|---|---|
| Patches, inverses, diff | five ops by address, all or nothing, `diff` correct | yes, one step a batch | planned | [format](format.md) | 78–83, 236 | every editing view |
| Undo and redo | every spec edit; previews end in a patch or a cancel | — | planned | [ux-foundations](ux-foundations.md) | 163, 164 | every editing view |
| Previews | a gesture shows its change as a transient overlay over the spec, ending in one patch or a cancel | its final patch | planned | [changes](changes.md) | 163 | every editing view |
| Snapshots | named patches over the base, a player | no: leaving a step is navigation | planned | [format](format.md) | 84, 85 | p53: a story |
| Patches on streamed things | a bookmark or a ref naming a dynamic node; waits while the node is not held | yes | planned | [changes](changes.md), [addressing](addressing.md) | — | p53: transcripts by region |
| Keep this region | a command that turns streamed data into spec content on purpose (an inline copy, a ref) | yes | planned | [changes](changes.md) | 175 | p53: the gene's region by region |
| Patch streams | patches from an agent, a collaborator or a host applied in order, all or nothing each; a stale base refused or rebased by the host | yes, each a step | planned | [changes](changes.md) | 81, 93 | an editing story with an agent's patches |
| Session state | the viewport, the selection, the hover, the cursor, open panels: transient, held by the viewer's models | no | planned | [changes](changes.md) | — | every view |
| Change kinds of commands | every command declares whether it patches the spec, changes session state or navigates; menus show what undo will undo | — | planned | [changes](changes.md) | — | the command palette |

### Streaming: the source protocol and the data store

| Feature | What | Undo | Status | Design | Cases | View |
|---|---|---|---|---|---|---|
| Self-describing sources | format, location, index, assembly, aliases, fields | no: data | planned | [sources](sources.md) | 185–191, 196–198 | p53: the gene's region from FASTA with its index |
| The source protocol | a source answers what it covers (systems, lengths, levels, nodes known up front) and a region at a level of detail; one protocol for inline and streamed data | no: data | planned | [sources](sources.md) | 169–175, 185–199 | p53: the gene's region by region |
| Answers carry nodes | tables plus data nodes (sequences, maps, annotations, nested sources), identity by (source, key) | no: data | planned | [dynamic-nodes](dynamic-nodes.md) | 175–180 | p53: transcripts by region |
| Levels of detail | the coarsest level that fits a device column; a format's own zoom levels used when it has them (bigWig) | no: data | planned | [performance](performance.md), [sources](sources.md) | 170, 171 | p53: conservation signal over the gene's region |
| The data store | entries by (source, region, level), versioned; pending, loaded, empty, failed; reactive per entry | no: data | planned | [data-store](data-store.md) | 169 | p53: the gene's region by region |
| Prefetch and cancellation | view, margin, prefetch; grace; retries | no: data | planned | [data-store](data-store.md) | 172–174 | p53: the gene's region by region |
| Eviction | byte-bounded, least recently used, never what is in view | no: data | planned | [data-store](data-store.md) | — | p53: the gene's region by region |
| Settle and open domains | open value domains computed when the view settles; legends state their domain | no: data | planned | [data-store](data-store.md) | 201–205 | p53: tracks |
| Export waits for data | `renderSvg` requests what its viewport needs and waits to settle, or says it is unsettled | no: data | planned | [data-store](data-store.md), [rendering](rendering.md) | 125 | every view (SVG) |
| Dynamic nodes and templates | data nodes that arrive with regions (item coordinate systems, nested sources); a track per node by template | no: data | planned | [dynamic-nodes](dynamic-nodes.md) | 175–180 | p53: transcripts by region |
| Item coordinate systems | an item's own system from its `map`, made on demand, keeping its id when it returns | no: data | planned | [dynamic-nodes](dynamic-nodes.md) | 176, 177 | p53: transcripts by region |
| Nested sources | requested on need, bounded depth | no: data | planned | [dynamic-nodes](dynamic-nodes.md) | 179 | p53: transcripts by region |
| Format runtimes | a registry entry per format, converting it onto the protocol; one shared conformance contract | no: data | planned | [sources](sources.md), [ports](ports.md) | — | every source |
| A genome format, as proof | bigWig read by region and zoom level, a format runtime like any other (the author: "at least one as proof of concept") | no: data | planned | [sources](sources.md) | 169–175, 186, 201 | p53: conservation signal over the gene's region |
| A synthetic asynchronous source | seeded, with a fake clock: the protocol's test source and the first slice's | no: data | planned | [version-1](version-1.md) | — | tests; a Storybook story |
| Range reads of assets | the asset manager's page cache for indexed formats (a base-library change) | no: data | planned | [base-library](base-library.md) | — | p53: conservation signal |

## Layers and themes

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Letters | cells and letters, density steps, letter caps | planned | [ux-foundations](ux-foundations.md) | 116, 117, 120–122 | p53: tracks |
| Ranges | lanes, labels, collapse, split, line styles, arrowheads | planned | — | 45, 104, 105, 128, 226, 232 | p53: tracks |
| Values | bars, line, area; shared domains; open domains at settle | planned | [data-store](data-store.md) | 52–57, 201, 202 | p53: tracks |
| Values for genomics | points, split by sample, log scale, reference line | deferred | — | 214–217 | trigger: the Genomes milestone ([version-1](version-1.md))
| Variants | sites colored by a field | planned | — | — | p53: variants |
| Variants by sample | samples as rows, long deletions solid | deferred | — | 222 | trigger: the Genomes milestone ([version-1](version-1.md))
| Pairs | arcs above and below, kinds, chords on rings | planned | — | 62, 63, 129 | tRNA: arcs |
| Matrix | positions × states, majority and maximum density | deferred | — | 60, 61, 220 | trigger: the Genomes milestone ([version-1](version-1.md))
| Grid and dot plots | a value per pair of positions; crosshair hover | planned | — | 58, 59, 109, 145 | tRNA: base-pair dot plot |
| Logo | stacked letters by bits or probability | planned | — | — | Globins: family |
| Secondary structure | helices and strands, cut ends square | planned | — | 67, 130 | p53: tracks |
| Ruler and numbering | tabular, skipped labels, absent or placed | planned | [ux-foundations](ux-foundations.md) | 234 | every view |
| Frames | outlined regions across tracks, open where cut | planned | — | 100 | Globins: ESPript figure |
| Rows layer | rows in tree or data order, filters, data columns | planned | — | 106, 107 | Globins: family |
| Fallbacks | every layer type declares one; readers degrade predictably | planned | [extensibility](extensibility.md) | 74, 238–244 | Nightingale reference |
| Themes | uniform, field palettes (binned, continuous, log), rules, schemes | planned | [ux-foundations](ux-foundations.md) | 133–137, 203–205, 215 | p53: themes |
| Residue schemes | Clustal X, Zappo, Taylor, hydrophobicity, nucleotides, ESPript | planned | [ux-foundations](ux-foundations.md) | 135 | Globins: family |
| Shared domains | domain groups across layers | planned | — | 56 | p53: tracks |
| Legends | ramps and steps, merged, painted vs vocabulary | planned | [ux-foundations](ux-foundations.md) | 132, 204 | p53: themes |
| Size channels | thickness, glyph, radius, width, opacity | planned | — | 137, 242 | p53: tracks |
| Layers across sections | ribbons between sections, a dot plot with its own y section | deferred | — | 209 | trigger: the Genomes milestone ([version-1](version-1.md))
| Contact maps | binned grids by region and level (Hi-C) | deferred | — | 223 | trigger: layers across sections |
| MAF over a genome | a multiple alignment on a genome axis | deferred | — | 213 | trigger: a dataset that needs it |

## Views and layouts

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Tracks layout | sections, tracks, groups, summaries | planned | [building-blocks](building-blocks.md) | 101–105 | p53: tracks |
| Text layout | wrapped groups of ten, lanes, inline labels | planned | — | 94–99 | p53: text figure |
| Alignment figures | rows in text, ESPript style, numbers per row | planned | — | 97 | Globins: ESPript figure |
| Blocks | side by side, stacked, many blocks | planned | — | 102, 115 | p53: structure chains |
| Two assemblies on one axis | blocks of two genomes on one axis or ring | deferred | — | 210 | trigger: the Genomes milestone ([version-1](version-1.md))
| Rings | tracks as rings, labels with leaders, a window | planned | — | 108, 114 | Plasmid: ring with linear detail |
| Navigators | strip, miniature, ideogram, vertical overview | planned | — | 110 | p53: the gene's region; Globins: 10,000 rows |
| Follow and scroll groups | sections follow through maps (scale 1, 3, 1/3) | planned | — | 113, 114, 211 | p53: protein on its gene |
| Viewport persistence | keep center and zoom across updates; pending intents | planned | — | 111, 112 | every view |
| Copies merged or split | `×N` blocks, hover for every member | planned | — | 65 | Globins: hemoglobin |

## Interaction and editing

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Configurable grammar | gestures bound to commands; presets; one per viewer | planned | [ux-foundations](ux-foundations.md) | 148–155 | every view; a preset story |
| Hover across linked sections | one owner, projected, emphasis outline | planned | [ux-foundations](ux-foundations.md) | 138, 141–147 | p53: orthologs |
| Selection | click, extend, toggle, snapped drag, projected | planned | [ux-foundations](ux-foundations.md) | 148–151 | every view |
| Keyboard cursor | crosses gaps, pieces, lazy rows | planned | [ux-foundations](ux-foundations.md) | 154, 155 | Globins: 10,000 rows |
| Focus and reveal | zoom to, or pan only when hidden | planned | [ux-foundations](ux-foundations.md) | 12, 158, 159, 168 | p53: search |
| Highlight channels | named, cleared independently | planned | [ux-foundations](ux-foundations.md) | 157 | p53: search |
| Search and results | matchers (exact, degenerate, regex, approximate); stepping; saved as sets | planned | — | 161 | p53: search |
| Query language | composable queries, printed and read as S-expressions | planned | [building-blocks](building-blocks.md) | 89, 229 | p53: queries |
| User-defined refs | named addresses in a `refs` container; bookmarks and location sets are refs; rewritten on move | planned | [addressing](addressing.md) | 162 | Antibody: CDR bookmarks |
| Track and group editing | hide, reorder, collapse, theme changes | planned | [ux-foundations](ux-foundations.md) | — | every dataset's editing view |
| Hover card and context menu | actions on the element; keyboard menu | planned | [ux-foundations](ux-foundations.md) | 156, 166 | p53: tracks |

## Rendering and output

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Resolved view and draw lists | public forms other renderers read | planned | [rendering](rendering.md) | 238–244 | Nightingale reference |
| Canvas2D backend | glyph atlas, hiDPI, band copy | planned | [rendering](rendering.md) | 120–123, 140 | every view |
| SVG and PNG | headless, same picture as the canvas | planned | [rendering](rendering.md) | 124, 125 | every view (side by side) |
| Render CLI | a spec file to SVG or PNG | planned | [rendering](rendering.md) | 126 | — |
| Text metrics port | the same layout in the browser and Node | planned | [rendering](rendering.md) | 124 | — |
| UX foundations | tokens, squircles, emphasis, motion, density steps | planned | [ux-foundations](ux-foundations.md) | 138, 139 | foundation stories |
| Writers | a selection as FASTA; a figure as PDF | deferred | [extensibility](extensibility.md) | — | trigger: the first extension of writers |

## Configuration and extensions

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Ports with default adapters | assets (`AssetManager`), text metrics, structure; formats, maps, backends as registry entries | planned | [ports](ports.md) | — | Extensibility page |
| Configurable layout and spacing | every foundation number a setting, in presets; layouts behind one internal interface | planned | [ux-foundations](ux-foundations.md) | — | foundation stories, a preset switcher |
| Typed configuration | one set of config items at composition (`ParamDefinition`) | planned | [extensibility](extensibility.md) | — | Extensibility page |
| Settings in layers | user settings over the configuration; effective values with origin | deferred | [extensibility](extensibility.md) | — | trigger: the typed configuration has users |
| Generated settings form | from the config items' `ParamDefinition` | deferred | [extensibility](extensibility.md) | — | trigger: settings in layers |
| Extension descriptors | passed at composition; conflicts are errors | planned | [extensibility](extensibility.md) | — | Extensibility page |
| Node kinds and custom fields | namespaced, carried through untouched | planned | [extensibility](extensibility.md) | 235 | Extensibility page |
| Registry entries | node kinds, layer types, themes, formats, map runtimes, backends, query functions, commands; layouts and chrome parts internal until an outside user | planned | [extensibility](extensibility.md) | — | Extensibility page |

## Integration

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Nightingale reference | a subset drawn by Nightingale from the resolved view | planned | [rendering](rendering.md) | 183 | Nightingale reference |
| Mol\* structure link, as a showcase | translations between a sequence's own system (often position-indexed: a lab's predictions per letter) and mmCIF `label` and `auth` numbering, insertion codes and unobserved residues; hover and selection both ways; per-position values colored in 3D. The structure runtime is application-layer code (the storybook package first; a package of its own when another host needs it) | planned | [maps](maps.md) | 3–5, 25, 118, 167, 181 | Storybook: predictions onto a structure |
| Addressing | local, hierarchical ids; addresses as `path/(query :key value)` (meso\*'s grammar) for nodes, locations in labels, elements and rects; commands take addresses; copy link | planned | [addressing](addressing.md) | 1–8, 182 | Storybook: addressing |
| Markdown links | links in Mol\*'s grammar driving the view; one href for Mol\* and Seq\* | planned | [addressing](addressing.md) | 182 | Storybook: addressing |
| Agent tools | read and write a spec through addresses and expressions; the skill | deferred (the author, 2026-11-17) | [addressing](addressing.md) | — | trigger: the language description and addressing done |
| React adapter | the default UI in React | deferred | [packages](packages.md) | — | trigger: a React host |
