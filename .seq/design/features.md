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
- **Before `.seq/` is deleted**, this list becomes the Storybook's Features page, each row linking its view.

Columns: the feature; what it does; its status; the design; the edge cases it covers (EC numbers, catalog); the view
that shows it (dataset: view).

## Language and spec

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| JSON spec, one node shape | `{ kind, id, params, children }`, addresses, refs | planned | [format](format.md) | 68–71, 77 | every view |
| Language as data | one description: validation, defaults, builder types, docs, agent skill | planned | [format](format.md) | 86–88 | — |
| Typed builder | the main way to write a spec; misuse is a type error | planned | [format](format.md) | — | every view |
| S-expression printing and reading | inspection, docs, tests, agent input; never stored | planned | [format](format.md) | — | every view (printed beside it) |
| Validation, every issue | stable codes, hints, counts, severity with one meaning | planned | [format](format.md) | 6, 76–78, 86–89 | Antibody: invalid labels |
| Versions and migrations | one semver, explicit migrations, lenient minors | planned | [format](format.md) | 72–74, 235 | — |
| Patches, inverses, diff | five ops by address, all or nothing, `diff` correct | planned | [format](format.md) | 78–83, 236 | every editing view |
| Snapshots | named patches over the base, a player | planned | [format](format.md) | 84, 85 | p53: a story |
| MessagePack encoding | the default for files and transfer; typed columns as extension types; interchangeable with JSON | planned | [format](format.md) | 53, 90 | every view (both encodings in its test) |
| Archives, through a port | a spec with its data files and an index, read and written by the archive port; a zip adapter out of the box | planned | [ports](ports.md) | 91 | every dataset, saved and opened as an archive |
| Published JSON Schema | derived from the description | deferred | [format](format.md) | — | trigger: an outside reader or agent tool needs it |

## Coordinates and maps

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Indices and labels | 0-based indices, labels converted at the edge | planned | [building-blocks](building-blocks.md) | 1–8 | Antibody: Kabat vs IMGT |
| Numbering schemes | offset, segments, labels, names; Kabat, IMGT | planned | [maps](maps.md) | 2–8 | Antibody: Kabat vs IMGT |
| Cyclic systems | origin, wrapped ranges, shortest arc | planned | [building-blocks](building-blocks.md) | 9–14 | Plasmid: across the origin |
| Map kinds | identity, offset, numbering, segments, pairwise, alignment, codon, structure | planned | [maps](maps.md) | 15–30 | p53: protein on its gene |
| Map sourcing | inline, from a source by region, computed by a runtime | planned | [maps](maps.md) | 212 | Genomes: liftover by chain file |
| Map encodings | implicit, runs, expanded by region; composed without expanding | planned | [maps](maps.md) | 27, 33, 34 | Globins: 10,000 rows |
| Direction on every run | inversions, reverse strands, reversed pairs rows | planned | [maps](maps.md) | 206 | Genomes: synteny with an inversion |
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
| Keyed tables and joins | rows keyed by id joined to items, rows, states | planned | — | 106, 224 | Genomes: transcripts by joined results |
| Gene hierarchies | gene, transcript, exon (several levels) | planned | — | 227 | Genomes: genes with transcripts |
| One source protocol | coverage and regions, levels of detail, pending vs empty | planned | [sources](sources.md) | 169–175, 185–199 | Genomes: by region |
| Self-describing sources | format, location, index, assembly, aliases, fields | planned | [sources](sources.md) | 185–191, 196–198 | Genomes: FASTA with its index |
| Formats in the core | FASTA, aligned FASTA, A3M, Stockholm, Clustal, Newick, TSV, CSV, JSON | planned | [sources](sources.md) | 38–40, 49, 50, 92 | every dataset |
| A genome format, as proof | bigWig read by region and zoom level, a format runtime like any other (the author: "at least one as proof of concept"; bigWig recommended: it exercises the index, levels of detail and values) | planned | [sources](sources.md) | 169–175, 186, 201 | Genomes: signal by region |
| More genome formats | 2bit, BigBed, VCF with tabix, PAF | deferred | [sources](sources.md) | 208 | trigger: the bigWig runtime done; each with the dataset that needs it |
| Services at the application layer | UniProt, PDBe (SIFTS), AlphaFold DB as runtimes an application composes, not in the core; recorded fixtures offline, live optional | planned | [sources](sources.md) | 51 | Storybook: an application composing services |
| Structure files at the application layer | mmCIF and BinaryCIF read (through Mol\*'s readers) for chains, `label` and `auth` numbering, insertion codes, observed residues | planned | [maps](maps.md) | 3–5, 118 | Storybook: predictions onto a structure |
| Prefetch and cancellation | view, margin, prefetch; grace; retries | planned | [sources](sources.md) | 172–174 | Genomes: by region |
| Dynamic nodes and templates | data nodes that arrive with regions; a track per node | deferred | — | 176–179 | trigger: a dataset with per-region nodes |
| Name index | open on a gene by name from a file read by region | planned | [sources](sources.md) | 228 | Genomes: open on a gene |

## Layers and themes

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Letters | cells and letters, density steps, letter caps | planned | [ux-foundations](ux-foundations.md) | 116, 117, 120–122 | p53: tracks |
| Ranges | lanes, labels, collapse, split, line styles, arrowheads | planned | — | 45, 104, 105, 128, 226, 232 | p53: tracks |
| Values | bars, line, area, points; split; log scale; reference line | planned | — | 52–57, 201, 202, 214–217 | Genomes: signal, Manhattan |
| Variants | sites, samples as rows, long deletions solid | planned | — | 222 | Genomes: samples' variants |
| Pairs | arcs above and below, kinds, chords on rings | planned | — | 62, 63, 129 | tRNA: arcs |
| Matrix | positions × states, majority and maximum density | planned | — | 60, 61, 220 | Genomes: methylation |
| Grid and dot plots | a value per pair of positions; crosshair hover | planned | — | 58, 59, 109, 145 | tRNA: base-pair dot plot |
| Logo | stacked letters by bits or probability | planned | — | — | Globins: family |
| Secondary structure | helices and strands, cut ends square | planned | — | 67, 130 | p53: tracks |
| Ruler and numbering | tabular, skipped labels, absent or placed | planned | [ux-foundations](ux-foundations.md) | 234 | every view |
| Frames | outlined regions across tracks, open where cut | planned | — | 100 | Globins: ESPript figure |
| Rows layer | rows in tree or data order, filters, data columns | planned | — | 106, 107 | Globins: family |
| Fallbacks | every layer type declares one; readers degrade predictably | planned | [extensibility](extensibility.md) | 74, 238–244 | Nightingale reference |
| Themes | uniform, field palettes (binned, continuous, log), rules, schemes | planned | [ux-foundations](ux-foundations.md) | 133–137, 203–205, 215 | p53: themes |
| Residue schemes | Clustal X, Zappo, Taylor, hydrophobicity, nucleotides, ESPript | planned | [ux-foundations](ux-foundations.md) | 135 | Globins: family |
| Shared domains | domain groups across layers | planned | — | 56 | Genomes: signal |
| Legends | ramps and steps, merged, painted vs vocabulary | planned | [ux-foundations](ux-foundations.md) | 132, 204 | p53: themes |
| Size channels | thickness, glyph, radius, width, opacity | planned | — | 137, 242 | p53: tracks |
| Layers across sections | ribbons between sections, a dot plot with its own y section | deferred | [rendering](rendering.md) | 209, 210 | trigger: synteny views; Genomes would carry them |
| Contact maps | binned grids by region and level (Hi-C) | deferred | — | 223 | trigger: layers across sections |
| MAF over a genome | a multiple alignment on a genome axis | deferred | — | 213 | trigger: a dataset that needs it |

## Views and layouts

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Tracks layout | sections, tracks, groups, summaries | planned | [building-blocks](building-blocks.md) | 101–105 | p53: tracks |
| Text layout | wrapped groups of ten, lanes, inline labels | planned | — | 94–99 | p53: text figure |
| Alignment figures | rows in text, ESPript style, numbers per row | planned | — | 97 | Globins: ESPript figure |
| Blocks | side by side, stacked, many blocks, two assemblies on one axis | planned | — | 102, 115, 210 | Genomes: both genomes |
| Rings | tracks as rings, labels with leaders, a window | planned | — | 108, 114 | Plasmid: ring with linear detail |
| Navigators | strip, miniature, ideogram, vertical overview | planned | — | 110 | Genomes: by region |
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
| Undo and redo | every spec edit; previews end in a patch or a cancel | planned | [ux-foundations](ux-foundations.md) | 163, 164 | every editing view |
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
| Ports with default adapters | fetch, archive, format and map runtimes, structure runtime, text metrics, backends, writers, host, compute, storage, clipboard, links, services; each with a suite | planned | [ports](ports.md) | — | Extensibility page |
| Configurable layout and spacing | every foundation number a setting, in presets; layouts as an extension point with a suite | planned | [ux-foundations](ux-foundations.md) | — | foundation stories, a preset switcher |
| Settings in layers | defaults, viewer configuration, user settings; effective values with origin | planned | [extensibility](extensibility.md) | — | Extensibility page |
| Generated settings form | from the settings' description | planned | [extensibility](extensibility.md) | — | Extensibility page |
| Extension descriptors | passed at composition; conflicts are errors | planned | [extensibility](extensibility.md) | — | Extensibility page |
| Node kinds and custom fields | namespaced, carried through untouched | planned | [extensibility](extensibility.md) | 235 | Extensibility page |
| Every registry open | layer types, themes, formats, map runtimes, queries, commands, backends, chrome parts | planned | [extensibility](extensibility.md) | — | Extensibility page |

## Integration

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Nightingale reference | a subset drawn by Nightingale from the resolved view | planned | [rendering](rendering.md) | 183 | Nightingale reference |
| Mol\* structure link, as a showcase | translations between a sequence's own system (often position-indexed: a lab's predictions per letter) and mmCIF `label` and `auth` numbering, insertion codes and unobserved residues; hover and selection both ways; per-position values colored in 3D. The structure runtime is application-layer code (the storybook package first; a package of its own when another host needs it) | planned | [maps](maps.md) | 3–5, 25, 118, 167, 181 | Storybook: predictions onto a structure |
| Addressing | local, hierarchical ids; addresses as `path/(query :key value)` (meso\*'s grammar) for nodes, locations in labels, elements and rects; commands take addresses; copy link | planned | [addressing](addressing.md) | 1–8, 182 | Storybook: addressing |
| Markdown links | links in Mol\*'s grammar driving the view; one href for Mol\* and Seq\* | planned | [addressing](addressing.md) | 182 | Storybook: addressing |
| Agent tools | read and write a spec through addresses and expressions; the skill | deferred (the author, 2026-11-17) | [addressing](addressing.md) | — | trigger: the language description and addressing done |
| React adapter | the default UI in React | deferred | [packages](packages.md) | — | trigger: a React host |
