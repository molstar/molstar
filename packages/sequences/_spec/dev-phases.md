# Development phase ids

Code comments and test names still carry the ids of the phases that built them (`ST4b`, `HV1`). This
list says what each id was, so those comments stay meaningful. An id is history: it says when and why
something was built, never what it must do (that is the spec). New code needs no phase id.

Some short ids were used by two plans; both are listed. Phases are grouped by the work they belonged to.

## The first viewer and format (2026-10-09 and 10)

| Id | What it was |
|---|---|
| A1, A1b | the aligners in `@molstar/seq-core/algorithms` (biomsa vendored); A1b, its review fixes |
| A2 | the `myers` aligner, edlib's bit-vector edit distance ported |
| B1 | the format's builder (in the format plan); also the fix plan's builder phase |
| B7–B10 | biomsa fixes: symmetric start-gap correction, typed minimizer arrays, short sequences, a size guard |
| G1, G1-1, G1-2, G1d | the index for genome scale (3M nucleotides); its findings (codon maps over the CDS only, Nightingale's dense values); G1d, transcripts, fragments and joins |
| I1 | the interactions extension (vocabulary, queries, fingerprints) |
| M1–M4, M2a, M2b | the model: path search and composed paths; dependents and input hygiene; build-once indexes and loci; text formats (fix plan); M2a parsers and alignments, M2b mappings and paths (format plan) |
| N1 | the Nightingale materializer |
| P1 | palettes and dictionaries copied from Mol\*, letter dictionaries, Clustal X |
| R0 | the rule "no allocation per frame, per pointer move or per lookup" (not a phase) |
| R1, R2 | the canvas engine: surfaces, frames, copies on pan; the kernel renderers |
| S1–S4 | the viewer core: state, commands, config, spec, registration; the transform graph and document layers; the managers; layout, hit testing and text measuring. In the fix plan, S1–S3 were the spec passes (format, viewer, projection) |
| T1 | themes, label providers, navigators and presets (viewer plan); test infrastructure (fix plan) |
| U1, U2 | the React UI; the DOM UI and `<seq-viewer>` |
| X1 | the interactions viewer half (fingerprint renderer, kind themes) |

## The fix plan (2026-10-10 review)

| Id | What it was |
|---|---|
| F0 | the coordinate-systems proposal (coordinate maps, projection, closure) |
| F1 | the domain-free abstractions wave (pieces, axis links, keyed sets, encodings, location sets) |
| V1–V7 | the format package: validator and patch shapes; unavailability; columns and tables; readers and transport; validation cost; registry, issue paths and reference; load and resolve caching |
| VS1–VS4 | the viewer's state: data mode and content keys; the snapshot layer and history; lineage; lifecycle and robustness |
| CP1a, CP1b, CP2, CP3 | marks (`MarkSet`); the projection core and closure; channels on marks; projection performance and Nightingale on the shared core |
| L1–L4 | ranged units' index spaces; themes' granularity and glyphs; features, rows and wrapped sizing; tiles |
| K1a, K2–K5 | zoomed-out keys and reductions; wrapped draws; marks and canvas chrome; indexes for scale; navigator and caches |
| UI1, UI2 | UI stability and the `<seq-viewer>` contract; menus and keyboard focus |

## The rebuild (2026-10-11 and 12)

| Id | What it was |
|---|---|
| RB0, RB0b | the rebuilt packages' entries and guards; the heap-scaling test helper |
| C1, C3 | checkpoints: the calibration of estimates; the streaming wave re-estimated |
| MD1, MD2, MD3 | the model: leaf modules and text formats; coordinate maps and mappings; annotations and column roles |
| VW1–VW5, VW4a | the format package: patches and validation; tables; loading and compile; resolve; projection routes and closure |
| VC1, VC2, VC3 | the viewer's state core (scheduler, states); the document layers and routing; the axis and viewports |
| VI1, VI2 | the rebuild's interfaces with the fix plan's phases (hover shape; streaming) |
| IN1–IN4 | marks and the projector; input dispatch and base events; hover, highlight and selection channels; the keyboard map and cursor |
| LY1–LY4, LY3a, LY3b | layout and hit testing; prepare and reductions; themes and the renderers' headless entries; the region cache |
| CS1, CS1a, CS1b, CS2 | the composer and surfaces: section marks, passes, drawing; the surface and hot paths |
| RN1, RN2, RN3 | renderers: letters, ranges, rulers and the rest; density, pairs, logo, matrix; data columns and the tree header |
| BK1 | the render backend's batches |
| SV1a, SV1b | the SVG backend; PNG export |
| XR1, XR1a, XR1b | the module layout and exports; the conformance kit |
| UX1, UX2, UX3 | UI slots and helpers; the DOM UI; the React UI |
| CM1, CM2 | generic commands; `@molstar/seq-core/links` |
| TR1 | `@molstar/seq-core/translate` |
| NG1 | Nightingale on the rebuilt stack |
| IX1 | the interactions extension rebuilt |
| BL1 | the rebuilt builder |
| FR1 | the coordinate-map rewrite of placements |
| SW1 | the switch-over: the old packages deleted, the new ones renamed |
| EN1 | `percent-identity` as a named `rules` table |
| LS1, LS2 | section highlights and bookmarks (location sets); the results path and the Search part |

## Streaming (2026-10-12)

| Id | What it was |
|---|---|
| ST1 | streamed values over a genome (hg38-sized) |
| ST1b-a, ST1b-b | dynamic nodes and the stable-id rule; nested sources and track templates |
| ST2 | `seqstar-tiles/1` as a stream provider outside the core |
| ST3a, ST3a-2, ST3b | streamed letters and items drawn; items regions as feature tables; progressive loading (skeletons, error bands, retry) |
| ST4a, ST4b | stream providers for indexed genome files (`@molstar/seq-core/stream/gmod`); the zoom gate |
| ST5 | the ideogram navigator |
| ST6 | prefetch (the lead, the next zoom level, snapshots' regions) |
| ST7 | sparse composed paths (`Paths`); streamed items and values asked and drawn through a path (a genome's tracks on a protein's axis); item coordinate systems from a streamed item's `map` |

## Features (2026-10-12 and 13)

| Id | What it was |
|---|---|
| AX1 | pieces: a block's range and direction |
| CY1 | cyclic coordinate systems |
| CG1a, CG1b, CG2, CG3 | rings (the circular geometry); their ruler and marks; circular features, chords and labels; the ring's window on its linear detail, labels hit as their items |
| LK1 | axis links: `follow` and `follow_anchor` |
| LM1 | letter maps |
| KS1 | keyed sets: order and filter over key tables, `row_stats`, `redundancy` |
| EN2b | size legends, the sashimi preset, more channels drawn |
| SR1 | the search library |
| SX1 | `.svsx` archives |
| NP1a, NP1b | guide trees (`nj`, `upgma`); `realign`, `codon_rows` |
| HV1 | hover emphasis: the hovered element outlined |
| HR1, HR2 | per-row marks on split rows over their own units; hits by drawn geometry |

## SeqQL and Mol\* (2026-10-12 and 13)

| Id | What it was |
|---|---|
| QL1, QL2a, QL2b, QL2c | `@molstar/seq-core/query`: the table; the reader, binder and printer; the checker, `sel`, the skill; `[SeqQL]` blocks |
| QR1, QR2a, QR2b | the runtime: the core; modifiers, `onto`, `partners`, `each`, filters; properties and scans |
| QD1 | the `query` computed provider |
| QE1a, QE1b | evals of how well models write SeqQL |
| QM1, QM1b | SeqQL with Mol\* (`molql`, the `structure` input); results kept on chains' units |
