# Features

## Intent

The scope of this version: what it promises, one row per feature. The designs say how things work, the edge-case
catalog ([edge-cases.md](edge-cases.md)) which behaviors must survive, `packages/seq/TODO.md` what is next. This says
what the version contains, so scope never grows or shrinks silently (the author, 2026-11-17: "ok" to the proposal).

## Rules

- **Status**: `planned` (in this version, not built), `building`, `done`, `deferred` (not in this version, with its
  trigger). A feature is `done` when its example view runs ([examples.md](examples.md)) and every edge case it covers
  is named by a test or has an outcome (checked: `.scaffolds/seq/_checks/ec-ids.ts`, §ec-tested-by-name). A `building` row says after a colon what is built so far (the first three slices, lanes named where it helps); the rest of its What is not.
- **Moved with the merges** (§features-move-with-merges): the lane that builds a feature moves its row here and its
  line on the Storybook's Design page "Features" (with its story's link) in its own change.
- **Changes are decisions**: moving a row to or from `deferred` is the author's call, with the reason in its row.
- **Sources**: the previous implementation's features page, the designs here, the JBrowse reply. A feature the
  previous implementation had that is not listed is an omission to fix, not a decision.
- **Provisional**: the first pass below is a proposal. `planned` rows are the recommended scope of this version; the
  `deferred` rows are the ones to review first.
- **At the scaffold's close**, this list becomes the Storybook's Features page, each row linking its view; the Cases
  column goes with the catalog it cites.
- **Version 1 and later milestones**: what version 1 contains and the order of work is [version-1.md](version-1.md);
  rows deferred to the Genomes milestone come back together.

Columns: the feature; what it does; (in Streaming and patches) what undo does to it; its status; the design; the edge cases it covers (EC numbers, catalog); the view
that shows it (dataset: view).

## Language and spec

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| One tree, one node shape | `{ kind, id, params, children }`, addresses, refs | done | [format](format.md) | 68–71, 77 | every view |
| Language as data | one description: validation, defaults, builder types, docs, agent skill | building: description, validation, defaults, builder types; docs and agent skill not | [format](format.md) | 86–88 | — |
| Typed builder | the main way to write a spec; misuse is a type error | done | [format](format.md) | — | every view; Features › Builder (six stories, the code beside the view) |
| S-expression printing and reading | inspection, docs, tests, agent input; never stored | done | [format](format.md) | — | every view (printed beside it) |
| Validation, every issue | stable codes, hints, counts, severity with one meaning | building: codes, hints, severity; the Antibody view exists (DA), its invalid labels not shown; EC-88 and EC-89's outcomes not recorded | [format](format.md) | 6, 76–78, 86–89 | Antibody: invalid labels |
| Versions and migrations | one semver, explicit migrations, lenient minors | building: the version gate only | [format](format.md) | 72–74, 235 | — |
| Encodings | JSON by default; MessagePack for large inline data (typed columns as extension types); columns declared, hashed by bytes | building: JSON, S-expression and MessagePack (core's codec, typed columns as extensions; records read from MessagePack, FS-A1); the views' fixtures round-tripped through MessagePack not; EC-53's outcome not recorded | [format](format.md) | 53, 90 | every view (both encodings in its test) |
| Archives, over the assets port | a spec with its data files and an index, read and written as zip, entries resolved as assets | building: zip read and written over the assets port, entries resolved as assets, EC-91's bounds (`ArchiveLimits`), the Plasmid and p53 saved and opened (AR; Features › Archives); every dataset as an archive not | [ports](ports.md) | 91 | every dataset, saved and opened as an archive |
| Published JSON Schema | derived from the description | deferred | [format](format.md) | — | trigger: an outside reader or agent tool needs it |

## Coordinates and maps

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Indices and labels | 0-based indices, labels converted at the edge | done: `index` and `ordinal`, labels in queries, Kabat and IMGT on the Antibody view (DA) | [building-blocks](building-blocks.md) | 1–8 | Antibody: Kabat vs IMGT |
| Numbering schemes | offset, segments, labels, names; Kabat, IMGT | building: the `numbering` kind with `labels` (inline or a records source's column), queries, addresses, the ruler and hover read it (NB), shown on the Antibody view (DA); `offset` and `segments` (MK) | [maps](maps.md) | 2–8 | Antibody: Kabat vs IMGT |
| Cyclic systems | origin, wrapped ranges, shortest arc | building: wrapped ranges, offset and codon maps over the origin, the Plasmid views (DP: ring, features across the origin, a rotated linear view); EC-13's outcome not recorded | [building-blocks](building-blocks.md) | 9–14 | Plasmid: across the origin |
| Map kinds | identity, offset, numbering, segments, pairwise, alignment, codon, structure | building: identity, offset, codon, numbering (`labels`, `offset`, `segments`), `segments`, `structure` with per-residue status (MK), an alignment's columns as a target (AC, a `segments` stopgap); `pairwise` and `alignment` not | [maps](maps.md) | 15–30 | p53: protein on its gene |
| Map sourcing | inline, from a source by region, computed by a runtime | building: inline and by region (item maps); computed not | [maps](maps.md) | 212 | p53: structure map computed by alignment |
| Map encodings | implicit, runs, expanded by region; composed without expanding | building: implicit and runs | [maps](maps.md) | 27, 33, 34 | Globins: 10,000 rows |
| Direction on every run | inversions, reverse strands, reversed pairs rows | building: codon maps on either strand; reversed pairs rows not | [maps](maps.md) | 206 | p53: the minus-strand gene; Plasmid: features on both strands |
| Paths and projection | bounded shortest paths, the copy rule, ambiguity reported; two maps of one pair chosen by a map's `priority` or a use's `via` | building: bounded paths, the copy rule, ambiguity reported; two maps of one pair by `priority` and `via`, a section's `via` over every path (MK, MX; p53: structure) | [maps](maps.md) | 19–28, 207 | p53: structure chains |
| Closure | a nucleotide lifts to its residue, one step | done | [building-blocks](building-blocks.md) | 30, 144 | p53: protein on its gene |
| Structure systems and maps | chain instances as units, SIFTS in abstract form (the `structure` kind kept) | building: units and copies (ML1, ML2), structure maps with priority and `via` (MK), the translation library `@molstar/seqx-core/structure` (ST1), the viewer's additions and `linkStructure` (ST2), the fake structure viewer and its stories (ST3; Features › Structure link), a map onto a source's unit pending until described (FX3) | [maps](maps.md), [structure-syncing](structure-syncing.md) | 65, 118 | p53: structure chains |
| Codon tables | NCBI tables, Sec, Pyl, stops | building: NCBI 1, 2 and 11 | [maps](maps.md) | 16–18 | p53: protein on its gene |
| Pieces and layout maps | ranges, reversed pieces, collapsed introns, letter maps | building: an axis of pieces, collapsed introns (pieces, `pad`, breaks), reversed pieces, letters through a map (AX, SL); EC-31 and EC-32's outcomes not recorded | [axis](axis.md) | 31–33, 233 | p53: introns collapsed |
| Letters through a map | another system's letters on a section's axis: a protein's on its gene, a residue a codon wide, a codon split by an intron drawn in each piece; a path that pairs positions one to one or a residue to a codon (others have no cell a position: an issue) | building: a protein a codon wide on its gene (AX; Features › Coordinates); EC-31's outcome not recorded | [axis](axis.md) | 15–18, 28–31 | p53: the protein and its gene on one axis |
| A section's direction | `genomic`, `reversed`, or `follow` (the leader's): a follower defaults to `follow`, a leader to `genomic`; a reversed axis counts down and says so; reversed DNA letters as their complement | building: `genomic`, `reversed`, `follow`, `{ along }` (AX; Features › Coordinates); EC-31's outcome not recorded | [axis](axis.md) | 15, 31, 206 | p53: the minus-strand gene moving with its protein |
| Collapsed introns | a section's `pieces` (`whole`, `coding`, ranges refs), `pad`; a break a zero-width hairline on every layout, hidden items a tick on it; `section.collapse` | building: zero width on every layout (AX, SL; p53: introns collapsed); EC-31 and EC-32's outcomes not recorded | [axis](axis.md) | 31–33, 117 | p53: introns collapsed |

## Data and sources

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Sequences, units, complexes, copies | molecules shared by units, copy groups | building: sequences and units, copies merged or split as blocks (RB, DC); the `complex` and `unit` kinds (ML1, ML2), a section of a complex drawn with block forms, molecule letters, modified marks and links (MR, [molecules](molecules.md)), a pair's ends on split copies by unit (CU), the examples on the new kinds (MX: 4HHB, 1EHZ, 1FIN, 3AVE, T-DM1; Features › Molecules); EC-36 and EC-66's outcomes not recorded | [building-blocks](building-blocks.md) | 35, 36, 65–67 | Globins: hemoglobin |
| Ligands | one position or one per atom; chips in text | building: chips of one position, one per heavy atom, arcs to the protein (DC, HB; Complex: Ligands, ATP per atom); the hemoglobin's hemes as molecules, 3AVE's Fc glycan, T-DM1's payloads (`UNL`) (MX); AlphaFill's transplants and chips flowing in text (EC-99) not | [layouts](layouts.md) | 64, 99 | Globins: hemoglobin with heme; Complex: ATP and AlphaFill's transplants |
| Alignments | gaps, insert columns, numbered rows, lazy rows | building: aligned FASTA and Stockholm read, rows, logo, column measures (DG; Globins: family); EC-42's outcome not recorded | [grids-and-rows](grids-and-rows.md) | 38–44 | Globins: family |
| Trees | given (Newick) or computed; clades collapse | building: Newick read, clades folded as session state, `clades.keep` (RW, CH; Globins: family by tree); groups of tracks (EC-103) not | [grids-and-rows](grids-and-rows.md) | 43, 100, 103 | Globins: family by tree |
| Annotations, six shapes | ranges, values, pairs, matrix, grid, rows | building: ranges, values, pairs, grid, rows; matrix (RG) | [building-blocks](building-blocks.md) | 45–63 | every dataset |
| Keyed tables and joins | rows keyed by id joined to items, rows, states | deferred | — | 106, 224 | trigger: the close-out slice ([version-1](version-1.md)) |
| Gene hierarchies | gene, transcript, exon (several levels) | deferred | — | 227 | trigger: the Genomes milestone ([version-1](version-1.md)) |
| Formats in the core | FASTA, aligned FASTA, A3M, Stockholm, Clustal, Newick, TSV, CSV, JSON | building: `seq`, FASTA with its index, aligned FASTA, Stockholm, Newick, records (TSV, CSV, JSON, MessagePack, GFF3; FS-A1, FS-A2); A3M and Clustal not | [sources](sources.md) | 38–40, 49, 50, 92 | every dataset |
| More genome formats | 2bit, BigBed, VCF with tabix, PAF, GenBank | building: VCF sites with tabix (VC: Features › Genome formats › VCF), by region past level 1 with exact counts and its genotype matrix (FR: Features › Genome formats › VCF › Genotypes); BED, bedGraph, GFF3 and GTF by region with levels (FR); GenBank read whole (GB: Features › Genome formats › GenBank, Services › Your own file); 2bit, BigBed, PAF designed only (genome-formats.md appendix) | [sources](sources.md) | 208 | trigger: the bigWig runtime done; each with the dataset that needs it |
| Services at the application layer | UniProt, PDBe (SIFTS), AlphaFold DB as runtimes an application composes, not in the core; recorded fixtures offline, live optional | building: UniProt, PDBe SIFTS and AlphaFold DB as runtimes from recorded fixtures with their licences, a dropped file (SV; Examples › Services); live not | [sources](sources.md) | 51 | Storybook: an application composing services |
| Structure files at the application layer | mmCIF and BinaryCIF read (through Mol\*'s readers) for chains, `label` and `auth` numbering, insertion codes, observed residues | building: an opt-in mmCIF runtime (`runtimes/mmcif`: CIF to a structure description to molecules' nodes; 1EHZ, 1FIN, 4HHB, 3SGJ; MC), `io/structure` (DP); BinaryCIF and Mol\*'s readers not | [maps](maps.md) | 3–5, 118 | Storybook: predictions onto a structure |
| Name index | open on a gene by name from a file read by region | deferred | — | 228 | trigger: the Genomes milestone ([version-1](version-1.md)) |
| Example datasets | real public entries under permissive licences, provenance and licence as data (`sources.json`), recorded by hand, large variants generated in code; the Stockholm and Newick readers | building: p53 (with its chromosome, CS), Globins, Plasmid, tRNA, Complex, Antibody, each with its index of sources and licences; EC-42's outcome not recorded | [examples](examples.md) | 38–43 | every dataset |

## Streaming and patches

How state changes ([changes.md](changes.md)): patches edit the spec and are the only undoable changes; streams bring
data into the store and are never undone or saved. Both are increments by address with versions. The Undo column says
what undo does to each.

### Patches

| Feature | What | Undo | Status | Design | Cases | View |
|---|---|---|---|---|---|---|
| Patches, inverses, diff | five ops by address, all or nothing, `diff` correct | yes, one step a batch | building: five ops by address, all or nothing, inverses, `diff` and `merge3` (`@molstar/core/repr/patch`), shown in p53: Edit and undo; EC-79, 83 and 236's outcomes not recorded | [format](format.md) | 78–83, 236 | every editing view |
| Undo and redo | every spec edit; previews end in a patch or a cancel | — | done: Seq\*'s history port, the host's stack (`@molstar/ui/history` adapters), one step a batch | [ux-foundations](ux-foundations.md) | 163, 164 | every editing view |
| Previews | a gesture shows its change as a transient overlay over the spec, ending in one patch or a cancel | its final patch | done: a patch over the root, applied again when the root changes, dropped when it no longer applies | [changes](changes.md) | 163 | every editing view |
| Snapshots | named patches over the base, a player | no: leaving a step is navigation | done: story steps, the player and a step's view (SY; p53: Story) | [format](format.md) | 84, 85 | p53: a story |
| Patches on streamed things | a bookmark or a ref naming a dynamic node; waits while the node is not held | yes | planned | [changes](changes.md), [addressing](addressing.md) | — | p53: transcripts by region |
| Keep this region | a command that turns streamed data into spec content on purpose (an inline copy, a ref) | yes | planned | [changes](changes.md) | 175 | p53: the gene's region by region |
| Patch streams | an agent's or a host's patches applied in order through the history, each a step, refused when made against another revision | yes, each a step | building: the revision check, a stale patch rebased through `merge3` on the host's opt-in, patches in order as one edit, an agent's turn as one step; the agent's editing story not | [changes](changes.md) | 81 | an editing story with an agent's patches |
| Concurrent writers | collaborators' patches rebased over local ones | — | deferred | [changes](changes.md) | 93 | trigger: a collaborative host |
| Session state | the viewport, the selection, the hover, the cursor, open panels: transient, held by the viewer's models | no | building: viewport, hover, follow, selection, highlights, cursor, the view's history of jumps; open panels not | [changes](changes.md) | — | every view |
| Change kinds of commands | every command declares whether it patches the spec, changes session state or navigates; menus show what undo will undo | — | building: the table of change kinds, each command of the default viewer tested against it; the command palette not | [changes](changes.md) | — | the command palette |

### Streaming: the source protocol and the data store

| Feature | What | Undo | Status | Design | Cases | View |
|---|---|---|---|---|---|---|
| Self-describing sources | format, location, index, assembly, aliases, fields | no: data | building: format, location, index, aliases, fields (FS-A1); assembly not | [sources](sources.md) | 185–191, 196–198 | p53: the gene's region from FASTA with its index |
| The source protocol | a source answers what it covers (systems, lengths, levels, nodes known up front) and a region at a level of detail; one protocol for inline and streamed data | no: data | done | [sources](sources.md) | 169–175, 185–199 | p53: the gene's region by region |
| Answers carry nodes | tables plus data nodes (sequences, maps, annotations, nested sources), identity by (source, key) | no: data | done | [dynamic-nodes](dynamic-nodes.md) | 175–180 | p53: transcripts by region |
| Levels of detail | the coarsest level that fits a device column; a format's own zoom levels used when it has them (bigWig) | no: data | building: bin-width levels and value summaries; bigWig's own levels not | [performance](performance.md), [sources](sources.md) | 170, 171 | p53: conservation signal over the gene's region |
| The data store | entries by (source, region, level), versioned; pending, loaded, empty, failed; reactive per entry | no: data | done | [data-store](data-store.md) | 169 | p53: the gene's region by region |
| Prefetch and cancellation | view, margin, prefetch; grace; retries | no: data | done | [data-store](data-store.md) | 172–174 | p53: the gene's region by region |
| Eviction | byte-bounded, least recently used, never what is in view | no: data | done | [data-store](data-store.md) | — | p53: the gene's region by region |
| Settle and open domains | open value domains computed when the view settles; legends state their domain | no: data | building: settle and open value domains; legends not | [data-store](data-store.md) | 201–205 | p53: tracks |
| Export waits for data | `renderSvg` requests what its viewport needs and waits to settle, or says it is unsettled | no: data | done | [data-store](data-store.md), [rendering](rendering.md) | 125 | every view (SVG) |
| Dynamic nodes and templates | data nodes that arrive with regions (item coordinate systems, nested sources); a track per node by template | no: data | done | [dynamic-nodes](dynamic-nodes.md) | 175–180 | p53: transcripts by region |
| Item coordinate systems | an item's own system from its `map`, made on demand, keeping its id when it returns | no: data | done | [dynamic-nodes](dynamic-nodes.md) | 176, 177 | p53: transcripts by region |
| Nested sources | requested on need, bounded depth | no: data | done | [dynamic-nodes](dynamic-nodes.md) | 179 | p53: transcripts by region |
| Format runtimes | a registry entry per format, converting it onto the protocol; one shared conformance contract | no: data | building: `seq`, FASTA, `records` (tsv, csv, json, msgpack, gff3), `alignment`, `aligned-fasta`, `stockholm`; the records runtimes' conformance suite (`data/_test/records`) | [sources](sources.md), [ports](ports.md) | — | every source |
| A genome format, as proof | bigWig read by region and zoom level, a format runtime like any other (the author: "at least one as proof of concept") | no: data | building: an opt-in bigWig runtime over range reads (BW); p53: Signal from a bigWig (FX1) | [sources](sources.md) | 169–175, 186, 201 | p53: conservation signal over the gene's region |
| A synthetic asynchronous source | seeded, with a fake clock: the protocol's test source and the first slice's | no: data | done | [version-1](version-1.md) | — | tests; a Storybook story |
| Range reads of assets | the asset manager's page cache for indexed formats (a base-library change) | no: data | building: range reads in `@molstar/core` (64 KiB pages, fail on change; BW) | [base-library](base-library.md) | — | p53: conservation signal |

## Layers and themes

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Letters | cells and letters, density steps, letter caps | building: cells, letters, density steps, letters for all columns or none (EC-116, 117); EC-120 to EC-122's outcomes not recorded | [ux-foundations](ux-foundations.md) | 116, 117, 120–122 | p53: tracks |
| Ranges | lanes, labels, collapse, split, line styles, arrowheads | building: lanes, labels, overflow, line styles, arrowheads (SL; EC-128, EC-232); split and collapse not | [themes](themes.md) | 45, 104, 105, 128, 226, 232 | p53: tracks |
| Values | bars, line, area; shared domains; open domains at settle | building: bars, line and area (SL, EC-57); EC-52 to EC-56, 201 and 202's outcomes not recorded | [themes](themes.md) | 52–57, 201, 202 | p53: tracks |
| Values for genomics | points, split by sample, log scale, reference line | deferred | — | 214–217 | trigger: the Genomes milestone ([version-1](version-1.md))
| Variants | sites colored by a field | done: sites colored by a field, heads stacked, a count strip (SL; p53: variants) | [themes](themes.md) | — | p53: variants |
| Variants by sample | samples as rows, long deletions solid | deferred | — | 222 | trigger: the Genomes milestone ([version-1](version-1.md))
| Pairs | arcs above and below, kinds, chords on rings | building: arcs above and below, chords on rings, crossing pairs, linked pairs (PG, LP; tRNA: arcs), on an alignment's columns (AC), links as arcs (MR), numeric extra columns (FX3); kinds at one cell as slices (EC-63, decided) | [grids-and-rows](grids-and-rows.md) | 62, 63, 129 | tRNA: arcs |
| Matrix | positions × states, majority and maximum density | building: the `matrix` layer by region (RG: Features › Region reads › Matrix); a VCF's genotypes as classes (FR: Features › Genome formats › VCF › Genotypes) | [region-reads](region-reads.md) | 60, 61, 220 | trigger: the Genomes milestone ([version-1](version-1.md))
| Grid and dot plots | a value per pair of positions; crosshair hover | building: a section of its own, rectangles, image tiles, navigation (PG, PZ; tRNA: dot plot; Complex: PAE); two grids on one block pair (EC-59) not | [grids-and-rows](grids-and-rows.md) | 58, 59, 109, 145 | tRNA: base-pair dot plot |
| Logo | stacked letters by bits or probability | done: letters by conservation, saturation (LG, LS; Globins: family) | [grids-and-rows](grids-and-rows.md) | — | Globins: family |
| Secondary structure | helices and strands, cut ends square | done: helices, strands with their arrowheads, turns (SL, EC-67, EC-130) | [themes](themes.md) | 67, 130 | p53: tracks |
| Ruler and numbering | tabular, skipped labels, absent or placed | building: steps of 1, 2, 5, insertion labels, a numbering a param, `unnumbered`, `units` (SL, CF, UN); EC-234's outcome not recorded | [themes](themes.md) | 234 | every view |
| Frames | outlined regions across tracks, open where cut | done: frames across tracks, open at a line break (TV, EC-100; Globins: ESPript figure) | [layouts](layouts.md) | 100 | Globins: ESPript figure |
| Rows layer | rows in tree or data order, filters, data columns | done: rows, columns beside the headers, pinned rows, order by a column (RW, RC, RX, EC-106, 107; Globins: family) | [grids-and-rows](grids-and-rows.md) | 106, 107 | Globins: family |
| Collapsed clades | a clade folded to a summary row (session state), bulk commands, reveal opens folds, `clades.keep` writes the opening set | building: folds as session state, a summary row, a light wedge (RW, CH, RX; Globins: folded); groups of tracks (EC-103) not | [grids-and-rows](grids-and-rows.md) | 100, 103 | Globins: collapsed clades |
| Lazy rows | 10,000 held and 10M lazy rows at one frame time; reductions by the theme's classes; the cursor over rows never held | building: 10,000 held, 10M lazy (RW, EC-41, 44; Globins: 10,000 rows); the cursor over rows never held (RK); the cursor across block gaps (EC-154) not | [grids-and-rows](grids-and-rows.md) | 41, 44, 154 | Globins: 10,000 rows |
| Fallbacks | every layer type declares one; readers degrade predictably | building: an extension's type drawn by its `fallback` chain, said once (EX, EC-235); built-in types' declared fallbacks not | [extensibility](extensibility.md) | 74, 238–244 | the test consumer |
| Themes | uniform, field palettes (binned, continuous, log), rules, schemes | building: tokens, the theme sublanguage (field palettes, rules, schemes, scales) (TH, GD; p53: themes); EC-133, 134, 136, 203, 205 and 215's outcomes not recorded | [themes](themes.md) | 133–137, 203–205, 215 | p53: themes |
| Residue schemes | Clustal X, Zappo, Taylor, hydrophobicity, nucleotides, ESPript | done: Clustal X, ESPript, `classes`, `classes-safe` (TH, GD; Globins: family) | [themes](themes.md) | 135 | Globins: family |
| Shared domains | domain groups across layers | building: `scale` nodes (TH; p53: themes); EC-56's outcome not recorded | [themes](themes.md) | 56 | p53: tracks |
| Legends | ramps and steps, merged, painted vs vocabulary; `omit`, a compact form, `orient`, a layer's `legend: none` | building: data in the core, a part a host places, `channels: reserve | collapse` (TH, LR, DS; p53: themes); `LegendSelection` (`of`, `omit`, `maxItems`), compact and full, vertical and horizontal, a layer's `legend: none` (LE; Features › Legends) | [themes](themes.md) | 132, 204 | p53: themes |
| Color-blind substitution | a reader's setting (`theme.cvdSafe`) drawing each unsafe scheme as its safe substitute; never in the spec | building: `theme.cvdSafe`, every substitute safe (TH, GD; p53: themes); EC-136's outcome not recorded | [themes](themes.md) | 136 | p53: themes |
| Column measures | share, identity, matches-reference per column, registered for themes (Clustal X, percent identity, ESPript) | done: share, identity, matches-reference, ESPript's gaps (RW, CF, EC-135; Globins: family) | [grids-and-rows](grids-and-rows.md) | 135 | Globins: family |
| Size channels | thickness, glyph, radius, width, opacity | building: a mark shrinks, never its lane (TH, EC-137); EC-242's outcome not recorded | [themes](themes.md) | 137, 242 | p53: tracks |
| Layers across sections | ribbons between sections, a dot plot with its own y section | deferred | — | 209 | trigger: the Genomes milestone ([version-1](version-1.md))
| Contact maps | binned grids by region and level (Hi-C) | done: region grids by tile and level (RG: Features › Region reads › Hi-C); over a genome GA2's, across sections LA3's; Genomes: contact map (GX) | [region-reads](region-reads.md) | 223 | trigger: layers across sections |
| MAF over a genome | a multiple alignment on a genome axis | deferred | — | 213 | trigger: a dataset that needs it |

## Views and layouts

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Tracks layout | sections, tracks, groups, summaries | building: sections, tracks, templates, heights by content (DB); groups and summaries (EC-103) not | [building-blocks](building-blocks.md) | 101–105 | p53: tracks |
| Text layout | wrapped groups of ten, lanes, inline labels | building: lines of groups, numbers, bands, chips as words (TX, TV; p53: text figure); EC-98's outcome not recorded | [layouts](layouts.md) | 94–99 | p53: text figure |
| Alignment figures | rows in text, ESPript style, numbers per row | done: every row on every line, columns, frames (TV, RC; Globins: ESPript figure) | [layouts](layouts.md) | 97 | Globins: ESPript figure |
| Blocks | side by side, stacked, many blocks | done: row and stack, chips, copies (RB, DC, EC-102, 115; Complex: blocks) | [layouts](layouts.md) | 102, 115 | p53: structure chains |
| Two assemblies on one axis | blocks of two genomes on one axis or ring | deferred | — | 210 | trigger: the Genomes milestone ([version-1](version-1.md))
| Rings | tracks as rings, labels with leaders, a window | done: a ring with its window and a linear detail, labels outside, groups opening on hover (RB, DS, EC-108; Plasmid: ring) | [layouts](layouts.md) | 108, 114 | Plasmid: ring with linear detail |
| Navigators | strip, miniature, ideogram, vertical overview; custom pictures and a host-drawn navigator (`NavigatorModel`) | building: the ruler rail, strips, ideograms, a vertical overview, a miniature (NV, IZ, CF; Features › Navigators); custom pictures as `navigator-picture` entries and a host-drawn navigator (`NavigatorModel`; CN); EC-110's outcome not recorded | [navigators](navigators.md) | 110 | p53: the gene's region; Globins: 10,000 rows |
| The ruler rail | every linear section's default navigator, on the ruler's bottom border (a slim rail without one); a drag on a positions ruler pans, on an alignment's ruler the lower 40% | building: on every linear section, grips, a zoom strip (NV, IZ, CR); EC-110's outcome not recorded | [navigators](navigators.md) | 110 | every view; Globins: family |
| The layout interface | strips, bands and warps behind one internal interface; placement and hit-testing as its inverse; one viewport record | building: tracks, text, blocks and rings behind one interface, a registry of layouts (MV); EC-123's outcome not recorded | [layouts](layouts.md) | 123 | every view |
| Follow and scroll groups | sections follow through maps (scale 1, 3, 1/3); a move asked of a follower moves its leader (`follow.followerMoves`, the second slice) | done | — | 113, 114, 211 | p53: protein on its gene |
| Viewport persistence | keep center and zoom across updates; pending intents | done | — | 111, 112 | every view |
| Copies merged or split | `×N` blocks, hover for every member | done: `copies: merge | each` (RB, DC, EC-65, 67; Globins: hemoglobin) | [layouts](layouts.md) | 65 | Globins: hemoglobin |

## Interaction and editing

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Configurable grammar | gestures bound to commands; presets; one per viewer | building: the binding table, the `default` and `genome-browser` presets, the view's keys as config | [ux-foundations](ux-foundations.md) | 148–155 | every view; a preset story |
| Hover across linked sections | one owner, projected, emphasis outline | building: owner, projection, codon closure, the projected outline drawn; the orthologs view not | [ux-foundations](ux-foundations.md) | 138, 141–147 | p53: orthologs |
| Selection | click, extend, toggle, snapped drag, projected | building: click, extend, toggle, several ranges (add, subtract), drag snapped to items, projected, drawn as one band a section (p53: Selection); the text layout's cases (TV); snapping to codons over nucleotides (EC-149) not | [ux-foundations](ux-foundations.md) | 148–151 | every view |
| Keyboard cursor | crosses gaps, pieces, lazy rows | building: across sections, through the codon map, extending and adding ranges, rows and grids' y (RK; Globins: rows); block gaps (EC-154) not | [ux-foundations](ux-foundations.md) | 154, 155 | Globins: 10,000 rows |
| Focus and reveal | zoom to, or pan only when hidden | building: focus, reveal, a double-click's focus of an item's whole extent, the view's history of jumps (back, forward), the zoom keys; the search view not | [ux-foundations](ux-foundations.md) | 12, 158, 159, 168 | p53: search |
| Highlight channels | named, cleared independently | building: named channels, drawn and projected, a link's hover on its own; the search view not | [ux-foundations](ux-foundations.md) | 157 | p53: search |
| Search and results | matchers (exact, degenerate, regex, approximate); stepping; saved as sets | planned | — | 161 | p53: search |
| Query language | composable queries, printed and read as S-expressions | building: range, at, all, item, node, extent, items, rect, row, chain; by role (fourth slice): `elements`, `groups`, `filter` over typed properties, `expand`, `linked`, `merge`/`intersect`/`except`, Keys (QE; Features › Query); extension forms and properties not | [building-blocks](building-blocks.md), [locations-and-queries](locations-and-queries.md), [properties](properties.md) | 89, 229 | p53: queries |
| User-defined refs | named addresses in a `refs` container; bookmarks and location sets are refs; rewritten on move | building: refs, bookmarks and sets, rewritten on move, dangling refused (p53: Bookmarks and links); refs in data params and a section's `of` and `follow`, the Antibody view not | [addressing](addressing.md) | 162 | Antibody: CDR bookmarks |
| Track and group editing | hide, reorder, collapse, theme changes | building: hide (tracks, layers), reorder, theme changes as patches (p53: Edit and undo); hidden tracks listed and restored by the host (HT); collapse and groups not | [ux-foundations](ux-foundations.md) | — | every dataset's editing view |
| Copy the selection's letters | what is shown, in reading order, rows as FASTA lines (`selection.copyLetters`, Cmd/Ctrl+C) | done: as shown (`selection.copyLetters`, Mod+C) | [layouts](layouts.md) | — | p53: text figure; Globins: family |
| Hover card and context menu | the place in words (every numbering, the item, a value, a grid cell, a rectangle's summary), a keys panel; actions on the element, a keyboard menu | building: the hover card and keys panel (HC, the third slice), its items capped (`hover.cardItems`, FX3), its track named (FX4); the context menu through track header parts (`on: 'element'` entries, HP2; flyouts, HT) | [ux-foundations](ux-foundations.md), [track-header-parts](track-header-parts.md) | 156, 166 | p53: tracks |
| Track header parts | actions as registry entries in fixed slots, count chips (hidden layers, filtered), hidden tracks listed by the host, the menu behind its port, drawn or a host's DOM parts over the header geometry, the secondary click a behavior param | building: core's header line (HP1: slots, state chips, the section's line, hits `header-slot`, `header-state`, `header-label`, `section-header`); the `header-action` registry kind, built-ins, the `headers` behavior, its menu and parts ports (HP2); the stories (HP3; Features › Track headers); the menu's own events (FX6); hidden tracks the host's, flyout submenus (HT) | [track-header-parts](track-header-parts.md) | 146, 155, 156 | Features/Track headers |

## Rendering and output

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Resolved view and draw lists | public forms other renderers read | building: public forms, a test-only second consumer (NG, EC-183) | [rendering](rendering.md) | 238–244 | the test consumer |
| Canvas2D backend | glyph atlas, hiDPI, band copy | building: hiDPI, three stacked canvases; glyph atlas and band copy not | [rendering](rendering.md) | 120–123, 140 | every view |
| SVG and PNG | headless, same picture as the canvas | building: SVG, PNG (OP); EC-124's outcome not recorded | [rendering](rendering.md) | 124, 125 | every view (side by side) |
| Render CLI | a spec file to SVG or PNG | done: `seqx-render`, a spec to SVG or PNG (OP, EC-126) | [rendering](rendering.md) | 126 | — |
| Text metrics port | the same layout in the browser and Node | done | [rendering](rendering.md) | 124 | — |
| UX foundations | tokens, squircles, emphasis, motion, density steps | building: tokens, squircles, wells, emphasis, density steps, motion (`view.motion`); EC-139's outcome not recorded | [ux-foundations](ux-foundations.md) | 138, 139 | foundation stories |
| Writers | a selection as FASTA; a figure as PDF | deferred | [extensibility](extensibility.md) | — | trigger: the first extension of writers |

## Configuration and extensions

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| Ports with default adapters | assets (`AssetManager`), text metrics; formats, maps, backends as registry entries | building: assets (with the policy and the manager's conformance suite), text metrics, history (a host's stack), animation, the `hoverCard` behavior's port, `trackSizing`, the header menu and header parts ports (HP2), the rasterizer (OP); a structure viewer is a library, not a port | [ports](ports.md) | — | Extensibility page |
| Security bounds | an asset policy by host, bounds on everything a spec sizes, escaped output | building: the asset policy passing credentials, redirect, signal and a streamed size to the manager (`checkAssetManager`), the spec's and answers' bounds, tile, template and nesting limits, archives' bounds (AR, EC-91); search bounds not | [security](security.md) | 76, 86, 91, 126 | tests |
| Configurable layout and spacing | every foundation number a setting, in presets; layouts behind one internal interface | building: render and viewer config items; presets not | [ux-foundations](ux-foundations.md) | — | foundation stories, a preset switcher |
| Typed configuration | one set of config items at composition (`ParamDefinition`) | building: plain default objects; `ParamDefinition`s not | [extensibility](extensibility.md) | — | Extensibility page |
| Settings in layers | user settings over the configuration; effective values with origin | deferred | [extensibility](extensibility.md) | — | trigger: the typed configuration has users |
| Generated settings form | from the config items' `ParamDefinition` | deferred | [extensibility](extensibility.md) | — | trigger: settings in layers |
| Extension descriptors | passed at composition; conflicts are errors; `<id>:<name>` types with a fallback, the spec's `extensions`, one `extension` entry as the contract | building: `SeqExtension`, `defineExtension`, `<id>:<name>` types with fallbacks, `spec.extensions`, `extend()` at run time with its undo (EX; Features › Extensibility) | [extensibility](extensibility.md), [extension-openings](extension-openings.md) | — | Extensibility page |
| Host data (seq data ports) | a host's in-memory structure served through the source protocol (a data adapter per format, the core's runtime the default), updates reopened in place by key, annotations given back as the host's | building: `changes` and `field` on the protocol, `withAdapter`, `io/structure`, the viewer reopening in place (DP), the fake host's adapter (HD; Features › Host data); sections over a source's complex not (TODO) | [seq-data-ports](seq-data-ports.md) | 176, 185 | Storybook: Host data |
| Node kinds and custom fields | namespaced, carried through untouched | building: an extension's params plain on its types, namespaced keys on others' nodes (EX); node kinds beyond layer types not | [extensibility](extensibility.md) | 235 | Extensibility page |
| Registry entries | node kinds, layer types, themes, formats, map runtimes, backends, query functions, commands; layouts and chrome parts internal until an outside user | building: format runtimes; layer types and layouts listed by the viewer's spec (MV); `header-action` (HP2), `scheme` and extensions' layer types (EX), `navigator-picture` (CN); others not | [extensibility](extensibility.md) | — | Extensibility page |

## Integration

| Feature | What | Status | Design | Cases | View |
|---|---|---|---|---|---|
| A second consumer of the resolved view | a subset drawn by a test-only consumer from the resolved view alone (Nightingale dropped as a dependency, 2026-11-19; a storybook showcase on it as a dev dependency built, NS, 2026-11-21) | building: the test consumer (NG: from the resolved view and the store alone, the field test EC-183) | [rendering](rendering.md) | 183 | — |
| Mol\* structure link, as a showcase | translations between a sequence's own system (often position-indexed: a lab's predictions per letter) and mmCIF `label` and `auth` numbering, insertion codes and unobserved residues; hover, selection and focus both ways; per-position values colored in 3D. Seq\* gives a library: translation on units (`@molstar/seqx-core/structure`) and `linkStructure` over the viewer's events; the vocabulary is the host's adapter's (the fake structure viewer's in the Storybook first; Mol\*'s with the new mol viewer) | building: hover, selection and focus both ways, chain maps by `priority` and `via`, a structure colored by a track, on a fake structure viewer (ST1–ST3; Features › Structure link), services onto a structure (SV; Examples › Services); Mol\*'s adapter with the new mol viewer not | [maps](maps.md), [structure-syncing](structure-syncing.md) | 3–5, 25, 118, 167, 181 | Storybook: predictions onto a structure |
| Addressing | local, hierarchical ids; addresses as `path/(query :key value)` (meso\*'s grammar) for nodes, locations in labels, elements and rects; commands take addresses; copy link | building: addresses and queries, refs, links, copy link (p53: Bookmarks and links); rects not | [addressing](addressing.md) | 1–8, 182 | Storybook: addressing |
| Markdown links | links in Mol\*'s grammar driving the view; one href for Mol\* and Seq\* | building: links in Mol\*'s grammar driving focus, reveal, select, highlight and hover (Storybook: links); one href driving Mol\* too not | [addressing](addressing.md) | 182 | Storybook: addressing |
| Agent tools | read and write a spec through addresses and expressions; the skill | deferred (the author, 2026-11-17) | [addressing](addressing.md) | — | trigger: the language description and addressing done |
| React adapter | the default UI in React | deferred | [packages](packages.md) | — | trigger: a React host |
