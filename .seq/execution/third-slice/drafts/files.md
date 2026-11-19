# Files as sources: tables, JSON and MessagePack, and the transforms they need

## Intent

The author, 2026-11-19: "we should have storybook example that load data from files too instead of just inlining the
values. This includes the annotation data with descriptions of the transforms that need to be applied. Similar to how
the base packages/sequences does that"; "Also json sources so we can import uniprot stuff directly? And possibly other
"well formatted" json?"; "since we support json, it should support msgpack too out of the box, and both the SoA and
AoS layouts for data"; "Maybe the msgpack stuff could be moved to core/data/??? And reexported in io?" (decided:
[base-library.md](base-library.md), the MessagePack row).

[sources.md](sources.md) says a source names its data in the data's own terms (a `format`, a `uri`, "a column-to-field
map") and that formats are converted at the edge by runtimes; the core ships "the small text formats (FASTA,
alignments, tables, Newick)". Tables have no runtime yet, the `source` node's `fields` param has no reader, and
every annotation the examples draw is inline or a `seq` file converted at recording. This design says how a file of
records (a TSV or CSV, a JSON or MessagePack document, GFF3) is read into the protocol's tables, how the spec says
what the file's columns mean (the transform), where that runs, and how every example reads its fixtures as files.

## Precedents

### Data

- **The current Seq\***: `annotation` by `uri` with `format`, or `from_source` (the file a sequence came from); `fields`
  remaps any column, positions and targets included (`{ beg: 'first', end: 'last' }`); a `numbering` the positions are
  in (`ordinal`, `kabat`); TSV, CSV (RFC 4180, cells trimmed, a header twice an error, EC-92), JSON and MessagePack
  tables as an array of row objects or an object of column arrays, read the same, columns dictionary-coded or not;
  cells text, typed by their column's role (`01` stays `01`, EC-49); `__proto__` an ordinary name (EC-50); a `uniprot`
  service converting an entry's `features[]` in code, a feature past the sequence's end left out and one running past
  it clipped, one note counting both (EC-51). Taken: the remap, both layouts, the roles, the clipping note, the numbering
  the positions are in. Not taken: the UniProt conversion as code behind a service name (the transform is written in
  the spec instead, so a reader without Seq\* sees it), `from_source` (one source node read by several data nodes is
  this codebase's shape already), a header row required.
- **MolViewSpec**: annotations `*_from_uri` with `format` (`cif`, `bcif`, `json`), a `schema` saying which columns are
  positions, `field_remapping`; JSON as an array of records. Taken: the remap as data beside the location. Not taken:
  schemas as a closed list (here the roles are the protocol's fields).
- **Vega-Lite's data**: `url` with `format: { type: 'json', property: 'features' }` (the records' path), `parse` (column
  types), then `transform: [{ filter }, { calculate }]`. Taken: a path to the records, declared types, a filter in the
  view's own expression language. Not taken: `calculate` (an expression language for new columns).
- **JBrowse**: one adapter per format (BED, GFF3, VCF, bigBed), each fixed by its format; columns mapped by the adapter's
  config, extra fields kept; jexl callbacks for anything else. Taken: presets fixed by their format (GFF3); "which file
  column is which field" stated in the document (EC-190). Not taken: code callbacks in a config.
- **IGV**: BED 0-based half-open, GFF3 1-based inclusive, `#` comment lines, a format by extension. Taken: the
  conventions as the presets' defaults. Not taken: a format guessed from the extension (the spec names it, EC-186).
- **pandas `read_csv`, d3-dsv**: `sep`, `comment`, `header=None` with `names`, `dtype`; d3's `autoType` guesses. Taken:
  the comment prefix, names for a headerless file, declared types. Not taken: guessing types ([format.md](format.md):
  typed columns are declared, not guessed).
- **Apache Arrow**: columns as typed buffers, read without a copy. Taken: a typed column already in the protocol's type
  and conventions is kept as given.
- **This codebase**: `@molstar/io`'s CSV reader (RFC 4180 over bytes, a delimiter, a comment character, a BOM skipped,
  headerless) and its MessagePack decoder; `@molstar/core`'s JSON and MessagePack tree codecs; the `seq` runtime that
  answers regions from request/answer pairs ([protocol.md](protocol.md)); `compileRowFilter` (the item form's filter,
  [themes.md](themes.md)); the assets port's Node adapter confined to `fixtures/` and the Storybook's Vite plugin serving
  the same folder ([examples-and-stories.md](examples-and-stories.md)).

### Interactions (what a user does with files in the precedents)

| Interaction | Precedent | Here |
|---|---|---|
| open a file by URL and pick its format | JBrowse, IGV | built: a `source` node (`uri`, `format`) is a patch like any edit |
| drop a local file onto the view | IGV, Jalview | refused this slice: a host's upload path ([examples.md](examples.md)'s standing tension, the fourth slice with services) |
| map columns in a dialog | JBrowse's add-track | refused: the application's; the spec's `fields` is the surface, and a host may build a form over the language description |
| features as one track a type | ProtVista, Nightingale | planned (FS-A): one annotation a type over one source (`where`, below) |
| hide features by type or source | ProtVista's filters | built: a layer's filter in `layer.data` ([themes.md](themes.md)), the "n filtered" chip |
| see the file as given beside the view | the current Seq\*'s format stories | planned (FS-A): the Files stories print the file's head and the source node beside the view |
| see why rows were left out | the current Seq\*'s notes | planned (FS-A): one issue a cause with a count and the first places (EC-86) |
| fix a misread file (0- against 1-based) | none (edit the file) | planned (FS-A): a patch of the source's transform, one undo step; the issue hints the other numbering (Decision, Issues) |
| re-read after the file changed | IGV's reload | refused: a source is read once a session; a host reloads by a patch of `uri` (a cache-busting query) |

## Decision

### A file of records is read whole, at open, into the protocol's tables

- **One runtime, the records runtime** (`RecordsRuntime`, formats `tsv`, `csv`, `json`, `msgpack`, `gff3`), a registry
  entry in its own module entry (`@molstar/seqx-core/data/records`), registered by the default viewer spec. Why one: the
  transform is the same whatever the encoding; only the step that yields records differs. Why its own entry: a reader of
  a spec never pays for it (the alignment readers' rule).
- **Read whole at open**, converted once into `seq` data (a coverage and one pair a system over the whole system at
  level 1), then answered by region by the `seq` runtime's join (`openSeq`). Why: these files are small (a fixture is
  under 1 MB; a protein's features a few hundred rows), "formats converted at the edge" ([sources.md](sources.md)), and
  one answering path is already tested. A file past the store's `maxRows` (1,000,000) or the policy's byte bound is
  refused at open (`seq.source.too_large`, with a hint: a format read by region). Files read by region (tabix, BGZF)
  wait for the Genomes milestone.
- **What a source serves**: `ranges` (rows with a start and an end, or a single position), `values` (a value a
  position), and `letters` when the source names them (`letters`, below). One source serves any of them; the reading
  node asks for its shape.
- **Its systems**: a `target` field names each row's system (BED's `chrom`, GFF3's `seqid`); without one the source has
  one system that answers to any name a reader asks by. Why: a file about one sequence (a protein's features, a table
  of scores) names it nowhere, and the reading node's target already says which; inventing a name param would be a
  second place for the same fact. Its length is the letters' when the source holds them, else unknown (the target's own
  length stands, [protocol.md](protocol.md)).
- **Rows sorted by start once**, each keeping its place in the file as `row` (its identity when it has no `id`), so an
  item cut across tiles stays one item ([protocol.md](protocol.md), Tables). Rows of one `id` with other starts are
  pieces of one item (a GFF3 CDS over several lines, an exon set).
- **Inline too**: a records source's `value` (an array of records or an object of columns, typed columns included) is
  read through the same transform. Why: a host holding parsed data in memory gives it without writing a file, and it is
  the case where typed columns are kept without a copy most often.

### The file's own syntax

- **TSV and CSV**: RFC 4180 over the UTF-8 text (a BOM skipped; LF, CRLF and CR line ends; a field starting with `"`
  quoted, `""` a quote, delimiters and line breaks inside; empty lines skipped), on `@molstar/io`'s CSV reader moved
  down to `@molstar/core` (Core additions). TSV quotes as CSV does (spreadsheets and pandas write quoted TSV). A header
  row names the columns unless `columns` names them (a headerless file); a header naming a column twice is malformed.
  `comment`: a prefix of lines read past (none by default: a `#` may start a header, BED's `#chrom`). Encodings: UTF-8
  only. Why: what UniProt, EBI, NCBI and every tool here write; a Latin-1 file reads its ASCII right and the rest as
  replacement characters, said by no check (a Tension).
- **Cells are text**, typed by their field's role or by `types`; a text cell is not trimmed (the RFC's rule, the moved
  reader's); a number or position parses past the spaces around it. EC-92 changes: the current Seq\* trimmed every
  unquoted cell; one CSV reader in the repository wins over the old rule, and only numbers needed it.
- **A row with another count of cells than the header** is read as the core reader reads it (missing cells empty, extra
  ones dropped) and said once with a count (`seq.source.records.cells`, a warning). Why: one ragged line in a
  hand-edited table should not blank a track, and the count tells the writer.
- **JSON**: `JSON.parse` over the text, a plain value (never the tree's codec: a member named `kind` would read as a
  node, a `$column` as a tagged column). **MessagePack**: `@molstar/core/data/msgpack`'s decoder with the column
  extensions of the tree's codec (`0x11`–`0x19` typed and category columns, `0x10` the numbers MessagePack drops), and
  nothing else of the codec: maps are records, never nodes. Why plain: "well formatted json" is anyone's JSON. Members
  live in null-prototype records: `__proto__` and `constructor` are ordinary names (EC-50). A document that does not
  parse is malformed with its parser's message; a JSON deeper than the engine's stack is caught as malformed, never a
  crash (EC-76).
- **GFF3**, a preset of the same reader (the format fixes what a user would otherwise write): tab-separated, nine
  columns (`seqid` the target, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`), `#` lines
  and `##` directives read past, everything after `##FASTA` ignored, `.` missing, 1-based inclusive; `attributes`
  (`key=value;…`, percent-decoded) become fields by their key (`ID` the `id`, `Name`, `Parent`, `product`, `gene`, …;
  a multi-valued one its text); rows of one `ID` are pieces of one item (GFF3's own rule for a CDS over several lines).
  `Parent` is kept as a field; no hierarchy is assembled (a transcript's exons under its gene are rows of their own,
  types filtered by `where`). Why GFF3 now: NCBI serves an INSDC record's features as GFF3, so a plasmid's features are
  read as published (Gaps). **BED** waits for a dataset that has one (the Genomes milestone): its preset would be the
  same reader with `numbering: index`, `ends: exclusive` and BED12's blocks as pieces.

### Records and columns: AoS and SoA, read the same

- `records`: the path from the document's root to its records; absent, the root. What the path reaches is either an
  **array of records** (AoS: each element an object; anything else malformed) or an **object of columns** (SoA: each
  member an array or a typed array, all of one length; else malformed). Both read the same. A table file is AoS by
  construction (a row a record).
- **A path** is segments joined by `.`: a member of an object, or an index into an array when the segment is all
  digits (`0.predicted_aligned_error`). A member whose name holds a `.` cannot be named. Why one dotted path and no
  query language (JSONPath, jq): every response the datasets read reaches its records and fields by members; a filter
  is `where`'s, not the path's. In AoS a field's path is walked in each record; in SoA a field's path is walked in the
  columns object and must end at a column.
- **Typed columns are kept without a copy** when a SoA column is already the protocol's type for its role (`Int32Array`
  for a start, an end or a position, `Float32Array` for a value, a category column for a category) and the transform
  changes nothing (the source in `index` numbering with `exclusive` ends, the protocol's own); otherwise converted once.
  Why: Arrow's lesson; a file written in the protocol's own conventions (by Seq\* itself, by a pipeline that knows it)
  costs no copy, and every other file costs exactly one.
- A field whose path reaches an object or a list is `null` in that row, said once (`seq.source.records.value`, info);
  no row is exploded into several. Why: no dataset needs it (the Proteins API's `alternativeSequence` is one string), and
  explosion is a transform with no second use.

### The transform: params of the `source` node

The transform is what the file's conventions are, written once beside its location: the source node is "the file as
read". Every param optional; a format's preset sets its own defaults. Within the language's `0.4.0` (Gaps).

| Param | What | Default |
|---|---|---|
| `fields` (exists) | the protocol's field (key) read from a column name or a path (value): `target`, `start`, `end`, `position` (a row of one position), `strand`, `id`, `value`; any other key a field of the ranges (`label`, `category`, `kind`, `score`, …), under that name | each role read from the column of its own name; every other column a field under its name |
| `records` | path to the records (JSON, MessagePack) | the root |
| `letters` | path from the root to the letters (a string) of the source's one system | none: the source serves no letters |
| `numbering` | the numbering positions are in: `ordinal` (1-based) or `index` (0-based) | `ordinal`; GFF3 `ordinal` |
| `ends` | whether a row's `end` is its last position (`inclusive`) or the one after (`exclusive`) | `inclusive` |
| `types` | a type a field (`int32`, `float32`, `float64`, `category`, `string`) | the role's type; else a text cell `category`, a JSON number `float64`, a JSON boolean `category` |
| `columns` | names of a headerless table's columns, in order | none: the first row is the header |
| `comment` | a prefix of lines read past (TSV, CSV) | none; GFF3 `#` |
| `attribution`, `licence` | the data's attribution line and licence id (decided in [examples.md](examples.md), the data's decisions 4) | none |

- **Why on the source and not on the reading node**: the conventions are the file's, whoever reads it; the store keys
  entries by source ([protocol.md](protocol.md)), so a transform on the reader would make one file several entries, and
  a runtime binds to the source, not to its readers ([sources.md](sources.md)). Several annotations from one file share
  one open, one parse and one set of entries.
- **Roles and their types are the protocol's**: `start`, `end`, `position` int32 after conversion; `strand` int8 from
  `+`, `-`, `.`, `?`, `1`, `-1`, `0` (anything else a row's strand 0, said once); `id` a category; `value` float32 (a
  JSON double rounded to float32, the protocol's column, EC-53); a position read from text or from a JSON string
  (`"175"`, the Proteins API's `begin`) parsed as an integer. Why roles fixed: the layers read them; a type param on a
  role would be a second spelling of the protocol.
- **Positions**: `ordinal` subtracts one; an `inclusive` end adds one; a `position` row is `[p, p + 1)`. Labels are
  converted once, at the edge, never mixed in one table (EC-1). The default is `ordinal` and `inclusive` because
  people's tables, UniProt, the Proteins API, GFF3, VCF and the default numbering of a reader's sequence count from one
  and name the last position (EC-195: the numbering is stated, not learned); BED and the protocol's own files say
  `index` and `exclusive`. **A declared numbering** (Kabat labels in a CDR table) waits for the language's `numbering`
  kind (Gaps): converting a label needs the target system's numbering, which a runtime does not see.
- **Values**: a `value` with a `position` a row (sparse: the positions not named are NaN, missing); a SoA or AoS with a
  `value` and no `position` is one value a row from the first position. A position named twice keeps the last, said once.
- **No expressions, no value recoding**: no computed column, no unit conversion, no string splitting, no map from one
  word to another. Why: each is a language of its own (Vega's `calculate`, JBrowse's jexl), none has a second use in the
  datasets, and a host with such a file adds a runtime (a registry entry) that converts at the edge. The one recoding the
  datasets need, UniProt's secondary-structure words, goes to the reader that already maps words to classes: the
  `secondary` layer's class table gains `Helix`, `Beta strand` and `Turn` (Gaps).

### Several annotations from one file: `where`

- `ranges` and `values` gain `where`: a form of the item filter's sublanguage (`item` by `:id` or a field, a list any
  of its values; `not`, `and`, `or`: what `layer.data`'s filter reads, [themes.md](themes.md)), written as its
  S-expression (`(item :category (Domain Region))`). The rows it rejects are not the annotation's: every reader of the
  annotation's rows (the shared item reader, the query model's items, legends and open domains, hits) reads through it,
  compiled once (`compileRowFilter`).
- **Why apart from the layer's filter**: `where` says what the annotation *is* (the writer's: "the domains of this
  entry"), a layer's filter what a user *hides* (counted by the "n filtered" chip, toggled). A track of UniProt's domains
  built from a layer filter would read "380 filtered" for rows that were never its own. Both apply: a layer's filter
  counts only what passes `where`.
- **One per feature type**: a source of an entry's features and an annotation a type, each `source` the same node, each
  its own `where`; one read of the file serves all.

### Where the transform runs, its key and its issues

- **In the runtime, once a source is opened** (`fields`, `records`, `letters`, `numbering`, `ends`, `types`, presets):
  its output is the protocol's tables, so nothing after the store knows a file was a TSV. **`where` in the row readers**
  (tier 3, compiled once a prepare): it belongs to the reading node, and the store's entries are the source's.
- **The key**: an opened source is keyed by its address and the hash of its params, so editing the transform (a patch:
  one undo step, previewable) opens the source again and its tiles are new entries; the old ones go with the detach
  ([data-store.md](data-store.md)). A preview of a transform's edit opens the patched source beside the live one, held
  by its own holder and released when the preview ends. `where` is in the reading node's hash: an edit of it re-reads
  no file.
- **Issues** (`seq.source.records.*`, each once a cause with a count and the first ten places, EC-86; severity as
  everywhere, EC-87):
  - the file cannot be read as its format (no header, a document that does not parse, a `records` path that reaches
    neither form, SoA columns of unequal lengths, a header naming a column twice): **malformed**, the source
    unavailable, its readers saying so (EC-72). Why: a structural misreading would draw rows in the wrong places.
  - a field's column or path not found in any record (`.field`, an error naming the field and listing the columns there
    are, or the first record's members): the source unavailable. Why: a missing `start` reads every row as nothing, and
    the list of columns is the fix (EC-6: errors must teach).
  - a row-level misfit, the row left out: a position or value not a number (`.value`), an end before its start
    (`.order`), a range starting past the letters' end (`.range`; one running past it is clipped and counted in the same
    issue, EC-51; checked where the source holds its letters, else by the target's length as an inline range is). A
    range one past the end when every row would fit one lower hints `numbering: index`; a first position 0 with
    `ordinal` hints the same. Why left out, not refused: one bad row in ten thousand should not blank a track, and the
    count says how many.
- **Bounds**: rows by the store's `maxRows`, bytes by the assets policy, a path's segments at most 32, a JSON's depth by
  the engine (a failure caught); a SoA column's length is checked against the others before any is read.

### Held data from files: trees and pairs read whole

A tree, pairs and a grid are held whole in the static model, inline today ([grids-and-rows.md](grids-and-rows.md)), so a
Newick file or a base-pair table cannot be a fixture read as a file.

- `tree` gains `source` (instead of `newick`: a `newick` source, the alignment readers' `readNewick`); `pairs` gains
  `source` (instead of `rows`: a records source with the pairs' fields, `a_start`, `a_end`, `b_start`, `b_end`, `type`,
  `id`, or one-position ends `a`, `b`; the same transform); `grid` the same rule, built when a dataset reads a grid from
  a file.
- **A whole read**: such a source is asked once for its whole shape (`OpenedSource.whole(shape, signal)`, a shape
  `tree`, `pairs` or `grid`), kept in the store as one entry (`<source>|whole|<shape>`), counted in bytes, never evicted
  while a shown layer reads it. The node's entry in the static model is pending until it arrives; the layers that read
  held data (rows and the tree, pairs, grid, the vertical overview) read it through one accessor, inline or the store's,
  at tier 3, and draw pending as pending (EC-169). Why not by region: these kinds are read whole by their layers
  (an order, an index on both ends), and region reads of pairs and grids are the Genomes milestone's (EC-208).
- Why not dynamic nodes: a node a source brings cannot carry the spec's own refs (`of`, `a`), and a layer would name
  `(node :key …)` addresses the static resolved view does not place.

### MessagePack everywhere a JSON file is read

- The `seq` format reads a MessagePack file too, told from JSON by its first byte ([format.md](format.md)): the medium
  lane's "every view in both encodings" then reads its fixtures both ways.
- **The assets port reads bytes**: `AssetReader` gains `readBytes(uri, signal)` beside `readText`; the policy adapter,
  the Node fixtures adapter, the export's refusing stub and the test adapters implement it. Why required, not optional:
  every adapter is the package's own, and a runtime should not test for a method.

### Examples and stories read files

- **Fixtures are read as published where the core reads the format** (amends [examples.md](examples.md)'s rule): TSV,
  CSV, JSON, MessagePack and GFF3 join FASTA, Stockholm and Newick. A JSON fixture is **trimmed, never reshaped**:
  members and records removed, none renamed or moved, so the spec's paths read the live response too and a showcase's
  switch to the service changes only `uri`. `seq` stays for what has no reader (mmCIF, GenBank, AlphaFill's ligands as
  systems) and for generated files read by region.
- **p53's entry** (FS-A): UniProt P04637's REST entry JSON (`sequence`, `features[]`; CC-BY-4.0, recorded by hand),
  trimmed to the members the views read; the p53 protein's letters, domains, secondary structure and sites read from it
  (one source, an annotation a feature type by `where`); the Proteins API's variants trimmed to UniProt's and
  ClinVar's sources as JSON (its `features[]`, `begin` and `end` strings parsed as positions). The generated gene stays
  `seq`.
- **Builders spell the transform**: a helper in `seqx-examples` writes a UniProt entry's source node and its
  annotations (`uniprotEntry(b, uri)`), so the views do not repeat it and the exported spec still holds it written out.
  No service name in the spec ([ports.md](ports.md): fetching is the showcases' registry entry).
- **Read the same in both hosts**: in Node through `node-assets.ts` (bytes too), in Storybook through the assets port
  and the Vite plugin serving `fixtures/`; tests never reach the network (EC-126's headless policy).
- **`Seq* (next)/Features/Files/<Story>`**: a table (TSV and CSV), JSON in both layouts, MessagePack in both layouts,
  a UniProt entry, GFF3, a headerless table, a misread file and its issues; each prints the file's head and the source
  node beside the canvas and the SVG, and has a headless twin reading the same files.

## The interaction matrix's cells

Files add no gesture: what a file holds is drawn by the layers that read it, as inline data is.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** | as its layer's: a row read from a file is an item (`id`, else `row`) like an inline one | as its layer's | as its layer's | as its layer's | as its layer's | — | — | no: a filter hides | no: data is not edited in the viewer |
| **a track or layer** | as built | as built | no: a track is added by `track.add` | as built | as built | — | — | a layer's filter in `layer.data` (the user's, counted); `where` is the writer's, not a gesture | the source's transform and the annotation's `where`: a patch each, one undo step, previewable, the source opened again |
| **a bookmark or link** | — | — | as built | as built | as built | — | — | — | no new cell: a bookmark of `(item :id X)` follows the id across a transform's edit; one of positions stays on positions (a Tension) |
| **a section, a follower, a range, the selection, the view, a navigator** | no change: positions are positions whatever file they came from | | | | | | | | |

## Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| the records runtime | registry entry (default viewer spec) | registered | out of the box | a host leaves it out |
| `source.fields` | spec param | each role from the column of its name; other columns fields | a file in the protocol's names needs no map | any column or path a role |
| `source.records` | spec param | the root | a table file's rows; a JSON array | a path |
| `source.letters` | spec param | none | most record files hold no letters | a path |
| `source.numbering` | spec param | `ordinal` | people, UniProt, GFF3, VCF count from one (EC-195) | `index`; a declared numbering later |
| `source.ends` | spec param | `inclusive` | the same files name their last position | `exclusive` (BED, the protocol) |
| `source.types` | spec param | the role's type; text `category`, JSON numbers `float64` | declared, not guessed | `int32`, `float32`, `float64`, `category`, `string` |
| `source.columns` | spec param | the header row | most tables have one | names for a headerless file |
| `source.comment` | spec param | none (GFF3 `#`) | a `#` may start a header | any prefix |
| quoting, line ends, encoding | fixed | RFC 4180 for CSV and TSV; LF, CRLF, CR; UTF-8, BOM skipped | the core reader's; what the services write | a host's runtime |
| `ranges.where`, `values.where` | spec param | none | what the annotation is | any item form |
| a path's segments | fixed | at most 32 | bounded work (security) | — |
| issues kept | fixed | one a cause, the first 10 places (EC-86) | the spec's own bound | — |
| rows a file | store setting (`maxRows`) | 1,000,000 | the protocol's bound | larger by the host |
| a whole entry's eviction | fixed | never while a shown layer reads it | a held kind has no coarser level | — |

## Tensions

- **The source holds the transform, so one file read two ways is two sources**: a TSV read once 1-based and once
  0-based (two writers' conventions in one file cannot happen; one file read for two systems by two `target` columns
  can) opens and parses twice. Accepted: the asset manager shares the download; a parsed-file cache keyed by bytes and
  format waits for a measurement that asks.
- **Declarative transforms against files that need computing**: a position in a string (`"12..45"`), a score in
  another unit, a value spread over columns. Refused here: a host adds a runtime. Standing: if datasets keep needing
  small computations, a minimal expression form is the next step, a language change of its own.
- **A bookmark of positions across a transform's edit**: switching `numbering` moves every row by one; a bookmark of
  positions stays where it was, now one off the feature it named. Options: refs by item (`(item :id X)`, follows),
  positions (stays). Recommendation: as the language does elsewhere (positions are positions); the issue hints at the
  edit, and the Files stories bookmark by item.
- **UTF-8 only**: a Latin-1 table misreads its accents silently. Options: an `encoding` param (a second use is not in
  the datasets), a check for replacement characters with an info issue. Recommendation: the info issue when the
  decoded text holds U+FFFD, no param.
- **GFF3 without its hierarchy**: a gene's mRNAs and exons are flat rows; a transcript drawn as one item needs the
  `Parent` chain assembled. Open until a dataset draws transcripts from GFF3 (the Genomes milestone).
- **Held data read whole is tier 3**: a tree from a file orders rows only once it arrives, so the rows layer draws
  pending first, where an inline tree orders at once. Accepted: the store's pending rule; the alternative (reading held
  files before the model) would make the static model wait on I/O.

## Self-check

The claim: **this design covers what a user does with files as sources.**

1. **Assume it holds.** A file of records in any of five formats, in either layout, becomes the protocol's tables by a
   transform written beside its location; several annotations share it; trees and pairs read whole; every dataset's
   fixtures are files read in both hosts; issues teach.
2. **What it rests on, negated.** (a) *The writer knows the file's conventions.* Negated: the writer does not know
   whether positions are 0- or 1-based, or which column is which. (b) *A file is read*. Its inverse: data written out
   to a file (an edited annotation, a selection's rows). (c) *A file is about one dataset, whole and small.* Its plural
   and its scale: one annotation over many files (a chromosome a file), and one file too large to read whole. (d) *The
   file is fixed.* Negated: it changes while viewed.
3. **The negated world.** (a) gives a misread file drawn one position off with no error, and a missing column read as
   nothing; it asks for the file as given beside the view, and issues that name the fix. (b) gives the spec exported
   with its sources and their transforms (the data stays in its file), but no way to write a table from Seq\*. (c) shows
   the source as one file and the read as whole: a genome-wide BED fails at `maxRows`. (d) shows no reload.
4. **Held both.** Integrated in the Decision: the off-by-one hints on `.range`, the missing-field error listing the
   columns there are, the Files stories printing the file's head and the source node, the misread-file story; the
   attribution and licence params carried into an exported spec. Planned: region reads (tabix, BGZF, bigBed) and BED
   with the Genomes milestone; a file written from Seq\* (an annotation as TSV or JSON, the inverse transform) with the
   fourth slice's writers, if a host asks. Standing: one annotation over many files; reload by a patch of `uri` only;
   the declarative transform's limit (Tensions).

## Consequences

- **`@molstar/core`** (base-library, each with its test and v6 changelog entry):
  - MessagePack moved down (decided): `data/msgpack` (encoder, decoder) and `repr/codecs/msgpack` (the tree's codec),
    `@molstar/io` re-exporting both at their old paths; the codec's column extensions exported for a plain value's
    decoding (a decoder hook beside the codec), a public addition.
  - **The CSV reader moved down** (a gap): `@molstar/io`'s RFC 4180 reader and its byte tokenizer to
    `@molstar/core/data/text` (a delimiter, a comment character, headerless), `@molstar/io` re-exporting it; non-breaking.
- **`seqx-core`**:
  - `data/records.ts` (its own entry): the runtime, the readers (delimited, JSON, MessagePack, GFF3), paths, the
    transform, conversion into `seq` data; `seq.source.records.*` issues; the `seq` runtime reading MessagePack;
    `AssetReader.readBytes` and the policy adapter's. Lane FS-B: `OpenedSource.whole`, whole entries in the store, the
    held accessor the rows, pairs and grid layers read through.
  - The language (`0.4.0`): `source`'s `records`, `letters`, `numbering`, `ends`, `types`, `columns`, `comment`,
    `attribution`, `licence` (`fields` gains its meaning); `ranges.where`, `values.where`; `tree.source`,
    `pairs.source` (FS-B), `newick` and `rows` then optional, one of the two required. The `secondary` class table
    gains UniProt's words.
  - Tests: the formats' conformance suite (coverage, region answers, empty and edge regions, cancellation; the format
    readers' shared contract, [ports.md](ports.md)); **one table in every encoding** (TSV, CSV, JSON and MessagePack in
    both layouts, inline) answering equal tables, a property over random tables; the transform's property (a random
    table written `ordinal`/`inclusive` and `index`/`exclusive` reads the same rows); a typed SoA column in protocol
    conventions kept by identity; each issue with its count, places and hint; fuzzed text and documents read without a
    throw (EC-40, EC-92); bounds refused before reading (a SoA of unequal columns, a huge row count, a deep JSON);
    `where` applied by every reader (items, query, legends, hits) while the chip counts only the layer's filter. A quick
    representative of the open's target ([performance.md](performance.md) gains the row): 10,000 TSV rows quick,
    100,000 under 300 ms in the heavy lane (`itWithin`, best of three); the path is the open, not a frame.
  - `_spec`: Data gains "Records" (formats, layouts, paths, the transform, issues, bounds) and the `seq` format's
    MessagePack; Protocol (`readBytes`, `whole`); Language (the params, Versions); Resolve (`where`, held from a source);
    Configuration.
- **`seqx-viewer`**: the default spec registers the records runtime; nothing else.
- **`seqx-examples`**: `node-assets.ts` reads bytes; p53's UniProt entry and variants as trimmed JSON with their
  `sources.json` rows; `uniprotEntry`; the p53 views read their protein, domains, secondary structure and variants from
  them (their records change once, reviewed as a diff). `scripts/record.ts` and the `sources.json` test move to FS-A
  (the first lane that records), DG adding its entries (§commit-often: the first to merge owns a shared helper). Tests:
  anchors on the read data (P04637 is 393 residues; R175, R248, R273 inside the DNA-binding domain), drawing against
  the data in each changed view.
- **Storybook**: `Features/Files` with twins; a Design page "Files and sources" (formats, the transform table, `where`,
  issues); the examples Design page's fixture rule amended.
- **Other designs**: [sources.md](sources.md) (the records runtime, the transform on the source),
  [protocol.md](protocol.md) (`readBytes`, `whole`), [examples.md](examples.md) (fixtures read as published where the
  core reads the format; JSON trimmed, never reshaped; p53's entry; DP's features as GFF3), [base-library.md](base-library.md)
  (the CSV row), [grids-and-rows.md](grids-and-rows.md) (held kinds from a source), [themes.md](themes.md) (the
  `secondary` class words).
- **Lanes**: **FS-A** (strong): the core moves first (MessagePack, CSV: hot `package.json` files through the
  orchestrator), then the runtime, the transform, `where`, `readBytes`, the `seq` format's MessagePack, p53's entry,
  the Files stories. **FS-B** (strong, after FS-A): whole reads, `tree.source` and `pairs.source`, the held accessor
  in the layers, Globins' tree read from its Newick fixture as the first use. The dataset lanes then record trimmed
  published files where the core reads them (DP's features as NCBI's GFF3, DG's tree as Newick, DT's pairs as a
  table), spell the transforms in their builders, and keep `seq` only where no reader exists.

## Gaps for the author

1. The transform's params (`source`'s, `where`, `tree.source`, `pairs.source`) land in `0.4.0`, this slice's one
   unreleased bump, not waiting for `0.5.0`? Recommendation: yes; they are additions a `0.4.0` reader without them
   reports as unknown params, and the slice's Versions line lists them.
2. `where` on `ranges` and `values` (the writer's row filter, no chip), apart from the layer's filter? Recommendation:
   yes; a feature-type track built from a layer filter shows its foreign rows as "filtered".
3. Defaults `numbering: ordinal` and `ends: inclusive` for TSV, CSV, JSON and MessagePack? Recommendation: yes (the
   files people and the services write); the protocol's own convention is two params away.
4. Move `@molstar/io`'s CSV reader (and its byte tokenizer) down to `@molstar/core`, re-exported by io, as MessagePack
   was, accepting its rules (text cells not trimmed: EC-92 changes)? Recommendation: yes; one RFC 4180 reader in the
   repository, no private copy in `seqx-core`.
5. GFF3 as a preset this slice, so DP reads pUC19c's features as NCBI publishes them; BED with the Genomes milestone?
   Recommendation: yes.
6. Held kinds from files (FS-B: `tree.source`, `pairs.source`, whole reads at tier 3) this slice? Recommendation: yes,
   Globins' tree as its first use; no to `grid.source` until a fixture needs it.
7. Positions in a declared numbering (Kabat labels in a file) wait for the language's `numbering` kind?
   Recommendation: yes; the Antibody's CDRs are written in code as [examples.md](examples.md) already says.
8. The `secondary` layer reads UniProt's words (`Helix`, `Beta strand`, `Turn`) instead of a recoding transform?
   Recommendation: yes; a recoding param would have this one use.
