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
what the file's columns mean (the transform), where that runs, and how the examples read their fixtures as files.

## Precedents

### Data

- **The current Seq\***: `annotation` by `uri` with `format`, or `from_source`; `fields` remaps any column, positions
  and targets included (`{ beg: 'first', end: 'last' }`); a `numbering` the positions are in (`ordinal`, `kabat`); TSV,
  CSV (RFC 4180, cells trimmed, a header twice an error, EC-92), JSON and MessagePack tables as an array of row objects
  or an object of column arrays, read the same; cells text, typed by their column's role (`01` stays `01`, EC-49);
  `__proto__` an ordinary name (EC-50); a `uniprot` service converting an entry's `features[]` in code, a feature past
  the sequence's end left out and one running past it clipped, one note counting both (EC-51). Taken: the remap, both
  layouts, the roles, the clipping note. Not taken: the UniProt conversion as code behind a service name (the transform
  is written in the spec, so a reader without Seq\* sees it), `from_source` (one source node read by several data nodes
  is this codebase's shape already), positions by a declared numbering (the `numbering` kind exists, NB; a records runtime does not convert its labels yet, Gaps).
- **MolViewSpec**: annotations `*_from_uri` with `format` (`cif`, `bcif`, `json`), a `schema` saying which columns are
  positions, `field_remapping`; JSON as an array of records. Taken: the remap as data beside the location. Not taken:
  schemas as a closed list (here the roles are the protocol's fields).
- **Vega-Lite's data**: `url` with `format: { type: 'json', property: 'features' }` (the records' path), `parse` (column
  types), then `transform: [{ filter }, { calculate }]`. Taken: a path to the records, declared types, a filter in the
  view's own expression language. Not taken: `calculate` (an expression language for new columns).
- **JBrowse**: one adapter per format (BED, GFF3, VCF, bigBed), each fixed by its format; columns mapped by the adapter's
  config; jexl callbacks for anything else. Taken: presets fixed by their format (GFF3); "which file column is which
  field" stated in the document (EC-190). Not taken: code callbacks in a config.
- **IGV**: BED 0-based half-open, GFF3 1-based inclusive, `#` comment lines, a format by extension. Taken: the
  conventions as the presets'. Not taken: a format guessed from the extension (the spec names it, EC-186).
- **pandas `read_csv`, d3-dsv**: `sep`, `comment`, `header=None` with `names`, `dtype`; d3's `autoType` guesses. Taken:
  declared types. Not taken: guessing types ([format.md](format.md): typed columns are declared, not guessed);
  `names` and `comment` as spec params (no second use: the GFF3 preset's internals).
- **Apache Arrow**: columns as typed buffers, read without a copy. Taken: a typed column already in the protocol's type
  and conventions is kept as given.
- **This codebase**: `@molstar/io`'s CSV reader (RFC 4180 over bytes, one-character delimiter and comment, a BOM
  skipped, headerless), built on io's CIF column model (`CifField`), its byte tokenizer, `ByteStringCache`, text-byte
  token columns, `ReaderResult` and core's `Task`; io's MessagePack decoder; core's JSON and MessagePack tree codecs;
  the `seq` runtime answering regions from request/answer pairs ([protocol.md](protocol.md)); `compileRowFilter` and
  the viewer's `dataFilter` (a layer's `data` address whose trailing `(except …)` is the user's filter; `except` is
  `not` renamed by the author, 2026-11-20, [locations-and-queries.md](locations-and-queries.md);
  [themes.md](themes.md)); the Node fixtures adapter and the Storybook's Vite plugin serving `fixtures/`
  ([examples-and-stories.md](examples-and-stories.md)).

### Interactions (what a user does with files in the precedents)

| Interaction | Precedent | Here |
|---|---|---|
| open a file by URL and pick its format | JBrowse, IGV | built: a `source` node (`uri`, `format`) is a patch like any edit |
| drop a local file onto the view | IGV, Jalview | refused this slice: a host's upload path ([examples.md](examples.md)'s standing tension, the fourth slice with services) |
| map columns in a dialog | JBrowse's add-track | refused: the application's; the spec's `fields` is the surface, and a host may build a form over the language description |
| features as one track a type | ProtVista, Nightingale | built (FS-A2): a layer's `data` with a leading item form (Decision, One track a feature type); `uniprotEntry` |
| hide features by type or source | ProtVista's filters | built: a layer's trailing `(except …)` filter, the "n filtered" chip |
| see the file as given beside the view | the current Seq\*'s format stories | built (FS-A2): the Files stories print the file's head and the source node beside the view |
| see why rows were left out | the current Seq\*'s notes | planned (FS-A1): one issue a cause with a count and the first places (EC-86) |
| fix a misread file (0- against 1-based) | none (edit the file) | planned (FS-A1): a patch of the source's transform, one undo step; the issue hints the other base |
| re-read after the file changed | IGV's reload | refused: a source is read once a session; a host reloads by a patch of `uri` (a cache-busting query) |

## Decision

### A file of records is read whole, at open, into the protocol's tables

- **One runtime, the records runtime** (`RecordsRuntime`, formats `tsv`, `csv`, `json`, `msgpack`, and `gff3` as a
  preset), a registry entry in its own module entry (`@molstar/seqx-core/data/records`), registered by the default
  viewer spec. Why one: the transform is the same whatever the encoding; only the step that yields records differs. Why
  its own entry: a reader of a spec never pays for it (the alignment readers' rule).
- **Read whole at open**, converted once into `seq` data (a coverage and one pair a system over the whole system at
  level 1), then answered by region by the `seq` runtime's join (`openSeq`). Why: these files are small (a fixture is
  under 1 MB; a protein's features a few hundred rows), "formats converted at the edge" ([sources.md](sources.md)), and
  one answering path is already tested. Files read by region (BGZF and tabix in `@molstar/io`) are the genome runtimes'
  (`runtimes/bed`, `runtimes/gff`, `runtimes/vcf`: [genome-formats](genome-formats.md), the fifth slice).
- **What a source serves**: `ranges` (rows with a start and an end, or one position), `values` (a value a position),
  and `letters` when the source names them (`letters`, below). One source serves any of them; the reading node asks for
  its shape.
- **Its systems, by name, as the protocol has them** (`Coverage.systems` lists fixed names; the store keys tiles by the
  system): a `target` field names each row's system (GFF3's `seqid`, a TSV's `chrom`); a row whose target is empty, and
  every row of a source without a `target` field, is on the source's own system, named by `system` (an optional param;
  absent, the source node's id). A reader names it as a FASTA record is named: its sequence's `name` param, else its id
  (a sequence `p53` reading `system: 'P04637'` writes `name: 'P04637'`). Why a name and not "any name": two readers'
  names would make two key families for one system's tiles, and the coverage would list no fixed name. Each system's
  length is its letters' when the source holds them (`letters` reads the source's own system), else unknown (the
  target's own length stands, [protocol.md](protocol.md)).
- **Rows sorted by start once**, each keeping its place in the file as `row` (its identity when it has no `id`), so an
  item cut across tiles stays one item ([protocol.md](protocol.md), Tables). Rows of one `id` with other starts are
  pieces of one item (a GFF3 CDS over several lines, a feature written in two pieces across a ring's origin).
- **Inline too**: a records source's `value` (an array of records or an object of columns, typed columns included) is
  read through the same transform; its positions default to `base: 0` and `ends: exclusive`, as `ranges.rows` beside it
  (the language's own convention inside a spec). Why: a host holding parsed data in memory gives it without writing a
  file, and one convention inside a spec.

### The file's own syntax

- **TSV and CSV**: RFC 4180 over the UTF-8 text (a BOM skipped; LF, CRLF and CR line ends; a field starting with `"`
  quoted, `""` a quote, delimiters and line breaks inside; empty lines skipped); TSV quotes as CSV does (spreadsheets
  and pandas write quoted TSV). The first row is the header; one naming a column twice is malformed. The reader is a
  delimited-record tokenizer with plain string cells, new in `@molstar/core` (Core additions; Gaps): io's CSV reader is
  not moved, its column model being io's CIF one.
- **Cells are text**, typed by their field's role or by `types`; a text cell is not trimmed (the RFC's rule); a number
  or position parses past the spaces around it. EC-92 changes: the current Seq\* trimmed every unquoted cell; only
  numbers needed it, and the RFC is the one rule a CSV writer is held to.
- **A row with another count of cells than the header**: missing cells empty, extra ones dropped, said once with a count
  (`seq.source.records.cells`, a warning). Why: one ragged line in a hand-edited table should not blank a track.
- **UTF-8 only**: text read as UTF-8; a decoded text holding U+FFFD says so once (`seq.source.records.encoding`, info:
  "not UTF-8? Latin-1 accents read as replacement characters"). Why: what UniProt, EBI, NCBI and every tool here write,
  and no dataset needs another encoding; the info issue tells the one who has one.
- **JSON**: `JSON.parse` over the text, a plain value (never the tree's codec: a member named `kind` would read as a
  node, a `$column` as a tagged column). **MessagePack**: `@molstar/core/data/msgpack`'s decoder with the tree codec's
  column extensions (`0x11`–`0x19` typed and category columns, `0x10` the numbers MessagePack drops), and nothing else of
  the codec: maps are records, never nodes. Why plain: "well formatted json" is anyone's JSON. Members live in
  null-prototype records (`__proto__`, `constructor` ordinary names, EC-50). A document that does not parse is
  malformed with its parser's message; a JSON or MessagePack document deeper than the stack (both parsers recurse) is
  caught as malformed, as the tree codec catches it, never a crash (EC-76).
- **GFF3**, a preset of the same reader (the format fixes its columns): tab-separated, nine columns (`seqid` the
  target, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`), headerless, `#` lines and `##`
  directives read past, everything after `##FASTA` ignored, `.` missing, base 1, ends inclusive; `attributes`
  (`key=value;…`, percent-decoded) become fields by their key (`ID` the `id`, `Name`, `Parent`, `product`, `gene`, …;
  a multi-valued one its text). `Parent` is kept as a field; no hierarchy is assembled. Why GFF3 now: NCBI serves an
  INSDC record's features as GFF3, so a plasmid's are read as published (Gaps). The genome runtimes (the fifth slice,
  [genome-formats](genome-formats.md)) add **GFF3 and GTF as transcripts** (`runtimes/gff`: the `Parent` chain
  assembled by `source.features`, `transcripts` by default) and **BED** (`runtimes/bed`: `base: 0`, `ends: exclusive`,
  BED12's blocks as pieces, extra columns named by `fields` with a column number, `{ stain: '4' }`; built by BE: BED's
  `name` read as the `label` field, `strand` as an int8 field and the role); `gff3` here stays
  the default records preset, the runtime taken by `withAdapter`.
- **A feature across the origin of a cyclic system** (pUC19c's, DP): GFF3 writes it with an end past the length, or as
  pieces sharing an `ID`. The runtime keeps the row as written; the row reader, which knows the target's length and
  `cyclic`, wraps a row ending past the length into two rows of one item (the language's wrapped range, EC-9), and on a
  linear system clips it with `.range`. Why not in the runtime: without letters it knows neither the length nor that
  the system is cyclic. Pieces sharing an `ID` are one item already.

### Records and columns: AoS and SoA, read the same

- `records`: the path from the document's root to its records; absent, the root. What it reaches is an **array of
  records** (AoS: each element an object; anything else malformed) or an **object of columns** (SoA: each member an
  array or a typed array, all of one length, checked before any is read; else malformed). Both read the same. A table
  file is AoS by construction.
- **A path** is segments joined by `.`: a member of an object, or an index into an array when the segment is all
  digits (`0.predicted_aligned_error`); at most 32 segments. A member whose name holds a `.` cannot be named. Why one
  dotted path and no query language (JSONPath, jq): every response the datasets read reaches its records and fields by
  members. In AoS a field's path is walked in each record; in SoA in the columns object, ending at a column.
- **Default fields are the scalar members**: without `fields`, each role is read from the member or column of its own
  name, and every other member whose values are scalars is a field under its name; objects and lists are not fields
  unless a path names one of their scalars. Why: UniProt's `location`, `evidences`, `featureCrossReferences` would
  otherwise be null fields on every read. A path in `fields` that reaches an object or a list is `null` in that row, said
  once (`.value`, info); no row is exploded into several (no dataset needs it).
- **Typed columns kept without a copy** when a column is already the protocol's type for its role (`Int32Array` for a
  start, an end or a position, `Float32Array` for a value, a category column for a category), aligned for that type,
  needs no conversion (`base: 0`, `ends: exclusive`, the protocol's own) and holds no value the role refuses (a NaN
  start, an end before its start: a row left out makes a new column); otherwise one copy. Why: a file written in the
  protocol's conventions costs nothing, every other costs exactly one copy.

### The transform: params of the `source` node

The transform is what the file's conventions are, written once beside its location: the source node is "the file as
read". Every param optional; the GFF3 preset fixes its own. Within the language's `0.4.0` (Gaps).

| Param | What | Default |
|---|---|---|
| `fields` (exists) | the protocol's field (key) read from a column name or a path (value): `target`, `start`, `end`, `position` (a row of one position), `strand`, `id`, `value`; any other key a field of the ranges (`label`, `category`, `start_modifier`, …), under that name | each role from the member or column of its name; other scalar members fields under their names |
| `records` | path to the records (JSON, MessagePack) | the root |
| `letters` | path from the root to the letters (a string) of the source's own system | none: the source serves no letters |
| `system` | the name of the source's own system (rows with no target) | the source node's id |
| `base` | the number of the first position: `1` or `0` | `1` for a file; `0` for an inline `value` |
| `ends` | whether a row's `end` is its last position (`inclusive`) or the one after (`exclusive`) | `inclusive` for a file; `exclusive` for an inline `value` |
| `types` | a type a non-role field (`int32`, `float32`, `float64`, `category`, `string`) | a text cell `category`, a JSON number `float64`, a JSON boolean `category` |
| `attribution`, `licence` | the data's attribution line and licence id (decided in [examples.md](examples.md), the data's decisions 4) | none |

- **Why on the source, not the reading node**: the conventions are the file's, whoever reads it; the store keys entries
  by source ([protocol.md](protocol.md)), so a transform on the reader would make one file several entries, and a
  runtime binds to the source, not to its readers ([sources.md](sources.md)). Several annotations from one file share
  one open, one parse and one set of entries.
- **Why `base` and not `numbering`**: `numbering` names a layer's labels today and a kind later
  ([molecules.md](molecules.md)); a file's first position is a fact of the file, not a numbering of a system.
- **Roles and their types are the protocol's**: `start`, `end`, `position` int32 after conversion; `strand` int8 from
  `+`, `-`, `.`, `?`, `1`, `-1`, `0` (anything else 0, said once); `id` a category; `value` float32 (a JSON double
  rounded to float32, EC-53); a position from text or a JSON string (`"175"`, the Proteins API's `begin`) parsed as an
  integer. Why roles fixed: the layers read them.
- **Positions**: `base: 1` subtracts one; an `inclusive` end adds one; a `position` row is `[p, p + 1)`. Converted
  once, at the edge, never mixed in one table (EC-1). The file default is `1` and `inclusive` because people's tables,
  UniProt, the Proteins API, GFF3 and VCF count from one and name the last position (EC-195: stated, not learned).
  **Labels of a declared numbering** (Kabat in a CDR table): the `numbering` kind exists (NB), but a records runtime
  does not convert a label against the target system's numbering, which it does not see (Gaps).
- **Values**: a `value` with a `position` a row (sparse: the positions not named are NaN); a `value` with no `position`
  is one value a row from the first position. A position named twice keeps the last, said once.
- **No expressions, no value recoding**: no computed column, no unit conversion, no string splitting, no map from one
  word to another. Why: each is a language of its own (Vega's `calculate`, JBrowse's jexl), none has a second use in the
  datasets, and a host with such a file adds a runtime that converts at the edge. The one recoding the datasets need,
  UniProt's secondary-structure words, goes to the reader that already maps words to classes: the `secondary` layer's
  class table gains `Helix`, `Beta strand` and `Turn` (Gaps).

### One track a feature type

- **Recommended: no new param.** A layer's `data` is already an address with query segments; a feature-type track
  writes its rows as `/data/features/(item :category Domain)`. The viewer's `dataFilter` takes only a **trailing
  `(except …)`** as the user's filter (the chip counts it, the theming commands toggle it, the base is everything before
  it); today it treats every segment that is not an id as the user's, so it changes to keep a leading item form in the
  base. The rows a leading form rejects are not counted as filtered. Why: the item filter sublanguage is built,
  `compileRowFilter` reads it, and the distinction the chip needs (what the track is against what the user hides) is
  the segment's place.
- **The alternative, `where`** on `ranges` (a row filter in the same sublanguage, read by every reader of the
  annotation): needed only where queries, legends or shared domains address the annotation node itself and should see
  the subset (a bookmark `(item :id X)` on "the domains", a scale over "the domains"). With the recommended form, those
  address `/data/features/(item :category Domain)/…`, which the query model already reads. Values have no items, so
  neither applies to them: a values file is one annotation (Gaps).
- **One per type from one source**: one `ranges` node over the entry's features, a layer a type, one read of the file.

### Where the transform runs, its key and its issues

- **In the runtime, once a source is opened** (`fields`, `records`, `letters`, `system`, `base`, `ends`, `types`, the
  preset): its output is the protocol's tables, so nothing after the store knows a file was a TSV. A track's item form
  in the row readers (tier 3, compiled once a prepare), as any layer filter. Origin wrapping in the row reader (above).
- **The key, as built**: the store keys a source by its address (`DataStore.attach`); the viewer's source sync opens a
  source again when its bound entry is a new object ([data-store.md](data-store.md)). An edit of the transform (a
  patch, one undo step) changes the source node, so its entry, so the source is reopened **in place**: detached, its
  tiles dropped, opened again. **No side-by-side preview of a source's params**: a preview keeps the live sources (one
  address, one key); a transform's edit shows when applied. Why: the store needs no params hash in its keys, and a
  transform is edited rarely, against a file whose head the stories print.
- **Issues** (`seq.source.records.*`, each once a cause with a count and the first ten places, EC-86; severity as
  everywhere, EC-87):
  - the file cannot be read as its format (no header, a document that does not parse, a `records` path that reaches
    neither form, SoA columns of unequal lengths, a header naming a column twice, too deep): **malformed**, the source
    unavailable, its readers saying so (EC-72). Why: a structural misreading would draw rows in the wrong places.
  - a field's column or path found in no record (`.field`, an error naming the field and listing the columns there
    are, or the first record's members): the source unavailable. Why: a missing `start` reads every row as nothing, and
    the list is the fix (EC-6: errors must teach).
  - a row-level misfit, the row left out: a position or value not a number (`.value`), a position missing (`.missing`:
    UniProt's `UNKNOWN` modifier gives a `null` value), an end before its start (`.order`), a range starting past the
    letters' end (`.range`; one running past it clipped and counted in the same issue, EC-51; checked where the source
    holds its letters, else by the row reader against the target's length). A range one past the end when every row
    would fit one lower hints `base: 0`; a first position 0 with `base: 1` hints the same. Why left out, not refused:
    one bad row in ten thousand should not blank a track, and the count says how many.
- **Bounds**, given to the runtime: `OpenContext` gains `limits` (`maxRows`, from the store's settings, 1,000,000;
  `maxBytes`, the records bound, 64 MB), passed by the viewer's source sync. A file past `maxBytes` is refused before it
  is parsed (`seq.source.too_large`, hinting a format read by region): `JSON.parse` makes several times a file's bytes,
  and the assets policy's 256 MB an asset is a download's bound, not a parse's. Rows past `maxRows` are refused the
  same way, counted before columns are built.

### Held data from files: trees and pairs

A tree, pairs and a grid are held whole in the static model, inline today ([grids-and-rows.md](grids-and-rows.md)): a
Newick file or a base-pair table cannot be read as a file. `tree` would gain `source` (instead of `newick`: a `newick`
source, `readNewick`) and `pairs` `source` (instead of `rows`: a records source with `a_start`, `a_end`, `b_start`,
`b_end`, `type`, `id`, or one-position ends `a`, `b`; the same transform); `grid` the same rule when a dataset reads
one. Two ways, for the author (Gaps):

- **(a) A whole-read tier**: `OpenedSource.whole(shape, signal)` for a shape `tree`, `pairs` or `grid`, a store key
  family (`<source>|whole|<shape>`), counted in bytes, never evicted while a shown layer reads it; the static entry
  pending until it arrives; the rows, pairs and grid layers and the vertical overview reading held data through one
  accessor (inline or the store's) at tier 3. Costs: the protocol, the store, every held layer.
- **(b) The loader resolves them first**: before the static model is built, the context reads each `tree.source` and
  `pairs.source` whole through its runtime and gives the model their content beside the spec (the model's key gains
  their hashes); until they arrive the nodes are pending and their layers say so; no store or layer change. Costs: the
  static model is built again when a held file lands, and a held file's bytes are outside the store's budget (bounded
  by `maxBytes`).
- Recommendation: **(b)** if built this slice (smaller, held kinds stay tier 1's); and build it only if DG's tree must
  be read from its Newick file this slice, else defer with the tree inline (generated, as now).
- **Pairs from files** (the fifth slice, [region-reads](region-reads.md)): `pairs.source` is read by region (pairs tiled
  in 1D by either end; a source indexed on one end or none read whole and indexed by both ends), so (a) is not needed
  for pairs; trees and inline grids keep this choice.

### MessagePack everywhere a JSON file is read

- The `seq` format reads a MessagePack file too, told from JSON by its first byte ([format.md](format.md)).
- **The assets port reads bytes**: `AssetReader` gains `readBytes(uri, signal)` beside `readText`; the policy adapter,
  the Node fixtures adapter, the export's refusing stub and the test adapters implement it. Why required: every adapter
  is the package's own, and a runtime should not test for a method.

### UniProt's entry as it is

- The entry JSON (`sequence.value` the letters, `features[]` the records, `location.start.value` and
  `location.end.value`, `type`, `description`, `featureId`) read by a builder helper in `seqx-examples`
  (`uniprotEntry(b, uri)`) that writes the source node and its layers; the exported spec holds the transform written
  out. No service name in the spec ([ports.md](ports.md): fetching is the showcases' registry entry).
- **Modifiers kept as fields**: `location.start.modifier`, `location.end.modifier` (`EXACT`, `OUTSIDE`, `UNSURE`) as
  `start_modifier`, `end_modifier`, so a theme or a label can show an uncertain end; `UNKNOWN` comes with a `null` value
  and the row is left out with `.missing`.
- **Isoform features**: a feature whose `location.sequence` names an isoform (`P04637-2`) is not on the canonical
  sequence; the helper maps `target` to `location.sequence`, so such rows land on that isoform's system (which nothing
  reads unless a view adds the isoform) and the empty target of every other row is the entry's own system
  (`system: 'P04637'`).
- The Files stories include one of each: an `OUTSIDE` end, an `UNKNOWN` one, an isoform's feature.

### Examples and stories read files

- **Fixtures read as published where the core reads the format** (amends [examples.md](examples.md)'s rule, Gaps):
  TSV, CSV, JSON, MessagePack and GFF3 join FASTA, Stockholm and Newick. A JSON fixture is **trimmed, never
  reshaped**: members and records removed, none renamed or moved, so the spec's paths read the live response too and a
  showcase's switch to the service changes only `uri`. `seq` stays for what has no reader (mmCIF, GenBank, AlphaFill's
  ligands as systems) and for generated files read by region.
- **p53's entry**: UniProt P04637's REST entry JSON (`sequence`, `features[]`; CC-BY-4.0, recorded by hand), trimmed;
  the Proteins API's variants trimmed to UniProt's and ClinVar's sources. FS-A2's Files stories read them; the p53 views
  switch to them in DG's commit with the `PCG` change, so p53's records change once (Gaps).
- **Read the same in both hosts**: in Node through `node-assets.ts` (bytes too), in Storybook through the assets port
  and the Vite plugin serving `fixtures/`; tests never reach the network (EC-126's headless policy).
- **`Seq* (next)/Features/Files/<Story>`**: a table (TSV and CSV), JSON in both layouts, MessagePack in both layouts,
  a UniProt entry (with its modifiers and an isoform's feature), GFF3 of a ring with a feature across the origin, a
  misread file and its issues; each prints the file's head and the source node beside the canvas and the SVG, and has
  a headless twin reading the same files.

## The interaction matrix's cells

Files add no gesture: what a file holds is drawn by the layers that read it, as inline data is.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** | as its layer's: a row read from a file is an item (`id`, else `row`) like an inline one; a feature across a ring's origin one item | as its layer's | as its layer's | as its layer's | as its layer's | — | — | no: a filter hides | no: data is not edited in the viewer |
| **a track or layer** | as built | as built | no: a track is added by `track.add` | as built | as built | — | — | the user's trailing `(except …)` filter (counted); a leading item form is what the track is, not a gesture | the source's transform: a patch, one undo step, the source reopened in place; not previewed (Decision, the key) |
| **a bookmark or link** | — | — | as built | as built | as built | — | — | — | no new cell: a bookmark of `(item :id X)` follows the id across a transform's edit; one of positions stays on positions (a Tension) |
| **a section, a follower, a range, the selection, the view, a navigator** | no change: positions are positions whatever file they came from | | | | | | | | |

## Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| the records runtime | registry entry (default viewer spec) | registered | out of the box | a host leaves it out |
| `source.fields` | spec param | each role from the member of its name; other scalar members fields | a file in the protocol's names needs no map | any column or path a role |
| `source.records` | spec param | the root | a table file's rows; a JSON array | a path |
| `source.letters` | spec param | none | most record files hold no letters | a path |
| `source.system` | spec param | the source node's id | a file about one sequence names it nowhere | any name |
| `source.base` | spec param | `1` (a file), `0` (inline) | people, UniProt, GFF3, VCF count from one (EC-195); inline as `ranges.rows` | the other |
| `source.ends` | spec param | `inclusive` (a file), `exclusive` (inline) | the same files name their last position | the other |
| `source.types` | spec param | text `category`, JSON numbers `float64` | declared, not guessed | `int32`, `float32`, `float64`, `category`, `string` |
| quoting, line ends, encoding | fixed | RFC 4180 for CSV and TSV; LF, CRLF, CR; UTF-8, BOM skipped, U+FFFD said | what the services write | a host's runtime |
| GFF3's comment and columns | fixed (the preset) | `#`, the nine columns | the format's | — |
| a path's segments | fixed | at most 32 | bounded work (security) | — |
| issues kept | fixed | one a cause, the first 10 places (EC-86) | the spec's own bound | — |
| `limits.maxRows`, `limits.maxBytes` | the store's settings, passed at open | 1,000,000; 64 MB | a parse's bound, before it allocates | larger by the host |

## In the specs

| Decision | Where |
|---|---|
| One records runtime (`tsv`, `csv`, `json`, `msgpack`, `gff3`), its own entry | core › Records |
| Read whole at open into `seq` data, answered by region by the `seq` runtime | core › Read whole at open |
| What a source serves, its systems by name, rows sorted once keeping `row` | core › Fields |
| Inline `value` too, its defaults `0` and `exclusive` (decided 3) | core › Positions |
| TSV and CSV by RFC 4180 on a core tokenizer, cells not trimmed (decided 4, EC-92) | core › The syntax |
| UTF-8 only; JSON a plain value; members own only | core › The syntax |
| GFF3 a preset of the same reader (decided 5) | core › GFF3 |
| A feature across a cyclic origin | core › A row past the system's end is the row reader's |
| Records and columns, AoS and SoA read the same; a path | core › Records and columns |
| Default fields the scalar members | core › Fields |
| Typed columns kept without a copy | core › Typed columns kept without a copy |
| The transform on the `source` node; why there, why `base` | core › A source's transform |
| Roles and their types are the protocol's; positions converted once (EC-195) | core › Positions |
| No expressions, no recoding (decided 10: UniProt's own words) | core › A source's transform |
| One track a feature type by a leading item form, no `where` (decided 2) | core › Filters |
| The transform runs in the runtime once opened; the store keys a source by address; an edit reopens in place | core › Reattach |
| Issues: one a cause, a count, ten places; one severity a code | core › Issues |
| Bounds: `limits` at open | core › limits |
| Trees from files (a whole-read tier, or the loader first; decided 6) | design-only: not built (the tree stays inline) |
| Pairs from files, read by region | core › Declared, read by later lanes |
| MessagePack wherever JSON is read; the assets port reads bytes | core › format runtime |
| UniProt's entry read as it is: modifiers kept as fields, isoform features on their own system | design-only: the examples' builder helper (`seqx-examples`, no spec the check reads) |
| Examples and stories read files as published (decided 7, 8) | design-only: the examples' rule (examples.md) |
| Labels of a declared numbering in a file (decided 9) | design-only: not built (waits for a converting runtime) |
| The params in `0.4.0` with no migration (decided 1) | core › Versions |
| The interaction matrix's cells | design-only: no cell beyond the transform's edit, a patch of the source (core › A source's transform) |
| The configuration table | core › Configuration |

## Tensions

- **The source holds the transform, so one file read two ways is two sources**: one file read for two systems by two
  conventions opens and parses twice. Accepted: the asset manager shares the download; a parsed-file cache keyed by
  bytes and format waits for a measurement that asks.
- **Declarative transforms against files that need computing**: a position in a string (`"12..45"`), a score in
  another unit, a value spread over columns. Refused here: a host adds a runtime. Standing: if datasets keep needing
  small computations, a minimal expression form is the next step, a language change of its own.
- **A bookmark of positions across a transform's edit**: switching `base` moves every row by one; a bookmark of
  positions stays where it was, one off the feature it named. Recommendation: positions are positions, as everywhere;
  the issue hints at the edit, and the Files stories bookmark by item.
- **No preview of a transform's edit**: a writer fixing `base` sees the result only once applied (undo is one step).
  The alternative, a params hash in the store's source key, reopens a preview beside the live source; it changes
  `keys.ts` and the store's attach, for an edit made rarely.
- **GFF3 without its hierarchy** (closed by [genome-formats](genome-formats.md), the fifth slice; built by GG): the records preset keeps
  a gene's mRNAs and exons as flat rows; `runtimes/gff` assembles the `Parent` chain into transcripts
  (`source.features`).
- **A track's subset lives in its layer, not in the data**: with a leading item form, "the domains" is no node of its
  own; two layers that show it repeat the form. Accepted while no query or scale needs the subset as a node; `where`
  is the alternative (Gaps).

## Self-check

The claim: **this design covers what a user does with files as sources.**

1. **Assume it holds.** A file of records in any of five formats, in either layout, becomes the protocol's tables by a
   transform written beside its location; several tracks share it; every dataset's fixtures are files read in both
   hosts; issues teach.
2. **What it rests on, negated.** (a) *The writer knows the file's conventions.* Negated: the writer does not know
   whether positions are 0- or 1-based, or which column is which. (b) *A file is read.* Its inverse: data written out
   to a file. (c) *A file is about one sequence.* Negated: an entry with its isoforms' features, a GFF3 of several
   `seqid`s. (d) *The target is linear.* Negated: a ring, with features across its origin. (e) *Whole is small.*
   Negated: a file within the 256 MB policy that a parse makes gigabytes of; and the plural, one annotation over many
   files. (f) *The file is fixed.* Negated: it changes while viewed.
3. **The negated world.** (a) gives a misread file drawn one position off with no error, a missing column read as
   nothing; it asks for the file as given beside the view and issues that name the fix. (b) gives the spec exported with
   its sources and transforms (the data stays in its file), but no table written from Seq\*. (c) gives rows on systems
   no reader names, and a target that must be a name, not "any". (d) gives a feature clipped at the origin, or drawn as
   two items. (e) gives an out-of-memory parse before any row count is seen, and a genome split into files. (f) shows
   no reload.
4. **Held both.** Integrated in the Decision: the base hints on `.range`, the missing-field error listing the columns,
   the Files stories printing the file's head and a misread-file story; `target` with `system` for the empty target, the
   isoform's rows on their own system; wrapping a row past a cyclic system's length in the row reader, pieces by `ID`,
   tested on pUC19c; `limits.maxBytes` before the parse; the attribution and licence params carried into an exported
   spec. Region reads (BGZF, tabix) and BED, GFF3, GTF and VCF: the fifth slice ([genome-formats](genome-formats.md); bigBed
   designed only); a file written from Seq\* (the
   inverse transform) with the fourth slice's writers, if a host asks. Standing: one annotation over many files; reload
   by a patch of `uri` only; the declarative transform's limit (Tensions).

## Consequences

- **`@molstar/core`** (base-library, each with its test and v6 changelog entry):
  - MessagePack moved down (decided): `data/msgpack` (encoder, decoder) and `repr/codecs/msgpack` (the tree's codec),
    `@molstar/io` re-exporting both at their old paths; the codec's column extensions exported as a decoder hook for a
    plain value, a public addition.
  - **A delimited-record tokenizer** (Gaps): `@molstar/core/data/text/delimited`, RFC 4180 records with plain string
    cells (a one-character delimiter, quote and comment, a BOM skipped, line ends), written from io's parser as the
    guide; io's CSV reader unchanged (its CIF column model stays io's), moving onto it a followup if wanted.
- **`seqx-core`**:
  - `data/records.ts` (its own entry): the runtime, the readers (delimited, JSON, MessagePack, the GFF3 preset), paths,
    the transform, conversion into `seq` data; `seq.source.records.*` issues; `OpenContext.limits`; the `seq` runtime
    reading MessagePack; `AssetReader.readBytes` and the policy adapter's. The row reader wraps a row past a cyclic
    system's length. Held kinds from files: (a) or (b), if built (Gaps).
  - The language (`0.4.0`): `source`'s `records`, `letters`, `system`, `base`, `ends`, `types`, `attribution`,
    `licence` (`fields` gains its meaning); `tree.source`, `pairs.source` if built, `newick` and `rows` then optional,
    one of the two required. The `secondary` class table gains UniProt's words.
  - Tests: the formats' conformance suite (coverage, region answers, empty and edge regions, cancellation; the format
    readers' shared contract, [ports.md](ports.md)); **one table in every encoding** (TSV, CSV, JSON and MessagePack in
    both layouts, inline) answering equal tables, a property over random tables; the transform's property (a random
    table written `1`/`inclusive` and `0`/`exclusive` reads the same rows); a typed aligned SoA column in protocol
    conventions kept by identity, an unaligned or NaN-holding one copied; each issue with its count, places and hint;
    fuzzed text and documents read without a throw (EC-40, EC-92); bounds refused before parsing (bytes, rows, a deep
    document, unequal SoA columns); pUC19c's GFF3 with a feature across the origin drawn as one item in canvas and SVG.
    A quick representative of the open's target ([performance.md](performance.md) gains the row): 10,000 TSV rows
    quick, 100,000 under 300 ms in the heavy lane (`itWithin`, best of three).
  - `_spec`: Data gains "Records" (formats, layouts, paths, the transform, issues, bounds) and the `seq` format's
    MessagePack; Protocol (`readBytes`, `limits`); Language (the params, Versions); Configuration.
- **`seqx-viewer`**: the default spec registers the records runtime; the source sync passes `limits`; `dataFilter`
  takes only a trailing `(except …)` as the user's filter (its test: a leading item form kept in the base, not counted).
- **`seqx-examples`**: `node-assets.ts` reads bytes; p53's UniProt entry and variants as trimmed JSON with their
  `sources.json` rows; `uniprotEntry`; pUC19c's GFF3 if DP takes it. Anchors on the read data (P04637 is 393 residues;
  R175, R248, R273 inside the DNA-binding domain).
- **Storybook**: `Features/Files` with twins; a Design page "Files and sources" (formats, the transform table, feature
  tracks, issues); the examples Design page's fixture rule amended.
- **Other designs**: [sources.md](sources.md) (the records runtime, the transform on the source),
  [protocol.md](protocol.md) (`readBytes`, `limits`; and `whole` if (a)), [examples.md](examples.md) (the fixture rule;
  p53's entry; DP's features as GFF3), [base-library.md](base-library.md) (the tokenizer row),
  [grids-and-rows.md](grids-and-rows.md) (held kinds from a source, if built), [themes.md](themes.md) (the `secondary`
  class words; `dataFilter`'s trailing rule).
- **Lanes**: **FS-A1** (strong): the MessagePack move and the tokenizer in core (hot `package.json` files through the
  orchestrator), `readBytes`, `limits`, the runtime for TSV, CSV, JSON and MessagePack with the transform, the `seq`
  format's MessagePack, the conformance suite. **FS-A2** (fast, after FS-A1): the GFF3 preset with the cyclic wrap
  (or GFF3 handed to DP), `dataFilter`'s trailing rule, p53's entry and `uniprotEntry`, the Files stories. **FS-B**
  only if DG's tree must be read from a file this slice (Gaps), else deferred with the Newick inline. The dataset lanes
  then record trimmed published files where the core reads them (DP's features as NCBI's GFF3), spell transforms in
  their builders, and keep `seq` only where no reader exists.

## Decided by the author

2026-11-19: "Ok as recommended" (every gap as recommended).

1. The params land in `0.4.0`, this slice's one unreleased bump? Recommendation: yes. `0.4.0` is pre-release: making
   `newick` and `rows` optional (one of them or `source` required), if FS-B is built, is acceptable there with no
   migration, as every `0.x` change is.
2. One track a feature type by a leading item form in `layer.data` (`dataFilter` takes only a trailing `(except …)` as
   the user's), instead of a `where` param on `ranges`? Recommendation: yes; `where` only if a query, legend or scale
   must address the subset as a node. No `values.where` either way (values have no items).
3. File defaults `base: 1`, `ends: inclusive`; inline `value` `0` and `exclusive` as `ranges.rows`? Recommendation:
   yes.
4. The CSV reader: (a) a new delimited-record tokenizer with plain string cells in `@molstar/core`, io's reader
   unchanged; (b) move io's whole CSV closure (its CIF column model, byte tokenizer, `ByteStringCache`, token columns,
   `ReaderResult`) to core; (c) `seqx-core` depending on `@molstar/io` (the plan rules it out). Recommendation: (a);
   text cells not trimmed (EC-92 changes).
5. GFF3 as a preset this slice, so DP reads pUC19c's features as NCBI publishes them (FS-A2, or DP itself); BED with
   the Genomes milestone? Recommendation: yes, in FS-A2.
6. Trees and pairs from files this slice? Recommendation: only if DG must read its tree from a file, and then by (b),
   the loader resolving `tree.source` and `pairs.source` before the static model, not (a), the whole-read tier; else
   no, the tree stays inline.
7. `scripts/record.ts` and the `sources.json` test move from DG to FS-A2 (the first lane that records)? Recommendation:
   yes, DG adding its entries; else FS-A2 records P04637 by hand with its `sources.json` rows and DG's script adopts
   them.
8. examples.md's fixture rule amended (read as published where the core reads the format; JSON trimmed, never
   reshaped), and p53's variants and secondary structure moving from `seq` to UniProt's JSON in DG's commit with the
   `PCG` change, so p53's records change once? Recommendation: yes.
9. Positions as labels of a declared numbering (Kabat in a file) wait for a runtime that converts them (the kind exists, NB)?
   Recommendation: yes; the Antibody's CDRs are in code as [examples.md](examples.md) says.
10. The `secondary` layer reads UniProt's words (`Helix`, `Beta strand`, `Turn`) instead of a recoding transform?
    Recommendation: yes.
