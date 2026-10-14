# Data by reference, streaming and archives

How data reaches a viewer when it is not in the document: sources read through the host's loader, data streamed by
region and level through the host's stream provider, the reference tile protocol, and the `.svsx` archive. The wire
shapes are TypeScript types in `core/src/format/contract.ts` (`Loaders`, `StreamSource`, `StreamProvider`,
`StreamRequest`, `Manifest`, `DataPatch`); this file says what they mean.

## Data by reference

A `source` node says where data comes from, by at most one of: `uri` (with a `format`, or a `protocol` for a streamed
source), a registered `service` with a record `id` (and a `resource`), a `host_file`, a `value` (a whole file in place),
or a `stream` with its `data_kind` (Streamed sources). A source is fetched once however many nodes read it.

- A source's `value` is a file where a data node's `value` is parsed content: it keeps what only the file holds (record
  ids, Stockholm tags, a structure). A binary format cannot be a `value`.
- An annotation may read a `uri` directly, or `from_source`: the file a sequence, alignment or complex came from.
- Outside an archive, a `uri` is an absolute `https:` URL (`svs.source.uri` otherwise).
- A source no registered format or service reads is unavailable (`svs.provider.unavailable`, a warning); one that fails
  to load (`svs.source.load`) or lacks its record (`svs.source.record`) makes its readers unavailable the same way: one
  issue, the readers left out quietly, nothing cascading. A document opens without the extension that reads a source.

## Formats and services

Formats and services are registered by name, one flat registry each (`seqspec:format#extensions`).

- **Kernel formats**: `fasta` (molecules, or an aligned alignment), `a3m`, `stockholm`, `clustal`, `newick`, `tsv`,
  `csv`, and `json` and `msgpack` (a table or an object). Text formats read LF, CRLF and CR line ends. A Stockholm tag
  reads by the reading annotation's geometry (`GS` rows; `GC` and `GR` ranges, values or base pairs in WUSS). A Clustal
  file (`.aln`; CLUSTAL W, Clustal Omega, MUSCLE and PROBCONS write it) is one aligned alignment: its header line, then
  blocks of `name residues [count]` lines, each block listing every row in the first block's order and its rows of one
  width; a row's blocks are joined in order, and counts and conservation lines are not read. A file without the header,
  a block naming other rows, or ragged rows is `svs.source.format`, saying which. TSV and CSV cells are text, typed by
  column role.
- **Kernel service**: `uniprot` (`sequence`, `features`, `variants`), which adds the numbering `uniprot:<id>` to a
  sequence it reads and writes features and variants as ranges. Where the response gives the protein's sequence, a
  feature starting past its end (a stop lost, a frameshift's new tail) is left out and one running past it is clipped,
  one `svs.source.skipped` note counting them.
- Structure formats and services (`mmcif`, `bcif`, `pdb`, `sdf`, `mol2`; `pdbe`, `afdb`) are the Mol\* adapter's
  (`seqspec:libraries#the-mol-adapter`). `afdb`'s `pae` resource is a model's predicted aligned error, AlphaFold DB's
  JSON (format `afdb_pae`: `predicted_aligned_error`, an N × N list, and `max_predicted_aligned_error`; the older
  `residue1`, `residue2`, `distance` lists read the same), which the `molstar:pae` provider reads into a `grid` over its
  complex's polymer units in order.
- A registered format that cannot read the bytes is `svs.source.format`.

## The loader

The host's I/O port, `Loaders`: bytes for a URI or a host file (`load`), optionally a parse the host shares (one Mol\*
model for the 3D view and the complex provider), and the stream provider for a streamed source (`stream`). The format
names no server; tests pass fixture loaders. With no loaders, a file load fails (`svs.source.load`) and a streamed
source finds no provider (`svs.provider.unavailable`).

## Streamed sources

A **streamed source** serves its data by region and level, never whole: a genome's letters, its transcripts and scores,
an alignment of a million rows. The format says what is streamed and in which shapes; how it is fetched (protocol,
endpoint, auth, batching, caching) is the host's **stream provider**'s. A source streams by `stream`, an id the provider
resolves, or by `protocol` with `uri`, the address of a provider the host registered for that protocol.

| `data_kind` | Read by                                                                                                                                             | Its coordinate systems                        |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| `letters`   | a `sequence` (a manifest key by `record`), or a `complex` with no provider: a unit `R/<key>` over a sequence `R/<key>/sequence` per manifest system | the reader's: lazy molecules                  |
| `alignment` | an `alignment`                                                                                                                                      | its columns and rows                          |
| `items`     | an annotation of geometry `ranges`                                                                                                                  | `coordinate_system`, else the reader's target |
| `values`    | an annotation of geometry `values`                                                                                                                  | as `items`                                    |

Requests are JSON base types and patches are data (their tables may hold typed-array columns), so a provider may run in
a worker or behind a server. A request's abort rejects that request alone: work it shares with other requests to the
provider (a file's open, a record, a page) goes on for them, and stops only when every request waiting on it has gone.
Rejecting with `StreamNotServed` means "not served": reported once per region, never retried. Any other rejection is a
failure, retried by the viewer: a region's or a tile's is `svs.stream.load`, a warning (it is retried, and the rest
stands), where a static source's load is `svs.source.load`, an error.

## Manifests

The **manifest** answers a source's first request and is all that loading reads of it (but the rows a streamed alignment
needs whole, Streamed alignments), so layout, palettes, legends and split rows are fixed before any data arrives. It
gives, per coordinate system, its key, its extent and its numbering; the summary **levels** served, finest first; for
`values` the range of its values over the whole source; for `items` the vocabulary of each categorical field, in order;
for an alignment its row and column counts, alphabet, insert columns and whether labels are served; and dynamic nodes
known up front.

- The manifest's extents are the readers' lengths; stated counts are checked against them.
- A palette over a streamed layer takes the vocabulary's order as its order of first occurrence; a value outside the
  vocabulary takes the missing color. A palette with no field reads the field an item's default does: its category, else
  its id, else its label (so a manifest that lists ids colors items with no category).

## Requests

A request names the source, a coordinate system (a manifest key; none asks for the manifest), a region (`start`, `stop`
in that system's indices; for an alignment also rows) and a level (an index into `levels`, or `-1` for the data itself).
It is always in the source's coordinate system, never the layout's: a view's window is taken back through the
placement's route to the source's system. Through a path (`seqspec:model#paths`) that is the runs of the source's system
the path takes the window back to: the stretches its maps pair (a protein's exons), never what lies between them.

- `priority` (`view`, `margin`, `prefetch`) lets a provider order or skip requests.
- `labels: true` asks only an alignment's row labels.
- A data request may carry `max_bytes`, the bytes a provider may read before it answers with bins instead, and
  `force: true`, read it whole.

## Data patches

A **data patch** answers a request, naming what it fills (source, system, region, level) and a `kind`: `manifest`;
`items` (a `ranges` table over the region, or dense values); `bins` (one row per bin: `min`, `max`, `mean` for values;
`count` for items; the majority residue `class` and its `share` for alignments); `letters` (the region's codes); `rows`
(an alignment's rows whole, or labels). Letters sources are asked only for `letters` patches, where a device pixel holds
at most `Stream.LettersPerPx` positions.

- **Idempotent.** A patch for a region and level already held changes nothing. Patches may arrive in any order, coarse
  before fine.
- **Pending, never empty.** Each region of a level is pending until a patch fills it, then loaded (empty when it has no
  rows), or failed. A pending region is never drawn or hit as empty.
- **Placement.** Streamed annotations (items and values) are placed on blocks of their own coordinate system and,
  through a path, on any block their system reaches (a genome's tracks on a protein's axis, by its codon map). Through a
  path an item is drawn as the parts the path takes to the block (an exon's coding part, a variant's residue), and is
  hovered and selected as their hull; a value over an element covers what reaches it (a residue's codon, both parts of
  one split by an intron). The copies of a merged block reached each by a path of its own are not placed
  (`svs.layer.streamed_merged`, a warning): their regions would land on one set of positions. A lazy molecule's letters
  are its own blocks'. Of the kernel's layer types, `ranges`, `points` and `variants` draw streamed items, region by
  region (density where too many are in view), and `values` streamed values.
- An item's own coordinate map (`map`) is read into an item coordinate system (below).
- **Streamed data is never document.** Patches and the nodes they add are never written into the document, its patches,
  snapshots, history or downloads.
- A viewer that cannot stream reports what it leaves out (`svs.stream.fallback`); Nightingale's way is
  `seqspec:libraries#nightingale`.

## Levels and the zoom gate

A viewer asks for the coarsest listed level whose bin holds at most the positions of a device column, or for the data
where the view is finer than every level. An alignment's manifest may list some of its levels: a level it does not list
is asked, and drawn, as the coarsest it lists no wider on either axis, else as the one coarser by the fewest levels.
Fidelity over budget: a listed level finer than `Stream.SampleBudget` would choose is asked as it is, so a manifest that
lists only fine levels costs, zoomed out, in proportion to its rows (and its tiles stay pinned while in view). For
values and alignments drawing never goes blank: where the view's level is not held, a coarser one held is drawn. Items
and letters draw a placeholder over a pending region; items take a coarser level when the view would need more than 16
of its regions.

A provider may answer a data request with a level's bins instead (a region holding too many items, or more bytes than
`max_bytes`), estimated from its index where it can. That region is **gated**: drawn as density, its data asked again
only for smaller regions where the view is fine enough, or whole with `force` (a reader's "Load anyway"). Bins should
come from the file's own summaries or index, never from reading every item of a chromosome.

## Dynamic nodes

A patch or a manifest may also add **data nodes**: molecules, units, mappings, annotations and nested sources. A
transcripts source adds, per transcript, its protein as a `sequence`, its `codon` mapping onto the genome, and a nested
source of the protein's domains.

- A node is `{kind, key, params}`; its ref is derived, `<source ref>/<key>`.
- **Identity** is (source, key): the same key in a later patch is the same node, and its coordinate system keeps its id
  while it comes and goes (`seqspec:model#identity-and-lineage`).
- A node's refs name the document's nodes, dynamic nodes held, or nodes the same patch adds; any other is
  `svs.stream.ref`, and the node is left out.
- A node is checked as a document node of its kind. A held key with equal params is a no-op; other params are
  `svs.stream.conflict`, and the held node is kept.
- A document may name a dynamic node before it arrives: a section's block naming it is left out quietly, a document
  mapping whose end names it and an annotation whose `target` names it join the set while what they name is held (an
  annotation's per-row targets cannot wait yet), and a layer whose data names only nodes not held is pending, never
  unplaced, and its track is not shown until one arrives. A pending layer is noted (`svs.ref.pending`, info) once per
  source it waits for, with a count and the first refs.
- **Lifetime** is the viewer's: a node stays while something holds it (a region that carried it, a shown placement, the
  selection, a highlight, a held node naming it) and is evicted with what it alone named, quietly. A view shows dynamic
  nodes through track templates (`seqspec:format#track-templates`).

## Item coordinate systems

An items row may carry a `map`: a `codon` or `segments` mapping's params (`cds`, `strand`, `from_range`, `segments`, …)
from the item's own positions (a transcript's residues) onto the source's coordinate system. An item whose rows carry
one has an **item coordinate system**, `<source ref>/<id>` (a molecule with a length and no letters, and its unit
`<source ref>/<id>/unit`), mapped onto the source's system by `<source ref>/<id>.map`
(`seqspec:model#coordinate-systems`). Providers that send `map` give it; the genome-file providers here do not yet.

- **Made on demand.** It exists only while something names the item: a section's block, a document mapping's end or an
  annotation's `target` naming `<source ref>/<id>` (or its unit), or a track template with
  `for_each: { kind: "item_system" }` over the source (every item). Until then the row's `map` is only data: an arrival
  adds nothing to the set.
- **Its region.** It is made when a region holding the item arrives: a layer in view must read the source (a genes track
  over the genome). Naming an item does not ask for its region, and a ref naming an item never seen waits, noted
  (`svs.ref.pending`) where a layer's data waits for it.
- The first row of an item carrying a `map` gives it. Its length is its `from_range`'s end, else its CDS (a multiple of
  3, each range `[start, stop)`) over 3 (`codon`), or its segments' furthest end (`segments`, each segment's length
  given). A map that gives none is `svs.stream.item_map` (a warning, once per source), and the item has none; any other
  fault is its mapping's, checked as a dynamic node's. Its mapping has no statuses: its side has no letters.
- **Lifetime** is the dynamic nodes' (above): it comes and goes with the regions that carry it, and a returning item
  keeps its coordinate system's id. While shown or selected it stays, its mapping with it.
- A protein mapped onto a transcript's item coordinate system shows its data along the transcript's exons in a genome
  view; data on the item's own block (the genome's, the protein's) is placed through the paths as on any block. Its own
  letters are not drawn yet.
- They are the viewer's, never a provider's (a patch's or a manifest's node of their kind is `svs.node.kind`), and a
  template over the source's sequences, units or mappings never lists them.
- Data streamed against item coordinate systems (a source's own `coordinate_system` naming one) is not read yet.

## Nested sources

A nested source streams too and states its `data_kind`. It is requested only when something needs its data: a template
that selects its nodes, a document ref under it, or a held node that reads it. Nesting stops at the document source's
`max_depth`; a deeper source is `svs.stream.depth` (a warning) and is not requested.

## Streamed alignments

Rows come in leaf order, so a tree over a streamed alignment needs `leaves: "order"`. A row is requested whole, at load,
only when a mapping endpoint or a positional annotation targets its residues, and so is the reference row. Row labels
arrive as rows are drawn. `column_stats` over a streamed alignment is the source's own `values` source (`<ref>#stats`),
with the default params only; a provider that serves none makes the statistic unavailable. Ordering by key tables,
`consensus` and `profile` are unavailable.

## The viewer's side

What a host's provider can rely on from the viewer (the numbers are config items,
[reference/viewer.md](reference/viewer.md)):

- Requests are coalesced to one pass a frame and run headless. For values and alignments the coarsest level first where
  nothing is held; then the view's level from the center out, then a margin of view widths on each side, then prefetch
  hints (the side the view last panned toward, the next zoom level, the snapshots next to the one shown).
- A miniature (a navigator's map, a `seqv:map` layer) draws each block whole: where it draws a streamed layer, the
  viewer also asks for the level whose bin is about a map pixel over each block, at most a handful of regions a block
  (else coarser), after the view's regions and before the margin; `settled()` waits for them.
- At most `Stream.MaxInFlight` requests at once; a region held or in flight is never asked twice; what an older view
  asked for and still lacks is aborted once no view has needed it for `Stream.CancelGraceMs`, or at once when its slot
  is wanted (the longest unneeded first). A tile computed in place that no view now needs pauses, its rows kept, and
  gives way to the tile the view needs.
- Failures are retried with backoff, then cool down; each outage is reported once. An answer of another shape than asked
  is a failure too. An outage ends when a region of its layer loads: the failed regions no view needs then are dropped
  (asked afresh when a view needs them again).
- The cache holds what drawing needs, bounded by `Stream.CacheBytes`, least recently needed first; an alignment's row
  labels are held and evicted by their chunks as regions are.
- Hover, selection and queries read loaded data only; a selection over a pending region lights it when it arrives.
- `settled()` resolves once every region in view is loaded, empty or failed (for export and tests).

## The tiles protocol

`seqstar-tiles/1` is the one `protocol` value the format knows (a source with `protocol` and `uri`); the model and the
viewer never name it. `@molstar/seq-core/stream/tiles` implements it: a client, the stream provider (`tilesStream`), and
an in-process reference server (`/server`, `tileServer`). It serves an alignment and `values` and `items` annotations
over its columns.

| Resource                               | Answer                                                                                                           |
| -------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `tiles/<ly>/<lx>/<ty>/<tx>`            | 256 × 256 samples at level `[ly, lx]` (each covering 2^ly rows × 2^lx columns): two byte planes, class and share |
| `labels/<a>/<b>`                       | the labels of rows `[a, b)`, a JSON list                                                                         |
| `rows/<k>`                             | row `k` whole: `{label, row}`                                                                                    |
| `stats/<measure>/<lx>/<tx>`            | optional: a column statistic as 1-D tiles                                                                        |
| `values/<lx>/<tx>`, `counts/<lx>/<tx>` | an annotation's 1-D tiles                                                                                        |

- **Class** is a byte: an upper-case residue code, `-` a gap, `0` no data; above level `[0, 0]` the most frequent
  residue among the sample's non-gap cells (ties to the lower code), `-` only where every cell is a gap. **Share**
  (0–255) is that class's count over all the sample's cells, gaps included, × 255. A server may estimate both from a
  subset of a large sample's cells (the reference server reads about 1,024 of them per sample, strided over rows and
  columns), so two servers may differ in the shares of one alignment.
- A 1-D tile has 256 samples of 2^lx positions: `values` and `stats` a float64 `[min, max]` pair each, `counts` a uint32
  each. Everything is little-endian.
- The provider builds its manifest from its readers' declared shape (row and column counts, length, target). 404 or 410
  on any resource is `StreamNotServed`, never retried (for a `stats` measure: the statistic is not served); any other
  failure is retried.

## Archives

An **`.svsx`** is a zip that carries a document and the files it reads, so a viewer, Nightingale or a headless
materializer opens it with no network (MolViewSpec's `.mvsx`, for this format).

- The main file is the document at the root, `index.svsm` (MessagePack, looked for first) or `index.svsj`. A reader may
  be given another main file's path. There is no manifest; the other files are bytes as their formats read them.
- In a document read from an archive, every `uri` may be a relative path, resolved against the main file's directory. A
  path that climbs above the root, a leading `/` or `\`, and any scheme but `https:` are refused (`svs.source.uri`).
  Absolute `https:` URIs stay valid.
- **Bundling** turns a document's absolute `https:` `uri`s into archive files: each is fetched once, stored under a
  stable relative path and the URI rewritten, in the tree and in nodes snapshots insert. Streamed (`protocol`) and
  `service` sources are not bundled.
- An archive is checked before it is unpacked: at most 4,096 entries and 1 GiB inflated by default (a host may set other
  bounds), each entry its declared size and CRC-32. A truncated, damaged or oversized archive is `svs.format.invalid`,
  saying which.
- `readSvsx` (in `@molstar/seq-core/format`) reads an archive; `archiveLoaders` answers relative URIs from it;
  `writeSvsx` and the bundle step are the builder's (`@molstar/seq-core/builder/svsx`).
- A viewer given an archive's loaders (`archiveLoaders`, which carry its base) opens the archive's document: the
  document, its patches and its snapshots are validated against that base, as loading reads them.
