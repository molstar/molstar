# The data model

What a sequence view is made of at runtime, as `@molstar/seq-core/model` holds it: molecules, units and complexes,
numberings, coordinate systems and the maps between them, alignments, annotations, locations. [format.md](format.md)
says how a document names them. The model is platform-free and has no dependencies (a MessagePack codec is vendored in);
its query runtime, `@molstar/seq-core/model/query`, uses `@molstar/seq-core/query` and
`@molstar/seq-core/algorithms/search`.

Its data is typed arrays, struct-of-arrays where things number in the thousands; things are referred to by integer ids,
`-1` meaning none; strings are interned; lookups take out-parameters or callbacks and allocate nothing. Objects are made
at the edges only (parsing, a loci a host keeps, a label).

## Indices and labels

- An **index** is a 0-based position. A range is `start`, `stop`, half-open. The model never says `end`.
- A **boundary** is a place between positions, an index `k` in `0..length`: the boundary before `k`.
- A **label** names a position in a numbering: `315`, `"27A"`, `"111.1"`, `"N13"`; an integer's suffix is letters (an
  insertion code, one or more: `"100AA"` past `"100Z"`) or `.n`. Label ranges (`beg`, `end`) are inclusive and are
  converted to indices once, at load.
- Numbers from data are kept losslessly: doubles as `Float64Array`, integers in integer arrays.

## Molecules

A **molecule** is one coordinate system that units instantiate, of two kinds:

- a **sequence**: a polymer's full residue sequence, unmodeled residues included. Its alphabet is protein, DNA, RNA or
  generic, declared or detected. Residues are ASCII, stored upper case; a character outside ASCII is an issue, never
  truncated.
- a **ligand**: indexed by its heavy atoms, each with its element. An ion is a ligand of one atom.

What differs between the kinds (letters, alphabets) is in `molecule.ts`. A ligand given to an aligner or to `realign` is
`svs.ligand.unsupported`; in a codon mapping, `model.mapping.codon`.

- **Detection**, when the source does not say: nucleotide only when every letter is an IUPAC nucleotide code and at
  least 90% are in `ACGTUN`; RNA when `U` is present and `T` absent; else protein. A sequence that looked like
  nucleotides but is read as protein is a warning.
- **A lazy molecule** has a length, an alphabet and numberings and no letters: a streamed source serves them by region.
  Nothing in the model reads its letters; providers over it are unavailable.

## Units and complexes

A **unit** is one instance of a molecule: a chain, a symmetry copy, a ligand instance, a pose. It uses its molecule's
indices over its **range** and never renumbers them. Its `observed` set marks what was modeled. Every `sequence` and
`ligand` node also has its own unit, `<ref>/unit`, over its full range, so a molecule on its own is shown and addressed
as any unit. Position and boundary loci name units; row, cell and column loci name an alignment; item loci an
annotation.

A **complex** is units shown together, in an order: a structure's chains and ligands, a set of poses. A unit belongs to
at most one. A complex is also a **composite coordinate system**: a coordinate in it is a unit and an index; its extents
are runs per member unit. It adds identity and grouping, never arithmetic: maps and conversions go through its units.

## Copies

A **copy group** is the live units with equal (molecule, range start, range stop). Copies are derived, never declared:
operator copies, same-entity chains and two units over one sequence alike. A molecule's own unit is always a group of
its own.

- **Copies never reach each other.** No path steps up from a unit to its molecule and back down to a sibling unit, so
  chain A's data never lands on chain B. Data meant for every copy targets the molecule, which reaches every unit over
  it.
- A view may merge a copy group into one block; a hover or selection on a merged block stands for every member.

## Identity and lineage

A **set** (`SeqSet`) is the model of one document or one data input. Building a set starts a **lineage**
(`SeqSets.build`); changing it (`SeqSets.withChanges`: upserts, removals, evictions by key) keeps it.

- Within a lineage, ids are append-only. An upsert of a held key keeps its id; a removed id is never given to another
  key; a key removed and upserted again gets its old id back. So tables holding ids are shared across versions, and
  locations stay valid.
- Across lineages, loci are remapped by key; what is gone is dropped.
- A change rebuilds or removes what depends on it. Removals are warnings (`model.unit.orphaned`, `.unit.stale`,
  `.mapping.orphaned`, `.annotation.orphaned`), except an own unit's, which is quiet; a complex drops a removed unit
  quietly. A full-range unit with no numberings of its own follows its molecule's new length.
- **Evicted** keys (a stream's dynamic nodes) go with their dependents, quietly; they keep their ids, so a location
  naming one names nothing meanwhile and resolves again when it returns.
- A version that changes only annotations or mappings costs its changes, not the set's size.

## Numberings

A **numbering** labels the positions of a molecule or a unit, by one of four methods: `offset` (label = index + first,
optionally skipping zero), `segments` (piecewise linear, holes allowed), `labels` (integers with sparse suffixes: Kabat,
IMGT, auth numbers) and `names` (a string per position: atom names).

- **Names**: `ordinal` (built in, first 1), `row` (an alignment row's, from its header) and `uniprot:<accession>` are
  known; any other name is opaque: a provider's (`auth`, `label`, `atom`) or an author's.
- A molecule's numberings hold for every unit over it; a unit's hold for it alone and win over its molecule's of the
  same name.
- **Order** is index order. A label absent from a numbering is placed by comparing its integer part, then its suffix:
  none first, then letters (shorter first, then alphabetically, case aside), then `.n` by number. A label range covers
  every position from `beg` to `end` in that order.
- A numbering's length equals its molecule's; a `segments` numbering's label ranges are disjoint.

## Coordinate systems

A **coordinate system** is a system whose elements are numbered by indices; its **domain** is the extent its indices may
take. Every index names its coordinate system, in its type or its name, and a `Coordinate` is a system and an index.

| Kind       | Elements                                                                                         | Domain                               |
| ---------- | ------------------------------------------------------------------------------------------------ | ------------------------------------ |
| `molecule` | a sequence's residues, a ligand's heavy atoms                                                    | `[0, length)`                        |
| `unit`     | its molecule's indices                                                                           | its range                            |
| `row`      | an alignment row's own residues, inserts included                                                | `[0, row length)`                    |
| `columns`  | an alignment's columns                                                                           | `[0, column count)`                  |
| `complex`  | (unit, index) pairs, composite                                                                   | its units' domains                   |
| `item`     | an annotation item as a whole: what a hit or a hover names when it names an item, not a position | none: its extent is its rows' ranges |
| `external` | a host's (a Mol\* chain's residues), registered through `translate`                              | as registered                        |

A streamed item's own positions (a transcript's residues) are an **item coordinate system**: a molecule with a length
and no letters and its own unit, keyed `<source ref>/<item id>`, made from the item's `map` and mapped onto its source's
system (`seqspec:sources#item-coordinate-systems`). It is named, placed on, hovered, selected and projected as any unit
is.

A coordinate system's id means something only within its set's lineage. In the viewer, each section also has a **layout
coordinate system**, its shown positions block by block, which the model does not have and loci never live in.

## Coordinate maps

A **coordinate map** relates the elements of one system (`from`) to the positions of another (`to`). Each `from` element
covers `span` `to` positions: 1, or 3 through a codon step. Its **support** is what it is defined over on each side, so
a map may be partial.

- **One shape.** Every map is a segment table, whatever made it: identity and a unit's range are one segment; a row onto
  its columns, `segments`, `pairwise` and `by_numbering` mappings one per ungapped run; a codon map one per CDS range; a
  path its steps composed. Maps are built from their method's input, never from per-position arrays.
- **One-to-many backward.** `from` ranges never overlap; `to` ranges may: a source used twice in a construct (tandem
  GFP-GFP) is one map whose segments repeat a `to` range, and a residue of the source reaches both copies. That is not
  ambiguous.
- **Direction.** A segment runs its `to` positions upward or, on a minus strand, downward.
- **Phase.** An element's place within its span (which nucleotide of its codon) is its phase.
- **Boundaries** map to the boundary before their element's first `to` position in reading order; with that element
  unmapped, or at the domain's end, to the boundary after element `k − 1`'s image, and nowhere when that one is unmapped
  too.
- A **mapping** is a coordinate map a document states between two systems (a molecule, unit or row on either side); each
  `from` element has a **status**: exact (equal letters, or a codon that translates to its residue), conflict, or none.
  Mappings are undirected for lookups.
- **Partial maps and chimeras.** A system may map to several others, each partially: a fusion writes one mapping per
  source, with several segments where its part is discontiguous. What no map covers maps nowhere, and that is never an
  issue. A mapping that pairs nothing at all is `model.mapping.empty`, a warning.
- A row and its alignment's columns are joined by an implicit edge (the row's gapped map), not a mapping.

## Paths

Data reaches a coordinate system through a **path**: steps of mappings (read either way), row edges and unit ↔ molecule
steps (the identity within the unit's range).

- Shortest paths only, at most 4 counted steps; unit ↔ molecule steps are free.
- At most one codon step in a path.
- The copy rule (Copies) is part of the search: a longer legal path wins over a shorter illegal one.
- A path that pairs nothing is no candidate. Several shortest candidates whose composed maps differ are **ambiguous**
  (`svs.mapping.ambiguous`), and the placement through them is skipped; candidates that compose to the same map are one.
  At most 16 candidates are compared.
- A path ending at a ranged unit is clipped to that unit's range.
- Paths are composed when a view is resolved, never on hover.
- **Sparse.** A composed path holds only the stretches its maps pair, never its systems' whole domains: a protein's path
  through a chromosome is the size of its CDS. What lies outside them maps nowhere. Read either way by ranges, a path
  takes a stretch of its `to` system back to the runs of its `from` system that reach it (a protein's residues to their
  exons), and the other way.

## Codon maps

A **codon map** pairs each residue of a protein (`from`) with the three nucleotides of its codon on the DNA or RNA that
codes it (`to`): span 3.

- **The CDS** is a list of half-open ranges on the nucleotides in transcript order (exons, spliced), or one range. On
  the minus strand the first range is the highest, each is read downward, complemented. Codons are read across the
  ranges, so an exon junction may split one.
- The codons code the residues of a `from` range (by default from residue 0). Codons past the last residue (a stop
  codon) pair with nothing; residues outside the codons have status none. The CDS's length is a multiple of 3.
- **Status** per residue: exact when the codon translates to it under the genetic code (NCBI's tables, by id, 1 the
  standard), conflict otherwise. `T` and `U` read alike; an ambiguity code translates to `X`. The first coded residue is
  also exact as `M` from a start codon; `U` is exact from `TGA` and `O` from `TAG`.
- A CDS cannot start or end mid-codon: a partial CDS is written from its first whole codon.
- A nucleotide converts to its residue (one); a residue to its three nucleotides.
- Numbering labels are never carried across a codon step: a residue's label does not name a nucleotide.

## Cyclic systems

A molecule may be **cyclic** (a plasmid, a circular genome): its domain is `[0, L)` and its indices are taken modulo
`L`. Its units and paths through it are cyclic too; rows and columns never are.

- **An extent never wraps.** A range across the origin is two runs, `[s, L)` and `[0, t)`. In data, a range written with
  `start > stop` wraps on a cyclic target in side a of `ranges`, `values` and `matrix`, and in a block's range; a pair's
  ends never wrap.
- A map segment that would cross the origin is split there when the map is built.
- Focus and scroll-to take the shortest arc holding an extent (`cyclicHull`), so a feature across the origin is one arc,
  not the whole system.
- A unit whose range wraps is not supported yet (`model.unit.range`); a block's range shows the wrap.

## Pieces and layout maps

A block shows a **piece** of its coordinate system: all of it or a range, laid out in a direction (1 or -1), on a cyclic
system from an origin. The piece's **layout map** takes it into the section's layout coordinate system; a section's
position filter is part of it. A collapsed intron, a reversed strand and a rotated plasmid are layout maps, never new
data. Two pieces of one unit are two blocks over one coordinate system: data, hover and selection reach both.

## Alignments and trees

An **alignment** is coordinate systems: its columns, and per row the row's own residues, with a gapped map from the row
onto the columns. Gapped strings are input syntax only.

- `-` and `.` are gaps. **Insert columns** stay columns and are flagged: from `#=GC RF`, or, without it in aligned
  input, columns whose residues are all lower case (input lower case throughout is read as letters,
  `model.alignment.inserts`). A3M's lower-case inserts map to no column.
- A row header `NAME/27-137` gives the row the numbering `row`, first 27.
- A **lazy alignment** knows its row and column counts and holds only the rows something needs whole (a mapping
  endpoint, an annotation on its residues, the reference row); everything else streams. Every operation on a row not
  held answers none; nothing throws.
- A row is bound to the units or molecules its own mappings reach (its **bindings**); anything further is a composed
  path.
- A **tree** has every row as exactly one leaf. A leaf matches the row with an equal label, else the row whose label
  before its last `/` is equal.

## Annotations

An **annotation** is a table of items over targets, with a **geometry**, a semantic type and provenance:

| Geometry | A row is                                                                                          |
| -------- | ------------------------------------------------------------------------------------------------- |
| `ranges` | a range (or a boundary) on its target                                                             |
| `values` | a value at a position or over a range; or dense, one value per position, no rows                  |
| `pairs`  | two ends, each a range on its own target                                                          |
| `matrix` | a value per (state, position)                                                                     |
| `grid`   | a value per pair of positions: a range on its target by a range on its `b_target`, values a-major |
| `rows`   | a record per alignment row, matrix state or dynamic node: a key table                             |

- Rows sharing an `id` are one **item**; an item's rows are read through an index, never by scanning. On a cyclic
  target, a row with `start > stop` is two rows of one item.
- **Column roles** are one table (`ColumnRoles`): a target or positional column the geometry does not list is an error
  (`model.columns.unknown`) and is left out; descriptive ones (`value`, `id`, `category`, `label`, …) are fields where
  it does not list them. Every other column is a **field**: numbers, strings or JSON.
- Where rows of one layer overlap and one value shows per position, the last covering row wins.
- A **pair**'s meaning is its `category`, an open string. The kernel computes no geometry.
- A **grid** row's b side is placed as a pair's b end. A side without a range (either) is its target's whole domain. Its
  values are one per pair of positions, a-major (row i of the a range holds the b range's values).
- `NaN` is a missing value, never zero.

## Locations and loci

A **location** points to one element; **loci** are many locations of one kind, stored compactly, and are what is kept
and passed (a selection, a highlight, an event).

| Element                    | Location           | Lives in                                              |
| -------------------------- | ------------------ | ----------------------------------------------------- |
| a position of a unit       | `SeqLocation`      | a unit system                                         |
| a boundary of a unit       | `BoundaryLocation` | a unit system, between positions                      |
| an alignment row's residue | `RowLocation`      | a row system                                          |
| a cell (a row at a column) | `AlignedLocation`  | the columns, qualified by a row; a gap has no residue |
| a column                   | `ColumnLocation`   | the columns                                           |
| an annotation item         | `ItemLocation`     | none: its extent is its rows' ranges                  |

- A band over many rows stays one rectangle.
- Conversions follow the maps: gaps drop going to residues, inserts going to columns, each counted; a molecule-level
  binding reaches every unit over the molecule; several bindings mark a result ambiguous.
- A ligand shown as one position is the locus of all its atoms.

## Extents

An **extent** is a set of disjoint half-open runs in one coordinate system. A domain, a map's support, a streamed region
and what a location covers elsewhere are extents. `Extents` is a reusable list of (system, start, stop) entries, in any
order and possibly overlapping, that stands for any loci kind without allocating (`Extents.addLoci`, `Extents.addItem`).

## Projection

**Projection** is the image of an extent through coordinate maps: a plain union of extents. It is pure: it depends on
the source and the routes, never on layers, the layout or another channel. For each entry, the first kind of route that
applies is taken, every route of that kind:

1. **Direct**: the block shows the entry's system (or is a merged copy group with it as a member); an entry on a
   molecule matches any block of a unit over it.
2. **Own molecule**: when the section shows the entry's molecule, only that molecule's blocks light.
3. **Links**: the view's links between sections whose other end shows the entry's system.
4. **Paths**: every shortest legal path (Paths), ambiguous ones skipped.

Through a map a run lands as the union of its elements' images: a residue on its codon, in pieces across an intron; a
source residue used twice on both copies. A row goes through its own row only: its binding's route; an unbound row, or
one bound more than once, lights nothing elsewhere. A column goes through every bound row. An item projects its rows.

- **A boundary** lights a **caret**, a place between two positions, where the consumer draws carets (the viewer's
  channels); elsewhere (a focus, a host's projection) the position after it.
- **An insert**, a stretch its map leaves unpaired between paired positions (a design's grafted loop, a wild-type
  stretch a design deletes, an intron), lights the caret at the boundary before the next pair on the other side, and so
  does a boundary inside one; a boundary that falls inside a codon lands after its residue. Data beyond what a map pairs
  (a fusion partner, a terminal extension) lights nothing there.

## Closure

**Closure** adds, to an element, the elements of other coordinate systems that contain it through a span map: a
nucleotide gains the residue whose codon holds it; a residue gains another protein's residue over the same three
nucleotides (an isoform in the same phase).

- One step: what closure adds is not closed again.
- Through a span-1 map it adds nothing: a cell never lifts to its column.
- It reads the span maps of the whole view, not of the section the element came from, so one element lights the same
  everywhere.
- The result stays in its elements' own systems.

Whether a channel closes is its own option, `closure: "element" | "containing"` (`seqspec:interaction#channels`).

## External systems

A host may own coordinate systems (Mol\*'s chain residues). It registers each with a descriptor and its maps to the
set's systems (`@molstar/seq-core/translate`, `seqspec:libraries#translate`); it joins the set's systems only through
those maps. A key keeps its id across unregistering and registering again, so coordinates held by key resolve again when
a structure is loaded again.

## Providers

Providers compute data from data, at load, cached by their inputs' content and their normalized options:

- **complex providers** read a source into molecules, units and numberings (`molstar:structure`);
- **annotation providers** yield an annotation of a kernel geometry (`coverage`, `column_stats`, `profile`, `paired`,
  `row_stats`, `redundancy`, `query`, an extension's);
- **node providers** yield a tree, an alignment or a mapping (`nj`, `upgma`, `realign`, `codon_rows`);
- **field providers** yield a field by name, usable wherever a field is named.

A provider's result shares nothing mutable with its inputs. Over a lazy molecule or a streamed alignment most providers
are unavailable (`svs.provider.lazy`); past a size bound, `svs.provider.size`; on a generic alphabet, statistics are
unavailable (`svs.provider.alphabet`). Translation is not a provider: a protein is a molecule, written with its codon
mapping.

The kernel's, where the result's meaning is not obvious from its name:

| Provider                  | Gives                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `coverage`                | ranges by category (`observed`, `unobserved`, `outside`, `conflict`); with a `reference`, on the reference, one item per unit                                                                                                                                                                                                                                                                                                                        |
| `column_stats`            | values per column over the counted symbols (an alphabet's standard letters, plus the gap with `gaps: "count"`, the default; others left out): `conservation`, 1 − H(p)/log2(k) over those k symbols (normalized Shannon information, not AMAS); `identity`, the share of counted symbols equal to the reference row's letter; `gap_fraction`; `consensus`, ranges labeled with the most frequent letter (ties to the lower code), its share as value |
| `profile`                 | a matrix of frequencies per column over the counted symbols, the alphabet's standard letters as states                                                                                                                                                                                                                                                                                                                                               |
| `paired`                  | with partners, the ends paired to them, marked on the counterpart units; without, every endpoint; `categories` filters                                                                                                                                                                                                                                                                                                                               |
| `row_stats`, `redundancy` | per-row identity, coverage, length and item counts; a greedy non-redundant subset by identity                                                                                                                                                                                                                                                                                                                                                        |
| `query`                   | a SeqQL query's groups as ranges (`seqspec:seqql#in-documents`)                                                                                                                                                                                                                                                                                                                                                                                      |
| `nj`, `upgma`             | a tree over an alignment's rows from p-distances, deterministic (ties by row order); `nj` midpoint-rooted                                                                                                                                                                                                                                                                                                                                            |
| `realign`                 | an alignment of the inputs by a registered aligner, each row bound to its input                                                                                                                                                                                                                                                                                                                                                                      |
| `codon_rows`              | a codon mapping from each protein row to its cDNA row                                                                                                                                                                                                                                                                                                                                                                                                |

## Custom data

Every model object has `custom`: its node's `custom` and runtime properties a provider attaches (not in the document or
its identity). An element's data is an annotation field. A UI reaches custom data through locations.

## Issues

Every package reports one shape: `{code, severity, path, message}`. Codes are stable, severities `error`, `warning` or
`info` (`seqspec:format#validation`), paths JSON pointers into the input, messages say what is wrong and what would be
valid. Model codes (`model.*`) are reported in a document as the format's (`svs.*`). Nothing is repaired silently.
