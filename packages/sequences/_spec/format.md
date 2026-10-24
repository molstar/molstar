# SeqViewSpec: what a document means

A SeqViewSpec (SVS) document says what a sequence view shows: data defined once, and views of it a reader switches
between. It never says how a viewer draws. This file states the format's meaning and invariants; every node's
parameters, the layer type and theme catalogs, letter maps, color schemes and issue codes are in
[reference/format.md](reference/format.md), generated from the tables in `core/src/format/tables/`. The model's terms
are [model.md](model.md)'s.

## The document

```json
{ "format": "seq_view_spec", "version": "1.0.0", "revision": "r12", "metadata": { "title": "p53" },
  "requires": [], "root": { "kind": "root", "children": [ … ] }, "snapshots": [ … ] }
```

- `version` is the format version the document was written for (Versioning); `revision` is the owner's name for this
  version of the document, compared for equality only, and named by patches as their base.
- `metadata` (title, description, authors, timestamp) is never read for meaning.
- `requires` lists extension namespaces needed to interpret the document (Extensions).
- A document nests at most 256 objects and lists deep, in either file form; a deeper one is not read
  (`svs.format.invalid`, at the path where it goes too deep). A dynamic node a patch adds is held to the same depth in
  its params (`svs.params.type`: the node is not read). A stored query's expression nests at most 64 deep, as its text
  may (`svs.query.invalid`, `seqql.read.depth`).
- Files are `.svsj` (JSON) or `.svsm` (MessagePack); their media types are `application/vnd.seqstar.svs+json` and
  `application/vnd.seqstar.svs+msgpack`. An archive with its data is `.svsx` (`seqspec:sources#archives`).

**Against MolViewSpec.** SVS keeps MVS's node shape, palettes, color names, `*_from_uri` and `*_from_source`, and
description format. It differs where identity and sequences need it: every node has a `ref`; `custom` is keyed by
namespace; snapshots are patches over one tree; positions are indices beside labels; data may be held in place
(`value`); a categorical palette repeats its list by default.

## Values

Every value in a document, and in the data it loads, is a JSON base type. A number means what its JSON text means:
`0.15` is the double nearest 0.15, compared, keyed and listed as that text.

- **Typed columns.** In memory a table's number column may be a typed array of any kind; it means the array of its
  numbers, NaN as missing. A `Float32Array`'s numbers mean their shortest decimals, so a column colors, sorts and lists
  as its JSON text does. Elsewhere lists are JSON arrays.
- **MessagePack** carries the same data model: integers within 32 bits as integers, other numbers as float64, non-finite
  numbers as nil. A typed column may be written as an extension type (`0x11`–`0x18`: int8, uint8, int16, uint16, int32,
  uint32, float32, float64, little-endian), decoding to exactly those numbers. A reader tells MessagePack from JSON by
  the first byte.

## Nodes and refs

Every node is `{kind, ref, params?, custom?, children?}`. The grammar (which kind under which, how many) is the tree
table's (reference: Nodes).

- **Every node has a `ref`**, unique in the document, `^[A-Za-z0-9][A-Za-z0-9_.-]*$`. `root` and `data` may leave theirs
  out. Refs are what patches address and what identity is.
- **References** between nodes are params typed `ref<kind>`. Nothing cascades from a parent: every default is stated in
  its table. A node's refs are exactly the ones its typed params name (a table's target column, a stored query's names,
  a theme's `column`), and they alone make its dependencies: what reruns when a node changes.
- **Derived refs.** A node that yields several things names each `<ref>/<suffix>` (`1brs/A`, `msa/3`, `p53/unit`, a
  streamed source's `hg38/chr17`). Author refs hold no `/`. Static validation checks a derived ref's syntax and that its
  prefix can yield it (a sequence or ligand only `/unit`, an alignment numeric rows, a provider's complex or a streamed
  source any suffix; else `svs.ref.wrong_kind`); whether it names something is known after loading
  (`svs.ref.derived_unknown`).
- A row is `R/<k>`, `k` its 0-based index; mappings, annotation targets and focus take it.
- A streamed item's own positions are `<source ref>/<item id>` (`genes_src/T1`, the item's `id` as its provider sends
  it), its unit `<source ref>/<item id>/unit`: an item coordinate system (`seqspec:sources#item-coordinate-systems`).
- Refs that form a cycle are `svs.ref.cycle`.

## Positions

Two forms, never mixed in one table:

- **Indices**, canonical: `start`, `stop` half-open, one position `index`, in the coordinate system the target names (a
  unit's indices are its molecule's: a unit over `[241, 493)` starts at 241).
- **Labels** in a named `numbering`: `beg`, `end` inclusive, one position `pos`; converted at load by the numbering's
  order (`seqspec:model#numberings`): a range's end the numbering lacks lands on the nearest position inside the range.
  A position label it lacks, or a range that covers none of its positions, leaves its row out (`svs.position.unknown`, a
  warning, naming the label and the labels the numbering has); a row with no position at all is `svs.params.missing`.
- **Insertion codes.** A label carries its code (`"52A"`), or a table gives it apart: `ins_code` (`b_ins_code` for side
  b) completes the author number in `pos` or `beg` (`b_beg`) of one residue, as MolViewSpec's `pdbx_PDB_ins_code` does
  with `auth_seq_id`: `pos 52, ins_code "A"` is 52A; `""`, `"?"` and `"."` are no code. A code is letters (else
  `svs.position.syntax`). A code on a range (`end` other than `beg`), on a label that has one, or in a table without
  that side's labels is `svs.position.form`, saying what to write: a range's ends with codes are whole labels
  (`beg "52A"`, `end "53"`), and a range by plain numbers holds the inserted residues between them.
- **Wrapping.** On a cyclic system a range with `start > stop` (or `beg` after `end`) runs across the origin; elsewhere
  it is `svs.position.range`, as is `start == stop`.
- Producers holding indices write indices; labels are for people and agents writing by hand.

## Data nodes

Data lives under `data`, defined once. A data node with content holds it in `value` or reads it from a `source`
(`source_ref`), exactly one; a source's `value` is a whole file in place (`seqspec:sources#data-by-reference`). A
numbering and a mapping have no content: their method's params say it.

- **`sequence`, `ligand`**: a molecule (`seqspec:model#molecules`); each has its own unit `<ref>/unit`. A `cyclic`
  sequence's indices wrap. A service may add numberings to what it reads.
- **`numbering`**: under a sequence, ligand or unit (`seqspec:model#numberings`). `ordinal` and `row` are reserved
  names.
- **`unit`**: an instance of a sequence or ligand over a `range`, with `observed` positions.
- **`complex`**: units shown together, from a registered **complex provider** reading a source (what it includes and how
  it names units are its options), by value (`unit` children), or from a streamed `letters` source (a unit per
  coordinate system of its manifest: a genome as a complex of chromosomes). Without its provider a complex is
  unavailable: its derived refs name nothing and nothing else fails.
- **`alignment`** and **`tree`**: by value, source or provider. An alignment is coordinate systems: its columns can be a
  section's axis, its rows map onto them (`seqspec:model#alignments-and-trees`). A computed alignment (`realign`) keeps
  its inputs' labels and order and binds each new row to its input by an identity mapping, so data on the inputs still
  lands.
- **`mapping`** (Mappings), **`annotation`** (Annotation tables), **`computed`** (Computed nodes).

**Load order** decides availability bottom-up: sources, complexes from providers or streams, molecules and alignments,
units, complexes by value, tables, rows and mappings, then providers.

## Mappings

A `mapping` is a coordinate map between two systems (a sequence, ligand, unit or row on either side), by a `method`:
`identity`, `segments`, `pairwise` (a CIGAR), `by_numbering` (positions with equal labels), `codon`, `align` (computed
at load by a registered aligner, stored as `pairwise`), or `provider`.

- Each element's status (exact, conflict, none) is computed at load from both sides' letters.
- **`codon`**: a protein's residues over their codons on DNA or RNA, by `cds` ranges, `strand` and `genetic_code`
  (`seqspec:model#codon-maps`). A feature's own strand field is not read by it.
- **`align`**: aligners are deterministic; the stored provenance names the aligner and its normalized options. The
  kernel's are `biomsa` (default), `needleman_wunsch` and `myers`. A ligand cannot be aligned
  (`svs.ligand.unsupported`). `needleman_wunsch` and `biomsa` work within 2^26 dynamic-programming cells; `myers` has no
  such bound. Past it, or on any failure, the mapping is `svs.mapping.align`, a warning: the sequences are valid and the
  mapping is a missing piece (for `realign`, `svs.provider.size`).
- **Partial and chimeric**: one mapping per source, several segments where needed; a source used twice is one mapping
  repeating a `to` range. Two mappings joining the same two systems make every path through them ambiguous.
- Mappings are undirected for lookups; `from` and `to` only fix how params read.

## Annotation tables

An `annotation` is a table of items over targets (`seqspec:model#annotations`), by `value`, `source_ref`, `uri` (with a
`format`), or `from_source` (the file a sequence, alignment or complex came from). It may sit under the layer that shows
it.

- A target that is a sequence or ligand reaches every unit over it; a unit, that unit; an alignment, its columns; a row,
  that row's residues.
- A **`Table`** is rows of named columns, written as an array of row objects or an object of column arrays (both read
  the same); a column may be dictionary-coded (`{dict, codes}`) or JSON (`{json}`).
- Which columns each geometry reads is the reference's Types: `Table<geometry>`. A target or positional column of
  another geometry is `svs.columns.unknown`; any other column is a field.
- `fields` remaps any column of a table read from a file, positions and targets included (a CIF category by author
  numbers: `{ pos: "auth_seq_id", ins_code: "pdbx_PDB_ins_code" }`); with an inline `value` it is
  `svs.params.exclusive`.
- A target value read from a file is a suffix of the node the file came from, or of the `target` param.
- A table read from a file holds text cells, typed by their column's role (an id `01` stays `01`).
- A **`grid`** row's `value` is a list of (a length) × (b length) numbers, a-major, `null` missing, or a `Float32Array`
  or `Float64Array` of them (NaN missing), kept as given: the host must not change a typed row's values once given (what
  is drawn from it would go stale); it gives a new array, and a new document or patch. A list of another length is
  `svs.columns.type` and leaves the row out, saying how the count is formed. Its sides never wrap a cyclic origin.

## Computed nodes

A `computed` node is an annotation a registered provider computes (`seqspec:model#providers`): always of a kernel
geometry, so every viewer can place it. A provider not registered leaves the node unavailable
(`svs.provider.unavailable`, a warning) and the layers over it unplaced; nothing cascades. A result is checked as a node
of its kind and cached by its inputs' content and normalized options.

## Views

A `view` is one way to show the data, picked by its `label`; `root.default_view` (else the first) opens. Views share the
data, so a switch reloads nothing. Its `layout` is `tracks` (blocks on an axis, tracks across them) or `text` (blocks
read as text, wrapping, no zoom).

## Sections

A `section` is one axis of a view with its tracks; a view's sections are shown one under another.

- **The axis** is blocks of units (side by side or stacked) or one alignment's columns, never both.
- **Linked sections**: hovering or selecting in one lights the counterparts in the others, through projection. A shared
  `scroll_group` also shares zoom and horizontal scroll; its sections' axes must be equal (the same units or members,
  ranges, directions, origins and shown positions).
- **Following** (`follow`): two sections set each other's window through a coordinate map, whatever their axes: the
  other centers the image of the leader's center element and takes the leader's cell times the route's scale (1, or 3 or
  1/3 through a codon step). Each keeps its vertical scroll. The anchor is the block under the center, or on an
  alignment axis `follow_anchor`, else the reference row. Cycles, other views and unjoined sections are
  `svs.section.follow`.
- **Numbering** (`numbering`): what rulers and labels use. Unset: per block the numbering marked default on its unit or
  molecule, else `ordinal`; on an alignment axis the column ordinal. A named numbering is looked up on the block's unit,
  its molecule, then along mapping paths (never across a codon step).
- **`positions`**: show only what these annotations cover (a filter, a layout map). Data unavailable: the filter is
  ignored.
- **`highlights`**: persistent highlights, each entry an annotation's items (a query node's included) in a style
  (`band`, `highlight`, `accent`), lit in that section only.
- **`style`**: a preset name (the viewer's) or an inline preset whose keys override that preset's sizes and switches for
  this section. The document's keys win over the preset's, which win over the host's; a host may still bound pixel
  values for its embedding. An unknown key is a warning and is ignored.

## Blocks

A `block` names one unit (or a sequence or ligand, through its own unit), a complex (a block per unit, in its order) or
an alignment. Without blocks a section takes every polymer unit in data order, then ligand units; with none, every
molecule's own unit; with none, the first alignment.

- **Pieces**: `range`, `direction` and `origin` make the block's piece (`seqspec:model#pieces-and-layout-maps`). A range
  makes no unit. A reverse complement is `direction: -1` with the `complement` letter map.
- **Size**: a block's piece lays out at most 2^31 - 1 positions (the axis holds them as Int32; a section's blocks
  together may hold more, a whole genome); a block past that is left out with `svs.position.range`, saying to show a
  `range` of it.
- **Copies** (`section.copies`): `one` merges the blocks of one copy group with equal range, origin and direction into
  one, labeled `<molecule> ×N (<units>)` unless the block has a label; `each` shows each unit. Copies never reach each
  other.
- **Ligand form**: a ligand block as one position labeled with its component, standing for its atoms, or one position
  per heavy atom. A viewer may flow one-position ligands after the chains in a `text` layout (the Seq\* viewer's chips:
  viewer.md, The text view).
- Blocks of an unavailable node are left out quietly.

## Tracks and layers

A `track` is one row of a section, in a `group` or not; groups expand and collapse at runtime from their `collapsed`
state at open. A `layer` is what a track shows: a type, its data and its channels' themes. A track's layers are
composited in order, later over earlier.

- A track `height` is a hint in rows. Under `row_sizing: "uniform"`, layers of `row` sizing share one row height. In a
  `text` layout whose style asks for it (`wrapped_tracks: "full"`; the default `"lanes"` keeps thin lanes) a track with
  a height is drawn in full under every line, that many of the text's rows (viewer.md, The text view); the style's
  `track_axes` and `track_labels` add its axis and its label there.
- A `summary` track shows under its group's header only while the group is collapsed.
- At most one `rows` layer per section: the first wins (`svs.layer.rows`).

## Placement

A layer is drawn in every block its data reaches. Each data slot is read in its own coordinate system and placed through
the map from it to the block's: directly, along a path, or onto an alignment's columns through rows
(`seqspec:model#projection` has the same route rules).

- **Direct placement wins.** Data on a unit or molecule whose molecule the section shows is placed only on that
  molecule's blocks, never through another path onto another molecule's block.
- **Onto columns, and through any partial map**, a range spans its first to last paired position, what lies unpaired
  between included (a gap; a wild-type domain on a design spans a grafted loop: one range). A range wholly inside an
  **insert**, unpaired between paired positions (a row's insert, a design's grafted loop, a wild-type stretch a design
  deletes), lands on the boundary after it (`svs.layer.insert`, counted).
- **Through a codon step**, a range covers the union of its residues' codons, split where the CDS is spliced;
  per-position types draw one mark per residue across its codon.
- **Through a partial map**, data beyond what it pairs (a fusion partner, a terminal extension) draws nothing; this is
  never reported.
- **Lanes.** In a track with a `rows` layer, other layers land on the lane of the row they target; a `frame` on the
  lanes of the rows track its span reaches.
- A layer without data is placed by its type (`letters` and `ruler` in every block).
- A layer that reaches no block is `svs.layer.unplaced` (a warning), unless its data has no rows or is pending.

## Track templates

A track with `for_each` stands for one track per dynamic node it selects (`seqspec:sources#dynamic-nodes`): by kind,
source, nesting depth, whether the node reaches a block of the section, with an order and filters over key tables keyed
by node. Its layers are written once over `"$node"`.

- The viewer adds an instance as a node arrives and removes it when the node is evicted; instances keep the template's
  place in the layout.
- An instance's ref is `<track ref>/<node ref>`; it is never `svs.layer.unplaced`.
- `"$node"` outside a template is `svs.ref.syntax`.

## Focus

A section's `focus` is the range it opens on: a target with a range (indices or labels), or an annotation's items (their
extent). Fit-to-all is the default. The last focus wins; one outside its target is `svs.focus.range`. On a cyclic system
a focus may wrap, and an annotation's extent is its shortest arc. Hover and selection are never stored.

## Layer types

A type is named by what it shows. Each declares the data it accepts, its params, where it is defined (layouts, axes),
its **sizing** (`row` or `free`), its channels and a **fallback**: the type a viewer draws when it cannot draw this one.
Every chain of fallbacks ends at `letters`, `ranges` or `omit`, and every viewer draws `letters` and `ranges`. Drawing a
fallback or omitting a layer is reported (`svs.render.*`). The catalog, with each type's params and fallback, is the
reference's Layer types.

- **Extension types** (`<namespace>:<name>`) have a registered table and a kernel `fallback`, which the document must
  state; their data must suit the fallback.
- **Split rows** (`ranges`, `points`, `variants`): one labeled row per value of a field (`split.by`), in a given order
  (a fixed grid, empty values shown or hidden). A track opens expanded or collapsed; a viewer that cannot split draws
  the layer collapsed. At most one split layer per track. A viewer may draw a faint guide between expanded rows (the
  Seq* viewer does, where a row is tall enough to show one: viewer.md, Render props).
- **Density** (`ranges`, `points`, `variants`): how a device column holding more items than it can draw shows them,
  `stacked` (each color's share) or `majority`. A layer draws items or density for all its blocks at once.
- **`letters` of other polymers** (`of`): a protein's residues over its coding DNA, placed as data is; `"mapped"` draws
  every protein whose codon mapping reaches the block.
- **`ruler`** may show another coordinate system's numbering (`of`), through the map that places it.
- **`secondary`**: a strand is an arrow toward its C-terminal end as shown (on a reversed piece, to the left), its head
  on the piece holding that end only; a strand split by an intron draws its other pieces as bodies. A helix is a
  cylinder, its ends rounded where the helix ends and square where a line, the view or an intron cuts it (it goes on).
  With `labels` each item's label (`α1`, `β2`) is drawn over its shape where it fits, whole or cut to three graphemes.
- **`rows`** are drawn in a text layout too, a text row per row on every line (viewer.md, Rows in text), except a
  streamed alignment's: text is a figure of rows held, so there they take their fallback (`svs.render.fallback`).
- **`frame`**: an outline per item and run over the cells it covers, across the tracks `span` names (from `from`'s band
  to `to`'s, in the frame's section), titled with the item's label; where a span's edge is a `rows` track, the rows the
  item targets (a row target each; on the alignment's columns, every row). A side is drawn only where the item ends, so
  a frame cut by a line break or the view stays open there. A span naming a track of another section is
  `svs.layer.span`, and the frame is drawn on its own track alone. Its fallback, `ranges`, shows the same items.
- **Dense states** (`matrix`): a band per state drawn, in its keyed order; bands thinner than a device pixel are reduced
  per device row by the layer's `density`: `majority` (the default) takes the plurality of the row's cells' colors (a
  tie: the color met first), opaque by the share of its cells that have a color (cells without one fade it); `maximum`
  takes the cell with the largest value (the first on a tie) in its own color, so a rare signal shows. A hit names the
  row's first state with the drawn color. A band's label shows only where the band is 9 px or taller. A matrix takes
  part in the row zoom: its height follows the row scale.
- **`pairs`** arcs stay within their track's band at every zoom, never lie flat along their ends, and move continuously
  as the view pans: an end beyond the window is taken at most a quarter of the window's width past its edge, and an arc
  is the parabola over its ends so taken, its height proportional to that span against the widest in view (or its
  `height` theme's share of the band), so nested arcs nest and the widest fills the band. A label shows only under the
  apex of an arc whose ends both lie within that reach, the apex in view.
- **`values`** are drawn as bars, a line or an area (`encoding`) in the layer's color theme at each position's value as
  given, not clamped to `domain` (which bounds the height only). A line takes each position's color over that position's
  cell, so it changes color halfway between two positions, and holds its end values over their cells to their edges (at
  a block's ends and around a missing value), as the area's cells reach them. Where a device column holds several
  positions, each column is drawn over its column in the color of its positions' mean. With no color theme a layer is
  drawn in one color. In a `text` layout a `values` layer is a strip under each line, a color per position (viewer.md,
  The text view).
- **`dot_plot`** has two axes over its section's layout coordinate system: x across, as every track's, and y down the
  track, showing the section's window (the positions the data area shows across, scaled to the track's height), so a
  zoom or a pan of the section moves both and a cell (x, y) pairs two positions in view. Each end is placed as a pair's
  end is (side a on x, side b on y, through paths); a ligand block takes no plot. A `pairs` row whose ends are both one
  position is a dot at every zoom; one with a longer end, a segment from (x, y) of its ends' first positions to the far
  corner of their last ones. A `grid` row is its cells, each in the color theme at its value (`domain`, and the numeric
  default color, as a matrix's), and where a device pixel holds several cells it shows their mean. A block pair's grid
  (the cells its rows place within the blocks' pieces) allocates at most 16,777,216 floats (`2^24`, 64 MB): its means at
  coarser zooms (about a third of its cells), plus a copy of its cells unless one `Float32Array` row is exactly the grid
  on unmapped ends, which is drawn from as given. So a single chain's 5,000 × 5,000 PAE given so fits (about 8.3M
  floats), and a grid copied from lists fits up to about 12.6M cells. Past it that pair is left out and reported
  (`svs.render.omitted`, with its cells, the floats and the limit). On a stacked section each lane's plot is its own
  block's, both ways. Its `aspect_ratio` (width ÷ height, 1 by default: square) sizes the track from the data area's
  width, at most what the section shows of its rows less a margin (a viewer's cap), but never flatter than the ratio
  stretched by a viewer's bound (its render props), which wins over that cap: a wide plot is taller than the box and its
  rows scroll. A ratio overrides the row zoom; `null` leaves it the track's own height, which follows it. The y axis is
  the section's own; another section's is not yet a choice.
- Drawing hints (`display`, `labels`, `residue_form`, `side`, `density`, `encoding`) never change meaning.

## Keyed sets

An alignment's rows, a matrix's states and a template's instances are members drawn in an order. A **key table** is a
`rows` annotation keyed by `row`, `state` or `node` that gives members fields (a document's, or computed: `row_stats`,
`redundancy`, a field provider's).

- `order`: a list of keys (`"data"`, `"tree"`, `{data, field, direction}`), sorted lexicographically, ties in data
  order, missing values last. A matrix has no tree: `null` is data order, and `tree` is reported (`svs.layer.data`) and
  ordered as data.
- A `rows` and a `matrix` layer take `order` and `filter` (their `filter` is a keyed set's, not a filter of the data's
  rows); a matrix's key tables are keyed by `state`, a row naming the first state of that name. A filter that keeps no
  state is noted (`svs.layer.empty`); the track is drawn empty.
- `filter`: every filter must keep a member: its field holds one of `values` and lies within `[min, max]`.
- A tree is drawn, and clades fold, only while the order is the tree's.
- A streamed alignment is ordered by tree or data only (`svs.provider.lazy`).

## Letter maps

A `letter_map` (on `letters` and `rows`) changes the letter shown for each stored one: a registered map by name
(`complement`, `complement_rna`, `lower`, an extension's) or an inline dict. A map may list the alphabets it is for;
letters of another alphabet are shown as stored, so a protein is never complemented. A map changes what is drawn and
what the `letter` theme reads, never data: the model, statuses and search keep the stored letters. Reversing never
complements by itself.

## Themes

A **theme** maps the data under a layer to the values of one **channel**: `color`, a size (`width`, `height`, `radius`,
`thickness`), `glyph` or `line_style`. Each channel is a layer param holding a `ThemeRef`, `{theme, params, fallback}`.
"Theme" means data-driven values only: the viewer's look is its skin.

- `params` are checked by the theme's table for the channel's kind; a theme given for a kind it does not serve is
  `svs.theme.channel`. Themes nest (a base theme as a param, a fallback), at most 4 levels counting the first
  (`svs.color.invalid`, `svs.glyph.invalid`, `svs.line_style.invalid`).
- An extension theme states a kernel `fallback` of the same kind.
- Size channels are `px` or `band` (a share of the element's row, lane or ring); sizes change how much of its slot an
  element fills, never the slot, so hits and marks never depend on sizes.
- **Line styles** (`ranges`' `line_style`; `null`: `solid`): how each item is drawn, one of `solid`, `dashed`, `dotted`,
  `wavy`, `double`, `box` and `highlight`; `uniform` with a `line_style` gives the whole layer one, `field_line_style`
  reads a column of the rows (`line_style` by default), as `field_color` reads colors. A style is drawn wherever the
  item is (tracks, the text's lanes, a track drawn in full, an export): solid is the bar or underline; dashed and dotted
  the line dashed, the pattern fixed to the content; wavy a wave of a fixed length; double two thin lines; box an
  outline around the item's positions; highlight a translucent fill behind them (in tracks the whole row; in text also
  behind the line's letters, as box draws there too). A value that is not a line style is drawn solid and said once
  (`svs.render.line_style`). A style never changes the item's slot, its hits or its marks. `points` take none: a glyph
  has no line.
- **Defaults**, when a channel is `null`: a provider's declared default color where every data node comes from one
  provider that declares one, else the type's (`letter` for letters, a categorical palette for features, a continuous
  one over the layer's domain for numbers).
- **One palette per layer.** Order of first occurrence and a normalized domain are taken over all rows of all the
  layer's data, whatever its filter keeps, so a value has one color in every block, section, filter and viewer. Where
  rows overlap, the last covering row wins.
- Legends are the viewer's; a value with no name lists nothing.

## Palettes

Palettes are MolViewSpec's, one implementation with Mol\*'s (`@molstar/core/util/color/value-palette`): the same kinds
(`categorical`, `discrete`, `continuous`), params, named lists and dictionaries, and evaluation, so a value gets one
color in SVS and in MVS. Two differences: a categorical list repeats by default, and a categorical `missing_color` also
colors a missing value. Letter dictionaries (Zappo, Taylor, …) are keyed by one-letter code; on ligand atoms, by element
symbol.

Streamed data has no rows to take keys from: its manifest's vocabulary of the field gives a categorical palette's keys
and their order (`seqspec:sources#manifests`). A value of a streamed item that the keys lack (the manifest lists no
vocabulary of the field, or one without it) takes, from a list palette, the color at the FNV-1a hash (32-bit, over the
value's UTF-8 bytes, case-folded under `case_insensitive`) modulo the list's length, so it has one color in every
region, view, viewer and reload; a dictionary palette leaves it the missing color. An item that sends no value of the
field (a stream with no such column) has the missing color. A legend cannot list hashed values: where the manifest lists
no vocabulary of the field, the layer has one note, `svs.legend.unlisted` (info), saying both.

## Rules

The `rules` theme is thresholds and composition as data, for any channel kind: the first rule whose clauses all hold
gives its value, or a nested theme's. A clause tests a **measure** per element: the layer's value, a field, a column
annotation's value, the share or rank of the cell's letter in its column, or equality with a row's letter. A measure
that does not apply is NaN, and a clause on NaN fails. Named color themes (`clustal_x`, `percent_identity`, …) are names
for `rules` tables, so a viewer implements `rules` alone. Where a measure cannot be computed (a lazy alignment's
columns), the viewer draws the fallback and reports `svs.render.theme`. `espript` is ESPript's shading as such a table:
a residue whose letter is every residue of its column (gaps not counted) in `identical` (black), else one whose letter's
group (`groups`) holds at least `similarity` (0.7) of the column in `similar` (grey), nothing elsewhere; an
approximation of ESPript's similarity score by groups.

## Shared scales

A layer's `domain_group` shares every normalized domain it resolves (a palette's, a size scale's, a `values` layer's y
range) with the view's layers in that group: the union per (group, field). A domain given explicitly wins in its layer,
and given ones form the group's union. A group of one is none.

## Color schemes

A **color scheme** is a named, curated way to color one kind of subject: residue letters (Clustal X, Zappo), numbers
(Viridis) or categories (Okabe–Ito). It is data, registered like any entry (`colorSchemes`), and becomes the `letter`,
`field_palette` or `rules` theme it stands for. A layer names one with the `scheme` theme (`name` or an inline `custom`,
with overrides). A scheme applies only to its subject and alphabets; one that does not fit is drawn as the layer's
default and reported. A scheme has a dark variant, given or derived; its legend lists the variant's colors, never the
other's. Picking a scheme at runtime is the viewer's (`seqspec:interaction#color-schemes-at-runtime`).

## Queries

A `query` computed node holds a SeqQL query (`seqspec:seqql#in-documents`) and yields a `ranges` annotation, one item
per group, so every place that takes an annotation takes it. A document's query is data only: it reads no host input and
does not search. Its issues are `svs.query.*`.

## Extensions

- A **namespace** is a short name or reverse-DNS (`lab`, `org.example`); a namespaced name is `<namespace>:<name>`.
  `svs`, `seqv` and `mvs` are reserved.
- A namespace owns its custom fields, providers, layer types, themes, letter maps, color schemes, commands and issue
  codes (`<namespace>.<code>`). Format and service names are one flat registry each.
- **`custom`** on any node holds extension fields by namespace. An extension adds meaning beside the kernel and never
  changes a kernel field's meaning.
- An extension ships **registry entries**, pure data (formats, services, providers, aligners, layer types, themes,
  letter maps, color schemes, custom tables, issue codes). A reader validates the namespaces it has registered and
  passes the others through; a namespace in `requires` it lacks rejects the document (`svs.requires.unknown`).
- **Missing pieces never fail a document**: an unregistered provider, format, service, type or theme makes its node
  unavailable or its layer fall back, reported once.
- **`mvs`** links units to a MolViewSpec state (`structure_ref`, `label_asym_id`, `auth_asym_id`, `instance_id`). A
  structure provider writes them; nothing in Seq\* reads them.

## What a change affects

A host changes a view by a new document or a patch; the format defines what is unchanged.

- **Identity** is a node's ref.
- **Modified**: a node whose kind, params, `custom` or list of child refs differ. Params compare normalized (defaults
  filled, tables as columns).
- **Affected**: a modified node, or one that refers to an affected node, transitively, through the typed references (ref
  params, `source_ref`, `from_source`, derived refs, `target` columns, a layer's `data` and `blocks`, a section's
  `positions`). A numbering affects its parent and a unit its complex. An ancestor is not affected by containing a
  modified node.
- A change to any data node re-resolves placement; a section keeps its viewport when its axis is equal.

## Patches

A **patch** changes a document by ref: `{format: "seq_view_spec_patch", base?, revision?, ops}`. Ops are `set` (a
shallow merge of `params` and `custom`; `null` removes a key, back to its default), `remove` (with its subtree),
`insert` (a subtree under a parent, after a sibling) and `move`.

- **Anchors**: `after` is a sibling's ref, `"first"` or `"last"`, never a number.
- **Format** is `seq_view_spec_patch`; another `format`, when given, rejects the patch (`svs.format.invalid`).
- **Base** is optional; when given, a document with another revision, or none, rejects the patch (`svs.patch.base`). The
  patch's `revision` becomes the document's; without one, the document has none. Revisions are the committer's.
- **All or nothing.** Shapes are checked first, then ops apply in order; only the result is validated. On any error the
  document is unchanged. An op's own errors are pathed `/ops/<k>`; the result's issues are pathed at the last op that
  touched their node (else `/ops`), naming the node's ref and path.
- **No cascade**: removing a node leaves refs to it dangling (`svs.ref.unresolved`), so the patch fails. `root` and
  `data` are never removed or moved. Ops never address derived refs. `set` never changes a kind.
- Applying a patch shares every untouched subtree and returns the **inverse**. `diff(a, b)` gives the patch from `a` to
  `b`, matching nodes by ref and kind.

## Snapshots

`snapshots` make a document a story: an ordered list, each a titled patch (`ops`) over the base tree, with a
description, a `duration_ms` and a `transition_ms`.

- Every snapshot applies to the base, never to another, so any one opens directly.
- A snapshot's `ref` is unique and is not `next` or `prev`, which are a link's steps (`svs.ref.syntax`).
- **Entering** a snapshot: when its ops set `root.default_view`, that view is shown; when they set none and the snapshot
  left set one that is still shown, the base's view returns (a view the reader picked stays). Every section whose focus
  changed moves to it over the transition; others keep their viewports.
- A ref a snapshot inserts is absent from the base and has one kind and parent in every snapshot.
- A base patch after which a snapshot no longer applies marks it invalid (`svs.snapshot.invalid`, a warning); a player
  skips it.

## Descriptions and invocations

Descriptions (the document's, every snapshot's) are text, `markdown` or `plaintext`. A viewer shows markdown only
through a host's renderer, sanitized; link text has no meaning in the format.

An **invocation** is the one encoding of a request to a viewer, plain JSON for a link, a button, an agent:
`{"command": "seqv:focus", "params": {"target": {"refs": ["1iep/A"], "args": {"beg": "315", "end": "315", "numbering": "auth"}}}}`.
`command` is a registered command's id; `params` are checked by its table when it runs. Refs are document refs, narrowed
by indices, labels or item ids. An invocation is never stored in the tree, and running one changes the document only
through a patch the viewer proposes.

## Interaction is the viewer's

Zoom, scroll, hover, selection, navigators, legends, level of detail and runtime expansion are the viewer's; the
document holds what a view starts with. A UI action that changes what is shown is a patch the viewer proposes to the
document's owner, who commits it or not.

## Builders

A document from any builder means the same: conveniences expand into ordinary nodes, refs are generated, defaults are
left out. Builders are not part of the format.

## Versioning

`version` is semver. Minor versions only add (optional params whose default keeps the old meaning, new kinds, types and
themes whose fallbacks end at an earlier one). A reader of `X.Y` given `X.Z`: `Z ≤ Y` strict; `Z > Y` lenient (each
unknown a warning, `svs.version.newer_minor`; an unknown type or theme drawn by its fallback, an unknown kind left out);
another major rejected. `upgrade` runs one step per older version. No JSON Schema is published: the reference is
generated from the tables.

## Validation

Every issue is `{code, severity, path, message}`: a stable code, a JSON pointer, and a message that says what is wrong
and what would be valid. Nothing is repaired silently. Issues come from five phases, one set of codes:

| Phase   | Checks                                                                                                       |
| ------- | ------------------------------------------------------------------------------------------------------------ |
| static  | the document alone: shape, grammar, params, refs, positions' syntax, inline tables, themes, patches, queries |
| load    | sources and their contents: formats, records, alphabets, numberings, mappings, providers                     |
| resolve | a view against the loaded set: axes, focus, links, follows, placement, split fields                          |
| stream  | each data patch and the dynamic nodes it adds                                                                |
| render  | what a viewer draws instead, or leaves out                                                                   |

An issue that repeats for one cause on one node (down a column, per record of a file, per row or tile fetched) is
reported once, with a count and its first examples; where each example has a path of its own (an annotation's rows, an
alignment's parents), the first ten are reported each at its path, then one issue counts the rest. Every code and its
severity is listed in the reference's Issue codes, and a code has that one severity wherever it is reported. A severity
means one thing: `error` is invalid input (what to write instead), `warning` a loss the person should know about
(something is left out or drawn otherwise, and the rest stands), and `info` a quiet note on how valid input was read (a
ref waiting for streamed data, inserts read from lower case, a file's last terminator missing). A note never changes
what is drawn and never fails a document, patch or commit; a viewer may hide it. Extensions report theirs in the same
shape, under their namespace, registered with their severity.
