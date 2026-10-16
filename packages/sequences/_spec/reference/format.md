# SeqViewSpec 1.0.0 reference

Generated from the parameter tables.

## Nodes

### `root`

the document. Parent: none (the root); ref: implied (`root`).

| Param | Type | Default | |
|---|---|---|---|
| `default_view` | `ref<view> \| null` | `null` | null: the first view |

### `data`

everything views show. Parent: `root`; exactly one under its parent; ref: implied (`data`).

No params.

### `source`

where data comes from, fetched once however many nodes read it. Parent: `data`; ref: required; exactly one of `uri`, `service`, `host_file`, `value`, `stream`.

| Param | Type | Default | |
|---|---|---|---|
| `uri` | `str \| null` | `null` | an absolute https:// URL |
| `service` | `str \| null` | `null` | a registered service |
| `id` | `str \| null` | `null` | the service's record |
| `resource` | `str \| null` | `null` | what to read from the service |
| `host_file` | `str \| null` | `null` | a name the host resolves |
| `value` | `str \| null` | `null` | a whole file in a registered text format |
| `format` | `str \| null` | `null` | a registered format; required with uri, host_file and value unless protocol is set |
| `protocol` | `"seqstar-tiles/1" \| null` | `null` | a streamed source by a registered protocol, at uri |
| `stream` | `str \| null` | `null` | a streamed source: an id the host's stream provider resolves |
| `data_kind` | `"letters" \| "items" \| "values" \| "alignment" \| null` | `null` | a streamed source's; required with stream |
| `coordinate_system` | `ref<sequence \| unit \| complex \| alignment \| row> \| null` | `null` | items, values: null, the target |
| `max_depth` | `int` | `4` | how deeply its patches may nest streamed sources |
| `max_bytes` | `int \| null` | `null` | a size bound; for value, null means 1 MiB |

### `sequence`

a polymer's residues. Parent: `data`; ref: required; exactly one of `value`, `source_ref`.

| Param | Type | Default | |
|---|---|---|---|
| `value` | `str \| null` | `null` | one-letter codes; no gaps |
| `alphabet` | `"protein" \| "dna" \| "rna" \| "generic" \| null` | `null` | null: detected |
| `cyclic` | `bool` | `false` | a cyclic coordinate system: positions modulo its length; ranges may wrap |
| `source_ref` | `ref<source> \| null` | `null` | read from a source node |
| `record` | `str \| int \| null` | `null` | a record id or 0-based index; required when the source has several |
| `chem_comp_ids` | `list<str \| null> \| null` | `null` | a chemical component per position |
| `label` | `str \| null` | `null` | null: from the source, else the ref |
| `identifiers` | `list<{namespace: str, value: str, version?: str \| null}>` | `[]` | uniprot:P04637, chembl:CHEMBL941; never imply a mapping |
| `provenance` | `{source: str, method?: str \| null, note?: str \| null, confidence?: float \| null, references?: list<str>, parameters?: object \| null} \| null` | `null` |  |

### `ligand`

a ligand's heavy atoms. Parent: `data`; ref: required; exactly one of `value`, `source_ref`.

| Param | Type | Default | |
|---|---|---|---|
| `value` | `list<str> \| null` | `null` | the element symbol of each heavy atom |
| `source_ref` | `ref<source> \| null` | `null` | read from a source node |
| `record` | `str \| int \| null` | `null` | a record id or 0-based index; required when the source has several |
| `chem_comp_ids` | `list<str \| null> \| null` | `null` | a chemical component per position |
| `label` | `str \| null` | `null` | null: from the source, else the ref |
| `identifiers` | `list<{namespace: str, value: str, version?: str \| null}>` | `[]` | uniprot:P04637, chembl:CHEMBL941; never imply a mapping |
| `provenance` | `{source: str, method?: str \| null, note?: str \| null, confidence?: float \| null, references?: list<str>, parameters?: object \| null} \| null` | `null` |  |

### `numbering`

a numbering of its parent. Parent: `sequence`, `ligand`, `unit`; ref: required.

A union on `method`: how labels are given.

`method: "offset"`:

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str` | required | unique among its parent's numberings; ordinal and row are reserved |
| `first` | `int` | `1` | the label of index 0 |
| `skip_zero` | `bool` | `false` | the numbering has no 0 |
| `default` | `bool` | `false` | the numbering rulers use when a section names none |

`method: "segments"`:

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str` | required |  |
| `segments` | `list<[int, int, int]>` | required | [index_start, label_start, length] |
| `default` | `bool` | `false` |  |

`method: "labels"`:

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str` | required |  |
| `labels` | `list<Label \| null>` | required | one per position; null: unlabeled |
| `default` | `bool` | `false` |  |

`method: "names"`:

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str` | required |  |
| `labels` | `list<Label \| null>` | required | one name per position; null: unlabeled |
| `default` | `bool` | `false` |  |

### `unit`

an instance of a sequence or a ligand. Parent: `data`, `complex`; ref: required; exactly one of `sequence`, `ligand`.

| Param | Type | Default | |
|---|---|---|---|
| `sequence` | `ref<sequence> \| null` | `null` | the polymer it instantiates |
| `ligand` | `ref<ligand> \| null` | `null` | the ligand it instantiates |
| `range` | `[int, int] \| null` | `null` | the shown positions; null: all |
| `observed` | `list<[int, int]> \| null` | `null` | observed positions; null: all of range |
| `label` | `str \| null` | `null` | the block's name |
| `origin` | `str \| null` | `null` | experimental, predicted, docked, designed, <ns>:<name>; sequence is reserved |
| `score` | `float \| null` | `null` | a score the instance carries |

### `complex`

units shown together: a provider's, or its unit children. Parent: `data`; ref: required.

| Param | Type | Default | |
|---|---|---|---|
| `provider` | `str \| null` | `null` | <namespace>:<name>, a registered complex provider; or unit children |
| `source_ref` | `ref<source> \| null` | `null` | required with provider; without, a streamed letters source |
| `options` | `object` | `{}` | the provider's, checked by its registered table |
| `origin` | `str \| null` | `null` | given to every unit; null: the provider's, or each unit's own |
| `label` | `str \| null` | `null` | null: the provider's, else the ref |

### `alignment`

an alignment, a coordinate system. Parent: `data`; ref: required; exactly one of `value`, `source_ref`, `provider`.

| Param | Type | Default | |
|---|---|---|---|
| `value` | `list<str> \| null` | `null` | the gapped rows |
| `flavor` | `"aligned" \| "a3m"` | `"aligned"` | a3m: lower case an insert, not a column |
| `labels` | `list<str> \| null` | `null` | value's row labels |
| `source_ref` | `ref<source> \| null` | `null` | read from a source node |
| `record` | `int \| null` | `null` | which alignment of a multi-alignment Stockholm file |
| `alphabet` | `"protein" \| "dna" \| "rna" \| "generic" \| null` | `null` | null: detected |
| `row_count` | `int \| null` | `null` | required for a protocol source |
| `column_count` | `int \| null` | `null` | required for a protocol source |
| `reference_row` | `int \| null` | `null` | the query row; null: row 0 |
| `parent_rows` | `list<int \| null> \| null` | `null` | per row: the row it is a variant of |
| `provider` | `str \| null` | `null` | a registered node provider (realign, or <namespace>:<name>): computed from inputs |
| `options` | `object` | `{}` | the provider's, checked by its registered table |
| `inputs` | `list<ref<sequence \| unit \| row>> \| ref<alignment> \| null` | `null` | with provider: what it aligns, in order: sequences, units or rows, or an alignment's rows |
| `label` | `str \| null` | `null` |  |

### `tree`

a tree over its alignment's rows. Parent: `alignment`; at most one under its parent; ref: required; exactly one of `value`, `source_ref`, `provider`.

| Param | Type | Default | |
|---|---|---|---|
| `value` | `str \| null` | `null` | Newick |
| `source_ref` | `ref<source> \| null` | `null` | read from a source node |
| `leaves` | `"labels" \| "order"` | `"labels"` | labels: by row label; order: leaf k is row k |
| `provider` | `str \| null` | `null` | a registered node provider (nj, upgma, or <namespace>:<name>): computed from the alignment |
| `options` | `object` | `{}` | the provider's, checked by its registered table |

### `mapping`

positions of one coordinate system paired with another's. Parent: `data`; ref: required.

A union on `method`: how positions are paired.

Every case:

| Param | Type | Default | |
|---|---|---|---|
| `provenance` | `{source: str, method?: str \| null, note?: str \| null, confidence?: float \| null, references?: list<str>, parameters?: object \| null} \| null` | `null` |  |

`method: "identity"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `from_range` | `[int, int] \| null` | `null` | the from indices it covers; null: all that fit |

`method: "segments"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `segments` | `list<[int, int, int \| null]>` | required | [from_start, to_start, length] |

`method: "pairwise"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `cigar` | `str` | required |  |

`method: "align"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `aligner` | `str` | `"biomsa"` | a registered aligner |
| `options` | `object` | `{}` | the aligner's params, checked by its table; {}: its defaults |

`method: "by_numbering"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `numbering` | `str` | required | pair equal labels |
| `to_numbering` | `str \| null` | `null` | null: the same name |

`method: "provider"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row \| alignment>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row \| alignment>` | required |  |
| `provider` | `str` | required | a registered node provider: codon_rows, or <namespace>:<name> |
| `options` | `object` | `{}` | the provider's, checked by its registered table |

`method: "codon"`:

| Param | Type | Default | |
|---|---|---|---|
| `from` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `to` | `ref<sequence \| ligand \| unit \| row>` | required |  |
| `from_range` | `[int, int] \| null` | `null` | the residues the CDS codes, the first with the first codon; null: from residue 0 |
| `cds` | `list<[int, int]> \| null` | `null` | the coding ranges on to (exons), in transcript order; null: one range from start |
| `start` | `int` | `0` | with cds null: the first index of the coding range |
| `strand` | `"+" \| "-"` | `"+"` | -: the ranges' reverse complement |
| `genetic_code` | `int` | `1` | an NCBI genetic code id (gc.prt; 1 the standard code) |

### `annotation`

a table of items; under a layer, that layer's own data. Parent: `data`, `layer`; ref: required; exactly one of `value`, `source_ref`, `uri`, `from_source`.

A union on `geometry`: the table's geometry.

Every case:

| Param | Type | Default | |
|---|---|---|---|
| `target` | `ref<sequence \| ligand \| unit \| alignment \| row> \| null` | `null` | every row's target without a target column |
| `numbering` | `str \| null` | `null` | required with label columns: their numbering |
| `source_ref` | `ref<source> \| null` | `null` | read from a source node |
| `uri` | `str \| null` | `null` | read from a URI directly; format required with it |
| `format` | `str \| null` | `null` | a registered format |
| `from_source` | `ref<sequence \| ligand \| complex \| alignment> \| null` | `null` | read from the file that node came from |
| `category_name` | `str \| null` | `null` | CIF: the category; null: the first |
| `block_header` | `str \| null` | `null` | CIF: the data block |
| `block_index` | `int` | `0` | CIF: the data block's index when block_header is null |
| `tag` | `str \| null` | `null` | Stockholm: GC:<tag>, GR:<tag>, GS:<tag>; FASTA: header |
| `fields` | `dict<str, str \| null>` | `{}` | column name -> name in the file, null to ignore one |
| `label` | `str \| null` | `null` |  |
| `semantic_type` | `str \| null` | `null` | open, namespaced: uniprot:Domain, so:0001067 |
| `length` | `int \| null` | `null` | required for a protocol source |
| `links` | `list<{label: str, url: str}>` | `[]` | per-row links |
| `provenance` | `{source: str, method?: str \| null, note?: str \| null, confidence?: float \| null, references?: list<str>, parameters?: object \| null} \| null` | `null` |  |

`geometry: "ranges"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<ranges> \| null` | `null` | the table |

`geometry: "values"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<values> \| null` | `null` | the table |

`geometry: "pairs"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<pairs> \| null` | `null` | the table |
| `b_numbering` | `str \| null` | `null` | side b's numbering; null: as numbering |

`geometry: "matrix"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<matrix> \| null` | `null` | the table |
| `states` | `list<str>` | required | the matrix's second axis |

`geometry: "grid"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<grid> \| null` | `null` | the table |
| `b_numbering` | `str \| null` | `null` | side b's numbering; null: as numbering |

`geometry: "rows"`:

| Param | Type | Default | |
|---|---|---|---|
| `value` | `Table<rows> \| null` | `null` | the table |

### `computed`

an annotation computed by a provider. Parent: `data`; ref: required.

A union on `provider`: a kernel provider or <namespace>:<name>.

`provider: "coverage"`:

| Param | Type | Default | |
|---|---|---|---|
| `of` | `ref<unit> \| null` | `null` | null: every unit |
| `reference` | `ref<sequence> \| null` | `null` | conflicts against it through a mapping |

`provider: "column_stats"`:

| Param | Type | Default | |
|---|---|---|---|
| `alignment` | `ref<alignment>` | required | the alignment |
| `row_subset` | `list<int> \| null` | `null` | a subgroup of rows |
| `gaps` | `"ignore" \| "count"` | `"count"` |  |
| `measure` | `"conservation" \| "consensus" \| "identity" \| "gap_fraction"` | required |  |
| `reference_row` | `int \| null` | `null` | null: the alignment's |

`provider: "profile"`:

| Param | Type | Default | |
|---|---|---|---|
| `alignment` | `ref<alignment>` | required | the alignment |
| `row_subset` | `list<int> \| null` | `null` | a subgroup of rows |
| `gaps` | `"ignore" \| "count"` | `"count"` |  |

`provider: "row_stats"`:

| Param | Type | Default | |
|---|---|---|---|
| `alignment` | `ref<alignment>` | required | the alignment |
| `reference_row` | `int \| null` | `null` | null: the alignment's |
| `items` | `ref<annotation \| computed> \| null` | `null` | count the items of this annotation on each row |
| `category` | `str \| null` | `null` | with items: only this category |
| `columns` | `[int, int] \| null` | `null` | with items: only within these columns, [start, stop) |
| `score` | `str \| null` | `null` | with items: the maximum of this field of the items on the row |

`provider: "redundancy"`:

| Param | Type | Default | |
|---|---|---|---|
| `alignment` | `ref<alignment>` | required | the alignment |
| `threshold` | `float (0..1)` | required | an identity from 0 to 1: a row joins the first representative this identical to it |
| `order` | `"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null` | `null` | the order rows are taken in; null: data order |

`provider: "paired"`:

| Param | Type | Default | |
|---|---|---|---|
| `pairs` | `list<ref<annotation>>` | required | pair annotations |
| `partner` | `list<ref<unit>> \| null` | `null` | null: every endpoint (an interface) |
| `categories` | `list<str> \| null` | `null` | only these pair categories |

`provider: "query"`:

| Param | Type | Default | |
|---|---|---|---|
| `query` | `SeqQLQuery` | required | a SeqQL query: { seqql, expr, text } |
| `label` | `str \| null` | `null` | null: the ref |

`provider: "<namespace>:<name>"`:

| Param | Type | Default | |
|---|---|---|---|
| `options` | `object` | `{}` | the provider's, checked by its registered table |

### `view`

one view. Parent: `root`; one or more under its parent; ref: required.

| Param | Type | Default | |
|---|---|---|---|
| `label` | `str` | required | what a reader picks it by |
| `layout` | `"tracks" \| "text"` | `"tracks"` | text: the blocks read as text, wrapping |

### `section`

one axis of the view with its tracks. Parent: `view`; one or more under its parent; ref: required.

| Param | Type | Default | |
|---|---|---|---|
| `label` | `str \| null` | `null` | shown when a view has several sections |
| `arrangement` | `"side" \| "stacked"` | `"side"` | blocks side by side, or one lane each |
| `copies` | `"one" \| "each"` | `"one"` | one block per copy group, or one per unit |
| `numbering` | `str \| null` | `null` | the numbering rulers and position labels use |
| `ligand_form` | `"position" \| "atoms"` | `"position"` | a ligand as one position or its atoms |
| `positions` | `ref<annotation \| computed> \| {data: ref<annotation \| computed>, filter: {field?: str \| null, values: list<str \| float>}} \| list<ref<annotation \| computed> \| {data: ref<annotation \| computed>, filter: {field?: str \| null, values: list<str \| float>}}> \| null` | `null` | show only what these cover |
| `scroll_group` | `str \| null` | `null` | sections of one group share horizontal scroll and zoom |
| `height` | `"content" \| "box"` | `"box"` | content: as tall as its rows |
| `height_rows` | `int \| null` | `null` | a hint for box |
| `row_sizing` | `"uniform" \| "natural"` | `"uniform"` |  |
| `style` | `str \| SectionStyle \| null` | `null` | a preset name (tracks, msa, strip, compact, wrapped, <ns>:<name>) or an inline preset |
| `direction` | `1 \| -1` | `1` | its blocks' direction, where a block gives none |
| `follow` | `ref<section> \| null` | `null` | this section and that one follow each other's window |
| `follow_anchor` | `ref<row \| unit> \| null` | `null` | an alignment axis's anchor; null: its reference_row |
| `highlights` | `list<{data: ref<annotation \| computed>, items?: list<str> \| null, style?: "band" \| "highlight" \| "accent"}>` | `[]` | persistent highlights, a channel each |
| `hidden` | `bool` | `false` |  |

### `block`

one block of the section's axis. Parent: `section`; ref: required; exactly one of `unit`, `sequence`, `ligand`, `complex`, `alignment`.

| Param | Type | Default | |
|---|---|---|---|
| `unit` | `ref<unit> \| null` | `null` |  |
| `sequence` | `ref<sequence> \| null` | `null` | shown through its own unit |
| `ligand` | `ref<ligand> \| null` | `null` | shown through its own unit |
| `complex` | `ref<complex> \| null` | `null` | its units in its order, a block each |
| `alignment` | `ref<alignment> \| null` | `null` | its columns |
| `label` | `str \| null` | `null` | null: the unit's, molecule's or alignment's |
| `numbering` | `str \| null` | `null` | overrides the section's |
| `range` | `[int, int] \| null` | `null` | with unit, sequence or ligand: the piece shown, [start, stop) of the unit (start > stop wraps a cyclic origin); null: all |
| `direction` | `1 \| -1 \| null` | `null` | -1: the piece from its highest index down; null: the section's |
| `origin` | `int` | `0` | a cyclic unit without range: the index its piece starts at |

### `group`

an expandable group of tracks. Parent: `section`, `group`; ref: required.

| Param | Type | Default | |
|---|---|---|---|
| `label` | `str` | required |  |
| `collapsed` | `bool` | `false` | the state at open |
| `hidden` | `bool` | `false` |  |

### `track`

one row of the section. Parent: `section`, `group`; ref: required.

| Param | Type | Default | |
|---|---|---|---|
| `label` | `str \| null` | `null` | null: the first layer's data label |
| `height` | `float \| null` | `null` | a hint, in rows |
| `hidden` | `bool` | `false` |  |
| `expanded` | `bool \| null` | `null` | a split layer's rows at open; null: expanded for variants, else collapsed |
| `summary` | `bool` | `false` | its group's summary: shown under the group's header only while the group is collapsed (a union of what it holds) |
| `for_each` | `{kind: "sequence" \| "ligand" \| "unit" \| "mapping" \| "annotation" \| "item_system", source: ref<source>, depth?: int \| null, reaching?: bool, order?: "data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null, filter?: list<{data: ref<annotation \| computed>, field: str, values?: list<str \| float \| bool> \| null, min?: float \| null, max?: float \| null}>} \| null` | `null` | a track template: one track per dynamic node it selects |

### `layer`

what a track shows. Parent: `track`; one or more under its parent; ref: required.

A union on `type`: a kernel layer type or <namespace>:<name>.

Every case:

| Param | Type | Default | |
|---|---|---|---|
| `data` | `ref<annotation \| computed> \| list<ref<annotation \| computed>> \| null` | `null` | annotations or computed nodes, drawn as one; absent with annotation children |
| `blocks` | `list<ref<block>> \| null` | `null` | only these blocks |
| `color` | `ThemeRef<color> \| null` | `null` | null: the type's default |
| `glyph` | `ThemeRef<glyph> \| null` | `null` | null: the type's default |
| `domain_group` | `str \| null` | `null` | layers of the view in one group share their normalized domains (the union, per field) |
| `opacity` | `float` | `1` |  |
| `filter` | `{field?: str \| null, values: list<str \| float>} \| null` | `null` | only the rows the filter keeps |
| `fallback` | `str \| null` | `null` | a kernel type or omit; null: the type's |
| `options` | `object` | `{}` | an extension type's params |

`type: "letters"`:

| Param | Type | Default | |
|---|---|---|---|
| `row` | `int \| null` | `null` | an alignment row's letters; null: the block's (the reference row's) |
| `of` | `"mapped" \| ref<sequence \| unit> \| null \| list<ref<sequence \| unit>>` | `null` | other polymers' letters where mapping paths place them: one, a list, or mapped (every codon mapping's protein) |
| `residue_form` | `"auto" \| "one" \| "three"` | `"auto"` | a protein's residues as M or Met; a hint |
| `letter_map` | `str \| dict<str, str> \| null` | `null` | the letter shown for each stored one: a registered map by name (complement, complement_rna, lower) or single letters mapped; null: as stored |

`type: "rows"`:

| Param | Type | Default | |
|---|---|---|---|
| `filter` | `list<{data: ref<annotation \| computed>, field: str, values?: list<str \| float \| bool> \| null, min?: float \| null, max?: float \| null}>` | `[]` | rows kept when every filter keeps them: {data, field, values, min, max} over a key table |
| `order` | `"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null` | `null` | tree, data, a key table's field ({data, field, direction}) or a list of them, lexicographic; null: tree when there is one |
| `show_tree` | `bool` | `true` | drawn while the order is the tree |
| `row_subset` | `list<int> \| null` | `null` | only these rows, applied before the filter |
| `data_columns` | `list<{data: ref<annotation \| computed>, field: str, label?: str \| null, type?: str}>` | `[]` | columns beside the rows |
| `letter_map` | `str \| dict<str, str> \| null` | `null` | the letter shown for each stored one: a registered map by name (complement, complement_rna, lower) or single letters mapped; null: as stored |

`type: "ranges"`:

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `line_style` | `ThemeRef<line_style> \| null` | `null` | a line style theme: how each item is drawn (solid, dashed, dotted, wavy, double, box, highlight); null: solid |
| `display` | `"lanes" \| "brackets" \| "underline" \| "union"` | `"lanes"` | a hint |
| `labels` | `bool` | `false` | each item's label inside its feature where it fits (on rings: outside the outermost ring, with a leader); a hint |
| `direction_field` | `str \| null` | `null` | a field of the rows giving each item's direction (+, 1, forward; -, -1, reverse): an arrowhead at its end |
| `collapse` | `str \| null` | `null` | a field of the rows: where the track's lanes do not fit, items sharing its value draw as one item in one lane (a gene's transcripts), until the track is expanded; streamed items; a hint |
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |

`type: "points"`:

| Param | Type | Default | |
|---|---|---|---|
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |

`type: "emphasis"`:

No params.

`type: "variants"`:

| Param | Type | Default | |
|---|---|---|---|
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |

`type: "secondary"`:

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `labels` | `bool` | `false` | each item's label over its shape where it fits (α1, β2); a hint |

`type: "values"`:

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `encoding` | `"bars" \| "line" \| "area"` | `"bars"` | a hint |
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |

`type: "color"`:

No params.

`type: "ruler"`:

| Param | Type | Default | |
|---|---|---|---|
| `of` | `ref<sequence \| unit \| row> \| null \| list<ref<sequence \| unit \| row>>` | `null` | the coordinate system whose numbering it shows, placed as data is; a list: the first that maps a position |
| `numbering` | `str \| null \| list<str>` | `null` | null: the section's or block's, or of's default; a list: the first that labels a position |

`type: "pairs"`:

| Param | Type | Default | |
|---|---|---|---|
| `width` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `height` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `display` | `"arcs" \| "ends" \| "both"` | `"both"` | a hint |
| `side` | `"above" \| "below" \| null` | `null` | null: toward the letters (the side facing the nearest letters track; above when there is none); a hint |
| `arcs` | `"always" \| "hover"` | `"always"` | a hint |
| `labels` | `bool` | `false` | each arc's label, else its width theme's field's value, inside it under its apex where it fits; a hint |

`type: "logo"`:

| Param | Type | Default | |
|---|---|---|---|
| `scale` | `"bits" \| "probability"` | `"bits"` | a hint |

`type: "matrix"`:

| Param | Type | Default | |
|---|---|---|---|
| `filter` | `list<{data: ref<annotation \| computed>, field: str, values?: list<str \| float \| bool> \| null, min?: float \| null, max?: float \| null}>` | `[]` | states kept when every filter keeps them: {data, field, values, min, max} over a table keyed by state |
| `height` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `density` | `"majority" \| "maximum"` | `"majority"` | states thinner than a device pixel, per device row: majority (the most common color, faded by the share of states with one) or maximum (the state with the largest value, so a rare signal shows); a hint |
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |
| `order` | `"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null` | `null` | the states' order: data, a key table's field ({data, field, direction}) over a table keyed by state, or a list of them; null: data |

`type: "dot_plot"`:

| Param | Type | Default | |
|---|---|---|---|
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |
| `aspect_ratio` | `float \| null` | `1` | the plot's width ÷ height (1: square), its height following the data area's width within what the section shows, never flatter than a viewer's stretch of it; null: the track's own height; a hint |

`type: "frame"`:

| Param | Type | Default | |
|---|---|---|---|
| `span` | `{from: ref<track>, to?: ref<track> \| null}` | required | the tracks each item's frame encloses, from the top of from's band to the bottom of to's (a rows track: the rows the item targets) |
| `line_style` | `ThemeRef<line_style> \| null` | `null` | a line style theme: solid, dashed or dotted (the others drawn solid); null: solid |
| `labels` | `bool` | `true` | each item's label as its frame's title, in the frame's own track |

`type: "<namespace>:<name>"`:

No params.

### `focus`

the section's initial range. Parent: `section`; ref: required; exactly one of `target`, `annotation`.

| Param | Type | Default | |
|---|---|---|---|
| `target` | `ref<unit \| sequence \| ligand \| alignment \| row> \| null` | `null` | with start, stop or beg, end |
| `start` | `int \| null` | `null` | an index; columns for an alignment |
| `stop` | `int \| null` | `null` | an index |
| `beg` | `Label \| null` | `null` | a label in numbering |
| `end` | `Label \| null` | `null` | a label in numbering |
| `numbering` | `str \| null` | `null` | the labels' numbering |
| `annotation` | `ref<annotation \| computed> \| null` | `null` | its items' extent |
| `items` | `list<str> \| null` | `null` | item ids; null: every item |

## Types

### `Table<geometry>`

Rows of named columns, written as an array of row objects or as an object of column arrays; a column may be `{dict, codes}` or `{json}`. The columns each geometry reads (a target or position column of another geometry is `svs.columns.unknown`; any other column is a field):

| Geometry | Columns |
|---|---|
| `ranges` | `target`, `start`, `stop`, `beg`, `end`, `ins_code`, `boundary`, `value`, `id`, `category`, `label`, `tooltip`, `alt` |
| `values` | `target`, `start`, `stop`, `beg`, `end`, `index`, `pos`, `ins_code`, `value`, `category`, `label`, `tooltip` |
| `pairs` | `target`, `b_target`, `start`, `stop`, `beg`, `end`, `ins_code`, `b_start`, `b_stop`, `b_beg`, `b_end`, `b_ins_code`, `id`, `category`, `label`, `tooltip`, `numbering`, `b_numbering`, `direction` |
| `matrix` | `target`, `index`, `pos`, `ins_code`, `value` |
| `grid` | `target`, `b_target`, `start`, `stop`, `beg`, `end`, `index`, `pos`, `ins_code`, `b_start`, `b_stop`, `b_beg`, `b_end`, `b_ins_code`, `value`, `id`, `label`, `tooltip`, `numbering`, `b_numbering` |
| `rows` | `row`, `state`, `node` |

### `Palette`

A union on `kind`: the palette's kind.

`kind: "categorical"`:

| Param | Type | Default | |
|---|---|---|---|
| `colors` | `str \| list<Color> \| dict<str, Color>` | `"Category10"` | a list name, a dictionary name, colors, or {value: Color} |
| `repeat_color_list` | `bool` | `true` | start the list over when it runs out (MVS: false) |
| `sort` | `"none" \| "lexical" \| "numeric"` | `"none"` | none: first occurrence |
| `sort_direction` | `"ascending" \| "descending"` | `"ascending"` |  |
| `case_insensitive` | `bool` | `false` |  |
| `missing_color` | `Color \| null` | `null` | null: not colored by this layer |

`kind: "discrete"`:

| Param | Type | Default | |
|---|---|---|---|
| `colors` | `str \| list<Color> \| list<[Color, float]> \| list<[Color \| null, float \| null, float \| null]>` | `"YlGn"` | a list name, colors, [Color, checkpoint] pairs, or [Color \| null, from \| null, to \| null] ranges |
| `reverse` | `bool` | `false` | reverse the list |
| `mode` | `"normalized" \| "absolute"` | `"normalized"` | normalized: over value_domain |
| `value_domain` | `[float \| null, float \| null]` | `[null,null]` | [min, max]; null: from the data |

`kind: "continuous"`:

| Param | Type | Default | |
|---|---|---|---|
| `colors` | `str \| list<Color> \| list<[Color, float]>` | `"YlGn"` | a list name, colors, or [Color, checkpoint] pairs |
| `reverse` | `bool` | `false` | reverse the list |
| `mode` | `"normalized" \| "absolute"` | `"normalized"` | normalized: over value_domain |
| `value_domain` | `[float \| null, float \| null]` | `[null,null]` | [min, max]; null: from the data |
| `underflow_color` | `"auto" \| Color \| null` | `null` | "auto", a color, or null: not colored |
| `overflow_color` | `"auto" \| Color \| null` | `null` | "auto", a color, or null: not colored |

### `Glyph`

A shape name (`circle`, `diamond`, `triangle`, `square`, `pin`, `star`, `cross`, `tick`) or an object with exactly one of `shape`, `char` (one grapheme) and `path` (SVG path data, at most 4 KB, with `box`: [w, h]), and optionally `label` and `fallback` (a shape).

### `LineStyle`

How a range's item is drawn: `solid`, `dashed`, `dotted`, `wavy`, `double`, `box`, `highlight`; `solid` is the plain bar or underline.

### `Color`

`#rgb`, `#rrggbb`, a Mol* color name, or a split color `a/b`. No alpha: a layer's `opacity` gives it.

### `Label`

A position label: an int (`315`), or a string without whitespace: `"27A"`, `"111.1"` in an integer numbering, a name (`"N13"`) in a `names` numbering.

### `ThemeRef`

`{theme, params, fallback}`: a registered theme with the params of its table (at most 4 deep); a `<namespace>:<name>` theme needs a kernel `fallback`.

## Layer types

| Type | Shows | Data | Where | Axes | Sizing | Fallback |
|---|---|---|---|---|---|---|
| `letters` | residues (element symbols for atoms), one alignment row's, or another polymer's through a mapping | none | tracks, text | blocks, columns | row | none |
| `rows` | an alignment's rows in tree order, with the tree and data columns | none | tracks, text | columns | free | letters |
| `ranges` | labeled ranges, discontinuous items, boundaries | ranges | tracks, text | blocks, columns | row | none |
| `points` | point features | ranges | tracks, text | blocks, columns | row | ranges |
| `emphasis` | positions emphasized in the letters or rows of its track | ranges | tracks, text | blocks, columns | row | points |
| `variants` | substitutions per position | ranges | tracks | blocks, columns | row | points |
| `secondary` | helices, strands, turns | ranges | tracks, text | blocks, columns | row | ranges |
| `values` | a number per position or column | values | tracks, text | blocks, columns | row | color |
| `color` | a color per position, from categories or values | ranges, values | tracks, text | blocks, columns | row | ranges |
| `ruler` | positions labeled in a numbering | none | tracks | blocks, columns | free | omit |
| `pairs` | pairs of positions, by category | pairs | tracks | blocks, columns | free | points |
| `logo` | letter frequencies per column | matrix | tracks | columns | free | values |
| `matrix` | a value per position and state | matrix | tracks | blocks, columns | free | values |
| `dot_plot` | positions against positions on two axes: pairs as segments, a grid's values as cells | pairs, grid | tracks | blocks | free | omit |
| `frame` | named regions outlined across tracks, titled | ranges | tracks, text | blocks, columns | row | ranges |

### `letters`

| Param | Type | Default | |
|---|---|---|---|
| `row` | `int \| null` | `null` | an alignment row's letters; null: the block's (the reference row's) |
| `of` | `"mapped" \| ref<sequence \| unit> \| null \| list<ref<sequence \| unit>>` | `null` | other polymers' letters where mapping paths place them: one, a list, or mapped (every codon mapping's protein) |
| `residue_form` | `"auto" \| "one" \| "three"` | `"auto"` | a protein's residues as M or Met; a hint |
| `letter_map` | `str \| dict<str, str> \| null` | `null` | the letter shown for each stored one: a registered map by name (complement, complement_rna, lower) or single letters mapped; null: as stored |

### `rows`

| Param | Type | Default | |
|---|---|---|---|
| `order` | `"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null` | `null` | tree, data, a key table's field ({data, field, direction}) or a list of them, lexicographic; null: tree when there is one |
| `filter` | `list<{data: ref<annotation \| computed>, field: str, values?: list<str \| float \| bool> \| null, min?: float \| null, max?: float \| null}>` | `[]` | rows kept when every filter keeps them: {data, field, values, min, max} over a key table |
| `show_tree` | `bool` | `true` | drawn while the order is the tree |
| `row_subset` | `list<int> \| null` | `null` | only these rows, applied before the filter |
| `data_columns` | `list<{data: ref<annotation \| computed>, field: str, label?: str \| null, type?: str}>` | `[]` | columns beside the rows |
| `letter_map` | `str \| dict<str, str> \| null` | `null` | the letter shown for each stored one: a registered map by name (complement, complement_rna, lower) or single letters mapped; null: as stored |

### `ranges`

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `line_style` | `ThemeRef<line_style> \| null` | `null` | a line style theme: how each item is drawn (solid, dashed, dotted, wavy, double, box, highlight); null: solid |
| `display` | `"lanes" \| "brackets" \| "underline" \| "union"` | `"lanes"` | a hint |
| `labels` | `bool` | `false` | each item's label inside its feature where it fits (on rings: outside the outermost ring, with a leader); a hint |
| `direction_field` | `str \| null` | `null` | a field of the rows giving each item's direction (+, 1, forward; -, -1, reverse): an arrowhead at its end |
| `collapse` | `str \| null` | `null` | a field of the rows: where the track's lanes do not fit, items sharing its value draw as one item in one lane (a gene's transcripts), until the track is expanded; streamed items; a hint |
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |

### `points`

| Param | Type | Default | |
|---|---|---|---|
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |

### `emphasis`

No params.

### `variants`

| Param | Type | Default | |
|---|---|---|---|
| `split` | `{by: str, order?: list<str> \| null, empty?: "show" \| "hide"} \| null` | `null` | one labeled row per value of the field; the track's expanded opens them |
| `density` | `"stacked" \| "majority"` | `"stacked"` | a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint |
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |

### `secondary`

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `labels` | `bool` | `false` | each item's label over its shape where it fits (α1, β2); a hint |

### `values`

| Param | Type | Default | |
|---|---|---|---|
| `thickness` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `encoding` | `"bars" \| "line" \| "area"` | `"bars"` | a hint |
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |

### `color`

No params.

### `ruler`

| Param | Type | Default | |
|---|---|---|---|
| `of` | `ref<sequence \| unit \| row> \| null \| list<ref<sequence \| unit \| row>>` | `null` | the coordinate system whose numbering it shows, placed as data is; a list: the first that maps a position |
| `numbering` | `str \| null \| list<str>` | `null` | null: the section's or block's, or of's default; a list: the first that labels a position |

### `pairs`

| Param | Type | Default | |
|---|---|---|---|
| `width` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `height` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `radius` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `display` | `"arcs" \| "ends" \| "both"` | `"both"` | a hint |
| `side` | `"above" \| "below" \| null` | `null` | null: toward the letters (the side facing the nearest letters track; above when there is none); a hint |
| `arcs` | `"always" \| "hover"` | `"always"` | a hint |
| `labels` | `bool` | `false` | each arc's label, else its width theme's field's value, inside it under its apex where it fits; a hint |

### `logo`

| Param | Type | Default | |
|---|---|---|---|
| `scale` | `"bits" \| "probability"` | `"bits"` | a hint |

### `matrix`

| Param | Type | Default | |
|---|---|---|---|
| `height` | `ThemeRef<size> \| null` | `null` | a size theme; null: the viewer's size |
| `density` | `"majority" \| "maximum"` | `"majority"` | states thinner than a device pixel, per device row: majority (the most common color, faded by the share of states with one) or maximum (the state with the largest value, so a rare signal shows); a hint |
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |
| `order` | `"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"} \| list<"data" \| "tree" \| {data: ref<annotation \| computed>, field: str, direction?: "asc" \| "desc"}> \| null` | `null` | the states' order: data, a key table's field ({data, field, direction}) over a table keyed by state, or a list of them; null: data |
| `filter` | `list<{data: ref<annotation \| computed>, field: str, values?: list<str \| float \| bool> \| null, min?: float \| null, max?: float \| null}>` | `[]` | states kept when every filter keeps them: {data, field, values, min, max} over a table keyed by state |

### `dot_plot`

| Param | Type | Default | |
|---|---|---|---|
| `domain` | `[float, float] \| null` | `null` | null: over all the layer's data; outside: clamped |
| `aspect_ratio` | `float \| null` | `1` | the plot's width ÷ height (1: square), its height following the data area's width within what the section shows, never flatter than a viewer's stretch of it; null: the track's own height; a hint |

### `frame`

| Param | Type | Default | |
|---|---|---|---|
| `span` | `{from: ref<track>, to?: ref<track> \| null}` | required | the tracks each item's frame encloses, from the top of from's band to the bottom of to's (a rows track: the rows the item targets) |
| `line_style` | `ThemeRef<line_style> \| null` | `null` | a line style theme: solid, dashed or dotted (the others drawn solid); null: solid |
| `labels` | `bool` | `true` | each item's label as its frame's title, in the frame's own track |

## Color themes

### `uniform`

one value

| Param | Type | Default | |
|---|---|---|---|
| `color` | `Color` | required | the color |

### `rules`

the first rule that holds

| Param | Type | Default | |
|---|---|---|---|
| `rules` | `list<{when?: list<{measure?: "value" \| {field: str} \| {column: ref<annotation \| computed>} \| {share: str} \| {rank: "letter"} \| {equals_row: int \| null} \| null, min?: float \| null, max?: float \| null, in?: list<str> \| null}>, color?: Color \| null, theme?: ThemeRef<color> \| null, alpha?: "value" \| {field: str} \| {column: ref<annotation \| computed>} \| {share: str} \| {rank: "letter"} \| {equals_row: int \| null} \| null, label?: str \| null}>` | required | the first that holds wins |
| `default` | `ThemeRef<color> \| null` | `null` | null: no value |
| `profile` | `ref<computed> \| null` | `null` | what share and rank read; null: the alignment's columns |

### `field_palette`

the row's field through the palette

| Param | Type | Default | |
|---|---|---|---|
| `field` | `str \| null` | `null` | null: value for values and matrix; else category, or the item (id, label) when there is no category column |
| `palette` | `Palette` | `{"kind":"categorical","colors":"Category10"}` |  |

### `field_color`

the field's values are Colors; others not colored

| Param | Type | Default | |
|---|---|---|---|
| `field` | `str` | `"color"` |  |

### `letter`

the residue letter, or a ligand atom's element, through a letter dictionary

| Param | Type | Default | |
|---|---|---|---|
| `colors` | `str \| dict<str, Color> \| null` | `null` | a dictionary name or {letter: Color}; null: the alphabet's |

### `clustal_x`

Clustal X by the column's composition (protein); a named rule table

No params.

### `percent_identity`

the consensus letter where more than 80, 60, 40% of the column agree (Jalview's); a named rule table

No params.

### `espript`

ESPript's shading: identical columns in one color, a letter whose group holds at least similarity of its column in another; none elsewhere; a named rule table

| Param | Type | Default | |
|---|---|---|---|
| `identical` | `Color` | `"#000000"` | where every row of the column has the letter |
| `similar` | `Color` | `"#bdbdbd"` | where the letter's group holds at least similarity of the column |
| `similarity` | `float` | `0.7` | the share of the column, from 0 to 1 |
| `groups` | `list<str>` | `["ILMV","FWY","KRH","DE","ST","NQ","AG"]` | letter groups read as similar |

### `scheme`

a color scheme by name, or an inline one: the letter, field_palette or rules theme it stands for, its dark variant on a dark skin

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str \| null` | `null` | a registered color scheme; exactly one of name and custom |
| `custom` | `SchemeDef \| null` | `null` | an inline scheme: letters, sequential, diverging or categorical |
| `overrides` | `dict<str, Color>` | `{}` | a letters scheme's letters recolored |
| `field` | `str \| null` | `null` | values and categories: the field, as field_palette's |
| `domain` | `[float \| null, float \| null] \| null` | `null` | values: [min, max]; null: the layer's |
| `reverse` | `bool` | `false` | values: the ramp reversed |

### `by_unit`

one color per unit, by unit ref in data order

| Param | Type | Default | |
|---|---|---|---|
| `palette` | `Palette` | `{"kind":"categorical","colors":"Category10"}` |  |

### `by_copy_group`

one color per copy group, by its first unit

| Param | Type | Default | |
|---|---|---|---|
| `palette` | `Palette` | `{"kind":"categorical","colors":"Category10"}` |  |

### `identity_to_row`

the base color where the letter equals that row's in the column; none elsewhere; a named rule table

| Param | Type | Default | |
|---|---|---|---|
| `row` | `int \| null` | `null` | null: the reference row |
| `base` | `ThemeRef<color>` | `{"theme":"letter"}` | the theme it shades |

### `conservation_shading`

the base color faded by the column's value where it is at least threshold; none below; a named rule table

| Param | Type | Default | |
|---|---|---|---|
| `data` | `ref<annotation \| computed>` | required | values over the alignment's columns |
| `threshold` | `float` | `0` |  |
| `base` | `ThemeRef<color>` | `{"theme":"letter"}` | the theme it shades |

## Glyph themes

### `uniform`

one value

| Param | Type | Default | |
|---|---|---|---|
| `glyph` | `Glyph` | required | the glyph |

### `rules`

the first rule that holds

| Param | Type | Default | |
|---|---|---|---|
| `rules` | `list<{when?: list<{measure?: "value" \| {field: str} \| {column: ref<annotation \| computed>} \| {share: str} \| {rank: "letter"} \| {equals_row: int \| null} \| null, min?: float \| null, max?: float \| null, in?: list<str> \| null}>, glyph?: Glyph \| null, theme?: ThemeRef<glyph> \| null, label?: str \| null}>` | required | the first that holds wins |
| `default` | `ThemeRef<glyph> \| null` | `null` | null: no value |
| `profile` | `ref<computed> \| null` | `null` | what share and rank read; null: the alignment's columns |

### `field_map`

the field's value through a map

| Param | Type | Default | |
|---|---|---|---|
| `field` | `str \| null` | `null` | null: category |
| `map` | `dict<str, Glyph>` | required | value -> glyph |
| `default` | `Glyph \| null` | `null` | for values the map lacks; null: none |

### `field_glyph`

the field's values are shape names or one character

| Param | Type | Default | |
|---|---|---|---|
| `field` | `str` | `"glyph"` |  |

## Line style themes

### `uniform`

one value

| Param | Type | Default | |
|---|---|---|---|
| `line_style` | `LineStyle` | required | the line style |

### `field_line_style`

the field's values are line styles; others solid, said once (svs.render.line_style)

| Param | Type | Default | |
|---|---|---|---|
| `field` | `str` | `"line_style"` |  |

## Letter maps

A `letter_map` (`letters`, `rows`): a name below, or a dict of single letters mapped; unlisted letters are kept.
A map listing alphabets leaves a layer with letters of another as stored.

- `complement` (dna, rna): A→T a→t T→A t→a U→A u→a C→G c→g G→C g→c R→Y r→y Y→R y→r K→M k→m M→K m→k B→V b→v V→B v→b D→H d→h H→D h→d
- `complement_rna` (dna, rna): A→U a→u U→A u→a T→A t→a C→G c→g G→C g→c R→Y r→y Y→R y→r K→M k→m M→K m→k B→V b→v V→B v→b D→H d→h H→D h→d
- `lower`: A→a B→b C→c D→d E→e F→f G→g H→h I→i J→j K→k L→l M→m N→n O→o P→p Q→q R→r S→s T→t U→u V→v W→w X→x Y→y Z→z

## Color schemes

A `scheme` color theme names one (`{"theme": "scheme", "params": {"name": "zappo"}}`); CB: colorblind-safe.

| Scheme | Label | Kind | CB | What it shows |
|---|---|---|---|---|
| `clustal_x` | Clustal X | theme: protein |  | residues by group where the column is conserved |
| `zappo` | Zappo | letters: protein |  | residues by physicochemical property |
| `taylor` | Taylor | letters: protein |  | each residue its own hue, similar residues near |
| `hydrophobicity` | Hydrophobicity | letters: protein |  | Kyte–Doolittle: red most hydrophobic, blue most hydrophilic |
| `percent_identity` | Percent identity | theme: protein, dna, rna, generic | yes | the consensus residue where more than 80, 60, 40% of the column agree |
| `espript` | ESPript | theme: protein | yes | identical columns black, similar residues grey (ESPript's shading, by residue groups) |
| `nucleotide` | Nucleotide | letters: dna, rna |  | A, C, G and T or U each a color |
| `purine_pyrimidine` | Purine/pyrimidine | letters: dna, rna |  | purines and pyrimidines |
| `nucleotide_cb` | Nucleotide, Okabe–Ito | letters: dna, rna | yes | A, C, G and T or U from the Okabe–Ito colors |
| `viridis` | Viridis | sequential | yes | perceptually uniform, dark low to yellow high |
| `cividis` | Cividis | sequential | yes | perceptually uniform for color vision deficiency, blue to yellow |
| `magma` | Magma | sequential | yes | perceptually uniform, black to cream |
| `blues` | Blues | sequential | yes | light low to dark blue high; the format's sequential default |
| `rd_bu` | Red–blue | diverging | yes | diverging at the middle: red low, blue high; the format's diverging default |
| `category10` | Category 10 | categorical |  | ten categories; the format's default |
| `okabe_ito` | Okabe–Ito | categorical | yes | seven categories distinguishable with color vision deficiency |
| `tableau10` | Tableau 10 | categorical |  | ten categories |

## Custom tables

### `mvs` on `complex`

| Param | Type | Default | |
|---|---|---|---|
| `structure_ref` | `str \| null` | `null` | the MVS ref of the structure node (the application writes it) |

### `mvs` on `unit`

| Param | Type | Default | |
|---|---|---|---|
| `structure_ref` | `str \| null` | `null` | the MVS ref of the structure node (the application writes it) |
| `label_asym_id` | `str \| null` | `null` | MVS's label_asym_id |
| `auth_asym_id` | `str \| null` | `null` | MVS's auth_asym_id |
| `instance_id` | `str \| null` | `null` | MVS's instance_id: the operator |

## Issue codes

| Code | Severity |
|---|---|
| `svs.format.invalid` | error |
| `svs.version.missing` | error |
| `svs.version.unsupported_major` | error |
| `svs.node.kind` | error |
| `svs.node.parent` | error |
| `svs.node.count` | error |
| `svs.params.missing` | error |
| `svs.params.unknown` | error |
| `svs.params.type` | error |
| `svs.params.exclusive` | error |
| `svs.ref.syntax` | error |
| `svs.ref.missing` | error |
| `svs.ref.duplicate` | error |
| `svs.ref.unresolved` | error |
| `svs.ref.wrong_kind` | error |
| `svs.ref.cycle` | error |
| `svs.ref.derived_unknown` | error |
| `svs.custom.namespace` | error |
| `svs.custom.invalid` | error |
| `svs.requires.unknown` | error |
| `svs.position.syntax` | error |
| `svs.position.form` | error |
| `svs.position.range` | error |
| `svs.columns.length` | error |
| `svs.columns.unknown` | error |
| `svs.columns.type` | error |
| `svs.color.invalid` | error |
| `svs.glyph.invalid` | error |
| `svs.line_style.invalid` | error |
| `svs.source.uri` | error |
| `svs.layer.data` | error |
| `svs.layer.fallback` | error |
| `svs.layer.annotation` | error |
| `svs.layer.split` | error |
| `svs.section.axis` | error |
| `svs.section.positions` | error |
| `svs.section.scroll_group` | error |
| `svs.tree.leaves` | error |
| `svs.patch.base` | error |
| `svs.patch.ref` | error |
| `svs.patch.anchor` | error |
| `svs.patch.cycle` | error |
| `svs.snapshot.ref` | error |
| `svs.source.load` | error |
| `svs.source.format` | error |
| `svs.source.record` | error |
| `svs.source.category` | error |
| `svs.source.size_bound` | error |
| `svs.sequence.alphabet` | error |
| `svs.ligand.unsupported` | error |
| `svs.numbering.length` | error |
| `svs.numbering.unknown` | error |
| `svs.alignment.ragged` | error |
| `svs.alignment.shape` | error |
| `svs.mapping.cigar` | error |
| `svs.mapping.range` | error |
| `svs.mapping.codon` | error |
| `svs.mapping.ambiguous` | error |
| `svs.focus.range` | error |
| `svs.theme.channel` | error |
| `svs.tree.streamed` | error |
| `svs.numbering.invalid` | error |
| `svs.section.follow` | error |
| `svs.scheme.unknown` | error |
| `svs.scheme.letter` | error |
| `svs.query.invalid` | error |
| `svs.query.version` | error |
| `svs.query.search` | error |
| `svs.query.input` | error |
| `svs.version.newer_minor` | warning |
| `svs.snapshot.invalid` | warning |
| `svs.position.unknown` | warning |
| `svs.numbering.descending` | warning |
| `svs.mapping.empty` | warning |
| `svs.numbering.fallback` | warning |
| `svs.provider.unavailable` | warning |
| `svs.layer.unplaced` | warning |
| `svs.layer.insert` | warning |
| `svs.render.fallback` | warning |
| `svs.render.omitted` | warning |
| `svs.section.style` | warning |
| `svs.render.theme` | warning |
| `svs.render.glyph` | warning |
| `svs.render.line_style` | warning |
| `svs.render.preset` | warning |
| `svs.render.direction` | warning |
| `svs.alphabet.detected` | warning |
| `svs.alignment.parent` | warning |
| `svs.annotation.rows` | warning |
| `svs.unit.orphaned` | warning |
| `svs.unit.stale` | warning |
| `svs.mapping.orphaned` | warning |
| `svs.annotation.orphaned` | warning |
| `svs.provider.unmapped` | warning |
| `svs.provider.size` | warning |
| `svs.provider.lazy` | warning |
| `svs.provider.alphabet` | warning |
| `svs.coverage.ambiguous` | warning |
| `svs.layer.rows` | warning |
| `svs.stream.ref` | warning |
| `svs.stream.conflict` | warning |
| `svs.stream.depth` | warning |
| `svs.stream.fallback` | warning |
| `svs.query.stale` | warning |
| `svs.query.warning` | warning |
| `svs.layer.streamed_merged` | warning |
| `svs.stream.item_map` | warning |
| `svs.mapping.align` | warning |
| `svs.source.id` | warning |
| `svs.stream.load` | warning |
| `svs.layer.span` | warning |
| `svs.alignment.inserts` | info |
| `svs.source.empty_record` | info |
| `svs.source.skipped` | info |
| `svs.source.unterminated` | info |
| `svs.layer.empty` | info |
| `svs.ref.pending` | info |
| `svs.query.note` | info |
| `svs.legend.unlisted` | info |
