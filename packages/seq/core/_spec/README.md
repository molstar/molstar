# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.0.0` until the first slice runs). Why: one language, one
  description, one number for readers to check.

## Queries

`src/query/` (`@molstar/seqx-core/query`): locations, element keys, the three result types, and the query forms as
`repr` form resolvers over a model of systems, numberings and annotations. Forms and addresses are one language: a form
is a segment of an address, `/data/p53/(range :from 94 :to 292)`.

### Result types

- **Nodes** (exact addresses of spec nodes and of nodes data makes), **Elements** (keys of what data makes: an
  annotation's items) and **Location** (positions on systems). A form gives one of them; they travel as `repr`'s `Named`
  (`{ kind: 'nodes' }`, `{ kind: 'other', type: 'elements' | 'location', value }`). A `Named` whose value lacks its
  type's shape is no query result (`fromNamed` checks it). Why: a hover, a click and a layer each hold the type they act
  on, and no form guesses a conversion.
- **Conversions are forms**: `(extent :of <elements>)` gives a Location, `(items :overlapping <location>)` Elements. A
  form given another type is `seq.query.form.type`.

### Location

- `{ parts: { system, ranges: Int32Array }[], rects: { a: part, b: part }[] }`: one part per system (sorted by address),
  its ranges sorted, merged (adjacent ones too), half-open `[start, end)` index pairs, none empty; a cyclic system's
  wrapped range is two ranges. Rects pair two parts (two axes), sorted, each once.
- **A value**: built normalized (`location`), frozen, never mutated. Equality is structural (`locationsEqual`): equal
  parts and an equal list of rects. Rects are never merged, so two rects covering the cells of one are another value.
  `hashLocation` is equal for equal Locations, cached per value.
- **Positions are int32**: `normalizeRanges` throws on an end that is not an int32 integer (a programmer error, never a
  wrap); the forms check data first and report `seq.query.position.invalid`.
- **Set operations** per system: union, intersection, subtraction. Rects intersect side by side over the same two
  systems; a rect is not subtracted, and its complement is refused (`seq.query.form.unsupported`).
- **Overlap** of an item `[start, end)`: a boundary (zero length) at `p` overlaps `[s, e)` when `s ≤ p < e`; one at the
  system's end (`p = length`) overlaps a range ending there, so `(all)` holds it. O(log n).

### Element keys

- `{ owner, key }`: the owner's canonical address and the data's own key, a string or a number. A number never equals a
  string (`1` is not `"1"`), as in an address. Elements are sorted (owner, numbers before strings) and each once;
  compared and hashed as values.

### Numberings and labels

- **Indices** are 0-based half-open inside; **labels** are inclusive names at the edge, converted once where a form is
  read.
- **Built-ins**: `index` (0-based) and `ordinal` (1-based), on every system; without `:numbering` a form reads
  `ordinal`. A declared numbering is looked up on the system, then on its molecule (`parent`, shifted by the unit's
  start, clipped to the unit): a unit's own beats its molecule's of the same id. Why: the author's numbering of a chain
  wins over the entity's. Not across a codon step: that waits for paths.
- **A label** is an integer with an optional suffix, letters (`94A`) or `.n` (`"100.1"`), or else a name. One total
  order: by integer; no suffix, then letters (shorter first, then alphabetic ignoring case), then `.n` by n; names after
  every number; spelling breaks the last ties. A `.n` label is written quoted: bare, `100.1` is the number 100.1 (the
  atoms' rule) and is refused with a hint to quote it; quoted, `"100.10"` and `"100.1"` are distinct labels. An integer
  outside int32 is refused (`seq.query.label.int32`). A repeated label names its first position.
- **Insertion codes**: `(at :label 94 :ins A)` is `94A`; `""`, `?`, `.` mean none; a code on a range, or on a label with
  a suffix, is `seq.query.label.insertion`.
- **Resolving a range**: each end is its label's position; a missing end snaps inward (a start to the first position
  whose label is at or after it, an end to the last at or before it, `seq.query.label.snapped`, info). Then one rule in
  index order for exact and snapped ends alike, never assuming labels rise with the index (IMGT): `from ≤ to` is the
  range; ends snapped past each other over positions with no label are empty (`seq.query.range.empty`); otherwise a
  cyclic system wraps as two ranges short of a whole lap (`seq.query.range.lap`), and any other is
  `seq.query.range.backwards`. A missing `:from` or `:to` is the system's own end.
- **Bounded label work**: a numbering may answer `snap` and `nearest` itself (`index` and `ordinal` do, in O(1));
  otherwise the positions are scanned, each a visit (signal, deadline). A hint lists the labels nearest an unknown one,
  kept a few at a time (no list of every label), and none past `HintScanLimit` positions.

### Forms

| Form                                           | Prefix            | Result                                                                                                          |
| ---------------------------------------------- | ----------------- | --------------------------------------------------------------------------------------------------------------- |
| `(range :from l :to l :numbering n)`           | systems           | Location                                                                                                        |
| `(at :label l :ins c :numbering n)`            | systems           | Location (one position)                                                                                         |
| `(all)`                                        | systems           | Location                                                                                                        |
| `(rect :a <location> :b <location>)`           | as its sides need | Location (one rect; each side one system; sides never wrap, `seq.query.range.wrap`)                             |
| `(item :id k)`, `(item :<field> v …)`          | annotations       | Elements (every param matches; a list param matches any of its atoms; a field the annotations lack is an issue) |
| `(items :overlapping <location>)`              | annotations       | Elements (every item without `:overlapping`)                                                                    |
| `(node :key k)`, `(chain :id k)`               | nodes             | Nodes, addressed by the form itself (`/data/1tup/(chain :id A)`), so location forms may follow a chain          |
| `(extent :of <elements>)`                      | any               | Location                                                                                                        |
| `(and :of […])`, `(or :of […])`, `(not :of x)` | any               | the operands' one type                                                                                          |

- **Combinators take `:of`**: a form's params are keywords only (the address grammar), so `and`, `or` and `not` read a
  list or one operand under `:of`. Operands of different types are `seq.query.form.type`.
- **`not`** complements within what the prefix names: its systems' positions, its annotations' items, its nodes'
  children with ids. Dynamic nodes are not in that universe yet.
- **Operands**: a nested form is evaluated against the same prefix (its resolver, else `repr`'s child selection); an
  address from the root. The location of `items :overlapping` is evaluated against the annotations' targets; a part on a
  system no target is reports `seq.query.location.unreachable` (projection is not built), a rect is refused. An operand
  naming a ref is refused (`seq.query.ref.unresolved`): refs are resolved by a later pass.
- **The table**: `queryResolvers(model, extra?)`; an extra kind named as a query form throws (a form's meaning is never
  replaced silently). `evaluateQuery` is `selectAddress` with the table, typed.

### The model the forms read

- An interface (`QueryModel`): `system(address)` (length, cyclic, numberings, parent), `annotation(address)` (target,
  fields, `items(over?)` with what is pending), `node(owner, key)` and `chain(owner, id)` for nodes data makes. The
  static model implements it; tests use a hand-written one.

### Staged and bounded

- **Staged**: over data read by region an evaluation gives `{ result, complete, pending }`, `pending` a Location of the
  regions not loaded; a form unions its operands' pending regions. Across two resolved segments of one address the
  selection keeps the last segment's (the union is `repr`'s, to come).
- **Bounded**: the caller's `limit` bounds the nodes and elements a form gives (more is
  `repr.tree.address.limit_exceeded`, never a truncation); the signal and the deadline are read every `CancelEvery`
  visits, one counter across nested forms and nested addresses (`repr.tree.address.cancelled`, `deadline_exceeded`). A
  Location's size is not bounded yet.

### Issues

- `seq.query.*` (`QueryIssues`): `form.prefix`, `form.param`, `form.type`, `form.unsupported`, `numbering.unknown`,
  `label.syntax`, `label.int32`, `label.unknown` (with the nearest labels), `label.insertion`, `label.snapped` (info),
  `range.backwards`, `range.lap`, `range.wrap`, `range.empty` (warning), `position.invalid`, `location.unreachable`
  (warning), `ref.unresolved`, `node.missing` (warning). A resolver's errors fail the selection; its warnings and infos
  are carried. A test checks the catalog against the codes the sources use.
