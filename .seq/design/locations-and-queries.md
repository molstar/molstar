# Locations, elements and queries

## Intent

The second review found the query language reduced to one paragraph, the meaning of address forms borrowed from
meso\* though Seq\*'s forms build locations over data (meso\*'s select child nodes), and the Location too thin for the
interaction grammar, which acts on items and rows, not only on ranges. This design settles the types and the first
forms.

## Settled

Accepted as recommended (the author, 2026-11-17, at the first slice's Gate 0): the three result types and labels 1-based
by default. The first forms grow with the examples.

## Decision

### Three result types

| Type | Holds | Made by | Used for |
|---|---|---|---|
| **Nodes** | exact addresses of spec nodes, data-made systems and held dynamic nodes | ids, kind forms (`(section :of …)`), `(node :key …)`, `(chain :id …)` | refs, layers' data, commands on tracks and views |
| **Elements** | keys of things data makes: items of an annotation, rows of an alignment, pairs | `(item …)`, `(row …)`, `(pair …)` | clicking, toggling, bookmarks on features, labels |
| **Location** | positions on systems: per system, sorted merged ranges; and rects over two axes | `(range …)`, `(at …)`, `(rect …)`, `(all)`, an element's extent | hover, selection, focus, highlights, layers' placement |

- **A Location** is `{ parts: { system, ranges: Int32Array }[], rects: { a: Part, b: Part }[] }`: ranges as sorted,
  merged, half-open index pairs on one system each; a cyclic system's wrapped range is two ranges (EC-9); rects pair
  two parts. Locations are values (immutable, comparable, hashable).
- **An element key** is `{ owner: address, key: string | number }` (an annotation's item id, a row's id, a chain's id);
  its extent is a Location, computed when asked.
- **Conversions are explicit**: `(extent <elements>)` gives a Location; `(items :overlapping <location>)` gives
  elements. A form that mixes types is an issue, not a guess.

### Addresses: meso\*'s syntax, a per-kind meaning

- **The syntax is meso\*'s**, extended with lists, strings and addresses as values ([addressing.md](addressing.md)).
- **The meaning is the form's**: `repr`'s `selectAddress` takes a table of form resolvers by form kind; meso\*'s child
  selection is the default resolver; Seq\*'s language registers its own (`range`, `at`, `rect`, `item`, `row`, `chain`,
  `node`, `all`, `extent`, `items`). A form resolver receives what the prefix names and returns one of the three types.
  Why: one grammar for addresses everywhere, each kind deciding what its forms mean.
- **Atoms** are disambiguated by the grammar: a token that parses as a JavaScript number is a number (`-1`, `1e5`);
  anything else is a symbol (`94A`, `HBA_HUMAN`); a string is quoted. `nil` is JSON's `null`.

### Numberings in forms

- **`:numbering <id>`** names a numbering of the system: one declared under its sequence (a `numbering` child: offset,
  segments, labels, a scheme), or a built-in: `index` (0-based positions) and `ordinal` (1-based). Without `:numbering`,
  labels are `ordinal` (1-based), stated once here (the JBrowse reply's guess, EC-195).
- **`.n` labels are quoted** (`(at :label "100.1")`): unquoted, `100.1` reads as a number.
- **A boundary** (a zero-length item) at p overlaps [s, e) when s ≤ p < e; a boundary at the system's end overlaps a
  range ending there, so `(all)` contains it.
- **Labels resolve at the edge** with the rules EC-1 to EC-8 (order, insertion codes, a missing end snapping inward, an
  unknown single label reported with the labels that exist).

### The first forms

| Form | Result | Example |
|---|---|---|
| `(range :from a :to b :numbering n)` | Location | `/data/p53/(range :from 94 :to 292)` |
| `(at :label l :numbering n)` | Location | `/data/1tup/(chain :id A)/(at :label 94A :numbering auth)` |
| `(rect :a <location> :b <location>)` | Location | a dot plot's rectangle |
| `(all)` | Location | a whole system |
| `(item :id k)`, `(item :<field> v)` | Elements | `/data/p53/domains/(item :label DNA-binding)` |
| `(chain :id k)` | Nodes (a system) | `/data/1tup/(chain :id A)` |
| `(row :id k)` | Elements | with alignments |
| `(node :key k)` | Nodes | `/data/genes/(node :key T1)` |
| `(extent :of <elements>)`, `(items :overlapping <location>)` | Location, Elements | conversions |
| `(and :of [ … ])`, `(or :of [ … ])`, `(not :of …)` | the type of their arguments | combinators over one type (keyword params: the address grammar has no positional ones) |

Later, with the examples that need them: `expand`, `onto` (projection as a query), `search` (motifs), `group`.

### Evaluation over data read by region

- **Staged**: a query over a source read by region evaluates over the regions loaded and says whether it is complete
  (`{ result, complete, pending: regions }`); a layer's data query is evaluated per frame over the regions in view, and
  a pending region is drawn as pending, never as empty. A query that needs the whole data (a motif over a genome) asks
  the source to evaluate it, when its runtime can; otherwise it is refused with an issue (later, with `search`).
- **Bounded**: every evaluation takes a limit, a depth and a deadline (meso\*'s bounds); more matches is an issue,
  never a truncation.

## Consequences

- Hover, selection, focus and highlights hold Locations and element keys ([viewer-models.md](viewer-models.md)).
- `repr` gains the form-resolver table in the tree core's first move; its default keeps meso\*'s meaning.
- Tests: every form's result type and examples; atoms; numbering cases EC-1 to EC-8; conversions; staged evaluation
  with a pending region.
