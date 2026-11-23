# Locations, keys and queries

## Intent

The second review found the query language reduced to one paragraph, the meaning of address forms borrowed from
meso\* though Seq\*'s forms build locations over data (meso\*'s select child nodes), and the Location too thin for the
interaction grammar, which acts on items and rows, not only on ranges. This design settles the types and the first
forms.

## Settled

Accepted as recommended (the author, 2026-11-17, at the first slice's Gate 0): the three result types and labels 1-based
by default. The first forms grow with the examples.

Decided by the author, 2026-11-20, for the structure forms ([The language](#the-language-forms-by-role)): "we should
not really have (atom), we should have (elements ....) and that gets "schema selection" etc .. similar to "atom
groups" combined with "schema item" in molQL .... could still have two generators"; on the two generators, "Accept
both"; on the result type *Elements*, "Rename to Keys"; on a Location as a prefix, "ok, yes, it can take a common
prefix"; "we should have a generic filter, and generic properties, but it should be easy to add own filter,
modified/combinator etc. ... you are mixing two things: property access and language design" (property access is
[properties.md](properties.md)); on the draft, "Accept, but allow positional operands"; on the combinators, "Rename".

## Decision

### Three result types

| Type | Holds | Made by | Used for |
|---|---|---|---|
| **Nodes** | exact addresses of spec nodes, data-made systems and held dynamic nodes | ids, kind forms (`(section :of …)`), `(node :key …)`, `(chain :id …)` | refs, layers' data, commands on tracks and views |
| **Keys** | keys of things data makes: items of an annotation, rows of an alignment, pairs, groups of a structure | `(item …)`, `(groups …)` (`(row …)` is a Location on the rows system, pairs are items) | clicking, toggling, bookmarks on features, labels |
| **Location** | positions on systems: per system, sorted merged ranges; and rects over two axes | `(range …)`, `(at …)`, `(rect …)`, `(all)`, an element's extent | hover, selection, focus, highlights, layers' placement |

- **A Location** is `{ parts: { system, ranges: Int32Array }[], rects: { a: Part, b: Part }[] }`: ranges as sorted,
  merged, half-open index pairs on one system each; a cyclic system's wrapped range is two ranges (EC-9); rects pair
  two parts. Locations are values (immutable, comparable, hashable).
- **A key** (`Key`; the type was *Elements*, renamed by the author, 2026-11-20) is `{ owner: address, key: string | number }` (an annotation's item id, a row's id, a chain's id);
  its extent is a Location, computed when asked.
- **Conversions are explicit**: `(extent <keys>)` gives a Location; `(items :overlapping <location>)` gives
  keys. A form that mixes types is an issue, not a guess.

### Addresses: meso\*'s syntax, a per-kind meaning

- **Pairs answer `items`** (decided, the author, 2026-11-19, for `linkedPairs`): a pairs annotation is queried as its
  items, a pair one item with a row per end, each row with its own system (`QueryItem.system`): `(items :overlapping
  <location>)` gives the pairs with an end in the location, on either system.
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

Forms compose monadically, as MolQL's: a form maps over every element its prefix yields and concatenates what each
gives; an element with nothing to give yields nothing, never a failure of the whole query (the author, 2026-11-20:
"it needs to work similar to how MolQL does. monadic composition.").

| Form | Result | Example |
|---|---|---|
| `(range :from a :to b :numbering n)` | Location | `/data/p53/(range :from 94 :to 292)` |
| `(at :label l :numbering n)` | Location | `/data/1tup/(chain :id A)/(at :label 94A :numbering auth)` |
| `(rect :a <location> :b <location>)` | Location | a dot plot's rectangle |
| `(all)` | Location | a whole system |
| `(item :id k)`, `(item :<field> v)` | Keys | `/data/p53/domains/(item :label DNA-binding)` |
| `(chain :id k)` | Nodes (a system) | `/data/1tup/(chain :id A)` |
| `(row :key k)`, or a list of keys | Location (on an alignment's rows system) | with alignments (amended 2026-11-23: the code's way, the author as recommended; no Keys row form, no `(pair …)`: pairs answer `items`) |
| `(node :key k)` | Nodes | `/data/genes/(node :key T1)` |
| `(extent :of <keys>)`, `(items :overlapping <location>)` | Location, Keys | conversions |
| `(merge :of [ … ])`, `(intersect :of [ … ])`, `(except :of [ … ])` | the type of their arguments | combinators over one type (renamed from `or`, `and`, `not`; `except` with one operand the complement within the prefix) |
| `(elements …)`, `(groups …)`, `(filter :test …)`, `(expand …)`, `(linked …)` | Location, Keys, the prefix's, Location, Location | [The language](#the-language-forms-by-role) |

Later, with the examples that need them: `onto` (projection as a query), `search` (motifs).

### The language: forms by role

MolQL's split (`struct.generator`, `struct.filter`, `struct.modifier`, `struct.combinator`, `core.*` operators), with
property access apart ([properties.md](properties.md)). Every form is monadic over a **common prefix**: nodes, a
Location, or Keys (read as their extents).

| Role | Built in | Result | MolQL |
|---|---|---|---|
| generator | `(elements <items>)` | Location | `elementSchema` |
| generator | `(groups :molecule-test p :unit-test p :residue-test p :atom-test p :by level)` | Keys | `atomGroups` |
| filter | `(filter :test p)` | the prefix's type, a subset | `pick`, per-element tests |
| modifier | `(expand :to residue \| unit \| molecule)` | Location | `wholeResidues` |
| modifier | `(linked :via <pairs> :type w)` | Location: the prefix and what links or a host's contacts join | `includeConnected` |
| combinator | `(merge :of [ … ])`, `(intersect :of [ … ])`, `(except :of [a b …])` | the operands' one type | `merge`, `intersect`, `exceptBy` |
| operator | `(= a b)`, `!=`, `<`, `<=`, `>`, `>=`, `(in-range x lo hi)`, `(in x [ … ])`, `(and p …)`, `(or p …)`, `(not p)` | a value | `core.rel`, `core.set.has`, `core.logic` |

- **`elements`** is the schema selection: a table of items ([addressing.md](addressing.md), Schema selections), each
  item a record of properties; its finest property's level decides the result: no residue or atom property, every
  matching unit's whole residue system; residue properties, those positions; an atom property, positions on the
  units' atoms systems (`<unit>/(atoms)`). Why: the record says how deep it goes, as MolViewSpec's does.
- **`elements` is sugar over `filter`, exactly**: a field `:f v` is `(= (prop f) v)`; `:beg b :end e` is
  `(in-range (prop seq) b e)` (`auth-beg`/`auth-end` over `auth-seq`, in the numbering's order, EC-1 to EC-8); a row is
  `and` of its fields; items are `or` of their rows; `(elements I)` is `(filter :test P(I))` over the prefix's elements
  at I's finest level. Only speed differs: rows are grouped by their unit fields and each unit's index walked once
  (MolQL's `elementSchema`): rows plus matches, not rows × units. Why one model: the record a structure viewer sends and
  a hand-written predicate mean the same, and a property test holds them equal.
- **`groups`**: each test reads its level's properties (a lower level's is `seq.query.property.level`); `:by` is
  `atom`, `residue` (default), `unit` or `molecule`; one key a group (owner the unit's or atoms system's address and a
  0-based position, or a complex's address and a unit or entity id). Why a second generator: clicking, labels and
  bookmarks act on keys.
- **`filter`** tests each element of its prefix: a position of a Location (a residue or an atom, by its system), a key
  of Keys (its level's element). One generic filter: an extension adds a property, not a filter, in most cases.
- **`expand`** uses the atoms-to-residues map and the unit and molecule steps; **`linked`** follows a `pairs` node
  without projecting through it ([molecules.md](molecules.md): links are never maps). Each has a use now (an atom
  selection to its residues; the glycan on an Asn). Not built: `include-surroundings` (no coordinates; a host's
  contacts are `pairs`, reached by `linked`).
- **Combinators renamed** (the author, 2026-11-20, "Rename"): `or`, `and`, `not` on results became `merge`,
  `intersect`, `except` (MolQL's; `except` with one operand is the complement within the prefix), and `and`, `or`,
  `not` are predicate logic. Why: one name one meaning; 0.5.0 is unreleased (no migration).
- **Positional operands in expressions only** (the author, 2026-11-20, "Accept, but allow positional operands"): the
  grammar reads positional values before keyword params in any form ([addressing.md](addressing.md)); only an
  expression form (an `operator`, `prop`, `value`) accepts them, and any other form with one is
  `seq.query.form.param`. Why the line there: a predicate nests deep and reads as MolQL's; a form over results names
  what it takes.
- **Types are checked when the query is read**: an operator over a word and a number, a test that is not boolean, a
  combinator over two types: `seq.query.form.type`.
- **Extensions add forms** as `query-form` registry entries ([extension-openings.md](extension-openings.md)):
  `{ kind: 'query-form', name: '<id>:<name>', role, params, result, resolve }`, `role` one of the five above, `resolve`
  the built-ins' signature (an operator's over values). Checked at composition (namespace, no built-in shadowed,
  conflicts across the descriptor); by role at evaluation (a filter's result a subset of its prefix, of its type; a
  combinator's operands one type; an operator's result a value); bounded through `ctx`. Printed by the grammar, params
  sorted by key, as every form is (one canonical spelling, so an address's text does not depend on whether the reader has
  the extension) (amended 2026-11-23: the code's way, the author as recommended). A namespaced form nobody registered is `seq.query.extension.missing`, naming the namespace; no
  fallback (a query cannot be approximated). A bare unknown kind keeps meso\*'s meaning (child selection by kind).
- **Replaced** (ML2's forms): `(atom :name O1G)` by `(elements :atom O1G)`; `(component :id MSE)` by
  `(elements :component MSE)`, its `:standard`, `:class`, `:annotations` by `filter` over the properties `modified`,
  `standard`, `component-class`; `(entity :id 3)` by `(elements :entity 3)` and the node `(molecule :entity 3)`.
  `(chain …)`, `(molecule …)`, `(atoms)` stay: Nodes, what `(section :of …)` takes.

Examples: `/data/4hhb/(elements :component HEM :atom FE)` (four irons); `/data/1fin/(elements :chain A :beg 10 :end 20
:atom CA)` (eleven Cα atoms); `/data/1fin/(elements :chain A :atom CA)/(filter :test (> (prop b-factor) 60))`;
`/data/af-p53/(elements)/(filter :test (>= (prop af:plddt) 70))`; `/data/glycan/(elements :atom C1)` over three
residues (three atoms); `/data/1ehz/(groups :residue-test (= (prop modified) true))` (a key a modified residue).

### Evaluation over data read by region

- **Staged**: a query over a source read by region evaluates over the regions loaded and says whether it is complete
  (`{ result, complete, pending: regions }`); a layer's data query is evaluated per frame over the regions in view, and
  a pending region is drawn as pending, never as empty. A query that needs the whole data (a motif over a genome) asks
  the source to evaluate it, when its runtime can; otherwise it is refused with an issue (later, with `search`).
- **Bounded**: every evaluation takes a limit, a depth and a deadline (meso\*'s bounds); more matches is an issue,
  never a truncation.

## Consequences

- Hover, selection, focus and highlights hold Locations and keys ([viewer-models.md](viewer-models.md)).
- `repr` gains the form-resolver table in the tree core's first move; its default keeps meso\*'s meaning.
- Tests: every form's result type and examples; `elements` equal to its `filter` (a property); SoA equal to AoS; atoms; numbering cases EC-1 to EC-8; conversions; staged evaluation
  with a pending region.

## In the specs

| Decision | Where |
|---|---|
| Three result types: Nodes, Keys, Location | core › Result types |
| A Location's shape; a value, comparable and hashable | core › Location |
| A key `{ owner, key }`, its extent | core › Keys |
| Conversions are forms, explicit; mixed types an issue | core › Conversions are forms |
| Pairs answer `items` | core › Pairs answer items |
| meso\*'s syntax with lists, strings, addresses | repr › Grammar |
| The meaning is the form's: a resolver table, child selection the default | repr › The form-resolver table |
| Seq\*'s resolvers registered as a table | core › The table |
| Atoms: a number if it parses as one, else a symbol | repr › Atoms' spelling |
| Labels 1-based by default (Settled) | core › Built-ins |
| `:numbering`, the built-ins `index` and `ordinal`, `ordinal` by default (EC-195) | core › Built-ins |
| `.n` labels quoted | core › A label |
| A boundary's overlap | core › Overlap |
| Labels resolve at the edge (EC-1 to EC-8) | core › Resolving a range |
| Monadic composition | core › Composition |
| The first forms | core › Forms |
| `(row :key k)` giving a Location on the rows system, no `(pair …)` form (pairs answer `items`) | core › Forms |
| `onto` and `search` | design-only: not built (later, with their examples) |
| The forms by role, MolQL's split, a common prefix | core › Language |
| `elements`: the schema selection, its depth its finest field | core › elements |
| `elements` sugar over `filter`, exactly; rows plus matches | core › Sugar, exactly |
| `groups` | core › groups |
| `filter` | core › filter |
| `expand` and `linked`; `include-surroundings` not built | core › expand |
| Combinators renamed to `merge`, `intersect`, `except` | core › Combinators take :of |
| Positional operands in expressions only | core › Expressions take positional operands |
| Types checked when the query is read | core › Types are checked when a test is compiled |
| Extensions add forms: namespaced, checked at composition and by role, no fallback | core › Composed |
| An extension form's params printed sorted by the grammar, like every form | core › Extensions |
| Replaced forms (`atom`, `component`, `entity`) | core › Molecule node forms |
| Staged evaluation over data read by region | core › Staged |
| A query needing the whole data asks the source, else refused | design-only: not built (with `search`) |
| Bounded: a limit, a depth, a deadline | core › Staged and bounded |
| Tests | core › Examples |
