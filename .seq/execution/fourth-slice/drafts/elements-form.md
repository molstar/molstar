# The `(elements …)` form: Seq\*'s structure query language and property access

## Intent

The author, 2026-11-20: "we should not really have (atom), we should have (elements ....) and that gets "schema
selection" etc .. similar to "atom groups" combined with "schema item" in molQL .... could still have two generators."
Forms compose monadically, as MolQL's (`(atom :name C1)` over a glycan's three residues gives three atoms).

ML2's branch built one form a level (`entity`, `chain`, `component`, `atom`, `atoms`, `molecule`). This design replaces
them with MolQL's split: **a language** of forms by role (generators, filters, modifiers, combinators, and the
operators a predicate uses), and, apart from it, **property access**: typed symbols read per element. The schema
selection a structure viewer speaks is sugar over the two. Each part is open to extensions through the registry.

## Decided by the author (2026-11-20)

1. Generators: "Accept both": `elements` and `groups`; `atom`, `entity`, `component` removed.
2. Field names: "support both short and long aliases"; "also make the language extensible so people could add
   generators, filters, etc on their own props (basically just adding a few symbols)".
3. The result type *Elements*: "Rename to Keys".
4. "ok, yes, it can take a common prefix. also remember the SoA and AoS interchangeable representations we have."
5. "we should have a generic filter, and generic properties, but it should be easy to add own filter,
   modified/combinator etc. likely won't need it many cases, but in general. you are mixing two things: property access
   and language design."
6. Records in an address as `:items [(fields …) …]`: "Yes".
7. The draft: "Accept, but allow positional operands": inside a predicate, operators take positional arguments
   (`(filter :test (> (prop b-factor) 60))`), an exception to the keyword-only rule for expressions only.
8. The combinators: "Rename": `merge`, `intersect`, `except` on results; `and`, `or`, `not` in predicates.

## Precedents

- **MolQL** (`packages/query-language`, `structure-query.ts`, `core.ts`): `struct.generator` (`atomGroups` with
  `entity-`, `chain-`, `residue-`, `atom-test` and `group-by`; `elementSchema` over `schemaItem` records),
  `struct.filter` (`pick :test`, `within`, `isConnectedTo`), `struct.modifier` (`wholeResidues`, `includeSurroundings`,
  `includeConnected`, `exceptBy`, `intersectBy`, `unionBy`), `struct.combinator` (`merge`, `intersect`), `core.rel`
  (`eq`, `lt`, `inRange`), `core.logic`, `core.set.has`; properties apart, as `struct.atomProperty.*` (`macromolecular`:
  `label_asym_id`, `auth_seq_id`, `label_comp_id`, `B_iso_or_equiv`, `occupancy`, …) read for the current element.
  Taken: the roles, the split, the names `merge`/`intersect`, the schema record as a table. Not taken: geometry
  (`within`, surroundings: a host's contacts instead), math beyond comparison (no use yet).
- **MolViewSpec `ComponentExpression`**: the round-trip record; its names are the long aliases.
- **Records files** ([files.md](../../../design/files.md)): AoS and SoA read the same. Taken for the selection's items.
- **Extension descriptor** ([extension-openings.md](../../../design/extension-openings.md)): entries `<id>:<name>`,
  registered whole. Taken: two entry kinds, `query-form` and `query-property`.
- **Interactions**: a typed or pasted selection (the query command, built), a link or bookmark (built), a structure
  viewer's selection coming in (ST2). No new gestures; the matrix's *selection* row gains these forms.

## Part 1: the language

### Forms by role

Forms take keyword params ([addressing.md](../../../design/addressing.md)); a predicate is a param (`:test`). **Inside a
predicate, operators take positional operands** (`(> (prop b-factor) 60)`), as MolQL's: the grammar reads positional
values before keyword params in any form, and only an expression form (an `operator`, `prop`, `value`) accepts them;
any other form with one is `seq.query.form.param`. Every form is
monadic over a **common prefix**: nodes, a Location, or Keys (read as their extents); each element its prefix yields
gives a part, one with nothing gives nothing, parts union.

| Role | Built in | Result | MolQL |
|---|---|---|---|
| generator | `(elements <items>)` | Location | `elementSchema` |
| generator | `(groups :molecule-test p :unit-test p :residue-test p :atom-test p :by level)` | Keys | `atomGroups` |
| filter | `(filter :test p)` | the prefix's type, a subset | `pick`, the per-element tests |
| modifier | `(expand :to residue \| unit \| molecule)` | Location | `wholeResidues` |
| modifier | `(linked :via <pairs> :type w)` | Location: the prefix plus what the links (or a host's contacts) join | `includeConnected` |
| combinator | `(merge :of [ … ])`, `(intersect :of [ … ])`, `(except :of [a b])` | the operands' one type | `merge`, `intersect`, `exceptBy` |
| operator | `(= a b)`, `!=`, `<`, `<=`, `>`, `>=`, `(in-range x lo hi)`, `(in x [ … ])`, `(and p …)`, `(or p …)`, `(not p)` | a value | `core.rel`, `core.set.has`, `core.logic` |

- **The finest level decides `elements`**: no residue or atom property among its fields: every matching unit's whole
  residue system; residue properties: those positions; an atom property: positions on the units' atoms systems
  (`<unit>/(atoms)`). Why: the record says how deep it goes, as MVS's does.
- **`groups`**: each test reads properties of its level (a lower level's property is `seq.query.property.level`); `:by`
  is `atom`, `residue` (default), `unit` or `molecule`; one key a group (owner the unit's or atoms system's address and a
  0-based position, or a complex's address and a unit or entity id), its extent its Location. Why a second generator:
  clicking, labels and bookmarks act on keys.
- **`filter`** evaluates `:test` per element of its prefix: a position of a Location (a residue or an atom, by its
  system), a key of Keys (its level's element). Why one generic filter: predicates over properties cover what the
  removed forms' params did (`:standard`, `:class`), and what extensions add.
- **`expand`** uses the atoms-to-residues map and the unit and molecule steps; **`linked`** follows a `pairs` node
  (links, a host's contacts) without projecting through it (molecules.md: links are never maps). Both have a use now (an
  atom selection to its residues; the glycan on an Asn, a payload on its Cys). `include-surroundings` is not built: no
  coordinates; a host's contacts are `pairs`, reached by `linked`.
- **Combinators renamed**: today's result combinators `and`, `or`, `not` become `intersect`, `merge`, `except`
  (MolQL's; `except` with one operand is today's `not`, the complement within the prefix), freeing `and`, `or`, `not`
  for predicate logic. Why: one name one meaning; 0.5.0 is unreleased.
- **Types are checked when the query is read**: an operator over a word and a number, a filter whose test is not
  boolean, a combinator over two types: `seq.query.form.type`.

### Items: the schema selection, AoS and SoA

- **One item** is `elements`' params (`(elements :chain A :seq 33)`); **columns** are lists, row *i* across them one
  item, a scalar the same in every row (`(elements :chain [A C] :auth-seq 87)`: A 87 and C 87); **records** are
  `:items [(fields :chain A :seq 33) (fields :chain B :beg 140 :end 146)]` (`fields` is a record, never evaluated). In
  JSON and the builder: an object (columns where its values are arrays) or an array of objects, as a records file.
  `nil` in a row is no test on that field.
- **Normalized to SoA** (typed columns; a field absent in every row dropped), **printed in the shorter layout** (one
  item as params; else columns or records, ties to columns), so `parse(print(a))` equals `a` as a value.
- **Bounds**: columns of unequal lengths, a `fields` that is not a record, `seq` with `beg`/`end` in one row, `beg` past
  `end`: form errors naming the row. Rows count against the evaluation's `limit`; more is an issue, never a truncation.
- **Sugar, exactly**: a field `:f v` is `(= (prop f) v)`; `:beg b :end e` is `(in-range (prop seq) b e)` (`auth-beg`/`auth-end` over `auth-seq`, through the numbering's order, EC-1 to EC-8); a row is `and` of
  its fields; items are `or` of their rows. So `(elements I)` is `(filter :test P(I))` over the prefix's elements at
  the finest level of I. Evaluation differs only in speed: it groups rows by their unit fields and walks each unit's
  index once (MolQL's `elementSchema`), so the cost is rows plus matches; a property test checks the walk equals the
  filter.

### Extensions add forms

```ts
| { kind: 'query-form'; name: string;   // `af:confident`
    role: 'generator' | 'filter' | 'modifier' | 'combinator' | 'operator';
    params: ParamDefs;                  // keyword params; `:of` for operands
    result: 'location' | 'keys' | 'nodes' | 'value';
    resolve: Resolve }                  // the built-ins' signature; an operator's is over values
```

- **Checked at composition** as every extension entry: names `<id>:<name>` (`seq.context.extension.namespace`), no
  built-in shadowed, conflicts across the descriptor whole. **Checked by role** at evaluation: a filter's result is a
  subset of its prefix, of its type; a combinator's operands one type; an operator's result a value. Bounds (`tick`,
  `limit`, deadline) come through `ctx`. `queryResolvers` becomes the built-ins plus the registry's entries.
- **Printing**: by the grammar (`(af:confident :min 70)`), params in declared order.
- **An unknown symbol**: a namespaced form nobody registered is an error, `seq.query.extension.missing`, naming the
  namespace; no fallback (a query cannot be approximated, as a format runtime has none). A bare unknown kind keeps
  today's meaning (child selection by kind).

## Part 2: property access

- **A property is a typed symbol read per element** of one level: `(prop chain)`, in a predicate or as an
  `elements` field. Types: `word`, `integer`, `number`, `boolean`, `words` (a set: `component` holds a position's
  alternatives; `=` on a set is membership).
- **Built in**, short name (printed) and long alias:

| Level | Properties |
|---|---|
| molecule | `entity` (`label_entity_id`), `molecule-class`, `identifier` (words) |
| unit | `chain` (`label_asym_id`), `auth-chain` (`auth_asym_id`), `operator` (`instance_id`) |
| residue | `seq` (`label_seq_id`: 1-based on the unit's molecule), `auth-seq` (`auth_seq_id`), `ins` (`pdbx_PDB_ins_code`), `index` (`residue_index`), `component` (`label_comp_id`, words), `component-class`, `standard`, `modified` (non-standard, or under a `modified`-type annotation item: `component`'s old `:annotations`) |
| atom | `atom` (`label_atom_id`), `element` (`type_symbol`), `atom-index` (`atom_index`), `b-factor` (`B_iso_or_equiv`), `occupancy` |
| any | `(value <values annotation> :field f)`: a values channel on the element's system |

  Both names read; the short one prints. Why: addresses are typed and read beside a sequence, and one printed form keeps
  equal selections equal as text; the long names are a structure viewer's, so inward translation is a copy. No
  `label_alt_id`, `auth_comp_id`, `auth_atom_id`, `atom_id`: Seq\* keeps one atom a name and no author names.
- **How a property is read**: a column per element on one system, `{ length, get(i) }` (typed where it can be), from a
  `QueryModel` hook: `property?(system, name): PropertyRead`, `PropertyRead = column | 'pending' | undefined`. Units,
  components, numberings and atoms answer the built-ins (`QueryUnit.entity` replaces ML2's key-of-`of` lookup;
  `QueryAtoms.element`, `residue`); a values channel reads the store (`values` annotations by region): `pending` stages
  the system, as a numbering not arrived; `b-factor` and `occupancy` are values channels the mmCIF runtime writes on a
  unit's atoms, lacking otherwise. A residue property on an atom reads its residue's (through the atoms map); the
  reverse is `seq.query.property.level`.
- **An element that lacks a property** (`chain` on a lone molecule, `auth-seq` with no `auth` numbering, `atom` without
  atoms, `entity` on an inline complex): its test is false and it yields nothing, with one warning
  `seq.query.property.absent` naming the property and the count. Monadic: never the whole query's failure.
- **Extensions add properties**, read by every generator test and filter, and usable as `elements` fields:

```ts
| { kind: 'query-property'; name: string;  // `af:plddt`
    level: 'molecule' | 'unit' | 'residue' | 'atom'; type: PropertyType; aliases?: readonly string[];
    read: ((model: QueryModel, system: string) => PropertyRead)
        | { values: string; field: string } }  // "a few symbols": a values annotation's field by its kind or path
```

  Namespaced and checked as forms; an unknown namespaced property is `seq.query.extension.missing`, an unknown bare one
  `seq.query.form.param` with the names in its hint (and the short name for a long one misspelled).
- **Outward**: no property reaches a structure viewer; the adapter translates the Location selected.

## Examples

- **4HHB**: `/data/4hhb/(elements :component HEM :atom FE)`: four irons; α's proximal His: `(elements :auth-chain [A C]
  :auth-seq 87)`.
- **1FIN**: `/data/1fin/(elements :chain A :beg 10 :end 20 :atom CA)`: eleven Cα atoms of the glycine loop;
  `/data/1fin/(elements :chain A :atom CA)/(filter :test (> (prop b-factor) 60))`: its mobile Cα atoms.
- **A pLDDT filter**: `/data/af-p53/(elements)/(filter :test (>= (prop af:plddt) 70))`, `af:plddt` an
  extension's property over a values annotation, or without one `(value /data/af-p53/plddt :field score)`.
- **A glycan** `[NAG NAG BMA]`: `/data/glycan/(elements :atom C1)`: three atoms; `<selection>/(elements :atom C1)`: the
  C1s of the residues it holds; `/data/1fin/(elements :chain A :seq 33)/(linked :via /data/1fin/links)`.
- **Modified residues**: `/data/1ehz/(groups :residue-test (= (prop modified) true))`: a key each (11);
  `(elements :component MSE :element SE)`: seleniums.

| ML2's form | Becomes |
|---|---|
| `(atom :name O1G)` | `(elements :atom O1G)` |
| `(component :id MSE)`, `:standard false`, `:class saccharide` | `(elements :component MSE)`, `(filter :test (= (prop modified) true))`, `(filter :test (= (prop component-class) saccharide))` |
| `(entity :id 3)` | `(elements :entity 3)`; the node: `(molecule :entity 3)` |
| `(chain …)`, `(molecule …)`, `(atoms)` | stay: Nodes, what `(section :of …)` takes |

## Configuration

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| property names | fixed (0.5.0) | short printed, long read | readable, one canonical text | — |
| item layouts | fixed | AoS and SoA read, SoA held, shorter printed | as records files | — |
| `groups :by` | form param | `residue` | what a user clicks and labels | `atom`, `unit`, `molecule` |
| an element lacking a property | fixed | false, one warning | monadic | a form error |
| forms and properties | registry entries | the built-ins | extensions add symbols | — |

## Tensions

- **Two ways to a range**: `(range :numbering auth …)` (any numbering, any system) and `elements`' `auth-beg`/`auth-end`
  (the record): kept, one numbering code under both.
- **Self-check** ("the language and properties cover what a structure user selects"): negated, a user selects by
  geometry (a pocket, a radius) or an altloc. Integrated: contacts a host computes, through `linked`; per-residue
  scores, through properties and `filter`. Standing: radius and altlocs are the structure viewer's; such a selection
  comes in as its residues.

## Consequences

- **`seqx-core`**: `query/` gains `elements`, `groups`, `filter`, `expand`, `linked`, the operators, `prop`, `value`, the
  items' normalization and a property table; loses `atom`, `entity`, `component`; renames `and`/`or`/`not`
  (combinators) and the result type `Keys` (`ElementKey` → `Key`); `QueryModel.property`. `context/`: `query-form`,
  `query-property`, their checks.
- **Tests**: each property and alias on the 1FIN and 4HHB fixtures; the examples; AoS, SoA and records equal; the sugar
  equal to its filter (a property test); the glycan's C1s; a lacking property; a pending values channel; an extension's
  property and filter, registered and missing; an MVS expression inward then outward through the fake's adapter.
- **Docs**: `_spec/README.md` "Queries" (Language, Properties, Keys, Extensions) and "Context"; molecules.md's Queries,
  locations-and-queries.md, extension-openings.md (what an extension may ship). ML2's branch takes the change before it
  merges. **Migration**: none; 0.5.0 is unreleased.

## Gaps for the author

All answered (2026-11-20): two generators, aliases, Keys, a common prefix, SoA/AoS, `:items [(fields …) …]`, a generic
filter and properties apart, positional operands in predicates, the combinators renamed. Agreed; folded into
locations-and-queries.md, properties.md, addressing.md, molecules.md, extension-openings.md and viewer-composition.md.
